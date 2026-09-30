! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! optional outputs requested on every call to the mock density matrix evaluating
    ! functions since the last reset, as bits in the order of their argument lists, so
    ! that tests can check that a routine only asks for what it needs; their number
    ! lets the mocks produce potential matrices which change between calls and thereby
    ! distinguish a cached quantity from one that was incorrectly recomputed from an
    ! unchanged input
    integer(ip), allocatable :: mock_requests(:)

    ! multiplier of the density matrix returned by the mock density matrix evaluating
    ! functions, which differs between the first and any subsequent call, further
    ! multipliers for extension-specific potentials follow the same convention
    real(rp), parameter :: mock_fock_factor(2) = [2.0_rp, 5.0_rp]

    ! names of the closed-shell and the open-shell case, indexed by the number of
    ! particle channels, for the failure messages of tests covering both
    character(len=12), parameter :: shell_names(2) = &
        [character(len=12) :: "closed-shell", "open-shell"]

contains

    function mock_factor(factors) result(factor)
        !
        ! this function returns the multiplier the mock density matrix evaluating
        ! functions apply on the current call
        !
        real(rp), intent(in) :: factors(:)
        real(rp) :: factor

        factor = factors(min(size(mock_requests), size(factors)))

    end function mock_factor

    subroutine record_mock_call(request)
        !
        ! this subroutine records a call to a mock density matrix evaluating function
        ! together with the optional outputs it was asked for
        !
        integer(ip), intent(in) :: request

        if (allocated(mock_requests)) then
            mock_requests = [mock_requests, request]
        else
            mock_requests = [request]
        end if

    end subroutine record_mock_call

    function identity_matrix(n) result(matrix)
        !
        ! this function returns the identity matrix
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        integer(ip) :: i

        matrix = 0.0_rp
        do i = 1, n
            matrix(i, i) = 1.0_rp
        end do

    end function identity_matrix

    function generate_random_symm_matrix(n) result(matrix)
        !
        ! this function generates a random symmetric matrix, corresponding to a valid
        ! density matrix displacement
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        call random_number(matrix)
        matrix = matrix + transpose(matrix)

    end function generate_random_symm_matrix

    function generate_random_orthogonal_matrix(n) result(matrix)
        !
        ! this function returns a random orthogonal matrix, obtained as the
        ! eigenvectors of a random symmetric matrix
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        integer(ip) :: lwork, info
        real(rp) :: eigvals(n)
        real(rp), allocatable :: work(:)
        external :: dsyev

        matrix = generate_random_symm_matrix(n)
        allocate(work(1))
        call dsyev("V", "U", n, matrix, n, eigvals, work, -1_ip, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        call dsyev("V", "U", n, matrix, n, eigvals, work, lwork, info)
        deallocate(work)

    end function generate_random_orthogonal_matrix

    function generate_random_density_matrix(n, n_occ) result(dm)
        !
        ! this function generates a random valid density matrix by first generating a
        ! random set of orthonormal basis vectors and then summing the outer products
        ! of these
        !
        use opentrustregion, only: numerical_zero

        integer(ip), intent(in) :: n, n_occ
        real(rp) :: dm(n, n)

        real(rp) :: vec(n), val, u(n, n_occ)
        integer(ip) :: i, j

        ! generate random orthonormal basis
        do j = 1, n_occ
            do
                call random_number(vec)
                do i = 1, j - 1
                    vec = vec - sum(vec * u(:, i)) * u(:, i)
                end do
                val = sqrt(sum(vec**2))
                if (val >= numerical_zero) then
                    u(:, j) = vec / val
                    exit
                end if
            end do
        end do

        ! construct valid density matrix by summing outer products of basis vectors
        dm = matmul(u, transpose(u))

    end function generate_random_density_matrix

    logical(c_bool) function test_matrix_exponential() bind(C)
        !
        ! this function tests the function which calculates the matrix exponential of a
        ! real antisymmetric matrix
        !
        use otr_common, only: matrix_exponential
        use opentrustregion, only: solver_settings_type
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n = 2
        real(rp), parameter :: angle = 0.3_rp

        real(rp) :: a(n, n), expected(n, n)
        real(rp), allocatable :: exp_a(:, :)
        type(solver_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_matrix_exponential = .true.

        ! set up the settings the routine logs through
        call setup_settings(settings)

        ! initialize antisymmetric matrix, whose exponential is the corresponding
        ! rotation matrix
        a = reshape([0.0_rp, -angle, angle, 0.0_rp], [n, n])

        ! initialize expected rotation matrix
        expected = reshape([cos(angle), -sin(angle), sin(angle), cos(angle)], [n, n])

        ! call routine and determine if dimensions and values of resulting matrix match
        exp_a = matrix_exponential(a, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_matrix_exponential failed: Produced error."
            test_matrix_exponential = .false.
        end if
        if (size(exp_a, 1) /= n .or. size(exp_a, 2) /= n) then
            write(stderr, *) "test_matrix_exponential failed: Incorrect matrix "// &
                "dimensions."
            test_matrix_exponential = .false.
            return
        end if
        if (norm2(exp_a - expected) > tol) then
            write(stderr, *) "test_matrix_exponential failed: Incorrect matrix values."
            test_matrix_exponential = .false.
        end if
        deallocate(exp_a)

        ! call routine for a vanishing matrix and determine if the identity matrix is
        ! returned
        a = 0.0_rp
        exp_a = matrix_exponential(a, settings, error)
        if (norm2(exp_a - identity_matrix(n)) > tol) then
            write(stderr, *) "test_matrix_exponential failed: Incorrect matrix "// &
                "values for vanishing matrix."
            test_matrix_exponential = .false.
        end if
        deallocate(exp_a)

    end function test_matrix_exponential

    logical(c_bool) function test_compute_sqrt_and_inv_sqrt() bind(C)
        !
        ! this function tests the subroutine which calculates the square root and
        ! inverse square root of a matrix
        !
        use otr_common, only: compute_sqrt_and_inv_sqrt
        use opentrustregion, only: solver_settings_type
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n = 3

        real(rp) :: a(n, n)
        real(rp), allocatable :: sqrt_a(:, :), inv_sqrt_a(:, :)
        type(solver_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_compute_sqrt_and_inv_sqrt = .true.

        ! set up the settings the routine logs through
        call setup_settings(settings)

        ! initialize a random symmetric positive definite matrix
        call random_number(a)
        a = matmul(a, transpose(a)) + real(n, kind=rp) * identity_matrix(n)

        ! call routine and determine if the square root squares to the matrix and if
        ! the inverse square root is its inverse
        call compute_sqrt_and_inv_sqrt(a, sqrt_a, inv_sqrt_a, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Produced error."
            test_compute_sqrt_and_inv_sqrt = .false.
            return
        end if
        if (norm2(matmul(sqrt_a, sqrt_a) - a) > tol) then
            write(stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Square root "// &
                "does not reproduce the matrix."
            test_compute_sqrt_and_inv_sqrt = .false.
        end if
        if (norm2(matmul(inv_sqrt_a, sqrt_a) - identity_matrix(n)) > tol) then
            write(stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Inverse "// &
                "square root is not the inverse of the square root."
            test_compute_sqrt_and_inv_sqrt = .false.
        end if
        if (norm2(sqrt_a - transpose(sqrt_a)) > tol) then
            write(stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Returned "// &
                "square root is not symmetric."
            test_compute_sqrt_and_inv_sqrt = .false.
        end if
        if (norm2(inv_sqrt_a - transpose(inv_sqrt_a)) > tol) then
            write(stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Returned "// &
                "inverse square root is not symmetric."
            test_compute_sqrt_and_inv_sqrt = .false.
        end if
        deallocate(sqrt_a, inv_sqrt_a)

    end function test_compute_sqrt_and_inv_sqrt

    logical(c_bool) function test_level_shifted_divisors() bind(C)
        !
        ! this function tests the function which returns the divisors of the
        ! level-shifted preconditioner
        !
        use otr_common, only: level_shifted_divisors
        use opentrustregion, only: precond_floor

        real(rp), parameter :: mu = 0.5_rp

        real(rp) :: eigval_pairs(5), expected(5)

        ! assume tests pass
        test_level_shifted_divisors = .true.

        ! initialize eigenvalue pairs whose shifted values are well separated from the
        ! floor on both sides of zero, and ones whose shifted values vanish or lie just
        ! below the floor
        eigval_pairs = [2.0_rp, -1.0_rp, mu, mu + 0.5_rp * precond_floor, &
                        mu - 0.5_rp * precond_floor]

        ! the shifted values are returned unless they are smaller in magnitude than the
        ! floor, in which case they are replaced by the positive floor
        expected = &
            [2.0_rp - mu, -1.0_rp - mu, precond_floor, precond_floor, precond_floor]

        ! call function and determine if the divisors are correct
        if (any(abs(level_shifted_divisors(eigval_pairs, mu) - expected) > &
                tol * abs(expected))) then
            write(stderr, *) "test_level_shifted_divisors failed: Returned "// &
                "divisors wrong."
            test_level_shifted_divisors = .false.
        end if

    end function test_level_shifted_divisors

    logical(c_bool) function test_positive_definite_divisors() bind(C)
        !
        ! this function tests the function which returns the divisors of the
        ! positive-definite preconditioner
        !
        use otr_common, only: positive_definite_divisors
        use opentrustregion, only: precond_floor, precond_rel_floor_factor

        real(rp), parameter :: largest = 4.0_rp

        real(rp) :: eigval_pairs(5), expected(5), rel_floor

        ! assume tests pass
        test_positive_definite_divisors = .true.

        ! floor relative to the largest eigenvalue pair magnitude
        rel_floor = precond_rel_floor_factor * largest

        ! initialize eigenvalue pairs above the relative floor with both signs, the
        ! largest magnitude belonging to a negative one, and ones below it with both
        ! signs, including a vanishing one
        eigval_pairs = [-largest, 0.5_rp * largest, 0.5_rp * rel_floor, &
                        -0.5_rp * rel_floor, 0.0_rp]

        ! the magnitudes are returned unless they lie below the relative floor, in
        ! which case they are replaced by it
        expected = [largest, 0.5_rp * largest, rel_floor, rel_floor, rel_floor]

        ! call function and determine if the divisors are correct
        if (any(abs(positive_definite_divisors(eigval_pairs) - expected) > &
                tol * abs(expected))) then
            write(stderr, *) "test_positive_definite_divisors failed: Returned "// &
                "divisors wrong for relative floor."
            test_positive_definite_divisors = .false.
        end if

        ! call function and determine if the divisors are correct, the absolute floor
        ! takes over once the relative floor falls below it
        eigval_pairs = &
            [0.5_rp * precond_floor, -0.1_rp * precond_floor, 0.0_rp, 0.0_rp, 0.0_rp]
        if (any(abs(positive_definite_divisors(eigval_pairs) - precond_floor) > &
                tol * precond_floor)) then
            write(stderr, *) "test_positive_definite_divisors failed: Returned "// &
                "divisors wrong for absolute floor."
            test_positive_definite_divisors = .false.
        end if

    end function test_positive_definite_divisors

end module otr_common_unit_tests
