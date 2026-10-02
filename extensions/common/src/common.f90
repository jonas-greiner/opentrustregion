! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common

    use opentrustregion, only: rp, ip, settings_type

    implicit none

    ! settings shared by the extensions parameterizing the orbitals in an orbital
    ! basis, which the orbital objects hold
    type, extends(settings_type) :: orbital_settings_type
    contains
        procedure :: init => init_orbital_settings
    end type

    type(orbital_settings_type), parameter :: default_orbital_settings = &
        orbital_settings_type(logger=null(), initialized=.true., verbose=0)

    ! orbital basis in which the orbital rotations are parameterized, holding the
    ! quantities every basis has, whose type-bound procedures perform the operations
    ! every basis provides: moving the current orbitals, the gradient and Hessian
    ! diagonal at the current density, and the eigendecomposition of the static part of
    ! the Hessian with the preconditioners and extra trial vectors built from it
    type, abstract :: orbital_basis_type
        integer(ip) :: n_ao, n_param, n_particle
        real(rp) :: energy
        real(rp), pointer, contiguous :: dm_ao(:, :, :) => null()
        real(rp), allocatable :: grad(:), h_diag(:)
        logical :: hess_eigen_stale = .true.
    contains
        procedure(rotate_orbitals_basis_type), deferred :: rotate_orbitals
        procedure(calculate_grad_h_diag_basis_type), deferred :: calculate_grad_h_diag
        procedure(refresh_hess_eigen_basis_type), deferred :: refresh_hess_eigen
        procedure(rotate_hess_eigenbasis_basis_type), deferred :: &
            rotate_to_hess_eigenbasis
        procedure(rotate_hess_eigenbasis_basis_type), deferred :: &
            rotate_from_hess_eigenbasis
        procedure(get_hess_eigval_pairs_basis_type), deferred :: get_hess_eigval_pairs
        procedure(get_extra_trial_vectors_basis_type), deferred :: &
            get_extra_trial_vectors
    end type

    abstract interface
        subroutine rotate_orbitals_basis_type(self, kappa, settings, error)
            import :: orbital_basis_type, settings_type, rp, ip

            class(orbital_basis_type), intent(inout) :: self
            real(rp), intent(in) :: kappa(:)
            class(settings_type), intent(in) :: settings
            integer(ip), intent(out) :: error
        end subroutine rotate_orbitals_basis_type

        subroutine calculate_grad_h_diag_basis_type(self, fock)
            import :: orbital_basis_type, rp

            class(orbital_basis_type), intent(inout) :: self
            real(rp), intent(in) :: fock(:, :, :)
        end subroutine calculate_grad_h_diag_basis_type

        subroutine refresh_hess_eigen_basis_type(self, settings, error)
            import :: orbital_basis_type, settings_type, ip

            class(orbital_basis_type), intent(inout) :: self
            class(settings_type), intent(in) :: settings
            integer(ip), intent(out) :: error
        end subroutine refresh_hess_eigen_basis_type

        function rotate_hess_eigenbasis_basis_type(self, vector) result(rotated)
            import :: orbital_basis_type, rp

            class(orbital_basis_type), intent(in) :: self
            real(rp), intent(in) :: vector(:)
            real(rp), allocatable :: rotated(:)
        end function rotate_hess_eigenbasis_basis_type

        function get_hess_eigval_pairs_basis_type(self) result(eigval_pairs)
            import :: orbital_basis_type, rp

            class(orbital_basis_type), intent(in) :: self
            real(rp), allocatable :: eigval_pairs(:)
        end function get_hess_eigval_pairs_basis_type

        subroutine get_extra_trial_vectors_basis_type(self, trial_vectors, settings, &
                                                      error)
            import :: orbital_basis_type, settings_type, rp, ip

            class(orbital_basis_type), intent(inout) :: self
            real(rp), intent(out) :: trial_vectors(:, :)
            class(settings_type), intent(in) :: settings
            integer(ip), intent(out) :: error
        end subroutine get_extra_trial_vectors_basis_type
    end interface

contains

    subroutine init_orbital_settings(self, error)
        !
        ! this subroutine initializes the settings shared by the extensions
        ! parameterizing the orbitals in an orbital basis
        !
        use opentrustregion, only: verbosity_error

        class(orbital_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        select type (settings => self)
        type is (orbital_settings_type)
            settings = default_orbital_settings
        class default
            call settings%log("Orbital settings could not be initialized because "// &
                              "initialization routine received the wrong type. The "// &
                              "type orbital_settings_type was likely subclassed "// &
                              "without providing an initialization routine.", &
                              verbosity_error, .true.)
            error = 1
        end select

    end subroutine init_orbital_settings

    subroutine compute_sqrt_and_inv_sqrt(A, sqrtA, inv_sqrtA, settings, error)
        !
        ! this subroutine calculates the square root and inverse square root of a
        ! symmetric positive definite matrix
        !
        use opentrustregion, only: verbosity_error

        real(rp), intent(in) :: A(:, :)
        real(rp), allocatable, intent(out) :: sqrtA(:, :), inv_sqrtA(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n, lwork, info, i
        real(rp), allocatable :: eigvecs(:, :), eigvals(:), work(:)
        character(len=300) :: msg
        external :: dsyev, dgemm

        ! initialize error flag
        error = 0

        ! get matrix dimension
        n = size(A, 1)

        ! allocate eigenvector and eigenvalue arrays
        allocate(eigvecs(n, n), eigvals(n))

        ! copy input because dsyev overwrites it
        eigvecs = A

        ! query optimal workspace size
        lwork = -1
        allocate(work(1))
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)

        ! deallocate work array
        deallocate(work)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! get square roots of eigenvalues
        eigvals = sqrt(eigvals)

        ! allocate and initialize output matrices
        allocate(sqrtA(n, n), inv_sqrtA(n, n))
        sqrtA = 0.0_rp
        inv_sqrtA = 0.0_rp

        ! construct the square root and inverse square root of A
        do i = 1, n
            call dgemm("N", "T", n, n, 1_ip, eigvals(i), eigvecs(:, i), n, &
                       eigvecs(:, i), n, 1.0_rp, sqrtA, n)
            call dgemm("N", "T", n, n, 1_ip, 1.0_rp / eigvals(i), eigvecs(:, i), n, &
                       eigvecs(:, i), n, 1.0_rp, inv_sqrtA, n)
        end do

        deallocate(eigvecs, eigvals)

    end subroutine compute_sqrt_and_inv_sqrt

    function matrix_exponential(A, settings, error) result(expA)
        !
        ! this function calculates the matrix exponential of a real antisymmetric
        ! matrix using the scaling and squaring method applied to the Taylor expansion
        ! of the exponential, the scale factor is derived from the Frobenius norm which
        ! is an upper bound for the spectral norm, convergence is tested against the
        ! last term of the expansion which works because the sum of the Frobenius norms
        ! of two matrices is larger than the Frobenius norm of the sum of both matrices
        !
        use opentrustregion, only: verbosity_error

        real(rp), intent(in) :: A(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        real(rp), allocatable :: expA(:, :)

        integer(ip) :: n, i, power
        real(rp) :: scale, fac, A_norm
        real(rp), allocatable :: An(:, :), tmp(:, :)
        external :: dgemm

        ! initialize error flag
        error = 0

        ! matrix size
        n = size(A, 1)

        ! workspace allocation
        allocate(An(n, n), tmp(n, n), expA(n, n))

        ! compute Frobenius norm of A
        A_norm = sqrt(sum(A**2))

        ! determine scale factor
        power = 3
        if (A_norm > 1.0_rp) then
            power = power + int(ceiling(log(A_norm) / log(2.0_rp)))
        end if
        scale = 2.0_rp**(-power)

        ! initialize exponential and product of matrices
        expA = 0.0_rp
        An = 0.0_rp
        do i = 1, n
            expA(i, i) = 1.0_rp
            An(i, i) = 1.0_rp
        end do

        ! perform Taylor expansion
        i = 1
        fac = 1.0_rp
        do while (A_norm > 1e-12_rp)
            ! get factorial
            fac = fac / real(i, kind=rp)

            ! multiply another matrix and change scale factor accordingly
            call dgemm("N", "N", n, n, n, scale, A, n, An, n, 0.0_rp, tmp, n)
            An = tmp

            ! add next expansion order
            expA = expA + fac * An

            ! convergence check for last expansion order
            A_norm = fac * sqrt(sum(An**2))

            ! check if maximum number of iterations is reached
            i = i + 1
            if (i > 100) then
                call settings%log("Maximum number of iterations for Taylor "// &
                                  "expansion of matrix exponential reached.", &
                                  verbosity_error, .true.)
                error = 1
                return
            end if
        end do
        deallocate(An)

        ! squaring step
        do i = 1, power
            call dgemm("N", "N", n, n, n, 1.0_rp, expA, n, expA, n, 0.0_rp, tmp, n)
            expA = tmp
        end do
        deallocate(tmp)

    end function matrix_exponential

    function level_shifted_divisors(eigval_pairs, mu) result(divisors)
        !
        ! this function returns the divisors of a level-shifted preconditioner acting
        ! in the eigenbasis of the static part of the Hessian: its eigenvalues shifted
        ! by the level shift, with divisors vanishing in magnitude replaced by a small
        ! floor
        !
        use opentrustregion, only: precond_floor

        real(rp), intent(in) :: eigval_pairs(:), mu
        real(rp) :: divisors(size(eigval_pairs))

        ! shift the eigenvalue pairs and floor vanishing divisors
        divisors = eigval_pairs - mu
        where (abs(divisors) < precond_floor) divisors = precond_floor

    end function level_shifted_divisors

    function positive_definite_divisors(eigval_pairs) result(divisors)
        !
        ! this function returns the divisors of a positive-definite preconditioner
        ! acting in the eigenbasis of the static part of the Hessian: the magnitudes of
        ! its eigenvalues, floored relative to the largest of them while guarding
        ! against vanishing eigenvalues
        !
        use opentrustregion, only: precond_floor, precond_rel_floor_factor

        real(rp), intent(in) :: eigval_pairs(:)
        real(rp) :: divisors(size(eigval_pairs))

        real(rp) :: floor_val

        ! floor relative to the largest eigenvalue pair magnitude
        floor_val = max(precond_rel_floor_factor * maxval(abs(eigval_pairs)), &
                        precond_floor)

        ! take magnitudes and floor small divisors
        divisors = abs(eigval_pairs)
        where (divisors < floor_val) divisors = floor_val

    end function positive_definite_divisors

    function channel_rows(lengths) result(rows)
        !
        ! this function returns the first and last row of every particle channel in a
        ! column which stacks the contributions of the channels, of the given lengths,
        ! one after another
        !
        integer(ip), intent(in) :: lengths(:)
        integer(ip) :: rows(2, size(lengths))

        integer(ip) :: i, offset

        ! stack the channels one after another
        offset = 0
        do i = 1, size(lengths, kind=ip)
            rows(1, i) = offset + 1
            rows(2, i) = offset + lengths(i)
            offset = rows(2, i)
        end do

    end function channel_rows

end module otr_common
