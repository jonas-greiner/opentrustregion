! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use otr_common, only: orbital_basis_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! optional outputs requested on every call to the mock density matrix evaluating
    ! functions since the last reset, as bits in the order of their argument lists, so
    ! that tests can check that a routine only asks for what it needs; their number
    ! lets the mocks produce potential matrices which change between calls and thereby
    ! distinguish a cached quantity from one that was incorrectly recomputed from an
    ! unchanged input
    integer(ip), allocatable :: mock_requests(:)

    ! sum of the density matrix passed in the latest call to the mock density matrix
    ! evaluating functions, so that tests can check which density matrix is evaluated
    real(rp) :: mock_dm_sum

    ! multiplier of the density matrix returned by the mock density matrix evaluating
    ! functions, which differs between the first and any subsequent call, further
    ! multipliers for extension-specific potentials follow the same convention
    real(rp), parameter :: mock_fock_factor(2) = [2.0_rp, 5.0_rp]

    ! multiplier of the density matrix returned by the mock response functions
    real(rp), parameter :: mock_response_factor = 2.0_rp

    ! names of the closed-shell and the open-shell case, indexed by the number of
    ! particle channels, for the failure messages of tests covering both
    character(len=12), parameter :: shell_names(2) = &
        [character(len=12) :: "closed-shell", "open-shell"]

    ! mock orbital basis for testing the operations every orbital basis shares, such as
    ! the preconditioners: it rotates vectors into and out of the eigenbasis of the
    ! static part of the Hessian with a fixed orthogonal matrix and returns cached
    ! eigenvalue pairs, which refreshing the eigendecomposition overwrites with those
    ! of the current static part; its orbital updating operations do nothing
    type, extends(orbital_basis_type) :: mock_orbital_basis_type
        real(rp), allocatable :: eigvecs(:, :), eigval_pairs(:), static_eigval_pairs(:)
    contains
        procedure :: rotate_orbitals => mock_rotate_orbitals
        procedure :: calculate_grad_h_diag => mock_calculate_grad_h_diag
        procedure :: refresh_hess_eigen => mock_refresh_hess_eigen
        procedure :: rotate_to_hess_eigenbasis => mock_rotate_to_hess_eigenbasis
        procedure :: rotate_from_hess_eigenbasis => mock_rotate_from_hess_eigenbasis
        procedure :: get_hess_eigval_pairs => mock_get_hess_eigval_pairs
        procedure :: get_extra_trial_vectors => mock_get_extra_trial_vectors
    end type

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

    subroutine mock_get_response_cs(dm, response, error)
        !
        ! this subroutine is a mock response function for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out), target, contiguous :: response(:, :)
        integer(ip), intent(out) :: error

        error = 0
        response = mock_response_factor * dm

    end subroutine mock_get_response_cs

    subroutine mock_get_response_os(dm, response, error)
        !
        ! this subroutine is a mock response function for the open-shell case
        !
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out), target :: response(:, :, :)
        integer(ip), intent(out) :: error

        error = 0
        response = mock_response_factor * dm

    end subroutine mock_get_response_os

    subroutine mock_evaluate_dm_cs(dm, energy, fock, get_response_funptr, error)
        !
        ! this subroutine is a mock density matrix evaluating function for the
        ! closed-shell case, which returns a multiple of the density matrix that
        ! changes between calls so that non-vanishing differences are produced
        !
        use otr_common, only: get_response_cs_type

        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target, contiguous :: fock(:, :)
        procedure(get_response_cs_type), intent(out), optional, pointer :: &
            get_response_funptr
        integer(ip), intent(out) :: error

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(get_response_funptr)))

        error = 0
        mock_dm_sum = sum(dm)
        energy = sum(dm)
        if (present(fock)) fock = mock_factor(mock_fock_factor) * dm
        if (present(get_response_funptr)) get_response_funptr => mock_get_response_cs

    end subroutine mock_evaluate_dm_cs

    subroutine mock_evaluate_dm_os(dm, energy, fock, get_response_funptr, error)
        !
        ! this subroutine is a mock density matrix evaluating function for the
        ! open-shell case, which returns a multiple of the density matrix that changes
        ! between calls so that non-vanishing differences are produced
        !
        use otr_common, only: get_response_os_type

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target :: fock(:, :, :)
        procedure(get_response_os_type), intent(out), optional, pointer :: &
            get_response_funptr
        integer(ip), intent(out) :: error

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(get_response_funptr)))

        error = 0
        mock_dm_sum = sum(dm)
        energy = sum(dm)
        if (present(fock)) fock = mock_factor(mock_fock_factor) * dm
        if (present(get_response_funptr)) get_response_funptr => mock_get_response_os

    end subroutine mock_evaluate_dm_os

    subroutine mock_evaluate_dm_failing_cs(dm, energy, fock, get_response_funptr, error)
        !
        ! this subroutine is a mock density matrix evaluating function for the
        ! closed-shell case, which fails
        !
        use otr_common, only: get_response_cs_type

        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target, contiguous :: fock(:, :)
        procedure(get_response_cs_type), intent(out), optional, pointer :: &
            get_response_funptr
        integer(ip), intent(out) :: error

        error = 1
        energy = sum(dm)
        if (present(fock)) fock = 0.0_rp
        if (present(get_response_funptr)) get_response_funptr => null()

    end subroutine mock_evaluate_dm_failing_cs

    subroutine mock_evaluate_dm_failing_os(dm, energy, fock, get_response_funptr, error)
        !
        ! this subroutine is a mock density matrix evaluating function for the
        ! open-shell case, which fails
        !
        use otr_common, only: get_response_os_type

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target :: fock(:, :, :)
        procedure(get_response_os_type), intent(out), optional, pointer :: &
            get_response_funptr
        integer(ip), intent(out) :: error

        error = 1
        energy = sum(dm)
        if (present(fock)) fock = 0.0_rp
        if (present(get_response_funptr)) get_response_funptr => null()

    end subroutine mock_evaluate_dm_failing_os

    subroutine mock_rotate_orbitals(self, kappa, settings, error)
        !
        ! this subroutine is a mock of moving the orbitals, which does nothing
        !
        use opentrustregion, only: settings_type

        class(mock_orbital_basis_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        error = 0

    end subroutine mock_rotate_orbitals

    subroutine mock_calculate_grad_h_diag(self, fock)
        !
        ! this subroutine is a mock of calculating the gradient and Hessian diagonal,
        ! which does nothing
        !
        class(mock_orbital_basis_type), intent(inout) :: self
        real(rp), intent(in) :: fock(:, :, :)

    end subroutine mock_calculate_grad_h_diag

    subroutine mock_refresh_hess_eigen(self, settings, error)
        !
        ! this subroutine is a mock of refreshing the eigendecomposition of the static
        ! part of the Hessian, which takes over the eigenvalue pairs of the current
        ! static part if the cached ones are stale
        !
        use opentrustregion, only: settings_type

        class(mock_orbital_basis_type), intent(inout) :: self
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        error = 0
        if (.not. self%hess_eigen_stale) return
        self%eigval_pairs = self%static_eigval_pairs
        self%hess_eigen_stale = .false.

    end subroutine mock_refresh_hess_eigen

    function mock_rotate_to_hess_eigenbasis(self, vector) result(rotated)
        !
        ! this function is a mock of rotating a vector into the eigenbasis of the
        ! static part of the Hessian with the eigenvector matrix
        !
        class(mock_orbital_basis_type), intent(in) :: self
        real(rp), intent(in) :: vector(:)
        real(rp), allocatable :: rotated(:)

        rotated = matmul(transpose(self%eigvecs), vector)

    end function mock_rotate_to_hess_eigenbasis

    function mock_rotate_from_hess_eigenbasis(self, vector) result(rotated)
        !
        ! this function is a mock of rotating a vector out of the eigenbasis of the
        ! static part of the Hessian with the eigenvector matrix
        !
        class(mock_orbital_basis_type), intent(in) :: self
        real(rp), intent(in) :: vector(:)
        real(rp), allocatable :: rotated(:)

        rotated = matmul(self%eigvecs, vector)

    end function mock_rotate_from_hess_eigenbasis

    function mock_get_hess_eigval_pairs(self) result(eigval_pairs)
        !
        ! this function is a mock of returning the cached eigenvalue pairs of the
        ! static part of the Hessian
        !
        class(mock_orbital_basis_type), intent(in) :: self
        real(rp), allocatable :: eigval_pairs(:)

        eigval_pairs = self%eigval_pairs

    end function mock_get_hess_eigval_pairs

    subroutine mock_get_extra_trial_vectors(self, trial_vectors, settings, error)
        !
        ! this subroutine is a mock of returning extra trial vectors, which contributes
        ! none
        !
        use opentrustregion, only: settings_type

        class(mock_orbital_basis_type), intent(inout) :: self
        real(rp), intent(out) :: trial_vectors(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        error = 0
        trial_vectors = 0.0_rp

    end subroutine mock_get_extra_trial_vectors

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

    subroutine setup_mock_orbital_basis(basis, n_param)
        !
        ! this subroutine sets up a mock orbital basis with a random orthogonal
        ! eigenvector matrix and random eigenvalue pairs of the current static part of
        ! the Hessian, whose cached ones are those of an earlier static part and marked
        ! stale, so that a routine using them has to refresh them first
        !
        type(mock_orbital_basis_type), intent(out) :: basis
        integer(ip), intent(in) :: n_param

        basis%n_param = n_param
        basis%eigvecs = generate_random_orthogonal_matrix(n_param)
        allocate(basis%static_eigval_pairs(n_param))
        call random_number(basis%static_eigval_pairs)
        basis%static_eigval_pairs = 2.0_rp * basis%static_eigval_pairs - 1.0_rp
        allocate(basis%eigval_pairs(n_param), source=1.0_rp)
        basis%hess_eigen_stale = .true.

    end subroutine setup_mock_orbital_basis

    logical(c_bool) function test_init_orbital_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the settings shared by
        ! the extensions parameterizing the orbitals in an orbital basis
        !
        use otr_common, only: orbital_settings_type, &
                              default_settings => default_orbital_settings
        use otr_common_test_reference, only: operator(==)

        type(orbital_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_init_orbital_settings = .true.

        ! initialize settings
        call settings%init(error)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_init_orbital_settings failed: Function raised error."
            test_init_orbital_settings = .false.
        end if

        ! check settings
        if (.not. (settings == default_settings)) then
            write(stderr, *) "test_init_orbital_settings failed: Settings not "// &
                "initialized correctly."
            test_init_orbital_settings = .false.
        end if

    end function test_init_orbital_settings

    logical(c_bool) function test_refresh_response_orbital_basis() bind(C)
        !
        ! this function tests the subroutine which rebuilds the response callbacks at
        ! the currently stored density matrix of an orbital basis
        !
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_ao, n_particle, n_occ

        type(mock_orbital_basis_type) :: basis
        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        integer(ip) :: i, i_shell, error
        logical :: response_set
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_refresh_response_orbital_basis = .true.

        ! set up the mock orbital basis and a density matrix for every particle channel
        call setup_settings(basis%settings)
        basis%n_ao = n_ao
        do i = 1, n_particle
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do

        ! loop over the closed-shell and the open-shell case
        do i_shell = 1, 2
            case_name = trim(shell_names(i_shell))

            ! store the density matrices of the particle channels of the case, set the
            ! mock density matrix evaluating function and mark the response as stale as
            ! ARH would after moving the density on its own
            basis%n_particle = i_shell
            basis%dm_ao => dm_ao(:, :, :i_shell)
            if (i_shell == 1) then
                basis%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                basis%evaluate_dm_cs => null()
                basis%evaluate_dm_os => mock_evaluate_dm_os
            end if
            basis%response_stale = .true.
            mock_requests = [integer(ip) :: ]

            ! call routine and determine if the density matrix evaluating function was
            ! called for the stored density matrix to rebuild only the response and if
            ! the flag was cleared
            call basis%refresh_response(error)
            if (error /= 0) then
                write(stderr, *) "test_refresh_response_orbital_basis failed: "// &
                    "Produced error for the "//case_name//" case."
                test_refresh_response_orbital_basis = .false.
            end if
            if (abs(mock_dm_sum - sum(dm_ao(:, :, :i_shell))) > tol) then
                write(stderr, *) "test_refresh_response_orbital_basis failed: "// &
                    "Stored density matrix not evaluated for the "//case_name//" case."
                test_refresh_response_orbital_basis = .false.
            end if
            if (size(mock_requests) /= 1) then
                write(stderr, *) "test_refresh_response_orbital_basis failed: "// &
                    "Density matrix evaluating function was not called for the "// &
                    case_name//" case."
                test_refresh_response_orbital_basis = .false.
            end if
            if (any(mock_requests /= 2)) then
                write(stderr, *) "test_refresh_response_orbital_basis failed: "// &
                    "Incorrect outputs requested from density matrix evaluating "// &
                    "function for the "//case_name//" case."
                test_refresh_response_orbital_basis = .false.
            end if
            if (i_shell == 1) then
                response_set = associated(basis%get_response_cs, mock_get_response_cs)
            else
                response_set = associated(basis%get_response_os, mock_get_response_os)
            end if
            if (.not. response_set) then
                write(stderr, *) "test_refresh_response_orbital_basis failed: "// &
                    "Response function not updated for the "//case_name//" case."
                test_refresh_response_orbital_basis = .false.
            end if
            if (basis%response_stale) then
                write(stderr, *) "test_refresh_response_orbital_basis failed: "// &
                    "Response still marked stale after being refreshed for the "// &
                    case_name//" case."
                test_refresh_response_orbital_basis = .false.
            end if
        end do

        ! call routine with a failing density matrix evaluating function and determine
        ! if the error is passed on and the response stays marked stale
        basis%evaluate_dm_os => mock_evaluate_dm_failing_os
        basis%response_stale = .true.
        call basis%refresh_response(error)
        if (error == 0) then
            write(stderr, *) "test_refresh_response_orbital_basis failed: Error of "// &
                "the density matrix evaluating function not passed on."
            test_refresh_response_orbital_basis = .false.
        end if
        if (.not. basis%response_stale) then
            write(stderr, *) "test_refresh_response_orbital_basis failed: Response "// &
                "not marked stale after a failed refresh."
            test_refresh_response_orbital_basis = .false.
        end if

    end function test_refresh_response_orbital_basis

    logical(c_bool) function test_precond_orbital_basis() bind(C)
        !
        ! this function tests the subroutine which applies a level-shifted
        ! preconditioner based on the eigendecomposition of the static part of the
        ! Hessian of an orbital basis
        !
        use otr_common, only: orbital_settings_type
        use opentrustregion, only: precond_floor
        use opentrustregion_unit_tests, only: setup_settings
        use test_reference, only: n_param

        type(mock_orbital_basis_type) :: basis
        type(orbital_settings_type) :: settings
        real(rp) :: residual(n_param), precond_residual(n_param), divisors(n_param), &
                    expected(n_param), mu
        integer(ip) :: error

        ! assume tests pass
        test_precond_orbital_basis = .true.

        ! setup settings object
        call setup_settings(settings)

        ! set up the mock orbital basis with a stale eigendecomposition and a random
        ! residual and level shift
        call setup_mock_orbital_basis(basis, n_param)
        call random_number(residual)
        mu = 0.3_rp

        ! construct the expected preconditioned residual from the eigenvalue pairs of
        ! the current static part
        divisors = basis%static_eigval_pairs - mu
        where (abs(divisors) < precond_floor) divisors = precond_floor
        expected = matmul(basis%eigvecs, &
                          matmul(transpose(basis%eigvecs), residual) / divisors)

        ! call routine and determine if values of the preconditioned residual match,
        ! which requires the eigendecomposition to be refreshed
        call basis%precond(residual, mu, precond_residual, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_precond_orbital_basis failed: Produced error."
            test_precond_orbital_basis = .false.
        end if
        if (norm2(precond_residual - expected) > tol) then
            write(stderr, *) "test_precond_orbital_basis failed: Incorrect "// &
                "preconditioned residual."
            test_precond_orbital_basis = .false.
        end if

    end function test_precond_orbital_basis

    logical(c_bool) function test_precond_pd_orbital_basis() bind(C)
        !
        ! this function tests the subroutine which applies the positive-definite
        ! preconditioner based on the eigendecomposition of the static part of the
        ! Hessian of an orbital basis
        !
        use otr_common, only: orbital_settings_type
        use opentrustregion, only: precond_floor, precond_rel_floor_factor
        use opentrustregion_unit_tests, only: setup_settings
        use test_reference, only: n_param

        type(mock_orbital_basis_type) :: basis
        type(orbital_settings_type) :: settings
        real(rp) :: residual(n_param), precond_residual(n_param), divisors(n_param), &
                    expected(n_param), floor_val
        integer(ip) :: error

        ! assume tests pass
        test_precond_pd_orbital_basis = .true.

        ! setup settings object
        call setup_settings(settings)

        ! set up the mock orbital basis with a stale eigendecomposition and a random
        ! residual
        call setup_mock_orbital_basis(basis, n_param)
        call random_number(residual)

        ! construct the expected preconditioned residual from the eigenvalue pairs of
        ! the current static part
        divisors = abs(basis%static_eigval_pairs)
        floor_val = max(precond_rel_floor_factor * maxval(divisors), precond_floor)
        where (divisors < floor_val) divisors = floor_val
        expected = matmul(basis%eigvecs, &
                          matmul(transpose(basis%eigvecs), residual) / divisors)

        ! call routine and determine if values of the preconditioned residual match,
        ! which requires the eigendecomposition to be refreshed
        call basis%precond_pd(residual, precond_residual, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_precond_pd_orbital_basis failed: Produced error."
            test_precond_pd_orbital_basis = .false.
        end if
        if (norm2(precond_residual - expected) > tol) then
            write(stderr, *) "test_precond_pd_orbital_basis failed: Incorrect "// &
                "preconditioned residual."
            test_precond_pd_orbital_basis = .false.
        end if

    end function test_precond_pd_orbital_basis

    logical(c_bool) function test_fill_extra_trial_vectors_orbital_basis() bind(C)
        !
        ! this function tests the subroutine which fills the extra trial vectors with
        ! the rotations belonging to the most negative eigenvalue pairs of the static
        ! part of the Hessian of an orbital basis
        !
        use test_reference, only: n_param

        real(rp), parameter :: eigval_pairs(n_param) = [-0.2_rp, 0.4_rp, -0.7_rp]
        integer(ip), parameter :: n_extra = 3

        type(mock_orbital_basis_type) :: basis
        real(rp) :: trial_vectors(n_param, n_extra)

        ! assume tests pass
        test_fill_extra_trial_vectors_orbital_basis = .true.

        ! set up the mock orbital basis with a random eigenvector matrix
        call setup_mock_orbital_basis(basis, n_param)

        ! call routine for one positive and two negative eigenvalue pairs, and
        ! determine if the first two vectors are the rotations out of the eigenbasis of
        ! the unit vectors along the negative pairs, the eigenvectors at their indices,
        ! in increasing order of the pairs, and if the slot of the positive pair
        ! vanishes
        call basis%fill_extra_trial_vectors(eigval_pairs, trial_vectors)
        if (norm2(trial_vectors(:, 1) - basis%eigvecs(:, 3)) > tol .or. &
            norm2(trial_vectors(:, 2) - basis%eigvecs(:, 1)) > tol) then
            write(stderr, *) "test_fill_extra_trial_vectors_orbital_basis failed: "// &
                "Incorrect extra trial vectors."
            test_fill_extra_trial_vectors_orbital_basis = .false.
        end if
        if (norm2(trial_vectors(:, 3)) > tol) then
            write(stderr, *) "test_fill_extra_trial_vectors_orbital_basis failed: "// &
                "Slot without a negative eigenvalue pair does not vanish."
            test_fill_extra_trial_vectors_orbital_basis = .false.
        end if

    end function test_fill_extra_trial_vectors_orbital_basis

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

    logical(c_bool) function test_channel_rows() bind(C)
        !
        ! this function tests the function which returns the first and last row of
        ! every particle channel in a column stacking the channels one after another
        !
        use otr_common, only: channel_rows

        integer(ip) :: rows(2, 3)

        ! assume tests pass
        test_channel_rows = .true.

        ! call function for channels of different lengths, including an empty one, and
        ! determine if the channels follow each other without gaps
        rows = channel_rows([4_ip, 0_ip, 3_ip])
        if (any(rows /= reshape([1_ip, 4_ip, &
                                 5_ip, 4_ip, &
                                 5_ip, 7_ip], [2, 3]))) then
            write(stderr, *) "test_channel_rows failed: Incorrect rows."
            test_channel_rows = .false.
        end if

    end function test_channel_rows

end module otr_common_unit_tests
