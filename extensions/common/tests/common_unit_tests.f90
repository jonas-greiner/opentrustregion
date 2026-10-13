! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use otr_common, only: extension_type, orbital_basis_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! name of the mock host callback function which the roles of an extension forwarded
    ! to last
    character(len=40) :: host_called = ""

    ! names of the levels of an extension object, that of the solver and that of its
    ! stability check
    character(len=15), parameter :: level_names(2) = &
        [character(len=15) :: "solver", "stability check"]

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

    ! mock extension for testing the routines which dispatch to an extension object:
    ! its orbital update and Hessian linear transformation record their name and return
    ! the results of the shared mocks with the error code the test sets, and the
    ! orbital update returns no Hessian linear transformation if the test disables it
    type, extends(extension_type) :: mock_extension_type
        character(len=40) :: called_role = ""
        logical :: return_hess_x = .true.
        integer(ip) :: mock_error = 0
    contains
        procedure :: update_orbs => update_orbs_mock_extension
        procedure :: hess_x => hess_x_mock_extension
    end type

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

    function record_host_call(name, context) result(error)
        !
        ! this function records the call of a mock host callback function, checks the
        ! host context it received and returns the error code the test set in it
        !
        use test_reference, only: host_context_type, check_host_context

        character(len=*), intent(in) :: name
        class(*), intent(in), pointer :: context
        integer(ip) :: error

        host_called = name
        call check_host_context(context)
        error = 0
        if (associated(context)) then
            select type (context)
            class is (host_context_type)
                error = int(context%mock_error, kind=ip)
            end select
        end if

    end function record_host_call

    subroutine mock_host_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                     context)
        !
        ! this subroutine is a mock host orbital update function, which records its
        ! call in addition to the shared mock
        !
        use opentrustregion, only: hess_x_type
        use otr_common_mock, only: mock_update_orbs

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, context)
        error = record_host_call("update_orbs", context)

    end subroutine mock_host_update_orbs

    function mock_host_obj_func(kappa, error, context) result(func)
        !
        ! this function is a mock host objective function, which records its call in
        ! addition to the shared mock
        !
        use otr_common_mock, only: mock_obj_func

        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        func = mock_obj_func(kappa, error, context)
        error = record_host_call("obj_func", context)

    end function mock_host_obj_func

    subroutine mock_host_hess_x(x, hess_x, error, context)
        !
        ! this subroutine is a mock host Hessian linear transformation, which records
        ! its call in addition to the shared mock
        !
        use otr_common_mock, only: mock_hess_x

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_hess_x(x, hess_x, error, context)
        error = record_host_call("hess_x", context)

    end subroutine mock_host_hess_x

    subroutine mock_host_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine is a mock host level-shifted preconditioner of the solver,
        ! which records its call in addition to the shared mock
        !
        use otr_common_mock, only: mock_precond

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_precond(residual, mu, precond_residual, error, context)
        error = record_host_call("precond", context)

    end subroutine mock_host_precond

    subroutine mock_host_stability_precond(residual, mu, precond_residual, error, &
                                           context)
        !
        ! this subroutine is a mock host level-shifted preconditioner of the stability
        ! check, which behaves like the solver's but is a different function
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_host_precond(residual, mu, precond_residual, error, context)
        host_called = "stability_precond"

    end subroutine mock_host_stability_precond

    subroutine mock_host_precond_pd(residual, precond_residual, error, context)
        !
        ! this subroutine is a mock host positive-definite preconditioner, which
        ! records its call in addition to the shared mock
        !
        use otr_common_mock, only: mock_precond_pd

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_precond_pd(residual, precond_residual, error, context)
        error = record_host_call("precond_pd", context)

    end subroutine mock_host_precond_pd

    subroutine mock_host_project(vector, error, context)
        !
        ! this subroutine is a mock host projection of the solver
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        error = record_host_call("project", context)
        vector = 2.0_rp * vector

    end subroutine mock_host_project

    subroutine mock_host_stability_project(vector, error, context)
        !
        ! this subroutine is a mock host projection of the stability check, which
        ! behaves like the solver's but is a different function
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_host_project(vector, error, context)
        host_called = "stability_project"

    end subroutine mock_host_stability_project

    subroutine mock_host_get_extra_trial_vectors(trial_vectors, error, context)
        !
        ! this subroutine is a mock host extra trial vector function of the solver,
        ! which records its call in addition to the shared mock
        !
        use otr_common_mock, only: &
            shared_mock_get_extra_trial_vectors => mock_get_extra_trial_vectors

        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call shared_mock_get_extra_trial_vectors(trial_vectors, error, context)
        error = record_host_call("get_extra_trial_vectors", context)

    end subroutine mock_host_get_extra_trial_vectors

    subroutine mock_host_stability_get_extra_trial_vectors(trial_vectors, error, &
                                                           context)
        !
        ! this subroutine is a mock host extra trial vector function of the stability
        ! check, which behaves like the solver's but is a different function
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_host_get_extra_trial_vectors(trial_vectors, error, context)
        host_called = "stability_get_extra_trial_vectors"

    end subroutine mock_host_stability_get_extra_trial_vectors

    function mock_host_conv_check(error, context) result(converged)
        !
        ! this function is a mock host convergence check of the solver
        !
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        error = record_host_call("conv_check", context)
        converged = .true.

    end function mock_host_conv_check

    function mock_host_conv_check_stability(residual, eigval, error, context) &
        result(converged)
        !
        ! this function is a mock host convergence check of the stability check, which
        ! reports convergence for a positive eigenvalue
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: eigval
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        error = record_host_call("conv_check_stability", context)
        converged = eigval > 0.0_rp .and. size(residual) > 0

    end function mock_host_conv_check_stability

    subroutine mock_host_stability_hess_x(x, hess_x, error, context)
        !
        ! this subroutine is a mock host Hessian linear transformation of the stability
        ! check, which behaves like the solver's but is a different function
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_host_hess_x(x, hess_x, error, context)
        host_called = "stability_hess_x"

    end subroutine mock_host_stability_hess_x

    subroutine mock_host_approx_hess_x(x, hess_x, error, context)
        !
        ! this subroutine is a mock host approximate Hessian linear transformation,
        ! which behaves like the Hessian linear transformation but is a different
        ! function
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call mock_host_hess_x(x, hess_x, error, context)
        host_called = "approx_hess_x"

    end subroutine mock_host_approx_hess_x

    subroutine mock_host_logger(message, context)
        !
        ! this subroutine is a mock host logging function of the solver, which appends
        ! the message to the log of the test context
        !
        use opentrustregion_unit_tests, only: logger

        character(len=*), intent(in) :: message
        class(*), intent(in), pointer :: context

        call logger(message, context)
        host_called = "logger"

    end subroutine mock_host_logger

    subroutine mock_host_stability_logger(message, context)
        !
        ! this subroutine is a mock host logging function of the stability check, which
        ! behaves like the solver's but is a different function
        !
        character(len=*), intent(in) :: message
        class(*), intent(in), pointer :: context

        call mock_host_logger(message, context)
        host_called = "stability_logger"

    end subroutine mock_host_stability_logger

    subroutine update_orbs_mock_extension(self, kappa, func, grad, h_diag, &
                                          hess_x_funptr, error)
        !
        ! this subroutine is the orbital update of the mock extension, which records
        ! its call in addition to the shared mock
        !
        use opentrustregion, only: hess_x_type
        use otr_common_mock, only: mock_update_orbs
        use test_reference, only: no_context

        class(mock_extension_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        self%called_role = "update_orbs"
        call mock_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                              no_context)
        if (.not. self%return_hess_x) hess_x_funptr => null()
        error = self%mock_error

    end subroutine update_orbs_mock_extension

    subroutine hess_x_mock_extension(self, x, hess_x, error)
        !
        ! this subroutine is the Hessian linear transformation of the mock extension,
        ! which records its call in addition to the shared mock
        !
        use otr_common_mock, only: mock_hess_x
        use test_reference, only: no_context

        class(mock_extension_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        self%called_role = "hess_x"
        call mock_hess_x(x, hess_x, error, no_context)
        error = self%mock_error

    end subroutine hess_x_mock_extension

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

    subroutine mock_get_response_cs(dm, response, error, context)
        !
        ! this subroutine is a mock response function for the closed-shell case
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out), target, contiguous :: response(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        error = 0
        response = mock_response_factor * dm

    end subroutine mock_get_response_cs

    subroutine mock_get_response_os(dm, response, error, context)
        !
        ! this subroutine is a mock response function for the open-shell case
        !
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out), target :: response(:, :, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        error = 0
        response = mock_response_factor * dm

    end subroutine mock_get_response_os

    subroutine mock_evaluate_dm_cs(dm, energy, fock, get_response_funptr, error, &
                                   context)
        !
        ! this subroutine is a mock density matrix evaluating function for the
        ! closed-shell case, which returns a multiple of the density matrix that
        ! changes between calls so that non-vanishing differences are produced
        !
        use otr_common, only: get_response_cs_type
        use test_reference, only: check_host_context

        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target, contiguous :: fock(:, :)
        procedure(get_response_cs_type), intent(out), optional, pointer :: &
            get_response_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(get_response_funptr)))

        error = 0
        mock_dm_sum = sum(dm)
        energy = sum(dm)
        if (present(fock)) fock = mock_factor(mock_fock_factor) * dm
        if (present(get_response_funptr)) get_response_funptr => mock_get_response_cs

    end subroutine mock_evaluate_dm_cs

    subroutine mock_evaluate_dm_os(dm, energy, fock, get_response_funptr, error, &
                                   context)
        !
        ! this subroutine is a mock density matrix evaluating function for the
        ! open-shell case, which returns a multiple of the density matrix that changes
        ! between calls so that non-vanishing differences are produced
        !
        use otr_common, only: get_response_os_type
        use test_reference, only: check_host_context

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target :: fock(:, :, :)
        procedure(get_response_os_type), intent(out), optional, pointer :: &
            get_response_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context
        call check_host_context(context)

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(get_response_funptr)))

        error = 0
        mock_dm_sum = sum(dm)
        energy = sum(dm)
        if (present(fock)) fock = mock_factor(mock_fock_factor) * dm
        if (present(get_response_funptr)) get_response_funptr => mock_get_response_os

    end subroutine mock_evaluate_dm_os

    subroutine mock_evaluate_dm_failing_cs(dm, energy, fock, get_response_funptr, &
                                           error, context)
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
        class(*), intent(in), pointer :: context

        error = 1
        energy = sum(dm)
        if (present(fock)) fock = 0.0_rp
        if (present(get_response_funptr)) get_response_funptr => null()

    end subroutine mock_evaluate_dm_failing_cs

    subroutine mock_evaluate_dm_failing_os(dm, energy, fock, get_response_funptr, &
                                           error, context)
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
        class(*), intent(in), pointer :: context

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

    subroutine setup_mock_extension(extension, context, stability_context)
        !
        ! this subroutine sets up a mock extension whose levels hold every mock host
        ! callback function, those of the stability check where the stability check has
        ! its own, and the given host contexts, so that a callback function forwarding
        ! to the wrong level is detected by the context it passes, and arms these
        !
        use opentrustregion_unit_tests, only: test_context_type
        use test_reference, only: reset_host_context, n_param

        type(mock_extension_type), intent(out), target :: extension
        type(test_context_type), intent(inout), target :: context, stability_context

        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        associate (host => extension%solver_level%host)
            host%context => context
            host%update_orbs => mock_host_update_orbs
            host%obj_func => mock_host_obj_func
            host%precond => mock_host_precond
            host%precond_pd => mock_host_precond_pd
            host%project => mock_host_project
            host%get_extra_trial_vectors => mock_host_get_extra_trial_vectors
            host%conv_check => mock_host_conv_check
            host%stability_hess_x => mock_host_stability_hess_x
            host%approx_hess_x => mock_host_approx_hess_x
            host%conv_check_stability => mock_host_conv_check_stability
            host%logger => mock_host_logger
        end associate
        associate (host => extension%stability_level%host)
            host%context => stability_context
            host%update_orbs => mock_host_update_orbs
            host%obj_func => mock_host_obj_func
            host%precond => mock_host_stability_precond
            host%precond_pd => mock_host_precond_pd
            host%project => mock_host_stability_project
            host%get_extra_trial_vectors => mock_host_stability_get_extra_trial_vectors
            host%conv_check => mock_host_conv_check
            host%stability_hess_x => mock_host_stability_hess_x
            host%approx_hess_x => mock_host_approx_hess_x
            host%conv_check_stability => mock_host_conv_check_stability
            host%logger => mock_host_stability_logger
        end associate

        ! clear the bookkeeping of the host contexts and arm them
        call reset_host_context()
        context%n_calls = 0
        context%mock_error = 0
        context%log_message = ""
        stability_context%n_calls = 0
        stability_context%mock_error = 0
        stability_context%log_message = ""
        host_called = ""

    end subroutine setup_mock_extension

    function level_context(extension, i_level) result(context)
        !
        ! this function returns the level of the solver (1) or of the stability check
        ! (2) of a mock extension as the context a callback function receives
        !
        type(mock_extension_type), intent(in), target :: extension
        integer(ip), intent(in) :: i_level
        class(*), pointer :: context

        if (i_level == 1) then
            context => extension%solver_level
        else
            context => extension%stability_level
        end if

    end function level_context

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

    logical function check_forwarded(test_name, i_level, host_callback, contexts)
        !
        ! this function checks that a callback function called with the given level of
        ! an extension forwarded to the given mock host callback function with the host
        ! context of that level and not with that of the other level
        !
        use opentrustregion_unit_tests, only: test_context_type
        use test_reference, only: host_context_reached

        character(len=*), intent(in) :: test_name, host_callback
        integer(ip), intent(in) :: i_level
        type(test_context_type), intent(in) :: contexts(2)

        ! assume tests pass
        check_forwarded = .true.

        if (host_called /= host_callback) then
            write(stderr, *) "test_"//test_name//" failed: Forwarded to wrong host "// &
                "callback function for the "//trim(level_names(i_level))//"."
            check_forwarded = .false.
        end if
        if (contexts(3 - i_level)%n_calls /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Forwarded with host "// &
                "context of the other level for the "//trim(level_names(i_level))//"."
            check_forwarded = .false.
        end if
        if (.not. host_context_reached(test_name, contexts(i_level))) &
            check_forwarded = .false.

    end function check_forwarded

    logical function check_role_called(test_name, extension, role)
        !
        ! this function checks that the given role of the mock extension was called
        !
        character(len=*), intent(in) :: test_name, role
        type(mock_extension_type), intent(in) :: extension

        ! assume tests pass
        check_role_called = .true.

        if (extension%called_role /= role) then
            write(stderr, *) "test_"//test_name//" failed: Role "//role// &
                " of the extension object not called."
            check_role_called = .false.
        end if

    end function check_role_called

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

    logical(c_bool) function test_wire_solver_settings_extension() bind(C)
        !
        ! this function tests the subroutine which saves the host's callback functions
        ! and host contexts of solver settings in the host records of the levels of an
        ! extension and wires the role callback functions into the settings
        !
        use opentrustregion, only: &
            solver_settings_type, stability_settings_type, update_orbs_type, &
            obj_func_type, default_solver_settings, default_stability_settings, &
            settings_uninitialized_warning_msg
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              setup_error_logging
        use otr_common, only: orbital_settings_type, host_record_type, &
                              precond_callback, precond_pd_callback, project_callback, &
                              get_extra_trial_vectors_callback, &
                              stability_hess_x_callback, approx_hess_x_callback, &
                              conv_check_stability_callback, logger_callback
        use otr_common_test_reference, only: same_callbacks

        type(mock_extension_type), target :: extension
        type(solver_settings_type) :: solver_settings, expected_settings
        type(stability_settings_type) :: expected_stability_settings
        type(host_record_type) :: expected_host, expected_stability_host
        type(orbital_settings_type) :: settings
        type(test_context_type), target :: context, stability_context, log_context
        type(test_context_type), pointer :: expected_stability_context
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        integer(ip) :: error, i_case, i_field
        character(len=40), parameter :: case_names(2) = &
            [character(len=40) :: "with nested settings of their own", &
             "with nested settings without their own"]
        character(len=41), parameter :: nested_fields(7) = &
            [character(len=41) :: "context", "preconditioner", "projection", &
             "logging function", "extra trial vector function", &
             "approximate Hessian linear transformation", "convergence check"]

        ! assume tests pass
        test_wire_solver_settings_extension = .true.

        ! set up the settings of the extension and the host's orbital update and
        ! objective function
        call setup_settings(settings)
        update_orbs_funptr => mock_host_update_orbs
        obj_func_funptr => mock_host_obj_func

        ! the host record of the solver has to hold the host's callback functions, of
        ! which the convergence check is unset
        expected_host = host_record_type( &
            update_orbs=mock_host_update_orbs, obj_func=mock_host_obj_func, &
            precond=mock_host_precond, precond_pd=mock_host_precond_pd, &
            project=mock_host_project, &
            get_extra_trial_vectors=mock_host_get_extra_trial_vectors, &
            stability_hess_x=mock_host_stability_hess_x, logger=mock_host_logger)

        ! wire solver settings that leave their convergence check unset and whose
        ! nested settings provide a context, a preconditioner, an extra trial vector
        ! function and a convergence check but no projection, logging function or
        ! approximate Hessian linear transformation of their own, which have to be
        ! inherited from the solver or left unset, and solver settings whose nested
        ! settings provide nothing of their own and have to be left untouched
        do i_case = 1, size(case_names)
            call setup_host_solver_settings()
            if (i_case == 1) then
                solver_settings%stability_settings%context => stability_context
                solver_settings%stability_settings%precond => &
                    mock_host_stability_precond
                solver_settings%stability_settings%get_extra_trial_vectors => &
                    mock_host_stability_get_extra_trial_vectors
                solver_settings%stability_settings%conv_check => &
                    mock_host_conv_check_stability
            end if

            ! the host record of the stability check has to hold the nested settings'
            ! own callback functions, the solver's preconditioner, projection and
            ! logging function where the nested settings provide none, as the solver
            ! lets them inherit these, and the Hessian linear transformation of the
            ! stability check the solver settings hold, and the solver settings have to
            ! hand the role callback functions to every role the host provided, as
            ! nested settings of their own do for the roles they provided or inherit
            expected_settings = default_solver_settings
            expected_settings%precond => precond_callback
            expected_settings%precond_pd => precond_pd_callback
            expected_settings%project => project_callback
            expected_settings%get_extra_trial_vectors => &
                get_extra_trial_vectors_callback
            expected_settings%stability_hess_x => stability_hess_x_callback
            expected_settings%logger => logger_callback
            if (i_case == 1) then
                expected_stability_context => stability_context
                expected_stability_host = host_record_type( &
                    precond=mock_host_stability_precond, project=mock_host_project, &
                    stability_hess_x=mock_host_stability_hess_x, &
                    conv_check_stability=mock_host_conv_check_stability, &
                    logger=mock_host_logger)
                expected_stability_host%get_extra_trial_vectors => &
                    mock_host_stability_get_extra_trial_vectors
                associate (nested => expected_settings%stability_settings)
                    nested%precond => precond_callback
                    nested%project => project_callback
                    nested%get_extra_trial_vectors => get_extra_trial_vectors_callback
                    nested%conv_check => conv_check_stability_callback
                    nested%logger => logger_callback
                end associate
            else
                expected_stability_context => context
                expected_stability_host = host_record_type( &
                    precond=mock_host_precond, project=mock_host_project, &
                    stability_hess_x=mock_host_stability_hess_x, &
                    logger=mock_host_logger)
            end if

            ! wire the settings
            call extension%wire_solver_settings(update_orbs_funptr, obj_func_funptr, &
                                                solver_settings, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_wire_solver_settings_extension failed: "// &
                    "Produced error "//trim(case_names(i_case))//"."
                test_wire_solver_settings_extension = .false.
            end if

            ! check that the host records hold the expected callback functions and the
            ! host context of their level, which the stability check inherits from the
            ! solver where the nested settings provide none
            if (.not. same_callbacks(extension%solver_level%host, expected_host) .or. &
                .not. associated(extension%solver_level%host%context, context)) then
                write(stderr, *) "test_wire_solver_settings_extension failed: Host "// &
                    "record of the solver wrong "//trim(case_names(i_case))//"."
                test_wire_solver_settings_extension = .false.
            end if
            if (.not. same_callbacks(extension%stability_level%host, &
                                     expected_stability_host) .or. &
                .not. associated(extension%stability_level%host%context, &
                                 expected_stability_context)) then
                write(stderr, *) "test_wire_solver_settings_extension failed: Host "// &
                    "record of the stability check wrong "//trim(case_names(i_case))// &
                    "."
                test_wire_solver_settings_extension = .false.
            end if

            ! check that the settings hold the expected role callback functions and the
            ! level of the solver as context, and nested settings of their own the
            ! level of the stability check, while nested settings without are left
            ! untouched to inherit the solver's
            if (.not. same_callbacks(solver_settings, expected_settings) .or. &
                .not. associated(solver_settings%context, extension%solver_level)) then
                write(stderr, *) "test_wire_solver_settings_extension failed: "// &
                    "Settings wired wrong "//trim(case_names(i_case))//"."
                test_wire_solver_settings_extension = .false.
            end if
            if ((associated(solver_settings%stability_settings%context) .neqv. &
                 i_case == 1) .or. (i_case == 1 .and. .not. associated( &
                     solver_settings%stability_settings%context, &
                     extension%stability_level))) then
                write(stderr, *) "test_wire_solver_settings_extension failed: "// &
                    "Context of the nested settings wrong "// &
                    trim(case_names(i_case))//"."
                test_wire_solver_settings_extension = .false.
            end if
        end do

        ! wire solver settings that were not initialized, which are initialized after a
        ! warning through their logging function, so that none of their callback
        ! functions or context is saved and no role is wired
        call setup_host_solver_settings()
        solver_settings%initialized = .false.
        call extension%wire_solver_settings(update_orbs_funptr, obj_func_funptr, &
                                            solver_settings, settings, error)
        if (index(context%log_message, " "//settings_uninitialized_warning_msg) == 0) &
            then
            write(stderr, *) "test_wire_solver_settings_extension failed: No "// &
                "warning printed for solver settings that were not initialized."
            test_wire_solver_settings_extension = .false.
        end if
        if (.not. solver_settings%initialized .or. &
            .not. same_callbacks(solver_settings, default_solver_settings) .or. &
            .not. associated(solver_settings%context, extension%solver_level) .or. &
            associated(solver_settings%stability_settings%context) .or. &
            .not. same_callbacks(extension%solver_level%host, host_record_type( &
                update_orbs=mock_host_update_orbs, obj_func=mock_host_obj_func)) .or. &
            associated(extension%solver_level%host%context)) then
            write(stderr, *) "test_wire_solver_settings_extension failed: Solver "// &
                "settings that were not initialized wired wrong."
            test_wire_solver_settings_extension = .false.
        end if

        ! wire solver settings whose nested settings provide a single callback function
        ! or context of their own, which has to suffice for them to be wired with the
        ! level of the stability check and the role callback functions for this one and
        ! for the preconditioner, projection and logging function they inherit, while
        ! the others stay unset
        do i_field = 1, size(nested_fields)
            call setup_host_solver_settings()
            expected_stability_settings = default_stability_settings
            expected_stability_settings%precond => precond_callback
            expected_stability_settings%project => project_callback
            expected_stability_settings%logger => logger_callback
            associate (nested => solver_settings%stability_settings)
                select case (trim(nested_fields(i_field)))
                case ("context")
                    nested%context => stability_context
                case ("preconditioner")
                    nested%precond => mock_host_stability_precond
                case ("projection")
                    nested%project => mock_host_stability_project
                case ("logging function")
                    nested%logger => mock_host_stability_logger
                case ("extra trial vector function")
                    nested%get_extra_trial_vectors => &
                        mock_host_stability_get_extra_trial_vectors
                    expected_stability_settings%get_extra_trial_vectors => &
                        get_extra_trial_vectors_callback
                case ("approximate Hessian linear transformation")
                    nested%approx_hess_x => mock_host_approx_hess_x
                    expected_stability_settings%approx_hess_x => approx_hess_x_callback
                case ("convergence check")
                    nested%conv_check => mock_host_conv_check_stability
                    expected_stability_settings%conv_check => &
                        conv_check_stability_callback
                end select
                call extension%wire_solver_settings(update_orbs_funptr, &
                                                    obj_func_funptr, solver_settings, &
                                                    settings, error)
                if (.not. same_callbacks(nested, expected_stability_settings) .or. &
                    .not. associated(nested%context, extension%stability_level)) then
                    write(stderr, *) "test_wire_solver_settings_extension failed: "// &
                        "Nested settings with only a "//trim(nested_fields(i_field))// &
                        " of their own wired wrong."
                    test_wire_solver_settings_extension = .false.
                end if
            end associate
        end do

        ! wire solver settings again whose solver or nested settings the extension
        ! already wired, which has to produce an error
        call setup_error_logging(settings, log_context)
        call setup_host_solver_settings()
        solver_settings%context => extension%solver_level
        call extension%wire_solver_settings(update_orbs_funptr, obj_func_funptr, &
                                            solver_settings, settings, error)
        if (error == 0) then
            write(stderr, *) "test_wire_solver_settings_extension failed: No error "// &
                "returned for solver settings wired by the extension."
            test_wire_solver_settings_extension = .false.
        end if
        if (log_context%log_message == "") then
            write(stderr, *) "test_wire_solver_settings_extension failed: No error "// &
                "message printed for solver settings wired by the extension."
            test_wire_solver_settings_extension = .false.
        end if
        log_context%log_message = ""
        call setup_host_solver_settings()
        solver_settings%stability_settings%context => extension%stability_level
        call extension%wire_solver_settings(update_orbs_funptr, obj_func_funptr, &
                                            solver_settings, settings, error)
        if (error == 0) then
            write(stderr, *) "test_wire_solver_settings_extension failed: No error "// &
                "returned for nested settings wired by the extension."
            test_wire_solver_settings_extension = .false.
        end if
        if (log_context%log_message == "") then
            write(stderr, *) "test_wire_solver_settings_extension failed: No error "// &
                "message printed for nested settings wired by the extension."
            test_wire_solver_settings_extension = .false.
        end if

    contains

        subroutine setup_host_solver_settings()
            !
            ! this subroutine sets up the host's solver settings with every mock host
            ! callback function of the solver except the convergence check and with the
            ! host context, and nested settings that provide nothing
            !
            call solver_settings%init(error)
            solver_settings%context => context
            solver_settings%precond => mock_host_precond
            solver_settings%precond_pd => mock_host_precond_pd
            solver_settings%project => mock_host_project
            solver_settings%get_extra_trial_vectors => mock_host_get_extra_trial_vectors
            solver_settings%stability_hess_x => mock_host_stability_hess_x
            solver_settings%logger => mock_host_logger
            context%log_message = ""

        end subroutine setup_host_solver_settings

    end function test_wire_solver_settings_extension

    logical(c_bool) function test_extension_level() bind(C)
        !
        ! this function tests the function which returns the level of an extension
        ! object a callback function receives as context
        !
        use otr_common, only: extension_level
        use test_reference, only: host_context

        type(mock_extension_type), target :: extension
        class(*), pointer :: context

        ! assume tests pass
        test_extension_level = .true.

        ! a level has to be returned
        context => extension%stability_level
        if (.not. associated(extension_level(context), extension%stability_level)) then
            write(stderr, *) "test_extension_level failed: Level not returned."
            test_extension_level = .false.
        end if

        ! a context of another type has to return none
        context => host_context
        if (associated(extension_level(context))) then
            write(stderr, *) "test_extension_level failed: Level returned for a "// &
                "context of another type."
            test_extension_level = .false.
        end if

        ! a missing context has to return none
        context => null()
        if (associated(extension_level(context))) then
            write(stderr, *) "test_extension_level failed: Level returned without "// &
                "a context."
            test_extension_level = .false.
        end if

    end function test_extension_level

    logical(c_bool) function test_update_orbs_callback() bind(C)
        !
        ! this function tests the orbital update, which dispatches to the extension
        ! object whose level it receives as context
        !
        use opentrustregion, only: hess_x_type
        use otr_common, only: update_orbs_callback
        use otr_common_mock, only: mock_hess_x
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        class(*), pointer :: context
        real(rp) :: kappa(n_param), func, grad(n_param), h_diag(n_param)
        procedure(hess_x_type), pointer :: hess_x_funptr
        integer(ip) :: error

        ! assume tests pass
        test_update_orbs_callback = .true.

        ! set up a step
        call random_number(kappa)

        ! call the orbital update with the level of the solver as context and determine
        ! if it calls that of the object and returns its results, Hessian linear
        ! transformation and error
        extension%mock_error = 3
        extension%solver_level%extension => extension
        context => extension%solver_level
        call update_orbs_callback(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                  context)
        if (.not. check_role_called("update_orbs_callback", extension, "update_orbs")) &
            test_update_orbs_callback = .false.
        if (error /= 3) then
            write(stderr, *) "test_update_orbs_callback failed: Error of the "// &
                "object not passed on."
            test_update_orbs_callback = .false.
        end if
        if (abs(func - sum(kappa)) > tol .or. &
            any(abs(grad - 2.0_rp * kappa) > tol) .or. &
            any(abs(h_diag - 3.0_rp * kappa) > tol)) then
            write(stderr, *) "test_update_orbs_callback failed: Results of the "// &
                "object not returned."
            test_update_orbs_callback = .false.
        end if
        if (.not. associated(hess_x_funptr, mock_hess_x)) then
            write(stderr, *) "test_update_orbs_callback failed: Hessian linear "// &
                "transformation of the object not returned."
            test_update_orbs_callback = .false.
        end if

        ! call the orbital update without a level, which has to produce an error and no
        ! Hessian linear transformation
        context => null()
        call update_orbs_callback(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                  context)
        if (error == 0) then
            write(stderr, *) "test_update_orbs_callback failed: No error produced "// &
                "without a level."
            test_update_orbs_callback = .false.
        end if
        if (associated(hess_x_funptr)) then
            write(stderr, *) "test_update_orbs_callback failed: Hessian linear "// &
                "transformation returned without a level."
            test_update_orbs_callback = .false.
        end if

    end function test_update_orbs_callback

    logical(c_bool) function test_obj_func_callback() bind(C)
        !
        ! this function tests the objective function, which forwards to the host's one
        ! saved in the level of an extension object it receives as context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: obj_func_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: kappa(n_param), func
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=8) :: "obj_func", "obj_func"]

        ! assume tests pass
        test_obj_func_callback = .true.

        ! set up a step
        call random_number(kappa)

        ! call the objective function with the levels of the solver and of the
        ! stability check as context and determine if it forwards to the host's one of
        ! the level with its host context and returns its results and error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            func = obj_func_callback(kappa, error, context)
            if (error /= 3) then
                write(stderr, *) "test_obj_func_callback failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_obj_func_callback = .false.
            end if
            if (abs(func - sum(kappa)) > tol) then
                write(stderr, *) "test_obj_func_callback failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_obj_func_callback = .false.
            end if
            if (.not. check_forwarded("obj_func_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_obj_func_callback = .false.
        end do

        ! call the objective function without a level, which has to produce an error
        context => null()
        func = obj_func_callback(kappa, error, context)
        if (error == 0) then
            write(stderr, *) "test_obj_func_callback failed: No error produced "// &
                "without a level."
            test_obj_func_callback = .false.
        end if

    end function test_obj_func_callback

    logical(c_bool) function test_hess_x_callback() bind(C)
        !
        ! this function tests the Hessian linear transformation, which dispatches to
        ! the extension object whose level it receives as context
        !
        use otr_common, only: hess_x_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        class(*), pointer :: context
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error, i_level

        ! assume tests pass
        test_hess_x_callback = .true.

        ! set up a vector
        call random_number(x)

        ! call the Hessian linear transformation with the levels of the solver and of
        ! the stability check as context, which both call it, and determine if it calls
        ! that of the object and returns its result and error
        extension%mock_error = 3
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        do i_level = 1, 2
            context => level_context(extension, i_level)
            extension%called_role = ""
            call hess_x_callback(x, hess_x, error, context)
            if (.not. check_role_called("hess_x_callback", extension, "hess_x")) &
                test_hess_x_callback = .false.
            if (error /= 3) then
                write(stderr, *) "test_hess_x_callback failed: Error of the object "// &
                    "not passed on for the "//trim(level_names(i_level))//"."
                test_hess_x_callback = .false.
            end if
            if (any(abs(hess_x - 4.0_rp * x) > tol)) then
                write(stderr, *) "test_hess_x_callback failed: Result of the "// &
                    "object not returned for the "//trim(level_names(i_level))//"."
                test_hess_x_callback = .false.
            end if
        end do

        ! call the Hessian linear transformation without a level, which has to produce
        ! an error
        context => null()
        call hess_x_callback(x, hess_x, error, context)
        if (error == 0) then
            write(stderr, *) "test_hess_x_callback failed: No error produced "// &
                "without a level."
            test_hess_x_callback = .false.
        end if

    end function test_hess_x_callback

    logical(c_bool) function test_precond_callback() bind(C)
        !
        ! this function tests the level-shifted preconditioner, which forwards to the
        ! host's one saved in the level of an extension object it receives as context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: precond_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: residual(n_param), precond_residual(n_param)
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=17) :: "precond", "stability_precond"]

        ! assume tests pass
        test_precond_callback = .true.

        ! set up a vector
        call random_number(residual)

        ! call the preconditioner with the levels of the solver and of the stability
        ! check as context and determine if it forwards to the host's one of the level
        ! with its host context and returns its results and error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            call precond_callback(residual, 2.0_rp, precond_residual, error, context)
            if (error /= 3) then
                write(stderr, *) "test_precond_callback failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_precond_callback = .false.
            end if
            if (any(abs(precond_residual - 2.0_rp * residual) > tol)) then
                write(stderr, *) "test_precond_callback failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_precond_callback = .false.
            end if
            if (.not. check_forwarded("precond_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_precond_callback = .false.
        end do

        ! call the preconditioner without a level, which has to produce an error
        context => null()
        call precond_callback(residual, 2.0_rp, precond_residual, error, context)
        if (error == 0) then
            write(stderr, *) "test_precond_callback failed: No error produced "// &
                "without a level."
            test_precond_callback = .false.
        end if

    end function test_precond_callback

    logical(c_bool) function test_precond_pd_callback() bind(C)
        !
        ! this function tests the positive-definite preconditioner, which forwards to
        ! the host's one saved in the level of an extension object it receives as
        ! context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: precond_pd_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: residual(n_param), precond_residual(n_param)
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=10) :: "precond_pd", "precond_pd"]

        ! assume tests pass
        test_precond_pd_callback = .true.

        ! set up a vector
        call random_number(residual)

        ! call the preconditioner with the levels of the solver and of the stability
        ! check as context and determine if it forwards to the host's one of the level
        ! with its host context and returns its results and error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            call precond_pd_callback(residual, precond_residual, error, context)
            if (error /= 3) then
                write(stderr, *) "test_precond_pd_callback failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_precond_pd_callback = .false.
            end if
            if (any(abs(precond_residual - 3.0_rp * residual) > tol)) then
                write(stderr, *) "test_precond_pd_callback failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_precond_pd_callback = .false.
            end if
            if (.not. check_forwarded("precond_pd_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_precond_pd_callback = .false.
        end do

        ! call the preconditioner without a level, which has to produce an error
        context => null()
        call precond_pd_callback(residual, precond_residual, error, context)
        if (error == 0) then
            write(stderr, *) "test_precond_pd_callback failed: No error produced "// &
                "without a level."
            test_precond_pd_callback = .false.
        end if

    end function test_precond_pd_callback

    logical(c_bool) function test_project_callback() bind(C)
        !
        ! this function tests the projection, which forwards to the host's one saved in
        ! the level of an extension object it receives as context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: project_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: vector(n_param), projected(n_param)
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=17) :: "project", "stability_project"]

        ! assume tests pass
        test_project_callback = .true.

        ! set up a vector
        call random_number(vector)

        ! call the projection with the levels of the solver and of the stability check
        ! as context and determine if it forwards to the host's one of the level with
        ! its host context and returns its results and error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            projected = vector
            call project_callback(projected, error, context)
            if (error /= 3) then
                write(stderr, *) "test_project_callback failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_project_callback = .false.
            end if
            if (any(abs(projected - 2.0_rp * vector) > tol)) then
                write(stderr, *) "test_project_callback failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_project_callback = .false.
            end if
            if (.not. check_forwarded("project_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_project_callback = .false.
        end do

        ! call the projection without a level, which has to produce an error
        context => null()
        projected = vector
        call project_callback(projected, error, context)
        if (error == 0) then
            write(stderr, *) "test_project_callback failed: No error produced "// &
                "without a level."
            test_project_callback = .false.
        end if

    end function test_project_callback

    logical(c_bool) function test_get_extra_trial_vectors_callback() bind(C)
        !
        ! this function tests the extra trial vector function, which forwards to the
        ! host's one saved in the level of an extension object it receives as context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: get_extra_trial_vectors_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: trial_vectors(n_param, 2)
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=33) :: "get_extra_trial_vectors", &
             "stability_get_extra_trial_vectors"]

        ! assume tests pass
        test_get_extra_trial_vectors_callback = .true.

        ! call the extra trial vector function with the levels of the solver and of the
        ! stability check as context and determine if it forwards to the host's one of
        ! the level with its host context and returns its results and error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            call get_extra_trial_vectors_callback(trial_vectors, error, context)
            if (error /= 3) then
                write(stderr, *) "test_get_extra_trial_vectors_callback failed: "// &
                    "Error of the host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_get_extra_trial_vectors_callback = .false.
            end if
            if (any(abs(trial_vectors(:, 2) - 2.0_rp) > tol)) then
                write(stderr, *) "test_get_extra_trial_vectors_callback failed: "// &
                    "Result of the host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_get_extra_trial_vectors_callback = .false.
            end if
            if (.not. check_forwarded("get_extra_trial_vectors_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_get_extra_trial_vectors_callback = .false.
        end do

        ! call the extra trial vector function without a level, which has to produce an
        ! error
        context => null()
        call get_extra_trial_vectors_callback(trial_vectors, error, context)
        if (error == 0) then
            write(stderr, *) "test_get_extra_trial_vectors_callback failed: No "// &
                "error produced without a level."
            test_get_extra_trial_vectors_callback = .false.
        end if

    end function test_get_extra_trial_vectors_callback

    logical(c_bool) function test_conv_check_callback() bind(C)
        !
        ! this function tests the convergence check of the solver, which forwards to
        ! the host's one saved in the level of an extension object it receives as
        ! context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: conv_check_callback

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        logical :: converged
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=10) :: "conv_check", "conv_check"]

        ! assume tests pass
        test_conv_check_callback = .true.

        ! call the convergence check with the levels of the solver and of the stability
        ! check as context and determine if it forwards to the host's one of the level
        ! with its host context and returns its results and error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            converged = conv_check_callback(error, context)
            if (error /= 3) then
                write(stderr, *) "test_conv_check_callback failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_conv_check_callback = .false.
            end if
            if (.not. converged) then
                write(stderr, *) "test_conv_check_callback failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_conv_check_callback = .false.
            end if
            if (.not. check_forwarded("conv_check_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_conv_check_callback = .false.
        end do

        ! call the convergence check without a level, which has to produce an error
        context => null()
        converged = conv_check_callback(error, context)
        if (error == 0) then
            write(stderr, *) "test_conv_check_callback failed: No error produced "// &
                "without a level."
            test_conv_check_callback = .false.
        end if

    end function test_conv_check_callback

    logical(c_bool) function test_stability_hess_x_callback() bind(C)
        !
        ! this function tests the Hessian linear transformation of the stability check,
        ! which forwards to the host's one saved in the level of an extension object it
        ! receives as context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: stability_hess_x_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=16) :: "stability_hess_x", "stability_hess_x"]

        ! assume tests pass
        test_stability_hess_x_callback = .true.

        ! set up a vector
        call random_number(x)

        ! call the Hessian linear transformation with the levels of the solver and of
        ! the stability check as context and determine if it forwards to the host's one
        ! of the level with its host context and returns its results and error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            call stability_hess_x_callback(x, hess_x, error, context)
            if (error /= 3) then
                write(stderr, *) "test_stability_hess_x_callback failed: Error of "// &
                    "the host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_stability_hess_x_callback = .false.
            end if
            if (any(abs(hess_x - 4.0_rp * x) > tol)) then
                write(stderr, *) "test_stability_hess_x_callback failed: Result of "// &
                    "the host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_stability_hess_x_callback = .false.
            end if
            if (.not. check_forwarded("stability_hess_x_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_stability_hess_x_callback = .false.
        end do

        ! call the Hessian linear transformation without a level, which has to produce
        ! an error
        context => null()
        call stability_hess_x_callback(x, hess_x, error, context)
        if (error == 0) then
            write(stderr, *) "test_stability_hess_x_callback failed: No error "// &
                "produced without a level."
            test_stability_hess_x_callback = .false.
        end if

    end function test_stability_hess_x_callback

    logical(c_bool) function test_approx_hess_x_callback() bind(C)
        !
        ! this function tests the approximate Hessian linear transformation, which
        ! forwards to the host's one saved in the level of an extension object it
        ! receives as context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: approx_hess_x_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=13) :: "approx_hess_x", "approx_hess_x"]

        ! assume tests pass
        test_approx_hess_x_callback = .true.

        ! set up a vector
        call random_number(x)

        ! call the approximate Hessian linear transformation with the levels of the
        ! solver and of the stability check as context and determine if it forwards to
        ! the host's one of the level with its host context and returns its results and
        ! error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            call approx_hess_x_callback(x, hess_x, error, context)
            if (error /= 3) then
                write(stderr, *) "test_approx_hess_x_callback failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_approx_hess_x_callback = .false.
            end if
            if (any(abs(hess_x - 4.0_rp * x) > tol)) then
                write(stderr, *) "test_approx_hess_x_callback failed: Result of "// &
                    "the host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_approx_hess_x_callback = .false.
            end if
            if (.not. check_forwarded("approx_hess_x_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_approx_hess_x_callback = .false.
        end do

        ! call the approximate Hessian linear transformation without a level, which has
        ! to produce an error
        context => null()
        call approx_hess_x_callback(x, hess_x, error, context)
        if (error == 0) then
            write(stderr, *) "test_approx_hess_x_callback failed: No error "// &
                "produced without a level."
            test_approx_hess_x_callback = .false.
        end if

    end function test_approx_hess_x_callback

    logical(c_bool) function test_conv_check_stability_callback() bind(C)
        !
        ! this function tests the convergence check of the stability check, which
        ! forwards to the host's one saved in the level of an extension object it
        ! receives as context
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: conv_check_stability_callback
        use test_reference, only: n_param

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: residual(n_param)
        logical :: converged
        integer(ip) :: error, i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=20) :: "conv_check_stability", "conv_check_stability"]

        ! assume tests pass
        test_conv_check_stability_callback = .true.

        ! set up a residual
        call random_number(residual)

        ! call the convergence check with the levels of the solver and of the stability
        ! check as context and determine if it forwards to the host's one of the level
        ! with its host context and returns its results and error
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            contexts(i_level)%mock_error = 3
            context => level_context(extension, i_level)
            converged = conv_check_stability_callback(residual, 1.0_rp, error, context)
            if (error /= 3) then
                write(stderr, *) "test_conv_check_stability_callback failed: Error "// &
                    "of the host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_conv_check_stability_callback = .false.
            end if
            if (.not. converged) then
                write(stderr, *) "test_conv_check_stability_callback failed: "// &
                    "Result of the host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_conv_check_stability_callback = .false.
            end if
            if (.not. check_forwarded("conv_check_stability_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_conv_check_stability_callback = .false.
        end do

        ! call the convergence check without a level, which has to produce an error
        context => null()
        converged = conv_check_stability_callback(residual, 1.0_rp, error, context)
        if (error == 0) then
            write(stderr, *) "test_conv_check_stability_callback failed: No error "// &
                "produced without a level."
            test_conv_check_stability_callback = .false.
        end if

    end function test_conv_check_stability_callback

    logical(c_bool) function test_logger_callback() bind(C)
        !
        ! this function tests the logging function, which forwards to the host's one
        ! saved in the level of an extension object it receives as context and discards
        ! the message without one
        !
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: logger_callback

        type(mock_extension_type), target :: extension
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        integer(ip) :: i_level
        character(len=*), parameter :: host_callbacks(2) = &
            [character(len=16) :: "logger", "stability_logger"]

        ! assume tests pass
        test_logger_callback = .true.

        ! log a message with the levels of the solver and of the stability check as
        ! context and determine if it reaches the host's logging function of the level
        ! with its host context
        do i_level = 1, 2
            call setup_mock_extension(extension, contexts(1), contexts(2))
            context => level_context(extension, i_level)
            call logger_callback("message", context)
            if (contexts(i_level)%log_message /= "message") then
                write(stderr, *) "test_logger_callback failed: Message not passed "// &
                    "to the host's logging function for the "// &
                    trim(level_names(i_level))//"."
                test_logger_callback = .false.
            end if
            if (.not. check_forwarded("logger_callback", i_level, &
                                      trim(host_callbacks(i_level)), contexts)) &
                test_logger_callback = .false.
        end do

        ! log a message without a level, which has to be discarded
        context => null()
        call logger_callback("message", context)

    end function test_logger_callback

    logical(c_bool) function test_refresh_response_orbital_basis() bind(C)
        !
        ! this function tests the subroutine which rebuilds the response callbacks at
        ! the currently stored density matrix of an orbital basis
        !
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_ao, n_particle, n_occ, shell_names
        use test_reference, only: arm_host_context, host_context_reached

        type(mock_orbital_basis_type) :: basis
        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        integer(ip) :: i, i_shell, error
        logical :: response_set
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_refresh_response_orbital_basis = .true.

        ! set up the mock orbital basis and a density matrix for every particle channel
        call setup_settings(basis%settings)
        call arm_host_context(basis%settings)
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

        ! check that the host context reached the callback functions unchanged
        if (.not. host_context_reached("refresh_response_orbital_basis")) &
            test_refresh_response_orbital_basis = .false.

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
        a = reshape([0.0_rp, -angle, &
                     angle, 0.0_rp], [n, n])

        ! initialize expected rotation matrix
        expected = reshape([cos(angle), -sin(angle), &
                            sin(angle), cos(angle)], [n, n])

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
