! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn

    use opentrustregion, only: rp, ip, kw_len, settings_type, update_orbs_type, &
                               hess_x_type, add_error_origin
    use otr_common, only: extension_type, provided_roles_type, error_qn

    implicit none

    ! interfaces for callback functions
    abstract interface
        subroutine transport_type(geodesic, tangent_vector, error, context)
            import :: rp, ip

            real(rp), intent(in), target :: geodesic(:)
            real(rp), intent(inout), target :: tangent_vector(:)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context

        end subroutine transport_type
    end interface

    abstract interface
        subroutine init_hess_type(vector, error, context)
            import :: rp, ip

            real(rp), intent(inout), target :: vector(:)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context

        end subroutine init_hess_type
    end interface

    ! derived type for settings
    type, extends(settings_type) :: qn_settings_type
        character(len=kw_len) :: hess_update_scheme
        integer(ip) :: max_points
    contains
        procedure :: init => init_qn_settings
    end type

    type(qn_settings_type), parameter :: default_qn_settings = &
        qn_settings_type(logger=null(), initialized=.true., verbose=0, max_points=10, &
                         hess_update_scheme="sr1")

    ! updating schemes
    character(len=kw_len), parameter :: hess_update_schemes(2) = &
        [character(len=kw_len) :: "sr1", "bfgs"]

    ! warning for settings that were not initialized
    character(len=*), parameter :: qn_settings_uninitialized_warning_msg = &
        "Quasi-Newton settings were not initialized. All quasi-Newton settings are "// &
        "set to default values"

    ! quasi-Newton object, holding the history of steps and gradient differences, the
    ! settings and the host's transport and initial Hessian functions, whose updating
    ! schemes extend it with their Hessian linear transformation
    type, abstract, extends(extension_type) :: updating_type
        real(rp), allocatable :: kappa_list(:, :), grad_diff_list(:, :), last_grad(:)
        integer(ip) :: n_points = 0
        type(qn_settings_type) :: settings
        procedure(transport_type), pointer, nopass :: transport => null()
        procedure(init_hess_type), pointer, nopass :: init_hess => null()
    contains
        procedure :: update_orbs => update_orbs_qn
        procedure :: init => init_qn
        procedure :: add => add_qn
        procedure(skip_add_interface), deferred :: skip_add
    end type

    type, extends(updating_type) :: sr1_updating_type
    contains
        procedure :: hess_x => hess_x_sr1
        procedure :: skip_add => sr1_skip_add
    end type

    type, extends(updating_type) :: bfgs_updating_type
    contains
        procedure :: hess_x => hess_x_bfgs
        procedure :: skip_add => bfgs_skip_add
    end type

    abstract interface
        function skip_add_interface(self, kappa, grad_diff, error) result(skip)
            import updating_type, rp, ip

            class(updating_type), intent(in) :: self
            real(rp), intent(in) :: kappa(:), grad_diff(:)
            logical :: skip
            integer(ip), intent(out) :: error
        end function skip_add_interface
    end interface

contains

    subroutine qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                          transport_funptr, init_hess_funptr, update_orbs_qn_funptr, &
                          obj_func_qn_funptr, solver_settings, solver_state, error, &
                          settings)
        !
        ! this subroutine returns the quasi-Newton orbital updating function, which
        ! wraps the host's, and the objective function, and wires the solver settings
        ! with a new quasi-Newton object, which they own and through which the host's
        ! callback functions keep receiving the host context, and which reads the state
        ! of the solver from the record the host also passes to the solver
        !
        use opentrustregion, only: obj_func_type, solver_state_type, &
                                   warn_uninitialized_settings
        use otr_common, only: extension_solver_settings_type, update_orbs_callback, &
                              obj_func_callback

        integer(ip), intent(in) :: n_param
        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        procedure(transport_type), intent(in), pointer :: transport_funptr
        procedure(init_hess_type), intent(in), pointer :: init_hess_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_qn_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_qn_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(qn_settings_type), intent(in), optional :: settings

        type(qn_settings_type) :: local_settings
        class(extension_type), allocatable, target :: object

        ! initialize error flag and the returned functions
        error = 0
        update_orbs_qn_funptr => null()
        obj_func_qn_funptr => null()

        ! work on a copy of the given settings, which is initialized if they were not,
        ! after printing a warning regardless of the verbosity, since the default
        ! settings would silence it, or, without any, of the default settings
        if (present(settings)) then
            local_settings = settings
            if (.not. local_settings%initialized) then
                call warn_uninitialized_settings(local_settings, &
                                                 qn_settings_uninitialized_warning_msg)
                call local_settings%init(error)
            end if
        else
            local_settings = default_qn_settings
        end if

        ! perform sanity check, which also converts the updating scheme to lowercase
        call qn_sanity_check(local_settings, n_param, error)
        call add_error_origin(error, error_qn, local_settings)
        if (error /= 0) return

        ! create the object of the updating scheme
        if (local_settings%hess_update_scheme == "sr1") then
            allocate(sr1_updating_type :: object)
        else
            allocate(bfgs_updating_type :: object)
        end if

        ! initialize the object with the host's transport and initial Hessian functions
        select type (object)
        class is (updating_type)
            call object%init(n_param, local_settings)
            object%transport => transport_funptr
            object%init_hess => init_hess_funptr
        end select

        ! wire the solver settings with the object, which reads the record of the state
        ! and saves the host's orbital update, objective function, callback functions
        ! and host contexts, and let the settings own it
        call solver_settings%wire(object, solver_state, local_settings, error, &
                                  update_orbs_funptr, obj_func_funptr)
        call add_error_origin(error, error_qn, local_settings)
        if (error /= 0) return

        ! get the orbital updating function and the objective function of the object
        update_orbs_qn_funptr => update_orbs_callback
        obj_func_qn_funptr => obj_func_callback

    end subroutine qn_factory

    subroutine qn_sanity_check(settings, n_param, error)
        !
        ! this subroutine performs a sanity check for quasi-Newton input parameters,
        ! after converting the updating scheme to lowercase
        !
        use opentrustregion, only: verbosity_error, string_to_lowercase

        type(qn_settings_type), intent(inout) :: settings
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! check that the number of parameters is positive
        if (n_param < 1) then
            call settings%log("Number of parameters should be larger than 0.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that the number of stored points is not negative, while zero stores no
        ! history
        if (settings%max_points < 0) then
            call settings%log("Maximum number of points should not be negative.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! convert strings to lowercase
        settings%hess_update_scheme = string_to_lowercase(settings%hess_update_scheme)

        ! check for character options
        if (.not. any(settings%hess_update_scheme == hess_update_schemes)) then
            call settings%log( &
                "Quasi-Newton updating scheme option unknown. Possible values are "// &
                """sr1"" (symmetric rank-1 update) and ""bfgs"" "// &
                "(Broyden-Fletcher-Goldfarb-Shanno update).", verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine qn_sanity_check

    recursive subroutine update_orbs_qn(self, kappa, func, grad, h_diag, &
                                        hess_x_funptr, error, context)
        !
        ! this subroutine is the quasi-Newton orbital updating function, which updates
        ! the orbitals with the host's orbital update and host context, replaces the
        ! Hessian linear transformation the host returns by the quasi-Newton one and
        ! adds the new point to the history
        !
        use otr_common, only: hess_x_callback

        class(updating_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! initialize error flag
        error = 0

        ! update orbitals
        call self%solver_level%host%update_orbs(kappa, func, grad, h_diag, &
                                                hess_x_funptr, error, &
                                                self%solver_level%host%context)
        if (error /= 0) return

        ! get pointer to the quasi-Newton Hessian linear transformation of the object
        hess_x_funptr => hess_x_callback

        ! add the new point, also the starting point the first call reaches without a
        ! step
        if (sum(abs(kappa)) > 0.0_rp .or. .not. allocated(self%last_grad)) then
            call self%add(kappa, grad, error)
            if (error /= 0) return
        end if

    end subroutine update_orbs_qn

    subroutine init_qn_settings(self, error)
        !
        ! this subroutine initializes the quasi-Newton settings
        !
        use opentrustregion, only: verbosity_error

        class(qn_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        select type (settings => self)
        type is (qn_settings_type)
            settings = default_qn_settings
        class default
            call settings%log("Quasi-Newton settings could not be initialized "// &
                              "because initialization routine received the wrong "// &
                              "type. The type qn_settings_type was likely "// &
                              "subclassed without providing an initialization "// &
                              "routine.", verbosity_error, .true.)
            error = 1
        end select

    end subroutine init_qn_settings

    subroutine add_qn(self, kappa, grad, error)
        !
        ! this subroutine adds a new rotation and gradient to the corresponding lists
        ! and transports the history accordingly; the gradient becomes the last one
        ! also if the skipping criterion rejects the pair, so that the next pair starts
        ! at the current point
        !
        class(updating_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), grad(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: grad_diff(:)
        integer(ip) :: n_param, i
        logical :: skip

        ! initialize error flag
        error = 0

        ! get number of parameters
        n_param = size(kappa)

        ! transport history
        do i = 1, self%n_points
            call self%transport(kappa, self%kappa_list(:, i), error, &
                                self%settings%context)
            if (error /= 0) return
            call self%transport(kappa, self%grad_diff_list(:, i), error, &
                                self%settings%context)
            if (error /= 0) return
        end do

        ! check if last gradient is available to calculate gradient difference
        skip = .false.
        if (allocated(self%last_grad)) then
            ! calculate gradient difference
            call self%transport(kappa, self%last_grad, error, self%settings%context)
            if (error /= 0) return
            grad_diff = grad - self%last_grad

            ! check if step addition should be skipped
            skip = self%skip_add(kappa, grad_diff, error)
            if (error /= 0) return
        end if

        ! the next pair starts at the current point, also if this one is skipped
        self%last_grad = grad
        if (skip) return

        ! check if point is added
        if (allocated(grad_diff) .and. self%settings%max_points > 0) then
            ! add new point
            if (self%n_points < self%settings%max_points) then
                self%n_points = self%n_points + 1
            ! replace old point
            else
                self%kappa_list(:, 1:self%n_points - 1) = &
                    self%kappa_list(:, 2:self%n_points)
                self%grad_diff_list(:, 1:self%n_points - 1) = &
                    self%grad_diff_list(:, 2:self%n_points)
            end if
            self%kappa_list(:, self%n_points) = kappa
            self%grad_diff_list(:, self%n_points) = grad_diff
            deallocate(grad_diff)
        end if

    end subroutine add_qn

    subroutine init_qn(self, n_param, settings)
        !
        ! this subroutine initializes an updating object with the given settings
        !
        class(updating_type), intent(inout) :: self
        integer(ip), intent(in) :: n_param
        type(qn_settings_type), intent(in) :: settings

        ! set the roles the object performs itself, the orbital update and the Hessian
        ! linear transformation, the settings and the number of parameters
        self%provides = provided_roles_type(update_orbs=.true., hess_x=.true.)
        self%error_origin = error_qn
        self%settings = settings
        self%n_param = n_param

        ! allocate arrays
        allocate(self%kappa_list(n_param, self%settings%max_points), &
                 self%grad_diff_list(n_param, self%settings%max_points))

    end subroutine init_qn

    function sr1_skip_add(self, kappa, grad_diff, error) result(skip)
        !
        ! this function implements the skipping criterion for the SR1 update
        !
        class(sr1_updating_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:), grad_diff(:)
        logical :: skip
        integer(ip), intent(out) :: error

        integer(ip) :: n_param
        real(rp) :: kappa_norm, residual_norm, curvature
        real(rp), allocatable :: residual(:)
        real(rp), parameter :: curvature_threshold = 1e-8_rp
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        skip = .false.
        n_param = size(kappa)
        kappa_norm = sqrt(ddot(n_param, kappa, 1_ip, kappa, 1_ip))
        residual = -kappa
        call self%init_hess(residual, error, self%settings%context)
        if (error /= 0) return
        residual = residual + grad_diff
        residual_norm = sqrt(ddot(n_param, residual, 1_ip, residual, 1_ip))
        curvature = abs(ddot(n_param, kappa, 1_ip, residual, 1_ip))
        if (curvature < curvature_threshold * kappa_norm * residual_norm) skip = .true.
        deallocate(residual)

    end function sr1_skip_add

    subroutine hess_x_sr1(self, x, hess_x, error, context)
        !
        ! this subroutine is the Hessian linear transformation of the SR1 updating
        ! scheme, which implements the L-SR1 version according to Byrd, R.H., Nocedal,
        ! J. & Schnabel, R.B. Mathematical Programming 63, 129–156 (1994).
        ! https://doi.org/10.1007/BF01582063
        !
        use opentrustregion, only: numerical_zero, verbosity_error

        class(sr1_updating_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: n_param, i, info
        real(rp), allocatable :: psi(:, :), U(:, :), q(:), z(:)
        integer(ip), allocatable :: ipiv(:)
        character(len=300) :: msg
        external :: dgemv, dgemm, dgetrf, dgetrs
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! get number of parameters
        n_param = size(x)

        ! check if trial vector vanishes
        if (abs(ddot(n_param, x, 1_ip, x, 1_ip)) < numerical_zero) then
            hess_x = 0.0_rp
            return
        end if

        ! initialize with Hessian diagonal
        hess_x = x
        call self%init_hess(hess_x, error, self%settings%context)
        if (error /= 0) return

        ! only current point is available
        if (self%n_points == 0) return

        ! build psi = grad_diff - h_start * kappa
        allocate(psi(n_param, self%n_points))
        do i = 1, self%n_points
            psi(:, i) = -self%kappa_list(:, i)
            call self%init_hess(psi(:, i), error, self%settings%context)
            if (error /= 0) return
            psi(:, i) = psi(:, i) + self%grad_diff_list(:, i)
        end do

        ! build U = psi^T * kappa
        allocate(U(self%n_points, self%n_points))
        call dgemm("T", "N", self%n_points, self%n_points, n_param, 1.0_rp, psi, &
                   n_param, self%kappa_list, n_param, 0.0_rp, U, self%n_points)

        ! symmetrize U
        U = 0.5_rp * (U + transpose(U))

        ! q = psi^T * x
        allocate(q(self%n_points))
        call dgemv("T", n_param, self%n_points, 1.0_rp, psi, n_param, x, 1_ip, 0.0_rp, &
                   q, 1_ip)

        ! LU factorize U
        allocate(ipiv(self%n_points))
        call dgetrf(self%n_points, self%n_points, U, self%n_points, ipiv, info)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Matrix inversion failed: Error in DGETRF, info = ", &
                info
            call self%settings%log(msg, verbosity_error, .true.)
            error = 1
            call add_error_origin(error, self%error_origin, self%settings)
            return
        end if

        ! solve U * z = q
        z = q
        deallocate(q)
        call dgetrs("N", self%n_points, 1_ip, U, self%n_points, ipiv, z, &
                    self%n_points, info)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Matrix solve failed: Error in DGETRS, info = ", info
            call self%settings%log(msg, verbosity_error, .true.)
            error = 1
            call add_error_origin(error, self%error_origin, self%settings)
            return
        end if
        deallocate(U, ipiv)

        ! w += psi * z
        call dgemv("N", n_param, self%n_points, 1.0_rp, psi, n_param, z, 1_ip, 1.0_rp, &
                   hess_x, 1_ip)

        ! deallocate arrays
        deallocate(psi, z)

    end subroutine hess_x_sr1

    function bfgs_skip_add(self, kappa, grad_diff, error) result(skip)
        !
        ! this function implements the skipping criterion for the BFGS update
        !
        class(bfgs_updating_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:), grad_diff(:)
        integer(ip), intent(out) :: error
        logical :: skip

        integer(ip) :: n_param
        real(rp) :: kappa_norm, grad_diff_norm, kappa_grad_dot
        real(rp), parameter :: curvature_threshold = 1e-14_rp
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        skip = .false.
        n_param = size(kappa)
        kappa_norm = sqrt(ddot(n_param, kappa, 1_ip, kappa, 1_ip))
        grad_diff_norm = sqrt(ddot(n_param, grad_diff, 1_ip, grad_diff, 1_ip))
        kappa_grad_dot = ddot(n_param, kappa, 1_ip, grad_diff, 1_ip)

        ! avoid zero near-zero curvature to avoid singular compact scaling matrix,
        ! negative curvature is not a problem since this is handled by the trust region
        ! solver
        if (abs(kappa_grad_dot) < curvature_threshold * kappa_norm * grad_diff_norm) &
            skip = .true.

    end function bfgs_skip_add

    subroutine hess_x_bfgs(self, x, hess_x, error, context)
        !
        ! this subroutine is the Hessian linear transformation of the BFGS updating
        ! scheme, which implements the L-BFGS version according to Byrd, R.H., Nocedal,
        ! J. & Schnabel, R.B. Mathematical Programming 63, 129–156 (1994).
        ! https://doi.org/10.1007/BF01582063
        !
        use opentrustregion, only: numerical_zero, verbosity_error

        class(bfgs_updating_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: n_param, i, j, info
        real(rp) :: kappa_grad_dot
        real(rp), allocatable :: compact_scal_mat(:, :), solution(:), init_hess_x(:)
        integer(ip), allocatable :: ipiv(:)
        character(len=300) :: msg
        external :: dgetrf, dgetrs
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(x)

        ! check if trial vector vanishes
        if (abs(ddot(n_param, x, 1_ip, x, 1_ip)) < numerical_zero) then
            hess_x = 0.0_rp
            return
        end if

        ! initialize with Hessian diagonal
        hess_x = x
        call self%init_hess(hess_x, error, self%settings%context)
        if (error /= 0) return

        ! only current point is available
        if (self%n_points == 0) return

        ! build compact scaling matrix
        ! [ kappa^T * h_start * kappa    L ]
        ! [             L^T            -D ]
        allocate(compact_scal_mat(2 * self%n_points, 2 * self%n_points))
        compact_scal_mat = 0.0_rp
        do j = 1, self%n_points
            init_hess_x = self%kappa_list(:, j)
            call self%init_hess(init_hess_x, error, self%settings%context)
            if (error /= 0) return
            do i = 1, self%n_points
                ! Weighted Top-Left: kappa_i^T * h_start * kappa_j
                compact_scal_mat(i, j) = ddot(n_param, self%kappa_list(:, i), 1_ip, &
                                              init_hess_x, 1_ip)

                ! L and D parts (do not depend on initial Hessian)
                kappa_grad_dot = ddot(n_param, self%kappa_list(:, i), 1_ip, &
                                      self%grad_diff_list(:, j), 1_ip)
                if (i > j) then
                    compact_scal_mat(i, j + self%n_points) = kappa_grad_dot
                    compact_scal_mat(j + self%n_points, i) = kappa_grad_dot
                else if (i == j) then
                    compact_scal_mat(i + self%n_points, j + self%n_points) = &
                        -kappa_grad_dot
                end if
            end do
        end do

        ! build right-hand side
        ! [ kappa^T * h_start * x ]
        ! [   grad_diff^T * x    ]
        allocate(solution(2 * self%n_points))
        do i = 1, self%n_points
            ! weighted dot product for the top half of RHS
            solution(i) = ddot(n_param, self%kappa_list(:, i), 1_ip, hess_x, 1_ip)
            solution(i + self%n_points) = ddot(n_param, self%grad_diff_list(:, i), &
                                               1_ip, x, 1_ip)
        end do

        ! LU factorize compact scaling matrix
        allocate(ipiv(2 * self%n_points))
        call dgetrf(2 * self%n_points, 2 * self%n_points, compact_scal_mat, &
                    2 * self%n_points, ipiv, info)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Matrix inversion failed: Error in DGETRF, info = ", &
                info
            call self%settings%log(msg, verbosity_error, .true.)
            error = 1
            call add_error_origin(error, self%error_origin, self%settings)
            return
        end if

        ! solve compact_scal_mat * solution = rhs
        call dgetrs("N", 2 * self%n_points, 1_ip, compact_scal_mat, 2 * self%n_points, &
                    ipiv, solution, 2 * self%n_points, info)
        deallocate(compact_scal_mat, ipiv)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Matrix inversion failed: Error in DGETRS, info = ", &
                info
            call self%settings%log(msg, verbosity_error, .true.)
            error = 1
            call add_error_origin(error, self%error_origin, self%settings)
            return
        end if

        ! construct Hessian linear transformation
        ! hess_x = hess_diag * x - ([hess_diag * kappa, grad_diff] * solution)
        do i = 1, self%n_points
            ! subtract rotation difference contributions
            init_hess_x = self%kappa_list(:, i)
            call self%init_hess(init_hess_x, error, self%settings%context)
            if (error /= 0) return
            hess_x = hess_x - (solution(i) * init_hess_x)
            ! subtract gradient difference contributions
            hess_x = hess_x - solution(i + self%n_points) * self%grad_diff_list(:, i)
        end do

        ! deallocate arrays
        deallocate(solution)

    end subroutine hess_x_bfgs

end module otr_qn
