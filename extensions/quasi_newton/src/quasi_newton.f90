! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn

    use opentrustregion, only: rp, ip, kw_len, settings_type, update_orbs_type, &
                               hess_x_type
    use otr_common, only: extension_type, provided_roles_type

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

    ! warning for settings that were not initialized
    character(len=*), parameter :: qn_settings_uninitialized_warning_msg = &
        "Quasi-Newton settings were not initialized. All quasi-Newton settings are "// &
        "set to default values"

    type, abstract, extends(extension_type) :: updating_type
        real(rp), allocatable :: kappa_list(:, :), grad_diff_list(:, :), last_grad(:)
        integer(ip) :: n_points = 0
        type(qn_settings_type) :: settings
    contains
        procedure :: update_orbs => update_orbs_qn
        procedure :: hess_x => hess_x_qn
        procedure :: init
        procedure :: add
        procedure(skip_add_interface), deferred :: skip_add
        procedure :: clear
    end type

    type, extends(updating_type) :: sr1_updating_type
    contains
        procedure :: skip_add => sr1_skip_add
    end type

    type, extends(updating_type) :: bfgs_updating_type
    contains
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

    ! global variables
    class(updating_type), pointer :: update_object => null()
    type(sr1_updating_type), target :: sr1_object
    type(bfgs_updating_type), target :: bfgs_object
    procedure(hess_x_type), pointer :: hess_x_qn_funptr => null()
    procedure(transport_type), pointer :: transport_funptr => null()
    procedure(init_hess_type), pointer :: init_hess_funptr => null()

contains

    subroutine qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                          transport_funptr_in, init_hess_funptr_in, &
                          update_orbs_qn_funptr, obj_func_qn_funptr, solver_settings, &
                          error, settings)
        !
        ! this subroutine returns the quasi-Newton orbital updating function, which
        ! wraps the host's, and the objective function, and wires the solver settings
        ! with the quasi-Newton object, through which the host's callback functions
        ! keep receiving the host context
        !
        use opentrustregion, only: solver_settings_type, obj_func_type, &
                                   verbosity_error, warn_uninitialized_settings
        use otr_common, only: update_orbs_callback, obj_func_callback

        integer(ip), intent(in) :: n_param
        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        procedure(transport_type), intent(in), pointer :: transport_funptr_in
        procedure(init_hess_type), intent(in), pointer :: init_hess_funptr_in
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_qn_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_qn_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(qn_settings_type), intent(in) :: settings

        type(qn_settings_type) :: local_settings
        class(updating_type), pointer :: object
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! initialize error flag and the returned functions
        error = 0
        update_orbs_qn_funptr => null()
        obj_func_qn_funptr => null()

        ! work on a copy of the settings, which is initialized if they were not, after
        ! printing a warning regardless of the verbosity, since the default settings
        ! would silence it
        local_settings = settings
        if (.not. local_settings%initialized) then
            call warn_uninitialized_settings(local_settings, &
                                             qn_settings_uninitialized_warning_msg)
            call local_settings%init(error)
            if (error /= 0) return
        end if

        ! get object for quasi-Newton and corresponding Hessian linear transformation
        ! function, which the module points to only once the object wired the solver
        ! settings, so that refused settings leave the object in use unchanged
        if (local_settings%hess_update_scheme == "sr1") then
            object => sr1_object
            hess_x_funptr => sr1_hess_x_fun
        else if (local_settings%hess_update_scheme == "bfgs") then
            object => bfgs_object
            hess_x_funptr => bfgs_hess_x_fun
        else
            call local_settings%log("Quasi-Newton updating scheme not implemented.", &
                                    verbosity_error, .true.)
            error = 1
            return
        end if

        ! refuse solver settings the object already wired before changing any state
        call object%check_unwired(solver_settings, local_settings, error)
        if (error /= 0) return

        ! initialize the object
        call object%init(n_param, local_settings)

        ! wire the solver settings with the object, which saves the host's orbital
        ! update, objective function, callback functions and host contexts
        call object%wire_solver_settings(solver_settings, local_settings, error, &
                                         update_orbs_funptr, obj_func_funptr)
        if (error /= 0) return
        update_object => object
        hess_x_qn_funptr => hess_x_funptr

        ! set pointer to transport function
        transport_funptr => transport_funptr_in

        ! set pointer to Hessian initialization function
        init_hess_funptr => init_hess_funptr_in

        ! get the orbital updating function and the objective function of the object
        update_orbs_qn_funptr => update_orbs_callback
        obj_func_qn_funptr => obj_func_callback

    end subroutine qn_factory

    subroutine update_orbs_qn(self, kappa, func, grad, h_diag, hess_x_funptr, error, &
                              context)
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

        ! add new step
        if (sum(abs(kappa)) > 0.0_rp) then
            call self%add(kappa, grad, error)
            if (error /= 0) return
        end if

    end subroutine update_orbs_qn

    subroutine hess_x_qn(self, x, hess_x, error, context)
        !
        ! this subroutine is the quasi-Newton Hessian linear transformation of the
        ! updating scheme of the object
        !
        class(updating_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        class(*), pointer :: no_context

        ! the Hessian linear transformation of the updating scheme reads the history of
        ! the scheme's object and needs no context
        no_context => null()
        call hess_x_qn_funptr(x, hess_x, error, no_context)

    end subroutine hess_x_qn

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

    subroutine update_orbs_qn_deconstructor()
        !
        ! this subroutine deallocates the quasi-Newton objects
        !
        if (associated(update_object)) then
            call update_object%clear()
            nullify(update_object)
        end if

    end subroutine update_orbs_qn_deconstructor

    subroutine add(self, kappa, grad, error)
        !
        ! this subroutine adds a new rotation and gradient to the corresponding lists
        ! and transports the history accordingly
        !
        class(updating_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), grad(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: grad_diff(:)
        integer(ip) :: n_param, i

        ! initialize error flag
        error = 0

        ! get number of parameters
        n_param = size(kappa)

        ! transport history
        do i = 1, self%n_points
            call transport_funptr(kappa, self%kappa_list(:, i), error, &
                                  self%settings%context)
            if (error /= 0) return
            call transport_funptr(kappa, self%grad_diff_list(:, i), error, &
                                  self%settings%context)
            if (error /= 0) return
        end do

        ! check if last gradient is available to calculate gradient difference
        if (allocated(self%last_grad)) then
            ! calculate gradient difference
            call transport_funptr(kappa, self%last_grad, error, self%settings%context)
            if (error /= 0) return
            grad_diff = grad - self%last_grad

            ! check if step addition should be skipped
            if (self%skip_add(kappa, grad_diff, error)) return
            if (error /= 0) return
        end if
        self%last_grad = grad

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

    end subroutine add

    subroutine init(self, n_param, settings)
        !
        ! this subroutine initializes an updating object with the given settings
        !
        class(updating_type), intent(inout) :: self
        integer(ip), intent(in) :: n_param
        type(qn_settings_type), intent(in) :: settings

        ! start from an empty object so that nothing of an earlier setup survives
        call self%clear()

        ! set the roles the object performs itself, the orbital update and the Hessian
        ! linear transformation, the settings and the number of parameters
        self%provides = provided_roles_type(update_orbs=.true., hess_x=.true.)
        self%settings = settings
        self%n_param = n_param

        ! allocate arrays
        allocate(self%kappa_list(n_param, self%settings%max_points), &
                 self%grad_diff_list(n_param, self%settings%max_points))

    end subroutine init

    subroutine clear(self)
        !
        ! this subroutine deallocates an updating object
        !
        class(updating_type), intent(inout) :: self

        if (allocated(self%kappa_list)) deallocate(self%kappa_list)
        if (allocated(self%grad_diff_list)) deallocate(self%grad_diff_list)
        if (allocated(self%last_grad)) deallocate(self%last_grad)
        self%n_points = 0

    end subroutine clear

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
        call init_hess_funptr(residual, error, self%settings%context)
        if (error /= 0) return
        residual = residual + grad_diff
        residual_norm = sqrt(ddot(n_param, residual, 1_ip, residual, 1_ip))
        curvature = abs(ddot(n_param, kappa, 1_ip, residual, 1_ip))
        if (curvature < curvature_threshold * kappa_norm * residual_norm) skip = .true.
        deallocate(residual)

    end function sr1_skip_add

    subroutine sr1_hess_x_fun(x, hess_x, error, context)
        !
        ! this implements the L-SR1 version according to Byrd, R.H., Nocedal, J. &
        ! Schnabel, R.B. Mathematical Programming 63, 129–156 (1994).
        ! https://doi.org/10.1007/BF01582063
        !
        use opentrustregion, only: numerical_zero, verbosity_error

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
        call init_hess_funptr(hess_x, error, update_object%settings%context)
        if (error /= 0) return

        ! only current point is available
        if (sr1_object%n_points == 0) return

        ! build psi = grad_diff - h_start * kappa
        allocate(psi(n_param, sr1_object%n_points))
        do i = 1, sr1_object%n_points
            psi(:, i) = -sr1_object%kappa_list(:, i)
            call init_hess_funptr(psi(:, i), error, update_object%settings%context)
            if (error /= 0) return
            psi(:, i) = psi(:, i) + sr1_object%grad_diff_list(:, i)
        end do

        ! build U = psi^T * kappa
        allocate(U(sr1_object%n_points, sr1_object%n_points))
        call dgemm("T", "N", sr1_object%n_points, sr1_object%n_points, n_param, &
                   1.0_rp, psi, n_param, sr1_object%kappa_list, n_param, 0.0_rp, U, &
                   sr1_object%n_points)

        ! symmetrize U
        U = 0.5_rp * (U + transpose(U))

        ! q = psi^T * x
        allocate(q(sr1_object%n_points))
        call dgemv("T", n_param, sr1_object%n_points, 1.0_rp, psi, n_param, x, 1_ip, &
                   0.0_rp, q, 1_ip)

        ! LU factorize U
        allocate(ipiv(sr1_object%n_points))
        call dgetrf(sr1_object%n_points, sr1_object%n_points, U, sr1_object%n_points, &
                    ipiv, info)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Matrix inversion failed: Error in DGETRF, info = ", &
                info
            call sr1_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! solve U * z = q
        z = q
        deallocate(q)
        call dgetrs("N", sr1_object%n_points, 1_ip, U, sr1_object%n_points, ipiv, z, &
                    sr1_object%n_points, info)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Matrix solve failed: Error in DGETRS, info = ", info
            call sr1_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if
        deallocate(U, ipiv)

        ! w += psi * z
        call dgemv("N", n_param, sr1_object%n_points, 1.0_rp, psi, n_param, z, 1_ip, &
                   1.0_rp, hess_x, 1_ip)

        ! deallocate arrays
        deallocate(psi, z)

    end subroutine sr1_hess_x_fun

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

    subroutine bfgs_hess_x_fun(x, hess_x, error, context)
        !
        ! this implements the L-BFGS version according to Byrd, R.H., Nocedal, J. &
        ! Schnabel, R.B. Mathematical Programming 63, 129–156 (1994).
        ! https://doi.org/10.1007/BF01582063
        !
        use opentrustregion, only: numerical_zero, verbosity_error

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
        call init_hess_funptr(hess_x, error, update_object%settings%context)
        if (error /= 0) return

        ! only current point is available
        if (bfgs_object%n_points == 0) return

        ! build compact scaling matrix
        ! [ kappa^T * h_start * kappa    L ]
        ! [             L^T            -D ]
        allocate(compact_scal_mat(2 * bfgs_object%n_points, 2 * bfgs_object%n_points))
        compact_scal_mat = 0.0_rp
        do j = 1, bfgs_object%n_points
            init_hess_x = bfgs_object%kappa_list(:, j)
            call init_hess_funptr(init_hess_x, error, update_object%settings%context)
            if (error /= 0) return
            do i = 1, bfgs_object%n_points
                ! Weighted Top-Left: kappa_i^T * h_start * kappa_j
                compact_scal_mat(i, j) = ddot(n_param, bfgs_object%kappa_list(:, i), &
                                              1_ip, init_hess_x, 1_ip)

                ! L and D parts (do not depend on initial Hessian)
                kappa_grad_dot = ddot(n_param, bfgs_object%kappa_list(:, i), 1_ip, &
                                      bfgs_object%grad_diff_list(:, j), 1_ip)
                if (i > j) then
                    compact_scal_mat(i, j + bfgs_object%n_points) = kappa_grad_dot
                    compact_scal_mat(j + bfgs_object%n_points, i) = kappa_grad_dot
                else if (i == j) then
                    compact_scal_mat(i + bfgs_object%n_points, &
                                     j + bfgs_object%n_points) = -kappa_grad_dot
                end if
            end do
        end do

        ! build right-hand side
        ! [ kappa^T * h_start * x ]
        ! [   grad_diff^T * x    ]
        allocate(solution(2 * bfgs_object%n_points))
        do i = 1, bfgs_object%n_points
            ! weighted dot product for the top half of RHS
            solution(i) = ddot(n_param, bfgs_object%kappa_list(:, i), 1_ip, hess_x, &
                               1_ip)
            solution(i + bfgs_object%n_points) = &
                ddot(n_param, bfgs_object%grad_diff_list(:, i), 1_ip, x, 1_ip)
        end do

        ! LU factorize compact scaling matrix
        allocate(ipiv(2 * bfgs_object%n_points))
        call dgetrf(2 * bfgs_object%n_points, 2 * bfgs_object%n_points, &
                    compact_scal_mat, 2 * bfgs_object%n_points, ipiv, info)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Matrix inversion failed: Error in DGETRF, info = ", &
                info
            call bfgs_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! solve compact_scal_mat * solution = rhs
        call dgetrs("N", 2 * bfgs_object%n_points, 1_ip, compact_scal_mat, &
                    2 * bfgs_object%n_points, ipiv, solution, &
                    2 * bfgs_object%n_points, info)
        deallocate(compact_scal_mat, ipiv)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Matrix inversion failed: Error in DGETRS, info = ", &
                info
            call bfgs_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! construct Hessian linear transformation
        ! hess_x = hess_diag * x - ([hess_diag * kappa, grad_diff] * solution)
        do i = 1, bfgs_object%n_points
            ! subtract rotation difference contributions
            init_hess_x = bfgs_object%kappa_list(:, i)
            call init_hess_funptr(init_hess_x, error, update_object%settings%context)
            if (error /= 0) return
            hess_x = hess_x - (solution(i) * init_hess_x)
            ! subtract gradient difference contributions
            hess_x = hess_x - solution(i + bfgs_object%n_points) * &
                     bfgs_object%grad_diff_list(:, i)
        end do

        ! deallocate arrays
        deallocate(solution)

    end subroutine bfgs_hess_x_fun

end module otr_qn
