! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_rotation_qn

    use opentrustregion_core, only: rp, ip, kw_len, settings_type, update_orbs_type, &
                                    hess_x_type, precond_type, add_error_origin
    use otr_common, only: extension_type, extension_level_type, provided_roles_type

    implicit none

    ! interfaces for callback functions
    abstract interface
        function obj_func_grad_type(kappa, grad, error, context) result(func)
            import :: rp, ip

            real(rp), intent(in), target :: kappa(:)
            real(rp), intent(out), target :: grad(:)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
            real(rp) :: func

        end function obj_func_grad_type
    end interface

    abstract interface
        subroutine transport_type(geodesic, tangent_vector, error, context)
            import :: rp, ip

            real(rp), intent(in), target :: geodesic(:)
            real(rp), intent(inout), target :: tangent_vector(:)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context

        end subroutine transport_type
    end interface

    ! derived type for settings, whose optional callback functions are the initial
    ! Hessian, the inverse of the level-shifted initial Hessian, the objective
    ! function that also returns the gradient at the trial point and the transport of
    ! tangent vectors along a geodesic
    type, extends(settings_type) :: rotation_qn_settings_type
        procedure(hess_x_type), pointer, nopass :: init_hess => null()
        procedure(precond_type), pointer, nopass :: init_precond => null()
        procedure(obj_func_grad_type), pointer, nopass :: obj_func_grad => null()
        procedure(transport_type), pointer, nopass :: transport => null()
        character(len=kw_len) :: qn_method = "sr1"
        integer(ip) :: max_points = 30
    end type

    type(rotation_qn_settings_type), parameter :: default_rotation_qn_settings = &
        rotation_qn_settings_type()

    ! updating schemes
    character(len=kw_len), parameter :: qn_methods(2) = &
        [character(len=kw_len) :: "sr1", "bfgs"]

    ! warning for settings that were not initialized
    character(len=*), parameter :: rotation_qn_settings_uninitialized_warning_msg = &
        "Rotation QN settings were not initialized. All rotation QN settings are "// &
        "set to default values"

    ! error origin of the rotation QN extension, with which its factory and object tag
    ! the errors they detect themselves
    integer(ip), parameter :: error_rotation_qn = 5400

    ! relative thresholds below which the updating schemes skip a pair of the history:
    ! the SR1 denominator s^T (y - B s) relative to |s| |y - B s| (Nocedal & Wright,
    ! Numerical Optimization, Eq. 6.26) and the BFGS curvatures s^T y relative to
    ! |s| |y| and s^T B s relative to |s| |B s|
    real(rp), parameter :: sr1_skip_thresh = 1e-8_rp, bfgs_skip_thresh = 1e-8_rp

    ! rotation QN object, holding the history of steps and gradient differences, the
    ! settings, the Hessian diagonal of the current point and the model the updating
    ! scheme builds from the history, the rank-1 terms by which it updates the initial
    ! Hessian, the initial Hessian of the settings or, without one, the Hessian
    ! diagonal:
    ! B x = B_0 x + sum_k coeffs(k) vectors(:, k) vectors(:, k)^T x
    type, abstract, extends(extension_type) :: updating_type
        real(rp), allocatable :: kappa_list(:, :), grad_diff_list(:, :), last_grad(:), &
                                 h_diag(:), vectors(:, :), coeffs(:)
        integer(ip) :: n_points = 0, n_terms = 0
        type(rotation_qn_settings_type) :: settings
    contains
        procedure :: update_orbs => update_orbs_rotation_qn
        procedure :: obj_func => obj_func_rotation_qn
        procedure :: hess_x => hess_x_rotation_qn
        procedure :: precond => precond_rotation_qn
        procedure :: init => init_rotation_qn
        procedure :: add => add_rotation_qn
        procedure :: add_pair => add_pair_rotation_qn
        procedure :: model_x => model_x_rotation_qn
        procedure :: init_precond => init_precond_rotation_qn
        procedure(build_model_interface), deferred :: build_model
    end type

    type, extends(updating_type) :: sr1_updating_type
    contains
        procedure :: build_model => build_model_sr1
    end type

    type, extends(updating_type) :: bfgs_updating_type
    contains
        procedure :: build_model => build_model_bfgs
    end type

    abstract interface
        subroutine build_model_interface(self, error)
            import updating_type, ip

            class(updating_type), intent(inout) :: self
            integer(ip), intent(out) :: error
        end subroutine build_model_interface
    end interface

contains

    subroutine rotation_qn_factory( &
        n_param, update_orbs_funptr, obj_func_funptr, update_orbs_rotation_qn_funptr, &
        obj_func_rotation_qn_funptr, solver_settings, solver_state, error, settings)
        !
        ! this subroutine returns the rotation QN orbital updating function, which
        ! wraps the host's, and the objective function, and wires the solver settings
        ! with a new rotation QN object, which they own and through which the host's
        ! callback functions keep receiving the host context, and which reads the state
        ! of the solver from the record the host also passes to the solver; the model
        ! of the object changes after a rejected step whose gradient the objective
        ! function of the settings returns, which the solver therefore refreshes
        !
        use opentrustregion_core, only: obj_func_type, solver_state_type
        use otr_common, only: extension_solver_settings_type, update_orbs_callback, &
                              obj_func_callback

        integer(ip), intent(in) :: n_param
        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        procedure(update_orbs_type), intent(out), pointer :: &
            update_orbs_rotation_qn_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_rotation_qn_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(rotation_qn_settings_type), intent(in), optional :: settings

        type(rotation_qn_settings_type) :: local_settings
        class(extension_type), allocatable, target :: object

        ! initialize error flag and the returned functions
        error = 0
        update_orbs_rotation_qn_funptr => null()
        obj_func_rotation_qn_funptr => null()

        ! work on a copy of the given settings or, without any, of the default settings
        if (present(settings)) then
            local_settings = settings
        else
            local_settings = default_rotation_qn_settings
        end if

        ! perform sanity check, which also converts the updating scheme to lowercase
        call rotation_qn_sanity_check(local_settings, n_param, error)
        call add_error_origin(error, error_rotation_qn, local_settings)
        if (error /= 0) return

        ! create the object of the updating scheme
        if (local_settings%qn_method == "sr1") then
            allocate(sr1_updating_type :: object)
        else
            allocate(bfgs_updating_type :: object)
        end if

        ! initialize the object
        select type (object)
        class is (updating_type)
            call object%init(n_param, local_settings)
        end select

        ! wire the solver settings with the object, which reads the record of the state
        ! and saves the host's orbital update, objective function, callback functions
        ! and host contexts, and let the settings own it
        call solver_settings%wire(object, solver_state, local_settings, error, &
                                  update_orbs_funptr, obj_func_funptr)
        call add_error_origin(error, error_rotation_qn, local_settings)
        if (error /= 0) return

        ! the model changes after a rejected step if the objective function returns
        ! the gradient at the trial point
        if (associated(local_settings%obj_func_grad)) &
            solver_settings%refresh_hess = .true.

        ! get the orbital updating function and the objective function of the object
        update_orbs_rotation_qn_funptr => update_orbs_callback
        obj_func_rotation_qn_funptr => obj_func_callback

    end subroutine rotation_qn_factory

    subroutine rotation_qn_sanity_check(settings, n_param, error)
        !
        ! this subroutine performs a sanity check for rotation QN input parameters,
        ! after converting the updating scheme to lowercase
        !
        use opentrustregion_core, only: verbosity_error, string_to_lowercase

        type(rotation_qn_settings_type), intent(inout) :: settings
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
        settings%qn_method = string_to_lowercase(settings%qn_method)

        ! check for character options
        if (.not. any(settings%qn_method == qn_methods)) then
            call settings%log("QN method option unknown. Possible values are "// &
                              """sr1"" (symmetric rank-1 update) and ""bfgs"" "// &
                              "(Broyden-Fletcher-Goldfarb-Shanno update).", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine rotation_qn_sanity_check

    subroutine invoke_init_hess(init_hess, x, init_hess_x, settings, error_origin, &
                                error)
        !
        ! this subroutine applies the initial Hessian of the settings with their host
        ! context and passes on its error, and returns an error of its own with the
        ! origin of the object if it returns values that are not finite
        !
        use opentrustregion_core, only: verbosity_error
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        procedure(hess_x_type), intent(in), pointer :: init_hess
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: init_hess_x(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(in) :: error_origin
        integer(ip), intent(out) :: error

        ! apply the initial Hessian
        call init_hess(x, init_hess_x, error, settings%context)
        if (error /= 0) return

        ! check the result
        if (.not. all(ieee_is_finite(init_hess_x))) then
            call settings%log("Initial Hessian returned values that are not finite.", &
                              verbosity_error, .true.)
            error = 1
            call add_error_origin(error, error_origin, settings)
        end if

    end subroutine invoke_init_hess

    subroutine invoke_init_precond(init_precond, x, mu, init_precond_x, settings, &
                                   error_origin, error)
        !
        ! this subroutine applies the preconditioner of the initial Hessian of the
        ! settings, shifted by a level shift, with their host context and passes on its
        ! error, and returns an error of its own with the origin of the object if it
        ! returns values that are not finite
        !
        use opentrustregion_core, only: verbosity_error
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        procedure(precond_type), intent(in), pointer :: init_precond
        real(rp), intent(in), target :: x(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: init_precond_x(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(in) :: error_origin
        integer(ip), intent(out) :: error

        ! apply the preconditioner of the initial Hessian
        call init_precond(x, mu, init_precond_x, error, settings%context)
        if (error /= 0) return

        ! check the result
        if (.not. all(ieee_is_finite(init_precond_x))) then
            call settings%log("Preconditioner of the initial Hessian returned "// &
                              "values that are not finite.", verbosity_error, .true.)
            error = 1
            call add_error_origin(error, error_origin, settings)
        end if

    end subroutine invoke_init_precond

    subroutine invoke_obj_func_grad(obj_func_grad, kappa, func, grad, settings, &
                                    error_origin, error)
        !
        ! this subroutine calls the objective function of the settings, which also
        ! returns the gradient at the trial point, with their host context and passes
        ! on its error; an energy that is not a number counts as infinite and leaves
        ! the gradient unchecked, as it does for an infinite energy, which rejects the
        ! step, while it returns an error of its own with the origin of the object for
        ! a finite energy whose gradient is not finite
        !
        use opentrustregion_core, only: verbosity_error
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, &
                                                 ieee_value, ieee_positive_inf

        procedure(obj_func_grad_type), intent(in), pointer :: obj_func_grad
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(in) :: error_origin
        integer(ip), intent(out) :: error

        ! call the objective function
        func = obj_func_grad(kappa, grad, error, settings%context)
        if (error /= 0) return

        ! an energy that is not a number counts as infinite, and the gradient of an
        ! infinite energy is not used
        if (ieee_is_nan(func)) func = ieee_value(func, ieee_positive_inf)
        if (.not. ieee_is_finite(func)) return

        ! check the gradient
        if (.not. all(ieee_is_finite(grad))) then
            call settings%log("Objective function returned a gradient that is not "// &
                              "finite.", verbosity_error, .true.)
            error = 1
            call add_error_origin(error, error_origin, settings)
        end if

    end subroutine invoke_obj_func_grad

    subroutine invoke_transport(transport, geodesic, tangent_vector, settings, &
                                error_origin, error)
        !
        ! this subroutine transports a tangent vector along a geodesic with the
        ! transport function of the settings and their host context and passes on its
        ! error, and returns an error of its own with the origin of the object if the
        ! function returns values that are not finite, while without a transport
        ! function it leaves the tangent vector unchanged, which transports it exactly
        ! in a flat space and in the frame that moves with the geodesic on the
        ! Grassmann manifold
        !
        use opentrustregion_core, only: verbosity_error
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        procedure(transport_type), intent(in), pointer :: transport
        real(rp), intent(in), target :: geodesic(:)
        real(rp), intent(inout), target :: tangent_vector(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(in) :: error_origin
        integer(ip), intent(out) :: error

        ! leave the tangent vector unchanged without a transport function
        error = 0
        if (.not. associated(transport)) return

        ! call the transport function
        call transport(geodesic, tangent_vector, error, settings%context)
        if (error /= 0) return

        ! check the result
        if (.not. all(ieee_is_finite(tangent_vector))) then
            call settings%log("Transport function returned values that are not "// &
                              "finite.", verbosity_error, .true.)
            error = 1
            call add_error_origin(error, error_origin, settings)
        end if

    end subroutine invoke_transport

    recursive subroutine update_orbs_rotation_qn(self, kappa, func, grad, h_diag, &
                                                 hess_x_funptr, error, context)
        !
        ! this subroutine is the rotation QN orbital updating function, which updates
        ! the orbitals with the host's orbital update and host context, keeps the
        ! Hessian diagonal, the initial Hessian and its preconditioner without those of
        ! the settings, returns the quasi-Newton Hessian linear transformation in place
        ! of the host's, which is not used, and adds the new point to the history
        !
        use otr_common, only: hess_x_callback

        class(updating_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        procedure(hess_x_type), pointer :: host_hess_x_funptr

        ! initialize error flag
        error = 0

        ! update orbitals
        host_hess_x_funptr => null()
        call self%solver_level%host%update_orbs(kappa, func, grad, h_diag, &
                                                host_hess_x_funptr, error, &
                                                self%solver_level%host%context)
        if (error /= 0) return

        ! keep the Hessian diagonal
        self%h_diag = h_diag

        ! get pointer to the quasi-Newton Hessian linear transformation of the object
        hess_x_funptr => hess_x_callback

        ! add the new point, also the starting point the first call reaches without a
        ! step
        if (sum(abs(kappa)) > 0.0_rp .or. .not. allocated(self%last_grad)) then
            call self%add(kappa, grad, error)
            if (error /= 0) return
        end if

    end subroutine update_orbs_rotation_qn

    recursive function obj_func_rotation_qn(self, kappa, error, context) result(func)
        !
        ! this function is the rotation QN objective function, which evaluates the
        ! objective function of the settings, if they have one, and adds the pair of
        ! the trial step and the gradient difference at the trial point, transported
        ! back to the current point, to the history, whose model it rebuilds, so that
        ! every evaluated trial step, also a rejected one, adds its pair, and otherwise
        ! evaluates the host's objective function with the host context of the level
        !
        use otr_common, only: obj_func_extension
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        class(updating_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        real(rp), allocatable :: grad(:), back(:)

        ! evaluate the host's objective function without one of the settings
        if (.not. associated(self%settings%obj_func_grad)) then
            func = obj_func_extension(self, kappa, error, context)
            return
        end if

        ! evaluate the energy and the gradient at the trial point
        allocate(grad(self%n_param))
        call invoke_obj_func_grad(self%settings%obj_func_grad, kappa, func, grad, &
                                  self%settings, self%error_origin, error)
        if (error /= 0) return

        ! a pair needs a finite energy, a step, the gradient of the current point and
        ! room in the history
        if (.not. ieee_is_finite(func) .or. .not. sum(abs(kappa)) > 0.0_rp .or. &
            .not. allocated(self%last_grad) .or. self%settings%max_points == 0) return

        ! transport the gradient back to the current point along the reversed step
        back = -kappa
        call invoke_transport(self%settings%transport, back, grad, self%settings, &
                              self%error_origin, error)
        if (error /= 0) return

        ! add the pair and rebuild the model
        call self%add_pair(kappa, grad - self%last_grad)
        call self%build_model(error)

    end function obj_func_rotation_qn

    subroutine add_rotation_qn(self, kappa, grad, error)
        !
        ! this subroutine transports the history to the new point, adds the pair of the
        ! step and the gradient difference to the point before, unless the settings
        ! have an objective function, which added the pairs of every trial step the
        ! solver evaluated, which leaves out a step the solver takes without evaluating
        ! it, as from a saddle point, and builds the model of the new point; the step
        ! needs no transport, since a geodesic transports its own direction to itself
        !
        class(updating_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), grad(:)
        integer(ip), intent(out) :: error

        integer(ip) :: i

        ! initialize error flag
        error = 0

        ! a history needs the gradient of the point before and room
        if (allocated(self%last_grad) .and. self%settings%max_points > 0) then
            ! transport history
            do i = 1, self%n_points
                call invoke_transport(self%settings%transport, kappa, &
                                      self%kappa_list(:, i), self%settings, &
                                      self%error_origin, error)
                if (error /= 0) return
                call invoke_transport(self%settings%transport, kappa, &
                                      self%grad_diff_list(:, i), self%settings, &
                                      self%error_origin, error)
                if (error /= 0) return
            end do

            ! add the pair of the step with the last gradient transported to the new
            ! point, unless the objective function of the settings added it
            if (.not. associated(self%settings%obj_func_grad)) then
                call invoke_transport(self%settings%transport, kappa, self%last_grad, &
                                      self%settings, self%error_origin, error)
                if (error /= 0) return
                call self%add_pair(kappa, grad - self%last_grad)
            end if
        end if

        ! the next pair starts at the new point
        self%last_grad = grad

        ! build the model of the new point, whose initial Hessian has changed
        call self%build_model(error)

    end subroutine add_rotation_qn

    subroutine add_pair_rotation_qn(self, kappa, grad_diff)
        !
        ! this subroutine appends a pair of a step and a gradient difference to the
        ! history, dropping the oldest pair once the history is full
        !
        class(updating_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), grad_diff(:)

        ! drop the oldest pair of a full history
        if (self%n_points == self%settings%max_points) then
            self%kappa_list(:, :self%n_points - 1) = self%kappa_list(:, 2:self%n_points)
            self%grad_diff_list(:, :self%n_points - 1) = &
                self%grad_diff_list(:, 2:self%n_points)
            self%n_points = self%n_points - 1
        end if

        ! append the new pair
        self%n_points = self%n_points + 1
        self%kappa_list(:, self%n_points) = kappa
        self%grad_diff_list(:, self%n_points) = grad_diff

    end subroutine add_pair_rotation_qn

    subroutine init_rotation_qn(self, n_param, settings)
        !
        ! this subroutine initializes an updating object with the given settings
        !
        class(updating_type), intent(inout) :: self
        integer(ip), intent(in) :: n_param
        type(rotation_qn_settings_type), intent(in) :: settings

        ! set the roles the object performs itself, the orbital update, the objective
        ! function, the Hessian linear transformation and the preconditioner, the
        ! settings and the number of parameters
        self%provides = provided_roles_type(update_orbs=.true., obj_func=.true., &
                                            hess_x=.true., precond=.true.)
        self%error_origin = error_rotation_qn
        self%settings = settings
        self%n_param = n_param

        ! allocate the history and the model, whose BFGS scheme adds two terms per
        ! pair
        allocate(self%kappa_list(n_param, self%settings%max_points), &
                 self%grad_diff_list(n_param, self%settings%max_points), &
                 self%vectors(n_param, 2 * self%settings%max_points), &
                 self%coeffs(2 * self%settings%max_points))

    end subroutine init_rotation_qn

    recursive subroutine model_x_rotation_qn(self, x, model_x, n_terms, error)
        !
        ! this subroutine applies the initial Hessian, updated by the first n_terms
        ! rank-1 terms of the model, to a vector
        !
        class(updating_type), intent(in) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: model_x(:)
        integer(ip), intent(in) :: n_terms
        integer(ip), intent(out) :: error

        real(rp) :: projections(n_terms)
        external :: dgemv

        ! apply the initial Hessian of the settings or the Hessian diagonal
        if (associated(self%settings%init_hess)) then
            call invoke_init_hess(self%settings%init_hess, x, model_x, self%settings, &
                                  self%error_origin, error)
            if (error /= 0) return
        else
            error = 0
            model_x = self%h_diag * x
        end if
        if (n_terms == 0) return

        ! add coeffs(k) vectors(:, k) vectors(:, k)^T x for every term
        call dgemv("T", self%n_param, n_terms, 1.0_rp, self%vectors, self%n_param, x, &
                   1_ip, 0.0_rp, projections, 1_ip)
        projections = self%coeffs(:n_terms) * projections
        call dgemv("N", self%n_param, n_terms, 1.0_rp, self%vectors, self%n_param, &
                   projections, 1_ip, 1.0_rp, model_x, 1_ip)

    end subroutine model_x_rotation_qn

    recursive subroutine hess_x_rotation_qn(self, x, hess_x, error, context)
        !
        ! this subroutine is the Hessian linear transformation of the rotation QN
        ! object, which applies its model
        !
        class(updating_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call self%model_x(x, hess_x, self%n_terms, error)

    end subroutine hess_x_rotation_qn

    recursive subroutine init_precond_rotation_qn(self, x, mu, init_precond_x, error, &
                                                  context)
        !
        ! this subroutine applies the inverse of the level-shifted initial Hessian, the
        ! preconditioner of the settings or, without one, the level-shifted diagonal
        ! preconditioner of the Hessian diagonal the solver applies without one,
        ! followed by the host's projection of the level it receives as context, and
        ! returns an error with the origin of the object if the projection returns
        ! values that are not finite
        !
        use opentrustregion_core, only: level_shifted_divisors, verbosity_error
        use otr_common, only: extension_level
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        class(updating_type), intent(in) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: init_precond_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        ! initialize error flag
        error = 0

        if (associated(self%settings%init_precond)) then
            call invoke_init_precond(self%settings%init_precond, x, mu, &
                                     init_precond_x, self%settings, self%error_origin, &
                                     error)
            return
        end if

        ! the floored diagonal keeps the result finite, unless the projection breaks it
        init_precond_x = x / level_shifted_divisors(self%h_diag, mu)
        level => extension_level(context)
        if (associated(level%host%project)) then
            call level%host%project(init_precond_x, error, level%host%context)
            if (error /= 0) return
            if (.not. all(ieee_is_finite(init_precond_x))) then
                call self%settings%log("Projection returned values that are not "// &
                                       "finite.", verbosity_error, .true.)
                error = 1
                call add_error_origin(error, self%error_origin, self%settings)
            end if
        end if

    end subroutine init_precond_rotation_qn

    recursive subroutine precond_rotation_qn(self, residual, mu, precond_residual, &
                                             error, context)
        !
        ! this subroutine is the preconditioner of the rotation QN object, which
        ! applies the inverse of its level-shifted model by the
        ! Sherman-Morrison-Woodbury formula,
        ! (A + V C V^T)^-1 = A^-1 - A^-1 V (C^-1 + V^T A^-1 V)^-1 V^T A^-1, with the
        ! inverse A^-1 of the level-shifted initial Hessian, or only that inverse where
        ! the level shift makes the model singular
        !
        class(updating_type), intent(inout) :: self
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp), allocatable :: init_precond_vectors(:, :), capacitance(:, :), &
                                 projections(:)
        integer(ip), allocatable :: pivots(:)
        integer(ip) :: k, info
        external :: dgemm, dgemv, dgesv

        ! apply the inverse of the level-shifted initial Hessian to the residual
        call self%init_precond(residual, mu, precond_residual, error, context)
        if (error /= 0 .or. self%n_terms == 0) return

        ! apply it to the vectors of the terms
        allocate(init_precond_vectors(self%n_param, self%n_terms))
        do k = 1, self%n_terms
            call self%init_precond(self%vectors(:, k), mu, init_precond_vectors(:, k), &
                                   error, context)
            if (error /= 0) return
        end do

        ! capacitance matrix C^-1 + V^T A^-1 V and projections V^T A^-1 r
        allocate(capacitance(self%n_terms, self%n_terms), projections(self%n_terms), &
                 pivots(self%n_terms))
        call dgemm("T", "N", self%n_terms, self%n_terms, self%n_param, 1.0_rp, &
                   self%vectors, self%n_param, init_precond_vectors, self%n_param, &
                   0.0_rp, capacitance, self%n_terms)
        do k = 1, self%n_terms
            capacitance(k, k) = capacitance(k, k) + 1.0_rp / self%coeffs(k)
        end do
        call dgemv("T", self%n_param, self%n_terms, 1.0_rp, self%vectors, &
                   self%n_param, precond_residual, 1_ip, 0.0_rp, projections, 1_ip)

        ! subtract A^-1 V (C^-1 + V^T A^-1 V)^-1 V^T A^-1 r, unless the level shift
        ! makes the model singular
        call dgesv(self%n_terms, 1_ip, capacitance, self%n_terms, pivots, projections, &
                   self%n_terms, info)
        if (info /= 0) return
        call dgemv("N", self%n_param, self%n_terms, -1.0_rp, init_precond_vectors, &
                   self%n_param, projections, 1_ip, 1.0_rp, precond_residual, 1_ip)

    end subroutine precond_rotation_qn

    subroutine keep_pair(self, i, n_kept)
        !
        ! this subroutine keeps the i-th pair of the history as the n_kept-th, so that
        ! the pairs the updating scheme skips are dropped from the history
        !
        class(updating_type), intent(inout) :: self
        integer(ip), intent(in) :: i
        integer(ip), intent(inout) :: n_kept

        n_kept = n_kept + 1
        self%kappa_list(:, n_kept) = self%kappa_list(:, i)
        self%grad_diff_list(:, n_kept) = self%grad_diff_list(:, i)

    end subroutine keep_pair

    subroutine build_model_sr1(self, error)
        !
        ! this subroutine builds the model of the SR1 updating scheme by applying the
        ! symmetric rank-1 update of every pair of the history, from the oldest to the
        ! newest, to the initial Hessian, B_j = B_{j-1} + r r^T / (s^T r) with
        ! r = y - B_{j-1} s, the recursive form of the limited-memory SR1 matrix of
        ! Byrd, R.H., Nocedal, J. & Schnabel, R.B. Mathematical Programming 63,
        ! 129–156 (1994), https://doi.org/10.1007/BF01582063, and drops every pair
        ! whose denominator is too small
        !
        class(sr1_updating_type), intent(inout) :: self
        integer(ip), intent(out) :: error

        integer(ip) :: i, n_kept
        real(rp) :: residual(self%n_param), denominator
        logical :: kept(self%n_points)
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        self%n_terms = 0
        kept = .false.
        do i = 1, self%n_points
            ! residual of the secant condition of the model so far
            call self%model_x(self%kappa_list(:, i), residual, self%n_terms, error)
            if (error /= 0) return
            residual = self%grad_diff_list(:, i) - residual

            ! skip the pair if the denominator is too small
            denominator = ddot(self%n_param, self%kappa_list(:, i), 1_ip, residual, &
                               1_ip)
            if (abs(denominator) <= &
                sr1_skip_thresh * dnrm2(self%n_param, self%kappa_list(:, i), 1_ip) * &
                dnrm2(self%n_param, residual, 1_ip)) cycle

            ! add the term of the pair
            kept(i) = .true.
            self%n_terms = self%n_terms + 1
            self%vectors(:, self%n_terms) = residual
            self%coeffs(self%n_terms) = 1.0_rp / denominator
        end do

        ! drop the skipped pairs from the history only now, so that an error leaves it
        ! whole
        n_kept = 0
        do i = 1, size(kept)
            if (kept(i)) call keep_pair(self, i, n_kept)
        end do
        self%n_points = n_kept

    end subroutine build_model_sr1

    subroutine build_model_bfgs(self, error)
        !
        ! this subroutine builds the model of the BFGS updating scheme by applying the
        ! BFGS update of every pair of the history, from the oldest to the newest, to
        ! the initial Hessian, B_j = B_{j-1} - B_{j-1} s s^T B_{j-1} / (s^T B_{j-1} s)
        ! + y y^T / (s^T y), the recursive form of the limited-memory BFGS matrix of
        ! Byrd, R.H., Nocedal, J. & Schnabel, R.B. Mathematical Programming 63,
        ! 129–156 (1994), https://doi.org/10.1007/BF01582063, and drops every pair
        ! whose curvature along its step, in the history or in the model so far, is
        ! not sufficiently positive, which keeps the model positive definite if the
        ! initial Hessian is
        !
        class(bfgs_updating_type), intent(inout) :: self
        integer(ip), intent(out) :: error

        integer(ip) :: i, n_kept
        real(rp) :: model_kappa(self%n_param), kappa_grad_diff, kappa_model_kappa
        logical :: kept(self%n_points)
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        self%n_terms = 0
        kept = .false.
        do i = 1, self%n_points
            ! skip the pair if its curvature is not sufficiently positive
            kappa_grad_diff = ddot(self%n_param, self%kappa_list(:, i), 1_ip, &
                                   self%grad_diff_list(:, i), 1_ip)
            if (kappa_grad_diff <= &
                bfgs_skip_thresh * dnrm2(self%n_param, self%kappa_list(:, i), 1_ip) * &
                dnrm2(self%n_param, self%grad_diff_list(:, i), 1_ip)) cycle

            ! skip the pair if the curvature of the model so far is not sufficiently
            ! positive
            call self%model_x(self%kappa_list(:, i), model_kappa, self%n_terms, error)
            if (error /= 0) return
            kappa_model_kappa = ddot(self%n_param, self%kappa_list(:, i), 1_ip, &
                                     model_kappa, 1_ip)
            if (kappa_model_kappa <= &
                bfgs_skip_thresh * dnrm2(self%n_param, self%kappa_list(:, i), 1_ip) * &
                dnrm2(self%n_param, model_kappa, 1_ip)) cycle

            ! add the two terms of the pair
            kept(i) = .true.
            self%vectors(:, self%n_terms + 1) = model_kappa
            self%coeffs(self%n_terms + 1) = -1.0_rp / kappa_model_kappa
            self%vectors(:, self%n_terms + 2) = self%grad_diff_list(:, i)
            self%coeffs(self%n_terms + 2) = 1.0_rp / kappa_grad_diff
            self%n_terms = self%n_terms + 2
        end do

        ! drop the skipped pairs from the history only now, so that an error leaves it
        ! whole
        n_kept = 0
        do i = 1, size(kept)
            if (kept(i)) call keep_pair(self, i, n_kept)
        end do
        self%n_points = n_kept

    end subroutine build_model_bfgs

end module otr_rotation_qn
