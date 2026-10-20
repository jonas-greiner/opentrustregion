! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common

    use opentrustregion, only: &
        rp, ip, settings_type, solver_settings_type, stability_settings_type, &
        solver_state_type, update_orbs_type, obj_func_type, hess_x_type, precond_type, &
        precond_pd_type, project_type, get_extra_trial_vectors_type, conv_check_type, &
        logger_type, add_error_origin

    implicit none

    ! settings shared by the extensions parameterizing the orbitals in an orbital
    ! basis, which the orbital objects hold
    type, extends(settings_type) :: orbital_settings_type
    end type

    ! callback functions and host context a factory found at one level of the solver
    ! settings before it wired its own, to which the callback functions it wired
    ! forward the roles its extension does not provide, together with the Hessian
    ! linear transformation the host's orbital update returned last
    type :: host_record_type
        class(*), pointer :: context => null()
        procedure(update_orbs_type), pointer, nopass :: update_orbs => null()
        procedure(obj_func_type), pointer, nopass :: obj_func => null()
        procedure(hess_x_type), pointer, nopass :: hess_x => null()
        procedure(precond_type), pointer, nopass :: precond => null()
        procedure(precond_pd_type), pointer, nopass :: precond_pd => null()
        procedure(project_type), pointer, nopass :: project => null()
        procedure(get_extra_trial_vectors_type), pointer, nopass :: &
            get_extra_trial_vectors => null()
        procedure(conv_check_type), pointer, nopass :: conv_check => null()
        procedure(hess_x_type), pointer, nopass :: stability_hess_x => null()
        procedure(hess_x_type), pointer, nopass :: approx_hess_x => null()
        procedure(logger_type), pointer, nopass :: logger => null()
    end type

    ! one level of an extension object, that of the solver or of its stability check,
    ! which the factory makes the context of the solver settings or their stability
    ! check settings, so that the callback functions it wires find the level that calls
    ! them, holding the object and the host record of the level
    type :: extension_level_type
        class(extension_type), pointer :: extension => null()
        type(host_record_type) :: host = host_record_type()
    end type

    ! roles an extension object performs itself, which have to be exactly those whose
    ! type-bound procedures its type overrides: the C wrappers of the other roles
    ! forward straight to the host's C functions, and the factory wires the role
    ! callback functions of the roles of the solver settings it performs also where the
    ! host provided none
    type :: provided_roles_type
        logical :: update_orbs = .false., obj_func = .false., hess_x = .false., &
                   precond = .false., precond_pd = .false., project = .false., &
                   get_extra_trial_vectors = .false.
    end type

    ! object of an extension, holding the number of parameters, the levels of the
    ! solver and of its stability check, the roles it performs itself, the origin with
    ! which it tags the errors it detects itself (none for an object no factory set up)
    ! and the state of the solver, which the host keeps current, whose type-bound
    ! procedures perform its roles, by default forwarding to the host record of the
    ! level they receive as context
    type, abstract :: extension_type
        integer(ip) :: n_param, error_origin = 0
        type(extension_level_type) :: solver_level, stability_level
        type(provided_roles_type) :: provides
        type(solver_state_type), pointer :: state => null()
    contains
        procedure :: update_orbs => update_orbs_extension
        procedure :: obj_func => obj_func_extension
        procedure :: hess_x => hess_x_extension
        procedure :: precond => precond_extension
        procedure :: precond_pd => precond_pd_extension
        procedure :: project => project_extension
        procedure :: get_extra_trial_vectors => get_extra_trial_vectors_extension
        procedure :: wire_solver_settings => wire_solver_settings_extension
    end type

    ! extension object owned by solver settings
    type :: owned_extension_type
        class(extension_type), allocatable :: object
    end type

    ! solver settings which own the extension objects wired into them, in the order of
    ! their wiring, so that the objects are freed with the settings, while the solver's
    ! copy of the settings takes only the parent part and therefore never copies them;
    ! only their wiring changes the list of objects, which hosts never modify
    type, extends(solver_settings_type) :: extension_solver_settings_type
        type(owned_extension_type), allocatable :: extensions(:)
    contains
        procedure :: wire => wire_extension_solver_settings
    end type

    ! warning for orbital bases set up without a density matrix evaluating function
    character(len=*), parameter :: evaluate_dm_missing_warning_msg = &
        "No density matrix evaluating function was given, so the orbital updating "// &
        "function, objective function and Hessian linear transformation of the "// &
        "orbital basis are not available and return an error when called."

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_callback_ptr => &
        update_orbs_callback
    procedure(obj_func_type), pointer :: obj_func_callback_ptr => obj_func_callback
    procedure(hess_x_type), pointer :: hess_x_callback_ptr => hess_x_callback
    procedure(precond_type), pointer :: precond_callback_ptr => precond_callback
    procedure(precond_pd_type), pointer :: precond_pd_callback_ptr => &
        precond_pd_callback
    procedure(project_type), pointer :: project_callback_ptr => project_callback
    procedure(get_extra_trial_vectors_type), pointer :: &
        get_extra_trial_vectors_callback_ptr => get_extra_trial_vectors_callback
    procedure(conv_check_type), pointer :: conv_check_callback_ptr => &
        conv_check_callback
    procedure(hess_x_type), pointer :: stability_hess_x_callback_ptr => &
        stability_hess_x_callback
    procedure(hess_x_type), pointer :: approx_hess_x_callback_ptr => &
        approx_hess_x_callback
    procedure(logger_type), pointer :: logger_callback_ptr => logger_callback

    ! density matrix evaluating functions and the response functions they return, with
    ! which the orbital bases evaluate the exact Hessian
    abstract interface
        subroutine get_response_cs_type(dm, response, error, context)
            import :: rp, ip

            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out), target, contiguous :: response(:, :)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine get_response_cs_type
    end interface

    abstract interface
        subroutine get_response_os_type(dm, response, error, context)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out), target :: response(:, :, :)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine get_response_os_type
    end interface

    abstract interface
        subroutine evaluate_dm_cs_type(dm, energy, fock, get_response_funptr, error, &
                                       context)
            import :: rp, ip

            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), optional, target, contiguous :: fock(:, :)
            procedure(get_response_cs_type), intent(out), optional, pointer :: &
                get_response_funptr
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine evaluate_dm_cs_type
    end interface

    abstract interface
        subroutine evaluate_dm_os_type(dm, energy, fock, get_response_funptr, error, &
                                       context)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), optional, target :: fock(:, :, :)
            procedure(get_response_os_type), intent(out), optional, pointer :: &
                get_response_funptr
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine evaluate_dm_os_type
    end interface

    ! orbital basis in which the orbital rotations are parameterized, an extension
    ! object holding the quantities every basis has together with the density matrix
    ! evaluating functions and the response functions they return, with which the basis
    ! evaluates the exact Hessian, whose type-bound procedures perform the operations
    ! every basis provides: moving the current orbitals, the gradient and Hessian
    ! diagonal at the current density, the response at the current density, and the
    ! eigendecomposition of the static part of the Hessian with the preconditioners and
    ! extra trial vectors built from it
    type, abstract, extends(extension_type) :: orbital_basis_type
        type(orbital_settings_type) :: settings
        integer(ip) :: n_ao, n_particle
        real(rp) :: energy
        real(rp), pointer, contiguous :: dm_ao(:, :, :) => null()
        real(rp), allocatable :: grad(:), h_diag(:)
        logical :: evaluation_stale = .true., response_stale = .true., &
                   hess_eigen_stale = .true.
        procedure(evaluate_dm_cs_type), pointer, nopass :: evaluate_dm_cs => null()
        procedure(get_response_cs_type), pointer, nopass :: get_response_cs => null()
        procedure(evaluate_dm_os_type), pointer, nopass :: evaluate_dm_os => null()
        procedure(get_response_os_type), pointer, nopass :: get_response_os => null()
    contains
        procedure(rotate_orbitals_basis_type), deferred :: rotate_orbitals
        procedure(calculate_grad_h_diag_basis_type), deferred :: calculate_grad_h_diag
        procedure(refresh_hess_eigen_basis_type), deferred :: refresh_hess_eigen
        procedure(rotate_hess_eigenbasis_basis_type), deferred :: &
            rotate_to_hess_eigenbasis
        procedure(rotate_hess_eigenbasis_basis_type), deferred :: &
            rotate_from_hess_eigenbasis
        procedure(get_hess_eigval_pairs_basis_type), deferred :: get_hess_eigval_pairs
        procedure :: check_evaluate_dm => check_evaluate_dm_orbital_basis
        procedure :: refresh_response => refresh_response_orbital_basis
        procedure :: precond => precond_orbital_basis
        procedure :: precond_pd => precond_pd_orbital_basis
        procedure :: fill_extra_trial_vectors => fill_extra_trial_vectors_orbital_basis
    end type

    abstract interface
        subroutine rotate_orbitals_basis_type(self, kappa, error)
            import :: orbital_basis_type, rp, ip

            class(orbital_basis_type), intent(inout) :: self
            real(rp), intent(in) :: kappa(:)
            integer(ip), intent(out) :: error
        end subroutine rotate_orbitals_basis_type

        subroutine calculate_grad_h_diag_basis_type(self, fock)
            import :: orbital_basis_type, rp

            class(orbital_basis_type), intent(inout) :: self
            real(rp), intent(in) :: fock(:, :, :)
        end subroutine calculate_grad_h_diag_basis_type

        subroutine refresh_hess_eigen_basis_type(self, error)
            import :: orbital_basis_type, ip

            class(orbital_basis_type), intent(inout) :: self
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
    end interface

contains

    subroutine wire_solver_settings_extension(self, solver_settings, settings, error, &
                                              update_orbs, obj_func)
        !
        ! this subroutine saves the callback functions and host contexts of the solver
        ! settings and their stability check settings and, for an object wrapping the
        ! host's, the host's orbital update and objective function in the host records
        ! of the levels of the object, and wires the role callback function of every
        ! role the host provided or the object performs into the settings, with the
        ! level as their context, after refusing an orbital update or objective
        ! function of the host for settings a factory wired, or one a factory returned
        ! for settings no factory wired
        !
        use opentrustregion, only: inherit_solver_callbacks, verbosity_error

        class(extension_type), intent(inout), target :: self
        type(solver_settings_type), intent(inout), target :: solver_settings
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        procedure(update_orbs_type), intent(in), pointer, optional :: update_orbs
        procedure(obj_func_type), intent(in), pointer, optional :: obj_func

        type(stability_settings_type), pointer :: nested
        type(solver_settings_type) :: inherited
        logical :: wired, wire_nested

        ! initialize error flag
        error = 0

        ! the orbital update and objective function the object wraps have to be those
        ! the factory that wired the settings last returned, which reach its object
        ! through the context of the settings, or the host's own for settings no
        ! factory wired, since the host's own would otherwise receive the level of that
        ! object as context and the returned ones the host's context
        wired = associated(extension_level(solver_settings%context))
        if (present(update_orbs)) then
            if (associated(update_orbs, update_orbs_callback) .neqv. wired) error = 1
        end if
        if (present(obj_func)) then
            if (associated(obj_func, obj_func_callback) .neqv. wired) error = 1
        end if
        if (error /= 0) then
            call settings%log("The orbital updating function and the objective "// &
                              "function have to be those the factory that wired "// &
                              "the solver settings last returned, or the host's "// &
                              "own for settings no factory wired.", verbosity_error, &
                              .true.)
            return
        end if

        ! save the host's callback functions and context of the solver in the level of
        ! the solver
        self%solver_level = extension_level_type()
        self%solver_level%extension => self
        associate (host => self%solver_level%host)
            host%context => solver_settings%context
            if (present(update_orbs)) host%update_orbs => update_orbs
            if (present(obj_func)) host%obj_func => obj_func
            host%precond => solver_settings%precond
            host%precond_pd => solver_settings%precond_pd
            host%project => solver_settings%project
            host%get_extra_trial_vectors => solver_settings%get_extra_trial_vectors
            host%conv_check => solver_settings%conv_check
            host%stability_hess_x => solver_settings%stability_hess_x
            host%logger => solver_settings%logger
        end associate

        ! save those of the stability check in its level, inheriting the solver's
        ! preconditioner, projection, logging function and context where the nested
        ! settings provide none, as the solver lets them, and the Hessian linear
        ! transformation of the stability check the solver settings hold, which the
        ! stability check calls with its own context
        inherited = solver_settings
        call inherit_solver_callbacks(inherited)
        self%stability_level = extension_level_type()
        self%stability_level%extension => self
        associate (host => self%stability_level%host, &
                   stability => inherited%stability_settings)
            host%context => stability%context
            host%precond => stability%precond
            host%project => stability%project
            host%logger => stability%logger
            host%get_extra_trial_vectors => stability%get_extra_trial_vectors
            host%stability_hess_x => solver_settings%stability_hess_x
            host%approx_hess_x => stability%approx_hess_x
            host%conv_check => stability%conv_check
        end associate

        ! wire the role callback functions of the roles the host provided or the object
        ! performs into the solver settings, leaving the others unset, since some change
        ! the algorithm by being set
        solver_settings%context => self%solver_level
        associate (host => self%solver_level%host, provides => self%provides)
            solver_settings%precond => null()
            if (associated(host%precond) .or. provides%precond) &
                solver_settings%precond => precond_callback
            solver_settings%precond_pd => null()
            if (associated(host%precond_pd) .or. provides%precond_pd) &
                solver_settings%precond_pd => precond_pd_callback
            solver_settings%project => null()
            if (associated(host%project) .or. provides%project) &
                solver_settings%project => project_callback
            solver_settings%get_extra_trial_vectors => null()
            if (associated(host%get_extra_trial_vectors) .or. &
                provides%get_extra_trial_vectors) &
                solver_settings%get_extra_trial_vectors => &
                get_extra_trial_vectors_callback
            solver_settings%conv_check => null()
            if (associated(host%conv_check)) &
                solver_settings%conv_check => conv_check_callback
            solver_settings%stability_hess_x => null()
            if (associated(host%stability_hess_x)) &
                solver_settings%stability_hess_x => stability_hess_x_callback
            solver_settings%logger => null()
            if (associated(host%logger)) solver_settings%logger => logger_callback
        end associate

        ! stability check settings which provide nothing of their own inherit the
        ! callback functions and context of the solver, which forward with the host
        ! context of the solver as before, and are left untouched unless the object
        ! performs a role of the stability check, since it does not inherit the extra
        ! trial vectors, while the others are wired like the solver settings with the
        ! level of the stability check as context
        nested => solver_settings%stability_settings
        wire_nested = associated(nested%context) .or. associated(nested%precond) .or. &
                      associated(nested%project) .or. associated(nested%logger) .or. &
                      associated(nested%get_extra_trial_vectors) .or. &
                      associated(nested%approx_hess_x) .or. &
                      associated(nested%conv_check) .or. self%provides%precond .or. &
                      self%provides%project .or. self%provides%get_extra_trial_vectors
        if (.not. wire_nested) return
        nested%context => self%stability_level
        associate (host => self%stability_level%host, provides => self%provides)
            nested%precond => null()
            if (associated(host%precond) .or. provides%precond) &
                nested%precond => precond_callback
            nested%project => null()
            if (associated(host%project) .or. provides%project) &
                nested%project => project_callback
            nested%get_extra_trial_vectors => null()
            if (associated(host%get_extra_trial_vectors) .or. &
                provides%get_extra_trial_vectors) &
                nested%get_extra_trial_vectors => get_extra_trial_vectors_callback
            nested%approx_hess_x => null()
            if (associated(host%approx_hess_x)) &
                nested%approx_hess_x => approx_hess_x_callback
            nested%conv_check => null()
            if (associated(host%conv_check)) nested%conv_check => conv_check_callback
            nested%logger => null()
            if (associated(host%logger)) nested%logger => logger_callback
        end associate

    end subroutine wire_solver_settings_extension

    subroutine wire_extension_solver_settings(self, object, solver_state, settings, &
                                              error, update_orbs, obj_func)
        !
        ! this subroutine points a new extension object at the record of the solver
        ! state, wires the solver settings with it and lets them take it over, after
        ! the objects they owned before, so that the objects whose role callback
        ! functions the latest one forwards to live as long as it does and no wired
        ! settings point into an object they do not own; the objects are moved, not
        ! copied, so that pointers into them stay valid, and an object whose wiring
        ! failed stays with the caller
        !
        class(extension_solver_settings_type), intent(inout), target :: self
        class(extension_type), allocatable, intent(inout), target :: object
        type(solver_state_type), intent(inout), target :: solver_state
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        procedure(update_orbs_type), intent(in), pointer, optional :: update_orbs
        procedure(obj_func_type), intent(in), pointer, optional :: obj_func

        type(owned_extension_type), allocatable, target :: extensions(:)
        integer(ip) :: n_owned, i

        ! point the object at the record of the state and wire the settings with it
        object%state => solver_state
        call object%wire_solver_settings(self%solver_settings_type, settings, error, &
                                         update_orbs, obj_func)
        if (error /= 0) return

        ! move the objects the settings own and the new one into a longer list
        n_owned = 0
        if (allocated(self%extensions)) n_owned = size(self%extensions)
        allocate(extensions(n_owned + 1))
        do i = 1, n_owned
            call move_alloc(self%extensions(i)%object, extensions(i)%object)
        end do
        call move_alloc(object, extensions(n_owned + 1)%object)
        call move_alloc(extensions, self%extensions)

    end subroutine wire_extension_solver_settings

    recursive subroutine update_orbs_extension(self, kappa, func, grad, h_diag, &
                                               hess_x_funptr, error, context)
        !
        ! this subroutine performs the host's orbital update saved in the level of the
        ! object it receives as context, with the host context of the level, for an
        ! extension that does not provide its own, keeps the Hessian linear
        ! transformation it returns in the host records of both levels, since the
        ! stability check calls it too, and returns the role callback function in its
        ! place, also if only the object provides one
        !
        class(extension_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level
        procedure(hess_x_type), pointer :: host_hess_x_funptr

        ! call the host's orbital update
        level => extension_level(context)
        host_hess_x_funptr => null()
        call level%host%update_orbs(kappa, func, grad, h_diag, host_hess_x_funptr, &
                                    error, level%host%context)

        ! keep the host's Hessian linear transformation and return the role callback
        ! function, through which the solver reaches it with the level as context
        self%solver_level%host%hess_x => host_hess_x_funptr
        self%stability_level%host%hess_x => host_hess_x_funptr
        hess_x_funptr => null()
        if (associated(host_hess_x_funptr) .or. self%provides%hess_x) &
            hess_x_funptr => hess_x_callback

    end subroutine update_orbs_extension

    recursive function obj_func_extension(self, kappa, error, context) result(func)
        !
        ! this function evaluates the host's objective function saved in the level of
        ! the object it receives as context, with the host context of the level, for an
        ! extension that does not provide its own
        !
        class(extension_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        func = level%host%obj_func(kappa, error, level%host%context)

    end function obj_func_extension

    recursive subroutine hess_x_extension(self, x, hess_x, error, context)
        !
        ! this subroutine performs the Hessian linear transformation the host's orbital
        ! update returned last, saved in the level of the object it receives as
        ! context, with the host context of the level, for an extension that does not
        ! provide its own
        !
        class(extension_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        call level%host%hess_x(x, hess_x, error, level%host%context)

    end subroutine hess_x_extension

    recursive subroutine precond_extension(self, residual, mu, precond_residual, &
                                           error, context)
        !
        ! this subroutine applies the host's level-shifted preconditioner saved in the
        ! level of the object it receives as context, with the host context of the
        ! level, for an extension that does not provide its own
        !
        class(extension_type), intent(inout) :: self
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        call level%host%precond(residual, mu, precond_residual, error, &
                                level%host%context)

    end subroutine precond_extension

    recursive subroutine precond_pd_extension(self, residual, precond_residual, error, &
                                              context)
        !
        ! this subroutine applies the host's positive-definite preconditioner saved in
        ! the level of the object it receives as context, with the host context of the
        ! level, for an extension that does not provide its own
        !
        class(extension_type), intent(inout) :: self
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        call level%host%precond_pd(residual, precond_residual, error, &
                                   level%host%context)

    end subroutine precond_pd_extension

    recursive subroutine project_extension(self, vector, error, context)
        !
        ! this subroutine applies the host's projection saved in the level of the
        ! object it receives as context, with the host context of the level, for an
        ! extension that does not provide its own
        !
        class(extension_type), intent(inout) :: self
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        call level%host%project(vector, error, level%host%context)

    end subroutine project_extension

    recursive subroutine get_extra_trial_vectors_extension(self, trial_vectors, error, &
                                                           context)
        !
        ! this subroutine returns the extra trial vectors of the host's function saved
        ! in the level of the object it receives as context, with the host context of
        ! the level, for an extension that does not provide its own
        !
        class(extension_type), intent(inout) :: self
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        call level%host%get_extra_trial_vectors(trial_vectors, error, &
                                                level%host%context)

    end subroutine get_extra_trial_vectors_extension

    function extension_level(context) result(level)
        !
        ! this function returns the level of an extension object a callback function
        ! receives as context, or a disassociated pointer for a missing context or one
        ! of another type
        !
        class(*), intent(in), pointer :: context
        type(extension_level_type), pointer :: level

        level => null()
        if (.not. associated(context)) return
        select type (context)
        type is (extension_level_type)
            level => context
        end select

    end function extension_level

    recursive subroutine update_orbs_callback(kappa, func, grad, h_diag, &
                                              hess_x_funptr, error, context)
        !
        ! this subroutine performs the orbital update of the extension object whose
        ! level it receives as context
        !
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            func = 0.0_rp
            hess_x_funptr => null()
            error = 1
            return
        end if
        call level%extension%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, &
                                         error, context)

    end subroutine update_orbs_callback

    recursive function obj_func_callback(kappa, error, context) result(func)
        !
        ! this function evaluates the objective function of the extension object whose
        ! level it receives as context
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            func = 0.0_rp
            error = 1
            return
        end if
        func = level%extension%obj_func(kappa, error, context)

    end function obj_func_callback

    recursive subroutine hess_x_callback(x, hess_x, error, context)
        !
        ! this subroutine performs the Hessian linear transformation of the extension
        ! object whose level it receives as context
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            error = 1
            return
        end if
        call level%extension%hess_x(x, hess_x, error, context)

    end subroutine hess_x_callback

    recursive subroutine precond_callback(residual, mu, precond_residual, error, &
                                          context)
        !
        ! this subroutine applies the level-shifted preconditioner of the extension
        ! object whose level it receives as context
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            error = 1
            return
        end if
        call level%extension%precond(residual, mu, precond_residual, error, context)

    end subroutine precond_callback

    recursive subroutine precond_pd_callback(residual, precond_residual, error, context)
        !
        ! this subroutine applies the positive-definite preconditioner of the extension
        ! object whose level it receives as context
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            error = 1
            return
        end if
        call level%extension%precond_pd(residual, precond_residual, error, context)

    end subroutine precond_pd_callback

    recursive subroutine project_callback(vector, error, context)
        !
        ! this subroutine applies the projection of the extension object whose level it
        ! receives as context
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            error = 1
            return
        end if
        call level%extension%project(vector, error, context)

    end subroutine project_callback

    recursive subroutine get_extra_trial_vectors_callback(trial_vectors, error, context)
        !
        ! this subroutine returns the extra trial vectors of the extension object whose
        ! level it receives as context
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            error = 1
            return
        end if
        call level%extension%get_extra_trial_vectors(trial_vectors, error, context)

    end subroutine get_extra_trial_vectors_callback

    recursive function conv_check_callback(error, context) result(converged)
        !
        ! this function checks the convergence of the solver or of its stability check
        ! with the host's convergence check saved in the level of an extension object
        ! it receives as context, with the host context of the level
        !
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            converged = .false.
            error = 1
            return
        end if
        converged = level%host%conv_check(error, level%host%context)

    end function conv_check_callback

    recursive subroutine stability_hess_x_callback(x, hess_x, error, context)
        !
        ! this subroutine performs the host's Hessian linear transformation of the
        ! stability check saved in the level of an extension object it receives as
        ! context, with the host context of the level
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            error = 1
            return
        end if
        call level%host%stability_hess_x(x, hess_x, error, level%host%context)

    end subroutine stability_hess_x_callback

    recursive subroutine approx_hess_x_callback(x, hess_x, error, context)
        !
        ! this subroutine performs the host's approximate Hessian linear transformation
        ! saved in the level of an extension object it receives as context, with the
        ! host context of the level
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (.not. associated(level)) then
            error = 1
            return
        end if
        call level%host%approx_hess_x(x, hess_x, error, level%host%context)

    end subroutine approx_hess_x_callback

    recursive subroutine logger_callback(message, context)
        !
        ! this subroutine logs a message with the host's logging function saved in the
        ! level of an extension object it receives as context, with the host context of
        ! the level, and discards it without a level
        !
        character(len=*), intent(in) :: message
        class(*), intent(in), pointer :: context

        type(extension_level_type), pointer :: level

        level => extension_level(context)
        if (associated(level)) call level%host%logger(message, level%host%context)

    end subroutine logger_callback

    subroutine check_evaluate_dm_orbital_basis(self, error)
        !
        ! this subroutine produces an error for an orbital basis set up without a
        ! density matrix evaluating function, whose orbital update, objective function
        ! and Hessian linear transformation need one
        !
        use opentrustregion, only: verbosity_error

        class(orbital_basis_type), intent(in) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! one of the density matrix evaluating functions is set if one was given
        if (associated(self%evaluate_dm_cs) .or. associated(self%evaluate_dm_os)) return
        call self%settings%log("The orbital basis was set up without a density "// &
                               "matrix evaluating function, which this function "// &
                               "needs.", verbosity_error, .true.)
        error = 1
        call add_error_origin(error, self%error_origin, self%settings)

    end subroutine check_evaluate_dm_orbital_basis

    subroutine refresh_response_orbital_basis(self, error)
        !
        ! this subroutine rebuilds the response callbacks at the currently stored
        ! density matrix of the orbital basis
        !
        class(orbital_basis_type), intent(inout) :: self
        integer(ip), intent(out) :: error

        real(rp) :: energy

        ! initialize error flag
        error = 0

        ! rebuild response, the energy is discarded since the caller already holds it
        ! for the current density
        if (associated(self%evaluate_dm_cs)) then
            call self%evaluate_dm_cs(self%dm_ao(:, :, 1), energy, &
                                     get_response_funptr=self%get_response_cs, &
                                     error=error, context=self%settings%context)
        else
            call self%evaluate_dm_os(self%dm_ao, energy, &
                                     get_response_funptr=self%get_response_os, &
                                     error=error, context=self%settings%context)
        end if
        if (error /= 0) return

        ! the response callbacks were just rebuilt at the current density
        self%response_stale = .false.

    end subroutine refresh_response_orbital_basis

    subroutine precond_orbital_basis(self, residual, mu, precond_residual, error, &
                                     context)
        !
        ! this subroutine applies a level-shifted preconditioner based on the exact
        ! eigendecomposition of the static part of the Hessian of the orbital basis
        !
        class(orbital_basis_type), intent(inout) :: self
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp), allocatable :: rotated_residual(:), eigval_pairs(:)

        ! initialize error flag
        error = 0

        ! refresh the eigendecomposition if the static Hessian part has changed
        call self%refresh_hess_eigen(error)
        if (error /= 0) return

        ! rotate residual into the eigenbasis of the static Hessian part
        rotated_residual = self%rotate_to_hess_eigenbasis(residual)

        ! get eigenvalue pairs of the static Hessian part
        eigval_pairs = self%get_hess_eigval_pairs()

        ! apply level-shifted diagonal scaling in the eigenbasis
        rotated_residual = rotated_residual / level_shifted_divisors(eigval_pairs, mu)

        ! rotate back to the original basis
        precond_residual = self%rotate_from_hess_eigenbasis(rotated_residual)

    end subroutine precond_orbital_basis

    subroutine precond_pd_orbital_basis(self, residual, precond_residual, error, &
                                        context)
        !
        ! this subroutine applies the positive-definite preconditioner based on the
        ! exact eigendecomposition of the static part of the Hessian of the orbital
        ! basis
        !
        class(orbital_basis_type), intent(inout) :: self
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp), allocatable :: rotated_residual(:), eigval_pairs(:)

        ! initialize error flag
        error = 0

        ! refresh the eigendecomposition if the static Hessian part has changed
        call self%refresh_hess_eigen(error)
        if (error /= 0) return

        ! rotate residual into the eigenbasis of the static Hessian part
        rotated_residual = self%rotate_to_hess_eigenbasis(residual)

        ! get eigenvalue pairs of the static Hessian part
        eigval_pairs = self%get_hess_eigval_pairs()

        ! apply positive-definite diagonal scaling in the eigenbasis
        rotated_residual = rotated_residual / positive_definite_divisors(eigval_pairs)

        ! rotate back to the original basis
        precond_residual = self%rotate_from_hess_eigenbasis(rotated_residual)

    end subroutine precond_pd_orbital_basis

    subroutine fill_extra_trial_vectors_orbital_basis(self, eigval_pairs, trial_vectors)
        !
        ! this subroutine fills the extra trial vectors with the rotations belonging to
        ! the most negative of the given eigenvalue pairs of the static part of the
        ! Hessian, in increasing order, rotated out of its eigenbasis
        !
        class(orbital_basis_type), intent(in) :: self
        real(rp), intent(in) :: eigval_pairs(:)
        real(rp), intent(out) :: trial_vectors(:, :)

        integer(ip) :: ivec, min_idx
        real(rp), allocatable :: remaining_pairs(:), unit_vector(:)

        ! a vanishing vector tells the solver that no direction is contributed for that
        ! slot, which is what is returned for every slot left unfilled below
        trial_vectors = 0.0_rp

        ! fill the requested slots with the rotations belonging to the most negative
        ! remaining eigenvalue pairs, stopping as soon as none is negative any more
        remaining_pairs = eigval_pairs
        allocate(unit_vector(size(eigval_pairs)))
        do ivec = 1, size(trial_vectors, 2, kind=ip)
            min_idx = minloc(remaining_pairs, dim=1)
            if (remaining_pairs(min_idx) >= 0.0_rp) exit
            unit_vector = 0.0_rp
            unit_vector(min_idx) = 1.0_rp
            trial_vectors(:, ivec) = self%rotate_from_hess_eigenbasis(unit_vector)
            remaining_pairs(min_idx) = huge(1.0_rp)
        end do

    end subroutine fill_extra_trial_vectors_orbital_basis

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
        ! of two matrices is larger than the Frobenius norm of the sum of both
        ! matrices; the scaling bounds the norm of the scaled matrix by 1/8, so that
        ! the terms fall below the threshold after about ten orders, while a matrix
        ! whose norm is not finite, for which no scale factor exists, is refused
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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

        ! compute Frobenius norm of A, which has to be finite
        A_norm = sqrt(sum(A**2))
        if (.not. ieee_is_finite(A_norm)) then
            call settings%log("Matrix exponential of a matrix whose norm is not "// &
                              "finite requested.", verbosity_error, .true.)
            error = 1
            return
        end if

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

            ! next expansion order
            i = i + 1
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
