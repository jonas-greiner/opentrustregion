! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface

    use opentrustregion_core, only: ip, rp, hess_x_type, solver_state_type
    use c_interface, only: &
        c_ip, c_rp, c_callbacks_type, solver_state_type_c, update_orbs_c_type, &
        obj_func_c_type, hess_x_c_type, precond_c_type, precond_pd_c_type, &
        project_c_type, get_extra_trial_vectors_c_type, conv_check_c_type, logger_c_type
    use otr_common, only: extension_level_type, extension_solver_settings_type, &
                          extension_level, update_orbs_callback, obj_func_callback, &
                          hess_x_callback, precond_callback, precond_pd_callback, &
                          project_callback, get_extra_trial_vectors_callback, &
                          evaluate_dm_cs_type, evaluate_dm_os_type, &
                          get_response_cs_type, get_response_os_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_ptr, c_funptr, c_funloc, &
                                           c_loc, c_f_pointer, c_f_procpointer, &
                                           c_associated, c_null_funptr

    implicit none

    ! C-interoperable interfaces for the density matrix evaluating functions and the
    ! response functions they return, which the orbital bases share
    abstract interface
        function evaluate_dm_c_type(dm_ao_c, energy_c, fock_c, get_response_c_funptr, &
                                    context_c) result(error_c) bind(C)
            import :: c_rp, c_ip, c_funptr, c_ptr

            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            real(c_rp), intent(out), optional :: fock_c(*)
            type(c_funptr), intent(out), optional :: get_response_c_funptr
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error_c
        end function evaluate_dm_c_type
    end interface

    abstract interface
        function get_response_c_type(dm_ao_c, response_c, context_c) result(error_c) &
            bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out), target :: response_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error_c
        end function get_response_c_type
    end interface

    ! bundle of the C callback functions of an orbital basis, which extends that of the
    ! core by the density matrix evaluating function and the response function it
    ! returned last
    type, extends(c_callbacks_type) :: orbital_c_callbacks_type
        procedure(evaluate_dm_c_type), pointer, nopass :: evaluate_dm => null()
        procedure(get_response_c_type), pointer, nopass :: get_response => null()
    end type

    ! one level of the holder of an extension object for C hosts, whose address is the
    ! context of the C wrappers wired into the C settings of the level, holding the
    ! level of the object, the bundle of the host's C callback functions and context of
    ! the level, for an object reading the state of the solver, the host's record of
    ! the state, from which the state the object reads is refreshed, and, for an object
    ! updating a copy of the host's orbitals since their real kind differs from C's,
    ! the copy and the host's array, to which the orbital update copies them back
    type :: holder_level_type
        type(extension_level_type), pointer :: level => null()
        type(c_callbacks_type) :: callbacks = c_callbacks_type()
        type(solver_state_type_c), pointer :: state_c => null()
        real(rp), pointer :: orbitals(:, :, :) => null()
        real(c_rp), pointer :: orbitals_c(:, :, :) => null()
    end type

    ! holder of an extension object for C hosts, with the levels of the solver and of
    ! its stability check, the Fortran solver settings the factory wires, which own the
    ! object of a holder a factory returns to the host, and, for such a holder, the
    ! bundle of the host's callback functions of the extension, which the context of
    ! the extension's settings points to, and the state the object reads; a holder
    ! holding an object must never be copied, since the copy would hold a copy of the
    ! object whose levels the levels of the copy do not hold
    type :: extension_holder_type
        type(holder_level_type) :: solver_level, stability_level
        type(extension_solver_settings_type) :: solver_settings
        class(c_callbacks_type), allocatable :: extension_callbacks
        type(solver_state_type) :: state
    end type

    ! interfaces for converting states from C to Fortran, with which the state an
    ! extension object reads is refreshed from the host's record
    interface assignment(=)
        module procedure assign_solver_state_f_c
        module procedure assign_stability_state_f_c
    end interface

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_c_type), pointer :: update_orbs_c_wrapper_ptr => &
        update_orbs_c_wrapper
    procedure(obj_func_c_type), pointer :: obj_func_c_wrapper_ptr => obj_func_c_wrapper
    procedure(hess_x_c_type), pointer :: hess_x_c_wrapper_ptr => hess_x_c_wrapper
    procedure(precond_c_type), pointer :: precond_c_wrapper_ptr => precond_c_wrapper
    procedure(precond_pd_c_type), pointer :: precond_pd_c_wrapper_ptr => &
        precond_pd_c_wrapper
    procedure(project_c_type), pointer :: project_c_wrapper_ptr => project_c_wrapper
    procedure(get_extra_trial_vectors_c_type), pointer :: &
        get_extra_trial_vectors_c_wrapper_ptr => get_extra_trial_vectors_c_wrapper
    procedure(conv_check_c_type), pointer :: conv_check_c_wrapper_ptr => &
        conv_check_c_wrapper
    procedure(hess_x_c_type), pointer :: stability_hess_x_c_wrapper_ptr => &
        stability_hess_x_c_wrapper
    procedure(hess_x_c_type), pointer :: approx_hess_x_c_wrapper_ptr => &
        approx_hess_x_c_wrapper
    procedure(logger_c_type), pointer :: logger_c_wrapper_ptr => logger_c_wrapper
    procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_f_wrapper_ptr => &
        evaluate_dm_cs_f_wrapper
    procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_f_wrapper_ptr => &
        evaluate_dm_os_f_wrapper
    procedure(get_response_cs_type), pointer :: get_response_cs_f_wrapper_ptr => &
        get_response_cs_f_wrapper
    procedure(get_response_os_type), pointer :: get_response_os_f_wrapper_ptr => &
        get_response_os_f_wrapper

contains

    function holder_level(context_c) result(held)
        !
        ! this function returns the level of a holder at the address a C wrapper
        ! receives as context, or a disassociated pointer for a null address
        !
        type(c_ptr), intent(in) :: context_c
        type(holder_level_type), pointer :: held

        held => null()
        if (c_associated(context_c)) call c_f_pointer(context_c, held)

    end function holder_level

    subroutine bundle_settings_c(callbacks, initialized, logger_c_funptr, context_c)
        !
        ! this subroutine hands the logging function and the host context of the C
        ! settings of an extension to the bundle in which the wrappers of the host's
        ! callback functions of the extension find them, only for initialized settings,
        ! since settings that were not initialized carry neither, as in the core
        !
        class(c_callbacks_type), intent(inout) :: callbacks
        logical(c_bool), intent(in) :: initialized
        type(c_funptr), intent(in) :: logger_c_funptr
        type(c_ptr), intent(in) :: context_c

        if (.not. initialized) return
        if (c_associated(logger_c_funptr)) &
            call c_f_procpointer(cptr=logger_c_funptr, fptr=callbacks%logger)
        callbacks%host_context = context_c

    end subroutine bundle_settings_c

    subroutine setup_extension_holder(holder, solver_settings_c, solver_state_c, &
                                      settings, error, update_orbs_c_funptr, &
                                      obj_func_c_funptr)
        !
        ! this subroutine sets up the holder of an extension object for a C factory by
        ! bundling the callback functions and host contexts of the host's C solver
        ! settings, which it converts to the Fortran solver settings of the holder the
        ! Fortran factory wires, and, for an object wrapping the host's, the host's
        ! orbital update and objective function, and by holding the host's record of
        ! the state, which starts from the default state as the solver's does and from
        ! which the state the object reads is refreshed
        !
        use opentrustregion_core, only: settings_type, verbosity_error
        use c_interface, only: solver_settings_type_c, convert_solver_settings_c, &
                               assignment(=)

        type(extension_holder_type), intent(inout), target :: holder
        type(solver_settings_type_c), intent(in) :: solver_settings_c
        type(c_ptr), intent(in) :: solver_state_c
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        type(c_funptr), intent(in), optional :: update_orbs_c_funptr, obj_func_c_funptr

        logical :: nested_own

        ! initialize error flag
        error = 0

        ! the object reads the state of the solver from the host's record, which is
        ! therefore required
        if (.not. c_associated(solver_state_c)) then
            call settings%log("A record of the solver state is required.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! bundle the host's callback functions and convert the solver settings, whose
        ! optional callback functions and host contexts join the bundles of the holder
        holder%solver_level = holder_level_type()
        holder%stability_level = holder_level_type()
        if (present(update_orbs_c_funptr)) call c_f_procpointer( &
            cptr=update_orbs_c_funptr, fptr=holder%solver_level%callbacks%update_orbs)
        if (present(obj_func_c_funptr)) call c_f_procpointer( &
            cptr=obj_func_c_funptr, fptr=holder%solver_level%callbacks%obj_func)
        call convert_solver_settings_c( &
            solver_settings_c, holder%solver_settings%solver_settings_type, &
            holder%solver_level%callbacks, holder%stability_level%callbacks)

        ! the conversion hands the nested bundle to the nested settings even if these
        ! provide nothing of their own, so their context is removed again for the
        ! factory to leave them inheriting the solver's, as for Fortran hosts
        associate (nested_c => solver_settings_c%stability_settings)
            nested_own = solver_settings_c%initialized .and. &
                         nested_c%initialized .and. &
                         (c_associated(nested_c%context) .or. &
                          c_associated(nested_c%precond) .or. &
                          c_associated(nested_c%project) .or. &
                          c_associated(nested_c%logger) .or. &
                          c_associated(nested_c%get_extra_trial_vectors) .or. &
                          c_associated(nested_c%approx_hess_x) .or. &
                          c_associated(nested_c%conv_check))
        end associate
        if (.not. nested_own) &
            holder%solver_settings%stability_settings%context => null()

        ! hold the host's record of the state at both levels and the state the object
        ! reads, both starting from the default state, since a C host's record need not
        ! be defined before the solver writes it
        call c_f_pointer(solver_state_c, holder%solver_level%state_c)
        holder%state = solver_state_type()
        holder%solver_level%state_c = holder%state
        holder%stability_level%state_c => holder%solver_level%state_c

    end subroutine setup_extension_holder

    subroutine refresh_state(held)
        !
        ! this subroutine refreshes the state the object of a level of a holder reads
        ! from the host's record the level holds, which the core copies the current
        ! state into before every callback function into C, before a C wrapper calls a
        ! role of the object
        !
        type(holder_level_type), intent(in) :: held

        if (.not. associated(held%state_c)) return
        if (associated(held%level%extension%state)) &
            held%level%extension%state = held%state_c

    end subroutine refresh_state

    subroutine assign_solver_state_f_c(state, state_c)
        !
        ! this subroutine converts a solver state from C to Fortran
        !
        type(solver_state_type), intent(out) :: state
        type(solver_state_type_c), intent(in) :: state_c

        ! convert logicals
        state%max_precision_reached = logical(state_c%max_precision_reached)

        ! convert reals
        state%func = real(state_c%func, kind=rp)
        state%grad_rms = real(state_c%grad_rms, kind=rp)
        state%trust_radius = real(state_c%trust_radius, kind=rp)
        state%relative_residual = real(state_c%relative_residual, kind=rp)
        state%level_shift = real(state_c%level_shift, kind=rp)
        state%step_size = real(state_c%step_size, kind=rp)

        ! convert integers
        state%macro_iteration = int(state_c%macro_iteration, kind=ip)
        state%micro_iteration = int(state_c%micro_iteration, kind=ip)
        state%jacobi_davidson_iteration = &
            int(state_c%jacobi_davidson_iteration, kind=ip)
        state%n_update_orbs = int(state_c%n_update_orbs, kind=ip)
        state%n_obj_func = int(state_c%n_obj_func, kind=ip)
        state%n_hess_x = int(state_c%n_hess_x, kind=ip)

        ! convert objects
        state%stability_state = state_c%stability_state

    end subroutine assign_solver_state_f_c

    subroutine assign_stability_state_f_c(state, state_c)
        !
        ! this subroutine converts a stability check state from C to Fortran
        !
        use opentrustregion_core, only: stability_state_type
        use c_interface, only: stability_state_type_c

        type(stability_state_type), intent(out) :: state
        type(stability_state_type_c), intent(in) :: state_c

        ! convert reals
        state%min_eigval = real(state_c%min_eigval, kind=rp)
        state%residual_rms = real(state_c%residual_rms, kind=rp)

        ! convert integers
        state%iteration = int(state_c%iteration, kind=ip)
        state%n_hess_x = int(state_c%n_hess_x, kind=ip)

    end subroutine assign_stability_state_f_c

    subroutine release_extension(extension_c) bind(C, name="release_extension")
        !
        ! this subroutine frees the holder of an extension object a C factory returned,
        ! together with the object
        !
        type(c_ptr), intent(in), value :: extension_c

        type(extension_holder_type), pointer :: holder

        if (.not. c_associated(extension_c)) return
        call c_f_pointer(extension_c, holder)
        deallocate(holder)

    end subroutine release_extension

    subroutine wire_solver_settings_c(holder, solver_settings_c)
        !
        ! this subroutine writes the Fortran solver settings of the holder, which the
        ! Fortran factory of an extension object wired, back into the host's C solver
        ! settings, with the C wrappers of the role callback functions in their place
        ! and the addresses of the levels of the holder as the contexts of the levels
        ! the factory wired, which also initializes C settings that were not
        ! initialized, as the conversion did, and holds the levels of the object, whose
        ! level of the solver the factory made the context of the settings
        !
        use c_interface, only: solver_settings_type_c, assignment(=)

        type(extension_holder_type), intent(inout), target :: holder
        type(solver_settings_type_c), intent(inout) :: solver_settings_c

        type(extension_level_type), pointer :: level

        ! hold the levels of the object, which the C wrappers find at the addresses of
        ! the levels of the holder
        level => extension_level(holder%solver_settings%context)
        holder%solver_level%level => level
        holder%stability_level%level => level%extension%stability_level

        associate (solver_settings => holder%solver_settings%solver_settings_type)
            ! convert the settings, which leaves the callback functions and contexts
            ! unset
            solver_settings_c = solver_settings

            ! the factory points the context of the solver settings at the level of the
            ! solver and wires the role callback function of every role that was set or
            ! that it provides
            solver_settings_c%context = c_loc(holder%solver_level)
            if (associated(solver_settings%precond)) &
                solver_settings_c%precond = c_funloc(precond_c_wrapper)
            if (associated(solver_settings%precond_pd)) &
                solver_settings_c%precond_pd = c_funloc(precond_pd_c_wrapper)
            if (associated(solver_settings%project)) &
                solver_settings_c%project = c_funloc(project_c_wrapper)
            if (associated(solver_settings%get_extra_trial_vectors)) &
                solver_settings_c%get_extra_trial_vectors = &
                c_funloc(get_extra_trial_vectors_c_wrapper)
            if (associated(solver_settings%conv_check)) &
                solver_settings_c%conv_check = c_funloc(conv_check_c_wrapper)
            if (associated(solver_settings%stability_hess_x)) &
                solver_settings_c%stability_hess_x = &
                c_funloc(stability_hess_x_c_wrapper)
            if (associated(solver_settings%logger)) &
                solver_settings_c%logger = c_funloc(logger_c_wrapper)

            ! the factory points the context of the nested settings at the level of the
            ! stability check only if it wired them, otherwise they provide nothing and
            ! inherit the solver's
            if (.not. associated(solver_settings%stability_settings%context)) return
            associate (nested => solver_settings%stability_settings, &
                       nested_c => solver_settings_c%stability_settings)
                nested_c%context = c_loc(holder%stability_level)
                if (associated(nested%precond)) &
                    nested_c%precond = c_funloc(precond_c_wrapper)
                if (associated(nested%project)) &
                    nested_c%project = c_funloc(project_c_wrapper)
                if (associated(nested%get_extra_trial_vectors)) &
                    nested_c%get_extra_trial_vectors = &
                    c_funloc(get_extra_trial_vectors_c_wrapper)
                if (associated(nested%approx_hess_x)) &
                    nested_c%approx_hess_x = c_funloc(approx_hess_x_c_wrapper)
                if (associated(nested%conv_check)) &
                    nested_c%conv_check = c_funloc(conv_check_c_wrapper)
                if (associated(nested%logger)) &
                    nested_c%logger = c_funloc(logger_c_wrapper)
            end associate
        end associate

    end subroutine wire_solver_settings_c

    subroutine finish_extension_holder(holder, error, solver_settings_c, &
                                       update_orbs_c_funptr, obj_func_c_funptr, &
                                       extension_c)
        !
        ! this subroutine finishes the factory call of a C host which created a holder
        ! for the new extension object, by writing the solver settings of the holder
        ! the Fortran factory wired back into the host's C solver settings and returning
        ! the C wrappers of the orbital update and the objective function and the
        ! holder as the extension, or by freeing the holder if the factory failed, which
        ! leaves the returned functions and the extension to the caller
        !
        use c_interface, only: solver_settings_type_c

        type(extension_holder_type), intent(inout), pointer :: holder
        integer(ip), intent(in) :: error
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_funptr), intent(inout) :: update_orbs_c_funptr, obj_func_c_funptr
        type(c_ptr), intent(inout) :: extension_c

        if (error == 0) then
            call wire_solver_settings_c(holder, solver_settings_c)
            update_orbs_c_funptr = c_funloc(update_orbs_c_wrapper)
            obj_func_c_funptr = c_funloc(obj_func_c_wrapper)
            extension_c = c_loc(holder)
        else
            deallocate(holder)
        end if

    end subroutine finish_extension_holder

    recursive function update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                             hess_x_c_funptr, context_c) &
        result(error_c) bind(C)
        !
        ! this function exposes the orbital update to C hosts, calling it with the
        ! level of an extension object held at the address it receives as context, or
        ! forwarding to the host's C function in the bundle of that level, with its
        ! host context, for an object that does not perform it, and returns the C
        ! wrapper of the Hessian linear transformation the orbital update returns
        !
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(update_orbs_c_type), pointer :: host_update_orbs
        type(c_funptr) :: host_hess_x_c_funptr
        class(*), pointer :: context
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp) :: func
        real(rp), pointer :: kappa(:), grad(:), h_diag(:)
        integer(ip) :: n_param, error

        ! find the level of the holder, without which there is nothing to update
        hess_x_c_funptr = c_null_funptr
        held => holder_level(context_c)
        if (.not. associated(held)) then
            func_c = 0.0_c_rp
            error_c = 1
            return
        end if
        context => held%level
        n_param = held%level%extension%n_param

        ! forward straight to the host's C function with its host context unless the
        ! object performs the role itself, keeping the Hessian linear transformation it
        ! returns in the bundles of both levels, since the stability check calls it too,
        ! and returning the C wrapper in its place, also if only the object provides one
        if (.not. held%level%extension%provides%update_orbs) then
            host_hess_x_c_funptr = c_null_funptr
            host_update_orbs => held%callbacks%update_orbs
            error_c = host_update_orbs(kappa_c, func_c, grad_c, h_diag_c, &
                                       host_hess_x_c_funptr, &
                                       held%callbacks%host_context)
            held%callbacks%hess_x => null()
            if (c_associated(host_hess_x_c_funptr)) call c_f_procpointer( &
                cptr=host_hess_x_c_funptr, fptr=held%callbacks%hess_x)
            if (associated(held%callbacks%stability)) &
                held%callbacks%stability%hess_x => held%callbacks%hess_x
            if (associated(held%callbacks%hess_x) .or. &
                held%level%extension%provides%hess_x) &
                hess_x_c_funptr = c_funloc(hess_x_c_wrapper)
            return
        end if

        ! hand the object the current state of the host's record
        call refresh_state(held)

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            kappa => kappa_c(:n_param)
            grad => grad_c(:n_param)
            h_diag => h_diag_c(:n_param)
        else
            allocate(kappa(n_param), grad(n_param), h_diag(n_param))
            kappa = real(kappa_c(:n_param), kind=rp)
        end if

        ! call the orbital update with the level, after disassociating the Hessian
        ! linear transformation so that a missing one is detected
        hess_x_funptr => null()
        call update_orbs_callback(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                  context)

        ! convert arguments to C kind, copy the orbitals back to the host's array if
        ! the object updates a copy of them and return the C wrapper of the Hessian
        ! linear transformation if the orbital update returned one
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            grad_c(:n_param) = real(grad, kind=c_rp)
            h_diag_c(:n_param) = real(h_diag, kind=c_rp)
            deallocate(kappa, grad, h_diag)
        end if
        if (associated(held%orbitals_c)) &
            held%orbitals_c = real(held%orbitals, kind=c_rp)
        if (associated(hess_x_funptr)) hess_x_c_funptr = c_funloc(hess_x_c_wrapper)

    end function update_orbs_c_wrapper

    recursive function obj_func_c_wrapper(kappa_c, func_c, context_c) result(error_c) &
        bind(C)
        !
        ! this function exposes the objective function to C hosts, calling it with the
        ! level of an extension object held at the address it receives as context, or
        ! forwarding to the host's C function in the bundle of that level, with its
        ! host context, for an object that does not perform it
        !
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(obj_func_c_type), pointer :: host_obj_func
        class(*), pointer :: context
        real(rp), pointer :: kappa(:)
        integer(ip) :: n_param, error

        ! find the level of the holder, without which there is nothing to call
        func_c = 0.0_c_rp
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if
        context => held%level
        n_param = held%level%extension%n_param

        ! forward straight to the host's C function with its host context unless the
        ! object performs the role itself
        if (.not. held%level%extension%provides%obj_func) then
            host_obj_func => held%callbacks%obj_func
            error_c = host_obj_func(kappa_c, func_c, held%callbacks%host_context)
            return
        end if

        ! hand the object the current state of the host's record
        call refresh_state(held)

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            kappa => kappa_c(:n_param)
        else
            allocate(kappa(n_param))
            kappa = real(kappa_c(:n_param), kind=rp)
        end if

        ! call the callback function with the level
        func_c = real(obj_func_callback(kappa, error, context), kind=c_rp)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) deallocate(kappa)

    end function obj_func_c_wrapper

    recursive function hess_x_c_wrapper(x_c, hess_x_c, context_c) result(error_c) &
        bind(C)
        !
        ! this function exposes the Hessian linear transformation to C hosts, calling
        ! it with the level of an extension object held at the address it receives as
        ! context, or forwarding to the host's C function in the bundle of that level,
        ! with its host context, for an object that does not perform it
        !
        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(hess_x_c_type), pointer :: host_hess_x
        class(*), pointer :: context
        real(rp), pointer :: x(:), hess_x(:)
        integer(ip) :: n_param, error

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if
        context => held%level
        n_param = held%level%extension%n_param

        ! forward straight to the host's C function with its host context unless the
        ! object performs the role itself
        if (.not. held%level%extension%provides%hess_x) then
            host_hess_x => held%callbacks%hess_x
            error_c = host_hess_x(x_c, hess_x_c, held%callbacks%host_context)
            return
        end if

        ! hand the object the current state of the host's record
        call refresh_state(held)

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            x => x_c(:n_param)
            hess_x => hess_x_c(:n_param)
        else
            allocate(x(n_param), hess_x(n_param))
            x = real(x_c(:n_param), kind=rp)
        end if

        ! call the callback function with the level
        call hess_x_callback(x, hess_x, error, context)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            hess_x_c(:n_param) = real(hess_x, kind=c_rp)
            deallocate(x, hess_x)
        end if

    end function hess_x_c_wrapper

    recursive function precond_c_wrapper(residual_c, mu_c, precond_residual_c, &
                                         context_c) result(error_c) bind(C)
        !
        ! this function exposes the level-shifted preconditioner to C hosts, calling it
        ! with the level of an extension object held at the address it receives as
        ! context, or forwarding to the host's C function in the bundle of that level,
        ! with its host context, for an object that does not perform it
        !
        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(in) :: mu_c
        real(c_rp), intent(out), target :: precond_residual_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(precond_c_type), pointer :: host_precond
        class(*), pointer :: context
        real(rp), pointer :: residual(:), precond_residual(:)
        integer(ip) :: n_param, error

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if
        context => held%level
        n_param = held%level%extension%n_param

        ! forward straight to the host's C function with its host context unless the
        ! object performs the role itself
        if (.not. held%level%extension%provides%precond) then
            host_precond => held%callbacks%precond
            error_c = host_precond(residual_c, mu_c, precond_residual_c, &
                                   held%callbacks%host_context)
            return
        end if

        ! hand the object the current state of the host's record
        call refresh_state(held)

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            residual => residual_c(:n_param)
            precond_residual => precond_residual_c(:n_param)
        else
            allocate(residual(n_param), precond_residual(n_param))
            residual = real(residual_c(:n_param), kind=rp)
        end if

        ! call the callback function with the level
        call precond_callback(residual, real(mu_c, kind=rp), precond_residual, error, &
                              context)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            precond_residual_c(:n_param) = real(precond_residual, kind=c_rp)
            deallocate(residual, precond_residual)
        end if

    end function precond_c_wrapper

    recursive function precond_pd_c_wrapper(residual_c, precond_residual_c, context_c) &
        result(error_c) bind(C)
        !
        ! this function exposes the positive-definite preconditioner to C hosts,
        ! calling it with the level of an extension object held at the address it
        ! receives as context, or forwarding to the host's C function in the bundle of
        ! that level, with its host context, for an object that does not perform it
        !
        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(out), target :: precond_residual_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(precond_pd_c_type), pointer :: host_precond_pd
        class(*), pointer :: context
        real(rp), pointer :: residual(:), precond_residual(:)
        integer(ip) :: n_param, error

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if
        context => held%level
        n_param = held%level%extension%n_param

        ! forward straight to the host's C function with its host context unless the
        ! object performs the role itself
        if (.not. held%level%extension%provides%precond_pd) then
            host_precond_pd => held%callbacks%precond_pd
            error_c = host_precond_pd(residual_c, precond_residual_c, &
                                      held%callbacks%host_context)
            return
        end if

        ! hand the object the current state of the host's record
        call refresh_state(held)

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            residual => residual_c(:n_param)
            precond_residual => precond_residual_c(:n_param)
        else
            allocate(residual(n_param), precond_residual(n_param))
            residual = real(residual_c(:n_param), kind=rp)
        end if

        ! call the callback function with the level
        call precond_pd_callback(residual, precond_residual, error, context)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            precond_residual_c(:n_param) = real(precond_residual, kind=c_rp)
            deallocate(residual, precond_residual)
        end if

    end function precond_pd_c_wrapper

    recursive function project_c_wrapper(vector_c, context_c) result(error_c) bind(C)
        !
        ! this function exposes the projection to C hosts, calling it with the level of
        ! an extension object held at the address it receives as context, or forwarding
        ! to the host's C function in the bundle of that level, with its host context,
        ! for an object that does not perform it
        !
        real(c_rp), intent(inout), target :: vector_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(project_c_type), pointer :: host_project
        class(*), pointer :: context
        real(rp), pointer :: vector(:)
        integer(ip) :: n_param, error

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if
        context => held%level
        n_param = held%level%extension%n_param

        ! forward straight to the host's C function with its host context unless the
        ! object performs the role itself
        if (.not. held%level%extension%provides%project) then
            host_project => held%callbacks%project
            error_c = host_project(vector_c, held%callbacks%host_context)
            return
        end if

        ! hand the object the current state of the host's record
        call refresh_state(held)

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            vector => vector_c(:n_param)
        else
            allocate(vector(n_param))
            vector = real(vector_c(:n_param), kind=rp)
        end if

        ! call the callback function with the level
        call project_callback(vector, error, context)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            vector_c(:n_param) = real(vector, kind=c_rp)
            deallocate(vector)
        end if

    end function project_c_wrapper

    recursive function get_extra_trial_vectors_c_wrapper( &
        trial_vectors_c, n_extra_trial_vectors_c, context_c) result(error_c) bind(C)
        !
        ! this function exposes the extra trial vector function to C hosts, calling it
        ! with the level of an extension object held at the address it receives as
        ! context, or forwarding to the host's C function in the bundle of that level,
        ! with its host context, for an object that does not perform it
        !
        real(c_rp), intent(out), target :: trial_vectors_c(*)
        integer(c_ip), intent(in), value :: n_extra_trial_vectors_c
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(get_extra_trial_vectors_c_type), pointer :: &
            host_get_extra_trial_vectors
        class(*), pointer :: context
        real(rp), pointer :: trial_vectors(:, :)
        integer(ip) :: n_param, n_extra_trial_vectors, error

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if
        context => held%level
        n_param = held%level%extension%n_param

        ! forward straight to the host's C function with its host context unless the
        ! object performs the role itself
        if (.not. held%level%extension%provides%get_extra_trial_vectors) then
            host_get_extra_trial_vectors => held%callbacks%get_extra_trial_vectors
            error_c = host_get_extra_trial_vectors( &
                trial_vectors_c, n_extra_trial_vectors_c, held%callbacks%host_context)
            return
        end if
        n_extra_trial_vectors = int(n_extra_trial_vectors_c, kind=ip)

        ! hand the object the current state of the host's record
        call refresh_state(held)

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            trial_vectors(1:n_param, 1:n_extra_trial_vectors) => &
                trial_vectors_c(:n_param * n_extra_trial_vectors)
        else
            allocate(trial_vectors(n_param, n_extra_trial_vectors))
        end if

        ! call the callback function with the level
        call get_extra_trial_vectors_callback(trial_vectors, error, context)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            trial_vectors_c(:n_param * n_extra_trial_vectors) = reshape( &
                real(trial_vectors, kind=c_rp), [n_param * n_extra_trial_vectors])
            deallocate(trial_vectors)
        end if

    end function get_extra_trial_vectors_c_wrapper

    recursive function conv_check_c_wrapper(converged_c, context_c) result(error_c) &
        bind(C)
        !
        ! this function forwards the convergence check of the solver or of its
        ! stability check for C hosts straight to the host's C function in the bundle
        ! of the level of a holder at the address it receives as context, with the host
        ! context of the level
        !
        logical(c_bool), intent(out) :: converged_c
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held

        ! find the level of the holder, without which there is nothing to call
        converged_c = .false._c_bool
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level
        error_c = held%callbacks%conv_check(converged_c, held%callbacks%host_context)

    end function conv_check_c_wrapper

    recursive function stability_hess_x_c_wrapper(x_c, hess_x_c, context_c) &
        result(error_c) bind(C)
        !
        ! this function forwards the Hessian linear transformation of the stability
        ! check for C hosts straight to the host's C function in the bundle of the
        ! level of a holder at the address it receives as context, with the host
        ! context of the level
        !
        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(hess_x_c_type), pointer :: host_stability_hess_x

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_stability_hess_x => held%callbacks%stability_hess_x
        error_c = host_stability_hess_x(x_c, hess_x_c, held%callbacks%host_context)

    end function stability_hess_x_c_wrapper

    recursive function approx_hess_x_c_wrapper(x_c, hess_x_c, context_c) &
        result(error_c) bind(C)
        !
        ! this function forwards the approximate Hessian linear transformation for C
        ! hosts straight to the host's C function in the bundle of the level of a
        ! holder at the address it receives as context, with the host context of the
        ! level
        !
        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(hess_x_c_type), pointer :: host_approx_hess_x

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_approx_hess_x => held%callbacks%approx_hess_x
        error_c = host_approx_hess_x(x_c, hess_x_c, held%callbacks%host_context)

    end function approx_hess_x_c_wrapper

    recursive subroutine logger_c_wrapper(message_c, context_c) bind(C)
        !
        ! this subroutine forwards the logging function for C hosts straight to the
        ! host's C function in the bundle of the level of a holder at the address it
        ! receives as context, with the host context of the level
        !
        character(kind=c_char), intent(in) :: message_c(*)
        type(c_ptr), intent(in), value :: context_c

        type(holder_level_type), pointer :: held
        procedure(logger_c_type), pointer :: host_logger

        ! find the level of the holder, without which the message is discarded
        held => holder_level(context_c)
        if (.not. associated(held)) return

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_logger => held%callbacks%logger
        call host_logger(message_c, held%callbacks%host_context)

    end subroutine logger_c_wrapper

    subroutine evaluate_dm_cs_f_wrapper(dm, energy, fock, get_response_cs_funptr, &
                                        error, context)
        !
        ! this subroutine wraps the density matrix evaluating C function of an orbital
        ! basis to convert Fortran variables to C variables for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target, contiguous :: fock(:, :)
        procedure(get_response_cs_type), intent(out), optional, pointer :: &
            get_response_cs_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp), pointer :: dm_3d(:, :, :), fock_3d(:, :, :)
        procedure(get_response_os_type), pointer :: get_response_funptr

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        nullify(fock_3d, get_response_funptr)
        if (present(fock)) fock_3d(1:size(fock, 1), 1:size(fock, 2), 1:1) => fock
        if (present(get_response_cs_funptr)) then
            call evaluate_dm_os_f_wrapper(dm_3d, energy, fock_3d, get_response_funptr, &
                                          error, context)
            get_response_cs_funptr => null()
            if (error == 0) get_response_cs_funptr => get_response_cs_f_wrapper
        else
            call evaluate_dm_os_f_wrapper(dm_3d, energy, fock_3d, error=error, &
                                          context=context)
        end if

    end subroutine evaluate_dm_cs_f_wrapper

    subroutine evaluate_dm_os_f_wrapper(dm, energy, fock, get_response_funptr, error, &
                                        context)
        !
        ! this subroutine wraps the density matrix evaluating C function of an orbital
        ! basis to convert Fortran variables to C variables for the open-shell case,
        ! finding the C function and the host context in the bundle it receives as
        ! context, in which it keeps the response C function the C function returns if
        ! one is requested
        !
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target :: fock(:, :, :)
        procedure(get_response_os_type), intent(out), optional, pointer :: &
            get_response_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :, :), fock_c(:, :, :)
        type(c_funptr) :: get_response_c_funptr
        integer(c_ip) :: error_c

        ! the bundle holds the C function
        energy = 0.0_rp
        error = 1
        if (present(get_response_funptr)) get_response_funptr => null()
        if (.not. associated(context)) return
        select type (callbacks => context)
        type is (orbital_c_callbacks_type)
            ! convert arguments to C kind
            nullify(fock_c)
            if (rp == c_rp) then
                dm_c => dm
                if (present(fock)) fock_c => fock
            else
                allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)))
                dm_c = real(dm, kind=c_rp)
                if (present(fock)) &
                    allocate(fock_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            end if

            ! call density matrix evaluating C function and keep the response C
            ! function it returns in the bundle
            if (present(get_response_funptr)) then
                error_c = callbacks%evaluate_dm(dm_c, energy_c, fock_c, &
                                                get_response_c_funptr, &
                                                callbacks%host_context)
                if (error_c == 0) call c_f_procpointer(cptr=get_response_c_funptr, &
                                                       fptr=callbacks%get_response)
            else
                error_c = callbacks%evaluate_dm(dm_c, energy_c, fock_c, &
                                                context_c=callbacks%host_context)
            end if

            ! convert arguments to Fortran kind
            energy = real(energy_c, kind=rp)
            error = int(error_c, kind=ip)
            if (rp /= c_rp) then
                if (present(fock)) then
                    fock = real(fock_c, kind=rp)
                    deallocate(fock_c)
                end if
                deallocate(dm_c)
            end if
        end select

        ! associate the Fortran procedure pointer with the wrapper of the response
        ! function
        if (present(get_response_funptr) .and. error == 0) &
            get_response_funptr => get_response_os_f_wrapper

    end subroutine evaluate_dm_os_f_wrapper

    subroutine get_response_cs_f_wrapper(dm, response, error, context)
        !
        ! this subroutine wraps the response C function of an orbital basis to convert
        ! Fortran variables to C variables for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out), target, contiguous :: response(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp), pointer :: dm_3d(:, :, :), response_3d(:, :, :)

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        response_3d(1:size(response, 1), 1:size(response, 2), 1:1) => response
        call get_response_os_f_wrapper(dm_3d, response_3d, error, context)

    end subroutine get_response_cs_f_wrapper

    subroutine get_response_os_f_wrapper(dm, response, error, context)
        !
        ! this subroutine wraps the response C function of an orbital basis to convert
        ! Fortran variables to C variables for the open-shell case, finding the C
        ! function the density matrix evaluating function returned last and the host
        ! context in the bundle it receives as context
        !
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out), target :: response(:, :, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp), pointer :: dm_c(:, :, :), response_c(:, :, :)
        integer(c_ip) :: error_c

        ! the bundle holds the C function
        error = 1
        if (.not. associated(context)) return
        select type (callbacks => context)
        type is (orbital_c_callbacks_type)
            ! convert arguments to C kind
            if (rp == c_rp) then
                dm_c => dm
                response_c => response
            else
                allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)), &
                         response_c(size(dm, 1), size(dm, 2), size(dm, 3)))
                dm_c = real(dm, kind=c_rp)
            end if

            ! call response C function
            error_c = callbacks%get_response(dm_c, response_c, callbacks%host_context)

            ! convert arguments to Fortran kind
            error = int(error_c, kind=ip)
            if (rp /= c_rp) then
                response = real(response_c, kind=rp)
                deallocate(dm_c, response_c)
            end if
        end select

    end subroutine get_response_os_f_wrapper

end module otr_common_c_interface
