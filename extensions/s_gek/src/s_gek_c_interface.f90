! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface

    use opentrustregion_core, only: ip, rp, kw_len
    use c_interface, only: c_ip, c_rp
    use otr_s_gek, only: change_reference_type, standard_s_gek_solver => s_gek_solver
    use otr_common_c_interface, only: obj_func_grad_c_callbacks_type
    use, intrinsic :: iso_c_binding, only: c_size_t, c_sizeof, c_bool, c_ptr, c_char, &
                                           c_funptr, c_f_procpointer, c_f_pointer, &
                                           c_associated, c_null_funptr, c_null_ptr

    implicit none

    ! C-interoperable interface for the chart transformation callback function, while
    ! otr_common_c_interface holds that of the objective function that also returns the
    ! gradient at a trial point
    abstract interface
        function change_reference_c_type(new_ref_c, n_points_c, kappa_list_c, &
                                         local_grad_list_c, grad_list_c, context_c) &
            result(error_c) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: new_ref_c(*)
            integer(c_ip), intent(in), value :: n_points_c
            real(c_rp), intent(inout), target :: kappa_list_c(*), &
                                                 local_grad_list_c(*), grad_list_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error_c
        end function change_reference_c_type
    end interface

    ! bundle of the C callback functions of the S-GEK settings, which holds the static
    ! Hessian in the slot of the core's Hessian linear transformation, whose Fortran
    ! wrapper the converted settings use, and extends the shared bundle with the
    ! objective function that also returns the gradient at a trial point by the chart
    ! transformation, reached through the context of the converted S-GEK settings
    type, extends(obj_func_grad_c_callbacks_type) :: s_gek_c_callbacks_type
        procedure(change_reference_c_type), pointer, nopass :: change_reference => &
            null()
    end type

    ! derived type for S-GEK settings
    type, bind(C) :: s_gek_settings_type_c
        type(c_funptr) :: obj_func_grad, init_hess, logger
        logical(c_bool) :: initialized, excess_dispersion
        integer(c_ip) :: verbose
        character(kind=c_char) :: inner_solver(kw_len + 1), prior_update(kw_len + 1), &
                                  acceptance(kw_len + 1)
        type(c_ptr) :: context
    end type

    procedure(standard_s_gek_solver), pointer :: s_gek_solver => standard_s_gek_solver

    ! create function pointers to ensure that routines comply with interface
    procedure(change_reference_type), pointer :: change_reference_f_wrapper_ptr => &
        change_reference_f_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_s_gek_f_c
        module procedure assign_s_gek_c_f
    end interface

contains

    function s_gek_solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                    change_reference_c_funptr, n_param_c, settings_c, &
                                    state_c_ptr, s_gek_settings_c) result(error_c) &
        bind(C, name="s_gek_solver")
        !
        ! this function exposes the S-GEK solver to C as the solver's entry point
        ! exposes the solver, with the chart transformation and the S-GEK settings,
        ! whose callback functions and host context a bundle of their own holds
        ! together with the chart transformation, which refers to the state like that
        ! of the solver settings
        !
        use opentrustregion_core, only: solver_settings_type, solver_state_type, &
                                        default_solver_settings, &
                                        warn_uninitialized_settings, verbosity_error
        use c_interface, only: solver_settings_type_c, c_callbacks_type, &
                               convert_solver_settings_c, copy_state_to_c, &
                               update_orbs_f_wrapper, obj_func_f_wrapper, assignment(=)
        use otr_s_gek, only: s_gek_settings_type, default_s_gek_settings, &
                             s_gek_settings_uninitialized_warning_msg, error_s_gek
        use otr_common_c_interface, only: bundle_settings_c

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             change_reference_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(solver_settings_type_c), intent(in), optional :: settings_c
        type(c_ptr), intent(in), value :: state_c_ptr
        type(s_gek_settings_type_c), intent(in), optional :: s_gek_settings_c
        integer(c_ip) :: error_c

        procedure(update_orbs_f_wrapper), pointer :: update_orbs
        procedure(obj_func_f_wrapper), pointer :: obj_func
        procedure(change_reference_type), pointer :: change_reference
        integer(ip) :: n_param, error
        type(solver_settings_type) :: settings
        type(s_gek_settings_type) :: s_gek_settings
        type(solver_settings_type_c) :: local_settings_c
        type(s_gek_settings_type_c) :: local_s_gek_settings_c
        type(solver_state_type), target :: state
        type(c_callbacks_type), target :: callbacks, stability_callbacks
        type(s_gek_c_callbacks_type), target :: s_gek_callbacks

        ! work on copies of the given settings or, without any, of the default settings
        if (present(settings_c)) then
            local_settings_c = settings_c
        else
            local_settings_c = default_solver_settings
        end if
        if (present(s_gek_settings_c)) then
            local_s_gek_settings_c = s_gek_settings_c
        else
            local_s_gek_settings_c = default_s_gek_settings
        end if

        ! bundle the C function pointers, the chart transformation with the callback
        ! functions of the S-GEK settings
        call c_f_procpointer(cptr=update_orbs_c_funptr, fptr=callbacks%update_orbs)
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=callbacks%obj_func)
        call c_f_procpointer(cptr=change_reference_c_funptr, &
                             fptr=s_gek_callbacks%change_reference)

        ! let the bundles refer to the state and to the host's record, if one is given,
        ! so that the state is copied into it before every callback function
        if (c_associated(state_c_ptr)) then
            callbacks%solver_state => state
            call c_f_pointer(cptr=state_c_ptr, fptr=callbacks%solver_state_c)
            s_gek_callbacks%solver_state => state
            s_gek_callbacks%solver_state_c => callbacks%solver_state_c
        end if

        ! associate procedure pointer to wrapper function
        update_orbs => update_orbs_f_wrapper
        obj_func => obj_func_f_wrapper
        change_reference => change_reference_f_wrapper

        ! convert dummy argument to Fortran kind
        n_param = int(n_param_c, kind=ip)

        ! convert solver settings, whose optional callback functions and host contexts
        ! join the bundles of the solver and of the internal stability check
        call convert_solver_settings_c(local_settings_c, settings, callbacks, &
                                       stability_callbacks)

        ! convert the S-GEK settings, which replaces settings that were not initialized
        ! by the default settings, and print the warning for these on stderr, since
        ! they provide no logging function
        s_gek_settings = local_s_gek_settings_c
        if (.not. local_s_gek_settings_c%initialized) &
            call warn_uninitialized_settings(s_gek_settings, &
                                             s_gek_settings_uninitialized_warning_msg)

        ! hand the objective function, the static Hessian and the logging function with
        ! the host context of initialized settings to their bundle, which the settings
        ! point to
        if (local_s_gek_settings_c%initialized) then
            if (c_associated(local_s_gek_settings_c%obj_func_grad)) &
                call c_f_procpointer(cptr=local_s_gek_settings_c%obj_func_grad, &
                                     fptr=s_gek_callbacks%obj_func_grad)
            if (c_associated(local_s_gek_settings_c%init_hess)) call c_f_procpointer( &
                cptr=local_s_gek_settings_c%init_hess, fptr=s_gek_callbacks%hess_x)
        end if
        call bundle_settings_c(s_gek_callbacks, local_s_gek_settings_c%initialized, &
                               local_s_gek_settings_c%logger, &
                               local_s_gek_settings_c%context)
        s_gek_settings%context => s_gek_callbacks

        ! the solver can not run without an orbital update, an objective function and a
        ! chart transformation
        if (.not. (c_associated(update_orbs_c_funptr) .and. &
                   c_associated(obj_func_c_funptr) .and. &
                   c_associated(change_reference_c_funptr))) then
            call settings%log("Orbital update, objective function and chart "// &
                              "transformation have to be given.", verbosity_error, &
                              .true.)
            call copy_state_to_c(callbacks)
            error_c = int(error_s_gek + 1, kind=c_ip)
            return
        end if

        ! call solver
        call s_gek_solver(update_orbs, obj_func, change_reference, n_param, error, &
                          settings, state, s_gek_settings)

        ! copy the state into the host's record once more, also after an error
        call copy_state_to_c(callbacks)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function s_gek_solver_c_wrapper

    subroutine change_reference_f_wrapper(new_ref, kappa_list, local_grad_list, &
                                          grad_list, error, context)
        !
        ! this subroutine exposes the chart transformation C function to Fortran,
        ! finding the C function and the host context in the bundle it receives as
        ! context, which it copies the state into first
        !
        use c_interface, only: copy_state_to_c

        real(rp), intent(in), target :: new_ref(:)
        real(rp), intent(inout), target :: kappa_list(:, :), local_grad_list(:, :), &
                                           grad_list(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(c_ip) :: error_c

        ! the bundle holds the C function
        error = 1
        if (.not. associated(context)) return
        select type (callbacks => context)
        type is (s_gek_c_callbacks_type)
            ! call chart transformation C function
            call copy_state_to_c(callbacks%c_callbacks_type)
            error_c = callbacks%change_reference( &
                new_ref, int(size(kappa_list, 2), kind=c_ip), kappa_list, &
                local_grad_list, grad_list, callbacks%host_context)

            ! convert error code to Fortran kind
            error = int(error_c, kind=ip)
        end select

    end subroutine change_reference_f_wrapper

    subroutine init_s_gek_settings_c(settings_c) bind(C, name="init_s_gek_settings")
        !
        ! this subroutine initializes the S-GEK settings
        !
        use otr_s_gek, only: default_s_gek_settings

        type(s_gek_settings_type_c), intent(out) :: settings_c

        settings_c = default_s_gek_settings

    end subroutine init_s_gek_settings_c

    integer(c_size_t) function s_gek_settings_size() bind(C)
        !
        ! this function returns the size of the C S-GEK settings, which the Python
        ! interface compares with its own to detect a library of another version and
        ! the C layout tests with the header's
        !
        type(s_gek_settings_type_c) :: settings_c

        s_gek_settings_size = c_sizeof(settings_c)

    end function s_gek_settings_size

    subroutine assign_s_gek_f_c(settings, settings_c)
        !
        ! this subroutine converts S-GEK settings from C to Fortran
        !
        use otr_s_gek, only: s_gek_settings_type, default_s_gek_settings
        use c_interface, only: character_from_c, logger_f_wrapper, hess_x_f_wrapper
        use otr_common_c_interface, only: obj_func_grad_f_wrapper

        type(s_gek_settings_type), intent(out) :: settings
        type(s_gek_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%obj_func_grad)) then
                settings%obj_func_grad => obj_func_grad_f_wrapper
            else
                settings%obj_func_grad => null()
            end if
            if (c_associated(settings_c%init_hess)) then
                settings%init_hess => hess_x_f_wrapper
            else
                settings%init_hess => null()
            end if
            if (c_associated(settings_c%logger)) then
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert logicals
            settings%excess_dispersion = logical(settings_c%excess_dispersion)

            ! convert integers
            settings%verbose = int(settings_c%verbose, kind=ip)

            ! convert characters
            settings%inner_solver = character_from_c(settings_c%inner_solver)
            settings%prior_update = character_from_c(settings_c%prior_update)
            settings%acceptance = character_from_c(settings_c%acceptance)
        else
            ! C settings that were not initialized hold no values and are converted as
            ! the default settings
            settings = default_s_gek_settings
        end if

    end subroutine assign_s_gek_f_c

    subroutine assign_s_gek_c_f(settings_c, settings)
        !
        ! this subroutine converts S-GEK settings from Fortran to C
        !
        use otr_s_gek, only: s_gek_settings_type
        use c_interface, only: character_to_c

        type(s_gek_settings_type_c), intent(out) :: settings_c
        type(s_gek_settings_type), intent(in) :: settings

        ! callback functions cannot be converted
        settings_c%obj_func_grad = c_null_funptr
        settings_c%init_hess = c_null_funptr
        settings_c%logger = c_null_funptr

        ! convert logicals
        settings_c%excess_dispersion = logical(settings%excess_dispersion, kind=c_bool)

        ! convert integers
        settings_c%verbose = int(settings%verbose, kind=c_ip)

        ! convert characters
        settings_c%inner_solver = character_to_c(settings%inner_solver)
        settings_c%prior_update = character_to_c(settings%prior_update)
        settings_c%acceptance = character_to_c(settings%acceptance)

        ! the host context cannot be converted
        settings_c%context = c_null_ptr

        ! Fortran settings are always initialized
        settings_c%initialized = .true._c_bool

    end subroutine assign_s_gek_c_f

end module otr_s_gek_c_interface
