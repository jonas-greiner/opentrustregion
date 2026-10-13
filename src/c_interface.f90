! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface

    use opentrustregion, only: &
        rp, ip, kw_len, standard_solver => solver, &
        standard_stability_check => stability_check, default_solver_settings, &
        default_stability_settings, update_orbs_type, hess_x_type, obj_func_type, &
        precond_type, precond_pd_type, project_type, conv_check_type, &
        get_extra_trial_vectors_type, conv_check_stability_type, logger_type, &
        solver_state_type, stability_state_type
    use, intrinsic :: iso_c_binding, only: &
        c_double, c_int64_t, c_int32_t, c_bool, c_ptr, c_funptr, c_f_pointer, &
        c_f_procpointer, c_associated, c_char, c_null_char, c_null_ptr, c_null_funptr

    implicit none

    integer, parameter :: c_rp = c_double
#ifdef USE_ILP64
    integer, parameter :: c_ip = c_int64_t  ! 64-bit integers
    logical(c_bool), bind(C) :: ilp64 = .true._c_bool
#else
    integer, parameter :: c_ip = c_int32_t  ! 32-bit integers
    logical(c_bool), bind(C) :: ilp64 = .false._c_bool
#endif

    ! maximum keyword length
    integer(c_ip), bind(C) :: kw_len_c = kw_len

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function update_orbs_c_type(kappa_c, func_c, grad_c, h_diag_c, &
                                    hess_x_c_funptr, context_c) result(error) bind(C)
            import :: c_rp, c_funptr, c_ip, c_ptr

            real(c_rp), intent(in), target :: kappa_c(*)
            real(c_rp), intent(out) :: func_c
            real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
            type(c_funptr), intent(inout) :: hess_x_c_funptr
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function update_orbs_c_type
    end interface

    abstract interface
        function hess_x_c_type(x_c, hess_x_c, context_c) result(error) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: x_c(*)
            real(c_rp), intent(out), target :: hess_x_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function hess_x_c_type
    end interface

    abstract interface
        function obj_func_c_type(kappa_c, func, context_c) result(error) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: kappa_c(*)
            real(c_rp), intent(out) :: func
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function obj_func_c_type
    end interface

    abstract interface
        function precond_c_type(residual_c, mu_c, precond_residual_c, context_c) &
            result(error) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: residual_c(*)
            real(c_rp), intent(in) :: mu_c
            real(c_rp), intent(out), target :: precond_residual_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function precond_c_type
    end interface

    abstract interface
        function precond_pd_c_type(residual_c, precond_residual_c, context_c) &
            result(error) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: residual_c(*)
            real(c_rp), intent(out), target :: precond_residual_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function precond_pd_c_type
    end interface

    abstract interface
        function project_c_type(vector_c, context_c) result(error) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(inout), target :: vector_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function project_c_type
    end interface

    abstract interface
        function conv_check_c_type(converged, context_c) result(error) bind(C)
            import :: c_bool, c_ip, c_ptr

            logical(c_bool), intent(out) :: converged
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function conv_check_c_type
    end interface

    abstract interface
        function get_extra_trial_vectors_c_type( &
            trial_vectors_c, n_extra_trial_vectors_c, context_c) result(error) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(out), target :: trial_vectors_c(*)
            integer(c_ip), intent(in), value :: n_extra_trial_vectors_c
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function get_extra_trial_vectors_c_type
    end interface

    abstract interface
        function conv_check_stability_c_type(residual_c, eigval_c, converged, &
                                             context_c) result(error) bind(C)
            import :: c_bool, c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: residual_c(*)
            real(c_rp), intent(in) :: eigval_c
            logical(c_bool), intent(out) :: converged
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error
        end function conv_check_stability_c_type
    end interface

    abstract interface
        subroutine logger_c_type(message, context_c) bind(C)
            import :: c_char, c_ptr

            character(kind=c_char), intent(in) :: message(*)
            type(c_ptr), intent(in), value :: context_c
        end subroutine logger_c_type
    end interface

    ! derived type for the stability check state
    type, bind(C) :: stability_state_type_c
        real(c_rp) :: min_eigval, residual_rms
        integer(c_ip) :: iteration, n_hess_x
    end type

    ! derived type for the solver state
    type, bind(C) :: solver_state_type_c
        logical(c_bool) :: max_precision_reached
        real(c_rp) :: func, grad_rms, trust_radius, relative_residual, level_shift, &
                      step_size
        integer(c_ip) :: macro_iteration, micro_iteration, jacobi_davidson_iteration, &
                         n_update_orbs, n_obj_func, n_hess_x
        type(stability_state_type_c) :: stability_state
    end type

    ! bundle of the C callback functions belonging to a single call, together with
    ! the host context and, if the host gave a record of the state, the state of the
    ! call and this record, into which the state is copied
    type :: c_callbacks_type
        procedure(update_orbs_c_type), pointer, nopass :: update_orbs => null()
        procedure(hess_x_c_type), pointer, nopass :: hess_x => null()
        procedure(obj_func_c_type), pointer, nopass :: obj_func => null()
        procedure(precond_c_type), pointer, nopass :: precond => null()
        procedure(precond_pd_c_type), pointer, nopass :: precond_pd => null()
        procedure(project_c_type), pointer, nopass :: project => null()
        procedure(get_extra_trial_vectors_c_type), pointer, nopass :: &
            get_extra_trial_vectors => null()
        procedure(conv_check_c_type), pointer, nopass :: conv_check => null()
        procedure(conv_check_stability_c_type), pointer, nopass :: &
            conv_check_stability => null()
        procedure(hess_x_c_type), pointer, nopass :: stability_hess_x => null()
        procedure(hess_x_c_type), pointer, nopass :: approx_hess_x => null()
        procedure(logger_c_type), pointer, nopass :: logger => null()
        type(c_ptr) :: host_context = c_null_ptr
        type(c_callbacks_type), pointer :: stability => null()
        type(solver_state_type), pointer :: solver_state => null()
        type(solver_state_type_c), pointer :: solver_state_c => null()
        type(stability_state_type), pointer :: stability_state => null()
        type(stability_state_type_c), pointer :: stability_state_c => null()
    end type

    ! derived type for stability check settings
    type, bind(C) :: stability_settings_type_c
        type(c_funptr) :: precond, project, approx_hess_x, get_extra_trial_vectors, &
                          conv_check, logger
        logical(c_bool) :: hess_symm, stop_on_instability, initialized
        real(c_rp) :: conv_tol
        integer(c_ip) :: n_random_trial_vectors, n_extra_trial_vectors, n_iter, &
                         jacobi_davidson_start, seed, verbose
        character(kind=c_char) :: diag_solver(kw_len + 1)
        type(c_ptr) :: context
    end type

    ! derived type for solver settings
    type, bind(C) :: solver_settings_type_c
        type(c_funptr) :: precond, precond_pd, project, get_extra_trial_vectors, &
                          conv_check, stability_hess_x, logger
        logical(c_bool) :: stability, line_search, refresh_hess, grad_noise, &
                           hess_symm, initialized
        real(c_rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(c_ip) :: n_random_trial_vectors, n_extra_trial_vectors, n_macro, &
                         n_micro, jacobi_davidson_start, seed, verbose
        character(kind=c_char) :: subsystem_solver(kw_len + 1)
        character(kind=c_char) :: trust_region_shape(kw_len + 1)
        type(stability_settings_type_c) :: stability_settings
        type(c_ptr) :: context
    end type

    procedure(standard_solver), pointer :: solver => standard_solver
    procedure(standard_stability_check), pointer :: stability_check => &
        standard_stability_check

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_f_wrapper_ptr => &
        update_orbs_f_wrapper
    procedure(hess_x_type), pointer :: hess_x_f_wrapper_ptr => hess_x_f_wrapper
    procedure(obj_func_type), pointer :: obj_func_f_wrapper_ptr => obj_func_f_wrapper
    procedure(precond_type), pointer :: precond_f_wrapper_ptr => precond_f_wrapper
    procedure(precond_pd_type), pointer :: precond_pd_f_wrapper_ptr => &
        precond_pd_f_wrapper
    procedure(project_type), pointer :: project_f_wrapper_ptr => project_f_wrapper
    procedure(conv_check_type), pointer :: conv_check_f_wrapper_ptr => &
        conv_check_f_wrapper
    procedure(hess_x_type), pointer :: stability_hess_x_f_wrapper_ptr => &
        stability_hess_x_f_wrapper
    procedure(hess_x_type), pointer :: approx_hess_x_f_wrapper_ptr => &
        approx_hess_x_f_wrapper
    procedure(get_extra_trial_vectors_type), pointer :: &
        get_extra_trial_vectors_f_wrapper_ptr => get_extra_trial_vectors_f_wrapper
    procedure(conv_check_stability_type), pointer :: &
        conv_check_stability_f_wrapper_ptr => conv_check_stability_f_wrapper
    procedure(logger_type), pointer :: logger_f_wrapper_ptr => logger_f_wrapper

    ! interfaces for converting settings between C and Fortran and states from
    ! Fortran to C
    interface assignment(=)
        module procedure assign_solver_f_c
        module procedure assign_stability_f_c
        module procedure assign_solver_c_f
        module procedure assign_stability_c_f
        module procedure assign_solver_state_c_f
        module procedure assign_stability_state_c_f
    end interface

contains

    function solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, n_param_c, &
                              settings_c, state_c_ptr) result(error_c) &
        bind(C, name="solver")
        !
        ! this function exposes a Fortran-implemented solver subroutine to C
        !
        use opentrustregion, only: solver_settings_type

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(solver_settings_type_c), intent(in), optional :: settings_c
        type(c_ptr), intent(in), value :: state_c_ptr
        integer(c_ip) :: error_c

        procedure(update_orbs_f_wrapper), pointer :: update_orbs
        procedure(obj_func_f_wrapper), pointer :: obj_func
        integer(ip) :: n_param, error
        type(solver_settings_type) :: settings
        type(solver_settings_type_c) :: local_settings_c
        type(solver_state_type), target :: state
        type(c_callbacks_type), target :: callbacks, stability_callbacks

        ! work on a copy of the given settings or, without any, of the default
        ! settings
        if (present(settings_c)) then
            local_settings_c = settings_c
        else
            local_settings_c = default_solver_settings
        end if

        ! bundle the C function pointers
        call c_f_procpointer(cptr=update_orbs_c_funptr, fptr=callbacks%update_orbs)
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=callbacks%obj_func)

        ! let the bundle refer to the state and to the host's record, if one is given,
        ! so that the state is copied into it before every callback function
        if (c_associated(state_c_ptr)) then
            callbacks%solver_state => state
            call c_f_pointer(cptr=state_c_ptr, fptr=callbacks%solver_state_c)
        end if

        ! associate procedure pointer to wrapper function
        update_orbs => update_orbs_f_wrapper
        obj_func => obj_func_f_wrapper

        ! convert dummy argument to Fortran kind
        n_param = int(n_param_c, kind=ip)

        ! convert settings, whose optional callback functions and host contexts join
        ! the bundles of the solver and of the internal stability check
        call convert_solver_settings_c(local_settings_c, settings, callbacks, &
                                       stability_callbacks)

        ! call solver
        call solver(update_orbs, obj_func, n_param, error, settings, state)

        ! copy the state into the host's record once more, also after an error
        call copy_state_to_c(callbacks)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function solver_c_wrapper

    function stability_check_c_wrapper( &
        h_diag_c, hess_x_c_funptr, n_param_c, stable_c, settings_c, state_c_ptr, &
        kappa_c_ptr, min_eigval_c_ptr) result(error_c) bind(C, name="stability_check")
        !
        ! this function exposes a Fortran-implemented stability check subroutine to C
        !
        use opentrustregion, only: stability_settings_type, stderr, &
                                   settings_uninitialized_warning_msg

        real(c_rp), intent(in), target :: h_diag_c(*)
        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: hess_x_c_funptr
        logical(c_bool), intent(out) :: stable_c
        type(stability_settings_type_c), intent(in), optional :: settings_c
        type(c_ptr), intent(in), value :: state_c_ptr, kappa_c_ptr, min_eigval_c_ptr
        integer(c_ip) :: error_c

        real(rp), pointer :: h_diag_ptr(:), kappa_ptr(:), min_eigval_ptr
        real(c_rp), pointer :: kappa_ptr_c(:), min_eigval_ptr_c
        logical :: stable
        integer(ip) :: error
        procedure(hess_x_f_wrapper), pointer :: hess_x
        type(stability_settings_type) :: settings
        type(stability_settings_type_c) :: local_settings_c
        type(stability_state_type), target :: state
        type(c_callbacks_type), target :: callbacks

        ! work on a copy of the given settings or, without any, of the default
        ! settings
        if (present(settings_c)) then
            local_settings_c = settings_c
        else
            local_settings_c = default_stability_settings
        end if

        ! bundle the C function pointers
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=callbacks%hess_x)
        call store_optional_c_callbacks( &
            callbacks, local_settings_c%initialized, local_settings_c%precond, &
            c_null_funptr, local_settings_c%project, &
            local_settings_c%get_extra_trial_vectors, c_null_funptr, &
            local_settings_c%conv_check, c_null_funptr, &
            local_settings_c%approx_hess_x, local_settings_c%logger, &
            local_settings_c%context)

        ! let the bundle refer to the state and to the host's record, if one is given,
        ! so that the state is copied into it before every callback function
        if (c_associated(state_c_ptr)) then
            callbacks%stability_state => state
            call c_f_pointer(cptr=state_c_ptr, fptr=callbacks%stability_state_c)
        end if

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_f_wrapper

        ! convert C pointers
        if (rp == c_rp) then
            h_diag_ptr => h_diag_c(:n_param_c)
        else
            allocate(h_diag_ptr(n_param_c))
            h_diag_ptr = real(h_diag_c(:n_param_c), kind=rp)
        end if
        nullify(kappa_ptr, min_eigval_ptr)
        if (c_associated(kappa_c_ptr)) then
            if (rp == c_rp) then
                call c_f_pointer(cptr=kappa_c_ptr, fptr=kappa_ptr, shape=[n_param_c])
            else
                call c_f_pointer(cptr=kappa_c_ptr, fptr=kappa_ptr_c, shape=[n_param_c])
                allocate(kappa_ptr(n_param_c))
            end if
        end if
        if (c_associated(min_eigval_c_ptr)) then
            if (rp == c_rp) then
                call c_f_pointer(cptr=min_eigval_c_ptr, fptr=min_eigval_ptr)
            else
                call c_f_pointer(cptr=min_eigval_c_ptr, fptr=min_eigval_ptr_c)
                allocate(min_eigval_ptr)
            end if
        end if

        ! convert settings and hand the bundle to the callback wrappers as opaque
        ! context
        settings = local_settings_c
        settings%context => callbacks

        ! the conversion replaces settings that were not initialized by the default
        ! settings, so the stability check cannot detect these anymore and the warning
        ! is printed here on stderr regardless of the verbosity since settings that
        ! were not initialized provide no logging function
        if (.not. local_settings_c%initialized) &
            write(stderr, '(A)') " "//settings_uninitialized_warning_msg

        ! call stability check
        call stability_check(h_diag_ptr, hess_x, stable, error, settings, state, &
                             kappa=kappa_ptr, min_eigval=min_eigval_ptr)

        ! copy the state into the host's record once more, also after an error
        call copy_state_to_c(callbacks)

        ! convert return arguments to C kind
        if (rp /= c_rp) then
            if (associated(kappa_ptr)) then
                kappa_ptr_c = real(kappa_ptr, kind=c_rp)
                deallocate(kappa_ptr)
            end if
            if (associated(min_eigval_ptr)) then
                min_eigval_ptr_c = real(min_eigval_ptr, kind=c_rp)
                deallocate(min_eigval_ptr)
            end if
            deallocate(h_diag_ptr)
        end if
        stable_c = logical(stable, kind=c_bool)
        error_c = int(error, kind=c_ip)

    end function stability_check_c_wrapper

    subroutine convert_solver_settings_c(settings_c, settings, callbacks, &
                                         stability_callbacks)
        !
        ! this subroutine converts C solver settings to Fortran, stores their optional
        ! callback functions and host contexts in the callback bundles of the solver
        ! and of its internal stability check, which the converted settings hand to the
        ! callback wrappers as opaque context, and prints the warning for settings that
        ! were not initialized, which the conversion replaces by the default settings
        !
        use opentrustregion, only: solver_settings_type, warn_uninitialized_settings, &
                                   settings_uninitialized_warning_msg, &
                                   stability_settings_uninitialized_warning_msg

        type(solver_settings_type_c), intent(in) :: settings_c
        type(solver_settings_type), intent(out) :: settings
        type(c_callbacks_type), intent(inout), target :: callbacks, stability_callbacks

        ! store the optional C function pointers and the host context in the bundle,
        ! which may already hold others
        call store_optional_c_callbacks( &
            callbacks, settings_c%initialized, settings_c%precond, &
            settings_c%precond_pd, settings_c%project, &
            settings_c%get_extra_trial_vectors, settings_c%conv_check, c_null_funptr, &
            settings_c%stability_hess_x, c_null_funptr, settings_c%logger, &
            settings_c%context)

        ! the solver lets the internal stability check inherit its optional callback
        ! functions when the nested settings do not provide their own, so the nested
        ! bundle starts as a copy of the solver's and only the provided ones override
        stability_callbacks = callbacks
        call store_optional_c_callbacks( &
            stability_callbacks, &
            settings_c%initialized .and. settings_c%stability_settings%initialized, &
            settings_c%stability_settings%precond, c_null_funptr, &
            settings_c%stability_settings%project, &
            settings_c%stability_settings%get_extra_trial_vectors, c_null_funptr, &
            settings_c%stability_settings%conv_check, c_null_funptr, &
            settings_c%stability_settings%approx_hess_x, &
            settings_c%stability_settings%logger, settings_c%stability_settings%context)
        callbacks%stability => stability_callbacks

        ! convert settings and hand the bundles to the callback wrappers as opaque
        ! context
        settings = settings_c
        settings%context => callbacks
        settings%stability_settings%context => stability_callbacks

        ! the conversion replaces settings that were not initialized by the default
        ! settings, so the solver cannot detect these anymore and the warning is
        ! printed here regardless of the verbosity, through the logging function of the
        ! converted settings or on stderr since settings that were not initialized
        ! provide none
        if (.not. settings_c%initialized) then
            call warn_uninitialized_settings(settings, &
                                             settings_uninitialized_warning_msg)
        else if (.not. settings_c%stability_settings%initialized) then
            call warn_uninitialized_settings( &
                settings, stability_settings_uninitialized_warning_msg)
        end if

    end subroutine convert_solver_settings_c

    subroutine store_optional_c_callbacks( &
        callbacks, initialized, precond_c_funptr, precond_pd_c_funptr, &
        project_c_funptr, get_extra_trial_vectors_c_funptr, conv_check_c_funptr, &
        conv_check_stability_c_funptr, stability_hess_x_c_funptr, &
        approx_hess_x_c_funptr, logger_c_funptr, context_c)
        !
        ! this subroutine stores the optional C callback functions and the host context
        ! of C settings in a callback bundle, skipping the ones that were not provided
        ! and all of them if the settings were not initialized
        !
        type(c_callbacks_type), intent(inout) :: callbacks
        logical(c_bool), intent(in) :: initialized
        type(c_funptr), intent(in) :: &
            precond_c_funptr, precond_pd_c_funptr, project_c_funptr, &
            get_extra_trial_vectors_c_funptr, conv_check_c_funptr, &
            conv_check_stability_c_funptr, stability_hess_x_c_funptr, &
            approx_hess_x_c_funptr, logger_c_funptr
        type(c_ptr), intent(in) :: context_c

        ! settings that were not initialized are replaced by the default settings and
        ! therefore provide neither optional callback functions nor a host context
        if (.not. initialized) return

        ! store the host context if one was provided
        if (c_associated(context_c)) callbacks%host_context = context_c

        ! associate the C pointers that were provided to Fortran procedure pointers
        if (c_associated(precond_c_funptr)) &
            call c_f_procpointer(cptr=precond_c_funptr, fptr=callbacks%precond)
        if (c_associated(precond_pd_c_funptr)) &
            call c_f_procpointer(cptr=precond_pd_c_funptr, fptr=callbacks%precond_pd)
        if (c_associated(project_c_funptr)) &
            call c_f_procpointer(cptr=project_c_funptr, fptr=callbacks%project)
        if (c_associated(get_extra_trial_vectors_c_funptr)) &
            call c_f_procpointer(cptr=get_extra_trial_vectors_c_funptr, &
                                 fptr=callbacks%get_extra_trial_vectors)
        if (c_associated(conv_check_c_funptr)) &
            call c_f_procpointer(cptr=conv_check_c_funptr, fptr=callbacks%conv_check)
        if (c_associated(conv_check_stability_c_funptr)) call c_f_procpointer( &
            cptr=conv_check_stability_c_funptr, fptr=callbacks%conv_check_stability)
        if (c_associated(stability_hess_x_c_funptr)) call c_f_procpointer( &
            cptr=stability_hess_x_c_funptr, fptr=callbacks%stability_hess_x)
        if (c_associated(approx_hess_x_c_funptr)) call c_f_procpointer( &
            cptr=approx_hess_x_c_funptr, fptr=callbacks%approx_hess_x)
        if (c_associated(logger_c_funptr)) &
            call c_f_procpointer(cptr=logger_c_funptr, fptr=callbacks%logger)

    end subroutine store_optional_c_callbacks

    subroutine copy_state_to_c(callbacks)
        !
        ! this subroutine copies the state of a call into the host's record if the
        ! host gave one, before a callback function into C and when the call returns
        !
        type(c_callbacks_type), intent(in) :: callbacks

        ! copy the solver and the stability check state into the host's records which
        ! were given
        if (associated(callbacks%solver_state_c)) &
            callbacks%solver_state_c = callbacks%solver_state
        if (associated(callbacks%stability_state_c)) &
            callbacks%stability_state_c = callbacks%stability_state

    end subroutine copy_state_to_c

    subroutine update_orbs_f_wrapper(kappa, func, grad, h_diag, hess_x, error, context)
        !
        ! this subroutine exposes a C-implemented orbital update function to Fortran
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp) :: func_c
        real(c_rp), pointer :: kappa_c(:), grad_c(:), h_diag_c(:)
        type(c_funptr) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            kappa_c => kappa
            grad_c => grad
            h_diag_c => h_diag
        else
            allocate(kappa_c(size(kappa)))
            allocate(grad_c(size(kappa)))
            allocate(h_diag_c(size(kappa)))
            kappa_c = real(kappa, kind=c_rp)
        end if

        ! call update_orbs C function and keep the Hessian linear transformation it
        ! returns in the context
        hess_x_c_funptr = c_null_funptr
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            error_c = callbacks%update_orbs(kappa_c, func_c, grad_c, h_diag_c, &
                                            hess_x_c_funptr, callbacks%host_context)
            if (c_associated(hess_x_c_funptr)) then
                call c_f_procpointer(cptr=hess_x_c_funptr, fptr=callbacks%hess_x)
                if (associated(callbacks%stability)) &
                    callbacks%stability%hess_x => callbacks%hess_x
            end if
        class default
            error_c = 1
        end select

        ! convert arguments to Fortran kind
        func = real(func_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            grad = real(grad_c, kind=rp)
            h_diag = real(h_diag_c, kind=rp)
            deallocate(kappa_c)
            deallocate(grad_c)
            deallocate(h_diag_c)
        end if

        ! associate procedure pointer to wrapper function if the C function provided a
        ! Hessian linear transformation, otherwise leave it disassociated so that the
        ! solver reports the missing one
        if (c_associated(hess_x_c_funptr)) then
            hess_x => hess_x_f_wrapper
        else
            hess_x => null()
        end if

    end subroutine update_orbs_f_wrapper

    subroutine hess_x_f_wrapper(x, hess_x, error, context)
        !
        ! this subroutine exposes a C-implemented Hessian linear transformation to
        ! Fortran
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! call C function
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            call hess_x_f_wrapper_impl(callbacks%hess_x, x, hess_x, error, &
                                       callbacks%host_context)
        class default
            error = 1
        end select

    end subroutine hess_x_f_wrapper

    subroutine hess_x_f_wrapper_impl(hess_x_funptr, x, hess_x, error, host_context)
        !
        ! this subroutine wraps the Hessian linear transformation subroutine to convert
        ! Fortran variables to C variables for a given function pointer, which receives
        ! the given host context
        !
        procedure(hess_x_c_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        type(c_ptr), intent(in) :: host_context

        real(c_rp), pointer :: x_c(:), hess_x_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            x_c => x
            hess_x_c => hess_x
        else
            allocate(x_c(size(x)))
            allocate(hess_x_c(size(x)))
            x_c = real(x, kind=c_rp)
        end if

        ! call C function
        error_c = hess_x_funptr(x_c, hess_x_c, host_context)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            hess_x = real(hess_x_c, kind=rp)
            deallocate(x_c)
            deallocate(hess_x_c)
        end if

    end subroutine hess_x_f_wrapper_impl

    function obj_func_f_wrapper(kappa, error, context) result(obj_func)
        !
        ! this function exposes a C-implemented objective function to Fortran
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: obj_func

        real(c_rp) :: obj_func_c
        real(c_rp), pointer :: kappa_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            kappa_c => kappa
        else
            allocate(kappa_c(size(kappa)))
            kappa_c = real(kappa, kind=c_rp)
        end if

        ! call obj_func C function
        obj_func_c = 0.0_c_rp
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            error_c = callbacks%obj_func(kappa_c, obj_func_c, callbacks%host_context)
        class default
            error_c = 1
        end select

        ! convert arguments to Fortran kind
        obj_func = real(obj_func_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            deallocate(kappa_c)
        end if

    end function obj_func_f_wrapper

    subroutine precond_f_wrapper(residual, mu, precond_residual, error, context)
        !
        ! this subroutine exposes a C-implemented preconditioner function to Fortran
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp) :: mu_c
        real(c_rp), pointer :: residual_c(:), precond_residual_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        mu_c = real(mu, kind=c_rp)
        if (rp == c_rp) then
            residual_c => residual
            precond_residual_c => precond_residual
        else
            allocate(residual_c(size(residual)))
            allocate(precond_residual_c(size(residual)))
            residual_c = real(residual, kind=c_rp)
        end if

        ! call precond C function
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            error_c = callbacks%precond(residual_c, mu_c, precond_residual_c, &
                                        callbacks%host_context)
        class default
            error_c = 1
        end select

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            precond_residual = real(precond_residual_c, kind=rp)
            deallocate(residual_c)
            deallocate(precond_residual_c)
        end if

    end subroutine precond_f_wrapper

    subroutine precond_pd_f_wrapper(residual, precond_residual, error, context)
        !
        ! this subroutine exposes a C-implemented positive-definite preconditioner
        ! function to Fortran
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp), pointer :: residual_c(:), precond_residual_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            residual_c => residual
            precond_residual_c => precond_residual
        else
            allocate(residual_c(size(residual)))
            allocate(precond_residual_c(size(residual)))
            residual_c = real(residual, kind=c_rp)
        end if

        ! call precond_pd C function
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            error_c = callbacks%precond_pd(residual_c, precond_residual_c, &
                                           callbacks%host_context)
        class default
            error_c = 1
        end select

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            precond_residual = real(precond_residual_c, kind=rp)
            deallocate(residual_c)
            deallocate(precond_residual_c)
        end if

    end subroutine precond_pd_f_wrapper

    subroutine project_f_wrapper(vector, error, context)
        !
        ! this subroutine exposes a C-implemented projection function to Fortran
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp), pointer :: vector_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            vector_c => vector
        else
            allocate(vector_c(size(vector)))
            vector_c = real(vector, kind=c_rp)
        end if

        ! call project C function
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            error_c = callbacks%project(vector_c, callbacks%host_context)
        class default
            error_c = 1
        end select

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            vector = real(vector_c, kind=rp)
            deallocate(vector_c)
        end if

    end subroutine project_f_wrapper

    function conv_check_f_wrapper(error, context) result(converged)
        !
        ! this function exposes a C-implemented convergence check function to Fortran
        !
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        integer(c_ip) :: error_c
        logical(c_bool) :: converged_c

        ! call conv_check C function
        converged_c = .false._c_bool
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            error_c = callbacks%conv_check(converged_c, callbacks%host_context)
        class default
            error_c = 1
        end select

        ! convert arguments to Fortran kind
        converged = logical(converged_c)
        error = int(error_c, kind=ip)

    end function conv_check_f_wrapper

    subroutine stability_hess_x_f_wrapper(x, hess_x, error, context)
        !
        ! this subroutine exposes a C-implemented stability check Hessian linear
        ! transformation to Fortran
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! call C function
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            call hess_x_f_wrapper_impl(callbacks%stability_hess_x, x, hess_x, error, &
                                       callbacks%host_context)
        class default
            error = 1
        end select

    end subroutine stability_hess_x_f_wrapper

    subroutine approx_hess_x_f_wrapper(x, hess_x, error, context)
        !
        ! this subroutine exposes a C-implemented approximate Hessian linear
        ! transformation to Fortran
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! call C function
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            call hess_x_f_wrapper_impl(callbacks%approx_hess_x, x, hess_x, error, &
                                       callbacks%host_context)
        class default
            error = 1
        end select

    end subroutine approx_hess_x_f_wrapper

    subroutine get_extra_trial_vectors_f_wrapper(trial_vectors, error, context)
        !
        ! this subroutine exposes a C-implemented extra trial vector function to Fortran
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp), pointer :: trial_vectors_c(:, :)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            trial_vectors_c => trial_vectors
        else
            allocate(trial_vectors_c(size(trial_vectors, 1), size(trial_vectors, 2)))
            trial_vectors_c = real(trial_vectors, kind=c_rp)
        end if

        ! call get_extra_trial_vectors C function
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            error_c = callbacks%get_extra_trial_vectors(trial_vectors_c, size( &
                trial_vectors, 2, kind=c_ip), callbacks%host_context)
        class default
            error_c = 1
        end select

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            trial_vectors = real(trial_vectors_c, kind=rp)
            deallocate(trial_vectors_c)
        end if

    end subroutine get_extra_trial_vectors_f_wrapper

    function conv_check_stability_f_wrapper(residual, eigval, error, context) &
        result(converged)
        !
        ! this function exposes a C-implemented convergence check function to Fortran
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: eigval
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        real(c_rp), pointer :: residual_c(:)
        real(c_rp) :: eigval_c
        logical(c_bool) :: converged_c
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            residual_c => residual
        else
            allocate(residual_c(size(residual)))
            residual_c = real(residual, kind=c_rp)
        end if
        eigval_c = real(eigval, kind=c_rp)

        ! call conv_check C function
        converged_c = .false._c_bool
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            error_c = callbacks%conv_check_stability( &
                residual_c, eigval_c, converged_c, callbacks%host_context)
        class default
            error_c = 1
        end select

        ! convert arguments to Fortran kind
        converged = logical(converged_c)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) deallocate(residual_c)

    end function conv_check_stability_f_wrapper

    subroutine logger_f_wrapper(message, context)
        !
        ! this subroutine exposes a C-implemented logger function to Fortran
        !
        character(len=*), intent(in) :: message
        class(*), intent(in), pointer :: context

        character(kind=c_char), allocatable :: message_c(:)
        integer(ip) :: message_len, i

        ! copy to C character
        message_len = len_trim(message)
        allocate(message_c(0:message_len))
        do i = 1, message_len
            message_c(i - 1) = message(i:i)
        end do

        ! append null terminator
        message_c(message_len) = c_null_char

        ! call logging C function
        select type (callbacks => context)
        class is (c_callbacks_type)
            call copy_state_to_c(callbacks)
            call callbacks%logger(message_c, callbacks%host_context)
        end select
        deallocate(message_c)

    end subroutine logger_f_wrapper

    subroutine warn_uninitialized_settings_c(logger_c_funptr, context_c, message)
        !
        ! this subroutine prints the warning that settings were not initialized through
        ! a C logging function with the host context or on stderr without one,
        ! regardless of the verbosity, since the default settings which replace
        ! settings that were not initialized would silence it
        !
        use opentrustregion, only: stderr

        type(c_funptr), intent(in) :: logger_c_funptr
        type(c_ptr), intent(in) :: context_c
        character(len=*), intent(in) :: message

        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        if (c_associated(logger_c_funptr)) then
            call c_f_procpointer(cptr=logger_c_funptr, fptr=callbacks%logger)
            callbacks%host_context = context_c
            context => callbacks
            call logger_f_wrapper(" "//message, context)
        else
            write(stderr, '(A)') " "//message
        end if

    end subroutine warn_uninitialized_settings_c

    subroutine init_solver_settings_c(settings) bind(C, name="init_solver_settings")
        !
        ! this subroutine initializes the C solver settings
        !
        type(solver_settings_type_c), intent(inout) :: settings

        settings = default_solver_settings

    end subroutine init_solver_settings_c

    subroutine init_stability_settings_c(settings) &
        bind(C, name="init_stability_settings")
        !
        ! this subroutine initializes the C stability settings
        !
        type(stability_settings_type_c), intent(inout) :: settings

        settings = default_stability_settings

    end subroutine init_stability_settings_c

    subroutine assign_solver_f_c(settings, settings_c)
        !
        ! this subroutine converts solver settings from C to Fortran
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(out) :: settings
        type(solver_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%precond)) then
                settings%precond => precond_f_wrapper
            else
                settings%precond => null()
            end if
            if (c_associated(settings_c%precond_pd)) then
                settings%precond_pd => precond_pd_f_wrapper
            else
                settings%precond_pd => null()
            end if
            if (c_associated(settings_c%project)) then
                settings%project => project_f_wrapper
            else
                settings%project => null()
            end if
            if (c_associated(settings_c%get_extra_trial_vectors)) then
                settings%get_extra_trial_vectors => get_extra_trial_vectors_f_wrapper
            else
                settings%get_extra_trial_vectors => null()
            end if
            if (c_associated(settings_c%conv_check)) then
                settings%conv_check => conv_check_f_wrapper
            else
                settings%conv_check => null()
            end if
            if (c_associated(settings_c%stability_hess_x)) then
                settings%stability_hess_x => stability_hess_x_f_wrapper
            else
                settings%stability_hess_x => null()
            end if
            if (c_associated(settings_c%logger)) then
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert logicals
            settings%stability = logical(settings_c%stability)
            settings%line_search = logical(settings_c%line_search)
            settings%refresh_hess = logical(settings_c%refresh_hess)
            settings%grad_noise = logical(settings_c%grad_noise)
            settings%hess_symm = logical(settings_c%hess_symm)

            ! convert reals
            settings%conv_tol = real(settings_c%conv_tol, kind=rp)
            settings%start_trust_radius = real(settings_c%start_trust_radius, kind=rp)
            settings%global_red_factor = real(settings_c%global_red_factor, kind=rp)
            settings%local_red_factor = real(settings_c%local_red_factor, kind=rp)

            ! convert integers
            settings%n_random_trial_vectors = &
                int(settings_c%n_random_trial_vectors, kind=ip)
            settings%n_extra_trial_vectors = &
                int(settings_c%n_extra_trial_vectors, kind=ip)
            settings%n_macro = int(settings_c%n_macro, kind=ip)
            settings%n_micro = int(settings_c%n_micro, kind=ip)
            settings%jacobi_davidson_start = &
                int(settings_c%jacobi_davidson_start, kind=ip)
            settings%seed = int(settings_c%seed, kind=ip)
            settings%verbose = int(settings_c%verbose, kind=ip)

            ! convert characters
            settings%subsystem_solver = character_from_c(settings_c%subsystem_solver)
            settings%trust_region_shape = &
                character_from_c(settings_c%trust_region_shape)

            ! convert objects
            settings%stability_settings = settings_c%stability_settings

            ! set settings to initialized
            settings%initialized = .true.
        else
            ! settings that were not initialized are set to the default values here
            ! rather than by the solver, whose initialization would also discard the
            ! callback bundles the C wrapper hands over as context
            settings = default_solver_settings
        end if

    end subroutine assign_solver_f_c

    subroutine assign_stability_f_c(settings, settings_c)
        !
        ! this subroutine converts stability check settings from C to Fortran
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type), intent(out) :: settings
        type(stability_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%precond)) then
                settings%precond => precond_f_wrapper
            else
                settings%precond => null()
            end if
            if (c_associated(settings_c%project)) then
                settings%project => project_f_wrapper
            else
                settings%project => null()
            end if
            if (c_associated(settings_c%approx_hess_x)) then
                settings%approx_hess_x => approx_hess_x_f_wrapper
            else
                settings%approx_hess_x => null()
            end if
            if (c_associated(settings_c%get_extra_trial_vectors)) then
                settings%get_extra_trial_vectors => get_extra_trial_vectors_f_wrapper
            else
                settings%get_extra_trial_vectors => null()
            end if
            if (c_associated(settings_c%conv_check)) then
                settings%conv_check => conv_check_stability_f_wrapper
            else
                settings%conv_check => null()
            end if
            if (c_associated(settings_c%logger)) then
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert logicals
            settings%hess_symm = logical(settings_c%hess_symm)
            settings%stop_on_instability = logical(settings_c%stop_on_instability)

            ! convert reals
            settings%conv_tol = real(settings_c%conv_tol, kind=rp)

            ! convert integers
            settings%n_random_trial_vectors = &
                int(settings_c%n_random_trial_vectors, kind=ip)
            settings%n_extra_trial_vectors = &
                int(settings_c%n_extra_trial_vectors, kind=ip)
            settings%n_iter = int(settings_c%n_iter, kind=ip)
            settings%jacobi_davidson_start = &
                int(settings_c%jacobi_davidson_start, kind=ip)
            settings%seed = int(settings_c%seed, kind=ip)
            settings%verbose = int(settings_c%verbose, kind=ip)

            ! convert characters
            settings%diag_solver = character_from_c(settings_c%diag_solver)

            ! set settings to initialized
            settings%initialized = .true.
        else
            ! settings that were not initialized are set to the default values here
            ! rather than by the stability check, whose initialization would also
            ! discard the callback bundle the C wrapper hands over as context
            settings = default_stability_settings
        end if

    end subroutine assign_stability_f_c

    subroutine assign_solver_c_f(settings_c, settings)
        !
        ! this subroutine converts solver settings from Fortran to C
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type_c), intent(out) :: settings_c
        type(solver_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%precond = c_null_funptr
            settings_c%precond_pd = c_null_funptr
            settings_c%project = c_null_funptr
            settings_c%get_extra_trial_vectors = c_null_funptr
            settings_c%conv_check = c_null_funptr
            settings_c%stability_hess_x = c_null_funptr
            settings_c%logger = c_null_funptr

            ! convert logicals
            settings_c%stability = logical(settings%stability, kind=c_bool)
            settings_c%line_search = logical(settings%line_search, kind=c_bool)
            settings_c%refresh_hess = logical(settings%refresh_hess, kind=c_bool)
            settings_c%grad_noise = logical(settings%grad_noise, kind=c_bool)
            settings_c%hess_symm = logical(settings%hess_symm, kind=c_bool)

            ! convert reals
            settings_c%conv_tol = real(settings%conv_tol, kind=c_rp)
            settings_c%start_trust_radius = real(settings%start_trust_radius, kind=c_rp)
            settings_c%global_red_factor = real(settings%global_red_factor, kind=c_rp)
            settings_c%local_red_factor = real(settings%local_red_factor, kind=c_rp)

            ! convert integers
            settings_c%n_random_trial_vectors = &
                int(settings%n_random_trial_vectors, kind=c_ip)
            settings_c%n_extra_trial_vectors = &
                int(settings%n_extra_trial_vectors, kind=c_ip)
            settings_c%n_macro = int(settings%n_macro, kind=c_ip)
            settings_c%n_micro = int(settings%n_micro, kind=c_ip)
            settings_c%jacobi_davidson_start = &
                int(settings%jacobi_davidson_start, kind=c_ip)
            settings_c%seed = int(settings%seed, kind=c_ip)
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! convert characters
            settings_c%subsystem_solver = character_to_c(settings%subsystem_solver)
            settings_c%trust_region_shape = character_to_c(settings%trust_region_shape)

            ! convert objects
            settings_c%stability_settings = settings%stability_settings

            ! the host context cannot be converted
            settings_c%context = c_null_ptr

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_solver_c_f

    subroutine assign_stability_c_f(settings_c, settings)
        !
        ! this subroutine converts stability settings from Fortran to C
        !
        use opentrustregion, only: stability_settings_type

        type(stability_settings_type_c), intent(out) :: settings_c
        type(stability_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%precond = c_null_funptr
            settings_c%project = c_null_funptr
            settings_c%approx_hess_x = c_null_funptr
            settings_c%get_extra_trial_vectors = c_null_funptr
            settings_c%conv_check = c_null_funptr
            settings_c%logger = c_null_funptr

            ! convert logicals
            settings_c%hess_symm = logical(settings%hess_symm, kind=c_bool)
            settings_c%stop_on_instability = &
                logical(settings%stop_on_instability, kind=c_bool)

            ! convert reals
            settings_c%conv_tol = real(settings%conv_tol, kind=c_rp)

            ! convert integers
            settings_c%n_random_trial_vectors = &
                int(settings%n_random_trial_vectors, kind=c_ip)
            settings_c%n_extra_trial_vectors = &
                int(settings%n_extra_trial_vectors, kind=c_ip)
            settings_c%n_iter = int(settings%n_iter, kind=c_ip)
            settings_c%jacobi_davidson_start = &
                int(settings%jacobi_davidson_start, kind=c_ip)
            settings_c%seed = int(settings%seed, kind=c_ip)
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! convert characters
            settings_c%diag_solver = character_to_c(settings%diag_solver)

            ! the host context cannot be converted
            settings_c%context = c_null_ptr

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_stability_c_f

    subroutine assign_solver_state_c_f(state_c, state)
        !
        ! this subroutine converts a solver state from Fortran to C
        !
        type(solver_state_type_c), intent(out) :: state_c
        type(solver_state_type), intent(in) :: state

        ! convert logicals
        state_c%max_precision_reached = &
            logical(state%max_precision_reached, kind=c_bool)

        ! convert reals
        state_c%func = real(state%func, kind=c_rp)
        state_c%grad_rms = real(state%grad_rms, kind=c_rp)
        state_c%trust_radius = real(state%trust_radius, kind=c_rp)
        state_c%relative_residual = real(state%relative_residual, kind=c_rp)
        state_c%level_shift = real(state%level_shift, kind=c_rp)
        state_c%step_size = real(state%step_size, kind=c_rp)

        ! convert integers
        state_c%macro_iteration = int(state%macro_iteration, kind=c_ip)
        state_c%micro_iteration = int(state%micro_iteration, kind=c_ip)
        state_c%jacobi_davidson_iteration = &
            int(state%jacobi_davidson_iteration, kind=c_ip)
        state_c%n_update_orbs = int(state%n_update_orbs, kind=c_ip)
        state_c%n_obj_func = int(state%n_obj_func, kind=c_ip)
        state_c%n_hess_x = int(state%n_hess_x, kind=c_ip)

        ! convert objects
        state_c%stability_state = state%stability_state

    end subroutine assign_solver_state_c_f

    subroutine assign_stability_state_c_f(state_c, state)
        !
        ! this subroutine converts a stability check state from Fortran to C
        !
        type(stability_state_type_c), intent(out) :: state_c
        type(stability_state_type), intent(in) :: state

        ! convert reals
        state_c%min_eigval = real(state%min_eigval, kind=c_rp)
        state_c%residual_rms = real(state%residual_rms, kind=c_rp)

        ! convert integers
        state_c%iteration = int(state%iteration, kind=c_ip)
        state_c%n_hess_x = int(state%n_hess_x, kind=c_ip)

    end subroutine assign_stability_state_c_f

    function character_from_c(char_c) result(char_f)
        !
        ! this function converts a C null-terminated character array to a Fortran
        ! character
        !
        character(kind=c_char), intent(in) :: char_c(*)
        character(len=:), allocatable :: char_f

        integer(ip) :: n

        ! allocate Fortran character
        n = 1
        do while (char_c(n) /= c_null_char)
            n = n + 1
        end do
        allocate(character(len=n - 1) :: char_f)

        ! copy and convert each character
        char_f = transfer(char_c(1:n - 1), char_f)

    end function character_from_c

    function character_to_c(char_f) result(char_c)
        !
        ! this function converts a Fortran keyword to a C null-terminated character
        ! array of the size of the keyword fields of the C settings, whose remainder
        ! is filled with null characters
        !
        character(len=*), intent(in) :: char_f
        character(kind=c_char) :: char_c(kw_len + 1)

        integer(ip) :: n

        ! fill with null characters, which terminate the copied characters
        char_c = c_null_char

        ! copy and convert each character
        n = len_trim(char_f)
        char_c(1:n) = transfer(char_f(1:n), char_c(1:n))

    end function character_to_c

end module c_interface
