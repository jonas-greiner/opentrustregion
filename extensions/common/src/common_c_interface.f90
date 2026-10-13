! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface

    use opentrustregion, only: ip, rp, hess_x_type
    use c_interface, only: c_ip, c_rp, c_callbacks_type, update_orbs_c_type, &
                           obj_func_c_type, hess_x_c_type, precond_c_type, &
                           precond_pd_c_type, project_c_type, &
                           get_extra_trial_vectors_c_type, conv_check_c_type, &
                           conv_check_stability_c_type, logger_c_type
    use otr_common, only: extension_type, extension_level_type, update_orbs_callback, &
                          hess_x_callback
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
    ! level of the object and the bundle of the host's C callback functions and context
    ! of the level
    type :: holder_level_type
        type(extension_level_type), pointer :: level => null()
        type(c_callbacks_type) :: callbacks
    end type

    ! holder of an extension object for C hosts, with the levels of the solver and of
    ! its stability check
    type :: extension_holder_type
        type(holder_level_type) :: solver_level, stability_level
    end type

    ! global variables
    integer(ip) :: n_param

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
    procedure(conv_check_stability_c_type), pointer :: &
        conv_check_stability_c_wrapper_ptr => conv_check_stability_c_wrapper
    procedure(logger_c_type), pointer :: logger_c_wrapper_ptr => logger_c_wrapper

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

    subroutine setup_extension_holder(holder, update_orbs_c_funptr, obj_func_c_funptr, &
                                      solver_settings_c, solver_settings, settings, &
                                      error)
        !
        ! this subroutine sets up the holder of an extension object for a C factory by
        ! bundling the host's orbital update and objective function with the callback
        ! functions and host contexts of the host's C solver settings, which it
        ! converts to the Fortran solver settings the Fortran factory wires, after
        ! checking that the holder was not already wired into them
        !
        use opentrustregion, only: solver_settings_type, settings_type, verbosity_error
        use c_interface, only: solver_settings_type_c, convert_solver_settings_c
        use otr_common, only: settings_wired_twice_error_msg

        type(extension_holder_type), intent(inout), target :: holder
        type(c_funptr), intent(in) :: update_orbs_c_funptr, obj_func_c_funptr
        type(solver_settings_type_c), intent(in) :: solver_settings_c
        type(solver_settings_type), intent(out) :: solver_settings
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        logical :: nested_own

        ! initialize error flag
        error = 0

        ! settings already wired by the holder would let the callback functions wired
        ! with its levels forward to themselves, which is detected before the holder is
        ! set up anew
        if (c_associated(solver_settings_c%context, c_loc(holder%solver_level)) .or. &
            c_associated(solver_settings_c%stability_settings%context, &
                         c_loc(holder%stability_level))) then
            call settings%log(settings_wired_twice_error_msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! bundle the host's callback functions and convert the solver settings, whose
        ! optional callback functions and host contexts join the bundles of the holder
        holder = extension_holder_type()
        call c_f_procpointer(cptr=update_orbs_c_funptr, &
                             fptr=holder%solver_level%callbacks%update_orbs)
        call c_f_procpointer(cptr=obj_func_c_funptr, &
                             fptr=holder%solver_level%callbacks%obj_func)
        call convert_solver_settings_c(solver_settings_c, solver_settings, &
                                       holder%solver_level%callbacks, &
                                       holder%stability_level%callbacks)

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
        if (.not. nested_own) solver_settings%stability_settings%context => null()

    end subroutine setup_extension_holder

    subroutine wire_solver_settings_c(holder, extension, solver_settings, &
                                      solver_settings_c)
        !
        ! this subroutine writes the solver settings the Fortran factory of an
        ! extension object wired back into the host's C solver settings, with the C
        ! wrappers of the role callback functions in their place and the addresses of
        ! the levels of the holder as the contexts of the levels the factory wired,
        ! which also initializes C settings that were not initialized, as the
        ! conversion and the factory did
        !
        use opentrustregion, only: solver_settings_type
        use c_interface, only: solver_settings_type_c, assignment(=)

        type(extension_holder_type), intent(inout), target :: holder
        class(extension_type), intent(in), target :: extension
        type(solver_settings_type), intent(in) :: solver_settings
        type(solver_settings_type_c), intent(inout) :: solver_settings_c

        ! hold the levels of the object, which the C wrappers find at the addresses of
        ! the levels of the holder
        holder%solver_level%level => extension%solver_level
        holder%stability_level%level => extension%stability_level

        ! convert the settings, which leaves the callback functions and contexts unset
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
            solver_settings_c%stability_hess_x = c_funloc(stability_hess_x_c_wrapper)
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
                nested_c%conv_check = c_funloc(conv_check_stability_c_wrapper)
            if (associated(nested%logger)) nested_c%logger = c_funloc(logger_c_wrapper)
        end associate

    end subroutine wire_solver_settings_c

    function update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, hess_x_c_funptr, &
                                   context_c) result(error_c) bind(C)
        !
        ! this function exposes the orbital update to C hosts, calling it with the
        ! level of an extension object held at the address it receives as context, and
        ! returns the C wrapper of the Hessian linear transformation the orbital update
        ! returns
        !
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
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

        ! convert arguments to C kind and return the C wrapper of the Hessian linear
        ! transformation if the orbital update returned one
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            grad_c(:n_param) = real(grad, kind=c_rp)
            h_diag_c(:n_param) = real(h_diag, kind=c_rp)
            deallocate(kappa, grad, h_diag)
        end if
        if (associated(hess_x_funptr)) hess_x_c_funptr = c_funloc(hess_x_c_wrapper)

    end function update_orbs_c_wrapper

    function obj_func_c_wrapper(kappa_c, func_c, context_c) result(error_c) bind(C)
        !
        ! this function forwards the objective function for C hosts straight to the
        ! host's C function in the bundle of the level of a holder at the address it
        ! receives as context, with the host context of the level
        !
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(obj_func_c_type), pointer :: host_obj_func

        ! find the level of the holder, without which there is nothing to call
        func_c = 0.0_c_rp
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_obj_func => held%callbacks%obj_func
        error_c = host_obj_func(kappa_c, func_c, held%callbacks%host_context)

    end function obj_func_c_wrapper

    function hess_x_c_wrapper(x_c, hess_x_c, context_c) result(error_c) bind(C)
        !
        ! this function exposes the Hessian linear transformation to C hosts, calling
        ! it with the level of an extension object held at the address it receives as
        ! context
        !
        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
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

    function precond_c_wrapper(residual_c, mu_c, precond_residual_c, context_c) &
        result(error_c) bind(C)
        !
        ! this function forwards the level-shifted preconditioner for C hosts straight
        ! to the host's C function in the bundle of the level of a holder at the
        ! address it receives as context, with the host context of the level
        !
        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(in) :: mu_c
        real(c_rp), intent(out), target :: precond_residual_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(precond_c_type), pointer :: host_precond

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_precond => held%callbacks%precond
        error_c = host_precond(residual_c, mu_c, precond_residual_c, &
                               held%callbacks%host_context)

    end function precond_c_wrapper

    function precond_pd_c_wrapper(residual_c, precond_residual_c, context_c) &
        result(error_c) bind(C)
        !
        ! this function forwards the positive-definite preconditioner for C hosts
        ! straight to the host's C function in the bundle of the level of a holder at
        ! the address it receives as context, with the host context of the level
        !
        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(out), target :: precond_residual_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(precond_pd_c_type), pointer :: host_precond_pd

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_precond_pd => held%callbacks%precond_pd
        error_c = &
            host_precond_pd(residual_c, precond_residual_c, held%callbacks%host_context)

    end function precond_pd_c_wrapper

    function project_c_wrapper(vector_c, context_c) result(error_c) bind(C)
        !
        ! this function forwards the projection for C hosts straight to the host's C
        ! function in the bundle of the level of a holder at the address it receives as
        ! context, with the host context of the level
        !
        real(c_rp), intent(inout), target :: vector_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(project_c_type), pointer :: host_project

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_project => held%callbacks%project
        error_c = host_project(vector_c, held%callbacks%host_context)

    end function project_c_wrapper

    function get_extra_trial_vectors_c_wrapper( &
        trial_vectors_c, n_extra_trial_vectors_c, context_c) result(error_c) bind(C)
        !
        ! this function forwards the extra trial vector function for C hosts straight
        ! to the host's C function in the bundle of the level of a holder at the
        ! address it receives as context, with the host context of the level
        !
        real(c_rp), intent(out), target :: trial_vectors_c(*)
        integer(c_ip), intent(in), value :: n_extra_trial_vectors_c
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(get_extra_trial_vectors_c_type), pointer :: &
            host_get_extra_trial_vectors

        ! find the level of the holder, without which there is nothing to call
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_get_extra_trial_vectors => held%callbacks%get_extra_trial_vectors
        error_c = host_get_extra_trial_vectors( &
            trial_vectors_c, n_extra_trial_vectors_c, held%callbacks%host_context)

    end function get_extra_trial_vectors_c_wrapper

    function conv_check_c_wrapper(converged_c, context_c) result(error_c) bind(C)
        !
        ! this function forwards the convergence check of the solver for C hosts
        ! straight to the host's C function in the bundle of the level of a holder at
        ! the address it receives as context, with the host context of the level
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

    function stability_hess_x_c_wrapper(x_c, hess_x_c, context_c) result(error_c) &
        bind(C)
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

    function approx_hess_x_c_wrapper(x_c, hess_x_c, context_c) result(error_c) bind(C)
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

    function conv_check_stability_c_wrapper(residual_c, eigval_c, converged_c, &
                                            context_c) result(error_c) bind(C)
        !
        ! this function forwards the convergence check of the stability check for C
        ! hosts straight to the host's C function in the bundle of the level of a
        ! holder at the address it receives as context, with the host context of the
        ! level
        !
        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(in) :: eigval_c
        logical(c_bool), intent(out) :: converged_c
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        type(holder_level_type), pointer :: held
        procedure(conv_check_stability_c_type), pointer :: host_conv_check_stability

        ! find the level of the holder, without which there is nothing to call
        converged_c = .false._c_bool
        held => holder_level(context_c)
        if (.not. associated(held)) then
            error_c = 1
            return
        end if

        ! call the host's C function with the host context of the level, through a
        ! local procedure pointer, since gfortran rejects assumed-size arrays as
        ! arguments of a procedure pointer component
        host_conv_check_stability => held%callbacks%conv_check_stability
        error_c = host_conv_check_stability(residual_c, eigval_c, converged_c, &
                                            held%callbacks%host_context)

    end function conv_check_stability_c_wrapper

    subroutine logger_c_wrapper(message_c, context_c) bind(C)
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

    function update_orbs_c_wrapper_impl( &
        update_orbs_before_wrapping, hess_x_before_wrapping_funptr, &
        hess_x_c_wrapper_funptr, kappa_c, func_c, grad_c, h_diag_c, hess_x_c_funptr, &
        context_c, update_orbs_orig_c) result(error_c)
        !
        ! this function wraps the orbital update subroutine to convert Fortran
        ! variables to C variables, handing the host context to the Fortran subroutine
        ! in a callback bundle, together with the host's original orbital update for
        ! an orbital update that wraps it
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use c_interface, only: hess_x_c_type, update_orbs_c_type

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_before_wrapping
        procedure(hess_x_type), intent(inout), pointer :: hess_x_before_wrapping_funptr
        procedure(hess_x_c_type) :: hess_x_c_wrapper_funptr
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in) :: context_c
        procedure(update_orbs_c_type), intent(in), pointer, optional :: &
            update_orbs_orig_c
        integer(c_ip) :: error_c

        real(rp) :: func
        real(rp), pointer :: kappa(:), grad(:), h_diag(:)
        integer(ip) :: error
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        ! hand the host context and the original orbital update over as context
        callbacks%host_context = context_c
        if (present(update_orbs_orig_c)) callbacks%update_orbs => update_orbs_orig_c
        context => callbacks

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            kappa => kappa_c(:n_param)
            grad => grad_c(:n_param)
            h_diag => h_diag_c(:n_param)
        else
            allocate(kappa(n_param))
            allocate(grad(n_param))
            allocate(h_diag(n_param))
            kappa = real(kappa_c(:n_param), kind=rp)
        end if

        ! call update_orbs Fortran function
        call update_orbs_before_wrapping(kappa, func, grad, h_diag, &
                                         hess_x_before_wrapping_funptr, error, context)

        ! convert arguments to Fortran kind
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            grad_c(:n_param) = real(grad, kind=c_rp)
            h_diag_c(:n_param) = real(h_diag, kind=c_rp)
            deallocate(kappa)
            deallocate(grad)
            deallocate(h_diag)
        end if

        ! get a C function pointer to the hess_x wrapper function
        hess_x_c_funptr = c_funloc(hess_x_c_wrapper_funptr)

    end function update_orbs_c_wrapper_impl

    function hess_x_c_wrapper_impl(hess_x_funptr, x_c, hess_x_c, context_c) &
        result(error_c)
        !
        ! this function wraps the Hessian linear transformation to convert Fortran
        ! variables to C variables, handing the host context to the Fortran subroutine
        ! in a callback bundle
        !
        use opentrustregion, only: hess_x_type

        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        type(c_ptr), intent(in) :: context_c
        integer(c_ip) :: error_c

        real(rp), pointer :: x(:), hess_x(:)
        integer(ip) :: error
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        ! hand the host context over as context
        callbacks%host_context = context_c
        context => callbacks

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            x => x_c(:n_param)
            hess_x => hess_x_c(:n_param)
        else
            allocate(x(n_param))
            allocate(hess_x(n_param))
            x = real(x_c(:n_param), kind=rp)
        end if

        ! call Fortran function
        call hess_x_funptr(x, hess_x, error, context)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            hess_x_c(:n_param) = real(hess_x, kind=c_rp)
            deallocate(x)
            deallocate(hess_x)
        end if

    end function hess_x_c_wrapper_impl

    function obj_func_c_wrapper_impl(obj_func_before_wrapping, kappa_c, func_c, &
                                     context_c) result(error_c)
        !
        ! this function wraps the objective function subroutine to convert Fortran
        ! variables to C variables, handing the host context to the Fortran function
        ! in a callback bundle
        !
        use opentrustregion, only: obj_func_type

        procedure(obj_func_type), intent(in), pointer :: obj_func_before_wrapping
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        type(c_ptr), intent(in) :: context_c
        integer(c_ip) :: error_c

        real(rp) :: func
        real(rp), pointer :: kappa(:)
        integer(ip) :: error
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        ! hand the host context over as context
        callbacks%host_context = context_c
        context => callbacks

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            kappa => kappa_c(:n_param)
        else
            allocate(kappa(n_param))
            kappa = real(kappa_c(:n_param), kind=rp)
        end if

        ! call obj_func Fortran function
        func = obj_func_before_wrapping(kappa, error, context)

        ! convert arguments to Fortran kind
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            deallocate(kappa)
        end if

    end function obj_func_c_wrapper_impl

    function precond_c_wrapper_impl(precond_before_wrapping, residual_c, mu_c, &
                                    precond_residual_c, context_c) result(error_c)
        !
        ! this function wraps the level-shifted preconditioner subroutine to convert
        ! Fortran variables to C variables, handing the host context to the Fortran
        ! subroutine in a callback bundle
        !
        use opentrustregion, only: precond_type

        procedure(precond_type), intent(in), pointer :: precond_before_wrapping
        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(in) :: mu_c
        real(c_rp), intent(out), target :: precond_residual_c(*)
        type(c_ptr), intent(in) :: context_c
        integer(c_ip) :: error_c

        real(rp) :: mu
        real(rp), pointer :: residual(:), precond_residual(:)
        integer(ip) :: error
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        ! hand the host context over as context
        callbacks%host_context = context_c
        context => callbacks

        ! convert arguments to Fortran kind
        mu = real(mu_c, kind=rp)
        if (rp == c_rp) then
            residual => residual_c(:n_param)
            precond_residual => precond_residual_c(:n_param)
        else
            allocate(residual(n_param), precond_residual(n_param))
            residual = real(residual_c(:n_param), kind=rp)
        end if

        ! call preconditioner Fortran subroutine
        call precond_before_wrapping(residual, mu, precond_residual, error, context)

        ! convert arguments to Fortran kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            precond_residual_c(:n_param) = real(precond_residual, kind=c_rp)
            deallocate(residual, precond_residual)
        end if

    end function precond_c_wrapper_impl

    function precond_pd_c_wrapper_impl(precond_pd_before_wrapping, residual_c, &
                                       precond_residual_c, context_c) result(error_c)
        !
        ! this function wraps the positive-definite preconditioner subroutine to
        ! convert Fortran variables to C variables, handing the host context to the
        ! Fortran subroutine in a callback bundle
        !
        use opentrustregion, only: precond_pd_type

        procedure(precond_pd_type), intent(in), pointer :: precond_pd_before_wrapping
        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(out), target :: precond_residual_c(*)
        type(c_ptr), intent(in) :: context_c
        integer(c_ip) :: error_c

        real(rp), pointer :: residual(:), precond_residual(:)
        integer(ip) :: error
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        ! hand the host context over as context
        callbacks%host_context = context_c
        context => callbacks

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            residual => residual_c(:n_param)
            precond_residual => precond_residual_c(:n_param)
        else
            allocate(residual(n_param), precond_residual(n_param))
            residual = real(residual_c(:n_param), kind=rp)
        end if

        ! call preconditioner Fortran subroutine
        call precond_pd_before_wrapping(residual, precond_residual, error, context)

        ! convert arguments to Fortran kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            precond_residual_c(:n_param) = real(precond_residual, kind=c_rp)
            deallocate(residual, precond_residual)
        end if

    end function precond_pd_c_wrapper_impl

    function get_extra_trial_vectors_c_wrapper_impl( &
        get_extra_trial_vectors_before_wrapping, trial_vectors_c, &
        n_extra_trial_vectors_c, context_c) result(error_c)
        !
        ! this function wraps the extra trial vector subroutine to convert Fortran
        ! variables to C variables, handing the host context to the Fortran subroutine
        ! in a callback bundle
        !
        use opentrustregion, only: get_extra_trial_vectors_type

        procedure(get_extra_trial_vectors_type), intent(in), pointer :: &
            get_extra_trial_vectors_before_wrapping
        real(c_rp), intent(out), target :: trial_vectors_c(*)
        integer(c_ip), intent(in) :: n_extra_trial_vectors_c
        type(c_ptr), intent(in) :: context_c
        integer(c_ip) :: error_c

        real(rp), pointer :: trial_vectors(:, :)
        integer(ip) :: n_extra_trial_vectors, error
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        ! hand the host context over as context
        callbacks%host_context = context_c
        context => callbacks

        ! convert arguments to Fortran kind
        n_extra_trial_vectors = int(n_extra_trial_vectors_c, kind=ip)
        if (rp == c_rp) then
            call c_f_pointer(c_loc(trial_vectors_c(1)), trial_vectors, &
                             [n_param, n_extra_trial_vectors])
        else
            allocate(trial_vectors(n_param, n_extra_trial_vectors))
        end if

        ! call extra trial vector Fortran subroutine
        call get_extra_trial_vectors_before_wrapping(trial_vectors, error, context)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            trial_vectors_c(:n_param * n_extra_trial_vectors) = real( &
                reshape(trial_vectors, [n_param * n_extra_trial_vectors]), kind=c_rp)
            deallocate(trial_vectors)
        end if

    end function get_extra_trial_vectors_c_wrapper_impl

    subroutine evaluate_dm_f_wrapper_impl(callbacks, dm, energy, get_response, error, &
                                          fock)
        !
        ! this subroutine wraps the density matrix evaluating C function of a callback
        ! bundle to convert Fortran variables to C variables, handing the host context
        ! of the bundle to it, and keeps the response C function it returns in the
        ! bundle if one is requested
        !
        type(orbital_c_callbacks_type), intent(inout) :: callbacks
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        logical, intent(in) :: get_response
        integer(ip), intent(out) :: error
        real(rp), intent(out), optional, target :: fock(:, :, :)

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :, :), fock_c(:, :, :)
        type(c_funptr) :: get_response_c_funptr
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        nullify(fock_c)
        if (rp == c_rp) then
            dm_c => dm
            if (present(fock)) fock_c => fock
        else
            allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            dm_c = real(dm, kind=c_rp)
            if (present(fock)) allocate(fock_c(size(dm, 1), size(dm, 2), size(dm, 3)))
        end if

        ! call density matrix evaluating C function and keep the response C function it
        ! returns in the bundle
        if (get_response) then
            error_c = callbacks%evaluate_dm( &
                dm_c, energy_c, fock_c, get_response_c_funptr, callbacks%host_context)
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

    end subroutine evaluate_dm_f_wrapper_impl

    subroutine get_response_f_wrapper_impl(callbacks, dm, response, error)
        !
        ! this subroutine wraps the response C function of a callback bundle to convert
        ! Fortran variables to C variables, handing the host context of the bundle to
        ! it
        !
        type(orbital_c_callbacks_type), intent(in) :: callbacks
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out), target :: response(:, :, :)
        integer(ip), intent(out) :: error

        real(c_rp), pointer :: dm_c(:, :, :), response_c(:, :, :)
        integer(c_ip) :: error_c

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

    end subroutine get_response_f_wrapper_impl

end module otr_common_c_interface
