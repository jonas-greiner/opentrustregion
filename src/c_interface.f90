! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface

    use opentrustregion, only: rp, ip, solver_type, stability_check_type, &
                               standard_solver => solver, &
                               standard_stability_check => stability_check
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_int32_t, c_bool, &
                                           c_ptr, c_funptr, c_f_pointer, &
                                           c_f_procpointer, c_associated, c_char, &
                                           c_null_char

    implicit none

    integer, parameter :: c_rp = c_double
#ifdef USE_ILP64
    integer, parameter :: c_ip = c_int64_t  ! 64-bit integers
#else
    integer, parameter :: c_ip = c_int32_t  ! 32-bit integers
#endif

    ! define procedure pointer which will point to the Fortran procedures
    procedure(update_orbs_c_type), pointer :: update_orbs_before_wrapping => null()
    procedure(hess_x_c_type), pointer :: hess_x_before_wrapping => null()
    procedure(obj_func_c_type), pointer :: obj_func_before_wrapping => null()
    procedure(precond_c_type), pointer :: precond_before_wrapping => null()
    procedure(conv_check_c_type), pointer :: conv_check_before_wrapping => null()
    procedure(logger_c_type), pointer :: logger_before_wrapping => null()

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function update_orbs_c_type(kappa_c, func_c, grad_c, h_diag_c, &
                                    hess_x_c_funptr) result(error) bind(C)
            import :: c_rp, c_funptr, c_ip

            real(c_rp), intent(in) :: kappa_c(*)
            real(c_rp), intent(out) :: func_c
            real(c_rp), intent(out) :: grad_c(*), h_diag_c(*)
            type(c_funptr), intent(out) :: hess_x_c_funptr
            integer(c_ip) :: error
        end function update_orbs_c_type
    end interface

    abstract interface
        function hess_x_c_type(x_c, hess_x_c) result(error) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in) :: x_c(*)
            real(c_rp), intent(out) :: hess_x_c(*)
            integer(c_ip) :: error
        end function hess_x_c_type
    end interface

    abstract interface
        function obj_func_c_type(kappa_c, func) result(error) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in) :: kappa_c(*)
            real(c_rp), intent(out) :: func
            integer(c_ip) :: error
        end function obj_func_c_type
    end interface

    abstract interface
        function precond_c_type(residual_c, mu_c, precond_residual_c) result(error) &
            bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in) :: residual_c(*), mu_c
            real(c_rp), intent(out) :: precond_residual_c(*)
            integer(c_ip) :: error
        end function precond_c_type
    end interface

    abstract interface
        function conv_check_c_type(converged) result(error) bind(C)
            import :: c_bool, c_ip

            logical(c_bool), intent(out) :: converged
            integer(c_ip) :: error
        end function conv_check_c_type
    end interface

    abstract interface
        subroutine logger_c_type(message) bind(C)
            import :: c_char

            character(kind=c_char), intent(in) :: message(*)
        end subroutine logger_c_type
    end interface

    ! interface for setting default values
    interface set_default_c_ptr
        module procedure set_default_c_ptr_integer, set_default_c_ptr_real, &
            set_default_c_ptr_logical
    end interface

    ! interfaces for solver and stability C wrapper subroutines
    interface
        function solver_c_wrapper_type(update_orbs_c_funptr, obj_func_c_funptr, &
                                       n_param_c, precond_c_funptr, &
                                       conv_check_c_funptr, stability_c_ptr, &
                                       line_search_c_ptr, davidson_c_ptr, &
                                       jacobi_davidson_c_ptr, &
                                       prefer_jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                       n_random_trial_vectors_c_ptr, &
                                       start_trust_radius_c_ptr, n_macro_c_ptr, &
                                       n_micro_c_ptr, global_red_factor_c_ptr, &
                                       local_red_factor_c_ptr, seed_c_ptr, &
                                       verbose_c_ptr, logger_c_funptr) result(error_c) &
            bind(C, name="solver")

        import :: c_funptr, c_ip, c_ptr      

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             precond_c_funptr, conv_check_c_funptr, &
                                             logger_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(c_ptr), intent(in), value :: stability_c_ptr, line_search_c_ptr, &
                                          davidson_c_ptr, jacobi_davidson_c_ptr, &
                                          prefer_jacobi_davidson_c_ptr, &
                                          conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, &
                                          start_trust_radius_c_ptr, n_macro_c_ptr, &
                                          n_micro_c_ptr, global_red_factor_c_ptr, &
                                          local_red_factor_c_ptr, seed_c_ptr, &
                                          verbose_c_ptr
        integer(c_ip) :: error_c

        end function solver_c_wrapper_type
    end interface

    interface
        function stability_check_c_wrapper_type(h_diag_c_ptr, hess_x_c_funptr, &
                                                n_param_c, stable_c, kappa_c_ptr, &
                                                precond_c_funptr, &
                                                jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                                n_random_trial_vectors_c_ptr, &
                                                n_iter_c_ptr, verbose_c_ptr, &
                                                logger_c_funptr) result(error_c) &
            bind(C, name="stability_check")

            import :: c_ptr, c_ip, c_funptr, c_bool

            type(c_ptr), intent(in), value :: h_diag_c_ptr, kappa_c_ptr, &
                                               jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                               n_random_trial_vectors_c_ptr, &
                                               n_iter_c_ptr, verbose_c_ptr
            integer(c_ip), intent(in), value :: n_param_c
            type(c_funptr), intent(in), value :: hess_x_c_funptr, precond_c_funptr, &
                                                 logger_c_funptr
            logical(c_bool), intent(out) :: stable_c
            integer(c_ip) :: error_c

        end function stability_check_c_wrapper_type
    end interface

    procedure(solver_type), pointer :: solver => standard_solver
    procedure(stability_check_type), pointer :: stability_check => &
        standard_stability_check

    ! create function pointers to ensure that routines comply with interface
    procedure(solver_c_wrapper_type), pointer :: solver_c_wrapper_ptr => &
        solver_c_wrapper
    procedure(stability_check_c_wrapper_type), pointer :: &
        stability_check_c_wrapper_ptr => stability_check_c_wrapper

contains

    function solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, n_param_c, &
                              precond_c_funptr, conv_check_c_funptr, stability_c_ptr, &
                              line_search_c_ptr, davidson_c_ptr, &
                              jacobi_davidson_c_ptr, prefer_jacobi_davidson_c_ptr, &
                              conv_tol_c_ptr, n_random_trial_vectors_c_ptr, &
                              start_trust_radius_c_ptr, n_macro_c_ptr, n_micro_c_ptr, &
                              global_red_factor_c_ptr, local_red_factor_c_ptr, &
                              seed_c_ptr, verbose_c_ptr, logger_c_funptr) &
        result(error_c) bind(C, name="solver")
        !
        ! this subroutine wraps the solver subroutine to convert C variables to Fortran
        ! variables
        !
        use opentrustregion, only: solver_stability_default, &
                                   solver_line_search_default, &
                                   solver_davidson_default, &
                                   solver_jacobi_davidson_default, &
                                   solver_prefer_jacobi_davidson_default, &
                                   solver_conv_tol_default, &
                                   solver_n_random_trial_vectors_default, &
                                   solver_start_trust_radius_default, &
                                   solver_n_macro_default, solver_n_micro_default, &
                                   solver_global_red_factor_default, &
                                   solver_local_red_factor_default, &
                                   solver_seed_default, solver_verbose_default
                                   

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             precond_c_funptr, conv_check_c_funptr, &
                                             logger_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(c_ptr), intent(in), value :: stability_c_ptr, line_search_c_ptr, &
                                          davidson_c_ptr, jacobi_davidson_c_ptr, &
                                          prefer_jacobi_davidson_c_ptr, &
                                          conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, &
                                          start_trust_radius_c_ptr, n_macro_c_ptr, &
                                          n_micro_c_ptr, global_red_factor_c_ptr, &
                                          local_red_factor_c_ptr, seed_c_ptr, &
                                          verbose_c_ptr
        integer(c_ip) :: error_c

        logical :: stability, line_search, davidson, jacobi_davidson, &
                   prefer_jacobi_davidson
        real(rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_param, n_random_trial_vectors, n_macro, n_micro, seed, &
                       verbose, error
        procedure(update_orbs_c_wrapper), pointer :: update_orbs
        procedure(obj_func_c_wrapper), pointer :: obj_func
        procedure(precond_c_wrapper), pointer :: precond
        procedure(conv_check_c_wrapper), pointer :: conv_check
        procedure(logger_c_wrapper), pointer :: logger

        ! set optional arguments to default values
        stability = set_default_c_ptr(stability_c_ptr, solver_stability_default)
        line_search = set_default_c_ptr(line_search_c_ptr, solver_line_search_default)
        davidson = set_default_c_ptr(davidson_c_ptr, solver_davidson_default)
        jacobi_davidson = set_default_c_ptr(jacobi_davidson_c_ptr, &
                                            solver_jacobi_davidson_default)
        prefer_jacobi_davidson = set_default_c_ptr(prefer_jacobi_davidson_c_ptr, &
                                                  solver_prefer_jacobi_davidson_default)
        conv_tol = set_default_c_ptr(conv_tol_c_ptr, solver_conv_tol_default)
        n_random_trial_vectors = set_default_c_ptr(n_random_trial_vectors_c_ptr, &
                                                  solver_n_random_trial_vectors_default)
        start_trust_radius = set_default_c_ptr(start_trust_radius_c_ptr, &
                                               solver_start_trust_radius_default)
        n_macro = set_default_c_ptr(n_macro_c_ptr, solver_n_macro_default)
        n_micro = set_default_c_ptr(n_micro_c_ptr, solver_n_micro_default)
        global_red_factor = set_default_c_ptr(global_red_factor_c_ptr, &
                                              solver_global_red_factor_default)
        local_red_factor = set_default_c_ptr(local_red_factor_c_ptr, &
                                             solver_local_red_factor_default)
        seed = set_default_c_ptr(seed_c_ptr, solver_seed_default)
        verbose = set_default_c_ptr(verbose_c_ptr, solver_verbose_default)

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=update_orbs_c_funptr, &
                             fptr=update_orbs_before_wrapping)
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_before_wrapping)

        ! associate procedure pointer to wrapper function
        update_orbs => update_orbs_c_wrapper
        obj_func => obj_func_c_wrapper

        ! convert dummy argument to Fortran kind
        n_param = int(n_param_c, kind=ip)

        ! call solver with or without optional procedures
        if (c_associated(precond_c_funptr)) then
            call c_f_procpointer(cptr=precond_c_funptr, fptr=precond_before_wrapping)
            precond => precond_c_wrapper
        else
            precond => null()
        end if
        if (c_associated(conv_check_c_funptr)) then
            call c_f_procpointer(cptr=conv_check_c_funptr, &
                                 fptr=conv_check_before_wrapping)
            conv_check => conv_check_c_wrapper
        else
            conv_check => null()
        end if
        if (c_associated(logger_c_funptr)) then
            call c_f_procpointer(cptr=logger_c_funptr, fptr=logger_before_wrapping)
            logger => logger_c_wrapper
        else
            logger => null()
        end if
        call solver(update_orbs, obj_func, n_param, error, precond, conv_check, &
                    stability, line_search, davidson, jacobi_davidson, &
                    prefer_jacobi_davidson, conv_tol, n_random_trial_vectors, &
                    start_trust_radius, n_macro, n_micro, global_red_factor, &
                    local_red_factor, seed, verbose, logger)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function solver_c_wrapper

    function stability_check_c_wrapper(h_diag_c_ptr, hess_x_c_funptr, n_param_c, &
                                       stable_c, kappa_c_ptr, precond_c_funptr, &
                                       jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                       n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                       verbose_c_ptr, logger_c_funptr) result(error_c) &
        bind(C, name="stability_check")
        !
        ! this subroutine wraps the stability check subroutine to convert C variables
        ! to Fortran variables
        !
        use opentrustregion, only: stability_jacobi_davidson_default, &
                                   stability_conv_tol_default, &
                                   stability_n_random_trial_vectors_default, &
                                   stability_n_iter_default, &
                                   stability_verbose_default

        type(c_ptr), intent(in), value :: h_diag_c_ptr, kappa_c_ptr, &
                                          jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                          verbose_c_ptr
        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: hess_x_c_funptr, precond_c_funptr, &
                                             logger_c_funptr
        logical(c_bool), intent(out) :: stable_c
        integer(c_ip) :: error_c

        real(rp) :: conv_tol
        real(rp), pointer :: h_diag_ptr(:), kappa_ptr(:)
        real(c_rp), pointer :: h_diag_ptr_c(:), kappa_ptr_c(:)
        logical :: stable, jacobi_davidson
        integer(ip) :: n_random_trial_vectors, n_iter, verbose, error
        procedure(hess_x_c_wrapper), pointer :: hess_x
        procedure(precond_c_wrapper), pointer :: precond
        procedure(logger_c_wrapper), pointer :: logger

        ! set optional arguments to default values
        jacobi_davidson = set_default_c_ptr(jacobi_davidson_c_ptr, &
                                            stability_jacobi_davidson_default)
        conv_tol = set_default_c_ptr(conv_tol_c_ptr, stability_conv_tol_default)
        n_random_trial_vectors = set_default_c_ptr(n_random_trial_vectors_c_ptr, &
                                               stability_n_random_trial_vectors_default)
        n_iter = set_default_c_ptr(n_iter_c_ptr, stability_n_iter_default)
        verbose = set_default_c_ptr(verbose_c_ptr, stability_verbose_default)

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_before_wrapping)

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_c_wrapper

        ! convert C pointers
        if (rp == c_rp) then
            call c_f_pointer(cptr=h_diag_c_ptr, fptr=h_diag_ptr, shape=[n_param_c])
        else
            call c_f_pointer(cptr=h_diag_c_ptr, fptr=h_diag_ptr_c, shape=[n_param_c])
            allocate(h_diag_ptr(n_param_c))
            h_diag_ptr = real(h_diag_ptr_c, kind=rp)
        end if
        if (c_associated(kappa_c_ptr)) then
            if (rp == c_rp) then
                call c_f_pointer(cptr=kappa_c_ptr, fptr=kappa_ptr, shape=[n_param_c])
            else
                call c_f_pointer(cptr=kappa_c_ptr, fptr=kappa_ptr_c, shape=[n_param_c])
                allocate(kappa_ptr(n_param_c))
            end if
        end if
        if (c_associated(precond_c_funptr)) then
            call c_f_procpointer(cptr=precond_c_funptr, fptr=precond_before_wrapping)
            precond => precond_c_wrapper
        else
            precond => null()
        end if
        if (c_associated(logger_c_funptr)) then
            call c_f_procpointer(cptr=logger_c_funptr, fptr=logger_before_wrapping)
            logger => logger_c_wrapper
        else
            logger => null()
        end if

        ! call stability check
        if (c_associated(kappa_c_ptr)) then
            call stability_check(h_diag_ptr, hess_x, stable, error, kappa_ptr, &
                                 precond, jacobi_davidson, conv_tol, &
                                 n_random_trial_vectors, n_iter, verbose, logger)
            if (rp /= c_rp) then
                call c_f_pointer(cptr=kappa_c_ptr, fptr=kappa_ptr_c, shape=[n_param_c])
                kappa_ptr_c = kappa_ptr
                deallocate(kappa_ptr)
            end if
        else
            call stability_check(h_diag_ptr, hess_x, stable, error, precond=precond, &
                                 jacobi_davidson=jacobi_davidson, conv_tol=conv_tol, &
                                 n_random_trial_vectors=n_random_trial_vectors, &
                                 n_iter=n_iter, verbose=verbose, logger=logger)
        end if
        if (rp /= c_rp) deallocate(h_diag_ptr)

        ! convert return arguments to C kind
        stable_c = logical(stable, kind=c_bool)
        error_c = int(error, kind=c_ip)

    end function stability_check_c_wrapper

    subroutine update_orbs_c_wrapper(kappa, func, grad, h_diag, hess_x, error)
        !
        ! this subroutine wraps the orbital update subroutine to convert C variables to
        ! Fortran variables
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x
        integer(ip), intent(out) :: error

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

        ! call update_orbs C function
        error_c = update_orbs_before_wrapping(kappa_c, func_c, grad_c, h_diag_c, &
                                              hess_x_c_funptr)

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

        ! associate the input C pointer to hess_x function to a Fortran procedure
        ! pointer
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_before_wrapping)

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_c_wrapper

    end subroutine update_orbs_c_wrapper

    subroutine hess_x_c_wrapper(x, hess_x, error)
        !
        ! this subroutine wraps the Hessian linear transformation to convert C variables
        ! to Fortran variables
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

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
        error_c = hess_x_before_wrapping(x_c, hess_x_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            hess_x = real(hess_x_c, kind=rp)
            deallocate(x_c)
            deallocate(hess_x_c)
        end if

    end subroutine hess_x_c_wrapper

    function obj_func_c_wrapper(kappa, error) result(obj_func)
        !
        ! this function wraps the objective function to convert C variables to Fortran
        ! variables
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
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
        error_c = obj_func_before_wrapping(kappa_c, obj_func_c)

        ! convert arguments to Fortran kind
        obj_func = real(obj_func_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            deallocate(kappa_c)
        end if

    end function obj_func_c_wrapper

    subroutine precond_c_wrapper(residual, mu, precond_residual, error)
        !
        ! this subroutine wraps the preconditioner function to convert C variables to 
        ! Fortran variables
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

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
        error_c = precond_before_wrapping(residual_c, mu_c, precond_residual_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            precond_residual = real(precond_residual_c, kind=rp)
            deallocate(residual_c)
            deallocate(precond_residual_c)
        end if

    end subroutine precond_c_wrapper

    function conv_check_c_wrapper(error) result(converged)
        !
        ! this function wraps the convergence check function to convert C variables to 
        ! Fortran variables
        !
        integer(ip), intent(out) :: error
        logical :: converged

        integer(c_ip) :: error_c
        logical(c_bool) :: converged_c

        ! call conv_check C function
        error_c = conv_check_before_wrapping(converged_c)

        ! convert arguments to Fortran kind
        converged = logical(converged_c)
        error = int(error_c, kind=ip)

    end function conv_check_c_wrapper

    subroutine logger_c_wrapper(message)
        !
        ! this function wraps the logging subroutine to convert C variables to Fortran 
        ! variables
        !
        character(*), intent(in) :: message

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
        call logger_before_wrapping(message_c)
        deallocate(message_c)

    end subroutine logger_c_wrapper

    function set_default_c_ptr_integer(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for reals
        !
        integer(ip), intent(in)  :: default_value
        type(c_ptr), intent(in) :: optional_value

        integer(ip) :: variable

        integer(c_ip), pointer :: variable_ptr

        if (c_associated(optional_value)) then
            call c_f_pointer(cptr=optional_value, fptr=variable_ptr)
            variable = int(variable_ptr, kind=ip)
        else
            variable = default_value
        end if

    end function set_default_c_ptr_integer

    function set_default_c_ptr_real(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for reals
        !
        real(rp), intent(in)  :: default_value
        type(c_ptr), intent(in) :: optional_value

        real(rp) :: variable

        real(c_rp), pointer :: variable_ptr

        if (c_associated(optional_value)) then
            call c_f_pointer(cptr=optional_value, fptr=variable_ptr)
            variable = real(variable_ptr, kind=rp)
        else
            variable = default_value
        end if

    end function set_default_c_ptr_real

    function set_default_c_ptr_logical(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for logicals
        !
        logical, intent(in)  :: default_value
        type(c_ptr), intent(in) :: optional_value

        logical :: variable

        logical(c_bool), pointer :: variable_ptr

        if (c_associated(optional_value)) then
            call c_f_pointer(cptr=optional_value, fptr=variable_ptr)
            variable = logical(variable_ptr)
        else
            variable = default_value
        end if

    end function set_default_c_ptr_logical

end module c_interface
