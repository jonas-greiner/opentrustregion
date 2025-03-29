module c_interface

    use opentrustregion, only: rp, ip, standard_solver => solver, &
                               standard_stability_check => stability_check
    use, intrinsic :: iso_c_binding, only: c_long, c_double, c_bool, c_ptr, c_funptr, &
                                              c_f_pointer, c_f_procpointer, c_associated

    implicit none

    ! define procedure pointer which will point to the Fortran procedures
    procedure(update_orbs_c_type), pointer :: update_orbs_before_wrapping => null()
    procedure(hess_x_c_type), pointer :: hess_x_before_wrapping => null()
    procedure(obj_func_c_type), pointer :: obj_func_before_wrapping => null()

    ! C-interoperable interfaces for the callback functions
    abstract interface
        subroutine update_orbs_c_type(kappa, func, grad, h_diag, hess_x) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr, c_funptr

            real(c_double), intent(in) :: kappa(:)

            real(c_double), intent(out) :: func
            type(c_ptr), intent(out) :: grad, h_diag
            type(c_funptr), intent(out) :: hess_x
        end subroutine update_orbs_c_type
    end interface

    abstract interface
        subroutine hess_x_c_type(x, hess_x) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr

            real(c_double), intent(in) :: x(:)

            type(c_ptr), intent(out) :: hess_x
        end subroutine hess_x_c_type
    end interface

    abstract interface
        function obj_func_c_type(kappa) result(func) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double

            real(c_double), intent(in) :: kappa(:)

            real(c_double) :: func
        end function obj_func_c_type
    end interface

    ! interface for setting default values
    interface set_default_c_ptr
        module procedure set_default_c_ptr_integer, set_default_c_ptr_real, &
            set_default_c_ptr_logical
    end interface

    ! interfaces for solver and stability subroutines, different routines can be
    ! injected, for example for testing
    interface
        subroutine solver_type(update_orbs, obj_func, n_param, stability, line_search, &
                               conv_tol, n_random_trial_vectors, start_trust_radius, &
                               n_macro, n_micro, global_red_factor, local_red_factor, &
                               verbose, seed)

            use opentrustregion, only: rp, ip, update_orbs_type, obj_func_type

            procedure(update_orbs_type), intent(in), pointer :: update_orbs
            procedure(obj_func_type), intent(in), pointer :: obj_func
            integer(ip), intent(in) :: n_param
            logical, intent(in), optional :: stability, line_search
            real(rp), intent(in), optional :: conv_tol, start_trust_radius, &
                                              global_red_factor, local_red_factor
            integer(ip), intent(in), optional :: n_random_trial_vectors, n_macro, &
                                                 n_micro, verbose, seed

        end subroutine solver_type
    end interface

    interface
        subroutine stability_check_type(grad, h_diag, hess_x, stable, kappa, conv_tol, &
                                        n_random_trial_vectors, n_iter, verbose)

            use opentrustregion, only: rp, ip, hess_x_type

            real(rp), intent(in) :: grad(:), h_diag(:)
            procedure(hess_x_type), pointer, intent(in) :: hess_x
            logical, intent(out) :: stable
            real(rp), intent(out) :: kappa(:)
            real(rp), intent(in), optional :: conv_tol
            integer(ip), intent(in), optional :: n_random_trial_vectors, n_iter, verbose

        end subroutine stability_check_type
    end interface

    procedure(solver_type), pointer :: solver => standard_solver
    procedure(stability_check_type), pointer :: stability_check => &
        standard_stability_check

contains

    subroutine solver_c_wrapper(update_orbs_c_ptr, obj_func_c_ptr, n_param_c, &
                                stability_c_ptr, line_search_c_ptr, conv_tol_c_ptr, &
                                n_random_trial_vectors_c_ptr, &
                                start_trust_radius_c_ptr, n_macro_c_ptr, &
                                n_micro_c_ptr, global_red_factor_c_ptr, &
                                local_red_factor_c_ptr, verbose_c_ptr, seed_c_ptr) &
        bind(C, name="solver")
        !
        ! this subroutine wraps the solver subroutine to convert C variables to Fortran
        ! variables
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, &
                                   solver_stability_default, &
                                   solver_line_search_default, &
                                   solver_conv_tol_default, &
                                   solver_n_random_trial_vectors_default, &
                                   solver_start_trust_radius_default, &
                                   solver_n_macro_default, solver_n_micro_default, &
                                   solver_global_red_factor_default, &
                                   solver_local_red_factor_default, &
                                   solver_verbose_default, solver_seed_default

        type(c_funptr), intent(in) :: update_orbs_c_ptr, obj_func_c_ptr
        integer(c_long), value, intent(in) :: n_param_c
        type(c_ptr), value, intent(in) :: stability_c_ptr, line_search_c_ptr, &
                                          conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, &
                                          start_trust_radius_c_ptr, n_macro_c_ptr, &
                                          n_micro_c_ptr, global_red_factor_c_ptr, &
                                          local_red_factor_c_ptr, verbose_c_ptr, &
                                          seed_c_ptr

        logical :: stability, line_search
        real(rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_param, n_random_trial_vectors, n_macro, n_micro, verbose, seed
        procedure(update_orbs_type), pointer :: update_orbs
        procedure(obj_func_type), pointer :: obj_func

        ! set optional arguments to default values
        stability = set_default_c_ptr(stability_c_ptr, solver_stability_default)
        line_search = set_default_c_ptr(line_search_c_ptr, solver_line_search_default)
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
        verbose = set_default_c_ptr(verbose_c_ptr, solver_verbose_default)
        seed = set_default_c_ptr(seed_c_ptr, solver_seed_default)

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=update_orbs_c_ptr, fptr=update_orbs_before_wrapping)
        call c_f_procpointer(cptr=obj_func_c_ptr, fptr=obj_func_before_wrapping)

        ! convert dummy argument to Fortran kind
        n_param = int(n_param_c, kind=ip)

        ! associate procedure pointer to wrapper function
        update_orbs => update_orbs_c_wrapper
        obj_func => obj_func_c_wrapper

        ! update orbitals
        call solver(update_orbs, obj_func, n_param, stability, line_search, conv_tol, &
                    n_random_trial_vectors, start_trust_radius, n_macro, n_micro, &
                    global_red_factor, local_red_factor, verbose, seed)

    end subroutine solver_c_wrapper

    subroutine stability_check_c_wrapper(grad_c, h_diag_c, hess_x_c_ptr, n_param_c, &
                                         stable_c, kappa_c, conv_tol_c_ptr, &
                                         n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                         verbose_c_ptr) bind(C, name="stability_check")
        !
        ! this subroutine wraps the stability check subroutine to convert C variables
        ! to Fortran variables
        !
        use opentrustregion, only: hess_x_type, stability_conv_tol_default, &
                                   stability_n_random_trial_vectors_default, &
                                   stability_n_iter_default, stability_verbose_default

        integer(c_long), value, intent(in) :: n_param_c
        real(c_double), intent(in), dimension(n_param_c) :: grad_c, h_diag_c
        type(c_funptr), intent(in) :: hess_x_c_ptr
        logical(c_bool), intent(out) :: stable_c
        real(c_double), intent(out) :: kappa_c(n_param_c)
        type(c_ptr), value, intent(in) :: conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                          verbose_c_ptr

        real(rp) :: conv_tol, grad(n_param_c), h_diag(n_param_c)
        real(rp) :: kappa(n_param_c)
        logical :: stable
        integer(ip) :: n_param, n_random_trial_vectors, n_iter, verbose
        procedure(hess_x_type), pointer :: hess_x

        ! set optional arguments to default values
        conv_tol = set_default_c_ptr(conv_tol_c_ptr, stability_conv_tol_default)
        n_random_trial_vectors = set_default_c_ptr(n_random_trial_vectors_c_ptr, &
                                               stability_n_random_trial_vectors_default)
        n_iter = set_default_c_ptr(n_iter_c_ptr, stability_n_iter_default)
        verbose = set_default_c_ptr(verbose_c_ptr, stability_verbose_default)

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=hess_x_c_ptr, fptr=hess_x_before_wrapping)

        ! convert dummy argument to Fortran kind
        n_param = int(n_param_c, kind=ip)

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_c_wrapper

        ! convert dummy argument to Fortran kind
        grad = grad_c
        h_diag = h_diag_c

        ! update orbitals
        call stability_check(grad, h_diag, hess_x, stable, kappa, conv_tol, &
                             n_random_trial_vectors, n_iter, verbose)

        ! convert return arguments to C kind
        stable_c = stable
        kappa_c = kappa

    end subroutine stability_check_c_wrapper

    subroutine update_orbs_c_wrapper(kappa, func, grad, h_diag, hess_x)
        !
        ! this subroutine wraps the orbital update subroutine to convert C variables to
        ! Fortran variables
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in) :: kappa(:)

        real(rp), intent(out) :: func, grad(:), h_diag(:)
        procedure(hess_x_type), pointer, intent(out) :: hess_x

        real(c_double) :: func_c
        real(c_double) :: kappa_c(size(kappa))
        type(c_ptr) :: grad_c_ptr, h_diag_c_ptr
        type(c_funptr) :: hess_x_c_ptr
        real(c_double), pointer :: grad_ptr(:), h_diag_ptr(:)

        ! convert dummy argument to C kind
        kappa_c = kappa

        ! call update_orbs C function
        call update_orbs_before_wrapping(kappa_c, func_c, grad_c_ptr, h_diag_c_ptr, &
                                         hess_x_c_ptr)

        ! associate C pointer to arrays with Fortran pointer
        call c_f_pointer(cptr=grad_c_ptr, fptr=grad_ptr, shape=[size(kappa)])
        call c_f_pointer(cptr=h_diag_c_ptr, fptr=h_diag_ptr, shape=[size(kappa)])

        ! associate the input C pointer to hess_x function to a Fortran procedure
        ! pointer
        call c_f_procpointer(cptr=hess_x_c_ptr, fptr=hess_x_before_wrapping)

        ! associate procedure pointer to wrapper function
        hess_x => hess_x_c_wrapper

        ! convert dummy argument to Fortran kind
        func = func_c
        grad = grad_ptr
        h_diag = h_diag_ptr

    end subroutine update_orbs_c_wrapper

    function hess_x_c_wrapper(x) result(hess_x)
        !
        ! this function wraps the Hessian linear transformation to convert C variables
        ! to Fortran variables
        !
        real(rp), intent(in) :: x(:)

        real(rp) :: hess_x(size(x))

        real(c_double) :: x_c(size(x))
        type(c_ptr) :: hess_x_c_ptr
        real(c_double), pointer :: hess_x_ptr(:)

        ! convert trial vector to C kind
        x_c = x

        ! call C function
        call hess_x_before_wrapping(x_c, hess_x_c_ptr)

        ! associate C pointer to arrays with Fortran pointer
        call c_f_pointer(cptr=hess_x_c_ptr, fptr=hess_x_ptr, shape=[size(x)])

        ! convert linear transformation to Fortran kind
        hess_x = hess_x_ptr

    end function hess_x_c_wrapper

    function obj_func_c_wrapper(kappa) result(obj_func)
        !
        ! this function wraps the objective function to convert C variables to Fortran
        ! variables
        !
        real(rp), intent(in) :: kappa(:)

        real(rp) :: obj_func

        real(c_double) :: func_c
        real(c_double) :: kappa_c(size(kappa))

        ! convert dummy argument to C kind
        kappa_c = kappa

        ! call update_orbs C function
        func_c = obj_func_before_wrapping(kappa_c)

        ! convert dummy argument to Fortran kind
        obj_func = func_c

    end function obj_func_c_wrapper

    function set_default_c_ptr_integer(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for reals
        !
        integer(ip), intent(in)  :: default_value
        type(c_ptr), intent(in) :: optional_value

        integer(ip) :: variable

        integer(c_long), pointer :: variable_ptr

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

        real(c_double), pointer :: variable_ptr

        if (c_associated(optional_value)) then
            call c_f_pointer(cptr=optional_value, fptr=variable_ptr)
            variable = variable_ptr
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

        logical, pointer :: variable_ptr

        if (c_associated(optional_value)) then
            call c_f_pointer(cptr=optional_value, fptr=variable_ptr)
            variable = variable_ptr
        else
            variable = default_value
        end if

    end function set_default_c_ptr_logical

end module c_interface
