module c_interface_mock

    use trustorbopt, only: stderr
    use c_interface, only: update_orbs_c_type, hess_x_c_type, obj_func_c_type
    use, intrinsic :: iso_c_binding, only: c_long, c_double, c_bool, c_ptr, c_funptr, &
                                              c_f_pointer, c_f_procpointer, c_associated

    implicit none

    logical :: solver_default = .true., stability_check_default = .true.
    logical(c_bool) :: test_solver_interface = .true., &
                       test_stability_check_interface = .true.
    real(c_double), parameter :: tol = 1.d-10

contains

    logical(c_bool) function test_solver_result() bind(C)
        !
        ! this function extracts whether the solver python interface test has passed
        !
        test_solver_result = test_solver_interface

    end function test_solver_result

    logical(c_bool) function test_stability_check_result() bind(C)
        !
        ! this function extracts whether the stability check python interface test has
        ! passed
        !
        test_stability_check_result = test_stability_check_interface

    end function test_stability_check_result

    subroutine mock_solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                     n_param_c, stability_c_ptr, line_search_c_ptr, &
                                     conv_tol_c_ptr, n_random_trial_vectors_c_ptr, &
                                     start_trust_radius_c_ptr, n_macro_c_ptr, &
                                     n_micro_c_ptr, global_red_factor_c_ptr, &
                                     local_red_factor_c_ptr, verbose_c_ptr, &
                                     seed_c_ptr) bind(C, name="mock_solver")
        !
        ! this subroutine is a mock routine for the solver C wrapper subroutine
        !
        type(c_funptr), intent(in) :: update_orbs_c_funptr, obj_func_c_funptr
        integer(c_long), value, intent(in) :: n_param_c
        type(c_ptr), value, intent(in) :: stability_c_ptr, line_search_c_ptr, &
                                          conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, &
                                          start_trust_radius_c_ptr, n_macro_c_ptr, &
                                          n_micro_c_ptr, global_red_factor_c_ptr, &
                                          local_red_factor_c_ptr, verbose_c_ptr, &
                                          seed_c_ptr

        type(c_funptr) :: hess_x_c_funptr
        procedure(update_orbs_c_type), pointer :: update_orbs_funptr
        procedure(hess_x_c_type), pointer :: hess_x_funptr
        procedure(obj_func_c_type), pointer :: obj_func_funptr
        real(c_double), pointer :: grad_ptr(:), h_diag_ptr(:), hess_x_ptr(:)
        real(c_double) :: func
        real(c_double), dimension(n_param_c) :: kappa, x
        type(c_ptr) :: grad_c_ptr, h_diag_c_ptr, hess_x_c_ptr

        logical, pointer :: stability_ptr, line_search_ptr
        real(c_double), pointer :: conv_tol_ptr, start_trust_radius_ptr, &
                                   global_red_factor_ptr, local_red_factor_ptr
        integer(c_long), pointer :: n_random_trial_vectors_ptr, n_macro_ptr, &
                                    n_micro_ptr, verbose_ptr, seed_ptr

        ! get Fortran pointer to passed orbital update routine and call it
        call c_f_procpointer(cptr=update_orbs_c_funptr, fptr=update_orbs_funptr)
        kappa = 1.0_c_double
        call update_orbs_funptr(kappa, func, grad_c_ptr, h_diag_c_ptr, hess_x_c_funptr)

        ! get Fortran pointer to Hessian linear transformation function and call it
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_funptr)
        x = 1.0_c_double
        call hess_x_funptr(x, hess_x_c_ptr)

        ! get Fortran pointers to generated arrays and check whether these are filled
        ! correctly
        call c_f_pointer(cptr=grad_c_ptr, fptr=grad_ptr, shape=[n_param_c])
        call c_f_pointer(cptr=h_diag_c_ptr, fptr=h_diag_ptr, shape=[n_param_c])
        call c_f_pointer(cptr=hess_x_c_ptr, fptr=hess_x_ptr, shape=[n_param_c])
        if (abs(func - 3.0_c_double) > tol) then
            write (stderr, *) "test_solver_py_interface failed: Returned function "// &
                "value of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(grad_ptr - 2.0_c_double) > tol)) then
            write (stderr, *) "test_solver_py_interface failed: Returned gradient "// &
                "of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(h_diag_ptr - 3.0_c_double) > tol)) then
            write (stderr, *) "test_solver_py_interface failed: Returned Hessian "// &
                "diagonal of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(hess_x_ptr - 4.0_c_double) > tol)) then
            write (stderr, *) "test_solver_py_interface failed: Returned Hessian "// &
                "linear transformation from Hessian linear transformation function "// &
                "of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if

        ! get Fortran pointer to passed objective function and call it
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_funptr)
        func = obj_func_funptr(kappa)

        ! check if generated function value is correct
        if (abs(func - 3.0_c_double) > tol) then
            write (stderr, *) "test_solver_py_interface failed: Returned function "// &
                "value of passed objective function wrong."
            test_solver_interface = .false.
        end if

        ! check if check if default arguments are correctly unassociated with values
        if (solver_default) then
            if (c_associated(stability_c_ptr) .or. c_associated(line_search_c_ptr) &
                .or. c_associated(conv_tol_c_ptr) .or. &
                c_associated(n_random_trial_vectors_c_ptr) .or. &
                c_associated(start_trust_radius_c_ptr) .or. &
                c_associated(n_macro_c_ptr) .or. c_associated(n_micro_c_ptr) .or. &
                c_associated(global_red_factor_c_ptr) .or. &
                c_associated(local_red_factor_c_ptr) .or. c_associated(verbose_c_ptr) &
                .or. c_associated(seed_c_ptr)) then
                write (stderr, *) "test_solver_py_interface failed: Passed "// &
                    "optional arguments associated with values."
                test_solver_interface = .false.
            end if
            ! check if check if optional arguments are associated with correct values
        else
            if (.not. (c_associated(stability_c_ptr) .and. &
                       c_associated(line_search_c_ptr) .and. &
                       c_associated(conv_tol_c_ptr) .and. &
                       c_associated(n_random_trial_vectors_c_ptr) .and. &
                       c_associated(start_trust_radius_c_ptr) .and. &
                       c_associated(n_macro_c_ptr) .and. &
                       c_associated(n_micro_c_ptr) .and. &
                       c_associated(global_red_factor_c_ptr) .and. &
                       c_associated(local_red_factor_c_ptr) .and. &
                       c_associated(verbose_c_ptr) .and. &
                       c_associated(seed_c_ptr))) then
                test_solver_interface = .false.
                write (stderr, *) "test_solver_py_interface failed: Passed "// &
                    "optional arguments not associated with values."
            end if
            call c_f_pointer(cptr=stability_c_ptr, fptr=stability_ptr)
            call c_f_pointer(cptr=line_search_c_ptr, fptr=line_search_ptr)
            call c_f_pointer(cptr=conv_tol_c_ptr, fptr=conv_tol_ptr)
            call c_f_pointer(cptr=n_random_trial_vectors_c_ptr, &
                             fptr=n_random_trial_vectors_ptr)
            call c_f_pointer(cptr=start_trust_radius_c_ptr, fptr=start_trust_radius_ptr)
            call c_f_pointer(cptr=n_macro_c_ptr, fptr=n_macro_ptr)
            call c_f_pointer(cptr=n_micro_c_ptr, fptr=n_micro_ptr)
            call c_f_pointer(cptr=global_red_factor_c_ptr, fptr=global_red_factor_ptr)
            call c_f_pointer(cptr=local_red_factor_c_ptr, fptr=local_red_factor_ptr)
            call c_f_pointer(cptr=verbose_c_ptr, fptr=verbose_ptr)
            call c_f_pointer(cptr=seed_c_ptr, fptr=seed_ptr)
            if (stability_ptr .or. .not. line_search_ptr .or. &
                abs(conv_tol_ptr - 1.e-3_c_double) > tol .or. &
                n_random_trial_vectors_ptr /= 5_c_long .or. &
                abs(start_trust_radius_ptr - 0.2_c_double) > tol .or. &
                n_macro_ptr /= 300_c_long .or. n_micro_ptr /= 200_c_long .or. &
                abs(global_red_factor_ptr - 1.e-2_c_double) > tol .or. &
                abs(local_red_factor_ptr - 1.e-3_c_double) > tol .or. &
                verbose_ptr /= 3_c_long .or. seed_ptr /= 33_c_long) then
                test_solver_interface = .false.
                write (stderr, *) "test_solver_py_interface failed: Passed "// &
                    "optional arguments associated with wrong values."
            end if
        end if

        ! move on to optional argument test for next test
        solver_default = .false.

    end subroutine mock_solver_c_wrapper

    subroutine mock_stability_check_c_wrapper(grad_c, h_diag_c, hess_x_c_funptr, &
                                              n_param_c, stable_c, kappa_c, &
                                              conv_tol_c_ptr, &
                                              n_random_trial_vectors_c_ptr, &
                                              n_iter_c_ptr, verbose_c_ptr) &
        bind(C, name="mock_stability_check")
        !
        ! this subroutine is a mock routine for the stability check C wrapper
        ! subroutine
        !
        integer(c_long), value, intent(in) :: n_param_c
        real(c_double), intent(in), dimension(n_param_c) :: grad_c, h_diag_c
        type(c_funptr), intent(in) :: hess_x_c_funptr
        logical(c_bool), intent(out) :: stable_c
        real(c_double), intent(out) :: kappa_c(n_param_c)
        type(c_ptr), value, intent(in) :: conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                          verbose_c_ptr

        procedure(hess_x_c_type), pointer :: hess_x_funptr
        type(c_ptr) :: hess_x_c_ptr
        real(c_double), pointer :: hess_x_ptr(:)
        real(c_double) :: x(n_param_c)
        real(c_double), pointer :: conv_tol_ptr
        integer(c_long), pointer :: n_random_trial_vectors_ptr, n_iter_ptr, verbose_ptr

        ! check if gradient is passed correctly
        if (any(abs(grad_c - 2.0_c_double) > tol)) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "gradient wrong."
            test_stability_check_interface = .false.
        end if

        ! check if Hessian diagonal is passed correctly
        if (any(abs(h_diag_c - 3.0_c_double) > tol)) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "Hessian diagonal wrong."
            test_stability_check_interface = .false.
        end if

        ! get Fortran pointer to Hessian linear transformation function and call it
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_funptr)
        x = 1.0_c_double
        call hess_x_funptr(x, hess_x_c_ptr)

        ! get Fortran pointer to generated array and check whether it is filled
        ! correctly
        call c_f_pointer(cptr=hess_x_c_ptr, fptr=hess_x_ptr, shape=[n_param_c])
        if (any(abs(hess_x_ptr - 4.0_c_double) > tol)) then
            write (stderr, *) "test_stability_check_py_interface failed: Returned "// &
                "Hessian linear transformation from passed Hessian linear "// &
                "transformation function wrong."
            test_stability_check_interface = .false.
        end if

        ! check if check if default arguments are correctly unassociated with values
        if (stability_check_default) then
            if (c_associated(conv_tol_c_ptr) .or. &
                c_associated(n_random_trial_vectors_c_ptr) .or. &
                c_associated(n_iter_c_ptr) .or. c_associated(verbose_c_ptr)) then
                write (stderr, *) "test_stability_check_py_interface failed: "// &
                    "Passed optional arguments associated with values."
                test_stability_check_interface = .false.
            end if
            ! check if check if optional arguments are associated with correct values
        else
            if (.not. (c_associated(conv_tol_c_ptr) .and. &
                       c_associated(n_random_trial_vectors_c_ptr) .and. &
                       c_associated(n_iter_c_ptr) .and. &
                       c_associated(verbose_c_ptr))) then
                test_stability_check_interface = .false.
                write (stderr, *) "test_stability_check_py_interface failed: "// &
                    "Passed optional arguments not associated with values."
            end if
            call c_f_pointer(cptr=conv_tol_c_ptr, fptr=conv_tol_ptr)
            call c_f_pointer(cptr=n_random_trial_vectors_c_ptr, &
                             fptr=n_random_trial_vectors_ptr)
            call c_f_pointer(cptr=n_iter_c_ptr, fptr=n_iter_ptr)
            call c_f_pointer(cptr=verbose_c_ptr, fptr=verbose_ptr)
            if (abs(conv_tol_ptr - 1.e-3_c_double) > tol .or. &
                n_random_trial_vectors_ptr /= 3_c_long .or. n_iter_ptr /= 50_c_long &
                .or. verbose_ptr /= 3_c_long) then
                test_solver_interface = .false.
                write (stderr, *) verbose_ptr /= 3_c_long
                write (stderr, *) "test_stability_check_py_interface failed: "// &
                    "Passed optional arguments associated with wrong values."
            end if
        end if

        ! set return arguments
        stable_c = .false.
        kappa_c = 1.0_c_double

        ! move on to optional argument test for next test
        stability_check_default = .false.

    end subroutine mock_stability_check_c_wrapper

end module c_interface_mock
