! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip, solver_c_wrapper_type, &
                           stability_check_c_wrapper_type, update_orbs_c_type, &
                           hess_x_c_type, obj_func_c_type, precond_c_type, &
                           conv_check_c_type, logger_c_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_f_pointer, &
                                           c_f_procpointer, c_associated, c_null_char

    implicit none

    logical :: solver_default = .true., stability_check_default = .true.
    logical(c_bool) :: test_solver_interface = .true., &
                       test_stability_check_interface = .true.
    real(c_rp), parameter :: tol = 1e-10_c_rp

    ! create function pointers to ensure that routines comply with interface
    procedure(solver_c_wrapper_type), pointer :: solver_c_wrapper_ptr => &
        mock_solver_c_wrapper
    procedure(stability_check_c_wrapper_type), pointer :: &
        stability_check_c_wrapper_ptr => mock_stability_check_c_wrapper

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

    function mock_solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                   n_param_c, precond_c_funptr, conv_check_c_funptr, &
                                   stability_c_ptr, line_search_c_ptr, davidson_c_ptr, &
                                   jacobi_davidson_c_ptr, &
                                   prefer_jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                   n_random_trial_vectors_c_ptr, &
                                   start_trust_radius_c_ptr, n_macro_c_ptr, &
                                   n_micro_c_ptr, global_red_factor_c_ptr, &
                                   local_red_factor_c_ptr, seed_c_ptr, verbose_c_ptr, &
                                   logger_c_funptr) result(error_c) &
        bind(C, name="mock_solver")
        !
        ! this subroutine is a mock routine for the solver C wrapper subroutine
        !
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

        type(c_funptr) :: hess_x_c_funptr
        procedure(update_orbs_c_type), pointer :: update_orbs_funptr
        procedure(hess_x_c_type), pointer :: hess_x_funptr
        procedure(obj_func_c_type), pointer :: obj_func_funptr
        real(c_rp) :: func
        real(c_rp), allocatable :: kappa(:), grad(:), h_diag(:), x(:), hess_x(:), &
                                   residual(:), precond_residual(:)
        procedure(precond_c_type), pointer :: precond_funptr
        procedure(conv_check_c_type), pointer :: conv_check_funptr
        logical, pointer :: stability_ptr, line_search_ptr, davidson_ptr, &
                            jacobi_davidson_ptr, prefer_jacobi_davidson_ptr
        real(c_rp), pointer :: conv_tol_ptr, start_trust_radius_ptr, &
                               global_red_factor_ptr, local_red_factor_ptr
        integer(c_ip), pointer :: n_random_trial_vectors_ptr, n_macro_ptr, &
                                  n_micro_ptr, seed_ptr, verbose_ptr
        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message
        logical(c_bool) :: converged
        integer(c_ip) :: error

        ! get Fortran pointer to passed orbital update routine and call it
        call c_f_procpointer(cptr=update_orbs_c_funptr, fptr=update_orbs_funptr)
        allocate(kappa(n_param_c), grad(n_param_c), h_diag(n_param_c))
        kappa = 1.0_c_rp
        error = update_orbs_funptr(kappa, func, grad, h_diag, hess_x_c_funptr)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_solver_py_interface failed: Passed orbital "// &
                "updating function returned error."
            test_solver_interface = .false.
        end if

        ! get Fortran pointer to Hessian linear transformation function and call it
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_funptr)
        allocate(x(n_param_c), hess_x(n_param_c))
        x = 1.0_c_rp
        error = hess_x_funptr(x, hess_x)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_solver_py_interface failed: Passed Hessian "// &
                "linear transformation function returned error."
            test_solver_interface = .false.
        end if

        ! get Fortran pointers to generated arrays and check whether these are filled
        ! correctly
        if (abs(func - 3.0_c_rp) > tol) then
            write (stderr, *) "test_solver_py_interface failed: Returned function "// &
                "value of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(grad - 2.0_c_rp) > tol)) then
            write (stderr, *) "test_solver_py_interface failed: Returned gradient "// &
                "of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(h_diag - 3.0_c_rp) > tol)) then
            write (stderr, *) "test_solver_py_interface failed: Returned Hessian "// &
                "diagonal of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        if (any(abs(hess_x - 4.0_c_rp) > tol)) then
            write (stderr, *) "test_solver_py_interface failed: Returned Hessian "// &
                "linear transformation from Hessian linear transformation function "// &
                "of passed orbital updating function wrong."
            test_solver_interface = .false.
        end if
        deallocate(grad, h_diag, x, hess_x)

        ! get Fortran pointer to passed objective function and call it
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_funptr)
        error = obj_func_funptr(kappa, func)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_solver_py_interface failed: Passed objective "// &
                "function returned error."
            test_solver_interface = .false.
        end if

        ! check if generated function value is correct
        if (abs(func - 3.0_c_rp) > tol) then
            write (stderr, *) "test_solver_py_interface failed: Returned function "// &
                "value of passed objective function wrong."
            test_solver_interface = .false.
        end if
        deallocate(kappa)

        ! check if check if default arguments are correctly unassociated with values
        if (solver_default) then
            if (c_associated(precond_c_funptr) .or. c_associated(conv_check_c_funptr) &
                .or. c_associated(stability_c_ptr) .or. &
                c_associated(line_search_c_ptr) .or. c_associated(davidson_c_ptr) .or. &
                c_associated(jacobi_davidson_c_ptr) .or. &
                c_associated(prefer_jacobi_davidson_c_ptr) .or. &
                c_associated(conv_tol_c_ptr) .or. &
                c_associated(n_random_trial_vectors_c_ptr) .or. &
                c_associated(start_trust_radius_c_ptr) .or. &
                c_associated(n_macro_c_ptr) .or. c_associated(n_micro_c_ptr) .or. &
                c_associated(global_red_factor_c_ptr) .or. &
                c_associated(local_red_factor_c_ptr) .or. c_associated(seed_c_ptr) &
                .or. c_associated(verbose_c_ptr) .or. c_associated(logger_c_funptr)) &
                then
                write (stderr, *) "test_solver_py_interface failed: Passed "// &
                    "optional arguments associated with values."
                test_solver_interface = .false.
            end if
            ! check if check if optional arguments are associated with correct values
        else
            if (.not. (c_associated(precond_c_funptr) .and. &
                       c_associated(conv_check_c_funptr) .and. &
                       c_associated(stability_c_ptr) .and. &
                       c_associated(line_search_c_ptr) .and. &
                       c_associated(davidson_c_ptr) .and. &
                       c_associated(jacobi_davidson_c_ptr) .and. &
                       c_associated(prefer_jacobi_davidson_c_ptr) .and. &
                       c_associated(conv_tol_c_ptr) .and. &
                       c_associated(n_random_trial_vectors_c_ptr) .and. &
                       c_associated(start_trust_radius_c_ptr) .and. &
                       c_associated(n_macro_c_ptr) .and. &
                       c_associated(n_micro_c_ptr) .and. &
                       c_associated(global_red_factor_c_ptr) .and. &
                       c_associated(local_red_factor_c_ptr) .and. &
                       c_associated(seed_c_ptr) .and. c_associated(verbose_c_ptr) &
                       .and. c_associated(logger_c_funptr))) then
                write (stderr, *) "test_solver_py_interface failed: Passed "// &
                    "optional arguments not associated with values."
                test_solver_interface = .false.
            end if

            ! get Fortran pointer to passed preconditioner function, call it and check
            ! result
            call c_f_procpointer(cptr=precond_c_funptr, fptr=precond_funptr)
            allocate(residual(n_param_c), precond_residual(n_param_c))
            residual = 1.0_c_rp
            error = precond_funptr(residual, 5.0_c_rp, precond_residual)
            if (error /= 0) then
                write (stderr, *) "test_solver_py_interface failed: Passed "// &
                    "preconditioner function returned error."
                test_solver_interface = .false.
            end if
            if (any(abs(precond_residual - 5.0_c_rp) > tol)) then
                write (stderr, *) "test_solver_py_interface failed: Returned "// &
                    "preconditioned residual from preconditioner function wrong."
                test_solver_interface = .false.
            end if
            deallocate(residual, precond_residual)

            ! get Fortran pointer to passed convergence check function, call it and 
            ! check result
            call c_f_procpointer(cptr=conv_check_c_funptr, fptr=conv_check_funptr)
            error = conv_check_funptr(converged)
            if (error /= 0) then
                write (stderr, *) "test_solver_py_interface failed: Passed "// &
                    "convergence check function returned error."
                test_solver_interface = .false.
            end if
            if (.not. converged) then
                write (stderr, *) "test_solver_py_interface failed: Returned "// &
                    "convergence logical of convergence check function wrong."
                test_solver_interface = .false.
            end if

            ! get Fortran pointers to optional arguments and check against reference 
            ! values
            call c_f_pointer(cptr=stability_c_ptr, fptr=stability_ptr)
            call c_f_pointer(cptr=line_search_c_ptr, fptr=line_search_ptr)
            call c_f_pointer(cptr=davidson_c_ptr, fptr=davidson_ptr)
            call c_f_pointer(cptr=jacobi_davidson_c_ptr, fptr=jacobi_davidson_ptr)
            call c_f_pointer(cptr=prefer_jacobi_davidson_c_ptr, &
                             fptr=prefer_jacobi_davidson_ptr)
            call c_f_pointer(cptr=conv_tol_c_ptr, fptr=conv_tol_ptr)
            call c_f_pointer(cptr=n_random_trial_vectors_c_ptr, &
                             fptr=n_random_trial_vectors_ptr)
            call c_f_pointer(cptr=start_trust_radius_c_ptr, fptr=start_trust_radius_ptr)
            call c_f_pointer(cptr=n_macro_c_ptr, fptr=n_macro_ptr)
            call c_f_pointer(cptr=n_micro_c_ptr, fptr=n_micro_ptr)
            call c_f_pointer(cptr=global_red_factor_c_ptr, fptr=global_red_factor_ptr)
            call c_f_pointer(cptr=local_red_factor_c_ptr, fptr=local_red_factor_ptr)
            call c_f_pointer(cptr=seed_c_ptr, fptr=seed_ptr)
            call c_f_pointer(cptr=verbose_c_ptr, fptr=verbose_ptr)
            if (stability_ptr .or. .not. line_search_ptr .or. davidson_ptr .or. &
                jacobi_davidson_ptr .or. .not. prefer_jacobi_davidson_ptr .or. &
                abs(conv_tol_ptr - 1e-3_c_rp) > tol .or. &
                n_random_trial_vectors_ptr /= 5_c_ip .or. &
                abs(start_trust_radius_ptr - 0.2_c_rp) > tol .or. &
                n_macro_ptr /= 300_c_ip .or. n_micro_ptr /= 200_c_ip .or. &
                abs(global_red_factor_ptr - 1e-2_c_rp) > tol .or. &
                abs(local_red_factor_ptr - 1e-3_c_rp) > tol .or. &
                seed_ptr /= 33_c_ip .or. verbose_ptr /= 3_c_ip) then
                write (stderr, *) "test_solver_py_interface failed: Passed "// &
                    "optional arguments associated with wrong values."
                test_solver_interface = .false.
            end if

            ! get Fortran pointer to passed logging function and call it
            message = "test" // c_null_char
            call c_f_procpointer(cptr=logger_c_funptr, fptr=logger_funptr)
            call logger_funptr(message)
        end if

        ! set return arguments
        error_c = 0

        ! move on to optional argument test for next test
        solver_default = .false.

    end function mock_solver_c_wrapper

    function mock_stability_check_c_wrapper(h_diag_c_ptr, hess_x_c_funptr, n_param_c, &
                                            stable_c, kappa_c_ptr, precond_c_funptr, &
                                            jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                            n_random_trial_vectors_c_ptr, &
                                            n_iter_c_ptr, verbose_c_ptr, &
                                            logger_c_funptr) result(error_c) &
        bind(C, name="mock_stability_check")
        !
        ! this subroutine is a mock routine for the stability check C wrapper
        ! subroutine
        !
        type(c_ptr), intent(in), value :: h_diag_c_ptr, kappa_c_ptr, &
                                          jacobi_davidson_c_ptr, conv_tol_c_ptr, &
                                          n_random_trial_vectors_c_ptr, n_iter_c_ptr, &
                                          verbose_c_ptr
        type(c_funptr), intent(in), value :: hess_x_c_funptr, precond_c_funptr, &
                                             logger_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        logical(c_bool), intent(out) :: stable_c
        integer(c_ip) :: error_c

        real(c_rp), pointer :: h_diag_ptr(:), kappa_ptr(:), conv_tol_ptr
        procedure(hess_x_c_type), pointer :: hess_x_funptr
        procedure(precond_c_type), pointer :: precond_funptr
        real(c_rp), allocatable :: x(:), hess_x(:), residual(:), precond_residual(:)
        logical, pointer :: jacobi_davidson_ptr
        integer(c_ip), pointer :: n_random_trial_vectors_ptr, n_iter_ptr, verbose_ptr
        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message
        integer(c_ip) :: error

        ! check if Hessian diagonal is passed correctly
        call c_f_pointer(h_diag_c_ptr, h_diag_ptr, [n_param_c])
        if (any(abs(h_diag_ptr - 3.0_c_rp) > tol)) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "Hessian diagonal wrong."
            test_stability_check_interface = .false.
        end if

        ! get Fortran pointer to Hessian linear transformation function, call it and 
        ! check result
        call c_f_procpointer(cptr=hess_x_c_funptr, fptr=hess_x_funptr)
        allocate(x(n_param_c), hess_x(n_param_c))
        x = 1.0_c_rp
        error = hess_x_funptr(x, hess_x)
        if (error /= 0) then
            write (stderr, *) "test_stability_check_py_interface failed: Passed "// &
                "Hessian linear transformation function returned error."
            test_stability_check_interface = .false.
        end if
        if (any(abs(hess_x - 4.0_c_rp) > tol)) then
            write (stderr, *) "test_stability_check_py_interface failed: Returned "// &
                "Hessian linear transformation from passed Hessian linear "// &
                "transformation function wrong."
            test_stability_check_interface = .false.
        end if
        deallocate(x, hess_x)

        ! check if check if default arguments are correctly unassociated with values
        if (stability_check_default) then
            if (c_associated(precond_c_funptr) .or. &
                c_associated(jacobi_davidson_c_ptr) .or. c_associated(conv_tol_c_ptr) &
                .or. c_associated(n_random_trial_vectors_c_ptr) .or. &
                c_associated(n_iter_c_ptr) .or. c_associated(verbose_c_ptr) .or. &
                c_associated(logger_c_funptr)) then
                write (stderr, *) "test_stability_check_py_interface failed: "// &
                    "Passed optional arguments associated with values."
                test_stability_check_interface = .false.
            end if
            ! check if check if optional arguments are associated with correct values
        else
            if (.not. (c_associated(precond_c_funptr) .and. &
                       c_associated(jacobi_davidson_c_ptr) .and. &
                       c_associated(conv_tol_c_ptr) .and. &
                       c_associated(n_random_trial_vectors_c_ptr) .and. &
                       c_associated(n_iter_c_ptr) .and. &
                       c_associated(verbose_c_ptr) .and. &
                       c_associated(logger_c_funptr))) then
                write (stderr, *) "test_stability_check_py_interface failed: "// &
                    "Passed optional arguments not associated with values."
                test_stability_check_interface = .false.
            end if

            ! get Fortran pointer to passed precond function, call it and check result
            call c_f_procpointer(cptr=precond_c_funptr, fptr=precond_funptr)
            allocate(residual(n_param_c), precond_residual(n_param_c))
            residual = 1.0_c_rp
            error = precond_funptr(residual, 5.0_c_rp, precond_residual)
            if (error /= 0) then
                write (stderr, *) "test_stability_check_py_interface failed: "// &
                    "Passed preconditioner function returned error."
                test_stability_check_interface = .false.
            end if
            if (any(abs(precond_residual - 5.0_c_rp) > tol)) then
                write (stderr, *) "test_stability_check_py_interface failed: "// &
                    "Returned preconditioned residual from preconditioner function "// &
                    "wrong."
                test_stability_check_interface = .false.
            end if
            deallocate(residual, precond_residual)

            ! get Fortran pointers to optional arguments and check against reference 
            ! values
            call c_f_pointer(cptr=conv_tol_c_ptr, fptr=conv_tol_ptr)
            call c_f_pointer(cptr=jacobi_davidson_c_ptr, fptr=jacobi_davidson_ptr)
            call c_f_pointer(cptr=n_random_trial_vectors_c_ptr, &
                             fptr=n_random_trial_vectors_ptr)
            call c_f_pointer(cptr=n_iter_c_ptr, fptr=n_iter_ptr)
            call c_f_pointer(cptr=verbose_c_ptr, fptr=verbose_ptr)
            if (jacobi_davidson_ptr .or. abs(conv_tol_ptr - 1e-3_c_rp) > tol .or. &
                n_random_trial_vectors_ptr /= 3_c_ip .or. n_iter_ptr /= 50_c_ip .or. &
                verbose_ptr /= 3_c_ip) then
                write (stderr, *) "test_stability_check_py_interface failed: "// &
                    "Passed optional arguments associated with wrong values."
                test_stability_check_interface = .false.
            end if

            ! get Fortran pointer to passed logging function and call it
            message = "test" // c_null_char
            call c_f_procpointer(cptr=logger_c_funptr, fptr=logger_funptr)
            call logger_funptr(message)
        end if

        ! set return arguments
        stable_c = .false.
        if (c_associated(kappa_c_ptr)) then
            call c_f_pointer(kappa_c_ptr, kappa_ptr, [n_param_c])
            kappa_ptr = 1.0_c_rp
        end if
        error_c = 0

        ! move on to optional argument test for next test
        stability_check_default = .false.

    end function mock_stability_check_c_wrapper

end module c_interface_mock
