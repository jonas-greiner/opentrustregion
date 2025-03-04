program testsuite

    use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
    use trustorbopt, only: ip
    use trustorbopt_unit_tests, only: test_solver, test_stability_check, &
                                      test_newton_step, test_bisection, test_bracket, &
                                      test_extend_symm_matrix, test_add_column, &
                                      test_symm_mat_min_eig, test_min_eigval, &
                                      test_init_rng, test_generate_trial_vectors, &
                                      test_gram_schmidt, test_init_solver_settings, &
                                      test_init_solver_settings, &
                                      test_init_stability_settings, test_set_default, &
                                      test_raise_error
    use c_interface_unit_tests, only: test_solver_c_wrapper, &
                                      test_stability_check_c_wrapper, &
                                      test_update_orbs_c_wrapper, &
                                      test_hess_x_c_wrapper, test_obj_func_c_wrapper, &
                                      test_set_default_c_ptr
    use trustorbopt_system_tests, only: test_h2o_atomic_fb, test_h2o_saddle_fb

    implicit none

    integer(ip) :: failures = 0, python_failure
#ifdef PYTHON_TEST_SCRIPT
    character(len=*), parameter :: cmd = PYTHON_TEST_SCRIPT
#endif

    ! run all tests and print results
    write (stdout, *) repeat("-", 50)
    write (stdout, *) "Running unit tests for trustorbopt..."
    write (stdout, *) repeat("-", 50)

    call run_test(test_raise_error, "raise_error")
    call run_test(test_set_default, "set_default")
    call run_test(test_init_solver_settings, "init_solver_settings")
    call run_test(test_init_stability_settings, "init_stability_settings")
    call run_test(test_init_rng, "init_rng")
    call run_test(test_gram_schmidt, "gram_schmidt")
    call run_test(test_generate_trial_vectors, "generate_trial_vectors")
    call run_test(test_add_column, "add_column")
    call run_test(test_extend_symm_matrix, "extend_symm_matrix")
    call run_test(test_min_eigval, "min_eigval")
    call run_test(test_symm_mat_min_eig, "symm_mat_min_eig")
    call run_test(test_bracket, "bracket")
    call run_test(test_newton_step, "newton_step")
    call run_test(test_bisection, "bisection")
    call run_test(test_solver, "solver")
    call run_test(test_stability_check, "stability_check")

    write (stdout, *) repeat("-", 50)
    write (stdout, *) "Running unit tests for C interface..."
    write (stdout, *) repeat("-", 50)

    call run_test(test_solver_c_wrapper, "solver_c_wrapper")
    call run_test(test_stability_check_c_wrapper, "stability_check_c_wrapper")
    call run_test(test_update_orbs_c_wrapper, "update_orbs_c_wrapper")
    call run_test(test_hess_x_c_wrapper, "hess_x_c_wrapper")
    call run_test(test_obj_func_c_wrapper, "obj_func_c_wrapper")
    call run_test(test_set_default_c_ptr, "set_default_c_ptr")

    ! only run python tests when python distribution and dependencies are present
#ifdef PYTHON_TEST_SCRIPT
    write (stdout, *) repeat("-", 50)
    write (stdout, *) "Running unit tests for Python interface..."
    write (stdout, *) repeat("-", 50)

    call execute_command_line(cmd, exitstat=python_failures)
    failures = failures + python_failures
#endif

    write (stdout, *) repeat("-", 50)
    write (stdout, *) "Running system tests for trustorbopt..."
    write (stdout, *) repeat("-", 50)

    call run_test(test_h2o_atomic_fb, "h2o_atomic_fb")
    call run_test(test_h2o_saddle_fb, "h2o_saddle_fb")

    write (stdout, *) repeat("-", 50)
    if (failures == 1) then
        write (stdout, '(I0, A)') failures, " test failed."
        error stop
    else if (failures > 1) then
        write (stdout, '(I0, A)') failures, " tests failed."
        error stop
    else
        write (stdout, *) "All tests passed!"
    end if
    write (stdout, *) repeat("-", 50)

contains

    subroutine run_test(test_func, test_name)
        !
        ! this subroutine runs tests and prints whether these have failed or passed.
        !
        logical :: test_func
        character(len=*), intent(in) :: test_name

        if (.not. test_func()) then
            failures = failures + 1
            write (stdout, *) trim(test_name), " FAILED"
        else
            write (stdout, *) trim(test_name), " PASSED"
        end if

    end subroutine run_test

end program testsuite
