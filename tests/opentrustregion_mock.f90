! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_mock

    use opentrustregion, only: rp, ip, stderr, solver, stability_check, &
                               solver_settings_type, stability_settings_type
    use c_interface, only: c_callbacks_type
    use test_reference, only: tol

    implicit none

    logical :: test_passed

    ! output values the mock solver and stability check return or write into the
    ! settings object, so that a test can check the C wrappers hand them back
    integer(ip), parameter :: mock_error = 7, mock_n_update_orbs = 3, &
                              mock_n_hess_x = 5, mock_stability_n_hess_x = 2, &
                              mock_stability_check_n_hess_x = 4

    ! stability the mock stability check returns, which a test sets
    logical :: mock_stable = .true.

    ! minimum eigenvalue the mock stability check returns
    real(rp), parameter :: mock_min_eigval = -1.0_rp

    ! whether the mock solver received a state record, so that a test can check
    ! whether the C wrapper passes one on
    logical :: solver_state_received

    ! settings and callback bundles the mock solver and stability check received,
    ! the bundle of the settings and that of the nested settings, so that a test can
    ! check what the C wrapper handed over
    type(solver_settings_type) :: received_solver_settings
    type(stability_settings_type) :: received_stability_settings
    type(c_callbacks_type) :: received_callbacks, received_stability_callbacks

    ! whether the mock stability check received the optional direction and minimum
    ! eigenvalue, so that a test can check that the C wrapper passes these on only
    ! where the caller provides memory for them
    logical :: kappa_received, min_eigval_received

    ! create function pointers to ensure that routines comply with interface
    procedure(solver), pointer :: mock_solver_ptr => mock_solver
    procedure(stability_check), pointer :: mock_stability_check_ptr => &
        mock_stability_check

contains

    subroutine mock_solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                           settings, state)
        !
        ! this subroutine is a mock routine for solver to test the C interface
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, solver_state_type
        use test_reference, only: check_update_orbs_funptr, check_obj_func_funptr, &
                                  check_hess_x_funptr, check_precond_funptr, &
                                  check_precond_pd_funptr, check_project_funptr, &
                                  check_get_extra_trial_vectors_funptr, &
                                  check_conv_check_funptr, n_param_ref => n_param

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        type(solver_settings_type), intent(inout) :: settings
        type(solver_state_type), intent(out), target, optional :: state

        ! assume tests pass
        test_passed = .true.

        ! record whether a state record was passed
        solver_state_received = present(state)

        ! test passed orbital update subroutine
        if (.not. check_update_orbs_funptr(update_orbs_funptr, "solver_c_wrapper", &
                                           " by given orbital updating subroutine", &
                                           settings%context)) test_passed = .false.

        ! test passed objective function
        if (.not. check_obj_func_funptr(obj_func_funptr, "solver_c_wrapper", &
                                        " by given objective function", &
                                        settings%context)) test_passed = .false.

        ! check number of parameters
        if (n_param /= n_param_ref) then
            write(stderr, *) "test_solver_c_wrapper failed: Passed number of "// &
                "parameters wrong."
            test_passed = .false.
        end if

        ! check the optional callback functions that were passed, which of them are
        ! passed the test checks
        if (associated(settings%precond)) then
            if (.not. check_precond_funptr(settings%precond, "solver_c_wrapper", &
                                           " by given preconditioner subroutine", &
                                           settings%context)) test_passed = .false.
        end if
        if (associated(settings%precond_pd)) then
            if (.not. check_precond_pd_funptr( &
                settings%precond_pd, "solver_c_wrapper", &
                " by given positive-definite preconditioner subroutine", &
                settings%context)) test_passed = .false.
        end if
        if (associated(settings%project)) then
            if (.not. check_project_funptr(settings%project, "solver_c_wrapper", &
                                           " by given projection subroutine", &
                                           settings%context)) test_passed = .false.
        end if
        if (associated(settings%get_extra_trial_vectors)) then
            if (.not. check_get_extra_trial_vectors_funptr( &
                settings%get_extra_trial_vectors, "solver_c_wrapper", &
                " by given extra trial vector subroutine", settings%context)) &
                test_passed = .false.
        end if
        if (associated(settings%stability_hess_x)) then
            if (.not. check_hess_x_funptr( &
                settings%stability_hess_x, "solver_c_wrapper", &
                " by given stability check Hessian linear transformation subroutine", &
                settings%context)) test_passed = .false.
        end if
        if (associated(settings%conv_check)) then
            if (.not. check_conv_check_funptr(settings%conv_check, "solver_c_wrapper", &
                                              " by given convergence check function", &
                                              settings%context)) test_passed = .false.
        end if
        if (associated(settings%logger)) call settings%logger("test", settings%context)

        ! record the settings and the callback bundles of the settings and the nested
        ! settings
        received_solver_settings = settings
        received_callbacks = c_callbacks_type()
        select type (callbacks => settings%context)
        type is (c_callbacks_type)
            received_callbacks = callbacks
        end select
        received_stability_callbacks = c_callbacks_type()
        select type (callbacks => settings%stability_settings%context)
        type is (c_callbacks_type)
            received_stability_callbacks = callbacks
        end select

        ! set output quantities and fields
        error = mock_error
        settings%max_precision_reached = .false.
        settings%n_update_orbs = mock_n_update_orbs
        settings%n_hess_x = mock_n_hess_x
        settings%stability_settings%n_hess_x = mock_stability_n_hess_x

    end subroutine mock_solver

    subroutine mock_stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                                    kappa, min_eigval)
        !
        ! this subroutine is a mock routine for the stability check to test the C
        ! interface
        !
        use opentrustregion, only: hess_x_type
        use test_reference, only: check_hess_x_funptr, check_precond_funptr, &
                                  check_project_funptr, &
                                  check_get_extra_trial_vectors_funptr, &
                                  check_conv_check_stability_funptr, n_param

        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable
        integer(ip), intent(out) :: error
        type(stability_settings_type), intent(inout) :: settings
        real(rp), intent(out), optional :: kappa(:), min_eigval

        ! assume tests pass
        test_passed = .true.

        ! check Hessian diagonal
        if (size(h_diag) /= n_param) then
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal has wrong size."
            test_passed = .false.
        else if (any(abs(h_diag - 3.0_rp) > tol)) then
            write(stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal wrong."
            test_passed = .false.
        end if

        ! test passed Hessian linear transformation subroutine
        if (.not. check_hess_x_funptr( &
            hess_x_funptr, "stability_check_c_wrapper", " by given Hessian linear "// &
            "transformation subroutine", settings%context)) test_passed = .false.

        ! check the optional callback functions that were passed, which of them are
        ! passed the test checks
        if (associated(settings%precond)) then
            if (.not. check_precond_funptr( &
                settings%precond, "stability_check_c_wrapper", " by given "// &
                "preconditioner subroutine", settings%context)) test_passed = .false.
        end if
        if (associated(settings%project)) then
            if (.not. check_project_funptr( &
                settings%project, "stability_check_c_wrapper", " by given "// &
                "projection subroutine", settings%context)) test_passed = .false.
        end if
        if (associated(settings%approx_hess_x)) then
            if (.not. check_hess_x_funptr( &
                settings%approx_hess_x, "stability_check_c_wrapper", &
                " by given approximate Hessian linear transformation subroutine", &
                settings%context)) test_passed = .false.
        end if
        if (associated(settings%get_extra_trial_vectors)) then
            if (.not. check_get_extra_trial_vectors_funptr( &
                settings%get_extra_trial_vectors, "stability_check_c_wrapper", &
                " by given extra trial vector subroutine", settings%context)) &
                test_passed = .false.
        end if
        if (associated(settings%conv_check)) then
            if (.not. check_conv_check_stability_funptr( &
                settings%conv_check, "stability_check_c_wrapper", " by given "// &
                "convergence check function", settings%context)) test_passed = .false.
        end if
        if (associated(settings%logger)) call settings%logger("test", settings%context)

        ! record the settings, their callback bundle and which optional outputs were
        ! requested
        received_stability_settings = settings
        received_callbacks = c_callbacks_type()
        select type (callbacks => settings%context)
        type is (c_callbacks_type)
            received_callbacks = callbacks
        end select
        kappa_received = present(kappa)
        min_eigval_received = present(min_eigval)

        ! set output quantities and fields
        stable = mock_stable
        if (present(kappa)) kappa = 1.0_rp
        if (present(min_eigval)) min_eigval = mock_min_eigval
        error = mock_error
        settings%n_hess_x = mock_stability_check_n_hess_x

    end subroutine mock_stability_check

end module opentrustregion_mock
