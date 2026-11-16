! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_mock

    use opentrustregion_core, only: ip, stderr, solver_settings_type
    use c_interface, only: c_callbacks_type
    use otr_s_gek, only: s_gek_solver, s_gek_settings_type
    use otr_s_gek_c_interface, only: s_gek_c_callbacks_type

    implicit none

    logical :: test_passed

    ! settings and callback bundles the mock S-GEK solver received, those of the solver
    ! settings and of the S-GEK settings, so that a test can check what the C wrapper
    ! handed over
    type(solver_settings_type) :: received_solver_settings
    type(s_gek_settings_type) :: received_s_gek_settings
    type(c_callbacks_type) :: received_callbacks
    type(s_gek_c_callbacks_type) :: received_s_gek_callbacks

    ! whether the callback bundles of the solver settings and of the S-GEK settings the
    ! mock S-GEK solver received refer to the state record it was given, so that a test
    ! can check that the C wrapper copies this record into the host's before every
    ! callback function of either settings
    logical :: solver_state_bundled, s_gek_state_bundled

    ! create function pointers to ensure that routines comply with interface
    procedure(s_gek_solver), pointer :: mock_s_gek_solver_ptr => mock_s_gek_solver

contains

    subroutine mock_s_gek_solver(update_orbs_funptr, obj_func_funptr, &
                                 change_reference_funptr, n_param, error, settings, &
                                 state, s_gek_settings)
        !
        ! this subroutine is a mock routine for the S-GEK solver to test the C
        ! interface, which returns the error code of the mock solver and writes its
        ! state
        !
        use opentrustregion_core, only: update_orbs_type, obj_func_type, &
                                        solver_state_type
        use opentrustregion_mock, only: mock_error, mock_solver_state
        use test_reference, only: check_update_orbs_funptr, check_obj_func_funptr, &
                                  check_hess_x_funptr, n_param_ref => n_param
        use otr_common_test_reference, only: check_obj_func_grad_funptr
        use otr_s_gek, only: change_reference_type
        use otr_s_gek_test_reference, only: check_change_reference_funptr

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        procedure(change_reference_type), intent(in), pointer :: change_reference_funptr
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        class(solver_settings_type), intent(in), optional :: settings
        type(solver_state_type), intent(out), target, optional :: state
        type(s_gek_settings_type), intent(in), optional :: s_gek_settings

        ! assume tests pass
        test_passed = .true.

        ! set error code and check that the C wrapper passes the settings it converted
        error = mock_error
        if (.not. (present(settings) .and. present(s_gek_settings))) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Settings not passed."
            test_passed = .false.
            return
        end if

        ! test passed orbital update subroutine and objective function, which find the
        ! C functions in the bundle of the solver settings
        if (.not. check_update_orbs_funptr( &
            update_orbs_funptr, "s_gek_solver_c_wrapper", " by given orbital "// &
            "updating subroutine", settings%context)) test_passed = .false.
        if (.not. check_obj_func_funptr(obj_func_funptr, "s_gek_solver_c_wrapper", &
                                        " by given objective function", &
                                        settings%context)) test_passed = .false.

        ! test passed chart transformation, which finds the C function in the bundle of
        ! the S-GEK settings
        if (.not. check_change_reference_funptr( &
            change_reference_funptr, "s_gek_solver_c_wrapper", " by given chart "// &
            "transformation", s_gek_settings%context)) test_passed = .false.

        ! check number of parameters
        if (n_param /= n_param_ref) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Passed number of "// &
                "parameters wrong."
            test_passed = .false.
        end if

        ! check the callback functions of the S-GEK settings that were passed, which
        ! find the C functions in the bundle of the S-GEK settings, which of them are
        ! passed the test checks
        if (associated(s_gek_settings%obj_func_grad)) then
            if (.not. check_obj_func_grad_funptr( &
                s_gek_settings%obj_func_grad, "s_gek_solver_c_wrapper", &
                " by the objective function of the S-GEK settings", &
                s_gek_settings%context)) test_passed = .false.
        end if
        if (associated(s_gek_settings%init_hess)) then
            if (.not. check_hess_x_funptr( &
                s_gek_settings%init_hess, "s_gek_solver_c_wrapper", &
                " by the initial Hessian of the S-GEK settings", &
                s_gek_settings%context)) test_passed = .false.
        end if
        if (associated(s_gek_settings%logger)) &
            call s_gek_settings%logger("test", s_gek_settings%context)

        ! record the settings and the callback bundles of the solver settings and of
        ! the S-GEK settings
        received_solver_settings = settings
        received_s_gek_settings = s_gek_settings
        received_callbacks = c_callbacks_type()
        select type (callbacks => settings%context)
        type is (c_callbacks_type)
            received_callbacks = callbacks
        end select
        received_s_gek_callbacks = s_gek_c_callbacks_type()
        select type (callbacks => s_gek_settings%context)
        type is (s_gek_c_callbacks_type)
            received_s_gek_callbacks = callbacks
        end select

        ! record whether both callback bundles refer to the state record
        solver_state_bundled = .false.
        s_gek_state_bundled = .false.
        if (present(state)) then
            solver_state_bundled = associated(received_callbacks%solver_state, state)
            s_gek_state_bundled = &
                associated(received_s_gek_callbacks%solver_state, state)
        end if

        ! set state
        if (present(state)) state = mock_solver_state

    end subroutine mock_s_gek_solver

end module otr_s_gek_mock
