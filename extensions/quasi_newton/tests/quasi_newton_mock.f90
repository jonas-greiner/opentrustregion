! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_mock

    use opentrustregion, only: ip, stderr
    use otr_qn, only: qn_factory
    use otr_qn_test_reference, only: ref_qn_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! error the mock factory returns
    integer(ip) :: mock_qn_factory_error = 0

    ! whether the mock factory expects the default settings without a logging function,
    ! which the C wrapper hands it for C settings that were not initialized and without
    ! any, instead of the reference settings
    logical :: mock_qn_factory_defaults = .false.

    ! create function pointers to ensure that routines comply with interface
    procedure(qn_factory), pointer :: mock_qn_factory_ptr => mock_qn_factory

contains

    subroutine mock_qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                               transport_funptr, init_hess_funptr, &
                               update_orbs_qn_funptr, obj_func_qn_funptr, &
                               solver_settings, solver_state, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function and objective function and wires the solver
        ! settings, which returns the error it is set to and otherwise creates a
        ! quasi-Newton object reading the given state, points the context of the solver
        ! settings at its level of the solver and lets the settings own it, as the
        ! factory does
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, solver_state_type
        use otr_common, only: extension_type, extension_solver_settings_type
        use otr_qn, only: transport_type, init_hess_type, qn_settings_type, &
                          sr1_updating_type, default_qn_settings
        use test_reference, only: check_update_orbs_funptr, check_obj_func_funptr, &
                                  n_param_ref => n_param
        use otr_qn_test_reference, only: check_transport_funptr, check_init_hess_funptr
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_common_unit_tests, only: check_mock_factory_logger, &
                                         setup_owning_solver_settings

        integer(ip), intent(in) :: n_param
        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        procedure(transport_type), intent(in), pointer :: transport_funptr
        procedure(init_hess_type), intent(in), pointer :: init_hess_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_qn_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_qn_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(qn_settings_type), intent(in), optional :: settings

        class(extension_type), allocatable, target :: object

        ! initialize logical
        test_passed = .true.

        ! check number of parameters
        if (n_param /= n_param_ref) then
            test_passed = .false.
            write(stderr, *) "test_qn_factory_c_wrapper failed: Passed number of "// &
                "parameters wrong."
        end if

        ! test passed orbital update and objective function, which call the host's C
        ! functions through the C wrappers of the core, so they find them with the host
        ! context in the callback bundle the converted solver settings point to, as the
        ! object does through the host record of its level of the solver
        if (.not. check_update_orbs_funptr( &
            update_orbs_funptr, "qn_factory_c_wrapper", " by given orbital "// &
            "updating function", solver_settings%context)) test_passed = .false.
        if (.not. check_obj_func_funptr(obj_func_funptr, "qn_factory_c_wrapper", &
                                        " by given objective function", &
                                        solver_settings%context)) test_passed = .false.

        ! the C wrapper always passes settings, since their context points at the
        ! bundle of the host's callback functions, which the checks below read
        if (.not. present(settings)) then
            test_passed = .false.
            write(stderr, *) "test_qn_factory_c_wrapper failed: Settings not passed."
            error = 1
            update_orbs_qn_funptr => null()
            obj_func_qn_funptr => null()
            return
        end if

        ! test passed transport function
        if (.not. check_transport_funptr(transport_funptr, "qn_factory_c_wrapper", &
                                         " by given transport function", &
                                         settings%context)) test_passed = .false.

        ! test passed initial Hessian function
        if (.not. check_init_hess_funptr(init_hess_funptr, "qn_factory_c_wrapper", &
                                         " by given initial Hessian function", &
                                         settings%context)) test_passed = .false.

        ! check if optional logging function is correctly passed
        if (.not. check_mock_factory_logger(settings, mock_qn_factory_defaults, &
                                            "test_qn_factory_c_wrapper")) &
            test_passed = .false.

        ! check if optional settings are correctly passed
        if (mock_qn_factory_defaults) then
            if (settings /= default_qn_settings) then
                test_passed = .false.
                write(stderr, *) "test_qn_factory_c_wrapper failed: Settings that "// &
                    "were not initialized not initialized to the default values."
            end if
        else if (settings /= ref_qn_settings) then
            test_passed = .false.
            write(stderr, *) "test_qn_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = mock_qn_factory_error
        update_orbs_qn_funptr => null()
        obj_func_qn_funptr => null()
        if (error /= 0) return
        allocate(sr1_updating_type :: object)
        object%n_param = n_param
        object%state => solver_state
        call setup_owning_solver_settings(solver_settings, object)
        update_orbs_qn_funptr => mock_update_orbs
        obj_func_qn_funptr => mock_obj_func

    end subroutine mock_qn_factory

end module otr_qn_mock
