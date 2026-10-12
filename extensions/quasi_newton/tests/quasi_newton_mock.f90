! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_mock

    use opentrustregion, only: ip, stderr
    use otr_qn, only: qn_factory, update_orbs_qn_deconstructor
    use otr_qn_test_reference, only: ref_qn_settings, operator(/=)
    use, intrinsic :: iso_c_binding, only: c_loc

    implicit none

    logical :: test_passed

    ! create function pointers to ensure that routines comply with interface
    procedure(qn_factory), pointer :: mock_qn_factory_ptr => mock_qn_factory
    procedure(update_orbs_qn_deconstructor), pointer :: &
        mock_update_orbs_qn_deconstructor_ptr => mock_update_orbs_qn_deconstructor

contains

    subroutine mock_qn_factory(update_orbs_orig_funptr, transport_funptr, &
                               init_hess_funptr, n_param, error, settings, &
                               update_orbs_qn_funptr)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function
        !
        use opentrustregion, only: update_orbs_type
        use otr_qn, only: transport_type, init_hess_type, qn_settings_type
        use test_reference, only: check_update_orbs_funptr, host_context, &
                                  n_param_ref => n_param
        use c_interface, only: c_callbacks_type
        use otr_qn_c_interface, only: update_orbs_orig_qn_before_wrapping
        use otr_qn_test_reference, only: test_transport_funptr, test_init_hess_funptr
        use otr_common_mock, only: mock_update_orbs

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_orig_funptr
        procedure(transport_type), intent(in), pointer :: transport_funptr
        procedure(init_hess_type), intent(in), pointer :: init_hess_funptr
        integer(ip), intent(in) :: n_param
        type(qn_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_qn_funptr

        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        ! initialize logical
        test_passed = .true.

        ! test passed orbital update function, which calls the original C function
        ! through the C wrapper of the core, so it finds it in the callback bundle the
        ! returned orbital update builds, here with the host context the test armed
        callbacks%update_orbs => update_orbs_orig_qn_before_wrapping
        callbacks%host_context = c_loc(host_context)
        context => callbacks
        if (.not. check_update_orbs_funptr( &
            update_orbs_orig_funptr, "qn_factory_c_wrapper", &
            " by given orbital updating function", context)) test_passed = .false.

        ! test passed transport function
        if (.not. test_transport_funptr(transport_funptr, "qn_factory_c_wrapper", &
                                        " by given transport function", &
                                        settings%context)) test_passed = .false.

        ! test passed initial Hessian function
        if (.not. test_init_hess_funptr(init_hess_funptr, "qn_factory_c_wrapper", &
                                        " by given initial Hessian function", &
                                        settings%context)) test_passed = .false.

        ! check number of parameters
        if (n_param /= n_param_ref) then
            test_passed = .false.
            write(stderr, *) "test_qn_factory_c_wrapper failed: Passed number of "// &
                "parameters wrong."
        end if

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_qn_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_qn_settings) then
            test_passed = .false.
            write(stderr, *) "test_qn_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        update_orbs_qn_funptr => mock_update_orbs

    end subroutine mock_qn_factory

    subroutine mock_update_orbs_qn_deconstructor()
        !
        ! this subroutine is a test subroutine for the quasi-Newton deconstructor
        !
        test_passed = .true.

    end subroutine mock_update_orbs_qn_deconstructor

end module otr_qn_mock
