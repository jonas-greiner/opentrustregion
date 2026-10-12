! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_mock

    use opentrustregion, only: ip, stderr
    use otr_s_gek, only: s_gek_factory, update_orbs_s_gek_deconstructor
    use otr_s_gek_test_reference, only: ref_s_gek_settings, operator(/=)
    use, intrinsic :: iso_c_binding, only: c_loc

    implicit none

    logical :: test_passed

    ! create function pointers to ensure that routines comply with interface
    procedure(s_gek_factory), pointer :: mock_s_gek_factory_ptr => mock_s_gek_factory
    procedure(update_orbs_s_gek_deconstructor), pointer :: &
        mock_update_orbs_s_gek_deconstructor_ptr => mock_update_orbs_s_gek_deconstructor

contains

    subroutine mock_s_gek_factory(update_orbs_orig_funptr, change_reference_funptr, &
                                  n_param, error, settings, update_orbs_s_gek_funptr)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function
        !
        use opentrustregion, only: update_orbs_type
        use otr_s_gek, only: change_reference_type, s_gek_settings_type
        use test_reference, only: check_update_orbs_funptr, host_context, &
                                  n_param_ref => n_param
        use c_interface, only: c_callbacks_type
        use otr_s_gek_c_interface, only: update_orbs_orig_s_gek_before_wrapping
        use otr_s_gek_test_reference, only: test_change_reference_funptr
        use otr_common_mock, only: mock_update_orbs

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_orig_funptr
        procedure(change_reference_type), intent(in), pointer :: change_reference_funptr
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        type(s_gek_settings_type), intent(in) :: settings
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_s_gek_funptr

        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context

        ! initialize logical
        test_passed = .true.

        ! test passed orbital update function, which calls the original C function
        ! through the C wrapper of the core, so it finds it in the callback bundle the
        ! returned orbital update builds, here with the host context the test armed
        callbacks%update_orbs => update_orbs_orig_s_gek_before_wrapping
        callbacks%host_context = c_loc(host_context)
        context => callbacks
        if (.not. check_update_orbs_funptr( &
            update_orbs_orig_funptr, "s_gek_factory_c_wrapper", &
            " by given orbital updating function", context)) test_passed = .false.

        ! test passed change reference function
        if (.not. test_change_reference_funptr( &
            change_reference_funptr, "s_gek_factory_c_wrapper", " by given change "// &
            "reference function", settings%context)) test_passed = .false.

        ! check number of parameters
        if (n_param /= n_param_ref) then
            test_passed = .false.
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Passed number "// &
                "of parameters wrong."
        end if

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_s_gek_settings) then
            test_passed = .false.
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        update_orbs_s_gek_funptr => mock_update_orbs

    end subroutine mock_s_gek_factory

    subroutine mock_update_orbs_s_gek_deconstructor()
        !
        ! this subroutine is a test subroutine for the S-GEK deconstructor
        !
        test_passed = .true.

    end subroutine mock_update_orbs_s_gek_deconstructor

end module otr_s_gek_mock
