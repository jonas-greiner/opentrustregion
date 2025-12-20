! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_mock

    use opentrustregion, only: rp, ip, stderr
    use otr_qn, only: update_orbs_qn_factory, update_orbs_qn_deconstructor
    use test_reference, only: tol
    use otr_qn_test_reference, only: ref_qn_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_qn_factory), pointer :: mock_update_orbs_qn_factory_ptr => &
        mock_update_orbs_qn_factory
    procedure(update_orbs_qn_deconstructor), pointer :: &
        mock_update_orbs_qn_deconstructor_ptr => mock_update_orbs_qn_deconstructor

contains

    subroutine mock_update_orbs_qn_factory(update_orbs_orig_funptr, n_param, settings, &
                                           error, update_orbs_qn_funptr)
        !
        ! this function returns is a test function for the function which returns a 
        ! modified orbital updating function
        !
        use opentrustregion, only: update_orbs_type
        use otr_qn, only: qn_settings_type
        use test_reference, only: test_update_orbs_funptr
        use otr_common_mock, only: mock_update_orbs

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_orig_funptr
        integer(ip), intent(in) :: n_param
        type(qn_settings_type), intent(inout) :: settings
        integer(ip), intent(out) :: error
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_qn_funptr

        ! initialize logical
        test_passed = .true.

        ! test passed orbital update function
        test_passed = test_passed .and. &
            test_update_orbs_funptr(update_orbs_orig_funptr, &
                                    "update_orbs_qn_factory_c_wrapper", &
                                    " by given orbital updating function")

        ! check number of parameters
        if (n_param /= 3) then
            test_passed = .false.
            write (stderr, *) "test_update_orbs_qn_factory_c_wrapper failed: "// &
                "Passed number of parameters wrong."
        end if

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write (stderr, *) "test_update_orbs_qn_factory_c_wrapper failed: "// &
                "Passed logging function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_qn_settings) then
            test_passed = .false.
            write (stderr, *) "test_update_orbs_qn_factory_c_wrapper failed: "// &
                "Passed optional settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        update_orbs_qn_funptr => mock_update_orbs

    end subroutine mock_update_orbs_qn_factory

    subroutine mock_update_orbs_qn_deconstructor()
        !
        ! this subroutine is a test subroutine for the quasi-Newton deconstructor
        !
        test_passed = .true.

    end subroutine mock_update_orbs_qn_deconstructor

end module otr_qn_mock
