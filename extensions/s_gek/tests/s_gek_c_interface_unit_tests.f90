! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol, tol_c, n_param, n_param_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_funloc, c_associated

    implicit none

contains

    logical(c_bool) function test_update_orbs_s_gek_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the S-GEK orbital update factory
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, &
                                         update_orbs_s_gek_factory, &
                                         update_orbs_s_gek_factory_c_wrapper
        use otr_s_gek_mock, only: mock_update_orbs_s_gek_factory, &
                                  test_passed_mock_update_orbs_s_gek_factory => &
                                  test_passed
        use otr_s_gek_test_reference, only: assignment(=), ref_s_gek_settings
        use c_interface_unit_tests, only: mock_update_orbs, mock_logger, test_logger
        use otr_common_c_interface_mock, only: mock_change_reference, test_name, &
                                               test_passed_mock_change_reference => &
                                               test_passed
        use test_reference, only: test_update_orbs_c_funptr

        type(c_funptr) :: update_orbs_orig_c_funptr, change_reference_c_funptr, &
                          update_orbs_s_gek_c_funptr
        type(s_gek_settings_type_c) :: settings
        integer(c_ip) :: error_c

        ! assume tests pass
        test_update_orbs_s_gek_factory_c_wrapper = .true.

        ! set test name
        test_name = "test_update_orbs_s_gek_factory_c_wrapper"

        ! inject mock function
        update_orbs_s_gek_factory => mock_update_orbs_s_gek_factory

        ! get C function pointers to Fortran functions
        update_orbs_orig_c_funptr = c_funloc(mock_update_orbs)
        change_reference_c_funptr = c_funloc(mock_change_reference)

        ! associate optional settings with values
        settings = ref_s_gek_settings
        settings%logger = c_funloc(mock_logger)

        ! initialize logger logical
        test_logger = .true.

        ! call S-GEK orbital updating factory C wrapper
        error_c = update_orbs_s_gek_factory_c_wrapper(update_orbs_orig_c_funptr, &
                                                      change_reference_c_funptr, &
                                                      n_param_c, settings, &
                                                      update_orbs_s_gek_c_funptr)

        ! check if logging subroutine was correctly called
        if (.not. test_logger) then
            test_update_orbs_s_gek_factory_c_wrapper = .false.
            write(stderr, *) "test_update_orbs_s_gek_factory_c_wrapper failed: "// &
                "Called logging subroutine wrong."
        end if

        ! check if output variables are as expected
        if (error_c /= 0) then
            test_update_orbs_s_gek_factory_c_wrapper = .false.
            write(stderr, *) "test_update_orbs_s_gek_factory_c_wrapper failed: "// &
                "Returned error code wrong."
        end if

        ! test returned orbital update function
        test_update_orbs_s_gek_factory_c_wrapper = &
            test_update_orbs_s_gek_factory_c_wrapper .and. &
            test_update_orbs_c_funptr(update_orbs_s_gek_c_funptr, &
                                      "update_orbs_s_gek_factory_c_wrapper", &
                                      " by returned orbital updating function") .and. &
            test_passed_mock_update_orbs_s_gek_factory .and. &
            test_passed_mock_change_reference

    end function test_update_orbs_s_gek_factory_c_wrapper

    logical(c_bool) function test_update_orbs_orig_s_gek_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the orbital update
        !
        use opentrustregion, only: update_orbs_type
        use otr_s_gek_c_interface, only: update_orbs_orig_s_gek_before_wrapping, &
                                         update_orbs_orig_s_gek_f_wrapper
        use c_interface_unit_tests, only: mock_update_orbs
        use test_reference, only: test_update_orbs_funptr

        procedure(update_orbs_type), pointer :: update_orbs_funptr

        ! inject mock subroutine
        update_orbs_orig_s_gek_before_wrapping => mock_update_orbs

        ! get pointer to subroutine
        update_orbs_funptr => update_orbs_orig_s_gek_f_wrapper

        ! test orbital update Fortran wrapper
        test_update_orbs_orig_s_gek_f_wrapper = &
            test_update_orbs_funptr(update_orbs_funptr, &
                                    "update_orbs_orig_s_gek_f_wrapper", "")

    end function test_update_orbs_orig_s_gek_f_wrapper

    logical(c_bool) function test_change_reference_s_gek_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the reference change
        !
        use otr_s_gek, only: change_reference_type
        use otr_s_gek_c_interface, only: change_reference_s_gek_before_wrapping, &
                                         change_reference_f_wrapper
        use otr_common_c_interface_mock, only: mock_change_reference, test_name, &
                                               test_passed
        use otr_common_test_reference, only: test_change_reference_funptr

        procedure(change_reference_type), pointer :: change_reference_funptr

        ! set test name
        test_name = "test_change_reference_s_gek_f_wrapper"

        ! inject mock subroutine
        change_reference_s_gek_before_wrapping => mock_change_reference

        ! get pointer to subroutine
        change_reference_funptr => change_reference_f_wrapper

        ! test change reference Fortran wrapper
        test_change_reference_s_gek_f_wrapper = &
            test_change_reference_funptr(change_reference_funptr, &
                                         "change_reference_s_gek_f_wrapper", "") .and. &
            test_passed

    end function test_change_reference_s_gek_f_wrapper

    logical(c_bool) function test_update_orbs_s_gek_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the S-GEK orbital update
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_s_gek_c_interface, only: update_orbs_s_gek_before_wrapping, &
                                         update_orbs_s_gek_c_wrapper
        use otr_common_mock, only: mock_update_orbs
        use test_reference, only: test_update_orbs_c_funptr

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        update_orbs_s_gek_before_wrapping => mock_update_orbs

        ! test orbital update C wrapper
        test_update_orbs_s_gek_c_wrapper = &
            test_update_orbs_c_funptr(c_funloc(update_orbs_s_gek_c_wrapper), &
                                      "update_orbs_s_gek_c_wrapper", "")

    end function test_update_orbs_s_gek_c_wrapper

    logical(c_bool) function test_hess_x_s_gek_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the S-GEK Hessian linear transformation
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_s_gek_c_interface, only: hess_x_s_gek_before_wrapping, &
                                         hess_x_s_gek_c_wrapper
        use otr_common_mock, only: mock_hess_x
        use test_reference, only: test_hess_x_c_funptr

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        hess_x_s_gek_before_wrapping => mock_hess_x

        ! test orbital update C wrapper
        test_hess_x_s_gek_c_wrapper = &
            test_hess_x_c_funptr(c_funloc(hess_x_s_gek_c_wrapper), &
                                 "hess_x_s_gek_c_wrapper", "")

    end function test_hess_x_s_gek_c_wrapper

    logical(c_bool) function test_init_s_gek_settings_c() bind(C)
        !
        ! this function tests that the S-GEK settings initialization routine correctly 
        ! initializes all settings to their default values
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, init_s_gek_settings_c
        use otr_s_gek, only: default_s_gek_settings
        use otr_s_gek_test_reference, only: operator(/=)

        type(s_gek_settings_type_c) :: settings

        ! assume test passes
        test_init_s_gek_settings_c = .true.

        ! initialize settings
        call init_s_gek_settings_c(settings)

        ! check function pointers
        if (c_associated(settings%logger)) then
            write(stderr, *) "test_init_s_gek_settings_c failed: Function pointers "// &
                "should not be initialized."
            test_init_s_gek_settings_c = .false.
        end if

        ! check settings
        if (settings /= default_s_gek_settings) then
            write(stderr, *) "test_init_s_gek_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_s_gek_settings_c = .false. 
        end if

    end function test_init_s_gek_settings_c

    logical(c_bool) function test_update_orbs_s_gek_deconstructor_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the S-GEK deconstructor
        !
        use otr_s_gek_c_interface, only: update_orbs_s_gek_deconstructor, &
                                         update_orbs_s_gek_deconstructor_c_wrapper
        use otr_s_gek_mock, only: mock_update_orbs_s_gek_deconstructor, test_passed

        ! assume tests pass
        test_update_orbs_s_gek_deconstructor_c_wrapper = .true.

        ! inject mock function
        update_orbs_s_gek_deconstructor => mock_update_orbs_s_gek_deconstructor

        ! initialize test logical
        test_passed = .false.

        ! call S-GEK orbital updating deconstructor C wrapper
        call update_orbs_s_gek_deconstructor_c_wrapper()

        ! check if test has passed
        test_update_orbs_s_gek_deconstructor_c_wrapper = test_passed

    end function test_update_orbs_s_gek_deconstructor_c_wrapper

    logical(c_bool) function test_assign_s_gek_f_c() bind(C)
        !
        ! this function tests that the function that converts S-GEK settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_settings_type
        use otr_s_gek_test_reference, only: assignment(=), ref_s_gek_settings, &
                                            operator(/=)
        use c_interface_unit_tests, only: mock_logger, test_logger

        type(s_gek_settings_type_c) :: settings_c
        type(s_gek_settings_type) :: settings

        ! assume test passes
        test_assign_s_gek_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_s_gek_settings
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings
        settings = settings_c

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_s_gek_f_c = .false.
            write(stderr, *) "test_assign_s_gek_f_c failed: Logging function not "// &
                "associated with value."
        else
            test_logger = .true.
            call settings%logger("test")
            if (.not. test_logger) then
                test_assign_s_gek_f_c = .false.
                write(stderr, *) "test_assign_s_gek_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
        end if

        ! check against reference values
        if (settings /= ref_s_gek_settings) then
            write(stderr, *) "test_assign_s_gek_f_c failed: Settings not converted "// &
                "correctly."
            test_assign_s_gek_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            write(stderr, *) "test_assign_s_gek_f_c failed: Settings not marked as "// &
                "initialized."
            test_assign_s_gek_f_c = .false.
        end if

    end function test_assign_s_gek_f_c

    logical(c_bool) function test_assign_s_gek_c_f() bind(C)
        !
        ! this function tests that the function that converts S-GEK settings from 
        ! Fortran to C correctly performs this conversion
        !
        use otr_s_gek, only: s_gek_settings_type
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek_test_reference, only: ref_s_gek_settings, assignment(=), &
                                            operator(/=)

        type(s_gek_settings_type)   :: settings
        type(s_gek_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_s_gek_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_s_gek_settings

        ! convert to C settings
        settings_c = settings

        ! check that callback function pointers are not associated
        if (c_associated(settings_c%logger)) then
            test_assign_s_gek_c_f = .false.
            write(stderr, *) "test_assign_s_gek_c_f failed: Logger function associated."
        end if

        ! check against reference values
        if (settings /= ref_s_gek_settings) then
            write(stderr, *) "test_assign_s_gek_c_f failed: Settings not converted "// &
                "correctly."
            test_assign_s_gek_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_s_gek_c_f = .false.
            write(stderr, *) "test_assign_s_gek_c_f failed: Settings not marked as "// &
                "initialized."
        end if

    end function test_assign_s_gek_c_f

end module otr_s_gek_c_interface_unit_tests
