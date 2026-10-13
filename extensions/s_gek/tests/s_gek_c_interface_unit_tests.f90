! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface_unit_tests

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol_c, n_param, n_param_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_associated, c_loc

    implicit none

    logical :: test_passed_mock_change_reference = .false.
    character(len=:), allocatable :: test_name

contains

    function mock_change_reference(new_ref_c, n_points_c, kappa_list_c, &
                                   local_grad_list_c, grad_list_c, context_c) &
        result(error) bind(C)
        !
        ! this subroutine is a mock subroutine for the reference change C function
        !
        use test_reference, only: check_host_context_c
        use otr_s_gek_test_reference, only: n_points

        real(c_rp), intent(in), target :: new_ref_c(*)
        integer(c_ip), intent(in), value :: n_points_c
        real(c_rp), intent(inout), target :: kappa_list_c(*), local_grad_list_c(*), &
                                             grad_list_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        ! initialize logical
        test_passed_mock_change_reference = .true.

        if (n_points_c /= n_points) then
            test_passed_mock_change_reference = .false.
            write(stderr, *) test_name//" failed: Number of points inside given "// &
                "change reference function wrong."
        end if

        if (any((new_ref_c(:n_param) - 1.0_c_rp) > tol_c)) then
            test_passed_mock_change_reference = .false.
            write(stderr, *) test_name//" failed: New reference parameters inside "// &
                "given change reference function wrong."
        end if

        kappa_list_c(:n_param * n_points_c) = 2.0_c_rp * &
                                              kappa_list_c(:n_param * n_points_c)
        local_grad_list_c(:n_param * n_points_c) = &
            3.0_c_rp * local_grad_list_c(:n_param * n_points_c)
        grad_list_c(:n_param * n_points_c) = 4.0_c_rp * &
                                             grad_list_c(:n_param * n_points_c)

        error = 0_c_ip

    end function mock_change_reference

    logical(c_bool) function test_s_gek_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the S-GEK orbital update factory
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_object
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, s_gek_factory, &
                                         s_gek_factory_c_wrapper, s_gek_callbacks, &
                                         s_gek_holder
        use otr_s_gek_mock, only: mock_s_gek_factory, mock_s_gek_factory_error, &
                                  test_passed_mock_s_gek_factory => test_passed
        use otr_s_gek_test_reference, only: assignment(=), ref_s_gek_settings
        use otr_common_c_interface_unit_tests, only: check_c_factory_wiring
        use c_interface_unit_tests, only: mock_update_orbs, mock_obj_func
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: host_context_type, arm_host_context_c, host_context, &
                                  host_context_reached

        type(c_funptr) :: update_orbs_s_gek_c_funptr, obj_func_s_gek_c_funptr
        type(solver_settings_type_c) :: solver_settings_c
        type(s_gek_settings_type_c) :: settings
        integer(c_ip) :: error_c
        type(host_context_type), target :: other_host_context

        ! assume tests pass
        test_s_gek_factory_c_wrapper = .true.

        ! set test name
        test_name = "test_s_gek_factory_c_wrapper"

        ! inject mock function
        s_gek_factory => mock_s_gek_factory

        ! associate optional settings with values
        settings = ref_s_gek_settings
        settings%logger = c_funloc(mock_logger_c)

        ! arm the host context of the settings, which the logging and change of
        ! reference functions receive, and of the solver settings, which the orbital
        ! update and objective function receive
        call arm_host_context_c(settings%context)
        solver_settings_c = default_solver_settings
        solver_settings_c%context = c_loc(host_context)

        ! call S-GEK orbital updating factory C wrapper
        test_passed_mock_s_gek_factory = .false.
        test_passed_mock_change_reference = .false.
        mock_s_gek_factory_error = 0
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, solver_settings_c, settings)

        ! check if the mock factory received the correct input and called the logging
        ! function, and if the given change of reference function received the correct
        ! input
        if (.not. test_passed_mock_s_gek_factory) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Mock factory "// &
                "received wrong input."
            test_s_gek_factory_c_wrapper = .false.
        end if
        if (.not. test_passed_mock_change_reference) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Given change of "// &
                "reference function received wrong input."
            test_s_gek_factory_c_wrapper = .false.
        end if
        if (.not. host_context%logger_called) then
            test_s_gek_factory_c_wrapper = .false.
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Called logging "// &
                "subroutine wrong."
        end if

        ! check that the host context reached the callback functions unchanged
        if (.not. host_context_reached("s_gek_factory_c_wrapper")) &
            test_s_gek_factory_c_wrapper = .false.

        ! check if output variables are as expected
        if (error_c /= 0) then
            test_s_gek_factory_c_wrapper = .false.
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Returned error "// &
                "code wrong."
        end if

        ! check if the callback slot points to the given function
        if (.not. associated(s_gek_callbacks%change_reference, mock_change_reference)) &
            then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Callback slot "// &
                "not set."
            test_s_gek_factory_c_wrapper = .false.
        end if

        ! check that the factory returned the C wrappers of the object and wired the
        ! solver settings with its holder
        if (.not. check_c_factory_wiring("s_gek_factory_c_wrapper", s_gek_holder, &
                                         s_gek_object, update_orbs_s_gek_c_funptr, &
                                         obj_func_s_gek_c_funptr, solver_settings_c)) &
            test_s_gek_factory_c_wrapper = .false.

        ! call the factory with the solver settings it wired and another host context,
        ! which has to produce an error without calling the factory, returning
        ! functions or replacing the bundle holding the callback functions of the
        ! object
        settings%context = c_loc(other_host_context)
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, solver_settings_c, settings)
        if (error_c == 0 .or. c_associated(update_orbs_s_gek_c_funptr) .or. &
            c_associated(obj_func_s_gek_c_funptr) .or. &
            .not. c_associated(s_gek_callbacks%host_context, c_loc(host_context))) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Solver settings "// &
                "the holder was wired into not refused before the object was changed."
            test_s_gek_factory_c_wrapper = .false.
        end if

        ! call the factory with a mock factory which produces an error, which has to be
        ! returned without functions or writing back the solver settings or changing
        ! the holder or the bundle of the object
        solver_settings_c = default_solver_settings
        update_orbs_s_gek_c_funptr = c_funloc(mock_logger_c)
        obj_func_s_gek_c_funptr = c_funloc(mock_logger_c)
        mock_s_gek_factory_error = 1
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, solver_settings_c, settings)
        mock_s_gek_factory_error = 0
        if (error_c /= 1 .or. c_associated(update_orbs_s_gek_c_funptr) .or. &
            c_associated(obj_func_s_gek_c_funptr) .or. &
            c_associated(solver_settings_c%context)) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Error of the "// &
                "factory not passed on."
            test_s_gek_factory_c_wrapper = .false.
        end if
        if (.not. associated(s_gek_holder%solver_level%level, &
                             s_gek_object%solver_level) .or. &
            .not. c_associated(s_gek_callbacks%host_context, c_loc(host_context))) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Holder or "// &
                "callback functions of the object not restored after an error of "// &
                "the factory."
            test_s_gek_factory_c_wrapper = .false.
        end if

    end function test_s_gek_factory_c_wrapper

    logical(c_bool) function test_change_reference_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the reference change
        !
        use otr_s_gek, only: change_reference_type
        use otr_s_gek_c_interface, only: s_gek_callbacks, change_reference_f_wrapper
        use otr_s_gek_test_reference, only: test_change_reference_funptr
        use test_reference, only: arm_host_context_c, host_context_reached

        procedure(change_reference_type), pointer :: change_reference_funptr
        class(*), pointer :: context

        ! assume tests pass
        test_change_reference_f_wrapper = .true.

        ! set test name
        test_name = "test_change_reference_f_wrapper"

        ! inject mock C function and hand the bundle to the wrapper as its context
        s_gek_callbacks%change_reference => mock_change_reference
        call arm_host_context_c(s_gek_callbacks%host_context)
        context => s_gek_callbacks

        ! get pointer to subroutine
        change_reference_funptr => change_reference_f_wrapper

        ! clear the flag of the mock
        test_passed_mock_change_reference = .false.

        ! test change reference Fortran wrapper
        if (.not. test_change_reference_funptr( &
            change_reference_funptr, "change_reference_f_wrapper", "", context)) &
            test_change_reference_f_wrapper = .false.
        if (.not. test_passed_mock_change_reference) &
            test_change_reference_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("change_reference_f_wrapper")) &
            test_change_reference_f_wrapper = .false.

    end function test_change_reference_f_wrapper

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

        if (.not. test_passed) then
            test_update_orbs_s_gek_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_update_orbs_s_gek_deconstructor_c_wrapper "// &
                "failed: Deconstructor called wrong."
        end if

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
        use c_interface, only: c_callbacks_type
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(s_gek_settings_type_c) :: settings_c
        type(s_gek_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_s_gek_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_s_gek_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the factory does
        settings = settings_c
        callbacks%logger => mock_logger_c
        settings%context => callbacks

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_s_gek_f_c = .false.
            write(stderr, *) "test_assign_s_gek_f_c failed: Logging function not "// &
                "associated with value."
        else
            call arm_host_context_c(callbacks%host_context)
            call settings%logger("test", settings%context)
            if (.not. host_context%logger_called) then
                test_assign_s_gek_f_c = .false.
                write(stderr, *) "test_assign_s_gek_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
            if (.not. host_context_reached("assign_s_gek_f_c")) &
                test_assign_s_gek_f_c = .false.
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

        ! convert C settings that were not initialized and check that they are left
        ! uninitialized without a logging function, so that the factory initializes
        ! them after its warning
        settings_c%initialized = .false.
        settings = settings_c
        if (settings%initialized .or. associated(settings%logger)) then
            write(stderr, *) "test_assign_s_gek_f_c failed: Settings that were not "// &
                "initialized converted."
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

        type(s_gek_settings_type) :: settings
        type(s_gek_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_s_gek_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_s_gek_settings

        ! convert to C settings
        settings_c = settings

        ! check that no callback function pointers and host contexts are associated
        if (c_associated(settings_c%logger)) then
            test_assign_s_gek_c_f = .false.
            write(stderr, *) "test_assign_s_gek_c_f failed: Logger function associated."
        end if
        if (c_associated(settings_c%context)) then
            write(stderr, *) "test_assign_s_gek_c_f failed: Host context associated."
            test_assign_s_gek_c_f = .false.
        end if

        ! check against reference values
        if (settings_c /= ref_s_gek_settings) then
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
