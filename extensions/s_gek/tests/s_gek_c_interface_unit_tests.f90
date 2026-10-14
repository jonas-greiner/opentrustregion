! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol_c, n_param, n_param_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_associated, c_loc, c_f_pointer, c_null_ptr

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
        use c_interface, only: solver_settings_type_c, solver_state_type_c, &
                               assignment(=)
        use otr_common, only: extension_level_type
        use otr_common_unit_tests, only: wired_level
        use otr_common_c_interface, only: extension_holder_type
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, &
                                         s_gek_c_callbacks_type, s_gek_factory, &
                                         s_gek_factory_c_wrapper
        use otr_s_gek_mock, only: mock_s_gek_factory, mock_s_gek_factory_error, &
                                  mock_s_gek_factory_defaults, &
                                  test_passed_mock_s_gek_factory => test_passed
        use otr_s_gek_test_reference, only: get_reference_s_gek_values
        use otr_common_c_interface_unit_tests, only: check_c_factory_wiring
        use c_interface_unit_tests, only: mock_update_orbs, mock_obj_func, mock_logger
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(c_funptr) :: update_orbs_s_gek_c_funptr, obj_func_s_gek_c_funptr
        type(solver_settings_type_c) :: solver_settings_c
        type(solver_state_type_c), target :: state_c
        type(s_gek_settings_type_c) :: settings
        type(c_ptr) :: extension_c, stacked_extension_c
        type(extension_holder_type), pointer :: holder
        type(extension_level_type), pointer :: level
        integer(c_ip) :: error_c
        integer(ip) :: i_case
        character(len=34), parameter :: default_cases(2) = &
            [character(len=34) :: "settings that were not initialized", "no settings"]

        ! assume tests pass
        test_s_gek_factory_c_wrapper = .true.

        ! set test name
        test_name = "test_s_gek_factory_c_wrapper"

        ! inject mock function
        s_gek_factory => mock_s_gek_factory

        ! associate optional settings with the reference values
        call get_reference_s_gek_values(settings)
        settings%logger = c_funloc(mock_logger)

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
            obj_func_s_gek_c_funptr, extension_c, solver_settings_c, c_loc(state_c), &
            settings)

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

        ! the returned extension is the holder, whose inspection would crash without
        ! one
        if (.not. c_associated(extension_c)) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: No extension "// &
                "returned."
            test_s_gek_factory_c_wrapper = .false.
            return
        end if
        call c_f_pointer(extension_c, holder)

        ! check that the settings of the holder, which own the object, were wired with
        ! it, that it reads the state of the holder and that the holder holds the
        ! bundle of the change of reference function and the host's record of the state
        level => wired_level(holder%solver_settings)
        if (.not. associated(level)) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Settings of the "// &
                "holder not wired with an object."
            test_s_gek_factory_c_wrapper = .false.
            return
        end if
        if (.not. associated(level%extension%state, holder%state)) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Object does not "// &
                "read the state of the holder."
            test_s_gek_factory_c_wrapper = .false.
        end if
        select type (callbacks => holder%extension_callbacks)
        type is (s_gek_c_callbacks_type)
            if (.not. associated(callbacks%change_reference, mock_change_reference)) &
                then
                write(stderr, *) "test_s_gek_factory_c_wrapper failed: Callback "// &
                    "slot not set."
                test_s_gek_factory_c_wrapper = .false.
            end if
        class default
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Bundle of the "// &
                "callback functions not held."
            test_s_gek_factory_c_wrapper = .false.
        end select
        if (.not. (associated(holder%solver_level%state_c, state_c) .and. &
                   associated(holder%stability_level%state_c, state_c))) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Host's record "// &
                "of the state not held."
            test_s_gek_factory_c_wrapper = .false.
        end if

        ! check that the factory returned the C wrappers of the object and wired the
        ! solver settings with its holder
        if (.not. check_c_factory_wiring("s_gek_factory_c_wrapper", holder, &
                                         level%extension, update_orbs_s_gek_c_funptr, &
                                         obj_func_s_gek_c_funptr, solver_settings_c)) &
            test_s_gek_factory_c_wrapper = .false.

        ! call the factory again with the solver settings it wired, which has to stack
        ! a new holder on the first
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, stacked_extension_c, solver_settings_c, &
            c_loc(state_c), settings)
        if (error_c /= 0 .or. .not. c_associated(stacked_extension_c) .or. &
            c_associated(stacked_extension_c, extension_c)) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Solver settings "// &
                "already wired not stacked with a new holder."
            test_s_gek_factory_c_wrapper = .false.
        else
            call c_f_pointer(stacked_extension_c, holder)
            deallocate(holder)
        end if
        if (c_associated(extension_c)) then
            call c_f_pointer(extension_c, holder)
            deallocate(holder)
        end if

        ! call the factory with settings that were not initialized and without any, for
        ! which the factory is handed the default settings with the context through
        ! which the change of reference function is found, which the mock checks, and
        ! the bundle keeps no logging function or host context
        do i_case = 1, size(default_cases)
            solver_settings_c = default_solver_settings
            mock_s_gek_factory_defaults = .true.
            test_passed_mock_s_gek_factory = .false.
            test_passed_mock_change_reference = .false.
            if (i_case == 1) then
                settings%initialized = .false.
                error_c = s_gek_factory_c_wrapper( &
                    n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
                    c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
                    obj_func_s_gek_c_funptr, extension_c, solver_settings_c, &
                    c_loc(state_c), settings)
                settings%initialized = .true.
            else
                error_c = s_gek_factory_c_wrapper( &
                    n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
                    c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
                    obj_func_s_gek_c_funptr, extension_c, solver_settings_c, &
                    c_loc(state_c))
            end if
            mock_s_gek_factory_defaults = .false.
            if (error_c /= 0 .or. .not. test_passed_mock_s_gek_factory .or. &
                .not. test_passed_mock_change_reference) then
                write(stderr, *) "test_s_gek_factory_c_wrapper failed: Default "// &
                    "settings not handed to the factory for "// &
                    trim(default_cases(i_case))//"."
                test_s_gek_factory_c_wrapper = .false.
            else
                call c_f_pointer(extension_c, holder)
                select type (callbacks => holder%extension_callbacks)
                type is (s_gek_c_callbacks_type)
                    if (associated(callbacks%logger) .or. &
                        c_associated(callbacks%host_context)) then
                        write(stderr, *) "test_s_gek_factory_c_wrapper failed: "// &
                            "Logging function or host context kept for "// &
                            trim(default_cases(i_case))//"."
                        test_s_gek_factory_c_wrapper = .false.
                    end if
                end select
            end if
            if (c_associated(extension_c)) then
                call c_f_pointer(extension_c, holder)
                deallocate(holder)
            end if
        end do

        ! call the factory without a record of the state, which has to produce an error
        ! without calling the factory or returning functions or an extension
        solver_settings_c = default_solver_settings
        test_passed_mock_s_gek_factory = .false.
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, extension_c, solver_settings_c, c_null_ptr, &
            settings)
        if (error_c == 0 .or. test_passed_mock_s_gek_factory .or. &
            c_associated(update_orbs_s_gek_c_funptr) .or. &
            c_associated(obj_func_s_gek_c_funptr) .or. c_associated(extension_c)) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: No error "// &
                "produced without a record of the state."
            test_s_gek_factory_c_wrapper = .false.
        end if

        ! call the factory with a mock factory which produces an error, which has to be
        ! returned without functions, an extension or writing back the solver settings
        solver_settings_c = default_solver_settings
        update_orbs_s_gek_c_funptr = c_funloc(mock_update_orbs)
        obj_func_s_gek_c_funptr = c_funloc(mock_obj_func)
        mock_s_gek_factory_error = 1
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, extension_c, solver_settings_c, c_loc(state_c), &
            settings)
        mock_s_gek_factory_error = 0
        if (error_c /= 1 .or. c_associated(update_orbs_s_gek_c_funptr) .or. &
            c_associated(obj_func_s_gek_c_funptr) .or. c_associated(extension_c) .or. &
            c_associated(solver_settings_c%context)) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Error of the "// &
                "factory not passed on."
            test_s_gek_factory_c_wrapper = .false.
        end if

    end function test_s_gek_factory_c_wrapper

    logical(c_bool) function test_change_reference_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the reference change, which finds
        ! the C function in the bundle it receives as context
        !
        use otr_s_gek, only: change_reference_type
        use otr_s_gek_c_interface, only: s_gek_c_callbacks_type, &
                                         change_reference_f_wrapper
        use otr_s_gek_test_reference, only: check_change_reference_funptr
        use test_reference, only: arm_host_context_c, host_context_reached

        type(s_gek_c_callbacks_type), target :: callbacks
        procedure(change_reference_type), pointer :: change_reference_funptr
        class(*), pointer :: context
        real(rp) :: vector(n_param), list(n_param, 1)
        integer(ip) :: error

        ! assume tests pass
        test_change_reference_f_wrapper = .true.

        ! set test name
        test_name = "test_change_reference_f_wrapper"

        ! inject mock C function into a bundle, which the wrapper receives as context
        callbacks%change_reference => mock_change_reference
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        change_reference_funptr => change_reference_f_wrapper

        ! clear the flag of the mock
        test_passed_mock_change_reference = .false.

        ! test change reference Fortran wrapper
        if (.not. check_change_reference_funptr( &
            change_reference_funptr, "change_reference_f_wrapper", "", context)) &
            test_change_reference_f_wrapper = .false.
        if (.not. test_passed_mock_change_reference) &
            test_change_reference_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("change_reference_f_wrapper")) &
            test_change_reference_f_wrapper = .false.

        ! a context that is no bundle has to produce an error
        vector = 1.0_rp
        list = 1.0_rp
        context => null()
        call change_reference_f_wrapper(vector, 1_ip, list, list, list, error, context)
        if (error == 0) then
            write(stderr, *) "test_change_reference_f_wrapper failed: No error "// &
                "produced without a bundle."
            test_change_reference_f_wrapper = .false.
        end if

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

    logical(c_bool) function test_assign_s_gek_f_c() bind(C)
        !
        ! this function tests that the function that converts S-GEK settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_settings_type
        use otr_s_gek_test_reference, only: get_reference_s_gek_values, &
                                            ref_s_gek_settings, operator(/=)
        use c_interface, only: c_callbacks_type
        use c_interface_unit_tests, only: mock_logger
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(s_gek_settings_type_c) :: settings_c
        type(s_gek_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_s_gek_f_c = .true.

        ! initialize the C settings with the reference values
        call get_reference_s_gek_values(settings_c)
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the factory does
        settings = settings_c
        callbacks%logger => mock_logger
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
        ! uninitialized without a logging function, so that the C factory initializes
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
        use otr_s_gek_test_reference, only: ref_s_gek_settings, operator(/=)

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
