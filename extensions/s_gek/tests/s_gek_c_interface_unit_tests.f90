! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: n_param, n_param_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_associated, c_loc, c_f_pointer, c_null_ptr

    implicit none

contains

    function mock_change_reference(new_ref_c, n_points_c, kappa_list_c, &
                                   local_grad_list_c, grad_list_c, context_c) &
        result(error) bind(C)
        !
        ! this subroutine is a mock subroutine for the reference change C function,
        ! which changes the stored points to the new reference in a flat parameter
        ! space, subtracting it from every point and taking the gradients over from the
        ! local gradients, triples the local gradients, as the change of reference
        ! checkers expect, and returns the error code set in the host context
        !
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(in), target :: new_ref_c(*)
        integer(c_ip), intent(in), value :: n_points_c
        real(c_rp), intent(inout), target :: kappa_list_c(*), local_grad_list_c(*), &
                                             grad_list_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: i_point, first, last

        ! check host context
        call check_host_context_c(context_c)

        ! change the reference of every point, whose entries follow one another
        do i_point = 1, n_points_c
            first = (i_point - 1) * n_param + 1
            last = i_point * n_param
            kappa_list_c(first:last) = kappa_list_c(first:last) - new_ref_c(:n_param)
            grad_list_c(first:last) = local_grad_list_c(first:last)
            local_grad_list_c(first:last) = 3.0_c_rp * local_grad_list_c(first:last)
        end do

        error = host_context_error(context_c)

    end function mock_change_reference

    logical(c_bool) function test_s_gek_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the S-GEK factory
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, solver_state_type_c, &
                               assignment(=)
        use otr_common, only: error_s_gek
        use otr_common_c_interface, only: extension_holder_type
        use otr_s_gek, only: standard_s_gek_factory => s_gek_factory
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

        type(c_funptr) :: update_orbs_s_gek_c_funptr, obj_func_s_gek_c_funptr, &
                          first_update_orbs_c_funptr, first_obj_func_c_funptr
        type(solver_settings_type_c) :: solver_settings_c
        type(solver_state_type_c), target :: state_c
        type(s_gek_settings_type_c) :: settings
        type(c_ptr) :: extension_c, stacked_extension_c
        type(extension_holder_type), pointer :: holder
        integer(c_ip) :: error_c
        integer(ip) :: i_case
        character(len=34), parameter :: default_cases(2) = &
            [character(len=34) :: "settings that were not initialized", "no settings"]

        ! assume tests pass
        test_s_gek_factory_c_wrapper = .true.

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

        ! call S-GEK factory C wrapper
        test_passed_mock_s_gek_factory = .false.
        mock_s_gek_factory_error = 0
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, extension_c, solver_settings_c, c_loc(state_c), &
            settings)

        ! check if the mock factory received the correct input, the change of reference
        ! function included, and called the logging function
        if (.not. test_passed_mock_s_gek_factory) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Mock factory "// &
                "received wrong input."
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
        else
            call c_f_pointer(extension_c, holder)

            ! check that the holder holds the bundle of the change of reference function
            select type (callbacks => holder%extension_callbacks)
            type is (s_gek_c_callbacks_type)
                if (.not. associated(callbacks%change_reference, &
                                     mock_change_reference)) then
                    write(stderr, *) "test_s_gek_factory_c_wrapper failed: "// &
                        "Callback slot not set."
                    test_s_gek_factory_c_wrapper = .false.
                end if
            class default
                write(stderr, *) "test_s_gek_factory_c_wrapper failed: Bundle of "// &
                    "the callback functions not held."
                test_s_gek_factory_c_wrapper = .false.
            end select

            ! check that the factory wired the settings of the holder with an object,
            ! which reads the state of the holder, that the holder holds the host's
            ! record of the state, and that the factory returned the C wrappers of the
            ! object and wired the solver settings with its holder
            if (.not. check_c_factory_wiring( &
                "s_gek_factory_c_wrapper", holder, state_c, &
                update_orbs_s_gek_c_funptr, obj_func_s_gek_c_funptr, &
                solver_settings_c)) test_s_gek_factory_c_wrapper = .false.
        end if

        ! call the factory again with the solver settings it wired and the functions it
        ! returned, as a host has to, which has to stack a new holder on the first,
        ! whose functions the mock factory checks through the C wrappers of the first
        ! holder, which forward to the host's functions with the armed host context
        call arm_host_context_c(settings%context)
        test_passed_mock_s_gek_factory = .false.
        first_update_orbs_c_funptr = update_orbs_s_gek_c_funptr
        first_obj_func_c_funptr = obj_func_s_gek_c_funptr
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, first_update_orbs_c_funptr, first_obj_func_c_funptr, &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, stacked_extension_c, solver_settings_c, &
            c_loc(state_c), settings)
        if (.not. test_passed_mock_s_gek_factory) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Mock factory "// &
                "received wrong input for solver settings already wired."
            test_s_gek_factory_c_wrapper = .false.
        end if
        if (.not. host_context_reached("s_gek_factory_c_wrapper")) &
            test_s_gek_factory_c_wrapper = .false.
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
            if (error_c /= 0 .or. .not. test_passed_mock_s_gek_factory) then
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
        ! without calling the factory or returning functions or an extension, which are
        ! set beforehand
        solver_settings_c = default_solver_settings
        update_orbs_s_gek_c_funptr = c_funloc(mock_update_orbs)
        obj_func_s_gek_c_funptr = c_funloc(mock_obj_func)
        extension_c = c_loc(state_c)
        test_passed_mock_s_gek_factory = .false.
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, extension_c, solver_settings_c, c_null_ptr, &
            settings)
        if (error_c /= error_s_gek + 1 .or. test_passed_mock_s_gek_factory .or. &
            c_associated(update_orbs_s_gek_c_funptr) .or. &
            c_associated(obj_func_s_gek_c_funptr) .or. c_associated(extension_c)) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: No error "// &
                "produced without a record of the state."
            test_s_gek_factory_c_wrapper = .false.
        end if

        ! call the factory with a mock factory which produces an error after its checks,
        ! with the armed host context, which has to be returned without functions, an
        ! extension or writing back the solver settings
        call arm_host_context_c(settings%context)
        solver_settings_c = default_solver_settings
        solver_settings_c%context = c_loc(host_context)
        update_orbs_s_gek_c_funptr = c_funloc(mock_update_orbs)
        obj_func_s_gek_c_funptr = c_funloc(mock_obj_func)
        extension_c = c_loc(state_c)
        test_passed_mock_s_gek_factory = .false.
        mock_s_gek_factory_error = 1
        error_c = s_gek_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_change_reference), update_orbs_s_gek_c_funptr, &
            obj_func_s_gek_c_funptr, extension_c, solver_settings_c, c_loc(state_c), &
            settings)
        mock_s_gek_factory_error = 0
        if (.not. test_passed_mock_s_gek_factory) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Mock factory "// &
                "received wrong input for a failing factory."
            test_s_gek_factory_c_wrapper = .false.
        end if
        if (.not. host_context_reached("s_gek_factory_c_wrapper")) &
            test_s_gek_factory_c_wrapper = .false.
        if (error_c /= 1 .or. c_associated(update_orbs_s_gek_c_funptr) .or. &
            c_associated(obj_func_s_gek_c_funptr) .or. c_associated(extension_c) .or. &
            .not. c_associated(solver_settings_c%context, c_loc(host_context))) then
            write(stderr, *) "test_s_gek_factory_c_wrapper failed: Error of the "// &
                "factory not passed on."
            test_s_gek_factory_c_wrapper = .false.
        end if

        ! restore the procedure pointer so later tests do not inherit the mock
        s_gek_factory => standard_s_gek_factory

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
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(s_gek_c_callbacks_type), target :: callbacks
        procedure(change_reference_type), pointer :: change_reference_funptr
        class(*), pointer :: context
        real(rp) :: new_ref(n_param), kappa_list(n_param, 1), &
                    local_grad_list(n_param, 1), grad_list(n_param, 1)
        integer(ip) :: error

        ! assume tests pass
        test_change_reference_f_wrapper = .true.

        ! inject mock C function into a bundle, which the wrapper receives as context
        callbacks%change_reference => mock_change_reference
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        change_reference_funptr => change_reference_f_wrapper

        ! test change reference Fortran wrapper
        if (.not. check_change_reference_funptr( &
            change_reference_funptr, "change_reference_f_wrapper", "", context)) &
            test_change_reference_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("change_reference_f_wrapper")) &
            test_change_reference_f_wrapper = .false.

        ! a C function that fails has to have its error returned
        call arm_host_context_c(callbacks%host_context)
        host_context%mock_error = 3
        call random_number(new_ref)
        call random_number(kappa_list)
        call random_number(local_grad_list)
        call random_number(grad_list)
        call change_reference_f_wrapper(new_ref, 1_ip, kappa_list, local_grad_list, &
                                        grad_list, error, context)
        if (error /= 3) then
            write(stderr, *) "test_change_reference_f_wrapper failed: Error of the "// &
                "C function not returned."
            test_change_reference_f_wrapper = .false.
        end if

        ! a context that is no bundle has to produce an error
        context => null()
        call change_reference_f_wrapper(new_ref, 1_ip, kappa_list, local_grad_list, &
                                        grad_list, error, context)
        if (error /= 1) then
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
        use otr_s_gek_test_reference, only: get_reference_s_gek_values, operator(/=)

        type(s_gek_settings_type_c) :: settings

        ! assume test passes
        test_init_s_gek_settings_c = .true.

        ! initialize settings which hold the reference values with a logging function
        ! and a host context, so that a field the initialization leaves alone keeps a
        ! value that differs from its default
        call get_reference_s_gek_values(settings)
        call init_s_gek_settings_c(settings)

        ! check function pointers and host context
        if (c_associated(settings%logger) .or. c_associated(settings%context)) then
            write(stderr, *) "test_init_s_gek_settings_c failed: Logging function "// &
                "or host context kept."
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
        use otr_s_gek, only: s_gek_settings_type, default_s_gek_settings
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

        ! convert C settings with a logging function that were not initialized, the
        ! custom values and logging function they still carry have to be replaced by the
        ! default settings
        settings_c%initialized = .false.
        settings = settings_c
        if (associated(settings%logger) .or. settings /= default_s_gek_settings) then
            write(stderr, *) "test_assign_s_gek_f_c failed: Settings that were not "// &
                "initialized not converted to default values."
            test_assign_s_gek_f_c = .false.
        end if

    end function test_assign_s_gek_f_c

    logical(c_bool) function test_assign_s_gek_c_f() bind(C)
        !
        ! this function tests that the function that converts S-GEK settings from
        ! Fortran to C correctly performs this conversion
        !
        use otr_s_gek, only: s_gek_settings_type, default_s_gek_settings
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek_test_reference, only: ref_s_gek_settings, &
                                            get_reference_s_gek_values, operator(/=)

        type(s_gek_settings_type) :: settings
        type(s_gek_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_s_gek_c_f = .true.

        ! convert Fortran settings with the reference values to C settings which hold a
        ! logging function and a host context that the conversion has to discard
        settings = ref_s_gek_settings
        call get_reference_s_gek_values(settings_c)
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

        ! convert the default settings into C settings which hold the reference values
        ! and are marked as not initialized, so that a field which is not converted,
        ! including the initialized flag, keeps a value that differs from the expected
        ! one
        call get_reference_s_gek_values(settings_c)
        settings_c%initialized = .false.
        settings = default_s_gek_settings
        settings_c = settings
        if (settings_c /= default_s_gek_settings) then
            write(stderr, *) "test_assign_s_gek_c_f failed: Default settings not "// &
                "converted correctly."
            test_assign_s_gek_c_f = .false.
        end if

    end function test_assign_s_gek_c_f

end module otr_s_gek_c_interface_unit_tests
