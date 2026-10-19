! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_rotation_qn_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: n_param, n_param_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_associated, c_loc, c_f_pointer, c_null_ptr

    implicit none

contains

    function mock_transport(geodesic_c, tangent_vector_c, context_c) result(error) &
        bind(C)
        !
        ! this subroutine is a mock subroutine for the transport C function, which
        ! returns the error code set in the host context
        !
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(in), target :: geodesic_c(*)
        real(c_rp), intent(inout), target :: tangent_vector_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        tangent_vector_c(:n_param) = geodesic_c(:n_param) * tangent_vector_c(:n_param)

        error = host_context_error(context_c)

    end function mock_transport

    function mock_init_hess(vector_c, context_c) result(error) bind(C)
        !
        ! this subroutine is a mock subroutine for the initial Hessian C function, which
        ! returns the error code set in the host context
        !
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(inout), target :: vector_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        vector_c(:n_param) = 2.0_c_rp * vector_c(:n_param)

        error = host_context_error(context_c)

    end function mock_init_hess

    logical(c_bool) function test_rotation_qn_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the rotation QN factory
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, solver_state_type_c, &
                               assignment(=)
        use otr_common, only: error_rotation_qn
        use otr_common_c_interface, only: extension_holder_type
        use otr_rotation_qn, only: standard_rotation_qn_factory => rotation_qn_factory
        use otr_rotation_qn_c_interface, only: &
            rotation_qn_settings_type_c, rotation_qn_c_callbacks_type, &
            rotation_qn_factory, rotation_qn_factory_c_wrapper
        use otr_rotation_qn_mock, only: &
            mock_rotation_qn_factory, mock_rotation_qn_factory_error, &
            mock_rotation_qn_factory_defaults, &
            test_passed_mock_rotation_qn_factory => test_passed
        use otr_rotation_qn_test_reference, only: get_reference_rotation_qn_values
        use otr_common_c_interface_unit_tests, only: check_c_factory_wiring
        use c_interface_unit_tests, only: mock_update_orbs, mock_obj_func, mock_logger
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(c_funptr) :: update_orbs_rotation_qn_c_funptr, &
                          obj_func_rotation_qn_c_funptr, first_update_orbs_c_funptr, &
                          first_obj_func_c_funptr
        type(solver_settings_type_c) :: solver_settings_c
        type(solver_state_type_c), target :: state_c
        type(rotation_qn_settings_type_c) :: settings
        type(c_ptr) :: extension_c, stacked_extension_c
        type(extension_holder_type), pointer :: holder
        integer(c_ip) :: error_c
        integer(ip) :: i_case
        character(len=34), parameter :: default_cases(2) = &
            [character(len=34) :: "settings that were not initialized", "no settings"]

        ! assume tests pass
        test_rotation_qn_factory_c_wrapper = .true.

        ! inject mock function
        rotation_qn_factory => mock_rotation_qn_factory

        ! associate optional settings with the reference values
        call get_reference_rotation_qn_values(settings)
        settings%logger = c_funloc(mock_logger)

        ! arm the host context of the settings, which the logging, transport and
        ! initial Hessian functions receive, and of the solver settings, which the
        ! orbital update and objective function receive
        call arm_host_context_c(settings%context)
        solver_settings_c = default_solver_settings
        solver_settings_c%context = c_loc(host_context)

        ! call rotation QN factory C wrapper
        test_passed_mock_rotation_qn_factory = .false.
        mock_rotation_qn_factory_error = 0
        error_c = rotation_qn_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_transport), c_funloc(mock_init_hess), &
            update_orbs_rotation_qn_c_funptr, obj_func_rotation_qn_c_funptr, &
            extension_c, solver_settings_c, c_loc(state_c), settings)

        ! check if the mock factory received the correct input and called the logging
        ! function
        if (.not. test_passed_mock_rotation_qn_factory) then
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: Mock "// &
                "factory received wrong input."
            test_rotation_qn_factory_c_wrapper = .false.
        end if
        if (.not. host_context%logger_called) then
            test_rotation_qn_factory_c_wrapper = .false.
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: Called "// &
                "logging subroutine wrong."
        end if

        ! check that the host context reached the callback functions unchanged
        if (.not. host_context_reached("rotation_qn_factory_c_wrapper")) &
            test_rotation_qn_factory_c_wrapper = .false.

        ! check if output variables are as expected
        if (error_c /= 0) then
            test_rotation_qn_factory_c_wrapper = .false.
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: Returned "// &
                "error code wrong."
        end if

        ! the returned extension is the holder, whose inspection would crash without
        ! one
        if (.not. c_associated(extension_c)) then
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: No "// &
                "extension returned."
            test_rotation_qn_factory_c_wrapper = .false.
        else
            call c_f_pointer(extension_c, holder)

            ! check that the holder holds the bundle of the transport and initial
            ! Hessian functions
            select type (callbacks => holder%extension_callbacks)
            type is (rotation_qn_c_callbacks_type)
                if (.not. (associated(callbacks%transport, mock_transport) .and. &
                           associated(callbacks%init_hess, mock_init_hess))) then
                    write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: "// &
                        "Callback slots not set."
                    test_rotation_qn_factory_c_wrapper = .false.
                end if
            class default
                write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: "// &
                    "Bundle of the callback functions not held."
                test_rotation_qn_factory_c_wrapper = .false.
            end select

            ! check that the factory wired the settings of the holder with an object,
            ! which reads the state of the holder, that the holder holds the host's
            ! record of the state, and that the factory returned the C wrappers of the
            ! object and wired the solver settings with its holder
            if (.not. check_c_factory_wiring( &
                "rotation_qn_factory_c_wrapper", holder, state_c, &
                update_orbs_rotation_qn_c_funptr, obj_func_rotation_qn_c_funptr, &
                solver_settings_c)) test_rotation_qn_factory_c_wrapper = .false.
        end if

        ! call the factory again with the solver settings it wired and the functions it
        ! returned, as a host has to, which has to stack a new holder on the first,
        ! whose functions the mock factory checks through the C wrappers of the first
        ! holder, which forward to the host's functions with the armed host context
        call arm_host_context_c(settings%context)
        test_passed_mock_rotation_qn_factory = .false.
        first_update_orbs_c_funptr = update_orbs_rotation_qn_c_funptr
        first_obj_func_c_funptr = obj_func_rotation_qn_c_funptr
        error_c = rotation_qn_factory_c_wrapper( &
            n_param_c, first_update_orbs_c_funptr, first_obj_func_c_funptr, &
            c_funloc(mock_transport), c_funloc(mock_init_hess), &
            update_orbs_rotation_qn_c_funptr, obj_func_rotation_qn_c_funptr, &
            stacked_extension_c, solver_settings_c, c_loc(state_c), settings)
        if (.not. test_passed_mock_rotation_qn_factory) then
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: Mock "// &
                "factory received wrong input for solver settings already wired."
            test_rotation_qn_factory_c_wrapper = .false.
        end if
        if (.not. host_context_reached("rotation_qn_factory_c_wrapper")) &
            test_rotation_qn_factory_c_wrapper = .false.
        if (error_c /= 0 .or. .not. c_associated(stacked_extension_c) .or. &
            c_associated(stacked_extension_c, extension_c)) then
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: Solver "// &
                "settings already wired not stacked with a new holder."
            test_rotation_qn_factory_c_wrapper = .false.
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
        ! which the transport and initial Hessian functions are found, which the mock
        ! checks, and the bundle keeps no logging function or host context
        do i_case = 1, size(default_cases)
            solver_settings_c = default_solver_settings
            mock_rotation_qn_factory_defaults = .true.
            test_passed_mock_rotation_qn_factory = .false.
            if (i_case == 1) then
                settings%initialized = .false.
                error_c = rotation_qn_factory_c_wrapper( &
                    n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
                    c_funloc(mock_transport), c_funloc(mock_init_hess), &
                    update_orbs_rotation_qn_c_funptr, obj_func_rotation_qn_c_funptr, &
                    extension_c, solver_settings_c, c_loc(state_c), settings)
                settings%initialized = .true.
            else
                error_c = rotation_qn_factory_c_wrapper( &
                    n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
                    c_funloc(mock_transport), c_funloc(mock_init_hess), &
                    update_orbs_rotation_qn_c_funptr, obj_func_rotation_qn_c_funptr, &
                    extension_c, solver_settings_c, c_loc(state_c))
            end if
            mock_rotation_qn_factory_defaults = .false.
            if (error_c /= 0 .or. .not. test_passed_mock_rotation_qn_factory) then
                write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: "// &
                    "Default settings not handed to the factory for "// &
                    trim(default_cases(i_case))//"."
                test_rotation_qn_factory_c_wrapper = .false.
            else
                call c_f_pointer(extension_c, holder)
                select type (callbacks => holder%extension_callbacks)
                type is (rotation_qn_c_callbacks_type)
                    if (associated(callbacks%logger) .or. &
                        c_associated(callbacks%host_context)) then
                        write(stderr, *) "test_rotation_qn_factory_c_wrapper "// &
                            "failed: Logging function or host context kept for "// &
                            trim(default_cases(i_case))//"."
                        test_rotation_qn_factory_c_wrapper = .false.
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
        update_orbs_rotation_qn_c_funptr = c_funloc(mock_update_orbs)
        obj_func_rotation_qn_c_funptr = c_funloc(mock_obj_func)
        extension_c = c_loc(state_c)
        test_passed_mock_rotation_qn_factory = .false.
        error_c = rotation_qn_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_transport), c_funloc(mock_init_hess), &
            update_orbs_rotation_qn_c_funptr, obj_func_rotation_qn_c_funptr, &
            extension_c, solver_settings_c, c_null_ptr, settings)
        if (error_c /= error_rotation_qn + 1 .or. &
            test_passed_mock_rotation_qn_factory .or. &
            c_associated(update_orbs_rotation_qn_c_funptr) .or. &
            c_associated(obj_func_rotation_qn_c_funptr) .or. &
            c_associated(extension_c)) then
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: No error "// &
                "produced without a record of the state."
            test_rotation_qn_factory_c_wrapper = .false.
        end if

        ! call the factory with a mock factory which produces an error after its checks,
        ! with the armed host context, which has to be returned without functions, an
        ! extension or writing back the solver settings
        call arm_host_context_c(settings%context)
        solver_settings_c = default_solver_settings
        solver_settings_c%context = c_loc(host_context)
        update_orbs_rotation_qn_c_funptr = c_funloc(mock_update_orbs)
        obj_func_rotation_qn_c_funptr = c_funloc(mock_obj_func)
        extension_c = c_loc(state_c)
        test_passed_mock_rotation_qn_factory = .false.
        mock_rotation_qn_factory_error = 1
        error_c = rotation_qn_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_transport), c_funloc(mock_init_hess), &
            update_orbs_rotation_qn_c_funptr, obj_func_rotation_qn_c_funptr, &
            extension_c, solver_settings_c, c_loc(state_c), settings)
        mock_rotation_qn_factory_error = 0
        if (.not. test_passed_mock_rotation_qn_factory) then
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: Mock "// &
                "factory received wrong input for a failing factory."
            test_rotation_qn_factory_c_wrapper = .false.
        end if
        if (.not. host_context_reached("rotation_qn_factory_c_wrapper")) &
            test_rotation_qn_factory_c_wrapper = .false.
        if (error_c /= 1 .or. c_associated(update_orbs_rotation_qn_c_funptr) .or. &
            c_associated(obj_func_rotation_qn_c_funptr) .or. &
            c_associated(extension_c) .or. &
            .not. c_associated(solver_settings_c%context, c_loc(host_context))) then
            write(stderr, *) "test_rotation_qn_factory_c_wrapper failed: Error of "// &
                "the factory not passed on."
            test_rotation_qn_factory_c_wrapper = .false.
        end if

        ! restore the procedure pointer so later tests do not inherit the mock
        rotation_qn_factory => standard_rotation_qn_factory

    end function test_rotation_qn_factory_c_wrapper

    logical(c_bool) function test_transport_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the transport, which finds the C
        ! function in the bundle it receives as context
        !
        use otr_rotation_qn, only: transport_type
        use otr_rotation_qn_c_interface, only: rotation_qn_c_callbacks_type, &
                                               transport_f_wrapper
        use otr_rotation_qn_test_reference, only: check_transport_funptr
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(rotation_qn_c_callbacks_type), target :: callbacks
        procedure(transport_type), pointer :: transport_funptr
        class(*), pointer :: context
        real(rp) :: geodesic(n_param), tangent_vector(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_transport_f_wrapper = .true.

        ! inject mock C function into a bundle, which the wrapper receives as context
        callbacks%transport => mock_transport
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        transport_funptr => transport_f_wrapper

        ! test transport Fortran wrapper
        if (.not. check_transport_funptr(transport_funptr, "transport_f_wrapper", "", &
                                         context)) test_transport_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("transport_f_wrapper")) &
            test_transport_f_wrapper = .false.

        ! a C function that fails has to have its error returned
        call arm_host_context_c(callbacks%host_context)
        host_context%mock_error = 3
        call random_number(geodesic)
        call random_number(tangent_vector)
        call transport_f_wrapper(geodesic, tangent_vector, error, context)
        if (error /= 3) then
            write(stderr, *) "test_transport_f_wrapper failed: Error of the C "// &
                "function not returned."
            test_transport_f_wrapper = .false.
        end if

        ! a context that is no bundle has to produce an error
        context => null()
        call transport_f_wrapper(geodesic, tangent_vector, error, context)
        if (error /= 1) then
            write(stderr, *) "test_transport_f_wrapper failed: No error produced "// &
                "without a bundle."
            test_transport_f_wrapper = .false.
        end if

    end function test_transport_f_wrapper

    logical(c_bool) function test_init_hess_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the initial Hessian, which finds
        ! the C function in the bundle it receives as context
        !
        use otr_rotation_qn, only: init_hess_type
        use otr_rotation_qn_c_interface, only: rotation_qn_c_callbacks_type, &
                                               init_hess_f_wrapper
        use otr_rotation_qn_test_reference, only: check_init_hess_funptr
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(rotation_qn_c_callbacks_type), target :: callbacks
        procedure(init_hess_type), pointer :: init_hess_funptr
        class(*), pointer :: context
        real(rp) :: vector(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_init_hess_f_wrapper = .true.

        ! inject mock C function into a bundle, which the wrapper receives as context
        callbacks%init_hess => mock_init_hess
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        init_hess_funptr => init_hess_f_wrapper

        ! test initial Hessian Fortran wrapper
        if (.not. check_init_hess_funptr(init_hess_funptr, "init_hess_f_wrapper", "", &
                                         context)) test_init_hess_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("init_hess_f_wrapper")) &
            test_init_hess_f_wrapper = .false.

        ! a C function that fails has to have its error returned
        call arm_host_context_c(callbacks%host_context)
        host_context%mock_error = 3
        call random_number(vector)
        call init_hess_f_wrapper(vector, error, context)
        if (error /= 3) then
            write(stderr, *) "test_init_hess_f_wrapper failed: Error of the C "// &
                "function not returned."
            test_init_hess_f_wrapper = .false.
        end if

        ! a context that is no bundle has to produce an error
        context => null()
        call init_hess_f_wrapper(vector, error, context)
        if (error /= 1) then
            write(stderr, *) "test_init_hess_f_wrapper failed: No error produced "// &
                "without a bundle."
            test_init_hess_f_wrapper = .false.
        end if

    end function test_init_hess_f_wrapper

    logical(c_bool) function test_init_rotation_qn_settings_c() bind(C)
        !
        ! this function tests that the rotation QN settings initialization routine
        ! correctly initializes all settings to their default values
        !
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c, &
                                               init_rotation_qn_settings_c
        use otr_rotation_qn, only: default_rotation_qn_settings
        use otr_rotation_qn_test_reference, only: get_reference_rotation_qn_values, &
                                                  operator(/=)

        type(rotation_qn_settings_type_c) :: settings

        ! assume test passes
        test_init_rotation_qn_settings_c = .true.

        ! initialize settings which hold the reference values with a logging function
        ! and a host context, so that a field the initialization leaves alone keeps a
        ! value that differs from its default
        call get_reference_rotation_qn_values(settings)
        call init_rotation_qn_settings_c(settings)

        ! check function pointers and host context
        if (c_associated(settings%logger) .or. c_associated(settings%context)) then
            write(stderr, *) "test_init_rotation_qn_settings_c failed: Logging "// &
                "function or host context kept."
            test_init_rotation_qn_settings_c = .false.
        end if

        ! check settings
        if (settings /= default_rotation_qn_settings) then
            write(stderr, *) "test_init_rotation_qn_settings_c failed: Settings "// &
                "not initialized correctly."
            test_init_rotation_qn_settings_c = .false.
        end if

    end function test_init_rotation_qn_settings_c

    logical(c_bool) function test_assign_rotation_qn_f_c() bind(C)
        !
        ! this function tests that the function that converts rotation QN settings from
        ! C to Fortran correctly perform this conversion
        !
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c, &
                                               assignment(=)
        use otr_rotation_qn, only: rotation_qn_settings_type, &
                                   default_rotation_qn_settings
        use otr_rotation_qn_test_reference, only: get_reference_rotation_qn_values, &
                                                  ref_rotation_qn_settings, operator(/=)
        use c_interface, only: c_callbacks_type
        use c_interface_unit_tests, only: mock_logger
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(rotation_qn_settings_type_c) :: settings_c
        type(rotation_qn_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_rotation_qn_f_c = .true.

        ! initialize the C settings with the reference values
        call get_reference_rotation_qn_values(settings_c)
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the factory does
        settings = settings_c
        callbacks%logger => mock_logger
        settings%context => callbacks

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_rotation_qn_f_c = .false.
            write(stderr, *) "test_assign_rotation_qn_f_c failed: Logging function "// &
                "not associated with value."
        else
            call arm_host_context_c(callbacks%host_context)
            call settings%logger("test", settings%context)
            if (.not. host_context%logger_called) then
                test_assign_rotation_qn_f_c = .false.
                write(stderr, *) "test_assign_rotation_qn_f_c failed: Called "// &
                    "logging subroutine wrong."
            end if
            if (.not. host_context_reached("assign_rotation_qn_f_c")) &
                test_assign_rotation_qn_f_c = .false.
        end if

        ! check against reference values
        if (settings /= ref_rotation_qn_settings) then
            write(stderr, *) "test_assign_rotation_qn_f_c failed: Settings not "// &
                "converted correctly."
            test_assign_rotation_qn_f_c = .false.
        end if

        ! convert C settings with a logging function that were not initialized, the
        ! custom values and logging function they still carry have to be replaced by the
        ! default settings
        settings_c%initialized = .false.
        settings = settings_c
        if (associated(settings%logger) .or. settings /= default_rotation_qn_settings) &
            then
            write(stderr, *) "test_assign_rotation_qn_f_c failed: Settings that "// &
                "were not initialized not converted to default values."
            test_assign_rotation_qn_f_c = .false.
        end if

    end function test_assign_rotation_qn_f_c

    logical(c_bool) function test_assign_rotation_qn_c_f() bind(C)
        !
        ! this function tests that the function that converts rotation QN settings from
        ! Fortran to C correctly performs this conversion
        !
        use otr_rotation_qn, only: rotation_qn_settings_type, &
                                   default_rotation_qn_settings
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c, &
                                               assignment(=)
        use otr_rotation_qn_test_reference, only: &
            ref_rotation_qn_settings, get_reference_rotation_qn_values, operator(/=)

        type(rotation_qn_settings_type) :: settings
        type(rotation_qn_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_rotation_qn_c_f = .true.

        ! convert Fortran settings with the reference values to C settings which hold a
        ! logging function and a host context that the conversion has to discard
        settings = ref_rotation_qn_settings
        call get_reference_rotation_qn_values(settings_c)
        settings_c = settings

        ! check that no callback function pointers and host contexts are associated
        if (c_associated(settings_c%logger)) then
            test_assign_rotation_qn_c_f = .false.
            write(stderr, *) "test_assign_rotation_qn_c_f failed: Logger function "// &
                "associated."
        end if
        if (c_associated(settings_c%context)) then
            write(stderr, *) "test_assign_rotation_qn_c_f failed: Host context "// &
                "associated."
            test_assign_rotation_qn_c_f = .false.
        end if

        ! check against reference values
        if (settings_c /= ref_rotation_qn_settings) then
            write(stderr, *) "test_assign_rotation_qn_c_f failed: Settings not "// &
                "converted correctly."
            test_assign_rotation_qn_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_rotation_qn_c_f = .false.
            write(stderr, *) "test_assign_rotation_qn_c_f failed: Settings not "// &
                "marked as initialized."
        end if

        ! convert the default settings into C settings which hold the reference values
        ! and are marked as not initialized, so that a field which is not converted,
        ! including the initialized flag, keeps a value that differs from the expected
        ! one
        call get_reference_rotation_qn_values(settings_c)
        settings_c%initialized = .false.
        settings = default_rotation_qn_settings
        settings_c = settings
        if (settings_c /= default_rotation_qn_settings) then
            write(stderr, *) "test_assign_rotation_qn_c_f failed: Default settings "// &
                "not converted correctly."
            test_assign_rotation_qn_c_f = .false.
        end if

    end function test_assign_rotation_qn_c_f

end module otr_rotation_qn_c_interface_unit_tests
