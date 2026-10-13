! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_c_interface_unit_tests

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: n_param, n_param_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_associated, c_loc

    implicit none

contains

    function mock_transport(geodesic_c, tangent_vector_c, context_c) result(error) &
        bind(C)
        !
        ! this subroutine is a mock subroutine for the transport C function
        !
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: geodesic_c(*)
        real(c_rp), intent(inout), target :: tangent_vector_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        tangent_vector_c(:n_param) = geodesic_c(:n_param) * tangent_vector_c(:n_param)

        error = 0_c_ip

    end function mock_transport

    function mock_init_hess(vector_c, context_c) result(error) bind(C)
        !
        ! this subroutine is a mock subroutine for the initial Hessian C function
        !
        use test_reference, only: check_host_context_c

        real(c_rp), intent(inout), target :: vector_c(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        vector_c(:n_param) = 2.0_c_rp * vector_c(:n_param)

        error = 0_c_ip

    end function mock_init_hess

    logical(c_bool) function test_qn_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the quasi-Newton orbital update factory
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, assignment(=)
        use otr_qn, only: update_object
        use otr_qn_c_interface, only: qn_settings_type_c, qn_factory, &
                                      qn_factory_c_wrapper, qn_callbacks, qn_holder
        use otr_qn_mock, only: mock_qn_factory, mock_qn_factory_error, &
                               test_passed_mock_qn_factory => test_passed
        use otr_qn_test_reference, only: assignment(=), ref_qn_settings
        use otr_common_c_interface_unit_tests, only: check_c_factory_wiring
        use c_interface_unit_tests, only: mock_update_orbs, mock_obj_func
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: host_context_type, arm_host_context_c, host_context, &
                                  host_context_reached

        type(c_funptr) :: update_orbs_qn_c_funptr, obj_func_qn_c_funptr
        type(solver_settings_type_c) :: solver_settings_c
        type(qn_settings_type_c) :: settings
        integer(c_ip) :: error_c
        type(host_context_type), target :: other_host_context

        ! assume tests pass
        test_qn_factory_c_wrapper = .true.

        ! inject mock function
        qn_factory => mock_qn_factory

        ! associate optional settings with values
        settings = ref_qn_settings
        settings%logger = c_funloc(mock_logger_c)

        ! arm the host context of the settings, which the logging, transport and
        ! initial Hessian functions receive, and of the solver settings, which the
        ! orbital update and objective function receive
        call arm_host_context_c(settings%context)
        solver_settings_c = default_solver_settings
        solver_settings_c%context = c_loc(host_context)

        ! call quasi-Newton orbital updating factory C wrapper
        test_passed_mock_qn_factory = .false.
        mock_qn_factory_error = 0
        error_c = qn_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_transport), c_funloc(mock_init_hess), &
            update_orbs_qn_c_funptr, obj_func_qn_c_funptr, solver_settings_c, settings)

        ! check if the mock factory received the correct input and called the logging
        ! function
        if (.not. test_passed_mock_qn_factory) then
            write(stderr, *) "test_qn_factory_c_wrapper failed: Mock factory "// &
                "received wrong input."
            test_qn_factory_c_wrapper = .false.
        end if
        if (.not. host_context%logger_called) then
            test_qn_factory_c_wrapper = .false.
            write(stderr, *) "test_qn_factory_c_wrapper failed: Called logging "// &
                "subroutine wrong."
        end if

        ! check that the host context reached the callback functions unchanged
        if (.not. host_context_reached("qn_factory_c_wrapper")) &
            test_qn_factory_c_wrapper = .false.

        ! check if output variables are as expected
        if (error_c /= 0) then
            test_qn_factory_c_wrapper = .false.
            write(stderr, *) "test_qn_factory_c_wrapper failed: Returned error "// &
                "code wrong."
        end if

        ! check if the callback slots point to the given functions
        if (.not. (associated(qn_callbacks%transport, mock_transport) .and. &
                   associated(qn_callbacks%init_hess, mock_init_hess))) then
            write(stderr, *) "test_qn_factory_c_wrapper failed: Callback slots not set."
            test_qn_factory_c_wrapper = .false.
        end if

        ! check that the factory returned the C wrappers of the object and wired the
        ! solver settings with its holder
        if (.not. check_c_factory_wiring("qn_factory_c_wrapper", qn_holder, &
                                         update_object, update_orbs_qn_c_funptr, &
                                         obj_func_qn_c_funptr, solver_settings_c)) &
            test_qn_factory_c_wrapper = .false.

        ! call the factory with the solver settings it wired and another host context,
        ! which has to produce an error without calling the factory, returning
        ! functions or replacing the bundle holding the callback functions of the
        ! object
        settings%context = c_loc(other_host_context)
        error_c = qn_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_transport), c_funloc(mock_init_hess), &
            update_orbs_qn_c_funptr, obj_func_qn_c_funptr, solver_settings_c, settings)
        if (error_c == 0 .or. c_associated(update_orbs_qn_c_funptr) .or. &
            c_associated(obj_func_qn_c_funptr) .or. &
            .not. c_associated(qn_callbacks%host_context, c_loc(host_context))) then
            write(stderr, *) "test_qn_factory_c_wrapper failed: Solver settings "// &
                "the holder was wired into not refused before the object was changed."
            test_qn_factory_c_wrapper = .false.
        end if

        ! call the factory with a mock factory which produces an error, which has to be
        ! returned without functions or writing back the solver settings or changing
        ! the holder or the bundle of the object
        solver_settings_c = default_solver_settings
        update_orbs_qn_c_funptr = c_funloc(mock_logger_c)
        obj_func_qn_c_funptr = c_funloc(mock_logger_c)
        mock_qn_factory_error = 1
        error_c = qn_factory_c_wrapper( &
            n_param_c, c_funloc(mock_update_orbs), c_funloc(mock_obj_func), &
            c_funloc(mock_transport), c_funloc(mock_init_hess), &
            update_orbs_qn_c_funptr, obj_func_qn_c_funptr, solver_settings_c, settings)
        mock_qn_factory_error = 0
        if (error_c /= 1 .or. c_associated(update_orbs_qn_c_funptr) .or. &
            c_associated(obj_func_qn_c_funptr) .or. &
            c_associated(solver_settings_c%context)) then
            write(stderr, *) "test_qn_factory_c_wrapper failed: Error of the "// &
                "factory not passed on."
            test_qn_factory_c_wrapper = .false.
        end if
        if (.not. associated(qn_holder%solver_level%level, &
                             update_object%solver_level) .or. &
            .not. c_associated(qn_callbacks%host_context, c_loc(host_context))) then
            write(stderr, *) "test_qn_factory_c_wrapper failed: Holder or callback "// &
                "functions of the object not restored after an error of the factory."
            test_qn_factory_c_wrapper = .false.
        end if

    end function test_qn_factory_c_wrapper

    logical(c_bool) function test_transport_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the transport
        !
        use otr_qn, only: transport_type
        use otr_qn_c_interface, only: qn_callbacks, transport_f_wrapper
        use otr_qn_test_reference, only: test_transport_funptr
        use test_reference, only: arm_host_context_c, host_context_reached

        procedure(transport_type), pointer :: transport_funptr
        class(*), pointer :: context

        ! assume tests pass
        test_transport_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        qn_callbacks%transport => mock_transport
        call arm_host_context_c(qn_callbacks%host_context)
        context => qn_callbacks

        ! get pointer to subroutine
        transport_funptr => transport_f_wrapper

        ! test transport Fortran wrapper
        if (.not. test_transport_funptr(transport_funptr, "transport_f_wrapper", "", &
                                        context)) test_transport_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("transport_f_wrapper")) &
            test_transport_f_wrapper = .false.

    end function test_transport_f_wrapper

    logical(c_bool) function test_init_hess_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the initial Hessian
        !
        use otr_qn, only: init_hess_type
        use otr_qn_c_interface, only: qn_callbacks, init_hess_f_wrapper
        use otr_qn_test_reference, only: test_init_hess_funptr
        use test_reference, only: arm_host_context_c, host_context_reached

        procedure(init_hess_type), pointer :: init_hess_funptr
        class(*), pointer :: context

        ! assume tests pass
        test_init_hess_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        qn_callbacks%init_hess => mock_init_hess
        call arm_host_context_c(qn_callbacks%host_context)
        context => qn_callbacks

        ! get pointer to subroutine
        init_hess_funptr => init_hess_f_wrapper

        ! test initial Hessian Fortran wrapper
        if (.not. test_init_hess_funptr(init_hess_funptr, "init_hess_f_wrapper", "", &
                                        context)) test_init_hess_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("init_hess_f_wrapper")) &
            test_init_hess_f_wrapper = .false.

    end function test_init_hess_f_wrapper

    logical(c_bool) function test_init_qn_settings_c() bind(C)
        !
        ! this function tests that the quasi-Newton settings initialization routine
        ! correctly initializes all settings to their default values
        !
        use otr_qn_c_interface, only: qn_settings_type_c, init_qn_settings_c
        use otr_qn, only: default_qn_settings
        use otr_qn_test_reference, only: operator(/=)

        type(qn_settings_type_c) :: settings

        ! assume test passes
        test_init_qn_settings_c = .true.

        ! initialize settings
        call init_qn_settings_c(settings)

        ! check function pointers
        if (c_associated(settings%logger)) then
            write(stderr, *) "test_init_qn_settings_c failed: Function pointers "// &
                "should not be initialized."
            test_init_qn_settings_c = .false.
        end if

        ! check settings
        if (settings /= default_qn_settings) then
            write(stderr, *) "test_init_qn_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_qn_settings_c = .false.
        end if

    end function test_init_qn_settings_c

    logical(c_bool) function test_update_orbs_qn_deconstructor_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the quasi-Newton deconstructor
        !
        use otr_qn_c_interface, only: update_orbs_qn_deconstructor, &
                                      update_orbs_qn_deconstructor_c_wrapper
        use otr_qn_mock, only: mock_update_orbs_qn_deconstructor, test_passed

        ! assume tests pass
        test_update_orbs_qn_deconstructor_c_wrapper = .true.

        ! inject mock function
        update_orbs_qn_deconstructor => mock_update_orbs_qn_deconstructor

        ! initialize test logical
        test_passed = .false.

        ! call quasi-Newton orbital updating deconstructor C wrapper
        call update_orbs_qn_deconstructor_c_wrapper()

        ! check if test has passed
        if (.not. test_passed) then
            test_update_orbs_qn_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_update_orbs_qn_deconstructor_c_wrapper failed: "// &
                "Deconstructor called wrong."
        end if

    end function test_update_orbs_qn_deconstructor_c_wrapper

    logical(c_bool) function test_assign_qn_f_c() bind(C)
        !
        ! this function tests that the function that converts quasi-Newton settings
        ! from C to Fortran correctly perform this conversion
        !
        use otr_qn_c_interface, only: qn_settings_type_c, assignment(=)
        use otr_qn, only: qn_settings_type
        use otr_qn_test_reference, only: assignment(=), ref_qn_settings, operator(/=)
        use c_interface, only: c_callbacks_type
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(qn_settings_type_c) :: settings_c
        type(qn_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_qn_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_qn_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the factory does
        settings = settings_c
        callbacks%logger => mock_logger_c
        settings%context => callbacks

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_qn_f_c = .false.
            write(stderr, *) "test_assign_qn_f_c failed: Logging function not "// &
                "associated with value."
        else
            call arm_host_context_c(callbacks%host_context)
            call settings%logger("test", settings%context)
            if (.not. host_context%logger_called) then
                test_assign_qn_f_c = .false.
                write(stderr, *) "test_assign_qn_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
            if (.not. host_context_reached("assign_qn_f_c")) &
                test_assign_qn_f_c = .false.
        end if

        ! check against reference values
        if (settings /= ref_qn_settings) then
            write(stderr, *) "test_assign_qn_f_c failed: Settings not converted "// &
                "correctly."
            test_assign_qn_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            write(stderr, *) "test_assign_qn_f_c failed: Settings not marked as "// &
                "initialized."
            test_assign_qn_f_c = .false.
        end if

        ! convert C settings that were not initialized and check that they are left
        ! uninitialized without a logging function, so that the factory initializes
        ! them after its warning
        settings_c%initialized = .false.
        settings = settings_c
        if (settings%initialized .or. associated(settings%logger)) then
            write(stderr, *) "test_assign_qn_f_c failed: Settings that were not "// &
                "initialized converted."
            test_assign_qn_f_c = .false.
        end if

    end function test_assign_qn_f_c

    logical(c_bool) function test_assign_qn_c_f() bind(C)
        !
        ! this function tests that the function that converts quasi-Newton settings
        ! from Fortran to C correctly performs this conversion
        !
        use otr_qn, only: qn_settings_type
        use otr_qn_c_interface, only: qn_settings_type_c, assignment(=)
        use otr_qn_test_reference, only: ref_qn_settings, assignment(=), operator(/=)

        type(qn_settings_type) :: settings
        type(qn_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_qn_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_qn_settings

        ! convert to C settings
        settings_c = settings

        ! check that no callback function pointers and host contexts are associated
        if (c_associated(settings_c%logger)) then
            test_assign_qn_c_f = .false.
            write(stderr, *) "test_assign_qn_c_f failed: Logger function associated."
        end if
        if (c_associated(settings_c%context)) then
            write(stderr, *) "test_assign_qn_c_f failed: Host context associated."
            test_assign_qn_c_f = .false.
        end if

        ! check against reference values
        if (settings_c /= ref_qn_settings) then
            write(stderr, *) "test_assign_qn_c_f failed: Settings not converted "// &
                "correctly."
            test_assign_qn_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_qn_c_f = .false.
            write(stderr, *) "test_assign_qn_c_f failed: Settings not marked as "// &
                "initialized."
        end if

    end function test_assign_qn_c_f

end module otr_qn_c_interface_unit_tests
