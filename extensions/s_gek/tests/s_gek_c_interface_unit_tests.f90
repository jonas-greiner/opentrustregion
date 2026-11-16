! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface_unit_tests

    use opentrustregion_core, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: n_param, n_param_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_associated, c_loc, c_null_ptr, &
                                           c_null_funptr

    implicit none

contains

    function mock_change_reference(new_ref_c, n_points_c, kappa_list_c, &
                                   local_grad_list_c, grad_list_c, context_c) &
        result(error) bind(C)
        !
        ! this function is a mock function for the chart transformation C function,
        ! which changes the points to the new reference in a flat parameter space,
        ! subtracting it from every point and taking the gradients over from the local
        ! gradients, triples the local gradients, as the chart transformation checkers
        ! expect, and returns the error code set in the host context
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

    logical(c_bool) function test_s_gek_solver_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the S-GEK solver
        !
        use opentrustregion_core, only: default_solver_settings, solver_state_type, &
                                        verbosity_error
        use c_interface, only: solver_settings_type_c, solver_state_type_c
        use otr_s_gek, only: standard_s_gek_solver => s_gek_solver, &
                             default_s_gek_settings, error_s_gek
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, s_gek_solver, &
                                         s_gek_solver_c_wrapper
        use otr_s_gek_mock, only: mock_s_gek_solver, test_passed, &
                                  received_solver_settings, received_s_gek_settings, &
                                  received_callbacks, received_s_gek_callbacks, &
                                  solver_state_bundled, s_gek_state_bundled
        use opentrustregion_mock, only: mock_error, mock_solver_state
        use otr_s_gek_test_reference, only: get_reference_s_gek_values, &
                                            ref_s_gek_settings, operator(/=)
        use otr_common_c_interface_unit_tests, only: mock_obj_func_grad
        use c_interface_unit_tests, only: mock_update_orbs, mock_obj_func, &
                                          mock_logger, mock_hess_x
        use test_reference, only: host_context_type, host_context, arm_host_context_c, &
                                  host_context_reached, get_reference_solver_values, &
                                  get_reference_solver_state_values, unset_callbacks, &
                                  ref_solver_settings, operator(/=), operator(==)

        type(c_funptr) :: update_orbs_c_funptr, obj_func_c_funptr, &
                          change_reference_c_funptr
        type(solver_settings_type_c) :: settings
        type(s_gek_settings_type_c) :: s_gek_settings
        type(solver_state_type_c), target :: state
        type(host_context_type), target :: s_gek_host_context
        integer(c_ip) :: error
        integer(ip) :: i_case
        character(len=20), parameter :: required_names(3) = &
            [character(len=20) :: "orbital update", "objective function", &
             "chart transformation"]
        character(len=*), parameter :: default_cases(2) = &
            [character(len=40) :: "no settings", &
             "S-GEK settings that were not initialized"]

        ! assume tests pass
        test_s_gek_solver_c_wrapper = .true.

        ! inject mock function
        s_gek_solver => mock_s_gek_solver

        ! get C function pointers to Fortran functions
        update_orbs_c_funptr = c_funloc(mock_update_orbs)
        obj_func_c_funptr = c_funloc(mock_obj_func)
        change_reference_c_funptr = c_funloc(mock_change_reference)

        ! associate the solver settings with the reference values without callback
        ! functions, whose conversion the test of the conversion routine covers, and
        ! the S-GEK settings with the reference values and every callback function, and
        ! arm the host context of both, which their callback functions receive
        call get_reference_solver_values(settings)
        call unset_callbacks(settings)
        call get_reference_s_gek_values(s_gek_settings)
        s_gek_settings%obj_func_grad = c_funloc(mock_obj_func_grad)
        s_gek_settings%init_hess = c_funloc(mock_hess_x)
        s_gek_settings%logger = c_funloc(mock_logger)
        call arm_host_context_c(s_gek_settings%context)
        settings%context = c_loc(host_context)

        ! clear the result of the mock so that a missing call is detected, and start
        ! the state from the reference values so that a state that is not copied is
        ! detected
        test_passed = .false.
        call get_reference_solver_state_values(state)

        ! call S-GEK solver
        error = s_gek_solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                       change_reference_c_funptr, n_param_c, settings, &
                                       c_loc(state), s_gek_settings)

        ! check if the logging subroutine of the S-GEK settings was called
        if (.not. host_context%logger_called) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Called logging "// &
                "subroutine wrong."
            test_s_gek_solver_c_wrapper = .false.
        end if

        ! check if the settings are correctly passed
        if (received_solver_settings /= ref_solver_settings) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Passed solver "// &
                "settings associated with wrong values."
            test_s_gek_solver_c_wrapper = .false.
        end if
        if (received_s_gek_settings /= ref_s_gek_settings) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Passed S-GEK "// &
                "settings associated with wrong values."
            test_s_gek_solver_c_wrapper = .false.
        end if

        ! check that the bundle of the S-GEK settings holds the chart transformation,
        ! their callback functions, the initial Hessian in the slot of the Hessian
        ! linear transformation, and their host context
        if (.not. associated(received_s_gek_callbacks%change_reference, &
                             mock_change_reference) .or. &
            .not. associated(received_s_gek_callbacks%obj_func_grad, &
                             mock_obj_func_grad) .or. &
            .not. associated(received_s_gek_callbacks%hess_x, mock_hess_x) .or. &
            .not. associated(received_s_gek_callbacks%logger, mock_logger) .or. &
            .not. c_associated(received_s_gek_callbacks%host_context, &
                               c_loc(host_context))) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Callback "// &
                "functions or host context of the S-GEK settings not bundled."
            test_s_gek_solver_c_wrapper = .false.
        end if

        ! check if output variables are as expected
        if (error /= int(mock_error, kind=c_ip)) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Produced error "// &
                "code wrong."
            test_s_gek_solver_c_wrapper = .false.
        end if

        ! check that the state the mock S-GEK solver set is copied into the given one,
        ! also after the error the mock returns
        if (state /= mock_solver_state) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Returned state wrong."
            test_s_gek_solver_c_wrapper = .false.
        end if

        ! check that the callback bundles of the solver settings and of the S-GEK
        ! settings both refer to the state record and the host's record, so that the
        ! state is copied before every callback function of either settings
        if (.not. solver_state_bundled .or. .not. s_gek_state_bundled .or. &
            .not. associated(received_callbacks%solver_state_c, state) .or. &
            .not. associated(received_s_gek_callbacks%solver_state_c, state)) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: State not handed "// &
                "to the callback bundles."
            test_s_gek_solver_c_wrapper = .false.
        end if

        ! check that the host context reached the callback functions unchanged
        if (.not. host_context_reached("s_gek_solver_c_wrapper")) &
            test_s_gek_solver_c_wrapper = .false.

        ! check if test has passed
        if (.not. test_passed) test_s_gek_solver_c_wrapper = .false.

        ! call the S-GEK solver without a state and with a host context of the S-GEK
        ! settings other than that of the solver settings, the callback bundles may not
        ! refer to a state and the bundle of the S-GEK settings has to hold their own
        ! host context, which their logging function receives
        s_gek_settings%context = c_loc(s_gek_host_context)
        test_passed = .false.
        error = s_gek_solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                       change_reference_c_funptr, n_param_c, settings, &
                                       c_null_ptr, s_gek_settings)
        if (solver_state_bundled .or. s_gek_state_bundled .or. &
            associated(received_callbacks%solver_state_c) .or. &
            associated(received_s_gek_callbacks%solver_state_c)) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: State handed to "// &
                "the callback bundles without a host's record."
            test_s_gek_solver_c_wrapper = .false.
        end if
        if (.not. c_associated(received_s_gek_callbacks%host_context, &
                               c_loc(s_gek_host_context)) .or. &
            .not. c_associated(received_callbacks%host_context, &
                               c_loc(host_context)) .or. &
            .not. s_gek_host_context%logger_called) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Host context of "// &
                "the S-GEK settings not kept apart from that of the solver settings."
            test_s_gek_solver_c_wrapper = .false.
        end if
        if (.not. test_passed) test_s_gek_solver_c_wrapper = .false.

        ! call the S-GEK solver without settings and with S-GEK settings that were not
        ! initialized, which still carry callback functions and a host context, the
        ! S-GEK solver has to receive the default settings and callback bundles that
        ! hold the given callback functions but neither callback functions of the
        ! settings nor a host context
        call get_reference_s_gek_values(s_gek_settings)
        s_gek_settings%obj_func_grad = c_funloc(mock_obj_func_grad)
        s_gek_settings%init_hess = c_funloc(mock_hess_x)
        s_gek_settings%logger = c_funloc(mock_logger)
        s_gek_settings%context = c_loc(host_context)
        s_gek_settings%initialized = .false.
        do i_case = 1, size(default_cases)
            test_passed = .false.
            if (i_case == 1) then
                error = s_gek_solver_c_wrapper( &
                    update_orbs_c_funptr, obj_func_c_funptr, &
                    change_reference_c_funptr, n_param_c, state_c_ptr=c_null_ptr)
            else
                error = s_gek_solver_c_wrapper( &
                    update_orbs_c_funptr, obj_func_c_funptr, &
                    change_reference_c_funptr, n_param_c, state_c_ptr=c_null_ptr, &
                    s_gek_settings_c=s_gek_settings)
            end if
            if (received_solver_settings /= default_solver_settings .or. &
                received_s_gek_settings /= default_s_gek_settings) then
                write(stderr, *) "test_s_gek_solver_c_wrapper failed: Default "// &
                    "settings not passed for "//trim(default_cases(i_case))//"."
                test_s_gek_solver_c_wrapper = .false.
            end if
            if (.not. associated(received_s_gek_callbacks%change_reference, &
                                 mock_change_reference)) then
                write(stderr, *) "test_s_gek_solver_c_wrapper failed: Chart "// &
                    "transformation not bundled for "//trim(default_cases(i_case))//"."
                test_s_gek_solver_c_wrapper = .false.
            end if
            if (associated(received_s_gek_settings%obj_func_grad) .or. &
                associated(received_s_gek_settings%init_hess) .or. &
                associated(received_s_gek_settings%logger) .or. &
                associated(received_s_gek_callbacks%obj_func_grad) .or. &
                associated(received_s_gek_callbacks%hess_x) .or. &
                associated(received_s_gek_callbacks%logger)) then
                write(stderr, *) "test_s_gek_solver_c_wrapper failed: Callback "// &
                    "functions of the S-GEK settings passed for "// &
                    trim(default_cases(i_case))//"."
                test_s_gek_solver_c_wrapper = .false.
            end if
            if (c_associated(received_callbacks%host_context) .or. &
                c_associated(received_s_gek_callbacks%host_context)) then
                write(stderr, *) "test_s_gek_solver_c_wrapper failed: Host context "// &
                    "passed for "//trim(default_cases(i_case))//"."
                test_s_gek_solver_c_wrapper = .false.
            end if
            if (.not. test_passed) test_s_gek_solver_c_wrapper = .false.
        end do

        ! call the S-GEK solver without an orbital update, without an objective
        ! function and without a chart transformation, the wrapper has to return an
        ! error with an error message, which is the only message passed to the logging
        ! function, without calling the S-GEK solver
        do i_case = 1, size(required_names)
            call get_reference_solver_values(settings)
            call unset_callbacks(settings)
            settings%logger = c_funloc(mock_logger)
            call arm_host_context_c(settings%context)
            settings%verbose = int(verbosity_error, kind=c_ip)
            test_passed = .false.
            error = s_gek_solver_c_wrapper( &
                merge(c_null_funptr, update_orbs_c_funptr, i_case == 1), &
                merge(c_null_funptr, obj_func_c_funptr, i_case == 2), &
                merge(c_null_funptr, change_reference_c_funptr, i_case == 3), &
                n_param_c, settings, c_null_ptr)
            if (error /= error_s_gek + 1) then
                write(stderr, *) "test_s_gek_solver_c_wrapper failed: No error "// &
                    "returned without "//trim(required_names(i_case))//"."
                test_s_gek_solver_c_wrapper = .false.
            end if
            if (host_context%n_messages == 0) then
                write(stderr, *) "test_s_gek_solver_c_wrapper failed: No error "// &
                    "message printed without "//trim(required_names(i_case))//"."
                test_s_gek_solver_c_wrapper = .false.
            end if
            if (test_passed) then
                write(stderr, *) "test_s_gek_solver_c_wrapper failed: S-GEK solver "// &
                    "called without "//trim(required_names(i_case))//"."
                test_s_gek_solver_c_wrapper = .false.
            end if
        end do

        ! call the S-GEK solver without an orbital update and without settings, whose
        ! logging function would copy the state into the host's record, the wrapper has
        ! to leave the default state there itself
        call get_reference_solver_state_values(state)
        error = s_gek_solver_c_wrapper(c_null_funptr, obj_func_c_funptr, &
                                       change_reference_c_funptr, n_param_c, &
                                       state_c_ptr=c_loc(state))
        if (.not. (state == solver_state_type())) then
            write(stderr, *) "test_s_gek_solver_c_wrapper failed: Default state "// &
                "not left in the host's record without orbital update."
            test_s_gek_solver_c_wrapper = .false.
        end if

        ! restore the procedure pointer so later tests do not inherit the mock
        s_gek_solver => standard_s_gek_solver

    end function test_s_gek_solver_c_wrapper

    logical(c_bool) function test_change_reference_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the chart transformation, which
        ! finds the C function in the bundle it receives as context and copies the
        ! state into the host's record first
        !
        use otr_s_gek, only: change_reference_type
        use otr_s_gek_c_interface, only: s_gek_c_callbacks_type, &
                                         change_reference_f_wrapper
        use otr_s_gek_test_reference, only: check_change_reference_funptr
        use test_reference, only: arm_host_context_c, arm_state_c, host_context, &
                                  host_context_reached

        type(s_gek_c_callbacks_type), target :: callbacks
        procedure(change_reference_type), pointer :: change_reference_funptr
        class(*), pointer :: context
        real(rp) :: new_ref(n_param), kappa_list(n_param, 1), &
                    local_grad_list(n_param, 1), grad_list(n_param, 1)
        integer(ip) :: error

        ! assume tests pass
        test_change_reference_f_wrapper = .true.

        ! inject mock C function into a bundle, which the wrapper receives as context
        ! and which refers to a state and a host's record the state has to be copied
        ! into
        callbacks%change_reference => mock_change_reference
        call arm_host_context_c(callbacks%host_context)
        call arm_state_c(callbacks%c_callbacks_type)
        context => callbacks

        ! get pointer to subroutine
        change_reference_funptr => change_reference_f_wrapper

        ! test chart transformation Fortran wrapper
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
        call change_reference_f_wrapper(new_ref, kappa_list, local_grad_list, &
                                        grad_list, error, context)
        if (error /= 3) then
            write(stderr, *) "test_change_reference_f_wrapper failed: Error of the "// &
                "C function not returned."
            test_change_reference_f_wrapper = .false.
        end if

        ! a context that is no bundle has to produce an error
        context => null()
        call change_reference_f_wrapper(new_ref, kappa_list, local_grad_list, &
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

        ! initialize settings which hold the reference values with callback functions
        ! and a host context, so that a field the initialization leaves alone keeps a
        ! value that differs from its default
        call get_reference_s_gek_values(settings)
        call init_s_gek_settings_c(settings)

        ! check function pointers and host context
        if (c_associated(settings%obj_func_grad) .or. &
            c_associated(settings%init_hess) .or. c_associated(settings%logger) .or. &
            c_associated(settings%context)) then
            write(stderr, *) "test_init_s_gek_settings_c failed: Callback function "// &
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
        ! Fortran correctly performs this conversion
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, &
                                         s_gek_c_callbacks_type, assignment(=)
        use otr_common_c_interface, only: obj_func_grad_f_wrapper
        use otr_common_c_interface_unit_tests, only: mock_obj_func_grad
        use c_interface, only: hess_x_f_wrapper
        use otr_s_gek, only: s_gek_settings_type, default_s_gek_settings
        use otr_s_gek_test_reference, only: get_reference_s_gek_values, &
                                            ref_s_gek_settings, operator(/=)
        use c_interface_unit_tests, only: mock_logger, mock_hess_x
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(s_gek_settings_type_c) :: settings_c
        type(s_gek_settings_type) :: settings
        type(s_gek_c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_s_gek_f_c = .true.

        ! initialize the C settings with the reference values
        call get_reference_s_gek_values(settings_c)
        settings_c%obj_func_grad = c_funloc(mock_obj_func_grad)
        settings_c%init_hess = c_funloc(mock_hess_x)
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the C wrapper of the S-GEK solver does
        settings = settings_c
        callbacks%logger => mock_logger
        settings%context => callbacks

        ! check that the objective function and the initial Hessian are converted to
        ! their Fortran wrappers
        if (.not. associated(settings%obj_func_grad, obj_func_grad_f_wrapper) .or. &
            .not. associated(settings%init_hess, hess_x_f_wrapper)) then
            test_assign_s_gek_f_c = .false.
            write(stderr, *) "test_assign_s_gek_f_c failed: Callback functions not "// &
                "converted to their Fortran wrappers."
        end if

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

        ! convert C settings with callback functions that were not initialized, the
        ! custom values and callback functions they still carry have to be replaced by
        ! the default settings
        call get_reference_s_gek_values(settings_c)
        settings_c%initialized = .false.
        settings = settings_c
        if (associated(settings%obj_func_grad) .or. &
            associated(settings%init_hess) .or. associated(settings%logger) .or. &
            settings /= default_s_gek_settings) then
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

        ! convert Fortran settings with the reference values to C settings which hold
        ! callback functions and a host context that the conversion has to discard
        settings = ref_s_gek_settings
        call get_reference_s_gek_values(settings_c)
        settings_c = settings

        ! check that no callback function pointers and host contexts are associated
        if (c_associated(settings_c%obj_func_grad) .or. &
            c_associated(settings_c%init_hess) .or. c_associated(settings_c%logger)) &
            then
            test_assign_s_gek_c_f = .false.
            write(stderr, *) "test_assign_s_gek_c_f failed: Callback function "// &
                "associated."
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
