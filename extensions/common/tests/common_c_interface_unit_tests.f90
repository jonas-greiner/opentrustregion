! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface_unit_tests

    use opentrustregion, only: ip, stderr
    use c_interface, only: c_rp, c_ip
    use otr_common_c_interface, only: evaluate_dm_c_type, get_response_c_type
    use test_reference, only: tol_c, n_param
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, c_loc, &
                                           c_associated, c_null_ptr, c_null_char

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(evaluate_dm_c_type), pointer :: &
        mock_evaluate_dm_cs_ptr => mock_evaluate_dm_cs, &
        mock_evaluate_dm_os_ptr => mock_evaluate_dm_os
    procedure(get_response_c_type), pointer :: &
        mock_get_response_cs_ptr => mock_get_response_cs, &
        mock_get_response_os_ptr => mock_get_response_os

contains

    function mock_evaluate_dm_cs(dm_ao, energy, fock, get_response_c_funptr, &
                                 context_c) result(error) bind(C)
        !
        ! this function is a test function for the density matrix evaluating C function
        ! for 2D density matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao
        use otr_common_unit_tests, only: record_mock_call
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*)
        type(c_funptr), intent(out), optional :: get_response_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(ip) :: flat_len = n_ao**2

        ! check host context
        call check_host_context_c(context_c)

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(get_response_c_funptr)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) fock(:flat_len) = evaluate_dm_factors(1) * dm_ao(:flat_len)
        if (present(get_response_c_funptr)) &
            get_response_c_funptr = c_funloc(mock_get_response_cs)

        error = 0_c_ip

    end function mock_evaluate_dm_cs

    function mock_evaluate_dm_os(dm_ao, energy, fock, get_response_c_funptr, &
                                 context_c) result(error) bind(C)
        !
        ! this function is a test function for the density matrix evaluating C function
        ! for 3D density matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao, n_particle
        use otr_common_unit_tests, only: record_mock_call
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*)
        type(c_funptr), intent(out), optional :: get_response_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(ip) :: flat_len = n_ao**2 * n_particle

        ! check host context
        call check_host_context_c(context_c)

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(get_response_c_funptr)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) fock(:flat_len) = evaluate_dm_factors(1) * dm_ao(:flat_len)
        if (present(get_response_c_funptr)) &
            get_response_c_funptr = c_funloc(mock_get_response_os)

        error = 0_c_ip

    end function mock_evaluate_dm_os

    function mock_get_response_cs(dm_ao, response, context_c) result(error) bind(C)
        !
        ! this function is a test function for the response C function for 2D density
        ! matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out), target :: response(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2

        ! check host context
        call check_host_context_c(context_c)

        response(:flat_len) = evaluate_dm_factors(2) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_get_response_cs

    function mock_get_response_os(dm_ao, response, context_c) result(error) bind(C)
        !
        ! this function is a test function for the response C function for 3D density
        ! matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao, n_particle
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out), target :: response(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2 * n_particle

        ! check host context
        call check_host_context_c(context_c)

        response(:flat_len) = evaluate_dm_factors(2) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_get_response_os

    function holder_level_context(holder, i_level) result(context_c)
        !
        ! this function returns the address of the level of the solver (1) or of the
        ! stability check (2) of a holder as the context a C wrapper receives
        !
        use otr_common_c_interface, only: extension_holder_type

        type(extension_holder_type), intent(in), target :: holder
        integer(ip), intent(in) :: i_level
        type(c_ptr) :: context_c

        if (i_level == 1) then
            context_c = c_loc(holder%solver_level)
        else
            context_c = c_loc(holder%stability_level)
        end if

    end function holder_level_context

    subroutine arm_holder_level(holder, i_level)
        !
        ! this subroutine sets up a holder whose level of the solver (1) or of the
        ! stability check (2) holds the host context, which it arms, in its bundle,
        ! while the bundle of the other level holds none, so that a C wrapper
        ! forwarding with the wrong level reaches the host's C function without the
        ! host context
        !
        use otr_common_c_interface, only: extension_holder_type
        use test_reference, only: arm_host_context_c

        type(extension_holder_type), intent(out), target :: holder
        integer(ip), intent(in) :: i_level

        if (i_level == 1) then
            call arm_host_context_c(holder%solver_level%callbacks%host_context)
        else
            call arm_host_context_c(holder%stability_level%callbacks%host_context)
        end if

    end subroutine arm_holder_level

    logical(c_bool) function test_holder_level() bind(C)
        !
        ! this function tests the function which returns the level of a holder at the
        ! address a C wrapper receives as context
        !
        use otr_common_c_interface, only: extension_holder_type, holder_level

        type(extension_holder_type), target :: holder

        ! assume tests pass
        test_holder_level = .true.

        ! the level of the holder has to be returned for its address
        if (.not. associated(holder_level(c_loc(holder%stability_level)), &
                             holder%stability_level)) then
            write(stderr, *) "test_holder_level failed: Level of the holder not "// &
                "returned."
            test_holder_level = .false.
        end if

        ! a null address has to return none
        if (associated(holder_level(c_null_ptr))) then
            write(stderr, *) "test_holder_level failed: Level returned for a null "// &
                "address."
            test_holder_level = .false.
        end if

    end function test_holder_level

    logical(c_bool) function test_setup_extension_holder() bind(C)
        !
        ! this function tests the subroutine which sets up the holder of an extension
        ! object for a C factory from the host's callback functions and C solver
        ! settings
        !
        use opentrustregion, only: solver_settings_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              setup_error_logging
        use c_interface, only: solver_settings_type_c
        use c_interface_unit_tests, only: &
            mock_update_orbs, mock_obj_func, mock_stability_precond, &
            mock_stability_project, mock_stability_logger, &
            mock_stability_get_extra_trial_vectors, mock_approx_hess_x, &
            mock_conv_check_stability
        use otr_common, only: orbital_settings_type
        use otr_common_c_interface, only: extension_holder_type, setup_extension_holder
        use otr_common_mock, only: mock_logger_c
        use otr_common_unit_tests, only: mock_extension_type
        use test_reference, only: host_context_type, get_reference_solver_values, &
                                  host_context, unset_callbacks, ref_solver_settings, &
                                  operator(/=)

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        type(solver_settings_type_c) :: settings_c
        type(solver_settings_type) :: solver_settings
        type(orbital_settings_type) :: settings
        type(host_context_type), target :: nested_host_context
        type(test_context_type), target :: log_context
        integer(ip) :: error, i_case, i_own
        character(len=55), parameter :: case_names(10) = &
            [character(len=55) :: "without nested callback functions or context", &
             "with a nested context", "with a nested preconditioner", &
             "with a nested projection", "with a nested logging function", "with a "// &
             "nested extra trial vector function", "with a nested approximate "// &
             "Hessian linear transformation", "with a nested convergence check", &
             "with uninitialized nested settings", "with uninitialized settings"]
        logical, parameter :: nested_own(10) = &
            [.false., (.true., i_own=1, 7), .false., .false.]

        ! assume tests pass
        test_setup_extension_holder = .true.

        ! set up the settings of the extension
        call setup_settings(settings)

        ! set up the holder with C solver settings whose nested settings provide
        ! nothing or a single context or callback function of their own, and with a
        ! context of nested settings or settings that were not initialized, of which
        ! only nested settings that were initialized and provide something of their own
        ! may keep the bundle of the holder as context, the holder has to be set up
        ! anew and its bundles have to hold the given callback functions and be handed
        ! to the converted settings
        do i_case = 1, size(case_names)
            call get_reference_solver_values(settings_c)
            call unset_callbacks(settings_c%stability_settings)
            settings_c%stability_settings%context = c_null_ptr
            settings_c%logger = c_funloc(mock_logger_c)
            settings_c%context = c_loc(host_context)
            associate (nested_c => settings_c%stability_settings)
                select case (i_case)
                case (2, 9, 10)
                    nested_c%context = c_loc(nested_host_context)
                case (3)
                    nested_c%precond = c_funloc(mock_stability_precond)
                case (4)
                    nested_c%project = c_funloc(mock_stability_project)
                case (5)
                    nested_c%logger = c_funloc(mock_stability_logger)
                case (6)
                    nested_c%get_extra_trial_vectors = &
                        c_funloc(mock_stability_get_extra_trial_vectors)
                case (7)
                    nested_c%approx_hess_x = c_funloc(mock_approx_hess_x)
                case (8)
                    nested_c%conv_check = c_funloc(mock_conv_check_stability)
                end select
                if (i_case == 9) nested_c%initialized = .false.
            end associate
            if (i_case == 10) settings_c%initialized = .false.
            holder%solver_level%level => extension%solver_level

            call setup_extension_holder(holder, c_funloc(mock_update_orbs), &
                                        c_funloc(mock_obj_func), settings_c, &
                                        solver_settings, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_setup_extension_holder failed: Produced "// &
                    "error "//trim(case_names(i_case))//"."
                test_setup_extension_holder = .false.
            end if
            if (associated(holder%solver_level%level)) then
                write(stderr, *) "test_setup_extension_holder failed: Holder not "// &
                    "set up anew "//trim(case_names(i_case))//"."
                test_setup_extension_holder = .false.
            end if
            if (.not. associated(holder%solver_level%callbacks%update_orbs, &
                                 mock_update_orbs) .or. &
                .not. associated(holder%solver_level%callbacks%obj_func, &
                                 mock_obj_func)) then
                write(stderr, *) "test_setup_extension_holder failed: Host's "// &
                    "callback functions not bundled "//trim(case_names(i_case))//"."
                test_setup_extension_holder = .false.
            end if
            if (.not. associated(solver_settings%context, &
                                 holder%solver_level%callbacks) .or. &
                .not. associated(holder%solver_level%callbacks%stability, &
                                 holder%stability_level%callbacks)) then
                write(stderr, *) "test_setup_extension_holder failed: Settings not "// &
                    "converted with the bundles of the holder "// &
                    trim(case_names(i_case))//"."
                test_setup_extension_holder = .false.
            end if
            if (i_case == 1 .and. solver_settings /= ref_solver_settings) then
                write(stderr, *) "test_setup_extension_holder failed: Given "// &
                    "settings not converted."
                test_setup_extension_holder = .false.
            end if
            if (associated(solver_settings%stability_settings%context, &
                           holder%stability_level%callbacks) .neqv. &
                nested_own(i_case)) then
                write(stderr, *) "test_setup_extension_holder failed: Context of "// &
                    "the nested settings wrong "//trim(case_names(i_case))//"."
                test_setup_extension_holder = .false.
            end if
        end do

        ! set up the holder with C solver settings or nested settings that the holder
        ! was already wired into, which has to produce an error and leave the holder
        ! untouched
        call setup_error_logging(settings, log_context)
        do i_case = 1, 2
            call get_reference_solver_values(settings_c)
            if (i_case == 1) then
                settings_c%context = c_loc(holder%solver_level)
            else
                settings_c%stability_settings%context = c_loc(holder%stability_level)
            end if
            holder%solver_level%level => extension%solver_level
            log_context%log_message = ""
            call setup_extension_holder(holder, c_funloc(mock_update_orbs), &
                                        c_funloc(mock_obj_func), settings_c, &
                                        solver_settings, settings, error)
            if (error == 0) then
                write(stderr, *) "test_setup_extension_holder failed: No error "// &
                    "returned for settings the holder was wired into."
                test_setup_extension_holder = .false.
            end if
            if (log_context%log_message == "") then
                write(stderr, *) "test_setup_extension_holder failed: No error "// &
                    "message printed for settings the holder was wired into."
                test_setup_extension_holder = .false.
            end if
            if (.not. associated(holder%solver_level%level, extension%solver_level)) &
                then
                write(stderr, *) "test_setup_extension_holder failed: Holder "// &
                    "changed for settings the holder was wired into."
                test_setup_extension_holder = .false.
            end if
        end do

    end function test_setup_extension_holder

    logical(c_bool) function test_wire_solver_settings_c() bind(C)
        !
        ! this function tests the subroutine which writes the solver settings a Fortran
        ! factory wired back into the host's C solver settings
        !
        use opentrustregion, only: solver_settings_type
        use c_interface, only: solver_settings_type_c
        use otr_common, only: precond_callback, project_callback, logger_callback, &
                              conv_check_stability_callback
        use otr_common_c_interface, only: &
            extension_holder_type, wire_solver_settings_c, precond_c_wrapper, &
            project_c_wrapper, logger_c_wrapper, conv_check_stability_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type
        use test_reference, only: get_reference_solver_values, unset_callbacks, &
                                  same_callbacks, ref_solver_settings, operator(/=)

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        type(solver_settings_type) :: solver_settings
        type(solver_settings_type_c) :: settings_c, expected_c
        integer(ip) :: i_case
        character(len=40), parameter :: case_names(2) = &
            [character(len=40) :: "with wired nested settings", &
             "with nested settings that were not wired"]

        ! assume tests pass
        test_wire_solver_settings_c = .true.

        ! write back solver settings with the reference values into which a factory
        ! wired a preconditioner, a projection and a logging function, and whose nested
        ! settings it wired with a preconditioner and a convergence check or left
        ! untouched, into C settings that were not initialized and hold other callback
        ! functions and contexts, which have to be replaced by the C wrappers of the
        ! wired roles and the addresses of the levels of the holder, which has to hold
        ! the levels of the extension object
        do i_case = 1, size(case_names)
            solver_settings = ref_solver_settings
            solver_settings%context => extension%solver_level
            solver_settings%precond => precond_callback
            solver_settings%project => project_callback
            solver_settings%logger => logger_callback
            call unset_callbacks(expected_c)
            expected_c%precond = c_funloc(precond_c_wrapper)
            expected_c%project = c_funloc(project_c_wrapper)
            expected_c%logger = c_funloc(logger_c_wrapper)
            if (i_case == 1) then
                solver_settings%stability_settings%context => extension%stability_level
                solver_settings%stability_settings%precond => precond_callback
                solver_settings%stability_settings%conv_check => &
                    conv_check_stability_callback
                expected_c%stability_settings%precond = c_funloc(precond_c_wrapper)
                expected_c%stability_settings%conv_check = &
                    c_funloc(conv_check_stability_c_wrapper)
            end if
            call get_reference_solver_values(settings_c)
            settings_c%initialized = .false.

            call wire_solver_settings_c(holder, extension, solver_settings, settings_c)
            if (.not. associated(holder%solver_level%level, &
                                 extension%solver_level) .or. &
                .not. associated(holder%stability_level%level, &
                                 extension%stability_level)) then
                write(stderr, *) "test_wire_solver_settings_c failed: Levels of "// &
                    "the extension object not held "//trim(case_names(i_case))//"."
                test_wire_solver_settings_c = .false.
            end if
            if (settings_c /= ref_solver_settings) then
                write(stderr, *) "test_wire_solver_settings_c failed: Settings not "// &
                    "converted "//trim(case_names(i_case))//"."
                test_wire_solver_settings_c = .false.
            end if
            if (.not. same_callbacks(settings_c, expected_c) .or. .not. c_associated( &
                settings_c%context, c_loc(holder%solver_level)) .or. &
                (c_associated(settings_c%stability_settings%context) .neqv. &
                 i_case == 1) .or. (i_case == 1 .and. .not. c_associated( &
                     settings_c%stability_settings%context, &
                     c_loc(holder%stability_level)))) then
                write(stderr, *) "test_wire_solver_settings_c failed: Settings "// &
                    "wired wrong "//trim(case_names(i_case))//"."
                test_wire_solver_settings_c = .false.
            end if
        end do

    end function test_wire_solver_settings_c

    logical(c_bool) function test_update_orbs_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the orbital update, which calls the
        ! callback function with the level of an extension object held at the address
        ! it receives as context and returns the C wrapper of the Hessian linear
        ! transformation
        !
        use otr_common_c_interface, only: extension_holder_type, &
                                          update_orbs_c_wrapper, hess_x_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type, check_role_called

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: kappa_c(n_param + 1), func_c, grad_c(n_param + 1), &
                      h_diag_c(n_param + 1)
        type(c_funptr) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        ! assume tests pass
        test_update_orbs_c_wrapper = .true.

        ! hold the level of the solver and set up a step, with arrays one element
        ! longer than the number of parameters, of which only these may be read and
        ! written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        holder%solver_level%level => extension%solver_level
        call random_number(kappa_c)
        grad_c = -1.0_c_rp
        h_diag_c = -1.0_c_rp

        ! call the orbital update and determine if it calls the orbital update of the
        ! object and returns its results and error and the C wrapper of the Hessian
        ! linear transformation
        extension%mock_error = 3
        error_c = update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                        hess_x_c_funptr, c_loc(holder%solver_level))
        if (.not. check_role_called("update_orbs_c_wrapper", extension, &
                                    "update_orbs")) test_update_orbs_c_wrapper = .false.
        if (error_c /= 3) then
            write(stderr, *) "test_update_orbs_c_wrapper failed: Error of the "// &
                "object not passed on."
            test_update_orbs_c_wrapper = .false.
        end if
        if (abs(func_c - sum(kappa_c(:n_param))) > tol_c .or. &
            any(abs(grad_c(:n_param) - 2.0_c_rp * kappa_c(:n_param)) > tol_c) .or. &
            any(abs(h_diag_c(:n_param) - 3.0_c_rp * kappa_c(:n_param)) > tol_c)) then
            write(stderr, *) "test_update_orbs_c_wrapper failed: Results of the "// &
                "object not returned."
            test_update_orbs_c_wrapper = .false.
        end if
        if (abs(grad_c(n_param + 1) + 1.0_c_rp) > tol_c .or. &
            abs(h_diag_c(n_param + 1) + 1.0_c_rp) > tol_c) then
            write(stderr, *) "test_update_orbs_c_wrapper failed: Arrays written "// &
                "beyond the number of parameters."
            test_update_orbs_c_wrapper = .false.
        end if
        if (.not. c_associated(hess_x_c_funptr, c_funloc(hess_x_c_wrapper))) then
            write(stderr, *) "test_update_orbs_c_wrapper failed: C wrapper of the "// &
                "Hessian linear transformation not returned."
            test_update_orbs_c_wrapper = .false.
        end if

        ! call the orbital update of an object which returns no Hessian linear
        ! transformation, no C wrapper may be returned
        extension%return_hess_x = .false.
        error_c = update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                        hess_x_c_funptr, c_loc(holder%solver_level))
        if (c_associated(hess_x_c_funptr)) then
            write(stderr, *) "test_update_orbs_c_wrapper failed: C wrapper of the "// &
                "Hessian linear transformation returned although the object "// &
                "returned none."
            test_update_orbs_c_wrapper = .false.
        end if

        ! call the orbital update without a holder, which has to produce an error and
        ! no Hessian linear transformation
        extension%return_hess_x = .true.
        error_c = update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                        hess_x_c_funptr, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_update_orbs_c_wrapper failed: No error produced "// &
                "without a holder."
            test_update_orbs_c_wrapper = .false.
        end if
        if (c_associated(hess_x_c_funptr)) then
            write(stderr, *) "test_update_orbs_c_wrapper failed: C wrapper of the "// &
                "Hessian linear transformation returned without a holder."
            test_update_orbs_c_wrapper = .false.
        end if

    end function test_update_orbs_c_wrapper

    logical(c_bool) function test_obj_func_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the objective function, which forwards
        ! to the host's C function in the bundle of the level of a holder at the
        ! address it receives as context
        !
        use c_interface_unit_tests, only: mock_obj_func
        use otr_common_c_interface, only: extension_holder_type, obj_func_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        real(c_rp) :: kappa_c(n_param), func_c
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_obj_func_c_wrapper = .true.

        ! set up a step
        call random_number(kappa_c)

        ! call the objective function with the addresses of the levels of the solver
        ! and of the stability check of a holder, whose bundles both hold the host's
        ! objective function but only the bundle of the level called holds the host
        ! context, and determine if it forwards to the host's objective function with
        ! that host context and returns its results and error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%obj_func => mock_obj_func
            holder%stability_level%callbacks%obj_func => mock_obj_func
            host_context%mock_error = 3
            error_c = obj_func_c_wrapper(kappa_c, func_c, &
                                         holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_obj_func_c_wrapper failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_obj_func_c_wrapper = .false.
            end if
            if (abs(func_c - sum(kappa_c)) > tol_c) then
                write(stderr, *) "test_obj_func_c_wrapper failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_obj_func_c_wrapper = .false.
            end if
            if (.not. host_context_reached("obj_func_c_wrapper")) &
                test_obj_func_c_wrapper = .false.
        end do

        ! call the objective function without a holder, which has to produce an error
        error_c = obj_func_c_wrapper(kappa_c, func_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_obj_func_c_wrapper failed: No error produced "// &
                "without a holder."
            test_obj_func_c_wrapper = .false.
        end if

    end function test_obj_func_c_wrapper

    logical(c_bool) function test_hess_x_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the Hessian linear transformation, which
        ! calls the callback function with the level of an extension object held at the
        ! address it receives as context
        !
        use otr_common_c_interface, only: extension_holder_type, hess_x_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type, check_role_called, &
                                         level_names

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: x_c(n_param + 1), hess_x_c(n_param + 1)
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_hess_x_c_wrapper = .true.

        ! hold the levels of the extension and set up arrays one element longer than
        ! the number of parameters, of which only these may be read and written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        holder%solver_level%level => extension%solver_level
        holder%stability_level%level => extension%stability_level
        call random_number(x_c)

        ! call the Hessian linear transformation with the addresses of the levels of
        ! the solver and of the stability check in the holder, which both call it, and
        ! determine if it calls that of the object and returns its result and error
        extension%mock_error = 3
        do i_level = 1, 2
            extension%called_role = ""
            hess_x_c = -1.0_c_rp
            error_c = hess_x_c_wrapper(x_c, hess_x_c, &
                                       holder_level_context(holder, i_level))
            if (.not. check_role_called("hess_x_c_wrapper", extension, "hess_x")) &
                test_hess_x_c_wrapper = .false.
            if (error_c /= 3) then
                write(stderr, *) "test_hess_x_c_wrapper failed: Error of the "// &
                    "object not passed on for the "//trim(level_names(i_level))//"."
                test_hess_x_c_wrapper = .false.
            end if
            if (any(abs(hess_x_c(:n_param) - 4.0_c_rp * x_c(:n_param)) > tol_c)) then
                write(stderr, *) "test_hess_x_c_wrapper failed: Result of the "// &
                    "object not returned for the "//trim(level_names(i_level))//"."
                test_hess_x_c_wrapper = .false.
            end if
            if (abs(hess_x_c(n_param + 1) + 1.0_c_rp) > tol_c) then
                write(stderr, *) "test_hess_x_c_wrapper failed: Array written "// &
                    "beyond the number of parameters for the "// &
                    trim(level_names(i_level))//"."
                test_hess_x_c_wrapper = .false.
            end if
        end do

        ! call the Hessian linear transformation without a holder, which has to produce
        ! an error
        error_c = hess_x_c_wrapper(x_c, hess_x_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_hess_x_c_wrapper failed: No error produced "// &
                "without a holder."
            test_hess_x_c_wrapper = .false.
        end if

    end function test_hess_x_c_wrapper

    logical(c_bool) function test_precond_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the level-shifted preconditioner, which
        ! forwards to the host's C function in the bundle of the level of a holder at
        ! the address it receives as context
        !
        use c_interface_unit_tests, only: mock_precond
        use otr_common_c_interface, only: extension_holder_type, precond_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        real(c_rp) :: residual_c(n_param), precond_residual_c(n_param)
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_precond_c_wrapper = .true.

        ! set up a residual
        call random_number(residual_c)

        ! call the preconditioner with the addresses of the levels of the solver and of
        ! the stability check of a holder, whose bundles both hold the host's
        ! preconditioner but only the bundle of the level called holds the host
        ! context, and determine if it forwards to the host's preconditioner with that
        ! host context and returns its results and error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%precond => mock_precond
            holder%stability_level%callbacks%precond => mock_precond
            host_context%mock_error = 3
            error_c = precond_c_wrapper(residual_c, 2.0_c_rp, precond_residual_c, &
                                        holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_precond_c_wrapper failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_precond_c_wrapper = .false.
            end if
            if (any(abs(precond_residual_c - 2.0_c_rp * residual_c) > tol_c)) then
                write(stderr, *) "test_precond_c_wrapper failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_precond_c_wrapper = .false.
            end if
            if (.not. host_context_reached("precond_c_wrapper")) &
                test_precond_c_wrapper = .false.
        end do

        ! call the preconditioner without a holder, which has to produce an error
        error_c = &
            precond_c_wrapper(residual_c, 2.0_c_rp, precond_residual_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_precond_c_wrapper failed: No error produced "// &
                "without a holder."
            test_precond_c_wrapper = .false.
        end if

    end function test_precond_c_wrapper

    logical(c_bool) function test_precond_pd_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the positive-definite preconditioner,
        ! which forwards to the host's C function in the bundle of the level of a
        ! holder at the address it receives as context
        !
        use c_interface_unit_tests, only: mock_precond_pd
        use otr_common_c_interface, only: extension_holder_type, precond_pd_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        real(c_rp) :: residual_c(n_param), precond_residual_c(n_param)
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_precond_pd_c_wrapper = .true.

        ! set up a residual
        call random_number(residual_c)

        ! call the preconditioner with the addresses of the levels of the solver and of
        ! the stability check of a holder, whose bundles both hold the host's
        ! preconditioner but only the bundle of the level called holds the host
        ! context, and determine if it forwards to the host's preconditioner with that
        ! host context and returns its results and error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%precond_pd => mock_precond_pd
            holder%stability_level%callbacks%precond_pd => mock_precond_pd
            host_context%mock_error = 3
            error_c = precond_pd_c_wrapper(residual_c, precond_residual_c, &
                                           holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_precond_pd_c_wrapper failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_precond_pd_c_wrapper = .false.
            end if
            if (any(abs(precond_residual_c - 3.0_c_rp * residual_c) > tol_c)) then
                write(stderr, *) "test_precond_pd_c_wrapper failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_precond_pd_c_wrapper = .false.
            end if
            if (.not. host_context_reached("precond_pd_c_wrapper")) &
                test_precond_pd_c_wrapper = .false.
        end do

        ! call the preconditioner without a holder, which has to produce an error
        error_c = precond_pd_c_wrapper(residual_c, precond_residual_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_precond_pd_c_wrapper failed: No error produced "// &
                "without a holder."
            test_precond_pd_c_wrapper = .false.
        end if

    end function test_precond_pd_c_wrapper

    logical(c_bool) function test_project_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the projection, which forwards to the
        ! host's C function in the bundle of the level of a holder at the address it
        ! receives as context
        !
        use c_interface_unit_tests, only: mock_project
        use otr_common_c_interface, only: extension_holder_type, project_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        real(c_rp) :: vector_c(n_param), projected_c(n_param)
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_project_c_wrapper = .true.

        ! set up a vector
        call random_number(vector_c)

        ! call the projection with the addresses of the levels of the solver and of the
        ! stability check of a holder, whose bundles both hold the host's projection
        ! but only the bundle of the level called holds the host context, and determine
        ! if it forwards to the host's projection with that host context and returns
        ! its results and error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%project => mock_project
            holder%stability_level%callbacks%project => mock_project
            host_context%mock_error = 3
            projected_c = vector_c
            error_c = project_c_wrapper(projected_c, &
                                        holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_project_c_wrapper failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_project_c_wrapper = .false.
            end if
            if (any(abs(projected_c - 2.0_c_rp * vector_c) > tol_c)) then
                write(stderr, *) "test_project_c_wrapper failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_project_c_wrapper = .false.
            end if
            if (.not. host_context_reached("project_c_wrapper")) &
                test_project_c_wrapper = .false.
        end do

        ! call the projection without a holder, which has to produce an error
        error_c = project_c_wrapper(projected_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_project_c_wrapper failed: No error produced "// &
                "without a holder."
            test_project_c_wrapper = .false.
        end if

    end function test_project_c_wrapper

    logical(c_bool) function test_get_extra_trial_vectors_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the extra trial vector function, which
        ! forwards to the host's C function in the bundle of the level of a holder at
        ! the address it receives as context
        !
        use c_interface_unit_tests, only: mock_get_extra_trial_vectors
        use otr_common_c_interface, only: extension_holder_type, &
                                          get_extra_trial_vectors_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        real(c_rp) :: trial_vectors_c(2 * n_param)
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_get_extra_trial_vectors_c_wrapper = .true.

        ! call the extra trial vector function with the addresses of the levels of the
        ! solver and of the stability check of a holder, whose bundles both hold the
        ! host's extra trial vector function but only the bundle of the level called
        ! holds the host context, and determine if it forwards to the host's extra
        ! trial vector function with that host context and returns its results and
        ! error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%get_extra_trial_vectors => &
                mock_get_extra_trial_vectors
            holder%stability_level%callbacks%get_extra_trial_vectors => &
                mock_get_extra_trial_vectors
            host_context%mock_error = 3
            error_c = get_extra_trial_vectors_c_wrapper( &
                trial_vectors_c, 2_c_ip, holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_get_extra_trial_vectors_c_wrapper failed: "// &
                    "Error of the host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_get_extra_trial_vectors_c_wrapper = .false.
            end if
            if (any(abs(trial_vectors_c(n_param + 1:) - 2.0_c_rp) > tol_c)) then
                write(stderr, *) "test_get_extra_trial_vectors_c_wrapper failed: "// &
                    "Result of the host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_get_extra_trial_vectors_c_wrapper = .false.
            end if
            if (.not. host_context_reached("get_extra_trial_vectors_c_wrapper")) &
                test_get_extra_trial_vectors_c_wrapper = .false.
        end do

        ! call the extra trial vector function without a holder, which has to produce
        ! an error
        error_c = get_extra_trial_vectors_c_wrapper(trial_vectors_c, 2_c_ip, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_get_extra_trial_vectors_c_wrapper failed: No "// &
                "error produced without a holder."
            test_get_extra_trial_vectors_c_wrapper = .false.
        end if

    end function test_get_extra_trial_vectors_c_wrapper

    logical(c_bool) function test_conv_check_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the convergence check of the solver,
        ! which forwards to the host's C function in the bundle of the level of a
        ! holder at the address it receives as context
        !
        use c_interface_unit_tests, only: mock_conv_check
        use otr_common_c_interface, only: extension_holder_type, conv_check_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        logical(c_bool) :: converged_c
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_conv_check_c_wrapper = .true.

        ! call the convergence check with the addresses of the levels of the solver and
        ! of the stability check of a holder, whose bundles both hold the host's
        ! convergence check but only the bundle of the level called holds the host
        ! context, and determine if it forwards to the host's convergence check with
        ! that host context and returns its results and error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%conv_check => mock_conv_check
            holder%stability_level%callbacks%conv_check => mock_conv_check
            host_context%mock_error = 3
            error_c = conv_check_c_wrapper(converged_c, &
                                           holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_conv_check_c_wrapper failed: Error of the "// &
                    "host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_conv_check_c_wrapper = .false.
            end if
            if (.not. converged_c) then
                write(stderr, *) "test_conv_check_c_wrapper failed: Result of the "// &
                    "host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_conv_check_c_wrapper = .false.
            end if
            if (.not. host_context_reached("conv_check_c_wrapper")) &
                test_conv_check_c_wrapper = .false.
        end do

        ! call the convergence check without a holder, which has to produce an error
        error_c = conv_check_c_wrapper(converged_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_conv_check_c_wrapper failed: No error produced "// &
                "without a holder."
            test_conv_check_c_wrapper = .false.
        end if

    end function test_conv_check_c_wrapper

    logical(c_bool) function test_stability_hess_x_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the Hessian linear transformation of the
        ! stability check, which forwards to the host's C function in the bundle of the
        ! level of a holder at the address it receives as context
        !
        use c_interface_unit_tests, only: mock_stability_hess_x
        use otr_common_c_interface, only: extension_holder_type, &
                                          stability_hess_x_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        real(c_rp) :: x_c(n_param), hess_x_c(n_param)
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_stability_hess_x_c_wrapper = .true.

        ! set up a vector
        call random_number(x_c)

        ! call the Hessian linear transformation with the addresses of the levels of
        ! the solver and of the stability check of a holder, whose bundles both hold
        ! the host's Hessian linear transformation but only the bundle of the level
        ! called holds the host context, and determine if it forwards to the host's
        ! Hessian linear transformation with that host context and returns its results
        ! and error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%stability_hess_x => mock_stability_hess_x
            holder%stability_level%callbacks%stability_hess_x => mock_stability_hess_x
            host_context%mock_error = 3
            error_c = stability_hess_x_c_wrapper(x_c, hess_x_c, &
                                                 holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_stability_hess_x_c_wrapper failed: Error of "// &
                    "the host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_stability_hess_x_c_wrapper = .false.
            end if
            if (any(abs(hess_x_c - 4.0_c_rp * x_c) > tol_c)) then
                write(stderr, *) "test_stability_hess_x_c_wrapper failed: Result "// &
                    "of the host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_stability_hess_x_c_wrapper = .false.
            end if
            if (.not. host_context_reached("stability_hess_x_c_wrapper")) &
                test_stability_hess_x_c_wrapper = .false.
        end do

        ! call the Hessian linear transformation without a holder, which has to produce
        ! an error
        error_c = stability_hess_x_c_wrapper(x_c, hess_x_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_stability_hess_x_c_wrapper failed: No error "// &
                "produced without a holder."
            test_stability_hess_x_c_wrapper = .false.
        end if

    end function test_stability_hess_x_c_wrapper

    logical(c_bool) function test_approx_hess_x_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the approximate Hessian linear
        ! transformation, which forwards to the host's C function in the bundle of the
        ! level of a holder at the address it receives as context
        !
        use c_interface_unit_tests, only: mock_approx_hess_x
        use otr_common_c_interface, only: extension_holder_type, approx_hess_x_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        real(c_rp) :: x_c(n_param), hess_x_c(n_param)
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_approx_hess_x_c_wrapper = .true.

        ! set up a vector
        call random_number(x_c)

        ! call the approximate Hessian linear transformation with the addresses of the
        ! levels of the solver and of the stability check of a holder, whose bundles
        ! both hold the host's approximate Hessian linear transformation but only the
        ! bundle of the level called holds the host context, and determine if it
        ! forwards to the host's approximate Hessian linear transformation with that
        ! host context and returns its results and error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%approx_hess_x => mock_approx_hess_x
            holder%stability_level%callbacks%approx_hess_x => mock_approx_hess_x
            host_context%mock_error = 3
            error_c = approx_hess_x_c_wrapper(x_c, hess_x_c, &
                                              holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_approx_hess_x_c_wrapper failed: Error of "// &
                    "the host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_approx_hess_x_c_wrapper = .false.
            end if
            if (any(abs(hess_x_c - 4.0_c_rp * x_c) > tol_c)) then
                write(stderr, *) "test_approx_hess_x_c_wrapper failed: Result of "// &
                    "the host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_approx_hess_x_c_wrapper = .false.
            end if
            if (.not. host_context_reached("approx_hess_x_c_wrapper")) &
                test_approx_hess_x_c_wrapper = .false.
        end do

        ! call the approximate Hessian linear transformation without a holder, which
        ! has to produce an error
        error_c = approx_hess_x_c_wrapper(x_c, hess_x_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_approx_hess_x_c_wrapper failed: No error "// &
                "produced without a holder."
            test_approx_hess_x_c_wrapper = .false.
        end if

    end function test_approx_hess_x_c_wrapper

    logical(c_bool) function test_conv_check_stability_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the convergence check of the stability
        ! check, which forwards to the host's C function in the bundle of the level of
        ! a holder at the address it receives as context
        !
        use c_interface_unit_tests, only: mock_conv_check_stability
        use otr_common_c_interface, only: extension_holder_type, &
                                          conv_check_stability_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        real(c_rp) :: residual_c(n_param)
        logical(c_bool) :: converged_c
        integer(c_ip) :: error_c
        integer(ip) :: i_level

        ! assume tests pass
        test_conv_check_stability_c_wrapper = .true.

        ! set up a residual, whose sum the host's convergence check compares to the
        ! eigenvalue
        call random_number(residual_c)

        ! call the convergence check with the addresses of the levels of the solver and
        ! of the stability check of a holder, whose bundles both hold the host's
        ! convergence check but only the bundle of the level called holds the host
        ! context, and determine if it forwards to the host's convergence check with
        ! that host context and returns its results and error
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%conv_check_stability => &
                mock_conv_check_stability
            holder%stability_level%callbacks%conv_check_stability => &
                mock_conv_check_stability
            host_context%mock_error = 3
            error_c = conv_check_stability_c_wrapper( &
                residual_c, sum(residual_c), converged_c, &
                holder_level_context(holder, i_level))
            if (error_c /= 3) then
                write(stderr, *) "test_conv_check_stability_c_wrapper failed: "// &
                    "Error of the host's callback function not passed on for the "// &
                    trim(level_names(i_level))//"."
                test_conv_check_stability_c_wrapper = .false.
            end if
            if (.not. converged_c) then
                write(stderr, *) "test_conv_check_stability_c_wrapper failed: "// &
                    "Result of the host's callback function not returned for the "// &
                    trim(level_names(i_level))//"."
                test_conv_check_stability_c_wrapper = .false.
            end if
            if (.not. host_context_reached("conv_check_stability_c_wrapper")) &
                test_conv_check_stability_c_wrapper = .false.
        end do

        ! call the convergence check without a holder, which has to produce an error
        error_c = conv_check_stability_c_wrapper(residual_c, sum(residual_c), &
                                                 converged_c, c_null_ptr)
        if (error_c == 0) then
            write(stderr, *) "test_conv_check_stability_c_wrapper failed: No error "// &
                "produced without a holder."
            test_conv_check_stability_c_wrapper = .false.
        end if

    end function test_conv_check_stability_c_wrapper

    logical(c_bool) function test_logger_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the logging function, which forwards to
        ! the host's C function in the bundle of the level of a holder at the address
        ! it receives as context
        !
        use c_interface_unit_tests, only: mock_logger
        use otr_common_c_interface, only: extension_holder_type, logger_c_wrapper
        use otr_common_unit_tests, only: level_names
        use test_reference, only: host_context, host_context_reached

        type(extension_holder_type), target :: holder
        integer(ip) :: i_level

        ! assume tests pass
        test_logger_c_wrapper = .true.

        ! call the logging function with the addresses of the levels of the solver and
        ! of the stability check of a holder, whose bundles both hold the host's
        ! logging function but only the bundle of the level called holds the host
        ! context, and determine if it forwards to the host's logging function with
        ! that host context
        do i_level = 1, 2
            call arm_holder_level(holder, i_level)
            holder%solver_level%callbacks%logger => mock_logger
            holder%stability_level%callbacks%logger => mock_logger
            call logger_c_wrapper("test"//c_null_char, &
                                  holder_level_context(holder, i_level))
            if (.not. host_context%logger_called) then
                write(stderr, *) "test_logger_c_wrapper failed: Message not passed "// &
                    "to the host's logging function for the "// &
                    trim(level_names(i_level))//"."
                test_logger_c_wrapper = .false.
            end if
            if (.not. host_context_reached("logger_c_wrapper")) &
                test_logger_c_wrapper = .false.
        end do

        ! call the logging function without a holder, which has to discard the message
        call logger_c_wrapper("test"//c_null_char, c_null_ptr)

    end function test_logger_c_wrapper

end module otr_common_c_interface_unit_tests
