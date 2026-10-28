! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface_unit_tests

    use opentrustregion_core, only: rp, ip, stderr, solver_state_type
    use c_interface, only: c_rp, c_ip, solver_state_type_c
    use otr_common_c_interface, only: evaluate_dm_c_type, get_response_c_type
    use test_reference, only: n_param, is_close
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, c_loc, &
                                           c_associated, c_null_ptr, c_null_funptr, &
                                           c_null_char

    implicit none

    ! whether the extension object of a test of the C wrappers of its roles performs
    ! the role itself or leaves it to the host
    character(len=*), parameter :: role_case_names(2) = &
        [character(len=23) :: "performing the role", "not performing the role"]

    ! state the object of the holder of a test of the C wrappers of the roles reads and
    ! the host's record of the state, which holds the reference values, from which the
    ! state has to be refreshed before a role of the object is called
    type(solver_state_type), target :: held_state
    type(solver_state_type_c), target :: held_state_c

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
        use test_reference, only: check_host_context_c, host_context_error

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

        error = host_context_error(context_c)

    end function mock_evaluate_dm_cs

    function mock_evaluate_dm_os(dm_ao, energy, fock, get_response_c_funptr, &
                                 context_c) result(error) bind(C)
        !
        ! this function is a test function for the density matrix evaluating C function
        ! for 3D density matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao, n_particle
        use otr_common_unit_tests, only: record_mock_call
        use test_reference, only: check_host_context_c, host_context_error

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

        error = host_context_error(context_c)

    end function mock_evaluate_dm_os

    function mock_get_response_cs(dm_ao, response, context_c) result(error) bind(C)
        !
        ! this function is a test function for the response C function for 2D density
        ! matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out), target :: response(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2

        ! check host context
        call check_host_context_c(context_c)

        response(:flat_len) = evaluate_dm_factors(2) * dm_ao(:flat_len)

        error = host_context_error(context_c)

    end function mock_get_response_cs

    function mock_get_response_os(dm_ao, response, context_c) result(error) bind(C)
        !
        ! this function is a test function for the response C function for 3D density
        ! matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao, n_particle
        use test_reference, only: check_host_context_c, host_context_error

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out), target :: response(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2 * n_particle

        ! check host context
        call check_host_context_c(context_c)

        response(:flat_len) = evaluate_dm_factors(2) * dm_ao(:flat_len)

        error = host_context_error(context_c)

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

    subroutine setup_mock_extension_holder(holder, extension, i_level)
        !
        ! this subroutine sets up a holder of the levels of a mock extension whose level
        ! of the solver (1) or of the stability check (2) holds the armed host context
        ! with an error code in its bundle, while the bundle of the other level holds
        ! none, and whose bundle of the solver refers to that of the stability check, so
        ! that a C wrapper calling the object's role or the host's C function with the
        ! wrong level is detected, whose levels refer to the host's record of the state
        ! and whose object reads its own state, and clears the role the object recorded
        !
        use otr_common_c_interface, only: extension_holder_type
        use otr_common_unit_tests, only: mock_extension_type
        use test_reference, only: host_context, get_reference_solver_state_values

        type(extension_holder_type), intent(out), target :: holder
        type(mock_extension_type), intent(inout), target :: extension
        integer(ip), intent(in) :: i_level

        call arm_holder_level(holder, i_level)
        host_context%mock_error = 3
        holder%solver_level%level => extension%solver_level
        holder%stability_level%level => extension%stability_level
        holder%solver_level%callbacks%stability => holder%stability_level%callbacks
        extension%called_role = ""

        ! both levels refer to the host's record holding the reference values, and the
        ! object reads the state with its default values
        held_state = solver_state_type()
        call get_reference_solver_state_values(held_state_c)
        holder%solver_level%state_c => held_state_c
        holder%stability_level%state_c => held_state_c
        extension%state => held_state

    end subroutine setup_mock_extension_holder

    function check_c_factory_wiring(test_name, holder, state_c, update_orbs_c_funptr, &
                                    obj_func_c_funptr, solver_settings_c, case_name) &
        result(test_passed)
        !
        ! this function checks that the C wrapper of a factory wired the settings of
        ! the holder it returned with an object, which reads the state of the holder,
        ! that the holder holds the host's record of the state at both levels, and that
        ! the wrapper returned the shared C wrappers of the orbital update and the
        ! objective function and wrote the solver settings back with the address of the
        ! level of the solver of its holder as context, which has to hold the levels of
        ! the object and the settings it wired with the object's level of the solver as
        ! context
        !
        use c_interface, only: solver_settings_type_c, solver_state_type_c
        use otr_common, only: extension_level_type
        use otr_common_c_interface, only: extension_holder_type, &
                                          update_orbs_c_wrapper, obj_func_c_wrapper
        use otr_common_unit_tests, only: wired_level

        character(len=*), intent(in) :: test_name
        type(extension_holder_type), intent(in), target :: holder
        type(solver_state_type_c), intent(in), target :: state_c
        type(c_funptr), intent(in) :: update_orbs_c_funptr, obj_func_c_funptr
        type(solver_settings_type_c), intent(in) :: solver_settings_c
        character(len=*), intent(in), optional :: case_name
        logical :: test_passed

        type(extension_level_type), pointer :: level
        character(len=:), allocatable :: suffix

        ! assume tests pass
        test_passed = .true.

        ! name the case in the messages
        suffix = ""
        if (present(case_name)) suffix = " for the "//case_name//" case"

        if (.not. c_associated(update_orbs_c_funptr, &
                               c_funloc(update_orbs_c_wrapper)) .or. &
            .not. c_associated(obj_func_c_funptr, c_funloc(obj_func_c_wrapper))) then
            write(stderr, *) "test_"//test_name// &
                " failed: C wrappers of the object not returned"//suffix//"."
            test_passed = .false.
        end if

        ! the settings of the holder, which own the object, have to be wired with it,
        ! whose inspection would crash without one
        level => wired_level(holder%solver_settings)
        if (.not. associated(level)) then
            write(stderr, *) "test_"//test_name// &
                " failed: Settings of the holder not wired with an object"//suffix//"."
            test_passed = .false.
            return
        end if
        if (.not. associated(level%extension%state, holder%state)) then
            write(stderr, *) "test_"//test_name// &
                " failed: Object does not read the state of the holder"//suffix//"."
            test_passed = .false.
        end if
        if (.not. (associated(holder%solver_level%state_c, state_c) .and. &
                   associated(holder%stability_level%state_c, state_c))) then
            write(stderr, *) "test_"//test_name// &
                " failed: Host's record of the state not held"//suffix//"."
            test_passed = .false.
        end if
        if (.not. c_associated(solver_settings_c%context, &
                               c_loc(holder%solver_level)) .or. &
            .not. associated(holder%solver_level%level, level) .or. .not. associated( &
                holder%stability_level%level, level%extension%stability_level) .or. &
            .not. associated(holder%solver_settings%context, level)) then
            write(stderr, *) "test_"//test_name//" failed: Solver settings not "// &
                "wired with the holder of the object"//suffix//"."
            test_passed = .false.
        end if

    end function check_c_factory_wiring

    function check_role_c_wrapper(test_name, extension, role, i_case, i_level, &
                                  error_c) result(passed)
        !
        ! this function checks that a C wrapper of a role, called with the address of
        ! the level of the solver (1) or of the stability check (2) in a holder of a
        ! mock extension, called the role of the object with the level, after
        ! refreshing the state of the object from the host's record, for an object
        ! performing it (case 1) and otherwise the host's C function with the host
        ! context of the level (case 2), and passed on the error of the one it called
        !
        use otr_common_unit_tests, only: mock_extension_type, check_role_called, &
                                         level_context, level_names
        use test_reference, only: host_context_reached, operator(/=)

        character(len=*), intent(in) :: test_name, role
        type(mock_extension_type), intent(in), target :: extension
        integer(ip), intent(in) :: i_case, i_level
        integer(c_ip), intent(in) :: error_c
        logical :: passed

        class(*), pointer :: level
        character(len=:), allocatable :: case_name

        ! assume tests pass
        passed = .true.

        ! name the level and the case in the messages
        case_name = "the "//trim(level_names(i_level))//" of an object "// &
                    trim(role_case_names(i_case))

        ! the object's role has to be called with the level if the object performs it
        ! and the host's C function with the host context of the level otherwise
        if (i_case == 1) then
            level => level_context(extension, i_level)
            if (.not. check_role_called(test_name, extension, role, level, case_name)) &
                passed = .false.
            if (held_state_c /= held_state) then
                write(stderr, *) "test_"//test_name//" failed: State not refreshed "// &
                    "from the host's record for "//case_name//"."
                passed = .false.
            end if
        else
            if (extension%called_role /= "") then
                write(stderr, *) "test_"//test_name//" failed: Role "//role// &
                    " of the object called for "//case_name//"."
                passed = .false.
            end if
            if (.not. host_context_reached(test_name)) then
                write(stderr, *) "test_"//test_name//" failed: Host's C function "// &
                    "not called with the host context for "//case_name//"."
                passed = .false.
            end if
        end if

        ! the error of the one called has to be passed on
        if (error_c /= 3) then
            write(stderr, *) "test_"//test_name//" failed: Error not passed on for "// &
                case_name//"."
            passed = .false.
        end if

    end function check_role_c_wrapper

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

    logical(c_bool) function test_bundle_settings_c() bind(C)
        !
        ! this function tests the subroutine which hands the logging function and the
        ! host context of the C settings of an extension to the bundle of the host's
        ! callback functions of the extension
        !
        use c_interface, only: c_callbacks_type
        use otr_common_c_interface, only: bundle_settings_c
        use c_interface_unit_tests, only: mock_logger
        use test_reference, only: host_context

        type(c_callbacks_type) :: callbacks
        type(c_funptr) :: logger_c_funptr
        integer(ip) :: i_case
        character(len=42), parameter :: case_names(3) = &
            [character(len=42) :: "for initialized settings", "for initialized "// &
             "settings without a logger", "for settings that were not initialized"]

        ! assume tests pass
        test_bundle_settings_c = .true.

        ! hand over the logging function and host context of initialized settings with
        ! and without a logging function and of settings that were not initialized, of
        ! which only the initialized ones may hand over their host context and only
        ! those with a logging function a logging function
        do i_case = 1, size(case_names)
            callbacks = c_callbacks_type()
            logger_c_funptr = c_funloc(mock_logger)
            if (i_case == 2) logger_c_funptr = c_null_funptr
            call bundle_settings_c(callbacks, logical(i_case /= 3, kind=c_bool), &
                                   logger_c_funptr, c_loc(host_context))
            if ((associated(callbacks%logger, mock_logger) .neqv. i_case == 1) .or. &
                (c_associated(callbacks%host_context, c_loc(host_context)) .neqv. &
                 i_case /= 3)) then
                write(stderr, *) "test_bundle_settings_c failed: Logging function "// &
                    "or host context handed over wrong "//trim(case_names(i_case))//"."
                test_bundle_settings_c = .false.
            end if
        end do

    end function test_bundle_settings_c

    logical(c_bool) function test_setup_extension_holder() bind(C)
        !
        ! this function tests the subroutine which sets up the holder of an extension
        ! object for a C factory from the host's callback functions and C solver
        ! settings
        !
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              setup_error_logging
        use c_interface, only: solver_settings_type_c
        use c_interface_unit_tests, only: &
            mock_logger, mock_update_orbs, mock_obj_func, mock_stability_precond, &
            mock_stability_project, mock_stability_logger, &
            mock_stability_get_extra_trial_vectors, mock_approx_hess_x, &
            mock_stability_conv_check
        use otr_common, only: orbital_settings_type
        use otr_common_c_interface, only: extension_holder_type, setup_extension_holder
        use otr_common_unit_tests, only: mock_extension_type
        use c_interface, only: c_callbacks_type
        use test_reference, only: host_context_type, get_reference_solver_values, &
                                  get_reference_solver_state_values, host_context, &
                                  unset_callbacks, ref_solver_settings, operator(/=)

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        type(solver_state_type_c), target :: state_c
        type(solver_settings_type_c) :: settings_c
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
        call setup_settings(settings, log_context)

        ! set up the holder with C solver settings whose nested settings provide
        ! nothing or a single context or callback function of their own, and with a
        ! context of nested settings or settings that were not initialized, of which
        ! only nested settings that were initialized and provide something of their own
        ! may keep the bundle of the holder as context, the holder has to be set up
        ! anew and its bundles have to hold the given callback functions and be handed
        ! to the settings converted into the holder
        do i_case = 1, size(case_names)
            call get_reference_solver_values(settings_c)
            call unset_callbacks(settings_c%stability_settings)
            settings_c%stability_settings%context = c_null_ptr
            settings_c%logger = c_funloc(mock_logger)
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
                    nested_c%conv_check = c_funloc(mock_stability_conv_check)
                end select
                if (i_case == 9) nested_c%initialized = .false.
            end associate
            if (i_case == 10) settings_c%initialized = .false.
            holder%solver_level%level => extension%solver_level

            call setup_extension_holder(holder, settings_c, c_loc(state_c), settings, &
                                        error, c_funloc(mock_update_orbs), &
                                        c_funloc(mock_obj_func))
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
            if (.not. associated(holder%solver_settings%context, &
                                 holder%solver_level%callbacks) .or. &
                .not. associated(holder%solver_level%callbacks%stability, &
                                 holder%stability_level%callbacks)) then
                write(stderr, *) "test_setup_extension_holder failed: Settings not "// &
                    "converted with the bundles of the holder "// &
                    trim(case_names(i_case))//"."
                test_setup_extension_holder = .false.
            end if
            if (i_case == 1 .and. &
                holder%solver_settings%solver_settings_type /= ref_solver_settings) then
                write(stderr, *) "test_setup_extension_holder failed: Given "// &
                    "settings not converted."
                test_setup_extension_holder = .false.
            end if
            if (associated(holder%solver_settings%stability_settings%context, &
                           holder%stability_level%callbacks) .neqv. &
                nested_own(i_case)) then
                write(stderr, *) "test_setup_extension_holder failed: Context of "// &
                    "the nested settings wrong "//trim(case_names(i_case))//"."
                test_setup_extension_holder = .false.
            end if
        end do

        ! set up the holder for an object which wraps no orbital update or objective
        ! function of the host, whose bundle has to hold neither
        call get_reference_solver_values(settings_c)
        call setup_extension_holder(holder, settings_c, c_loc(state_c), settings, error)
        if (error /= 0 .or. associated(holder%solver_level%callbacks%update_orbs) .or. &
            associated(holder%solver_level%callbacks%obj_func)) then
            write(stderr, *) "test_setup_extension_holder failed: Holder set up "// &
                "wrong without the host's orbital update and objective function."
            test_setup_extension_holder = .false.
        end if

        ! set up the holder with the host's record of the state, which its levels have
        ! to hold and which, as the state of the holder, has to be reset to the default
        ! state from the values it held before, and which has to keep the bundle the
        ! extension's settings point to
        call get_reference_solver_values(settings_c)
        call get_reference_solver_state_values(state_c)
        holder%state%max_precision_reached = .true.
        allocate(c_callbacks_type :: holder%extension_callbacks)
        call setup_extension_holder(holder, settings_c, c_loc(state_c), settings, error)
        if (error /= 0 .or. .not. allocated(holder%extension_callbacks)) then
            write(stderr, *) "test_setup_extension_holder failed: Holder set up "// &
                "wrong with a record of the state."
            test_setup_extension_holder = .false.
        end if
        if (.not. (associated(holder%solver_level%state_c, state_c) .and. &
                   associated(holder%stability_level%state_c, state_c))) then
            write(stderr, *) "test_setup_extension_holder failed: Record of the "// &
                "state not held by both levels."
            test_setup_extension_holder = .false.
        else if (state_c /= solver_state_type() .or. state_c /= holder%state) then
            write(stderr, *) "test_setup_extension_holder failed: Record of the "// &
                "state or state of the holder not reset to the default state."
            test_setup_extension_holder = .false.
        end if

        ! set up the holder without a record of the state, which has to produce an
        ! error and leave the holder untouched
        call setup_error_logging(settings, log_context)
        call get_reference_solver_values(settings_c)
        holder%solver_level%level => extension%solver_level
        log_context%log_message = ""
        call setup_extension_holder(holder, settings_c, c_null_ptr, settings, error)
        if (error == 0 .or. log_context%log_message == "" .or. &
            .not. associated(holder%solver_level%level, extension%solver_level)) then
            write(stderr, *) "test_setup_extension_holder failed: No error "// &
                "returned and printed or holder changed without a record of the state."
            test_setup_extension_holder = .false.
        end if

    end function test_setup_extension_holder

    logical(c_bool) function test_refresh_state() bind(C)
        !
        ! this function tests the subroutine which refreshes the state the object of a
        ! level of a holder reads from the host's record the level holds
        !
        use otr_common_c_interface, only: holder_level_type, refresh_state
        use otr_common_unit_tests, only: mock_extension_type
        use test_reference, only: get_reference_solver_state_values, operator(/=)

        type(mock_extension_type), target :: extension
        type(holder_level_type) :: held
        type(solver_state_type), target :: state
        type(solver_state_type_c), target :: state_c

        ! assume tests pass
        test_refresh_state = .true.

        ! a level holding the host's record, which holds the reference values, has to
        ! refresh the state the object reads from it
        call get_reference_solver_state_values(state_c)
        state = solver_state_type()
        extension%state => state
        extension%solver_level%extension => extension
        held%level => extension%solver_level
        held%state_c => state_c
        call refresh_state(held)
        if (state_c /= state) then
            write(stderr, *) "test_refresh_state failed: State not refreshed from "// &
                "the host's record."
            test_refresh_state = .false.
        end if

        ! a level without a record has to leave the state the object reads unchanged
        state = solver_state_type()
        held%state_c => null()
        call refresh_state(held)
        if (.not. (state_c /= state)) then
            write(stderr, *) "test_refresh_state failed: State changed without a "// &
                "record."
            test_refresh_state = .false.
        end if

        ! an object without a state has to be left alone
        extension%state => null()
        held%state_c => state_c
        call refresh_state(held)
        if (associated(extension%state) .or. .not. (state_c /= state)) then
            write(stderr, *) "test_refresh_state failed: State of an object "// &
                "without one changed."
            test_refresh_state = .false.
        end if

    end function test_refresh_state

    logical(c_bool) function test_assign_solver_state_f_c() bind(C)
        !
        ! this function tests that the function that converts a solver state from C to
        ! Fortran correctly performs this conversion
        !
        use otr_common_c_interface, only: assignment(=)
        use test_reference, only: get_reference_solver_state_values, operator(/=)

        type(solver_state_type) :: state
        type(solver_state_type_c) :: state_c

        ! assume tests pass
        test_assign_solver_state_f_c = .true.

        ! convert the reference state, which is set field by field without the
        ! conversions under test and whose fields all differ from each other, and check
        ! every field, including those of the nested stability check state
        call get_reference_solver_state_values(state_c)
        state = state_c
        if (state_c /= state) then
            write(stderr, *) "test_assign_solver_state_f_c failed: State not "// &
                "converted correctly."
            test_assign_solver_state_f_c = .false.
        end if

        ! convert the reference state with the opposite maximum precision flag, so that
        ! a conversion which writes a constant flag is detected
        state_c%max_precision_reached = .not. state_c%max_precision_reached
        state = state_c
        if (state_c /= state) then
            write(stderr, *) "test_assign_solver_state_f_c failed: State with the "// &
                "opposite maximum precision flag not converted correctly."
            test_assign_solver_state_f_c = .false.
        end if

    end function test_assign_solver_state_f_c

    logical(c_bool) function test_assign_stability_state_f_c() bind(C)
        !
        ! this function tests that the function that converts a stability check state
        ! from C to Fortran correctly performs this conversion
        !
        use opentrustregion_core, only: stability_state_type
        use c_interface, only: stability_state_type_c
        use otr_common_c_interface, only: assignment(=)
        use test_reference, only: get_reference_stability_state_values, operator(/=)

        type(stability_state_type) :: state
        type(stability_state_type_c) :: state_c

        ! assume tests pass
        test_assign_stability_state_f_c = .true.

        ! convert the reference state, which is set field by field without the
        ! conversions under test and whose fields all differ from each other, and check
        ! every field
        call get_reference_stability_state_values(state_c)
        state = state_c
        if (state_c /= state) then
            write(stderr, *) "test_assign_stability_state_f_c failed: State not "// &
                "converted correctly."
            test_assign_stability_state_f_c = .false.
        end if

    end function test_assign_stability_state_f_c

    logical(c_bool) function test_release_extension() bind(C)
        !
        ! this function tests the subroutine which frees the holder of an extension
        ! object a C factory returned, together with the object
        !
        use c_interface, only: c_callbacks_type
        use otr_common, only: extension_type
        use otr_common_c_interface, only: extension_holder_type, release_extension
        use otr_common_unit_tests, only: setup_owned_mock_extension, freed_extension_ids

        type(extension_holder_type), pointer :: holder
        class(extension_type), allocatable, target :: object

        ! assume tests pass
        test_release_extension = .true.

        ! release a holder whose settings own an object and which holds a bundle, which
        ! frees the object
        allocate(holder)
        call setup_owned_mock_extension(object, 1_ip)
        allocate(holder%solver_settings%extensions(1))
        call move_alloc(object, holder%solver_settings%extensions(1)%object)
        allocate(c_callbacks_type :: holder%extension_callbacks)
        call release_extension(c_loc(holder))
        if (size(freed_extension_ids) /= 1) then
            write(stderr, *) "test_release_extension failed: Object not freed."
            test_release_extension = .false.
        end if

        ! releasing no extension does nothing
        call release_extension(c_null_ptr)

    end function test_release_extension

    logical(c_bool) function test_wire_solver_settings_c() bind(C)
        !
        ! this function tests the subroutine which writes the solver settings a Fortran
        ! factory wired back into the host's C solver settings
        !
        use opentrustregion_core, only: solver_settings_type
        use c_interface, only: solver_settings_type_c
        use otr_common, only: precond_callback, precond_pd_callback, project_callback, &
                              get_extra_trial_vectors_callback, conv_check_callback, &
                              stability_hess_x_callback, approx_hess_x_callback, &
                              logger_callback
        use otr_common_c_interface, only: &
            extension_holder_type, wire_solver_settings_c, precond_c_wrapper, &
            precond_pd_c_wrapper, project_c_wrapper, &
            get_extra_trial_vectors_c_wrapper, conv_check_c_wrapper, &
            stability_hess_x_c_wrapper, approx_hess_x_c_wrapper, logger_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type
        use test_reference, only: get_reference_solver_values, unset_callbacks, &
                                  same_callbacks, ref_solver_settings, operator(/=)

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        type(solver_settings_type) :: solver_settings
        type(solver_settings_type_c) :: settings_c, expected_c
        integer(ip) :: i_case
        character(len=48), parameter :: case_names(2) = &
            [character(len=48) :: "with every role wired at both levels", &
             "with nested settings that were not wired"]

        ! assume tests pass
        test_wire_solver_settings_c = .true.

        ! the level of the solver of the extension object points back at it, as after
        ! its wiring
        extension%solver_level%extension => extension

        ! write back solver settings with the reference values into which a factory
        ! wired every role of both levels, or a preconditioner, a projection and a
        ! logging function while leaving the nested settings untouched, into C settings
        ! that were not initialized and hold other callback functions and contexts,
        ! which have to be replaced by the C wrappers of the wired roles and the
        ! addresses of the levels of the holder, while the roles that were not wired
        ! stay unset, and the holder has to hold the levels of the extension object
        ! whose level of the solver the settings of the holder hold as context
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
                solver_settings%precond_pd => precond_pd_callback
                solver_settings%get_extra_trial_vectors => &
                    get_extra_trial_vectors_callback
                solver_settings%conv_check => conv_check_callback
                solver_settings%stability_hess_x => stability_hess_x_callback
                expected_c%precond_pd = c_funloc(precond_pd_c_wrapper)
                expected_c%get_extra_trial_vectors = &
                    c_funloc(get_extra_trial_vectors_c_wrapper)
                expected_c%conv_check = c_funloc(conv_check_c_wrapper)
                expected_c%stability_hess_x = c_funloc(stability_hess_x_c_wrapper)
                associate (nested => solver_settings%stability_settings, &
                           nested_c => expected_c%stability_settings)
                    nested%context => extension%stability_level
                    nested%precond => precond_callback
                    nested%project => project_callback
                    nested%get_extra_trial_vectors => get_extra_trial_vectors_callback
                    nested%approx_hess_x => approx_hess_x_callback
                    nested%conv_check => conv_check_callback
                    nested%logger => logger_callback
                    nested_c%precond = c_funloc(precond_c_wrapper)
                    nested_c%project = c_funloc(project_c_wrapper)
                    nested_c%get_extra_trial_vectors = &
                        c_funloc(get_extra_trial_vectors_c_wrapper)
                    nested_c%approx_hess_x = c_funloc(approx_hess_x_c_wrapper)
                    nested_c%conv_check = c_funloc(conv_check_c_wrapper)
                    nested_c%logger = c_funloc(logger_c_wrapper)
                end associate
            end if
            call get_reference_solver_values(settings_c)
            settings_c%initialized = .false.

            holder%solver_settings%solver_settings_type = solver_settings
            call wire_solver_settings_c(holder, settings_c)
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

    logical(c_bool) function test_finish_extension_holder() bind(C)
        !
        ! this function tests the subroutine which finishes the factory call of a C host
        ! by writing back the solver settings and returning the C wrappers and the
        ! holder, or by freeing the holder if the factory failed
        !
        use c_interface, only: solver_settings_type_c
        use otr_common_c_interface, only: extension_holder_type, &
                                          finish_extension_holder, &
                                          update_orbs_c_wrapper, obj_func_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type
        use test_reference, only: get_reference_solver_values

        type(mock_extension_type), target :: extension
        type(extension_holder_type), pointer :: holder
        type(solver_settings_type_c) :: settings_c
        type(c_funptr) :: update_orbs_c_funptr, obj_func_c_funptr
        type(c_ptr) :: extension_c

        ! assume tests pass
        test_finish_extension_holder = .true.

        ! finish a successful factory call with a holder whose settings the factory
        ! wired with an extension object, which has to write them back with the address
        ! of the level of the solver of the holder as context and return the C wrappers
        ! of the orbital update and the objective function and the holder
        extension%solver_level%extension => extension
        allocate(holder)
        holder%solver_settings%context => extension%solver_level
        call get_reference_solver_values(settings_c)
        update_orbs_c_funptr = c_null_funptr
        obj_func_c_funptr = c_null_funptr
        extension_c = c_null_ptr
        call finish_extension_holder(holder, 0_ip, settings_c, update_orbs_c_funptr, &
                                     obj_func_c_funptr, extension_c)
        if (.not. associated(holder)) then
            write(stderr, *) "test_finish_extension_holder failed: Holder freed "// &
                "after a successful factory call."
            test_finish_extension_holder = .false.
            return
        end if
        if (.not. c_associated(settings_c%context, c_loc(holder%solver_level)) .or. &
            .not. associated(holder%solver_level%level, extension%solver_level)) then
            write(stderr, *) "test_finish_extension_holder failed: Solver settings "// &
                "not written back with the holder."
            test_finish_extension_holder = .false.
        end if
        if (.not. c_associated(update_orbs_c_funptr, &
                               c_funloc(update_orbs_c_wrapper)) .or. &
            .not. c_associated(obj_func_c_funptr, c_funloc(obj_func_c_wrapper)) .or. &
            .not. c_associated(extension_c, c_loc(holder))) then
            write(stderr, *) "test_finish_extension_holder failed: C wrappers or "// &
                "holder not returned."
            test_finish_extension_holder = .false.
        end if
        deallocate(holder)

        ! finish a failed factory call, which has to free the holder and leave the
        ! functions and the extension unset
        allocate(holder)
        update_orbs_c_funptr = c_null_funptr
        obj_func_c_funptr = c_null_funptr
        extension_c = c_null_ptr
        call finish_extension_holder(holder, 1_ip, settings_c, update_orbs_c_funptr, &
                                     obj_func_c_funptr, extension_c)
        if (associated(holder) .or. c_associated(update_orbs_c_funptr) .or. &
            c_associated(obj_func_c_funptr) .or. c_associated(extension_c)) then
            write(stderr, *) "test_finish_extension_holder failed: Holder not "// &
                "freed or functions or extension returned after a failed factory call."
            test_finish_extension_holder = .false.
        end if

    end function test_finish_extension_holder

    logical(c_bool) function test_update_orbs_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the orbital update, which calls the
        ! callback function with the level of an extension object held at the address
        ! it receives as context, or the host's C function for an object that does not
        ! perform it, whose Hessian linear transformation it keeps in the bundles of
        ! both levels, and returns the C wrapper of the Hessian linear transformation
        !
        use c_interface_unit_tests, only: mock_update_orbs, &
                                          mock_update_orbs_no_hess_x, mock_hess_x
        use otr_common_c_interface, only: extension_holder_type, &
                                          update_orbs_c_wrapper, hess_x_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type
        use test_reference, only: host_context_reached

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: kappa_c(n_param + 1), func_c, grad_c(n_param + 1), &
                      h_diag_c(n_param + 1)
        type(c_funptr) :: hess_x_c_funptr
        integer(c_ip) :: error_c
        integer(ip) :: i_case
        character(len=*), parameter :: no_hess_x_case_names(3) = &
            [character(len=53) :: "an object returning none", "a host returning none", &
             "a host returning none to an object performing its own"]

        ! assume tests pass
        test_update_orbs_c_wrapper = .true.

        ! set up the extension and a step, with arrays one element longer than the
        ! number of parameters, of which only these may be read and written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        extension%mock_error = 3
        call random_number(kappa_c)

        ! call the orbital update of an object that performs it and of one that does
        ! not, whose holder bundles the host's, with the address of the level of the
        ! solver, the only level that calls it, and determine if it calls the object's
        ! or the host's, returns its results and error and the C wrapper of the Hessian
        ! linear transformation and keeps the host's in the bundles of both levels
        do i_case = 1, 2
            extension%provides%update_orbs = i_case == 1
            call setup_mock_extension_holder(holder, extension, 1_ip)
            holder%solver_level%callbacks%update_orbs => mock_update_orbs
            grad_c = -1.0_c_rp
            h_diag_c = -1.0_c_rp
            error_c = update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                            hess_x_c_funptr, c_loc(holder%solver_level))
            if (.not. check_role_c_wrapper("update_orbs_c_wrapper", extension, &
                                           "update_orbs", i_case, 1_ip, error_c)) &
                test_update_orbs_c_wrapper = .false.
            if (.not. is_close(func_c, sum(kappa_c(:n_param))) .or. &
                .not. all(is_close(grad_c(:n_param), &
                                   2.0_c_rp * kappa_c(:n_param))) .or. &
                .not. all(is_close(h_diag_c(:n_param), 3.0_c_rp * kappa_c(:n_param)))) &
                then
                write(stderr, *) "test_update_orbs_c_wrapper failed: Results not "// &
                    "returned for an object "//trim(role_case_names(i_case))//"."
                test_update_orbs_c_wrapper = .false.
            end if
            if (.not. is_close(grad_c(n_param + 1), -1.0_c_rp) .or. &
                .not. is_close(h_diag_c(n_param + 1), -1.0_c_rp)) then
                write(stderr, *) "test_update_orbs_c_wrapper failed: Arrays "// &
                    "written beyond the number of parameters for an object "// &
                    trim(role_case_names(i_case))//"."
                test_update_orbs_c_wrapper = .false.
            end if
            if (.not. c_associated(hess_x_c_funptr, c_funloc(hess_x_c_wrapper))) then
                write(stderr, *) "test_update_orbs_c_wrapper failed: C wrapper of "// &
                    "the Hessian linear transformation not returned for an object "// &
                    trim(role_case_names(i_case))//"."
                test_update_orbs_c_wrapper = .false.
            end if
            if (i_case == 2 .and. .not. ( &
                associated(holder%solver_level%callbacks%hess_x, mock_hess_x) .and. &
                associated(holder%stability_level%callbacks%hess_x, mock_hess_x))) then
                write(stderr, *) "test_update_orbs_c_wrapper failed: Hessian "// &
                    "linear transformation of the host not kept in the bundles of "// &
                    "both levels."
                test_update_orbs_c_wrapper = .false.
            end if
        end do

        ! call the orbital update of an object that performs it and returns no Hessian
        ! linear transformation, and of one that does not, whose host returns none,
        ! without and with a Hessian linear transformation of the object itself, after
        ! leaving one of an earlier orbital update in the bundles, and determine if a C
        ! wrapper is returned only where a Hessian linear transformation can be called
        ! and the host's earlier one discarded, and that the host's orbital update is
        ! called with the host context and its error passed on
        do i_case = 1, size(no_hess_x_case_names)
            extension%provides%update_orbs = i_case == 1
            extension%provides%hess_x = i_case == 3
            extension%return_hess_x = i_case /= 1
            call setup_mock_extension_holder(holder, extension, 1_ip)
            holder%solver_level%callbacks%update_orbs => mock_update_orbs_no_hess_x
            holder%solver_level%callbacks%hess_x => mock_hess_x
            holder%stability_level%callbacks%hess_x => mock_hess_x
            error_c = update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                            hess_x_c_funptr, c_loc(holder%solver_level))
            if (error_c /= 3) then
                write(stderr, *) "test_update_orbs_c_wrapper failed: Error not "// &
                    "passed on for "//trim(no_hess_x_case_names(i_case))//"."
                test_update_orbs_c_wrapper = .false.
            end if
            if (i_case > 1) then
                if (.not. host_context_reached("update_orbs_c_wrapper")) then
                    write(stderr, *) "test_update_orbs_c_wrapper failed: Host's C "// &
                        "function not called with the host context for "// &
                        trim(no_hess_x_case_names(i_case))//"."
                    test_update_orbs_c_wrapper = .false.
                end if
            end if
            if (c_associated(hess_x_c_funptr) .neqv. i_case == 3) then
                write(stderr, *) "test_update_orbs_c_wrapper failed: C wrapper of "// &
                    "the Hessian linear transformation returned wrongly for "// &
                    trim(no_hess_x_case_names(i_case))//"."
                test_update_orbs_c_wrapper = .false.
            end if
            if (i_case > 1 .and. &
                (associated(holder%solver_level%callbacks%hess_x) .or. &
                 associated(holder%stability_level%callbacks%hess_x))) then
                write(stderr, *) "test_update_orbs_c_wrapper failed: Earlier "// &
                    "Hessian linear transformation of the host kept for "// &
                    trim(no_hess_x_case_names(i_case))//"."
                test_update_orbs_c_wrapper = .false.
            end if
        end do

        ! call the orbital update without a holder, which has to produce an error and
        ! no Hessian linear transformation
        error_c = update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                        hess_x_c_funptr, c_null_ptr)
        if (error_c /= 1) then
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
        ! this function tests the C wrapper of the objective function, which calls the
        ! callback function with the level of an extension object held at the address
        ! it receives as context, or the host's C function for an object that does not
        ! perform it
        !
        use c_interface_unit_tests, only: mock_obj_func
        use otr_common_c_interface, only: extension_holder_type, obj_func_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type, level_names

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: kappa_c(n_param + 1), func_c
        integer(c_ip) :: error_c
        integer(ip) :: i_case, i_level

        ! assume tests pass
        test_obj_func_c_wrapper = .true.

        ! set up the extension and arrays one element longer than the number of
        ! parameters, of which only these may be read and written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        extension%mock_error = 3
        call random_number(kappa_c)

        ! call the objective function of an object that performs it and of one that
        ! does not, whose holder bundles the host's, with the addresses of the levels
        ! of the solver and of the stability check in the holder and determine if it
        ! calls the object's with the level or the host's with the host context of the
        ! level and returns its result and error
        do i_case = 1, 2
            extension%provides%obj_func = i_case == 1
            do i_level = 1, 2
                call setup_mock_extension_holder(holder, extension, i_level)
                holder%solver_level%callbacks%obj_func => mock_obj_func
                holder%stability_level%callbacks%obj_func => mock_obj_func
                error_c = obj_func_c_wrapper(kappa_c, func_c, &
                                             holder_level_context(holder, i_level))
                if (.not. check_role_c_wrapper("obj_func_c_wrapper", extension, &
                                               "obj_func", i_case, i_level, error_c)) &
                    test_obj_func_c_wrapper = .false.
                if (.not. is_close(func_c, sum(kappa_c(:n_param)))) then
                    write(stderr, *) "test_obj_func_c_wrapper failed: Result not "// &
                        "returned for the "//trim(level_names(i_level))// &
                        " of an object "//trim(role_case_names(i_case))//"."
                    test_obj_func_c_wrapper = .false.
                end if
            end do
        end do

        ! call the objective function without a holder, which has to produce an error
        error_c = obj_func_c_wrapper(kappa_c, func_c, c_null_ptr)
        if (error_c /= 1) then
            write(stderr, *) "test_obj_func_c_wrapper failed: No error produced "// &
                "without a holder."
            test_obj_func_c_wrapper = .false.
        end if

    end function test_obj_func_c_wrapper

    logical(c_bool) function test_hess_x_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the Hessian linear transformation, which
        ! calls the callback function with the level of an extension object held at the
        ! address it receives as context, or the host's C function for an object that
        ! does not perform it
        !
        use c_interface_unit_tests, only: mock_hess_x
        use otr_common_c_interface, only: extension_holder_type, hess_x_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type, level_names

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: x_c(n_param + 1), hess_x_c(n_param + 1)
        integer(c_ip) :: error_c
        integer(ip) :: i_case, i_level

        ! assume tests pass
        test_hess_x_c_wrapper = .true.

        ! set up the extension and arrays one element longer than the number of
        ! parameters, of which only these may be read and written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        extension%mock_error = 3
        call random_number(x_c)

        ! call the Hessian linear transformation of an object that performs it and of
        ! one that does not, whose holder bundles the host's, with the addresses of the
        ! levels of the solver and of the stability check in the holder and determine
        ! if it calls the object's with the level or the host's with the host context
        ! of the level and returns its result and error
        do i_case = 1, 2
            extension%provides%hess_x = i_case == 1
            do i_level = 1, 2
                call setup_mock_extension_holder(holder, extension, i_level)
                holder%solver_level%callbacks%hess_x => mock_hess_x
                holder%stability_level%callbacks%hess_x => mock_hess_x
                hess_x_c = -1.0_c_rp
                error_c = hess_x_c_wrapper(x_c, hess_x_c, &
                                           holder_level_context(holder, i_level))
                if (.not. check_role_c_wrapper("hess_x_c_wrapper", extension, &
                                               "hess_x", i_case, i_level, error_c)) &
                    test_hess_x_c_wrapper = .false.
                if (.not. all(is_close(hess_x_c(:n_param), 4.0_c_rp * x_c(:n_param)))) &
                    then
                    write(stderr, *) "test_hess_x_c_wrapper failed: Result not "// &
                        "returned for the "//trim(level_names(i_level))// &
                        " of an object "//trim(role_case_names(i_case))//"."
                    test_hess_x_c_wrapper = .false.
                end if
                if (.not. is_close(hess_x_c(n_param + 1), -1.0_c_rp)) then
                    write(stderr, *) "test_hess_x_c_wrapper failed: Array written "// &
                        "beyond the number of parameters for the "// &
                        trim(level_names(i_level))//" of an object "// &
                        trim(role_case_names(i_case))//"."
                    test_hess_x_c_wrapper = .false.
                end if
            end do
        end do

        ! call the Hessian linear transformation without a holder, which has to produce
        ! an error
        error_c = hess_x_c_wrapper(x_c, hess_x_c, c_null_ptr)
        if (error_c /= 1) then
            write(stderr, *) "test_hess_x_c_wrapper failed: No error produced "// &
                "without a holder."
            test_hess_x_c_wrapper = .false.
        end if

    end function test_hess_x_c_wrapper

    logical(c_bool) function test_precond_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the level-shifted preconditioner, which
        ! calls the callback function with the level of an extension object held at the
        ! address it receives as context, or the host's C function for an object that
        ! does not perform it
        !
        use c_interface_unit_tests, only: mock_precond
        use otr_common_c_interface, only: extension_holder_type, precond_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type, level_names

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: residual_c(n_param + 1), precond_residual_c(n_param + 1)
        integer(c_ip) :: error_c
        integer(ip) :: i_case, i_level

        ! assume tests pass
        test_precond_c_wrapper = .true.

        ! set up the extension and arrays one element longer than the number of
        ! parameters, of which only these may be read and written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        extension%mock_error = 3
        call random_number(residual_c)

        ! call the level-shifted preconditioner of an object that performs it and of
        ! one that does not, whose holder bundles the host's, with the addresses of the
        ! levels of the solver and of the stability check in the holder and determine
        ! if it calls the object's with the level or the host's with the host context
        ! of the level and returns its result and error
        do i_case = 1, 2
            extension%provides%precond = i_case == 1
            do i_level = 1, 2
                call setup_mock_extension_holder(holder, extension, i_level)
                holder%solver_level%callbacks%precond => mock_precond
                holder%stability_level%callbacks%precond => mock_precond
                precond_residual_c = -1.0_c_rp
                error_c = precond_c_wrapper(residual_c, 2.0_c_rp, precond_residual_c, &
                                            holder_level_context(holder, i_level))
                if (.not. check_role_c_wrapper("precond_c_wrapper", extension, &
                                               "precond", i_case, i_level, error_c)) &
                    test_precond_c_wrapper = .false.
                if (.not. all(is_close(precond_residual_c(:n_param), &
                                       2.0_c_rp * residual_c(:n_param)))) then
                    write(stderr, *) "test_precond_c_wrapper failed: Result not "// &
                        "returned for the "//trim(level_names(i_level))// &
                        " of an object "//trim(role_case_names(i_case))//"."
                    test_precond_c_wrapper = .false.
                end if
                if (.not. is_close(precond_residual_c(n_param + 1), -1.0_c_rp)) then
                    write(stderr, *) "test_precond_c_wrapper failed: Array written "// &
                        "beyond the number of parameters for the "// &
                        trim(level_names(i_level))//" of an object "// &
                        trim(role_case_names(i_case))//"."
                    test_precond_c_wrapper = .false.
                end if
            end do
        end do

        ! call the level-shifted preconditioner without a holder, which has to produce
        ! an error
        error_c = &
            precond_c_wrapper(residual_c, 2.0_c_rp, precond_residual_c, c_null_ptr)
        if (error_c /= 1) then
            write(stderr, *) "test_precond_c_wrapper failed: No error produced "// &
                "without a holder."
            test_precond_c_wrapper = .false.
        end if

    end function test_precond_c_wrapper

    logical(c_bool) function test_precond_pd_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the positive-definite preconditioner,
        ! which calls the callback function with the level of an extension object held
        ! at the address it receives as context, or the host's C function for an object
        ! that does not perform it
        !
        use c_interface_unit_tests, only: mock_precond_pd
        use otr_common_c_interface, only: extension_holder_type, precond_pd_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type, level_names

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: residual_c(n_param + 1), precond_residual_c(n_param + 1)
        integer(c_ip) :: error_c
        integer(ip) :: i_case, i_level

        ! assume tests pass
        test_precond_pd_c_wrapper = .true.

        ! set up the extension and arrays one element longer than the number of
        ! parameters, of which only these may be read and written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        extension%mock_error = 3
        call random_number(residual_c)

        ! call the positive-definite preconditioner of an object that performs it and
        ! of one that does not, whose holder bundles the host's, with the addresses of
        ! the levels of the solver and of the stability check in the holder and
        ! determine if it calls the object's with the level or the host's with the host
        ! context of the level and returns its result and error
        do i_case = 1, 2
            extension%provides%precond_pd = i_case == 1
            do i_level = 1, 2
                call setup_mock_extension_holder(holder, extension, i_level)
                holder%solver_level%callbacks%precond_pd => mock_precond_pd
                holder%stability_level%callbacks%precond_pd => mock_precond_pd
                precond_residual_c = -1.0_c_rp
                error_c = precond_pd_c_wrapper(residual_c, precond_residual_c, &
                                               holder_level_context(holder, i_level))
                if (.not. check_role_c_wrapper( &
                    "precond_pd_c_wrapper", extension, "precond_pd", i_case, i_level, &
                    error_c)) test_precond_pd_c_wrapper = .false.
                if (.not. all(is_close(precond_residual_c(:n_param), &
                                       3.0_c_rp * residual_c(:n_param)))) then
                    write(stderr, *) "test_precond_pd_c_wrapper failed: Result not "// &
                        "returned for the "//trim(level_names(i_level))// &
                        " of an object "//trim(role_case_names(i_case))//"."
                    test_precond_pd_c_wrapper = .false.
                end if
                if (.not. is_close(precond_residual_c(n_param + 1), -1.0_c_rp)) then
                    write(stderr, *) "test_precond_pd_c_wrapper failed: Array "// &
                        "written beyond the number of parameters for the "// &
                        trim(level_names(i_level))//" of an object "// &
                        trim(role_case_names(i_case))//"."
                    test_precond_pd_c_wrapper = .false.
                end if
            end do
        end do

        ! call the positive-definite preconditioner without a holder, which has to
        ! produce an error
        error_c = precond_pd_c_wrapper(residual_c, precond_residual_c, c_null_ptr)
        if (error_c /= 1) then
            write(stderr, *) "test_precond_pd_c_wrapper failed: No error produced "// &
                "without a holder."
            test_precond_pd_c_wrapper = .false.
        end if

    end function test_precond_pd_c_wrapper

    logical(c_bool) function test_project_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the projection, which calls the callback
        ! function with the level of an extension object held at the address it receives
        ! as context, or the host's C function for an object that does not perform it
        !
        use c_interface_unit_tests, only: mock_project
        use otr_common_c_interface, only: extension_holder_type, project_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type, level_names

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: vector_c(n_param + 1), projected_c(n_param + 1)
        integer(c_ip) :: error_c
        integer(ip) :: i_case, i_level

        ! assume tests pass
        test_project_c_wrapper = .true.

        ! set up the extension and arrays one element longer than the number of
        ! parameters, of which only these may be read and written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        extension%mock_error = 3
        call random_number(vector_c)

        ! call the projection of an object that performs it and of one that does not,
        ! whose holder bundles the host's, with the addresses of the levels of the
        ! solver and of the stability check in the holder and determine if it calls the
        ! object's with the level or the host's with the host context of the level and
        ! returns its result and error
        do i_case = 1, 2
            extension%provides%project = i_case == 1
            do i_level = 1, 2
                call setup_mock_extension_holder(holder, extension, i_level)
                holder%solver_level%callbacks%project => mock_project
                holder%stability_level%callbacks%project => mock_project
                projected_c = vector_c
                error_c = project_c_wrapper(projected_c, &
                                            holder_level_context(holder, i_level))
                if (.not. check_role_c_wrapper("project_c_wrapper", extension, &
                                               "project", i_case, i_level, error_c)) &
                    test_project_c_wrapper = .false.
                if (.not. all(is_close(projected_c(:n_param), &
                                       2.0_c_rp * vector_c(:n_param)))) then
                    write(stderr, *) "test_project_c_wrapper failed: Result not "// &
                        "returned for the "//trim(level_names(i_level))// &
                        " of an object "//trim(role_case_names(i_case))//"."
                    test_project_c_wrapper = .false.
                end if
                if (.not. is_close(projected_c(n_param + 1), vector_c(n_param + 1))) &
                    then
                    write(stderr, *) "test_project_c_wrapper failed: Array written "// &
                        "beyond the number of parameters for the "// &
                        trim(level_names(i_level))//" of an object "// &
                        trim(role_case_names(i_case))//"."
                    test_project_c_wrapper = .false.
                end if
            end do
        end do

        ! call the projection without a holder, which has to produce an error
        error_c = project_c_wrapper(projected_c, c_null_ptr)
        if (error_c /= 1) then
            write(stderr, *) "test_project_c_wrapper failed: No error produced "// &
                "without a holder."
            test_project_c_wrapper = .false.
        end if

    end function test_project_c_wrapper

    logical(c_bool) function test_get_extra_trial_vectors_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the extra trial vector function, which
        ! calls the callback function with the level of an extension object held at the
        ! address it receives as context, or the host's C function for an object that
        ! does not perform it
        !
        use c_interface_unit_tests, only: mock_get_extra_trial_vectors
        use otr_common_c_interface, only: extension_holder_type, &
                                          get_extra_trial_vectors_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type, level_names

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(c_rp) :: trial_vectors_c(2 * n_param + 1)
        integer(c_ip) :: error_c
        integer(ip) :: i_case, i_level

        ! assume tests pass
        test_get_extra_trial_vectors_c_wrapper = .true.

        ! set up the extension and arrays one element longer than the number of
        ! parameters, of which only these may be read and written
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        extension%mock_error = 3

        ! call the extra trial vector function of an object that performs it and of one
        ! that does not, whose holder bundles the host's, with the addresses of the
        ! levels of the solver and of the stability check in the holder and determine
        ! if it calls the object's with the level or the host's with the host context
        ! of the level and returns its result and error
        do i_case = 1, 2
            extension%provides%get_extra_trial_vectors = i_case == 1
            do i_level = 1, 2
                call setup_mock_extension_holder(holder, extension, i_level)
                holder%solver_level%callbacks%get_extra_trial_vectors => &
                    mock_get_extra_trial_vectors
                holder%stability_level%callbacks%get_extra_trial_vectors => &
                    mock_get_extra_trial_vectors
                trial_vectors_c = -1.0_c_rp
                error_c = get_extra_trial_vectors_c_wrapper( &
                    trial_vectors_c, 2_c_ip, holder_level_context(holder, i_level))
                if (.not. check_role_c_wrapper("get_extra_trial_vectors_c_wrapper", &
                                               extension, "get_extra_trial_vectors", &
                                               i_case, i_level, error_c)) &
                    test_get_extra_trial_vectors_c_wrapper = .false.
                if (.not. all(is_close(trial_vectors_c(n_param + 1:2 * n_param), &
                                       2.0_c_rp))) then
                    write(stderr, *) "test_get_extra_trial_vectors_c_wrapper "// &
                        "failed: Result not returned for the "// &
                        trim(level_names(i_level))//" of an object "// &
                        trim(role_case_names(i_case))//"."
                    test_get_extra_trial_vectors_c_wrapper = .false.
                end if
                if (.not. is_close(trial_vectors_c(2 * n_param + 1), -1.0_c_rp)) then
                    write(stderr, *) "test_get_extra_trial_vectors_c_wrapper "// &
                        "failed: Array written beyond the number of parameters for "// &
                        "the "//trim(level_names(i_level))//" of an object "// &
                        trim(role_case_names(i_case))//"."
                    test_get_extra_trial_vectors_c_wrapper = .false.
                end if
            end do
        end do

        ! call the extra trial vector function without a holder, which has to produce
        ! an error
        error_c = get_extra_trial_vectors_c_wrapper(trial_vectors_c, 2_c_ip, c_null_ptr)
        if (error_c /= 1) then
            write(stderr, *) "test_get_extra_trial_vectors_c_wrapper failed: No "// &
                "error produced without a holder."
            test_get_extra_trial_vectors_c_wrapper = .false.
        end if

    end function test_get_extra_trial_vectors_c_wrapper

    logical(c_bool) function test_conv_check_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the convergence check of the solver or
        ! of its stability check, which forwards to the host's C function in the bundle
        ! of the level of a holder at the address it receives as context
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
        if (error_c /= 1) then
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
            if (.not. all(is_close(hess_x_c, 5.0_c_rp * x_c))) then
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
        if (error_c /= 1) then
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
            if (.not. all(is_close(hess_x_c, 6.0_c_rp * x_c))) then
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
        if (error_c /= 1) then
            write(stderr, *) "test_approx_hess_x_c_wrapper failed: No error "// &
                "produced without a holder."
            test_approx_hess_x_c_wrapper = .false.
        end if

    end function test_approx_hess_x_c_wrapper

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

    logical(c_bool) function test_evaluate_dm_cs_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating
        ! function of an orbital basis for the closed-shell case, which hands the C
        ! function in the bundle it receives as context the density matrix as that of a
        ! single particle channel
        !
        use otr_common, only: evaluate_dm_cs_type, get_response_cs_type
        use otr_common_c_interface, only: orbital_c_callbacks_type, &
                                          evaluate_dm_cs_f_wrapper
        use otr_common_test_reference, only: check_evaluate_dm_cs_funptr, n_ao
        use otr_common_unit_tests, only: mock_requests
        use test_reference, only: arm_host_context_c, host_context_reached, host_context

        type(orbital_c_callbacks_type), target :: callbacks
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(get_response_cs_type), pointer :: get_response_cs_funptr
        class(*), pointer :: context
        real(rp) :: dm(n_ao, n_ao), energy
        integer(ip) :: request, error

        ! assume tests pass
        test_evaluate_dm_cs_f_wrapper = .true.

        ! inject mock C function into a bundle, which the wrapper receives as context
        callbacks%evaluate_dm => mock_evaluate_dm_cs
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        if (.not. check_evaluate_dm_cs_funptr( &
            evaluate_dm_cs_funptr, "evaluate_dm_cs_f_wrapper", "", context)) &
            test_evaluate_dm_cs_f_wrapper = .false.
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request=0, 3)])) then
            write(stderr, *) "test_evaluate_dm_cs_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function wrong."
            test_evaluate_dm_cs_f_wrapper = .false.
        end if

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("evaluate_dm_cs_f_wrapper")) &
            test_evaluate_dm_cs_f_wrapper = .false.

        ! a C function that fails has to have its error returned and neither return nor
        ! keep the response function it produced
        call arm_host_context_c(callbacks%host_context)
        host_context%mock_error = 3
        callbacks%get_response => null()
        dm = 0.0_rp
        call evaluate_dm_cs_f_wrapper(dm, energy, &
                                      get_response_cs_funptr=get_response_cs_funptr, &
                                      error=error, context=context)
        if (error /= 3) then
            write(stderr, *) "test_evaluate_dm_cs_f_wrapper failed: Error of the C "// &
                "function not returned."
            test_evaluate_dm_cs_f_wrapper = .false.
        end if
        if (associated(get_response_cs_funptr) .or. &
            associated(callbacks%get_response)) then
            write(stderr, *) "test_evaluate_dm_cs_f_wrapper failed: Response "// &
                "function returned or kept after an error."
            test_evaluate_dm_cs_f_wrapper = .false.
        end if

    end function test_evaluate_dm_cs_f_wrapper

    logical(c_bool) function test_evaluate_dm_os_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating
        ! function of an orbital basis for the open-shell case, which finds the C
        ! function in the bundle it receives as context and keeps the response C
        ! function it returns there
        !
        use otr_common, only: evaluate_dm_os_type, get_response_os_type
        use otr_common_c_interface, only: orbital_c_callbacks_type, &
                                          evaluate_dm_os_f_wrapper
        use otr_common_test_reference, only: check_evaluate_dm_os_funptr, n_ao, &
                                             n_particle
        use otr_common_unit_tests, only: mock_requests
        use test_reference, only: arm_host_context_c, host_context_reached, host_context

        type(orbital_c_callbacks_type), target :: callbacks
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(get_response_os_type), pointer :: get_response_os_funptr
        class(*), pointer :: context
        real(rp) :: dm(n_ao, n_ao, n_particle), energy
        integer(ip) :: request, error

        ! assume tests pass
        test_evaluate_dm_os_f_wrapper = .true.

        ! inject mock C function into a bundle, which the wrapper receives as context
        callbacks%evaluate_dm => mock_evaluate_dm_os
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        if (.not. check_evaluate_dm_os_funptr( &
            evaluate_dm_os_funptr, "evaluate_dm_os_f_wrapper", "", context)) &
            test_evaluate_dm_os_f_wrapper = .false.
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request=0, 3)])) then
            write(stderr, *) "test_evaluate_dm_os_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function wrong."
            test_evaluate_dm_os_f_wrapper = .false.
        end if

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("evaluate_dm_os_f_wrapper")) &
            test_evaluate_dm_os_f_wrapper = .false.

        ! a C function that fails has to have its error returned and neither return nor
        ! keep the response function it produced
        call arm_host_context_c(callbacks%host_context)
        host_context%mock_error = 3
        callbacks%get_response => null()
        dm = 0.0_rp
        call evaluate_dm_os_f_wrapper(dm, energy, &
                                      get_response_funptr=get_response_os_funptr, &
                                      error=error, context=context)
        if (error /= 3) then
            write(stderr, *) "test_evaluate_dm_os_f_wrapper failed: Error of the C "// &
                "function not returned."
            test_evaluate_dm_os_f_wrapper = .false.
        end if
        if (associated(get_response_os_funptr) .or. &
            associated(callbacks%get_response)) then
            write(stderr, *) "test_evaluate_dm_os_f_wrapper failed: Response "// &
                "function returned or kept after an error."
            test_evaluate_dm_os_f_wrapper = .false.
        end if

        ! a context that is no bundle has to produce an error
        context => null()
        call evaluate_dm_os_f_wrapper(dm, energy, error=error, context=context)
        if (error /= 1) then
            write(stderr, *) "test_evaluate_dm_os_f_wrapper failed: No error "// &
                "produced without a bundle."
            test_evaluate_dm_os_f_wrapper = .false.
        end if

    end function test_evaluate_dm_os_f_wrapper

    logical(c_bool) function test_get_response_cs_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the response function of an
        ! orbital basis for the closed-shell case, which hands the C function the
        ! density matrix evaluating function returned last the density matrix as that
        ! of a single particle channel
        !
        use otr_common, only: get_response_cs_type
        use otr_common_c_interface, only: orbital_c_callbacks_type, &
                                          get_response_cs_f_wrapper
        use otr_common_test_reference, only: check_get_response_cs_funptr, n_ao
        use test_reference, only: arm_host_context_c, host_context_reached, host_context

        type(orbital_c_callbacks_type), target :: callbacks
        procedure(get_response_cs_type), pointer :: get_response_cs_funptr
        class(*), pointer :: context
        real(rp) :: dm(n_ao, n_ao), response(n_ao, n_ao)
        integer(ip) :: error

        ! assume tests pass
        test_get_response_cs_f_wrapper = .true.

        ! inject mock C function into a bundle, which the wrapper receives as context
        callbacks%get_response => mock_get_response_cs
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        get_response_cs_funptr => get_response_cs_f_wrapper

        ! test response function wrapper
        if (.not. check_get_response_cs_funptr( &
            get_response_cs_funptr, "get_response_cs_f_wrapper", "", context)) &
            test_get_response_cs_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("get_response_cs_f_wrapper")) &
            test_get_response_cs_f_wrapper = .false.

        ! a C function that fails has to have its error returned
        call arm_host_context_c(callbacks%host_context)
        host_context%mock_error = 3
        dm = 0.0_rp
        call get_response_cs_f_wrapper(dm, response, error, context)
        if (error /= 3) then
            write(stderr, *) "test_get_response_cs_f_wrapper failed: Error of the "// &
                "C function not returned."
            test_get_response_cs_f_wrapper = .false.
        end if

    end function test_get_response_cs_f_wrapper

    logical(c_bool) function test_get_response_os_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the response function of an
        ! orbital basis for the open-shell case, which finds the C function the density
        ! matrix evaluating function returned last in the bundle it receives as context
        !
        use otr_common, only: get_response_os_type
        use otr_common_c_interface, only: orbital_c_callbacks_type, &
                                          get_response_os_f_wrapper
        use otr_common_test_reference, only: check_get_response_os_funptr, n_ao, &
                                             n_particle
        use test_reference, only: arm_host_context_c, host_context_reached, host_context

        type(orbital_c_callbacks_type), target :: callbacks
        procedure(get_response_os_type), pointer :: get_response_os_funptr
        class(*), pointer :: context
        real(rp) :: dm(n_ao, n_ao, n_particle), response(n_ao, n_ao, n_particle)
        integer(ip) :: error

        ! assume tests pass
        test_get_response_os_f_wrapper = .true.

        ! inject mock C function into a bundle, which the wrapper receives as context
        callbacks%get_response => mock_get_response_os
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        get_response_os_funptr => get_response_os_f_wrapper

        ! test response function wrapper
        if (.not. check_get_response_os_funptr( &
            get_response_os_funptr, "get_response_os_f_wrapper", "", context)) &
            test_get_response_os_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("get_response_os_f_wrapper")) &
            test_get_response_os_f_wrapper = .false.

        ! a C function that fails has to have its error returned
        call arm_host_context_c(callbacks%host_context)
        host_context%mock_error = 3
        dm = 0.0_rp
        call get_response_os_f_wrapper(dm, response, error, context)
        if (error /= 3) then
            write(stderr, *) "test_get_response_os_f_wrapper failed: Error of the "// &
                "C function not returned."
            test_get_response_os_f_wrapper = .false.
        end if

        ! a context that is no bundle has to produce an error
        context => null()
        call get_response_os_f_wrapper(dm, response, error, context)
        if (error /= 1) then
            write(stderr, *) "test_get_response_os_f_wrapper failed: No error "// &
                "produced without a bundle."
            test_get_response_os_f_wrapper = .false.
        end if

    end function test_get_response_os_f_wrapper

end module otr_common_c_interface_unit_tests
