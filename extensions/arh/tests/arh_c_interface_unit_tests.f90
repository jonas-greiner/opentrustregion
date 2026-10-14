! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use otr_arh_c_interface, only: evaluate_dm_c_type
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_ptr, c_funptr, &
                                           c_funloc, c_loc, c_f_pointer, c_null_ptr

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(evaluate_dm_c_type), pointer :: mock_arh_evaluate_dm_cs_ptr => &
        mock_arh_evaluate_dm_cs
    procedure(evaluate_dm_c_type), pointer :: mock_arh_evaluate_dm_os_ptr => &
        mock_arh_evaluate_dm_os

contains

    function mock_arh_evaluate_dm_cs(dm_ao, energy, fock, v_coulomb, v_exchange, &
                                     v_nonlinear, context_c) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function with separate Coulomb, exact-exchange and non-linear potential
        ! contributions for the closed-shell case
        !
        use otr_common_test_reference, only: n_ao
        use otr_arh_test_reference, only: arh_evaluate_dm_factors
        use otr_common_unit_tests, only: record_mock_call
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*), v_coulomb(*), v_exchange(*), &
                                             v_nonlinear(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2

        ! check host context
        call check_host_context_c(context_c)

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + merge( &
            2_ip, 0_ip, present(v_coulomb)) + merge(4_ip, 0_ip, present(v_exchange)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) &
            fock(:flat_len) = arh_evaluate_dm_factors(1) * dm_ao(:flat_len)
        if (present(v_coulomb)) &
            v_coulomb(:flat_len) = arh_evaluate_dm_factors(2) * dm_ao(:flat_len)
        if (present(v_exchange)) &
            v_exchange(:flat_len) = arh_evaluate_dm_factors(3) * dm_ao(:flat_len)
        if (present(v_nonlinear)) &
            v_nonlinear(:flat_len) = arh_evaluate_dm_factors(4) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_arh_evaluate_dm_cs

    function mock_arh_evaluate_dm_os(dm_ao, energy, fock, v_coulomb, v_exchange, &
                                     v_nonlinear, context_c) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function with separate Coulomb, exact-exchange and non-linear potential
        ! contributions for the open-shell case
        !
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_arh_test_reference, only: arh_evaluate_dm_factors
        use otr_common_unit_tests, only: record_mock_call
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*), v_coulomb(*), v_exchange(*), &
                                             v_nonlinear(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2 * n_particle

        ! check host context
        call check_host_context_c(context_c)

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + merge( &
            2_ip, 0_ip, present(v_coulomb)) + merge(4_ip, 0_ip, present(v_exchange)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) &
            fock(:flat_len) = arh_evaluate_dm_factors(1) * dm_ao(:flat_len)
        if (present(v_coulomb)) &
            v_coulomb(:flat_len) = arh_evaluate_dm_factors(2) * dm_ao(:flat_len)
        if (present(v_exchange)) &
            v_exchange(:flat_len) = arh_evaluate_dm_factors(3) * dm_ao(:flat_len)
        if (present(v_nonlinear)) &
            v_nonlinear(:flat_len) = arh_evaluate_dm_factors(4) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_arh_evaluate_dm_os

    logical(c_bool) function test_arh_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH factory on a closed-shell MO
        ! object and an open-shell OAO object, held by holders whose level of the
        ! solver is the context of the C solver settings
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, solver_state_type_c, &
                               assignment(=)
        use otr_common, only: extension_type, extension_level_type, orbital_basis_type
        use otr_common_unit_tests, only: wired_level
        use otr_common_c_interface, only: extension_holder_type
        use otr_mo, only: mo_type
        use otr_oao, only: oao_type
        use otr_arh_c_interface, only: arh_settings_type_c, arh_c_callbacks_type, &
                                       arh_factory_common, arh_factory_c_wrapper
        use otr_arh_mock, only: mock_arh_factory_common, mock_arh_factory_error, &
                                mock_arh_factory_defaults, received_arh_level, &
                                test_passed
        use otr_arh_test_reference, only: assignment(=), ref_arh_settings
        use otr_common_c_interface_unit_tests, only: check_c_factory_wiring
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(extension_holder_type), target :: orbitals_holder
        class(extension_type), allocatable, target :: orbitals
        type(c_funptr) :: evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr
        type(arh_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        type(solver_state_type_c), target :: state_c
        type(c_ptr) :: extension_c
        type(extension_holder_type), pointer :: holder
        type(extension_level_type), pointer :: level
        real(rp), target :: copy(1, 1, 1)
        real(c_rp), target :: host_array(1, 1, 1)
        integer(ip) :: n_particle, i_case
        integer(c_ip) :: error_c
        character(len=:), allocatable :: case_name
        character(len=34), parameter :: default_cases(2) = &
            [character(len=34) :: "settings that were not initialized", "no settings"]

        ! assume tests pass
        test_arh_factory_c_wrapper = .true.

        ! inject mock function
        arh_factory_common => mock_arh_factory_common

        ! associate optional settings with values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! the factory builds on a closed-shell MO object and on an open-shell OAO
        ! object, whose holders hold the copy of the orbitals the orbital update copies
        ! back
        do n_particle = 1, 2
            if (n_particle == 1) then
                case_name = "closed-shell MO"
                allocate(mo_type :: orbitals)
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_cs)
            else
                case_name = "open-shell OAO"
                allocate(oao_type :: orbitals)
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_os)
            end if
            select type (orbitals)
            class is (orbital_basis_type)
                orbitals%n_particle = n_particle
            end select
            orbitals%solver_level%extension => orbitals
            orbitals_holder%solver_level%level => orbitals%solver_level
            orbitals_holder%solver_level%orbitals => copy
            orbitals_holder%solver_level%orbitals_c => host_array

            ! arm the host context of the ARH settings, which the logging and density
            ! matrix evaluating functions receive
            call arm_host_context_c(settings_c%context)
            solver_settings_c = default_solver_settings
            solver_settings_c%context = c_loc(orbitals_holder%solver_level)

            ! call ARH factory C wrapper
            test_passed = .false.
            error_c = arh_factory_c_wrapper( &
                evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr, &
                extension_c, solver_settings_c, c_loc(state_c), settings_c)

            ! check if the mock factory received the correct input, which includes the
            ! level of the orbital basis object, and called the logging function
            if (.not. test_passed .or. &
                .not. associated(received_arh_level, orbitals%solver_level)) then
                test_arh_factory_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_c_wrapper failed: Mock factory "// &
                    "received wrong input for the "//case_name//" case."
            end if
            if (.not. host_context%logger_called) then
                test_arh_factory_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_c_wrapper failed: Called "// &
                    "logging subroutine wrong for the "//case_name//" case."
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("arh_factory_c_wrapper")) &
                test_arh_factory_c_wrapper = .false.

            ! check if output variables are as expected, the returned extension being
            ! the holder, whose inspection would crash without one
            if (error_c /= 0 .or. .not. c_associated(extension_c)) then
                write(stderr, *) "test_arh_factory_c_wrapper failed: Error "// &
                    "returned or no extension returned for the "//case_name//" case."
                test_arh_factory_c_wrapper = .false.
                return
            end if
            call c_f_pointer(extension_c, holder)

            ! check that the settings of the holder, which own the object, were wired
            ! with it, that it reads the state of the holder and that the holder holds
            ! the bundle of the density matrix evaluating function, the host's record
            ! of the state and the copy of the orbitals of the orbital basis object
            level => wired_level(holder%solver_settings)
            if (.not. associated(level)) then
                write(stderr, *) "test_arh_factory_c_wrapper failed: Settings of "// &
                    "the holder not wired with an object for the "//case_name//" case."
                test_arh_factory_c_wrapper = .false.
                return
            end if
            if (.not. associated(level%extension%state, holder%state)) then
                write(stderr, *) "test_arh_factory_c_wrapper failed: Object does "// &
                    "not read the state of the holder for the "//case_name//" case."
                test_arh_factory_c_wrapper = .false.
            end if
            select type (callbacks => holder%extension_callbacks)
            type is (arh_c_callbacks_type)
                if (.not. ( &
                    associated(callbacks%evaluate_dm, mock_arh_evaluate_dm_cs) .or. &
                    associated(callbacks%evaluate_dm, mock_arh_evaluate_dm_os))) then
                    write(stderr, *) "test_arh_factory_c_wrapper failed: Density "// &
                        "matrix evaluating function not held for the "//case_name// &
                        " case."
                    test_arh_factory_c_wrapper = .false.
                end if
            class default
                write(stderr, *) "test_arh_factory_c_wrapper failed: Bundle of the "// &
                    "callback functions not held for the "//case_name//" case."
                test_arh_factory_c_wrapper = .false.
            end select
            if (.not. (associated(holder%solver_level%state_c, state_c) .and. &
                       associated(holder%stability_level%state_c, state_c))) then
                write(stderr, *) "test_arh_factory_c_wrapper failed: Host's record "// &
                    "of the state not held for the "//case_name//" case."
                test_arh_factory_c_wrapper = .false.
            end if
            if (.not. (associated(holder%solver_level%orbitals, copy) .and. &
                       associated(holder%solver_level%orbitals_c, host_array))) then
                write(stderr, *) "test_arh_factory_c_wrapper failed: Copy of the "// &
                    "orbitals of the orbital basis object not held for the "// &
                    case_name//" case."
                test_arh_factory_c_wrapper = .false.
            end if

            ! check that the factory returned the C wrappers of the object and wired the
            ! solver settings with its holder
            if (.not. check_c_factory_wiring("arh_factory_c_wrapper", holder, &
                                             level%extension, update_orbs_c_funptr, &
                                             obj_func_c_funptr, solver_settings_c)) &
                test_arh_factory_c_wrapper = .false.
            if (c_associated(extension_c)) then
                call c_f_pointer(extension_c, holder)
                deallocate(holder)
            end if

            ! call the factory with a mock factory which produces an error, which has
            ! to be returned without functions, an extension or writing back the solver
            ! settings
            solver_settings_c = default_solver_settings
            solver_settings_c%context = c_loc(orbitals_holder%solver_level)
            update_orbs_c_funptr = c_funloc(mock_logger_c)
            obj_func_c_funptr = c_funloc(mock_logger_c)
            mock_arh_factory_error = 1
            error_c = arh_factory_c_wrapper( &
                evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr, &
                extension_c, solver_settings_c, c_loc(state_c), settings_c)
            mock_arh_factory_error = 0
            if (error_c /= 1 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. c_associated(extension_c) .or. &
                .not. c_associated(solver_settings_c%context, &
                                   c_loc(orbitals_holder%solver_level))) then
                test_arh_factory_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_c_wrapper failed: Error of the "// &
                    "factory not passed on for the "//case_name//" case."
            end if
            if (n_particle == 1) deallocate(orbitals)
        end do

        ! call the factory with settings that were not initialized and without any, for
        ! which the factory is handed the default settings with the context through
        ! which the density matrix evaluating function is found, which the mock checks,
        ! and the bundle keeps no logging function or host context
        do i_case = 1, size(default_cases)
            solver_settings_c = default_solver_settings
            solver_settings_c%context = c_loc(orbitals_holder%solver_level)
            mock_arh_factory_defaults = .true.
            test_passed = .false.
            if (i_case == 1) then
                settings_c%initialized = .false.
                error_c = arh_factory_c_wrapper( &
                    evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr, &
                    extension_c, solver_settings_c, c_loc(state_c), settings_c)
                settings_c%initialized = .true.
            else
                error_c = arh_factory_c_wrapper( &
                    evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr, &
                    extension_c, solver_settings_c, c_loc(state_c))
            end if
            mock_arh_factory_defaults = .false.
            if (error_c /= 0 .or. .not. test_passed) then
                write(stderr, *) "test_arh_factory_c_wrapper failed: Default "// &
                    "settings not handed to the factory for "// &
                    trim(default_cases(i_case))//"."
                test_arh_factory_c_wrapper = .false.
            else
                call c_f_pointer(extension_c, holder)
                select type (callbacks => holder%extension_callbacks)
                type is (arh_c_callbacks_type)
                    if (associated(callbacks%logger) .or. &
                        c_associated(callbacks%host_context)) then
                        write(stderr, *) "test_arh_factory_c_wrapper failed: "// &
                            "Logging function or host context kept for "// &
                            trim(default_cases(i_case))//"."
                        test_arh_factory_c_wrapper = .false.
                    end if
                end select
            end if
            if (c_associated(extension_c)) then
                call c_f_pointer(extension_c, holder)
                deallocate(holder)
            end if
        end do

        ! call the factory for C solver settings without a context, which no factory
        ! wired, for which the factory has to be handed no level
        solver_settings_c = default_solver_settings
        error_c = arh_factory_c_wrapper(evaluate_dm_c_funptr, update_orbs_c_funptr, &
                                        obj_func_c_funptr, extension_c, &
                                        solver_settings_c, c_loc(state_c), settings_c)
        if (associated(received_arh_level)) then
            write(stderr, *) "test_arh_factory_c_wrapper failed: Level handed to "// &
                "the factory for C solver settings without a context."
            test_arh_factory_c_wrapper = .false.
        end if
        if (c_associated(extension_c)) then
            call c_f_pointer(extension_c, holder)
            deallocate(holder)
        end if

        ! call the factory without a record of the state, which has to produce an error
        ! without calling the factory or returning functions or an extension
        solver_settings_c = default_solver_settings
        solver_settings_c%context = c_loc(orbitals_holder%solver_level)
        test_passed = .false.
        error_c = arh_factory_c_wrapper(evaluate_dm_c_funptr, update_orbs_c_funptr, &
                                        obj_func_c_funptr, extension_c, &
                                        solver_settings_c, c_null_ptr, settings_c)
        if (error_c == 0 .or. test_passed .or. c_associated(update_orbs_c_funptr) .or. &
            c_associated(obj_func_c_funptr) .or. c_associated(extension_c)) then
            write(stderr, *) "test_arh_factory_c_wrapper failed: No error produced "// &
                "without a record of the state."
            test_arh_factory_c_wrapper = .false.
        end if
        deallocate(orbitals)

    end function test_arh_factory_c_wrapper

    logical(c_bool) function test_evaluate_dm_cs_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating C
        ! function with Coulomb, exact-exchange and non-linear potential contributions
        ! for the closed-shell case
        !
        use otr_arh, only: evaluate_dm_cs_type
        use otr_arh_c_interface, only: arh_c_callbacks_type, evaluate_dm_cs_f_wrapper
        use otr_arh_test_reference, only: test_evaluate_dm_cs_funptr
        use otr_common_unit_tests, only: mock_requests
        use test_reference, only: arm_host_context_c, host_context_reached

        type(arh_c_callbacks_type), target :: callbacks
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        class(*), pointer :: context
        integer(ip) :: request

        ! assume tests pass
        test_evaluate_dm_cs_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%evaluate_dm => mock_arh_evaluate_dm_cs
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        if (.not. test_evaluate_dm_cs_funptr(evaluate_dm_cs_funptr, &
                                             "evaluate_dm_cs_f_wrapper", "", context)) &
            test_evaluate_dm_cs_f_wrapper = .false.
        if (size(mock_requests) /= 16 .or. &
            .not. all([(any(mock_requests == request), request=0, 15)])) then
            write(stderr, *) "test_evaluate_dm_cs_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for closed-shell case "// &
                "wrong."
            test_evaluate_dm_cs_f_wrapper = .false.
        end if

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("evaluate_dm_cs_f_wrapper")) &
            test_evaluate_dm_cs_f_wrapper = .false.

    end function test_evaluate_dm_cs_f_wrapper

    logical(c_bool) function test_evaluate_dm_os_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating C
        ! function with Coulomb, exact-exchange and non-linear potential contributions
        ! for the open-shell case
        !
        use otr_arh, only: evaluate_dm_os_type
        use otr_arh_c_interface, only: arh_c_callbacks_type, evaluate_dm_os_f_wrapper
        use otr_arh_test_reference, only: test_evaluate_dm_os_funptr
        use otr_common_unit_tests, only: mock_requests
        use test_reference, only: arm_host_context_c, host_context_reached

        type(arh_c_callbacks_type), target :: callbacks
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        class(*), pointer :: context
        real(rp) :: dm(1, 1, 1), energy
        integer(ip) :: request, error

        ! assume tests pass
        test_evaluate_dm_os_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%evaluate_dm => mock_arh_evaluate_dm_os
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        if (.not. test_evaluate_dm_os_funptr(evaluate_dm_os_funptr, &
                                             "evaluate_dm_os_f_wrapper", "", context)) &
            test_evaluate_dm_os_f_wrapper = .false.
        if (size(mock_requests) /= 16 .or. &
            .not. all([(any(mock_requests == request), request=0, 15)])) then
            write(stderr, *) "test_evaluate_dm_os_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for open-shell case wrong."
            test_evaluate_dm_os_f_wrapper = .false.
        end if

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("evaluate_dm_os_f_wrapper")) &
            test_evaluate_dm_os_f_wrapper = .false.

        ! a context that is no bundle has to produce an error
        dm = 0.0_rp
        context => null()
        call evaluate_dm_os_f_wrapper(dm, energy, error=error, context=context)
        if (error == 0) then
            write(stderr, *) "test_evaluate_dm_os_f_wrapper failed: No error "// &
                "produced without a bundle."
            test_evaluate_dm_os_f_wrapper = .false.
        end if

    end function test_evaluate_dm_os_f_wrapper

    logical(c_bool) function test_init_arh_settings_c() bind(C)
        !
        ! this function tests that the ARH settings initialization routine correctly
        ! initializes all settings to their default values
        !
        use otr_arh_c_interface, only: arh_settings_type_c, init_arh_settings_c
        use otr_arh, only: default_arh_settings
        use otr_arh_test_reference, only: operator(/=)

        type(arh_settings_type_c) :: settings

        ! assume test passes
        test_init_arh_settings_c = .true.

        ! initialize settings
        call init_arh_settings_c(settings)

        ! check settings
        if (settings /= default_arh_settings) then
            write(stderr, *) "test_init_arh_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_arh_settings_c = .false.
        end if

    end function test_init_arh_settings_c

    logical(c_bool) function test_assign_arh_f_c() bind(C)
        !
        ! this function tests that the function that converts ARH settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh, only: arh_settings_type
        use otr_arh_test_reference, only: assignment(=), operator(/=), ref_arh_settings
        use c_interface, only: c_callbacks_type
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(arh_settings_type_c) :: settings_c
        type(arh_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_arh_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the factory does
        settings = settings_c
        callbacks%logger => mock_logger_c
        settings%context => callbacks

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_arh_f_c = .false.
            write(stderr, *) "test_assign_arh_f_c failed: Logging function not "// &
                "associated with value."
        else
            call arm_host_context_c(callbacks%host_context)
            call settings%logger("test", settings%context)
            if (.not. host_context%logger_called) then
                test_assign_arh_f_c = .false.
                write(stderr, *) "test_assign_arh_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
            if (.not. host_context_reached("assign_arh_f_c")) &
                test_assign_arh_f_c = .false.
        end if

        ! check against reference values
        if (settings /= ref_arh_settings) then
            write(stderr, *) "test_assign_arh_f_c failed: Settings not converted "// &
                "correctly."
            test_assign_arh_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            write(stderr, *) "test_assign_arh_f_c failed: Settings not marked as "// &
                "initialized."
            test_assign_arh_f_c = .false.
        end if

        ! convert C settings that were not initialized and check that they are left
        ! uninitialized without a logging function, so that the factory initializes
        ! them after its warning
        settings_c%initialized = .false.
        settings = settings_c
        if (settings%initialized .or. associated(settings%logger)) then
            write(stderr, *) "test_assign_arh_f_c failed: Settings that were not "// &
                "initialized converted."
            test_assign_arh_f_c = .false.
        end if

    end function test_assign_arh_f_c

    logical(c_bool) function test_assign_arh_c_f() bind(C)
        !
        ! this function tests that the function that converts ARH settings from Fortran
        ! to C correctly performs this conversion
        !
        use otr_arh, only: arh_settings_type
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh_test_reference, only: assignment(=), operator(/=), ref_arh_settings

        type(arh_settings_type) :: settings
        type(arh_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_arh_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_arh_settings

        ! convert to C settings
        settings_c = settings

        ! check that no callback function pointers and host contexts are associated
        if (c_associated(settings_c%logger)) then
            test_assign_arh_c_f = .false.
            write(stderr, *) "test_assign_arh_c_f failed: Logger function associated."
        end if
        if (c_associated(settings_c%context)) then
            write(stderr, *) "test_assign_arh_c_f failed: Host context associated."
            test_assign_arh_c_f = .false.
        end if

        ! check against reference values
        if (settings_c /= ref_arh_settings) then
            write(stderr, *) "test_assign_arh_c_f failed: Settings not converted "// &
                "correctly."
            test_assign_arh_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_arh_c_f = .false.
            write(stderr, *) "test_assign_arh_c_f failed: Settings not marked as "// &
                "initialized."
        end if

    end function test_assign_arh_c_f

end module otr_arh_c_interface_unit_tests
