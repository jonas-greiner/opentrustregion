! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_c_interface_unit_tests

    use opentrustregion_core, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: is_close
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_associated, c_loc, c_f_pointer, &
                                           c_null_ptr, c_null_funptr

    implicit none

contains

    logical(c_bool) function test_oao_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO factory for the closed- and the
        ! open-shell case
        !
        use opentrustregion_core, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, solver_state_type_c, &
                               assignment(=)
        use otr_common_c_interface, only: extension_holder_type, &
                                          orbital_c_callbacks_type
        use otr_oao, only: standard_oao_factory_cs => oao_factory_cs, &
                           standard_oao_factory_os => oao_factory_os, error_oao
        use otr_oao_c_interface, only: oao_settings_type_c, oao_factory_cs, &
                                       oao_factory_os, oao_factory_c_wrapper
        use otr_oao_mock, only: mock_oao_factory_cs, mock_oao_factory_os, &
                                mock_oao_factory_error, mock_oao_factory_defaults, &
                                mock_oao_factory_no_evaluate_dm, test_passed
        use otr_oao_test_reference, only: get_reference_oao_values, density_pattern
        use otr_common_test_reference, only: n_particle, n_ao_c, shell_names, &
                                             ao_overlap_pattern
        use otr_common_c_interface_unit_tests, only: &
            mock_evaluate_dm_cs, mock_evaluate_dm_os, check_c_factory_wiring
        use c_interface_unit_tests, only: mock_logger, mock_update_orbs, mock_obj_func
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        real(c_rp), allocatable :: ao_overlap_c(:, :), dm_ao_c(:, :, :)
        type(c_funptr) :: evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr
        type(oao_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        type(solver_state_type_c), target :: state_c
        type(c_ptr) :: extension_c
        type(extension_holder_type), pointer :: holder
        integer(c_ip) :: n_particle_c, error_c
        integer(ip) :: i_case
        character(len=:), allocatable :: case_name
        character(len=34), parameter :: default_cases(2) = &
            [character(len=34) :: "settings that were not initialized", "no settings"]

        ! assume tests pass
        test_oao_factory_c_wrapper = .true.

        ! inject mock functions
        oao_factory_cs => mock_oao_factory_cs
        oao_factory_os => mock_oao_factory_os

        ! initialize AO overlap matrix with values encoding the pairs of its indices
        ao_overlap_c = ao_overlap_pattern(n_ao_c)

        ! associate optional settings with the reference values
        call get_reference_oao_values(settings_c)
        settings_c%logger = c_funloc(mock_logger)

        ! both spin cases pass through the same C wrapper
        do n_particle_c = 1, int(n_particle, kind=c_ip)
            case_name = trim(shell_names(n_particle_c))

            ! initialize density matrix with values encoding the particle channel and
            ! the pairs of its indices and get C function pointer to the density matrix
            ! evaluating function of the spin case
            dm_ao_c = density_pattern(n_ao_c, n_particle_c, 0.0_c_rp)
            if (n_particle_c == 1) then
                evaluate_dm_c_funptr = c_funloc(mock_evaluate_dm_cs)
            else
                evaluate_dm_c_funptr = c_funloc(mock_evaluate_dm_os)
            end if

            ! arm the host context of the OAO settings, which the logging and density
            ! matrix evaluating functions receive
            call arm_host_context_c(settings_c%context)
            solver_settings_c = default_solver_settings

            ! call OAO factory C wrapper
            test_passed = .false.
            error_c = oao_factory_c_wrapper( &
                dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
                update_orbs_c_funptr, obj_func_c_funptr, extension_c, &
                solver_settings_c, c_loc(state_c), settings_c)

            ! check if the mock factory received the correct input and called the
            ! logging function
            if (.not. test_passed) then
                test_oao_factory_c_wrapper = .false.
                write(stderr, *) "test_oao_factory_c_wrapper failed: Mock factory "// &
                    "received wrong input for the "//case_name//" case."
            end if
            if (.not. host_context%logger_called) then
                test_oao_factory_c_wrapper = .false.
                write(stderr, *) "test_oao_factory_c_wrapper failed: Called "// &
                    "logging subroutine wrong for the "//case_name//" case."
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("oao_factory_c_wrapper")) &
                test_oao_factory_c_wrapper = .false.

            ! check if output variables are as expected
            if (error_c /= 0) then
                test_oao_factory_c_wrapper = .false.
                write(stderr, *) "test_oao_factory_c_wrapper failed: Returned "// &
                    "error code wrong for the "//case_name//" case."
            end if

            ! check that the factory received the density matrix of the caller in place
            if (.not. all(is_close(dm_ao_c, 2.0_c_rp))) then
                test_oao_factory_c_wrapper = .false.
                write(stderr, *) "test_oao_factory_c_wrapper failed: Density "// &
                    "matrix of the caller not passed in place for the "//case_name// &
                    " case."
            end if

            ! the returned extension is the holder, whose inspection would crash
            ! without one
            if (.not. c_associated(extension_c)) then
                write(stderr, *) "test_oao_factory_c_wrapper failed: No extension "// &
                    "returned for the "//case_name//" case."
                test_oao_factory_c_wrapper = .false.
                cycle
            end if
            call c_f_pointer(extension_c, holder)

            ! check that the holder holds the bundle of the density matrix evaluating
            ! function
            select type (callbacks => holder%extension_callbacks)
            type is (orbital_c_callbacks_type)
                if (n_particle_c == 1) then
                    if (.not. associated(callbacks%evaluate_dm, mock_evaluate_dm_cs)) &
                        then
                        write(stderr, *) "test_oao_factory_c_wrapper failed: "// &
                            "Density matrix evaluating function not held for the "// &
                            case_name//" case."
                        test_oao_factory_c_wrapper = .false.
                    end if
                else if (.not. associated(callbacks%evaluate_dm, mock_evaluate_dm_os)) &
                    then
                    write(stderr, *) "test_oao_factory_c_wrapper failed: Density "// &
                        "matrix evaluating function not held for the "//case_name// &
                        " case."
                    test_oao_factory_c_wrapper = .false.
                end if
            class default
                write(stderr, *) "test_oao_factory_c_wrapper failed: Bundle of the "// &
                    "callback functions not held for the "//case_name//" case."
                test_oao_factory_c_wrapper = .false.
            end select

            ! check that the factory wired the settings of the holder with an object,
            ! which reads the state of the holder, that the holder holds the host's
            ! record of the state, and that the factory returned the C wrappers of the
            ! object and wired the solver settings with its holder
            if (.not. check_c_factory_wiring("oao_factory_c_wrapper", holder, state_c, &
                                             update_orbs_c_funptr, obj_func_c_funptr, &
                                             solver_settings_c, case_name)) &
                test_oao_factory_c_wrapper = .false.
            if (c_associated(extension_c)) then
                call c_f_pointer(extension_c, holder)
                deallocate(holder)
            end if

            ! call the factory without a density matrix evaluating function, which the
            ! factory has to be handed as none and the bundle must not hold
            dm_ao_c = density_pattern(n_ao_c, n_particle_c, 0.0_c_rp)
            solver_settings_c = default_solver_settings
            mock_oao_factory_no_evaluate_dm = .true.
            test_passed = .false.
            error_c = oao_factory_c_wrapper( &
                dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, c_null_funptr, &
                update_orbs_c_funptr, obj_func_c_funptr, extension_c, &
                solver_settings_c, c_loc(state_c), settings_c)
            mock_oao_factory_no_evaluate_dm = .false.
            if (error_c /= 0 .or. .not. test_passed) then
                write(stderr, *) "test_oao_factory_c_wrapper failed: Missing "// &
                    "density matrix evaluating function not handed to the factory "// &
                    "for the "//case_name//" case."
                test_oao_factory_c_wrapper = .false.
            else
                call c_f_pointer(extension_c, holder)
                select type (callbacks => holder%extension_callbacks)
                type is (orbital_c_callbacks_type)
                    if (associated(callbacks%evaluate_dm)) then
                        write(stderr, *) "test_oao_factory_c_wrapper failed: "// &
                            "Density matrix evaluating function held without one "// &
                            "for the "//case_name//" case."
                        test_oao_factory_c_wrapper = .false.
                    end if
                end select
            end if
            if (c_associated(extension_c)) then
                call c_f_pointer(extension_c, holder)
                deallocate(holder)
            end if

            ! call the factory for the initial density matrix with a mock factory which
            ! produces an error, which has to be returned without functions, an
            ! extension or writing back the solver settings
            dm_ao_c = density_pattern(n_ao_c, n_particle_c, 0.0_c_rp)
            solver_settings_c = default_solver_settings
            update_orbs_c_funptr = c_funloc(mock_update_orbs)
            obj_func_c_funptr = c_funloc(mock_obj_func)
            extension_c = c_loc(state_c)
            mock_oao_factory_error = 1
            error_c = oao_factory_c_wrapper( &
                dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
                update_orbs_c_funptr, obj_func_c_funptr, extension_c, &
                solver_settings_c, c_loc(state_c), settings_c)
            mock_oao_factory_error = 0
            if (error_c /= 1 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. c_associated(extension_c) .or. &
                c_associated(solver_settings_c%context)) then
                test_oao_factory_c_wrapper = .false.
                write(stderr, *) "test_oao_factory_c_wrapper failed: Error of the "// &
                    "factory not passed on for the "//case_name//" case."
            end if
        end do

        ! call the factory with settings that were not initialized and without any, for
        ! which the factory is handed the default settings with the context through
        ! which the density matrix evaluating function is found, which the mock checks,
        ! and the bundle keeps no logging function or host context
        evaluate_dm_c_funptr = c_funloc(mock_evaluate_dm_cs)
        do i_case = 1, size(default_cases)
            dm_ao_c = density_pattern(n_ao_c, 1_c_ip, 0.0_c_rp)
            solver_settings_c = default_solver_settings
            mock_oao_factory_defaults = .true.
            test_passed = .false.
            if (i_case == 1) then
                settings_c%initialized = .false.
                error_c = oao_factory_c_wrapper( &
                    dm_ao_c, ao_overlap_c, 1_c_ip, n_ao_c, evaluate_dm_c_funptr, &
                    update_orbs_c_funptr, obj_func_c_funptr, extension_c, &
                    solver_settings_c, c_loc(state_c), settings_c)
                settings_c%initialized = .true.
            else
                error_c = oao_factory_c_wrapper( &
                    dm_ao_c, ao_overlap_c, 1_c_ip, n_ao_c, evaluate_dm_c_funptr, &
                    update_orbs_c_funptr, obj_func_c_funptr, extension_c, &
                    solver_settings_c, c_loc(state_c))
            end if
            mock_oao_factory_defaults = .false.
            if (error_c /= 0 .or. .not. test_passed) then
                write(stderr, *) "test_oao_factory_c_wrapper failed: Default "// &
                    "settings not handed to the factory for "// &
                    trim(default_cases(i_case))//"."
                test_oao_factory_c_wrapper = .false.
            else
                call c_f_pointer(extension_c, holder)
                select type (callbacks => holder%extension_callbacks)
                type is (orbital_c_callbacks_type)
                    if (associated(callbacks%logger) .or. &
                        c_associated(callbacks%host_context)) then
                        write(stderr, *) "test_oao_factory_c_wrapper failed: "// &
                            "Logging function or host context kept for "// &
                            trim(default_cases(i_case))//"."
                        test_oao_factory_c_wrapper = .false.
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
        update_orbs_c_funptr = c_funloc(mock_update_orbs)
        obj_func_c_funptr = c_funloc(mock_obj_func)
        extension_c = c_loc(state_c)
        test_passed = .false.
        error_c = oao_factory_c_wrapper(dm_ao_c, ao_overlap_c, 1_c_ip, n_ao_c, &
                                        evaluate_dm_c_funptr, update_orbs_c_funptr, &
                                        obj_func_c_funptr, extension_c, &
                                        solver_settings_c, c_null_ptr, settings_c)
        if (error_c /= error_oao + 1 .or. test_passed .or. &
            c_associated(update_orbs_c_funptr) .or. &
            c_associated(obj_func_c_funptr) .or. c_associated(extension_c)) then
            write(stderr, *) "test_oao_factory_c_wrapper failed: No error produced "// &
                "without a record of the state."
            test_oao_factory_c_wrapper = .false.
        end if

        ! restore the procedure pointers so later tests do not inherit the mocks
        oao_factory_cs => standard_oao_factory_cs
        oao_factory_os => standard_oao_factory_os

    end function test_oao_factory_c_wrapper

    logical(c_bool) function test_init_oao_settings_c() bind(C)
        !
        ! this function tests that the OAO settings initialization routine correctly
        ! initializes all settings to their default values
        !
        use otr_oao_c_interface, only: oao_settings_type_c, init_oao_settings_c
        use otr_oao, only: default_oao_settings
        use otr_oao_test_reference, only: operator(/=)

        type(oao_settings_type_c) :: settings

        ! assume test passes
        test_init_oao_settings_c = .true.

        ! initialize settings
        call init_oao_settings_c(settings)

        ! check function pointers
        if (c_associated(settings%logger)) then
            write(stderr, *) "test_init_oao_settings_c failed: Function pointers "// &
                "should not be initialized."
            test_init_oao_settings_c = .false.
        end if

        ! check settings
        if (settings /= default_oao_settings) then
            write(stderr, *) "test_init_oao_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_oao_settings_c = .false.
        end if

    end function test_init_oao_settings_c

    logical(c_bool) function test_assign_oao_f_c() bind(C)
        !
        ! this function tests that the function that converts OAO settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao, only: oao_settings_type, default_oao_settings
        use otr_oao_test_reference, only: get_reference_oao_values, ref_oao_settings, &
                                          operator(/=)
        use c_interface, only: c_callbacks_type
        use c_interface_unit_tests, only: mock_logger
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(oao_settings_type_c) :: settings_c
        type(oao_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_oao_f_c = .true.

        ! initialize the C settings with the reference values
        call get_reference_oao_values(settings_c)
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the factory does
        settings = settings_c
        callbacks%logger => mock_logger
        settings%context => callbacks

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_oao_f_c = .false.
            write(stderr, *) "test_assign_oao_f_c failed: Logging function not "// &
                "associated with value."
        else
            call arm_host_context_c(callbacks%host_context)
            call settings%logger("test", settings%context)
            if (.not. host_context%logger_called) then
                test_assign_oao_f_c = .false.
                write(stderr, *) "test_assign_oao_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
            if (.not. host_context_reached("assign_oao_f_c")) &
                test_assign_oao_f_c = .false.
        end if

        ! check against reference values
        if (settings /= ref_oao_settings) then
            write(stderr, *) "test_assign_oao_f_c failed: Settings not converted "// &
                "correctly."
            test_assign_oao_f_c = .false.
        end if

        ! convert C settings with a logging function that were not initialized, the
        ! custom values and logging function they still carry have to be replaced by the
        ! default settings
        settings_c%initialized = .false.
        settings = settings_c
        if (associated(settings%logger) .or. settings /= default_oao_settings) then
            write(stderr, *) "test_assign_oao_f_c failed: Settings that were not "// &
                "initialized not converted to default values."
            test_assign_oao_f_c = .false.
        end if

    end function test_assign_oao_f_c

    logical(c_bool) function test_assign_oao_c_f() bind(C)
        !
        ! this function tests that the function that converts OAO settings from Fortran
        ! to C correctly performs this conversion
        !
        use otr_oao, only: oao_settings_type, default_oao_settings
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao_test_reference, only: ref_oao_settings, get_reference_oao_values, &
                                          operator(/=)

        type(oao_settings_type) :: settings
        type(oao_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_oao_c_f = .true.

        ! convert Fortran settings with the reference values to C settings which hold a
        ! logging function and a host context that the conversion has to discard
        settings = ref_oao_settings
        call get_reference_oao_values(settings_c)
        settings_c = settings

        ! check that no callback function pointers and host contexts are associated
        if (c_associated(settings_c%logger)) then
            test_assign_oao_c_f = .false.
            write(stderr, *) "test_assign_oao_c_f failed: Logger function associated."
        end if
        if (c_associated(settings_c%context)) then
            write(stderr, *) "test_assign_oao_c_f failed: Host context associated."
            test_assign_oao_c_f = .false.
        end if

        ! check against reference values
        if (settings_c /= ref_oao_settings) then
            write(stderr, *) "test_assign_oao_c_f failed: Settings not converted "// &
                "correctly."
            test_assign_oao_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_oao_c_f = .false.
            write(stderr, *) "test_assign_oao_c_f failed: Settings not marked as "// &
                "initialized."
        end if

        ! convert the default settings into C settings which hold the reference values
        ! and are marked as not initialized, so that a field which is not converted,
        ! including the initialized flag, keeps a value that differs from the expected
        ! one
        call get_reference_oao_values(settings_c)
        settings_c%initialized = .false.
        settings = default_oao_settings
        settings_c = settings
        if (settings_c /= default_oao_settings) then
            write(stderr, *) "test_assign_oao_c_f failed: Default settings not "// &
                "converted correctly."
            test_assign_oao_c_f = .false.
        end if

    end function test_assign_oao_c_f

end module otr_oao_c_interface_unit_tests
