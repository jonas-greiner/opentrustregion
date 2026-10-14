! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_associated, c_loc, c_f_pointer, &
                                           c_null_ptr, c_null_funptr

    implicit none

contains

    logical(c_bool) function test_mo_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO factory for the closed- and the
        ! open-shell case, without and with irreps of the MOs
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, solver_state_type_c, &
                               assignment(=)
        use otr_common, only: extension_level_type
        use otr_common_unit_tests, only: wired_level
        use otr_common_c_interface, only: extension_holder_type, &
                                          orbital_c_callbacks_type
        use otr_mo_c_interface, only: mo_settings_type_c, mo_factory_cs, &
                                      mo_factory_os, mo_factory_c_wrapper
        use otr_mo_mock, only: mock_mo_factory_cs, mock_mo_factory_os, &
                               mock_mo_factory_error, mock_mo_factory_defaults, &
                               mock_mo_factory_no_evaluate_dm, test_passed, &
                               orbsym_passed
        use otr_mo_test_reference, only: assignment(=), ref_mo_settings, n_mo, &
                                         mo_coeff_pattern, case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_occ, n_particle, n_ao_c, &
                                             shell_names
        use otr_common_c_interface_unit_tests, only: &
            mock_evaluate_dm_cs, mock_evaluate_dm_os, check_c_factory_wiring
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        real(c_rp), allocatable :: ao_overlap_c(:, :), mo_coeff_c(:, :, :)
        integer(c_ip), allocatable :: n_occ_c(:), orbsym_c(:)
        integer(ip) :: irreps(n_mo, n_particle), i_call, i_sym, i_case
        logical :: with_irreps
        type(c_funptr) :: evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr
        type(mo_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        type(solver_state_type_c), target :: state_c
        type(c_ptr) :: extension_c
        type(extension_holder_type), pointer :: holder
        type(extension_level_type), pointer :: level
        integer(c_ip) :: n_particle_c, error_c
        character(len=:), allocatable :: case_name
        character(len=34), parameter :: default_cases(2) = &
            [character(len=34) :: "settings that were not initialized", "no settings"]

        ! assume tests pass
        test_mo_factory_c_wrapper = .true.

        ! inject mock functions
        mo_factory_cs => mock_mo_factory_cs
        mo_factory_os => mock_mo_factory_os

        ! initialize AO overlap matrix
        allocate(ao_overlap_c(n_ao, n_ao))
        ao_overlap_c = 2.0_c_rp

        ! associate optional settings with values
        settings_c = ref_mo_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! both spin cases pass through the same C wrapper, first without and then with
        ! the irreps of the MOs of the symmetric occupation case of the shell
        do i_call = 1, 2 * n_particle
            n_particle_c = int(1 + mod(i_call - 1, n_particle), kind=c_ip)
            with_irreps = i_call > n_particle
            case_name = trim(shell_names(n_particle_c))
            if (with_irreps) case_name = case_name//" with symmetry"

            ! initialize MO coefficients with values encoding their indices and
            ! occupations
            mo_coeff_c = mo_coeff_pattern(n_ao_c, int(n_mo, kind=c_ip), n_particle_c, &
                                          0.0_c_rp)
            n_occ_c = int(n_occ(:n_particle_c), kind=c_ip)

            ! get C function pointers to Fortran functions
            if (n_particle_c == 1) then
                evaluate_dm_c_funptr = c_funloc(mock_evaluate_dm_cs)
            else
                evaluate_dm_c_funptr = c_funloc(mock_evaluate_dm_os)
            end if

            ! arm the host context of the MO settings, which the logging and density
            ! matrix evaluating functions receive
            call arm_host_context_c(settings_c%context)
            solver_settings_c = default_solver_settings

            ! set the irreps of the MOs passed to the wrapper, which an unallocated
            ! array passes as absent
            irreps = 0
            if (allocated(orbsym_c)) deallocate(orbsym_c)
            if (with_irreps) then
                i_sym = findloc(case_names, case_name, dim=1, kind=ip)
                irreps(:, :n_particle_c) = case_irreps(:, :n_particle_c, i_sym)
                orbsym_c = int( &
                    reshape(irreps(:, :n_particle_c), [n_mo * n_particle_c]), kind=c_ip)
            end if

            ! call MO factory C wrapper
            test_passed = .false.
            error_c = mo_factory_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, &
                int(n_mo, kind=c_ip), evaluate_dm_c_funptr, update_orbs_c_funptr, &
                obj_func_c_funptr, extension_c, solver_settings_c, c_loc(state_c), &
                settings_c, orbsym_c)

            ! check if the mock factory received the correct input and called the
            ! logging function
            if (.not. test_passed) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Mock factory "// &
                    "received wrong input for the "//case_name//" case."
            end if
            if (.not. host_context%logger_called) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Called logging "// &
                    "subroutine wrong for the "//case_name//" case."
            end if
            if (orbsym_passed .neqv. with_irreps) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Irreps of the "// &
                    "MOs passed on wrongly for the "//case_name//" case."
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("mo_factory_c_wrapper")) &
                test_mo_factory_c_wrapper = .false.

            ! check if output variables are as expected
            if (error_c /= 0) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Returned error "// &
                    "code wrong for the "//case_name//" case."
            end if

            ! check that the factory received the MO coefficients of the caller in
            ! place, unless the kinds differ and the object updates a copy
            if (rp == c_rp .and. any(abs(mo_coeff_c - 2.0_c_rp) > tol)) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: MO "// &
                    "coefficients of the caller not passed in place for the "// &
                    case_name//" case."
            end if

            ! the returned extension is the holder, whose inspection would crash
            ! without one
            if (.not. c_associated(extension_c)) then
                write(stderr, *) "test_mo_factory_c_wrapper failed: No extension "// &
                    "returned for the "//case_name//" case."
                test_mo_factory_c_wrapper = .false.
                return
            end if
            call c_f_pointer(extension_c, holder)

            ! check that the settings of the holder, which own the object, were wired
            ! with it, that it reads the state of the holder and that the holder holds
            ! the bundle of the density matrix evaluating function and the host's
            ! record of the state
            level => wired_level(holder%solver_settings)
            if (.not. associated(level)) then
                write(stderr, *) "test_mo_factory_c_wrapper failed: Settings of "// &
                    "the holder not wired with an object for the "//case_name//" case."
                test_mo_factory_c_wrapper = .false.
                return
            end if
            if (.not. associated(level%extension%state, holder%state)) then
                write(stderr, *) "test_mo_factory_c_wrapper failed: Object does "// &
                    "not read the state of the holder for the "//case_name//" case."
                test_mo_factory_c_wrapper = .false.
            end if
            select type (callbacks => holder%extension_callbacks)
            type is (orbital_c_callbacks_type)
                if (n_particle_c == 1) then
                    if (.not. associated(callbacks%evaluate_dm, mock_evaluate_dm_cs)) &
                        then
                        write(stderr, *) "test_mo_factory_c_wrapper failed: "// &
                            "Density matrix evaluating function not held for the "// &
                            case_name//" case."
                        test_mo_factory_c_wrapper = .false.
                    end if
                else if (.not. associated(callbacks%evaluate_dm, mock_evaluate_dm_os)) &
                    then
                    write(stderr, *) "test_mo_factory_c_wrapper failed: Density "// &
                        "matrix evaluating function not held for the "//case_name// &
                        " case."
                    test_mo_factory_c_wrapper = .false.
                end if
            class default
                write(stderr, *) "test_mo_factory_c_wrapper failed: Bundle of the "// &
                    "callback functions not held for the "//case_name//" case."
                test_mo_factory_c_wrapper = .false.
            end select
            if (.not. (associated(holder%solver_level%state_c, state_c) .and. &
                       associated(holder%stability_level%state_c, state_c))) then
                write(stderr, *) "test_mo_factory_c_wrapper failed: Host's record "// &
                    "of the state not held for the "//case_name//" case."
                test_mo_factory_c_wrapper = .false.
            end if

            ! check that the factory returned the C wrappers of the object and wired the
            ! solver settings with its holder
            if (.not. check_c_factory_wiring("mo_factory_c_wrapper", holder, &
                                             level%extension, update_orbs_c_funptr, &
                                             obj_func_c_funptr, solver_settings_c)) &
                test_mo_factory_c_wrapper = .false.
            if (c_associated(extension_c)) then
                call c_f_pointer(extension_c, holder)
                deallocate(holder)
            end if

            ! call the factory without a density matrix evaluating function, which the
            ! factory has to be handed as none and the bundle must not hold
            mo_coeff_c = mo_coeff_pattern(n_ao_c, int(n_mo, kind=c_ip), n_particle_c, &
                                          0.0_c_rp)
            solver_settings_c = default_solver_settings
            mock_mo_factory_no_evaluate_dm = .true.
            test_passed = .false.
            error_c = mo_factory_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, &
                int(n_mo, kind=c_ip), c_null_funptr, update_orbs_c_funptr, &
                obj_func_c_funptr, extension_c, solver_settings_c, c_loc(state_c), &
                settings_c, orbsym_c)
            mock_mo_factory_no_evaluate_dm = .false.
            if (error_c /= 0 .or. .not. test_passed) then
                write(stderr, *) "test_mo_factory_c_wrapper failed: Missing "// &
                    "density matrix evaluating function not handed to the factory "// &
                    "for the "//case_name//" case."
                test_mo_factory_c_wrapper = .false.
            else
                call c_f_pointer(extension_c, holder)
                select type (callbacks => holder%extension_callbacks)
                type is (orbital_c_callbacks_type)
                    if (associated(callbacks%evaluate_dm)) then
                        write(stderr, *) "test_mo_factory_c_wrapper failed: "// &
                            "Density matrix evaluating function held without one "// &
                            "for the "//case_name//" case."
                        test_mo_factory_c_wrapper = .false.
                    end if
                end select
            end if
            if (c_associated(extension_c)) then
                call c_f_pointer(extension_c, holder)
                deallocate(holder)
            end if

            ! call the factory for the initial MO coefficients with a mock factory
            ! which produces an error, which has to be returned without functions, an
            ! extension or writing back the solver settings
            mo_coeff_c = mo_coeff_pattern(n_ao_c, int(n_mo, kind=c_ip), n_particle_c, &
                                          0.0_c_rp)
            solver_settings_c = default_solver_settings
            update_orbs_c_funptr = c_funloc(mock_logger_c)
            obj_func_c_funptr = c_funloc(mock_logger_c)
            mock_mo_factory_error = 1
            error_c = mo_factory_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, &
                int(n_mo, kind=c_ip), evaluate_dm_c_funptr, update_orbs_c_funptr, &
                obj_func_c_funptr, extension_c, solver_settings_c, c_loc(state_c), &
                settings_c, orbsym_c)
            mock_mo_factory_error = 0
            if (error_c /= 1 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. c_associated(extension_c) .or. &
                c_associated(solver_settings_c%context)) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Error of the "// &
                    "factory not passed on for the "//case_name//" case."
            end if
        end do

        ! call the factory with settings that were not initialized and without any, for
        ! which the factory is handed the default settings with the context through
        ! which the density matrix evaluating function is found, which the mock checks,
        ! and the bundle keeps no logging function or host context
        n_occ_c = int(n_occ(:1), kind=c_ip)
        evaluate_dm_c_funptr = c_funloc(mock_evaluate_dm_cs)
        do i_case = 1, size(default_cases)
            mo_coeff_c = mo_coeff_pattern(n_ao_c, int(n_mo, kind=c_ip), 1_c_ip, &
                                          0.0_c_rp)
            solver_settings_c = default_solver_settings
            mock_mo_factory_defaults = .true.
            test_passed = .false.
            if (i_case == 1) then
                settings_c%initialized = .false.
                error_c = mo_factory_c_wrapper( &
                    mo_coeff_c, ao_overlap_c, n_occ_c, 1_c_ip, n_ao_c, &
                    int(n_mo, kind=c_ip), evaluate_dm_c_funptr, update_orbs_c_funptr, &
                    obj_func_c_funptr, extension_c, solver_settings_c, c_loc(state_c), &
                    settings_c)
                settings_c%initialized = .true.
            else
                error_c = mo_factory_c_wrapper( &
                    mo_coeff_c, ao_overlap_c, n_occ_c, 1_c_ip, n_ao_c, &
                    int(n_mo, kind=c_ip), evaluate_dm_c_funptr, update_orbs_c_funptr, &
                    obj_func_c_funptr, extension_c, solver_settings_c, c_loc(state_c))
            end if
            mock_mo_factory_defaults = .false.
            if (error_c /= 0 .or. .not. test_passed) then
                write(stderr, *) "test_mo_factory_c_wrapper failed: Default "// &
                    "settings not handed to the factory for "// &
                    trim(default_cases(i_case))//"."
                test_mo_factory_c_wrapper = .false.
            else
                call c_f_pointer(extension_c, holder)
                select type (callbacks => holder%extension_callbacks)
                type is (orbital_c_callbacks_type)
                    if (associated(callbacks%logger) .or. &
                        c_associated(callbacks%host_context)) then
                        write(stderr, *) "test_mo_factory_c_wrapper failed: "// &
                            "Logging function or host context kept for "// &
                            trim(default_cases(i_case))//"."
                        test_mo_factory_c_wrapper = .false.
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
        test_passed = .false.
        error_c = mo_factory_c_wrapper( &
            mo_coeff_c, ao_overlap_c, n_occ_c, 1_c_ip, n_ao_c, int(n_mo, kind=c_ip), &
            evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr, &
            extension_c, solver_settings_c, c_null_ptr, settings_c)
        if (error_c == 0 .or. test_passed .or. c_associated(update_orbs_c_funptr) .or. &
            c_associated(obj_func_c_funptr) .or. c_associated(extension_c)) then
            write(stderr, *) "test_mo_factory_c_wrapper failed: No error produced "// &
                "without a record of the state."
            test_mo_factory_c_wrapper = .false.
        end if

    end function test_mo_factory_c_wrapper

    logical(c_bool) function test_init_mo_settings_c() bind(C)
        !
        ! this function tests that the MO settings initialization routine correctly
        ! initializes all settings to their default values
        !
        use otr_mo_c_interface, only: mo_settings_type_c, init_mo_settings_c
        use otr_mo, only: default_mo_settings
        use otr_mo_test_reference, only: operator(/=)

        type(mo_settings_type_c) :: settings

        ! assume test passes
        test_init_mo_settings_c = .true.

        ! initialize settings
        call init_mo_settings_c(settings)

        ! check function pointers
        if (c_associated(settings%logger)) then
            write(stderr, *) "test_init_mo_settings_c failed: Function pointers "// &
                "should not be initialized."
            test_init_mo_settings_c = .false.
        end if

        ! check settings
        if (settings /= default_mo_settings) then
            write(stderr, *) "test_init_mo_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_mo_settings_c = .false.
        end if

    end function test_init_mo_settings_c

    logical(c_bool) function test_assign_mo_f_c() bind(C)
        !
        ! this function tests that the function that converts MO settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_mo_c_interface, only: mo_settings_type_c, assignment(=)
        use otr_mo, only: mo_settings_type
        use otr_mo_test_reference, only: assignment(=), ref_mo_settings, operator(/=)
        use c_interface, only: c_callbacks_type
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(mo_settings_type_c) :: settings_c
        type(mo_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_mo_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_mo_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the factory does
        settings = settings_c
        callbacks%logger => mock_logger_c
        settings%context => callbacks

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_mo_f_c = .false.
            write(stderr, *) "test_assign_mo_f_c failed: Logging function not "// &
                "associated with value."
        else
            call arm_host_context_c(callbacks%host_context)
            call settings%logger("test", settings%context)
            if (.not. host_context%logger_called) then
                test_assign_mo_f_c = .false.
                write(stderr, *) "test_assign_mo_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
            if (.not. host_context_reached("assign_mo_f_c")) &
                test_assign_mo_f_c = .false.
        end if

        ! check against reference values
        if (settings /= ref_mo_settings) then
            write(stderr, *) "test_assign_mo_f_c failed: Settings not converted "// &
                "correctly."
            test_assign_mo_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            write(stderr, *) "test_assign_mo_f_c failed: Settings not marked as "// &
                "initialized."
            test_assign_mo_f_c = .false.
        end if

        ! convert C settings that were not initialized and check that they are left
        ! uninitialized without a logging function, so that the factory initializes
        ! them after its warning
        settings_c%initialized = .false.
        settings = settings_c
        if (settings%initialized .or. associated(settings%logger)) then
            write(stderr, *) "test_assign_mo_f_c failed: Settings that were not "// &
                "initialized converted."
            test_assign_mo_f_c = .false.
        end if

    end function test_assign_mo_f_c

    logical(c_bool) function test_assign_mo_c_f() bind(C)
        !
        ! this function tests that the function that converts MO settings from Fortran
        ! to C correctly performs this conversion
        !
        use otr_mo, only: mo_settings_type
        use otr_mo_c_interface, only: mo_settings_type_c, assignment(=)
        use otr_mo_test_reference, only: ref_mo_settings, assignment(=), operator(/=)

        type(mo_settings_type) :: settings
        type(mo_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_mo_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_mo_settings

        ! convert to C settings
        settings_c = settings

        ! check that no callback function pointers and host contexts are associated
        if (c_associated(settings_c%logger)) then
            test_assign_mo_c_f = .false.
            write(stderr, *) "test_assign_mo_c_f failed: Logger function associated."
        end if
        if (c_associated(settings_c%context)) then
            write(stderr, *) "test_assign_mo_c_f failed: Host context associated."
            test_assign_mo_c_f = .false.
        end if

        ! check against reference values
        if (settings_c /= ref_mo_settings) then
            write(stderr, *) "test_assign_mo_c_f failed: Settings not converted "// &
                "correctly."
            test_assign_mo_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_mo_c_f = .false.
            write(stderr, *) "test_assign_mo_c_f failed: Settings not marked as "// &
                "initialized."
        end if

    end function test_assign_mo_c_f

end module otr_mo_c_interface_unit_tests
