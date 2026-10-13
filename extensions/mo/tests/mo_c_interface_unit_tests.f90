! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_funloc, c_associated, &
                                           c_loc

    implicit none

contains

    logical(c_bool) function test_mo_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO factory for the closed- and the
        ! open-shell case, without and with irreps of the MOs
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, assignment(=)
        use otr_mo, only: mo_object
        use otr_mo_c_interface, only: update_orbs_mo_c_wrapper, mo_callbacks, &
                                      mo_settings_type_c, mo_factory_cs, &
                                      mo_factory_os, mo_factory_c_wrapper, mo_holder
        use otr_mo_mock, only: mock_mo_factory_cs, mock_mo_factory_os, &
                               mock_mo_factory_error, test_passed, orbsym_passed
        use otr_mo_test_reference, only: assignment(=), ref_mo_settings, n_mo, &
                                         mo_coeff_pattern, case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_occ, n_particle, n_ao_c, &
                                             shell_names
        use otr_common_c_interface_unit_tests, only: &
            mock_evaluate_dm_cs, mock_evaluate_dm_os, check_c_factory_wiring
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: host_context_type, arm_host_context_c, host_context, &
                                  host_context_reached

        real(c_rp), allocatable :: ao_overlap_c(:, :), mo_coeff_c(:, :, :)
        integer(c_ip), allocatable :: n_occ_c(:), orbsym_c(:)
        integer(ip) :: irreps(n_mo, n_particle), i_call, i_sym
        logical :: with_irreps
        type(c_funptr) :: evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr
        type(mo_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        integer(c_ip) :: n_particle_c, error_c
        character(len=:), allocatable :: case_name
        type(host_context_type), target :: other_host_context

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
            mock_mo_factory_error = 0
            error_c = mo_factory_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, &
                int(n_mo, kind=c_ip), evaluate_dm_c_funptr, update_orbs_c_funptr, &
                obj_func_c_funptr, solver_settings_c, settings_c, orbsym_c)

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

            ! check that the factory returned the C wrappers of the object and wired the
            ! solver settings with its holder
            if (.not. check_c_factory_wiring( &
                "mo_factory_c_wrapper", mo_holder, mo_object, update_orbs_c_funptr, &
                obj_func_c_funptr, solver_settings_c, &
                c_funloc(update_orbs_mo_c_wrapper))) test_mo_factory_c_wrapper = .false.

            ! call the factory with the solver settings it wired and another host
            ! context, which has to produce an error without calling the factory,
            ! returning functions or replacing the bundle holding the callback functions
            ! of the object
            settings_c%context = c_loc(other_host_context)
            error_c = mo_factory_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, &
                int(n_mo, kind=c_ip), evaluate_dm_c_funptr, update_orbs_c_funptr, &
                obj_func_c_funptr, solver_settings_c, settings_c, orbsym_c)
            if (error_c == 0 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. &
                .not. c_associated(mo_callbacks%host_context, c_loc(host_context))) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Solver "// &
                    "settings the holder was wired into not refused before the "// &
                    "object was changed for the "//case_name//" case."
            end if

            ! call the factory for the initial MO coefficients with a mock factory
            ! which produces an error, which has to be returned without functions or
            ! writing back the solver settings or changing the holder or the bundle of
            ! the object
            mo_coeff_c = mo_coeff_pattern(n_ao_c, int(n_mo, kind=c_ip), n_particle_c, &
                                          0.0_c_rp)
            solver_settings_c = default_solver_settings
            update_orbs_c_funptr = c_funloc(mock_logger_c)
            obj_func_c_funptr = c_funloc(mock_logger_c)
            mock_mo_factory_error = 1
            error_c = mo_factory_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, &
                int(n_mo, kind=c_ip), evaluate_dm_c_funptr, update_orbs_c_funptr, &
                obj_func_c_funptr, solver_settings_c, settings_c, orbsym_c)
            mock_mo_factory_error = 0
            if (error_c /= 1 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. &
                c_associated(solver_settings_c%context)) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Error of the "// &
                    "factory not passed on for the "//case_name//" case."
            end if
            if (.not. associated(mo_holder%solver_level%level, &
                                 mo_object%solver_level) .or. &
                .not. c_associated(mo_callbacks%host_context, c_loc(host_context))) then
                write(stderr, *) "test_mo_factory_c_wrapper failed: Holder or "// &
                    "callback functions of the object not restored after an error "// &
                    "of the factory for the "//case_name//" case."
                test_mo_factory_c_wrapper = .false.
            end if
        end do
        deallocate(ao_overlap_c)
        if (allocated(mo_object)) deallocate(mo_object)

    end function test_mo_factory_c_wrapper

    logical(c_bool) function test_update_orbs_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper of the orbital update of the MO object,
        ! which passes on what the shared C wrapper of the orbital update returns and
        ! copies the MO coefficients the object rotated back to the array of the C host
        ! only if the object updates a copy of it
        !
        use otr_mo, only: mo_object
        use otr_mo_c_interface, only: update_orbs_mo_c_wrapper, mo_coeff_3d_c
        use otr_common_c_interface, only: extension_holder_type, hess_x_c_wrapper
        use otr_common_unit_tests, only: mock_extension_type
        use otr_common_c_interface_unit_tests, only: setup_mock_extension_holder
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_mo_test_reference, only: n_mo
        use test_reference, only: tol_c, n_param

        type(mock_extension_type), target :: extension
        type(extension_holder_type), target :: holder
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(c_rp), target :: mo_coeff_c(n_ao, n_mo, n_particle)
        real(c_rp) :: kappa_c(n_param), func_c, grad_c(n_param), h_diag_c(n_param)
        type(c_funptr) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        ! assume tests pass
        test_update_orbs_mo_c_wrapper = .true.

        ! set up an extension object performing the orbital update, which returns an
        ! error, in a holder, and a step
        extension%n_param = n_param
        extension%solver_level%extension => extension
        extension%stability_level%extension => extension
        extension%provides%update_orbs = .true.
        extension%mock_error = 3
        call setup_mock_extension_holder(holder, extension, 1_ip)
        call random_number(kappa_c)

        ! set up the MO coefficients of the MO object and a different array of the C
        ! host, which the module points to
        allocate(mo_object)
        call random_number(mo_coeff)
        mo_object%mo_coeff => mo_coeff
        mo_coeff_c = 0.0_c_rp
        mo_coeff_3d_c => mo_coeff_c

        ! call the orbital update with the address of the level of the solver in the
        ! holder and determine if it passes on the results and the error of the
        ! object's orbital update and the C wrapper of the Hessian linear transformation
        ! the shared C wrapper returns
        error_c = update_orbs_mo_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                           hess_x_c_funptr, c_loc(holder%solver_level))
        if (error_c /= 3) then
            write(stderr, *) "test_update_orbs_mo_c_wrapper failed: Error not "// &
                "passed on."
            test_update_orbs_mo_c_wrapper = .false.
        end if
        if (abs(func_c - sum(kappa_c)) > tol_c .or. &
            any(abs(grad_c - 2.0_c_rp * kappa_c) > tol_c) .or. &
            any(abs(h_diag_c - 3.0_c_rp * kappa_c) > tol_c)) then
            write(stderr, *) "test_update_orbs_mo_c_wrapper failed: Results not "// &
                "passed on."
            test_update_orbs_mo_c_wrapper = .false.
        end if
        if (.not. c_associated(hess_x_c_funptr, c_funloc(hess_x_c_wrapper))) then
            write(stderr, *) "test_update_orbs_mo_c_wrapper failed: C wrapper of "// &
                "the Hessian linear transformation not passed on."
            test_update_orbs_mo_c_wrapper = .false.
        end if

        ! determine if the MO coefficients of the object is copied to the host's array
        ! only if the object updates a copy of it, since it otherwise updates the
        ! host's array itself
        if (rp == c_rp) then
            if (any(abs(mo_coeff_c) > tol_c)) then
                write(stderr, *) "test_update_orbs_mo_c_wrapper failed: Orbitals "// &
                    "copied although the object updates the host's array in place."
                test_update_orbs_mo_c_wrapper = .false.
            end if
        else if (any(abs(mo_coeff_c - real(mo_coeff, kind=c_rp)) > tol_c)) then
            write(stderr, *) "test_update_orbs_mo_c_wrapper failed: Orbitals not "// &
                "copied back to the host's array."
            test_update_orbs_mo_c_wrapper = .false.
        end if

        ! release the module's pointer and the object
        mo_coeff_3d_c => null()
        deallocate(mo_object)

    end function test_update_orbs_mo_c_wrapper

    logical(c_bool) function test_evaluate_dm_mo_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating
        ! function
        !
        use otr_common, only: evaluate_dm_cs_type, evaluate_dm_os_type
        use otr_mo_c_interface, only: mo_callbacks, evaluate_dm_mo_cs_f_wrapper, &
                                      evaluate_dm_mo_os_f_wrapper
        use otr_common_test_reference, only: test_evaluate_dm_cs_funptr, &
                                             test_evaluate_dm_os_funptr
        use otr_common_unit_tests, only: mock_requests
        use otr_common_c_interface_unit_tests, only: mock_evaluate_dm_cs, &
                                                     mock_evaluate_dm_os
        use test_reference, only: arm_host_context_c, host_context_reached

        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        class(*), pointer :: context
        integer(ip) :: request

        ! assume tests pass
        test_evaluate_dm_mo_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        mo_callbacks%evaluate_dm => mock_evaluate_dm_cs
        call arm_host_context_c(mo_callbacks%host_context)
        context => mo_callbacks

        ! get pointer to subroutine
        evaluate_dm_cs_funptr => evaluate_dm_mo_cs_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        if (.not. test_evaluate_dm_cs_funptr(evaluate_dm_cs_funptr, &
                                             "evaluate_dm_mo_f_wrapper", "", context)) &
            test_evaluate_dm_mo_f_wrapper = .false.
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request=0, 3)])) then
            write(stderr, *) "test_evaluate_dm_mo_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for closed-shell case "// &
                "wrong."
            test_evaluate_dm_mo_f_wrapper = .false.
        end if

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("evaluate_dm_mo_f_wrapper")) &
            test_evaluate_dm_mo_f_wrapper = .false.

        ! inject mock C function and arm the host context again
        mo_callbacks%evaluate_dm => mock_evaluate_dm_os
        call arm_host_context_c(mo_callbacks%host_context)

        ! get pointer to subroutine
        evaluate_dm_os_funptr => evaluate_dm_mo_os_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        if (.not. test_evaluate_dm_os_funptr(evaluate_dm_os_funptr, &
                                             "evaluate_dm_mo_f_wrapper", "", context)) &
            test_evaluate_dm_mo_f_wrapper = .false.
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request=0, 3)])) then
            write(stderr, *) "test_evaluate_dm_mo_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for open-shell case wrong."
            test_evaluate_dm_mo_f_wrapper = .false.
        end if

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("evaluate_dm_mo_f_wrapper")) &
            test_evaluate_dm_mo_f_wrapper = .false.

    end function test_evaluate_dm_mo_f_wrapper

    logical(c_bool) function test_get_response_mo_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the response function
        !
        use otr_common, only: get_response_cs_type, get_response_os_type
        use otr_mo_c_interface, only: mo_callbacks, get_response_mo_cs_f_wrapper, &
                                      get_response_mo_os_f_wrapper
        use otr_common_test_reference, only: test_get_response_cs_funptr, &
                                             test_get_response_os_funptr
        use otr_common_c_interface_unit_tests, only: mock_get_response_cs, &
                                                     mock_get_response_os
        use test_reference, only: arm_host_context_c, host_context_reached

        procedure(get_response_cs_type), pointer :: get_response_cs_funptr
        procedure(get_response_os_type), pointer :: get_response_os_funptr
        class(*), pointer :: context

        ! assume tests pass
        test_get_response_mo_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        mo_callbacks%get_response => mock_get_response_cs
        call arm_host_context_c(mo_callbacks%host_context)
        context => mo_callbacks

        ! get pointer to subroutine
        get_response_cs_funptr => get_response_mo_cs_f_wrapper

        ! test response function wrapper
        if (.not. test_get_response_cs_funptr( &
            get_response_cs_funptr, "get_response_mo_f_wrapper", "", context)) &
            test_get_response_mo_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("get_response_mo_f_wrapper")) &
            test_get_response_mo_f_wrapper = .false.

        ! inject mock C function and arm the host context again
        mo_callbacks%get_response => mock_get_response_os
        call arm_host_context_c(mo_callbacks%host_context)

        ! get pointer to subroutine
        get_response_os_funptr => get_response_mo_os_f_wrapper

        ! test response function wrapper
        if (.not. test_get_response_os_funptr( &
            get_response_os_funptr, "get_response_mo_f_wrapper", "", context)) &
            test_get_response_mo_f_wrapper = .false.

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("get_response_mo_f_wrapper")) &
            test_get_response_mo_f_wrapper = .false.

    end function test_get_response_mo_f_wrapper

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

    logical(c_bool) function test_mo_deconstructor_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO deconstructor
        !
        use otr_mo_c_interface, only: mo_deconstructor, mo_deconstructor_c_wrapper
        use otr_mo_mock, only: mock_mo_deconstructor, test_passed

        ! assume tests pass
        test_mo_deconstructor_c_wrapper = .true.

        ! inject mock functions
        mo_deconstructor => mock_mo_deconstructor

        ! initialize test logical
        test_passed = .false.

        ! call MO orbital updating deconstructor C wrapper
        call mo_deconstructor_c_wrapper()

        ! check if test has passed
        if (.not. test_passed) then
            test_mo_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_mo_deconstructor_c_wrapper failed: "// &
                "Deconstructor called wrong."
        end if

    end function test_mo_deconstructor_c_wrapper

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
