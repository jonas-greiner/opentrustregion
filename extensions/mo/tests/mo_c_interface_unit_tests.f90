! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_funloc, c_associated

    implicit none

contains

    logical(c_bool) function test_mo_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO factory for the closed- and the
        ! open-shell case, without and with irreps of the MOs
        !
        use otr_mo_c_interface, only: &
            mo_settings_type_c, mo_factory_cs, mo_factory_os, mo_factory_c_wrapper, &
            obj_func_mo_before_wrapping, update_orbs_mo_before_wrapping, &
            precond_mo_before_wrapping, precond_pd_mo_before_wrapping, &
            get_extra_trial_vectors_mo_before_wrapping
        use otr_mo_mock, only: mock_mo_factory_cs, mock_mo_factory_os, test_passed, &
                               mo_coeff_3d, mock_update_orbs, orbsym_passed
        use otr_mo_test_reference, only: assignment(=), ref_mo_settings, n_mo, &
                                         mo_coeff_pattern, case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_occ, n_particle, n_ao_c, &
                                             shell_names
        use c_interface_unit_tests, only: mock_logger, test_logger, mock_project
        use test_reference, only: test_obj_func_c_funptr, test_update_orbs_c_funptr, &
                                  test_precond_c_funptr, test_precond_pd_c_funptr, &
                                  test_get_extra_trial_vectors_c_funptr, ref_settings, &
                                  assignment(=), operator(/=), n_param_ref => n_param
        use otr_common_mock, only: mock_obj_func, mock_precond, mock_precond_pd, &
                                   mock_get_extra_trial_vectors
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_mo_unit_tests, only: ref_count_mo_params
        use c_interface, only: solver_settings_type_c
        use otr_common_c_interface_unit_tests, only: mock_evaluate_dm_cs, &
                                                     mock_evaluate_dm_os

        real(c_rp), allocatable :: ao_overlap_c(:, :), mo_coeff_c(:, :, :)
        integer(c_ip), allocatable :: n_occ_c(:), orbsym_c(:)
        integer(ip) :: irreps(n_mo, n_particle), n_param_expected, i_call, i_sym
        logical :: with_irreps
        type(c_funptr) :: evaluate_dm_c_funptr, obj_func_c_funptr, update_orbs_c_funptr
        type(mo_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        integer(c_ip) :: n_particle_c, error_c
        character(len=:), allocatable :: case_name

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
        settings_c%logger = c_funloc(mock_logger)

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

            ! initialize logger logical, global number of parameters and solver
            ! settings, uninitialized for the closed-shell case and set to reference
            ! values with a projection of the caller, which have to be kept, for the
            ! open-shell case
            test_logger = .false.
            n_param_global = -1
            if (n_particle_c == 1) then
                solver_settings_c%initialized = .false._c_bool
            else
                solver_settings_c = ref_settings
                solver_settings_c%project = c_funloc(mock_project)
                solver_settings_c%stability_settings%project = c_funloc(mock_project)
            end if

            ! set the irreps of the MOs passed to the wrapper, which an unallocated
            ! array passes as absent, and the expected number of parameters, the
            ! occupied-virtual pairs of the same irrep, with all MOs in one irrep if no
            ! irreps are passed
            irreps = 0
            if (allocated(orbsym_c)) deallocate(orbsym_c)
            if (with_irreps) then
                i_sym = findloc(case_names, case_name, dim=1, kind=ip)
                irreps(:, :n_particle_c) = case_irreps(:, :n_particle_c, i_sym)
                orbsym_c = int( &
                    reshape(irreps(:, :n_particle_c), [n_mo * n_particle_c]), kind=c_ip)
            end if
            n_param_expected = ref_count_mo_params(n_occ(:n_particle_c), &
                                                   irreps(:, :n_particle_c))

            ! clear the callback slots, which only the factory call may set
            nullify(obj_func_mo_before_wrapping, update_orbs_mo_before_wrapping, &
                    precond_mo_before_wrapping, precond_pd_mo_before_wrapping, &
                    get_extra_trial_vectors_mo_before_wrapping)

            ! call MO factory C wrapper
            error_c = mo_factory_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, &
                int(n_mo, kind=c_ip), evaluate_dm_c_funptr, obj_func_c_funptr, &
                update_orbs_c_funptr, solver_settings_c, settings_c, orbsym_c)

            ! check if the mock factory received the correct input and called the
            ! logging function
            if (.not. test_passed) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Mock factory "// &
                    "received wrong input for the "//case_name//" case."
            end if
            if (.not. test_logger) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Called logging "// &
                    "subroutine wrong for the "//case_name//" case."
            end if
            if (orbsym_passed .neqv. with_irreps) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Irreps of the "// &
                    "MOs passed on wrongly for the "//case_name//" case."
            end if

            ! check if output variables are as expected
            if (error_c /= 0) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Returned error "// &
                    "code wrong for the "//case_name//" case."
            end if

            ! check if the callback slots point to the returned and wired functions,
            ! without which calling these below would crash
            if (.not. ( &
                associated(obj_func_mo_before_wrapping, mock_obj_func) .and. &
                associated(update_orbs_mo_before_wrapping, mock_update_orbs) .and. &
                associated(precond_mo_before_wrapping, mock_precond) .and. &
                associated(precond_pd_mo_before_wrapping, mock_precond_pd) .and. &
                associated(get_extra_trial_vectors_mo_before_wrapping, &
                           mock_get_extra_trial_vectors))) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Callback slots "// &
                    "not set for the "//case_name//" case."
                nullify(mo_coeff_3d)
                return
            end if

            ! determine if the number of parameters of the MO basis was set
            if (n_param_global /= n_param_expected) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Number of "// &
                    "parameters not set for the "//case_name//" case."
            end if

            ! determine if the solver settings are initialized and set as the factory
            ! sets them, if initialized solver settings and the projection of the
            ! caller were kept and if no projection was wired
            if (.not. solver_settings_c%initialized) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Solver "// &
                    "settings not initialized for the "//case_name//" case."
            end if
            if (n_particle_c == 2) then
                if (solver_settings_c /= ref_settings) then
                    test_mo_factory_c_wrapper = .false.
                    write(stderr, *) "test_mo_factory_c_wrapper failed: "// &
                        "Initialized solver settings not kept for the "//case_name// &
                        " case."
                end if
            end if
            if ((c_associated(solver_settings_c%project) .or. &
                 c_associated(solver_settings_c%stability_settings%project)) .neqv. &
                n_particle_c == 2) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Projection "// &
                    "wired or projection supplied by the caller not kept for the "// &
                    case_name//" case."
            end if
            if (n_particle_c == 2 .and. .not. ( &
                c_associated(solver_settings_c%project, c_funloc(mock_project)) .and. &
                c_associated(solver_settings_c%stability_settings%project, &
                             c_funloc(mock_project)))) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Projection "// &
                    "supplied by the caller replaced for the "//case_name//" case."
            end if

            ! test the returned and wired functions for the reference number of
            ! parameters of the wrapper tests
            n_param_global = n_param_ref
            test_mo_factory_c_wrapper = &
                test_mo_factory_c_wrapper .and. test_obj_func_c_funptr( &
                    obj_func_c_funptr, "mo_factory_c_wrapper", &
                    " by returned objective function for the "//case_name//" case")
            test_mo_factory_c_wrapper = &
                test_mo_factory_c_wrapper .and. test_update_orbs_c_funptr( &
                    update_orbs_c_funptr, "mo_factory_c_wrapper", " by returned "// &
                    "orbital updating function for the "//case_name//" case")

            ! check if the MO coefficients were updated by the returned orbital
            ! updating function
            if (any(abs(mo_coeff_c - 2.0_c_rp) > tol)) then
                test_mo_factory_c_wrapper = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: MO "// &
                    "coefficients not updated correctly by returned orbital "// &
                    "updating function for the "//case_name//" case."
            end if
            deallocate(mo_coeff_c)

            ! test functions wired into the solver and stability check settings
            test_mo_factory_c_wrapper = &
                test_mo_factory_c_wrapper .and. test_precond_c_funptr( &
                    solver_settings_c%precond, "mo_factory_c_wrapper", &
                    " by wired level-shifted preconditioner function for the "// &
                    case_name//" case")
            test_mo_factory_c_wrapper = &
                test_mo_factory_c_wrapper .and. test_precond_pd_c_funptr( &
                    solver_settings_c%precond_pd, "mo_factory_c_wrapper", &
                    " by wired positive-definite preconditioner function for the "// &
                    case_name//" case")
            test_mo_factory_c_wrapper = &
                test_mo_factory_c_wrapper .and. test_get_extra_trial_vectors_c_funptr( &
                    solver_settings_c%get_extra_trial_vectors, "mo_factory_c_wrapper", &
                    " by wired extra trial vector function for the "//case_name// &
                    " case")
            test_mo_factory_c_wrapper = &
                test_mo_factory_c_wrapper .and. &
                test_precond_c_funptr(solver_settings_c%stability_settings%precond, &
                                      "mo_factory_c_wrapper", " by wired stability "// &
                                      "check level-shifted preconditioner function "// &
                                      "for the "//case_name//" case")
            test_mo_factory_c_wrapper = &
                test_mo_factory_c_wrapper .and. test_get_extra_trial_vectors_c_funptr( &
                    solver_settings_c%stability_settings%get_extra_trial_vectors, &
                    "mo_factory_c_wrapper", " by wired stability check extra trial "// &
                    "vector function for the "//case_name//" case")
        end do
        deallocate(ao_overlap_c)
        nullify(mo_coeff_3d)

    end function test_mo_factory_c_wrapper

    logical(c_bool) function test_evaluate_dm_mo_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating
        ! function
        !
        use otr_common, only: evaluate_dm_cs_type, evaluate_dm_os_type
        use otr_mo_c_interface, only: evaluate_dm_mo_before_wrapping, &
                                      evaluate_dm_mo_cs_f_wrapper, &
                                      evaluate_dm_mo_os_f_wrapper
        use otr_common_test_reference, only: test_evaluate_dm_cs_funptr, &
                                             test_evaluate_dm_os_funptr
        use otr_common_unit_tests, only: mock_requests
        use otr_common_c_interface_unit_tests, only: mock_evaluate_dm_cs, &
                                                     mock_evaluate_dm_os

        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        integer(ip) :: request

        ! inject mock subroutine
        evaluate_dm_mo_before_wrapping => mock_evaluate_dm_cs

        ! get pointer to subroutine
        evaluate_dm_cs_funptr => evaluate_dm_mo_cs_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        test_evaluate_dm_mo_f_wrapper = test_evaluate_dm_cs_funptr( &
            evaluate_dm_cs_funptr, "evaluate_dm_mo_f_wrapper", "")
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request=0, 3)])) then
            write(stderr, *) "test_evaluate_dm_mo_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for closed-shell case "// &
                "wrong."
            test_evaluate_dm_mo_f_wrapper = .false.
        end if

        ! inject mock subroutine
        evaluate_dm_mo_before_wrapping => mock_evaluate_dm_os

        ! get pointer to subroutine
        evaluate_dm_os_funptr => evaluate_dm_mo_os_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        test_evaluate_dm_mo_f_wrapper = &
            test_evaluate_dm_mo_f_wrapper .and. test_evaluate_dm_os_funptr( &
                evaluate_dm_os_funptr, "evaluate_dm_mo_f_wrapper", "")
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request=0, 3)])) then
            write(stderr, *) "test_evaluate_dm_mo_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for open-shell case wrong."
            test_evaluate_dm_mo_f_wrapper = .false.
        end if

    end function test_evaluate_dm_mo_f_wrapper

    logical(c_bool) function test_get_response_mo_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the response function
        !
        use otr_common, only: get_response_cs_type, get_response_os_type
        use otr_mo_c_interface, only: get_response_mo_before_wrapping, &
                                      get_response_mo_cs_f_wrapper, &
                                      get_response_mo_os_f_wrapper
        use otr_common_test_reference, only: test_get_response_cs_funptr, &
                                             test_get_response_os_funptr
        use otr_common_c_interface_unit_tests, only: mock_get_response_cs, &
                                                     mock_get_response_os

        procedure(get_response_cs_type), pointer :: get_response_cs_funptr
        procedure(get_response_os_type), pointer :: get_response_os_funptr

        ! inject mock subroutine
        get_response_mo_before_wrapping => mock_get_response_cs

        ! get pointer to subroutine
        get_response_cs_funptr => get_response_mo_cs_f_wrapper

        ! test response function wrapper
        test_get_response_mo_f_wrapper = test_get_response_cs_funptr( &
            get_response_cs_funptr, "get_response_mo_f_wrapper", "")

        ! inject mock subroutine
        get_response_mo_before_wrapping => mock_get_response_os

        ! get pointer to subroutine
        get_response_os_funptr => get_response_mo_os_f_wrapper

        ! test response function wrapper
        test_get_response_mo_f_wrapper = &
            test_get_response_mo_f_wrapper .and. test_get_response_os_funptr( &
                get_response_os_funptr, "get_response_mo_f_wrapper", "")

    end function test_get_response_mo_f_wrapper

    logical(c_bool) function test_obj_func_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO objective function
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_mo_c_interface, only: obj_func_mo_before_wrapping, obj_func_mo_c_wrapper
        use otr_common_mock, only: mock_obj_func
        use test_reference, only: test_obj_func_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        obj_func_mo_before_wrapping => mock_obj_func

        ! test objective function
        test_obj_func_mo_c_wrapper = test_obj_func_c_funptr( &
            c_funloc(obj_func_mo_c_wrapper), "obj_func_mo_c_wrapper", "")

    end function test_obj_func_mo_c_wrapper

    logical(c_bool) function test_update_orbs_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO orbital update
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_mo_c_interface, only: update_orbs_mo_before_wrapping, &
                                      update_orbs_mo_c_wrapper
        use otr_common_mock, only: mock_update_orbs
        use test_reference, only: test_update_orbs_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        update_orbs_mo_before_wrapping => mock_update_orbs

        ! test orbital updating function
        test_update_orbs_mo_c_wrapper = test_update_orbs_c_funptr( &
            c_funloc(update_orbs_mo_c_wrapper), "update_orbs_mo_c_wrapper", "")

    end function test_update_orbs_mo_c_wrapper

    logical(c_bool) function test_hess_x_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO Hessian linear transformation
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_mo_c_interface, only: hess_x_mo_before_wrapping, hess_x_mo_c_wrapper
        use otr_common_mock, only: mock_hess_x
        use test_reference, only: test_hess_x_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        hess_x_mo_before_wrapping => mock_hess_x

        ! test orbital updating function
        test_hess_x_mo_c_wrapper = test_hess_x_c_funptr(c_funloc(hess_x_mo_c_wrapper), &
                                                        "hess_x_mo_c_wrapper", "")

    end function test_hess_x_mo_c_wrapper

    logical(c_bool) function test_precond_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO level-shifted preconditioner
        ! subroutine
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_mo_c_interface, only: precond_mo_before_wrapping, precond_mo_c_wrapper
        use otr_common_mock, only: mock_precond
        use test_reference, only: test_precond_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        precond_mo_before_wrapping => mock_precond

        ! test preconditioner function
        test_precond_mo_c_wrapper = test_precond_c_funptr( &
            c_funloc(precond_mo_c_wrapper), "precond_mo_c_wrapper", "")

    end function test_precond_mo_c_wrapper

    logical(c_bool) function test_precond_pd_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO positive-definite preconditioner
        ! subroutine
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_mo_c_interface, only: precond_pd_mo_before_wrapping, &
                                      precond_pd_mo_c_wrapper
        use otr_common_mock, only: mock_precond_pd
        use test_reference, only: test_precond_pd_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        precond_pd_mo_before_wrapping => mock_precond_pd

        ! test preconditioner function
        test_precond_pd_mo_c_wrapper = test_precond_pd_c_funptr( &
            c_funloc(precond_pd_mo_c_wrapper), "precond_pd_mo_c_wrapper", "")

    end function test_precond_pd_mo_c_wrapper

    logical(c_bool) function test_get_extra_trial_vectors_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the MO extra trial vector subroutine
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_mo_c_interface, only: get_extra_trial_vectors_mo_before_wrapping, &
                                      get_extra_trial_vectors_mo_c_wrapper
        use otr_common_mock, only: mock_get_extra_trial_vectors
        use test_reference, only: test_get_extra_trial_vectors_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        get_extra_trial_vectors_mo_before_wrapping => mock_get_extra_trial_vectors

        ! test extra trial vector function
        test_get_extra_trial_vectors_mo_c_wrapper = &
            test_get_extra_trial_vectors_c_funptr( &
                c_funloc(get_extra_trial_vectors_mo_c_wrapper), &
                "get_extra_trial_vectors_mo_c_wrapper", "")

    end function test_get_extra_trial_vectors_mo_c_wrapper

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
        use c_interface_unit_tests, only: mock_logger, test_logger

        type(mo_settings_type_c) :: settings_c
        type(mo_settings_type) :: settings

        ! assume test passes
        test_assign_mo_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_mo_settings
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings
        settings = settings_c

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_mo_f_c = .false.
            write(stderr, *) "test_assign_mo_f_c failed: Logging function not "// &
                "associated with value."
        else
            test_logger = .true.
            call settings%logger("test")
            if (.not. test_logger) then
                test_assign_mo_f_c = .false.
                write(stderr, *) "test_assign_mo_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
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

        ! check that callback function pointers are not associated
        if (c_associated(settings_c%logger)) then
            test_assign_mo_c_f = .false.
            write(stderr, *) "test_assign_mo_c_f failed: Logger function associated."
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
