! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol, tol_c
    use otr_arh_c_interface, only: evaluate_dm_os_c_type, evaluate_dm_cs_c_type
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_funptr, c_funloc, &
                                           c_f_procpointer

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(evaluate_dm_cs_c_type), pointer :: mock_arh_evaluate_dm_cs_ptr => &
        mock_arh_evaluate_dm_cs
    procedure(evaluate_dm_os_c_type), pointer :: mock_arh_evaluate_dm_os_ptr => &
        mock_arh_evaluate_dm_os

contains

    function mock_arh_evaluate_dm_cs(dm_ao, energy, fock, v_nonlinear) result(error) &
        bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function with a separate non-linear potential contribution for the
        ! closed-shell case
        !
        use otr_common_test_reference, only: n_ao
        use otr_arh_test_reference, only: evaluate_dm_cs_factors
        use otr_common_unit_tests, only: record_mock_call

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*), v_nonlinear(*)
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(v_nonlinear)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) &
            fock(:flat_len) = evaluate_dm_cs_factors(1) * dm_ao(:flat_len)
        if (present(v_nonlinear)) &
            v_nonlinear(:flat_len) = evaluate_dm_cs_factors(2) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_arh_evaluate_dm_cs

    function mock_arh_evaluate_dm_os(dm_ao, energy, fock, v_same_spin, &
                                     v_opposite_spin, v_nonlinear) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function with separate same- and opposite-spin potential contributions and a
        ! non-linear potential contribution for the open-shell case
        !
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_arh_test_reference, only: evaluate_dm_os_factors
        use otr_common_unit_tests, only: record_mock_call

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*), v_same_spin(*), &
                                             v_opposite_spin(*), v_nonlinear(*)
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2 * n_particle

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(v_same_spin)) + &
                              merge(4_ip, 0_ip, present(v_opposite_spin)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) &
            fock(:flat_len) = evaluate_dm_os_factors(1) * dm_ao(:flat_len)
        if (present(v_same_spin)) &
            v_same_spin(:flat_len) = evaluate_dm_os_factors(2) * dm_ao(:flat_len)
        if (present(v_opposite_spin)) &
            v_opposite_spin(:flat_len) = evaluate_dm_os_factors(3) * dm_ao(:flat_len)
        if (present(v_nonlinear)) &
            v_nonlinear(:flat_len) = evaluate_dm_os_factors(4) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_arh_evaluate_dm_os

    logical(c_bool) function test_arh_factory_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH factory for orbitals
        ! parameterized in the MO basis for the closed- and the open-shell case
        !
        use otr_arh_c_interface, only: arh_settings_type_c, arh_factory_mo_cs, &
                                       arh_factory_mo_os, arh_factory_mo_c_wrapper, &
                                       update_orbs_arh_before_wrapping
        use otr_arh_mock, only: mock_arh_factory_mo_cs, mock_arh_factory_mo_os, &
                                test_passed
        use otr_mo_mock, only: mo_coeff_3d, mock_update_orbs
        use otr_arh, only: arh_n_micro
        use otr_arh_test_reference, only: assignment(=), ref_arh_settings
        use otr_mo_test_reference, only: mo_coeff_pattern, n_mo
        use otr_common_test_reference, only: n_ao, n_particle, n_occ, n_ao_c
        use otr_common_unit_tests, only: shell_names
        use c_interface_unit_tests, only: mock_logger, test_logger, mock_project
        use test_reference, only: test_obj_func_c_funptr, test_update_orbs_c_funptr, &
                                  test_precond_c_funptr, test_precond_pd_c_funptr, &
                                  test_get_extra_trial_vectors_c_funptr, ref_settings, &
                                  assignment(=), operator(/=), n_param_ref => n_param
        use otr_mo_c_interface, only: &
            obj_func_mo_before_wrapping, precond_mo_before_wrapping, &
            precond_pd_mo_before_wrapping, get_extra_trial_vectors_mo_before_wrapping
        use otr_common_mock, only: mock_obj_func, mock_precond, mock_precond_pd, &
                                   mock_get_extra_trial_vectors
        use otr_common_c_interface, only: n_param_global => n_param
        use c_interface, only: solver_settings_type_c

        real(c_rp), allocatable :: ao_overlap_c(:, :), mo_coeff_c(:, :, :)
        type(c_funptr) :: evaluate_dm_c_funptr, obj_func_c_funptr, update_orbs_c_funptr
        type(arh_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        integer(c_ip) :: n_particle_c, error_c, n_mo_c, n_occ_c(n_particle)
        character(:), allocatable :: case_name

        ! assume tests pass
        test_arh_factory_mo_c_wrapper = .true.

        ! inject mock functions
        arh_factory_mo_cs => mock_arh_factory_mo_cs
        arh_factory_mo_os => mock_arh_factory_mo_os

        ! dimensions differ between AOs and MOs and occupations between the spin
        ! channels, so that the shapes and the number of parameters are checked
        n_mo_c = int(n_mo, kind=c_ip)
        n_occ_c = int(n_occ, kind=c_ip)

        ! initialize AO overlap matrix
        allocate(ao_overlap_c(n_ao, n_ao))
        ao_overlap_c = 2.0_c_rp

        ! associate optional settings with values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger)

        ! both spin cases pass through the same C wrapper
        do n_particle_c = 1, n_particle
            case_name = trim(shell_names(n_particle_c))

            ! initialize MO coefficients with values encoding their indices
            mo_coeff_c = mo_coeff_pattern(n_ao_c, n_mo_c, n_particle_c, 0.0_c_rp)

            ! get C function pointers to Fortran functions
            if (n_particle_c == 1) then
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_cs)
            else
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_os)
            end if

            ! initialize logger logical, global number of parameters and solver
            ! settings, uninitialized for the closed-shell case and, for the open-shell
            ! case, set to reference values with a projection supplied by the caller,
            ! which have to be kept
            test_logger = .false.
            n_param_global = -1
            if (n_particle_c == 1) then
                solver_settings_c%initialized = .false._c_bool
            else
                solver_settings_c = ref_settings
                solver_settings_c%project = c_funloc(mock_project)
                solver_settings_c%stability_settings%project = c_funloc(mock_project)
            end if

            ! clear the callback slots, which only the factory call may set
            nullify(obj_func_mo_before_wrapping, update_orbs_arh_before_wrapping, &
                    precond_mo_before_wrapping, precond_pd_mo_before_wrapping, &
                    get_extra_trial_vectors_mo_before_wrapping)

            ! call ARH MO factory C wrapper
            error_c = arh_factory_mo_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
                evaluate_dm_c_funptr, obj_func_c_funptr, update_orbs_c_funptr, &
                solver_settings_c, settings_c)

            ! check if the mock factory received the correct input and called the
            ! logging function
            if (.not. test_passed) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Mock "// &
                    "factory received wrong input for the "//case_name//" case."
            end if
            if (.not. test_logger) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Called "// &
                    "logging subroutine wrong for the "//case_name//" case."
            end if

            ! check if output variables are as expected
            if (error_c /= 0) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Returned "// &
                    "error code wrong for the "//case_name//" case."
            end if

            ! check if the callback slots point to the returned and wired functions,
            ! without which calling these below would crash
            if (.not. ( &
                associated(obj_func_mo_before_wrapping, mock_obj_func) .and. &
                associated(update_orbs_arh_before_wrapping, mock_update_orbs) .and. &
                associated(precond_mo_before_wrapping, mock_precond) .and. &
                associated(precond_pd_mo_before_wrapping, mock_precond_pd) .and. &
                associated(get_extra_trial_vectors_mo_before_wrapping, &
                           mock_get_extra_trial_vectors))) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Callback "// &
                    "slots not set for the "//case_name//" case."
                nullify(mo_coeff_3d)
                return
            end if

            ! determine if the number of parameters of the MO basis, the number of
            ! occupied-virtual pairs of all particle channels, was set
            if (n_param_global /= &
                sum(n_occ_c(:n_particle_c) * (n_mo_c - n_occ_c(:n_particle_c)))) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Number of "// &
                    "parameters not set for the "//case_name//" case."
            end if

            ! determine if the solver settings are initialized and set as the factory
            ! sets them, without a projection wired but with the projection supplied by
            ! the caller, and if initialized solver settings apart from the micro
            ! iteration limit were kept
            if (.not. solver_settings_c%initialized) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Solver "// &
                    "settings not initialized for the "//case_name//" case."
            end if
            if (.not. solver_settings_c%refresh_hess) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Hessian "// &
                    "refresh not requested for the "//case_name//" case."
            end if
            if (solver_settings_c%hess_symm) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Symmetry "// &
                    "of the approximate Hessian not passed on for the "//case_name// &
                    " case."
            end if
            if (solver_settings_c%n_micro /= arh_n_micro) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: Micro "// &
                    "iteration limit not passed on for the "//case_name//" case."
            end if
            if ((c_associated(solver_settings_c%project) .or. &
                 c_associated(solver_settings_c%stability_settings%project)) .neqv. &
                n_particle_c == 2) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: "// &
                    "Projection wired or projection supplied by the caller not "// &
                    "kept for the "//case_name//" case."
            end if
            if (n_particle_c == 2) then
                solver_settings_c%n_micro = int(ref_settings%n_micro, kind=c_ip)
                if (solver_settings_c /= ref_settings) then
                    test_arh_factory_mo_c_wrapper = .false.
                    write (stderr, *) "test_arh_factory_mo_c_wrapper failed: "// &
                        "Initialized solver settings not kept for the "//case_name// &
                        " case."
                end if
            end if

            ! test the returned and wired functions for the reference number of
            ! parameters of the wrapper tests
            n_param_global = n_param_ref
            test_arh_factory_mo_c_wrapper = &
                test_arh_factory_mo_c_wrapper .and. test_obj_func_c_funptr( &
                    obj_func_c_funptr, "arh_factory_mo_c_wrapper", &
                    " by returned objective function for the "//case_name//" case")
            test_arh_factory_mo_c_wrapper = &
                test_arh_factory_mo_c_wrapper .and. test_update_orbs_c_funptr( &
                    update_orbs_c_funptr, "arh_factory_mo_c_wrapper", " by "// &
                    "returned orbital updating function for the "//case_name//" case")

            ! check if the MO coefficients passed from C were rotated in place by the
            ! returned orbital updating function
            if (any(abs(mo_coeff_c - 2.0_c_rp) > tol_c)) then
                test_arh_factory_mo_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_mo_c_wrapper failed: MO "// &
                    "coefficients not updated correctly by returned orbital "// &
                    "updating function for the "//case_name//" case."
            end if
            deallocate(mo_coeff_c)

            ! test functions wired into the solver and stability check settings
            test_arh_factory_mo_c_wrapper = &
                test_arh_factory_mo_c_wrapper .and. test_precond_c_funptr( &
                    solver_settings_c%precond, "arh_factory_mo_c_wrapper", &
                    " by wired level-shifted preconditioner function for the "// &
                    case_name//" case")
            test_arh_factory_mo_c_wrapper = &
                test_arh_factory_mo_c_wrapper .and. test_precond_pd_c_funptr( &
                    solver_settings_c%precond_pd, "arh_factory_mo_c_wrapper", &
                    " by wired positive-definite preconditioner function for the "// &
                    case_name//" case")
            test_arh_factory_mo_c_wrapper = &
                test_arh_factory_mo_c_wrapper .and. &
                test_get_extra_trial_vectors_c_funptr( &
                    solver_settings_c%get_extra_trial_vectors, &
                    "arh_factory_mo_c_wrapper", " by wired extra trial vector "// &
                    "function for the "//case_name//" case")
            test_arh_factory_mo_c_wrapper = &
                test_arh_factory_mo_c_wrapper .and. test_precond_c_funptr( &
                    solver_settings_c%stability_settings%precond, &
                    "arh_factory_mo_c_wrapper", " by wired stability check "// &
                    "level-shifted preconditioner function for "// &
                    "the "//case_name//" case")
            test_arh_factory_mo_c_wrapper = &
                test_arh_factory_mo_c_wrapper .and. &
                test_get_extra_trial_vectors_c_funptr( &
                    solver_settings_c%stability_settings%get_extra_trial_vectors, &
                    "arh_factory_mo_c_wrapper", " by wired stability check extra "// &
                    "trial vector function for the "//case_name//" case")
        end do
        deallocate(ao_overlap_c)
        nullify(mo_coeff_3d)

    end function test_arh_factory_mo_c_wrapper

    logical(c_bool) function test_arh_factory_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH factory for orbitals
        ! parameterized in the OAO basis for the closed- and the open-shell case
        !
        use otr_arh_c_interface, only: arh_settings_type_c, arh_factory_oao_cs, &
                                       arh_factory_oao_os, arh_factory_oao_c_wrapper, &
                                       update_orbs_arh_before_wrapping
        use otr_arh_mock, only: mock_arh_factory_oao_cs, mock_arh_factory_oao_os, &
                                test_passed
        use otr_oao_mock, only: dm_ao_3d, mock_update_orbs, mock_project_oao
        use otr_arh, only: arh_n_micro
        use otr_arh_test_reference, only: assignment(=), ref_arh_settings
        use otr_common_test_reference, only: n_ao, n_particle, n_ao_c
        use otr_common_unit_tests, only: shell_names
        use c_interface_unit_tests, only: mock_logger, test_logger
        use test_reference, only: test_obj_func_c_funptr, test_update_orbs_c_funptr, &
                                  test_precond_c_funptr, test_precond_pd_c_funptr, &
                                  test_project_c_funptr, &
                                  test_get_extra_trial_vectors_c_funptr, ref_settings, &
                                  assignment(=), operator(/=), n_param_ref => n_param
        use otr_oao_c_interface, only: &
            obj_func_oao_before_wrapping, precond_oao_before_wrapping, &
            precond_pd_oao_before_wrapping, &
            get_extra_trial_vectors_oao_before_wrapping, project_oao_before_wrapping
        use otr_common_mock, only: mock_obj_func, mock_precond, mock_precond_pd, &
                                   mock_get_extra_trial_vectors
        use otr_common_c_interface, only: n_param_global => n_param
        use c_interface, only: solver_settings_type_c

        real(c_rp), allocatable :: ao_overlap_c(:, :), dm_ao_c(:, :, :)
        type(c_funptr) :: evaluate_dm_c_funptr, obj_func_c_funptr, update_orbs_c_funptr
        type(arh_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        integer(c_ip) :: n_particle_c, error_c
        character(:), allocatable :: case_name

        ! assume tests pass
        test_arh_factory_oao_c_wrapper = .true.

        ! inject mock functions
        arh_factory_oao_cs => mock_arh_factory_oao_cs
        arh_factory_oao_os => mock_arh_factory_oao_os

        ! initialize AO overlap matrix
        allocate(ao_overlap_c(n_ao, n_ao))
        ao_overlap_c = 2.0_c_rp

        ! associate optional settings with values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger)

        ! both spin cases pass through the same C wrapper
        do n_particle_c = 1, n_particle
            case_name = trim(shell_names(n_particle_c))

            ! initialize density matrix
            allocate(dm_ao_c(n_ao, n_ao, n_particle_c))
            dm_ao_c = 1.0_c_rp

            ! get C function pointers to Fortran functions
            if (n_particle_c == 1) then
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_cs)
            else
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_os)
            end if

            ! initialize logger logical, global number of parameters and solver
            ! settings, uninitialized for the closed-shell case and set to reference
            ! values, which have to be kept, for the open-shell case
            test_logger = .false.
            n_param_global = -1
            if (n_particle_c == 1) then
                solver_settings_c%initialized = .false._c_bool
            else
                solver_settings_c = ref_settings
            end if

            ! clear the callback slots, which only the factory call may set
            nullify(obj_func_oao_before_wrapping, update_orbs_arh_before_wrapping, &
                    precond_oao_before_wrapping, precond_pd_oao_before_wrapping, &
                    get_extra_trial_vectors_oao_before_wrapping, &
                    project_oao_before_wrapping)

            ! call ARH OAO factory C wrapper
            error_c = arh_factory_oao_c_wrapper( &
                dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
                obj_func_c_funptr, update_orbs_c_funptr, solver_settings_c, settings_c)

            ! check if the mock factory received the correct input and called the
            ! logging function
            if (.not. test_passed) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Mock "// &
                    "factory received wrong input for the "//case_name//" case."
            end if
            if (.not. test_logger) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Called "// &
                    "logging subroutine wrong for the "//case_name//" case."
            end if

            ! check if output variables are as expected
            if (error_c /= 0) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Returned "// &
                    "error code wrong for the "//case_name//" case."
            end if

            ! check if the callback slots point to the returned and wired functions,
            ! without which calling these below would crash
            if (.not. ( &
                associated(obj_func_oao_before_wrapping, mock_obj_func) .and. &
                associated(update_orbs_arh_before_wrapping, mock_update_orbs) .and. &
                associated(precond_oao_before_wrapping, mock_precond) .and. &
                associated(precond_pd_oao_before_wrapping, mock_precond_pd) .and. &
                associated(get_extra_trial_vectors_oao_before_wrapping, &
                           mock_get_extra_trial_vectors) .and. &
                associated(project_oao_before_wrapping, mock_project_oao))) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Callback "// &
                    "slots not set for the "//case_name//" case."
                nullify(dm_ao_3d)
                return
            end if

            ! determine if the number of parameters of the OAO basis was set
            if (n_param_global /= n_particle_c * n_ao * (n_ao - 1) / 2) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Number "// &
                    "of parameters not set for the "//case_name//" case."
            end if

            ! determine if the solver settings are initialized and set as the factory
            ! sets them, and if initialized solver settings apart from the micro
            ! iteration limit were kept
            if (.not. solver_settings_c%initialized) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Solver "// &
                    "settings not initialized for the "//case_name//" case."
            end if
            if (.not. solver_settings_c%refresh_hess) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Hessian "// &
                    "refresh not requested for the "//case_name//" case."
            end if
            if (solver_settings_c%hess_symm) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Symmetry "// &
                    "of the approximate Hessian not passed on for the "//case_name// &
                    " case."
            end if
            if (solver_settings_c%n_micro /= arh_n_micro) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Micro "// &
                    "iteration limit not passed on for the "//case_name//" case."
            end if
            if (n_particle_c == 2) then
                solver_settings_c%n_micro = int(ref_settings%n_micro, kind=c_ip)
                if (solver_settings_c /= ref_settings) then
                    test_arh_factory_oao_c_wrapper = .false.
                    write (stderr, *) "test_arh_factory_oao_c_wrapper failed: "// &
                        "Initialized solver settings not kept for the "//case_name// &
                        " case."
                end if
            end if

            ! test the returned and wired functions for the reference number of
            ! parameters of the wrapper tests
            n_param_global = n_param_ref
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. test_obj_func_c_funptr( &
                    obj_func_c_funptr, "arh_factory_oao_c_wrapper", &
                    " by returned objective function for the "//case_name//" case")
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. test_update_orbs_c_funptr( &
                    update_orbs_c_funptr, "arh_factory_oao_c_wrapper", " by "// &
                    "returned orbital updating function for the "//case_name//" case")

            ! check if the density matrix was updated by the returned orbital updating
            ! function
            if (any(abs(dm_ao_c - 2.0_c_rp) > tol)) then
                test_arh_factory_oao_c_wrapper = .false.
                write (stderr, *) "test_arh_factory_oao_c_wrapper failed: Density "// &
                    "matrix not updated correctly by returned orbital updating "// &
                    "function for the "//case_name//" case."
            end if
            deallocate(dm_ao_c)

            ! test functions wired into the solver and stability check settings
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. test_precond_c_funptr( &
                    solver_settings_c%precond, "arh_factory_oao_c_wrapper", &
                    " by wired level-shifted preconditioner function for the "// &
                    case_name//" case")
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. test_precond_pd_c_funptr( &
                    solver_settings_c%precond_pd, "arh_factory_oao_c_wrapper", &
                    " by wired positive-definite preconditioner function for the "// &
                    case_name//" case")
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. test_project_c_funptr( &
                    solver_settings_c%project, "arh_factory_oao_c_wrapper", &
                    " by wired projection function for the "//case_name//" case")
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. &
                test_get_extra_trial_vectors_c_funptr( &
                    solver_settings_c%get_extra_trial_vectors, &
                    "arh_factory_oao_c_wrapper", " by wired extra trial vector "// &
                    "function for the "//case_name//" case")
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. test_precond_c_funptr( &
                    solver_settings_c%stability_settings%precond, &
                    "arh_factory_oao_c_wrapper", " by wired stability check "// &
                    "level-shifted preconditioner function for "// &
                    "the "//case_name//" case")
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. test_project_c_funptr( &
                    solver_settings_c%stability_settings%project, &
                    "arh_factory_oao_c_wrapper", " by wired stability check "// &
                    "projection function for the "//case_name//" case")
            test_arh_factory_oao_c_wrapper = &
                test_arh_factory_oao_c_wrapper .and. &
                test_get_extra_trial_vectors_c_funptr( &
                    solver_settings_c%stability_settings%get_extra_trial_vectors, &
                    "arh_factory_oao_c_wrapper", " by wired stability check extra "// &
                    "trial vector function for the "//case_name//" case")
        end do
        deallocate(ao_overlap_c)
        nullify(dm_ao_3d)

    end function test_arh_factory_oao_c_wrapper

    logical(c_bool) function test_evaluate_dm_cs_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating C
        ! function with a separate non-linear potential contribution for the
        ! closed-shell case
        !
        use otr_arh, only: evaluate_dm_cs_type
        use otr_arh_c_interface, only: evaluate_dm_cs_before_wrapping, &
                                       evaluate_dm_cs_f_wrapper
        use otr_arh_test_reference, only: test_evaluate_dm_cs_funptr
        use otr_common_unit_tests, only: mock_requests

        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        integer(ip) :: request

        ! inject mock subroutine
        evaluate_dm_cs_before_wrapping => mock_arh_evaluate_dm_cs

        ! get pointer to subroutine
        evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) ::]
        test_evaluate_dm_cs_f_wrapper = test_evaluate_dm_cs_funptr( &
            evaluate_dm_cs_funptr, "evaluate_dm_cs_f_wrapper", "")
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request = 0, 3)])) then
            write (stderr, *) "test_evaluate_dm_cs_f_wrapper failed: Outputs "// &
                "passed on to density matrix evaluating C function for "// &
                "closed-shell case wrong."
            test_evaluate_dm_cs_f_wrapper = .false.
        end if

    end function test_evaluate_dm_cs_f_wrapper

    logical(c_bool) function test_evaluate_dm_os_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating C
        ! function with same- and opposite-spin potential contributions for the
        ! open-shell case
        !
        use otr_arh, only: evaluate_dm_os_type
        use otr_arh_c_interface, only: evaluate_dm_os_before_wrapping, &
                                       evaluate_dm_os_f_wrapper
        use otr_arh_test_reference, only: test_evaluate_dm_os_funptr
        use otr_common_unit_tests, only: mock_requests

        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        integer(ip) :: request

        ! inject mock subroutine
        evaluate_dm_os_before_wrapping => mock_arh_evaluate_dm_os

        ! get pointer to subroutine
        evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) ::]
        test_evaluate_dm_os_f_wrapper = test_evaluate_dm_os_funptr( &
            evaluate_dm_os_funptr, "evaluate_dm_os_f_wrapper", "")
        if (size(mock_requests) /= 16 .or. &
            .not. all([(any(mock_requests == request), request = 0, 15)])) then
            write (stderr, *) "test_evaluate_dm_os_f_wrapper failed: Outputs "// &
                "passed on to density matrix evaluating C function for open-shell "// &
                "case wrong."
            test_evaluate_dm_os_f_wrapper = .false.
        end if

    end function test_evaluate_dm_os_f_wrapper

    logical(c_bool) function test_update_orbs_arh_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH orbital update
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_arh_c_interface, only: update_orbs_arh_before_wrapping, &
                                       update_orbs_arh_c_wrapper
        use otr_common_mock, only: mock_update_orbs
        use test_reference, only: test_update_orbs_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        update_orbs_arh_before_wrapping => mock_update_orbs

        ! test orbital updating function
        test_update_orbs_arh_c_wrapper = test_update_orbs_c_funptr( &
            c_funloc(update_orbs_arh_c_wrapper), "update_orbs_arh_c_wrapper", "")

    end function test_update_orbs_arh_c_wrapper

    logical(c_bool) function test_hess_x_arh_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH Hessian linear transformation
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_arh_c_interface, only: hess_x_arh_before_wrapping, hess_x_arh_c_wrapper
        use otr_common_mock, only: mock_hess_x
        use test_reference, only: test_hess_x_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        hess_x_arh_before_wrapping => mock_hess_x

        ! test orbital updating function
        test_hess_x_arh_c_wrapper = test_hess_x_c_funptr( &
            c_funloc(hess_x_arh_c_wrapper), "hess_x_arh_c_wrapper", "")

    end function test_hess_x_arh_c_wrapper

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

    logical(c_bool) function test_arh_deconstructor_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH deconstructor
        !
        use otr_arh_c_interface, only: arh_deconstructor, arh_deconstructor_c_wrapper
        use otr_arh_mock, only: mock_arh_deconstructor, test_passed

        ! assume tests pass
        test_arh_deconstructor_c_wrapper = .true.

        ! inject mock functions
        arh_deconstructor => mock_arh_deconstructor

        ! initialize test logical
        test_passed = .false.

        ! call ARH orbital updating deconstructor C wrapper
        call arh_deconstructor_c_wrapper()

        ! check if test has passed
        test_arh_deconstructor_c_wrapper = test_passed

        ! check if test has passed
        if (.not. test_passed) then
            test_arh_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_arh_deconstructor_c_wrapper failed: "// &
                "Deconstructor called wrong."
        end if

    end function test_arh_deconstructor_c_wrapper

    logical(c_bool) function test_assign_arh_f_c() bind(C)
        !
        ! this function tests that the function that converts ARH settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh, only: arh_settings_type
        use otr_arh_test_reference, only: assignment(=), operator(/=), ref_arh_settings
        use c_interface_unit_tests, only: mock_logger, test_logger

        type(arh_settings_type_c) :: settings_c
        type(arh_settings_type) :: settings

        ! assume test passes
        test_assign_arh_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings
        settings = settings_c

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_arh_f_c = .false.
            write(stderr, *) "test_assign_arh_f_c failed: Logging function not "// &
                "associated with value."
        else
            test_logger = .true.
            call settings%logger("test")
            if (.not. test_logger) then
                test_assign_arh_f_c = .false.
                write(stderr, *) "test_assign_arh_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
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

    end function test_assign_arh_f_c

    logical(c_bool) function test_assign_arh_c_f() bind(C)
        !
        ! this function tests that the function that converts ARH settings from Fortran
        ! to C correctly performs this conversion
        !
        use otr_arh, only: arh_settings_type
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh_test_reference, only: assignment(=), operator(/=), ref_arh_settings

        type(arh_settings_type)   :: settings
        type(arh_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_arh_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_arh_settings

        ! convert to C settings
        settings_c = settings

        ! check logging function
        if (c_associated(settings_c%logger)) then
            test_assign_arh_c_f = .false.
            write(stderr, *) "test_assign_arh_c_f failed: Logger function associated."
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
