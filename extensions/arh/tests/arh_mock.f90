! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_mock

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_ip, c_rp
    use otr_arh, only: arh_factory_mo_cs, arh_factory_mo_os, arh_factory_oao_cs, &
                       arh_factory_oao_os, arh_deconstructor
    use test_reference, only: tol
    use otr_arh_test_reference, only: ref_arh_settings

    implicit none

    logical :: test_passed

    ! MO coefficients passed to the mock MO factories, which the mock orbital updating
    ! function for the MO basis overwrites to test that they are rotated in place
    real(rp), pointer, contiguous :: mo_coeff_3d(:, :, :) => null()

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory_mo_cs), pointer :: mock_arh_factory_mo_cs_ptr => &
        mock_arh_factory_mo_cs
    procedure(arh_factory_mo_os), pointer :: mock_arh_factory_mo_os_ptr => &
        mock_arh_factory_mo_os
    procedure(arh_factory_oao_cs), pointer :: mock_arh_factory_oao_cs_ptr => &
        mock_arh_factory_oao_cs
    procedure(arh_factory_oao_os), pointer :: mock_arh_factory_oao_os_ptr => &
        mock_arh_factory_oao_os
    procedure(arh_deconstructor), pointer :: mock_arh_deconstructor_ptr => &
        mock_arh_deconstructor

contains

    subroutine mock_arh_set_solver_settings(solver_settings, project)
        !
        ! this subroutine wires the OAO mock routines and, if given, a projection into
        ! the solver settings and asks for the Hessian refresh and the raised micro
        ! iteration limit, as the ARH factories do
        !
        use opentrustregion, only: solver_settings_type, project_type
        use otr_arh, only: arh_n_micro
        use otr_common_mock, only: mock_precond, mock_precond_pd, &
                                   mock_get_extra_trial_vectors

        type(solver_settings_type), intent(inout) :: solver_settings
        procedure(project_type), optional :: project

        integer(ip) :: error

        if (.not. solver_settings%initialized) call solver_settings%init(error)
        solver_settings%precond => mock_precond
        solver_settings%precond_pd => mock_precond_pd
        solver_settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
        solver_settings%stability_settings%precond => mock_precond
        solver_settings%stability_settings%get_extra_trial_vectors => &
            mock_get_extra_trial_vectors
        if (present(project)) then
            solver_settings%project => project
            solver_settings%stability_settings%project => project
        end if
        solver_settings%refresh_hess = .true.
        solver_settings%hess_symm = .false.
        solver_settings%n_micro = arh_n_micro

    end subroutine mock_arh_set_solver_settings

    subroutine mock_update_orbs_mo(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this subroutine is a test subroutine for the orbital update function for
        ! orbitals parameterized in the MO basis
        !
        use opentrustregion, only: hess_x_type
        use otr_common_mock, only: orig_mock_update_orbs => mock_update_orbs

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        call orig_mock_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)

        mo_coeff_3d = 2.0_rp

    end subroutine mock_update_orbs_mo

    subroutine check_factory_mo_input(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, &
                                      n_mo, settings)
        !
        ! this subroutine checks the input the mock MO factories are passed by the C
        ! wrapper, which is the same for both spin cases apart from the occupations
        !
        use otr_arh, only: arh_settings_type
        use otr_arh_test_reference, only: operator(/=)
        use otr_mo_test_reference, only: n_mo_ref => n_mo, mo_coeff_pattern
        use otr_common_test_reference, only: n_ao_ref => n_ao, n_occ_ref => n_occ

        real(rp), intent(in) :: mo_coeff(:, :, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        type(arh_settings_type), intent(inout) :: settings

        real(rp), allocatable :: pattern(:, :, :)

        ! expected MO coefficients, whose values encode their indices
        pattern = real(mo_coeff_pattern(int(n_ao_ref, kind=c_ip), int( &
            n_mo_ref, kind=c_ip), size(n_occ, kind=c_ip), 0.0_c_rp), kind=rp)

        ! check passed arrays
        if (any(shape(mo_coeff) /= shape(pattern))) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed MO "// &
                "coefficients have wrong shape."
        else if (any(abs(mo_coeff - pattern) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed MO "// &
                "coefficients wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed AO "// &
                "overlap matrix wrong."
        end if

        ! check dimensions
        if (n_particle /= size(n_occ)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed number "// &
                "of particles wrong."
        end if
        if (n_ao /= n_ao_ref .or. n_mo /= n_mo_ref) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed number "// &
                "of AOs or MOs wrong."
        end if
        if (any(n_occ /= n_occ_ref(:size(n_occ)))) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed number "// &
                "of occupied orbitals wrong."
        end if

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

    end subroutine check_factory_mo_input

    subroutine mock_arh_factory_mo_cs( &
        mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm_funptr, &
        obj_func_arh_funptr, update_orbs_arh_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the closed-shell case with the orbitals
        ! parameterized in the MO basis
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_settings_type, evaluate_dm_cs_type
        use otr_arh_test_reference, only: test_evaluate_dm_cs_funptr
        use otr_common_mock, only: mock_obj_func

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ, n_particle, n_ao, n_mo
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed input
        call check_factory_mo_input(reshape( &
            mo_coeff, [size(mo_coeff, 1, kind=ip), size(mo_coeff, 2, kind=ip), 1_ip]), &
            ao_overlap, [n_occ], n_particle, n_ao, n_mo, settings)

        ! test passed density matrix evaluating function
        test_passed = test_passed .and. test_evaluate_dm_cs_funptr( &
            evaluate_dm_funptr, "arh_factory_mo_c_wrapper", " by given density "// &
            "matrix evaluating function with non-linear potential contribution")

        ! set output quantities
        error = 0
        obj_func_arh_funptr => mock_obj_func
        update_orbs_arh_funptr => mock_update_orbs_mo
        call mock_arh_set_solver_settings(solver_settings)
        mo_coeff_3d(1:n_ao, 1:n_mo, 1:1) => mo_coeff

    end subroutine mock_arh_factory_mo_cs

    subroutine mock_arh_factory_mo_os( &
        mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm_funptr, &
        obj_func_arh_funptr, update_orbs_arh_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the open-shell case with the orbitals
        ! parameterized in the MO basis
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_settings_type, evaluate_dm_os_type
        use otr_arh_test_reference, only: test_evaluate_dm_os_funptr
        use otr_common_mock, only: mock_obj_func

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed input
        call check_factory_mo_input(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, &
                                    n_mo, settings)

        ! test passed density matrix evaluating function
        test_passed = test_passed .and. test_evaluate_dm_os_funptr( &
            evaluate_dm_funptr, "arh_factory_mo_c_wrapper", " by given density "// &
            "matrix evaluating function with same- and opposite-spin potential "// &
            "contributions")

        ! set output quantities
        error = 0
        obj_func_arh_funptr => mock_obj_func
        update_orbs_arh_funptr => mock_update_orbs_mo
        call mock_arh_set_solver_settings(solver_settings)
        mo_coeff_3d => mo_coeff

    end subroutine mock_arh_factory_mo_os

    subroutine mock_arh_factory_oao_cs( &
        dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, obj_func_arh_funptr, &
        update_orbs_arh_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the closed-shell case
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_settings_type, evaluate_dm_cs_type
        use otr_arh_test_reference, only: test_evaluate_dm_cs_funptr, operator(/=)
        use otr_common_mock, only: mock_obj_func
        use otr_oao_mock, only: mock_update_orbs, mock_project_oao, dm_ao_3d
        use otr_common_test_reference, only: n_ao_ref => n_ao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed AO "// &
                "density matrix for closed-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed AO "// &
                "overlap matrix wrong."
        end if

        ! check number of particles
        if (n_particle /= 1) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed number "// &
                "of particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= n_ao_ref) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed number "// &
                "of AOs wrong."
        end if

        ! test passed density matrix evaluating function
        test_passed = test_passed .and. test_evaluate_dm_cs_funptr( &
            evaluate_dm_funptr, "arh_factory_oao_c_wrapper", " by given density "// &
            "matrix evaluating function with non-linear potential contribution")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed "// &
                "logging function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_arh_funptr => mock_obj_func
        update_orbs_arh_funptr => mock_update_orbs
        call mock_arh_set_solver_settings(solver_settings, mock_project_oao)
        dm_ao_3d(1:n_ao, 1:n_ao, 1:1) => dm_ao

    end subroutine mock_arh_factory_oao_cs

    subroutine mock_arh_factory_oao_os( &
        dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_os_funptr, &
        obj_func_arh_funptr, update_orbs_arh_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the open-shell case
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: evaluate_dm_os_type, arh_settings_type
        use otr_arh_test_reference, only: test_evaluate_dm_os_funptr, operator(/=)
        use otr_common_mock, only: mock_obj_func
        use otr_oao_mock, only: mock_update_orbs, mock_project_oao, dm_ao_3d
        use otr_common_test_reference, only: n_ao_ref => n_ao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed AO "// &
                "density matrix for open-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed AO "// &
                "overlap matrix wrong."
        end if

        ! check number of particles
        if (n_particle /= 2) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed number "// &
                "of particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= n_ao_ref) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed number "// &
                "of AOs wrong."
        end if

        ! test passed density matrix evaluating function
        test_passed = test_passed .and. test_evaluate_dm_os_funptr( &
            evaluate_dm_os_funptr, "arh_factory_oao_c_wrapper", " by given density "// &
            "matrix evaluating function with same- and opposite-spin potential "// &
            "contributions")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed "// &
                "logging function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_arh_funptr => mock_obj_func
        update_orbs_arh_funptr => mock_update_orbs
        call mock_arh_set_solver_settings(solver_settings, mock_project_oao)
        dm_ao_3d => dm_ao

    end subroutine mock_arh_factory_oao_os

    subroutine mock_arh_deconstructor()
        !
        ! this subroutine is a test function for the ARH deconstructor
        !
        test_passed = .true.

    end subroutine mock_arh_deconstructor

end module otr_arh_mock
