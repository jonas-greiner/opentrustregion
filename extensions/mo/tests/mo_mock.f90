! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_mock

    use opentrustregion, only: rp, ip, stderr, obj_func_type
    use c_interface, only: c_ip, c_rp
    use otr_mo, only: mo_factory_cs, mo_factory_os, mo_deconstructor
    use test_reference, only: tol
    use otr_mo_test_reference, only: ref_mo_settings, operator(/=)

    implicit none

    logical :: test_passed
    real(rp), pointer, contiguous :: mo_coeff_3d(:, :, :)

    ! create function pointers to ensure that routines comply with interface
    procedure(mo_factory_cs), pointer :: mock_mo_factory_cs_ptr => mock_mo_factory_cs
    procedure(mo_factory_os), pointer :: mock_mo_factory_os_ptr => mock_mo_factory_os
    procedure(mo_deconstructor), pointer :: mock_mo_deconstructor_ptr => &
        mock_mo_deconstructor

contains

    subroutine mock_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this subroutine is a test subroutine for the orbital update function
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

    end subroutine mock_update_orbs

    subroutine mock_mo_set_solver_settings(solver_settings)
        !
        ! this subroutine wires the MO mock routines into the solver settings
        !
        use opentrustregion, only: solver_settings_type
        use otr_common_mock, only: mock_precond, mock_precond_pd, &
                                   mock_get_extra_trial_vectors

        type(solver_settings_type), intent(inout) :: solver_settings

        integer(ip) :: error

        if (.not. solver_settings%initialized) call solver_settings%init(error)
        solver_settings%precond => mock_precond
        solver_settings%precond_pd => mock_precond_pd
        solver_settings%get_extra_trial_vectors => mock_get_extra_trial_vectors
        solver_settings%stability_settings%precond => mock_precond
        solver_settings%stability_settings%get_extra_trial_vectors => &
            mock_get_extra_trial_vectors

    end subroutine mock_mo_set_solver_settings

    subroutine mock_mo_factory_cs( &
        mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm_funptr, &
        obj_func_mo_funptr, update_orbs_mo_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the closed-shell case
        !
        use opentrustregion, only: update_orbs_type, solver_settings_type
        use otr_common_mock, only: mock_obj_func
        use otr_mo, only: mo_settings_type
        use otr_common, only: evaluate_dm_cs_type
        use otr_common_test_reference, only: test_evaluate_dm_cs_funptr, &
                                             n_ao_ref => n_ao, n_occ_ref => n_occ
        use otr_mo_test_reference, only: n_mo_ref => n_mo, mo_coeff_pattern

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ, n_particle, n_ao, n_mo
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_mo_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_mo_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(mo_settings_type), intent(inout) :: settings

        real(rp), allocatable :: pattern(:, :, :)

        ! initialize logical
        test_passed = .true.

        ! expected MO coefficients, whose values encode their indices
        pattern = real(mo_coeff_pattern(int(n_ao_ref, kind=c_ip), int( &
            n_mo_ref, kind=c_ip), 1_c_ip, 0.0_c_rp), kind=rp)

        ! check passed arrays
        if (any(shape(mo_coeff) /= shape(pattern(:, :, 1)))) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed MO "// &
                "coefficients for closed-shell case have wrong shape."
        else if (any(abs(mo_coeff - pattern(:, :, 1)) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed MO "// &
                "coefficients for closed-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed AO overlap "// &
                "matrix wrong."
        end if

        ! check number of occupied orbitals
        if (n_occ /= n_occ_ref(1)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed number of "// &
                "occupied orbitals wrong."
        end if

        ! check number of particles
        if (n_particle /= 1) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= n_ao_ref) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! check number of MOs
        if (n_mo /= n_mo_ref) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed number of "// &
                "MOs wrong."
        end if

        ! test passed density matrix evaluating function
        test_passed = test_passed .and. test_evaluate_dm_cs_funptr( &
            evaluate_dm_funptr, "mo_factory_c_wrapper", &
            " by given density matrix evaluating function")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_mo_settings) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_mo_funptr => mock_obj_func
        update_orbs_mo_funptr => mock_update_orbs
        call mock_mo_set_solver_settings(solver_settings)
        mo_coeff_3d(1:n_ao, 1:n_mo, 1:1) => mo_coeff

    end subroutine mock_mo_factory_cs

    subroutine mock_mo_factory_os( &
        mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, evaluate_dm_funptr, &
        obj_func_mo_funptr, update_orbs_mo_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the open-shell case
        !
        use opentrustregion, only: update_orbs_type, solver_settings_type
        use otr_common_mock, only: mock_obj_func
        use otr_mo, only: mo_settings_type
        use otr_common, only: evaluate_dm_os_type
        use otr_common_test_reference, only: test_evaluate_dm_os_funptr, &
                                             n_ao_ref => n_ao, n_occ_ref => n_occ
        use otr_mo_test_reference, only: n_mo_ref => n_mo, mo_coeff_pattern

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_mo_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_mo_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(mo_settings_type), intent(inout) :: settings

        real(rp), allocatable :: pattern(:, :, :)

        ! initialize logical
        test_passed = .true.

        ! expected MO coefficients, whose values encode their indices
        pattern = real(mo_coeff_pattern(int(n_ao_ref, kind=c_ip), int( &
            n_mo_ref, kind=c_ip), 2_c_ip, 0.0_c_rp), kind=rp)

        ! check passed arrays
        if (any(shape(mo_coeff) /= shape(pattern))) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed MO "// &
                "coefficients for open-shell case have wrong shape."
        else if (any(abs(mo_coeff - pattern) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed MO "// &
                "coefficients for open-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed AO overlap "// &
                "matrix wrong."
        end if

        ! check number of occupied orbitals
        if (any(n_occ /= n_occ_ref)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed number of "// &
                "occupied orbitals wrong."
        end if

        ! check number of particles
        if (n_particle /= 2) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= n_ao_ref) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! check number of MOs
        if (n_mo /= n_mo_ref) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed number of "// &
                "MOs wrong."
        end if

        ! test passed density matrix evaluating function
        test_passed = test_passed .and. test_evaluate_dm_os_funptr( &
            evaluate_dm_funptr, "mo_factory_c_wrapper", &
            " by given density matrix evaluating function")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_mo_settings) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_mo_funptr => mock_obj_func
        update_orbs_mo_funptr => mock_update_orbs
        call mock_mo_set_solver_settings(solver_settings)
        mo_coeff_3d => mo_coeff

    end subroutine mock_mo_factory_os

    subroutine mock_mo_deconstructor()
        !
        ! this subroutine is a test function for the MO deconstructor
        !
        test_passed = .true.

    end subroutine mock_mo_deconstructor

end module otr_mo_mock
