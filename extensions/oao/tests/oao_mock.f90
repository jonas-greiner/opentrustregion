! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_mock

    use opentrustregion, only: rp, ip, stderr, obj_func_type, precond_type, &
                               precond_pd_type, project_type, &
                               get_extra_trial_vectors_type
    use otr_oao, only: oao_factory_cs, oao_factory_os, oao_deconstructor
    use test_reference, only: tol
    use otr_oao_test_reference, only: ref_oao_settings, operator(/=)

    implicit none

    logical :: test_passed
    real(rp), pointer, contiguous :: dm_ao_3d(:, :, :)

    ! create function pointers to ensure that routines comply with interface
    procedure(oao_factory_cs), pointer :: mock_oao_factory_cs_ptr => mock_oao_factory_cs
    procedure(oao_factory_os), pointer :: mock_oao_factory_os_ptr => mock_oao_factory_os
    procedure(oao_deconstructor), pointer :: mock_oao_deconstructor_ptr => &
        mock_oao_deconstructor
    procedure(obj_func_type), pointer :: mock_obj_func_oao_ptr => mock_obj_func_oao
    procedure(precond_type), pointer :: mock_precond_oao_ptr => mock_precond_oao
    procedure(precond_pd_type), pointer :: mock_precond_pd_oao_ptr => &
        mock_precond_pd_oao
    procedure(project_type), pointer ::  mock_project_oao_ptr => mock_project_oao
    procedure(get_extra_trial_vectors_type), pointer :: &
        mock_get_extra_trial_vectors_oao_ptr => mock_get_extra_trial_vectors_oao

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

        dm_ao_3d = 2.0_rp

    end subroutine mock_update_orbs

    subroutine mock_oao_set_solver_settings(solver_settings)
        !
        ! this subroutine wires the OAO mock routines into the solver settings
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type), intent(inout) :: solver_settings

        integer(ip) :: error

        if (.not. solver_settings%initialized) call solver_settings%init(error)
        solver_settings%precond => mock_precond_oao
        solver_settings%precond_pd => mock_precond_pd_oao
        solver_settings%project => mock_project_oao
        solver_settings%get_extra_trial_vectors => mock_get_extra_trial_vectors_oao
        solver_settings%stability_settings%precond => mock_precond_oao
        solver_settings%stability_settings%project => mock_project_oao
        solver_settings%stability_settings%get_extra_trial_vectors => &
            mock_get_extra_trial_vectors_oao

    end subroutine mock_oao_set_solver_settings

    subroutine mock_oao_factory_cs( &
        dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, obj_func_oao_funptr, &
        update_orbs_oao_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the closed-shell case
        !
        use opentrustregion, only: update_orbs_type, solver_settings_type
        use otr_oao, only: evaluate_dm_cs_type, oao_settings_type
        use otr_oao_test_reference, only: test_evaluate_dm_cs_funptr

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed AO "// &
                "density matrix for closed-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed AO "// &
                "overlap matrix wrong."
        end if

        ! check number of particles
        if (n_particle /= 1) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= 3) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! test passed density matrix evaluating function
        test_passed = test_passed .and. test_evaluate_dm_cs_funptr( &
            evaluate_dm_funptr, "oao_factory_c_wrapper", &
            " by given density matrix evaluating function")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_oao_settings) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_oao_funptr => mock_obj_func_oao
        update_orbs_oao_funptr => mock_update_orbs
        call mock_oao_set_solver_settings(solver_settings)
        dm_ao_3d(1:n_ao, 1:n_ao, 1:1) => dm_ao

    end subroutine mock_oao_factory_cs

    subroutine mock_oao_factory_os( &
        dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, obj_func_oao_funptr, &
        update_orbs_oao_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the open-shell case
        !
        use opentrustregion, only: update_orbs_type, solver_settings_type
        use otr_oao, only: evaluate_dm_os_type, oao_settings_type
        use otr_oao_test_reference, only: test_evaluate_dm_os_funptr

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed AO "// &
                "density matrix for open-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed AO "// &
                "overlap matrix wrong."
        end if

        ! check number of particles
        if (n_particle /= 2) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= 3) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! test passed density matrix evaluating function
        test_passed = test_passed .and. test_evaluate_dm_os_funptr( &
            evaluate_dm_funptr, "oao_factory_c_wrapper", &
            " by given density matrix evaluating function")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_oao_settings) then
            test_passed = .false.
            write (stderr, *) "test_oao_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_oao_funptr => mock_obj_func_oao
        update_orbs_oao_funptr => mock_update_orbs
        call mock_oao_set_solver_settings(solver_settings)
        dm_ao_3d => dm_ao

    end subroutine mock_oao_factory_os

    subroutine mock_oao_deconstructor()
        !
        ! this subroutine is a test function for the OAO deconstructor
        !
        test_passed = .true.

    end subroutine mock_oao_deconstructor

    function mock_obj_func_oao(kappa, error) result(func)
        !
        ! this function is a test function for the OAO objective function
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: func

        func = sum(kappa)

        error = 0

    end function mock_obj_func_oao

    subroutine mock_precond_oao(residual, mu, precond_residual, error)
        !
        ! this function is a test function for the OAO level-shifted preconditioner
        ! function
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        precond_residual = mu * residual

        error = 0

    end subroutine mock_precond_oao

    subroutine mock_precond_pd_oao(residual, precond_residual, error)
        !
        ! this function is a test function for the OAO positive-definite preconditioner 
        ! function
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        precond_residual = 3.0_rp * residual

        error = 0

    end subroutine mock_precond_pd_oao

    subroutine mock_project_oao(vector, error)
        !
        ! this function is a test function for the OAO projection function
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error

        vector = 2.0_rp * vector

        error = 0

    end subroutine mock_project_oao

    subroutine mock_get_extra_trial_vectors_oao(trial_vectors, error)
        !
        ! this function is a test function for the OAO extra trial vector function
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error

        integer(ip) :: i

        do i = 1, size(trial_vectors, 2)
            trial_vectors(:, i) = real(i, kind=rp)
        end do

        error = 0

    end subroutine mock_get_extra_trial_vectors_oao

end module otr_oao_mock
