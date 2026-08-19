! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_mock

    use opentrustregion, only: rp, ip, stderr, obj_func_type, project_type
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
    procedure(project_type), pointer ::  mock_project_oao_ptr => mock_project_oao

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

    subroutine mock_oao_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, &
                                   get_energy_funptr, update_dm_funptr, &
                                   obj_func_oao_funptr, update_orbs_oao_funptr, &
                                   project_oao_funptr, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the closed-shell case
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_oao, only: get_energy_cs_type, update_dm_cs_type, oao_settings_type
        use otr_oao_test_reference, only: test_get_energy_cs_funptr, &
                                          test_update_dm_cs_funptr

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_cs_type), intent(in), pointer :: get_energy_funptr
        procedure(update_dm_cs_type), intent(in), pointer :: update_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(project_type), intent(out), pointer :: project_oao_funptr
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

        ! test passed energy function
        test_passed = test_passed .and. &
            test_get_energy_cs_funptr(get_energy_funptr, "oao_factory_c_wrapper", &
                                      " by given energy function")

        ! test passed density matrix updating function
        test_passed = test_passed .and. &
            test_update_dm_cs_funptr(update_dm_funptr, "oao_factory_c_wrapper", &
                                     " by given density matrix updating function")

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
        project_oao_funptr => mock_project_oao
        dm_ao_3d(1:n_ao, 1:n_ao, 1:1) => dm_ao

    end subroutine mock_oao_factory_cs

    subroutine mock_oao_factory_os(dm_ao, ao_overlap, n_particle, n_ao, &
                                   get_energy_funptr, update_dm_funptr, &
                                   obj_func_oao_funptr, update_orbs_oao_funptr, &
                                   project_oao_funptr, error, settings)          
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the open-shell case
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_oao, only: get_energy_os_type, update_dm_os_type, oao_settings_type
        use otr_oao_test_reference, only: test_get_energy_os_funptr, &
                                          test_update_dm_os_funptr

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_os_type), intent(in), pointer :: get_energy_funptr
        procedure(update_dm_os_type), intent(in), pointer :: update_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(project_type), intent(out), pointer :: project_oao_funptr
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

        ! test passed energy function
        test_passed = test_passed .and. &
            test_get_energy_os_funptr(get_energy_funptr, "oao_factory_c_wrapper", &
                                      " by given energy function")

        ! test passed density matrix updating function
        test_passed = test_passed .and. &
            test_update_dm_os_funptr(update_dm_funptr, "oao_factory_c_wrapper", &
                                     " by given density matrix updating function")

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
        project_oao_funptr => mock_project_oao
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

    subroutine mock_project_oao(vector, error)
        !
        ! this function is a test function for the OAO projection function
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error

        vector = 2.0_rp * vector

        error = 0

    end subroutine mock_project_oao

end module otr_oao_mock
