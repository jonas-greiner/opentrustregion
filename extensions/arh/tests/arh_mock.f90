! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_mock

    use opentrustregion, only: rp, ip, stderr
    use otr_arh, only: arh_factory_closed_shell, arh_factory_open_shell, &
                       arh_deconstructor
    use test_reference, only: tol
    use otr_arh_test_reference, only: ref_arh_settings

    implicit none

    logical :: test_passed

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory_closed_shell), pointer :: mock_arh_factory_closed_shell_ptr &
        => mock_arh_factory_closed_shell
    procedure(arh_factory_open_shell), pointer :: mock_arh_factory_open_shell_ptr => &
        mock_arh_factory_open_shell
    procedure(arh_deconstructor), pointer :: mock_arh_deconstructor_ptr => &
        mock_arh_deconstructor

contains

    subroutine mock_arh_factory_closed_shell(dm_ao, ao_overlap, n_particle, n_ao, &
                                             get_energy_funptr, update_dm_funptr, &
                                             obj_func_arh_funptr, &
                                             update_orbs_arh_funptr, &
                                             project_arh_funptr, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the closed-shell case
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, hess_x_type, &
                                   project_type
        use otr_oao, only: get_energy_2d_type, update_dm_2d_type
        use otr_arh, only: arh_settings_type
        use otr_oao_test_reference, only: test_get_energy_2d_funptr, &
                                          test_update_dm_2d_funptr, operator(/=)
        use otr_oao_mock, only: mock_obj_func_oao, mock_update_orbs, mock_project_oao, &
                                dm_ao_3d

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_2d_type), intent(in), pointer :: &
            get_energy_funptr
        procedure(update_dm_2d_type), intent(in), pointer :: update_dm_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed AO "// &
                "density matrix for closed-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed AO "// &
                "overlap matrix wrong."
        end if

        ! check number of particles
        if (n_particle /= 1) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= 3) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! test passed energy function
        test_passed = test_passed .and. &
            test_get_energy_2d_funptr(get_energy_funptr, "arh_factory_c_wrapper", &
                                      " by given energy function")

        ! test passed density matrix updating function
        test_passed = test_passed .and. &
            test_update_dm_2d_funptr(update_dm_funptr, "arh_factory_c_wrapper", &
                                     " by given density matrix updating function")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_arh_funptr => mock_obj_func_oao
        update_orbs_arh_funptr => mock_update_orbs
        project_arh_funptr => mock_project_oao
        dm_ao_3d(1:n_ao, 1:n_ao, 1:1) => dm_ao

    end subroutine mock_arh_factory_closed_shell

    subroutine mock_arh_factory_open_shell(dm_ao, ao_overlap, n_particle, n_ao, &
                                           get_energy_funptr, update_dm_spin_funptr, &
                                           obj_func_arh_funptr, &
                                           update_orbs_arh_funptr, project_arh_funptr, &
                                           error, settings)          
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the open-shell case
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, hess_x_type, &
                                   project_type
        use otr_oao, only: get_energy_3d_type
        use otr_arh, only: update_dm_spin_type, arh_settings_type
        use otr_oao_test_reference, only: test_get_energy_3d_funptr, operator(/=)
        use otr_arh_test_reference, only: test_update_dm_spin_funptr
        use otr_oao_mock, only: mock_obj_func_oao, mock_update_orbs, mock_project_oao, &
                                dm_ao_3d

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_3d_type), intent(in), pointer :: get_energy_funptr
        procedure(update_dm_spin_type), intent(in), pointer :: update_dm_spin_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed AO "// &
                "density matrix for open-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed AO "// &
                "overlap matrix wrong."
        end if

        ! check number of particles
        if (n_particle /= 2) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= 3) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! test passed energy function
        test_passed = test_passed .and. &
            test_get_energy_3d_funptr(get_energy_funptr, "arh_factory_c_wrapper", &
                                      " by given energy function")

        ! test passed density matrix updating function
        test_passed = test_passed .and. &
            test_update_dm_spin_funptr(update_dm_spin_funptr, "arh_factory_c_wrapper", &
                                       " by given density matrix updating function "// &
                                       "with same- and opposite-spin potential "// &
                                       "contributions")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_arh_funptr => mock_obj_func_oao
        update_orbs_arh_funptr => mock_update_orbs
        project_arh_funptr => mock_project_oao
        dm_ao_3d => dm_ao

    end subroutine mock_arh_factory_open_shell

    subroutine mock_arh_deconstructor()
        !
        ! this subroutine is a test function for the ARH deconstructor
        !
        test_passed = .true.

    end subroutine mock_arh_deconstructor

end module otr_arh_mock
