! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_mock

    use opentrustregion, only: rp, ip, stderr, obj_func_type, project_type
    use otr_arh, only: arh_factory_closed_shell, arh_factory_open_shell, &
                       arh_deconstructor_closed_shell, arh_deconstructor_open_shell
    use test_reference, only: tol
    use otr_arh_test_reference, only: ref_arh_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory_closed_shell), pointer :: mock_arh_factory_closed_shell_ptr &
        => mock_arh_factory_closed_shell
    procedure(arh_factory_open_shell), pointer :: mock_arh_factory_open_shell_ptr => &
        mock_arh_factory_open_shell
    procedure(arh_deconstructor_closed_shell), pointer :: &
        mock_arh_deconstructor_closed_shell_ptr => mock_arh_deconstructor_closed_shell
    procedure(arh_deconstructor_open_shell), pointer :: &
        mock_arh_deconstructor_open_shell_ptr => mock_arh_deconstructor_open_shell
    procedure(obj_func_type), pointer :: mock_obj_func_arh_ptr => mock_obj_func_arh
    procedure(project_type), pointer ::  mock_project_arh_ptr => mock_project_arh

contains

    subroutine mock_arh_factory_closed_shell(dm_ao, ao_overlap, n_particle, n_ao, &
                                             get_energy_funptr, get_fock_funptr, &
                                             obj_func_arh_funptr, &
                                             update_orbs_arh_funptr, &
                                             project_arh_funptr, error, settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the closed-shell case
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_arh, only: get_energy_closed_shell_type, get_fock_type, &
                           arh_settings_type
        use otr_arh_test_reference, only: test_get_energy_closed_shell_funptr, &
                                          test_get_fock_funptr
        use otr_common_mock, only: mock_update_orbs

        real(rp), intent(in) :: dm_ao(:, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_closed_shell_type), intent(in), pointer :: &
            get_energy_funptr
        procedure(get_fock_type), intent(in), pointer :: get_fock_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

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
            test_get_energy_closed_shell_funptr(get_energy_funptr, &
                                                "arh_factory_c_wrapper", &
                                                " by given energy function")

        ! test passed Fock matrix function
        test_passed = test_passed .and. &
            test_get_fock_funptr(get_fock_funptr, "arh_factory_c_wrapper", &
                                 " by given Fock matrix function")

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
        obj_func_arh_funptr => mock_obj_func_arh
        update_orbs_arh_funptr => mock_update_orbs
        project_arh_funptr => mock_project_arh

    end subroutine mock_arh_factory_closed_shell

    subroutine mock_arh_factory_open_shell(dm_ao, ao_overlap, n_particle, n_ao, &
                                           get_energy_funptr, get_fock_jk_funptr, &
                                           obj_func_arh_funptr, &
                                           update_orbs_arh_funptr, project_arh_funptr, &
                                           error, settings)          
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function for the open-shell case
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_arh, only: get_energy_open_shell_type, get_fock_jk_type, &
                           arh_settings_type
        use otr_arh_test_reference, only: test_get_energy_open_shell_funptr, &
                                          test_get_fock_jk_funptr
        use otr_common_mock, only: mock_update_orbs

        real(rp), intent(in) :: dm_ao(:, :, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_open_shell_type), intent(in), pointer :: get_energy_funptr
        procedure(get_fock_jk_type), intent(in), pointer :: get_fock_jk_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

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
            test_get_energy_open_shell_funptr(get_energy_funptr, &
                                              "arh_factory_c_wrapper", &
                                              " by given energy function")

        ! test passed Fock, Coulomb, and exchange matrix function
        test_passed = test_passed .and. &
            test_get_fock_jk_funptr(get_fock_jk_funptr, "arh_factory_c_wrapper", &
                                    " by given Fock, Coulomb, and exchange matrix "// &
                                    "function")

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
        obj_func_arh_funptr => mock_obj_func_arh
        update_orbs_arh_funptr => mock_update_orbs
        project_arh_funptr => mock_project_arh

    end subroutine mock_arh_factory_open_shell

    subroutine mock_arh_deconstructor_closed_shell(dm_ao, error)
        !
        ! this subroutine is a test function for the ARH deconstructor for the 
        ! closed-shell case
        !
        real(rp), intent(out) :: dm_ao(:, :)
        integer(ip), intent(out) :: error

        dm_ao = 1.0_rp

        error = 0

    end subroutine mock_arh_deconstructor_closed_shell

    subroutine mock_arh_deconstructor_open_shell(dm_ao, error)
        !
        ! this subroutine is a test function for the ARH deconstructor for the 
        ! open-shell case
        !
        real(rp), intent(out) :: dm_ao(:, :, :)
        integer(ip), intent(out) :: error

        dm_ao = 1.0_rp

        error = 0

    end subroutine mock_arh_deconstructor_open_shell

    function mock_obj_func_arh(kappa, error) result(func)
        !
        ! this function is a test function for the ARH C objective function
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: func

        func = sum(kappa)

        error = 0

    end function mock_obj_func_arh

    subroutine mock_project_arh(vector, error)
        !
        ! this function is a test function for the ARH C projection function
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error

        vector = 2.0_rp * vector

        error = 0

    end subroutine mock_project_arh

end module otr_arh_mock
