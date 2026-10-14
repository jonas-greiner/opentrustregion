! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_mock

    use opentrustregion, only: ip, stderr
    use otr_arh, only: arh_factory_common
    use otr_common, only: extension_level_type
    use otr_arh_test_reference, only: ref_arh_settings

    implicit none

    logical :: test_passed

    ! error the mock factory returns
    integer(ip) :: mock_arh_factory_error = 0

    ! whether the mock factory expects the default settings without a logging function,
    ! which the C wrapper hands it for C settings that were not initialized and without
    ! any, instead of the reference settings
    logical :: mock_arh_factory_defaults = .false.

    ! level of the orbital basis object the mock factory received
    type(extension_level_type), pointer :: received_arh_level => null()

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory_common), pointer :: mock_arh_factory_common_ptr => &
        mock_arh_factory_common

contains

    subroutine mock_arh_factory_common( &
        level, evaluate_dm_cs_funptr, evaluate_dm_os_funptr, update_orbs_arh_funptr, &
        obj_func_arh_funptr, solver_settings, solver_state, error, settings)
        !
        ! this function is a test function for the common setup of the ARH factories,
        ! which records the level it receives, checks that it receives the density
        ! matrix evaluating function of the shell of the orbitals of an orbital basis
        ! object of the level, returns the error it is set to and otherwise creates an
        ! ARH object reading the given state, points the context of the solver settings
        ! at its level of the solver, lets the settings own it and adjusts them, as the
        ! factory does
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_state_type
        use otr_arh, only: arh_settings_type, arh_mo_type, default_arh_settings, &
                           arh_evaluate_dm_cs_type, arh_evaluate_dm_os_type, arh_n_micro
        use otr_common, only: extension_type, extension_solver_settings_type, &
                              orbital_basis_type
        use otr_arh_test_reference, only: check_arh_evaluate_dm_cs_funptr, &
                                          check_arh_evaluate_dm_os_funptr, operator(/=)
        use otr_common_mock, only: mock_update_orbs, mock_obj_func

        type(extension_level_type), intent(in), pointer :: level
        procedure(arh_evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_cs_funptr
        procedure(arh_evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in), optional :: settings

        class(extension_type), allocatable, target :: object
        integer(ip) :: n_particle

        ! initialize logical and record the level
        test_passed = .true.
        received_arh_level => level
        update_orbs_arh_funptr => null()
        obj_func_arh_funptr => null()

        ! the C wrapper always passes settings, since their context points at the
        ! bundle of the host's callback functions
        if (.not. present(settings)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_c_wrapper failed: Settings not passed."
            error = 1
            return
        end if

        ! test passed density matrix evaluating function of the shell of the orbitals
        n_particle = 0
        if (associated(level)) then
            select type (orbitals => level%extension)
            class is (orbital_basis_type)
                n_particle = orbitals%n_particle
            end select
        end if
        if (n_particle == 1) then
            if (associated(evaluate_dm_os_funptr) .or. &
                .not. check_arh_evaluate_dm_cs_funptr( &
                    evaluate_dm_cs_funptr, "arh_factory_c_wrapper", &
                    " by given closed-shell density matrix evaluating function", &
                    settings%context)) test_passed = .false.
        else if (n_particle == 2) then
            if (associated(evaluate_dm_cs_funptr) .or. &
                .not. check_arh_evaluate_dm_os_funptr( &
                    evaluate_dm_os_funptr, "arh_factory_c_wrapper", &
                    " by given open-shell density matrix evaluating function", &
                    settings%context)) test_passed = .false.
        end if

        ! check if optional logging function is correctly passed
        if (mock_arh_factory_defaults) then
            if (associated(settings%logger)) then
                test_passed = .false.
                write(stderr, *) "test_arh_factory_c_wrapper failed: Logging "// &
                    "function passed for settings that were not initialized."
            end if
        else if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (mock_arh_factory_defaults) then
            if (settings /= default_arh_settings) then
                test_passed = .false.
                write(stderr, *) "test_arh_factory_c_wrapper failed: Settings that "// &
                    "were not initialized not initialized to the default values."
            end if
        else if (settings /= ref_arh_settings) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = mock_arh_factory_error
        if (error /= 0) return
        allocate(arh_mo_type :: object)
        object%state => solver_state
        object%solver_level%extension => object
        solver_settings%context => object%solver_level
        allocate(solver_settings%extensions(1))
        call move_alloc(object, solver_settings%extensions(1)%object)
        solver_settings%refresh_hess = .true.
        solver_settings%hess_symm = .false.
        solver_settings%n_micro = arh_n_micro
        update_orbs_arh_funptr => mock_update_orbs
        obj_func_arh_funptr => mock_obj_func

    end subroutine mock_arh_factory_common

end module otr_arh_mock
