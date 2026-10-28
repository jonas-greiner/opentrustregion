! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_mock

    use opentrustregion_core, only: rp, ip, stderr, obj_func_type
    use c_interface, only: c_ip, c_rp
    use otr_oao, only: oao_factory_cs, oao_factory_os
    use test_reference, only: tol
    use otr_oao_test_reference, only: ref_oao_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! error the mock factories return
    integer(ip) :: mock_oao_factory_error = 0

    ! whether the mock factories expect the default settings without a logging
    ! function, which the C wrapper hands them for C settings that were not initialized
    ! and without any, instead of the reference settings, and whether they expect no
    ! density matrix evaluating function
    logical :: mock_oao_factory_defaults = .false., &
               mock_oao_factory_no_evaluate_dm = .false.

    ! create function pointers to ensure that routines comply with interface
    procedure(oao_factory_cs), pointer :: mock_oao_factory_cs_ptr => mock_oao_factory_cs
    procedure(oao_factory_os), pointer :: mock_oao_factory_os_ptr => mock_oao_factory_os

contains

    subroutine check_and_wire_oao_object(ao_overlap, n_ao, update_orbs_oao_funptr, &
                                         obj_func_oao_funptr, solver_settings, &
                                         solver_state, error, settings)
        !
        ! this subroutine checks the inputs both mock factories share, returns the
        ! error the mock factories are set to and otherwise creates an OAO object
        ! reading the given state, points the context of the solver settings at its
        ! level of the solver and lets the settings own it, as the factory does
        !
        use opentrustregion_core, only: update_orbs_type, solver_state_type
        use otr_common, only: extension_type, extension_solver_settings_type
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_common_unit_tests, only: check_mock_factory_logger, &
                                         check_mock_factory_ao_overlap, &
                                         setup_owning_solver_settings
        use otr_oao, only: oao_type, oao_settings_type, default_oao_settings

        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_ao
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(in), optional :: settings

        class(extension_type), allocatable, target :: object

        ! check shape and values of passed overlap matrix and number of AOs
        if (.not. check_mock_factory_ao_overlap( &
            ao_overlap, n_ao, "test_oao_factory_c_wrapper")) test_passed = .false.

        ! the C wrapper always passes settings, since their context points at the
        ! bundle of the host's callback functions
        update_orbs_oao_funptr => null()
        obj_func_oao_funptr => null()
        if (.not. present(settings)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Settings not passed."
            error = 1
            return
        end if

        ! check if optional logging function is correctly passed
        if (.not. check_mock_factory_logger(settings, mock_oao_factory_defaults, &
                                            "test_oao_factory_c_wrapper")) &
            test_passed = .false.

        ! check if optional settings are correctly passed
        if (mock_oao_factory_defaults) then
            if (settings /= default_oao_settings) then
                test_passed = .false.
                write(stderr, *) "test_oao_factory_c_wrapper failed: Settings that "// &
                    "were not initialized not initialized to the default values."
            end if
        else if (settings /= ref_oao_settings) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! set output quantities
        error = mock_oao_factory_error
        if (error /= 0) return
        allocate(oao_type :: object)
        object%state => solver_state
        call setup_owning_solver_settings(solver_settings, object)
        update_orbs_oao_funptr => mock_update_orbs
        obj_func_oao_funptr => mock_obj_func

    end subroutine check_and_wire_oao_object

    subroutine mock_oao_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, &
                                   evaluate_dm_funptr, update_orbs_oao_funptr, &
                                   obj_func_oao_funptr, solver_settings, solver_state, &
                                   error, settings)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the closed-shell case and wires
        ! the solver settings, which marks the density matrix of the caller it receives
        ! in place and otherwise behaves like the factory
        !
        use opentrustregion_core, only: update_orbs_type, solver_state_type
        use otr_common, only: evaluate_dm_cs_type, extension_solver_settings_type
        use otr_oao, only: oao_settings_type
        use otr_common_test_reference, only: check_evaluate_dm_cs_funptr, &
                                             n_ao_ref => n_ao
        use otr_oao_test_reference, only: density_pattern

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(in), optional :: settings

        ! initialize logical
        test_passed = .true.

        ! check shape and values of passed density matrix and number of particles
        if (any(shape(dm_ao) /= [n_ao_ref, n_ao_ref])) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed AO density "// &
                "matrix for closed-shell case shaped wrong."
        else if (any(abs(dm_ao - reshape( &
            real(density_pattern(int(n_ao_ref, kind=c_ip), 1_c_ip, 0.0_c_rp), &
                 kind=rp), [n_ao_ref, n_ao_ref])) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed AO density "// &
                "matrix for closed-shell case wrong."
        end if
        if (n_particle /= 1) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! test passed density matrix evaluating function, if one is expected
        if (mock_oao_factory_no_evaluate_dm) then
            if (associated(evaluate_dm_funptr)) then
                test_passed = .false.
                write(stderr, *) "test_oao_factory_c_wrapper failed: Density "// &
                    "matrix evaluating function passed without one."
            end if
        else if (present(settings)) then
            if (.not. check_evaluate_dm_cs_funptr( &
                evaluate_dm_funptr, "oao_factory_c_wrapper", " by given density "// &
                "matrix evaluating function", settings%context)) test_passed = .false.
        end if

        ! mark the density matrix of the caller, check the inputs both mock factories
        ! share and wire an object, as the factory does
        dm_ao = 2.0_rp
        call check_and_wire_oao_object(ao_overlap, n_ao, update_orbs_oao_funptr, &
                                       obj_func_oao_funptr, solver_settings, &
                                       solver_state, error, settings)

    end subroutine mock_oao_factory_cs

    subroutine mock_oao_factory_os(dm_ao, ao_overlap, n_particle, n_ao, &
                                   evaluate_dm_funptr, update_orbs_oao_funptr, &
                                   obj_func_oao_funptr, solver_settings, solver_state, &
                                   error, settings)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the open-shell case and wires
        ! the solver settings, which marks the density matrix of the caller it receives
        ! in place and otherwise behaves like the factory
        !
        use opentrustregion_core, only: update_orbs_type, solver_state_type
        use otr_common, only: evaluate_dm_os_type, extension_solver_settings_type
        use otr_oao, only: oao_settings_type
        use otr_common_test_reference, only: check_evaluate_dm_os_funptr, &
                                             n_ao_ref => n_ao
        use otr_oao_test_reference, only: density_pattern

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(in), optional :: settings

        ! initialize logical
        test_passed = .true.

        ! check shape and values of passed density matrix and number of particles
        if (any(shape(dm_ao) /= [n_ao_ref, n_ao_ref, 2_ip])) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed AO density "// &
                "matrix for open-shell case shaped wrong."
        else if (any(abs(dm_ao - real(density_pattern( &
            int(n_ao_ref, kind=c_ip), 2_c_ip, 0.0_c_rp), kind=rp)) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed AO density "// &
                "matrix for open-shell case wrong."
        end if
        if (n_particle /= 2) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! test passed density matrix evaluating function, if one is expected
        if (mock_oao_factory_no_evaluate_dm) then
            if (associated(evaluate_dm_funptr)) then
                test_passed = .false.
                write(stderr, *) "test_oao_factory_c_wrapper failed: Density "// &
                    "matrix evaluating function passed without one."
            end if
        else if (present(settings)) then
            if (.not. check_evaluate_dm_os_funptr( &
                evaluate_dm_funptr, "oao_factory_c_wrapper", " by given density "// &
                "matrix evaluating function", settings%context)) test_passed = .false.
        end if

        ! mark the density matrix of the caller, check the inputs both mock factories
        ! share and wire an object, as the factory does
        dm_ao = 2.0_rp
        call check_and_wire_oao_object(ao_overlap, n_ao, update_orbs_oao_funptr, &
                                       obj_func_oao_funptr, solver_settings, &
                                       solver_state, error, settings)

    end subroutine mock_oao_factory_os

end module otr_oao_mock
