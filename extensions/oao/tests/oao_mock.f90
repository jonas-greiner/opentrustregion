! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_mock

    use opentrustregion, only: rp, ip, stderr, obj_func_type
    use otr_oao, only: oao_factory_cs, oao_factory_os, oao_deconstructor
    use test_reference, only: tol
    use otr_oao_test_reference, only: ref_oao_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! error the mock factories return
    integer(ip) :: mock_oao_factory_error = 0

    ! create function pointers to ensure that routines comply with interface
    procedure(oao_factory_cs), pointer :: mock_oao_factory_cs_ptr => mock_oao_factory_cs
    procedure(oao_factory_os), pointer :: mock_oao_factory_os_ptr => mock_oao_factory_os
    procedure(oao_deconstructor), pointer :: mock_oao_deconstructor_ptr => &
        mock_oao_deconstructor

contains

    subroutine mock_oao_factory_cs( &
        dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
        update_orbs_oao_funptr, obj_func_oao_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the closed-shell case and wires
        ! the solver settings, which marks the density matrix of the caller it receives
        ! in place, returns the error it is set to and otherwise points the context of
        ! the solver settings at the level of the solver of the OAO object, as the
        ! factory does
        !
        use opentrustregion, only: update_orbs_type, solver_settings_type
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_oao, only: oao_settings_type, oao_object
        use otr_common, only: evaluate_dm_cs_type
        use otr_common_test_reference, only: test_evaluate_dm_cs_funptr, &
                                             n_ao_ref => n_ao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(in) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed AO density "// &
                "matrix for closed-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed AO overlap "// &
                "matrix wrong."
        end if

        ! check number of particles
        if (n_particle /= 1) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= n_ao_ref) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! test passed density matrix evaluating function
        if (.not. test_evaluate_dm_cs_funptr( &
            evaluate_dm_funptr, "oao_factory_c_wrapper", " by given density matrix "// &
            "evaluating function", settings%context)) test_passed = .false.

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_oao_settings) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! mark the density matrix of the caller and set output quantities
        dm_ao = 2.0_rp
        error = mock_oao_factory_error
        update_orbs_oao_funptr => null()
        obj_func_oao_funptr => null()
        if (error /= 0) return
        if (.not. allocated(oao_object)) allocate(oao_object)
        oao_object%solver_level%extension => oao_object
        solver_settings%context => oao_object%solver_level
        update_orbs_oao_funptr => mock_update_orbs
        obj_func_oao_funptr => mock_obj_func

    end subroutine mock_oao_factory_cs

    subroutine mock_oao_factory_os( &
        dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
        update_orbs_oao_funptr, obj_func_oao_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the open-shell case and wires the
        ! solver settings, which marks the density matrix of the caller it receives in
        ! place, returns the error it is set to and otherwise points the context of the
        ! solver settings at the level of the solver of the OAO object, as the factory
        ! does
        !
        use opentrustregion, only: update_orbs_type, solver_settings_type
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_oao, only: oao_settings_type, oao_object
        use otr_common, only: evaluate_dm_os_type
        use otr_common_test_reference, only: test_evaluate_dm_os_funptr, &
                                             n_ao_ref => n_ao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(in) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed AO density "// &
                "matrix for open-shell case wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed AO overlap "// &
                "matrix wrong."
        end if

        ! check number of particles
        if (n_particle /= 2) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "particles wrong."
        end if

        ! check number of AOs
        if (n_ao /= n_ao_ref) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! test passed density matrix evaluating function
        if (.not. test_evaluate_dm_os_funptr( &
            evaluate_dm_funptr, "oao_factory_c_wrapper", " by given density matrix "// &
            "evaluating function", settings%context)) test_passed = .false.

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_oao_settings) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! mark the density matrix of the caller and set output quantities
        dm_ao = 2.0_rp
        error = mock_oao_factory_error
        update_orbs_oao_funptr => null()
        obj_func_oao_funptr => null()
        if (error /= 0) return
        if (.not. allocated(oao_object)) allocate(oao_object)
        oao_object%solver_level%extension => oao_object
        solver_settings%context => oao_object%solver_level
        update_orbs_oao_funptr => mock_update_orbs
        obj_func_oao_funptr => mock_obj_func

    end subroutine mock_oao_factory_os

    subroutine mock_oao_deconstructor()
        !
        ! this subroutine is a test function for the OAO deconstructor
        !
        test_passed = .true.

    end subroutine mock_oao_deconstructor

end module otr_oao_mock
