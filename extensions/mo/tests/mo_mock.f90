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

    ! error the mock factories return
    integer(ip) :: mock_mo_factory_error = 0

    ! whether the last call of a mock factory for orbitals parameterized in the MO
    ! basis was passed irreps of the MOs
    logical :: orbsym_passed = .false.

    ! create function pointers to ensure that routines comply with interface
    procedure(mo_factory_cs), pointer :: mock_mo_factory_cs_ptr => mock_mo_factory_cs
    procedure(mo_factory_os), pointer :: mock_mo_factory_os_ptr => mock_mo_factory_os
    procedure(mo_deconstructor), pointer :: mock_mo_deconstructor_ptr => &
        mock_mo_deconstructor

contains

    subroutine mock_mo_factory_cs(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                                  evaluate_dm_funptr, update_orbs_mo_funptr, &
                                  obj_func_mo_funptr, solver_settings, error, &
                                  settings, orbsym)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the closed-shell case and wires
        ! the solver settings, which marks the MO coefficients of the caller it
        ! receives in place, returns the error it is set to and otherwise points the
        ! context of the solver settings at the level of the solver of the MO object,
        ! as the factory does
        !
        use opentrustregion, only: update_orbs_type, solver_settings_type
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_mo, only: mo_settings_type, mo_object
        use otr_common, only: evaluate_dm_cs_type
        use otr_common_test_reference, only: test_evaluate_dm_cs_funptr, &
                                             n_ao_ref => n_ao, n_occ_ref => n_occ
        use otr_mo_test_reference, only: n_mo_ref => n_mo, mo_coeff_pattern, &
                                         case_irreps, case_names

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ, n_particle, n_ao, n_mo
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_mo_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(mo_settings_type), intent(in) :: settings
        integer(ip), intent(in), optional :: orbsym(:)

        real(rp), allocatable :: pattern(:, :, :)
        integer(ip) :: i_sym

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
        if (.not. test_evaluate_dm_cs_funptr( &
            evaluate_dm_funptr, "mo_factory_c_wrapper", " by given density matrix "// &
            "evaluating function", settings%context)) test_passed = .false.

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_mo_settings) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! record whether irreps of the MOs are passed and check them against those of
        ! the closed-shell occupation case with symmetry
        orbsym_passed = present(orbsym)
        if (present(orbsym)) then
            i_sym = findloc(case_names, "closed-shell with symmetry", dim=1, kind=ip)
            if (any(shape(orbsym) /= [n_mo_ref])) then
                test_passed = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Passed irreps "// &
                    "of the MOs wrongly shaped."
            else if (any(orbsym /= case_irreps(:, 1, i_sym))) then
                test_passed = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Passed irreps "// &
                    "of the MOs wrong."
            end if
        end if

        ! mark the MO coefficients of the caller and set output quantities
        mo_coeff = 2.0_rp
        error = mock_mo_factory_error
        update_orbs_mo_funptr => null()
        obj_func_mo_funptr => null()
        if (error /= 0) return
        if (.not. allocated(mo_object)) allocate(mo_object)
        mo_object%solver_level%extension => mo_object
        solver_settings%context => mo_object%solver_level
        update_orbs_mo_funptr => mock_update_orbs
        obj_func_mo_funptr => mock_obj_func

    end subroutine mock_mo_factory_cs

    subroutine mock_mo_factory_os(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                                  evaluate_dm_funptr, update_orbs_mo_funptr, &
                                  obj_func_mo_funptr, solver_settings, error, &
                                  settings, orbsym)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the open-shell case and wires the
        ! solver settings, which marks the MO coefficients of the caller it receives in
        ! place, returns the error it is set to and otherwise points the context of the
        ! solver settings at the level of the solver of the MO object, as the factory
        ! does
        !
        use opentrustregion, only: update_orbs_type, solver_settings_type
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_mo, only: mo_settings_type, mo_object
        use otr_common, only: evaluate_dm_os_type
        use otr_common_test_reference, only: test_evaluate_dm_os_funptr, &
                                             n_ao_ref => n_ao, n_occ_ref => n_occ
        use otr_mo_test_reference, only: n_mo_ref => n_mo, mo_coeff_pattern, &
                                         case_irreps, case_names

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_mo_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(mo_settings_type), intent(in) :: settings
        integer(ip), intent(in), optional :: orbsym(:, :)

        real(rp), allocatable :: pattern(:, :, :)
        integer(ip) :: i_sym

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
        if (.not. test_evaluate_dm_os_funptr( &
            evaluate_dm_funptr, "mo_factory_c_wrapper", " by given density matrix "// &
            "evaluating function", settings%context)) test_passed = .false.

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_mo_settings) then
            test_passed = .false.
            write(stderr, *) "test_mo_factory_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

        ! record whether irreps of the MOs are passed and check them against those of
        ! the open-shell occupation case with symmetry
        orbsym_passed = present(orbsym)
        if (present(orbsym)) then
            i_sym = findloc(case_names, "open-shell with symmetry", dim=1, kind=ip)
            if (any(shape(orbsym) /= [n_mo_ref, 2_ip])) then
                test_passed = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Passed irreps "// &
                    "of the MOs wrongly shaped."
            else if (any(orbsym /= case_irreps(:, :, i_sym))) then
                test_passed = .false.
                write(stderr, *) "test_mo_factory_c_wrapper failed: Passed irreps "// &
                    "of the MOs wrong."
            end if
        end if

        ! mark the MO coefficients of the caller and set output quantities
        mo_coeff = 2.0_rp
        error = mock_mo_factory_error
        update_orbs_mo_funptr => null()
        obj_func_mo_funptr => null()
        if (error /= 0) return
        if (.not. allocated(mo_object)) allocate(mo_object)
        mo_object%solver_level%extension => mo_object
        solver_settings%context => mo_object%solver_level
        update_orbs_mo_funptr => mock_update_orbs
        obj_func_mo_funptr => mock_obj_func

    end subroutine mock_mo_factory_os

    subroutine mock_mo_deconstructor()
        !
        ! this subroutine is a test function for the MO deconstructor
        !
        test_passed = .true.

    end subroutine mock_mo_deconstructor

end module otr_mo_mock
