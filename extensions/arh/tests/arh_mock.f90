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

    ! error the mock factories return
    integer(ip) :: mock_arh_factory_error = 0

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

    subroutine mock_arh_wire_solver_settings(solver_settings, mo_basis)
        !
        ! this subroutine points the context of the solver settings at the level of the
        ! solver of the ARH object for the given orbital basis, which it allocates if
        ! necessary, and asks for the Hessian refresh and the raised micro iteration
        ! limit, as the ARH factories do
        !
        use opentrustregion, only: solver_settings_type
        use otr_arh, only: arh_object, arh_mo_type, arh_oao_type, arh_n_micro

        type(solver_settings_type), intent(inout), target :: solver_settings
        logical, intent(in) :: mo_basis

        if (.not. allocated(arh_object)) then
            if (mo_basis) then
                allocate(arh_mo_type :: arh_object)
            else
                allocate(arh_oao_type :: arh_object)
            end if
        end if
        arh_object%solver_level%extension => arh_object
        solver_settings%context => arh_object%solver_level
        solver_settings%refresh_hess = .true.
        solver_settings%hess_symm = .false.
        solver_settings%n_micro = arh_n_micro

    end subroutine mock_arh_wire_solver_settings

    subroutine check_factory_mo_input(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, &
                                      n_mo, settings, orbsym)
        !
        ! this subroutine checks the input the mock MO factories are passed by the C
        ! wrapper, which is the same for both spin cases apart from the occupations,
        ! and the irreps of the MOs of every particle channel if they are passed
        !
        use otr_arh, only: arh_settings_type
        use otr_arh_test_reference, only: operator(/=)
        use otr_mo_test_reference, only: n_mo_ref => n_mo, mo_coeff_pattern, &
                                         case_irreps, case_names
        use otr_common_test_reference, only: n_ao_ref => n_ao, n_occ_ref => n_occ
        use otr_mo_mock, only: orbsym_passed

        real(rp), intent(in) :: mo_coeff(:, :, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(in), optional :: orbsym(:, :)

        real(rp), allocatable :: pattern(:, :, :)
        integer(ip) :: i_sym

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
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

        ! record whether irreps of the MOs are passed and check them against those of
        ! the occupation case with symmetry of the shell
        orbsym_passed = present(orbsym)
        if (present(orbsym)) then
            if (n_particle == 1) then
                i_sym = &
                    findloc(case_names, "closed-shell with symmetry", dim=1, kind=ip)
            else
                i_sym = findloc(case_names, "open-shell with symmetry", dim=1, kind=ip)
            end if
            if (any(shape(orbsym) /= [n_mo_ref, n_particle])) then
                test_passed = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed "// &
                    "irreps of the MOs wrongly shaped."
            else if (any(orbsym /= case_irreps(:, :n_particle, i_sym))) then
                test_passed = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Passed "// &
                    "irreps of the MOs wrong."
            end if
        end if

    end subroutine check_factory_mo_input

    subroutine mock_arh_factory_mo_cs(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, &
                                      n_mo, evaluate_dm_funptr, &
                                      update_orbs_arh_funptr, obj_func_arh_funptr, &
                                      solver_settings, error, settings, orbsym)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the closed-shell case with the
        ! orbitals parameterized in the MO basis and wires the solver settings, which
        ! marks the MO coefficients of the caller it receives in place, returns the
        ! error it is set to and otherwise points the context of the solver settings at
        ! the level of the solver of the ARH object and adjusts them, as the factory
        ! does
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_settings_type, evaluate_dm_cs_type
        use otr_arh_test_reference, only: test_evaluate_dm_cs_funptr
        use otr_common_mock, only: mock_update_orbs, mock_obj_func

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ, n_particle, n_ao, n_mo
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(in), optional :: orbsym(:)

        integer(ip), allocatable :: orbsym_2d(:, :)

        ! initialize logical
        test_passed = .true.

        ! check passed input, with the irreps of the MOs reshaped to a single particle
        ! channel
        if (present(orbsym)) orbsym_2d = reshape(orbsym, [size(orbsym, kind=ip), 1_ip])
        call check_factory_mo_input(reshape( &
            mo_coeff, [size(mo_coeff, 1, kind=ip), size(mo_coeff, 2, kind=ip), 1_ip]), &
            ao_overlap, [n_occ], n_particle, n_ao, n_mo, settings, orbsym_2d)

        ! test passed density matrix evaluating function
        if (.not. test_evaluate_dm_cs_funptr( &
            evaluate_dm_funptr, "arh_factory_mo_c_wrapper", " by given density "// &
            "matrix evaluating function with Coulomb, exact-exchange and "// &
            "non-linear potential contributions", settings%context)) &
            test_passed = .false.

        ! mark the MO coefficients of the caller and set output quantities
        mo_coeff = 2.0_rp
        error = mock_arh_factory_error
        update_orbs_arh_funptr => null()
        obj_func_arh_funptr => null()
        if (error /= 0) return
        call mock_arh_wire_solver_settings(solver_settings, .true.)
        update_orbs_arh_funptr => mock_update_orbs
        obj_func_arh_funptr => mock_obj_func

    end subroutine mock_arh_factory_mo_cs

    subroutine mock_arh_factory_mo_os(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, &
                                      n_mo, evaluate_dm_funptr, &
                                      update_orbs_arh_funptr, obj_func_arh_funptr, &
                                      solver_settings, error, settings, orbsym)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the open-shell case with the
        ! orbitals parameterized in the MO basis and wires the solver settings, which
        ! marks the MO coefficients of the caller it receives in place, returns the
        ! error it is set to and otherwise points the context of the solver settings at
        ! the level of the solver of the ARH object and adjusts them, as the factory
        ! does
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_settings_type, evaluate_dm_os_type
        use otr_arh_test_reference, only: test_evaluate_dm_os_funptr
        use otr_common_mock, only: mock_update_orbs, mock_obj_func

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(in), optional :: orbsym(:, :)

        ! initialize logical
        test_passed = .true.

        ! check passed input
        call check_factory_mo_input(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, &
                                    n_mo, settings, orbsym)

        ! test passed density matrix evaluating function
        if (.not. test_evaluate_dm_os_funptr( &
            evaluate_dm_funptr, "arh_factory_mo_c_wrapper", " by given density "// &
            "matrix evaluating function with Coulomb, exact-exchange and "// &
            "non-linear potential contributions", settings%context)) &
            test_passed = .false.

        ! mark the MO coefficients of the caller and set output quantities
        mo_coeff = 2.0_rp
        error = mock_arh_factory_error
        update_orbs_arh_funptr => null()
        obj_func_arh_funptr => null()
        if (error /= 0) return
        call mock_arh_wire_solver_settings(solver_settings, .true.)
        update_orbs_arh_funptr => mock_update_orbs
        obj_func_arh_funptr => mock_obj_func

    end subroutine mock_arh_factory_mo_os

    subroutine mock_arh_factory_oao_cs( &
        dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
        update_orbs_arh_funptr, obj_func_arh_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the closed-shell case with the
        ! orbitals parameterized in the OAO basis and wires the solver settings, which
        ! marks the density matrix of the caller it receives in place, returns the
        ! error it is set to and otherwise points the context of the solver settings at
        ! the level of the solver of the ARH object and adjusts them, as the factory
        ! does
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_settings_type, evaluate_dm_cs_type
        use otr_arh_test_reference, only: test_evaluate_dm_cs_funptr, operator(/=)
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_common_test_reference, only: n_ao_ref => n_ao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

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
        if (.not. test_evaluate_dm_cs_funptr( &
            evaluate_dm_funptr, "arh_factory_oao_c_wrapper", " by given density "// &
            "matrix evaluating function with Coulomb, exact-exchange and "// &
            "non-linear potential contributions", settings%context)) &
            test_passed = .false.

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed "// &
                "logging function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

        ! mark the density matrix of the caller and set output quantities
        dm_ao = 2.0_rp
        error = mock_arh_factory_error
        update_orbs_arh_funptr => null()
        obj_func_arh_funptr => null()
        if (error /= 0) return
        call mock_arh_wire_solver_settings(solver_settings, .false.)
        update_orbs_arh_funptr => mock_update_orbs
        obj_func_arh_funptr => mock_obj_func

    end subroutine mock_arh_factory_oao_cs

    subroutine mock_arh_factory_oao_os( &
        dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_os_funptr, &
        update_orbs_arh_funptr, obj_func_arh_funptr, solver_settings, error, settings)
        !
        ! this function is a test function for the function which returns the orbital
        ! updating function and objective function for the open-shell case with the
        ! orbitals parameterized in the OAO basis and wires the solver settings, which
        ! marks the density matrix of the caller it receives in place, returns the
        ! error it is set to and otherwise points the context of the solver settings at
        ! the level of the solver of the ARH object and adjusts them, as the factory
        ! does
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: evaluate_dm_os_type, arh_settings_type
        use otr_arh_test_reference, only: test_evaluate_dm_os_funptr, operator(/=)
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_common_test_reference, only: n_ao_ref => n_ao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        type(solver_settings_type), intent(inout), target :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

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
        if (.not. test_evaluate_dm_os_funptr( &
            evaluate_dm_os_funptr, "arh_factory_oao_c_wrapper", " by given density "// &
            "matrix evaluating function with Coulomb, exact-exchange and "// &
            "non-linear potential contributions", settings%context)) &
            test_passed = .false.

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed "// &
                "logging function not associated with value."
        else
            call settings%logger("test", settings%context)
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

        ! mark the density matrix of the caller and set output quantities
        dm_ao = 2.0_rp
        error = mock_arh_factory_error
        update_orbs_arh_funptr => null()
        obj_func_arh_funptr => null()
        if (error /= 0) return
        call mock_arh_wire_solver_settings(solver_settings, .false.)
        update_orbs_arh_funptr => mock_update_orbs
        obj_func_arh_funptr => mock_obj_func

    end subroutine mock_arh_factory_oao_os

    subroutine mock_arh_deconstructor()
        !
        ! this subroutine is a test function for the ARH deconstructor
        !
        test_passed = .true.

    end subroutine mock_arh_deconstructor

end module otr_arh_mock
