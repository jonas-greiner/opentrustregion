! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_mock

    use opentrustregion, only: ip, stderr
    use c_interface, only: c_rp, c_ip
    use otr_arh_test_reference, only: ref_arh_settings
    use otr_arh_c_interface, only: arh_factory_mo_c_wrapper, &
                                   arh_factory_oao_c_wrapper, init_arh_settings_c, &
                                   arh_deconstructor_c_wrapper
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char, c_f_pointer, c_loc

    implicit none

    logical(c_bool), bind(C) :: test_arh_factory_mo_interface = .true._c_bool, &
                                test_arh_factory_oao_interface = .true._c_bool, &
                                test_arh_deconstructor_interface = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory_mo_c_wrapper), pointer :: &
        mock_arh_factory_mo_c_wrapper_ptr => mock_arh_factory_mo_c_wrapper
    procedure(arh_factory_oao_c_wrapper), pointer :: &
        mock_arh_factory_oao_c_wrapper_ptr => mock_arh_factory_oao_c_wrapper
    procedure(init_arh_settings_c), pointer :: mock_init_arh_settings_c_ptr => &
        mock_init_arh_settings_c
    procedure(arh_deconstructor_c_wrapper), pointer :: &
        mock_arh_deconstructor_c_wrapper_ptr => mock_arh_deconstructor_c_wrapper

contains

    function mock_arh_factory_mo_c_wrapper( &
        mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
        evaluate_dm_c_funptr, update_orbs_arh_c_funptr, obj_func_arh_c_funptr, &
        solver_settings_c, settings_c, orbsym_c) result(error_c) &
        bind(C, name="mock_arh_factory_mo")
        !
        ! this subroutine is a mock routine for the C wrapper of the ARH factory for
        ! orbitals parameterized in the MO basis
        !
        use, intrinsic :: iso_c_binding, only: c_null_ptr
        use otr_arh, only: arh_n_micro
        use otr_arh_c_interface, only: arh_settings_type_c
        use c_interface, only: logger_c_type, solver_settings_type_c
        use test_reference, only: tol_c
        use otr_arh_test_reference, only: test_evaluate_dm_c_funptr, operator(/=)
        use otr_mo_test_reference, only: n_mo_c_ref => n_mo_c, n_occ_c_ref => n_occ_c, &
                                         orbsym_c_ref => orbsym_c, mo_coeff_pattern
        use otr_mo_c_interface_mock, only: mock_update_orbs_mo, orbsym_passed, &
                                           mo_coeff_3d_c
        use otr_common_test_reference, only: n_ao_c_ref => n_ao_c
        use otr_common_c_interface_mock, only: mock_level_obj_func, &
                                               mock_wire_solver_settings_c

        real(c_rp), intent(inout), target :: mo_coeff_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in) :: n_occ_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c, n_mo_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_arh_c_funptr, obj_func_arh_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(arh_settings_type_c), intent(in) :: settings_c
        integer(c_ip), intent(in), optional :: orbsym_c(*)
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! check passed dimensions
        if (n_ao_c /= n_ao_c_ref .or. n_mo_c /= n_mo_c_ref) then
            write(stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
                "number of AOs or MOs wrong."
            test_arh_factory_mo_interface = .false.
            error_c = 1
            return
        end if
        if (n_particle_c < 1 .or. n_particle_c > 2) then
            write(stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
                "number of particles wrong."
            test_arh_factory_mo_interface = .false.
            error_c = 1
            return
        end if
        if (any(n_occ_c(:n_particle_c) /= n_occ_c_ref(:n_particle_c))) then
            write(stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
                "number of occupied orbitals wrong."
            test_arh_factory_mo_interface = .false.
        end if

        ! set global pointer to MO coefficients so that they can be accessed in the
        ! mock orbital updating function
        call c_f_pointer(c_loc(mo_coeff_c(1)), mo_coeff_3d_c, &
                         [n_ao_c, n_mo_c, n_particle_c])

        ! check passed arrays, whose index pattern verifies column-major order
        if (any(abs(mo_coeff_3d_c - mo_coeff_pattern(n_ao_c, n_mo_c, n_particle_c, &
                                                     0.0_c_rp)) > tol_c)) then
            write(stderr, *) "test_arh_factory_mo_py_interface failed: Passed MO "// &
                "coefficients wrong."
            test_arh_factory_mo_interface = .false.
        end if
        if (any(abs(ao_overlap_c(:n_ao_c**2) - 2.0_c_rp) > tol_c)) then
            write(stderr, *) "test_arh_factory_mo_py_interface failed: Passed AO "// &
                "overlap matrix wrong."
            test_arh_factory_mo_interface = .false.
        end if

        ! test passed density matrix evaluating function
        if (.not. test_evaluate_dm_c_funptr( &
            evaluate_dm_c_funptr, int(n_particle_c, kind=ip), &
            "arh_factory_mo_py_interface", " by given density matrix evaluating "// &
            "function with Coulomb, exact-exchange and non-linear potential "// &
            "contributions", c_null_ptr)) test_arh_factory_mo_interface = .false.

        ! get Fortran pointer to passed logging function and call it
        message = "test"//c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message, c_null_ptr)

        ! check optional settings against reference values
        if (settings_c /= ref_arh_settings) then
            write(stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
                "settings associated with wrong values."
            test_arh_factory_mo_interface = .false.
        end if

        ! record whether irreps of the MOs are passed and check them against the
        ! reference values of the passed particle channels
        orbsym_passed = present(orbsym_c)
        if (present(orbsym_c)) then
            if (any(orbsym_c(:n_mo_c * n_particle_c) /= reshape( &
                orbsym_c_ref(:, :n_particle_c), [n_mo_c * n_particle_c]))) then
                write(stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
                    "irreps of the MOs wrong."
                test_arh_factory_mo_interface = .false.
            end if
        end if

        ! set function pointers to the mock functions of the level of the solver, and
        ! check the passed solver settings, wire them and adjust them as the factory
        ! does
        update_orbs_arh_c_funptr = c_funloc(mock_update_orbs_mo)
        obj_func_arh_c_funptr = c_funloc(mock_level_obj_func)
        if (.not. mock_wire_solver_settings_c(solver_settings_c, &
                                              "arh_factory_mo_py_interface")) &
            test_arh_factory_mo_interface = .false.
        solver_settings_c%refresh_hess = .true._c_bool
        solver_settings_c%n_micro = int(arh_n_micro, kind=c_ip)
        solver_settings_c%hess_symm = .false._c_bool

        ! set return arguments
        error_c = 0

    end function mock_arh_factory_mo_c_wrapper

    function mock_arh_factory_oao_c_wrapper( &
        dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
        update_orbs_arh_c_funptr, obj_func_arh_c_funptr, solver_settings_c, &
        settings_c) result(error_c) bind(C, name="mock_arh_factory_oao")
        !
        ! this subroutine is a mock routine for the ARH orbital updating factory C
        ! wrapper subroutine
        !
        use, intrinsic :: iso_c_binding, only: c_null_ptr
        use otr_arh, only: arh_n_micro
        use otr_arh_c_interface, only: arh_settings_type_c
        use c_interface, only: logger_c_type, solver_settings_type_c
        use test_reference, only: tol_c
        use otr_common_test_reference, only: n_ao_c_ref => n_ao_c
        use otr_arh_test_reference, only: test_evaluate_dm_c_funptr, operator(/=)
        use otr_common_c_interface_mock, only: &
            mock_level_update_orbs, mock_level_obj_func, mock_wire_solver_settings_c

        real(c_rp), intent(inout), target :: dm_ao_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_arh_c_funptr, obj_func_arh_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(arh_settings_type_c), intent(in) :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! closed-shell case
        if (n_particle_c == 1) then
            ! check passed arrays
            if (any(abs(dm_ao_c(:n_ao_c**2) - 1.0_c_rp) > tol_c)) then
                write(stderr, *) "test_arh_factory_oao_py_interface failed: Passed "// &
                    "AO density matrix for closed-shell case wrong."
                test_arh_factory_oao_interface = .false.
            end if
            if (any(abs(ao_overlap_c(:n_ao_c**2) - 2.0_c_rp) > tol_c)) then
                write(stderr, *) "test_arh_factory_oao_py_interface failed: Passed "// &
                    "AO overlap matrix wrong."
                test_arh_factory_oao_interface = .false.
            end if

            ! test passed density matrix evaluating function
            if (.not. test_evaluate_dm_c_funptr( &
                evaluate_dm_c_funptr, 1_ip, "arh_factory_oao_py_interface", " by "// &
                "given density matrix evaluating function with Coulomb, "// &
                "exact-exchange and non-linear potential contributions", c_null_ptr)) &
                test_arh_factory_oao_interface = .false.

            ! check if passed number of AOs is correct
            if (n_ao_c /= n_ao_c_ref) then
                write(stderr, *) "test_arh_factory_oao_py_interface failed: Passed "// &
                    "number of AOs wrong."
                test_arh_factory_oao_interface = .false.
            end if

            ! get Fortran pointer to passed logging function and call it
            message = "test"//c_null_char
            call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
            call logger_funptr(message, c_null_ptr)

            ! check optional settings against reference values
            if (settings_c /= ref_arh_settings) then
                write(stderr, *) "test_arh_factory_oao_py_interface failed: Passed "// &
                    "settings associated with wrong values."
                test_arh_factory_oao_interface = .false.
            end if

        ! open-shell case
        else if (n_particle_c == 2) then
            ! check passed arrays
            if (any(abs(dm_ao_c(:n_ao_c**2 * n_particle_c) - 1.0_c_rp) > tol_c)) then
                write(stderr, *) "test_arh_factory_oao_py_interface failed: Passed "// &
                    "AO density matrix for open-shell case wrong."
                test_arh_factory_oao_interface = .false.
            end if

            ! test passed density matrix evaluating function
            if (.not. test_evaluate_dm_c_funptr( &
                evaluate_dm_c_funptr, 2_ip, "arh_factory_oao_py_interface", " by "// &
                "given density matrix evaluating function with Coulomb, "// &
                "exact-exchange and non-linear potential contributions", c_null_ptr)) &
                test_arh_factory_oao_interface = .false.

        ! number of particles is not correct
        else
            write(stderr, *) "test_arh_factory_oao_py_interface failed: Passed "// &
                "number of particles wrong."
            test_arh_factory_oao_interface = .false.

        end if

        ! mark the density matrix of the caller, which the factory receives in place,
        ! set function pointers to the mock functions of the level of the solver, and
        ! check the passed solver settings, wire them and adjust them as the factory
        ! does
        dm_ao_c(:n_ao_c**2 * n_particle_c) = 2.0_c_rp
        update_orbs_arh_c_funptr = c_funloc(mock_level_update_orbs)
        obj_func_arh_c_funptr = c_funloc(mock_level_obj_func)
        if (.not. mock_wire_solver_settings_c(solver_settings_c, &
                                              "arh_factory_oao_py_interface")) &
            test_arh_factory_oao_interface = .false.
        solver_settings_c%refresh_hess = .true._c_bool
        solver_settings_c%hess_symm = .false._c_bool
        solver_settings_c%n_micro = int(arh_n_micro, kind=c_ip)

        ! set return arguments
        error_c = 0

    end function mock_arh_factory_oao_c_wrapper

    subroutine mock_init_arh_settings_c(settings) bind(C, name="mock_init_arh_settings")
        !
        ! this subroutine is a mock routine for the C ARH setting initialization
        ! subroutine
        !
        use otr_arh_c_interface, only: arh_settings_type_c
        use otr_arh_test_reference, only: assignment(=)

        type(arh_settings_type_c), intent(inout) :: settings

        ! set reference values
        settings = ref_arh_settings

    end subroutine mock_init_arh_settings_c

    subroutine mock_arh_deconstructor_c_wrapper() bind(C, name="mock_arh_deconstructor")
        !
        ! this subroutine is a mock routine for the C ARH deconstructor subroutine
        !
        test_arh_deconstructor_interface = .true._c_bool

    end subroutine mock_arh_deconstructor_c_wrapper

end module otr_arh_c_interface_mock
