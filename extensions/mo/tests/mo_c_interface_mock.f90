! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip, update_orbs_c_type
    use otr_mo_c_interface, only: mo_factory_c_wrapper, init_mo_settings_c, &
                                  mo_deconstructor_c_wrapper
    use otr_mo_test_reference, only: ref_mo_settings
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char, c_f_pointer, c_loc

    implicit none

    logical(c_bool), bind(C) :: test_mo_factory_interface = .true._c_bool, &
                                test_mo_deconstructor_interface = .false._c_bool

    ! whether the last call of a mock factory for orbitals parameterized in the MO
    ! basis was passed irreps of the MOs
    logical(c_bool), bind(C, name="test_orbsym_passed") :: orbsym_passed = &
        .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_c_type), pointer :: mock_update_orbs_mo_ptr => &
        mock_update_orbs_mo
    procedure(mo_factory_c_wrapper), pointer :: mock_mo_factory_c_wrapper_ptr => &
        mock_mo_factory_c_wrapper
    procedure(init_mo_settings_c), pointer :: mock_init_mo_settings_c_ptr => &
        mock_init_mo_settings_c
    procedure(mo_deconstructor_c_wrapper), pointer :: &
        mock_mo_deconstructor_c_wrapper_ptr => mock_mo_deconstructor_c_wrapper

contains

    function mock_update_orbs_mo(kappa, func, grad, h_diag, hess_x_c_funptr) &
        result(error) bind(C)
        !
        ! this function is a test function for the orbital update C function for
        ! orbitals parameterized in the MO basis, which overwrites the MO coefficients
        ! with a pattern encoding their indices
        !
        use c_interface_unit_tests, only: mock_update_orbs_orig => mock_update_orbs
        use otr_mo_c_interface, only: mo_coeff_3d_c
        use otr_mo_test_reference, only: mo_coeff_pattern

        real(c_rp), intent(in), target :: kappa(*)
        real(c_rp), intent(out) :: func
        real(c_rp), intent(out), target :: grad(*), h_diag(*)
        type(c_funptr), intent(out) :: hess_x_c_funptr
        integer(c_ip) :: error

        error = mock_update_orbs_orig(kappa, func, grad, h_diag, hess_x_c_funptr)

        mo_coeff_3d_c = mo_coeff_pattern(size(mo_coeff_3d_c, 1, kind=c_ip), &
                                         size(mo_coeff_3d_c, 2, kind=c_ip), &
                                         size(mo_coeff_3d_c, 3, kind=c_ip), 1000.0_c_rp)

    end function mock_update_orbs_mo

    function mock_mo_factory_c_wrapper( &
        mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
        evaluate_dm_c_funptr, obj_func_mo_c_funptr, update_orbs_mo_c_funptr, &
        solver_settings_c, settings_c, orbsym_c) result(error_c) &
        bind(C, name="mock_mo_factory")
        !
        ! this subroutine is a mock routine for the MO orbital updating factory C
        ! wrapper subroutine
        !
        use opentrustregion, only: default_solver_settings
        use otr_mo_c_interface, only: mo_settings_type_c, mo_coeff_3d_c
        use c_interface, only: logger_c_type, solver_settings_type_c, assignment(=)
        use test_reference, only: tol_c
        use otr_mo_test_reference, only: n_mo_c_ref => n_mo_c, n_occ_c_ref => n_occ_c, &
                                         orbsym_c_ref => orbsym_c, mo_coeff_pattern, &
                                         operator(/=)
        use otr_common_test_reference, only: n_ao_c_ref => n_ao_c, &
                                             test_evaluate_dm_cs_c_funptr, &
                                             test_evaluate_dm_os_c_funptr
        use c_interface_unit_tests, only: mock_obj_func, mock_precond, &
                                          mock_precond_pd, mock_get_extra_trial_vectors

        real(c_rp), intent(inout), target :: mo_coeff_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in) :: n_occ_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c, n_mo_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: obj_func_mo_c_funptr, update_orbs_mo_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(mo_settings_type_c), intent(inout) :: settings_c
        integer(c_ip), intent(in), optional :: orbsym_c(*)
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! check passed dimensions
        if (n_ao_c /= n_ao_c_ref .or. n_mo_c /= n_mo_c_ref) then
            write(stderr, *) "test_mo_factory_py_interface failed: Passed number "// &
                "of AOs or MOs wrong."
            test_mo_factory_interface = .false.
            error_c = 1
            return
        end if
        if (n_particle_c < 1 .or. n_particle_c > 2) then
            write(stderr, *) "test_mo_factory_py_interface failed: Passed number "// &
                "of particles wrong."
            test_mo_factory_interface = .false.
            error_c = 1
            return
        end if
        if (any(n_occ_c(:n_particle_c) /= n_occ_c_ref(:n_particle_c))) then
            write(stderr, *) "test_mo_factory_py_interface failed: Passed number "// &
                "of occupied orbitals wrong."
            test_mo_factory_interface = .false.
        end if

        ! set global pointer to MO coefficients so that they can be accessed in the
        ! mock orbital updating function
        call c_f_pointer(c_loc(mo_coeff_c(1)), mo_coeff_3d_c, &
                         [n_ao_c, n_mo_c, n_particle_c])

        ! check passed arrays, whose index pattern verifies column-major order
        if (any(abs(mo_coeff_3d_c - mo_coeff_pattern(n_ao_c, n_mo_c, n_particle_c, &
                                                     0.0_c_rp)) > tol_c)) then
            write(stderr, *) "test_mo_factory_py_interface failed: Passed MO "// &
                "coefficients wrong."
            test_mo_factory_interface = .false.
        end if
        if (any(abs(ao_overlap_c(:n_ao_c**2) - 2.0_c_rp) > tol_c)) then
            write(stderr, *) "test_mo_factory_py_interface failed: Passed AO "// &
                "overlap matrix wrong."
            test_mo_factory_interface = .false.
        end if

        ! test passed density matrix evaluating function
        if (n_particle_c == 1) then
            test_mo_factory_interface = &
                test_mo_factory_interface .and. test_evaluate_dm_cs_c_funptr( &
                    evaluate_dm_c_funptr, "mo_factory_py_interface", &
                    " by given density matrix evaluating function")
        else
            test_mo_factory_interface = &
                test_mo_factory_interface .and. test_evaluate_dm_os_c_funptr( &
                    evaluate_dm_c_funptr, "mo_factory_py_interface", &
                    " by given density matrix evaluating function")
        end if

        ! get Fortran pointer to passed logging function and call it
        message = "test"//c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message)

        ! check optional settings against reference values
        if (settings_c /= ref_mo_settings) then
            write(stderr, *) "test_mo_factory_py_interface failed: Passed settings "// &
                "associated with wrong values."
            test_mo_factory_interface = .false.
        end if

        ! record whether irreps of the MOs are passed and check them against the
        ! reference values of the passed particle channels
        orbsym_passed = present(orbsym_c)
        if (present(orbsym_c)) then
            if (any(orbsym_c(:n_mo_c * n_particle_c) /= reshape( &
                orbsym_c_ref(:, :n_particle_c), [n_mo_c * n_particle_c]))) then
                write(stderr, *) "test_mo_factory_py_interface failed: Passed "// &
                    "irreps of the MOs wrong."
                test_mo_factory_interface = .false.
            end if
        end if

        ! set function pointers to mock to MO mock functions, and set solver settings
        ! without wiring a projection, since the parameters are non-redundant, so that
        ! the projection is left to the caller
        obj_func_mo_c_funptr = c_funloc(mock_obj_func)
        update_orbs_mo_c_funptr = c_funloc(mock_update_orbs_mo)
        if (.not. solver_settings_c%initialized) &
            solver_settings_c = default_solver_settings
        solver_settings_c%precond = c_funloc(mock_precond)
        solver_settings_c%precond_pd = c_funloc(mock_precond_pd)
        solver_settings_c%get_extra_trial_vectors = &
            c_funloc(mock_get_extra_trial_vectors)
        solver_settings_c%stability_settings%precond = c_funloc(mock_precond)
        solver_settings_c%stability_settings%get_extra_trial_vectors = &
            c_funloc(mock_get_extra_trial_vectors)

        ! set return arguments
        error_c = 0

    end function mock_mo_factory_c_wrapper

    subroutine mock_init_mo_settings_c(settings) bind(C, name="mock_init_mo_settings")
        !
        ! this subroutine is a mock routine for the C MO setting initialization
        ! subroutine
        !
        use otr_mo_c_interface, only: mo_settings_type_c
        use otr_mo_test_reference, only: assignment(=)

        type(mo_settings_type_c), intent(inout) :: settings

        ! set reference values
        settings = ref_mo_settings

    end subroutine mock_init_mo_settings_c

    subroutine mock_mo_deconstructor_c_wrapper() bind(C, name="mock_mo_deconstructor")
        !
        ! this subroutine is a mock routine for the C MO deconstructor subroutine
        !
        test_mo_deconstructor_interface = .true._c_bool

    end subroutine mock_mo_deconstructor_c_wrapper

end module otr_mo_c_interface_mock
