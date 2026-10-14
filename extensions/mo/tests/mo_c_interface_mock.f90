! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_c_interface_mock

    use opentrustregion, only: ip, stderr
    use c_interface, only: c_rp, c_ip, update_orbs_c_type
    use otr_mo_c_interface, only: mo_factory_c_wrapper, init_mo_settings_c
    use otr_mo_test_reference, only: ref_mo_settings
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char, c_f_pointer, c_loc, &
                                           c_null_funptr, c_null_ptr

    implicit none

    logical(c_bool), bind(C) :: test_mo_factory_interface = .true._c_bool

    ! whether the mock factory received settings and a density matrix evaluating
    ! function, which the Python test reads
    logical(c_bool), bind(C) :: mo_settings_passed = .false._c_bool, &
                                mo_evaluate_dm_passed = .false._c_bool

    ! whether the last call of a mock factory for orbitals parameterized in the MO
    ! basis was passed irreps of the MOs
    logical(c_bool), bind(C, name="test_orbsym_passed") :: orbsym_passed = &
        .false._c_bool

    ! MO coefficients passed to the mock factory, which the mock orbital update
    ! overwrites to test that the Python interface copies them back
    real(c_rp), pointer :: mo_coeff_3d_c(:, :, :) => null()

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_c_type), pointer :: mock_update_orbs_mo_ptr => &
        mock_update_orbs_mo
    procedure(mo_factory_c_wrapper), pointer :: mock_mo_factory_c_wrapper_ptr => &
        mock_mo_factory_c_wrapper
    procedure(init_mo_settings_c), pointer :: mock_init_mo_settings_c_ptr => &
        mock_init_mo_settings_c

contains

    function mock_update_orbs_mo(kappa, func, grad, h_diag, hess_x_c_funptr, &
                                 context_c) result(error) bind(C)
        !
        ! this function is a test function for the orbital update C function for
        ! orbitals parameterized in the MO basis, which overwrites the MO coefficients
        ! with a pattern encoding their indices in addition to the mock orbital update
        ! of the level of the solver
        !
        use otr_common_c_interface_mock, only: mock_level_update_orbs
        use otr_mo_test_reference, only: mo_coeff_pattern

        real(c_rp), intent(in), target :: kappa(*)
        real(c_rp), intent(out) :: func
        real(c_rp), intent(out), target :: grad(*), h_diag(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_level_update_orbs(kappa, func, grad, h_diag, hess_x_c_funptr, &
                                       context_c)

        mo_coeff_3d_c = mo_coeff_pattern(size(mo_coeff_3d_c, 1, kind=c_ip), &
                                         size(mo_coeff_3d_c, 2, kind=c_ip), &
                                         size(mo_coeff_3d_c, 3, kind=c_ip), 1000.0_c_rp)

    end function mock_update_orbs_mo

    function mock_mo_factory_c_wrapper( &
        mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
        evaluate_dm_c_funptr, update_orbs_mo_c_funptr, obj_func_mo_c_funptr, &
        extension_c, solver_settings_c, solver_state_c, settings_c, orbsym_c) &
        result(error_c) bind(C, name="mock_mo_factory")
        !
        ! this subroutine is a mock routine for the MO orbital updating factory C
        ! wrapper subroutine
        !
        use, intrinsic :: iso_c_binding, only: c_null_ptr, c_associated
        use otr_mo_c_interface, only: mo_settings_type_c
        use c_interface, only: logger_c_type, solver_settings_type_c
        use test_reference, only: tol_c
        use otr_mo_test_reference, only: n_mo_c_ref => n_mo_c, n_occ_c_ref => n_occ_c, &
                                         orbsym_c_ref => orbsym_c, mo_coeff_pattern, &
                                         operator(/=)
        use otr_common_test_reference, only: &
            n_ao_c_ref => n_ao_c, check_evaluate_dm_c_funptr, ao_overlap_pattern
        use otr_common_c_interface_mock, only: &
            mock_level_obj_func, check_and_wire_solver_settings_c, mock_extension

        real(c_rp), intent(inout), target :: mo_coeff_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in) :: n_occ_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c, n_mo_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_mo_c_funptr, obj_func_mo_c_funptr
        type(c_ptr), intent(out) :: extension_c
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_ptr), intent(in), value :: solver_state_c
        type(mo_settings_type_c), intent(in), optional :: settings_c
        integer(c_ip), intent(in), optional :: orbsym_c(*)
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! check passed dimensions
        extension_c = c_null_ptr
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
        if (any(abs(ao_overlap_c(:n_ao_c**2) - &
                    reshape(ao_overlap_pattern(n_ao_c), [n_ao_c**2])) > tol_c)) then
            write(stderr, *) "test_mo_factory_py_interface failed: Passed AO "// &
                "overlap matrix wrong."
            test_mo_factory_interface = .false.
        end if

        ! record whether a density matrix evaluating function was passed and, if so,
        ! test it for the number of particles
        mo_evaluate_dm_passed = c_associated(evaluate_dm_c_funptr)
        if (mo_evaluate_dm_passed) then
            if (.not. check_evaluate_dm_c_funptr( &
                evaluate_dm_c_funptr, int(n_particle_c, kind=ip), &
                "mo_factory_py_interface", " by given density matrix evaluating "// &
                "function", c_null_ptr)) test_mo_factory_interface = .false.
        end if

        ! record whether settings were passed and, if so, call their logging function
        ! and check them against the reference values
        mo_settings_passed = present(settings_c)
        if (present(settings_c)) then
            ! get Fortran pointer to passed logging function and call it
            message = "test"//c_null_char
            call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
            call logger_funptr(message, c_null_ptr)

            ! check optional settings against reference values
            if (settings_c /= ref_mo_settings) then
                write(stderr, *) "test_mo_factory_py_interface failed: Passed "// &
                    "settings associated with wrong values."
                test_mo_factory_interface = .false.
            end if
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

        ! check the passed solver settings and wire them, set function pointers to the
        ! mock functions of the level of the solver and return an extension, recording
        ! the record of the state
        if (.not. check_and_wire_solver_settings_c(solver_settings_c, &
                                                   "mo_factory_py_interface")) &
            test_mo_factory_interface = .false.
        update_orbs_mo_c_funptr = c_funloc(mock_update_orbs_mo)
        obj_func_mo_c_funptr = c_funloc(mock_level_obj_func)
        extension_c = mock_extension(solver_state_c)

        ! set return arguments
        error_c = 0

    end function mock_mo_factory_c_wrapper

    subroutine mock_init_mo_settings_c(settings) bind(C, name="mock_init_mo_settings")
        !
        ! this subroutine is a mock routine for the C MO setting initialization
        ! subroutine
        !
        use otr_mo_c_interface, only: mo_settings_type_c
        use otr_mo_test_reference, only: get_reference_mo_values

        type(mo_settings_type_c), intent(inout) :: settings

        ! set reference values without logging function and host context
        call get_reference_mo_values(settings)
        settings%logger = c_null_funptr
        settings%context = c_null_ptr

    end subroutine mock_init_mo_settings_c

end module otr_mo_c_interface_mock
