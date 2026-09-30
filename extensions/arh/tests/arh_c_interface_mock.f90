! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip
    use otr_arh_test_reference, only: ref_arh_settings
    use otr_arh_c_interface, only: arh_factory_mo_c_wrapper, &
                                   arh_factory_oao_c_wrapper, init_arh_settings_c, &
                                   arh_deconstructor_c_wrapper
    use c_interface, only: update_orbs_c_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char, c_f_pointer, c_loc, &
                                           c_null_funptr

    implicit none

    logical(c_bool), bind(C) :: test_arh_factory_mo_interface = .true._c_bool, &
                                test_arh_factory_oao_interface = .true._c_bool, &
                                test_arh_deconstructor_interface = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_c_type), pointer :: mock_update_orbs_mo_c_ptr => &
        mock_update_orbs_mo_c
    procedure(arh_factory_mo_c_wrapper), pointer :: &
        mock_arh_factory_mo_c_wrapper_ptr => mock_arh_factory_mo_c_wrapper
    procedure(arh_factory_oao_c_wrapper), pointer :: &
        mock_arh_factory_oao_c_wrapper_ptr => mock_arh_factory_oao_c_wrapper
    procedure(init_arh_settings_c), pointer :: mock_init_arh_settings_c_ptr => &
        mock_init_arh_settings_c
    procedure(arh_deconstructor_c_wrapper), pointer :: &
        mock_arh_deconstructor_c_wrapper_ptr => mock_arh_deconstructor_c_wrapper

contains

    function mock_update_orbs_mo_c(kappa, func, grad, h_diag, hess_x_c_funptr) &
        result(error) bind(C)
        !
        ! this function is a test function for the orbital update C function for
        ! orbitals parameterized in the MO basis, which overwrites the MO coefficients
        ! with a pattern encoding their indices
        !
        use c_interface_unit_tests, only: mock_update_orbs_orig => mock_update_orbs
        use otr_arh_c_interface, only: mo_coeff_3d_c
        use otr_arh_test_reference, only: mo_coeff_pattern

        real(c_rp), intent(in), target :: kappa(*)
        real(c_rp), intent(out) :: func
        real(c_rp), intent(out), target :: grad(*), h_diag(*)
        type(c_funptr), intent(out) :: hess_x_c_funptr
        integer(c_ip) :: error

        error = mock_update_orbs_orig(kappa, func, grad, h_diag, hess_x_c_funptr)

        mo_coeff_3d_c = mo_coeff_pattern(size(mo_coeff_3d_c, 1, kind=c_ip), &
                                         size(mo_coeff_3d_c, 2, kind=c_ip), &
                                         size(mo_coeff_3d_c, 3, kind=c_ip), 1000.0_c_rp)

    end function mock_update_orbs_mo_c

    function mock_arh_factory_mo_c_wrapper( &
        mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
        evaluate_dm_c_funptr, obj_func_arh_c_funptr, update_orbs_arh_c_funptr, &
        solver_settings_c, settings_c) result(error_c) &
        bind(C, name="mock_arh_factory_mo")
        !
        ! this subroutine is a mock routine for the C wrapper of the ARH factory for
        ! orbitals parameterized in the MO basis
        !
        use opentrustregion, only: default_solver_settings
        use otr_arh, only: arh_n_micro
        use otr_arh_c_interface, only: arh_settings_type_c, mo_coeff_3d_c
        use c_interface, only: logger_c_type, solver_settings_type_c, assignment(=)
        use test_reference, only: tol_c
        use otr_arh_test_reference, only: &
            test_evaluate_dm_os_c_funptr, test_evaluate_dm_cs_c_funptr, &
            n_mo_c_ref => n_mo_c, n_occ_c_ref => n_occ_c, mo_coeff_pattern, operator(/=)
        use otr_common_test_reference, only: n_ao_c_ref => n_ao_c
        use c_interface_unit_tests, only: mock_obj_func, mock_precond, &
                                          mock_precond_pd, mock_get_extra_trial_vectors

        real(c_rp), intent(inout), target :: mo_coeff_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in) :: n_occ_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c, n_mo_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: obj_func_arh_c_funptr, update_orbs_arh_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(arh_settings_type_c), intent(inout) :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message

        ! check passed dimensions
        if (n_ao_c /= n_ao_c_ref .or. n_mo_c /= n_mo_c_ref) then
            write (stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
                "number of AOs or MOs wrong."
            test_arh_factory_mo_interface = .false.
            error_c = 1
            return
        end if
        if (n_particle_c < 1 .or. n_particle_c > 2) then
            write (stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
                "number of particles wrong."
            test_arh_factory_mo_interface = .false.
            error_c = 1
            return
        end if
        if (any(n_occ_c(:n_particle_c) /= n_occ_c_ref(:n_particle_c))) then
            write (stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
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
            write (stderr, *) "test_arh_factory_mo_py_interface failed: Passed MO "// &
                "coefficients wrong."
            test_arh_factory_mo_interface = .false.
        end if
        if (any(abs(ao_overlap_c(:n_ao_c**2) - 2.0_c_rp) > tol_c)) then
            write (stderr, *) "test_arh_factory_mo_py_interface failed: Passed AO "// &
                "overlap matrix wrong."
            test_arh_factory_mo_interface = .false.
        end if

        ! test passed density matrix evaluating function
        if (n_particle_c == 1) then
            test_arh_factory_mo_interface = &
                test_arh_factory_mo_interface .and. test_evaluate_dm_cs_c_funptr( &
                    evaluate_dm_c_funptr, "arh_factory_mo_py_interface", " by "// &
                    "given density matrix evaluating function with non-linear "// &
                    "potential contribution")
        else
            test_arh_factory_mo_interface = &
                test_arh_factory_mo_interface .and. test_evaluate_dm_os_c_funptr( &
                    evaluate_dm_c_funptr, "arh_factory_mo_py_interface", " by "// &
                    "given density matrix evaluating function with separate same- "// &
                    "and opposite-spin potential contributions")
        end if

        ! get Fortran pointer to passed logging function and call it
        message = "test" // c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message)

        ! check optional settings against reference values
        if (settings_c /= ref_arh_settings) then
            write (stderr, *) "test_arh_factory_mo_py_interface failed: Passed "// &
                "settings associated with wrong values."
            test_arh_factory_mo_interface = .false.
        end if

        ! set function pointers to mock to ARH mock functions, and set solver settings
        ! without wiring a projection, since the parameters are non-redundant, so that
        ! the projection is left to the caller
        obj_func_arh_c_funptr = c_funloc(mock_obj_func)
        update_orbs_arh_c_funptr = c_funloc(mock_update_orbs_mo_c)
        if (.not. solver_settings_c%initialized) &
            solver_settings_c = default_solver_settings
        solver_settings_c%precond = c_funloc(mock_precond)
        solver_settings_c%precond_pd = c_funloc(mock_precond_pd)
        solver_settings_c%get_extra_trial_vectors = &
            c_funloc(mock_get_extra_trial_vectors)
        solver_settings_c%stability_settings%precond = c_funloc(mock_precond)
        solver_settings_c%stability_settings%get_extra_trial_vectors = &
            c_funloc(mock_get_extra_trial_vectors)
        solver_settings_c%refresh_hess = .true._c_bool
        solver_settings_c%n_micro = int(arh_n_micro, kind=c_ip)
        solver_settings_c%hess_symm = .false._c_bool

        ! set return arguments
        error_c = 0

    end function mock_arh_factory_mo_c_wrapper

    function mock_arh_factory_oao_c_wrapper( &
        dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
        obj_func_arh_c_funptr, update_orbs_arh_c_funptr, solver_settings_c, &
        settings_c) result(error_c) bind(C, name="mock_arh_factory_oao")
        !
        ! this subroutine is a mock routine for the ARH orbital updating factory C
        ! wrapper subroutine
        !
        use opentrustregion, only: default_solver_settings
        use otr_arh, only: arh_n_micro
        use otr_arh_c_interface, only: arh_settings_type_c
        use c_interface, only: obj_func_c_type, update_orbs_c_type, precond_c_type, &
                               precond_pd_c_type, project_c_type, logger_c_type, &
                               solver_settings_type_c, assignment(=)
        use otr_oao_c_interface, only: dm_ao_3d_c
        use test_reference, only: tol_c
        use otr_common_test_reference, only: n_ao_c_ref => n_ao_c
        use otr_arh_test_reference, only: test_evaluate_dm_os_c_funptr, &
                                          test_evaluate_dm_cs_c_funptr, operator(/=)
        use c_interface_unit_tests, only: mock_obj_func, mock_precond, &
                                          mock_precond_pd, mock_project, &
                                          mock_get_extra_trial_vectors
        use otr_oao_c_interface_mock, only: mock_update_orbs_oao

        real(c_rp), intent(in), target :: dm_ao_c(*), ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: obj_func_arh_c_funptr, update_orbs_arh_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(arh_settings_type_c), intent(inout) :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message
        procedure(obj_func_c_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_c_type), pointer :: update_orbs_arh_funptr
        procedure(precond_c_type), pointer :: precond_arh_funptr
        procedure(precond_pd_c_type), pointer :: precond_pd_arh_funptr
        procedure(project_c_type), pointer :: project_arh_funptr

        ! set global pointer to density matrix so that it can be accessed in the mock
        ! orbital updating function
        call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_3d_c, [n_ao_c, n_ao_c, n_particle_c])

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
            test_arh_factory_oao_interface = &
                test_arh_factory_oao_interface .and. test_evaluate_dm_cs_c_funptr( &
                    evaluate_dm_c_funptr, "arh_factory_oao_py_interface", " by "// &
                    "given density matrix evaluating function with non-linear "// &
                    "potential contribution")

            ! check if passed number of AOs is correct
            if (n_ao_c /= n_ao_c_ref) then
                write (stderr, *) "test_arh_factory_oao_py_interface failed: "// &
                    "Passed number of AOs wrong."
                test_arh_factory_oao_interface = .false.
            end if

            ! get Fortran pointer to passed logging function and call it
            message = "test" // c_null_char
            call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
            call logger_funptr(message)

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
            test_arh_factory_oao_interface = &
                test_arh_factory_oao_interface .and. test_evaluate_dm_os_c_funptr( &
                    evaluate_dm_c_funptr, "arh_factory_oao_py_interface", " by "// &
                    "given density matrix evaluating function with separate same- "// &
                    "and opposite-spin potential contributions")

        ! number of particles is not correct
        else
            write (stderr, *) "test_arh_factory_oao_py_interface failed: Passed "// &
                "number of particles wrong."
            test_arh_factory_oao_interface = .false.

        end if

        ! set function pointers to mock to ARH mock functions, and set solver settings
        obj_func_arh_c_funptr = c_funloc(mock_obj_func)
        update_orbs_arh_c_funptr = c_funloc(mock_update_orbs_oao)
        if (.not. solver_settings_c%initialized) &
            solver_settings_c = default_solver_settings
        solver_settings_c%precond = c_funloc(mock_precond)
        solver_settings_c%precond_pd = c_funloc(mock_precond_pd)
        solver_settings_c%project = c_funloc(mock_project)
        solver_settings_c%get_extra_trial_vectors = &
            c_funloc(mock_get_extra_trial_vectors)
        solver_settings_c%stability_settings%precond = c_funloc(mock_precond)
        solver_settings_c%stability_settings%project = c_funloc(mock_project)
        solver_settings_c%stability_settings%get_extra_trial_vectors = &
            c_funloc(mock_get_extra_trial_vectors)
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
