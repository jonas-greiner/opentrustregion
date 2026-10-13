! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_c_interface_mock

    use opentrustregion, only: ip, stderr
    use c_interface, only: c_rp, c_ip
    use otr_oao_c_interface, only: oao_factory_c_wrapper, init_oao_settings_c, &
                                   oao_deconstructor_c_wrapper
    use otr_oao_test_reference, only: ref_oao_settings
    use otr_common_test_reference, only: n_particle, n_ao
    use test_reference, only: n_param
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char

    implicit none

    logical(c_bool), bind(C) :: test_oao_factory_interface = .true._c_bool, &
                                test_oao_deconstructor_interface = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(oao_factory_c_wrapper), pointer :: mock_oao_factory_c_wrapper_ptr => &
        mock_oao_factory_c_wrapper
    procedure(init_oao_settings_c), pointer :: mock_init_oao_settings_c_ptr => &
        mock_init_oao_settings_c
    procedure(oao_deconstructor_c_wrapper), pointer :: &
        mock_oao_deconstructor_c_wrapper_ptr => mock_oao_deconstructor_c_wrapper

contains

    function mock_oao_factory_c_wrapper( &
        dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
        update_orbs_oao_c_funptr, obj_func_oao_c_funptr, solver_settings_c, &
        settings_c) result(error_c) bind(C, name="mock_oao_factory")
        !
        ! this subroutine is a mock routine for the OAO orbital updating factory C
        ! wrapper subroutine
        !
        use, intrinsic :: iso_c_binding, only: c_null_ptr
        use otr_oao_c_interface, only: oao_settings_type_c
        use c_interface, only: logger_c_type, solver_settings_type_c
        use test_reference, only: tol_c
        use otr_oao_test_reference, only: operator(/=)
        use otr_common_test_reference, only: test_evaluate_dm_c_funptr
        use otr_common_c_interface_mock, only: &
            mock_level_update_orbs, mock_level_obj_func, mock_wire_solver_settings_c

        real(c_rp), intent(inout), target :: dm_ao_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_oao_c_funptr, obj_func_oao_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(oao_settings_type_c), intent(in) :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! closed-shell case
        if (n_particle_c == 1) then
            ! check passed arrays
            if (any(abs(dm_ao_c(:n_ao_c**2) - 1.0_c_rp) > tol_c)) then
                write(stderr, *) "test_oao_factory_py_interface failed: Passed AO "// &
                    "density matrix for closed-shell case wrong."
                test_oao_factory_interface = .false.
            end if
            if (any(abs(ao_overlap_c(:n_ao_c**2) - 2.0_c_rp) > tol_c)) then
                write(stderr, *) "test_oao_factory_py_interface failed: Passed AO "// &
                    "overlap matrix wrong."
                test_oao_factory_interface = .false.
            end if

            ! test passed density matrix evaluating function
            if (.not. test_evaluate_dm_c_funptr( &
                evaluate_dm_c_funptr, 1_ip, "oao_factory_py_interface", &
                " by given density matrix evaluating function", c_null_ptr)) &
                test_oao_factory_interface = .false.

            ! check if passed number of AOs is correct
            if (n_ao_c /= n_ao) then
                write(stderr, *) "test_oao_factory_py_interface failed: Passed "// &
                    "number of AOs wrong."
                test_oao_factory_interface = .false.
            end if

            ! get Fortran pointer to passed logging function and call it
            message = "test"//c_null_char
            call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
            call logger_funptr(message, c_null_ptr)

            ! check optional settings against reference values
            if (settings_c /= ref_oao_settings) then
                write(stderr, *) "test_oao_factory_py_interface failed: Passed "// &
                    "settings associated with wrong values."
                test_oao_factory_interface = .false.
            end if

        ! open-shell case
        else if (n_particle_c == 2) then
            ! check passed arrays
            if (any(abs(dm_ao_c(:n_ao_c**2 * n_particle_c) - 1.0_c_rp) > tol_c)) then
                write(stderr, *) "test_oao_factory_py_interface failed: Passed AO "// &
                    "density matrix for open-shell case wrong."
                test_oao_factory_interface = .false.
            end if

            ! test passed density matrix evaluating function
            if (.not. test_evaluate_dm_c_funptr( &
                evaluate_dm_c_funptr, 2_ip, "oao_factory_py_interface", &
                " by given density matrix evaluating function", c_null_ptr)) &
                test_oao_factory_interface = .false.

        ! number of particles is not correct
        else
            write(stderr, *) "test_oao_factory_py_interface failed: Passed number "// &
                "of particles wrong."
            test_oao_factory_interface = .false.

        end if

        ! mark the density matrix of the caller, which the factory receives in place,
        ! check the passed solver settings and wire them, and set function pointers to
        ! the mock functions of the level of the solver
        dm_ao_c(:n_ao_c**2 * n_particle_c) = 2.0_c_rp
        if (.not. mock_wire_solver_settings_c(solver_settings_c, &
                                              "oao_factory_py_interface")) &
            test_oao_factory_interface = .false.
        update_orbs_oao_c_funptr = c_funloc(mock_level_update_orbs)
        obj_func_oao_c_funptr = c_funloc(mock_level_obj_func)

        ! set return arguments
        error_c = 0

    end function mock_oao_factory_c_wrapper

    subroutine mock_init_oao_settings_c(settings) bind(C, name="mock_init_oao_settings")
        !
        ! this subroutine is a mock routine for the C OAO setting initialization
        ! subroutine
        !
        use otr_oao_c_interface, only: oao_settings_type_c
        use otr_oao_test_reference, only: assignment(=)

        type(oao_settings_type_c), intent(inout) :: settings

        ! set reference values
        settings = ref_oao_settings

    end subroutine mock_init_oao_settings_c

    subroutine mock_oao_deconstructor_c_wrapper() bind(C, name="mock_oao_deconstructor")
        !
        ! this subroutine is a mock routine for the C OAO deconstructor subroutine
        !
        test_oao_deconstructor_interface = .true._c_bool

    end subroutine mock_oao_deconstructor_c_wrapper

end module otr_oao_c_interface_mock
