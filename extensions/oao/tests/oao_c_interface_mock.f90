! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_c_interface_mock

    use opentrustregion, only: ip, stderr
    use c_interface, only: c_rp, c_ip
    use otr_oao_c_interface, only: oao_factory_c_wrapper, init_oao_settings_c
    use otr_oao_test_reference, only: ref_oao_settings
    use otr_common_test_reference, only: n_particle, n_ao
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char, c_null_funptr, &
                                           c_null_ptr

    implicit none

    logical(c_bool), bind(C) :: test_oao_factory_interface = .true._c_bool

    ! whether the mock factory received settings and a density matrix evaluating
    ! function, which the Python test reads
    logical(c_bool), bind(C) :: oao_settings_passed = .false._c_bool, &
                                oao_evaluate_dm_passed = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(oao_factory_c_wrapper), pointer :: mock_oao_factory_c_wrapper_ptr => &
        mock_oao_factory_c_wrapper
    procedure(init_oao_settings_c), pointer :: mock_init_oao_settings_c_ptr => &
        mock_init_oao_settings_c

contains

    function mock_oao_factory_c_wrapper( &
        dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
        update_orbs_oao_c_funptr, obj_func_oao_c_funptr, extension_c, &
        solver_settings_c, solver_state_c, settings_c) result(error_c) &
        bind(C, name="mock_oao_factory")
        !
        ! this subroutine is a mock routine for the OAO orbital updating factory C
        ! wrapper subroutine
        !
        use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated
        use otr_oao_c_interface, only: oao_settings_type_c
        use c_interface, only: logger_c_type, solver_settings_type_c
        use test_reference, only: tol_c
        use otr_oao_test_reference, only: operator(/=)
        use otr_common_test_reference, only: check_evaluate_dm_c_funptr
        use otr_common_c_interface_mock, only: &
            mock_level_update_orbs, mock_level_obj_func, &
            check_and_wire_solver_settings_c, mock_extension

        real(c_rp), intent(inout), target :: dm_ao_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_oao_c_funptr, obj_func_oao_c_funptr
        type(c_ptr), intent(out) :: extension_c
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_ptr), intent(in), value :: solver_state_c
        type(oao_settings_type_c), intent(in), optional :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! check passed arrays and number of AOs
        if (n_particle_c < 1 .or. n_particle_c > n_particle) then
            write(stderr, *) "test_oao_factory_py_interface failed: Passed number "// &
                "of particles wrong."
            test_oao_factory_interface = .false.
        else if (any(abs(dm_ao_c(:n_ao_c**2 * n_particle_c) - 1.0_c_rp) > tol_c)) then
            write(stderr, *) "test_oao_factory_py_interface failed: Passed AO "// &
                "density matrix wrong."
            test_oao_factory_interface = .false.
        end if
        if (any(abs(ao_overlap_c(:n_ao_c**2) - 2.0_c_rp) > tol_c)) then
            write(stderr, *) "test_oao_factory_py_interface failed: Passed AO "// &
                "overlap matrix wrong."
            test_oao_factory_interface = .false.
        end if
        if (n_ao_c /= n_ao) then
            write(stderr, *) "test_oao_factory_py_interface failed: Passed number "// &
                "of AOs wrong."
            test_oao_factory_interface = .false.
        end if

        ! record whether a density matrix evaluating function was passed and, if so,
        ! test it for the number of particles
        oao_evaluate_dm_passed = c_associated(evaluate_dm_c_funptr)
        if (oao_evaluate_dm_passed) then
            if (.not. check_evaluate_dm_c_funptr( &
                evaluate_dm_c_funptr, int(n_particle_c, kind=ip), &
                "oao_factory_py_interface", " by given density matrix evaluating "// &
                "function", c_null_ptr)) test_oao_factory_interface = .false.
        end if

        ! record whether settings were passed and, if so, call their logging function
        ! and check them against the reference values
        oao_settings_passed = present(settings_c)
        if (present(settings_c)) then
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
        end if

        ! mark the density matrix of the caller, which the factory receives in place,
        ! and check the passed solver settings and wire them
        dm_ao_c(:n_ao_c**2 * n_particle_c) = 2.0_c_rp
        if (.not. check_and_wire_solver_settings_c(solver_settings_c, &
                                                   "oao_factory_py_interface")) &
            test_oao_factory_interface = .false.

        ! set function pointers to the mock functions of the level of the solver and
        ! return an extension, recording the record of the state
        update_orbs_oao_c_funptr = c_funloc(mock_level_update_orbs)
        obj_func_oao_c_funptr = c_funloc(mock_level_obj_func)
        extension_c = mock_extension(solver_state_c)

        ! set return arguments
        error_c = 0

    end function mock_oao_factory_c_wrapper

    subroutine mock_init_oao_settings_c(settings) bind(C, name="mock_init_oao_settings")
        !
        ! this subroutine is a mock routine for the C OAO setting initialization
        ! subroutine
        !
        use otr_oao_c_interface, only: oao_settings_type_c
        use otr_oao_test_reference, only: get_reference_oao_values

        type(oao_settings_type_c), intent(inout) :: settings

        ! set reference values without logging function and host context
        call get_reference_oao_values(settings)
        settings%logger = c_null_funptr
        settings%context = c_null_ptr

    end subroutine mock_init_oao_settings_c

end module otr_oao_c_interface_mock
