! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_mock

    use opentrustregion, only: ip, stderr
    use c_interface, only: c_ip
    use otr_arh_test_reference, only: ref_arh_settings
    use otr_arh_c_interface, only: arh_factory_c_wrapper, init_arh_settings_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char

    implicit none

    logical(c_bool), bind(C) :: test_arh_factory_interface = .true._c_bool

    ! whether the mock factory received settings, which the Python test reads, and the
    ! number of particles of the orbitals of the solver settings the Python test
    ! passes, for which the mock factory tests the density matrix evaluating function
    logical(c_bool), bind(C) :: arh_settings_passed = .false._c_bool
    integer(c_ip), bind(C) :: arh_n_particle = 1

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory_c_wrapper), pointer :: mock_arh_factory_c_wrapper_ptr => &
        mock_arh_factory_c_wrapper
    procedure(init_arh_settings_c), pointer :: mock_init_arh_settings_c_ptr => &
        mock_init_arh_settings_c

contains

    function mock_arh_factory_c_wrapper( &
        evaluate_dm_c_funptr, update_orbs_arh_c_funptr, obj_func_arh_c_funptr, &
        extension_c, solver_settings_c, solver_state_c, settings_c) result(error_c) &
        bind(C, name="mock_arh_factory")
        !
        ! this subroutine is a mock routine for the C wrapper of the ARH factory
        !
        use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
        use otr_arh, only: arh_n_micro
        use otr_arh_c_interface, only: arh_settings_type_c
        use c_interface, only: logger_c_type, solver_settings_type_c
        use otr_arh_test_reference, only: check_arh_evaluate_dm_c_funptr, operator(/=)
        use otr_common_c_interface_mock, only: &
            mock_level_update_orbs, mock_level_obj_func, mock_wire_solver_settings_c, &
            mock_extension

        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_arh_c_funptr, obj_func_arh_c_funptr
        type(c_ptr), intent(out) :: extension_c
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_ptr), intent(in), value :: solver_state_c
        type(arh_settings_type_c), intent(in), optional :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! test passed density matrix evaluating function for the number of particles of
        ! the orbitals
        if (.not. check_arh_evaluate_dm_c_funptr( &
            evaluate_dm_c_funptr, int(arh_n_particle, kind=ip), &
            "arh_factory_py_interface", " by given density matrix evaluating "// &
            "function with Coulomb, exact-exchange and non-linear potential "// &
            "contributions", c_null_ptr)) test_arh_factory_interface = .false.

        ! record whether settings were passed and, if so, call their logging function
        ! and check them against the reference values
        arh_settings_passed = present(settings_c)
        if (present(settings_c)) then
            ! get Fortran pointer to passed logging function and call it
            message = "test"//c_null_char
            call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
            call logger_funptr(message, c_null_ptr)

            ! check optional settings against reference values
            if (settings_c /= ref_arh_settings) then
                write(stderr, *) "test_arh_factory_py_interface failed: Passed "// &
                    "settings associated with wrong values."
                test_arh_factory_interface = .false.
            end if
        end if

        ! check the passed solver settings, wire them and adjust them as the factory
        ! does
        if (.not. mock_wire_solver_settings_c(solver_settings_c, &
                                              "arh_factory_py_interface")) &
            test_arh_factory_interface = .false.
        solver_settings_c%refresh_hess = .true._c_bool
        solver_settings_c%n_micro = int(arh_n_micro, kind=c_ip)
        solver_settings_c%hess_symm = .false._c_bool

        ! set function pointers to the mock functions of the level of the solver and
        ! return an extension, recording the record of the state
        update_orbs_arh_c_funptr = c_funloc(mock_level_update_orbs)
        obj_func_arh_c_funptr = c_funloc(mock_level_obj_func)
        extension_c = mock_extension(solver_state_c)

        ! set return arguments
        error_c = 0

    end function mock_arh_factory_c_wrapper

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

end module otr_arh_c_interface_mock
