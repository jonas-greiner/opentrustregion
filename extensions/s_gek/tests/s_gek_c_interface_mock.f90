! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_ip
    use otr_s_gek_c_interface, only: s_gek_factory_c_wrapper, init_s_gek_settings_c
    use otr_s_gek_test_reference, only: ref_s_gek_settings
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char, c_null_funptr, &
                                           c_null_ptr

    implicit none

    logical(c_bool), bind(C) :: test_s_gek_factory_interface = .true._c_bool

    ! whether the mock factory received settings, which the Python test reads
    logical(c_bool), bind(C) :: s_gek_settings_passed = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(s_gek_factory_c_wrapper), pointer :: mock_s_gek_factory_c_wrapper_ptr => &
        mock_s_gek_factory_c_wrapper
    procedure(init_s_gek_settings_c), pointer :: mock_init_s_gek_settings_c_ptr => &
        mock_init_s_gek_settings_c

contains

    function mock_s_gek_factory_c_wrapper( &
        n_param_c, update_orbs_c_funptr, obj_func_c_funptr, change_reference_c_funptr, &
        update_orbs_s_gek_c_funptr, obj_func_s_gek_c_funptr, extension_c, &
        solver_settings_c, solver_state_c, settings_c) result(error_c) &
        bind(C, name="mock_s_gek_factory")
        !
        ! this subroutine is a mock routine for the S-GEK orbital updating factory C
        ! wrapper subroutine
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c
        use c_interface, only: logger_c_type, solver_settings_type_c
        use test_reference, only: check_update_orbs_c_funptr, check_obj_func_c_funptr, &
                                  n_param
        use otr_s_gek_test_reference, only: check_change_reference_c_funptr, &
                                            operator(/=)
        use otr_common_c_interface_mock, only: &
            mock_level_update_orbs, mock_level_obj_func, &
            check_and_wire_solver_settings_c, mock_extension

        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             change_reference_c_funptr
        type(c_funptr), intent(out) :: update_orbs_s_gek_c_funptr, &
                                       obj_func_s_gek_c_funptr
        type(c_ptr), intent(out) :: extension_c
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_ptr), intent(in), value :: solver_state_c
        type(s_gek_settings_type_c), intent(in), optional :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! check if passed number of parameters is correct
        if (n_param_c /= n_param) then
            write(stderr, *) "test_s_gek_factory_py_interface failed: Passed "// &
                "number of parameters wrong."
            test_s_gek_factory_interface = .false.
        end if

        ! test passed orbital update function
        if (.not. check_update_orbs_c_funptr( &
            update_orbs_c_funptr, "s_gek_factory_py_interface", " by given orbital "// &
            "updating function", c_null_ptr)) test_s_gek_factory_interface = .false.

        ! test passed objective function
        if (.not. check_obj_func_c_funptr( &
            obj_func_c_funptr, "s_gek_factory_py_interface", " by given objective "// &
            "function", c_null_ptr)) test_s_gek_factory_interface = .false.

        ! test passed change reference function
        if (.not. check_change_reference_c_funptr( &
            change_reference_c_funptr, "s_gek_factory_py_interface", &
            " by given change reference function", c_null_ptr)) &
            test_s_gek_factory_interface = .false.

        ! record whether settings were passed and, if so, call their logging function
        ! and check them against the reference values
        s_gek_settings_passed = present(settings_c)
        if (present(settings_c)) then
            ! get Fortran pointer to passed logging function and call it
            message = "test"//c_null_char
            call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
            call logger_funptr(message, c_null_ptr)

            ! check optional settings against reference values
            if (settings_c /= ref_s_gek_settings) then
                write(stderr, *) "test_s_gek_factory_py_interface failed: Passed "// &
                    "settings associated with wrong values."
                test_s_gek_factory_interface = .false.
            end if
        end if

        ! check the passed solver settings and wire them
        if (.not. check_and_wire_solver_settings_c(solver_settings_c, &
                                                   "s_gek_factory_py_interface")) &
            test_s_gek_factory_interface = .false.

        ! set function pointers to the mock functions of the level of the solver and
        ! return an extension, recording the record of the state
        update_orbs_s_gek_c_funptr = c_funloc(mock_level_update_orbs)
        obj_func_s_gek_c_funptr = c_funloc(mock_level_obj_func)
        extension_c = mock_extension(solver_state_c)

        ! set return arguments
        error_c = 0

    end function mock_s_gek_factory_c_wrapper

    subroutine mock_init_s_gek_settings_c(settings) &
        bind(C, name="mock_init_s_gek_settings")
        !
        ! this subroutine is a mock routine for the C S-GEK setting initialization
        ! subroutine
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c
        use otr_s_gek_test_reference, only: get_reference_s_gek_values

        type(s_gek_settings_type_c), intent(inout) :: settings

        ! set reference values without logging function and host context
        call get_reference_s_gek_values(settings)
        settings%logger = c_null_funptr
        settings%context = c_null_ptr

    end subroutine mock_init_s_gek_settings_c

end module otr_s_gek_c_interface_mock
