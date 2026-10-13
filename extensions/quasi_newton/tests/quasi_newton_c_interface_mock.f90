! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_ip
    use otr_qn_c_interface, only: qn_factory_c_wrapper, init_qn_settings_c, &
                                  update_orbs_qn_deconstructor_c_wrapper
    use otr_qn_test_reference, only: ref_qn_settings
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char

    implicit none

    logical(c_bool), bind(C) :: &
        test_qn_factory_interface = .true._c_bool, &
        test_update_orbs_qn_deconstructor_interface = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(qn_factory_c_wrapper), pointer :: mock_qn_factory_c_wrapper_ptr => &
        mock_qn_factory_c_wrapper
    procedure(init_qn_settings_c), pointer :: mock_init_qn_settings_c_ptr => &
        mock_init_qn_settings_c
    procedure(update_orbs_qn_deconstructor_c_wrapper), pointer :: &
        mock_update_orbs_qn_deconstructor_c_wrapper_ptr => &
        mock_update_orbs_qn_deconstructor_c_wrapper

contains

    function mock_qn_factory_c_wrapper( &
        n_param_c, update_orbs_c_funptr, obj_func_c_funptr, transport_c_funptr, &
        init_hess_c_funptr, update_orbs_qn_c_funptr, obj_func_qn_c_funptr, &
        solver_settings_c, settings_c) result(error_c) bind(C, name="mock_qn_factory")
        !
        ! this subroutine is a mock routine for the quasi-Newton orbital updating
        ! factory C wrapper subroutine
        !
        use, intrinsic :: iso_c_binding, only: c_null_ptr
        use otr_qn_c_interface, only: qn_settings_type_c
        use c_interface, only: logger_c_type, solver_settings_type_c
        use test_reference, only: check_update_orbs_c_funptr, check_obj_func_c_funptr, &
                                  n_param
        use otr_qn_test_reference, only: test_transport_c_funptr, &
                                         test_init_hess_c_funptr, operator(/=)
        use otr_common_c_interface_mock, only: &
            mock_level_update_orbs, mock_level_obj_func, mock_wire_solver_settings_c

        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             transport_c_funptr, init_hess_c_funptr
        type(c_funptr), intent(out) :: update_orbs_qn_c_funptr, obj_func_qn_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(qn_settings_type_c), intent(in) :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(len=:), allocatable, target :: message

        ! check if passed number of parameters is correct
        if (n_param_c /= n_param) then
            write(stderr, *) "test_qn_factory_py_interface failed: Passed number "// &
                "of parameters wrong."
            test_qn_factory_interface = .false.
        end if

        ! test passed orbital update function
        if (.not. check_update_orbs_c_funptr( &
            update_orbs_c_funptr, "qn_factory_py_interface", " by given orbital "// &
            "updating function", c_null_ptr)) test_qn_factory_interface = .false.

        ! test passed objective function
        if (.not. check_obj_func_c_funptr( &
            obj_func_c_funptr, "qn_factory_py_interface", " by given objective "// &
            "function", c_null_ptr)) test_qn_factory_interface = .false.

        ! test passed transport function
        if (.not. test_transport_c_funptr( &
            transport_c_funptr, "qn_factory_py_interface", " by given transport "// &
            "function", c_null_ptr)) test_qn_factory_interface = .false.

        ! test passed initial Hessian function
        if (.not. test_init_hess_c_funptr( &
            init_hess_c_funptr, "qn_factory_py_interface", " by given initial "// &
            "Hessian function", c_null_ptr)) test_qn_factory_interface = .false.

        ! get Fortran pointer to passed logging function and call it
        message = "test"//c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message, c_null_ptr)

        ! check optional settings against reference values
        if (settings_c /= ref_qn_settings) then
            write(stderr, *) "test_qn_factory_py_interface failed: Passed settings "// &
                "associated with wrong values."
            test_qn_factory_interface = .false.
        end if

        ! check the passed solver settings and wire them
        if (.not. mock_wire_solver_settings_c(solver_settings_c, &
                                              "qn_factory_py_interface")) &
            test_qn_factory_interface = .false.

        ! set function pointers to the mock functions of the level of the solver
        update_orbs_qn_c_funptr = c_funloc(mock_level_update_orbs)
        obj_func_qn_c_funptr = c_funloc(mock_level_obj_func)

        ! set return arguments
        error_c = 0

    end function mock_qn_factory_c_wrapper

    subroutine mock_init_qn_settings_c(settings) bind(C, name="mock_init_qn_settings")
        !
        ! this subroutine is a mock routine for the C quasi-Newton setting
        ! initialization subroutine
        !
        use otr_qn_c_interface, only: qn_settings_type_c
        use otr_qn_test_reference, only: assignment(=)

        type(qn_settings_type_c), intent(inout) :: settings

        ! set reference values
        settings = ref_qn_settings

    end subroutine mock_init_qn_settings_c

    subroutine mock_update_orbs_qn_deconstructor_c_wrapper() &
        bind(C, name="mock_update_orbs_qn_deconstructor")
        !
        ! this subroutine is a mock routine for the C quasi-Newton deconstructor
        ! subroutine
        !
        test_update_orbs_qn_deconstructor_interface = .true._c_bool

    end subroutine mock_update_orbs_qn_deconstructor_c_wrapper

end module otr_qn_c_interface_mock
