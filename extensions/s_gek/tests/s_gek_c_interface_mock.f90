! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface_mock

    use opentrustregion_core, only: stderr
    use c_interface, only: c_rp, c_ip
    use otr_s_gek_c_interface, only: s_gek_solver_c_wrapper, init_s_gek_settings_c
    use otr_s_gek_test_reference, only: ref_s_gek_settings
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_f_pointer, &
                                           c_f_procpointer, c_associated, c_null_char, &
                                           c_null_funptr, c_null_ptr

    implicit none

    ! result of the mock S-GEK solver, which the Python test clears before and reads
    ! after every call whose arguments it checks, and whether the mock received solver
    ! settings and S-GEK settings, which the Python test reads
    logical(c_bool), bind(C) :: test_s_gek_solver_interface = .false._c_bool, &
                                s_gek_solver_settings_passed = .false._c_bool, &
                                s_gek_settings_passed = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(s_gek_solver_c_wrapper), pointer :: mock_s_gek_solver_c_wrapper_ptr => &
        mock_s_gek_solver_c_wrapper
    procedure(init_s_gek_settings_c), pointer :: mock_init_s_gek_settings_c_ptr => &
        mock_init_s_gek_settings_c

contains

    function mock_s_gek_solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                         change_reference_c_funptr, n_param_c, &
                                         settings_c, state_c_ptr, s_gek_settings_c) &
        result(error_c) bind(C, name="mock_s_gek_solver")
        !
        ! this function is a mock routine for the S-GEK solver C wrapper function
        !
        use c_interface, only: solver_settings_type_c, solver_state_type_c, &
                               logger_c_type, obj_func_c_type
        use otr_s_gek_c_interface, only: s_gek_settings_type_c
        use test_reference, only: check_update_orbs_c_funptr, check_obj_func_c_funptr, &
                                  check_hess_x_c_funptr, &
                                  get_reference_solver_state_values, n_param
        use otr_common_test_reference, only: check_obj_func_grad_c_funptr
        use otr_s_gek_test_reference, only: check_change_reference_c_funptr, &
                                            operator(/=)

        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             change_reference_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(solver_settings_type_c), intent(in), optional :: settings_c
        type(c_ptr), intent(in), value :: state_c_ptr
        type(s_gek_settings_type_c), intent(in), optional :: s_gek_settings_c
        integer(c_ip) :: error_c

        type(solver_state_type_c), pointer :: state_c
        procedure(logger_c_type), pointer :: logger_funptr
        procedure(obj_func_c_type), pointer :: obj_func_funptr
        character(len=:), allocatable, target :: message
        real(c_rp) :: kappa(n_param), func

        ! assume tests pass
        test_s_gek_solver_interface = .true.

        ! record whether solver settings and S-GEK settings were passed
        s_gek_solver_settings_passed = present(settings_c)
        s_gek_settings_passed = present(s_gek_settings_c)

        ! evaluate the objective function once and return its error as the S-GEK solver
        ! would, so that a test can check how an error of a callback function is
        ! reported
        kappa = 0.0_c_rp
        call c_f_procpointer(cptr=obj_func_c_funptr, fptr=obj_func_funptr)
        error_c = obj_func_funptr(kappa, func, c_null_ptr)
        if (error_c /= 0) return

        ! test passed orbital update function
        if (.not. check_update_orbs_c_funptr( &
            update_orbs_c_funptr, "s_gek_solver_py_interface", " by given orbital "// &
            "updating function", c_null_ptr)) test_s_gek_solver_interface = .false.

        ! test passed objective function
        if (.not. check_obj_func_c_funptr( &
            obj_func_c_funptr, "s_gek_solver_py_interface", " by given objective "// &
            "function", c_null_ptr)) test_s_gek_solver_interface = .false.

        ! test passed chart transformation
        if (.not. check_change_reference_c_funptr( &
            change_reference_c_funptr, "s_gek_solver_py_interface", " by given "// &
            "chart transformation", c_null_ptr)) test_s_gek_solver_interface = .false.

        ! check if passed number of parameters is correct
        if (n_param_c /= n_param) then
            write(stderr, *) "test_s_gek_solver_py_interface failed: Passed number "// &
                "of parameters wrong."
            test_s_gek_solver_interface = .false.
        end if

        ! check the callback functions and the values of the S-GEK settings where these
        ! were passed and call their logging function
        if (present(s_gek_settings_c)) then
            if (.not. check_obj_func_grad_c_funptr( &
                s_gek_settings_c%obj_func_grad, "s_gek_solver_py_interface", &
                " by the objective function of the S-GEK settings", c_null_ptr)) &
                test_s_gek_solver_interface = .false.
            if (.not. check_hess_x_c_funptr( &
                s_gek_settings_c%init_hess, "s_gek_solver_py_interface", &
                " by the initial Hessian of the S-GEK settings", c_null_ptr)) &
                test_s_gek_solver_interface = .false.

            ! get Fortran pointer to passed logging function and call it
            message = "test"//c_null_char
            call c_f_procpointer(cptr=s_gek_settings_c%logger, fptr=logger_funptr)
            call logger_funptr(message, c_null_ptr)

            ! check S-GEK settings against reference values
            if (s_gek_settings_c /= ref_s_gek_settings) then
                write(stderr, *) "test_s_gek_solver_py_interface failed: Passed "// &
                    "S-GEK settings associated with wrong values."
                test_s_gek_solver_interface = .false.
            end if
        end if

        ! fill a given state with the reference values, so that a test can check that
        ! the Python interface passes its state
        if (c_associated(state_c_ptr)) then
            call c_f_pointer(state_c_ptr, state_c)
            call get_reference_solver_state_values(state_c)
        end if

        ! set return arguments
        error_c = 0

    end function mock_s_gek_solver_c_wrapper

    subroutine mock_init_s_gek_settings_c(settings) &
        bind(C, name="mock_init_s_gek_settings")
        !
        ! this subroutine is a mock routine for the C S-GEK setting initialization
        ! subroutine
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c
        use otr_s_gek_test_reference, only: get_reference_s_gek_values

        type(s_gek_settings_type_c), intent(out) :: settings

        ! set reference values without callback functions and host context
        call get_reference_s_gek_values(settings)
        settings%obj_func_grad = c_null_funptr
        settings%init_hess = c_null_funptr
        settings%logger = c_null_funptr
        settings%context = c_null_ptr

    end subroutine mock_init_s_gek_settings_c

end module otr_s_gek_c_interface_mock
