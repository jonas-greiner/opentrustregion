! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip, update_orbs_c_type, obj_func_c_type, &
                           hess_x_c_type, precond_c_type
    use, intrinsic :: iso_c_binding, only: c_ptr, c_funptr, c_funloc, c_loc, &
                                           c_f_pointer, c_associated, c_null_ptr

    implicit none

    ! levels whose addresses the mock factories write into the contexts of the C solver
    ! settings and of their nested settings, with which the Python interface has to
    ! call the functions they return or wire
    logical, target :: mock_solver_level, mock_stability_level

    ! record of the solver state the latest mock factory received, extension it
    ! returned, extension the latest mock release received and number of mock releases,
    ! which the Python tests compare with the record they passed and the extensions the
    ! factories returned
    type(c_ptr), bind(C) :: received_solver_state = c_null_ptr, &
                            returned_extension = c_null_ptr, &
                            released_extension = c_null_ptr
    integer(c_ip), bind(C) :: n_released_extensions = 0

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_c_type), pointer :: mock_level_update_orbs_ptr => &
        mock_level_update_orbs
    procedure(hess_x_c_type), pointer :: mock_level_hess_x_ptr => mock_level_hess_x
    procedure(obj_func_c_type), pointer :: mock_level_obj_func_ptr => &
        mock_level_obj_func
    procedure(precond_c_type), pointer :: &
        mock_level_precond_ptr => mock_level_precond, &
        mock_level_stability_precond_ptr => mock_level_stability_precond

contains

    function mock_level_update_orbs(kappa, func, grad, h_diag, hess_x_c_funptr, &
                                    context_c) result(error) bind(C)
        !
        ! this function is a mock of the orbital update a factory returns, which
        ! returns the mock Hessian linear transformation of the level of the solver and
        ! produces an error unless it is called with that level
        !
        use c_interface_unit_tests, only: mock_update_orbs

        real(c_rp), intent(in), target :: kappa(*)
        real(c_rp), intent(out) :: func
        real(c_rp), intent(out), target :: grad(*), h_diag(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_update_orbs(kappa, func, grad, h_diag, hess_x_c_funptr, c_null_ptr)
        hess_x_c_funptr = c_funloc(mock_level_hess_x)
        if (.not. c_associated(context_c, c_loc(mock_solver_level))) error = 1

    end function mock_level_update_orbs

    function mock_level_hess_x(x, hess_x, context_c) result(error) bind(C)
        !
        ! this function is a mock of the Hessian linear transformation the orbital
        ! update of a factory returns, which produces an error unless it is called with
        ! the level of the solver
        !
        use c_interface_unit_tests, only: mock_hess_x

        real(c_rp), intent(in), target :: x(*)
        real(c_rp), intent(out), target :: hess_x(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_hess_x(x, hess_x, c_null_ptr)
        if (.not. c_associated(context_c, c_loc(mock_solver_level))) error = 1

    end function mock_level_hess_x

    function mock_level_obj_func(kappa, func, context_c) result(error) bind(C)
        !
        ! this function is a mock of the objective function a factory returns, which
        ! produces an error unless it is called with the level of the solver
        !
        use c_interface_unit_tests, only: mock_obj_func

        real(c_rp), intent(in), target :: kappa(*)
        real(c_rp), intent(out) :: func
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_obj_func(kappa, func, c_null_ptr)
        if (.not. c_associated(context_c, c_loc(mock_solver_level))) error = 1

    end function mock_level_obj_func

    function mock_level_precond(residual, mu, precond_residual, context_c) &
        result(error) bind(C)
        !
        ! this function is a mock of the preconditioner a factory wires into the solver
        ! settings, which produces an error unless it is called with the level of the
        ! solver
        !
        use c_interface_unit_tests, only: mock_precond

        real(c_rp), intent(in), target :: residual(*)
        real(c_rp), intent(in) :: mu
        real(c_rp), intent(out), target :: precond_residual(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_precond(residual, mu, precond_residual, c_null_ptr)
        if (.not. c_associated(context_c, c_loc(mock_solver_level))) error = 1

    end function mock_level_precond

    function mock_level_stability_precond(residual, mu, precond_residual, context_c) &
        result(error) bind(C)
        !
        ! this function is a mock of the preconditioner a factory wires into the nested
        ! stability check settings, which produces an error unless it is called with
        ! the level of the stability check
        !
        use c_interface_unit_tests, only: mock_precond

        real(c_rp), intent(in), target :: residual(*)
        real(c_rp), intent(in) :: mu
        real(c_rp), intent(out), target :: precond_residual(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_precond(residual, mu, precond_residual, c_null_ptr)
        if (.not. c_associated(context_c, c_loc(mock_stability_level))) error = 1

    end function mock_level_stability_precond

    function mock_wire_solver_settings_c(solver_settings_c, test_name) &
        result(test_passed)
        !
        ! this function checks that a factory was handed the C solver settings of the
        ! Python test with the callback functions of its Python functions, the
        ! projection and logging function of the solver and the convergence check of
        ! the stability check, and without a context, and wires them as a factory does,
        ! with the addresses of the levels as contexts and the mock preconditioners of
        ! both levels
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        character(len=*), intent(in) :: test_name
        logical :: test_passed

        ! assume tests pass
        test_passed = .true.

        ! check the callback functions and the context of the solver settings
        if (.not. (c_associated(solver_settings_c%project) .and. &
                   c_associated(solver_settings_c%logger) .and. &
                   c_associated(solver_settings_c%stability_settings%conv_check))) then
            write(stderr, *) "test_"//test_name// &
                " failed: Callback functions of Python functions not passed."
            test_passed = .false.
        end if
        if (c_associated(solver_settings_c%context)) then
            write(stderr, *) "test_"//test_name// &
                " failed: Context of Python settings passed."
            test_passed = .false.
        end if

        ! wire the solver settings
        solver_settings_c%context = c_loc(mock_solver_level)
        solver_settings_c%precond = c_funloc(mock_level_precond)
        solver_settings_c%stability_settings%context = c_loc(mock_stability_level)
        solver_settings_c%stability_settings%precond = &
            c_funloc(mock_level_stability_precond)

    end function mock_wire_solver_settings_c

    function mock_extension(solver_state_c) result(extension_c) &
        bind(C, name="mock_extension")
        !
        ! this function records the record of the solver state a mock factory received
        ! and returns and records an empty holder as the extension it returns, which
        ! the real release frees, so that the Python owner of the extension can release
        ! it at any time
        !
        use otr_common_c_interface, only: extension_holder_type

        type(c_ptr), intent(in), value :: solver_state_c
        type(c_ptr) :: extension_c

        type(extension_holder_type), pointer :: holder

        received_solver_state = solver_state_c
        allocate(holder)
        extension_c = c_loc(holder)
        returned_extension = extension_c

    end function mock_extension

    subroutine mock_release_extension(extension_c) &
        bind(C, name="mock_release_extension")
        !
        ! this subroutine is a mock routine for the release of an extension, which
        ! records the extension it received and frees the holder the mock factory
        ! allocated
        !
        use otr_common_c_interface, only: extension_holder_type

        type(c_ptr), intent(in), value :: extension_c

        type(extension_holder_type), pointer :: holder

        released_extension = extension_c
        n_released_extensions = n_released_extensions + 1
        if (c_associated(extension_c)) then
            call c_f_pointer(extension_c, holder)
            deallocate(holder)
        end if

    end subroutine mock_release_extension

end module otr_common_c_interface_mock
