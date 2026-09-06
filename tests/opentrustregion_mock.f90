! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_mock

    use opentrustregion, only: rp, ip, stderr, solver, stability_check, &
                               update_orbs_type, hess_x_type, obj_func_type
    use test_reference, only: tol, ref_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! create function pointers to ensure that routines comply with interface
    procedure(solver), pointer :: mock_solver_ptr => mock_solver
    procedure(stability_check), pointer :: mock_stability_check_ptr => &
        mock_stability_check

contains

    subroutine mock_solver(update_orbs_funptr, obj_func_funptr, n_param, error, &
                           settings)
        !
        ! this subroutine is a mock routine for solver to test the C interface
        !
        use opentrustregion, only: solver_settings_type
        use test_reference, only: test_update_orbs_funptr, test_obj_func_funptr, &
                                  test_mock_solver_funptr, test_mock_stability_funptr

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        type(solver_settings_type), intent(inout) :: settings

        ! initialize logical
        test_passed = .true.

        ! test passed orbital update subroutine
        test_passed = test_passed .and. &
            test_update_orbs_funptr(update_orbs_funptr, "solver_c_wrapper", &
                                    " by given orbital updating subroutine")

        ! test passed objective function
        test_passed = test_passed .and. &
            test_obj_func_funptr(obj_func_funptr, "solver_c_wrapper", &
                                 " by given objective function")

        ! check number of parameters
        if (n_param /= 3) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed number of "// &
                "parameters wrong."
        end if

        ! set output quantities
        error = 0

        ! check optional function pointers
        if (.not. test_mock_solver_funptr(settings, "solver_c_wrapper", " given")) &
            test_passed = .false.

        ! check if optional settings are correctly passed
        if (settings /= ref_settings) then
            test_passed = .false.
            write (stderr, *) "test_solver_c_wrapper failed: Passed optional "// &
                "settings associated with wrong values."
        end if

    end subroutine mock_solver

    subroutine mock_stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                                    kappa, min_eigval)
        !
        ! this subroutine is a mock routine for the stability check to test the C 
        ! interface
        !
        use opentrustregion, only: stability_settings_type
        use test_reference, only: test_hess_x_funptr, test_mock_stability_funptr

        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable
        integer(ip), intent(out) :: error
        type(stability_settings_type), intent(inout) :: settings
        real(rp), intent(out), optional :: kappa(:), min_eigval

        ! initialize logical
        test_passed = .true.

        ! check Hessian diagonal
        if (any(abs(h_diag - 3.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "Hessian diagonal wrong."
        end if

        ! test passed Hessian linear transformation subroutine
        test_passed = test_passed .and. &
            test_hess_x_funptr(hess_x_funptr, "stability_check_c_wrapper", &
                               " by given Hessian linear transformation subroutine")

        ! set output quantities
        stable = .false.
        if (present(kappa)) kappa = 1.0_rp
        if (present(min_eigval)) min_eigval = -1.0_rp
        error = 0

        ! check optional function pointers
        if (.not. test_mock_stability_funptr(settings, "stability_check_c_wrapper", &
                                             " given")) test_passed = .false.

        ! check if optional settings are correctly passed
        if (settings /= ref_settings) then
            test_passed = .false.
            write (stderr, *) "test_stability_check_c_wrapper failed: Passed "// &
                "optional settings associated with wrong values."
        end if

    end subroutine mock_stability_check

end module opentrustregion_mock
