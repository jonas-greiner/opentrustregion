! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_mock

    use opentrustregion, only: rp, ip, stderr, obj_func_type, update_orbs_type, &
                               hess_x_type, precond_type, precond_pd_type, &
                               get_extra_trial_vectors_type
    use c_interface, only: logger_c_type
    use, intrinsic :: iso_c_binding, only: c_char, c_ptr

    implicit none

    ! whether the C logging function was called with the test message, which is
    ! recorded here rather than in a host context since the C settings of the
    ! extensions do not hand one to their logging function
    logical :: test_logger = .false.

    ! create function pointers to ensure that routines comply with interface
    procedure(obj_func_type), pointer :: mock_obj_func_ptr => mock_obj_func
    procedure(update_orbs_type), pointer :: mock_update_orbs_ptr => mock_update_orbs
    procedure(hess_x_type), pointer :: mock_hess_x_ptr => mock_hess_x
    procedure(precond_type), pointer :: mock_precond_ptr => mock_precond
    procedure(precond_pd_type), pointer :: mock_precond_pd_ptr => mock_precond_pd
    procedure(get_extra_trial_vectors_type), pointer :: &
        mock_get_extra_trial_vectors_ptr => mock_get_extra_trial_vectors
    procedure(logger_c_type), pointer :: mock_logger_c_ptr => mock_logger_c

contains

    function mock_obj_func(kappa, error, context) result(func)
        !
        ! this function is a test function for the objective function
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        func = sum(kappa)

        error = 0

    end function mock_obj_func

    subroutine mock_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                context)
        !
        ! this subroutine is a test subroutine for the orbital update function
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        func = sum(kappa)

        grad = 2 * kappa

        h_diag = 3 * kappa

        hess_x_funptr => mock_hess_x

        error = 0

    end subroutine mock_update_orbs

    subroutine mock_hess_x(x, hess_x, error, context)
        !
        ! this subroutine is a test subroutine for the Hessian linear transformation
        ! function
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        hess_x = 4 * x

        error = 0

    end subroutine mock_hess_x

    subroutine mock_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine is a test subroutine for the level-shifted preconditioner
        ! function
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        precond_residual = mu * residual

        error = 0

    end subroutine mock_precond

    subroutine mock_precond_pd(residual, precond_residual, error, context)
        !
        ! this subroutine is a test subroutine for the positive-definite preconditioner
        ! function
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        precond_residual = 3.0_rp * residual

        error = 0

    end subroutine mock_precond_pd

    subroutine mock_get_extra_trial_vectors(trial_vectors, error, context)
        !
        ! this subroutine is a test subroutine for the extra trial vector function
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: i

        do i = 1, size(trial_vectors, 2)
            trial_vectors(:, i) = real(i, kind=rp)
        end do

        error = 0

    end subroutine mock_get_extra_trial_vectors

    subroutine mock_logger_c(message_c, context_c) bind(C)
        !
        ! this subroutine is a test subroutine for the C logging function, which
        ! records a call with the test message
        !
        use test_reference, only: check_host_context_c

        character(kind=c_char), intent(in) :: message_c(*)
        type(c_ptr), intent(in), value :: context_c

        character(len=4) :: message

        ! check host context
        call check_host_context_c(context_c)

        ! record call with test message
        message = transfer(message_c(1:4), message)
        if (message == "test") test_logger = .true.

    end subroutine mock_logger_c

end module otr_common_mock
