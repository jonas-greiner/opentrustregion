! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_mock

    use opentrustregion, only: rp, ip, obj_func_type, update_orbs_type, hess_x_type, &
                               precond_type, precond_pd_type, &
                               get_extra_trial_vectors_type

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(obj_func_type), pointer :: mock_obj_func_ptr => mock_obj_func
    procedure(update_orbs_type), pointer :: mock_update_orbs_ptr => mock_update_orbs
    procedure(hess_x_type), pointer :: mock_hess_x_ptr => mock_hess_x
    procedure(precond_type), pointer :: mock_precond_ptr => mock_precond
    procedure(precond_pd_type), pointer :: mock_precond_pd_ptr => mock_precond_pd
    procedure(get_extra_trial_vectors_type), pointer :: &
        mock_get_extra_trial_vectors_ptr => mock_get_extra_trial_vectors

    ! name of the mock host callback function called last, so that tests can check to
    ! which host callback function a callback function forwarded
    character(len=40) :: host_called = ""

contains

    function record_host_call(name, context) result(error)
        !
        ! this function records the call of a mock host callback function, checks the
        ! host context it received and returns the error code the test set in it, for
        ! the call the test chose or every call
        !
        use test_reference, only: host_context_type, check_host_context

        character(len=*), intent(in) :: name
        class(*), intent(in), pointer :: context
        integer(ip) :: error

        host_called = name
        call check_host_context(context)
        error = 0
        if (associated(context)) then
            select type (context)
            class is (host_context_type)
                if (context%failing_call == 0 .or. &
                    context%n_calls == context%failing_call) &
                    error = int(context%mock_error, kind=ip)
            end select
        end if

    end function record_host_call

    function mock_obj_func(kappa, error, context) result(func)
        !
        ! this function is a test function for the objective function, which records
        ! its call, checks the host context it receives and returns the error code set
        ! in it
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        func = sum(kappa)

        error = record_host_call("obj_func", context)

    end function mock_obj_func

    subroutine mock_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                context)
        !
        ! this subroutine is a test subroutine for the orbital update function, which
        ! records its call, checks the host context it receives and returns the error
        ! code set in it
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

        error = record_host_call("update_orbs", context)

    end subroutine mock_update_orbs

    subroutine mock_hess_x(x, hess_x, error, context)
        !
        ! this subroutine is a test subroutine for the Hessian linear transformation
        ! function, which records its call, checks the host context it receives and
        ! returns the error code set in it
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        hess_x = 4 * x

        error = record_host_call("hess_x", context)

    end subroutine mock_hess_x

    subroutine mock_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine is a test subroutine for the level-shifted preconditioner
        ! function, which records its call, checks the host context it receives and
        ! returns the error code set in it
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        precond_residual = mu * residual

        error = record_host_call("precond", context)

    end subroutine mock_precond

    subroutine mock_precond_pd(residual, precond_residual, error, context)
        !
        ! this subroutine is a test subroutine for the positive-definite preconditioner
        ! function, which records its call, checks the host context it receives and
        ! returns the error code set in it
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        precond_residual = 3.0_rp * residual

        error = record_host_call("precond_pd", context)

    end subroutine mock_precond_pd

    subroutine mock_get_extra_trial_vectors(trial_vectors, error, context)
        !
        ! this subroutine is a test subroutine for the extra trial vector function,
        ! which records its call, checks the host context it receives and returns the
        ! error code set in it
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: i

        do i = 1, size(trial_vectors, 2)
            trial_vectors(:, i) = real(i, kind=rp)
        end do

        error = record_host_call("get_extra_trial_vectors", context)

    end subroutine mock_get_extra_trial_vectors

end module otr_common_mock
