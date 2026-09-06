! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_mock

    use opentrustregion, only: rp, ip, stderr, update_orbs_type, hess_x_type

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: mock_update_orbs_ptr => mock_update_orbs
    procedure(hess_x_type), pointer :: mock_hess_x_ptr => mock_hess_x

contains

    subroutine mock_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this subroutine is a test subroutine for the orbital update function
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        func = sum(kappa)

        grad = 2 * kappa

        h_diag = 3 * kappa

        hess_x_funptr => mock_hess_x

        error = 0

    end subroutine mock_update_orbs

    subroutine mock_hess_x(x, hess_x, error)
        !
        ! this subroutine is a test subroutine for the Hessian linear transformation 
        ! function
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        hess_x = 4 * x

        error = 0

    end subroutine mock_hess_x

end module otr_common_mock
