! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common

    use opentrustregion, only: ip, rp

    implicit none

    abstract interface
        subroutine change_reference_type(new_ref, n_points, kappa_list, &
                                         local_grad_list, grad_list, error)
            import :: rp, ip

            real(rp), intent(in), target :: new_ref(:)
            integer(ip), intent(in) :: n_points
            real(rp), intent(inout), target :: kappa_list(:, :), &
                                               local_grad_list(:, :), grad_list(:, :)
            integer(ip), intent(out) :: error

            real(rp) :: energy
        end subroutine change_reference_type
    end interface

contains

end module otr_common
