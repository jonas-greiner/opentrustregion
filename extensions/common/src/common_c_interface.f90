! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface

    use opentrustregion, only: ip, rp
    use c_interface, only: c_ip, c_rp
    use, intrinsic :: iso_c_binding, only: c_funptr, c_funloc

    implicit none

    ! global variables
    integer(ip) :: n_param

contains

    function update_orbs_c_wrapper_impl(update_orbs_before_wrapping, &
                                        hess_x_before_wrapping_funptr, &
                                        hess_x_c_wrapper_funptr, kappa_c, func_c, &
                                        grad_c, h_diag_c, hess_x_c_funptr) &
        result(error_c)
        !
        ! this function wraps the orbital update subroutine to convert Fortran 
        ! variables to C variables
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use c_interface, only: hess_x_c_type

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_before_wrapping
        procedure(hess_x_type), intent(out), pointer :: &
            hess_x_before_wrapping_funptr
        procedure(hess_x_c_type) :: hess_x_c_wrapper_funptr
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(out) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        real(rp) :: func
        real(rp), pointer :: kappa(:), grad(:), h_diag(:)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            kappa => kappa_c(:n_param)
            grad => grad_c(:n_param)
            h_diag => h_diag_c(:n_param)
        else
            allocate(kappa(n_param))
            allocate(grad(n_param))
            allocate(h_diag(n_param))
            kappa = real(kappa_c(:n_param), kind=rp)
        end if

        ! call update_orbs Fortran function
        call update_orbs_before_wrapping(kappa, func, grad, h_diag, &
                                         hess_x_before_wrapping_funptr, error)

        ! convert arguments to Fortran kind
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            grad_c(:n_param) = real(grad, kind=c_rp)
            h_diag_c(:n_param) = real(h_diag, kind=c_rp)
            deallocate(kappa)
            deallocate(grad)
            deallocate(h_diag)
        end if

        ! get a C function pointer to the hess_x wrapper function
        hess_x_c_funptr = c_funloc(hess_x_c_wrapper_funptr)

    end function update_orbs_c_wrapper_impl

    function hess_x_c_wrapper_impl(hess_x_funptr, x_c, hess_x_c) result(error_c)
        !
        ! this function wraps the Hessian linear transformation to convert Fortran 
        ! variables to C variables
        !
        use opentrustregion, only: hess_x_type

        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        integer(c_ip) :: error_c

        real(rp), pointer :: x(:), hess_x(:)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            x => x_c(:n_param)
            hess_x => hess_x_c(:n_param)
        else
            allocate(x(n_param))
            allocate(hess_x(n_param))
            x = real(x_c(:n_param), kind=rp)
        end if

        ! call Fortran function
        call hess_x_funptr(x, hess_x, error)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            hess_x_c(:n_param) = real(hess_x, kind=c_rp)
            deallocate(x)
            deallocate(hess_x)
        end if

    end function hess_x_c_wrapper_impl

end module otr_common_c_interface
