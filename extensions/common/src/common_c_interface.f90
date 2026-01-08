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

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function change_reference_c_type(new_ref_c, n_points_c, kappa_list_c, &
                                         local_grad_list_c, grad_list_c) &
            result(error_c) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in), target :: new_ref_c(*)
            integer(c_ip), intent(in), value :: n_points_c
            real(c_rp), intent(inout), target :: kappa_list_c(*), &
                                                 local_grad_list_c(*), grad_list_c(*)
            integer(c_ip) :: error_c
        end function change_reference_c_type
    end interface

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

    subroutine change_reference_f_wrapper_impl(change_reference_funptr, new_ref, &
                                               n_points, kappa_list, local_grad_list, &
                                               grad_list, error)
        !
        ! this subroutine wraps the change reference subroutine to convert Fortran 
        ! variables to C variables for a given function pointer
        !
        procedure(change_reference_c_type), intent(in), pointer :: &
            change_reference_funptr
        real(rp), intent(in), target :: new_ref(:)
        integer(ip), intent(in) :: n_points
        real(rp), intent(inout), target :: kappa_list(:, :), local_grad_list(:, :), &
                                           grad_list(:, :)
        integer(ip), intent(out) :: error

        real(c_rp), pointer :: new_ref_c(:), kappa_list_c(:, :), &
                               local_grad_list_c(:, :), grad_list_c(:, :)
        integer(c_ip) :: n_points_c, error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            new_ref_c => new_ref
            n_points_c = n_points
            kappa_list_c => kappa_list
            local_grad_list_c => local_grad_list
            grad_list_c => grad_list
        else
            allocate(new_ref_c(size(new_ref, 1)))
            allocate(kappa_list_c(size(kappa_list, 1), size(kappa_list, 2)))
            allocate(local_grad_list_c(size(local_grad_list, 1), &
                                       size(local_grad_list, 2)))
            allocate(grad_list_c(size(grad_list, 1), size(grad_list, 2)))
            new_ref_c = real(new_ref, kind=c_rp)
            n_points_c = int(n_points, kind=c_ip)
            kappa_list_c = real(kappa_list, kind=c_rp)
            local_grad_list_c = real(local_grad_list, kind=c_rp)
            grad_list_c = real(grad_list, kind=c_rp)
        end if

        ! call change reference C function
        error_c = change_reference_funptr(new_ref_c, n_points_c, kappa_list_c, &
                                          local_grad_list_c, grad_list_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            kappa_list = real(kappa_list_c, kind=rp)
            local_grad_list = real(local_grad_list_c, kind=rp)
            grad_list = real(grad_list_c, kind=rp)
            deallocate(new_ref_c)
            deallocate(kappa_list_c)
            deallocate(local_grad_list_c)
            deallocate(grad_list_c)
        end if

    end subroutine change_reference_f_wrapper_impl

end module otr_common_c_interface
