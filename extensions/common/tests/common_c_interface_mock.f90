! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol_c, n_param

    implicit none

    logical :: test_passed
    character(:), allocatable :: test_name

contains

    function mock_change_reference(new_ref_c, n_points_c, kappa_list_c, &
                                   local_grad_list_c, grad_list_c) result(error) bind(C)
        !
        ! this subroutine is a mock subroutine for the reference change C function
        !
        real(c_rp), intent(in), target :: new_ref_c(*)
        integer(c_ip), intent(in), value :: n_points_c
        real(c_rp), intent(inout), target :: kappa_list_c(*), local_grad_list_c(*), &
                                             grad_list_c(*)
        integer(c_ip) :: error

        ! initialize logical
        test_passed = .true.

        if (n_points_c /= 4_c_ip) then
            test_passed = .false.
            write (stderr, *) test_name//" failed: Number of points inside given "// &
                "change reference function wrong."
        end if

        if (any((new_ref_c(:n_param) - 1.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) test_name//" failed: New reference parameters inside "// &
                "given change reference function wrong."
        end if

        kappa_list_c(:n_param * n_points_c) = 2.0_c_rp * &
            kappa_list_c(:n_param * n_points_c)
        local_grad_list_c(:n_param * n_points_c) = 3.0_c_rp * &
            local_grad_list_c(:n_param * n_points_c)
        grad_list_c(:n_param * n_points_c) = 4.0_c_rp * &
            grad_list_c(:n_param * n_points_c)

        error = 0_c_ip

    end function mock_change_reference

end module otr_common_c_interface_mock
