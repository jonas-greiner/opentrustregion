! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_test_reference

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: n_param, tol, tol_c
    use, intrinsic :: iso_c_binding, only: c_funptr, c_f_procpointer

    implicit none

contains

    function test_change_reference_funptr(change_reference_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided change reference function pointer
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_common, only: change_reference_type

        procedure(change_reference_type), intent(in), pointer :: change_reference_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: new_ref(:), kappa_list(:, :), local_grad_list(:, :), &
                                 grad_list(:, :)
        integer(ip) :: n_points, error

        ! assume tests pass
        test_passed = .true.

        ! initialize number of points
        n_points = 4

        ! allocate arrays
        allocate(new_ref(n_param), kappa_list(n_param, n_points), &
                 local_grad_list(n_param, n_points), grad_list(n_param, n_points))

        ! initialize input arrays
        new_ref = 1.0_rp
        kappa_list = 2.0_rp
        local_grad_list = 3.0_rp
        grad_list = 4.0_rp

        ! call change reference function pointer
        call change_reference_funptr(new_ref, n_points, kappa_list, local_grad_list, &
                                     grad_list, error)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check kappa list
        if (any((kappa_list - 4.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Parameter list returned"// &
                message//" wrong."
        end if

        ! check local gradient list
        if (any((local_grad_list - 9.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Local gradient list "// &
                "returned"//message//" wrong."
        end if

        ! check gradient list
        if (any((grad_list - 16.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Gradient list returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(new_ref, kappa_list, local_grad_list, grad_list)

    end function test_change_reference_funptr

    function test_change_reference_c_funptr(change_reference_c_funptr, test_name, &
                                            message) result(test_passed)
        !
        ! this function tests a provided change reference C function pointer
        !
        use otr_common_c_interface, only: change_reference_c_type

        type(c_funptr), intent(in) :: change_reference_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(change_reference_c_type), pointer :: change_reference_funptr
        real(c_rp), allocatable :: new_ref(:), kappa_list(:, :), &
                                   local_grad_list(:, :), grad_list(:, :)
        integer(c_ip) :: n_points, error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=change_reference_c_funptr, &
                             fptr=change_reference_funptr)

        ! initialize number of points
        n_points = 4_c_ip

        ! allocate arrays
        allocate(new_ref(n_param), kappa_list(n_param, n_points), &
                 local_grad_list(n_param, n_points), grad_list(n_param, n_points))

        ! initialize new reference
        new_ref = 1.0_c_rp
        kappa_list = 2.0_c_rp
        local_grad_list = 3.0_c_rp
        grad_list = 4.0_c_rp

        ! call change reference function pointer
        error = change_reference_funptr(new_ref, n_points, kappa_list, &
                                        local_grad_list, grad_list)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check kappa list
        if (any((kappa_list - 4.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Parameter list returned"// &
                message//" wrong."
        end if

        ! check local gradient list
        if (any((local_grad_list - 9.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Local gradient list "// &
                "returned"//message//" wrong."
        end if

        ! check gradient list
        if (any((grad_list - 16.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Gradient list returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(new_ref, kappa_list, local_grad_list, grad_list)

    end function test_change_reference_c_funptr

end module otr_common_test_reference
