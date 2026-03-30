! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_test_reference

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
     use test_reference, only: n_param, tol, tol_c
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer

    implicit none

    ! derived types for S-GEK settings
    type ref_s_gek_settings_type
        logical :: use_subspace
        integer(ip) :: verbose, max_points
    end type

    type, bind(C) :: ref_s_gek_settings_type_c
        logical(c_bool) :: use_subspace
        integer(c_ip) :: verbose, max_points
    end type

    ! general reference parameters
    type(ref_s_gek_settings_type) :: ref_s_gek_settings = &
        ref_s_gek_settings_type(use_subspace = .false., verbose = 3, max_points = 10)

    interface assignment(=)
        module procedure assign_ref_to_s_gek
        module procedure assign_ref_to_s_gek_c
        module procedure assign_ref_to_ref_c
    end interface

    interface operator(==)
        module procedure equal_s_gek_to_ref
        module procedure equal_s_gek_c_to_ref
        module procedure equal_s_gek
        module procedure equal_s_gek_c
    end interface

    interface operator(/=)
        module procedure not_equal_s_gek_to_ref
        module procedure not_equal_s_gek_c_to_ref
        module procedure not_equal_s_gek
        module procedure not_equal_s_gek_c
    end interface

contains

    function test_change_reference_funptr(change_reference_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided change reference function pointer
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_s_gek, only: change_reference_type

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
        use otr_s_gek_c_interface, only: change_reference_c_type

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

    subroutine get_reference_s_gek_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the S-GEK reference values for tests
        !
        type(ref_s_gek_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_s_gek_settings

    end subroutine get_reference_s_gek_values

    subroutine assign_ref_to_s_gek(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set S-GEK settings to 
        ! reference values
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(out) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%logger => null()

        ! set reference values
        lhs%use_subspace = rhs%use_subspace
        lhs%verbose = rhs%verbose
        lhs%max_points = rhs%max_points

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_s_gek

    subroutine assign_ref_to_s_gek_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C S-GEK settings to 
        ! reference values
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(out) :: lhs_c
        type(ref_s_gek_settings_type), intent(in) :: rhs

        type(s_gek_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_s_gek_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        use c_interface, only: character_to_c

        type(ref_s_gek_settings_type_c), intent(out) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs

        lhs%use_subspace = logical(rhs%use_subspace, kind=c_bool)
        lhs%verbose = int(rhs%verbose, kind=c_ip)
        lhs%max_points = int(rhs%max_points, kind=c_ip)

    end subroutine assign_ref_to_ref_c

    logical function equal_s_gek_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to 
        ! reference values
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(in) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs

        equal_s_gek_to_ref = (lhs%use_subspace .eqv. rhs%use_subspace) .and. &
                             lhs%verbose == rhs%verbose .and. lhs%max_points == &
                             rhs%max_points

    end function equal_s_gek_to_ref

    logical function not_equal_s_gek_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK 
        ! settings to reference values
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(in) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs

        not_equal_s_gek_to_ref = .not. (lhs == rhs)

    end function not_equal_s_gek_to_ref

    logical function equal_s_gek_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to
        ! reference values
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(ref_s_gek_settings_type), intent(in) :: rhs

        type(s_gek_settings_type) :: lhs

        lhs = lhs_c
        equal_s_gek_c_to_ref = lhs == rhs

    end function equal_s_gek_c_to_ref

    logical function not_equal_s_gek_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK
        ! settings to reference values
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c

        type(s_gek_settings_type_c), intent(in) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs
        
        not_equal_s_gek_c_to_ref = .not. (lhs == rhs)

    end function not_equal_s_gek_c_to_ref

    logical function equal_s_gek(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to 
        ! different S-GEK settings
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(in) :: lhs, rhs
        
        equal_s_gek = (lhs%use_subspace .eqv. rhs%use_subspace) .and. lhs%verbose == &
                      rhs%verbose .and. lhs%max_points == rhs%max_points

    end function equal_s_gek

    logical function not_equal_s_gek(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK
        ! settings to different S-GEK settings
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(in) :: lhs, rhs
        
        not_equal_s_gek = .not. (lhs == rhs)

    end function not_equal_s_gek

    logical function equal_s_gek_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to
        ! settings to different S-GEK settings
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(s_gek_settings_type), intent(in) :: rhs
        
        type(s_gek_settings_type) :: lhs

        lhs = lhs_c
        equal_s_gek_c = lhs == rhs

    end function equal_s_gek_c

    logical function not_equal_s_gek_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK
        ! settings to different S-GEK settings
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(s_gek_settings_type), intent(in) :: rhs
        
        not_equal_s_gek_c = .not. (lhs_c == rhs)

    end function not_equal_s_gek_c

end module otr_s_gek_test_reference
