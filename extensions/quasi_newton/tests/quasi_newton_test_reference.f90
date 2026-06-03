! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_test_reference

    use opentrustregion, only: rp, ip, kw_len, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: n_param, tol, tol_c
    use, intrinsic :: iso_c_binding, only: c_char, c_funptr, c_f_procpointer

    implicit none

    ! derived types for quasi-Newton settings
    type ref_qn_settings_type
        integer(ip) :: verbose, max_points
        character(kw_len, c_char) :: hess_update_scheme
    end type

    type, bind(C) :: ref_qn_settings_type_c
        integer(c_ip) :: verbose, max_points
        character(c_char) :: hess_update_scheme(kw_len + 1)
    end type

    ! general reference parameters
    type(ref_qn_settings_type) :: ref_qn_settings = &
        ref_qn_settings_type(verbose = 3, max_points = 100, hess_update_scheme = "bfgs")

    interface assignment(=)
        module procedure assign_ref_to_qn
        module procedure assign_ref_to_qn_c
        module procedure assign_ref_to_ref_c
    end interface

    interface operator(==)
        module procedure equal_qn_to_ref
        module procedure equal_qn_c_to_ref
        module procedure equal_qn
        module procedure equal_qn_c
    end interface

    interface operator(/=)
        module procedure not_equal_qn_to_ref
        module procedure not_equal_qn_c_to_ref
        module procedure not_equal_qn
        module procedure not_equal_qn_c
    end interface

contains

    function test_transport_funptr(transport_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided transport function pointer
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_qn, only: transport_type

        procedure(transport_type), intent(in), pointer :: transport_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: geodesic(:), tangent_vector(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(geodesic(n_param), tangent_vector(n_param))

        ! initialize input arrays
        geodesic = 2.0_rp
        tangent_vector = 3.0_rp

        ! call transport function pointer
        call transport_funptr(geodesic, tangent_vector, error)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check tangent vector
        if (any((tangent_vector - 6.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Tangent vector returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(geodesic, tangent_vector)

    end function test_transport_funptr

    function test_transport_c_funptr(transport_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided transport C function pointer
        !
        use otr_qn_c_interface, only: transport_c_type

        type(c_funptr), intent(in) :: transport_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(transport_c_type), pointer :: transport_funptr
        real(c_rp), allocatable :: geodesic(:), tangent_vector(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=transport_c_funptr, fptr=transport_funptr)

        ! allocate arrays
        allocate(geodesic(n_param), tangent_vector(n_param))

        ! initialize new reference
        geodesic = 2.0_c_rp
        tangent_vector = 3.0_c_rp

        ! call transport function pointer
        error = transport_funptr(geodesic, tangent_vector)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check tangent vector
        if (any((tangent_vector - 6.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Tangent vector returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(geodesic, tangent_vector)

    end function test_transport_c_funptr

    function test_init_hess_funptr(init_hess_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided initial Hessian function pointer
        !
        use otr_qn, only: init_hess_type

        procedure(init_hess_type), intent(in), pointer :: init_hess_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: vector(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(vector(n_param))

        ! initialize input arrays
        vector = 2.0_rp

        ! call initial Hessian function pointer
        call init_hess_funptr(vector, error)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check vector
        if (any((vector - 4.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Vector returned"//message &
                //" wrong."
        end if

        ! deallocate arrays
        deallocate(vector)

    end function test_init_hess_funptr

    function test_init_hess_c_funptr(init_hess_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided initial Hessian C function pointer
        !
        use otr_qn_c_interface, only: init_hess_c_type

        type(c_funptr), intent(in) :: init_hess_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(init_hess_c_type), pointer :: init_hess_funptr
        real(c_rp), allocatable :: vector(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=init_hess_c_funptr, fptr=init_hess_funptr)

        ! allocate arrays
        allocate(vector(n_param))

        ! initialize new reference
        vector = 2.0_c_rp

        ! call initial Hessian function pointer
        error = init_hess_funptr(vector)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
        end if

        ! check vector
        if (any((vector - 4.0_c_rp) > tol_c)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Vector returned"//message &
                //" wrong."
        end if

        ! deallocate arrays
        deallocate(vector)

    end function test_init_hess_c_funptr

    subroutine get_reference_qn_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the quasi-Newton reference values for tests
        !
        type(ref_qn_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_qn_settings

    end subroutine get_reference_qn_values

    subroutine assign_ref_to_qn(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set quasi-Newton 
        ! settings to reference values
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(out) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%logger => null()

        ! set reference values
        lhs%verbose = rhs%verbose
        lhs%max_points = rhs%max_points
        lhs%hess_update_scheme = rhs%hess_update_scheme

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_qn

    subroutine assign_ref_to_qn_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C quasi-Newton 
        ! settings to reference values
        !
        use otr_qn_c_interface, only: qn_settings_type_c, assignment(=)
        use otr_qn, only: qn_settings_type

        type(qn_settings_type_c), intent(out) :: lhs_c
        type(ref_qn_settings_type), intent(in) :: rhs

        type(qn_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_qn_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        use c_interface, only: character_to_c

        type(ref_qn_settings_type_c), intent(out) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs

        lhs%verbose = int(rhs%verbose, kind=c_ip)
        lhs%max_points = int(rhs%max_points, kind=c_ip)
        lhs%hess_update_scheme = character_to_c(rhs%hess_update_scheme)

    end subroutine assign_ref_to_ref_c

    logical function equal_qn_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare quasi-Newton 
        ! settings to reference values
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(in) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs

        equal_qn_to_ref = lhs%verbose == rhs%verbose .and. &
                          lhs%max_points == rhs%max_points .and. &
                          lhs%hess_update_scheme == rhs%hess_update_scheme

    end function equal_qn_to_ref

    logical function not_equal_qn_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare 
        ! quasi-Newton settings to reference values
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(in) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs

        not_equal_qn_to_ref = .not. (lhs == rhs)

    end function not_equal_qn_to_ref

    logical function equal_qn_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare quasi-Newton 
        ! settings to reference values
        !
        use otr_qn_c_interface, only: qn_settings_type_c, assignment(=)
        use otr_qn, only: qn_settings_type

        type(qn_settings_type_c), intent(in) :: lhs_c
        type(ref_qn_settings_type), intent(in) :: rhs

        type(qn_settings_type) :: lhs

        lhs = lhs_c
        equal_qn_c_to_ref = lhs == rhs

    end function equal_qn_c_to_ref

    logical function not_equal_qn_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare 
        ! quasi-Newton settings to reference values
        !
        use otr_qn_c_interface, only: qn_settings_type_c

        type(qn_settings_type_c), intent(in) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs
        
        not_equal_qn_c_to_ref = .not. (lhs == rhs)

    end function not_equal_qn_c_to_ref

    logical function equal_qn(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare quasi-Newton 
        ! settings to different quasi-Newton settings
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(in) :: lhs, rhs
        
        equal_qn = lhs%verbose == rhs%verbose .and. &
                   lhs%max_points == rhs%max_points .and. &
                   lhs%hess_update_scheme == rhs%hess_update_scheme

    end function equal_qn

    logical function not_equal_qn(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare 
        ! quasi-Newton settings to different quasi-Newton settings
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(in) :: lhs, rhs
        
        not_equal_qn = .not. (lhs == rhs)

    end function not_equal_qn

    logical function equal_qn_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare quasi-Newton 
        ! settings to different quasi-Newton settings
        !
        use otr_qn_c_interface, only: qn_settings_type_c, assignment(=)
        use otr_qn, only: qn_settings_type

        type(qn_settings_type_c), intent(in) :: lhs_c
        type(qn_settings_type), intent(in) :: rhs
        
        type(qn_settings_type) :: lhs

        lhs = lhs_c
        equal_qn_c = lhs == rhs

    end function equal_qn_c

    logical function not_equal_qn_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare 
        ! quasi-Newton settings to different quasi-Newton settings
        !
        use otr_qn_c_interface, only: qn_settings_type_c
        use otr_qn, only: qn_settings_type

        type(qn_settings_type_c), intent(in) :: lhs_c
        type(qn_settings_type), intent(in) :: rhs
        
        not_equal_qn_c = .not. (lhs_c == rhs)

    end function not_equal_qn_c

end module otr_qn_test_reference
