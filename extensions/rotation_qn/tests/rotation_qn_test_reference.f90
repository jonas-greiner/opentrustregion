! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_rotation_qn_test_reference

    use opentrustregion, only: rp, ip, kw_len, stderr
    use c_interface, only: c_rp, c_ip
    use otr_rotation_qn, only: rotation_qn_settings_type
    use test_reference, only: n_param, tol, tol_c
    use, intrinsic :: iso_c_binding, only: &
        c_size_t, c_sizeof, c_bool, c_char, c_ptr, c_funptr, c_f_procpointer, &
        c_associated, c_intptr_t, c_null_char, c_null_funptr, c_null_ptr

    implicit none

    ! reference rotation QN settings, every field but the logging function and the host
    ! context differing from its default and every integer from the other
    type(rotation_qn_settings_type), parameter :: ref_rotation_qn_settings = &
        rotation_qn_settings_type(logger=null(), verbose=3, max_points=100, &
                                  qn_method="bfgs")

    interface operator(==)
        module procedure equal_rotation_qn
        module procedure equal_rotation_qn_c
    end interface

    interface operator(/=)
        module procedure not_equal_rotation_qn
        module procedure not_equal_rotation_qn_c
    end interface

contains

    function check_transport_funptr(transport_funptr, test_name, message, context) &
        result(test_passed)
        !
        ! this function tests a provided transport function pointer
        !
        use otr_rotation_qn, only: transport_type

        procedure(transport_type), intent(in), pointer :: transport_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), allocatable :: geodesic(:), tangent_vector(:), expected(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(transport_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Transport function "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(geodesic(n_param), tangent_vector(n_param))

        ! generate random input arrays, which the mock transport multiplies
        call random_number(geodesic)
        call random_number(tangent_vector)
        expected = geodesic * tangent_vector

        ! call transport function pointer
        call transport_funptr(geodesic, tangent_vector, error, context)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check tangent vector
        if (any(abs(tangent_vector - expected) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Tangent vector returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(geodesic, tangent_vector, expected)

    end function check_transport_funptr

    function check_transport_c_funptr(transport_c_funptr, test_name, message, &
                                      context_c) result(test_passed)
        !
        ! this function tests a provided transport C function pointer
        !
        use otr_rotation_qn_c_interface, only: transport_c_type

        type(c_funptr), intent(in) :: transport_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(transport_c_type), pointer :: transport_funptr
        real(c_rp), allocatable :: geodesic(:), tangent_vector(:), expected(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(transport_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Transport function "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=transport_c_funptr, fptr=transport_funptr)

        ! allocate arrays
        allocate(geodesic(n_param), tangent_vector(n_param))

        ! generate random input arrays, which the mock transport multiplies
        call random_number(geodesic)
        call random_number(tangent_vector)
        expected = geodesic * tangent_vector

        ! call transport function pointer
        error = transport_funptr(geodesic, tangent_vector, context_c)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check tangent vector
        if (any(abs(tangent_vector - expected) > tol_c)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Tangent vector returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(geodesic, tangent_vector, expected)

    end function check_transport_c_funptr

    function check_init_hess_funptr(init_hess_funptr, test_name, message, context) &
        result(test_passed)
        !
        ! this function tests a provided initial Hessian function pointer
        !
        use otr_rotation_qn, only: init_hess_type

        procedure(init_hess_type), intent(in), pointer :: init_hess_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), allocatable :: vector(:), expected(:)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(init_hess_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Initial Hessian "// &
                "function provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(vector(n_param))

        ! generate a random input array, which the mock initial Hessian doubles
        call random_number(vector)
        expected = 2.0_rp * vector

        ! call initial Hessian function pointer
        call init_hess_funptr(vector, error, context)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check vector
        if (any(abs(vector - expected) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Vector returned"//message// &
                " wrong."
        end if

        ! deallocate arrays
        deallocate(vector, expected)

    end function check_init_hess_funptr

    function check_init_hess_c_funptr(init_hess_c_funptr, test_name, message, &
                                      context_c) result(test_passed)
        !
        ! this function tests a provided initial Hessian C function pointer
        !
        use otr_rotation_qn_c_interface, only: init_hess_c_type

        type(c_funptr), intent(in) :: init_hess_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(init_hess_c_type), pointer :: init_hess_funptr
        real(c_rp), allocatable :: vector(:), expected(:)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(init_hess_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Initial Hessian "// &
                "function provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=init_hess_c_funptr, fptr=init_hess_funptr)

        ! allocate arrays
        allocate(vector(n_param))

        ! generate a random input array, which the mock initial Hessian doubles
        call random_number(vector)
        expected = 2.0_c_rp * vector

        ! call initial Hessian function pointer
        error = init_hess_funptr(vector, context_c)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check vector
        if (any(abs(vector - expected) > tol_c)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Vector returned"//message// &
                " wrong."
        end if

        ! deallocate arrays
        deallocate(vector, expected)

    end function check_init_hess_c_funptr

    subroutine get_reference_rotation_qn_values(values_out) bind(C)
        !
        ! this subroutine exports the reference rotation QN settings as C settings
        ! whose fields are set one by one without the conversion routines under test,
        ! the logging function and the host context are set to distinct addresses
        !
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c
        use test_reference, only: ref_character_to_c

        type(rotation_qn_settings_type_c), intent(out) :: values_out

        ! set logging function and host context to distinct addresses
        values_out%logger = transfer(1_c_intptr_t, c_null_funptr)
        values_out%context = transfer(2_c_intptr_t, c_null_ptr)

        ! set logical, integers and keyword
        values_out%initialized = .true.
        values_out%verbose = int(ref_rotation_qn_settings%verbose, kind=c_ip)
        values_out%max_points = int(ref_rotation_qn_settings%max_points, kind=c_ip)
        values_out%qn_method = ref_character_to_c(ref_rotation_qn_settings%qn_method)

    end subroutine get_reference_rotation_qn_values

    subroutine reference_rotation_qn_field(name_c, field_value, keyword_c) bind(C)
        !
        ! this subroutine returns the value of a field of the C rotation QN settings
        ! exported by get_reference_rotation_qn_values by its name, a numeric, logical
        ! (as 1 or 0) or pointer (as its address) field in field_value and a keyword
        ! field as a null-terminated character array in keyword_c, so that the Python
        ! tests need no values of their own
        !
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c
        use test_reference, only: ref_character_from_c

        character(kind=c_char), intent(in) :: name_c(*)
        real(c_rp), intent(out) :: field_value
        character(kind=c_char), intent(out) :: keyword_c(kw_len + 1)

        type(rotation_qn_settings_type_c) :: values

        ! export reference values and select field by name
        field_value = 0.0_c_rp
        keyword_c = c_null_char
        call get_reference_rotation_qn_values(values)
        select case (ref_character_from_c(name_c))
        case ("logger")
            field_value = real(transfer(values%logger, 0_c_intptr_t), kind=c_rp)
        case ("context")
            field_value = real(transfer(values%context, 0_c_intptr_t), kind=c_rp)
        case ("initialized")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%initialized)
        case ("verbose")
            field_value = values%verbose
        case ("max_points")
            field_value = values%max_points
        case ("qn_method")
            keyword_c = values%qn_method
        end select

    end subroutine reference_rotation_qn_field

    logical(c_bool) function is_default_rotation_qn_settings(settings_c) bind(C)
        !
        ! this function checks whether C rotation QN settings hold the default values
        ! without logging function and host context, so that a C test can compare
        ! against them without values of its own
        !
        use otr_rotation_qn, only: default_rotation_qn_settings
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c

        type(rotation_qn_settings_type_c), intent(in) :: settings_c

        is_default_rotation_qn_settings = &
            settings_c == default_rotation_qn_settings .and. &
            .not. c_associated(settings_c%logger) .and. &
            .not. c_associated(settings_c%context)

    end function is_default_rotation_qn_settings

    integer(c_size_t) function rotation_qn_settings_size() bind(C)
        !
        ! this function returns the size of the C rotation QN settings, so that the
        ! layout tests detect a field missing from the C or Python settings
        !
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c

        type(rotation_qn_settings_type_c) :: settings_c

        rotation_qn_settings_size = c_sizeof(settings_c)

    end function rotation_qn_settings_size

    logical function equal_rotation_qn(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare rotation QN
        ! settings to different rotation QN settings
        !
        type(rotation_qn_settings_type), intent(in) :: lhs, rhs

        equal_rotation_qn = lhs%verbose == rhs%verbose .and. &
                            lhs%max_points == rhs%max_points .and. &
                            lhs%qn_method == rhs%qn_method

    end function equal_rotation_qn

    logical function not_equal_rotation_qn(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare rotation
        ! QN settings to different rotation QN settings
        !
        type(rotation_qn_settings_type), intent(in) :: lhs, rhs

        not_equal_rotation_qn = .not. (lhs == rhs)

    end function not_equal_rotation_qn

    logical function equal_rotation_qn_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare C rotation QN
        ! settings to Fortran rotation QN settings field by field, independently of the
        ! conversion routines under test
        !
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c
        use test_reference, only: ref_character_from_c

        type(rotation_qn_settings_type_c), intent(in) :: lhs_c
        type(rotation_qn_settings_type), intent(in) :: rhs

        equal_rotation_qn_c = &
            lhs_c%initialized .and. lhs_c%verbose == rhs%verbose .and. &
            lhs_c%max_points == rhs%max_points .and. &
            ref_character_from_c(lhs_c%qn_method) == trim(rhs%qn_method)

    end function equal_rotation_qn_c

    logical function not_equal_rotation_qn_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare C rotation
        ! QN settings to Fortran rotation QN settings
        !
        use otr_rotation_qn_c_interface, only: rotation_qn_settings_type_c

        type(rotation_qn_settings_type_c), intent(in) :: lhs_c
        type(rotation_qn_settings_type), intent(in) :: rhs

        not_equal_rotation_qn_c = .not. (lhs_c == rhs)

    end function not_equal_rotation_qn_c

end module otr_rotation_qn_test_reference
