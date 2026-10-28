! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_test_reference

    use opentrustregion_core, only: rp, ip, kw_len, stderr
    use c_interface, only: c_rp, c_ip
    use otr_s_gek, only: s_gek_settings_type
    use test_reference, only: n_param, tol, tol_c
    use, intrinsic :: iso_c_binding, only: &
        c_size_t, c_sizeof, c_bool, c_char, c_ptr, c_funptr, c_f_procpointer, &
        c_associated, c_intptr_t, c_null_char, c_null_funptr, c_null_ptr

    implicit none

    ! number of stored points the change of reference checkers pass
    integer(ip), parameter :: n_points = 4_ip

    ! reference S-GEK settings, every field but the logging function and the host
    ! context differing from its default and every integer from the other, the use of
    ! the subspace false unlike the initialization flag of the C settings, so that the
    ! two logicals cannot be swapped unnoticed
    type(s_gek_settings_type), parameter :: ref_s_gek_settings = s_gek_settings_type( &
        logger=null(), use_subspace=.false., verbose=3, max_points=10)

    interface operator(==)
        module procedure equal_s_gek
        module procedure equal_s_gek_c
    end interface

    interface operator(/=)
        module procedure not_equal_s_gek
        module procedure not_equal_s_gek_c
    end interface

contains

    function check_change_reference_funptr(change_reference_funptr, test_name, &
                                           message, context) result(test_passed)
        !
        ! this function tests a provided change reference function pointer
        !
        use otr_s_gek, only: change_reference_type

        procedure(change_reference_type), intent(in), pointer :: change_reference_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), allocatable :: new_ref(:), kappa_list(:, :), local_grad_list(:, :), &
                                 grad_list(:, :), expected_lists(:, :, :)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(change_reference_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Change reference "// &
                "function provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(new_ref(n_param), kappa_list(n_param, n_points), &
                 local_grad_list(n_param, n_points), grad_list(n_param, n_points))

        ! generate a random new reference and random lists, which the mock change of
        ! reference changes to the new reference in a flat parameter space, subtracting
        ! it from every stored point and taking the gradients over from the local
        ! gradients, and whose local gradients it triples, so that a transposed or
        ! shifted view of the lists or the new reference is detected
        call random_number(new_ref)
        call random_number(kappa_list)
        call random_number(local_grad_list)
        call random_number(grad_list)
        expected_lists = reshape([kappa_list - spread(new_ref, 2, n_points), &
                                  3.0_rp * local_grad_list, local_grad_list], &
                                 [n_param, n_points, 3_ip])

        ! call change reference function pointer
        call change_reference_funptr(new_ref, n_points, kappa_list, local_grad_list, &
                                     grad_list, error, context)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check kappa list
        if (any(abs(kappa_list - expected_lists(:, :, 1)) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Parameter list returned"// &
                message//" wrong."
        end if

        ! check local gradient list
        if (any(abs(local_grad_list - expected_lists(:, :, 2)) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Local gradient list returned"//message//" wrong."
        end if

        ! check gradient list
        if (any(abs(grad_list - expected_lists(:, :, 3)) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Gradient list returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(new_ref, kappa_list, local_grad_list, grad_list, expected_lists)

    end function check_change_reference_funptr

    function check_change_reference_c_funptr(change_reference_c_funptr, test_name, &
                                             message, context_c) result(test_passed)
        !
        ! this function tests a provided change reference C function pointer
        !
        use otr_s_gek_c_interface, only: change_reference_c_type

        type(c_funptr), intent(in) :: change_reference_c_funptr
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(change_reference_c_type), pointer :: change_reference_funptr
        real(c_rp), allocatable :: new_ref(:), kappa_list(:, :), &
                                   local_grad_list(:, :), grad_list(:, :), &
                                   expected_lists(:, :, :)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(change_reference_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Change reference "// &
                "function provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=change_reference_c_funptr, &
                             fptr=change_reference_funptr)

        ! allocate arrays
        allocate(new_ref(n_param), kappa_list(n_param, n_points), &
                 local_grad_list(n_param, n_points), grad_list(n_param, n_points))

        ! generate a random new reference and random lists, which the mock change of
        ! reference changes to the new reference in a flat parameter space, subtracting
        ! it from every stored point and taking the gradients over from the local
        ! gradients, and whose local gradients it triples, so that a transposed or
        ! shifted view of the lists or the new reference is detected
        call random_number(new_ref)
        call random_number(kappa_list)
        call random_number(local_grad_list)
        call random_number(grad_list)
        expected_lists = reshape([kappa_list - spread(new_ref, 2, n_points), &
                                  3.0_c_rp * local_grad_list, local_grad_list], &
                                 [n_param, n_points, 3_ip])

        ! call change reference function pointer
        error = change_reference_funptr(new_ref, int(n_points, kind=c_ip), kappa_list, &
                                        local_grad_list, grad_list, context_c)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check kappa list
        if (any(abs(kappa_list - expected_lists(:, :, 1)) > tol_c)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Parameter list returned"// &
                message//" wrong."
        end if

        ! check local gradient list
        if (any(abs(local_grad_list - expected_lists(:, :, 2)) > tol_c)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Local gradient list returned"//message//" wrong."
        end if

        ! check gradient list
        if (any(abs(grad_list - expected_lists(:, :, 3)) > tol_c)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Gradient list returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(new_ref, kappa_list, local_grad_list, grad_list, expected_lists)

    end function check_change_reference_c_funptr

    subroutine get_reference_s_gek_values(values_out) bind(C)
        !
        ! this subroutine exports the reference S-GEK settings as C settings whose
        ! fields are set one by one without the conversion routines under test, the
        ! logging function and the host context are set to distinct addresses
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c

        type(s_gek_settings_type_c), intent(out) :: values_out

        ! set logging function and host context to distinct addresses
        values_out%logger = transfer(1_c_intptr_t, c_null_funptr)
        values_out%context = transfer(2_c_intptr_t, c_null_ptr)

        ! set logicals and integers
        values_out%initialized = .true.
        values_out%use_subspace = ref_s_gek_settings%use_subspace
        values_out%verbose = int(ref_s_gek_settings%verbose, kind=c_ip)
        values_out%max_points = int(ref_s_gek_settings%max_points, kind=c_ip)

    end subroutine get_reference_s_gek_values

    subroutine reference_s_gek_field(name_c, field_value, keyword_c) bind(C)
        !
        ! this subroutine returns the value of a field of the C S-GEK settings exported
        ! by get_reference_s_gek_values by its name, a numeric, logical (as 1 or 0) or
        ! pointer (as its address) field in field_value and a keyword field as a
        ! null-terminated character array in keyword_c, so that the Python tests need no
        ! values of their own
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c
        use test_reference, only: ref_character_from_c

        character(kind=c_char), intent(in) :: name_c(*)
        real(c_rp), intent(out) :: field_value
        character(kind=c_char), intent(out) :: keyword_c(kw_len + 1)

        type(s_gek_settings_type_c) :: values

        ! export reference values and select field by name
        field_value = 0.0_c_rp
        keyword_c = c_null_char
        call get_reference_s_gek_values(values)
        select case (ref_character_from_c(name_c))
        case ("logger")
            field_value = real(transfer(values%logger, 0_c_intptr_t), kind=c_rp)
        case ("context")
            field_value = real(transfer(values%context, 0_c_intptr_t), kind=c_rp)
        case ("initialized")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%initialized)
        case ("use_subspace")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%use_subspace)
        case ("verbose")
            field_value = values%verbose
        case ("max_points")
            field_value = values%max_points
        end select

    end subroutine reference_s_gek_field

    logical(c_bool) function is_default_s_gek_settings(settings_c) bind(C)
        !
        ! this function checks whether C S-GEK settings hold the default values without
        ! logging function and host context, so that a C test can compare against them
        ! without values of its own
        !
        use otr_s_gek, only: default_s_gek_settings
        use otr_s_gek_c_interface, only: s_gek_settings_type_c

        type(s_gek_settings_type_c), intent(in) :: settings_c

        is_default_s_gek_settings = settings_c == default_s_gek_settings .and. &
                                    .not. c_associated(settings_c%logger) .and. &
                                    .not. c_associated(settings_c%context)

    end function is_default_s_gek_settings

    integer(c_size_t) function s_gek_settings_size() bind(C)
        !
        ! this function returns the size of the C S-GEK settings, so that the layout
        ! tests detect a field missing from the C or Python settings
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c

        type(s_gek_settings_type_c) :: settings_c

        s_gek_settings_size = c_sizeof(settings_c)

    end function s_gek_settings_size

    logical function equal_s_gek(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to
        ! different S-GEK settings
        !
        type(s_gek_settings_type), intent(in) :: lhs, rhs

        equal_s_gek = (lhs%use_subspace .eqv. rhs%use_subspace) .and. &
                      lhs%verbose == rhs%verbose .and. lhs%max_points == rhs%max_points

    end function equal_s_gek

    logical function not_equal_s_gek(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK
        ! settings to different S-GEK settings
        !
        type(s_gek_settings_type), intent(in) :: lhs, rhs

        not_equal_s_gek = .not. (lhs == rhs)

    end function not_equal_s_gek

    logical function equal_s_gek_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare C S-GEK settings
        ! to Fortran S-GEK settings field by field, independently of the conversion
        ! routines under test
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(s_gek_settings_type), intent(in) :: rhs

        equal_s_gek_c = &
            lhs_c%initialized .and. &
            (lhs_c%use_subspace .eqv. logical(rhs%use_subspace, kind=c_bool)) .and. &
            lhs_c%verbose == rhs%verbose .and. lhs_c%max_points == rhs%max_points

    end function equal_s_gek_c

    logical function not_equal_s_gek_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare C S-GEK
        ! settings to Fortran S-GEK settings
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(s_gek_settings_type), intent(in) :: rhs

        not_equal_s_gek_c = .not. (lhs_c == rhs)

    end function not_equal_s_gek_c

end module otr_s_gek_test_reference
