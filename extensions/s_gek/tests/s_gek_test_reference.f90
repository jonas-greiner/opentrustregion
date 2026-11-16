! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_test_reference

    use opentrustregion_core, only: rp, ip, kw_len, stderr
    use c_interface, only: c_rp, c_ip
    use otr_s_gek, only: s_gek_settings_type, s_gek_settings_uninitialized_warning_msg
    use test_reference, only: n_param, is_close
    use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_ptr, c_funptr, &
                                           c_f_procpointer, c_associated, c_intptr_t, &
                                           c_null_char, c_null_funptr, c_null_ptr

    implicit none

    ! number of points the chart transformation checkers pass, which differs from the
    ! number of parameters so that the two dimensions of the lists cannot be mixed up
    integer(ip), parameter :: n_points = 4_ip

    ! reference S-GEK settings, every field but the callback functions and the host
    ! context differing from its default and every keyword from the others, the update
    ! of the prior and the acceptance by the case of their common value "none", which
    ! the S-GEK solver does not distinguish
    type(s_gek_settings_type), parameter :: ref_s_gek_settings = &
        s_gek_settings_type(logger=null(), obj_func_grad=null(), init_hess=null(), &
                            verbose=3, inner_solver="rs-rfo", prior_update="none", &
                            acceptance="NONE", excess_dispersion=.false.)

    ! warning the C entry point prints for S-GEK settings that were not initialized, as
    ! a null-terminated string for the C integration tests
    integer(ip), parameter :: settings_warning_len = &
        len(s_gek_settings_uninitialized_warning_msg) + 1
    character(kind=c_char), protected, &
        bind(C, name="test_s_gek_settings_uninitialized_warning") :: &
        settings_uninitialized_warning_c(settings_warning_len) = &
        transfer(s_gek_settings_uninitialized_warning_msg//c_null_char, "a", &
                 settings_warning_len)

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
        ! this function tests a provided chart transformation function pointer
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
            write(stderr, *) "test_"//test_name//" failed: Chart transformation "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(new_ref(n_param), kappa_list(n_param, n_points), &
                 local_grad_list(n_param, n_points), grad_list(n_param, n_points))

        ! generate a random new reference and random lists, which the mock chart
        ! transformation changes to the new reference in a flat parameter space,
        ! subtracting it from every point and taking the gradients over from the local
        ! gradients, and whose local gradients it triples, so that a transposed or
        ! shifted view of the lists or the new reference is detected
        call random_number(new_ref)
        call random_number(kappa_list)
        call random_number(local_grad_list)
        call random_number(grad_list)
        expected_lists = reshape([kappa_list - spread(new_ref, 2, n_points), &
                                  3.0_rp * local_grad_list, local_grad_list], &
                                 [n_param, n_points, 3_ip])

        ! call chart transformation function pointer
        call change_reference_funptr(new_ref, kappa_list, local_grad_list, grad_list, &
                                     error, context)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check positions
        if (.not. all(is_close(kappa_list, expected_lists(:, :, 1)))) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Positions returned"// &
                message//" wrong."
        end if

        ! check local gradients
        if (.not. all(is_close(local_grad_list, expected_lists(:, :, 2)))) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Local gradients returned"// &
                message//" wrong."
        end if

        ! check gradients
        if (.not. all(is_close(grad_list, expected_lists(:, :, 3)))) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Gradients returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(new_ref, kappa_list, local_grad_list, grad_list, expected_lists)

    end function check_change_reference_funptr

    function check_change_reference_c_funptr(change_reference_c_funptr, test_name, &
                                             message, context_c) result(test_passed)
        !
        ! this function tests a provided chart transformation C function pointer
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
            write(stderr, *) "test_"//test_name//" failed: Chart transformation "// &
                "provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=change_reference_c_funptr, &
                             fptr=change_reference_funptr)

        ! allocate arrays
        allocate(new_ref(n_param), kappa_list(n_param, n_points), &
                 local_grad_list(n_param, n_points), grad_list(n_param, n_points))

        ! generate a random new reference and random lists, which the mock chart
        ! transformation changes to the new reference in a flat parameter space,
        ! subtracting it from every point and taking the gradients over from the local
        ! gradients, and whose local gradients it triples, so that a transposed or
        ! shifted view of the lists or the new reference is detected
        call random_number(new_ref)
        call random_number(kappa_list)
        call random_number(local_grad_list)
        call random_number(grad_list)
        expected_lists = reshape([kappa_list - spread(new_ref, 2, n_points), &
                                  3.0_c_rp * local_grad_list, local_grad_list], &
                                 [n_param, n_points, 3_ip])

        ! call chart transformation function pointer
        error = change_reference_funptr(new_ref, int(n_points, kind=c_ip), kappa_list, &
                                        local_grad_list, grad_list, context_c)

        ! check for error
        if (error /= 0) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message//"."
        end if

        ! check positions
        if (.not. all(is_close(kappa_list, expected_lists(:, :, 1)))) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Positions returned"// &
                message//" wrong."
        end if

        ! check local gradients
        if (.not. all(is_close(local_grad_list, expected_lists(:, :, 2)))) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Local gradients returned"// &
                message//" wrong."
        end if

        ! check gradients
        if (.not. all(is_close(grad_list, expected_lists(:, :, 3)))) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Gradients returned"// &
                message//" wrong."
        end if

        ! deallocate arrays
        deallocate(new_ref, kappa_list, local_grad_list, grad_list, expected_lists)

    end function check_change_reference_c_funptr

    subroutine get_reference_s_gek_values(values_out) bind(C)
        !
        ! this subroutine exports the reference S-GEK settings as C settings whose
        ! fields are set one by one without the conversion routines under test, the
        ! callback functions and the host context are set to distinct addresses
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c
        use test_reference, only: ref_character_to_c

        type(s_gek_settings_type_c), intent(out) :: values_out

        ! set callback functions and host context to distinct addresses
        values_out%obj_func_grad = transfer(1_c_intptr_t, c_null_funptr)
        values_out%init_hess = transfer(2_c_intptr_t, c_null_funptr)
        values_out%logger = transfer(3_c_intptr_t, c_null_funptr)
        values_out%context = transfer(4_c_intptr_t, c_null_ptr)

        ! set logicals, integer and keywords
        values_out%initialized = .true.
        values_out%excess_dispersion = ref_s_gek_settings%excess_dispersion
        values_out%verbose = int(ref_s_gek_settings%verbose, kind=c_ip)
        values_out%inner_solver = ref_character_to_c(ref_s_gek_settings%inner_solver)
        values_out%prior_update = ref_character_to_c(ref_s_gek_settings%prior_update)
        values_out%acceptance = ref_character_to_c(ref_s_gek_settings%acceptance)

    end subroutine get_reference_s_gek_values

    subroutine reference_s_gek_field(name_c, field_value, keyword_c) bind(C)
        !
        ! this subroutine returns the value of a field of the C S-GEK settings exported
        ! by get_reference_s_gek_values by its name, a numeric, logical (as 1 or 0) or
        ! pointer (as its address) field in field_value and a keyword field as a
        ! null-terminated character array in keyword_c, so that the Python tests need
        ! no values of their own
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
        case ("obj_func_grad")
            field_value = real(transfer(values%obj_func_grad, 0_c_intptr_t), kind=c_rp)
        case ("init_hess")
            field_value = real(transfer(values%init_hess, 0_c_intptr_t), kind=c_rp)
        case ("logger")
            field_value = real(transfer(values%logger, 0_c_intptr_t), kind=c_rp)
        case ("context")
            field_value = real(transfer(values%context, 0_c_intptr_t), kind=c_rp)
        case ("initialized")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%initialized)
        case ("excess_dispersion")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%excess_dispersion)
        case ("verbose")
            field_value = values%verbose
        case ("inner_solver")
            keyword_c = values%inner_solver
        case ("prior_update")
            keyword_c = values%prior_update
        case ("acceptance")
            keyword_c = values%acceptance
        end select

    end subroutine reference_s_gek_field

    logical(c_bool) function is_default_s_gek_settings(settings_c) bind(C)
        !
        ! this function checks whether C S-GEK settings hold the default values without
        ! callback functions and host context, so that a C test can compare against
        ! them without values of its own
        !
        use otr_s_gek, only: default_s_gek_settings
        use otr_s_gek_c_interface, only: s_gek_settings_type_c

        type(s_gek_settings_type_c), intent(in) :: settings_c

        is_default_s_gek_settings = settings_c == default_s_gek_settings .and. &
                                    .not. c_associated(settings_c%obj_func_grad) .and. &
                                    .not. c_associated(settings_c%init_hess) .and. &
                                    .not. c_associated(settings_c%logger) .and. &
                                    .not. c_associated(settings_c%context)

    end function is_default_s_gek_settings

    logical function equal_s_gek(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to
        ! different S-GEK settings
        !
        type(s_gek_settings_type), intent(in) :: lhs, rhs

        equal_s_gek = lhs%verbose == rhs%verbose .and. &
                      lhs%inner_solver == rhs%inner_solver .and. &
                      lhs%prior_update == rhs%prior_update .and. &
                      lhs%acceptance == rhs%acceptance .and. &
                      (lhs%excess_dispersion .eqv. rhs%excess_dispersion)

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
        use test_reference, only: ref_character_from_c

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(s_gek_settings_type), intent(in) :: rhs

        equal_s_gek_c = &
            lhs_c%initialized .and. lhs_c%verbose == rhs%verbose .and. &
            ref_character_from_c(lhs_c%inner_solver) == trim(rhs%inner_solver) .and. &
            ref_character_from_c(lhs_c%prior_update) == trim(rhs%prior_update) .and. &
            ref_character_from_c(lhs_c%acceptance) == trim(rhs%acceptance) .and. &
            (logical(lhs_c%excess_dispersion) .eqv. rhs%excess_dispersion)

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
