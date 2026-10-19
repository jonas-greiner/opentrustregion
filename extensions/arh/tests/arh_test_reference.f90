! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_density_qn_test_reference

    use opentrustregion, only: ip, rp, kw_len, stderr
    use c_interface, only: c_ip, c_rp
    use otr_density_qn, only: density_qn_settings_type, &
                              density_qn_settings_uninitialized_warning_msg
    use otr_common_test_reference, only: ref_orbital_settings
    use, intrinsic :: iso_c_binding, only: &
        c_size_t, c_sizeof, c_bool, c_char, c_ptr, c_funptr, c_f_procpointer, &
        c_associated, c_intptr_t, c_null_char, c_null_funptr, c_null_ptr

    implicit none

    ! reference density QN settings, those shared by the orbital bases and a QN method
    ! other than the default
    type(density_qn_settings_type), parameter :: ref_density_qn_settings = &
        density_qn_settings_type(orbital_settings_type=ref_orbital_settings, &
                                 qn_method="symm_arh")

    ! warning the C factory prints for density QN settings that were not initialized,
    ! as a null-terminated string for the C integration tests
    integer(ip), parameter :: settings_warning_len = &
        len(density_qn_settings_uninitialized_warning_msg) + 1
    character(kind=c_char), protected, &
        bind(C, name="test_density_qn_settings_uninitialized_warning") :: &
        settings_uninitialized_warning_c(settings_warning_len) = &
        transfer(density_qn_settings_uninitialized_warning_msg//c_null_char, "a", &
                 settings_warning_len)

    ! multiples of the density matrix the mock density matrix evaluating functions of
    ! both shells return for each optional output, in the order of
    ! density_qn_evaluate_dm_outputs
    real(c_rp), protected, bind(C, name="test_density_qn_evaluate_dm_factors") :: &
        density_qn_evaluate_dm_factors(4) = [2.0_c_rp, 3.0_c_rp, 4.0_c_rp, 5.0_c_rp]

    ! optional outputs of the density matrix evaluating functions of both shells in the
    ! order of their argument lists
    character(len=24), parameter :: density_qn_evaluate_dm_outputs(4) = &
        [character(len=24) :: "Fock matrix", "Coulomb potential", &
         "exact-exchange potential", "non-linear potential"]

    interface operator(==)
        module procedure equal_density_qn
        module procedure equal_density_qn_c
    end interface

    interface operator(/=)
        module procedure not_equal_density_qn
        module procedure not_equal_density_qn_c
    end interface

contains

    function capitalized(text) result(capital)
        !
        ! this function returns a text with its first letter capitalized, so that the
        ! name of an output can start a failure message
        !
        character(len=*), intent(in) :: text
        character(len=len(text)) :: capital

        capital = text
        if (lge(text(1:1), "a") .and. lle(text(1:1), "z")) &
            capital(1:1) = achar(iachar(text(1:1)) - iachar("a") + iachar("A"))

    end function capitalized

    function check_density_qn_evaluate_dm_cs_funptr( &
        evaluate_dm_funptr, test_name, message, context) result(test_passed)
        !
        ! this function tests a provided density matrix evaluating function pointer
        ! with Coulomb, exact-exchange and non-linear potential contributions for the
        ! closed-shell case for every combination of requested outputs
        !
        use otr_density_qn, only: density_qn_evaluate_dm_cs_type
        use otr_common_test_reference, only: n_ao, request_label
        use test_reference, only: tol

        procedure(density_qn_evaluate_dm_cs_type), intent(in), pointer :: &
            evaluate_dm_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), target :: dm_ao(n_ao, n_ao), &
                            outputs(n_ao, n_ao, size(density_qn_evaluate_dm_outputs))
        real(rp), pointer, contiguous :: fock(:, :), v_coulomb(:, :), &
                                         v_exchange(:, :), v_nonlinear(:, :)
        real(rp) :: energy
        character(len=:), allocatable :: requested
        integer(ip) :: request, i, error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(evaluate_dm_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function with Coulomb, exact-exchange and non-linear "// &
                "potential contributions for closed-shell case provided"//message// &
                " not associated with value."
            return
        end if

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating subroutine for every combination of requested
        ! outputs
        do request = 0, 15
            requested = request_label(density_qn_evaluate_dm_outputs, request)
            nullify(fock, v_coulomb, v_exchange, v_nonlinear)
            if (btest(request, 0)) fock => outputs(:, :, 1)
            if (btest(request, 1)) v_coulomb => outputs(:, :, 2)
            if (btest(request, 2)) v_exchange => outputs(:, :, 3)
            if (btest(request, 3)) v_nonlinear => outputs(:, :, 4)
            energy = 0.0_rp
            outputs = 0.0_rp
            call evaluate_dm_funptr(dm_ao, energy, fock, v_coulomb, v_exchange, &
                                    v_nonlinear, error, context)

            ! check for error
            if (error /= 0) then
                write(stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for closed-shell case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol) then
                write(stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for closed-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check requested outputs
            do i = 1, size(density_qn_evaluate_dm_outputs)
                if (btest(request, i - 1) .and. &
                    any(abs(outputs(:, :, i) - &
                            density_qn_evaluate_dm_factors(i) * dm_ao) > tol)) then
                    write(stderr, *) "test_"//test_name//" failed: "// &
                        trim(capitalized(density_qn_evaluate_dm_outputs(i)))// &
                        " returned"//message//" for closed-shell case wrong"// &
                        requested//"."
                    test_passed = .false.
                end if
            end do
        end do

    end function check_density_qn_evaluate_dm_cs_funptr

    function check_density_qn_evaluate_dm_os_funptr( &
        evaluate_dm_funptr, test_name, message, context) result(test_passed)
        !
        ! this function tests a provided density matrix evaluating function pointer
        ! with Coulomb, exact-exchange and non-linear potential contributions for the
        ! open-shell case for every combination of requested outputs
        !
        use otr_density_qn, only: density_qn_evaluate_dm_os_type
        use otr_common_test_reference, only: n_ao, n_particle, request_label
        use test_reference, only: tol

        procedure(density_qn_evaluate_dm_os_type), intent(in), pointer :: &
            evaluate_dm_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
        logical :: test_passed

        real(rp), target :: &
            dm_ao(n_ao, n_ao, n_particle), &
            outputs(n_ao, n_ao, n_particle, size(density_qn_evaluate_dm_outputs))
        real(rp), pointer :: fock(:, :, :), v_coulomb(:, :, :), v_exchange(:, :, :), &
                             v_nonlinear(:, :, :)
        real(rp) :: energy
        character(len=:), allocatable :: requested
        integer(ip) :: request, i, error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(evaluate_dm_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function with Coulomb, exact-exchange and non-linear "// &
                "potential contributions for open-shell case provided"//message// &
                " not associated with value."
            return
        end if

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating subroutine for every combination of requested
        ! outputs
        do request = 0, 15
            requested = request_label(density_qn_evaluate_dm_outputs, request)
            nullify(fock, v_coulomb, v_exchange, v_nonlinear)
            if (btest(request, 0)) fock => outputs(:, :, :, 1)
            if (btest(request, 1)) v_coulomb => outputs(:, :, :, 2)
            if (btest(request, 2)) v_exchange => outputs(:, :, :, 3)
            if (btest(request, 3)) v_nonlinear => outputs(:, :, :, 4)
            energy = 0.0_rp
            outputs = 0.0_rp
            call evaluate_dm_funptr(dm_ao, energy, fock, v_coulomb, v_exchange, &
                                    v_nonlinear, error, context)

            ! check for error
            if (error /= 0) then
                write(stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for open-shell case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol) then
                write(stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for open-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check requested outputs
            do i = 1, size(density_qn_evaluate_dm_outputs)
                if (btest(request, i - 1) .and. &
                    any(abs(outputs(:, :, :, i) - &
                            density_qn_evaluate_dm_factors(i) * dm_ao) > tol)) then
                    write(stderr, *) "test_"//test_name//" failed: "// &
                        trim(capitalized(density_qn_evaluate_dm_outputs(i)))// &
                        " returned"//message//" for open-shell case wrong"// &
                        requested//"."
                    test_passed = .false.
                end if
            end do
        end do

    end function check_density_qn_evaluate_dm_os_funptr

    function check_density_qn_evaluate_dm_c_funptr(evaluate_dm_c_funptr, n_particle, &
                                                   test_name, message, context_c) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating C function pointer
        ! with Coulomb, exact-exchange and non-linear potential contributions for the
        ! closed- (n_particle = 1) or open-shell (n_particle = 2) case for every
        ! combination of requested outputs
        !
        use otr_density_qn_c_interface, only: density_qn_evaluate_dm_c_type
        use otr_common_test_reference, only: n_ao, request_label, shell_names
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: evaluate_dm_c_funptr
        integer(ip), intent(in) :: n_particle
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(density_qn_evaluate_dm_c_type), pointer :: evaluate_dm_funptr
        real(c_rp), target :: &
            dm_ao(n_ao, n_ao, n_particle), &
            outputs(n_ao, n_ao, n_particle, size(density_qn_evaluate_dm_outputs))
        real(c_rp), pointer :: fock(:, :, :), v_coulomb(:, :, :), v_exchange(:, :, :), &
                               v_nonlinear(:, :, :)
        real(c_rp) :: energy
        character(len=:), allocatable :: requested, shell
        integer(ip) :: request, i
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! name of the shell for the messages
        shell = trim(shell_names(n_particle))

        ! check if function pointer is associated
        if (.not. c_associated(evaluate_dm_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function with Coulomb, exact-exchange and non-linear "// &
                "potential contributions for "//shell//" case provided"//message// &
                " not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=evaluate_dm_funptr)

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating function for every combination of requested
        ! outputs
        do request = 0, 15
            requested = request_label(density_qn_evaluate_dm_outputs, request)
            nullify(fock, v_coulomb, v_exchange, v_nonlinear)
            if (btest(request, 0)) fock => outputs(:, :, :, 1)
            if (btest(request, 1)) v_coulomb => outputs(:, :, :, 2)
            if (btest(request, 2)) v_exchange => outputs(:, :, :, 3)
            if (btest(request, 3)) v_nonlinear => outputs(:, :, :, 4)
            energy = 0.0_c_rp
            outputs = 0.0_c_rp
            error = evaluate_dm_funptr(dm_ao, energy, fock, v_coulomb, v_exchange, &
                                       v_nonlinear, context_c)

            ! check for error
            if (error /= 0) then
                write(stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for "//shell//" case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol_c) then
                write(stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for "//shell//" case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check requested outputs
            do i = 1, size(density_qn_evaluate_dm_outputs)
                if (btest(request, i - 1) .and. &
                    any(abs(outputs(:, :, :, i) - &
                            density_qn_evaluate_dm_factors(i) * dm_ao) > tol_c)) then
                    write(stderr, *) "test_"//test_name//" failed: "// &
                        trim(capitalized(density_qn_evaluate_dm_outputs(i)))// &
                        " returned"//message//" for "//shell//" case wrong"// &
                        requested//"."
                    test_passed = .false.
                end if
            end do
        end do

    end function check_density_qn_evaluate_dm_c_funptr

    subroutine get_reference_density_qn_values(values_out) bind(C)
        !
        ! this subroutine exports the reference density QN settings as C settings whose
        ! fields are set one by one without the conversion routines under test, the
        ! logging function and the host context are set to distinct addresses
        !
        use otr_density_qn_c_interface, only: density_qn_settings_type_c
        use test_reference, only: ref_character_to_c

        type(density_qn_settings_type_c), intent(out) :: values_out

        ! set logging function and host context to distinct addresses
        values_out%logger = transfer(1_c_intptr_t, c_null_funptr)
        values_out%context = transfer(2_c_intptr_t, c_null_ptr)

        ! set logical, integer and keyword
        values_out%initialized = .true.
        values_out%verbose = int(ref_density_qn_settings%verbose, kind=c_ip)
        values_out%qn_method = ref_character_to_c(ref_density_qn_settings%qn_method)

    end subroutine get_reference_density_qn_values

    subroutine reference_density_qn_field(name_c, field_value, keyword_c) bind(C)
        !
        ! this subroutine returns the value of a field of the C density QN settings
        ! exported by get_reference_density_qn_values by its name, a numeric, logical
        ! (as 1 or 0) or pointer (as its address) field in field_value and a keyword
        ! field as a null-terminated character array in keyword_c, so that the Python
        ! tests need no values of their own
        !
        use otr_density_qn_c_interface, only: density_qn_settings_type_c
        use test_reference, only: ref_character_from_c

        character(kind=c_char), intent(in) :: name_c(*)
        real(c_rp), intent(out) :: field_value
        character(kind=c_char), intent(out) :: keyword_c(kw_len + 1)

        type(density_qn_settings_type_c) :: values

        ! export reference values and select field by name
        field_value = 0.0_c_rp
        keyword_c = c_null_char
        call get_reference_density_qn_values(values)
        select case (ref_character_from_c(name_c))
        case ("logger")
            field_value = real(transfer(values%logger, 0_c_intptr_t), kind=c_rp)
        case ("context")
            field_value = real(transfer(values%context, 0_c_intptr_t), kind=c_rp)
        case ("initialized")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%initialized)
        case ("verbose")
            field_value = values%verbose
        case ("qn_method")
            keyword_c = values%qn_method
        end select

    end subroutine reference_density_qn_field

    logical(c_bool) function is_default_density_qn_settings(settings_c) bind(C)
        !
        ! this function checks whether C density QN settings hold the default values
        ! without logging function and host context, so that a C test can compare
        ! against them without values of its own
        !
        use otr_density_qn, only: default_density_qn_settings
        use otr_density_qn_c_interface, only: density_qn_settings_type_c

        type(density_qn_settings_type_c), intent(in) :: settings_c

        is_default_density_qn_settings = &
            settings_c == default_density_qn_settings .and. &
            .not. c_associated(settings_c%logger) .and. &
            .not. c_associated(settings_c%context)

    end function is_default_density_qn_settings

    integer(c_size_t) function density_qn_settings_size() bind(C)
        !
        ! this function returns the size of the C density QN settings, so that the
        ! layout tests detect a field missing from the C or Python settings
        !
        use otr_density_qn_c_interface, only: density_qn_settings_type_c

        type(density_qn_settings_type_c) :: settings_c

        density_qn_settings_size = c_sizeof(settings_c)

    end function density_qn_settings_size

    logical function equal_density_qn(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare density QN
        ! settings to different density QN settings
        !
        use otr_common_test_reference, only: operator(==)

        type(density_qn_settings_type), intent(in) :: lhs, rhs

        equal_density_qn = &
            lhs%orbital_settings_type == rhs%orbital_settings_type .and. &
            lhs%qn_method == rhs%qn_method

    end function equal_density_qn

    logical function not_equal_density_qn(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare density QN
        ! settings to different density QN settings
        !
        type(density_qn_settings_type), intent(in) :: lhs, rhs

        not_equal_density_qn = .not. (lhs == rhs)

    end function not_equal_density_qn

    logical function equal_density_qn_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare C density QN
        ! settings to Fortran density QN settings field by field, independently of the
        ! conversion routines under test
        !
        use otr_density_qn_c_interface, only: density_qn_settings_type_c
        use test_reference, only: ref_character_from_c

        type(density_qn_settings_type_c), intent(in) :: lhs_c
        type(density_qn_settings_type), intent(in) :: rhs

        equal_density_qn_c = &
            lhs_c%initialized .and. lhs_c%verbose == rhs%verbose .and. &
            ref_character_from_c(lhs_c%qn_method) == trim(rhs%qn_method)

    end function equal_density_qn_c

    logical function not_equal_density_qn_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare C density
        ! QN settings to Fortran density QN settings
        !
        use otr_density_qn_c_interface, only: density_qn_settings_type_c

        type(density_qn_settings_type_c), intent(in) :: lhs_c
        type(density_qn_settings_type), intent(in) :: rhs

        not_equal_density_qn_c = .not. (lhs_c == rhs)

    end function not_equal_density_qn_c

end module otr_density_qn_test_reference
