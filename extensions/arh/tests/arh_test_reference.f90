! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_test_reference

    use opentrustregion, only: ip, rp, kw_len, stderr
    use c_interface, only: c_ip, c_rp
    use otr_oao_test_reference, only: n_particle, n_ao, ref_oao_settings_type, &
                                      ref_oao_settings
    use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_funptr, c_f_procpointer, &
                                           c_associated

    implicit none

    ! derived types for ARH settings
    type, extends(ref_oao_settings_type) :: ref_arh_settings_type
        character(kw_len, c_char) :: arh_type
    end type

    type, bind(C) :: ref_arh_settings_type_c
        integer(c_ip) :: verbose
        character(c_char) :: arh_type(kw_len + 1)
    end type

    ! general reference parameters
    type(ref_arh_settings_type), parameter :: ref_arh_settings = &
        ref_arh_settings_type(ref_oao_settings_type = ref_oao_settings, &
                              arh_type = "symm_arh")

    ! multiples of the density matrix the mock density matrix evaluating functions
    ! return for each optional output, in the order of evaluate_dm_*_outputs
    real(c_rp), protected, bind(C, name="test_evaluate_dm_cs_factors") :: &
        evaluate_dm_cs_factors(2) = [2.0_c_rp, 3.0_c_rp]
    real(c_rp), protected, bind(C, name="test_evaluate_dm_os_factors") :: &
        evaluate_dm_os_factors(4) = [2.0_c_rp, 3.0_c_rp, 4.0_c_rp, 5.0_c_rp]

    ! optional outputs of the density matrix evaluating functions in the order of their
    ! argument lists
    character(20), parameter :: evaluate_dm_cs_outputs(2) = &
        [character(20) :: "Fock matrix", "non-linear potential"]
    character(23), parameter :: evaluate_dm_os_outputs(4) = &
        [character(23) :: "Fock matrix", "same-spin potential", &
         "opposite-spin potential", "non-linear potential"]

    interface assignment(=)
        module procedure assign_ref_to_arh
        module procedure assign_ref_to_arh_c
        module procedure assign_ref_to_ref_c
    end interface

    interface operator(==)
        module procedure equal_arh_to_ref
        module procedure equal_arh_c_to_ref
        module procedure equal_arh
        module procedure equal_arh_c
    end interface

    interface operator(/=)
        module procedure not_equal_arh_to_ref
        module procedure not_equal_arh_c_to_ref
        module procedure not_equal_arh
        module procedure not_equal_arh_c
    end interface

contains

    function test_evaluate_dm_cs_funptr(evaluate_dm_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating function pointer
        ! with a separate non-linear potential contribution for the closed-shell case
        ! for every combination of requested outputs
        !
        use otr_arh, only: evaluate_dm_cs_type
        use otr_oao_test_reference, only: request_label, capitalized
        use test_reference, only: tol

        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), target :: dm_ao(n_ao, n_ao), outputs(n_ao, n_ao, 2)
        real(rp), pointer, contiguous :: fock(:, :), v_nonlinear(:, :)
        real(rp) :: energy
        character(:), allocatable :: requested
        integer(ip) :: request, i, error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(evaluate_dm_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function with non-linear potential contribution for "// &
                "closed-shell case provided"//message//" not associated with value."
            return
        end if

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating subroutine for every combination of requested
        ! outputs
        do request = 0, 3
            requested = request_label(evaluate_dm_cs_outputs, request)
            nullify(fock, v_nonlinear)
            if (btest(request, 0)) fock => outputs(:, :, 1)
            if (btest(request, 1)) v_nonlinear => outputs(:, :, 2)
            energy = 0.0_rp
            outputs = 0.0_rp
            call evaluate_dm_funptr(dm_ao, energy, fock, v_nonlinear, error)

            ! check for error
            if (error /= 0) then
                write (stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for closed-shell case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol) then
                write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for closed-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check requested outputs
            do i = 1, size(evaluate_dm_cs_outputs)
                if (btest(request, i - 1) .and. any( &
                    abs(outputs(:, :, i) - evaluate_dm_cs_factors(i) * dm_ao) > tol)) &
                    then
                    write (stderr, *) "test_"//test_name//" failed: "// &
                        trim(capitalized(evaluate_dm_cs_outputs(i)))//" returned"// &
                        message//" for closed-shell case wrong"//requested//"."
                    test_passed = .false.
                end if
            end do
        end do

    end function test_evaluate_dm_cs_funptr

    function test_evaluate_dm_cs_c_funptr(evaluate_dm_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating C function pointer
        ! with a separate non-linear potential contribution for the closed-shell case
        ! for every combination of requested outputs
        !
        use otr_arh_c_interface, only: evaluate_dm_cs_c_type
        use otr_oao_test_reference, only: request_label, capitalized
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: evaluate_dm_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed


        procedure(evaluate_dm_cs_c_type), pointer :: evaluate_dm_funptr
        real(c_rp), target :: dm_ao(n_ao, n_ao), outputs(n_ao, n_ao, 2)
        real(c_rp), pointer :: fock(:, :), v_nonlinear(:, :)
        real(c_rp) :: energy
        character(:), allocatable :: requested
        integer(ip) :: request, i
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(evaluate_dm_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function with non-linear potential contribution for "// &
                "closed-shell case provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=evaluate_dm_funptr)

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating function for every combination of requested
        ! outputs
        do request = 0, 3
            requested = request_label(evaluate_dm_cs_outputs, request)
            nullify(fock, v_nonlinear)
            if (btest(request, 0)) fock => outputs(:, :, 1)
            if (btest(request, 1)) v_nonlinear => outputs(:, :, 2)
            energy = 0.0_c_rp
            outputs = 0.0_c_rp
            error = evaluate_dm_funptr(dm_ao, energy, fock, v_nonlinear)

            ! check for error
            if (error /= 0) then
                write (stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for closed-shell case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol_c) then
                write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for closed-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check requested outputs
            do i = 1, size(evaluate_dm_cs_outputs)
                if (btest(request, i - 1) .and. &
                    any(abs(outputs(:, :, i) - evaluate_dm_cs_factors(i) * dm_ao) > &
                        tol_c)) then
                    write (stderr, *) "test_"//test_name//" failed: "// &
                        trim(capitalized(evaluate_dm_cs_outputs(i)))//" returned"// &
                        message//" for closed-shell case wrong"//requested//"."
                    test_passed = .false.
                end if
            end do
        end do

    end function test_evaluate_dm_cs_c_funptr

    function test_evaluate_dm_os_funptr(evaluate_dm_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating function pointer
        ! with same- and opposite-spin potential contributions for every combination of
        ! requested outputs
        !
        use otr_arh, only: evaluate_dm_os_type
        use otr_oao_test_reference, only: request_label, capitalized
        use test_reference, only: tol

        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed


        real(rp), target :: dm_ao(n_ao, n_ao, n_particle), &
                            outputs(n_ao, n_ao, n_particle, 4)
        real(rp), pointer :: fock(:, :, :), v_same_spin(:, :, :), &
                             v_opposite_spin(:, :, :), v_nonlinear(:, :, :)
        real(rp) :: energy
        character(:), allocatable :: requested
        integer(ip) :: request, i, error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(evaluate_dm_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function with same- and opposite-spin potential "// &
                "contributions provided"//message//" not associated with value."
            return
        end if

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating subroutine for every combination of requested
        ! outputs
        do request = 0, 15
            requested = request_label(evaluate_dm_os_outputs, request)
            nullify(fock, v_same_spin, v_opposite_spin, v_nonlinear)
            if (btest(request, 0)) fock => outputs(:, :, :, 1)
            if (btest(request, 1)) v_same_spin => outputs(:, :, :, 2)
            if (btest(request, 2)) v_opposite_spin => outputs(:, :, :, 3)
            if (btest(request, 3)) v_nonlinear => outputs(:, :, :, 4)
            energy = 0.0_rp
            outputs = 0.0_rp
            call evaluate_dm_funptr(dm_ao, energy, fock, v_same_spin, v_opposite_spin, &
                                    v_nonlinear, error)

            ! check for error
            if (error /= 0) then
                write (stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for open-shell case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol) then
                write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for open-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check requested outputs
            do i = 1, size(evaluate_dm_os_outputs)
                if (btest(request, i - 1) .and. &
                    any(abs(outputs(:, :, :, i) - evaluate_dm_os_factors(i) * dm_ao) > &
                        tol)) then
                    write (stderr, *) "test_"//test_name//" failed: "// &
                        trim(capitalized(evaluate_dm_os_outputs(i)))//" returned"// &
                        message//" for open-shell case wrong"//requested//"."
                    test_passed = .false.
                end if
            end do
        end do

    end function test_evaluate_dm_os_funptr

    function test_evaluate_dm_os_c_funptr(evaluate_dm_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating C function pointer
        ! with same- and opposite-spin potential contributions for every combination of
        ! requested outputs
        !
        use otr_arh_c_interface, only: evaluate_dm_os_c_type
        use otr_oao_test_reference, only: request_label, capitalized
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: evaluate_dm_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed


        procedure(evaluate_dm_os_c_type), pointer :: evaluate_dm_funptr
        real(c_rp), target :: dm_ao(n_ao, n_ao, n_particle), &
                              outputs(n_ao, n_ao, n_particle, 4)
        real(c_rp), pointer :: fock(:, :, :), v_same_spin(:, :, :), &
                               v_opposite_spin(:, :, :), v_nonlinear(:, :, :)
        real(c_rp) :: energy
        character(:), allocatable :: requested
        integer(ip) :: request, i
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(evaluate_dm_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function with same- and opposite-spin potential "// &
                "contributions provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=evaluate_dm_funptr)

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating function for every combination of requested
        ! outputs
        do request = 0, 15
            requested = request_label(evaluate_dm_os_outputs, request)
            nullify(fock, v_same_spin, v_opposite_spin, v_nonlinear)
            if (btest(request, 0)) fock => outputs(:, :, :, 1)
            if (btest(request, 1)) v_same_spin => outputs(:, :, :, 2)
            if (btest(request, 2)) v_opposite_spin => outputs(:, :, :, 3)
            if (btest(request, 3)) v_nonlinear => outputs(:, :, :, 4)
            energy = 0.0_c_rp
            outputs = 0.0_c_rp
            error = evaluate_dm_funptr(dm_ao, energy, fock, v_same_spin, &
                                       v_opposite_spin, v_nonlinear)

            ! check for error
            if (error /= 0) then
                write (stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for open-shell case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol_c) then
                write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for open-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check requested outputs
            do i = 1, size(evaluate_dm_os_outputs)
                if (btest(request, i - 1) .and. &
                    any(abs(outputs(:, :, :, i) - evaluate_dm_os_factors(i) * dm_ao) > &
                        tol_c)) then
                    write (stderr, *) "test_"//test_name//" failed: "// &
                        trim(capitalized(evaluate_dm_os_outputs(i)))//" returned"// &
                        message//" for open-shell case wrong"//requested//"."
                    test_passed = .false.
                end if
            end do
        end do

    end function test_evaluate_dm_os_c_funptr

    subroutine get_reference_arh_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the ARH reference values for tests
        !
        type(ref_arh_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_arh_settings

    end subroutine get_reference_arh_values

    subroutine get_default_arh_values(default_values_out) bind(C)
        !
        ! this subroutine exports the default values for tests
        !
        use otr_arh, only: default_arh_settings
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)

        type(arh_settings_type_c), intent(out) :: default_values_out

        default_values_out = default_arh_settings

    end subroutine get_default_arh_values

    subroutine assign_ref_to_arh(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set ARH settings to
        ! reference values
        !
        use otr_arh, only: arh_settings_type
        use otr_oao_test_reference, only: assignment(=)

        type(arh_settings_type), intent(out) :: lhs
        type(ref_arh_settings_type), intent(in) :: rhs

        ! set OAO settings using type extension
        lhs%oao_settings_type = rhs%ref_oao_settings_type

        ! set reference values
        lhs%arh_type = rhs%arh_type

    end subroutine assign_ref_to_arh

    subroutine assign_ref_to_arh_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C ARH settings to
        ! reference values
        !
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh, only: arh_settings_type

        type(arh_settings_type_c), intent(out) :: lhs_c
        type(ref_arh_settings_type), intent(in) :: rhs

        type(arh_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_arh_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        use c_interface, only: character_to_c

        type(ref_arh_settings_type_c), intent(out) :: lhs
        type(ref_arh_settings_type), intent(in) :: rhs

        lhs%verbose = int(rhs%verbose, kind=c_ip)
        lhs%arh_type = character_to_c(rhs%arh_type)

    end subroutine assign_ref_to_ref_c

    logical function equal_arh_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare ARH settings to
        ! reference values
        !
        use otr_arh, only: arh_settings_type
        use otr_oao_test_reference, only: operator(==)

        type(arh_settings_type), intent(in) :: lhs
        type(ref_arh_settings_type), intent(in) :: rhs

        equal_arh_to_ref = lhs%oao_settings_type == rhs%ref_oao_settings_type .and. &
                           (lhs%arh_type == rhs%arh_type)

    end function equal_arh_to_ref

    logical function not_equal_arh_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare ARH
        ! settings to reference values
        !
        use otr_arh, only: arh_settings_type

        type(arh_settings_type), intent(in) :: lhs
        type(ref_arh_settings_type), intent(in) :: rhs

        not_equal_arh_to_ref = .not. (lhs == rhs)

    end function not_equal_arh_to_ref

    logical function equal_arh_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare ARH settings to
        ! reference values
        !
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh, only: arh_settings_type

        type(arh_settings_type_c), intent(in) :: lhs_c
        type(ref_arh_settings_type), intent(in) :: rhs

        type(arh_settings_type) :: lhs

        lhs = lhs_c
        equal_arh_c_to_ref = lhs == rhs

    end function equal_arh_c_to_ref

    logical function not_equal_arh_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare ARH
        ! settings to reference values
        !
        use otr_arh_c_interface, only: arh_settings_type_c

        type(arh_settings_type_c), intent(in) :: lhs
        type(ref_arh_settings_type), intent(in) :: rhs
        
        not_equal_arh_c_to_ref = .not. (lhs == rhs)

    end function not_equal_arh_c_to_ref

    logical function equal_arh(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare ARH settings to
        ! different ARH settings
        !
        use otr_arh, only: arh_settings_type
        use otr_oao_test_reference, only: operator(==)

        type(arh_settings_type), intent(in) :: lhs, rhs
        
        equal_arh = lhs%oao_settings_type == rhs%oao_settings_type .and. &
                    (lhs%arh_type == rhs%arh_type)

    end function equal_arh

    logical function not_equal_arh(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare ARH
        ! settings to different ARH settings
        !
        use otr_arh, only: arh_settings_type

        type(arh_settings_type), intent(in) :: lhs, rhs

        not_equal_arh = .not. (lhs == rhs)

    end function not_equal_arh

    logical function equal_arh_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare ARH settings to
        ! different ARH settings
        !
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh, only: arh_settings_type

        type(arh_settings_type_c), intent(in) :: lhs_c
        type(arh_settings_type), intent(in) :: rhs
        
        type(arh_settings_type) :: lhs

        lhs = lhs_c
        equal_arh_c = lhs == rhs

    end function equal_arh_c

    logical function not_equal_arh_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare ARH
        ! settings to different ARH settings
        !
        use otr_arh_c_interface, only: arh_settings_type_c
        use otr_arh, only: arh_settings_type

        type(arh_settings_type_c), intent(in) :: lhs_c
        type(arh_settings_type), intent(in) :: rhs
        
        not_equal_arh_c = .not. (lhs_c == rhs)

    end function not_equal_arh_c

end module otr_arh_test_reference
