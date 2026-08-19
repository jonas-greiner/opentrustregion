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

        function test_update_dm_cs_funptr(update_dm_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix updating function pointer with 
        ! a separate non-linear potential contribution for the closed-shell case
        !
        use otr_arh, only: update_dm_cs_type
        use test_reference, only: tol

        procedure(update_dm_cs_type), intent(in), pointer :: update_dm_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :), fock(:, :), v_nonlinear(:, :)
        real(rp) :: energy
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(update_dm_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "updating function with non-linear potential contribution for "// &
                "closed-shell case not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), fock(n_ao, n_ao), v_nonlinear(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call density matrix updating subroutine
        call update_dm_funptr(dm_ao, energy, fock, v_nonlinear, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for closed-shell case."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 9.0_rp) > tol) then
            write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! check Fock matrix
        if (any(abs(fock - 2.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! check non-linear potential
        if (any(abs(v_nonlinear - 3.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Non-linear potential "// &
                "returned"//message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, fock, v_nonlinear)

    end function test_update_dm_cs_funptr

    function test_update_dm_cs_c_funptr(update_dm_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix updating C function pointer 
        ! with a separate non-linear potential contribution for the closed-shell case
        !
        use otr_arh_c_interface, only: update_dm_cs_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: update_dm_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(update_dm_cs_c_type), pointer :: update_dm_funptr
        real(c_rp), allocatable :: dm_ao(:, :), fock(:, :), v_nonlinear(:, :)
        real(c_rp) :: energy
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(update_dm_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "updating function with non-linear potential contribution for "// &
                "closed-shell case not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=update_dm_c_funptr, fptr=update_dm_funptr)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), fock(n_ao, n_ao), v_nonlinear(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call density matrix updating function
        error = update_dm_funptr(dm_ao, energy, fock, v_nonlinear)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for closed-shell case."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 9.0_c_rp) > tol_c) then
            write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! check Fock matrix
        if (any(abs(fock - 2.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! check non-linear potential
        if (any(abs(v_nonlinear - 3.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Non-linear potential "// &
                "returned"//message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, fock, v_nonlinear)

    end function test_update_dm_cs_c_funptr

    function test_update_dm_os_funptr(update_dm_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix function pointer with same- and 
        ! opposite-spin potential contributions
        !
        use otr_arh, only: update_dm_os_type
        use test_reference, only: tol

        procedure(update_dm_os_type), intent(in), pointer :: update_dm_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :, :), fock(:, :, :), v_same_spin(:, :, :), &
                                 v_opposite_spin(:, :, :), v_nonlinear(:, :, :)
        real(rp) :: energy
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(update_dm_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "updating function with same- and opposite-spin potential "// &
                "contributions not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle), &
                 v_same_spin(n_ao, n_ao, n_particle), &
                 v_opposite_spin(n_ao, n_ao, n_particle), &
                 v_nonlinear(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call density matrix updating subroutine
        call update_dm_funptr(dm_ao, energy, fock, v_same_spin, v_opposite_spin, &
                              v_nonlinear, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 18.0_rp) > tol) then
            write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! check Fock matrix
        if (any(abs(fock - 2.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! check same-spin potential
        if (any(abs(v_same_spin - 3.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Same-spin potential "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

        ! check opposite-spin potential
        if (any(abs(v_opposite_spin - 4.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Opposite-spin "// &
                "potential returned"//message//" wrong."
            test_passed = .false.
        end if

        ! check non-linear potential
        if (any(abs(v_nonlinear - 5.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Non-linear potential "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, fock, v_same_spin, v_opposite_spin, v_nonlinear)

    end function test_update_dm_os_funptr

    function test_update_dm_os_c_funptr(update_dm_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix updating C function pointer
        ! with same- and opposite-spin potential contributions
        !
        use otr_arh_c_interface, only: update_dm_os_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: update_dm_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(update_dm_os_c_type), pointer :: update_dm_funptr
        real(c_rp), allocatable :: dm_ao(:, :, :), fock(:, :, :), &
                                   v_same_spin(:, :, :), v_opposite_spin(:, :, :), &
                                   v_nonlinear(:, :, :)
        real(c_rp) :: energy
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(update_dm_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "updating function with same- and opposite-spin potential "// &
                "contributions not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=update_dm_c_funptr, fptr=update_dm_funptr)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle), &
                 v_same_spin(n_ao, n_ao, n_particle), &
                 v_opposite_spin(n_ao, n_ao, n_particle), &
                 v_nonlinear(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call density matrix updating function
        error = update_dm_funptr(dm_ao, energy, fock, v_same_spin, v_opposite_spin, &
                                 v_nonlinear)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 18.0_c_rp) > tol_c) then
            write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! check Fock matrix
        if (any(abs(fock - 2.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! check same-spin potential
        if (any(abs(v_same_spin - 3.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Same-spin potential "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

        ! check opposite-spin potential
        if (any(abs(v_opposite_spin - 4.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Opposite-spin "// &
                "potential returned"//message//" wrong."
            test_passed = .false.
        end if

        ! check non-linear potential
        if (any(abs(v_nonlinear - 5.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Non-linear potential "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, fock, v_same_spin, v_opposite_spin, v_nonlinear)

    end function test_update_dm_os_c_funptr

    subroutine get_reference_arh_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the ARH reference values for tests
        !
        type(ref_arh_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_arh_settings

    end subroutine get_reference_arh_values

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
