! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_test_reference

    use opentrustregion, only: ip, rp, stderr
    use c_interface, only: c_ip, c_rp
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer

    implicit none

    ! number of AOs
    integer(ip), parameter :: n_particle = 2_ip, n_ao = 3_ip
    integer(c_ip), parameter :: n_particle_c = int(n_particle, kind=c_ip), &
                                n_ao_c = int(n_ao, kind=c_ip)

    ! derived types for ARH settings
    type ref_arh_settings_type
        logical :: restricted
        integer(ip) :: verbose
    end type

    type, bind(C) :: ref_arh_settings_type_c
        logical(c_bool) :: restricted
        integer(c_ip) :: verbose
    end type

    ! general reference parameters
    type(ref_arh_settings_type) :: ref_arh_settings = &
        ref_arh_settings_type(restricted = .true., verbose = 3)

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

    function test_get_energy_closed_shell_funptr(get_energy_funptr, test_name, &
                                                 message) result(test_passed)
        !
        ! this function tests a provided energy function pointer for the closed-shell 
        ! case
        !
        use otr_arh, only: get_energy_closed_shell_type
        use test_reference, only: tol

        procedure(get_energy_closed_shell_type), intent(in), pointer :: &
            get_energy_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :)
        real(rp) :: energy
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call energy function
        energy = get_energy_funptr(dm_ao, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 9.0_rp) > tol) then
            write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_energy_closed_shell_funptr

    function test_get_energy_open_shell_funptr(get_energy_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided energy function pointer for the open-shell case
        !
        use otr_arh, only: get_energy_open_shell_type
        use test_reference, only: tol

        procedure(get_energy_open_shell_type), intent(in), pointer :: get_energy_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :, :)
        real(rp) :: energy
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call energy function
        energy = get_energy_funptr(dm_ao, error)

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

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_energy_open_shell_funptr

    function test_get_fock_funptr(get_fock_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided Fock matrix function pointer
        !
        use otr_arh, only: get_fock_type
        use test_reference, only: tol

        procedure(get_fock_type), intent(in), pointer :: get_fock_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :), fock(:, :)
        real(rp) :: energy
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), fock(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call Fock matrix subroutine
        call get_fock_funptr(dm_ao, energy, fock, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 9.0_rp) > tol) then
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

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_fock_funptr

    function test_get_fock_jk_funptr(get_fock_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided Fock, Coulomb, and exchange matrix function 
        ! pointer
        !
        use otr_arh, only: get_fock_jk_type
        use test_reference, only: tol

        procedure(get_fock_jk_type), intent(in), pointer :: get_fock_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :, :), fock(:, :, :), coulomb(:, :, :), &
                                 exchange(:, :, :)
        real(rp) :: energy
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle), &
                 coulomb(n_ao, n_ao, n_particle), exchange(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call Fock matrix subroutine
        call get_fock_funptr(dm_ao, energy, fock, coulomb, exchange, error)

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

        ! check Coulomb matrix
        if (any(abs(coulomb - 3.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Coulomb matrix returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! check exchange matrix
        if (any(abs(exchange - 4.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Exchange matrix "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_fock_jk_funptr

    function test_get_energy_closed_shell_c_funptr(get_energy_c_funptr, test_name, &
                                                 message) result(test_passed)
        !
        ! this function tests a provided energy C function pointer for the closed-shell 
        ! case
        !
        use otr_arh_c_interface, only: get_energy_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_energy_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(get_energy_c_type), pointer :: get_energy_funptr
        real(c_rp), allocatable :: dm_ao(:, :)
        real(c_rp) :: energy
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_energy_c_funptr, fptr=get_energy_funptr)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call energy function
        error = get_energy_funptr(dm_ao, energy)

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

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_energy_closed_shell_c_funptr

    function test_get_energy_open_shell_c_funptr(get_energy_c_funptr, test_name, &
                                                   message) result(test_passed)
        !
        ! this function tests a provided energy C function pointer for the open-shell  
        ! case
        !
        use otr_arh_c_interface, only: get_energy_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_energy_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(get_energy_c_type), pointer :: get_energy_funptr
        real(c_rp), allocatable :: dm_ao(:, :, :)
        real(c_rp) :: energy
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_energy_c_funptr, fptr=get_energy_funptr)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle))

        ! initialize orbital update
        dm_ao = 1.0_c_rp

        ! call energy function
        error = get_energy_funptr(dm_ao, energy)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for open-shell case."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 18.0_c_rp) > tol_c) then
            write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_energy_open_shell_c_funptr

    function test_get_fock_c_funptr(get_fock_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided Fock matrix C function pointer
        !
        use otr_arh_c_interface, only: get_fock_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_fock_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(get_fock_c_type), pointer :: get_fock_funptr
        real(c_rp), allocatable :: dm_ao(:, :), fock(:, :)
        real(c_rp) :: energy
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_fock_c_funptr, fptr=get_fock_funptr)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), fock(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call Fock matrix function
        error = get_fock_funptr(dm_ao, energy, fock)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                "."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 9.0_c_rp) > tol_c) then
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

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_fock_c_funptr

    function test_get_fock_jk_c_funptr(get_fock_jk_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided Fock matrix C function pointer
        !
        use otr_arh_c_interface, only: get_fock_jk_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_fock_jk_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(get_fock_jk_c_type), pointer :: get_fock_jk_funptr
        real(c_rp), allocatable :: dm_ao(:, :, :), fock(:, :, :), coulomb(:, :, :), &
                                   exchange(:, :, :)
        real(c_rp) :: energy
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_fock_jk_c_funptr, fptr=get_fock_jk_funptr)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle), &
                 coulomb(n_ao, n_ao, n_particle), exchange(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call Fock matrix function
        error = get_fock_jk_funptr(dm_ao, energy, fock, coulomb, exchange)

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

        ! check Coulomb matrix
        if (any(abs(coulomb - 3.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Coulomb matrix returned"// &
                message//" wrong."
            test_passed = .false.
        end if

        ! check exchange matrix
        if (any(abs(exchange - 4.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Exchange matrix "// &
                "returned"//message//" wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, fock, coulomb, exchange)

    end function test_get_fock_jk_c_funptr

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

        type(arh_settings_type), intent(out) :: lhs
        type(ref_arh_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%logger => null()

        ! set reference values
        lhs%restricted = rhs%restricted
        lhs%verbose = rhs%verbose

        ! set initialization logical
        lhs%initialized = .true.

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

        lhs%restricted = logical(rhs%restricted, kind=c_bool)
        lhs%verbose = int(rhs%verbose, kind=c_ip)

    end subroutine assign_ref_to_ref_c

    logical function equal_arh_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare ARH settings to 
        ! reference values
        !
        use otr_arh, only: arh_settings_type

        type(arh_settings_type), intent(in) :: lhs
        type(ref_arh_settings_type), intent(in) :: rhs

        equal_arh_to_ref = (lhs%restricted .eqv. rhs%restricted) .and. &
                           lhs%verbose == rhs%verbose

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

        type(arh_settings_type), intent(in) :: lhs, rhs
        
        equal_arh = (lhs%restricted .eqv. rhs%restricted) .and. &
                    lhs%verbose == rhs%verbose

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
