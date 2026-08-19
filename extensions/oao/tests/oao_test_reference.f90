! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_test_reference

    use opentrustregion, only: ip, rp, stderr
    use c_interface, only: c_ip, c_rp
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_associated

    implicit none

    ! number of particles and AOs
    integer(ip), parameter :: n_particle = 2_ip, n_ao = 3_ip
    integer(c_ip), parameter :: n_particle_c = int(n_particle, kind=c_ip)
    integer(c_ip), protected, bind(C, name="test_n_ao") :: n_ao_c = int(n_ao, kind=c_ip)

    ! derived types for OAO settings
    type :: ref_oao_settings_type
        integer(ip) :: verbose
    end type

    type, bind(C) :: ref_oao_settings_type_c
        integer(c_ip) :: verbose
    end type

    ! general reference parameters
    type(ref_oao_settings_type), parameter :: ref_oao_settings = &
        ref_oao_settings_type(verbose = 3)

    interface assignment(=)
        module procedure assign_ref_to_ref_c
        module procedure assign_ref_to_oao
        module procedure assign_ref_to_oao_c
    end interface

    interface operator(==)
        module procedure equal_oao_to_ref
        module procedure equal_oao_c_to_ref
        module procedure equal_oao
        module procedure equal_oao_c
    end interface

    interface operator(/=)
        module procedure not_equal_oao_to_ref
        module procedure not_equal_oao_c_to_ref
        module procedure not_equal_oao
        module procedure not_equal_oao_c
    end interface

contains

    function test_get_energy_cs_funptr(get_energy_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided energy function pointer for the closed-shell 
        ! case
        !
        use otr_oao, only: get_energy_cs_type
        use test_reference, only: tol

        procedure(get_energy_cs_type), intent(in), pointer :: get_energy_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :)
        real(rp) :: energy
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(get_energy_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Energy function for "// &
                "closed-shell case not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call energy function
        energy = get_energy_funptr(dm_ao, error)

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

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_energy_cs_funptr

    function test_get_energy_cs_c_funptr(get_energy_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided energy C function pointer for the closed-shell 
        ! case
        !
        use otr_oao_c_interface, only: get_energy_c_type
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

        ! check if function pointer is associated
        if (.not. c_associated(get_energy_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Energy function for "// &
                "closed-shell case not associated with value."
            return
        end if

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

    end function test_get_energy_cs_c_funptr

    function test_get_energy_os_funptr(get_energy_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided energy function pointer for the open-shell case
        !
        use otr_oao, only: get_energy_os_type
        use test_reference, only: tol

        procedure(get_energy_os_type), intent(in), pointer :: get_energy_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :, :)
        real(rp) :: energy
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(get_energy_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Energy function for "// &
                "open-shell case not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call energy function
        energy = get_energy_funptr(dm_ao, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for open-shell case."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 18.0_rp) > tol) then
            write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao)

    end function test_get_energy_os_funptr

    function test_get_energy_os_c_funptr(get_energy_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided energy C function pointer for the open-shell  
        ! case
        !
        use otr_oao_c_interface, only: get_energy_c_type
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

        ! check if function pointer is associated
        if (.not. c_associated(get_energy_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Energy function for "// &
                "open-shell case not associated with value."
            return
        end if

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

    end function test_get_energy_os_c_funptr

    function test_update_dm_cs_funptr(update_dm_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix updating function pointer for 
        ! the closed-shell case
        !
        use otr_oao, only: update_dm_cs_type, get_response_cs_type
        use test_reference, only: tol

        procedure(update_dm_cs_type), intent(in), pointer :: update_dm_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :), fock(:, :)
        real(rp) :: energy
        procedure(get_response_cs_type), pointer :: get_response_funptr
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(update_dm_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "updating function for closed-shell case not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), fock(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call density matrix updating subroutine
        call update_dm_funptr(dm_ao, energy, fock, get_response_funptr, error)

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

        ! deallocate arrays
        deallocate(dm_ao, fock)

        ! test returned response function
        test_passed = test_passed .and. &
            test_get_response_cs_funptr(get_response_funptr, test_name, " by "// &
                                        "response function returned"//message)

    end function test_update_dm_cs_funptr

    function test_update_dm_cs_c_funptr(update_dm_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix updating C function pointer for 
        ! the closed-shell case
        !
        use otr_oao_c_interface, only: update_dm_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: update_dm_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(update_dm_c_type), pointer :: update_dm_funptr
        real(c_rp), allocatable :: dm_ao(:, :), fock(:, :)
        real(c_rp) :: energy
        type(c_funptr) :: get_response_c_funptr
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(update_dm_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "updating function for closed-shell case not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=update_dm_c_funptr, fptr=update_dm_funptr)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), fock(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call density matrix updating function
        error = update_dm_funptr(dm_ao, energy, fock, get_response_c_funptr)

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

        ! deallocate arrays
        deallocate(dm_ao, fock)

        ! test returned response function
        test_passed = test_passed .and. &
            test_get_response_cs_c_funptr(get_response_c_funptr, test_name, " by "// &
                                          "response function returned"//message)

    end function test_update_dm_cs_c_funptr

    function test_update_dm_os_funptr(update_dm_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix updating function pointer for 
        ! the open-shell case
        !
        use otr_oao, only: update_dm_os_type, get_response_os_type
        use test_reference, only: tol

        procedure(update_dm_os_type), intent(in), pointer :: update_dm_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :, :), fock(:, :, :)
        real(rp) :: energy
        procedure(get_response_os_type), pointer :: get_response_funptr
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(update_dm_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "updating function for open-shell case not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call density matrix updating subroutine
        call update_dm_funptr(dm_ao, energy, fock, get_response_funptr, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for open-shell case."
            test_passed = .false.
        end if

        ! check energy
        if (abs(energy - 18.0_rp) > tol) then
            write (stderr, *) "test_"//test_name//" failed: Energy returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! check Fock matrix
        if (any(abs(fock - 2.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, fock)

        ! test returned response function
        test_passed = test_passed .and. &
            test_get_response_os_funptr(get_response_funptr, test_name, " by "// &
                                        "response function returned"//message)

    end function test_update_dm_os_funptr

    function test_update_dm_os_c_funptr(update_dm_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix updating C function pointer for 
        ! the open-shell case
        !
        use otr_oao_c_interface, only: update_dm_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: update_dm_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(update_dm_c_type), pointer :: update_dm_funptr
        real(c_rp), allocatable :: dm_ao(:, :, :), fock(:, :, :)
        real(c_rp) :: energy
        type(c_funptr) :: get_response_c_funptr
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(update_dm_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "updating function for open-shell case not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=update_dm_c_funptr, fptr=update_dm_funptr)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call density matrix updating function
        error = update_dm_funptr(dm_ao, energy, fock, get_response_c_funptr)

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

        ! check Fock matrix
        if (any(abs(fock - 2.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, fock)

        ! test returned response function
        test_passed = test_passed .and. &
            test_get_response_os_c_funptr(get_response_c_funptr, test_name, " by "// &
                                          "response function returned"//message)

    end function test_update_dm_os_c_funptr

    function test_get_response_cs_funptr(get_response_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided response function pointer for the closed-shell 
        ! case
        !
        use otr_oao, only: get_response_cs_type
        use test_reference, only: tol

        procedure(get_response_cs_type), intent(in), pointer :: get_response_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :), response(:, :)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(get_response_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Response function for "// &
                "closed-shell case not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), response(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call response subroutine
        call get_response_funptr(dm_ao, response, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for closed-shell case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - 2.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Response returned"// &
                message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, response)

    end function test_get_response_cs_funptr

    function test_get_response_cs_c_funptr(get_response_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided response C function pointer for the 
        ! closed-shell case
        !
        use otr_oao_c_interface, only: get_response_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_response_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(get_response_c_type), pointer :: get_response_funptr_c
        real(c_rp), allocatable :: dm_ao(:, :), response(:, :)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(get_response_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Response function for "// &
                "closed-shell case not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_response_c_funptr, fptr=get_response_funptr_c)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), response(n_ao, n_ao))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call response function
        error = get_response_funptr_c(dm_ao, response)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for closed-shell case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - 2.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Response returned"// &
                message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, response)

    end function test_get_response_cs_c_funptr

    function test_get_response_os_funptr(get_response_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided response function pointer for the 
        ! open-shell case
        !
        use otr_oao, only: get_response_os_type
        use test_reference, only: tol

        procedure(get_response_os_type), intent(in), pointer :: get_response_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :, :), response(:, :, :)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(get_response_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Response function for "// &
                "open-shell case not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), response(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_rp

        ! call response subroutine
        call get_response_funptr(dm_ao, response, error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for open-shell case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - 2.0_rp) > tol)) then
            write (stderr, *) "test_"//test_name//" failed: Response returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, response)

    end function test_get_response_os_funptr

    function test_get_response_os_c_funptr(get_response_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided response C function pointer for the 
        ! open-shell case
        !
        use otr_oao_c_interface, only: get_response_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_response_c_funptr
        character(*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(get_response_c_type), pointer :: get_response_funptr_c
        real(c_rp), allocatable :: dm_ao(:, :, :), response(:, :, :)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(get_response_c_funptr)) then
            test_passed = .false.
            write (stderr, *) "test_"//test_name//" failed: Response function for "// &
                "open-shell case not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_response_c_funptr, fptr=get_response_funptr_c)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), response(n_ao, n_ao, n_particle))

        ! initialize density matrix
        dm_ao = 1.0_c_rp

        ! call response function
        error = get_response_funptr_c(dm_ao, response)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for open-shell case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - 2.0_c_rp) > tol_c)) then
            write (stderr, *) "test_"//test_name//" failed: Response returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, response)

    end function test_get_response_os_c_funptr

    subroutine get_reference_oao_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the OAO reference values for tests
        !
        type(ref_oao_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_oao_settings

    end subroutine get_reference_oao_values

    subroutine assign_ref_to_oao(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set OAO settings to
        ! reference values
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type), intent(out) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%logger => null()

        ! set reference values
        lhs%verbose = rhs%verbose

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_oao

    subroutine assign_ref_to_oao_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C OAO settings to 
        ! reference values
        !
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(out) :: lhs_c
        type(ref_oao_settings_type), intent(in) :: rhs

        type(oao_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_oao_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        use c_interface, only: character_to_c

        type(ref_oao_settings_type_c), intent(out) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs

        lhs%verbose = int(rhs%verbose, kind=c_ip)

    end subroutine assign_ref_to_ref_c

    logical function equal_oao_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to 
        ! reference values
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type), intent(in) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs

        equal_oao_to_ref = lhs%verbose == rhs%verbose

    end function equal_oao_to_ref

    logical function not_equal_oao_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare OAO
        ! settings to reference values
        !
        use otr_oao, only: oao_settings_type

        class(oao_settings_type), intent(in) :: lhs
        class(ref_oao_settings_type), intent(in) :: rhs

        not_equal_oao_to_ref = .not. (lhs == rhs)

    end function not_equal_oao_to_ref

    logical function equal_oao_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to 
        ! reference values
        !
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(in) :: lhs_c
        type(ref_oao_settings_type), intent(in) :: rhs

        type(oao_settings_type) :: lhs

        lhs = lhs_c
        equal_oao_c_to_ref = lhs == rhs

    end function equal_oao_c_to_ref

    logical function not_equal_oao_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare OAO 
        ! settings to reference values
        !
        use otr_oao_c_interface, only: oao_settings_type_c

        type(oao_settings_type_c), intent(in) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs
        
        not_equal_oao_c_to_ref = .not. (lhs == rhs)

    end function not_equal_oao_c_to_ref

    logical function equal_oao(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to 
        ! different OAO settings
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type), intent(in) :: lhs, rhs
        
        equal_oao = lhs%verbose == rhs%verbose

    end function equal_oao

    logical function not_equal_oao(lhs, rhs) 
        !
        ! this function overloads the negated comparison operator to compare OAO 
        ! settings to different OAO settings
        !
        use otr_oao, only: oao_settings_type

        class(oao_settings_type), intent(in) :: lhs, rhs
        
        not_equal_oao = .not. (lhs == rhs)

    end function not_equal_oao

    logical function equal_oao_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to 
        ! different OAO settings
        !
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(in) :: lhs_c
        type(oao_settings_type), intent(in) :: rhs
        
        type(oao_settings_type) :: lhs

        lhs = lhs_c
        equal_oao_c = lhs == rhs

    end function equal_oao_c

    logical function not_equal_oao_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare OAO 
        ! settings to different OAO settings
        !
        use otr_oao_c_interface, only: oao_settings_type_c
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(in) :: lhs_c
        type(oao_settings_type), intent(in) :: rhs
        
        not_equal_oao_c = .not. (lhs_c == rhs)

    end function not_equal_oao_c

end module otr_oao_test_reference
