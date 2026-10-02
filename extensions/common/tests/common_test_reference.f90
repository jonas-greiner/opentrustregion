! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_test_reference

    use opentrustregion, only: ip, rp, stderr
    use c_interface, only: c_ip, c_rp
    use, intrinsic :: iso_c_binding, only: c_funptr, c_f_procpointer, c_associated, &
                                           c_null_funptr

    implicit none

    ! numbers of particle channels, number of AOs and number of the occupied orbitals
    ! of every channel, which differ between the channels
    integer(ip), parameter :: n_particle = 2_ip, n_ao = 5_ip, &
                              n_occ(n_particle) = [2_ip, 1_ip]
    integer(c_ip), protected, bind(C, name="test_n_ao") :: n_ao_c = int(n_ao, kind=c_ip)

    ! derived type for the settings shared by the extensions parameterizing the
    ! orbitals in an orbital basis
    type :: ref_orbital_settings_type
        integer(ip) :: verbose
    end type

    ! general reference parameters
    type(ref_orbital_settings_type), parameter :: ref_orbital_settings = &
        ref_orbital_settings_type(verbose=3)

    ! multiples of the density matrix the mock density matrix evaluating functions
    ! return for each optional output, in the order of evaluate_dm_outputs
    real(c_rp), protected, bind(C, name="test_evaluate_dm_factors") :: &
        evaluate_dm_factors(2) = [2.0_c_rp, 3.0_c_rp]

    ! optional outputs of the density matrix evaluating function in the order of its
    ! argument list
    character(len=17), parameter :: evaluate_dm_outputs(2) = &
        [character(len=17) :: "Fock matrix", "response function"]

    interface operator(==)
        module procedure equal_orbital
    end interface

contains

    function test_evaluate_dm_cs_funptr(evaluate_dm_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating function pointer for
        ! the closed-shell case for every combination of requested outputs
        !
        use otr_common, only: evaluate_dm_cs_type, get_response_cs_type
        use test_reference, only: tol

        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        character(len=*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), target :: dm_ao(n_ao, n_ao), fock(n_ao, n_ao)
        real(rp), pointer, contiguous :: fock_ptr(:, :)
        real(rp) :: energy
        procedure(get_response_cs_type), pointer :: get_response_funptr
        character(len=:), allocatable :: requested
        integer(ip) :: request, error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(evaluate_dm_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function for closed-shell case provided"//message// &
                " not associated with value."
            return
        end if

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating subroutine for every combination of requested
        ! outputs
        do request = 0, 3
            requested = request_label(evaluate_dm_outputs, request)
            nullify(fock_ptr)
            get_response_funptr => null()
            if (btest(request, 0)) fock_ptr => fock
            energy = 0.0_rp
            fock = 0.0_rp
            if (btest(request, 1)) then
                call evaluate_dm_funptr(dm_ao, energy, fock_ptr, get_response_funptr, &
                                        error)
            else
                call evaluate_dm_funptr(dm_ao, energy, fock_ptr, error=error)
            end if

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

            ! check Fock matrix
            if (btest(request, 0) .and. &
                any(abs(fock - evaluate_dm_factors(1) * dm_ao) > tol)) then
                write(stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                    message//" for closed-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! test returned response function
            if (btest(request, 1)) &
                test_passed = test_passed .and. test_get_response_cs_funptr( &
                    get_response_funptr, test_name, &
                    " by response function returned"//message//requested)
        end do

    end function test_evaluate_dm_cs_funptr

    function test_evaluate_dm_cs_c_funptr(evaluate_dm_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating C function pointer
        ! for the closed-shell case for every combination of requested outputs
        !
        use otr_common_c_interface, only: evaluate_dm_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: evaluate_dm_c_funptr
        character(len=*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(evaluate_dm_c_type), pointer :: evaluate_dm_funptr
        real(c_rp), target :: dm_ao(n_ao, n_ao), fock(n_ao, n_ao)
        real(c_rp), pointer :: fock_ptr(:, :)
        real(c_rp) :: energy
        type(c_funptr) :: get_response_c_funptr
        character(len=:), allocatable :: requested
        integer(ip) :: request
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(evaluate_dm_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function for closed-shell case provided"//message// &
                " not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=evaluate_dm_funptr)

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating function for every combination of requested
        ! outputs
        do request = 0, 3
            requested = request_label(evaluate_dm_outputs, request)
            nullify(fock_ptr)
            get_response_c_funptr = c_null_funptr
            if (btest(request, 0)) fock_ptr => fock
            energy = 0.0_c_rp
            fock = 0.0_c_rp
            if (btest(request, 1)) then
                error = &
                    evaluate_dm_funptr(dm_ao, energy, fock_ptr, get_response_c_funptr)
            else
                error = evaluate_dm_funptr(dm_ao, energy, fock_ptr)
            end if

            ! check for error
            if (error /= 0) then
                write(stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for closed-shell case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol_c) then
                write(stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for closed-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check Fock matrix
            if (btest(request, 0) .and. &
                any(abs(fock - evaluate_dm_factors(1) * dm_ao) > tol_c)) then
                write(stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                    message//" for closed-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! test returned response function
            if (btest(request, 1)) &
                test_passed = test_passed .and. test_get_response_cs_c_funptr( &
                    get_response_c_funptr, test_name, &
                    " by response function returned"//message//requested)
        end do

    end function test_evaluate_dm_cs_c_funptr

    function test_evaluate_dm_os_funptr(evaluate_dm_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating function pointer for
        ! the open-shell case for every combination of requested outputs
        !
        use otr_common, only: evaluate_dm_os_type, get_response_os_type
        use test_reference, only: tol

        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        character(len=*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle)
        real(rp), pointer :: fock_ptr(:, :, :)
        real(rp) :: energy
        procedure(get_response_os_type), pointer :: get_response_funptr
        character(len=:), allocatable :: requested
        integer(ip) :: request, error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(evaluate_dm_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function for open-shell case provided"//message// &
                " not associated with value."
            return
        end if

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating subroutine for every combination of requested
        ! outputs
        do request = 0, 3
            requested = request_label(evaluate_dm_outputs, request)
            nullify(fock_ptr)
            get_response_funptr => null()
            if (btest(request, 0)) fock_ptr => fock
            energy = 0.0_rp
            fock = 0.0_rp
            if (btest(request, 1)) then
                call evaluate_dm_funptr(dm_ao, energy, fock_ptr, get_response_funptr, &
                                        error)
            else
                call evaluate_dm_funptr(dm_ao, energy, fock_ptr, error=error)
            end if

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

            ! check Fock matrix
            if (btest(request, 0) .and. &
                any(abs(fock - evaluate_dm_factors(1) * dm_ao) > tol)) then
                write(stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                    message//" for open-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! test returned response function
            if (btest(request, 1)) &
                test_passed = test_passed .and. test_get_response_os_funptr( &
                    get_response_funptr, test_name, &
                    " by response function returned"//message//requested)
        end do

    end function test_evaluate_dm_os_funptr

    function test_evaluate_dm_os_c_funptr(evaluate_dm_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided density matrix evaluating C function pointer
        ! for the open-shell case for every combination of requested outputs
        !
        use otr_common_c_interface, only: evaluate_dm_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: evaluate_dm_c_funptr
        character(len=*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(evaluate_dm_c_type), pointer :: evaluate_dm_funptr
        real(c_rp), target :: dm_ao(n_ao, n_ao, n_particle), &
                              fock(n_ao, n_ao, n_particle)
        real(c_rp), pointer :: fock_ptr(:, :, :)
        real(c_rp) :: energy
        type(c_funptr) :: get_response_c_funptr
        character(len=:), allocatable :: requested
        integer(ip) :: request
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(evaluate_dm_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Density matrix "// &
                "evaluating function for open-shell case provided"//message// &
                " not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=evaluate_dm_funptr)

        ! generate random density matrix
        call random_number(dm_ao)

        ! call density matrix evaluating function for every combination of requested
        ! outputs
        do request = 0, 3
            requested = request_label(evaluate_dm_outputs, request)
            nullify(fock_ptr)
            get_response_c_funptr = c_null_funptr
            if (btest(request, 0)) fock_ptr => fock
            energy = 0.0_c_rp
            fock = 0.0_c_rp
            if (btest(request, 1)) then
                error = &
                    evaluate_dm_funptr(dm_ao, energy, fock_ptr, get_response_c_funptr)
            else
                error = evaluate_dm_funptr(dm_ao, energy, fock_ptr)
            end if

            ! check for error
            if (error /= 0) then
                write(stderr, *) "test_"//test_name//" failed: Error produced"// &
                    message//" for open-shell case"//requested//"."
                test_passed = .false.
                cycle
            end if

            ! check energy
            if (abs(energy - sum(dm_ao)) > tol_c) then
                write(stderr, *) "test_"//test_name//" failed: Energy returned"// &
                    message//" for open-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! check Fock matrix
            if (btest(request, 0) .and. &
                any(abs(fock - evaluate_dm_factors(1) * dm_ao) > tol_c)) then
                write(stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                    message//" for open-shell case wrong"//requested//"."
                test_passed = .false.
            end if

            ! test returned response function
            if (btest(request, 1)) &
                test_passed = test_passed .and. test_get_response_os_c_funptr( &
                    get_response_c_funptr, test_name, &
                    " by response function returned"//message//requested)
        end do

    end function test_evaluate_dm_os_c_funptr

    function test_get_response_cs_funptr(get_response_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided response function pointer for the closed-shell
        ! case
        !
        use otr_common, only: get_response_cs_type
        use test_reference, only: tol

        procedure(get_response_cs_type), intent(in), pointer :: get_response_funptr
        character(len=*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :), response(:, :)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(get_response_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Response function for "// &
                "closed-shell case provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), response(n_ao, n_ao))

        ! generate random density matrix
        call random_number(dm_ao)

        ! call response subroutine
        call get_response_funptr(dm_ao, response, error)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for closed-shell case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - evaluate_dm_factors(2) * dm_ao) > tol)) then
            write(stderr, *) "test_"//test_name//" failed: Response returned"// &
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
        use otr_common_c_interface, only: get_response_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_response_c_funptr
        character(len=*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(get_response_c_type), pointer :: get_response_funptr_c
        real(c_rp), allocatable :: dm_ao(:, :), response(:, :)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(get_response_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Response function for "// &
                "closed-shell case provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_response_c_funptr, fptr=get_response_funptr_c)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao), response(n_ao, n_ao))

        ! generate random density matrix
        call random_number(dm_ao)

        ! call response function
        error = get_response_funptr_c(dm_ao, response)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for closed-shell case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - evaluate_dm_factors(2) * dm_ao) > tol_c)) then
            write(stderr, *) "test_"//test_name//" failed: Response returned"// &
                message//" for closed-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, response)

    end function test_get_response_cs_c_funptr

    function test_get_response_os_funptr(get_response_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided response function pointer for the open-shell
        ! case
        !
        use otr_common, only: get_response_os_type
        use test_reference, only: tol

        procedure(get_response_os_type), intent(in), pointer :: get_response_funptr
        character(len=*), intent(in) :: test_name, message
        logical :: test_passed

        real(rp), allocatable :: dm_ao(:, :, :), response(:, :, :)
        integer(ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. associated(get_response_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Response function for "// &
                "open-shell case provided"//message//" not associated with value."
            return
        end if

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), response(n_ao, n_ao, n_particle))

        ! generate random density matrix
        call random_number(dm_ao)

        ! call response subroutine
        call get_response_funptr(dm_ao, response, error)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for open-shell case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - evaluate_dm_factors(2) * dm_ao) > tol)) then
            write(stderr, *) "test_"//test_name//" failed: Response returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, response)

    end function test_get_response_os_funptr

    function test_get_response_os_c_funptr(get_response_c_funptr, test_name, message) &
        result(test_passed)
        !
        ! this function tests a provided response C function pointer for the open-shell
        ! case
        !
        use otr_common_c_interface, only: get_response_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_response_c_funptr
        character(len=*), intent(in) :: test_name, message
        logical :: test_passed

        procedure(get_response_c_type), pointer :: get_response_funptr_c
        real(c_rp), allocatable :: dm_ao(:, :, :), response(:, :, :)
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! check if function pointer is associated
        if (.not. c_associated(get_response_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Response function for "// &
                "open-shell case provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_response_c_funptr, fptr=get_response_funptr_c)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), response(n_ao, n_ao, n_particle))

        ! generate random density matrix
        call random_number(dm_ao)

        ! call response function
        error = get_response_funptr_c(dm_ao, response)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for open-shell case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - evaluate_dm_factors(2) * dm_ao) > tol_c)) then
            write(stderr, *) "test_"//test_name//" failed: Response returned"// &
                message//" for open-shell case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, response)

    end function test_get_response_os_c_funptr

    function request_label(outputs, request) result(label)
        !
        ! this function describes which optional outputs of a density matrix evaluating
        ! function a test requests besides the energy, given as the bits of an integer
        ! in the order of the outputs, for the failure messages of tests which go
        ! through every combination of them
        !
        character(len=*), intent(in) :: outputs(:)
        integer(ip), intent(in) :: request
        character(len=:), allocatable :: label

        integer(ip) :: i, n_requested, n_listed

        ! list the energy and every requested output
        n_requested = count([(btest(request, i - 1), i=1, size(outputs))])
        n_listed = 0
        label = " when requesting the energy"
        do i = 1, size(outputs)
            if (.not. btest(request, i - 1)) cycle
            n_listed = n_listed + 1
            if (n_listed == n_requested) then
                label = label//" and the "//trim(outputs(i))
            else
                label = label//", the "//trim(outputs(i))
            end if
        end do

    end function request_label

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

    logical function equal_orbital(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare the settings
        ! shared by the extensions parameterizing the orbitals in an orbital basis
        !
        use otr_common, only: orbital_settings_type

        type(orbital_settings_type), intent(in) :: lhs, rhs

        equal_orbital = lhs%verbose == rhs%verbose

    end function equal_orbital

end module otr_common_test_reference
