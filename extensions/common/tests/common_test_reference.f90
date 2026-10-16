! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_test_reference

    use opentrustregion, only: ip, rp, stderr
    use c_interface, only: c_ip, c_rp
    use otr_common, only: orbital_settings_type, evaluate_dm_missing_warning_msg
    use test_reference, only: same_callbacks
    use, intrinsic :: iso_c_binding, only: c_ptr, c_funptr, c_f_procpointer, &
                                           c_associated, c_null_funptr, c_char, &
                                           c_null_char

    implicit none

    ! numbers of particle channels, number of AOs and number of the occupied orbitals
    ! of every channel, which differ between the channels
    integer(ip), parameter :: n_particle = 2_ip, n_ao = 5_ip, &
                              n_occ(n_particle) = [2_ip, 1_ip]
    integer(c_ip), protected, bind(C, name="test_n_ao") :: n_ao_c = int(n_ao, kind=c_ip)

    ! reference values of the settings shared by the extensions parameterizing the
    ! orbitals in an orbital basis, every field but the logical, the logging function
    ! and the host context differing from its default
    type(orbital_settings_type), parameter :: ref_orbital_settings = &
        orbital_settings_type(logger=null(), initialized=.true., verbose=3)

    ! names of the closed-shell and the open-shell case, indexed by the number of
    ! particle channels, for the failure messages of tests covering both
    character(len=12), parameter :: shell_names(2) = &
        [character(len=12) :: "closed-shell", "open-shell"]

    ! multiples of the density matrix the mock density matrix evaluating functions
    ! return for each optional output, in the order of evaluate_dm_outputs
    real(c_rp), protected, bind(C, name="test_evaluate_dm_factors") :: &
        evaluate_dm_factors(2) = [2.0_c_rp, 3.0_c_rp]

    ! optional outputs of the density matrix evaluating function in the order of its
    ! argument list
    character(len=17), parameter :: evaluate_dm_outputs(2) = &
        [character(len=17) :: "Fock matrix", "response function"]

    ! warning the factories of the orbital bases print without a density matrix
    ! evaluating function, as a null-terminated string for the C integration tests
    integer(ip), parameter :: evaluate_dm_warning_len = &
        len(evaluate_dm_missing_warning_msg) + 1
    character(kind=c_char), protected, &
        bind(C, name="test_evaluate_dm_missing_warning") :: &
        evaluate_dm_missing_warning_c(evaluate_dm_warning_len) = transfer( &
            evaluate_dm_missing_warning_msg//c_null_char, "a", evaluate_dm_warning_len)

    interface operator(==)
        module procedure equal_orbital
    end interface

    ! the comparison of the callback functions of settings and callback bundles,
    ! extended to the host records of the levels of extension objects
    interface same_callbacks
        module procedure same_host_record_callbacks
    end interface

contains

    function ao_overlap_pattern(n_rows) result(pattern)
        !
        ! this function returns an AO overlap matrix with the given number of AOs whose
        ! values encode the pair of its AO indices (counted from one) regardless of
        ! their order, so that it is symmetric, as the Python interface expects, while
        ! its values reveal a matrix of the wrong shape
        !
        integer(c_ip), intent(in) :: n_rows
        real(c_rp) :: pattern(n_rows, n_rows)

        integer(c_ip) :: i, j

        do j = 1, n_rows
            do i = 1, n_rows
                pattern(i, j) = real(10 * max(i, j) + min(i, j), kind=c_rp)
            end do
        end do

    end function ao_overlap_pattern

    subroutine get_toy_model(h_c, s_c, g_c, nonlinear_factor_c, mo_coeff_c, n_occ_c) &
        bind(C)
        !
        ! this subroutine returns the toy model which the integration tests of the
        ! orbital bases and ARH solve through the C and Python interfaces, a mean-field
        ! energy of the shared numbers of AOs and occupied orbitals in a non-orthogonal
        ! AO basis: the core Hamiltonian, the AO overlap matrix, the positive
        ! semidefinite two-electron kernel g_pqrs = sum_k M^k_pq M^k_rs of a few
        ! symmetric matrices M^k, the factor of the non-linear term and starting
        ! orbitals orthonormal in the AO metric, the inverse square root of the overlap
        ! matrix rotated so that they are not stationary, all deterministic so that C
        ! and Python solve the same model
        !
        use opentrustregion_unit_tests, only: identity_matrix

        real(c_rp), intent(out) :: h_c(n_ao, n_ao), s_c(n_ao, n_ao), &
                                   g_c(n_ao, n_ao, n_ao, n_ao), nonlinear_factor_c, &
                                   mo_coeff_c(n_ao, n_ao)
        integer(c_ip), intent(out) :: n_occ_c(n_particle)

        integer(ip), parameter :: &
            n_kernel = 3, n_rotations = 3, &
            rotation_pairs(2, n_rotations) = reshape([1, 4, &
                                                      2, 5, &
                                                      1, 3], [2, n_rotations])
        real(rp), parameter :: rotation_angles(n_rotations) = [0.6_rp, -0.6_rp, 0.3_rp]

        real(rp) :: s(n_ao, n_ao), m(n_ao, n_ao, n_kernel), eigvals(n_ao), &
                    eigvecs(n_ao, n_ao), rotation(n_ao, n_ao), givens(n_ao, n_ao)
        real(rp), allocatable :: work(:)
        integer(ip) :: p, q, k, i, lwork, info
        external :: dsyev

        ! core Hamiltonian with diagonal elements -2 to 2 coupled by small off-diagonal
        ! elements, and AO overlap matrix close to the identity
        do q = 1, n_ao
            do p = 1, n_ao
                if (p == q) then
                    h_c(p, q) = real(p - 3, kind=c_rp)
                    s(p, q) = 1.0_rp
                else
                    h_c(p, q) = real(0.2_rp / real(p + q, kind=rp), kind=c_rp)
                    s(p, q) = 0.1_rp / real(p + q, kind=rp)
                end if
            end do
        end do
        s_c = real(s, kind=c_rp)

        ! symmetric matrices M^k and the two-electron kernel they span, a Gram matrix
        ! of AO pairs and therefore positive semidefinite
        do k = 1, n_kernel
            do q = 1, n_ao
                do p = 1, n_ao
                    m(p, q, k) = 0.4_rp * cos(real(k * p * q, kind=rp)) / &
                                 real(k, kind=rp)
                end do
            end do
        end do
        g_c = real(reshape(matmul(reshape(m, [n_ao**2, n_kernel]), &
                                  transpose(reshape(m, [n_ao**2, n_kernel]))), &
                           [n_ao, n_ao, n_ao, n_ao]), kind=c_rp)

        ! factor of the non-linear term
        nonlinear_factor_c = 0.05_c_rp

        ! inverse square root of the overlap matrix
        eigvecs = s
        allocate(work(1))
        call dsyev("V", "U", n_ao, eigvecs, n_ao, eigvals, work, -1_ip, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))
        call dsyev("V", "U", n_ao, eigvecs, n_ao, eigvals, work, lwork, info)
        deallocate(work)

        ! rotate it by plane rotations of a few pairs of orbitals
        rotation = identity_matrix(n_ao)
        do i = 1, n_rotations
            p = rotation_pairs(1, i)
            q = rotation_pairs(2, i)
            givens = identity_matrix(n_ao)
            givens(p, p) = cos(rotation_angles(i))
            givens(q, q) = cos(rotation_angles(i))
            givens(p, q) = -sin(rotation_angles(i))
            givens(q, p) = sin(rotation_angles(i))
            rotation = matmul(rotation, givens)
        end do
        mo_coeff_c = real(matmul(matmul(eigvecs * spread( &
            1.0_rp / sqrt(eigvals), 1, n_ao), transpose(eigvecs)), rotation), kind=c_rp)

        ! occupations of the particle channels
        n_occ_c = int(n_occ, kind=c_ip)

    end subroutine get_toy_model

    function check_evaluate_dm_cs_funptr(evaluate_dm_funptr, test_name, message, &
                                         context) result(test_passed)
        !
        ! this function tests a provided density matrix evaluating function pointer for
        ! the closed-shell case for every combination of requested outputs
        !
        use otr_common, only: evaluate_dm_cs_type, get_response_cs_type
        use test_reference, only: tol

        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
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
                                        error, context)
            else
                call evaluate_dm_funptr(dm_ao, energy, fock_ptr, error=error, &
                                        context=context)
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
            if (btest(request, 1)) then
                if (.not. check_get_response_cs_funptr( &
                    get_response_funptr, test_name, " by response function returned"// &
                    message//requested, context)) test_passed = .false.
            end if
        end do

    end function check_evaluate_dm_cs_funptr

    function check_evaluate_dm_os_funptr(evaluate_dm_funptr, test_name, message, &
                                         context) result(test_passed)
        !
        ! this function tests a provided density matrix evaluating function pointer for
        ! the open-shell case for every combination of requested outputs
        !
        use otr_common, only: evaluate_dm_os_type, get_response_os_type
        use test_reference, only: tol

        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
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
                                        error, context)
            else
                call evaluate_dm_funptr(dm_ao, energy, fock_ptr, error=error, &
                                        context=context)
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
            if (btest(request, 1)) then
                if (.not. check_get_response_os_funptr( &
                    get_response_funptr, test_name, " by response function returned"// &
                    message//requested, context)) test_passed = .false.
            end if
        end do

    end function check_evaluate_dm_os_funptr

    function check_evaluate_dm_c_funptr(evaluate_dm_c_funptr, n_particle, test_name, &
                                        message, context_c) result(test_passed)
        !
        ! this function tests a provided density matrix evaluating C function pointer
        ! for the closed- (n_particle = 1) or open-shell (n_particle = 2) case for every
        ! combination of requested outputs
        !
        use otr_common_c_interface, only: evaluate_dm_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: evaluate_dm_c_funptr
        integer(ip), intent(in) :: n_particle
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(evaluate_dm_c_type), pointer :: evaluate_dm_funptr
        real(c_rp), target :: dm_ao(n_ao, n_ao, n_particle), &
                              fock(n_ao, n_ao, n_particle)
        real(c_rp), pointer :: fock_ptr(:, :, :)
        real(c_rp) :: energy
        type(c_funptr) :: get_response_c_funptr
        character(len=:), allocatable :: requested, shell
        integer(ip) :: request
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! name of the shell for the messages
        shell = trim(shell_names(n_particle))

        ! check if function pointer is associated
        if (.not. c_associated(evaluate_dm_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name// &
                " failed: Density matrix evaluating function for "//shell// &
                " case provided"//message//" not associated with value."
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
                error = evaluate_dm_funptr(dm_ao, energy, fock_ptr, &
                                           get_response_c_funptr, context_c)
            else
                error = evaluate_dm_funptr(dm_ao, energy, fock_ptr, context_c=context_c)
            end if

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

            ! check Fock matrix
            if (btest(request, 0) .and. &
                any(abs(fock - evaluate_dm_factors(1) * dm_ao) > tol_c)) then
                write(stderr, *) "test_"//test_name//" failed: Fock matrix returned"// &
                    message//" for "//shell//" case wrong"//requested//"."
                test_passed = .false.
            end if

            ! test returned response function
            if (btest(request, 1)) then
                if (.not. check_get_response_c_funptr( &
                    get_response_c_funptr, n_particle, test_name, &
                    " by response function returned"//message//requested, context_c)) &
                    test_passed = .false.
            end if
        end do

    end function check_evaluate_dm_c_funptr

    function check_get_response_cs_funptr(get_response_funptr, test_name, message, &
                                          context) result(test_passed)
        !
        ! this function tests a provided response function pointer for the closed-shell
        ! case
        !
        use otr_common, only: get_response_cs_type
        use test_reference, only: tol

        procedure(get_response_cs_type), intent(in), pointer :: get_response_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
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
        call get_response_funptr(dm_ao, response, error, context)

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

    end function check_get_response_cs_funptr

    function check_get_response_os_funptr(get_response_funptr, test_name, message, &
                                          context) result(test_passed)
        !
        ! this function tests a provided response function pointer for the open-shell
        ! case
        !
        use otr_common, only: get_response_os_type
        use test_reference, only: tol

        procedure(get_response_os_type), intent(in), pointer :: get_response_funptr
        character(len=*), intent(in) :: test_name, message
        class(*), intent(in), pointer :: context
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
        call get_response_funptr(dm_ao, response, error, context)

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

    end function check_get_response_os_funptr

    function check_get_response_c_funptr(get_response_c_funptr, n_particle, test_name, &
                                         message, context_c) result(test_passed)
        !
        ! this function tests a provided response C function pointer for the closed-
        ! (n_particle = 1) or open-shell (n_particle = 2) case
        !
        use otr_common_c_interface, only: get_response_c_type
        use test_reference, only: tol_c

        type(c_funptr), intent(in) :: get_response_c_funptr
        integer(ip), intent(in) :: n_particle
        character(len=*), intent(in) :: test_name, message
        type(c_ptr), intent(in) :: context_c
        logical :: test_passed

        procedure(get_response_c_type), pointer :: get_response_funptr_c
        real(c_rp), allocatable :: dm_ao(:, :, :), response(:, :, :)
        character(len=:), allocatable :: shell
        integer(c_ip) :: error

        ! assume tests pass
        test_passed = .true.

        ! name of the shell for the messages
        shell = trim(shell_names(n_particle))

        ! check if function pointer is associated
        if (.not. c_associated(get_response_c_funptr)) then
            test_passed = .false.
            write(stderr, *) "test_"//test_name//" failed: Response function for "// &
                shell//" case provided"//message//" not associated with value."
            return
        end if

        ! convert to Fortran function pointer
        call c_f_procpointer(cptr=get_response_c_funptr, fptr=get_response_funptr_c)

        ! allocate arrays
        allocate(dm_ao(n_ao, n_ao, n_particle), response(n_ao, n_ao, n_particle))

        ! generate random density matrix
        call random_number(dm_ao)

        ! call response function
        error = get_response_funptr_c(dm_ao, response, context_c)

        ! check for error
        if (error /= 0) then
            write(stderr, *) "test_"//test_name//" failed: Error produced"//message// &
                " for "//shell//" case."
            test_passed = .false.
        end if

        ! check response
        if (any(abs(response - evaluate_dm_factors(2) * dm_ao) > tol_c)) then
            write(stderr, *) "test_"//test_name//" failed: Response returned"// &
                message//" for "//shell//" case wrong."
            test_passed = .false.
        end if

        ! deallocate arrays
        deallocate(dm_ao, response)

    end function check_get_response_c_funptr

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

    logical function equal_orbital(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare the settings
        ! shared by the extensions parameterizing the orbitals in an orbital basis
        !
        type(orbital_settings_type), intent(in) :: lhs, rhs

        equal_orbital = (lhs%initialized .eqv. rhs%initialized) .and. &
                        lhs%verbose == rhs%verbose

    end function equal_orbital

    logical function same_host_record_callbacks(host, expected)
        !
        ! this function checks that every callback function of the host record of a
        ! level of an extension object is that of the expected record, unset where this
        ! has none
        !
        use otr_common, only: host_record_type

        type(host_record_type), intent(in) :: host, expected

        same_host_record_callbacks = &
            (associated(host%update_orbs, expected%update_orbs) .or. &
             .not. (associated(host%update_orbs) .or. &
                    associated(expected%update_orbs))) .and. &
            (associated(host%obj_func, expected%obj_func) .or. &
             .not. (associated(host%obj_func) .or. &
                    associated(expected%obj_func))) .and. &
            (associated(host%hess_x, expected%hess_x) .or. &
             .not. (associated(host%hess_x) .or. associated(expected%hess_x))) .and. &
            (associated(host%precond, expected%precond) .or. &
             .not. (associated(host%precond) .or. associated(expected%precond))) .and. &
            (associated(host%precond_pd, expected%precond_pd) .or. &
             .not. (associated(host%precond_pd) .or. &
                    associated(expected%precond_pd))) .and. &
            (associated(host%project, expected%project) .or. &
             .not. (associated(host%project) .or. associated(expected%project))) .and. &
            (associated(host%get_extra_trial_vectors, &
                        expected%get_extra_trial_vectors) .or. &
             .not. (associated(host%get_extra_trial_vectors) .or. &
                    associated(expected%get_extra_trial_vectors))) .and. &
            (associated(host%conv_check, expected%conv_check) .or. &
             .not. (associated(host%conv_check) .or. &
                    associated(expected%conv_check))) .and. &
            (associated(host%stability_hess_x, expected%stability_hess_x) .or. &
             .not. (associated(host%stability_hess_x) .or. &
                    associated(expected%stability_hess_x))) .and. &
            (associated(host%approx_hess_x, expected%approx_hess_x) .or. &
             .not. (associated(host%approx_hess_x) .or. &
                    associated(expected%approx_hess_x))) .and. &
            (associated(host%conv_check_stability, expected%conv_check_stability) .or. &
             .not. (associated(host%conv_check_stability) .or. &
                    associated(expected%conv_check_stability))) .and. &
            (associated(host%logger, expected%logger) .or. &
             .not. (associated(host%logger) .or. associated(expected%logger)))

    end function same_host_record_callbacks

    logical function same_logger_and_context(settings, expected)
        !
        ! this function checks that the logging function and the host context of the
        ! settings of an extension are those of the expected settings, unset where these
        ! have none, which the comparisons of the settings values leave out
        !
        use opentrustregion, only: settings_type

        class(settings_type), intent(in) :: settings, expected

        same_logger_and_context = &
            (associated(settings%logger, expected%logger) .or. &
             .not. (associated(settings%logger) .or. &
                    associated(expected%logger))) .and. &
            (associated(settings%context, expected%context) .or. &
             .not. (associated(settings%context) .or. associated(expected%context)))

    end function same_logger_and_context

end module otr_common_test_reference
