! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! global number of calls to the mock density matrix updating functions, so that
    ! tests can produce potential matrices which change between calls and thereby
    ! distinguish a cached quantity from one that was incorrectly recomputed from an
    ! unchanged input
    integer(ip) :: n_mock_calls = 0

    ! multiplier of the density matrix returned by the mock density matrix updating
    ! functions, which differs between the first and any subsequent call, further
    ! multipliers for extension-specific potentials follow the same convention
    real(rp), parameter :: mock_fock_factor(2) = [2.0_rp, 5.0_rp]

contains

    function identity_matrix(n) result(matrix)
        !
        ! this function returns the identity matrix
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        integer(ip) :: i

        matrix = 0.0_rp
        do i = 1, n
            matrix(i, i) = 1.0_rp
        end do

    end function identity_matrix

    function generate_random_symm_matrix(n) result(matrix)
        !
        ! this function generates a random symmetric matrix, corresponding to a valid
        ! density matrix displacement
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        call random_number(matrix)
        matrix = matrix + transpose(matrix)

    end function generate_random_symm_matrix

    function generate_random_density_matrix(n, n_electrons) result(dm)
        !
        ! this function generates a random valid density matrix by first generating a
        ! random set of orthonormal basis vectors and then summing the outer products
        ! of these
        !
        use opentrustregion, only: numerical_zero

        integer(ip), intent(in) :: n, n_electrons
        real(rp) :: dm(n, n)

        real(rp) :: vec(n), val, u(n, n_electrons)
        integer(ip) :: i, j

        ! generate random orthonormal basis
        do j = 1, n_electrons
            do
                call random_number(vec)
                do i = 1, j - 1
                    vec = vec - sum(vec * u(:, i)) * u(:, i)
                end do
                val = sqrt(sum(vec**2))
                if (val >= numerical_zero) then
                    u(:, j) = vec / val
                    exit
                end if
            end do
        end do

        ! construct valid density matrix by summing outer products of basis vectors
        dm = matmul(u, transpose(u))

    end function generate_random_density_matrix

    function mock_factor(factors) result(factor)
        !
        ! this function returns the multiplier the mock density matrix updating
        ! functions apply on the current call
        !
        real(rp), intent(in) :: factors(:)
        real(rp) :: factor

        factor = factors(min(n_mock_calls, size(factors)))

    end function mock_factor

    function ref_unpack_asymm(matrix_nonred, n_particle_in, n_ao_in) result(matrix)
        !
        ! this function unpacks an antisymmetric matrix, reproducing the corresponding
        ! OAO routine so that tests of routines which unpack internally do not depend
        ! on it
        !
        real(rp), intent(in) :: matrix_nonred(:)
        integer(ip), intent(in) :: n_particle_in, n_ao_in
        real(rp) :: matrix(n_ao_in, n_ao_in, n_particle_in)

        integer(ip) :: i, j, k, idx

        matrix = 0.0_rp
        idx = 1
        do k = 1, n_particle_in
            do j = 1, n_ao_in
                do i = 1, j - 1
                    matrix(i, j, k) = matrix_nonred(idx)
                    matrix(j, i, k) = -matrix_nonred(idx)
                    idx = idx + 1
                end do
            end do
        end do

    end function ref_unpack_asymm

    function ref_pack_asymm(matrix, n_param) result(matrix_nonred)
        !
        ! this function packs an antisymmetric matrix, reproducing the corresponding
        ! OAO routine so that tests of routines which pack internally do not depend on
        ! it
        !
        real(rp), intent(in) :: matrix(:, :, :)
        integer(ip), intent(in) :: n_param
        real(rp) :: matrix_nonred(n_param)

        integer(ip) :: i, j, k, idx

        idx = 1
        do k = 1, size(matrix, 3)
            do j = 1, size(matrix, 2)
                do i = 1, j - 1
                    matrix_nonred(idx) = matrix(i, j, k)
                    idx = idx + 1
                end do
            end do
        end do

    end function ref_pack_asymm

    function ref_project_asymm(matrix, dm_oao) result(projected_matrix)
        !
        ! this function retains only the occupied-virtual and virtual-occupied
        ! contributions to a matrix in antisymmetric form, reproducing the
        ! corresponding OAO routine so that tests of routines which project internally 
        ! do not depend on it
        !
        real(rp), intent(in) :: matrix(:, :, :), dm_oao(:, :, :)
        real(rp) :: projected_matrix(size(matrix, 1), size(matrix, 2), size(matrix, 3))

        integer(ip) :: i
        real(rp) :: proj_v(size(matrix, 1), size(matrix, 1))

        do i = 1, size(matrix, 3)
            ! construct projection matrix on virtual space
            proj_v = identity_matrix(size(matrix, 1, kind=ip)) - dm_oao(:, :, i)

            ! construct virtual-occupied and occupied-virtual contributions
            projected_matrix(:, :, i) = matmul(dm_oao(:, :, i), &
                                               matmul(matrix(:, :, i), proj_v))
            projected_matrix(:, :, i) = projected_matrix(:, :, i) - &
                                        transpose(projected_matrix(:, :, i))
        end do

    end function ref_project_asymm

    function ref_project_symm(x_full, dm_oao) result(delta_dm)
        !
        ! this function retains only the occupied-virtual and virtual-occupied
        ! contributions to a matrix in symmetric form, reproducing the corresponding 
        ! OAO routine so that tests of routines which project internally do not depend 
        ! on it
        !
        real(rp), intent(in) :: x_full(:, :, :), dm_oao(:, :, :)
        real(rp) :: delta_dm(size(x_full, 1), size(x_full, 2), size(x_full, 3))

        integer(ip) :: i

        do i = 1, size(x_full, 3)
            delta_dm(:, :, i) = matmul(dm_oao(:, :, i), x_full(:, :, i))
            delta_dm(:, :, i) = delta_dm(:, :, i) + transpose(delta_dm(:, :, i))
        end do

    end function ref_project_symm

    function ref_hess_x(x_full, response, dm_oao, fock_oo, fock_vv, n_param) &
        result(hess_x)
        !
        ! this function assembles a Hessian linear transformation in the OAO basis from
        ! a given response contribution
        !
        real(rp), intent(in) :: x_full(:, :, :), response(:, :, :), dm_oao(:, :, :), &
                                fock_oo(:, :, :), fock_vv(:, :, :)
        integer(ip), intent(in) :: n_param
        real(rp) :: hess_x(n_param)

        integer(ip) :: i
        real(rp) :: hess_x_full(size(x_full, 1), size(x_full, 2), size(x_full, 3))

        ! get static part
        do i = 1, size(x_full, 3)
            hess_x_full(:, :, i) = matmul(fock_vv(:, :, i) - fock_oo(:, :, i), &
                                          x_full(:, :, i))
            hess_x_full(:, :, i) = hess_x_full(:, :, i) - &
                                   transpose(hess_x_full(:, :, i))
        end do

        ! project the combined static and response contributions onto the
        ! occupied-virtual and virtual-occupied subspace, matching the production code
        hess_x_full = ref_project_asymm(hess_x_full + response, dm_oao)

        ! pack Hessian linear transformation
        hess_x = merge(4.0_rp, 2.0_rp, size(x_full, 3) == 1) * &
                 ref_pack_asymm(hess_x_full, n_param)

    end function ref_hess_x

    subroutine ref_diagonalize_static_part(fock_oo, fock_vv, eigvecs, eigvals)
        !
        ! this subroutine independently diagonalizes the static part of the Hessian
        ! for each particle channel
        !
        real(rp), intent(in) :: fock_oo(:, :, :), fock_vv(:, :, :)
        real(rp), intent(out) :: eigvecs(:, :, :), eigvals(:, :)

        integer(ip) :: n_ao, n_particle, i, lwork, info
        real(rp), allocatable :: a(:, :), work(:)
        external :: dsyev

        n_ao = size(fock_oo, 1)
        n_particle = size(fock_oo, 3)

        allocate(a(n_ao, n_ao))
        do i = 1, n_particle
            a = fock_vv(:, :, i) - fock_oo(:, :, i)
            allocate(work(1))
            call dsyev("V", "U", n_ao, a, n_ao, eigvals(:, i), work, -1_ip, info)
            lwork = int(work(1))
            deallocate(work)
            allocate(work(lwork))
            call dsyev("V", "U", n_ao, a, n_ao, eigvals(:, i), work, lwork, info)
            deallocate(work)
            eigvecs(:, :, i) = a
        end do
        deallocate(a)

    end subroutine ref_diagonalize_static_part

    function ref_rotate_to_eigenbasis(vector, eigvecs, n_particle, n_ao) &
        result(rotated)
        !
        ! this function independently rotates a packed antisymmetric vector into a
        ! given eigenbasis, one particle channel at a time, reproducing the 
        ! corresponding OAO routine so that tests of routines which project internally 
        ! do not depend on it
        !
        real(rp), intent(in) :: vector(:), eigvecs(:, :, :)
        integer(ip), intent(in) :: n_particle, n_ao
        real(rp) :: rotated(size(vector))

        real(rp) :: full(n_ao, n_ao, n_particle), rotated_full(n_ao, n_ao, n_particle)
        integer(ip) :: i

        full = ref_unpack_asymm(vector, n_particle, n_ao)
        do i = 1, n_particle
            rotated_full(:, :, i) = matmul(transpose(eigvecs(:, :, i)), &
                                           matmul(full(:, :, i), eigvecs(:, :, i)))
        end do

        rotated = ref_pack_asymm(rotated_full, size(vector, kind=ip))

    end function ref_rotate_to_eigenbasis

    function ref_rotate_from_eigenbasis(vector, eigvecs, n_particle, n_ao) &
        result(rotated)
        !
        ! this function independently rotates a packed antisymmetric vector out of a
        ! given eigenbasis, one particle channel at a time, reproducing the 
        ! corresponding OAO routine so that tests of routines which project internally 
        ! do not depend on it
        !
        real(rp), intent(in) :: vector(:), eigvecs(:, :, :)
        integer(ip), intent(in) :: n_particle, n_ao
        real(rp) :: rotated(size(vector))

        real(rp) :: full(n_ao, n_ao, n_particle), rotated_full(n_ao, n_ao, n_particle)
        integer(ip) :: i

        full = ref_unpack_asymm(vector, n_particle, n_ao)
        do i = 1, n_particle
            rotated_full(:, :, i) = matmul(eigvecs(:, :, i), &
                                           matmul(full(:, :, i), &
                                                  transpose(eigvecs(:, :, i))))
        end do

        rotated = ref_pack_asymm(rotated_full, size(vector, kind=ip))

    end function ref_rotate_from_eigenbasis

    function ref_hess_eigval_pairs(eigvals, n_particle, n_ao, n_param) &
        result(eigval_pairs)
        !
        ! this function independently constructs the pairwise sums of the static
        ! Hessian part eigenvalues, reproducing the corresponding OAO routine so that 
        ! tests of routines which project internally do not depend on it
        !
        real(rp), intent(in) :: eigvals(:, :)
        integer(ip), intent(in) :: n_particle, n_ao, n_param
        real(rp) :: eigval_pairs(n_param)

        integer(ip) :: i, j, k, idx

        idx = 1
        do k = 1, n_particle
            do j = 1, n_ao
                do i = 1, j - 1
                    eigval_pairs(idx) = eigvals(i, k) + eigvals(j, k)
                    idx = idx + 1
                end do
            end do
        end do
        eigval_pairs = merge(4.0_rp, 2.0_rp, n_particle == 1) * eigval_pairs

    end function ref_hess_eigval_pairs

    function mock_get_energy_cs(dm, error) result(energy)
        !
        ! this function is a mock energy function for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        error = 0
        energy = sum(dm)

    end function mock_get_energy_cs

    function mock_get_energy_os(dm, error) result(energy)
        !
        ! this function is a mock energy function for the open-shell case
        !
        real(rp), intent(in), target :: dm(:, :, :)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        error = 0
        energy = sum(dm)

    end function mock_get_energy_os

    subroutine mock_get_response_cs(dm, response, error)
        !
        ! this subroutine is a mock response function for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out), target, contiguous :: response(:, :)
        integer(ip), intent(out) :: error

        error = 0
        response = 2.0_rp * dm

    end subroutine mock_get_response_cs

    subroutine mock_get_response_os(dm, response, error)
        !
        ! this subroutine is a mock response function for the open-shell case
        !
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out), target :: response(:, :, :)
        integer(ip), intent(out) :: error

        error = 0
        response = 2.0_rp * dm

    end subroutine mock_get_response_os

    subroutine mock_update_dm_cs(dm, energy, fock, get_response_funptr, error)
        !
        ! this subroutine is a mock density matrix updating function for the
        ! closed-shell case, which returns a multiple of the density matrix that
        ! changes between calls so that non-vanishing differences are produced
        !
        use otr_oao, only: get_response_cs_type

        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target, contiguous :: fock(:, :)
        procedure(get_response_cs_type), intent(out), pointer :: get_response_funptr
        integer(ip), intent(out) :: error

        n_mock_calls = n_mock_calls + 1

        error = 0
        energy = sum(dm)
        fock = mock_factor(mock_fock_factor) * dm
        get_response_funptr => mock_get_response_cs

    end subroutine mock_update_dm_cs

    subroutine mock_update_dm_os(dm, energy, fock, get_response_funptr, error)
        !
        ! this subroutine is a mock density matrix updating function for the open-shell
        ! case, which returns a multiple of the density matrix that changes between
        ! calls so that non-vanishing differences are produced
        !
        use otr_oao, only: get_response_os_type

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target :: fock(:, :, :)
        procedure(get_response_os_type), intent(out), pointer :: get_response_funptr
        integer(ip), intent(out) :: error

        n_mock_calls = n_mock_calls + 1

        error = 0
        energy = sum(dm)
        fock = mock_factor(mock_fock_factor) * dm
        get_response_funptr => mock_get_response_os

    end subroutine mock_update_dm_os

    logical(c_bool) function test_init_oao_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the OAO settings
        !
        use otr_oao, only: oao_settings_type, default_settings => default_oao_settings
        use otr_oao_test_reference, only: operator(==)

        type(oao_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_init_oao_settings = .true.

        ! initialize settings
        call settings%init(error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_init_oao_settings failed: Function raised error."
            test_init_oao_settings = .false.
        end if

        ! check settings
        if (.not. (settings == default_settings)) then
            write (stderr, *) "test_init_oao_settings failed: Settings not "// &
                "initialized correctly."
            test_init_oao_settings = .false.
        end if

    end function test_init_oao_settings

    logical(c_bool) function test_oao_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the OAO
        ! input parameters
        !
        use otr_oao, only: oao_settings_type, oao_sanity_check
        use opentrustregion_unit_tests, only: setup_settings

        type(oao_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_oao_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if a positive number of AOs is accepted
        call oao_sanity_check(settings, 1_ip, error)
        if (error /= 0) then
            write (stderr, *) "test_oao_sanity_check failed: Error thrown for "// &
                "positive number of AOs."
            test_oao_sanity_check = .false.
        end if

        ! check if a vanishing number of AOs is rejected
        call oao_sanity_check(settings, 0_ip, error)
        if (error == 0) then
            write (stderr, *) "test_oao_sanity_check failed: Error not thrown for "// &
                "vanishing number of AOs."
            test_oao_sanity_check = .false.
        end if

    end function test_oao_sanity_check

    logical(c_bool) function test_oao_deconstructor() bind(C)
        !
        ! this function tests the subroutine which deallocates the OAO objects
        !
        use otr_oao, only: oao_deconstructor, oao_object

        ! assume tests pass
        test_oao_deconstructor = .true.

        ! allocate OAO object
        if (.not. allocated(oao_object)) allocate(oao_object)

        ! call routine and determine if the object is deallocated
        call oao_deconstructor()
        if (allocated(oao_object)) then
            write (stderr, *) "test_oao_deconstructor failed: Object not deallocated."
            test_oao_deconstructor = .false.
        end if

        ! call routine again and determine if an already deallocated object is handled
        call oao_deconstructor()
        if (allocated(oao_object)) then
            write (stderr, *) "test_oao_deconstructor failed: Already deallocated "// &
                "object not handled."
            test_oao_deconstructor = .false.
        end if

    end function test_oao_deconstructor

    logical(c_bool) function test_unpack_asymm() bind(C)
        !
        ! this function tests the function which unpacks an antisymmetric matrix
        !
        use otr_oao, only: unpack_asymm
        use otr_oao_test_reference, only: n_ao, n_particle

        real(rp) :: expected(n_ao, n_ao, n_particle)
        real(rp), allocatable :: matrix(:, :, :)

        ! assume tests pass
        test_unpack_asymm = .true.

        ! initialize expected matrices
        expected(:, :, 1) = reshape([0.0_rp, -1.0_rp, -2.0_rp, &
                                     1.0_rp, 0.0_rp, -3.0_rp, &
                                     2.0_rp, 3.0_rp, 0.0_rp], [n_ao, n_ao])
        expected(:, :, 2) = reshape([0.0_rp, -4.0_rp, -5.0_rp, &
                                     4.0_rp, 0.0_rp, -6.0_rp, &
                                     5.0_rp, 6.0_rp, 0.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrices
        ! match
        matrix = unpack_asymm([1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp], &
                              n_particle, n_ao)
        if (size(matrix, 1) /= n_ao .or. size(matrix, 2) /= n_ao .or. &
            size(matrix, 3) /= n_particle) then
            write (stderr, *) "test_unpack_asymm failed: Incorrect matrix "// &
                "dimensions."
            test_unpack_asymm = .false.
            return
        end if
        if (norm2(matrix - expected) > tol) then
            write (stderr, *) "test_unpack_asymm failed: Incorrect matrix values."
            test_unpack_asymm = .false.
        end if
        deallocate(matrix)

    end function test_unpack_asymm

    logical(c_bool) function test_pack_asymm() bind(C)
        !
        ! this function tests the function which packs an antisymmetric matrix
        !
        use otr_oao, only: pack_asymm
        use otr_oao_test_reference, only: n_ao, n_particle, n_param

        real(rp) :: matrix(n_ao, n_ao, n_particle)
        real(rp), allocatable :: matrix_nonred(:)

        ! assume tests pass
        test_pack_asymm = .true.

        ! initialize antisymmetric matrices
        matrix(:, :, 1) = reshape([0.0_rp, -1.0_rp, -2.0_rp, &
                                   1.0_rp, 0.0_rp, -3.0_rp, &
                                   2.0_rp, 3.0_rp, 0.0_rp], [n_ao, n_ao])
        matrix(:, :, 2) = reshape([0.0_rp, -3.0_rp, -4.0_rp, &
                                   3.0_rp, 0.0_rp, -5.0_rp, &
                                   4.0_rp, 5.0_rp, 0.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and
        ! values of resulting vector match, where each spin is packed separately
        matrix_nonred = pack_asymm(matrix, n_param)
        if (size(matrix_nonred) /= n_param) then
            write (stderr, *) "test_pack_asymm failed: Incorrect vector dimension."
            test_pack_asymm = .false.
            return
        end if
        if (norm2(matrix_nonred - [1.0_rp, 2.0_rp, 3.0_rp, 3.0_rp, 4.0_rp, 5.0_rp]) > &
            tol) then
            write (stderr, *) "test_pack_asymm failed: Incorrect vector values."
            test_pack_asymm = .false.
        end if
        deallocate(matrix_nonred)

    end function test_pack_asymm

    logical(c_bool) function test_project_asymm() bind(C)
        !
        ! this function tests the function which retains only the occupied-virtual and
        ! virtual-occupied contributions to a matrix in antisymmetric form
        !
        use otr_oao, only: project_asymm
        use otr_oao_test_reference, only: n_ao, n_particle

        real(rp) :: matrix(n_ao, n_ao, n_particle), dm_oao(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle)
        real(rp), allocatable :: projected_matrix(:, :, :)

        ! assume tests pass
        test_project_asymm = .true.

        ! initialize matrices and density matrices, each occupying a single, distinct
        ! orbital only
        matrix(:, :, 1) = reshape([1.0_rp, 2.0_rp, 3.0_rp, &
                                   4.0_rp, 5.0_rp, 6.0_rp, &
                                   7.0_rp, 8.0_rp, 9.0_rp], [n_ao, n_ao])
        matrix(:, :, 2) = reshape([9.0_rp, 8.0_rp, 7.0_rp, &
                                   6.0_rp, 5.0_rp, 4.0_rp, &
                                   3.0_rp, 2.0_rp, 1.0_rp], [n_ao, n_ao])
        dm_oao = 0.0_rp
        dm_oao(1, 1, 1) = 1.0_rp
        dm_oao(2, 2, 2) = 1.0_rp

        ! initialize expected matrices, where only the antisymmetrized occupied-virtual
        ! elements survive
        expected(:, :, 1) = reshape([0.0_rp, -4.0_rp, -7.0_rp, &
                                     4.0_rp, 0.0_rp, 0.0_rp, &
                                     7.0_rp, 0.0_rp, 0.0_rp], [n_ao, n_ao])
        expected(:, :, 2) = reshape([0.0_rp, 8.0_rp, 0.0_rp, &
                                     -8.0_rp, 0.0_rp, -2.0_rp, &
                                     0.0_rp, 2.0_rp, 0.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrix match
        projected_matrix = project_asymm(matrix, dm_oao)
        if (size(projected_matrix, 1) /= n_ao .or. &
            size(projected_matrix, 3) /= n_particle) then
            write (stderr, *) "test_project_asymm failed: Incorrect matrix dimensions."
            test_project_asymm = .false.
            return
        end if
        if (norm2(projected_matrix - expected) > tol) then
            write (stderr, *) "test_project_asymm failed: Incorrect matrix values."
            test_project_asymm = .false.
        end if
        deallocate(projected_matrix)

    end function test_project_asymm

    logical(c_bool) function test_project_symm() bind(C)
        !
        ! this function tests the function which retains only the occupied-virtual and
        ! virtual-occupied contributions to a matrix in symmetric form
        !
        use otr_oao, only: project_symm
        use otr_oao_test_reference, only: n_ao, n_particle

        real(rp) :: x_full(n_ao, n_ao, n_particle), dm_oao(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle)
        real(rp), allocatable :: projected_matrix(:, :, :)

        ! assume tests pass
        test_project_symm = .true.

        ! initialize antisymmetric trial vectors and density matrices, each occupying
        ! a single, distinct orbital only
        x_full(:, :, 1) = reshape([0.0_rp, -1.0_rp, -2.0_rp, &
                                   1.0_rp, 0.0_rp, -3.0_rp, &
                                   2.0_rp, 3.0_rp, 0.0_rp], [n_ao, n_ao])
        x_full(:, :, 2) = reshape([0.0_rp, -4.0_rp, -5.0_rp, &
                                   4.0_rp, 0.0_rp, -6.0_rp, &
                                   5.0_rp, 6.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_oao = 0.0_rp
        dm_oao(1, 1, 1) = 1.0_rp
        dm_oao(2, 2, 2) = 1.0_rp

        ! initialize expected matrices, where only the symmetrized occupied-virtual
        ! elements survive
        expected(:, :, 1) = reshape([0.0_rp, 1.0_rp, 2.0_rp, &
                                     1.0_rp, 0.0_rp, 0.0_rp, &
                                     2.0_rp, 0.0_rp, 0.0_rp], [n_ao, n_ao])
        expected(:, :, 2) = reshape([0.0_rp, -4.0_rp, 0.0_rp, &
                                     -4.0_rp, 0.0_rp, 6.0_rp, &
                                     0.0_rp, 6.0_rp, 0.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrix match
        projected_matrix = project_symm(x_full, dm_oao)
        if (size(projected_matrix, 1) /= n_ao .or. &
            size(projected_matrix, 3) /= n_particle) then
            write (stderr, *) "test_project_symm failed: Incorrect matrix dimensions."
            test_project_symm = .false.
            return
        end if
        if (norm2(projected_matrix - expected) > tol) then
            write (stderr, *) "test_project_symm failed: Incorrect matrix values."
            test_project_symm = .false.
        end if
        deallocate(projected_matrix)

    end function test_project_symm

    logical(c_bool) function test_purify() bind(C)
        !
        ! this function tests the subroutine which purifies a density matrix
        !
        use otr_oao, only: purify
        use otr_oao_test_reference, only: n_ao, n_particle

        integer(ip), parameter :: n_electrons = 2
        real(rp) :: dm(n_ao, n_ao, n_particle), expected(n_ao, n_ao, n_particle), &
                    idempotent_dm(n_ao, n_ao, n_particle), &
                    purified_dm(n_ao, n_ao, n_particle)
        integer(ip) :: i

        ! assume tests pass
        test_purify = .true.

        ! initialize density matrices with fractional occupations
        dm(:, :, 1) = reshape([0.6_rp, 0.0_rp, 0.0_rp, &
                               0.0_rp, 0.2_rp, 0.0_rp, &
                               0.0_rp, 0.0_rp, 0.3_rp], [n_ao, n_ao])
        dm(:, :, 2) = reshape([0.7_rp, 0.0_rp, 0.0_rp, &
                               0.0_rp, 0.4_rp, 0.0_rp, &
                               0.0_rp, 0.0_rp, 0.5_rp], [n_ao, n_ao])

        ! initialize expected density matrices, where the occupations are driven
        ! towards zero and one
        expected(:, :, 1) = reshape([0.648_rp, 0.0_rp, 0.0_rp, &
                                     0.0_rp, 0.104_rp, 0.0_rp, &
                                     0.0_rp, 0.0_rp, 0.216_rp], [n_ao, n_ao])
        expected(:, :, 2) = reshape([0.784_rp, 0.0_rp, 0.0_rp, &
                                     0.0_rp, 0.352_rp, 0.0_rp, &
                                     0.0_rp, 0.0_rp, 0.5_rp], [n_ao, n_ao])

        ! call routine and determine if values of resulting density matrices match
        call purify(dm)
        if (norm2(dm - expected) > tol) then
            write (stderr, *) "test_purify failed: Incorrect density matrix values."
            test_purify = .false.
        end if

        ! initialize idempotent density matrices
        do i = 1, n_particle
            idempotent_dm(:, :, i) = generate_random_density_matrix(n_ao, n_electrons)
        end do
        purified_dm = idempotent_dm

        ! call routine and determine if idempotent density matrices are left unchanged
        call purify(purified_dm)
        if (norm2(purified_dm - idempotent_dm) > tol) then
            write (stderr, *) "test_purify failed: Idempotent density matrices not "// &
                "left unchanged."
            test_purify = .false.
        end if

    end function test_purify

    logical(c_bool) function test_symmetric_transformation() bind(C)
        !
        ! this function tests the function which performs a symmetric transformation
        !
        use otr_oao, only: symmetric_transformation

        integer(ip), parameter :: n_ao = 2, n_particle = 1
        real(rp) :: trans_matrix(n_ao, n_ao), matrix(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle)
        real(rp), allocatable :: matrix_transformed(:, :, :)

        ! assume tests pass
        test_symmetric_transformation = .true.

        ! initialize transformation matrix and matrix to be transformed
        trans_matrix = reshape([1.0_rp, 0.0_rp, &
                                2.0_rp, 1.0_rp], [n_ao, n_ao])
        matrix(:, :, 1) = reshape([1.0_rp, 0.0_rp, &
                                   0.0_rp, 2.0_rp], [n_ao, n_ao])

        ! initialize expected matrix
        expected(:, :, 1) = reshape([1.0_rp, 0.0_rp, &
                                     6.0_rp, 2.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrix match
        matrix_transformed = symmetric_transformation(trans_matrix, matrix)
        if (size(matrix_transformed, 1) /= n_ao .or. &
            size(matrix_transformed, 3) /= n_particle) then
            write (stderr, *) "test_symmetric_transformation failed: Incorrect "// &
                "matrix dimensions."
            test_symmetric_transformation = .false.
            return
        end if
        if (norm2(matrix_transformed - expected) > tol) then
            write (stderr, *) "test_symmetric_transformation failed: Incorrect "// &
                "matrix values."
            test_symmetric_transformation = .false.
        end if
        deallocate(matrix_transformed)

    end function test_symmetric_transformation

    logical(c_bool) function test_matrix_exponential() bind(C)
        !
        ! this function tests the function which calculates the matrix exponential of a
        ! real antisymmetric matrix
        !
        use otr_oao, only: matrix_exponential

        integer(ip), parameter :: n_ao = 2
        real(rp), parameter :: angle = 0.3_rp

        real(rp) :: a(n_ao, n_ao), expected(n_ao, n_ao)
        real(rp), allocatable :: exp_a(:, :)
        integer(ip) :: error

        ! assume tests pass
        test_matrix_exponential = .true.

        ! initialize antisymmetric matrix, whose exponential is the corresponding
        ! rotation matrix
        a = reshape([0.0_rp, -angle, &
                     angle, 0.0_rp], [n_ao, n_ao])

        ! initialize expected rotation matrix
        expected = reshape([cos(angle), -sin(angle), &
                            sin(angle), cos(angle)], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrix match
        exp_a = matrix_exponential(a, error)
        if (error /= 0) then
            write (stderr, *) "test_matrix_exponential failed: Produced error."
            test_matrix_exponential = .false.
        end if
        if (size(exp_a, 1) /= n_ao .or. size(exp_a, 2) /= n_ao) then
            write (stderr, *) "test_matrix_exponential failed: Incorrect matrix "// &
                "dimensions."
            test_matrix_exponential = .false.
            return
        end if
        if (norm2(exp_a - expected) > tol) then
            write (stderr, *) "test_matrix_exponential failed: Incorrect matrix "// &
                "values."
            test_matrix_exponential = .false.
        end if
        deallocate(exp_a)

        ! call routine for a vanishing matrix and determine if the identity matrix is
        ! returned
        a = 0.0_rp
        exp_a = matrix_exponential(a, error)
        if (norm2(exp_a - identity_matrix(n_ao)) > tol) then
            write (stderr, *) "test_matrix_exponential failed: Incorrect matrix "// &
                "values for vanishing matrix."
            test_matrix_exponential = .false.
        end if
        deallocate(exp_a)

    end function test_matrix_exponential

    logical(c_bool) function test_compute_sqrt_and_inv_sqrt() bind(C)
        !
        ! this function tests the subroutine which calculates the square root and
        ! inverse square root of a matrix
        !
        use otr_oao, only: compute_sqrt_and_inv_sqrt, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao

        real(rp) :: a(n_ao, n_ao)
        real(rp), allocatable :: sqrt_a(:, :), inv_sqrt_a(:, :)
        integer(ip) :: error

        ! assume tests pass
        test_compute_sqrt_and_inv_sqrt = .true.

        ! set up the OAO object, whose settings the routine logs through
        allocate(oao_object)
        call setup_settings(oao_object%settings)

        ! initialize a diagonally dominant, and therefore symmetric positive definite,
        ! matrix
        a = reshape([4.0_rp, 1.0_rp, 0.0_rp, &
                     1.0_rp, 3.0_rp, 1.0_rp, &
                     0.0_rp, 1.0_rp, 2.0_rp], [n_ao, n_ao])

        ! call routine and determine if the square root squares to the matrix and if
        ! the inverse square root is its inverse
        call compute_sqrt_and_inv_sqrt(a, sqrt_a, inv_sqrt_a, error)
        if (error /= 0) then
            write (stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Produced error."
            test_compute_sqrt_and_inv_sqrt = .false.
            deallocate(oao_object)
            return
        end if
        if (norm2(matmul(sqrt_a, sqrt_a) - a) > tol) then
            write (stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Square root "// &
                "does not reproduce the matrix."
            test_compute_sqrt_and_inv_sqrt = .false.
        end if
        if (norm2(matmul(inv_sqrt_a, sqrt_a) - identity_matrix(n_ao)) > tol) then
            write (stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Inverse "// &
                "square root is not the inverse of the square root."
            test_compute_sqrt_and_inv_sqrt = .false.
        end if
        if (norm2(sqrt_a - transpose(sqrt_a)) > tol) then
            write (stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Returned "// &
                "square root is not symmetric."
            test_compute_sqrt_and_inv_sqrt = .false.
        end if
        if (norm2(inv_sqrt_a - transpose(inv_sqrt_a)) > tol) then
            write (stderr, *) "test_compute_sqrt_and_inv_sqrt failed: Returned "// &
                "inverse square root is not symmetric."
            test_compute_sqrt_and_inv_sqrt = .false.
        end if
        deallocate(sqrt_a, inv_sqrt_a, oao_object)

    end function test_compute_sqrt_and_inv_sqrt

    logical(c_bool) function test_rotate_dm_ao() bind(C)
        !
        ! this function tests the subroutine which returns the rotated density matrix
        ! in the AO basis
        !
        use otr_oao, only: rotate_dm_ao, oao_object
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_ao = 2, n_particle = 1, n_electrons = 1
        real(rp), parameter :: angle = 0.3_rp

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), s_inv_sqrt(n_ao, n_ao), &
                    rotation(n_ao, n_ao), expected_dm_oao(n_ao, n_ao, n_particle), &
                    expected_dm_ao(n_ao, n_ao, n_particle), &
                    rot_dm_ao(n_ao, n_ao, n_particle), &
                    rot_dm_oao(n_ao, n_ao, n_particle)
        integer(ip) :: error

        ! assume tests pass
        test_rotate_dm_ao = .true.

        ! set up the OAO object with a non-trivial transformation to the AO basis
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_electrons)
        oao_object%dm_oao = dm_oao
        s_inv_sqrt = reshape([2.0_rp, 0.0_rp, &
                              0.0_rp, 3.0_rp], [n_ao, n_ao])
        oao_object%s_inv_sqrt = s_inv_sqrt

        ! initialize expected density matrices, where the rotation matrix is the
        ! exponential of the antisymmetric matrix the rotation is unpacked into and the
        ! rotated density matrix stays idempotent so that purification leaves it
        ! unchanged
        rotation = reshape([cos(angle), -sin(angle), &
                            sin(angle), cos(angle)], [n_ao, n_ao])
        expected_dm_oao(:, :, 1) = matmul(transpose(rotation), &
                                          matmul(dm_oao(:, :, 1), rotation))
        expected_dm_ao(:, :, 1) = matmul(s_inv_sqrt, &
                                         matmul(expected_dm_oao(:, :, 1), s_inv_sqrt))

        ! call routine and determine if values of the resulting density matrices match
        call rotate_dm_ao([angle], n_particle, n_ao, rot_dm_ao, error, rot_dm_oao)
        if (error /= 0) then
            write (stderr, *) "test_rotate_dm_ao failed: Produced error."
            test_rotate_dm_ao = .false.
            deallocate(oao_object)
            return
        end if
        if (norm2(rot_dm_oao - expected_dm_oao) > tol) then
            write (stderr, *) "test_rotate_dm_ao failed: Incorrect rotated density "// &
                "matrix in the OAO basis."
            test_rotate_dm_ao = .false.
        end if
        if (norm2(rot_dm_ao - expected_dm_ao) > tol) then
            write (stderr, *) "test_rotate_dm_ao failed: Incorrect rotated density "// &
                "matrix in the AO basis."
            test_rotate_dm_ao = .false.
        end if

        ! call routine without the optional rotated density matrix in the OAO basis and
        ! determine if values of the resulting density matrix match
        rot_dm_ao = 0.0_rp
        call rotate_dm_ao([angle], n_particle, n_ao, rot_dm_ao, error)
        if (norm2(rot_dm_ao - expected_dm_ao) > tol) then
            write (stderr, *) "test_rotate_dm_ao failed: Incorrect rotated density "// &
                "matrix in the AO basis without the optional argument."
            test_rotate_dm_ao = .false.
        end if

        ! call routine for a vanishing rotation and determine if the density matrix is
        ! left unchanged
        call rotate_dm_ao([0.0_rp], n_particle, n_ao, rot_dm_ao, error, rot_dm_oao)
        if (norm2(rot_dm_oao - dm_oao) > tol) then
            write (stderr, *) "test_rotate_dm_ao failed: Density matrix not left "// &
                "unchanged for vanishing rotation."
            test_rotate_dm_ao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_rotate_dm_ao

    logical(c_bool) function test_calculate_grad_h_diag() bind(C)
        !
        ! this function tests the subroutine which calculates the gradient, Hessian
        ! diagonal and the occupied-occupied and virtual-virtual parts of the Fock
        ! matrix in the OAO basis
        !
        use otr_oao, only: calculate_grad_h_diag
        use otr_oao_test_reference, only: n_ao

        integer(ip), parameter :: n_electrons = 2

        integer(ip) :: n_particle, n_param
        real(rp) :: dm_oao(n_ao, n_ao, 2), fock_oao(n_ao, n_ao, 2), &
                    proj_v(n_ao, n_ao), expected_fock_oo(n_ao, n_ao, 2), &
                    expected_fock_vv(n_ao, n_ao, 2), fock_ov(n_ao, n_ao, 2), &
                    grad_full(n_ao, n_ao, 2), fock_oo(n_ao, n_ao, 2), &
                    fock_vv(n_ao, n_ao, 2)
        real(rp), allocatable :: grad(:), h_diag(:), expected_grad(:), &
                                 expected_h_diag(:)
        integer(ip) :: i, j, k, idx

        ! assume tests pass
        test_calculate_grad_h_diag = .true.

        ! initialize density and Fock matrices for both particle slots, since the
        ! closed- and open-shell cases below both read from these shared arrays
        do i = 1, size(dm_oao, 3)
            dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_electrons)
            fock_oao(:, :, i) = generate_random_symm_matrix(n_ao)
        end do

        ! initialize expected occupancy-resolved parts of the Fock matrix
        do i = 1, size(dm_oao, 3)
            proj_v = identity_matrix(n_ao) - dm_oao(:, :, i)
            expected_fock_oo(:, :, i) = matmul(dm_oao(:, :, i), &
                                               matmul(fock_oao(:, :, i), &
                                                      dm_oao(:, :, i)))
            expected_fock_vv(:, :, i) = matmul(proj_v, matmul(fock_oao(:, :, i), &
                                                              proj_v))
            fock_ov(:, :, i) = matmul(dm_oao(:, :, i), matmul(fock_oao(:, :, i), &
                                                              proj_v))
        end do

        ! initialize expected gradient and Hessian diagonal for the closed-shell case
        n_particle = 1
        n_param = n_particle * n_ao * (n_ao - 1) / 2
        do i = 1, n_particle
            grad_full(:, :, i) = 4.0_rp * (fock_ov(:, :, i) - &
                                           transpose(fock_ov(:, :, i)))
        end do
        expected_grad = ref_pack_asymm(grad_full(:, :, 1:n_particle), n_param)
        allocate(expected_h_diag(n_param))
        idx = 1
        do k = 1, n_particle
            do j = 1, n_ao
                do i = 1, j - 1
                    expected_h_diag(idx) = &
                        4.0_rp * (expected_fock_vv(i, i, k) + &
                                  expected_fock_vv(j, j, k) - &
                                  expected_fock_oo(i, i, k) - expected_fock_oo(j, j, k))
                    idx = idx + 1
                end do
            end do
        end do

        ! call routine and determine if values of resulting quantities match
        allocate(grad(n_param), h_diag(n_param))
        call calculate_grad_h_diag(dm_oao, fock_oao, n_particle, n_ao, grad, h_diag, &
                                   fock_oo, fock_vv)
        if (norm2(grad - expected_grad) > tol) then
            write (stderr, *) "test_calculate_grad_h_diag failed: Incorrect "// &
                "gradient for closed-shell case."
            test_calculate_grad_h_diag = .false.
        end if
        if (norm2(h_diag - &
                  expected_h_diag) > tol) then
            write (stderr, *) "test_calculate_grad_h_diag failed: Incorrect "// &
                "Hessian diagonal for closed-shell case."
            test_calculate_grad_h_diag = .false.
        end if
        deallocate(expected_h_diag, grad, h_diag)

        ! initialize expected gradient and Hessian diagonal for the open-shell case
        n_particle = 2
        n_param = n_particle * n_ao * (n_ao - 1) / 2
        do i = 1, n_particle
            grad_full(:, :, i) = 2.0_rp * (fock_ov(:, :, i) - &
                                           transpose(fock_ov(:, :, i)))
        end do
        expected_grad = ref_pack_asymm(grad_full, n_param)
        allocate(expected_h_diag(n_param))
        idx = 1
        do k = 1, n_particle
            do j = 1, n_ao
                do i = 1, j - 1
                    expected_h_diag(idx) = &
                        2.0_rp * (expected_fock_vv(i, i, k) + &
                                  expected_fock_vv(j, j, k) - &
                                  expected_fock_oo(i, i, k) - &
                                  expected_fock_oo(j, j, k))
                    idx = idx + 1
                end do
            end do
        end do

        ! call routine and determine if values of resulting quantities match
        allocate(grad(n_param), h_diag(n_param))
        call calculate_grad_h_diag(dm_oao, fock_oao, n_particle, n_ao, grad, h_diag, &
                                   fock_oo, fock_vv)
        if (norm2(fock_oo - expected_fock_oo) > tol) then
            write (stderr, *) "test_calculate_grad_h_diag failed: Incorrect "// &
                "occupied-occupied part of the Fock matrix."
            test_calculate_grad_h_diag = .false.
        end if
        if (norm2(fock_vv - expected_fock_vv) > tol) then
            write (stderr, *) "test_calculate_grad_h_diag failed: Incorrect "// &
                "virtual-virtual part of the Fock matrix."
            test_calculate_grad_h_diag = .false.
        end if
        if (norm2(grad - expected_grad) > tol) then
            write (stderr, *) "test_calculate_grad_h_diag failed: Incorrect "// &
                "gradient for open-shell case."
            test_calculate_grad_h_diag = .false.
        end if
        if (norm2(h_diag - expected_h_diag) > tol) then
            write (stderr, *) "test_calculate_grad_h_diag failed: Incorrect "// &
                "Hessian diagonal for open-shell case."
            test_calculate_grad_h_diag = .false.
        end if
        deallocate(expected_h_diag, grad, h_diag)

    end function test_calculate_grad_h_diag

    logical(c_bool) function test_refresh_oao_response() bind(C)
        !
        ! this function tests the subroutine which rebuilds the OAO response callbacks 
        ! at the currently stored density matrix
        !
        use otr_oao, only: refresh_oao_response, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        integer(ip) :: i, error

        ! assume tests pass
        test_refresh_oao_response = .true.

        ! set up the OAO object with a density matrix and a mock density matrix
        ! updating function, and mark the response as stale as ARH would after moving
        ! the density on its own
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        do i = 1, n_particle
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, 2_ip)
        end do
        oao_object%dm_ao => dm_ao
        oao_object%update_dm_cs => mock_update_dm_cs
        oao_object%response_stale = .true.

        ! reset mock density matrix updating function call count
        n_mock_calls = 0

        ! call routine and determine if the density matrix updating function was
        ! called to rebuild the response and if the flag was cleared
        call refresh_oao_response(error)
        if (error /= 0) then
            write (stderr, *) "test_refresh_oao_response failed: Produced error "// &
                "for the closed-shell case."
            test_refresh_oao_response = .false.
            deallocate(oao_object)
            return
        end if
        if (n_mock_calls /= 1) then
            write (stderr, *) "test_refresh_oao_response failed: Density matrix "// &
                "updating function was not called for the closed-shell case."
            test_refresh_oao_response = .false.
        end if
        if (.not. associated(oao_object%get_response_cs, mock_get_response_cs)) then
            write (stderr, *) "test_refresh_oao_response failed: Response function "// &
                "not updated for the closed-shell case."
            test_refresh_oao_response = .false.
        end if
        if (oao_object%response_stale) then
            write (stderr, *) "test_refresh_oao_response failed: Response still "// &
                "marked stale after being refreshed for the closed-shell case."
            test_refresh_oao_response = .false.
        end if

        ! repeat for the open-shell case
        oao_object%update_dm_cs => null()
        oao_object%update_dm_os => mock_update_dm_os
        oao_object%response_stale = .true.
        n_mock_calls = 0
        call refresh_oao_response(error)
        if (error /= 0) then
            write (stderr, *) "test_refresh_oao_response failed: Produced error "// &
                "for the open-shell case."
            test_refresh_oao_response = .false.
            deallocate(oao_object)
            return
        end if
        if (n_mock_calls /= 1) then
            write (stderr, *) "test_refresh_oao_response failed: Density matrix "// &
                "updating function was not called for the open-shell case."
            test_refresh_oao_response = .false.
        end if
        if (.not. associated(oao_object%get_response_os, mock_get_response_os)) then
            write (stderr, *) "test_refresh_oao_response failed: Response function "// &
                "not updated for the open-shell case."
            test_refresh_oao_response = .false.
        end if
        if (oao_object%response_stale) then
            write (stderr, *) "test_refresh_oao_response failed: Response still "// &
                "marked stale after being refreshed for the open-shell case."
            test_refresh_oao_response = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_refresh_oao_response

    logical(c_bool) function test_obj_func_oao() bind(C)
        !
        ! this function tests the function which defines the energy evaluation in the
        ! OAO basis
        !
        use otr_oao, only: obj_func_oao, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle, n_param

        integer(ip), parameter :: n_electrons = 2

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), kappa(n_param), energy
        integer(ip) :: i, error

        ! assume tests pass
        test_obj_func_oao = .true.

        ! set up the OAO object with an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        do i = 1, n_particle
            dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_electrons)
        end do
        oao_object%dm_oao = dm_oao
        oao_object%s_inv_sqrt = identity_matrix(n_ao)

        ! call routine without an orbital rotation for the closed-shell case and
        ! determine if the energy of the energy function is returned, where only the
        ! first spin channel contributes
        oao_object%get_energy_cs => mock_get_energy_cs
        energy = obj_func_oao(kappa, error)
        if (error /= 0) then
            write (stderr, *) "test_obj_func_oao failed: Produced error for "// &
                "closed-shell case."
            test_obj_func_oao = .false.
        end if
        if (abs(energy - sum(dm_oao(:, :, 1))) > tol) then
            write (stderr, *) "test_obj_func_oao failed: Incorrect energy for "// &
                "closed-shell case."
            test_obj_func_oao = .false.
        end if

        ! call routine without an orbital rotation for the open-shell case and
        ! determine if the energy of the energy function is returned
        oao_object%get_energy_cs => null()
        oao_object%get_energy_os => mock_get_energy_os
        kappa = 0.0_rp
        energy = obj_func_oao(kappa, error)
        if (error /= 0) then
            write (stderr, *) "test_obj_func_oao failed: Produced error for "// &
                "open-shell case."
            test_obj_func_oao = .false.
        end if
        if (abs(energy - sum(dm_oao)) > tol) then
            write (stderr, *) "test_obj_func_oao failed: Incorrect energy for "// &
                "open-shell case."
            test_obj_func_oao = .false.
        end if
        deallocate(oao_object)

    end function test_obj_func_oao

    logical(c_bool) function test_update_orbs_oao() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation in the OAO basis
        !
        use otr_oao, only: update_orbs_oao, oao_object, hess_x_oao_ptr
        use opentrustregion, only: hess_x_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle, n_param

        integer(ip), parameter :: n_electrons = 2

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: kappa(n_param), grad(n_param), h_diag(n_param), func
        integer(ip) :: i, error
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_update_orbs_oao = .true.

        ! set up the OAO object with an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        do i = 1, n_particle
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_electrons)
        end do
        oao_object%dm_ao => dm_ao
        oao_object%dm_oao = dm_ao
        oao_object%s_inv_sqrt = identity_matrix(n_ao)
        allocate(oao_object%fock_oo(n_ao, n_ao, n_particle), &
                 oao_object%fock_vv(n_ao, n_ao, n_particle))
        oao_object%update_dm_os => mock_update_dm_os

        ! reset mock density matrix updating function
        n_mock_calls = 0

        ! call routine without an orbital rotation for an uninitialized object and
        ! determine if the energy, gradient, Hessian diagonal, and Hessian linear 
        ! transformation are correct and if the response function is set
        kappa = 0.0_rp
        oao_object%hess_eigen_stale = .false.
        call update_orbs_oao(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_oao failed: Produced error."
            test_update_orbs_oao = .false.
            deallocate(oao_object)
            return
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_oao failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after the static Hessian part was rebuilt."
            test_update_orbs_oao = .false.
        end if
        if (abs(func - sum(dm_ao)) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Incorrect energy."
            test_update_orbs_oao = .false.
        end if
        if (norm2(grad - oao_object%grad) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Gradient not returned."
            test_update_orbs_oao = .false.
        end if
        if (norm2(h_diag - oao_object%h_diag) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Hessian diagonal not "// &
                "returned."
            test_update_orbs_oao = .false.
        end if
        if (.not. associated(oao_object%get_response_os, mock_get_response_os)) then
            write (stderr, *) "test_update_orbs_oao failed: Response function not "// &
                "stored correctly."
            test_update_orbs_oao = .false.
        end if
        if (.not. associated(hess_x_funptr, hess_x_oao_ptr)) then
            write (stderr, *) "test_update_orbs_oao failed: Returned Hessian "// &
                "linear transformation is wrong."
            test_update_orbs_oao = .false.
        end if

        ! call routine again without an orbital rotation and determine if the
        ! quantities of the already initialized object are reused, including the
        ! cached eigendecomposition of the static Hessian part, which should remain
        ! valid since it was not rebuilt
        oao_object%hess_eigen_stale = .false.
        call update_orbs_oao(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (n_mock_calls /= 1) then
            write (stderr, *) "test_update_orbs_oao failed: Quantities recomputed "// &
                "without an orbital rotation."
            test_update_orbs_oao = .false.
        end if
        if (oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_oao failed: Cached "// &
                "eigendecomposition of the static Hessian part marked stale even "// &
                "though the static Hessian part was not rebuilt."
            test_update_orbs_oao = .false.
        end if

        ! mark the response as stale (as an approximate-Hessian extension such as ARH
        ! would after moving the density without going through this routine) and call
        ! again without an orbital rotation, and determine if this still forces a
        ! recompute and clears the flag
        oao_object%response_stale = .true.
        call update_orbs_oao(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (n_mock_calls /= 2) then
            write (stderr, *) "test_update_orbs_oao failed: Quantities not "// &
                "recomputed for a stale response."
            test_update_orbs_oao = .false.
        end if
        if (oao_object%response_stale) then
            write (stderr, *) "test_update_orbs_oao failed: Response still "// &
                "marked stale after being recomputed."
            test_update_orbs_oao = .false.
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_oao failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after the static Hessian part was rebuilt for a stale response."
            test_update_orbs_oao = .false.
        end if

        ! call routine with an orbital rotation and determine if the energy, gradient,
        ! Hessian diagonal, and Hessian linear transformation are recomputed and
        ! correct and if the response function is still set
        kappa = 0.1_rp
        oao_object%hess_eigen_stale = .false.
        call update_orbs_oao(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_oao failed: Produced error after "// &
                "an orbital rotation."
            test_update_orbs_oao = .false.
        end if
        if (n_mock_calls /= 3) then
            write (stderr, *) "test_update_orbs_oao failed: Quantities not "// &
                "recomputed after an orbital rotation."
            test_update_orbs_oao = .false.
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_oao failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after an orbital rotation."
            test_update_orbs_oao = .false.
        end if
        if (abs(func - sum(oao_object%dm_oao)) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Incorrect energy after "// &
                "an orbital rotation."
            test_update_orbs_oao = .false.
        end if
        if (norm2(grad - oao_object%grad) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Gradient not returned "// &
                "after an orbital rotation."
            test_update_orbs_oao = .false.
        end if
        if (norm2(h_diag - oao_object%h_diag) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Hessian diagonal not "// &
                "returned after an orbital rotation."
            test_update_orbs_oao = .false.
        end if
        if (.not. associated(oao_object%get_response_os, mock_get_response_os)) then
            write (stderr, *) "test_update_orbs_oao failed: Response function not "// &
                "stored correctly after an orbital rotation."
            test_update_orbs_oao = .false.
        end if
        if (.not. associated(hess_x_funptr, hess_x_oao_ptr)) then
            write (stderr, *) "test_update_orbs_oao failed: Returned Hessian "// &
                "linear transformation is wrong after an orbital rotation."
            test_update_orbs_oao = .false.
        end if

        ! call routine for the closed-shell case and determine if the energy is correct 
        ! and if the energy, gradient, Hessian diagonal, response function, and Hessian 
        ! linear transformation are correct
        oao_object%update_dm_os => null()
        oao_object%update_dm_cs => mock_update_dm_cs
        oao_object%hess_eigen_stale = .false.
        call update_orbs_oao(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_oao failed: Produced error for the "// &
                "closed-shell case."
            test_update_orbs_oao = .false.
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_oao failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "for the closed-shell case."
            test_update_orbs_oao = .false.
        end if
        if (abs(func - sum(oao_object%dm_oao(:, :, 1))) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Incorrect energy for "// &
                "closed-shell case."
            test_update_orbs_oao = .false.
        end if
        if (norm2(grad - oao_object%grad) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Gradient not returned "// &
                "for the closed-shell case."
            test_update_orbs_oao = .false.
        end if
        if (norm2(h_diag - oao_object%h_diag) > tol) then
            write (stderr, *) "test_update_orbs_oao failed: Hessian diagonal not "// &
                "returned for the closed-shell case."
            test_update_orbs_oao = .false.
        end if
        if (.not. associated(oao_object%get_response_cs, mock_get_response_cs)) then
            write (stderr, *) "test_update_orbs_oao failed: Closed-shell response "// &
                "function not stored correctly."
            test_update_orbs_oao = .false.
        end if
        if (.not. associated(hess_x_funptr, hess_x_oao_ptr)) then
            write (stderr, *) "test_update_orbs_oao failed: Returned Hessian "// &
                "linear transformation is wrong for the closed-shell case."
            test_update_orbs_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_update_orbs_oao

    logical(c_bool) function test_hess_x_oao() bind(C)
        !
        ! this function tests the subroutine which defines the Hessian linear
        ! transformation in the OAO basis
        !
        use otr_oao, only: hess_x_oao, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao

        ! multiplier the mock response function applies to the density matrix response
        real(rp), parameter :: mock_response_factor = 2.0_rp

        integer(ip) :: n_particle, n_param
        real(rp), target :: dm_oao(n_ao, n_ao, 2), fock_oo(n_ao, n_ao, 2), &
                            fock_vv(n_ao, n_ao, 2)
        real(rp), allocatable :: x(:)
        real(rp), allocatable :: x_full(:, :, :), delta_dm(:, :, :), hess_x(:), &
                                 expected_hess_x(:)
        integer(ip) :: j, error

        ! assume tests pass
        test_hess_x_oao = .true.

        ! generate random density matrices, Fock matrix contributions and trial vector
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, 2_ip)
        dm_oao(:, :, 2) = generate_random_density_matrix(n_ao, 1_ip)
        do j = 1, 2
            fock_oo(:, :, j) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, j) = generate_random_symm_matrix(n_ao)
        end do
        call random_number(x)

        ! set up the OAO object with an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%dm_oao = dm_oao
        oao_object%fock_oo = fock_oo
        oao_object%fock_vv = fock_vv
        oao_object%s_inv_sqrt = identity_matrix(n_ao)

        ! set up the closed-shell case
        n_particle = 1
        n_param = n_ao * (n_ao - 1) / 2
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%get_response_cs => mock_get_response_cs

        ! the response is the mock response of the density matrix displacement
        allocate(x(n_param))
        x_full = ref_unpack_asymm(x, n_particle, n_ao)
        delta_dm = ref_project_symm(x_full, dm_oao(:, :, 1:1))
        expected_hess_x = ref_hess_x(x_full, mock_response_factor * delta_dm, &
                                     dm_oao(:, :, 1:1), fock_oo(:, :, 1:1), &
                                     fock_vv(:, :, 1:1), n_param)

        ! call routine and determine if values of resulting Hessian linear
        ! transformation match
        allocate(hess_x(n_param))
        call hess_x_oao(x, hess_x, error)
        if (error /= 0) then
            write (stderr, *) "test_hess_x_oao failed: Produced error for "// &
                "closed-shell case."
            test_hess_x_oao = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write (stderr, *) "test_hess_x_oao failed: Incorrect Hessian linear "// &
                "transformation for closed-shell case."
            test_hess_x_oao = .false.
        end if
        deallocate(x, hess_x)

        ! set up the open-shell case
        n_particle = 2
        n_param = n_particle * n_param
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%get_response_cs => null()
        oao_object%get_response_os => mock_get_response_os

        ! the response is the mock response of the density matrix displacements
        allocate(x(n_param))
        x_full = ref_unpack_asymm(x, n_particle, n_ao)
        delta_dm = ref_project_symm(x_full, dm_oao)
        expected_hess_x = ref_hess_x(x_full, mock_response_factor * delta_dm, dm_oao, &
                                     fock_oo, fock_vv, n_param)

        ! call routine and determine if values of resulting Hessian linear
        ! transformation match
        allocate(hess_x(n_param))
        call hess_x_oao(x, hess_x, error)
        if (error /= 0) then
            write (stderr, *) "test_hess_x_oao failed: Produced error for "// &
                "open-shell case."
            test_hess_x_oao = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write (stderr, *) "test_hess_x_oao failed: Incorrect Hessian linear "// &
                "transformation for open-shell case."
            test_hess_x_oao = .false.
        end if
        deallocate(x, hess_x)

        ! mark the response as stale (as ARH would after updating the density without 
        ! going through update_orbs_oao) and call again, and determine if this triggers 
        ! the density matrix updating function to refresh the response and clears the 
        ! flag
        oao_object%dm_ao => dm_oao
        oao_object%update_dm_os => mock_update_dm_os
        oao_object%response_stale = .true.
        n_mock_calls = 0
        allocate(x(n_param), hess_x(n_param))
        call random_number(x)
        call hess_x_oao(x, hess_x, error)
        if (error /= 0) then
            write (stderr, *) "test_hess_x_oao failed: Produced error while "// &
                "refreshing a stale response."
            test_hess_x_oao = .false.
        end if
        if (n_mock_calls /= 1) then
            write (stderr, *) "test_hess_x_oao failed: Stale response was not "// &
                "refreshed."
            test_hess_x_oao = .false.
        end if
        if (oao_object%response_stale) then
            write (stderr, *) "test_hess_x_oao failed: Response still marked "// &
                "stale after being refreshed."
            test_hess_x_oao = .false.
        end if
        deallocate(x, hess_x)

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_hess_x_oao

    logical(c_bool) function test_project_oao() bind(C)
        !
        ! this function tests the subroutine which discards the redundant
        ! occupied-occupied and virtual-virtual rotations from a vector, retaining
        ! only its occupied-virtual and virtual-occupied contributions
        !
        use otr_oao, only: project_oao, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle, n_param


        real(rp) :: dm_oao(n_ao, n_ao, n_particle), vector(n_param), expected(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_project_oao = .true.

        ! generate random density matrices and vector
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, 2_ip)
        dm_oao(:, :, 2) = generate_random_density_matrix(n_ao, 1_ip)
        call random_number(vector)

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%dm_oao = dm_oao

        ! initialize expected vector
        oao_object%n_particle = n_particle
        expected = ref_pack_asymm(ref_project_asymm(ref_unpack_asymm( &
            vector, n_particle, n_ao), dm_oao), n_param)

        ! call routine and determine if values of resulting vector match
        call project_oao(vector, error)
        if (error /= 0) then
            write (stderr, *) "test_project_oao failed: Produced error."
            test_project_oao = .false.
        end if
        if (norm2(vector - expected) > tol) then
            write (stderr, *) "test_project_oao failed: Incorrect vector."
            test_project_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_project_oao

    logical(c_bool) function test_precond_oao() bind(C)
        !
        ! this function tests the subroutine that defines a level-shifted
        ! preconditioner based on the exact eigendecomposition of the static part of
        ! the Hessian
        !
        use otr_oao, only: precond_oao, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use opentrustregion, only: precond_floor
        use otr_oao_test_reference, only: n_ao, n_particle, n_param


        real(rp) :: fock_oo(n_ao, n_ao, n_particle), fock_vv(n_ao, n_ao, n_particle), &
                    eigvecs(n_ao, n_ao, n_particle), eigvals(n_ao, n_particle), &
                    eigval_pairs(n_param), residual(n_param), &
                    precond_residual(n_param), rotated_residual(n_param), &
                    expected(n_param), mu
        integer(ip) :: j, error

        ! assume tests pass
        test_precond_oao = .true.

        ! generate random Fock matrix contributions, ensuring a non-diagonal static
        ! part, and a residual vector and level shift
        do j = 1, n_particle
            fock_oo(:, :, j) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, j) = generate_random_symm_matrix(n_ao)
        end do
        call random_number(residual)
        mu = 0.3_rp

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%fock_oo = fock_oo
        oao_object%fock_vv = fock_vv

        ! independently diagonalize the static part and construct the expected
        ! preconditioned residual
        call ref_diagonalize_static_part(fock_oo, fock_vv, eigvecs, eigvals)
        eigval_pairs = ref_hess_eigval_pairs(eigvals, n_particle, n_ao, n_param)
        rotated_residual = ref_rotate_to_eigenbasis(residual, eigvecs, n_particle, n_ao)
        eigval_pairs = eigval_pairs - mu
        where (abs(eigval_pairs) < precond_floor) eigval_pairs = precond_floor
        rotated_residual = rotated_residual / eigval_pairs
        expected = ref_rotate_from_eigenbasis(rotated_residual, eigvecs, n_particle, &
                                              n_ao)

        ! call routine and determine if values of the preconditioned residual match
        call precond_oao(residual, mu, precond_residual, error)
        if (error /= 0) then
            write (stderr, *) "test_precond_oao failed: Produced error."
            test_precond_oao = .false.
        end if
        if (norm2(precond_residual - expected) > tol) then
            write (stderr, *) "test_precond_oao failed: Incorrect preconditioned "// &
                "residual."
            test_precond_oao = .false.
        end if
        if (oao_object%hess_eigen_stale) then
            write (stderr, *) "test_precond_oao failed: Eigendecomposition still "// &
                "marked stale after being refreshed."
            test_precond_oao = .false.
        end if

        ! mutate the Fock matrix without marking the cache stale, and confirm the
        ! routine reuses the cached eigendecomposition rather than recomputing it
        oao_object%fock_oo = fock_oo + 1.0_rp
        oao_object%fock_vv = fock_vv + 1.0_rp
        call precond_oao(residual, mu, precond_residual, error)
        if (error /= 0 .or. norm2(precond_residual - expected) > tol) then
            write (stderr, *) "test_precond_oao failed: Recomputed the "// &
                "eigendecomposition despite the cache not being marked stale."
            test_precond_oao = .false.
        end if

        ! mark the cache stale and confirm the eigendecomposition is refreshed to
        ! reflect the mutated Fock matrix
        oao_object%hess_eigen_stale = .true.
        call ref_diagonalize_static_part(fock_oo + 1.0_rp, fock_vv + 1.0_rp, eigvecs, &
                                         eigvals)
        eigval_pairs = ref_hess_eigval_pairs(eigvals, n_particle, n_ao, n_param)
        rotated_residual = ref_rotate_to_eigenbasis(residual, eigvecs, n_particle, n_ao)
        eigval_pairs = eigval_pairs - mu
        where (abs(eigval_pairs) < precond_floor) eigval_pairs = precond_floor
        rotated_residual = rotated_residual / eigval_pairs
        expected = ref_rotate_from_eigenbasis(rotated_residual, eigvecs, n_particle, &
                                              n_ao)
        call precond_oao(residual, mu, precond_residual, error)
        if (error /= 0 .or. norm2(precond_residual - expected) > tol) then
            write (stderr, *) "test_precond_oao failed: Did not refresh the "// &
                "eigendecomposition after being marked stale."
            test_precond_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_precond_oao

    logical(c_bool) function test_precond_pd_oao() bind(C)
        !
        ! this function tests the subroutine that defines the positive-definite
        ! preconditioner based on the exact eigendecomposition of the static part of
        ! the Hessian
        !
        use otr_oao, only: precond_pd_oao, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use opentrustregion, only: precond_floor, precond_rel_floor_factor
        use otr_oao_test_reference, only: n_ao, n_particle, n_param


        real(rp) :: fock_oo(n_ao, n_ao, n_particle), fock_vv(n_ao, n_ao, n_particle), &
                    eigvecs(n_ao, n_ao, n_particle), eigvals(n_ao, n_particle), &
                    eigval_pairs(n_param), residual(n_param), &
                    precond_residual(n_param), rotated_residual(n_param), &
                    expected(n_param), floor_val
        integer(ip) :: j, error

        ! assume tests pass
        test_precond_pd_oao = .true.

        ! generate random Fock matrix contributions, ensuring a non-diagonal static
        ! part, and a residual vector
        do j = 1, n_particle
            fock_oo(:, :, j) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, j) = generate_random_symm_matrix(n_ao)
        end do
        call random_number(residual)

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%fock_oo = fock_oo
        oao_object%fock_vv = fock_vv

        ! independently diagonalize the static part and construct the expected
        ! preconditioned residual
        call ref_diagonalize_static_part(fock_oo, fock_vv, eigvecs, eigvals)
        eigval_pairs = ref_hess_eigval_pairs(eigvals, n_particle, n_ao, n_param)
        rotated_residual = ref_rotate_to_eigenbasis(residual, eigvecs, n_particle, n_ao)
        floor_val = max(precond_rel_floor_factor * maxval(abs(eigval_pairs)), &
                        precond_floor)
        eigval_pairs = abs(eigval_pairs)
        where (eigval_pairs < floor_val) eigval_pairs = floor_val
        rotated_residual = rotated_residual / eigval_pairs
        expected = ref_rotate_from_eigenbasis(rotated_residual, eigvecs, n_particle, &
                                              n_ao)

        ! call routine and determine if values of the preconditioned residual match
        call precond_pd_oao(residual, precond_residual, error)
        if (error /= 0) then
            write (stderr, *) "test_precond_pd_oao failed: Produced error."
            test_precond_pd_oao = .false.
        end if
        if (norm2(precond_residual - expected) > tol) then
            write (stderr, *) "test_precond_pd_oao failed: Incorrect "// &
                "preconditioned residual."
            test_precond_pd_oao = .false.
        end if
        if (oao_object%hess_eigen_stale) then
            write (stderr, *) "test_precond_pd_oao failed: Eigendecomposition "// &
                "still marked stale after being refreshed."
            test_precond_pd_oao = .false.
        end if

        ! mutate the Fock matrix without marking the cache stale, and confirm the
        ! routine reuses the cached eigendecomposition rather than recomputing it
        oao_object%fock_oo = fock_oo + 1.0_rp
        oao_object%fock_vv = fock_vv + 1.0_rp
        call precond_pd_oao(residual, precond_residual, error)
        if (error /= 0 .or. norm2(precond_residual - expected) > tol) then
            write (stderr, *) "test_precond_pd_oao failed: Recomputed the "// &
                "eigendecomposition despite the cache not being marked stale."
            test_precond_pd_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_precond_pd_oao

    logical(c_bool) function test_refresh_hess_eigen() bind(C)
        !
        ! this function tests the subroutine that diagonalizes the static part of the
        ! Hessian and caches the result
        !
        use otr_oao, only: refresh_hess_eigen, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle

        real(rp) :: fock_oo(n_ao, n_ao, n_particle), fock_vv(n_ao, n_ao, n_particle), &
                    expected_eigvecs(n_ao, n_ao, n_particle), &
                    expected_eigvals(n_ao, n_particle)
        integer(ip) :: j, error

        ! assume tests pass
        test_refresh_hess_eigen = .true.

        ! generate random Fock matrix contributions
        do j = 1, n_particle
            fock_oo(:, :, j) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, j) = generate_random_symm_matrix(n_ao)
        end do

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%fock_oo = fock_oo
        oao_object%fock_vv = fock_vv

        ! independently diagonalize the static part
        call ref_diagonalize_static_part(fock_oo, fock_vv, expected_eigvecs, &
                                         expected_eigvals)

        ! call routine and determine if the cached eigendecomposition matches
        call refresh_hess_eigen(error)
        if (error /= 0) then
            write (stderr, *) "test_refresh_hess_eigen failed: Produced error."
            test_refresh_hess_eigen = .false.
        end if
        if (norm2(oao_object%hess_eigvecs - expected_eigvecs) > tol) then
            write (stderr, *) "test_refresh_hess_eigen failed: Incorrect eigenvectors."
            test_refresh_hess_eigen = .false.
        end if
        if (norm2(oao_object%hess_eigvals - expected_eigvals) > tol) then
            write (stderr, *) "test_refresh_hess_eigen failed: Incorrect eigenvalues."
            test_refresh_hess_eigen = .false.
        end if
        if (oao_object%hess_eigen_stale) then
            write (stderr, *) "test_refresh_hess_eigen failed: Eigendecomposition "// &
                "still marked stale after being refreshed."
            test_refresh_hess_eigen = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_refresh_hess_eigen

    logical(c_bool) function test_rotate_to_hess_eigenbasis() bind(C)
        !
        ! this function tests the function that rotates a packed antisymmetric vector 
        ! into the eigenbasis of the cached static Hessian part
        !
        use otr_oao, only: rotate_to_hess_eigenbasis, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle, n_param


        real(rp) :: eigvecs(n_ao, n_ao, n_particle), vector(n_param), &
                    expected(n_param), rotated(n_param)
        integer(ip) :: i, j

        ! assume tests pass
        test_rotate_to_hess_eigenbasis = .true.

        ! use a cyclic permutation matrix as a simple, exactly orthogonal, non-identity 
        ! rotation
        eigvecs = 0.0_rp
        do j = 1, n_particle
            do i = 1, n_ao
                eigvecs(i, mod(i, n_ao) + 1, j) = 1.0_rp
            end do
        end do
        call random_number(vector)

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%hess_eigvecs = eigvecs

        ! independently construct the expected rotated vector
        expected = ref_rotate_to_eigenbasis(vector, eigvecs, n_particle, n_ao)

        ! call routine and determine if values match
        rotated = rotate_to_hess_eigenbasis(vector)
        if (norm2(rotated - expected) > tol) then
            write (stderr, *) "test_rotate_to_hess_eigenbasis failed: Incorrect "// &
                "vector when rotating into the eigenbasis."
            test_rotate_to_hess_eigenbasis = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_rotate_to_hess_eigenbasis

    logical(c_bool) function test_rotate_from_hess_eigenbasis() bind(C)
        !
        ! this function tests the function that rotates a packed antisymmetric
        ! vector out of the eigenbasis of the cached static Hessian part
        !
        use otr_oao, only: rotate_from_hess_eigenbasis, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle, n_param


        real(rp) :: eigvecs(n_ao, n_ao, n_particle), vector(n_param), &
                    expected(n_param), rotated(n_param)
        integer(ip) :: i, j

        ! assume tests pass
        test_rotate_from_hess_eigenbasis = .true.

        ! use a cyclic permutation matrix as a simple, exactly orthogonal, non-identity 
        ! rotation
        eigvecs = 0.0_rp
        do j = 1, n_particle
            do i = 1, n_ao
                eigvecs(i, mod(i, n_ao) + 1, j) = 1.0_rp
            end do
        end do
        call random_number(vector)

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%hess_eigvecs = eigvecs

        ! independently construct the expected rotated vector
        expected = ref_rotate_from_eigenbasis(vector, eigvecs, n_particle, n_ao)

        ! call routine and determine if values match
        rotated = rotate_from_hess_eigenbasis(vector)
        if (norm2(rotated - expected) > tol) then
            write (stderr, *) "test_rotate_from_hess_eigenbasis failed: Incorrect "// &
                "vector when rotating out of the eigenbasis."
            test_rotate_from_hess_eigenbasis = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_rotate_from_hess_eigenbasis

    logical(c_bool) function test_get_hess_eigval_pairs() bind(C)
        !
        ! this function tests the function that returns the pairwise sums of the cached 
        ! static Hessian part eigenvalues; the closed-shell and open-shell cases apply 
        ! different scaling factors
        !
        use otr_oao, only: get_hess_eigval_pairs, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle

        real(rp), allocatable :: eigvals(:, :), expected(:), eigval_pairs(:)

        ! assume tests pass
        test_get_hess_eigval_pairs = .true.

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao

        ! closed-shell case
        oao_object%n_particle = 1
        oao_object%n_param = oao_object%n_particle * n_ao * (n_ao - 1) / 2
        allocate(eigvals(n_ao, oao_object%n_particle))
        call random_number(eigvals)
        oao_object%hess_eigvals = eigvals

        ! independently construct the expected pairwise sums
        expected = ref_hess_eigval_pairs(eigvals, oao_object%n_particle, n_ao, &
                                         oao_object%n_param)

        ! call routine and determine if values match
        eigval_pairs = get_hess_eigval_pairs()
        if (norm2(eigval_pairs - expected) > tol) then
            write (stderr, *) "test_get_hess_eigval_pairs failed: Incorrect "// &
                "pairwise eigenvalue sums for the closed-shell case."
            test_get_hess_eigval_pairs = .false.
        end if
        deallocate(eigvals)

        ! open-shell case
        oao_object%n_particle = n_particle
        oao_object%n_param = oao_object%n_particle * n_ao * (n_ao - 1) / 2
        allocate(eigvals(n_ao, oao_object%n_particle))
        call random_number(eigvals)
        oao_object%hess_eigvals = eigvals

        ! independently construct the expected pairwise sums
        expected = ref_hess_eigval_pairs(eigvals, oao_object%n_particle, n_ao, &
                                         oao_object%n_param)

        ! call routine and determine if values match
        eigval_pairs = get_hess_eigval_pairs()
        if (norm2(eigval_pairs - expected) > tol) then
            write (stderr, *) "test_get_hess_eigval_pairs failed: Incorrect "// &
                "pairwise eigenvalue sums for the open-shell case."
            test_get_hess_eigval_pairs = .false.
        end if
        deallocate(eigvals)

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_get_hess_eigval_pairs

    logical(c_bool) function test_oao_factory_cs() bind(C)
        !
        ! this function tests the subroutine which returns the modified OAO orbital
        ! updating function for the closed-shell case
        !
        use otr_oao, only: oao_factory_cs, oao_object, oao_settings_type, &
                           get_energy_cs_type, update_dm_cs_type, obj_func_oao_ptr, &
                           update_orbs_oao_ptr, precond_oao_ptr, precond_pd_oao_ptr, &
                           project_oao_ptr
        use opentrustregion, only: obj_func_type, update_orbs_type, precond_type, &
                                   precond_pd_type, project_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, operator(==)

        integer(ip), parameter :: n_particle = 1, n_electrons = 2, &
                                  n_param = n_ao * (n_ao - 1) / 2

        real(rp), target :: dm_ao(n_ao, n_ao)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: error
        type(oao_settings_type) :: settings
        procedure(get_energy_cs_type), pointer :: get_energy_funptr
        procedure(update_dm_cs_type), pointer :: update_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        procedure(precond_type), pointer :: precond_oao_funptr
        procedure(precond_pd_type), pointer :: precond_pd_oao_funptr
        procedure(project_type), pointer :: project_oao_funptr

        ! assume tests pass
        test_oao_factory_cs = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize density matrix and an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        dm_ao = generate_random_density_matrix(n_ao, n_electrons)
        ao_overlap = identity_matrix(n_ao)

        ! initialize callback function pointers
        get_energy_funptr => mock_get_energy_cs
        update_dm_funptr => mock_update_dm_cs

        ! call routine and determine if an error is produced
        call oao_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, &
                                      get_energy_funptr, update_dm_funptr, &
                                      obj_func_oao_funptr, update_orbs_oao_funptr, &
                                      precond_oao_funptr, precond_pd_oao_funptr, &
                                      project_oao_funptr, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_oao_factory_cs failed: Produced error."
            test_oao_factory_cs = .false.
            return
        end if

        ! determine if OAO object is set up correctly
        if (.not. allocated(oao_object)) then
            write (stderr, *) "test_oao_factory_cs failed: OAO object not allocated."
            test_oao_factory_cs = .false.
            return
        end if
        if (.not. (oao_object%settings == settings)) then
            write (stderr, *) "test_oao_factory_cs failed: Settings not stored "// &
                "correctly."
            test_oao_factory_cs = .false.
        end if
        if (oao_object%n_ao /= n_ao) then
            write (stderr, *) "test_oao_factory_cs failed: Number of AOs not "// &
                "stored correctly."
            test_oao_factory_cs = .false.
        end if
        if (oao_object%n_particle /= n_particle) then
            write (stderr, *) "test_oao_factory_cs failed: Number of particles not "// &
                "stored correctly."
            test_oao_factory_cs = .false.
        end if
        if (oao_object%n_param /= n_param) then
            write (stderr, *) "test_oao_factory_cs failed: Number of parameters "// &
                "not stored correctly."
            test_oao_factory_cs = .false.
        end if
        if (norm2(oao_object%dm_oao(:, :, 1) - dm_ao) > tol) then
            write (stderr, *) "test_oao_factory_cs failed: Density matrix not set "// &
                "up correctly."
            test_oao_factory_cs = .false.
        end if
        if (norm2(oao_object%s_sqrt - identity_matrix(n_ao)) > tol) then
            write (stderr, *) "test_oao_factory_cs failed: Overlap matrix square "// &
                "root not set up correctly."
            test_oao_factory_cs = .false.
        end if
        if (norm2(oao_object%s_inv_sqrt - identity_matrix(n_ao)) > tol) then
            write (stderr, *) "test_oao_factory_cs failed: Overlap matrix inverse "// &
                "square root not set up correctly."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(oao_object%get_energy_cs, mock_get_energy_cs)) then
            write (stderr, *) "test_oao_factory_cs failed: Energy function not "// &
                "stored correctly."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(oao_object%update_dm_cs, mock_update_dm_cs)) then
            write (stderr, *) "test_oao_factory_cs failed: Density matrix updating "// &
                "function not stored correctly."
            test_oao_factory_cs = .false.
        end if

        ! determine if returned function pointers point to the correct routines
        if (.not. associated(obj_func_oao_funptr, obj_func_oao_ptr)) then
            write (stderr, *) "test_oao_factory_cs failed: Returned objective "// &
                "function is wrong."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(update_orbs_oao_funptr, update_orbs_oao_ptr)) then
            write (stderr, *) "test_oao_factory_cs failed: Returned orbital "// &
                "updating function is wrong."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(precond_oao_funptr, precond_oao_ptr)) then
            write (stderr, *) "test_oao_factory_cs failed: Returned level-shifted "// &
                "preconditioner function is wrong."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(precond_pd_oao_funptr, precond_pd_oao_ptr)) then
            write (stderr, *) "test_oao_factory_cs failed: Returned "// &
                "positive-definite preconditioner function is wrong."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(project_oao_funptr, project_oao_ptr)) then
            write (stderr, *) "test_oao_factory_cs failed: Returned projection "// &
                "function is wrong."
            test_oao_factory_cs = .false.
        end if
        deallocate(oao_object)

    end function test_oao_factory_cs

    logical(c_bool) function test_oao_factory_os() bind(C)
        !
        ! this function tests the subroutine which returns the modified OAO orbital
        ! updating function for the open-shell case
        !
        use otr_oao, only: oao_factory_os, oao_object, oao_settings_type, &
                           get_energy_os_type, update_dm_os_type, obj_func_oao_ptr, &
                           update_orbs_oao_ptr, precond_oao_ptr, precond_pd_oao_ptr, &
                           project_oao_ptr
        use opentrustregion, only: obj_func_type, update_orbs_type, precond_type, &
                                   precond_pd_type, project_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_ao, n_particle, n_param, operator(==)

        integer(ip), parameter :: n_electrons = 2

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: i, error
        type(oao_settings_type) :: settings
        procedure(get_energy_os_type), pointer :: get_energy_funptr
        procedure(update_dm_os_type), pointer :: update_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        procedure(precond_type), pointer :: precond_oao_funptr
        procedure(precond_pd_type), pointer :: precond_pd_oao_funptr
        procedure(project_type), pointer :: project_oao_funptr

        ! assume tests pass
        test_oao_factory_os = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize density matrices and an orthonormal AO basis, so that the AO and
        ! the OAO basis coincide
        do i = 1, n_particle
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_electrons)
        end do
        ao_overlap = identity_matrix(n_ao)

        ! initialize callback function pointers
        get_energy_funptr => mock_get_energy_os
        update_dm_funptr => mock_update_dm_os

        ! call routine and determine if an error is produced
        call oao_factory_os(dm_ao, ao_overlap, n_particle, n_ao, get_energy_funptr, &
                            update_dm_funptr, obj_func_oao_funptr, &
                            update_orbs_oao_funptr, precond_oao_funptr, &
                            precond_pd_oao_funptr, project_oao_funptr, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_oao_factory_os failed: Produced error."
            test_oao_factory_os = .false.
            return
        end if

        ! determine if OAO object is set up correctly
        if (.not. allocated(oao_object)) then
            write (stderr, *) "test_oao_factory_os failed: OAO object not allocated."
            test_oao_factory_os = .false.
            return
        end if
        if (.not. (oao_object%settings == settings)) then
            write (stderr, *) "test_oao_factory_os failed: Settings not stored "// &
                "correctly."
            test_oao_factory_os = .false.
        end if
        if (oao_object%n_ao /= n_ao) then
            write (stderr, *) "test_oao_factory_os failed: Number of AOs not "// &
                "stored correctly."
            test_oao_factory_os = .false.
        end if
        if (oao_object%n_particle /= n_particle) then
            write (stderr, *) "test_oao_factory_os failed: Number of particles not "// &
                "stored correctly."
            test_oao_factory_os = .false.
        end if
        if (oao_object%n_param /= n_param) then
            write (stderr, *) "test_oao_factory_os failed: Number of parameters "// &
                "not stored correctly."
            test_oao_factory_os = .false.
        end if
        if (norm2(oao_object%dm_oao - dm_ao) > tol) then
            write (stderr, *) "test_oao_factory_os failed: Density matrices not "// &
                "set up correctly."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(oao_object%get_energy_os, mock_get_energy_os)) then
            write (stderr, *) "test_oao_factory_os failed: Energy function not "// &
                "stored correctly."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(oao_object%update_dm_os, mock_update_dm_os)) then
            write (stderr, *) "test_oao_factory_os failed: Density matrix updating "// &
                "function not stored correctly."
            test_oao_factory_os = .false.
        end if

        ! determine if returned function pointers point to the correct routines
        if (.not. associated(obj_func_oao_funptr, obj_func_oao_ptr)) then
            write (stderr, *) "test_oao_factory_os failed: Returned objective "// &
                "function is wrong."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(update_orbs_oao_funptr, update_orbs_oao_ptr)) then
            write (stderr, *) "test_oao_factory_os failed: Returned orbital "// &
                "updating function is wrong."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(precond_oao_funptr, precond_oao_ptr)) then
            write (stderr, *) "test_oao_factory_os failed: Returned level-shifted "// &
                "preconditioner function is wrong."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(precond_pd_oao_funptr, precond_pd_oao_ptr)) then
            write (stderr, *) "test_oao_factory_os failed: Returned "// &
                "positive-definite preconditioner function is wrong."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(project_oao_funptr, project_oao_ptr)) then
            write (stderr, *) "test_oao_factory_os failed: Returned projection "// &
                "function is wrong."
            test_oao_factory_os = .false.
        end if
        deallocate(oao_object)

    end function test_oao_factory_os

end module otr_oao_unit_tests
