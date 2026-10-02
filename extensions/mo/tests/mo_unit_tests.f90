! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

contains

    function diagonal_matrix(values) result(matrix)
        !
        ! this function returns the diagonal matrix of the given values
        !
        real(rp), intent(in) :: values(:)
        real(rp) :: matrix(size(values), size(values))

        integer(ip) :: i

        matrix = 0.0_rp
        do i = 1, size(values, kind=ip)
            matrix(i, i) = values(i)
        end do

    end function diagonal_matrix

    function generate_random_ao_overlap(n) result(ao_overlap)
        !
        ! this function generates a random symmetric positive definite AO overlap
        ! matrix with a unit diagonal perturbed by random off-diagonal couplings
        !
        integer(ip), intent(in) :: n
        real(rp) :: ao_overlap(n, n)

        real(rp) :: a(n, n)
        integer(ip) :: i

        call random_number(a)
        ao_overlap = 0.2_rp * matmul(a, transpose(a)) / real(n, kind=rp)
        do i = 1, n
            ao_overlap(i, i) = ao_overlap(i, i) + 1.0_rp
        end do

    end function generate_random_ao_overlap

    function generate_random_mo_coeff(ao_overlap, n_mo, n_particle) result(mo_coeff)
        !
        ! this function generates random MO coefficients for every particle channel,
        ! orthonormal with respect to the AO overlap matrix, as the leading columns of
        ! a random orthogonal matrix transformed with the inverse square root of the
        ! overlap matrix
        !
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix

        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_mo, n_particle
        real(rp) :: mo_coeff(size(ao_overlap, 1), n_mo, n_particle)

        integer(ip) :: n, lwork, info, i
        real(rp) :: eigvecs(size(ao_overlap, 1), size(ao_overlap, 1)), &
                    eigvals(size(ao_overlap, 1)), &
                    inv_sqrt(size(ao_overlap, 1), size(ao_overlap, 1)), &
                    orth(size(ao_overlap, 1), size(ao_overlap, 1))
        real(rp), allocatable :: work(:)
        external :: dsyev

        ! inverse square root of the overlap matrix
        n = size(ao_overlap, 1, kind=ip)
        eigvecs = ao_overlap
        allocate(work(1))
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, -1_ip, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
        deallocate(work)
        inv_sqrt = matmul(eigvecs * spread(1.0_rp / sqrt(eigvals), 1, n), &
                          transpose(eigvecs))

        ! orthonormal orbitals of every particle channel
        do i = 1, n_particle
            orth = generate_random_orthogonal_matrix(n)
            mo_coeff(:, :, i) = matmul(inv_sqrt, orth(:, :n_mo))
        end do

    end function generate_random_mo_coeff

    subroutine setup_mo_object(mo_coeff, ao_overlap, n_occ)
        !
        ! this subroutine sets up the module-global MO object the way the MO factory
        ! would, pointing to the test-local MO coefficients
        !
        use otr_mo, only: mo_object

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:)

        integer(ip) :: n_ao, n_mo, n_particle, i

        ! dimensions
        n_ao = size(mo_coeff, 1, kind=ip)
        n_mo = size(mo_coeff, 2, kind=ip)
        n_particle = size(n_occ, kind=ip)

        ! set up the MO object
        allocate(mo_object)
        mo_object%n_ao = n_ao
        mo_object%n_mo = n_mo
        mo_object%n_particle = n_particle
        mo_object%n_param = sum(n_occ * (n_mo - n_occ))
        mo_object%mo_coeff => mo_coeff
        mo_object%ao_overlap = ao_overlap
        allocate(mo_object%mo_channels(n_particle))
        mo_object%mo_channels%n_occ = n_occ
        mo_object%mo_channels%n_virt = n_mo - n_occ
        allocate(mo_object%dm_ao(n_ao, n_ao, n_particle), &
                 mo_object%grad(mo_object%n_param), mo_object%h_diag(mo_object%n_param))
        do i = 1, n_particle
            mo_object%dm_ao(:, :, i) = matmul(mo_coeff(:, :n_occ(i), i), &
                                              transpose(mo_coeff(:, :n_occ(i), i)))
        end do

    end subroutine setup_mo_object

    subroutine setup_random_mo_channels(n_occ, n_mo)
        !
        ! this subroutine fills the occupied-occupied and virtual-virtual Fock matrix
        ! blocks of the MO channels of the MO object with random symmetric matrices
        ! together with their eigendecompositions
        !
        use otr_mo, only: mo_object
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix

        integer(ip), intent(in) :: n_occ(:), n_mo

        integer(ip) :: k

        do k = 1, size(n_occ, kind=ip)
            associate(channel => mo_object%mo_channels(k))
                call generate_random_symm_eigen(n_occ(k), channel%occ_eigvals, &
                                                channel%occ_eigvecs, channel%fock_oo)
                call generate_random_symm_eigen(n_mo - n_occ(k), channel%virt_eigvals, &
                                                channel%virt_eigvecs, channel%fock_vv)
            end associate
        end do
        mo_object%hess_eigen_stale = .false.

    contains

        subroutine generate_random_symm_eigen(n, eigvals, eigvecs, matrix)
            !
            ! this subroutine returns a random symmetric matrix of dimension n together
            ! with its eigendecomposition
            !
            integer(ip), intent(in) :: n
            real(rp), allocatable, intent(out) :: eigvals(:), eigvecs(:, :), &
                                                  matrix(:, :)

            ! an empty block has an empty eigendecomposition
            if (n == 0) then
                allocate(eigvals(0), eigvecs(0, 0), matrix(0, 0))
                return
            end if

            allocate(eigvals(n))
            call random_number(eigvals)
            eigvals = 2.0_rp * eigvals - 1.0_rp

            ! non-symmetric eigenvectors, unlike those LAPACK returns for small
            ! matrices, so that a transposed rotation cannot go unnoticed
            eigvecs = matmul(generate_random_orthogonal_matrix(n), &
                             generate_random_orthogonal_matrix(n))
            matrix = matmul(eigvecs * spread(eigvals, 1, n), transpose(eigvecs))

        end subroutine generate_random_symm_eigen

    end subroutine setup_random_mo_channels

    subroutine setup_minimal_mo_object(n_occ)
        !
        ! this subroutine sets up the module-global MO object with only the dimensions
        ! and the particle channels of the given occupations, whose Fock matrix blocks
        ! and eigendecompositions are random, which is all that the routines acting on
        ! the static part of the Hessian read
        !
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo

        integer(ip), intent(in) :: n_occ(:)

        allocate(mo_object)
        mo_object%n_mo = n_mo
        mo_object%n_particle = size(n_occ, kind=ip)
        mo_object%n_param = sum(n_occ * (n_mo - n_occ))
        allocate(mo_object%mo_channels(size(n_occ)))
        mo_object%mo_channels%n_occ = n_occ
        mo_object%mo_channels%n_virt = n_mo - n_occ
        call setup_random_mo_channels(n_occ, n_mo)

    end subroutine setup_minimal_mo_object

    subroutine setup_identity_mo_eigenbasis(occ_eigval, virt_eigval)
        !
        ! this subroutine sets up an up-to-date eigendecomposition of the static part
        ! of the Hessian in the MO basis itself for the particle channels of the MO
        ! object, with the given eigenvalues of the occupied-occupied and the
        ! virtual-virtual blocks of the Fock matrix, so that all eigenvalue pairs are
        ! the same and routines using the eigendecomposition reduce to a scaling
        !
        use otr_mo, only: mo_object
        use otr_common_unit_tests, only: identity_matrix

        real(rp), intent(in) :: occ_eigval, virt_eigval

        integer(ip) :: k

        do k = 1, mo_object%n_particle
            associate(channel => mo_object%mo_channels(k))
                channel%occ_eigvecs = identity_matrix(channel%n_occ)
                channel%occ_eigvals = spread(occ_eigval, 1, channel%n_occ)
                channel%virt_eigvecs = identity_matrix(channel%n_virt)
                channel%virt_eigvals = spread(virt_eigval, 1, channel%n_virt)
            end associate
        end do
        mo_object%hess_eigen_stale = .false.

    end subroutine setup_identity_mo_eigenbasis

    function ref_mo_transform(coeff, matrix) result(transformed)
        !
        ! this function independently reproduces the transformation C^T M C of a matrix
        ! given in the AO basis
        !
        real(rp), intent(in) :: coeff(:, :), matrix(:, :)
        real(rp) :: transformed(size(coeff, 2), size(coeff, 2))

        transformed = matmul(transpose(coeff), matmul(matrix, coeff))

    end function ref_mo_transform

    function ref_pack_ov(matrices, n_occ) result(packed)
        !
        ! this function independently packs the occupied-virtual blocks of a set of
        ! matrices in the MO basis, one per particle channel, with the occupied index
        ! running fastest
        !
        real(rp), intent(in) :: matrices(:, :, :)
        integer(ip), intent(in) :: n_occ(:)
        real(rp), allocatable :: packed(:)

        integer(ip) :: n_mo, i, j, k

        n_mo = size(matrices, 1, kind=ip)
        allocate(packed(0))
        do k = 1, size(n_occ, kind=ip)
            do j = n_occ(k) + 1, n_mo
                do i = 1, n_occ(k)
                    packed = [packed, matrices(i, j, k)]
                end do
            end do
        end do

    end function ref_pack_ov

    function ref_unpack_ov(packed, n_occ, n_mo) result(matrices)
        !
        ! this function independently unpacks a parameter vector in the MO basis into
        ! the occupied-virtual blocks of a set of otherwise vanishing matrices, one per
        ! particle channel
        !
        real(rp), intent(in) :: packed(:)
        integer(ip), intent(in) :: n_occ(:), n_mo
        real(rp) :: matrices(n_mo, n_mo, size(n_occ))

        integer(ip) :: idx, i, j, k

        matrices = 0.0_rp
        idx = 0
        do k = 1, size(n_occ, kind=ip)
            do j = n_occ(k) + 1, n_mo
                do i = 1, n_occ(k)
                    idx = idx + 1
                    matrices(i, j, k) = packed(idx)
                end do
            end do
        end do

    end function ref_unpack_ov

    function ref_hess_x_static_mo(x, channels) result(hess_x)
        !
        ! this function independently reproduces the static part of the Hessian in the
        ! MO basis, X F_vv - F_oo X for the occupied-virtual block X of every particle
        ! channel, scaled by 4 for closed-shell and 2 for open-shell systems
        !
        use otr_mo, only: mo_channel_type

        real(rp), intent(in) :: x(:)
        type(mo_channel_type), intent(in) :: channels(:)
        real(rp) :: hess_x(size(x))

        integer(ip) :: offset, n_occ, n_virt, k
        real(rp) :: shell_scale
        real(rp), allocatable :: ov_block(:, :)

        shell_scale = merge(4.0_rp, 2.0_rp, size(channels) == 1)
        offset = 0
        do k = 1, size(channels, kind=ip)
            n_occ = channels(k)%n_occ
            n_virt = channels(k)%n_virt
            ov_block = reshape(x(offset + 1:offset + n_occ * n_virt), [n_occ, n_virt])
            hess_x(offset + 1:offset + n_occ * n_virt) = shell_scale * reshape( &
                matmul(ov_block, channels(k)%fock_vv) - &
                matmul(channels(k)%fock_oo, ov_block), [n_occ * n_virt])
            offset = offset + n_occ * n_virt
        end do

    end function ref_hess_x_static_mo

    function ref_hess_x_mo(x, channels, mo_coeff, response_factor) result(hess_x)
        !
        ! this function independently reproduces the Hessian linear transformation in
        ! the MO basis for a response function returning a multiple of the density
        ! matrix displacement, by displacing the density matrix of every particle
        ! channel in the MO basis by the symmetrized occupied-virtual block of the
        ! trial vector, transforming the displacement to the AO basis as C dD C^T and
        ! the resulting response back to the MO basis, whose occupied-virtual block,
        ! scaled like the static part, is added to the static part
        !
        use otr_mo, only: mo_channel_type

        real(rp), intent(in) :: x(:), mo_coeff(:, :, :), response_factor
        type(mo_channel_type), intent(in) :: channels(:)
        real(rp) :: hess_x(size(x))

        real(rp) :: x_full(size(mo_coeff, 2), size(mo_coeff, 2), size(channels)), &
                    response_mo(size(mo_coeff, 2), size(mo_coeff, 2), size(channels)), &
                    shell_scale
        integer(ip) :: k

        shell_scale = merge(4.0_rp, 2.0_rp, size(channels) == 1)
        x_full = ref_unpack_ov(x, channels%n_occ, size(mo_coeff, 2, kind=ip))
        do k = 1, size(channels, kind=ip)
            response_mo(:, :, k) = ref_mo_transform( &
                mo_coeff(:, :, k), response_factor * matmul(mo_coeff(:, :, k), matmul( &
                    x_full(:, :, k) + transpose(x_full(:, :, k)), &
                    transpose(mo_coeff(:, :, k)))))
        end do
        hess_x = ref_hess_x_static_mo(x, channels) + &
                 shell_scale * ref_pack_ov(response_mo, channels%n_occ)

    end function ref_hess_x_mo

    function ref_rotate_eigenbasis_mo(x, channels, n_mo, to_eigenbasis) result(rotated)
        !
        ! this function independently rotates a parameter vector in the MO basis into
        ! or out of the eigenbasis of the static part of the Hessian given by the
        ! cached eigenvectors of the Fock matrix blocks of every particle channel, by
        ! rotating the occupied-virtual block X of every channel as U_o^T X U_v or as
        ! U_o X U_v^T
        !
        use otr_mo, only: mo_channel_type

        real(rp), intent(in) :: x(:)
        type(mo_channel_type), intent(in) :: channels(:)
        integer(ip), intent(in) :: n_mo
        logical, intent(in) :: to_eigenbasis
        real(rp), allocatable :: rotated(:)

        real(rp) :: x_full(n_mo, n_mo, size(channels)), &
                    rotated_full(n_mo, n_mo, size(channels))
        integer(ip) :: n_occ, k

        x_full = ref_unpack_ov(x, channels%n_occ, n_mo)
        rotated_full = 0.0_rp
        do k = 1, size(channels, kind=ip)
            n_occ = channels(k)%n_occ
            if (to_eigenbasis) then
                rotated_full(:n_occ, n_occ + 1:, k) = matmul( &
                    transpose(channels(k)%occ_eigvecs), &
                    matmul(x_full(:n_occ, n_occ + 1:, k), channels(k)%virt_eigvecs))
            else
                rotated_full(:n_occ, n_occ + 1:, k) = matmul( &
                    channels(k)%occ_eigvecs, &
                    matmul(x_full(:n_occ, n_occ + 1:, k), &
                           transpose(channels(k)%virt_eigvecs)))
            end if
        end do
        rotated = ref_pack_ov(rotated_full, channels%n_occ)

    end function ref_rotate_eigenbasis_mo

    function ref_hess_eigval_pairs_mo(channels) result(eigval_pairs)
        !
        ! this function independently constructs the eigenvalues of the static part of
        ! the Hessian in the MO basis, the scaled differences of the cached virtual and
        ! occupied Fock matrix block eigenvalues of every particle channel
        !
        use otr_mo, only: mo_channel_type

        type(mo_channel_type), intent(in) :: channels(:)
        real(rp), allocatable :: eigval_pairs(:)

        real(rp) :: shell_scale
        integer(ip) :: i, j, k

        shell_scale = merge(4.0_rp, 2.0_rp, size(channels) == 1)
        eigval_pairs = [(( &
            (shell_scale * (channels(k)%virt_eigvals(j) - channels(k)%occ_eigvals(i)), &
             i = 1, channels(k)%n_occ), j = 1, channels(k)%n_virt), k = 1, &
            size(channels, kind=ip))]

    end function ref_hess_eigval_pairs_mo

    function ref_rotate_mo_coeff(kappa, mo_coeff, n_occ) result(rot_mo_coeff)
        !
        ! this function independently reproduces the rotation of orthonormal MO
        ! coefficients of every particle channel as C exp(K)^T with K the antisymmetric
        ! matrix whose occupied-virtual block is given by kappa
        !
        real(rp), intent(in) :: kappa(:), mo_coeff(:, :, :)
        integer(ip), intent(in) :: n_occ(:)
        real(rp) :: &
            rot_mo_coeff(size(mo_coeff, 1), size(mo_coeff, 2), size(mo_coeff, 3))

        real(rp) :: kappa_full(size(mo_coeff, 2), size(mo_coeff, 2), size(n_occ))
        integer(ip) :: k

        kappa_full = ref_unpack_ov(kappa, n_occ, size(mo_coeff, 2, kind=ip))
        do k = 1, size(n_occ, kind=ip)
            rot_mo_coeff(:, :, k) = matmul(mo_coeff(:, :, k), transpose( &
                ref_expm(kappa_full(:, :, k) - transpose(kappa_full(:, :, k)))))
        end do

    contains

        function ref_expm(a) result(exp_a)
            !
            ! this function independently reproduces the exponential of a small matrix
            ! by summing its Taylor series until the terms vanish at working precision
            !
            use otr_common_unit_tests, only: identity_matrix

            real(rp), intent(in) :: a(:, :)
            real(rp) :: exp_a(size(a, 1), size(a, 1))

            real(rp) :: term(size(a, 1), size(a, 1))
            integer(ip) :: order

            exp_a = identity_matrix(size(a, 1, kind=ip))
            term = exp_a
            do order = 1, 100
                term = matmul(term, a) / real(order, kind=rp)
                exp_a = exp_a + term
                if (maxval(abs(term)) < epsilon(1.0_rp)) exit
            end do

        end function ref_expm

    end function ref_rotate_mo_coeff

    function check_rotate_hess_eigenbasis_mo(to_eigenbasis) result(passed)
        !
        ! this function checks the rotation of a random parameter vector into or out of
        ! the eigenbasis of the static part of the Hessian for the MO basis in every
        ! occupation case
        !
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_names

        logical, intent(in) :: to_eigenbasis
        logical :: passed

        real(rp), allocatable :: x(:), rotated(:)
        integer(ip) :: i_case
        character(:), allocatable :: test_name

        ! assume test passes
        passed = .true.
        test_name = trim(merge("test_rotate_to_hess_eigenbasis_mo  ", &
                               "test_rotate_from_hess_eigenbasis_mo", to_eigenbasis))

        ! call routine for every occupation case and determine if the rotation matches
        do i_case = 1, n_cases
            call setup_minimal_mo_object(case_n_occ(:case_n_particle(i_case), i_case))
            allocate(x(mo_object%n_param))
            call random_number(x)
            if (to_eigenbasis) then
                rotated = mo_object%rotate_to_hess_eigenbasis(x)
            else
                rotated = mo_object%rotate_from_hess_eigenbasis(x)
            end if
            if (norm2(rotated - ref_rotate_eigenbasis_mo( &
                x, mo_object%mo_channels, n_mo, to_eigenbasis)) > tol) then
                write (stderr, *) test_name//" failed: Incorrect rotation for the "// &
                    trim(case_names(i_case))//" case."
                passed = .false.
            end if
            deallocate(x, mo_object)
        end do

    end function check_rotate_hess_eigenbasis_mo

    logical(c_bool) function test_mo_factory_cs() bind(C)
        !
        ! this function tests the subroutine which returns the modified MO orbital
        ! updating function for the closed-shell case
        !
        use otr_mo, only: mo_factory_cs, mo_object, mo_settings_type, &
                          obj_func_mo_callback_ptr, update_orbs_mo_callback_ptr, &
                          precond_mo_callback_ptr
        use otr_common, only: evaluate_dm_cs_type
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ
        use otr_common_unit_tests, only: mock_evaluate_dm_cs, mock_evaluate_dm_os

        integer(ip), parameter :: n_particle = 1

        real(rp), target :: mo_coeff(n_ao, n_mo)
        real(rp) :: ao_overlap(n_ao, n_ao), mo_coeff_3d(n_ao, n_mo, n_particle)
        integer(ip) :: error
        type(mo_settings_type) :: settings
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_mo_funptr
        procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
        type(solver_settings_type) :: solver_settings

        ! assume tests pass
        test_mo_factory_cs = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff_3d = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        mo_coeff = mo_coeff_3d(:, :, 1)

        ! initialize callback function pointers
        evaluate_dm_funptr => mock_evaluate_dm_cs

        ! leave an MO object from an open-shell calculation behind, whose density
        ! matrix evaluating function has to be cleared
        allocate(mo_object)
        mo_object%evaluate_dm_os => mock_evaluate_dm_os

        ! call routine and determine if an error is produced
        call mo_factory_cs(mo_coeff, ao_overlap, n_occ(1), n_particle, n_ao, n_mo, &
                           evaluate_dm_funptr, obj_func_mo_funptr, &
                           update_orbs_mo_funptr, solver_settings, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_mo_factory_cs failed: Produced error."
            test_mo_factory_cs = .false.
            if (allocated(mo_object)) deallocate(mo_object)
            return
        end if

        ! determine if the MO object points to the MO coefficients of the caller as its
        ! only particle channel, the remaining common setup is covered by the test of
        ! the common setup
        if (.not. associated(mo_object%mo_coeff)) then
            write (stderr, *) "test_mo_factory_cs failed: MO coefficients not "// &
                "associated."
            test_mo_factory_cs = .false.
            deallocate(mo_object)
            return
        end if
        if (any(shape(mo_object%mo_coeff) /= [n_ao, n_mo, n_particle])) then
            write (stderr, *) "test_mo_factory_cs failed: MO coefficients not "// &
                "associated as a single particle channel."
            test_mo_factory_cs = .false.
            deallocate(mo_object)
            return
        end if
        mo_coeff(1, 1) = mo_coeff(1, 1) + 1.0_rp
        if (norm2(mo_object%mo_coeff(:, :, 1) - mo_coeff) > tol) then
            write (stderr, *) "test_mo_factory_cs failed: MO coefficients of the "// &
                "caller not taken over."
            test_mo_factory_cs = .false.
        end if
        if (.not. associated(mo_object%evaluate_dm_cs, mock_evaluate_dm_cs)) then
            write (stderr, *) "test_mo_factory_cs failed: Density matrix "// &
                "evaluating function not stored correctly."
            test_mo_factory_cs = .false.
        end if
        if (associated(mo_object%evaluate_dm_os)) then
            write (stderr, *) "test_mo_factory_cs failed: Open-shell density "// &
                "matrix evaluating function of a previous calculation kept."
            test_mo_factory_cs = .false.
        end if

        ! determine if returned function pointers point to the correct routines
        if (.not. associated(obj_func_mo_funptr, obj_func_mo_callback_ptr)) then
            write (stderr, *) "test_mo_factory_cs failed: Returned objective "// &
                "function is wrong."
            test_mo_factory_cs = .false.
        end if
        if (.not. associated(update_orbs_mo_funptr, update_orbs_mo_callback_ptr)) then
            write (stderr, *) "test_mo_factory_cs failed: Returned orbital "// &
                "updating function is wrong."
            test_mo_factory_cs = .false.
        end if

        ! determine if the MO routines are wired into the solver settings
        if (.not. associated(solver_settings%precond, precond_mo_callback_ptr)) then
            write (stderr, *) "test_mo_factory_cs failed: MO routines not wired "// &
                "into solver settings."
            test_mo_factory_cs = .false.
        end if
        deallocate(mo_object)

    end function test_mo_factory_cs

    logical(c_bool) function test_mo_factory_os() bind(C)
        !
        ! this function tests the subroutine which returns the modified MO orbital
        ! updating function for the open-shell case
        !
        use otr_mo, only: mo_factory_os, mo_object, mo_settings_type, &
                          obj_func_mo_callback_ptr, update_orbs_mo_callback_ptr, &
                          precond_mo_callback_ptr
        use otr_common, only: evaluate_dm_os_type
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_particle, n_occ
        use otr_common_unit_tests, only: mock_evaluate_dm_cs, mock_evaluate_dm_os

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: error
        type(mo_settings_type) :: settings
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_mo_funptr
        procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
        type(solver_settings_type) :: solver_settings

        ! assume tests pass
        test_mo_factory_os = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)

        ! initialize callback function pointers
        evaluate_dm_funptr => mock_evaluate_dm_os

        ! leave an MO object from a closed-shell calculation behind, whose density
        ! matrix evaluating function has to be cleared
        allocate(mo_object)
        mo_object%evaluate_dm_cs => mock_evaluate_dm_cs

        ! call routine and determine if an error is produced
        call mo_factory_os(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                           evaluate_dm_funptr, obj_func_mo_funptr, &
                           update_orbs_mo_funptr, solver_settings, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_mo_factory_os failed: Produced error."
            test_mo_factory_os = .false.
            if (allocated(mo_object)) deallocate(mo_object)
            return
        end if

        ! determine if the MO object points to the MO coefficients of the caller, the
        ! remaining common setup is covered by the test of the common setup
        if (.not. associated(mo_object%mo_coeff, mo_coeff)) then
            write (stderr, *) "test_mo_factory_os failed: MO coefficients of the "// &
                "caller not taken over."
            test_mo_factory_os = .false.
        end if
        if (.not. associated(mo_object%evaluate_dm_os, mock_evaluate_dm_os)) then
            write (stderr, *) "test_mo_factory_os failed: Density matrix "// &
                "evaluating function not stored correctly."
            test_mo_factory_os = .false.
        end if
        if (associated(mo_object%evaluate_dm_cs)) then
            write (stderr, *) "test_mo_factory_os failed: Closed-shell density "// &
                "matrix evaluating function of a previous calculation kept."
            test_mo_factory_os = .false.
        end if

        ! determine if returned function pointers point to the correct routines
        if (.not. associated(obj_func_mo_funptr, obj_func_mo_callback_ptr)) then
            write (stderr, *) "test_mo_factory_os failed: Returned objective "// &
                "function is wrong."
            test_mo_factory_os = .false.
        end if
        if (.not. associated(update_orbs_mo_funptr, update_orbs_mo_callback_ptr)) then
            write (stderr, *) "test_mo_factory_os failed: Returned orbital "// &
                "updating function is wrong."
            test_mo_factory_os = .false.
        end if

        ! determine if the MO routines are wired into the solver settings
        if (.not. associated(solver_settings%precond, precond_mo_callback_ptr)) then
            write (stderr, *) "test_mo_factory_os failed: MO routines not wired "// &
                "into solver settings."
            test_mo_factory_os = .false.
        end if
        deallocate(mo_object)

    end function test_mo_factory_os

    logical(c_bool) function test_mo_factory_common() bind(C)
        !
        ! this function tests the subroutine which performs the common MO
        ! initialization operations, which sets the MO object up anew unless it was
        ! already set up for the same dimensions and occupations, for the closed- and
        ! the open-shell case
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: mo_factory_common, mo_object
        use otr_mo_test_reference, only: n_mo, n_param_cs, n_param_os
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, operator(==)
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_unit_tests, only: shell_names, mock_get_response_cs, &
                                         mock_get_response_os, mock_evaluate_dm_cs, &
                                         mock_evaluate_dm_os

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref), &
                            mo_coeff_new(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: n_particle, leftover, n_ao_old, n_mo_old, error
        integer(ip), allocatable :: n_occ_old(:)
        type(orbital_settings_type) :: settings
        character(:), allocatable :: shell, leftover_case
        character(39), parameter :: leftover_names(4) = [ &
            character(39) :: "a different number of AOs", "a different number of MOs", &
            "a different number of particle channels", "different occupations"]

        ! assume tests pass
        test_mo_factory_common = .true.

        ! setup settings object
        call setup_settings(settings)

        do n_particle = 1, n_particle_ref
            shell = trim(shell_names(n_particle))

            ! initialize random orthonormal MO coefficients
            ao_overlap = generate_random_ao_overlap(n_ao)
            mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
            mo_coeff_new = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)

            ! leave MO objects from calculations which differ only in the number of
            ! AOs, of MOs or of particle channels or in the occupations, which are
            ! lowered for the closed and swapped for the open shell, behind and
            ! determine if each is set up anew
            do leftover = 1, size(leftover_names, kind=ip)
                n_ao_old = merge(n_ao + 1, n_ao, leftover == 1)
                n_mo_old = merge(n_mo - 1, n_mo, leftover == 2)
                if (leftover == 3) then
                    n_occ_old = n_occ(:3 - n_particle)
                else if (leftover == 4) then
                    n_occ_old = merge(n_occ(:n_particle) - 1, n_occ(n_particle:1:-1), &
                                      n_particle == 1)
                else
                    n_occ_old = n_occ(:n_particle)
                end if
                leftover_case = shell//" object left from a calculation with "// &
                                trim(leftover_names(leftover))
                if (allocated(mo_object)) deallocate(mo_object)
                allocate(mo_object)
                mo_object%n_ao = n_ao_old
                mo_object%n_mo = n_mo_old
                allocate(mo_object%mo_channels(size(n_occ_old)), mo_object%grad(1), &
                         mo_object%h_diag(1), &
                         mo_object%dm_ao(n_ao_old, n_ao_old, size(n_occ_old)))
                mo_object%mo_channels%n_occ = n_occ_old
                call mo_factory_common(mo_coeff(:, :, :n_particle), ao_overlap, &
                                       n_occ(:n_particle), n_particle, n_ao, n_mo, &
                                       error, settings)
                if (error /= 0) then
                    write (stderr, *) "test_mo_factory_common failed: Produced "// &
                        "error for the "//leftover_case//"."
                    test_mo_factory_common = .false.
                end if
                if (.not. check_mo_object(mo_coeff(:, :, :n_particle), leftover_case)) &
                    test_mo_factory_common = .false.
            end do

            ! leave evaluated quantities with their response functions, the density
            ! matrix evaluating function and only the gradient allocated behind and
            ! call routine again with new settings and starting orbitals, and determine
            ! if the object is kept while the new settings and orbitals are taken over,
            ! its evaluated state is discarded, the density matrix evaluating function,
            ! which only the factories set, is kept and the Hessian diagonal is
            ! allocated
            mo_object%evaluation_stale = .false.
            mo_object%response_stale = .false.
            mo_object%hess_eigen_stale = .false.
            mo_object%get_response_cs => mock_get_response_cs
            mo_object%get_response_os => mock_get_response_os
            if (n_particle == 1) then
                mo_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                mo_object%evaluate_dm_os => mock_evaluate_dm_os
            end if
            mo_object%mo_channels(1)%fock_oo = reshape([7.0_rp], [1, 1])
            deallocate(mo_object%h_diag)
            settings%verbose = settings%verbose + 1
            call mo_factory_common(mo_coeff_new(:, :, :n_particle), ao_overlap, &
                                   n_occ(:n_particle), n_particle, n_ao, n_mo, error, &
                                   settings)
            if (error /= 0) then
                write (stderr, *) "test_mo_factory_common failed: Produced error "// &
                    "for the reused "//shell//" object."
                test_mo_factory_common = .false.
            end if
            if (.not. check_mo_object(mo_coeff_new(:, :, :n_particle), &
                                      "reused "//shell//" object")) &
                test_mo_factory_common = .false.
            if (.not. allocated(mo_object%mo_channels(1)%fock_oo)) then
                write (stderr, *) "test_mo_factory_common failed: MO object not "// &
                    "kept for the same dimensions and occupations for the "//shell// &
                    " object."
                test_mo_factory_common = .false.
            end if
            if (.not. (associated(mo_object%evaluate_dm_cs, mock_evaluate_dm_cs) .or. &
                       associated(mo_object%evaluate_dm_os, mock_evaluate_dm_os))) then
                write (stderr, *) "test_mo_factory_common failed: Density matrix "// &
                    "evaluating function not kept for the reused "//shell//" object."
                test_mo_factory_common = .false.
            end if

            ! leave an MO object behind whose previous setup failed before its particle
            ! channels were set up and determine if it is set up anew although the
            ! dimensions match
            deallocate(mo_object%mo_channels)
            call mo_factory_common(mo_coeff(:, :, :n_particle), ao_overlap, &
                                   n_occ(:n_particle), n_particle, n_ao, n_mo, error, &
                                   settings)
            if (error /= 0) then
                write (stderr, *) "test_mo_factory_common failed: Produced error "// &
                    "for the "//shell//" object whose previous setup failed."
                test_mo_factory_common = .false.
            end if
            if (.not. check_mo_object(mo_coeff(:, :, :n_particle), &
                                      shell//" object whose previous setup failed")) &
                test_mo_factory_common = .false.

            ! call routine for dimensions which do not match the MO coefficients and
            ! determine if the sanity check rejects them before the object is changed
            mo_object%hess_eigen_stale = .false.
            call mo_factory_common(mo_coeff(:, :, :n_particle), ao_overlap, &
                                   n_occ(:n_particle), n_particle, n_ao, n_mo + 1, &
                                   error, settings)
            if (error == 0) then
                write (stderr, *) "test_mo_factory_common failed: Error not thrown "// &
                    "for mismatching dimensions for the "//shell//" object."
                test_mo_factory_common = .false.
            end if
            if (mo_object%hess_eigen_stale) then
                write (stderr, *) "test_mo_factory_common failed: MO object "// &
                    "changed for mismatching dimensions for the "//shell//" object."
                test_mo_factory_common = .false.
            end if

            ! deallocate MO object
            deallocate(mo_object)
        end do

    contains

        function check_mo_object(caller_mo_coeff, case_name) result(passed)
            !
            ! this function checks if the MO object is set up for the given MO
            ! coefficients of the caller and the AO overlap matrix and occupations of
            ! the current shell, with the density matrix constructed from the occupied
            ! orbitals, the evaluated state with its response functions discarded and
            ! the settings stored
            !
            real(rp), intent(in), target :: caller_mo_coeff(:, :, :)
            character(*), intent(in) :: case_name
            logical :: passed

            integer(ip) :: n_param, k

            ! assume test passes
            passed = .true.

            ! the checks below compare arrays of these dimensions
            n_param = merge(n_param_cs, n_param_os, n_particle == 1)
            if (any([mo_object%n_ao, mo_object%n_mo, mo_object%n_particle, &
                     mo_object%n_param, size(mo_object%mo_channels, kind=ip)] /= &
                    [n_ao, n_mo, n_particle, n_param, n_particle])) then
                write (stderr, *) "test_mo_factory_common failed: Incorrect "// &
                    "dimensions for the "//case_name//"."
                passed = .false.
                return
            end if

            ! determine if the occupations, the AO overlap matrix, the MO coefficients
            ! of the caller, the density matrix, the gradient and Hessian diagonal and
            ! the discarded evaluated state are set up
            if (any(mo_object%mo_channels%n_occ /= n_occ(:n_particle)) .or. &
                any(mo_object%mo_channels%n_virt /= n_mo - n_occ(:n_particle))) then
                write (stderr, *) "test_mo_factory_common failed: Incorrect "// &
                    "occupations for the "//case_name//"."
                passed = .false.
            end if
            if (norm2(mo_object%ao_overlap - ao_overlap) > tol) then
                write (stderr, *) "test_mo_factory_common failed: AO overlap "// &
                    "matrix not stored for the "//case_name//"."
                passed = .false.
            end if
            if (.not. associated(mo_object%mo_coeff, caller_mo_coeff)) then
                write (stderr, *) "test_mo_factory_common failed: MO coefficients "// &
                    "of the caller not taken over for the "//case_name//"."
                passed = .false.
            end if
            do k = 1, n_particle
                if (norm2(mo_object%dm_ao(:, :, k) - matmul( &
                    caller_mo_coeff(:, :n_occ(k), k), &
                    transpose(caller_mo_coeff(:, :n_occ(k), k)))) > tol) then
                    write (stderr, *) "test_mo_factory_common failed: Density "// &
                        "matrix not constructed from the occupied orbitals for the "// &
                        case_name//"."
                    passed = .false.
                end if
            end do
            if (.not. (allocated(mo_object%grad) .and. allocated(mo_object%h_diag))) &
                then
                write (stderr, *) "test_mo_factory_common failed: Gradient or "// &
                    "Hessian diagonal not allocated for the "//case_name//"."
                passed = .false.
            else if (any([size(mo_object%grad), size(mo_object%h_diag)] /= n_param)) &
                then
                write (stderr, *) "test_mo_factory_common failed: Gradient or "// &
                    "Hessian diagonal not allocated with the number of parameters "// &
                    "for the "//case_name//"."
                passed = .false.
            end if
            if (.not. (mo_object%evaluation_stale .and. mo_object%response_stale .and. &
                       mo_object%hess_eigen_stale)) then
                write (stderr, *) "test_mo_factory_common failed: Evaluated state "// &
                    "not discarded for the "//case_name//"."
                passed = .false.
            end if
            if (associated(mo_object%get_response_cs) .or. &
                associated(mo_object%get_response_os)) then
                write (stderr, *) "test_mo_factory_common failed: Response "// &
                    "function kept for the "//case_name//"."
                passed = .false.
            end if
            if (.not. (mo_object%settings == settings)) then
                write (stderr, *) "test_mo_factory_common failed: Settings not "// &
                    "stored for the "//case_name//"."
                passed = .false.
            end if

        end function check_mo_object

    end function test_mo_factory_common

    logical(c_bool) function test_mo_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! dimensions of orbitals parameterized in the MO basis
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: mo_sanity_check
        use opentrustregion_unit_tests, only: setup_settings, log_message

        type :: mo_sanity_case_type
            integer(ip) :: mo_coeff_shape(3), ao_overlap_shape(2), n_occ(2), &
                           n_occ_passed, n_particle, n_ao, n_mo
            logical :: valid
            character(54) :: name
            character(80) :: message
        end type mo_sanity_case_type

        type(mo_sanity_case_type), parameter :: cases(*) = [ &
            ! mo_coeff_shape, ao_overlap_shape, n_occ, n_occ_passed, n_particle, n_ao,
            ! n_mo, valid, then name and message
            mo_sanity_case_type( &
                [4, 3, 1], [4, 4], [1, 0], 1, 1, 4, 3, .true., &
                "valid closed-shell dimensions with fewer MOs than AOs", ""), &
            mo_sanity_case_type( &
                [3, 3, 2], [3, 3], [2, 0], 2, 2, 3, 3, .true., &
                "valid open-shell dimensions with as many MOs as AOs", ""), &
            mo_sanity_case_type( &
                [4, 3, 2], [4, 4], [3, 1], 2, 2, 4, 3, .true., &
                "valid open-shell dimensions with fewer MOs than AOs", ""), &
            mo_sanity_case_type( &
                [0, 1, 1], [0, 0], [1, 0], 1, 1, 0, 1, .false., &
                "vanishing number of AOs", "Number of AOs should be larger than 0."), &
            mo_sanity_case_type( &
                [3, 0, 1], [3, 3], [0, 0], 1, 1, 3, 0, .false., &
                "vanishing number of MOs", &
                "Number of MOs should be larger than 0 and not larger than the "// &
                "number of AOs."), &
            mo_sanity_case_type( &
                [3, 4, 1], [3, 3], [1, 0], 1, 1, 3, 4, .false., &
                "more MOs than AOs", ""), &
            mo_sanity_case_type( &
                [3, 3, 0], [3, 3], [1, 0], 0, 0, 3, 3, .false., &
                "vanishing number of particle channels", &
                "Number of particles should be 1 or 2."), &
            mo_sanity_case_type( &
                [3, 3, 3], [3, 3], [1, 1], 2, 3, 3, 3, .false., &
                "three particle channels", "Number of particles should be 1 or 2."), &
            mo_sanity_case_type( &
                [3, 3, 2], [3, 3], [1, 0], 1, 2, 3, 3, .false., &
                "single occupation for two particle channels", ""), &
            mo_sanity_case_type( &
                [3, 3, 1], [3, 3], [1, 1], 2, 1, 3, 3, .false., &
                "two occupations for a single particle channel", ""), &
            mo_sanity_case_type( &
                [2, 3, 1], [3, 3], [1, 0], 1, 1, 3, 3, .false., &
                "MO coefficients with wrong number of AOs", ""), &
            mo_sanity_case_type( &
                [3, 2, 1], [3, 3], [1, 0], 1, 1, 3, 3, .false., &
                "MO coefficients with wrong number of MOs", ""), &
            mo_sanity_case_type( &
                [3, 3, 2], [3, 3], [1, 0], 1, 1, 3, 3, .false., &
                "MO coefficients with wrong number of particle channels", ""), &
            mo_sanity_case_type( &
                [3, 3, 1], [3, 2], [1, 0], 1, 1, 3, 3, .false., &
                "wrongly shaped AO overlap matrix", ""), &
            mo_sanity_case_type( &
                [6, 6, 2], [6, 6], [2, -1], 2, 2, 6, 6, .false., &
                "negative number of occupied orbitals", ""), &
            mo_sanity_case_type( &
                [6, 6, 2], [6, 6], [3, 7], 2, 2, 6, 6, .false., &
                "more occupied orbitals than MOs", ""), &
            mo_sanity_case_type( &
                [3, 3, 2], [3, 3], [3, 0], 2, 2, 3, 3, .false., &
                "missing occupied-virtual rotations", "")]

        type(mo_sanity_case_type) :: c
        type(orbital_settings_type) :: settings
        integer(ip) :: i_case, error

        ! assume tests pass
        test_mo_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! call routine for every case and determine if valid dimensions are accepted,
        ! invalid ones rejected and, where the case expects one, the right check fires
        do i_case = 1, size(cases)
            c = cases(i_case)
            log_message = ""
            call mo_sanity_check(settings, c%mo_coeff_shape, c%ao_overlap_shape, &
                                 c%n_occ(:c%n_occ_passed), c%n_particle, c%n_ao, &
                                 c%n_mo, error)
            if ((error == 0) .neqv. c%valid) then
                write (stderr, *) "test_mo_sanity_check failed: Incorrect result "// &
                    "for "//trim(c%name)//"."
                test_mo_sanity_check = .false.
            end if
            if (len_trim(c%message) > 0 .and. adjustl(log_message) /= c%message) then
                write (stderr, *) "test_mo_sanity_check failed: Incorrect log "// &
                    "message for "//trim(c%name)//"."
                test_mo_sanity_check = .false.
            end if
        end do

    end function test_mo_sanity_check

    logical(c_bool) function test_mo_set_solver_settings() bind(C)
        !
        ! this function tests the subroutine which wires the MO preconditioners and
        ! extra trial vectors into the solver settings
        !
        use otr_mo, only: mo_set_solver_settings, precond_mo_callback_ptr, &
                          precond_pd_mo_callback_ptr, &
                          get_extra_trial_vectors_mo_callback_ptr
        use opentrustregion, only: solver_settings_type, default_solver_settings
        use test_reference, only: ref_settings, assignment(=), operator(/=)

        type(solver_settings_type) :: solver_settings
        integer(ip) :: i_case, error
        character(:), allocatable :: case_name

        ! assume tests pass
        test_mo_set_solver_settings = .true.

        ! wire the MO routines into uninitialized settings, which have to be
        ! initialized to their defaults, and into settings initialized to the reference
        ! values with a projection of the caller, which have to be kept
        do i_case = 1, 2
            if (i_case == 1) then
                case_name = "for uninitialized settings"
            else
                case_name = "for initialized settings"
                solver_settings = ref_settings
                solver_settings%project => mock_project
                solver_settings%stability_settings%project => mock_project
            end if
            call mo_set_solver_settings(solver_settings, error)
            if (error /= 0) then
                write (stderr, *) "test_mo_set_solver_settings failed: Produced "// &
                    "error " // case_name // "."
                test_mo_set_solver_settings = .false.
            end if
            if (.not. solver_settings%initialized) then
                write (stderr, *) "test_mo_set_solver_settings failed: Settings "// &
                    "not initialized " // case_name // "."
                test_mo_set_solver_settings = .false.
            end if
            if (i_case == 1 .and. solver_settings /= default_solver_settings) then
                write (stderr, *) "test_mo_set_solver_settings failed: Settings "// &
                    "not set to their defaults " // case_name // "."
                test_mo_set_solver_settings = .false.
            end if
            if (i_case == 2 .and. solver_settings /= ref_settings) then
                write (stderr, *) "test_mo_set_solver_settings failed: Settings "// &
                    "not kept " // case_name // "."
                test_mo_set_solver_settings = .false.
            end if
            if (.not. (associated( &
                solver_settings%precond, precond_mo_callback_ptr) .and. associated( &
                    solver_settings%precond_pd, precond_pd_mo_callback_ptr) .and. &
                associated(solver_settings%get_extra_trial_vectors, &
                           get_extra_trial_vectors_mo_callback_ptr))) then
                write (stderr, *) "test_mo_set_solver_settings failed: MO routines "// &
                    "not wired into solver settings " // case_name // "."
                test_mo_set_solver_settings = .false.
            end if
            if (.not. ( &
                associated(solver_settings%stability_settings%precond, &
                           precond_mo_callback_ptr) .and. &
                associated(solver_settings%stability_settings%get_extra_trial_vectors, &
                           get_extra_trial_vectors_mo_callback_ptr))) then
                write (stderr, *) "test_mo_set_solver_settings failed: MO routines "// &
                    "not wired into stability check settings " // case_name // "."
                test_mo_set_solver_settings = .false.
            end if
            if (i_case == 1 .and. &
                (associated(solver_settings%project) .or. &
                 associated(solver_settings%stability_settings%project))) then
                write (stderr, *) "test_mo_set_solver_settings failed: Projection "// &
                    "wired into settings " // case_name // "."
                test_mo_set_solver_settings = .false.
            end if
            if (i_case == 2 .and. .not. ( &
                associated(solver_settings%project, mock_project) .and. associated( &
                    solver_settings%stability_settings%project, mock_project))) then
                write (stderr, *) "test_mo_set_solver_settings failed: Projection "// &
                    "of the caller not kept " // case_name // "."
                test_mo_set_solver_settings = .false.
            end if
        end do

    contains

        subroutine mock_project(vector, error_out)
            !
            ! this subroutine is a mock projection of the caller, which leaves the
            ! vector unchanged
            !
            real(rp), intent(inout), target :: vector(:)
            integer(ip), intent(out) :: error_out

            error_out = 0
            vector = vector

        end subroutine mock_project

    end function test_mo_set_solver_settings

    logical(c_bool) function test_obj_func_mo_callback() bind(C)
        !
        ! this function tests the function which defines the energy evaluation in the
        ! MO basis
        !
        use otr_mo, only: obj_func_mo_callback, mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle_ref => n_particle
        use otr_common_unit_tests, only: mock_requests, shell_names, &
                                         mock_evaluate_dm_cs, mock_evaluate_dm_os

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), &
                    mo_coeff_start(n_ao, n_mo, n_particle_ref), energy, expected_energy
        real(rp), allocatable :: kappa(:), rot_mo_coeff(:, :, :), dm_ao_start(:, :, :)
        integer(ip) :: n_particle, i, error
        character(:), allocatable :: case_name

        ! assume tests pass
        test_obj_func_mo_callback = .true.

        ! initialize random orthonormal MO coefficients in a non-orthonormal AO basis
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
        mo_coeff_start = mo_coeff

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))

            ! set up the MO object
            call setup_mo_object(mo_coeff(:, :, :n_particle), ao_overlap, &
                                 n_occ(:n_particle))
            call setup_settings(mo_object%settings)
            if (n_particle == 1) then
                mo_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                mo_object%evaluate_dm_os => mock_evaluate_dm_os
            end if
            allocate(kappa(mo_object%n_param))

            ! call routine without an orbital rotation and determine if only the energy
            ! of the density matrix evaluating function is requested and returned for
            ! the current density matrix
            kappa = 0.0_rp
            mock_requests = [integer(ip) ::]
            energy = obj_func_mo_callback(kappa, error)
            if (error /= 0) then
                write (stderr, *) "test_obj_func_mo_callback failed: Produced "// &
                    "error without an orbital rotation for the "//case_name//" case."
                test_obj_func_mo_callback = .false.
            end if
            if (abs(energy - sum(mo_object%dm_ao)) > tol) then
                write (stderr, *) "test_obj_func_mo_callback failed: Incorrect "// &
                    "energy without an orbital rotation for the "//case_name//" case."
                test_obj_func_mo_callback = .false.
            end if
            if (size(mock_requests) /= 1) then
                write (stderr, *) "test_obj_func_mo_callback failed: Density "// &
                    "matrix evaluating function not called exactly once for the "// &
                    case_name//" case."
                test_obj_func_mo_callback = .false.
            end if
            if (any(mock_requests /= 0)) then
                write (stderr, *) "test_obj_func_mo_callback failed: Incorrect "// &
                    "outputs requested from density matrix evaluating function for "// &
                    "the "//case_name//" case."
                test_obj_func_mo_callback = .false.
            end if

            ! call routine with an orbital rotation and determine if only the energy of
            ! the density matrix of the rotated orbitals is requested and returned
            ! while the current orbitals and density matrix are left untouched
            call random_number(kappa)
            kappa = 0.2_rp * (kappa - 0.5_rp)
            rot_mo_coeff = ref_rotate_mo_coeff( &
                kappa, mo_coeff_start(:, :, :n_particle), n_occ(:n_particle))
            expected_energy = 0.0_rp
            do i = 1, n_particle
                expected_energy = expected_energy + &
                                  sum(matmul(rot_mo_coeff(:, :n_occ(i), i), &
                                             transpose(rot_mo_coeff(:, :n_occ(i), i))))
            end do
            dm_ao_start = mo_object%dm_ao
            mock_requests = [integer(ip) ::]
            energy = obj_func_mo_callback(kappa, error)
            if (error /= 0) then
                write (stderr, *) "test_obj_func_mo_callback failed: Produced "// &
                    "error for an orbital rotation for the "//case_name//" case."
                test_obj_func_mo_callback = .false.
            end if
            if (abs(energy - expected_energy) > tol) then
                write (stderr, *) "test_obj_func_mo_callback failed: Incorrect "// &
                    "energy for an orbital rotation for the "//case_name//" case."
                test_obj_func_mo_callback = .false.
            end if
            if (size(mock_requests) /= 1 .or. any(mock_requests /= 0)) then
                write (stderr, *) "test_obj_func_mo_callback failed: Incorrect "// &
                    "outputs requested from density matrix evaluating function for "// &
                    "an orbital rotation for the "//case_name//" case."
                test_obj_func_mo_callback = .false.
            end if
            if (norm2(mo_coeff - mo_coeff_start) > tol) then
                write (stderr, *) "test_obj_func_mo_callback failed: Current "// &
                    "orbitals changed by an orbital rotation for the "//case_name// &
                    " case."
                test_obj_func_mo_callback = .false.
            end if
            if (norm2(mo_object%dm_ao - dm_ao_start) > tol) then
                write (stderr, *) "test_obj_func_mo_callback failed: Current "// &
                    "density matrix changed by an orbital rotation for the "// &
                    case_name//" case."
                test_obj_func_mo_callback = .false.
            end if

            ! deallocate MO object
            deallocate(mo_object, kappa)
        end do

    end function test_obj_func_mo_callback

    logical(c_bool) function test_update_orbs_mo_callback() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation in the MO basis, for the closed-shell and the
        ! open-shell case
        !
        use otr_mo, only: update_orbs_mo_callback, mo_object, hess_x_mo_callback_ptr
        use opentrustregion, only: hess_x_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle_ref => n_particle
        use otr_common_unit_tests, only: &
            mock_fock_factor, mock_requests, shell_names, mock_get_response_cs, &
            mock_get_response_os, mock_evaluate_dm_cs, mock_evaluate_dm_os, &
            mock_evaluate_dm_failing_cs, mock_evaluate_dm_failing_os

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), fock_mo(n_mo, n_mo), func
        real(rp), allocatable :: mo_coeff_start(:, :, :), kappa(:), grad(:), h_diag(:)
        integer(ip) :: n_particle, i, error
        logical :: response_set
        character(:), allocatable :: case_name
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_update_orbs_mo_callback = .true.

        ! initialize a non-orthonormal AO basis
        ao_overlap = generate_random_ao_overlap(n_ao)

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))

            ! set up the MO object with random orthonormal MO coefficients
            mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
            call setup_mo_object(mo_coeff(:, :, :n_particle), ao_overlap, &
                                 n_occ(:n_particle))
            call setup_settings(mo_object%settings)
            allocate(kappa(mo_object%n_param), grad(mo_object%n_param), &
                     h_diag(mo_object%n_param))
            if (n_particle == 1) then
                mo_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                mo_object%evaluate_dm_os => mock_evaluate_dm_os
            end if
            mock_requests = [integer(ip) ::]

            ! call routine without an orbital rotation for starting orbitals which have
            ! not been evaluated yet, with the response marked current so that only the
            ! evaluation flag forces the evaluation, and determine if the energy of the
            ! density matrix, the static part of the Hessian from the Fock matrix
            ! transformed to the MO basis, the gradient, the Hessian diagonal, the
            ! response function and the Hessian linear transformation are obtained
            kappa = 0.0_rp
            mo_object%response_stale = .false.
            if (.not. check_update_orbs_mo_stage(1_ip, "starting orbitals")) &
                test_update_orbs_mo_callback = .false.
            if (abs(func - sum(mo_object%dm_ao)) > tol) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Incorrect "// &
                    "energy for the "//case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if
            do i = 1, n_particle
                fock_mo = ref_mo_transform( &
                    mo_coeff(:, :, i), mock_fock_factor(1) * mo_object%dm_ao(:, :, i))
                if (norm2(mo_object%mo_channels(i)%fock_oo - &
                          fock_mo(:n_occ(i), :n_occ(i))) > tol .or. &
                    norm2(mo_object%mo_channels(i)%fock_vv - &
                          fock_mo(n_occ(i) + 1:, n_occ(i) + 1:)) > tol) then
                    write (stderr, *) "test_update_orbs_mo_callback failed: Static "// &
                        "Hessian part not built from the Fock matrix in the MO "// &
                        "basis for the "//case_name//" case."
                    test_update_orbs_mo_callback = .false.
                end if
            end do
            if (norm2(grad - mo_object%grad) > tol) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Gradient "// &
                    "not returned for the "//case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if
            if (norm2(h_diag - mo_object%h_diag) > tol) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Hessian "// &
                    "diagonal not returned for the "//case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if
            if (n_particle == 1) then
                response_set = &
                    associated(mo_object%get_response_cs, mock_get_response_cs)
            else
                response_set = &
                    associated(mo_object%get_response_os, mock_get_response_os)
            end if
            if (.not. response_set) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Response "// &
                    "function not stored for the "//case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if
            if (.not. associated(hess_x_funptr, hess_x_mo_callback_ptr)) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Returned "// &
                    "Hessian linear transformation is wrong for the "//case_name// &
                    " case."
                test_update_orbs_mo_callback = .false.
            end if
            if (mo_object%evaluation_stale) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Evaluation "// &
                    "still marked stale after the quantities were computed for the "// &
                    case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if

            ! call routine again without an orbital rotation after clearing the outputs
            ! and setting a vanishing stored energy, and determine if the quantities of
            ! the evaluated orbitals are returned without an evaluation, even for the
            ! vanishing energy
            mo_object%energy = 0.0_rp
            grad = 0.0_rp
            h_diag = 0.0_rp
            if (.not. check_update_orbs_mo_stage(1_ip, "evaluated orbitals")) &
                test_update_orbs_mo_callback = .false.
            if (abs(func) > tol) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Stored "// &
                    "energy not returned for evaluated orbitals for the "//case_name// &
                    " case."
                test_update_orbs_mo_callback = .false.
            end if
            if (norm2(grad - mo_object%grad) > tol) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Stored "// &
                    "gradient not returned for evaluated orbitals for the "// &
                    case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if
            if (norm2(h_diag - mo_object%h_diag) > tol) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Stored "// &
                    "Hessian diagonal not returned for evaluated orbitals for the "// &
                    case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if

            ! mark the response as stale, as an approximate-Hessian extension such as
            ! ARH would after moving the orbitals without going through this routine,
            ! and call again without an orbital rotation, and determine if this forces
            ! an evaluation and clears the flag
            mo_object%response_stale = .true.
            if (.not. check_update_orbs_mo_stage(2_ip, "a stale response")) &
                test_update_orbs_mo_callback = .false.
            if (mo_object%response_stale) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Response "// &
                    "still marked stale after being recomputed for the "//case_name// &
                    " case."
                test_update_orbs_mo_callback = .false.
            end if

            ! mark the evaluation as stale, as the factory does for new starting
            ! orbitals, and call again without an orbital rotation, and determine if
            ! this forces an evaluation and clears the flag
            mo_object%evaluation_stale = .true.
            if (.not. check_update_orbs_mo_stage(3_ip, "a stale evaluation")) &
                test_update_orbs_mo_callback = .false.
            if (mo_object%evaluation_stale) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Evaluation "// &
                    "still marked stale after being recomputed for the "//case_name// &
                    " case."
                test_update_orbs_mo_callback = .false.
            end if

            ! call routine with an orbital rotation whose evaluation fails and
            ! determine if the error is passed on and the rotated orbitals are marked
            ! as not evaluated
            if (n_particle == 1) then
                mo_object%evaluate_dm_cs => mock_evaluate_dm_failing_cs
            else
                mo_object%evaluate_dm_os => mock_evaluate_dm_failing_os
            end if
            kappa = 0.1_rp
            call update_orbs_mo_callback(kappa, func, grad, h_diag, hess_x_funptr, &
                                         error)
            if (error == 0) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Error of "// &
                    "the density matrix evaluating function not passed on for the "// &
                    case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if
            if (.not. mo_object%evaluation_stale) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Evaluation "// &
                    "not marked stale after a failed evaluation for the "//case_name// &
                    " case."
                test_update_orbs_mo_callback = .false.
            end if
            if (n_particle == 1) then
                mo_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                mo_object%evaluate_dm_os => mock_evaluate_dm_os
            end if

            ! call routine with an orbital rotation, keeping the MO coefficients it
            ! starts from, and determine if the orbitals are rotated, consistently with
            ! the density matrix, and evaluated, and if every evaluation requested the
            ! Fock matrix and the response function
            mo_coeff_start = mo_object%mo_coeff
            if (.not. check_update_orbs_mo_stage(4_ip, "an orbital rotation")) &
                test_update_orbs_mo_callback = .false.
            if (norm2(mo_object%mo_coeff - ref_rotate_mo_coeff( &
                kappa, mo_coeff_start, n_occ(:n_particle))) > tol) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Orbitals "// &
                    "not rotated correctly by an orbital rotation for the "// &
                    case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if
            do i = 1, n_particle
                if (norm2(mo_object%dm_ao(:, :, i) - matmul( &
                    mo_object%mo_coeff(:, :n_occ(i), i), &
                    transpose(mo_object%mo_coeff(:, :n_occ(i), i)))) > tol) then
                    write (stderr, *) "test_update_orbs_mo_callback failed: "// &
                        "Density matrix not moved consistently with the orbitals "// &
                        "for the "//case_name//" case."
                    test_update_orbs_mo_callback = .false.
                end if
            end do
            if (any(mock_requests /= 3)) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Incorrect "// &
                    "outputs requested from density matrix evaluating function for "// &
                    "the "//case_name//" case."
                test_update_orbs_mo_callback = .false.
            end if

            ! deallocate MO object and case arrays
            deallocate(mo_object, kappa, grad, h_diag)
        end do

    contains

        function check_update_orbs_mo_stage(n_calls, stage) result(passed)
            !
            ! this function calls the MO orbital updating function at the current
            ! rotation for one stage of its test and determines if it produces no error
            ! and if the density matrix evaluating function has been called the
            ! expected number of times in total
            !
            integer(ip), intent(in) :: n_calls
            character(*), intent(in) :: stage
            logical :: passed

            ! assume the stage passes
            passed = .true.

            ! call routine and determine if it produces an error and if the density
            ! matrix evaluating function was called the expected number of times
            call update_orbs_mo_callback(kappa, func, grad, h_diag, hess_x_funptr, &
                                         error)
            if (error /= 0) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Produced "// &
                    "error for "//stage//" for the "//case_name//" case."
                passed = .false.
            end if
            if (size(mock_requests) /= n_calls) then
                write (stderr, *) "test_update_orbs_mo_callback failed: Density "// &
                    "matrix evaluating function not called the expected number of "// &
                    "times for "//stage//" for the "//case_name//" case."
                passed = .false.
            end if

        end function check_update_orbs_mo_stage

    end function test_update_orbs_mo_callback

    logical(c_bool) function test_hess_x_mo_callback() bind(C)
        !
        ! this function tests the subroutine which defines the Hessian linear
        ! transformation in the MO basis
        !
        use otr_mo, only: hess_x_mo_callback, mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_names
        use otr_common_test_reference, only: n_ao, n_occ, n_particle_ref => n_particle
        use otr_common_unit_tests, only: &
            generate_random_symm_matrix, mock_requests, mock_response_factor, &
            shell_names, mock_get_response_cs, mock_get_response_os, mock_evaluate_dm_os

        ! step of the finite difference and its tolerance relative to the second
        ! derivative it approximates
        real(rp), parameter :: step = 1e-3_rp, fd_tol = 1e-5_rp

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), model_h(n_ao, n_ao), fock_mo(n_mo, n_mo), &
                    model_scale, second_diff
        real(rp), allocatable :: x(:), hess_x(:), expected_hess_x(:)
        integer(ip) :: n_particle, i_case, i, error
        character(:), allocatable :: case_name

        ! assume tests pass
        test_hess_x_mo_callback = .true.

        ! generate random orthonormal MO coefficients in a non-orthonormal AO basis, so
        ! that the MO coefficients differ from their inverse and the density matrix
        ! response and the Fock matrix response have to be transformed with the right
        ! one
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)

        ! loop over every occupation case
        do i_case = 1, n_cases
            case_name = trim(case_names(i_case))
            n_particle = case_n_particle(i_case)

            ! set up the MO object with random Fock matrix blocks and the mock response
            ! function
            call setup_mo_object(mo_coeff(:, :, :n_particle), ao_overlap, &
                                 case_n_occ(:n_particle, i_case))
            call setup_settings(mo_object%settings)
            call setup_random_mo_channels(case_n_occ(:n_particle, i_case), n_mo)
            if (n_particle == 1) then
                mo_object%get_response_cs => mock_get_response_cs
            else
                mo_object%get_response_os => mock_get_response_os
            end if
            mo_object%response_stale = .false.

            ! call routine for a random trial vector and determine if values of
            ! resulting Hessian linear transformation match
            allocate(x(mo_object%n_param), hess_x(mo_object%n_param))
            call random_number(x)
            expected_hess_x = ref_hess_x_mo(x, mo_object%mo_channels, &
                                            mo_coeff(:, :, :n_particle), &
                                            mock_response_factor)
            call hess_x_mo_callback(x, hess_x, error)
            if (error /= 0) then
                write (stderr, *) "test_hess_x_mo_callback failed: Produced error "// &
                    "for the "//case_name//" case."
                test_hess_x_mo_callback = .false.
            end if
            if (norm2(hess_x - expected_hess_x) > tol) then
                write (stderr, *) "test_hess_x_mo_callback failed: Incorrect "// &
                    "Hessian linear transformation for the "//case_name//" case."
                test_hess_x_mo_callback = .false.
            end if
            deallocate(mo_object, x, hess_x)
        end do

        ! set up the open-shell MO object again, mark the response as stale, as ARH
        ! would after moving the orbitals without going through the orbital updating
        ! routine, replace the response function by one of outdated orbitals, and
        ! determine if this triggers the density matrix evaluating function to refresh
        ! the response, which is then used
        call setup_mo_object(mo_coeff, ao_overlap, n_occ)
        call setup_settings(mo_object%settings)
        call setup_random_mo_channels(n_occ, n_mo)
        mo_object%evaluate_dm_os => mock_evaluate_dm_os
        mo_object%get_response_os => mock_get_stale_response_os
        mo_object%response_stale = .true.
        allocate(x(mo_object%n_param), hess_x(mo_object%n_param))
        call random_number(x)
        expected_hess_x = &
            ref_hess_x_mo(x, mo_object%mo_channels, mo_coeff, mock_response_factor)
        mock_requests = [integer(ip) ::]
        call hess_x_mo_callback(x, hess_x, error)
        if (error /= 0) then
            write (stderr, *) "test_hess_x_mo_callback failed: Produced error "// &
                "while refreshing a stale response."
            test_hess_x_mo_callback = .false.
        end if
        if (size(mock_requests) /= 1) then
            write (stderr, *) "test_hess_x_mo_callback failed: Stale response was "// &
                "not refreshed."
            test_hess_x_mo_callback = .false.
        end if
        if (any(mock_requests /= 2)) then
            write (stderr, *) "test_hess_x_mo_callback failed: Incorrect outputs "// &
                "requested from density matrix evaluating function while "// &
                "refreshing a stale response."
            test_hess_x_mo_callback = .false.
        end if
        if (mo_object%response_stale) then
            write (stderr, *) "test_hess_x_mo_callback failed: Response still "// &
                "marked stale after being refreshed."
            test_hess_x_mo_callback = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write (stderr, *) "test_hess_x_mo_callback failed: Incorrect Hessian "// &
                "linear transformation with a refreshed response."
            test_hess_x_mo_callback = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object, x, hess_x)

        ! set up the MO object with the Fock matrix and the response of a quadratic
        ! model energy at its current orbitals, and determine if the quadratic form of
        ! the Hessian linear transformation along a random trial vector matches the
        ! second finite difference of the model energy along it, which checks the
        ! static and the response part together with their scaling independently of how
        ! the routine assembles them
        model_h = generate_random_symm_matrix(n_ao)
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            model_scale = merge(1.0_rp, 0.5_rp, n_particle == 1)
            call setup_mo_object(mo_coeff(:, :, :n_particle), ao_overlap, &
                                 n_occ(:n_particle))
            call setup_settings(mo_object%settings)
            do i = 1, n_particle
                fock_mo = ref_mo_transform( &
                    mo_coeff(:, :, i), model_h + model_scale * mo_object%dm_ao(:, :, i))
                mo_object%mo_channels(i)%fock_oo = fock_mo(:n_occ(i), :n_occ(i))
                mo_object%mo_channels(i)%fock_vv = fock_mo(n_occ(i) + 1:, n_occ(i) + 1:)
            end do
            if (n_particle == 1) then
                mo_object%get_response_cs => model_response_cs
            else
                mo_object%get_response_os => model_response_os
            end if
            mo_object%response_stale = .false.
            allocate(x(mo_object%n_param), hess_x(mo_object%n_param))
            call random_number(x)
            x = x - 0.5_rp
            call hess_x_mo_callback(x, hess_x, error)
            if (error /= 0) then
                write (stderr, *) "test_hess_x_mo_callback failed: Produced error "// &
                    "for the quadratic model energy for the "//case_name//" case."
                test_hess_x_mo_callback = .false.
            end if
            second_diff = (model_energy(step * x) - 2.0_rp * &
                           model_energy(0.0_rp * x) + model_energy(-step * x)) / step**2
            if (abs(dot_product(x, hess_x) - second_diff) > &
                fd_tol * max(1.0_rp, abs(second_diff))) then
                write (stderr, *) "test_hess_x_mo_callback failed: Hessian linear "// &
                    "transformation does not match the finite difference of the "// &
                    "quadratic model energy for the "//case_name//" case."
                test_hess_x_mo_callback = .false.
            end if
            deallocate(mo_object, x, hess_x)
        end do

    contains

        function model_energy(kappa) result(energy)
            !
            ! this function evaluates the quadratic model energy at the current
            ! orbitals of the shell rotated by kappa: for the closed-shell density
            ! matrix D of a single spin it is 2 tr(H D) + tr(D D), whose Fock matrix,
            ! half its derivative as in the orbital bases, is H + D, and for the
            ! open-shell one it is the sum over the spins of tr(H D) + tr(D D) / 4,
            ! whose Fock matrix is H + D / 2
            !
            real(rp), intent(in) :: kappa(:)
            real(rp) :: energy

            real(rp) :: rot_mo_coeff(n_ao, n_mo, n_particle), dm(n_ao, n_ao)
            integer(ip) :: k

            rot_mo_coeff = ref_rotate_mo_coeff(kappa, mo_coeff(:, :, :n_particle), &
                                               n_occ(:n_particle))
            energy = 0.0_rp
            do k = 1, n_particle
                dm = matmul(rot_mo_coeff(:, :n_occ(k), k), &
                            transpose(rot_mo_coeff(:, :n_occ(k), k)))
                energy = energy + merge(2.0_rp, 1.0_rp, n_particle == 1) * &
                         (sum(model_h * dm) + 0.5_rp * model_scale * sum(dm * dm))
            end do

        end function model_energy

        subroutine model_response_cs(dm, response_out, error_out)
            !
            ! this subroutine is the response function of the closed-shell quadratic
            ! model energy
            !
            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out), target, contiguous :: response_out(:, :)
            integer(ip), intent(out) :: error_out

            error_out = 0
            response_out = model_scale * dm

        end subroutine model_response_cs

        subroutine model_response_os(dm, response_out, error_out)
            !
            ! this subroutine is the response function of the open-shell quadratic
            ! model energy
            !
            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out), target :: response_out(:, :, :)
            integer(ip), intent(out) :: error_out

            error_out = 0
            response_out = model_scale * dm

        end subroutine model_response_os

        subroutine mock_get_stale_response_os(dm, response_out, error_out)
            !
            ! this subroutine is a mock response function for the open-shell case which
            ! belongs to outdated orbitals and returns a vanishing response
            !
            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out), target :: response_out(:, :, :)
            integer(ip), intent(out) :: error_out

            error_out = 0
            response_out = 0.0_rp * dm

        end subroutine mock_get_stale_response_os

    end function test_hess_x_mo_callback

    logical(c_bool) function test_precond_mo_callback() bind(C)
        !
        ! this function tests the subroutine which defines a level-shifted
        ! preconditioner of the MO basis, which applies the one of the orbital basis to
        ! the MO object
        !
        use otr_mo, only: precond_mo_callback, mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_occ

        real(rp), parameter :: mu = 0.3_rp

        real(rp), allocatable :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! assume tests pass
        test_precond_mo_callback = .true.

        ! set up the MO object with an eigendecomposition whose eigenvalue pairs, the
        ! open-shell differences 2 (e_v - e_o), shifted by the level shift are all 2,
        ! so that the preconditioner of the orbital basis halves the residual
        call setup_minimal_mo_object(n_occ)
        call setup_settings(mo_object%settings)
        call setup_identity_mo_eigenbasis(0.0_rp, (2.0_rp + mu) / 2.0_rp)
        allocate(residual(mo_object%n_param), precond_residual(mo_object%n_param))
        call random_number(residual)

        ! call routine and determine if the residual is halved
        call precond_mo_callback(residual, mu, precond_residual, error)
        if (error /= 0) then
            write (stderr, *) "test_precond_mo_callback failed: Produced error."
            test_precond_mo_callback = .false.
        end if
        if (norm2(precond_residual - 0.5_rp * residual) > tol) then
            write (stderr, *) "test_precond_mo_callback failed: Incorrect "// &
                "preconditioned residual."
            test_precond_mo_callback = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object)

    end function test_precond_mo_callback

    logical(c_bool) function test_precond_pd_mo_callback() bind(C)
        !
        ! this function tests the subroutine which defines the positive-definite
        ! preconditioner of the MO basis, which applies the one of the orbital basis to
        ! the MO object
        !
        use otr_mo, only: precond_pd_mo_callback, mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_occ

        real(rp), allocatable :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! assume tests pass
        test_precond_pd_mo_callback = .true.

        ! set up the MO object with an eigendecomposition whose eigenvalue pairs, the
        ! open-shell differences 2 (e_v - e_o), are all -2, so that the
        ! positive-definite preconditioner of the orbital basis, which divides by their
        ! magnitudes, halves the residual, unlike the level-shifted one
        call setup_minimal_mo_object(n_occ)
        call setup_settings(mo_object%settings)
        call setup_identity_mo_eigenbasis(1.0_rp, 0.0_rp)
        allocate(residual(mo_object%n_param), precond_residual(mo_object%n_param))
        call random_number(residual)

        ! call routine and determine if the residual is halved
        call precond_pd_mo_callback(residual, precond_residual, error)
        if (error /= 0) then
            write (stderr, *) "test_precond_pd_mo_callback failed: Produced error."
            test_precond_pd_mo_callback = .false.
        end if
        if (norm2(precond_residual - 0.5_rp * residual) > tol) then
            write (stderr, *) "test_precond_pd_mo_callback failed: Incorrect "// &
                "preconditioned residual."
            test_precond_pd_mo_callback = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object)

    end function test_precond_pd_mo_callback

    logical(c_bool) function test_get_extra_trial_vectors_mo_callback() bind(C)
        !
        ! this function tests the subroutine which returns the extra trial vectors of
        ! the MO basis for the solver's initial trial space
        !
        use otr_mo, only: get_extra_trial_vectors_mo_callback, mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_occ

        integer(ip), parameter :: n_extra = 2

        real(rp), allocatable :: trial_vectors(:, :), expected(:, :)
        integer(ip) :: error

        ! assume tests pass
        test_get_extra_trial_vectors_mo_callback = .true.

        ! set up the MO object with an eigendecomposition in the MO basis itself whose
        ! only negative eigenvalue pair belongs to the first occupied and the first
        ! virtual orbital of the first particle channel, at packed index 1, so that the
        ! first extra trial vector is the unit vector along it and the second vanishes
        call setup_minimal_mo_object(n_occ)
        call setup_settings(mo_object%settings)
        call setup_identity_mo_eigenbasis(0.0_rp, 1.0_rp)
        mo_object%mo_channels(1)%occ_eigvals(1) = 1.5_rp
        mo_object%mo_channels(1)%virt_eigvals(2) = 2.0_rp
        allocate(trial_vectors(mo_object%n_param, n_extra), &
                 expected(mo_object%n_param, n_extra))
        expected = 0.0_rp
        expected(1, 1) = 1.0_rp

        ! call routine and determine if the extra trial vectors of the MO basis are
        ! returned
        call get_extra_trial_vectors_mo_callback(trial_vectors, error)
        if (error /= 0) then
            write (stderr, *) "test_get_extra_trial_vectors_mo_callback failed: "// &
                "Produced error."
            test_get_extra_trial_vectors_mo_callback = .false.
        end if
        if (norm2(trial_vectors - expected) > tol) then
            write (stderr, *) "test_get_extra_trial_vectors_mo_callback failed: "// &
                "Incorrect extra trial vectors."
            test_get_extra_trial_vectors_mo_callback = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object)

    end function test_get_extra_trial_vectors_mo_callback

    logical(c_bool) function test_init_mo_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the MO settings
        !
        use otr_mo, only: mo_settings_type, default_settings => default_mo_settings
        use otr_mo_test_reference, only: operator(==)

        type(mo_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_init_mo_settings = .true.

        ! initialize settings
        call settings%init(error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_init_mo_settings failed: Function raised error."
            test_init_mo_settings = .false.
        end if

        ! check settings
        if (.not. (settings == default_settings)) then
            write (stderr, *) "test_init_mo_settings failed: Settings not "// &
                "initialized correctly."
            test_init_mo_settings = .false.
        end if

    end function test_init_mo_settings

    logical(c_bool) function test_mo_deconstructor() bind(C)
        !
        ! this function tests the subroutine which deallocates the MO objects
        !
        use otr_mo, only: mo_deconstructor, mo_object

        ! assume tests pass
        test_mo_deconstructor = .true.

        ! allocate MO object
        if (.not. allocated(mo_object)) allocate(mo_object)

        ! call routine and determine if the object is deallocated
        call mo_deconstructor()
        if (allocated(mo_object)) then
            write (stderr, *) "test_mo_deconstructor failed: Object not deallocated."
            test_mo_deconstructor = .false.
        end if

        ! call routine again and determine if an already deallocated object is handled
        call mo_deconstructor()
        if (allocated(mo_object)) then
            write (stderr, *) "test_mo_deconstructor failed: Already deallocated "// &
                "object not handled."
            test_mo_deconstructor = .false.
        end if

    end function test_mo_deconstructor

    logical(c_bool) function test_rotate_orbitals_mo() bind(C)
        !
        ! this function tests the subroutine which moves the current orbitals by an
        ! orbital rotation in the MO basis
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_mo, n_param_os
        use otr_common_test_reference, only: n_ao, n_occ, n_particle

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), expected(n_ao, n_mo, n_particle), &
                    kappa(n_param_os)
        integer(ip) :: error, k
        type(orbital_settings_type) :: settings

        ! assume tests pass
        test_rotate_orbitals_mo = .true.

        ! setup settings object
        call setup_settings(settings)

        ! set up the MO object with random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_mo_object(mo_coeff, ao_overlap, n_occ)

        ! call routine and determine if the MO coefficients of the caller and the
        ! density matrix are rotated and the response, which was not rebuilt, is marked
        ! stale
        call random_number(kappa)
        kappa = 0.2_rp * (kappa - 0.5_rp)
        expected = ref_rotate_mo_coeff(kappa, mo_coeff, n_occ)
        mo_object%response_stale = .false.
        call mo_object%rotate_orbitals(kappa, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_rotate_orbitals_mo failed: Produced error."
            test_rotate_orbitals_mo = .false.
        end if
        if (.not. mo_object%response_stale) then
            write (stderr, *) "test_rotate_orbitals_mo failed: Response not marked "// &
                "stale."
            test_rotate_orbitals_mo = .false.
        end if
        if (norm2(mo_coeff - expected) > tol) then
            write (stderr, *) "test_rotate_orbitals_mo failed: MO coefficients not "// &
                "rotated in place."
            test_rotate_orbitals_mo = .false.
        end if
        do k = 1, n_particle
            if (norm2(mo_object%dm_ao(:, :, k) - matmul( &
                expected(:, :n_occ(k), k), transpose(expected(:, :n_occ(k), k)))) > &
                tol) then
                write (stderr, *) "test_rotate_orbitals_mo failed: Density matrix "// &
                    "not updated."
                test_rotate_orbitals_mo = .false.
            end if
        end do

        ! deallocate MO object
        deallocate(mo_object)

    end function test_rotate_orbitals_mo

    logical(c_bool) function test_calculate_grad_h_diag_mo() bind(C)
        !
        ! this function tests the subroutine which calculates the gradient, Hessian
        ! diagonal and static part of the Hessian at the current density for the MO
        ! basis
        !
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo, n_param_cs, n_cases, case_n_particle, &
                                         case_n_occ, case_names
        use otr_common_test_reference, only: n_ao, n_occ_ref => n_occ, &
                                             n_particle_ref => n_particle
        use otr_common_unit_tests, only: generate_random_symm_matrix

        real(rp), parameter :: fd_step = 1e-4_rp

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), fock_ao(n_ao, n_ao, n_particle_ref), &
                    fock_mo(n_mo, n_mo, n_particle_ref), shell_scale, &
                    kappa(n_param_cs), rot_mo_coeff(n_ao, n_mo, 1), energy, &
                    fd_grad(n_param_cs)
        real(rp), allocatable :: expected_h_diag(:)
        integer(ip), allocatable :: n_occ(:)
        integer(ip) :: n_particle, n_occ_k, i, j, k, i_case, i_sign
        character(:), allocatable :: case_name

        ! assume tests pass
        test_calculate_grad_h_diag_mo = .true.

        ! loop over every occupation case
        do i_case = 1, n_cases
            case_name = trim(case_names(i_case))
            n_particle = case_n_particle(i_case)
            n_occ = case_n_occ(:n_particle, i_case)
            shell_scale = merge(4.0_rp, 2.0_rp, n_particle == 1)

            ! set up the MO object with random orthonormal MO coefficients and
            ! independently transform random Fock matrices to the MO basis
            ao_overlap = generate_random_ao_overlap(n_ao)
            mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
            call setup_mo_object(mo_coeff(:, :, :n_particle), ao_overlap, n_occ)
            do k = 1, n_particle
                fock_ao(:, :, k) = generate_random_symm_matrix(n_ao)
                fock_mo(:, :, k) = ref_mo_transform(mo_coeff(:, :, k), fock_ao(:, :, k))
            end do
            mo_object%hess_eigen_stale = .false.

            ! expected Hessian diagonal from the diagonal Fock matrix elements
            expected_h_diag = [(((shell_scale * (fock_mo(j, j, k) - fock_mo(i, i, k)), &
                                  i = 1, n_occ(k)), j = n_occ(k) + 1, n_mo), k = 1, &
                                n_particle)]

            ! call routine and determine if the gradient, Hessian diagonal and Fock
            ! matrix blocks are obtained in the MO basis and the eigendecomposition is
            ! marked stale
            call mo_object%calculate_grad_h_diag(fock_ao(:, :, :n_particle))
            if (norm2(mo_object%grad - shell_scale * &
                      ref_pack_ov(fock_mo(:, :, :n_particle), n_occ)) > tol) then
                write (stderr, *) "test_calculate_grad_h_diag_mo failed: Incorrect "// &
                    "gradient for the "//case_name//" case."
                test_calculate_grad_h_diag_mo = .false.
            end if
            if (norm2(mo_object%h_diag - expected_h_diag) > tol) then
                write (stderr, *) "test_calculate_grad_h_diag_mo failed: Incorrect "// &
                    "Hessian diagonal for the "//case_name//" case."
                test_calculate_grad_h_diag_mo = .false.
            end if
            do k = 1, n_particle
                n_occ_k = n_occ(k)
                if (norm2(mo_object%mo_channels(k)%fock_oo - &
                          fock_mo(:n_occ_k, :n_occ_k, k)) > tol) then
                    write (stderr, *) "test_calculate_grad_h_diag_mo failed: "// &
                        "Incorrect occupied-occupied Fock matrix block for the "// &
                        case_name//" case."
                    test_calculate_grad_h_diag_mo = .false.
                end if
                if (norm2(mo_object%mo_channels(k)%fock_vv - &
                          fock_mo(n_occ_k + 1:, n_occ_k + 1:, k)) > tol) then
                    write (stderr, *) "test_calculate_grad_h_diag_mo failed: "// &
                        "Incorrect virtual-virtual Fock matrix block for the "// &
                        case_name//" case."
                    test_calculate_grad_h_diag_mo = .false.
                end if
            end do
            if (.not. mo_object%hess_eigen_stale) then
                write (stderr, *) "test_calculate_grad_h_diag_mo failed: "// &
                    "Eigendecomposition not marked stale for the "//case_name//" case."
                test_calculate_grad_h_diag_mo = .false.
            end if

            ! deallocate MO object
            deallocate(mo_object)
        end do

        ! call routine for a random symmetric closed-shell Fock matrix H, which, unlike
        ! that of the mock, is consistent with an energy, the linear energy 2 tr(H D),
        ! and determine if the gradient is the finite-difference derivative of the
        ! energy along the reference orbital rotation, which pins the rotation
        ! convention and parameter order against the gradient
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
        call setup_mo_object(mo_coeff(:, :, :1), ao_overlap, n_occ_ref(:1))
        fock_ao(:, :, 1) = generate_random_symm_matrix(n_ao)
        call mo_object%calculate_grad_h_diag(fock_ao(:, :, :1))
        do k = 1, n_param_cs
            fd_grad(k) = 0.0_rp
            do i_sign = -1, 1, 2
                kappa = 0.0_rp
                kappa(k) = i_sign * fd_step
                rot_mo_coeff = ref_rotate_mo_coeff(kappa, mo_coeff(:, :, :1), &
                                                   n_occ_ref(:1))
                energy = 2.0_rp * sum(fock_ao(:, :, 1) * matmul( &
                    rot_mo_coeff(:, :n_occ_ref(1), 1), &
                    transpose(rot_mo_coeff(:, :n_occ_ref(1), 1))))
                fd_grad(k) = fd_grad(k) + i_sign * energy / (2.0_rp * fd_step)
            end do
        end do
        if (norm2(mo_object%grad - fd_grad) > 1e-6_rp * norm2(fd_grad)) then
            write (stderr, *) "test_calculate_grad_h_diag_mo failed: Gradient not "// &
                "the derivative of the energy."
            test_calculate_grad_h_diag_mo = .false.
        end if
        deallocate(mo_object)

    end function test_calculate_grad_h_diag_mo

    logical(c_bool) function test_refresh_hess_eigen_mo() bind(C)
        !
        ! this function tests the subroutine which refreshes the eigendecomposition of
        ! the static part of the Hessian in the MO basis, which diagonalizes the
        ! occupied-occupied and virtual-virtual blocks of the Fock matrix, if it is
        ! stale
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_cases, case_n_particle, case_n_occ, &
                                         case_names
        use otr_common_unit_tests, only: identity_matrix, generate_random_symm_matrix

        real(rp), allocatable :: eigvals_before(:)
        integer(ip) :: error, k, n_occ, n_virt, i_case, n_particle
        character(:), allocatable :: case_name
        type(orbital_settings_type) :: settings

        ! assume tests pass
        test_refresh_hess_eigen_mo = .true.

        ! setup settings object
        call setup_settings(settings)

        ! loop over every occupation case: both Fock matrix blocks of every channel,
        ! which are replaced so that they no longer match the cached
        ! eigendecomposition, are diagonalized
        do i_case = 1, n_cases
            n_particle = case_n_particle(i_case)
            case_name = trim(case_names(i_case))
            call setup_minimal_mo_object(case_n_occ(:n_particle, i_case))
            do k = 1, n_particle
                associate(channel => mo_object%mo_channels(k))
                    channel%fock_oo = generate_random_symm_matrix(channel%n_occ)
                    channel%fock_vv = generate_random_symm_matrix(channel%n_virt)
                end associate
            end do
            mo_object%hess_eigen_stale = .true.
            call mo_object%refresh_hess_eigen(settings, error)
            if (error /= 0) then
                write (stderr, *) "test_refresh_hess_eigen_mo failed: Produced "// &
                    "error for the "//case_name//" case."
                test_refresh_hess_eigen_mo = .false.
            end if
            if (mo_object%hess_eigen_stale) then
                write (stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                    "Eigendecomposition still stale for the "//case_name//" case."
                test_refresh_hess_eigen_mo = .false.
            end if
            do k = 1, n_particle
                associate(channel => mo_object%mo_channels(k))
                    n_occ = channel%n_occ
                    n_virt = channel%n_virt
                    if (norm2(matmul(transpose(channel%occ_eigvecs), &
                                     matmul(channel%fock_oo, channel%occ_eigvecs)) - &
                              diagonal_matrix(channel%occ_eigvals)) > tol) then
                        write (stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Occupied-occupied block not diagonalized for the "// &
                            case_name//" case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                    if (norm2(matmul(transpose(channel%occ_eigvecs), &
                                     channel%occ_eigvecs) - identity_matrix(n_occ)) > &
                        tol) then
                        write (stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Occupied eigenvectors not orthonormal for the "// &
                            case_name//" case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                    if (norm2(matmul(transpose(channel%virt_eigvecs), &
                                     matmul(channel%fock_vv, channel%virt_eigvecs)) - &
                              diagonal_matrix(channel%virt_eigvals)) > tol) then
                        write (stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Virtual-virtual block not diagonalized for the "// &
                            case_name//" case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                    if (norm2(matmul(transpose(channel%virt_eigvecs), &
                                     channel%virt_eigvecs) - &
                              identity_matrix(n_virt)) > tol) then
                        write (stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Virtual eigenvectors not orthonormal for the "// &
                            case_name//" case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                end associate
            end do

            ! double the occupied block of the first channel, which is occupied in
            ! every case, and determine if an up-to-date eigendecomposition is not
            ! recomputed
            eigvals_before = mo_object%mo_channels(1)%occ_eigvals
            mo_object%mo_channels(1)%fock_oo = 2.0_rp * mo_object%mo_channels(1)%fock_oo
            call mo_object%refresh_hess_eigen(settings, error)
            if (error /= 0) then
                write (stderr, *) "test_refresh_hess_eigen_mo failed: Produced "// &
                    "error without stale eigendecomposition for the "//case_name// &
                    " case."
                test_refresh_hess_eigen_mo = .false.
            end if
            if (norm2(mo_object%mo_channels(1)%occ_eigvals - eigvals_before) > tol) then
                write (stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                    "Eigendecomposition recomputed although not stale for the "// &
                    case_name//" case."
                test_refresh_hess_eigen_mo = .false.
            end if
            deallocate(mo_object)
        end do

    end function test_refresh_hess_eigen_mo

    logical(c_bool) function test_rotate_to_hess_eigenbasis_mo() bind(C)
        !
        ! this function tests the function which rotates a parameter vector into the
        ! eigenbasis of the static part of the Hessian for the MO basis
        !
        test_rotate_to_hess_eigenbasis_mo = check_rotate_hess_eigenbasis_mo(.true.)

    end function test_rotate_to_hess_eigenbasis_mo

    logical(c_bool) function test_rotate_from_hess_eigenbasis_mo() bind(C)
        !
        ! this function tests the function which rotates a parameter vector out of the
        ! eigenbasis of the static part of the Hessian for the MO basis
        !
        test_rotate_from_hess_eigenbasis_mo = check_rotate_hess_eigenbasis_mo(.false.)

    end function test_rotate_from_hess_eigenbasis_mo

    logical(c_bool) function test_get_hess_eigval_pairs_mo() bind(C)
        !
        ! this function tests the function which returns the eigenvalues of the static
        ! part of the Hessian for the MO basis
        !
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_cases, case_n_particle, case_n_occ, &
                                         case_names

        integer(ip) :: i_case

        ! assume tests pass
        test_get_hess_eigval_pairs_mo = .true.

        ! loop over every occupation case: differences of virtual and occupied
        ! eigenvalues
        do i_case = 1, n_cases
            call setup_minimal_mo_object(case_n_occ(:case_n_particle(i_case), i_case))
            if (norm2(mo_object%get_hess_eigval_pairs() - &
                      ref_hess_eigval_pairs_mo(mo_object%mo_channels)) > tol) then
                write (stderr, *) "test_get_hess_eigval_pairs_mo failed: Incorrect "// &
                    "eigenvalue differences for the "//trim(case_names(i_case))// &
                    " case."
                test_get_hess_eigval_pairs_mo = .false.
            end if
            deallocate(mo_object)
        end do

    end function test_get_hess_eigval_pairs_mo

    logical(c_bool) function test_get_extra_trial_vectors_mo() bind(C)
        !
        ! this function tests the subroutine which returns extra trial vectors along
        ! the most negative eigenvalue pairs of the static part of the Hessian for the
        ! MO basis
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_occ
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_extra = 3

        type(orbital_settings_type) :: settings
        real(rp), allocatable :: trial_vectors(:, :), expected(:, :), pairs(:), &
                                 unit_vector(:)
        integer(ip) :: error, i, j, k, ivec, min_idx

        ! assume tests pass
        test_get_extra_trial_vectors_mo = .true.

        ! setup settings object
        call setup_settings(settings)

        ! set orbital energies with exactly two negative pairs in the first channel and
        ! none in the second, so that the last requested vector has to vanish
        call setup_minimal_mo_object(n_occ)
        mo_object%mo_channels(1)%occ_eigvals = [-1.0_rp, 2.0_rp]
        mo_object%mo_channels(1)%virt_eigvals = [0.5_rp, 1.5_rp]
        mo_object%mo_channels(2)%occ_eigvals = [-3.0_rp]
        mo_object%mo_channels(2)%virt_eigvals = [0.5_rp, 1.0_rp, 2.0_rp]

        ! construct the expected trial vectors, the rotations along the eigenvalue
        ! pairs in increasing order of the negative ones, rotated out of the
        ! eigenbasis, with vanishing slots once no pair is negative any more
        pairs = ref_hess_eigval_pairs_mo(mo_object%mo_channels)
        allocate(expected(size(pairs), n_extra), unit_vector(size(pairs)))
        expected = 0.0_rp
        do ivec = 1, n_extra
            min_idx = minloc(pairs, dim=1)
            if (pairs(min_idx) >= 0.0_rp) exit
            unit_vector = 0.0_rp
            unit_vector(min_idx) = 1.0_rp
            expected(:, ivec) = ref_rotate_eigenbasis_mo( &
                unit_vector, mo_object%mo_channels, n_mo, .false.)
            pairs(min_idx) = huge(1.0_rp)
        end do
        if (norm2(expected(:, n_extra)) > tol .or. norm2(expected(:, 2)) < tol) then
            write (stderr, *) "test_get_extra_trial_vectors_mo failed: Test "// &
                "fixture does not have exactly two negative eigenvalue pairs."
            test_get_extra_trial_vectors_mo = .false.
        end if

        ! call routine and determine if the trial vectors are correct
        allocate(trial_vectors(mo_object%n_param, n_extra))
        call mo_object%get_extra_trial_vectors(trial_vectors, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_extra_trial_vectors_mo failed: Produced error."
            test_get_extra_trial_vectors_mo = .false.
        end if
        if (norm2(trial_vectors - expected) > tol) then
            write (stderr, *) "test_get_extra_trial_vectors_mo failed: Incorrect "// &
                "trial vectors."
            test_get_extra_trial_vectors_mo = .false.
        end if

        ! mark the eigendecomposition stale for diagonal Fock matrix blocks with the
        ! same orbital energies, whose eigenvectors are unit vectors since the
        ! diagonals are ascending, and determine if it is refreshed before the trial
        ! vectors, then unit vectors up to their sign, are constructed
        do k = 1, 2
            associate(channel => mo_object%mo_channels(k))
                channel%fock_oo = diagonal_matrix(channel%occ_eigvals)
                channel%fock_vv = diagonal_matrix(channel%virt_eigvals)
            end associate
        end do
        mo_object%hess_eigen_stale = .true.
        pairs = [(((mo_object%mo_channels(k)%virt_eigvals(j) - &
                    mo_object%mo_channels(k)%occ_eigvals(i), i = 1, n_occ(k)), j = 1, &
                   n_mo - n_occ(k)), k = 1, 2)]
        expected = 0.0_rp
        do ivec = 1, 2
            min_idx = minloc(pairs, dim=1)
            expected(min_idx, ivec) = 1.0_rp
            pairs(min_idx) = huge(1.0_rp)
        end do
        call mo_object%get_extra_trial_vectors(trial_vectors, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_extra_trial_vectors_mo failed: Produced "// &
                "error for a stale eigendecomposition."
            test_get_extra_trial_vectors_mo = .false.
        end if
        if (norm2(abs(trial_vectors) - expected) > tol) then
            write (stderr, *) "test_get_extra_trial_vectors_mo failed: Incorrect "// &
                "trial vectors for a stale eigendecomposition."
            test_get_extra_trial_vectors_mo = .false.
        end if
        deallocate(mo_object)

    end function test_get_extra_trial_vectors_mo

    logical(c_bool) function test_finalize_mo() bind(C)
        !
        ! this function tests the subroutine which frees the density matrix the MO
        ! object allocated
        !
        use otr_mo, only: mo_type, finalize_mo

        type(mo_type) :: mo

        ! assume tests pass
        test_finalize_mo = .true.

        ! call routine for an MO object with an allocated density matrix and determine
        ! if it is freed
        allocate(mo%dm_ao(2, 2, 1))
        call finalize_mo(mo)
        if (associated(mo%dm_ao)) then
            write (stderr, *) "test_finalize_mo failed: Density matrix not freed."
            test_finalize_mo = .false.
        end if

        ! call routine for an MO object without a density matrix and determine if it is
        ! left without one
        call finalize_mo(mo)
        if (associated(mo%dm_ao)) then
            write (stderr, *) "test_finalize_mo failed: Density matrix associated "// &
                "for an object without one."
            test_finalize_mo = .false.
        end if

    end function test_finalize_mo

    logical(c_bool) function test_rotate_mo_coeff() bind(C)
        !
        ! this function tests the subroutine which rotates MO coefficients by an
        ! orbital rotation, reorthonormalizes them and returns the resulting density
        ! matrix
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: rotate_mo_coeff, mo_channel_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_names
        use otr_common_test_reference, only: n_ao, n_occ_ref => n_occ, &
                                             n_particle_ref => n_particle

        real(rp) :: mo_coeff(n_ao, n_mo, n_particle_ref), ao_overlap(n_ao, n_ao), &
                    rot_mo_coeff(n_ao, n_mo, n_particle_ref), &
                    rot_dm_ao(n_ao, n_ao, n_particle_ref), &
                    expected(n_ao, n_mo, n_particle_ref)
        real(rp), allocatable :: kappa(:)
        integer(ip), allocatable :: n_occ(:)
        integer(ip) :: n_particle, error, k, i_case
        character(:), allocatable :: case_name
        type(mo_channel_type), allocatable :: mo_channels(:)
        type(orbital_settings_type) :: settings

        ! assume tests pass
        test_rotate_mo_coeff = .true.

        ! setup settings object
        call setup_settings(settings)

        ! loop over every occupation case
        do i_case = 1, n_cases
            case_name = trim(case_names(i_case))
            n_particle = case_n_particle(i_case)
            n_occ = case_n_occ(:n_particle, i_case)

            ! generate random orthonormal MO coefficients and the occupations of their
            ! particle channels
            ao_overlap = generate_random_ao_overlap(n_ao)
            mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
            allocate(mo_channels(n_particle))
            mo_channels%n_occ = n_occ
            mo_channels%n_virt = n_mo - n_occ

            ! call routine for a random orbital rotation and determine if the rotated
            ! orbitals and density matrix are correct
            allocate(kappa(sum(n_occ * (n_mo - n_occ))))
            call random_number(kappa)
            kappa = 0.2_rp * (kappa - 0.5_rp)
            expected(:, :, :n_particle) = &
                ref_rotate_mo_coeff(kappa, mo_coeff(:, :, :n_particle), n_occ)
            call rotate_mo_coeff(kappa, mo_coeff(:, :, :n_particle), ao_overlap, &
                                 mo_channels, rot_mo_coeff(:, :, :n_particle), &
                                 rot_dm_ao(:, :, :n_particle), settings, error)
            if (error /= 0) then
                write (stderr, *) "test_rotate_mo_coeff failed: Produced error for "// &
                    "the "//case_name//" case."
                test_rotate_mo_coeff = .false.
            end if
            if (norm2(rot_mo_coeff(:, :, :n_particle) - expected(:, :, :n_particle)) > &
                tol) then
                write (stderr, *) "test_rotate_mo_coeff failed: Incorrect rotated "// &
                    "MO coefficients for the "//case_name//" case."
                test_rotate_mo_coeff = .false.
            end if
            do k = 1, n_particle
                if (norm2(rot_dm_ao(:, :, k) - &
                          matmul(expected(:, :n_occ(k), k), &
                                 transpose(expected(:, :n_occ(k), k)))) > tol) then
                    write (stderr, *) "test_rotate_mo_coeff failed: Incorrect "// &
                        "rotated density matrix for the "//case_name//" case."
                    test_rotate_mo_coeff = .false.
                end if
            end do
            deallocate(kappa, mo_channels)
        end do

        ! call routine without rotation for MO coefficients which are off
        ! orthonormality and determine if they are symmetrically reorthonormalized
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
        call random_number(rot_mo_coeff)
        mo_coeff = mo_coeff + 1e-2_rp * rot_mo_coeff
        allocate(mo_channels(1))
        mo_channels%n_occ = n_occ_ref(1)
        mo_channels%n_virt = n_mo - n_occ_ref(1)
        allocate(kappa(n_occ_ref(1) * (n_mo - n_occ_ref(1))))
        kappa = 0.0_rp
        call rotate_mo_coeff(kappa, mo_coeff(:, :, :1), ao_overlap, mo_channels, &
                             rot_mo_coeff(:, :, :1), rot_dm_ao(:, :, :1), settings, &
                             error)
        if (error /= 0) then
            write (stderr, *) "test_rotate_mo_coeff failed: Produced error for MO "// &
                "coefficients off orthonormality."
            test_rotate_mo_coeff = .false.
        end if
        if (norm2(rot_mo_coeff(:, :, 1) - &
                  ref_lowdin_orthonormalize(mo_coeff(:, :, 1))) > tol) then
            write (stderr, *) "test_rotate_mo_coeff failed: MO coefficients not "// &
                "symmetrically reorthonormalized."
            test_rotate_mo_coeff = .false.
        end if

    contains

        function ref_lowdin_orthonormalize(coeff) result(orth_coeff)
            !
            ! this function independently reproduces the symmetric reorthonormalization
            ! C (C^T S C)^(-1/2) of MO coefficients with respect to the AO overlap
            ! matrix of the test
            !
            real(rp), intent(in) :: coeff(:, :)
            real(rp) :: orth_coeff(n_ao, n_mo)

            integer(ip) :: lwork, info
            real(rp) :: eigvecs(n_mo, n_mo), eigvals(n_mo)
            real(rp), allocatable :: work(:)
            external :: dsyev

            ! eigendecomposition of the metric of the MO coefficients
            eigvecs = matmul(transpose(coeff), matmul(ao_overlap, coeff))
            allocate(work(1))
            call dsyev("V", "U", n_mo, eigvecs, n_mo, eigvals, work, -1_ip, info)
            lwork = int(work(1))
            deallocate(work)
            allocate(work(lwork))
            call dsyev("V", "U", n_mo, eigvecs, n_mo, eigvals, work, lwork, info)
            deallocate(work)

            ! apply the inverse square root of the metric
            orth_coeff = matmul(coeff, matmul( &
                eigvecs * spread(1.0_rp / sqrt(eigvals), 1, n_mo), transpose(eigvecs)))

        end function ref_lowdin_orthonormalize

    end function test_rotate_mo_coeff

    logical(c_bool) function test_hess_x_static_mo() bind(C)
        !
        ! this function tests the function which applies the static part of the Hessian
        ! in the MO basis to a trial vector
        !
        use otr_mo, only: hess_x_static_mo, mo_object
        use otr_mo_test_reference, only: n_cases, case_n_particle, case_n_occ, &
                                         case_names

        real(rp), allocatable :: x(:)
        integer(ip) :: i_case

        ! assume tests pass
        test_hess_x_static_mo = .true.

        ! loop over every occupation case
        do i_case = 1, n_cases
            call setup_minimal_mo_object(case_n_occ(:case_n_particle(i_case), i_case))
            allocate(x(mo_object%n_param))
            call random_number(x)
            if (norm2(hess_x_static_mo(x, mo_object%mo_channels) - &
                      ref_hess_x_static_mo(x, mo_object%mo_channels)) > tol) then
                write (stderr, *) "test_hess_x_static_mo failed: Incorrect static "// &
                    "part for the "//trim(case_names(i_case))//" case."
                test_hess_x_static_mo = .false.
            end if
            deallocate(x, mo_object)
        end do

    end function test_hess_x_static_mo

    logical(c_bool) function test_mo_transform() bind(C)
        !
        ! this function tests the function which transforms a matrix from the AO basis
        ! with a set of orbital coefficients
        !
        use otr_mo, only: mo_transform
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao

        real(rp) :: coeff(n_ao, n_mo), matrix(n_ao, n_ao)
        real(rp), allocatable :: transformed(:, :)

        ! assume tests pass
        test_mo_transform = .true.

        ! random rectangular coefficients and a random, generically non-symmetric
        ! matrix, so that a transposed matrix cannot go unnoticed
        call random_number(coeff)
        call random_number(matrix)

        ! call function and determine if the transformed matrix is correct
        transformed = mo_transform(coeff, matrix)
        if (size(transformed, 1) /= n_mo .or. size(transformed, 2) /= n_mo) then
            write (stderr, *) "test_mo_transform failed: Incorrect dimensions."
            test_mo_transform = .false.
        else if (norm2(transformed - ref_mo_transform(coeff, matrix)) > tol) then
            write (stderr, *) "test_mo_transform failed: Incorrect transformed matrix."
            test_mo_transform = .false.
        end if

    end function test_mo_transform

end module otr_mo_unit_tests
