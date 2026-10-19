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

    function masked_mo_channels(n_occ, n_mo, irreps) result(mo_channels)
        !
        ! this function returns the particle channels of the given occupations with the
        ! irreps of their MOs, all of the same irrep if none are given, and the mask of
        ! their occupied-virtual pairs of the same irrep, the way the MO factory sets
        ! them up
        !
        use otr_mo, only: mo_channel_type

        integer(ip), intent(in) :: n_occ(:), n_mo
        integer(ip), intent(in), optional :: irreps(:, :)
        type(mo_channel_type) :: mo_channels(size(n_occ))

        integer(ip) :: j, k

        mo_channels%n_occ = n_occ
        mo_channels%n_virt = n_mo - n_occ
        do k = 1, size(n_occ, kind=ip)
            associate (channel => mo_channels(k))
                if (present(irreps)) then
                    channel%irreps = irreps(:, k)
                else
                    allocate(channel%irreps(n_mo), source=0_ip)
                end if
                allocate(channel%param_mask(n_occ(k), n_mo - n_occ(k)))
                do j = 1, n_mo - n_occ(k)
                    channel%param_mask(:, j) = &
                        channel%irreps(:n_occ(k)) == channel%irreps(n_occ(k) + j)
                end do
            end associate
        end do

    end function masked_mo_channels

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

    subroutine setup_mo_object(mo, mo_coeff, ao_overlap, n_occ, irreps)
        !
        ! this subroutine sets up the given MO object the way the MO factory would,
        ! pointing to the test-local MO coefficients
        !
        use otr_common, only: error_mo
        use otr_mo, only: mo_type

        type(mo_type), allocatable, intent(out), target :: mo
        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:)
        integer(ip), intent(in), optional :: irreps(:, :)

        integer(ip) :: n_ao, n_mo, n_particle, i

        ! dimensions
        n_ao = size(mo_coeff, 1, kind=ip)
        n_mo = size(mo_coeff, 2, kind=ip)
        n_particle = size(n_occ, kind=ip)

        ! set up the MO object
        allocate(mo)
        mo%error_origin = error_mo
        mo%n_ao = n_ao
        mo%n_mo = n_mo
        mo%n_particle = n_particle
        mo%mo_coeff => mo_coeff
        mo%ao_overlap = ao_overlap
        mo%mo_channels = masked_mo_channels(n_occ, n_mo, irreps)
        mo%n_param = &
            sum([(count(mo%mo_channels(i)%param_mask, kind=ip), i=1, n_particle)])
        allocate(mo%own_dm_ao(n_ao, n_ao, n_particle), mo%grad(mo%n_param), &
                 mo%h_diag(mo%n_param))
        mo%dm_ao => mo%own_dm_ao
        do i = 1, n_particle
            mo%dm_ao(:, :, i) = matmul(mo_coeff(:, :n_occ(i), i), &
                                       transpose(mo_coeff(:, :n_occ(i), i)))
        end do

    end subroutine setup_mo_object

    subroutine setup_random_mo_channels(mo, n_occ, n_mo)
        !
        ! this subroutine fills the occupied-occupied and virtual-virtual Fock matrix
        ! blocks of the MO channels of the MO object with random symmetric matrices,
        ! which couple only orbitals of the same irrep, together with their
        ! eigendecompositions, whose eigenvectors sit at the orbitals of their irrep
        !
        use otr_mo, only: mo_type
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix

        type(mo_type), intent(inout) :: mo
        integer(ip), intent(in) :: n_occ(:), n_mo

        integer(ip) :: k

        do k = 1, size(n_occ, kind=ip)
            associate (channel => mo%mo_channels(k))
                call generate_random_symm_eigen(channel%irreps(:n_occ(k)), &
                                                channel%occ_eigvals, &
                                                channel%occ_eigvecs, channel%fock_oo)
                call generate_random_symm_eigen(channel%irreps(n_occ(k) + 1:n_mo), &
                                                channel%virt_eigvals, &
                                                channel%virt_eigvecs, channel%fock_vv)
            end associate
        end do
        mo%hess_eigen_stale = .false.

    contains

        subroutine generate_random_symm_eigen(block_irreps, eigvals, eigvecs, matrix)
            !
            ! this subroutine returns a random symmetric matrix between orbitals of the
            ! given irreps, which couples only orbitals of the same irrep, together
            ! with its eigendecomposition
            !
            integer(ip), intent(in) :: block_irreps(:)
            real(rp), allocatable, intent(out) :: eigvals(:), eigvecs(:, :), &
                                                  matrix(:, :)

            integer(ip) :: n, i, j
            integer(ip), allocatable :: idx(:)
            logical, allocatable :: done(:)

            ! dimension of the matrix
            n = size(block_irreps, kind=ip)

            ! an empty block has an empty eigendecomposition
            if (n == 0) then
                allocate(eigvals(0), eigvecs(0, 0), matrix(0, 0))
                return
            end if

            allocate(eigvals(n))
            call random_number(eigvals)
            eigvals = 2.0_rp * eigvals - 1.0_rp

            ! non-symmetric eigenvectors within every irrep, unlike those LAPACK
            ! returns for small matrices, so that a transposed rotation cannot go
            ! unnoticed
            allocate(eigvecs(n, n), source=0.0_rp)
            allocate(done(n), source=.false.)
            do i = 1, n
                if (done(i)) cycle
                idx = pack([(j, j=1, n)], block_irreps == block_irreps(i))
                eigvecs(idx, idx) = matmul( &
                    generate_random_orthogonal_matrix(size(idx, kind=ip)), &
                    generate_random_orthogonal_matrix(size(idx, kind=ip)))
                done(idx) = .true.
            end do
            matrix = matmul(eigvecs * spread(eigvals, 1, n), transpose(eigvecs))

        end subroutine generate_random_symm_eigen

    end subroutine setup_random_mo_channels

    subroutine setup_minimal_mo_object(mo, n_occ, irreps)
        !
        ! this subroutine sets up the given MO object with only the dimensions and the
        ! particle channels of the given occupations and, if given, irreps of their
        ! MOs, whose Fock matrix blocks and eigendecompositions are random, which is
        ! all that the routines acting on the static part of the Hessian read
        !
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo

        type(mo_type), allocatable, intent(out) :: mo
        integer(ip), intent(in) :: n_occ(:)
        integer(ip), intent(in), optional :: irreps(:, :)

        integer(ip) :: i

        allocate(mo)
        mo%n_mo = n_mo
        mo%n_particle = size(n_occ, kind=ip)
        mo%mo_channels = masked_mo_channels(n_occ, n_mo, irreps)
        mo%n_param = &
            sum([(count(mo%mo_channels(i)%param_mask, kind=ip), i=1, mo%n_particle)])
        call setup_random_mo_channels(mo, n_occ, n_mo)

    end subroutine setup_minimal_mo_object

    function ref_mo_transform(coeff, matrix) result(transformed)
        !
        ! this function independently reproduces the transformation C^T M C of a matrix
        ! given in the AO basis
        !
        real(rp), intent(in) :: coeff(:, :), matrix(:, :)
        real(rp) :: transformed(size(coeff, 2), size(coeff, 2))

        transformed = matmul(transpose(coeff), matmul(matrix, coeff))

    end function ref_mo_transform

    function ref_count_mo_params(n_occ, irreps) result(n_param)
        !
        ! this function independently counts the parameters in the MO basis, the
        ! occupied-virtual pairs between orbitals of the same irrep of every particle
        ! channel
        !
        integer(ip), intent(in) :: n_occ(:), irreps(:, :)
        integer(ip) :: n_param

        integer(ip) :: i, j, k

        n_param = 0
        do k = 1, size(n_occ, kind=ip)
            do j = n_occ(k) + 1, size(irreps, 1, kind=ip)
                do i = 1, n_occ(k)
                    if (irreps(i, k) == irreps(j, k)) n_param = n_param + 1
                end do
            end do
        end do

    end function ref_count_mo_params

    function ref_pack_ov(matrices, channels) result(packed)
        !
        ! this function independently packs the occupied-virtual blocks of a set of
        ! matrices in the MO basis, one per particle channel, keeping only the entries
        ! between orbitals of the same irrep, with the occupied index running fastest
        !
        use otr_mo, only: mo_channel_type

        real(rp), intent(in) :: matrices(:, :, :)
        type(mo_channel_type), intent(in) :: channels(:)
        real(rp), allocatable :: packed(:)

        integer(ip) :: n_mo, i, j, k

        n_mo = size(matrices, 1, kind=ip)
        allocate(packed(0))
        do k = 1, size(channels, kind=ip)
            do j = channels(k)%n_occ + 1, n_mo
                do i = 1, channels(k)%n_occ
                    if (channels(k)%irreps(i) == channels(k)%irreps(j)) &
                        packed = [packed, matrices(i, j, k)]
                end do
            end do
        end do

    end function ref_pack_ov

    function ref_unpack_ov(packed, channels) result(matrices)
        !
        ! this function independently unpacks a parameter vector in the MO basis into
        ! the entries between orbitals of the same irrep of the occupied-virtual blocks
        ! of a set of otherwise vanishing matrices, one per particle channel
        !
        use otr_mo, only: mo_channel_type

        real(rp), intent(in) :: packed(:)
        type(mo_channel_type), intent(in) :: channels(:)
        real(rp) :: &
            matrices(size(channels(1)%irreps), size(channels(1)%irreps), size(channels))

        integer(ip) :: idx, i, j, k

        matrices = 0.0_rp
        idx = 0
        do k = 1, size(channels, kind=ip)
            do j = channels(k)%n_occ + 1, size(channels(k)%irreps, kind=ip)
                do i = 1, channels(k)%n_occ
                    if (channels(k)%irreps(i) /= channels(k)%irreps(j)) cycle
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

        real(rp) :: x_full(size(channels(1)%irreps), size(channels(1)%irreps), &
                           size(channels)), &
                    static_part(size(channels(1)%irreps), size(channels(1)%irreps), &
                                size(channels)), shell_scale
        integer(ip) :: n_occ, k

        shell_scale = merge(4.0_rp, 2.0_rp, size(channels) == 1)
        x_full = ref_unpack_ov(x, channels)
        static_part = 0.0_rp
        do k = 1, size(channels, kind=ip)
            n_occ = channels(k)%n_occ
            static_part(:n_occ, n_occ + 1:, k) = &
                matmul(x_full(:n_occ, n_occ + 1:, k), channels(k)%fock_vv) - &
                matmul(channels(k)%fock_oo, x_full(:n_occ, n_occ + 1:, k))
        end do
        hess_x = shell_scale * ref_pack_ov(static_part, channels)

    end function ref_hess_x_static_mo

    function ref_hess_x_mo(x, channels, mo_coeff, response_factor, response_coupling) &
        result(hess_x)
        !
        ! this function independently reproduces the Hessian linear transformation in
        ! the MO basis for a response function returning a multiple of the density
        ! matrix displacement of every particle channel plus a multiple of those of the
        ! other channels, by displacing the density matrix of every particle channel in
        ! the MO basis by the symmetrized occupied-virtual block of the trial vector,
        ! transforming the displacements to the AO basis as C dD C^T and the resulting
        ! response back to the MO basis, whose occupied-virtual block, scaled like the
        ! static part, is added to the static part
        !
        use otr_mo, only: mo_channel_type

        real(rp), intent(in) :: x(:), mo_coeff(:, :, :), response_factor, &
                                response_coupling
        type(mo_channel_type), intent(in) :: channels(:)
        real(rp) :: hess_x(size(x))

        real(rp) :: x_full(size(mo_coeff, 2), size(mo_coeff, 2), size(channels)), &
                    dm_response(size(mo_coeff, 1), size(mo_coeff, 1), size(channels)), &
                    response_mo(size(mo_coeff, 2), size(mo_coeff, 2), size(channels)), &
                    shell_scale
        integer(ip) :: k

        shell_scale = merge(4.0_rp, 2.0_rp, size(channels) == 1)
        x_full = ref_unpack_ov(x, channels)
        do k = 1, size(channels, kind=ip)
            dm_response(:, :, k) = matmul(mo_coeff(:, :, k), matmul( &
                x_full(:, :, k) + transpose(x_full(:, :, k)), &
                transpose(mo_coeff(:, :, k))))
        end do
        do k = 1, size(channels, kind=ip)
            response_mo(:, :, k) = ref_mo_transform( &
                mo_coeff(:, :, k), response_factor * dm_response(:, :, k) + &
                response_coupling * (sum(dm_response, dim=3) - dm_response(:, :, k)))
        end do
        hess_x = ref_hess_x_static_mo(x, channels) + &
                 shell_scale * ref_pack_ov(response_mo, channels)

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

        x_full = ref_unpack_ov(x, channels)
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
        rotated = ref_pack_ov(rotated_full, channels)

    end function ref_rotate_eigenbasis_mo

    function ref_hess_eigval_pairs_mo(channels) result(eigval_pairs)
        !
        ! this function independently constructs the eigenvalues of the static part of
        ! the Hessian in the MO basis, the scaled differences of the cached virtual and
        ! occupied Fock matrix block eigenvalues of every particle channel, for the
        ! pairs of eigenvectors at orbitals of the same irrep
        !
        use otr_mo, only: mo_channel_type

        type(mo_channel_type), intent(in) :: channels(:)
        real(rp), allocatable :: eigval_pairs(:)

        real(rp) :: shell_scale
        integer(ip) :: n_occ, i, j, k

        shell_scale = merge(4.0_rp, 2.0_rp, size(channels) == 1)
        allocate(eigval_pairs(0))
        do k = 1, size(channels, kind=ip)
            n_occ = channels(k)%n_occ
            do j = 1, channels(k)%n_virt
                do i = 1, n_occ
                    if (channels(k)%irreps(i) /= channels(k)%irreps(n_occ + j)) cycle
                    eigval_pairs = [eigval_pairs, &
                                    shell_scale * (channels(k)%virt_eigvals(j) - &
                                                   channels(k)%occ_eigvals(i))]
                end do
            end do
        end do

    end function ref_hess_eigval_pairs_mo

    function ref_rotate_mo_coeff(kappa, mo_coeff, channels) result(rot_mo_coeff)
        !
        ! this function independently reproduces the rotation of orthonormal MO
        ! coefficients of every particle channel as C exp(K)^T with K the antisymmetric
        ! matrix whose occupied-virtual block is given by kappa
        !
        use otr_mo, only: mo_channel_type

        real(rp), intent(in) :: kappa(:), mo_coeff(:, :, :)
        type(mo_channel_type), intent(in) :: channels(:)
        real(rp) :: &
            rot_mo_coeff(size(mo_coeff, 1), size(mo_coeff, 2), size(mo_coeff, 3))

        real(rp) :: kappa_full(size(mo_coeff, 2), size(mo_coeff, 2), size(channels))
        integer(ip) :: k

        kappa_full = ref_unpack_ov(kappa, channels)
        do k = 1, size(channels, kind=ip)
            rot_mo_coeff(:, :, k) = matmul(mo_coeff(:, :, k), transpose( &
                ref_expm(kappa_full(:, :, k) - transpose(kappa_full(:, :, k)))))
        end do

    contains

        function ref_expm(a) result(exp_a)
            !
            ! this function independently reproduces the exponential of a small matrix
            ! by summing its Taylor series until the terms vanish at working precision
            !
            use opentrustregion_unit_tests, only: identity_matrix

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
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names

        type(mo_type), allocatable, target :: mo
        logical, intent(in) :: to_eigenbasis
        logical :: passed

        real(rp), allocatable :: x(:), rotated(:)
        integer(ip) :: i_case
        character(len=:), allocatable :: test_name

        ! assume test passes
        passed = .true.
        test_name = trim(merge("test_rotate_to_hess_eigenbasis_mo  ", &
                               "test_rotate_from_hess_eigenbasis_mo", to_eigenbasis))

        ! call routine for every occupation case and determine if the rotation matches
        do i_case = 1, n_cases
            call setup_minimal_mo_object( &
                mo, case_n_occ(:case_n_particle(i_case), i_case), &
                case_irreps(:, :case_n_particle(i_case), i_case))
            allocate(x(mo%n_param))
            call random_number(x)
            if (to_eigenbasis) then
                rotated = mo%rotate_to_hess_eigenbasis(x)
            else
                rotated = mo%rotate_from_hess_eigenbasis(x)
            end if
            if (norm2(rotated - ref_rotate_eigenbasis_mo(x, mo%mo_channels, n_mo, &
                                                         to_eigenbasis)) > tol) then
                write(stderr, *) test_name//" failed: Incorrect rotation for the "// &
                    trim(case_names(i_case))//" case."
                passed = .false.
            end if
            deallocate(x, mo)
        end do

    end function check_rotate_hess_eigenbasis_mo

    logical(c_bool) function test_mo_factory_cs() bind(C)
        !
        ! this function tests the subroutine which returns the MO orbital updating
        ! function and objective function for the closed-shell case, which hands the MO
        ! coefficients of the caller as the only particle channel, the irreps of their
        ! MOs, if given, as those of that channel and the closed-shell density matrix
        ! evaluating function to the common setup
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_state_type
        use otr_mo, only: mo_factory_cs, mo_type, mo_settings_type
        use otr_common, only: evaluate_dm_cs_type, extension_solver_settings_type, &
                              extension_level_type, update_orbs_callback, &
                              obj_func_callback, error_mo
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_mo_test_reference, only: n_mo, case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_occ, operator(==)
        use otr_common_unit_tests, only: mock_evaluate_dm_cs, wired_level

        type(test_context_type), target :: context
        integer(ip), parameter :: n_particle = 1

        real(rp), target :: mo_coeff(n_ao, n_mo)
        real(rp) :: ao_overlap(n_ao, n_ao), mo_coeff_3d(n_ao, n_mo, n_particle)
        integer(ip) :: i_sym, error
        type(mo_settings_type) :: settings
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), pointer :: obj_func_mo_funptr
        type(extension_solver_settings_type), target :: solver_settings
        type(solver_state_type), target :: state
        type(extension_level_type), pointer :: level

        ! assume tests pass
        test_mo_factory_cs = .true.

        ! setup settings objects
        call setup_settings(settings, context)

        ! initialize random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff_3d = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        mo_coeff = mo_coeff_3d(:, :, 1)

        ! initialize callback function pointers
        evaluate_dm_funptr => mock_evaluate_dm_cs

        ! call routine with the irreps of the MOs of the symmetric closed-shell
        ! occupation case and determine if the role callback functions are returned and
        ! the solver settings wired with an MO object, the remaining setup is covered
        ! by the test of the common setup
        i_sym = findloc(case_names, "closed-shell with symmetry", dim=1, kind=ip)
        call mo_factory_cs(mo_coeff, ao_overlap, n_occ(1), n_particle, n_ao, n_mo, &
                           evaluate_dm_funptr, update_orbs_mo_funptr, &
                           obj_func_mo_funptr, solver_settings, state, error, &
                           settings, case_irreps(:, 1, i_sym))
        if (error /= 0) then
            write(stderr, *) "test_mo_factory_cs failed: Produced error."
            test_mo_factory_cs = .false.
        end if
        if (.not. associated(update_orbs_mo_funptr, update_orbs_callback) .or. &
            .not. associated(obj_func_mo_funptr, obj_func_callback)) then
            write(stderr, *) "test_mo_factory_cs failed: Role callback functions "// &
                "not returned."
            test_mo_factory_cs = .false.
        end if
        level => wired_level(solver_settings)
        if (.not. associated(level)) then
            write(stderr, *) "test_mo_factory_cs failed: Solver settings not wired "// &
                "with an object."
            test_mo_factory_cs = .false.
        else
            ! determine if the MO object points to the MO coefficients of the caller as
            ! its only particle channel, whose occupation and irreps are those given,
            ! reads the record of the state and holds the given settings and only the
            ! closed-shell density matrix evaluating function
            select type (mo => level%extension)
            type is (mo_type)
                if (.not. associated(mo%mo_coeff)) then
                    write(stderr, *) "test_mo_factory_cs failed: MO coefficients "// &
                        "not associated."
                    test_mo_factory_cs = .false.
                else if (any(shape(mo%mo_coeff) /= [n_ao, n_mo, n_particle])) then
                    write(stderr, *) "test_mo_factory_cs failed: MO coefficients "// &
                        "not associated as a single particle channel."
                    test_mo_factory_cs = .false.
                else
                    mo_coeff(1, 1) = mo_coeff(1, 1) + 1.0_rp
                    if (norm2(mo%mo_coeff(:, :, 1) - mo_coeff) > tol) then
                        write(stderr, *) "test_mo_factory_cs failed: MO "// &
                            "coefficients of the caller not taken over."
                        test_mo_factory_cs = .false.
                    end if
                end if
                if (size(mo%mo_channels) /= n_particle) then
                    write(stderr, *) "test_mo_factory_cs failed: Not a single "// &
                        "particle channel."
                    test_mo_factory_cs = .false.
                else
                    if (mo%mo_channels(1)%n_occ /= n_occ(1)) then
                        write(stderr, *) "test_mo_factory_cs failed: Occupation "// &
                            "not passed on."
                        test_mo_factory_cs = .false.
                    end if
                    if (any(mo%mo_channels(1)%irreps /= case_irreps(:, 1, i_sym))) then
                        write(stderr, *) "test_mo_factory_cs failed: Irreps not "// &
                            "passed on."
                        test_mo_factory_cs = .false.
                    end if
                end if
                if (.not. associated(mo%state, state)) then
                    write(stderr, *) "test_mo_factory_cs failed: Record of the "// &
                        "state not passed on."
                    test_mo_factory_cs = .false.
                end if
                if (.not. (mo%settings == settings%orbital_settings_type)) then
                    write(stderr, *) "test_mo_factory_cs failed: Settings not "// &
                        "passed on."
                    test_mo_factory_cs = .false.
                end if
                if (.not. associated(mo%evaluate_dm_cs, mock_evaluate_dm_cs) .or. &
                    associated(mo%evaluate_dm_os)) then
                    write(stderr, *) "test_mo_factory_cs failed: Closed-shell "// &
                        "density matrix evaluating function not handed over."
                    test_mo_factory_cs = .false.
                end if
            class default
                write(stderr, *) "test_mo_factory_cs failed: Solver settings not "// &
                    "wired with an MO object."
                test_mo_factory_cs = .false.
            end select
        end if

        ! call routine for dimensions which do not match the MO coefficients and
        ! determine if the error of the common setup is passed on
        call mo_factory_cs(mo_coeff, ao_overlap, n_occ(1), n_particle, n_ao, n_mo + 1, &
                           evaluate_dm_funptr, update_orbs_mo_funptr, &
                           obj_func_mo_funptr, solver_settings, state, error, settings)
        if (error /= error_mo + 1) then
            write(stderr, *) "test_mo_factory_cs failed: Error of the common setup "// &
                "not passed on."
            test_mo_factory_cs = .false.
        end if

    end function test_mo_factory_cs

    logical(c_bool) function test_mo_factory_os() bind(C)
        !
        ! this function tests the subroutine which returns the MO orbital updating
        ! function and objective function for the open-shell case, which hands the MO
        ! coefficients of the caller, the irreps of their MOs, if given, and the
        ! open-shell density matrix evaluating function to the common setup
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_state_type
        use otr_mo, only: mo_factory_os, mo_type, mo_settings_type
        use otr_common, only: evaluate_dm_os_type, extension_solver_settings_type, &
                              extension_level_type, update_orbs_callback, &
                              obj_func_callback, error_mo
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_mo_test_reference, only: n_mo, case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_particle, n_occ, operator(==)
        use otr_common_unit_tests, only: mock_evaluate_dm_os, wired_level

        type(test_context_type), target :: context
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: i_sym, error
        type(mo_settings_type) :: settings
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), pointer :: obj_func_mo_funptr
        type(extension_solver_settings_type), target :: solver_settings
        type(solver_state_type), target :: state
        type(extension_level_type), pointer :: level

        ! assume tests pass
        test_mo_factory_os = .true.

        ! setup settings objects
        call setup_settings(settings, context)

        ! initialize random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)

        ! initialize callback function pointers
        evaluate_dm_funptr => mock_evaluate_dm_os

        ! call routine with the irreps of the MOs of the symmetric open-shell
        ! occupation case and determine if the role callback functions are returned and
        ! the solver settings wired with an MO object, the remaining setup is covered
        ! by the test of the common setup
        i_sym = findloc(case_names, "open-shell with symmetry", dim=1, kind=ip)
        call mo_factory_os(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                           evaluate_dm_funptr, update_orbs_mo_funptr, &
                           obj_func_mo_funptr, solver_settings, state, error, &
                           settings, case_irreps(:, :, i_sym))
        if (error /= 0) then
            write(stderr, *) "test_mo_factory_os failed: Produced error."
            test_mo_factory_os = .false.
        end if
        if (.not. associated(update_orbs_mo_funptr, update_orbs_callback) .or. &
            .not. associated(obj_func_mo_funptr, obj_func_callback)) then
            write(stderr, *) "test_mo_factory_os failed: Role callback functions "// &
                "not returned."
            test_mo_factory_os = .false.
        end if
        level => wired_level(solver_settings)
        if (.not. associated(level)) then
            write(stderr, *) "test_mo_factory_os failed: Solver settings not wired "// &
                "with an object."
            test_mo_factory_os = .false.
        else
            ! determine if the MO object points to the MO coefficients of the caller,
            ! whose occupations and irreps are those given, reads the record of the
            ! state and holds the given settings and only the open-shell density matrix
            ! evaluating function
            select type (mo => level%extension)
            type is (mo_type)
                if (.not. associated(mo%mo_coeff, mo_coeff)) then
                    write(stderr, *) "test_mo_factory_os failed: MO coefficients "// &
                        "of the caller not taken over."
                    test_mo_factory_os = .false.
                end if
                if (size(mo%mo_channels) /= n_particle) then
                    write(stderr, *) "test_mo_factory_os failed: Wrong number of "// &
                        "particle channels."
                    test_mo_factory_os = .false.
                else
                    if (any(mo%mo_channels%n_occ /= n_occ)) then
                        write(stderr, *) "test_mo_factory_os failed: Occupations "// &
                            "not passed on."
                        test_mo_factory_os = .false.
                    end if
                    if (any(mo%mo_channels(1)%irreps /= case_irreps(:, 1, i_sym)) .or. &
                        any(mo%mo_channels(2)%irreps /= case_irreps(:, 2, i_sym))) then
                        write(stderr, *) "test_mo_factory_os failed: Irreps not "// &
                            "passed on."
                        test_mo_factory_os = .false.
                    end if
                end if
                if (.not. associated(mo%state, state)) then
                    write(stderr, *) "test_mo_factory_os failed: Record of the "// &
                        "state not passed on."
                    test_mo_factory_os = .false.
                end if
                if (.not. (mo%settings == settings%orbital_settings_type)) then
                    write(stderr, *) "test_mo_factory_os failed: Settings not "// &
                        "passed on."
                    test_mo_factory_os = .false.
                end if
                if (.not. associated(mo%evaluate_dm_os, mock_evaluate_dm_os) .or. &
                    associated(mo%evaluate_dm_cs)) then
                    write(stderr, *) "test_mo_factory_os failed: Open-shell "// &
                        "density matrix evaluating function not handed over."
                    test_mo_factory_os = .false.
                end if
            class default
                write(stderr, *) "test_mo_factory_os failed: Solver settings not "// &
                    "wired with an MO object."
                test_mo_factory_os = .false.
            end select
        end if

        ! call routine for dimensions which do not match the MO coefficients and
        ! determine if the error of the common setup is passed on
        call mo_factory_os(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo + 1, &
                           evaluate_dm_funptr, update_orbs_mo_funptr, &
                           obj_func_mo_funptr, solver_settings, state, error, settings)
        if (error /= error_mo + 1) then
            write(stderr, *) "test_mo_factory_os failed: Error of the common setup "// &
                "not passed on."
            test_mo_factory_os = .false.
        end if

    end function test_mo_factory_os

    logical(c_bool) function test_mo_factory_common() bind(C)
        !
        ! this function tests the subroutine which sets up a new MO object for the MO
        ! coefficients of the caller and, if given, the irreps of their MOs with the
        ! density matrix evaluating function of either shell, if one is given, after
        ! the sanity check of its settings, wires the solver settings with it and the
        ! record of the state and returns its orbital updating function and objective
        ! function
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_state_type, &
                                   verbosity_warning
        use otr_common, only: evaluate_dm_cs_type, evaluate_dm_os_type, &
                              extension_solver_settings_type, extension_level_type, &
                              update_orbs_callback, obj_func_callback, &
                              evaluate_dm_missing_warning_msg, error_mo
        use otr_mo, only: mo_factory_common, mo_type, mo_settings_type, &
                          default_mo_settings
        use otr_mo_test_reference, only: n_mo, case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, &
                                             operator(==), same_logger_and_context
        use opentrustregion_unit_tests, only: setup_settings, test_context_type
        use otr_common_unit_tests, only: mock_evaluate_dm_cs, mock_evaluate_dm_os, &
                                         wired_level

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao)
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), pointer :: obj_func_mo_funptr
        type(extension_solver_settings_type), allocatable, target :: solver_settings, &
                                                                     other_settings
        type(solver_state_type), target :: state, other_state
        type(mo_settings_type) :: settings, expected_settings
        type(test_context_type), target :: context
        type(extension_level_type), pointer :: level, other_level
        integer(ip) :: n_particle, n_param, i_sym, i_case, k, error
        integer(ip) :: irreps(n_mo, n_particle_ref)
        logical :: owned, stored
        logical, allocatable :: same_irrep(:, :)
        character(len=45), parameter :: case_names_factory(4) = &
            [character(len=45) :: "for the closed-shell case", &
             "for the open-shell case with symmetry", "without settings", &
             "without a density matrix evaluating function"]
        character(len=59), parameter :: case_names_rejected(2) = &
            [character(len=59) :: "mismatching dimensions", &
             "irreps without occupied-virtual pairs of the same irrep"]

        ! assume tests pass
        test_mo_factory_common = .true.

        ! initialize random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)

        ! create objects for the closed shell without irreps, for the open shell with
        ! the irreps of the symmetric open-shell occupation case, for no settings and
        ! without a density matrix evaluating function, which prints a warning at the
        ! verbosity of warnings, each of which
        ! has to be an MO object set up for the MO coefficients of the caller and the
        ! irreps of their MOs with the settings and the density matrix evaluating
        ! function and to wire the solver settings, which own it, with the record of
        ! the state
        i_sym = findloc(case_names, "open-shell with symmetry", dim=1, kind=ip)
        do i_case = 1, size(case_names_factory)
            allocate(solver_settings)
            n_particle = merge(2_ip, 1_ip, i_case == 2)
            call setup_settings(settings, context)
            if (i_case == 4) settings%verbose = verbosity_warning
            expected_settings = settings
            if (i_case == 3) expected_settings = default_mo_settings
            evaluate_dm_cs_funptr => null()
            evaluate_dm_os_funptr => null()
            if (i_case == 2) then
                evaluate_dm_os_funptr => mock_evaluate_dm_os
            else if (i_case /= 4) then
                evaluate_dm_cs_funptr => mock_evaluate_dm_cs
            end if
            irreps = 0
            if (i_case == 2) then
                irreps = case_irreps(:, :, i_sym)
                call mo_factory_common( &
                    mo_coeff(:, :, :n_particle), ao_overlap, n_occ(:n_particle), &
                    n_particle, n_ao, n_mo, evaluate_dm_cs_funptr, &
                    evaluate_dm_os_funptr, update_orbs_mo_funptr, obj_func_mo_funptr, &
                    solver_settings, state, error, settings, irreps(:, :n_particle))
            else if (i_case == 3) then
                call mo_factory_common(mo_coeff(:, :, :n_particle), ao_overlap, &
                                       n_occ(:n_particle), n_particle, n_ao, n_mo, &
                                       evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                                       update_orbs_mo_funptr, obj_func_mo_funptr, &
                                       solver_settings, state, error)
            else
                call mo_factory_common(mo_coeff(:, :, :n_particle), ao_overlap, &
                                       n_occ(:n_particle), n_particle, n_ao, n_mo, &
                                       evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                                       update_orbs_mo_funptr, obj_func_mo_funptr, &
                                       solver_settings, state, error, settings)
            end if
            if (error /= 0) then
                write(stderr, *) "test_mo_factory_common failed: Produced error "// &
                    trim(case_names_factory(i_case))//"."
                test_mo_factory_common = .false.
                deallocate(solver_settings)
                cycle
            end if
            if (.not. associated(update_orbs_mo_funptr, update_orbs_callback) .or. &
                .not. associated(obj_func_mo_funptr, obj_func_callback)) then
                write(stderr, *) "test_mo_factory_common failed: Role callback "// &
                    "functions not returned "//trim(case_names_factory(i_case))//"."
                test_mo_factory_common = .false.
            end if
            level => wired_level(solver_settings)
            if (.not. associated(level)) then
                write(stderr, *) "test_mo_factory_common failed: Solver settings "// &
                    "not wired with an object "//trim(case_names_factory(i_case))//"."
                test_mo_factory_common = .false.
                deallocate(solver_settings)
                cycle
            end if
            owned = .false.
            if (allocated(solver_settings%extensions)) then
                if (size(solver_settings%extensions) == 1) owned = &
                    associated(level%extension, solver_settings%extensions(1)%object)
            end if
            if (.not. owned .or. .not. associated(level%extension%state, state)) then
                write(stderr, *) "test_mo_factory_common failed: Object not owned "// &
                    "by the solver settings or not reading the record of the state "// &
                    trim(case_names_factory(i_case))//"."
                test_mo_factory_common = .false.
            end if
            select type (mo => level%extension)
            type is (mo_type)
                ! determine if the object has the dimensions of the case, which the
                ! remaining checks compare arrays of, the number of parameters counting
                ! the occupied-virtual pairs of the same irrep
                n_param = ref_count_mo_params(n_occ(:n_particle), &
                                              irreps(:, :n_particle))
                if (any([mo%n_ao, mo%n_mo, mo%n_particle, mo%n_param, &
                         size(mo%mo_channels, kind=ip)] /= &
                        [n_ao, n_mo, n_particle, n_param, n_particle])) then
                    write(stderr, *) "test_mo_factory_common failed: Incorrect "// &
                        "dimensions "//trim(case_names_factory(i_case))//"."
                    test_mo_factory_common = .false.
                else
                    ! determine if the occupations, the irreps and parameters, the AO
                    ! overlap matrix, the MO coefficients of the caller, the density
                    ! matrix, the gradient and Hessian diagonal, the roles, the
                    ! unevaluated state and the settings are set up
                    if (any(mo%mo_channels%n_occ /= n_occ(:n_particle)) .or. &
                        any(mo%mo_channels%n_virt /= n_mo - n_occ(:n_particle))) then
                        write(stderr, *) "test_mo_factory_common failed: Incorrect "// &
                            "occupations "//trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                    do k = 1, n_particle
                        same_irrep = &
                            spread(irreps(:n_occ(k), k), 2, n_mo - n_occ(k)) == &
                            spread(irreps(n_occ(k) + 1:, k), 1, n_occ(k))
                        if (any(mo%mo_channels(k)%irreps /= irreps(:, k))) then
                            write(stderr, *) "test_mo_factory_common failed: "// &
                                "Incorrect irreps "// &
                                trim(case_names_factory(i_case))//"."
                            test_mo_factory_common = .false.
                        else if (any(mo%mo_channels(k)%param_mask .neqv. same_irrep)) &
                            then
                            write(stderr, *) "test_mo_factory_common failed: "// &
                                "Incorrect parameters "// &
                                trim(case_names_factory(i_case))//"."
                            test_mo_factory_common = .false.
                        end if
                    end do
                    if (norm2(mo%ao_overlap - ao_overlap) > tol) then
                        write(stderr, *) "test_mo_factory_common failed: AO "// &
                            "overlap matrix not stored "// &
                            trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                    if (.not. associated(mo%mo_coeff, mo_coeff(:, :, :n_particle))) then
                        write(stderr, *) "test_mo_factory_common failed: MO "// &
                            "coefficients of the caller not taken over "// &
                            trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                    if (.not. associated(mo%dm_ao, mo%own_dm_ao)) then
                        write(stderr, *) "test_mo_factory_common failed: Density "// &
                            "matrix not held by the object "// &
                            trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                    do k = 1, n_particle
                        if (norm2(mo%dm_ao(:, :, k) - matmul( &
                            mo_coeff(:, :n_occ(k), k), &
                            transpose(mo_coeff(:, :n_occ(k), k)))) > tol) then
                            write(stderr, *) "test_mo_factory_common failed: "// &
                                "Density matrix not constructed from the occupied "// &
                                "orbitals "//trim(case_names_factory(i_case))//"."
                            test_mo_factory_common = .false.
                        end if
                    end do
                    if (.not. (allocated(mo%grad) .and. allocated(mo%h_diag))) then
                        write(stderr, *) "test_mo_factory_common failed: Gradient "// &
                            "or Hessian diagonal not allocated "// &
                            trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    else if (any([size(mo%grad), size(mo%h_diag)] /= n_param)) then
                        write(stderr, *) "test_mo_factory_common failed: Gradient "// &
                            "or Hessian diagonal not allocated with the number of "// &
                            "parameters "//trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                    if (.not. all([mo%provides%update_orbs, mo%provides%obj_func, &
                                   mo%provides%hess_x, mo%provides%precond, &
                                   mo%provides%precond_pd, &
                                   mo%provides%get_extra_trial_vectors]) .or. &
                        mo%provides%project) then
                        write(stderr, *) "test_mo_factory_common failed: Roles the "// &
                            "object performs itself wrong "// &
                            trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                    if (mo%error_origin /= error_mo) then
                        write(stderr, *) "test_mo_factory_common failed: Error "// &
                            "origin of the object not set "// &
                            trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                    if (.not. (mo%evaluation_stale .and. mo%response_stale .and. &
                               mo%hess_eigen_stale)) then
                        write(stderr, *) "test_mo_factory_common failed: Object "// &
                            "not marked unevaluated "// &
                            trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                    if (.not. ( &
                        mo%settings == expected_settings%orbital_settings_type .and. &
                        same_logger_and_context(mo%settings, expected_settings))) then
                        write(stderr, *) "test_mo_factory_common failed: Settings "// &
                            "not stored "//trim(case_names_factory(i_case))//"."
                        test_mo_factory_common = .false.
                    end if
                end if
                if (i_case == 2) then
                    stored = associated(mo%evaluate_dm_os, mock_evaluate_dm_os) .and. &
                             .not. associated(mo%evaluate_dm_cs)
                else if (i_case == 4) then
                    stored = .not. (associated(mo%evaluate_dm_cs) .or. &
                                    associated(mo%evaluate_dm_os))
                else
                    stored = associated(mo%evaluate_dm_cs, mock_evaluate_dm_cs) .and. &
                             .not. associated(mo%evaluate_dm_os)
                end if
                if (.not. stored) then
                    write(stderr, *) "test_mo_factory_common failed: Density "// &
                        "matrix evaluating function not stored "// &
                        trim(case_names_factory(i_case))//"."
                    test_mo_factory_common = .false.
                end if
            class default
                write(stderr, *) "test_mo_factory_common failed: Solver settings "// &
                    "not wired with an MO object "//trim(case_names_factory(i_case))// &
                    "."
                test_mo_factory_common = .false.
            end select
            if ((i_case == 4) .neqv. (index( &
                context%log_message, " "//evaluate_dm_missing_warning_msg) /= 0)) then
                write(stderr, *) "test_mo_factory_common failed: Warning for a "// &
                    "missing density matrix evaluating function printed wrongly "// &
                    trim(case_names_factory(i_case))//"."
                test_mo_factory_common = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create objects for dimensions which do not match the MO coefficients and for
        ! irreps without an occupied-virtual pair of the same irrep, which the sanity
        ! check rejects, each of which has to return its error, return no functions in
        ! place of those returned before and leave the solver settings untouched
        evaluate_dm_cs_funptr => mock_evaluate_dm_cs
        evaluate_dm_os_funptr => null()
        irreps = 0
        irreps(n_occ(1) + 1:, 1) = 1
        do i_case = 1, size(case_names_rejected)
            allocate(solver_settings)
            call setup_settings(settings, context)
            update_orbs_mo_funptr => update_orbs_callback
            obj_func_mo_funptr => obj_func_callback
            if (i_case == 1) then
                call mo_factory_common(mo_coeff(:, :, :1), ao_overlap, n_occ(:1), &
                                       1_ip, n_ao, n_mo + 1, evaluate_dm_cs_funptr, &
                                       evaluate_dm_os_funptr, update_orbs_mo_funptr, &
                                       obj_func_mo_funptr, solver_settings, state, &
                                       error, settings)
            else
                call mo_factory_common(mo_coeff(:, :, :1), ao_overlap, n_occ(:1), &
                                       1_ip, n_ao, n_mo, evaluate_dm_cs_funptr, &
                                       evaluate_dm_os_funptr, update_orbs_mo_funptr, &
                                       obj_func_mo_funptr, solver_settings, state, &
                                       error, settings, irreps(:, :1))
            end if
            if (error /= error_mo + 1) then
                write(stderr, *) "test_mo_factory_common failed: Error not "// &
                    "returned for "//trim(case_names_rejected(i_case))//"."
                test_mo_factory_common = .false.
            end if
            if (associated(update_orbs_mo_funptr) .or. associated(obj_func_mo_funptr)) &
                then
                write(stderr, *) "test_mo_factory_common failed: Functions "// &
                    "returned for "//trim(case_names_rejected(i_case))//"."
                test_mo_factory_common = .false.
            end if
            if (associated(solver_settings%context) .or. &
                allocated(solver_settings%extensions)) then
                write(stderr, *) "test_mo_factory_common failed: Solver settings "// &
                    "changed for "//trim(case_names_rejected(i_case))//"."
                test_mo_factory_common = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create objects for two separate solver settings, which have to be distinct
        ! objects reading their own record of the state, so that calculations with
        ! separate settings can be used alternately
        allocate(solver_settings, other_settings)
        call mo_factory_common(mo_coeff(:, :, :1), ao_overlap, n_occ(:1), 1_ip, n_ao, &
                               n_mo, evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                               update_orbs_mo_funptr, obj_func_mo_funptr, &
                               solver_settings, state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_mo_factory_common failed: Produced error for "// &
                "the first of two separate solver settings."
            test_mo_factory_common = .false.
        end if
        call mo_factory_common(mo_coeff(:, :, :1), ao_overlap, n_occ(:1), 1_ip, n_ao, &
                               n_mo, evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                               update_orbs_mo_funptr, obj_func_mo_funptr, &
                               other_settings, other_state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_mo_factory_common failed: Produced error for "// &
                "the second of two separate solver settings."
            test_mo_factory_common = .false.
        end if
        level => wired_level(solver_settings)
        other_level => wired_level(other_settings)
        if (.not. (associated(level) .and. associated(other_level))) then
            write(stderr, *) "test_mo_factory_common failed: Separate solver "// &
                "settings not wired."
            test_mo_factory_common = .false.
        else if (associated(level, other_level) .or. &
                 .not. associated(level%extension%state, state) .or. &
                 .not. associated(other_level%extension%state, other_state)) then
            write(stderr, *) "test_mo_factory_common failed: Separate solver "// &
                "settings not wired with distinct objects reading their own state."
            test_mo_factory_common = .false.
        end if
        deallocate(solver_settings, other_settings)

    end function test_mo_factory_common

    logical(c_bool) function test_mo_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! dimensions of orbitals parameterized in the MO basis and of the irreps of
        ! their MOs
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: mo_sanity_check
        use opentrustregion_unit_tests, only: setup_error_logging, test_context_type

        type :: mo_sanity_case_type
            integer(ip) :: mo_coeff_shape(3), ao_overlap_shape(2), n_occ(3), &
                           n_occ_passed, n_particle, n_ao, n_mo
            logical :: valid
            character(len=54) :: name
        end type

        ! every invalid case fails a single check, except where a failed check implies
        ! another: no MO fits into a vanishing number of AOs, and vanishing numbers of
        ! MOs or particle channels leave no occupied-virtual rotation
        type(mo_sanity_case_type), parameter :: cases(*) = [ &
            ! mo_coeff_shape, ao_overlap_shape, n_occ, n_occ_passed, n_particle, n_ao,
            ! n_mo, valid, then name
            mo_sanity_case_type( &
                [4, 3, 1], [4, 4], [1, 0, 0], 1, 1, 4, 3, .true., &
                "valid closed-shell dimensions with fewer MOs than AOs"), &
            mo_sanity_case_type( &
                [3, 3, 2], [3, 3], [2, 0, 0], 2, 2, 3, 3, .true., &
                "valid open-shell dimensions with as many MOs as AOs"), &
            mo_sanity_case_type( &
                [4, 3, 2], [4, 4], [3, 1, 0], 2, 2, 4, 3, .true., &
                "valid open-shell dimensions with fewer MOs than AOs"), &
            mo_sanity_case_type( &
                [0, 1, 1], [0, 0], [1, 0, 0], 1, 1, 0, 1, .false., &
                "vanishing number of AOs"), &
            mo_sanity_case_type( &
                [3, 0, 1], [3, 3], [0, 0, 0], 1, 1, 3, 0, .false., &
                "vanishing number of MOs"), &
            mo_sanity_case_type( &
                [3, 4, 1], [3, 3], [1, 0, 0], 1, 1, 3, 4, .false., &
                "more MOs than AOs"), &
            mo_sanity_case_type( &
                [3, 3, 0], [3, 3], [1, 0, 0], 0, 0, 3, 3, .false., &
                "vanishing number of particle channels"), &
            mo_sanity_case_type( &
                [3, 3, 3], [3, 3], [1, 1, 1], 3, 3, 3, 3, .false., &
                "three particle channels"), &
            mo_sanity_case_type( &
                [3, 3, 2], [3, 3], [1, 0, 0], 1, 2, 3, 3, .false., &
                "single occupation for two particle channels"), &
            mo_sanity_case_type( &
                [3, 3, 1], [3, 3], [1, 1, 0], 2, 1, 3, 3, .false., &
                "two occupations for a single particle channel"), &
            mo_sanity_case_type( &
                [2, 3, 1], [3, 3], [1, 0, 0], 1, 1, 3, 3, .false., &
                "MO coefficients with wrong number of AOs"), &
            mo_sanity_case_type( &
                [3, 2, 1], [3, 3], [1, 0, 0], 1, 1, 3, 3, .false., &
                "MO coefficients with wrong number of MOs"), &
            mo_sanity_case_type( &
                [3, 3, 2], [3, 3], [1, 0, 0], 1, 1, 3, 3, .false., &
                "MO coefficients with wrong number of particle channels"), &
            mo_sanity_case_type( &
                [3, 3, 1], [3, 2], [1, 0, 0], 1, 1, 3, 3, .false., &
                "wrongly shaped AO overlap matrix"), &
            mo_sanity_case_type( &
                [6, 6, 2], [6, 6], [2, -1, 0], 2, 2, 6, 6, .false., &
                "negative number of occupied orbitals"), &
            mo_sanity_case_type( &
                [6, 6, 2], [6, 6], [3, 7, 0], 2, 2, 6, 6, .false., &
                "more occupied orbitals than MOs"), &
            mo_sanity_case_type( &
                [3, 3, 2], [3, 3], [3, 0, 0], 2, 2, 3, 3, .false., &
                "missing occupied-virtual rotations")]

        type(mo_sanity_case_type) :: c
        type(orbital_settings_type) :: settings
        type(test_context_type), target :: context
        integer(ip) :: i_case, error

        ! assume tests pass
        test_mo_sanity_check = .true.

        ! call routine for every case and determine if valid dimensions are accepted
        ! without a message and invalid ones rejected with an error message
        do i_case = 1, size(cases)
            c = cases(i_case)
            call setup_error_logging(settings, context)
            call mo_sanity_check(settings, c%mo_coeff_shape, c%ao_overlap_shape, &
                                 c%n_occ(:c%n_occ_passed), c%n_particle, c%n_ao, &
                                 c%n_mo, error)
            if (.not. check_result(c%valid, trim(c%name))) &
                test_mo_sanity_check = .false.
        end do

        ! call routine for valid dimensions with irreps of the MOs and determine if
        ! irreps with an occupied-virtual pair of the same irrep in any channel are
        ! accepted, and wrongly shaped ones, ones without such a pair and ones given
        ! for occupations which do not fit into the MOs rejected
        if (.not. check_irreps_case([1_ip], reshape([0_ip, 1_ip, 0_ip], [3, 1]), &
                                    .true., "valid closed-shell irreps")) &
            test_mo_sanity_check = .false.
        if (.not. check_irreps_case( &
            [1_ip, 1_ip], reshape([0_ip, 1_ip, 0_ip, &
                                   0_ip, 1_ip, 1_ip], [3, 2]), .true., &
            "valid open-shell irreps with a pair in only one channel")) &
            test_mo_sanity_check = .false.
        if (.not. check_irreps_case([1_ip], reshape([0_ip, 1_ip, 0_ip, 0_ip], [4, 1]), &
                                    .false., "irreps of four MOs")) &
            test_mo_sanity_check = .false.
        if (.not. check_irreps_case([1_ip], reshape([0_ip, 1_ip, 0_ip, &
                                                     0_ip, 1_ip, 0_ip], [3, 2]), &
                                    .false., "irreps of two particle channels")) &
            test_mo_sanity_check = .false.
        if (.not. check_irreps_case( &
            [1_ip], reshape([0_ip, 1_ip, 1_ip], [3, 1]), .false., &
            "closed-shell irreps without occupied-virtual pairs of the same irrep")) &
            test_mo_sanity_check = .false.
        if (.not. check_irreps_case( &
            [1_ip, 1_ip], reshape([0_ip, 1_ip, 1_ip, &
                                   1_ip, 0_ip, 0_ip], [3, 2]), .false., &
            "open-shell irreps without occupied-virtual pairs of the same irrep")) &
            test_mo_sanity_check = .false.
        if (.not. check_irreps_case([4_ip], reshape([0_ip, 1_ip, 0_ip], [3, 1]), &
                                    .false., "irreps for more occupied orbitals "// &
                                    "than MOs")) test_mo_sanity_check = .false.

    contains

        function check_irreps_case(n_occ, orbsym, valid, name) result(passed)
            !
            ! this function calls the sanity check for valid dimensions with three MOs
            ! and AOs, the given occupations and the given irreps and determines if
            ! they are accepted or rejected as expected
            !
            integer(ip), intent(in) :: n_occ(:), orbsym(:, :)
            logical, intent(in) :: valid
            character(len=*), intent(in) :: name
            logical :: passed

            call setup_error_logging(settings, context)
            call mo_sanity_check(settings, [3_ip, 3_ip, size(n_occ, kind=ip)], &
                                 [3_ip, 3_ip], n_occ, size(n_occ, kind=ip), 3_ip, &
                                 3_ip, error, orbsym)
            passed = check_result(valid, name)

        end function check_irreps_case

        function check_result(valid, name) result(passed)
            !
            ! this function determines if the sanity check accepted a valid case
            ! without a message or rejected an invalid one with an error message
            !
            logical, intent(in) :: valid
            character(len=*), intent(in) :: name
            logical :: passed

            ! assume the case passes
            passed = .true.

            ! check the error code and the printed message
            if (valid) then
                if (error /= 0 .or. context%log_message /= "") then
                    write(stderr, *) "test_mo_sanity_check failed: Valid input "// &
                        "rejected for "//name//"."
                    passed = .false.
                end if
            else
                if (error /= 1) then
                    write(stderr, *) "test_mo_sanity_check failed: No error "// &
                        "returned for "//name//"."
                    passed = .false.
                end if
                if (context%log_message == "") then
                    write(stderr, *) "test_mo_sanity_check failed: No error "// &
                        "message printed for "//name//"."
                    passed = .false.
                end if
            end if

        end function check_result

    end function test_mo_sanity_check

    logical(c_bool) function test_obj_func_mo() bind(C)
        !
        ! this function tests the function which defines the energy evaluation in the
        ! MO basis
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use otr_common, only: error_mo
        use test_reference, only: arm_host_context, host_context_reached
        use otr_mo, only: mo_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, shell_names
        use otr_common_unit_tests, only: level_context, mock_requests, &
                                         mock_evaluate_dm_cs, mock_evaluate_dm_os

        type(test_context_type), target :: context
        type(mo_type), allocatable, target :: mo
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), &
                    mo_coeff_start(n_ao, n_mo, n_particle_ref), energy, expected_energy
        real(rp), allocatable :: kappa(:), rot_mo_coeff(:, :, :), dm_ao_start(:, :, :)
        integer(ip) :: n_particle, i, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_obj_func_mo = .true.

        ! initialize random orthonormal MO coefficients in a non-orthonormal AO basis
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
        mo_coeff_start = mo_coeff

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))

            ! set up the MO object
            call setup_mo_object(mo, mo_coeff(:, :, :n_particle), ao_overlap, &
                                 n_occ(:n_particle))
            call setup_settings(mo%settings, context)
            call arm_host_context(mo%settings, context)
            if (n_particle == 1) then
                mo%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                mo%evaluate_dm_os => mock_evaluate_dm_os
            end if
            allocate(kappa(mo%n_param))

            ! call routine without an orbital rotation and determine if only the energy
            ! of the density matrix evaluating function is requested and returned for
            ! the current density matrix
            kappa = 0.0_rp
            mock_requests = [integer(ip) :: ]
            energy = mo%obj_func(kappa, error, level_context(mo, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_obj_func_mo failed: Produced error without "// &
                    "an orbital rotation for the "//case_name//" case."
                test_obj_func_mo = .false.
            end if
            if (abs(energy - sum(mo%dm_ao)) > tol) then
                write(stderr, *) "test_obj_func_mo failed: Incorrect energy "// &
                    "without an orbital rotation for the "//case_name//" case."
                test_obj_func_mo = .false.
            end if
            if (size(mock_requests) /= 1) then
                write(stderr, *) "test_obj_func_mo failed: Density matrix "// &
                    "evaluating function not called exactly once for the "// &
                    case_name//" case."
                test_obj_func_mo = .false.
            end if
            if (any(mock_requests /= 0)) then
                write(stderr, *) "test_obj_func_mo failed: Incorrect outputs "// &
                    "requested from density matrix evaluating function for the "// &
                    case_name//" case."
                test_obj_func_mo = .false.
            end if

            ! call routine with an orbital rotation and determine if only the energy of
            ! the density matrix of the rotated orbitals is requested and returned
            ! while the current orbitals and density matrix are left untouched
            call random_number(kappa)
            kappa = 0.2_rp * (kappa - 0.5_rp)
            rot_mo_coeff = ref_rotate_mo_coeff( &
                kappa, mo_coeff_start(:, :, :n_particle), mo%mo_channels)
            expected_energy = 0.0_rp
            do i = 1, n_particle
                expected_energy = expected_energy + &
                                  sum(matmul(rot_mo_coeff(:, :n_occ(i), i), &
                                             transpose(rot_mo_coeff(:, :n_occ(i), i))))
            end do
            dm_ao_start = mo%dm_ao
            mock_requests = [integer(ip) :: ]
            energy = mo%obj_func(kappa, error, level_context(mo, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_obj_func_mo failed: Produced error for an "// &
                    "orbital rotation for the "//case_name//" case."
                test_obj_func_mo = .false.
            end if
            if (abs(energy - expected_energy) > tol) then
                write(stderr, *) "test_obj_func_mo failed: Incorrect energy for an "// &
                    "orbital rotation for the "//case_name//" case."
                test_obj_func_mo = .false.
            end if
            if (size(mock_requests) /= 1 .or. any(mock_requests /= 0)) then
                write(stderr, *) "test_obj_func_mo failed: Incorrect outputs "// &
                    "requested from density matrix evaluating function for an "// &
                    "orbital rotation for the "//case_name//" case."
                test_obj_func_mo = .false.
            end if
            if (norm2(mo_coeff - mo_coeff_start) > tol) then
                write(stderr, *) "test_obj_func_mo failed: Current orbitals "// &
                    "changed by an orbital rotation for the "//case_name//" case."
                test_obj_func_mo = .false.
            end if
            if (norm2(mo%dm_ao - dm_ao_start) > tol) then
                write(stderr, *) "test_obj_func_mo failed: Current density matrix "// &
                    "changed by an orbital rotation for the "//case_name//" case."
                test_obj_func_mo = .false.
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("obj_func_mo", context)) &
                test_obj_func_mo = .false.

            ! call routine with a density matrix evaluating function that fails and
            ! determine if its error is returned unchanged
            context%mock_error = 3
            energy = mo%obj_func(kappa, error, level_context(mo, 1_ip))
            if (error /= 3) then
                write(stderr, *) "test_obj_func_mo failed: Error of the density "// &
                    "matrix evaluating function not returned for the "//case_name// &
                    " case."
                test_obj_func_mo = .false.
            end if
            context%mock_error = 0

            ! call routine with a rotation that is not finite, which the matrix
            ! exponential refuses, and determine if its error is returned with the
            ! origin of the object
            kappa(1) = ieee_value(1.0_rp, ieee_quiet_nan)
            energy = mo%obj_func(kappa, error, level_context(mo, 1_ip))
            if (error /= error_mo + 1) then
                write(stderr, *) "test_obj_func_mo failed: Error for a rotation "// &
                    "that is not finite not returned with the origin of the object "// &
                    "for the "//case_name//" case."
                test_obj_func_mo = .false.
            end if

            ! call routine without a density matrix evaluating function and determine
            ! if the error of its check is returned
            mo%evaluate_dm_cs => null()
            mo%evaluate_dm_os => null()
            energy = mo%obj_func(kappa, error, level_context(mo, 1_ip))
            if (error /= error_mo + 1) then
                write(stderr, *) "test_obj_func_mo failed: No error returned "// &
                    "without a density matrix evaluating function for the "// &
                    case_name//" case."
                test_obj_func_mo = .false.
            end if

            ! deallocate MO object
            deallocate(mo, kappa)
        end do

    end function test_obj_func_mo

    logical(c_bool) function test_update_orbs_mo() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation in the MO basis, for the closed-shell and the
        ! open-shell case
        !
        use test_reference, only: arm_host_context, host_context_reached
        use otr_mo, only: mo_type
        use otr_common, only: hess_x_callback, error_mo
        use opentrustregion, only: hess_x_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, shell_names
        use otr_common_unit_tests, only: &
            level_context, mock_fock_factor, mock_requests, mock_get_response_cs, &
            mock_get_response_os, mock_evaluate_dm_cs, mock_evaluate_dm_os

        type(test_context_type), target :: context
        type(mo_type), allocatable, target :: mo
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), fock_mo(n_mo, n_mo), func
        real(rp), allocatable :: mo_coeff_start(:, :, :), kappa(:), grad(:), h_diag(:)
        integer(ip) :: n_particle, i, error
        logical :: response_set
        character(len=:), allocatable :: case_name
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_update_orbs_mo = .true.

        ! initialize a non-orthonormal AO basis
        ao_overlap = generate_random_ao_overlap(n_ao)

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))

            ! set up the MO object with random orthonormal MO coefficients
            mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
            call setup_mo_object(mo, mo_coeff(:, :, :n_particle), ao_overlap, &
                                 n_occ(:n_particle))
            call setup_settings(mo%settings, context)
            call arm_host_context(mo%settings, context)
            allocate(kappa(mo%n_param), grad(mo%n_param), h_diag(mo%n_param))
            if (n_particle == 1) then
                mo%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                mo%evaluate_dm_os => mock_evaluate_dm_os
            end if
            mock_requests = [integer(ip) :: ]

            ! call routine without an orbital rotation for starting orbitals which have
            ! not been evaluated yet, with the response marked current so that only the
            ! evaluation flag forces the evaluation, and determine if the energy of the
            ! density matrix, the static part of the Hessian from the Fock matrix
            ! transformed to the MO basis, the gradient, the Hessian diagonal, the
            ! response function and the Hessian linear transformation are obtained
            kappa = 0.0_rp
            mo%response_stale = .false.
            if (.not. check_update_orbs_mo_stage(1_ip, "starting orbitals")) &
                test_update_orbs_mo = .false.
            if (abs(func - sum(mo%dm_ao)) > tol) then
                write(stderr, *) "test_update_orbs_mo failed: Incorrect energy for "// &
                    "the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if
            do i = 1, n_particle
                fock_mo = ref_mo_transform(mo_coeff(:, :, i), &
                                           mock_fock_factor(1) * mo%dm_ao(:, :, i))
                if (norm2(mo%mo_channels(i)%fock_oo - fock_mo(:n_occ(i), :n_occ(i))) > &
                    tol .or. norm2(mo%mo_channels(i)%fock_vv - &
                                   fock_mo(n_occ(i) + 1:, n_occ(i) + 1:)) > tol) then
                    write(stderr, *) "test_update_orbs_mo failed: Static Hessian "// &
                        "part not built from the Fock matrix in the MO basis for "// &
                        "the "//case_name//" case."
                    test_update_orbs_mo = .false.
                end if
            end do
            if (norm2(grad - mo%grad) > tol) then
                write(stderr, *) "test_update_orbs_mo failed: Gradient not "// &
                    "returned for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if
            if (norm2(h_diag - mo%h_diag) > tol) then
                write(stderr, *) "test_update_orbs_mo failed: Hessian diagonal not "// &
                    "returned for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if
            if (n_particle == 1) then
                response_set = associated(mo%get_response_cs, mock_get_response_cs)
            else
                response_set = associated(mo%get_response_os, mock_get_response_os)
            end if
            if (.not. response_set) then
                write(stderr, *) "test_update_orbs_mo failed: Response function "// &
                    "not stored for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if
            if (.not. associated(hess_x_funptr, hess_x_callback)) then
                write(stderr, *) "test_update_orbs_mo failed: Returned Hessian "// &
                    "linear transformation is wrong for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if
            if (mo%evaluation_stale) then
                write(stderr, *) "test_update_orbs_mo failed: Evaluation still "// &
                    "marked stale after the quantities were computed for the "// &
                    case_name//" case."
                test_update_orbs_mo = .false.
            end if

            ! call routine again without an orbital rotation after clearing the outputs
            ! and setting a vanishing stored energy, and determine if the quantities of
            ! the evaluated orbitals are returned without an evaluation, even for the
            ! vanishing energy
            mo%energy = 0.0_rp
            grad = 0.0_rp
            h_diag = 0.0_rp
            if (.not. check_update_orbs_mo_stage(1_ip, "evaluated orbitals")) &
                test_update_orbs_mo = .false.
            if (abs(func) > tol) then
                write(stderr, *) "test_update_orbs_mo failed: Stored energy not "// &
                    "returned for evaluated orbitals for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if
            if (norm2(grad - mo%grad) > tol) then
                write(stderr, *) "test_update_orbs_mo failed: Stored gradient not "// &
                    "returned for evaluated orbitals for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if
            if (norm2(h_diag - mo%h_diag) > tol) then
                write(stderr, *) "test_update_orbs_mo failed: Stored Hessian "// &
                    "diagonal not returned for evaluated orbitals for the "// &
                    case_name//" case."
                test_update_orbs_mo = .false.
            end if

            ! mark the response as stale, as an approximate-Hessian extension such as
            ! ARH would after moving the orbitals without going through this routine,
            ! and call again without an orbital rotation, and determine if this forces
            ! an evaluation and clears the flag
            mo%response_stale = .true.
            if (.not. check_update_orbs_mo_stage(2_ip, "a stale response")) &
                test_update_orbs_mo = .false.
            if (mo%response_stale) then
                write(stderr, *) "test_update_orbs_mo failed: Response still "// &
                    "marked stale after being recomputed for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if

            ! mark the evaluation as stale, as the factory does for new starting
            ! orbitals, and call again without an orbital rotation, and determine if
            ! this forces an evaluation and clears the flag
            mo%evaluation_stale = .true.
            if (.not. check_update_orbs_mo_stage(3_ip, "a stale evaluation")) &
                test_update_orbs_mo = .false.
            if (mo%evaluation_stale) then
                write(stderr, *) "test_update_orbs_mo failed: Evaluation still "// &
                    "marked stale after being recomputed for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if

            ! call routine with a random orbital rotation of unequal elements of either
            ! sign, whose sum is negative, at evaluated orbitals, so that only the
            ! rotation itself calls for an evaluation, keeping the MO coefficients it
            ! starts from, and determine if the orbitals are rotated, consistently with
            ! the density matrix, and evaluated, and if every evaluation requested the
            ! Fock matrix and the response function
            call random_number(kappa)
            kappa = 0.2_rp * (kappa - 0.5_rp)
            kappa = kappa - sum(kappa) / size(kappa) - 0.01_rp
            mo_coeff_start = mo%mo_coeff
            if (.not. check_update_orbs_mo_stage(4_ip, "an orbital rotation")) &
                test_update_orbs_mo = .false.
            if (norm2(mo%mo_coeff - ref_rotate_mo_coeff(kappa, mo_coeff_start, &
                                                        mo%mo_channels)) > tol) then
                write(stderr, *) "test_update_orbs_mo failed: Orbitals not rotated "// &
                    "correctly by an orbital rotation for the "//case_name//" case."
                test_update_orbs_mo = .false.
            end if
            do i = 1, n_particle
                if (norm2(mo%dm_ao(:, :, i) - &
                          matmul(mo%mo_coeff(:, :n_occ(i), i), &
                                 transpose(mo%mo_coeff(:, :n_occ(i), i)))) > tol) then
                    write(stderr, *) "test_update_orbs_mo failed: Density matrix "// &
                        "not moved consistently with the orbitals for the "// &
                        case_name//" case."
                    test_update_orbs_mo = .false.
                end if
            end do
            if (any(mock_requests /= 3)) then
                write(stderr, *) "test_update_orbs_mo failed: Incorrect outputs "// &
                    "requested from density matrix evaluating function for the "// &
                    case_name//" case."
                test_update_orbs_mo = .false.
            end if

            ! call routine with the orbital rotation again, whose evaluation fails, and
            ! determine if the error is passed on and the rotated orbitals are marked
            ! as not evaluated
            context%mock_error = 1
            call mo%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                level_context(mo, 1_ip))
            if (error /= 1) then
                write(stderr, *) "test_update_orbs_mo failed: Error of the density "// &
                    "matrix evaluating function not passed on for the "//case_name// &
                    " case."
                test_update_orbs_mo = .false.
            end if
            if (.not. mo%evaluation_stale) then
                write(stderr, *) "test_update_orbs_mo failed: Evaluation not "// &
                    "marked stale after a failed evaluation for the "//case_name// &
                    " case."
                test_update_orbs_mo = .false.
            end if
            context%mock_error = 0

            ! call routine without a density matrix evaluating function and determine
            ! if the error of its check is returned
            mo%evaluate_dm_cs => null()
            mo%evaluate_dm_os => null()
            call mo%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                level_context(mo, 1_ip))
            if (error /= error_mo + 1) then
                write(stderr, *) "test_update_orbs_mo failed: No error returned "// &
                    "without a density matrix evaluating function for the "// &
                    case_name//" case."
                test_update_orbs_mo = .false.
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("update_orbs_mo", context)) &
                test_update_orbs_mo = .false.

            ! deallocate MO object and case arrays
            deallocate(mo, kappa, grad, h_diag)
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
            character(len=*), intent(in) :: stage
            logical :: passed

            ! assume the stage passes
            passed = .true.

            ! call routine and determine if it produces an error and if the density
            ! matrix evaluating function was called the expected number of times
            call mo%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                level_context(mo, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_update_orbs_mo failed: Produced error for "// &
                    stage//" for the "//case_name//" case."
                passed = .false.
            end if
            if (size(mock_requests) /= n_calls) then
                write(stderr, *) "test_update_orbs_mo failed: Density matrix "// &
                    "evaluating function not called the expected number of times "// &
                    "for "//stage//" for the "//case_name//" case."
                passed = .false.
            end if

        end function check_update_orbs_mo_stage

    end function test_update_orbs_mo

    logical(c_bool) function test_hess_x_mo() bind(C)
        !
        ! this function tests the subroutine which defines the Hessian linear
        ! transformation in the MO basis
        !
        use otr_common, only: error_mo
        use test_reference, only: arm_host_context, host_context_reached
        use otr_mo, only: mo_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              generate_random_symm_matrix
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, shell_names
        use otr_common_unit_tests, only: &
            level_context, mock_requests, mock_response_factor, &
            mock_response_coupling, model_response_scale, mock_get_response_cs, &
            mock_get_coupled_response_os, mock_get_stale_response_os, &
            model_response_cs, model_response_os, mock_evaluate_dm_cs, &
            mock_evaluate_dm_os

        type(test_context_type), target :: context
        type(mo_type), allocatable, target :: mo
        ! step of the finite difference and its tolerance relative to the second
        ! derivative it approximates
        real(rp), parameter :: step = 1e-3_rp, fd_tol = 1e-5_rp

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), model_h(n_ao, n_ao), fock_mo(n_mo, n_mo), &
                    model_scale, second_diff
        real(rp), allocatable :: x(:), hess_x(:), expected_hess_x(:)
        integer(ip) :: n_particle, i_case, i, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_hess_x_mo = .true.

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
            ! function, which couples the particle channels in the open-shell cases
            call setup_mo_object(mo, mo_coeff(:, :, :n_particle), ao_overlap, &
                                 case_n_occ(:n_particle, i_case), &
                                 case_irreps(:, :n_particle, i_case))
            call setup_settings(mo%settings, context)
            call arm_host_context(mo%settings, context)
            call setup_random_mo_channels(mo, case_n_occ(:n_particle, i_case), n_mo)
            if (n_particle == 1) then
                mo%evaluate_dm_cs => mock_evaluate_dm_cs
                mo%get_response_cs => mock_get_response_cs
            else
                mo%evaluate_dm_os => mock_evaluate_dm_os
                mo%get_response_os => mock_get_coupled_response_os
            end if
            mo%response_stale = .false.

            ! call routine for a random trial vector and determine if values of
            ! resulting Hessian linear transformation match
            allocate(x(mo%n_param), hess_x(mo%n_param))
            call random_number(x)
            expected_hess_x = ref_hess_x_mo( &
                x, mo%mo_channels, mo_coeff(:, :, :n_particle), mock_response_factor, &
                mock_response_coupling)
            call mo%hess_x(x, hess_x, error, level_context(mo, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_hess_x_mo failed: Produced error for the "// &
                    case_name//" case."
                test_hess_x_mo = .false.
            end if
            if (norm2(hess_x - expected_hess_x) > tol) then
                write(stderr, *) "test_hess_x_mo failed: Incorrect Hessian linear "// &
                    "transformation for the "//case_name//" case."
                test_hess_x_mo = .false.
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("hess_x_mo", context)) &
                test_hess_x_mo = .false.

            ! call routine with a response function that fails and determine if its
            ! error is returned unchanged
            context%mock_error = 3
            call mo%hess_x(x, hess_x, error, level_context(mo, 1_ip))
            if (error /= 3) then
                write(stderr, *) "test_hess_x_mo failed: Error of the response "// &
                    "function not returned for the "//case_name//" case."
                test_hess_x_mo = .false.
            end if
            context%mock_error = 0

            deallocate(mo, x, hess_x)
        end do

        ! set up the open-shell MO object again, mark the response as stale, as ARH
        ! would after moving the orbitals without going through the orbital updating
        ! routine, replace the response function by one of outdated orbitals, and
        ! determine if this triggers the density matrix evaluating function to refresh
        ! the response, which is then used
        call setup_mo_object(mo, mo_coeff, ao_overlap, n_occ)
        call setup_settings(mo%settings, context)
        call arm_host_context(mo%settings, context)
        call setup_random_mo_channels(mo, n_occ, n_mo)
        mo%evaluate_dm_os => mock_evaluate_dm_os
        mo%get_response_os => mock_get_stale_response_os
        mo%response_stale = .true.
        allocate(x(mo%n_param), hess_x(mo%n_param))
        call random_number(x)
        expected_hess_x = &
            ref_hess_x_mo(x, mo%mo_channels, mo_coeff, mock_response_factor, 0.0_rp)
        mock_requests = [integer(ip) :: ]
        call mo%hess_x(x, hess_x, error, level_context(mo, 1_ip))
        if (error /= 0) then
            write(stderr, *) "test_hess_x_mo failed: Produced error while "// &
                "refreshing a stale response."
            test_hess_x_mo = .false.
        end if
        if (size(mock_requests) /= 1) then
            write(stderr, *) "test_hess_x_mo failed: Stale response was not refreshed."
            test_hess_x_mo = .false.
        end if
        if (any(mock_requests /= 2)) then
            write(stderr, *) "test_hess_x_mo failed: Incorrect outputs requested "// &
                "from density matrix evaluating function while refreshing a stale "// &
                "response."
            test_hess_x_mo = .false.
        end if
        if (mo%response_stale) then
            write(stderr, *) "test_hess_x_mo failed: Response still marked stale "// &
                "after being refreshed."
            test_hess_x_mo = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write(stderr, *) "test_hess_x_mo failed: Incorrect Hessian linear "// &
                "transformation with a refreshed response."
            test_hess_x_mo = .false.
        end if

        ! check that the host context reached the callback functions unchanged
        if (.not. host_context_reached("hess_x_mo", context)) test_hess_x_mo = .false.

        ! call routine without a density matrix evaluating function and determine if
        ! the error of its check is returned
        mo%evaluate_dm_cs => null()
        mo%evaluate_dm_os => null()
        call mo%hess_x(x, hess_x, error, level_context(mo, 1_ip))
        if (error /= error_mo + 1) then
            write(stderr, *) "test_hess_x_mo failed: No error returned without a "// &
                "density matrix evaluating function."
            test_hess_x_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo, x, hess_x)

        ! set up the MO object with the Fock matrix and the response of a quadratic
        ! model energy at its current orbitals, and determine if the quadratic form of
        ! the Hessian linear transformation along a random trial vector matches the
        ! second finite difference of the model energy along it, which checks the
        ! static and the response part together with their scaling independently of how
        ! the routine assembles them
        model_h = generate_random_symm_matrix(n_ao)
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            model_scale = model_response_scale(n_particle)
            call setup_mo_object(mo, mo_coeff(:, :, :n_particle), ao_overlap, &
                                 n_occ(:n_particle))
            call setup_settings(mo%settings, context)
            do i = 1, n_particle
                fock_mo = ref_mo_transform(mo_coeff(:, :, i), &
                                           model_h + model_scale * mo%dm_ao(:, :, i))
                mo%mo_channels(i)%fock_oo = fock_mo(:n_occ(i), :n_occ(i))
                mo%mo_channels(i)%fock_vv = fock_mo(n_occ(i) + 1:, n_occ(i) + 1:)
            end do
            if (n_particle == 1) then
                mo%evaluate_dm_cs => mock_evaluate_dm_cs
                mo%get_response_cs => model_response_cs
            else
                mo%evaluate_dm_os => mock_evaluate_dm_os
                mo%get_response_os => model_response_os
            end if
            mo%response_stale = .false.
            allocate(x(mo%n_param), hess_x(mo%n_param))
            call random_number(x)
            x = x - 0.5_rp
            call mo%hess_x(x, hess_x, error, level_context(mo, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_hess_x_mo failed: Produced error for the "// &
                    "quadratic model energy for the "//case_name//" case."
                test_hess_x_mo = .false.
            end if
            second_diff = (model_energy(step * x) - 2.0_rp * &
                           model_energy(0.0_rp * x) + model_energy(-step * x)) / step**2
            if (abs(dot_product(x, hess_x) - second_diff) > &
                fd_tol * max(1.0_rp, abs(second_diff))) then
                write(stderr, *) "test_hess_x_mo failed: Hessian linear "// &
                    "transformation does not match the finite difference of the "// &
                    "quadratic model energy for the "//case_name//" case."
                test_hess_x_mo = .false.
            end if
            deallocate(mo, x, hess_x)
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
                                               mo%mo_channels)
            energy = 0.0_rp
            do k = 1, n_particle
                dm = matmul(rot_mo_coeff(:, :n_occ(k), k), &
                            transpose(rot_mo_coeff(:, :n_occ(k), k)))
                energy = energy + merge(2.0_rp, 1.0_rp, n_particle == 1) * &
                         (sum(model_h * dm) + 0.5_rp * model_scale * sum(dm * dm))
            end do

        end function model_energy

    end function test_hess_x_mo

    logical(c_bool) function test_rotate_orbitals_mo() bind(C)
        !
        ! this function tests the subroutine which moves the current orbitals by an
        ! orbital rotation in the MO basis
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use otr_common, only: error_mo
        use otr_mo, only: mo_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_mo_test_reference, only: n_mo, n_param
        use otr_common_test_reference, only: n_ao, n_occ, n_particle

        type(test_context_type), target :: context
        type(mo_type), allocatable, target :: mo
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), expected(n_ao, n_mo, n_particle), &
                    kappa(n_param)
        integer(ip) :: error, k

        ! assume tests pass
        test_rotate_orbitals_mo = .true.

        ! set up the MO object with random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_mo_object(mo, mo_coeff, ao_overlap, n_occ)
        call setup_settings(mo%settings, context)

        ! call routine and determine if the MO coefficients of the caller and the
        ! density matrix are rotated and the response, which was not rebuilt, is marked
        ! stale
        call random_number(kappa)
        kappa = 0.2_rp * (kappa - 0.5_rp)
        expected = ref_rotate_mo_coeff(kappa, mo_coeff, mo%mo_channels)
        mo%response_stale = .false.
        call mo%rotate_orbitals(kappa, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_orbitals_mo failed: Produced error."
            test_rotate_orbitals_mo = .false.
        end if
        if (.not. mo%response_stale) then
            write(stderr, *) "test_rotate_orbitals_mo failed: Response not marked "// &
                "stale."
            test_rotate_orbitals_mo = .false.
        end if
        if (norm2(mo_coeff - expected) > tol) then
            write(stderr, *) "test_rotate_orbitals_mo failed: MO coefficients not "// &
                "rotated in place."
            test_rotate_orbitals_mo = .false.
        end if
        do k = 1, n_particle
            if (norm2(mo%dm_ao(:, :, k) - matmul( &
                expected(:, :n_occ(k), k), transpose(expected(:, :n_occ(k), k)))) > &
                tol) then
                write(stderr, *) "test_rotate_orbitals_mo failed: Density matrix "// &
                    "not updated."
                test_rotate_orbitals_mo = .false.
            end if
        end do

        ! call routine with a rotation that is not finite, which the matrix exponential
        ! refuses, and determine if its error is returned with the origin of the object
        ! while the orbitals stay in place
        expected = mo_coeff
        kappa(1) = ieee_value(1.0_rp, ieee_quiet_nan)
        call mo%rotate_orbitals(kappa, error)
        if (error /= error_mo + 1) then
            write(stderr, *) "test_rotate_orbitals_mo failed: Error for a rotation "// &
                "that is not finite not returned with the origin of the object."
            test_rotate_orbitals_mo = .false.
        end if
        if (norm2(mo_coeff - expected) > tol) then
            write(stderr, *) "test_rotate_orbitals_mo failed: Orbitals moved by a "// &
                "rotation that is not finite."
            test_rotate_orbitals_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo)

    end function test_rotate_orbitals_mo

    logical(c_bool) function test_calculate_grad_h_diag_mo() bind(C)
        !
        ! this function tests the subroutine which calculates the gradient, Hessian
        ! diagonal and static part of the Hessian at the current density for the MO
        ! basis
        !
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_occ_ref => n_occ, &
                                             n_particle_ref => n_particle
        use opentrustregion_unit_tests, only: generate_random_symm_matrix

        real(rp), parameter :: fd_step = 1e-4_rp

        type(mo_type), allocatable, target :: mo
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), fock_ao(n_ao, n_ao, n_particle_ref), &
                    fock_mo(n_mo, n_mo, n_particle_ref), &
                    diag_diff(n_mo, n_mo, n_particle_ref), shell_scale, &
                    rot_mo_coeff(n_ao, n_mo, 1), energy
        real(rp), allocatable :: expected_h_diag(:), kappa(:), fd_grad(:)
        integer(ip), allocatable :: n_occ(:)
        integer(ip) :: n_particle, n_occ_k, i, j, k, i_case, i_sign
        character(len=:), allocatable :: case_name

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
            call setup_mo_object(mo, mo_coeff(:, :, :n_particle), ao_overlap, n_occ, &
                                 case_irreps(:, :n_particle, i_case))
            do k = 1, n_particle
                fock_ao(:, :, k) = generate_random_symm_matrix(n_ao)
                fock_mo(:, :, k) = ref_mo_transform(mo_coeff(:, :, k), fock_ao(:, :, k))
            end do
            mo%hess_eigen_stale = .false.

            ! expected Hessian diagonal from the diagonal Fock matrix elements
            diag_diff = 0.0_rp
            do k = 1, n_particle
                do j = n_occ(k) + 1, n_mo
                    do i = 1, n_occ(k)
                        diag_diff(i, j, k) = shell_scale * &
                                             (fock_mo(j, j, k) - fock_mo(i, i, k))
                    end do
                end do
            end do
            expected_h_diag = ref_pack_ov(diag_diff(:, :, :n_particle), mo%mo_channels)

            ! call routine and determine if the gradient, Hessian diagonal and Fock
            ! matrix blocks are obtained in the MO basis and the eigendecomposition is
            ! marked stale
            call mo%calculate_grad_h_diag(fock_ao(:, :, :n_particle))
            if (norm2(mo%grad - shell_scale * ref_pack_ov(fock_mo(:, :, :n_particle), &
                                                          mo%mo_channels)) > tol) then
                write(stderr, *) "test_calculate_grad_h_diag_mo failed: Incorrect "// &
                    "gradient for the "//case_name//" case."
                test_calculate_grad_h_diag_mo = .false.
            end if
            if (norm2(mo%h_diag - expected_h_diag) > tol) then
                write(stderr, *) "test_calculate_grad_h_diag_mo failed: Incorrect "// &
                    "Hessian diagonal for the "//case_name//" case."
                test_calculate_grad_h_diag_mo = .false.
            end if
            do k = 1, n_particle
                n_occ_k = n_occ(k)
                if (norm2(mo%mo_channels(k)%fock_oo - &
                          fock_mo(:n_occ_k, :n_occ_k, k)) > tol) then
                    write(stderr, *) "test_calculate_grad_h_diag_mo failed: "// &
                        "Incorrect occupied-occupied Fock matrix block for the "// &
                        case_name//" case."
                    test_calculate_grad_h_diag_mo = .false.
                end if
                if (norm2(mo%mo_channels(k)%fock_vv - &
                          fock_mo(n_occ_k + 1:, n_occ_k + 1:, k)) > tol) then
                    write(stderr, *) "test_calculate_grad_h_diag_mo failed: "// &
                        "Incorrect virtual-virtual Fock matrix block for the "// &
                        case_name//" case."
                    test_calculate_grad_h_diag_mo = .false.
                end if
            end do
            if (.not. mo%hess_eigen_stale) then
                write(stderr, *) "test_calculate_grad_h_diag_mo failed: "// &
                    "Eigendecomposition not marked stale for the "//case_name//" case."
                test_calculate_grad_h_diag_mo = .false.
            end if

            ! deallocate MO object
            deallocate(mo)
        end do

        ! call routine for a random symmetric closed-shell Fock matrix H, which, unlike
        ! that of the mock, is consistent with an energy, the linear energy 2 tr(H D),
        ! and determine if the gradient is the finite-difference derivative of the
        ! energy along the reference orbital rotation, which pins the rotation
        ! convention and parameter order against the gradient
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
        call setup_mo_object(mo, mo_coeff(:, :, :1), ao_overlap, n_occ_ref(:1))
        fock_ao(:, :, 1) = generate_random_symm_matrix(n_ao)
        call mo%calculate_grad_h_diag(fock_ao(:, :, :1))
        allocate(kappa(mo%n_param), fd_grad(mo%n_param))
        do k = 1, mo%n_param
            fd_grad(k) = 0.0_rp
            do i_sign = -1, 1, 2
                kappa = 0.0_rp
                kappa(k) = i_sign * fd_step
                rot_mo_coeff = ref_rotate_mo_coeff(kappa, mo_coeff(:, :, :1), &
                                                   mo%mo_channels)
                energy = 2.0_rp * sum(fock_ao(:, :, 1) * matmul( &
                    rot_mo_coeff(:, :n_occ_ref(1), 1), &
                    transpose(rot_mo_coeff(:, :n_occ_ref(1), 1))))
                fd_grad(k) = fd_grad(k) + i_sign * energy / (2.0_rp * fd_step)
            end do
        end do
        if (norm2(mo%grad - fd_grad) > 1e-6_rp * norm2(fd_grad)) then
            write(stderr, *) "test_calculate_grad_h_diag_mo failed: Gradient not "// &
                "the derivative of the energy."
            test_calculate_grad_h_diag_mo = .false.
        end if
        deallocate(mo)

    end function test_calculate_grad_h_diag_mo

    logical(c_bool) function test_refresh_hess_eigen_mo() bind(C)
        !
        ! this function tests the subroutine which refreshes the eigendecomposition of
        ! the static part of the Hessian in the MO basis, which diagonalizes the
        ! occupied-occupied and virtual-virtual blocks of the Fock matrix per irrep, if
        ! it is stale
        !
        use otr_mo, only: mo_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix
        use otr_mo_test_reference, only: n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names

        type(test_context_type), target :: context
        type(mo_type), allocatable, target :: mo
        real(rp), allocatable :: eigvals_before(:)
        integer(ip) :: error, k, n_occ, n_virt, i_case, n_particle
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_refresh_hess_eigen_mo = .true.

        ! loop over every occupation case, replace both Fock matrix blocks of every
        ! channel by random ones which also couple orbitals of different irreps,
        ! discard the eigendecomposition, which an object the factory set up does not
        ! have yet, mark it stale and determine if the part of the blocks between
        ! orbitals of the same irrep is diagonalized by orthonormal eigenvectors which
        ! do not mix irreps, ignoring the coupling between irreps
        do i_case = 1, n_cases
            n_particle = case_n_particle(i_case)
            case_name = trim(case_names(i_case))
            call setup_minimal_mo_object(mo, case_n_occ(:n_particle, i_case), &
                                         case_irreps(:, :n_particle, i_case))
            call setup_settings(mo%settings, context)
            do k = 1, n_particle
                associate (channel => mo%mo_channels(k))
                    channel%fock_oo = generate_random_symm_matrix(channel%n_occ)
                    channel%fock_vv = generate_random_symm_matrix(channel%n_virt)
                    deallocate(channel%occ_eigvecs, channel%occ_eigvals, &
                               channel%virt_eigvecs, channel%virt_eigvals)
                end associate
            end do
            mo%hess_eigen_stale = .true.
            call mo%refresh_hess_eigen(error)
            if (error /= 0) then
                write(stderr, *) "test_refresh_hess_eigen_mo failed: Produced "// &
                    "error for the "//case_name//" case."
                test_refresh_hess_eigen_mo = .false.
            end if
            if (mo%hess_eigen_stale) then
                write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                    "Eigendecomposition still stale for the "//case_name//" case."
                test_refresh_hess_eigen_mo = .false.
            end if
            do k = 1, n_particle
                associate (channel => mo%mo_channels(k))
                    n_occ = channel%n_occ
                    n_virt = channel%n_virt
                    if (norm2(matmul(transpose(channel%occ_eigvecs), matmul( &
                        same_irrep_part(channel%fock_oo, channel%irreps(:n_occ)), &
                        channel%occ_eigvecs)) - &
                              diagonal_matrix(channel%occ_eigvals)) > tol) then
                        write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Occupied-occupied block not diagonalized for the "// &
                            case_name//" case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                    if (norm2(matmul(transpose(channel%occ_eigvecs), &
                                     channel%occ_eigvecs) - identity_matrix(n_occ)) > &
                        tol) then
                        write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Occupied eigenvectors not orthonormal for the "// &
                            case_name//" case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                    if (mixes_irreps(channel%occ_eigvecs, channel%irreps(:n_occ))) then
                        write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Occupied eigenvectors mix irreps for the "//case_name// &
                            " case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                    if (norm2(matmul(transpose(channel%virt_eigvecs), matmul( &
                        same_irrep_part(channel%fock_vv, channel%irreps(n_occ + 1:)), &
                        channel%virt_eigvecs)) - &
                              diagonal_matrix(channel%virt_eigvals)) > tol) then
                        write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Virtual-virtual block not diagonalized for the "// &
                            case_name//" case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                    if (norm2(matmul(transpose(channel%virt_eigvecs), &
                                     channel%virt_eigvecs) - &
                              identity_matrix(n_virt)) > tol) then
                        write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Virtual eigenvectors not orthonormal for the "// &
                            case_name//" case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                    if (mixes_irreps(channel%virt_eigvecs, &
                                     channel%irreps(n_occ + 1:))) then
                        write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                            "Virtual eigenvectors mix irreps for the "//case_name// &
                            " case."
                        test_refresh_hess_eigen_mo = .false.
                    end if
                end associate
            end do

            ! double the occupied block of the first channel, which is occupied in
            ! every case, and determine if an up-to-date eigendecomposition is not
            ! recomputed
            eigvals_before = mo%mo_channels(1)%occ_eigvals
            mo%mo_channels(1)%fock_oo = 2.0_rp * mo%mo_channels(1)%fock_oo
            call mo%refresh_hess_eigen(error)
            if (error /= 0) then
                write(stderr, *) "test_refresh_hess_eigen_mo failed: Produced "// &
                    "error without stale eigendecomposition for the "//case_name// &
                    " case."
                test_refresh_hess_eigen_mo = .false.
            end if
            if (norm2(mo%mo_channels(1)%occ_eigvals - eigvals_before) > tol) then
                write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                    "Eigendecomposition recomputed although not stale for the "// &
                    case_name//" case."
                test_refresh_hess_eigen_mo = .false.
            end if

            ! mark the eigendecomposition stale and determine if the one computed
            ! before is replaced by that of the doubled block, whose eigenvalues double
            mo%hess_eigen_stale = .true.
            call mo%refresh_hess_eigen(error)
            if (error /= 0) then
                write(stderr, *) "test_refresh_hess_eigen_mo failed: Produced "// &
                    "error when replacing an eigendecomposition for the "//case_name// &
                    " case."
                test_refresh_hess_eigen_mo = .false.
            end if
            if (norm2(mo%mo_channels(1)%occ_eigvals - 2.0_rp * eigvals_before) > tol) &
                then
                write(stderr, *) "test_refresh_hess_eigen_mo failed: "// &
                    "Eigendecomposition not replaced for the "//case_name//" case."
                test_refresh_hess_eigen_mo = .false.
            end if
            deallocate(mo)
        end do

    contains

        function same_irrep_part(matrix, block_irreps) result(coupled)
            !
            ! this function returns the part of a matrix between orbitals of the given
            ! irreps which couples only orbitals of the same irrep
            !
            real(rp), intent(in) :: matrix(:, :)
            integer(ip), intent(in) :: block_irreps(:)
            real(rp) :: coupled(size(matrix, 1), size(matrix, 2))

            coupled = merge(matrix, 0.0_rp, spread(block_irreps, 1, size( &
                matrix, 1)) == spread(block_irreps, 2, size(matrix, 1)))

        end function same_irrep_part

        logical function mixes_irreps(eigvecs, block_irreps)
            !
            ! this function determines if any eigenvector has a component at an orbital
            ! of another irrep than the orbital it sits at
            !
            real(rp), intent(in) :: eigvecs(:, :)
            integer(ip), intent(in) :: block_irreps(:)

            mixes_irreps = any(abs(eigvecs) > tol .and. spread(block_irreps, 1, size( &
                eigvecs, 1)) /= spread(block_irreps, 2, size(eigvecs, 1)))

        end function mixes_irreps

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
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names

        type(mo_type), allocatable, target :: mo
        integer(ip) :: i_case

        ! assume tests pass
        test_get_hess_eigval_pairs_mo = .true.

        ! loop over every occupation case: differences of virtual and occupied
        ! eigenvalues
        do i_case = 1, n_cases
            call setup_minimal_mo_object( &
                mo, case_n_occ(:case_n_particle(i_case), i_case), &
                case_irreps(:, :case_n_particle(i_case), i_case))
            if (norm2(mo%get_hess_eigval_pairs() - &
                      ref_hess_eigval_pairs_mo(mo%mo_channels)) > tol) then
                write(stderr, *) "test_get_hess_eigval_pairs_mo failed: Incorrect "// &
                    "eigenvalue differences for the "//trim(case_names(i_case))// &
                    " case."
                test_get_hess_eigval_pairs_mo = .false.
            end if
            deallocate(mo)
        end do

    end function test_get_hess_eigval_pairs_mo

    logical(c_bool) function test_get_extra_trial_vectors_mo() bind(C)
        !
        ! this function tests the subroutine which returns extra trial vectors along
        ! the most negative eigenvalue pairs of the static part of the Hessian for the
        ! MO basis
        !
        use otr_common_unit_tests, only: level_context
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_occ
        use opentrustregion_unit_tests, only: test_context_type, setup_settings

        type(test_context_type), target :: context
        integer(ip), parameter :: n_extra = 3

        type(mo_type), allocatable, target :: mo
        real(rp), allocatable :: trial_vectors(:, :), expected(:, :), pairs(:), &
                                 unit_vector(:)
        integer(ip) :: error, i, j, k, ivec, min_idx

        ! assume tests pass
        test_get_extra_trial_vectors_mo = .true.

        ! set orbital energies with exactly two negative pairs in the first channel and
        ! none in the second, so that the last requested vector has to vanish
        call setup_minimal_mo_object(mo, n_occ)
        call setup_settings(mo%settings, context)
        mo%mo_channels(1)%occ_eigvals = [-1.0_rp, 2.0_rp]
        mo%mo_channels(1)%virt_eigvals = [0.5_rp, 1.5_rp]
        mo%mo_channels(2)%occ_eigvals = [-3.0_rp]
        mo%mo_channels(2)%virt_eigvals = [0.5_rp, 1.0_rp, 2.0_rp]

        ! construct the expected trial vectors, the rotations along the eigenvalue
        ! pairs in increasing order of the negative ones, rotated out of the
        ! eigenbasis, with vanishing slots once no pair is negative any more
        pairs = ref_hess_eigval_pairs_mo(mo%mo_channels)
        allocate(expected(size(pairs), n_extra), unit_vector(size(pairs)))
        expected = 0.0_rp
        do ivec = 1, n_extra
            min_idx = minloc(pairs, dim=1)
            if (pairs(min_idx) >= 0.0_rp) exit
            unit_vector = 0.0_rp
            unit_vector(min_idx) = 1.0_rp
            expected(:, ivec) = &
                ref_rotate_eigenbasis_mo(unit_vector, mo%mo_channels, n_mo, .false.)
            pairs(min_idx) = huge(1.0_rp)
        end do
        if (norm2(expected(:, n_extra)) > tol .or. norm2(expected(:, 2)) < tol) then
            write(stderr, *) "test_get_extra_trial_vectors_mo failed: Test fixture "// &
                "does not have exactly two negative eigenvalue pairs."
            test_get_extra_trial_vectors_mo = .false.
        end if

        ! call routine and determine if the trial vectors are correct
        allocate(trial_vectors(mo%n_param, n_extra))
        call mo%get_extra_trial_vectors(trial_vectors, error, level_context(mo, 1_ip))
        if (error /= 0) then
            write(stderr, *) "test_get_extra_trial_vectors_mo failed: Produced error."
            test_get_extra_trial_vectors_mo = .false.
        end if
        if (norm2(trial_vectors - expected) > tol) then
            write(stderr, *) "test_get_extra_trial_vectors_mo failed: Incorrect "// &
                "trial vectors."
            test_get_extra_trial_vectors_mo = .false.
        end if

        ! mark the eigendecomposition stale for diagonal Fock matrix blocks with the
        ! same orbital energies, whose eigenvectors are unit vectors since the
        ! diagonals are ascending, and determine if it is refreshed before the trial
        ! vectors, then unit vectors up to their sign, are constructed
        do k = 1, 2
            associate (channel => mo%mo_channels(k))
                channel%fock_oo = diagonal_matrix(channel%occ_eigvals)
                channel%fock_vv = diagonal_matrix(channel%virt_eigvals)
            end associate
        end do
        mo%hess_eigen_stale = .true.
        pairs = [( &
            ((mo%mo_channels(k)%virt_eigvals(j) - mo%mo_channels(k)%occ_eigvals(i), &
              i=1, n_occ(k)), j=1, n_mo - n_occ(k)), k=1, 2)]
        expected = 0.0_rp
        do ivec = 1, 2
            min_idx = minloc(pairs, dim=1)
            expected(min_idx, ivec) = 1.0_rp
            pairs(min_idx) = huge(1.0_rp)
        end do
        call mo%get_extra_trial_vectors(trial_vectors, error, level_context(mo, 1_ip))
        if (error /= 0) then
            write(stderr, *) "test_get_extra_trial_vectors_mo failed: Produced "// &
                "error for a stale eigendecomposition."
            test_get_extra_trial_vectors_mo = .false.
        end if
        if (norm2(abs(trial_vectors) - expected) > tol) then
            write(stderr, *) "test_get_extra_trial_vectors_mo failed: Incorrect "// &
                "trial vectors for a stale eigendecomposition."
            test_get_extra_trial_vectors_mo = .false.
        end if
        deallocate(mo)

    end function test_get_extra_trial_vectors_mo

    logical(c_bool) function test_rotate_mo_coeff() bind(C)
        !
        ! this function tests the subroutine which rotates MO coefficients by an
        ! orbital rotation, reorthonormalizes them and returns the resulting density
        ! matrix
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: rotate_mo_coeff, mo_channel_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_occ_ref => n_occ, &
                                             n_particle_ref => n_particle

        type(test_context_type), target :: context
        real(rp) :: mo_coeff(n_ao, n_mo, n_particle_ref), ao_overlap(n_ao, n_ao), &
                    rot_mo_coeff(n_ao, n_mo, n_particle_ref), &
                    rot_dm_ao(n_ao, n_ao, n_particle_ref), &
                    expected(n_ao, n_mo, n_particle_ref)
        real(rp), allocatable :: kappa(:)
        integer(ip), allocatable :: n_occ(:)
        integer(ip) :: n_particle, error, k, i_case
        character(len=:), allocatable :: case_name
        type(mo_channel_type), allocatable :: mo_channels(:)
        type(orbital_settings_type) :: settings

        ! assume tests pass
        test_rotate_mo_coeff = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! loop over every occupation case
        do i_case = 1, n_cases
            case_name = trim(case_names(i_case))
            n_particle = case_n_particle(i_case)
            n_occ = case_n_occ(:n_particle, i_case)

            ! generate random orthonormal MO coefficients and the occupations and
            ! irreps of their particle channels, whose parameters are the
            ! occupied-virtual pairs of the same irrep
            ao_overlap = generate_random_ao_overlap(n_ao)
            mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle_ref)
            mo_channels = masked_mo_channels(n_occ, n_mo, &
                                             case_irreps(:, :n_particle, i_case))

            ! call routine for a random orbital rotation and determine if the rotated
            ! orbitals and density matrix are correct
            allocate(kappa( &
                sum([(count(mo_channels(k)%param_mask, kind=ip), k=1, n_particle)])))
            call random_number(kappa)
            kappa = 0.2_rp * (kappa - 0.5_rp)
            expected(:, :, :n_particle) = &
                ref_rotate_mo_coeff(kappa, mo_coeff(:, :, :n_particle), mo_channels)
            call rotate_mo_coeff(kappa, mo_coeff(:, :, :n_particle), ao_overlap, &
                                 mo_channels, rot_mo_coeff(:, :, :n_particle), &
                                 rot_dm_ao(:, :, :n_particle), settings, error)
            if (error /= 0) then
                write(stderr, *) "test_rotate_mo_coeff failed: Produced error for "// &
                    "the "//case_name//" case."
                test_rotate_mo_coeff = .false.
            end if
            if (norm2(rot_mo_coeff(:, :, :n_particle) - expected(:, :, :n_particle)) > &
                tol) then
                write(stderr, *) "test_rotate_mo_coeff failed: Incorrect rotated "// &
                    "MO coefficients for the "//case_name//" case."
                test_rotate_mo_coeff = .false.
            end if
            do k = 1, n_particle
                if (norm2(rot_dm_ao(:, :, k) - &
                          matmul(expected(:, :n_occ(k), k), &
                                 transpose(expected(:, :n_occ(k), k)))) > tol) then
                    write(stderr, *) "test_rotate_mo_coeff failed: Incorrect "// &
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
        mo_channels = masked_mo_channels(n_occ_ref(:1), n_mo)
        allocate(kappa(n_occ_ref(1) * (n_mo - n_occ_ref(1))))
        kappa = 0.0_rp
        call rotate_mo_coeff(kappa, mo_coeff(:, :, :1), ao_overlap, mo_channels, &
                             rot_mo_coeff(:, :, :1), rot_dm_ao(:, :, :1), settings, &
                             error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_mo_coeff failed: Produced error for MO "// &
                "coefficients off orthonormality."
            test_rotate_mo_coeff = .false.
        end if
        if (norm2(rot_mo_coeff(:, :, 1) - &
                  ref_lowdin_orthonormalize(mo_coeff(:, :, 1))) > tol) then
            write(stderr, *) "test_rotate_mo_coeff failed: MO coefficients not "// &
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
        use otr_mo, only: hess_x_static_mo, mo_type
        use otr_mo_test_reference, only: n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names

        type(mo_type), allocatable, target :: mo
        real(rp), allocatable :: x(:)
        integer(ip) :: i_case

        ! assume tests pass
        test_hess_x_static_mo = .true.

        ! loop over every occupation case
        do i_case = 1, n_cases
            call setup_minimal_mo_object( &
                mo, case_n_occ(:case_n_particle(i_case), i_case), &
                case_irreps(:, :case_n_particle(i_case), i_case))
            allocate(x(mo%n_param))
            call random_number(x)
            if (norm2(hess_x_static_mo(x, mo%mo_channels) - &
                      ref_hess_x_static_mo(x, mo%mo_channels)) > tol) then
                write(stderr, *) "test_hess_x_static_mo failed: Incorrect static "// &
                    "part for the "//trim(case_names(i_case))//" case."
                test_hess_x_static_mo = .false.
            end if
            deallocate(x, mo)
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
            write(stderr, *) "test_mo_transform failed: Incorrect dimensions."
            test_mo_transform = .false.
        else if (norm2(transformed - ref_mo_transform(coeff, matrix)) > tol) then
            write(stderr, *) "test_mo_transform failed: Incorrect transformed matrix."
            test_mo_transform = .false.
        end if

    end function test_mo_transform

    logical(c_bool) function test_count_mo_params() bind(C)
        !
        ! this function tests the function which returns the number of parameters of
        ! orbitals parameterized in the MO basis, with and without irreps of the MOs
        !
        use otr_mo, only: count_mo_params
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names

        integer(ip) :: n_particle, i_case

        ! assume tests pass
        test_count_mo_params = .true.

        ! loop over every occupation case and determine if all occupied-virtual pairs
        ! are counted without irreps and only those of the same irrep with irreps
        do i_case = 1, n_cases
            n_particle = case_n_particle(i_case)
            associate (n_occ => case_n_occ(:n_particle, i_case), &
                       irreps => case_irreps(:, :n_particle, i_case))
                if (count_mo_params(n_occ, n_mo) /= sum(n_occ * (n_mo - n_occ))) then
                    write(stderr, *) "test_count_mo_params failed: Incorrect "// &
                        "number of parameters without irreps for the "// &
                        trim(case_names(i_case))//" case."
                    test_count_mo_params = .false.
                end if
                if (count_mo_params(n_occ, n_mo, irreps) /= &
                    ref_count_mo_params(n_occ, irreps)) then
                    write(stderr, *) "test_count_mo_params failed: Incorrect "// &
                        "number of parameters with irreps for the "// &
                        trim(case_names(i_case))//" case."
                    test_count_mo_params = .false.
                end if
            end associate
        end do

    end function test_count_mo_params

    logical(c_bool) function test_mo_param_rows() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in the parameter vector of the MO basis
        !
        use otr_mo, only: mo_param_rows, mo_channel_type
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names

        integer(ip) :: n_particle, n_rows, i_case, k
        integer(ip), allocatable :: rows(:, :), expected(:, :)
        type(mo_channel_type), allocatable :: mo_channels(:)

        ! assume tests pass
        test_mo_param_rows = .true.

        ! loop over every occupation case and determine if the channels occupy the rows
        ! of their occupied-virtual pairs of the same irrep one after another
        do i_case = 1, n_cases
            n_particle = case_n_particle(i_case)
            mo_channels = masked_mo_channels(case_n_occ(:n_particle, i_case), n_mo, &
                                             case_irreps(:, :n_particle, i_case))
            allocate(expected(2, n_particle))
            n_rows = 0
            do k = 1, n_particle
                expected(1, k) = n_rows + 1
                n_rows = n_rows + ref_count_mo_params(case_n_occ(k:k, i_case), &
                                                      case_irreps(:, k:k, i_case))
                expected(2, k) = n_rows
            end do
            rows = mo_param_rows(mo_channels)
            if (any(rows /= expected)) then
                write(stderr, *) "test_mo_param_rows failed: Incorrect rows for "// &
                    "the "//trim(case_names(i_case))//" case."
                test_mo_param_rows = .false.
            end if
            deallocate(expected)
        end do

    end function test_mo_param_rows

    logical(c_bool) function test_pack_ov() bind(C)
        !
        ! this function tests the function which packs the parameters of a particle
        ! channel from its occupied-virtual block
        !
        use otr_mo, only: pack_ov, mo_channel_type
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names

        real(rp) :: matrix(n_mo, n_mo, 1)
        integer(ip) :: n_particle, n_occ, i_case, k
        type(mo_channel_type), allocatable :: mo_channels(:)

        ! assume tests pass
        test_pack_ov = .true.

        ! loop over every occupation case and particle channel and determine if a
        ! random occupied-virtual block is packed
        do i_case = 1, n_cases
            n_particle = case_n_particle(i_case)
            mo_channels = masked_mo_channels(case_n_occ(:n_particle, i_case), n_mo, &
                                             case_irreps(:, :n_particle, i_case))
            do k = 1, n_particle
                n_occ = case_n_occ(k, i_case)
                call random_number(matrix)
                if (norm2(pack_ov(matrix(:n_occ, n_occ + 1:, 1), mo_channels(k)) - &
                          ref_pack_ov(matrix, mo_channels(k:k))) > tol) then
                    write(stderr, *) "test_pack_ov failed: Incorrect parameters "// &
                        "for the "//trim(case_names(i_case))//" case."
                    test_pack_ov = .false.
                end if
            end do
        end do

    end function test_pack_ov

    logical(c_bool) function test_unpack_ov() bind(C)
        !
        ! this function tests the function which unpacks the parameters of a particle
        ! channel into its occupied-virtual block
        !
        use otr_mo, only: unpack_ov, mo_channel_type
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names

        real(rp) :: expected(n_mo, n_mo, 1)
        real(rp), allocatable :: packed(:)
        integer(ip) :: n_particle, n_occ, i_case, k
        type(mo_channel_type), allocatable :: mo_channels(:)

        ! assume tests pass
        test_unpack_ov = .true.

        ! loop over every occupation case and particle channel and determine if random
        ! parameters are unpacked into the occupied-virtual block, which vanishes
        ! between orbitals of different irreps
        do i_case = 1, n_cases
            n_particle = case_n_particle(i_case)
            mo_channels = masked_mo_channels(case_n_occ(:n_particle, i_case), n_mo, &
                                             case_irreps(:, :n_particle, i_case))
            do k = 1, n_particle
                n_occ = case_n_occ(k, i_case)
                allocate(packed(count(mo_channels(k)%param_mask)))
                call random_number(packed)
                expected = ref_unpack_ov(packed, mo_channels(k:k))
                if (norm2(unpack_ov(packed, mo_channels(k)) - &
                          expected(:n_occ, n_occ + 1:, 1)) > tol) then
                    write(stderr, *) "test_unpack_ov failed: Incorrect "// &
                        "occupied-virtual block for the "//trim(case_names(i_case))// &
                        " case."
                    test_unpack_ov = .false.
                end if
                deallocate(packed)
            end do
        end do

    end function test_unpack_ov

    logical(c_bool) function test_diagonalize_per_irrep() bind(C)
        !
        ! this function tests the subroutine which diagonalizes a symmetric matrix
        ! between orbitals of given irreps separately for every irrep
        !
        use otr_common, only: orbital_settings_type
        use otr_mo, only: diagonalize_per_irrep
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix

        type(test_context_type), target :: context
        integer(ip), parameter :: n = 4, n_matrices = 3
        integer(ip), parameter :: matrix_irreps(n, n_matrices) = &
            reshape([0_ip, 0_ip, 0_ip, 0_ip, &
                     0_ip, 1_ip, 0_ip, 1_ip, &
                     0_ip, 1_ip, 0_ip, 1_ip], [n, n_matrices])
        character(len=33), parameter :: matrix_names(n_matrices) = &
            [character(len=33) :: "a single irrep", &
             "two irreps with degenerate blocks", "two irreps coupled to each other"]

        real(rp) :: matrix(n, n), same_irrep_part(n, n), eigvals(n), eigvecs(n, n), &
                    block(2, 2)
        integer(ip) :: error, i_matrix
        logical :: same_irrep(n, n)
        character(len=:), allocatable :: name
        type(orbital_settings_type) :: settings

        ! assume tests pass
        test_diagonalize_per_irrep = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! loop over a random matrix of a single irrep, one of two irreps whose blocks
        ! are the same so that every eigenvalue is degenerate between the irreps, and a
        ! random matrix coupling two irreps, and determine if the part of the matrix
        ! between orbitals of the same irrep is diagonalized by orthonormal
        ! eigenvectors at the orbitals of their irrep, ignoring any coupling between
        ! the irreps
        do i_matrix = 1, n_matrices
            name = trim(matrix_names(i_matrix))
            if (i_matrix == 2) then
                block = generate_random_symm_matrix(2_ip)
                matrix = 0.0_rp
                matrix([1, 3], [1, 3]) = block
                matrix([2, 4], [2, 4]) = block
            else
                matrix = generate_random_symm_matrix(n)
            end if
            same_irrep = spread(matrix_irreps(:, i_matrix), 1, n) == &
                         spread(matrix_irreps(:, i_matrix), 2, n)
            same_irrep_part = merge(matrix, 0.0_rp, same_irrep)
            call diagonalize_per_irrep(matrix, matrix_irreps(:, i_matrix), eigvals, &
                                       eigvecs, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_diagonalize_per_irrep failed: Produced "// &
                    "error for "//name//"."
                test_diagonalize_per_irrep = .false.
            end if
            if (norm2(matmul(transpose(eigvecs), matmul(same_irrep_part, eigvecs)) - &
                      diagonal_matrix(eigvals)) > tol) then
                write(stderr, *) "test_diagonalize_per_irrep failed: Matrix not "// &
                    "diagonalized for "//name//"."
                test_diagonalize_per_irrep = .false.
            end if
            if (norm2(matmul(transpose(eigvecs), eigvecs) - identity_matrix(n)) > tol) &
                then
                write(stderr, *) "test_diagonalize_per_irrep failed: Eigenvectors "// &
                    "not orthonormal for "//name//"."
                test_diagonalize_per_irrep = .false.
            end if
            if (any(abs(eigvecs) > tol .and. .not. same_irrep)) then
                write(stderr, *) "test_diagonalize_per_irrep failed: Eigenvectors "// &
                    "mix irreps for "//name//"."
                test_diagonalize_per_irrep = .false.
            end if
        end do

    end function test_diagonalize_per_irrep

end module otr_mo_unit_tests
