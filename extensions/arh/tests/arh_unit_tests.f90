! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_unit_tests

    use opentrustregion, only: rp, ip, kw_len, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! multipliers of the density matrix returned by the mock density matrix updating
    ! functions, which differ between the first and any subsequent call, following the
    ! same convention as the shared multiplier for the Fock matrix
    real(rp), parameter :: mock_v_same_spin_factor(2) = [3.0_rp, 7.0_rp], &
                           mock_v_opposite_spin_factor(2) = [4.0_rp, 9.0_rp], &
                           mock_v_nonlinear_factor(2) = [6.0_rp, 11.0_rp]

    ! dispatch the mock potential by density matrix rank
    interface mock_potential
        module procedure mock_potential_cs, mock_potential_os
    end interface mock_potential

contains

    function mock_potential_cs(factor, dm) result(v)
        !
        ! this function evaluates a mock potential for a given density matrix as a
        ! multiple of it plus the anticommutator with a fixed symmetric matrix; the 
        ! anticommutator makes sure the potential does not commute with the density 
        ! matrix and is therefore not annihilated by the occupied-virtual projection, 
        ! which makes testing of projected quantities possible
        !
        real(rp), intent(in) :: factor, dm(:, :)
        real(rp) :: v(size(dm, 1), size(dm, 1))

        integer(ip) :: i, j
        real(rp) :: coupling(size(dm, 1), size(dm, 1))

        do j = 1, size(dm, 1)
            do i = 1, size(dm, 1)
                coupling(i, j) = 1.0_rp / real(i + j, kind=rp)
            end do
        end do
        v = factor * dm + 0.25_rp * (matmul(coupling, dm) + matmul(dm, coupling))

    end function mock_potential_cs

    function mock_potential_os(factor, dm) result(v)
        !
        ! this function applies mock_potential_cs to every spin channel of an
        ! open-shell density matrix
        !
        real(rp), intent(in) :: factor, dm(:, :, :)
        real(rp) :: v(size(dm, 1), size(dm, 1), size(dm, 3))

        integer(ip) :: i

        do i = 1, size(dm, 3)
            v(:, :, i) = mock_potential(factor, dm(:, :, i))
        end do

    end function mock_potential_os

    subroutine mock_update_dm_cs(dm, energy, fock, v_nonlinear, error)
        !
        ! this subroutine is a mock density matrix updating function with a separate
        ! non-linear potential contribution for the closed-shell case, which returns
        ! multiples of the density matrix that change between calls so that
        ! non-vanishing differences are produced
        !
        use otr_oao_unit_tests, only: n_mock_calls, mock_factor, mock_fock_factor

        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target, contiguous :: fock(:, :), v_nonlinear(:, :)
        integer(ip), intent(out) :: error

        n_mock_calls = n_mock_calls + 1

        error = 0
        energy = sum(dm)
        fock = mock_potential(mock_factor(mock_fock_factor), dm)
        v_nonlinear = mock_potential(mock_factor(mock_v_nonlinear_factor), dm)

    end subroutine mock_update_dm_cs

    subroutine mock_update_dm_os(dm, energy, fock, v_same_spin, v_opposite_spin, &
                                 v_nonlinear, error)
        !
        ! this subroutine is a mock density matrix updating function with spin-resolved
        ! and non-linear potential contributions for the open-shell case, which returns
        ! multiples of the density matrix that change between calls so that
        ! non-vanishing differences are produced
        !
        use otr_oao_unit_tests, only: n_mock_calls, mock_factor, mock_fock_factor

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target :: fock(:, :, :), v_same_spin(:, :, :), &
                                         v_opposite_spin(:, :, :), v_nonlinear(:, :, :)
        integer(ip), intent(out) :: error

        integer(ip) :: i

        n_mock_calls = n_mock_calls + 1

        error = 0
        energy = sum(dm)
        do i = 1, size(dm, 3)
            fock(:, :, i) = mock_potential(mock_factor(mock_fock_factor), dm(:, :, i))
            v_same_spin(:, :, i) = &
                mock_potential(mock_factor(mock_v_same_spin_factor), dm(:, :, i))
            v_opposite_spin(:, :, i) = &
                mock_potential(mock_factor(mock_v_opposite_spin_factor), dm(:, :, i))
            v_nonlinear(:, :, i) = &
                mock_potential(mock_factor(mock_v_nonlinear_factor), dm(:, :, i))
        end do

    end subroutine mock_update_dm_os

        function embed_channel(v, channel, n_ao, n_particle) result(embedded)
        !
        ! this function embeds a single (n_ao, n_ao) matrix into one channel of a
        ! full (n_ao, n_ao, n_particle) array, zeroing the other channel(s)
        !
        real(rp), intent(in) :: v(:, :)
        integer(ip), intent(in) :: channel, n_ao, n_particle
        real(rp) :: embedded(n_ao, n_ao, n_particle)

        embedded = 0.0_rp
        embedded(:, :, channel) = v

    end function embed_channel

    function generate_random_dm_diff(dm, n) result(dm_diff)
        !
        ! this function generates a valid density matrix difference as the commutator
        ! of a random antisymmetric matrix with a density matrix
        !
        real(rp), intent(in) :: dm(:, :)
        integer(ip), intent(in) :: n
        real(rp) :: dm_diff(n, n)

        real(rp) :: x(n, n)

        call random_number(x)
        x = x - transpose(x)
        dm_diff = matmul(dm, x)
        dm_diff = dm_diff + transpose(dm_diff)

    end function generate_random_dm_diff

    function generate_random_symm_hessian(n) result(hess)
        !
        ! this function generates a random mock Hessian operator with the full
        ! permutational symmetry of the two-electron integrals (the array is only ever
        ! read for a first index pair whose first index is smaller than or equal to its
        ! second index)
        !
        integer(ip), intent(in) :: n
        real(rp) :: hess(n, n, n, n)

        integer(ip) :: i, j, k, l
        real(rp) :: val

        do l = 1, n
            do k = 1, l
                do j = 1, l
                    do i = 1, merge(k, j, j == l)
                        call random_number(val)
                        hess(i, j, k, l) = val
                        hess(i, j, l, k) = val
                        hess(k, l, i, j) = val
                        hess(k, l, j, i) = val
                    end do
                end do
            end do
        end do

    end function generate_random_symm_hessian

    function contract_symm_hessian(hess, matrix) result(contracted)
        !
        ! this function contracts a mock Hessian operator with a symmetric matrix to
        ! produce the corresponding symmetric potential matrix
        !
        real(rp), intent(in) :: hess(:, :, :, :), matrix(:, :)
        real(rp) :: contracted(size(matrix, 1), size(matrix, 2))

        integer(ip) :: i, j

        do j = 1, size(matrix, 2)
            do i = 1, j
                contracted(i, j) = sum(hess(i, j, :, :) * matrix)
                contracted(j, i) = contracted(i, j)
            end do
        end do

    end function contract_symm_hessian

    subroutine generate_fock_partition(n_ao, n_particle, n_electrons, dm_oao, fock_oo, &
                                       fock_vv)
        !
        ! this subroutine generates a random density matrix per particle/spin channel
        ! and the corresponding occupied-occupied and virtual-virtual Fock matrix
        ! blocks
        !
        use otr_oao_unit_tests, only: generate_random_density_matrix, &
                                      generate_random_symm_matrix, identity_matrix

        integer(ip), intent(in) :: n_ao, n_particle, n_electrons
        real(rp), intent(out) :: dm_oao(:, :, :), fock_oo(:, :, :), fock_vv(:, :, :)

        real(rp) :: full_fock(n_ao, n_ao), complement(n_ao, n_ao)
        integer(ip) :: j

        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)
            full_fock = generate_random_symm_matrix(n_ao)
            complement = identity_matrix(n_ao) - dm_oao(:, :, j)
            fock_oo(:, :, j) = matmul(dm_oao(:, :, j), matmul(full_fock, &
                                                              dm_oao(:, :, j)))
            fock_vv(:, :, j) = matmul(complement, matmul(full_fock, complement))
        end do

    end subroutine generate_fock_partition

    function ref_cache_dirs(v_diff, dm_oao, n_param) result(dirs)
        !
        ! this function independently projects and packs every column of a
        ! history-difference array into the non-redundant subspace
        !
        use otr_oao_unit_tests, only: ref_project_asymm, ref_pack_asymm

        real(rp), intent(in) :: v_diff(:, :, :, :), dm_oao(:, :, :)
        integer(ip), intent(in) :: n_param
        real(rp) :: dirs(n_param, size(v_diff, 4))

        integer(ip) :: k

        do k = 1, size(v_diff, 4, kind=ip)
            dirs(:, k) = ref_pack_asymm(ref_project_asymm(v_diff(:, :, :, k), dm_oao), &
                                        n_param)
        end do

    end function ref_cache_dirs

    subroutine setup_arh_and_oao_objects(arh_type, dm_oao, fock_oo, fock_vv, &
                                         n_ao_target, n_particle_target, n_param_target)
        !
        ! this subroutine sets up the module-global OAO and ARH objects for the
        ! preconditioner tests; the allocatable caches the individual preconditioner 
        ! branches need are filled in by the caller
        !
        use otr_arh, only: arh_object
        use otr_oao, only: oao_object
        use opentrustregion_unit_tests, only: setup_settings

        character(*), intent(in) :: arh_type
        real(rp), target :: dm_oao(:, :, :), fock_oo(:, :, :), fock_vv(:, :, :)
        integer(ip), target :: n_ao_target, n_particle_target, n_param_target

        ! set up the OAO object so that the static-Hessian eigendecomposition used by
        ! the preconditioner can be refreshed
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao_target
        oao_object%n_particle = n_particle_target
        oao_object%n_param = n_param_target
        oao_object%fock_oo = fock_oo
        oao_object%fock_vv = fock_vv
        oao_object%hess_eigen_stale = .true.

        ! set up the ARH object
        allocate(arh_object)
        call setup_settings(arh_object%settings)
        arh_object%settings%arh_type = arh_type
        arh_object%n_ao => n_ao_target
        arh_object%n_particle => n_particle_target
        arh_object%n_param => n_param_target
        arh_object%dm_oao => dm_oao
        arh_object%fock_oo => fock_oo
        arh_object%fock_vv => fock_vv

    end subroutine setup_arh_and_oao_objects

    function generate_nonredundant_vector(n_param, n_particle, n_ao, dm_oao) &
        result(vector)
        !
        ! this function generates a random vector confined to the non-redundant
        ! subspace
        !
        use otr_oao_unit_tests, only: ref_unpack_asymm, ref_project_asymm, &
                                      ref_pack_asymm

        integer(ip), intent(in) :: n_param, n_particle, n_ao
        real(rp), intent(in) :: dm_oao(:, :, :)
        real(rp) :: vector(n_param)

        call random_number(vector)
        vector = ref_pack_asymm(ref_project_asymm(ref_unpack_asymm( &
            vector, n_particle, n_ao), dm_oao), n_param)

    end function generate_nonredundant_vector

    subroutine check_inv_hess_x_arh(hess, vector, mu, test_name, passed)
        !
        ! this subroutine checks the exact-inverse of the ARH Woodbury preconditioner 
        ! against a reference Hessian by applying the reference (B - mu*I) to the 
        ! preconditioned vector which has to recover the vector, both with and without 
        ! the level shift
        !
        use otr_arh, only: inv_hess_x_arh

        real(rp), intent(in) :: hess(:, :), vector(:), mu
        character(*), intent(in) :: test_name
        logical, intent(inout) :: passed

        real(rp), allocatable :: actual(:)
        integer(ip) :: error

        allocate(actual(size(vector)))

        ! shifted entry point: applying the reference (B - mu*I) to the preconditioned 
        ! vector has to recover the vector
        call inv_hess_x_arh(vector, actual, error, mu)
        if (error /= 0) then
            write (stderr, *) test_name//" failed: Produced error for the shifted "// &
                "inverse."
            passed = .false.
        end if
        if (norm2(matmul(hess, actual) - mu * actual - vector) > tol * &
            (1.0_rp + norm2(hess) * norm2(actual))) then
            write (stderr, *) test_name//" failed: Shifted inverse is not inverse "// &
                "of shifted reference Hessian."
            passed = .false.
        end if

        ! the unshifted entry point has to reproduce the same round trip against the
        ! reference Hessian itself, without a level shift
        call inv_hess_x_arh(vector, actual, error)
        if (error /= 0) then
            write (stderr, *) test_name//" failed: Produced error for the "// &
                "unshifted inverse."
            passed = .false.
        end if
        if (norm2(matmul(hess, actual) - vector) > tol * &
            (1.0_rp + norm2(hess) * norm2(actual))) then
            write (stderr, *) test_name//" failed: Inverse is not inverse "// &
                "of reference Hessian."
            passed = .false.
        end if

    end subroutine check_inv_hess_x_arh

    function ref_response_cs(arh_type, delta_dm, dm_diff, fock_diff, v_linear_diff, &
                             v_nonlinear_diff, metric_inv, a_sym, a_inv, a_inv_comb, &
                             n_ao, n_diff) result(response)
        !
        ! this function independently reimplements the closed-shell ARH response
        ! contribution for every ARH types
        !
        character(*), intent(in) :: arh_type
        real(rp), intent(in) :: delta_dm(:, :, :), dm_diff(:, :, :, :), &
                                fock_diff(:, :, :, :), v_linear_diff(:, :, :, :), &
                                v_nonlinear_diff(:, :, :, :), metric_inv(:, :), &
                                a_sym(:, :), a_inv(:, :), a_inv_comb(:, :)
        integer(ip), intent(in) :: n_ao, n_diff
        real(rp) :: response(n_ao, n_ao, 1)

        integer(ip) :: i
        real(rp) :: s_proj(n_diff), alpha_s(n_diff), y_proj(n_diff), alpha_y(n_diff), &
                   sy(n_diff), alpha_sy(n_diff), lin_proj(n_diff), alpha_lin(n_diff), &
                   nl_proj(n_diff), alpha_nl(n_diff)

        response = 0.0_rp

        ! MS-SR1 keeps the linear and non-linear potential differences on their own, 
        ! independent systems
        if (arh_type == "ms_sr1") then
            do i = 1, n_diff
                lin_proj(i) = sum(v_linear_diff(:, :, 1, i) * delta_dm(:, :, 1))
            end do
            alpha_lin = matmul(a_inv, lin_proj)
            do i = 1, n_diff
                response(:, :, 1) = response(:, :, 1) + alpha_lin(i) * &
                                    v_linear_diff(:, :, 1, i)
            end do
            do i = 1, n_diff
                nl_proj(i) = sum(v_nonlinear_diff(:, :, 1, i) * delta_dm(:, :, 1))
            end do
            alpha_nl = matmul(a_inv_comb, nl_proj)
            do i = 1, n_diff
                response(:, :, 1) = response(:, :, 1) + alpha_nl(i) * &
                                    v_nonlinear_diff(:, :, 1, i)
            end do
            return
        end if

        ! every other type starts by contracting the density-matrix-difference
        ! history
        do i = 1, n_diff
            s_proj(i) = sum(dm_diff(:, :, 1, i) * delta_dm(:, :, 1))
        end do
        alpha_s = matmul(metric_inv, s_proj)

        ! the shared output direction of standard ARH
        do i = 1, n_diff
            response(:, :, 1) = response(:, :, 1) + alpha_s(i) * fock_diff(:, :, 1, i)
        end do

        ! standard ARH stops here
        if (arh_type == "arh") return

        ! the transposed contraction shared by symmetrized ARH and MS-PSB
        do i = 1, n_diff
            y_proj(i) = sum(fock_diff(:, :, 1, i) * delta_dm(:, :, 1))
        end do
        alpha_y = matmul(metric_inv, y_proj)

        if (arh_type == "symm_arh") then
            ! symmetrized ARH averages the two contractions
            response = 0.5_rp * response
            do i = 1, n_diff
                response(:, :, 1) = response(:, :, 1) + &
                                    0.5_rp * alpha_y(i) * dm_diff(:, :, 1, i)
            end do
            return
        end if

        ! the subspace-projected term shared by the MS-SP and MS-PSB methods
        sy = matmul(a_sym, alpha_s)
        alpha_sy = matmul(metric_inv, sy)

        if (arh_type == "ms_sp") then
            ! the MS-SP method keeps only that term
            response = 0.0_rp
            do i = 1, n_diff
                response(:, :, 1) = response(:, :, 1) + alpha_sy(i) * &
                                    dm_diff(:, :, 1, i)
            end do
        else
            ! MS-PSB adds the transposed contraction and subtracts the
            ! subspace-projected term
            do i = 1, n_diff
                response(:, :, 1) = response(:, :, 1) + &
                                    (alpha_y(i) - alpha_sy(i)) * dm_diff(:, :, 1, i)
            end do
        end if

    end function ref_response_cs

    function ref_response_os(arh_type, delta_dm, dm_diff, v_same_eff, v_opp, v_nl, &
                             metric_inv, a_block, a_inv_lin, a_inv_comb, n_ao, n_diff) &
        result(response)
        !
        ! this function independently reimplements the open-shell ARH response
        ! contribution for every ARH type
        !
        character(*), intent(in) :: arh_type
        real(rp), intent(in) :: delta_dm(:, :, :), dm_diff(:, :, :, :), &
                                v_same_eff(:, :, :, :), v_opp(:, :, :, :), &
                                v_nl(:, :, :, :), metric_inv(:, :), a_block(:, :), &
                                a_inv_lin(:, :), a_inv_comb(:, :)
        integer(ip), intent(in) :: n_ao, n_diff
        real(rp) :: response(n_ao, n_ao, 2)

        integer(ip) :: i, j
        real(rp) :: s_proj(2*n_diff), alpha_s(2*n_diff), y_proj(2*n_diff), &
                    alpha_y(2*n_diff), sy(2*n_diff), alpha_sy(2*n_diff), &
                    nl_proj(n_diff), alpha_nl(n_diff)

        response = 0.0_rp

        ! MS-SR1 keeps the linear part on a joint spin-separated system and the 
        ! non-linear part on a dedicated spin-combined one
        if (arh_type == "ms_sr1") then
            do i = 1, n_diff
                y_proj(i) = sum(v_same_eff(:, :, 1, i) * delta_dm(:, :, 1)) + &
                            sum(v_opp(:, :, 2, i) * delta_dm(:, :, 2))
                y_proj(n_diff + i) = sum(v_opp(:, :, 1, i) * delta_dm(:, :, 1)) + &
                                     sum(v_same_eff(:, :, 2, i) * delta_dm(:, :, 2))
            end do
            alpha_y = matmul(a_inv_lin, y_proj)
            do i = 1, n_diff
                response(:, :, 1) = response(:, :, 1) + &
                                    alpha_y(i) * v_same_eff(:, :, 1, i) + &
                                    alpha_y(n_diff + i) * v_opp(:, :, 1, i)
                response(:, :, 2) = response(:, :, 2) + &
                                    alpha_y(i) * v_opp(:, :, 2, i) + &
                                    alpha_y(n_diff + i) * v_same_eff(:, :, 2, i)
            end do

            ! non-linear part, contracted jointly over both spin channels
            do i = 1, n_diff
                nl_proj(i) = sum(v_nl(:, :, :, i) * delta_dm)
            end do
            alpha_nl = matmul(a_inv_comb, nl_proj)
            do i = 1, n_diff
                do j = 1, 2
                    response(:, :, j) = response(:, :, j) + alpha_nl(i) * &
                                        v_nl(:, :, j, i)
                end do
            end do
            return
        end if

        ! ARH and symmetric variants contract the density matrix history per spin
        ! channel, with the coefficients of one channel driving that channel's
        ! same-spin direction and the other channel's opposite-spin direction
        do j = 1, 2
            do i = 1, n_diff
                s_proj((j - 1)*n_diff + i) = sum(dm_diff(:, :, j, i) * &
                                                 delta_dm(:, :, j))
            end do
        end do
        alpha_s = matmul(metric_inv, s_proj)

        ! the shared output direction of every type below
        do j = 1, 2
            do i = 1, n_diff
                response(:, :, j) = response(:, :, j) + &
                                    alpha_s((j - 1)*n_diff + i) * &
                                    v_same_eff(:, :, j, i) + &
                                    alpha_s((2 - j)*n_diff + i) * v_opp(:, :, j, i)
            end do
        end do

        ! standard ARH stops here
        if (arh_type == "arh") return

        ! the transposed contraction shared by symmetrized ARH and MS-PSB
        do j = 1, 2
            do i = 1, n_diff
                y_proj((j - 1)*n_diff + i) = &
                    sum(v_same_eff(:, :, j, i) * delta_dm(:, :, j)) + &
                    sum(v_opp(:, :, 3 - j, i) * delta_dm(:, :, 3 - j))
            end do
        end do
        alpha_y = matmul(metric_inv, y_proj)

        if (arh_type == "symm_arh") then
            ! symmetrized ARH averages the two contractions
            response = 0.5_rp * response
            do j = 1, 2
                do i = 1, n_diff
                    response(:, :, j) = response(:, :, j) + 0.5_rp * &
                                        alpha_y((j - 1)*n_diff + i) * &
                                        dm_diff(:, :, j, i)
                end do
            end do
            return
        end if

        ! the subspace-projected term shared by MS-SP and MS-PSB
        sy = matmul(a_block, alpha_s)
        alpha_sy = matmul(metric_inv, sy)

        if (arh_type == "ms_sp") then
            ! the MS-SP method keeps only that term
            response = 0.0_rp
            do j = 1, 2
                do i = 1, n_diff
                    response(:, :, j) = response(:, :, j) + &
                                        alpha_sy((j - 1)*n_diff + i) * &
                                        dm_diff(:, :, j, i)
                end do
            end do
        else
            ! MS-PSB adds the transposed contraction and subtracts the
            ! subspace-projected term
            do j = 1, 2
                do i = 1, n_diff
                    response(:, :, j) = response(:, :, j) + &
                                        (alpha_y((j - 1)*n_diff + i) - &
                                         alpha_sy((j - 1)*n_diff + i)) * &
                                         dm_diff(:, :, j, i)
                end do
            end do
        end if

    end function ref_response_os

    function check_inv_hess_x_arh_cs(arh_type, test_name) result(passed)
        !
        ! this function performs the exact-inverse check of the closed-shell ARH 
        ! Woodbury preconditioner for a single ARH type
        !
        use otr_arh, only: arh_object
        use otr_oao, only: oao_object
        use otr_oao_test_reference, only: n_ao
        use otr_oao_unit_tests, only: ref_unpack_asymm, ref_project_symm, ref_hess_x, &
                                      generate_random_symm_matrix

        character(*), intent(in) :: arh_type, test_name
        logical :: passed

        integer(ip), parameter :: n_diff = 2, n_electrons = 2
        real(rp), parameter :: mu = 0.1_rp

        integer(ip), target :: n_ao_target, n_particle_target, n_param_target
        real(rp), target :: dm_oao(n_ao, n_ao, 1), fock_oo(n_ao, n_ao, 1), &
                            fock_vv(n_ao, n_ao, 1)
        real(rp) :: dm_diff(n_ao, n_ao, 1, n_diff), fock_diff(n_ao, n_ao, 1, n_diff), &
                    v_linear_diff(n_ao, n_ao, 1, n_diff), &
                    v_nonlinear_diff(n_ao, n_ao, 1, n_diff), &
                    metric_inv(n_diff, n_diff), a_sym(n_diff, n_diff), &
                    a_inv(n_diff, n_diff), a_inv_comb(n_diff, n_diff)
        integer(ip) :: n_param, i, k
        real(rp), allocatable :: x_full(:, :, :), delta_dm(:, :, :), &
                                 response(:, :, :), hess(:, :), e_i(:), vector(:)

        ! assume test passes
        passed = .true.

        ! get number of parameters
        n_param = n_ao * (n_ao - 1) / 2

        ! generate a random density matrix, the corresponding occupied-occupied and
        ! virtual-virtual Fock matrix blocks, and history differences
        call generate_fock_partition(n_ao, 1_ip, n_electrons, dm_oao, fock_oo, &
                                     fock_vv)
        do k = 1, n_diff
            dm_diff(:, :, 1, k) = generate_random_symm_matrix(n_ao)
            fock_diff(:, :, 1, k) = generate_random_symm_matrix(n_ao)
            v_linear_diff(:, :, 1, k) = generate_random_symm_matrix(n_ao)
            v_nonlinear_diff(:, :, 1, k) = generate_random_symm_matrix(n_ao)
        end do

        ! randomly generated (symmetric) metric pseudoinverse, A_sym, and multisecant 
        ! SR1 pseudoinverse cores
        metric_inv = generate_random_symm_matrix(n_diff)
        a_sym = generate_random_symm_matrix(n_diff)
        a_inv = generate_random_symm_matrix(n_diff)
        a_inv_comb = generate_random_symm_matrix(n_diff)

        ! set up the ARH and OAO objects
        n_ao_target = n_ao
        n_particle_target = 1
        n_param_target = n_param
        call setup_arh_and_oao_objects(arh_type, dm_oao, fock_oo, fock_vv, &
                                       n_ao_target, n_particle_target, n_param_target)

        ! populate the cached history projections
        arh_object%dm_dirs = ref_cache_dirs(dm_diff, dm_oao, n_param)
        arh_object%potential_dirs = ref_cache_dirs(fock_diff, dm_oao, n_param)
        arh_object%linear_potential_dirs = ref_cache_dirs(v_linear_diff, dm_oao, &
                                                          n_param)
        arh_object%nonlinear_potential_dirs = ref_cache_dirs(v_nonlinear_diff, dm_oao, &
                                                             n_param)

        ! assemble the low-rank factors
        select case (arh_type)
        case ("ms_sr1")
            allocate(arh_object%expansion_dirs(n_param, 2 * n_diff), &
                     arh_object%coupling_matrix(2 * n_diff, 2 * n_diff))
            arh_object%expansion_dirs(:, :n_diff) = arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, n_diff + 1:) = &
                arh_object%nonlinear_potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_diff, :n_diff) = 8.0_rp * a_inv
            arh_object%coupling_matrix(n_diff + 1:, n_diff + 1:) = 8.0_rp * a_inv_comb
        case ("ms_sp")
            arh_object%expansion_dirs = arh_object%dm_dirs
            arh_object%coupling_matrix = 8.0_rp * matmul(metric_inv, &
                                                         matmul(a_sym, metric_inv))
        case ("symm_arh")
            allocate(arh_object%expansion_dirs(n_param, 2 * n_diff), &
                     arh_object%coupling_matrix(2 * n_diff, 2 * n_diff))
            arh_object%expansion_dirs(:, :n_diff) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_diff + 1:) = arh_object%potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_diff, n_diff + 1:) = 4.0_rp * metric_inv
            arh_object%coupling_matrix(n_diff + 1:, :n_diff) = 4.0_rp * metric_inv
        case ("ms_psb")
            allocate(arh_object%expansion_dirs(n_param, 2 * n_diff), &
                     arh_object%coupling_matrix(2 * n_diff, 2 * n_diff))
            arh_object%expansion_dirs(:, :n_diff) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_diff + 1:) = arh_object%potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_diff, :n_diff) = &
                -8.0_rp * matmul(metric_inv, matmul(a_sym, metric_inv))
            arh_object%coupling_matrix(:n_diff, n_diff + 1:) = 8.0_rp * metric_inv
            arh_object%coupling_matrix(n_diff + 1:, :n_diff) = 8.0_rp * metric_inv
        case default
            arh_object%expansion_dirs = arh_object%potential_dirs
            arh_object%projection_dirs = arh_object%dm_dirs
            arh_object%coupling_matrix = 8.0_rp * metric_inv
        end select
        if (.not. allocated(arh_object%projection_dirs)) &
            arh_object%projection_dirs = arh_object%expansion_dirs

        ! build the dense reference Hessian
        allocate(hess(n_param, n_param), e_i(n_param))
        do i = 1, n_param
            e_i = 0.0_rp
            e_i(i) = 1.0_rp
            x_full = ref_unpack_asymm(e_i, 1_ip, n_ao)
            delta_dm = ref_project_symm(x_full, dm_oao)
            allocate(response(n_ao, n_ao, 1))
            response = ref_response_cs(arh_type, delta_dm, dm_diff, fock_diff, &
                                       v_linear_diff, v_nonlinear_diff, metric_inv, &
                                       a_sym, a_inv, a_inv_comb, n_ao, n_diff)
            hess(:, i) = ref_hess_x(x_full, response, dm_oao, fock_oo, fock_vv, n_param)
            deallocate(response, x_full, delta_dm)
        end do

        ! generate a vector confined to the non-redundant subspace and check inverse 
        ! Hessian
        vector = generate_nonredundant_vector(n_param, 1_ip, n_ao, dm_oao)
        call check_inv_hess_x_arh(hess, vector, mu, test_name, passed)

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    end function check_inv_hess_x_arh_cs

    function check_inv_hess_x_arh_os(arh_type, test_name) result(passed)
        !
        ! this function performs the exact-inverse check of the open-shell ARH Woodbury 
        ! preconditioner for a single ARH type
        !
        use otr_arh, only: arh_object
        use otr_oao, only: oao_object
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao_unit_tests, only: ref_unpack_asymm, ref_project_symm, &
                                      ref_project_asymm, ref_pack_asymm, ref_hess_x, &
                                      generate_random_symm_matrix

        character(*), intent(in) :: arh_type, test_name
        logical :: passed

        integer(ip), parameter :: n_diff = 2, n_electrons = 1
        real(rp), parameter :: mu = 0.1_rp

        integer(ip), target :: n_ao_target, n_particle_target, n_param_target
        real(rp), target :: dm_oao(n_ao, n_ao, n_particle), &
                            fock_oo(n_ao, n_ao, n_particle), &
                            fock_vv(n_ao, n_ao, n_particle)
        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_same_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_opp_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_nl_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_same_eff(n_ao, n_ao, n_particle, n_diff), &
                    metric_inv(2 * n_diff, 2 * n_diff), &
                    a_block(2 * n_diff, 2 * n_diff), &
                    a_inv_lin(2 * n_diff, 2 * n_diff), a_inv_comb(n_diff, n_diff)
        integer(ip) :: i, j, k
        real(rp), allocatable :: x_full(:, :, :), delta_dm(:, :, :), &
                                 response(:, :, :), hess(:, :), e_i(:), vector(:)

        ! assume test passes
        passed = .true.

        ! generate a random density matrix per spin channel and the corresponding
        ! occupied-occupied and virtual-virtual Fock matrix blocks
        call generate_fock_partition(n_ao, n_particle, n_electrons, dm_oao, fock_oo, &
                                     fock_vv)

        ! generate spin-resolved history differences
        do k = 1, n_diff
            do j = 1, n_particle
                dm_diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
                v_same_diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
                v_opp_diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
                v_nl_diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! multisecant SR1 keeps the non-linear potential on its own spin-combined
        ! system, while every other type folds it into the same-spin potential
        if (arh_type == "ms_sr1") then
            v_same_eff = v_same_diff
        else
            v_same_eff = v_same_diff + v_nl_diff
        end if

        ! randomly generated (symmetric) metric pseudoinverse, A_sym, and multisecant 
        ! SR1 pseudoinverse cores
        metric_inv = 0.0_rp
        metric_inv(1:n_diff, 1:n_diff) = generate_random_symm_matrix(n_diff)
        metric_inv(n_diff + 1:2*n_diff, n_diff + 1:2*n_diff) = &
            generate_random_symm_matrix(n_diff)
        a_block = generate_random_symm_matrix(2*n_diff)
        a_inv_lin = generate_random_symm_matrix(2*n_diff)
        a_inv_comb = generate_random_symm_matrix(n_diff)

        ! set up the ARH and OAO objects
        n_ao_target = n_ao
        n_particle_target = n_particle
        n_param_target = n_param
        call setup_arh_and_oao_objects(arh_type, dm_oao, fock_oo, fock_vv, &
                                       n_ao_target, n_particle_target, n_param_target)

        ! populate the cached history projections
        allocate(arh_object%dm_dirs(n_param, 2 * n_diff), &
                 arh_object%potential_dirs(n_param, 2 * n_diff), &
                 arh_object%linear_potential_dirs(n_param, 2 * n_diff), &
                 arh_object%nonlinear_potential_dirs(n_param, n_diff))
        do k = 1, n_diff
            do j = 1, n_particle
                arh_object%dm_dirs(:, (j - 1) * n_diff + k) = ref_pack_asymm( &
                    ref_project_asymm(embed_channel(dm_diff(:, :, j, k), j, n_ao, &
                                                    n_particle), dm_oao), n_param)
            end do
            arh_object%potential_dirs(:, k) = ref_pack_asymm(ref_project_asymm( &
                embed_channel(v_same_eff(:, :, 1, k), 1_ip, n_ao, n_particle) + &
                embed_channel(v_opp_diff(:, :, 2, k), 2_ip, n_ao, n_particle), &
                dm_oao), n_param)
            arh_object%potential_dirs(:, n_diff + k) = &
                ref_pack_asymm(ref_project_asymm( &
                    embed_channel(v_same_eff(:, :, 2, k), 2_ip, n_ao, n_particle) + &
                    embed_channel(v_opp_diff(:, :, 1, k), 1_ip, n_ao, n_particle), &
                    dm_oao), n_param)
            arh_object%nonlinear_potential_dirs(:, k) = &
                ref_pack_asymm(ref_project_asymm(v_nl_diff(:, :, :, k), dm_oao), &
                               n_param)
        end do
        arh_object%linear_potential_dirs = arh_object%potential_dirs

        ! assemble the low-rank factors
        select case (arh_type)
        case ("ms_sr1")
            allocate(arh_object%expansion_dirs(n_param, 3 * n_diff), &
                     arh_object%coupling_matrix(3 * n_diff, 3 * n_diff))
            arh_object%expansion_dirs(:, :2 * n_diff) = arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, 2 * n_diff + 1:) = &
                arh_object%nonlinear_potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:2 * n_diff, :2 * n_diff) = 4.0_rp * a_inv_lin
            arh_object%coupling_matrix(2 * n_diff + 1:, 2 * n_diff + 1:) = 4.0_rp * &
                                                                           a_inv_comb
        case ("ms_sp")
            arh_object%expansion_dirs = arh_object%dm_dirs
            arh_object%coupling_matrix = 4.0_rp * &
                                         matmul(metric_inv, matmul(a_block, metric_inv))
        case ("symm_arh")
            allocate(arh_object%expansion_dirs(n_param, 4 * n_diff), &
                     arh_object%coupling_matrix(4 * n_diff, 4 * n_diff))
            arh_object%expansion_dirs(:, :2 * n_diff) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, 2 * n_diff + 1:) = arh_object%potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:2 * n_diff, 2 * n_diff + 1:) = 2.0_rp * &
                                                                       metric_inv
            arh_object%coupling_matrix(2 * n_diff + 1:, :2 * n_diff) = 2.0_rp * &
                                                                       metric_inv
        case ("ms_psb")
            allocate(arh_object%expansion_dirs(n_param, 4 * n_diff), &
                     arh_object%coupling_matrix(4 * n_diff, 4 * n_diff))
            arh_object%expansion_dirs(:, :2 * n_diff) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, 2 * n_diff + 1:) = arh_object%potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:2 * n_diff, :2 * n_diff) = &
                -4.0_rp * matmul(metric_inv, matmul(a_block, metric_inv))
            arh_object%coupling_matrix(:2 * n_diff, 2 * n_diff + 1:) = 4.0_rp * &
                                                                       metric_inv
            arh_object%coupling_matrix(2 * n_diff + 1:, :2 * n_diff) = 4.0_rp * &
                                                                       metric_inv
        case default
            arh_object%expansion_dirs = arh_object%potential_dirs
            arh_object%projection_dirs = arh_object%dm_dirs
            arh_object%coupling_matrix = 4.0_rp * metric_inv
        end select
        if (.not. allocated(arh_object%projection_dirs)) &
            arh_object%projection_dirs = arh_object%expansion_dirs

        ! build the dense reference Hessian one column at a time
        allocate(hess(n_param, n_param), e_i(n_param))
        do i = 1, n_param
            e_i = 0.0_rp
            e_i(i) = 1.0_rp
            x_full = ref_unpack_asymm(e_i, n_particle, n_ao)
            delta_dm = ref_project_symm(x_full, dm_oao)
            allocate(response(n_ao, n_ao, n_particle))
            response = ref_response_os(arh_type, delta_dm, dm_diff, v_same_eff, &
                                       v_opp_diff, v_nl_diff, metric_inv, a_block, &
                                       a_inv_lin, a_inv_comb, n_ao, n_diff)
            hess(:, i) = ref_hess_x(x_full, response, dm_oao, fock_oo, fock_vv, &
                                    n_param)
            deallocate(response, x_full, delta_dm)
        end do

        ! generate a vector confined to the non-redundant subspace and check inverse 
        ! Hessian
        vector = generate_nonredundant_vector(n_param, n_particle, n_ao, dm_oao)
        call check_inv_hess_x_arh(hess, vector, mu, test_name, passed)

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    end function check_inv_hess_x_arh_os

    logical(c_bool) function test_arh_factory_cs() bind(C)
        !
        ! this function tests the subroutine which returns the modified ARH orbital
        ! updating function for the closed-shell case
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, precond_type, &
                                   precond_pd_type, project_type
        use otr_arh, only: arh_factory, arh_object, arh_settings_type, &
                           update_dm_cs_type, update_orbs_arh_cs_ptr, precond_arh_ptr
        use otr_oao_test_reference, only: n_ao
        use otr_arh_test_reference, only: operator(==)
        use otr_oao, only: oao_object, get_energy_cs_type, obj_func_oao_ptr, &
                           precond_pd_oao_ptr, project_oao_ptr
        use otr_oao_unit_tests, only: mock_get_energy_cs, identity_matrix, &
                                      generate_random_density_matrix
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_particle = 1, n_electrons = 2, &
                                  n_param = n_ao * (n_ao - 1) / 2

        real(rp), target :: dm_ao(n_ao, n_ao)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: error
        type(arh_settings_type) :: settings
        procedure(get_energy_cs_type), pointer :: get_energy_funptr
        procedure(update_dm_cs_type), pointer :: update_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        procedure(precond_type), pointer :: precond_arh_funptr
        procedure(precond_pd_type), pointer :: precond_pd_arh_funptr
        procedure(project_type), pointer :: project_arh_funptr

        ! assume tests pass
        test_arh_factory_cs = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%arh_type = "ms_psb"

        ! initialize density matrix and an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        dm_ao = generate_random_density_matrix(n_ao, n_electrons)
        ao_overlap = identity_matrix(n_ao)

        ! initialize callback function pointers
        get_energy_funptr => mock_get_energy_cs
        update_dm_funptr => mock_update_dm_cs

        ! call routine and determine if an error is produced
        call arh_factory(dm_ao, ao_overlap, n_particle, n_ao, get_energy_funptr, &
                         update_dm_funptr, obj_func_arh_funptr, &
                         update_orbs_arh_funptr, precond_arh_funptr, &
                         precond_pd_arh_funptr, project_arh_funptr, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_arh_factory_cs failed: Produced error."
            test_arh_factory_cs = .false.
            return
        end if

        ! determine if ARH object is set up correctly
        if (.not. allocated(arh_object)) then
            write (stderr, *) "test_arh_factory_cs failed: ARH object not allocated."
            test_arh_factory_cs = .false.
            return
        end if
        if (.not. (arh_object%settings == settings)) then
            write (stderr, *) "test_arh_factory_cs failed: Settings not stored "// &
                "correctly."
            test_arh_factory_cs = .false.
        end if
        if (arh_object%n_ao /= n_ao) then
            write (stderr, *) "test_arh_factory_cs failed: Number of AOs not "// &
                "associated correctly."
            test_arh_factory_cs = .false.
        end if
        if (arh_object%n_particle /= n_particle) then
            write (stderr, *) "test_arh_factory_cs failed: Number of particles not "// &
                "associated correctly."
            test_arh_factory_cs = .false.
        end if
        if (arh_object%n_param /= n_param) then
            write (stderr, *) "test_arh_factory_cs failed: Number of parameters "// &
                "not associated correctly."
            test_arh_factory_cs = .false.
        end if
        if (norm2(arh_object%dm_oao(:, :, 1) - dm_ao) > tol) then
            write (stderr, *) "test_arh_factory_cs failed: Density matrix not "// &
                "associated correctly."
            test_arh_factory_cs = .false.
        end if
        if (.not. associated(oao_object%get_energy_cs, mock_get_energy_cs)) then
            write (stderr, *) "test_arh_factory_cs failed: Energy function not "// &
                "stored correctly."
            test_arh_factory_cs = .false.
        end if
        if (.not. associated(arh_object%update_dm_cs, mock_update_dm_cs)) then
            write (stderr, *) "test_arh_factory_cs failed: Density matrix updating "// &
                "function not stored correctly."
            test_arh_factory_cs = .false.
        end if

        ! determine if returned function pointers point to the correct routines
        if (.not. associated(obj_func_arh_funptr, obj_func_oao_ptr)) then
            write (stderr, *) "test_arh_factory_cs failed: Returned objective "// &
                "function is wrong."
            test_arh_factory_cs = .false.
        end if
        if (.not. associated(update_orbs_arh_funptr, update_orbs_arh_cs_ptr)) then
            write (stderr, *) "test_arh_factory_cs failed: Returned orbital "// &
                "updating function is wrong."
            test_arh_factory_cs = .false.
        end if
        if (.not. associated(precond_arh_funptr, precond_arh_ptr)) then
            write (stderr, *) "test_arh_factory_cs failed: Returned level-shifted "// &
                "preconditioner function is wrong."
            test_arh_factory_cs = .false.
        end if
        if (.not. associated(precond_pd_arh_funptr, precond_pd_oao_ptr)) then
            write (stderr, *) "test_arh_factory_cs failed: Returned "// &
                "positive-definite preconditioner function is wrong."
            test_arh_factory_cs = .false.
        end if
        if (.not. associated(project_arh_funptr, project_oao_ptr)) then
            write (stderr, *) "test_arh_factory_cs failed: Returned projection "// &
                "function is wrong."
            test_arh_factory_cs = .false.
        end if

        ! call routine with an unknown ARH type and determine if the sanity check
        ! rejects it
        settings%arh_type = "unknown"
        call arh_factory(dm_ao, ao_overlap, n_particle, n_ao, get_energy_funptr, &
                         update_dm_funptr, obj_func_arh_funptr, &
                         update_orbs_arh_funptr, precond_arh_funptr, &
                         precond_pd_arh_funptr, project_arh_funptr, error, settings)
        if (error == 0) then
            write (stderr, *) "test_arh_factory_cs failed: Error not thrown "// &
                "for unknown ARH type."
            test_arh_factory_cs = .false.
        end if

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    end function test_arh_factory_cs

    logical(c_bool) function test_arh_factory_os() bind(C)
        !
        ! this function tests the subroutine which returns the modified ARH orbital
        ! updating function for the open-shell case
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, precond_type, &
                                   precond_pd_type, project_type
        use otr_arh, only: arh_factory, arh_object, arh_settings_type, &
                           update_dm_os_type, update_orbs_arh_os_ptr, precond_arh_ptr
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_arh_test_reference, only: operator(==)
        use otr_oao, only: oao_object, get_energy_os_type, obj_func_oao_ptr, &
                           precond_pd_oao_ptr, project_oao_ptr
        use otr_oao_unit_tests, only: mock_get_energy_os, identity_matrix, &
                                      generate_random_density_matrix
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_electrons = 2

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: i, error
        type(arh_settings_type) :: settings
        procedure(get_energy_os_type), pointer :: get_energy_funptr
        procedure(update_dm_os_type), pointer :: update_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        procedure(precond_type), pointer :: precond_arh_funptr
        procedure(precond_pd_type), pointer :: precond_pd_arh_funptr
        procedure(project_type), pointer :: project_arh_funptr

        ! assume tests pass
        test_arh_factory_os = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%arh_type = "ms_psb"

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
        call arh_factory(dm_ao, ao_overlap, n_particle, n_ao, get_energy_funptr, &
                         update_dm_funptr, obj_func_arh_funptr, &
                         update_orbs_arh_funptr, precond_arh_funptr, &
                         precond_pd_arh_funptr, project_arh_funptr, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_arh_factory_os failed: Produced error."
            test_arh_factory_os = .false.
            return
        end if

        ! determine if ARH object is set up correctly
        if (.not. allocated(arh_object)) then
            write (stderr, *) "test_arh_factory_os failed: ARH object not allocated."
            test_arh_factory_os = .false.
            return
        end if
        if (.not. (arh_object%settings == settings)) then
            write (stderr, *) "test_arh_factory_os failed: Settings not stored "// &
                "correctly."
            test_arh_factory_os = .false.
        end if
        if (arh_object%n_ao /= n_ao) then
            write (stderr, *) "test_arh_factory_os failed: Number of AOs not "// &
                "associated correctly."
            test_arh_factory_os = .false.
        end if
        if (arh_object%n_particle /= n_particle) then
            write (stderr, *) "test_arh_factory_os failed: Number of particles not "// &
                "associated correctly."
            test_arh_factory_os = .false.
        end if
        if (arh_object%n_param /= n_param) then
            write (stderr, *) "test_arh_factory_os failed: Number of parameters "// &
                "not associated correctly."
            test_arh_factory_os = .false.
        end if
        if (norm2(arh_object%dm_oao - dm_ao) > tol) then
            write (stderr, *) "test_arh_factory_os failed: Density matrices not "// &
                "associated correctly."
            test_arh_factory_os = .false.
        end if
        if (.not. associated(oao_object%get_energy_os, mock_get_energy_os)) then
            write (stderr, *) "test_arh_factory_os failed: Energy function not "// &
                "stored correctly."
            test_arh_factory_os = .false.
        end if
        if (.not. associated(arh_object%update_dm_os, mock_update_dm_os)) then
            write (stderr, *) "test_arh_factory_os failed: Density matrix updating "// &
                "function not stored correctly."
            test_arh_factory_os = .false.
        end if

        ! determine if returned function pointers point to the correct routines
        if (.not. associated(obj_func_arh_funptr, obj_func_oao_ptr)) then
            write (stderr, *) "test_arh_factory_os failed: Returned objective "// &
                "function is wrong."
            test_arh_factory_os = .false.
        end if
        if (.not. associated(update_orbs_arh_funptr, update_orbs_arh_os_ptr)) then
            write (stderr, *) "test_arh_factory_os failed: Returned orbital "// &
                "updating function is wrong."
            test_arh_factory_os = .false.
        end if
        if (.not. associated(precond_arh_funptr, precond_arh_ptr)) then
            write (stderr, *) "test_arh_factory_os failed: Returned level-shifted "// &
                "preconditioner function is wrong."
            test_arh_factory_os = .false.
        end if
        if (.not. associated(precond_pd_arh_funptr, precond_pd_oao_ptr)) then
            write (stderr, *) "test_arh_factory_os failed: Returned "// &
                "positive-definite preconditioner function is wrong."
            test_arh_factory_os = .false.
        end if
        if (.not. associated(project_arh_funptr, project_oao_ptr)) then
            write (stderr, *) "test_arh_factory_os failed: Returned projection "// &
                "function is wrong."
            test_arh_factory_os = .false.
        end if

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    end function test_arh_factory_os

    logical(c_bool) function test_arh_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the ARH
        ! input parameters
        !
        use otr_arh, only: arh_settings_type, arh_sanity_check, arh_types
        use opentrustregion_unit_tests, only: setup_settings

        type(arh_settings_type) :: settings
        integer(ip) :: i, error

        ! assume tests pass
        test_arh_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if all available ARH types are accepted
        do i = 1, size(arh_types)
            settings%arh_type = arh_types(i)
            call arh_sanity_check(settings, error)
            if (error /= 0) then
                write (stderr, *) "test_arh_sanity_check failed: Error thrown for "// &
                    trim(arh_types(i))//" ARH type."
                test_arh_sanity_check = .false.
            end if
        end do

        ! check if ARH type is converted to lowercase
        settings%arh_type = "MS_PSB"
        call arh_sanity_check(settings, error)
        if (settings%arh_type /= "ms_psb") then
            write (stderr, *) "test_arh_sanity_check failed: ARH type not "// &
                "converted to lowercase."
            test_arh_sanity_check = .false.
        end if

        ! check if unknown ARH type is rejected
        settings%arh_type = "unknown"
        call arh_sanity_check(settings, error)
        if (error == 0) then
            write (stderr, *) "test_arh_sanity_check failed: Error not thrown for "// &
                "unknown ARH type."
            test_arh_sanity_check = .false.
        end if

    end function test_arh_sanity_check

    logical(c_bool) function test_update_orbs_arh_cs() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation in the OAO basis for the closed-shell case
        !
        use opentrustregion, only: hess_x_type
        use otr_arh, only: update_orbs_arh_cs, arh_object, hess_x_arh_ptr
        use otr_oao_test_reference, only: n_ao
        use otr_oao, only: oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_unit_tests, only: n_mock_calls, mock_fock_factor, identity_matrix, &
                                      generate_random_density_matrix, &
                                      ref_pack_asymm, ref_project_asymm

        integer(ip), parameter :: n_particle = 1, n_electrons = 2, &
                                  n_param = n_ao * (n_ao - 1) / 2

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: dm_saved(n_ao, n_ao, n_particle), &
                    fock_saved(n_ao, n_ao, n_particle), &
                    v_nonlinear_saved(n_ao, n_ao, n_particle), &
                    dm_saved_2(n_ao, n_ao, n_particle), &
                    fock_saved_2(n_ao, n_ao, n_particle), &
                    v_nonlinear_saved_2(n_ao, n_ao, n_particle), metric(2, 2), &
                    kappa(n_param), grad(n_param), h_diag(n_param), func
        real(rp), allocatable :: a_linear(:, :)
        integer(ip) :: i, j, n_diff, error
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_update_orbs_arh_cs = .true.

        ! set up the OAO object with an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%s_inv_sqrt = identity_matrix(n_ao)
        dm_ao(:, :, 1) = generate_random_density_matrix(n_ao, n_electrons)
        oao_object%dm_ao => dm_ao
        oao_object%dm_oao = dm_ao
        allocate(oao_object%fock_oo(n_ao, n_ao, n_particle), &
                 oao_object%fock_vv(n_ao, n_ao, n_particle))

        ! set up the ARH object the way the ARH factory would
        allocate(arh_object)
        call setup_settings(arh_object%settings)
        arh_object%settings%arh_type = "ms_psb"
        arh_object%n_ao => oao_object%n_ao
        arh_object%n_param => oao_object%n_param
        arh_object%n_particle => oao_object%n_particle
        arh_object%dm_ao => oao_object%dm_ao
        arh_object%s_inv_sqrt => oao_object%s_inv_sqrt
        arh_object%dm_oao => oao_object%dm_oao
        arh_object%fock_oo => oao_object%fock_oo
        arh_object%fock_vv => oao_object%fock_vv
        arh_object%energy => oao_object%energy
        arh_object%update_dm_cs => mock_update_dm_cs

        ! reset mock density matrix updating function
        n_mock_calls = 0

        ! call routine without an orbital rotation for an uninitialized object and
        ! determine if an error is produced
        kappa = 0.0_rp
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_cs(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Produced error."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after the static Hessian part was rebuilt."
            test_update_orbs_arh_cs = .false.
        end if

        ! determine if the energy, Fock matrix and non-linear potential of the density
        ! matrix updating function are picked up
        if (abs(func - sum(dm_ao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect energy."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%fock_oao - mock_potential(mock_fock_factor(1), dm_ao)) > &
            tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect Fock matrix."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%v_nonlinear_oao - &
                  mock_potential(mock_v_nonlinear_factor(1), dm_ao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect "// &
                "non-linear potential."
            test_update_orbs_arh_cs = .false.
        end if

        ! determine if the correct gradient, Hessian diagonal and Hessian linear 
        ! transformation are returned
        if (norm2(grad - arh_object%grad) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Gradient not returned."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(h_diag - arh_object%h_diag) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Hessian diagonal "// &
                "not returned."
            test_update_orbs_arh_cs = .false.
        end if
        if (.not. associated(hess_x_funptr, hess_x_arh_ptr)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Returned Hessian "// &
                "linear transformation is wrong."
            test_update_orbs_arh_cs = .false.
        end if

        ! determine if the history is initialized empty
        if (size(arh_object%dm_list, 4) /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Density matrix "// &
                "history not initialized empty."
            test_update_orbs_arh_cs = .false.
        end if
        if (size(arh_object%dm_dirs, 2) /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Difference "// &
                "directions not initialized empty."
            test_update_orbs_arh_cs = .false.
        end if

        ! call routine again without an orbital rotation and determine if the
        ! quantities of the already initialized object are reused, including the
        ! cached eigendecomposition of the static Hessian part, which should remain
        ! valid since it was not rebuilt
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_cs(kappa, func, grad, h_diag, hess_x_funptr, &
                                          error)
        if (n_mock_calls /= 1) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Quantities "// &
                "recomputed without an orbital rotation."
            test_update_orbs_arh_cs = .false.
        end if
        if (size(arh_object%dm_list, 4) /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: History extended "// &
                "without an orbital rotation."
            test_update_orbs_arh_cs = .false.
        end if
        if (oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Cached "// &
                "eigendecomposition of the static Hessian part marked stale even "// &
                "though the static Hessian part was not rebuilt."
            test_update_orbs_arh_cs = .false.
        end if

        ! save the current quantities, which the history has to retain, and rotate
        ! twice in a row
        dm_saved = arh_object%dm_oao
        fock_saved = arh_object%fock_oao
        v_nonlinear_saved = arh_object%v_nonlinear_oao
        kappa = 0.1_rp
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_cs(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Produced error "// &
                "after the first orbital rotation."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after an orbital rotation."
            test_update_orbs_arh_cs = .false.
        end if
        dm_saved_2 = arh_object%dm_oao
        fock_saved_2 = arh_object%fock_oao
        v_nonlinear_saved_2 = arh_object%v_nonlinear_oao
        call update_orbs_arh_cs(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Produced error "// &
                "after the second orbital rotation."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (size(arh_object%dm_list, 4) /= 2) then
            write (stderr, *) "test_update_orbs_arh_cs failed: History not "// &
                "extended to two entries."
            test_update_orbs_arh_cs = .false.
            return
        end if

        ! determine if the history retains the raw quantities at both columns
        if (norm2(arh_object%dm_list(:, :, :, 1) - dm_saved_2) > tol .or. &
            norm2(arh_object%dm_list(:, :, :, 2) - dm_saved) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect density "// &
                "matrix history."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%fock_list(:, :, :, 1) - fock_saved_2) > tol .or. &
            norm2(arh_object%fock_list(:, :, :, 2) - fock_saved) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect Fock "// &
                "matrix history."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%v_nonlinear_list(:, :, :, 1) - v_nonlinear_saved_2) > &
            tol .or. &
            norm2(arh_object%v_nonlinear_list(:, :, :, 2) - v_nonlinear_saved) > tol) &
        then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect "// &
                "non-linear potential history."
            test_update_orbs_arh_cs = .false.
        end if

        ! determine if the quantities the density matrix and Fock matrix differences
        ! feed into are built from the history and the current quantities
        if (.not. allocated(arh_object%metric_inv)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: ARH metric "// &
                "pseudoinverse not constructed."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (.not. allocated(arh_object%a_sym)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Symmetrized A "// &
                "matrix not constructed."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (.not. allocated(arh_object%dm_dirs)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Density matrix "// &
                "difference directions not constructed."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (.not. allocated(arh_object%potential_dirs)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Potential "// &
                "difference directions not constructed."
            test_update_orbs_arh_cs = .false.
            return
        end if
        metric(1, 1) = sum((dm_saved_2 - arh_object%dm_oao)**2)
        metric(2, 2) = sum((dm_saved - arh_object%dm_oao)**2)
        metric(1, 2) = sum((dm_saved_2 - arh_object%dm_oao) * &
                           (dm_saved - arh_object%dm_oao))
        metric(2, 1) = metric(1, 2)
        if (norm2(matmul(matmul(arh_object%metric_inv, metric), &
                         arh_object%metric_inv) - arh_object%metric_inv) > tol * &
            (1.0_rp + norm2(metric) * norm2(arh_object%metric_inv))) then
            write (stderr, *) "test_update_orbs_arh_cs failed: ARH metric "// &
                "pseudoinverse is not the pseudoinverse of the metric."
            test_update_orbs_arh_cs = .false.
        end if
        if (abs(arh_object%a_sym(1, 1) - sum((dm_saved_2 - arh_object%dm_oao) * &
                                             (fock_saved_2 - arh_object%fock_oao))) > &
            tol .or. &
            abs(arh_object%a_sym(2, 2) - sum((dm_saved - arh_object%dm_oao) * &
                                             (fock_saved - arh_object%fock_oao))) > &
            tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect "// &
                "diagonal block of the symmetrized A matrix."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%dm_dirs(:, 1) - &
                  ref_pack_asymm(ref_project_asymm(dm_saved_2 - arh_object%dm_oao, &
                                                   arh_object%dm_oao), n_param)) > &
            tol .or. &
            norm2(arh_object%dm_dirs(:, 2) - &
                  ref_pack_asymm(ref_project_asymm(dm_saved - arh_object%dm_oao, &
                                                   arh_object%dm_oao), n_param)) > &
            tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect density "// &
                "matrix difference directions."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%potential_dirs(:, 1) - ref_pack_asymm(ref_project_asymm( &
                      fock_saved_2 - arh_object%fock_oao, arh_object%dm_oao), &
                      n_param)) > tol .or. &
            norm2(arh_object%potential_dirs(:, 2) - ref_pack_asymm(ref_project_asymm( &
                      fock_saved - arh_object%fock_oao, arh_object%dm_oao), n_param)) &
            > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect potential "// &
                "difference directions."
            test_update_orbs_arh_cs = .false.
        end if

        ! call routine for multisecant SR1 and determine if the separately regularized
        ! multisecant SR1 systems are constructed
        arh_object%settings%arh_type = "ms_sr1"
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_cs(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Produced error for "// &
                "multisecant SR1."
            test_update_orbs_arh_cs = .false.
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "for multisecant SR1."
            test_update_orbs_arh_cs = .false.
        end if
        if (.not. allocated(arh_object%a_inv)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Multisecant SR1 "// &
                "system not constructed."
            test_update_orbs_arh_cs = .false.
        end if
        if (.not. allocated(arh_object%a_inv_comb)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Spin-combined "// &
                "multisecant SR1 system not constructed."
            test_update_orbs_arh_cs = .false.
        end if
        if (.not. allocated(arh_object%linear_potential_dirs)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Linear potential "// &
                "difference directions not constructed."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (.not. allocated(arh_object%nonlinear_potential_dirs)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Non-linear "// &
                "potential difference directions not constructed."
            test_update_orbs_arh_cs = .false.
            return
        end if

        ! determine if the quantities the density matrix and potential differences
        ! feed into are built from the history and the current quantities
        n_diff = size(arh_object%dm_list, 4)
        allocate(a_linear(n_diff, n_diff))
        do i = 1, n_diff
            do j = 1, n_diff
                a_linear(i, j) = sum((arh_object%dm_list(:, :, :, i) - &
                                      arh_object%dm_oao) * &
                                     (arh_object%fock_list(:, :, :, j) - &
                                      arh_object%fock_oao - &
                                      (arh_object%v_nonlinear_list(:, :, :, j) - &
                                       arh_object%v_nonlinear_oao)))
            end do
        end do
        a_linear = 0.5_rp * (a_linear + transpose(a_linear))
        if (norm2(matmul(matmul(arh_object%a_inv, a_linear), arh_object%a_inv) - &
                  arh_object%a_inv) > tol .or. norm2(arh_object%a_inv) < tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Multisecant SR1 "// &
                "pseudoinverse is not the pseudoinverse of the linear system built "// &
                "from the history differences."
            test_update_orbs_arh_cs = .false.
        end if
        deallocate(a_linear)
        do i = 1, size(arh_object%dm_list, 4)
            if (norm2(arh_object%linear_potential_dirs(:, i) - ref_pack_asymm( &
                      ref_project_asymm( &
                          arh_object%fock_list(:, :, :, i) - arh_object%fock_oao - &
                          (arh_object%v_nonlinear_list(:, :, :, i) - &
                           arh_object%v_nonlinear_oao), arh_object%dm_oao), n_param)) &
                > tol) then
                write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect "// &
                    "linear potential difference directions."
                test_update_orbs_arh_cs = .false.
            end if
            if (norm2(arh_object%nonlinear_potential_dirs(:, i) - ref_pack_asymm( &
                      ref_project_asymm( &
                          arh_object%v_nonlinear_list(:, :, :, i) - &
                          arh_object%v_nonlinear_oao, arh_object%dm_oao), n_param)) > &
                tol) then
                write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect "// &
                    "non-linear potential difference directions."
                test_update_orbs_arh_cs = .false.
            end if
        end do

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    end function test_update_orbs_arh_cs

    logical(c_bool) function test_update_orbs_arh_os() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation in the OAO basis for the open-shell case
        !
        use opentrustregion, only: hess_x_type
        use otr_arh, only: update_orbs_arh_os, arh_object, hess_x_arh_ptr
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao, only: oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_unit_tests, only: n_mock_calls, identity_matrix, &
                                      generate_random_density_matrix, &
                                      ref_pack_asymm, ref_project_asymm

        integer(ip), parameter :: n_electrons = 2

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: dm_saved(n_ao, n_ao, n_particle), &
                    v_same_spin_saved(n_ao, n_ao, n_particle), &
                    v_opposite_spin_saved(n_ao, n_ao, n_particle), &
                    v_nonlinear_saved(n_ao, n_ao, n_particle), &
                    dm_saved_2(n_ao, n_ao, n_particle), &
                    v_same_spin_saved_2(n_ao, n_ao, n_particle), &
                    v_opposite_spin_saved_2(n_ao, n_ao, n_particle), &
                    v_nonlinear_saved_2(n_ao, n_ao, n_particle), kappa(n_param), &
                    grad(n_param), h_diag(n_param), func
        integer(ip) :: i, k, col, n_diff, error
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_update_orbs_arh_os = .true.

        ! set up the OAO object with an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%s_inv_sqrt = identity_matrix(n_ao)
        do i = 1, n_particle
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_electrons)
        end do
        oao_object%dm_ao => dm_ao
        oao_object%dm_oao = dm_ao
        allocate(oao_object%fock_oo(n_ao, n_ao, n_particle), &
                 oao_object%fock_vv(n_ao, n_ao, n_particle))

        ! set up the ARH object the way the ARH factory would
        allocate(arh_object)
        call setup_settings(arh_object%settings)
        arh_object%settings%arh_type = "ms_psb"
        arh_object%n_ao => oao_object%n_ao
        arh_object%n_param => oao_object%n_param
        arh_object%n_particle => oao_object%n_particle
        arh_object%dm_ao => oao_object%dm_ao
        arh_object%s_inv_sqrt => oao_object%s_inv_sqrt
        arh_object%dm_oao => oao_object%dm_oao
        arh_object%fock_oo => oao_object%fock_oo
        arh_object%fock_vv => oao_object%fock_vv
        arh_object%energy => oao_object%energy
        arh_object%update_dm_os => mock_update_dm_os

        ! reset mock density matrix updating function
        n_mock_calls = 0

        ! call routine without an orbital rotation for an uninitialized object and
        ! determine if an error is produced
        kappa = 0.0_rp
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_os(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Produced error."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_os failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after the static Hessian part was rebuilt."
            test_update_orbs_arh_os = .false.
        end if

        ! determine if the energy and the potentials of the density matrix updating
        ! function are picked up
        if (abs(func - sum(dm_ao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect energy."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_same_spin_oao - &
                  mock_potential(mock_v_same_spin_factor(1), dm_ao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect same-spin "// &
                "potential."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_opposite_spin_oao - &
                  mock_potential(mock_v_opposite_spin_factor(1), dm_ao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "opposite-spin potential."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_nonlinear_oao - &
                  mock_potential(mock_v_nonlinear_factor(1), dm_ao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "non-linear potential."
            test_update_orbs_arh_os = .false.
        end if

        ! determine if the gradient and Hessian diagonal are returned
        if (norm2(grad - arh_object%grad) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Gradient not returned."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(h_diag - arh_object%h_diag) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Hessian diagonal "// &
                "not returned."
            test_update_orbs_arh_os = .false.
        end if

        ! determine if the history is initialized empty
        if (size(arh_object%dm_list, 4) /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Density matrix "// &
                "history not initialized empty."
            test_update_orbs_arh_os = .false.
        end if
        if (size(arh_object%dm_dirs, 2) /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Difference "// &
                "directions not initialized empty."
            test_update_orbs_arh_os = .false.
        end if

        ! determine if the returned Hessian linear transformation is correct
        if (.not. associated(hess_x_funptr, hess_x_arh_ptr)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Returned Hessian "// &
                "linear transformation is wrong."
            test_update_orbs_arh_os = .false.
        end if

        ! call routine again without an orbital rotation and determine if the
        ! quantities of the already initialized object are reused, including the
        ! cached eigendecomposition of the static Hessian part, which should remain
        ! valid since it was not rebuilt
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_os(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (n_mock_calls /= 1) then
            write (stderr, *) "test_update_orbs_arh_os failed: Quantities "// &
                "recomputed without an orbital rotation."
            test_update_orbs_arh_os = .false.
        end if
        if (size(arh_object%dm_list, 4) /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: History extended "// &
                "without an orbital rotation."
            test_update_orbs_arh_os = .false.
        end if
        if (oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_os failed: Cached "// &
                "eigendecomposition of the static Hessian part marked stale even "// &
                "though the static Hessian part was not rebuilt."
            test_update_orbs_arh_os = .false.
        end if

        ! save the current quantities, which the history has to retain, and rotate
        ! twice in a row
        dm_saved = arh_object%dm_oao
        v_same_spin_saved = arh_object%v_same_spin_oao
        v_opposite_spin_saved = arh_object%v_opposite_spin_oao
        v_nonlinear_saved = arh_object%v_nonlinear_oao
        kappa = 0.1_rp
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_os(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Produced error "// &
                "after the first orbital rotation."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_os failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after an orbital rotation."
            test_update_orbs_arh_os = .false.
        end if
        dm_saved_2 = arh_object%dm_oao
        v_same_spin_saved_2 = arh_object%v_same_spin_oao
        v_opposite_spin_saved_2 = arh_object%v_opposite_spin_oao
        v_nonlinear_saved_2 = arh_object%v_nonlinear_oao
        call update_orbs_arh_os(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Produced error "// &
                "after the second orbital rotation."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (size(arh_object%dm_list, 4) /= 2) then
            write (stderr, *) "test_update_orbs_arh_os failed: History not "// &
                "extended to two entries."
            test_update_orbs_arh_os = .false.
            return
        end if
        n_diff = size(arh_object%dm_list, 4)
        if (norm2(arh_object%dm_list(:, :, :, 1) - dm_saved_2) > tol .or. &
            norm2(arh_object%dm_list(:, :, :, 2) - dm_saved) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect density "// &
                "matrix history."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_same_spin_list(:, :, :, 1) - v_same_spin_saved_2) > tol &
            .or. norm2(arh_object%v_same_spin_list(:, :, :, 2) - v_same_spin_saved) > &
            tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect same-spin "// &
                "potential history."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_opposite_spin_list(:, :, :, 1) - &
                  v_opposite_spin_saved_2) > tol .or. &
            norm2(arh_object%v_opposite_spin_list(:, :, :, 2) - v_opposite_spin_saved) &
                  > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "opposite-spin potential history."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_nonlinear_list(:, :, :, 1) - v_nonlinear_saved_2) > &
            tol .or. &
            norm2(arh_object%v_nonlinear_list(:, :, :, 2) - v_nonlinear_saved) > tol) &
        then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "non-linear potential history."
            test_update_orbs_arh_os = .false.
        end if

        ! determine if the quantities the density matrix and potential differences
        ! feed into are built from the history and the current quantities
        if (.not. allocated(arh_object%metric_inv)) then
            write (stderr, *) "test_update_orbs_arh_os failed: ARH metric "// &
                "pseudoinverse not constructed."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (.not. allocated(arh_object%a_sym)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Symmetrized A "// &
                "matrix not constructed."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (.not. allocated(arh_object%dm_dirs)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Density matrix "// &
                "difference directions not constructed."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (.not. allocated(arh_object%potential_dirs)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Potential "// &
                "difference directions not constructed."
            test_update_orbs_arh_os = .false.
            return
        end if
        do i = 1, n_particle
            do k = 1, n_diff
                col = (i - 1) * n_diff + k
                if (k == 1) then
                    if (abs(arh_object%a_sym(col, col) - &
                            sum((dm_saved_2(:, :, i) - arh_object%dm_oao(:, :, i)) * &
                                (v_same_spin_saved_2(:, :, i) - &
                                 arh_object%v_same_spin_oao(:, :, i) + &
                                 v_nonlinear_saved_2(:, :, i) - &
                                 arh_object%v_nonlinear_oao(:, :, i)))) > tol) then
                        write (stderr, *) "test_update_orbs_arh_os failed: "// &
                            "Incorrect diagonal block of the symmetrized A matrix."
                        test_update_orbs_arh_os = .false.
                    end if
                    if (norm2(arh_object%dm_dirs(:, col) - ref_pack_asymm( &
                            ref_project_asymm(embed_channel( &
                                dm_saved_2(:, :, i) - arh_object%dm_oao(:, :, i), i, &
                                n_ao, n_particle), arh_object%dm_oao), n_param)) > &
                        tol) then
                        write (stderr, *) "test_update_orbs_arh_os failed: "// &
                            "Incorrect density matrix difference directions."
                        test_update_orbs_arh_os = .false.
                    end if
                    if (norm2(arh_object%potential_dirs(:, col) - ref_pack_asymm( &
                            ref_project_asymm(embed_channel( &
                                v_same_spin_saved_2(:, :, i) - &
                                arh_object%v_same_spin_oao(:, :, i) + &
                                v_nonlinear_saved_2(:, :, i) - &
                                arh_object%v_nonlinear_oao(:, :, i), i, n_ao, &
                                n_particle) + embed_channel( &
                                v_opposite_spin_saved_2(:, :, 3 - i) - &
                                arh_object%v_opposite_spin_oao(:, :, 3 - i), 3 - i, &
                                n_ao, n_particle), arh_object%dm_oao), n_param)) > &
                        tol) then
                        write (stderr, *) "test_update_orbs_arh_os failed: "// &
                            "Incorrect potential difference directions."
                        test_update_orbs_arh_os = .false.
                    end if
                else
                    if (abs(arh_object%a_sym(col, col) - &
                            sum((dm_saved(:, :, i) - arh_object%dm_oao(:, :, i)) * &
                                (v_same_spin_saved(:, :, i) - &
                                 arh_object%v_same_spin_oao(:, :, i) + &
                                 v_nonlinear_saved(:, :, i) - &
                                 arh_object%v_nonlinear_oao(:, :, i)))) > tol) then
                        write (stderr, *) "test_update_orbs_arh_os failed: "// &
                            "Incorrect diagonal block of the symmetrized A matrix."
                        test_update_orbs_arh_os = .false.
                    end if
                    if (norm2(arh_object%dm_dirs(:, col) - ref_pack_asymm( &
                            ref_project_asymm(embed_channel( &
                                dm_saved(:, :, i) - arh_object%dm_oao(:, :, i), i, &
                                n_ao, n_particle), arh_object%dm_oao), n_param)) > &
                        tol) then
                        write (stderr, *) "test_update_orbs_arh_os failed: "// &
                            "Incorrect density matrix difference directions."
                        test_update_orbs_arh_os = .false.
                    end if
                    if (norm2(arh_object%potential_dirs(:, col) - ref_pack_asymm( &
                            ref_project_asymm(embed_channel( &
                                v_same_spin_saved(:, :, i) - &
                                arh_object%v_same_spin_oao(:, :, i) + &
                                v_nonlinear_saved(:, :, i) - &
                                arh_object%v_nonlinear_oao(:, :, i), i, n_ao, &
                                n_particle) + embed_channel( &
                                v_opposite_spin_saved(:, :, 3 - i) - &
                                arh_object%v_opposite_spin_oao(:, :, 3 - i), 3 - i, &
                                n_ao, n_particle), arh_object%dm_oao), n_param)) > &
                        tol) then
                        write (stderr, *) "test_update_orbs_arh_os failed: "// &
                            "Incorrect potential difference directions."
                        test_update_orbs_arh_os = .false.
                    end if
                end if
            end do
        end do
        if (abs(arh_object%a_sym(1, n_diff + 1) - &
                (norm2(dm_saved_2(:, :, 2) - arh_object%dm_oao(:, :, 2)) * &
                 sum((dm_saved_2(:, :, 2) - arh_object%dm_oao(:, :, 2)) * &
                     (v_opposite_spin_saved_2(:, :, 2) - &
                      arh_object%v_opposite_spin_oao(:, :, 2))) + &
                 norm2(dm_saved_2(:, :, 1) - arh_object%dm_oao(:, :, 1)) * &
                 sum((dm_saved_2(:, :, 1) - arh_object%dm_oao(:, :, 1)) * &
                     (v_opposite_spin_saved_2(:, :, 1) - &
                      arh_object%v_opposite_spin_oao(:, :, 1)))) / &
                (norm2(dm_saved_2(:, :, 1) - arh_object%dm_oao(:, :, 1)) + &
                 norm2(dm_saved_2(:, :, 2) - arh_object%dm_oao(:, :, 2)))) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "cross-symmetrized off-diagonal block of the A matrix."
            test_update_orbs_arh_os = .false.
        end if

        ! call routine for multisecant SR1 and determine if the spin-separated and
        ! spin-combined multisecant SR1 systems are constructed
        arh_object%settings%arh_type = "ms_sr1"
        call update_orbs_arh_os(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Produced error for "// &
                "multisecant SR1."
            test_update_orbs_arh_os = .false.
        end if
        if (.not. allocated(arh_object%a_inv)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Spin-separated "// &
                "multisecant SR1 system not constructed."
            test_update_orbs_arh_os = .false.
        end if
        if (.not. allocated(arh_object%a_inv_comb)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Spin-combined "// &
                "multisecant SR1 system not constructed."
            test_update_orbs_arh_os = .false.
        end if
        if (.not. allocated(arh_object%linear_potential_dirs)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Spin-separated "// &
                "potential difference directions not constructed."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (.not. allocated(arh_object%nonlinear_potential_dirs)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Non-linear "// &
                "potential difference directions not constructed."
            test_update_orbs_arh_os = .false.
            return
        end if

        ! determine that the potential difference directions multisecant SR1 keeps
        ! spin-separated are built from the history and the current quantities
        do i = 1, size(arh_object%dm_list, 4)
            if (norm2(arh_object%nonlinear_potential_dirs(:, i) - &
                      ref_pack_asymm(ref_project_asymm( &
                          arh_object%v_nonlinear_list(:, :, :, i) - &
                          arh_object%v_nonlinear_oao, arh_object%dm_oao), n_param)) > &
                tol) then
                write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                    "non-linear potential difference directions."
                test_update_orbs_arh_os = .false.
            end if
            if (norm2(arh_object%linear_potential_dirs(:, i) - &
                      ref_pack_asymm(ref_project_asymm( &
                          embed_channel(arh_object%v_same_spin_list(:, :, 1, i) - &
                                        arh_object%v_same_spin_oao(:, :, 1), 1_ip, &
                                        n_ao, n_particle) + &
                          embed_channel(arh_object%v_opposite_spin_list(:, :, 2, i) - &
                                        arh_object%v_opposite_spin_oao(:, :, 2), 2_ip, &
                                        n_ao, n_particle), arh_object%dm_oao), &
                          n_param)) > tol) then
                write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                    "spin-separated potential difference directions."
                test_update_orbs_arh_os = .false.
            end if
        end do

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    end function test_update_orbs_arh_os

    logical(c_bool) function test_hess_x_arh() bind(C)
        !
        ! this function tests the subroutine which defines the Hessian linear
        ! transformation on the basis of augmented Roothaan-Hall and related methods
        !
        use otr_arh, only: hess_x_arh, arh_object
        use otr_oao_test_reference, only: n_ao, n_particle
        use otr_oao_unit_tests, only: ref_unpack_asymm, ref_pack_asymm, &
                                      ref_project_asymm, ref_project_symm, ref_hess_x, &
                                      generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_diff = 2

        integer(ip), target :: n_particle_target, n_param, n_ao_target
        real(rp), target :: dm_oao(n_ao, n_ao, n_particle), &
                            fock_oo(n_ao, n_ao, n_particle), &
                            fock_vv(n_ao, n_ao, n_particle)
        real(rp) :: expansion_history(n_ao, n_ao, n_particle, n_diff), &
                    projection_history(n_ao, n_ao, n_particle, n_diff), &
                    coupling(n_diff, n_diff), projected(n_diff), coupled(n_diff)
        integer(ip) :: i, j, error
        real(rp), allocatable :: x(:), x_full(:, :, :), delta_dm(:, :, :), &
                                 response(:, :, :), hess_x(:), expected_hess_x(:)

        ! assume tests pass
        test_hess_x_arh = .true.

        ! generate random density matrices, Fock matrix contributions and history
        ! matrices the low-rank directions are built from
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, 2_ip)
        dm_oao(:, :, 2) = generate_random_density_matrix(n_ao, 1_ip)
        do j = 1, n_particle
            fock_oo(:, :, j) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, j) = generate_random_symm_matrix(n_ao)
            do i = 1, n_diff
                expansion_history(:, :, j, i) = generate_random_symm_matrix(n_ao)
                projection_history(:, :, j, i) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! asymmetric coupling matrix, so that a swap of the expansion and projection
        ! directions would not go unnoticed
        coupling = reshape([2.0_rp, -1.0_rp, 3.0_rp, 0.5_rp], [n_diff, n_diff])

        ! set up the ARH object with the quantities the Hessian linear transformation
        ! requires
        allocate(arh_object)
        n_ao_target = n_ao
        arh_object%n_ao => n_ao_target
        arh_object%n_particle => n_particle_target
        arh_object%n_param => n_param
        arh_object%dm_oao => dm_oao
        arh_object%fock_oo => fock_oo
        arh_object%fock_vv => fock_vv
        arh_object%coupling_matrix = coupling

        ! set up the closed-shell case
        n_particle_target = 1
        n_param = n_ao * (n_ao - 1) / 2
        allocate(arh_object%expansion_dirs(n_param, n_diff), &
                 arh_object%projection_dirs(n_param, n_diff))
        do i = 1, n_diff
            arh_object%expansion_dirs(:, i) = ref_pack_asymm( &
                ref_project_asymm(expansion_history(:, :, 1:1, i), dm_oao(:, :, 1:1)), &
                n_param)
            arh_object%projection_dirs(:, i) = ref_pack_asymm( &
                ref_project_asymm(projection_history(:, :, 1:1, i), &
                                  dm_oao(:, :, 1:1)), n_param)
        end do
        allocate(x(n_param))
        call random_number(x)
        x_full = ref_unpack_asymm(x, n_particle_target, n_ao)
        delta_dm = ref_project_symm(x_full, dm_oao(:, :, 1:1))
        do i = 1, n_diff
            projected(i) = 0.5_rp * sum(projection_history(:, :, 1:1, i) * delta_dm)
        end do
        coupled = matmul(coupling, projected)
        allocate(response(n_ao, n_ao, n_particle_target))
        response = 0.0_rp
        do i = 1, n_diff
            response = response + coupled(i) / 4.0_rp * expansion_history(:, :, 1:1, i)
        end do
        expected_hess_x = ref_hess_x(x_full, response, dm_oao(:, :, 1:1), &
                                     fock_oo(:, :, 1:1), fock_vv(:, :, 1:1), n_param)

        ! call routine and determine if values of resulting Hessian linear
        ! transformation match
        allocate(hess_x(n_param))
        call hess_x_arh(x, hess_x, error)
        if (error /= 0) then
            write (stderr, *) "test_hess_x_arh failed: Produced error for "// &
                "closed-shell case."
            test_hess_x_arh = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write (stderr, *) "test_hess_x_arh failed: Incorrect Hessian linear "// &
                "transformation for closed-shell case."
            test_hess_x_arh = .false.
        end if

        ! test whether absent low-rank part leaves the static part alone 
        deallocate(arh_object%expansion_dirs, arh_object%projection_dirs, &
                  arh_object%coupling_matrix)
        response = 0.0_rp
        expected_hess_x = ref_hess_x(x_full, response, dm_oao(:, :, 1:1), &
                                     fock_oo(:, :, 1:1), fock_vv(:, :, 1:1), n_param)
        call hess_x_arh(x, hess_x, error)
        if (error /= 0) then
            write (stderr, *) "test_hess_x_arh failed: Produced error for empty "// &
                "history."
            test_hess_x_arh = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write (stderr, *) "test_hess_x_arh failed: Incorrect Hessian linear "// &
                "transformation for empty history."
            test_hess_x_arh = .false.
        end if
        deallocate(x, hess_x, response, x_full, delta_dm)

        ! set up the open-shell case
        n_particle_target = n_particle
        n_param = n_particle_target * n_ao * (n_ao - 1) / 2
        arh_object%coupling_matrix = coupling
        allocate(arh_object%expansion_dirs(n_param, n_diff), &
                 arh_object%projection_dirs(n_param, n_diff))
        do i = 1, n_diff
            arh_object%expansion_dirs(:, i) = ref_pack_asymm( &
                ref_project_asymm(expansion_history(:, :, :, i), dm_oao), n_param)
            arh_object%projection_dirs(:, i) = ref_pack_asymm( &
                ref_project_asymm(projection_history(:, :, :, i), dm_oao), n_param)
        end do
        allocate(x(n_param))
        call random_number(x)
        x_full = ref_unpack_asymm(x, n_particle_target, n_ao)
        delta_dm = ref_project_symm(x_full, dm_oao)
        do i = 1, n_diff
            projected(i) = 0.5_rp * sum(projection_history(:, :, :, i) * delta_dm)
        end do
        coupled = matmul(coupling, projected)
        allocate(response(n_ao, n_ao, n_particle_target))
        response = 0.0_rp
        do i = 1, n_diff
            response = response + coupled(i) / 2.0_rp * expansion_history(:, :, :, i)
        end do
        expected_hess_x = ref_hess_x(x_full, response, dm_oao, fock_oo, fock_vv, &
                                     n_param)

        ! call routine and determine if values of resulting Hessian linear
        ! transformation match
        allocate(hess_x(n_param))
        call hess_x_arh(x, hess_x, error)
        if (error /= 0) then
            write (stderr, *) "test_hess_x_arh failed: Produced error for "// &
                "open-shell case."
            test_hess_x_arh = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write (stderr, *) "test_hess_x_arh failed: Incorrect Hessian linear "// &
                "transformation for open-shell case."
            test_hess_x_arh = .false.
        end if
        deallocate(x, hess_x, response)

        ! deallocate ARH object
        deallocate(arh_object)

    end function test_hess_x_arh

    logical(c_bool) function test_inv_hess_x_arh() bind(C)
        !
        ! this function tests the exact, optionally level-shifted inverse of the
        ! approximate Hessian for every ARH type in both the closed- and the
        ! open-shell case
        !
        use otr_arh, only: arh_types

        integer(ip) :: i

        ! assume tests pass
        test_inv_hess_x_arh = .true.

        ! closed- and open-shell cases are driven through the shared per-shell checker, 
        ! one case per ARH type
        do i = 1, size(arh_types)
            if (.not. check_inv_hess_x_arh_cs(arh_types(i), &
                    "test_inv_hess_x_arh failed for "//trim(arh_types(i))//"_cs")) &
                test_inv_hess_x_arh = .false.
            if (.not. check_inv_hess_x_arh_os(arh_types(i), &
                    "test_inv_hess_x_arh failed for "//trim(arh_types(i))//"_os")) &
                test_inv_hess_x_arh = .false.
        end do

    end function test_inv_hess_x_arh

    logical(c_bool) function test_precond_arh() bind(C)
        !
        ! this function tests the preconditioner entry point, which merely hands the
        ! level shift to inv_hess_x_arh, so it only checks that rather than
        ! re-deriving the inverse the inv_hess_x_arh tests already cover
        !
        use otr_arh, only: precond_arh, inv_hess_x_arh, arh_object
        use otr_oao, only: oao_object
        use otr_oao_test_reference, only: n_ao

        real(rp), parameter :: mu = 0.1_rp

        integer(ip), target :: n_ao_target, n_particle_target, n_param_target
        real(rp), target :: dm_oao(n_ao, n_ao, 1), fock_oo(n_ao, n_ao, 1), &
                            fock_vv(n_ao, n_ao, 1)
        integer(ip) :: n_param, error
        real(rp), allocatable :: residual(:), preconditioned(:), inverted(:)

        ! assume test passes
        test_precond_arh = .true.

        ! get number of parameters
        n_param = n_ao * (n_ao - 1) / 2

        ! generate a random density matrix and the oo- and vv-blocks of the Fock matrix
        call generate_fock_partition(n_ao, 1_ip, 2_ip, dm_oao, fock_oo, fock_vv)

        ! set up the ARH and OAO objects
        n_ao_target = n_ao
        n_particle_target = 1
        n_param_target = n_param
        call setup_arh_and_oao_objects("arh", dm_oao, fock_oo, fock_vv, n_ao_target, &
                                       n_particle_target, n_param_target)

        ! call both entry points with the same level shift and determine if the
        ! preconditioner reproduces the level-shifted inverse
        allocate(residual(n_param), preconditioned(n_param), inverted(n_param))
        call random_number(residual)
        call precond_arh(residual, mu, preconditioned, error)
        if (error /= 0) then
            write (stderr, *) "test_precond_arh failed: Produced error."
            test_precond_arh = .false.
        end if
        call inv_hess_x_arh(residual, inverted, error, mu)
        if (norm2(preconditioned - inverted) > tol) then
            write (stderr, *) "test_precond_arh failed: Preconditioner does not "// &
                "reproduce the level-shifted inverse."
            test_precond_arh = .false.
        end if

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    end function test_precond_arh

    logical(c_bool) function test_init_arh_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the ARH settings
        !
        use otr_arh, only: arh_settings_type, default_settings => default_arh_settings
        use otr_arh_test_reference, only: operator(==)

        type(arh_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_init_arh_settings = .true.

        ! initialize settings
        call settings%init(error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_init_arh_settings failed: Function raised error."
            test_init_arh_settings = .false.
        end if

        ! check settings
        if (.not. (settings == default_settings)) then
            write (stderr, *) "test_init_arh_settings failed: Settings not "// &
                "initialized correctly."
            test_init_arh_settings = .false.
        end if

    end function test_init_arh_settings

    logical(c_bool) function test_arh_deconstructor() bind(C)
        !
        ! this function tests the subroutine which deallocates the ARH objects
        !
        use otr_arh, only: arh_deconstructor, arh_object
        use otr_oao, only: oao_object

        ! assume tests pass
        test_arh_deconstructor = .true.

        ! allocate ARH and OAO objects
        if (.not. allocated(arh_object)) allocate(arh_object)
        if (.not. allocated(oao_object)) allocate(oao_object)

        ! call routine and determine if both objects are deallocated
        call arh_deconstructor()
        if (allocated(arh_object)) then
            write (stderr, *) "test_arh_deconstructor failed: ARH object not "// &
                "deallocated."
            test_arh_deconstructor = .false.
        end if
        if (allocated(oao_object)) then
            write (stderr, *) "test_arh_deconstructor failed: OAO object not "// &
                "deallocated."
            test_arh_deconstructor = .false.
        end if

        ! call routine again and determine if already deallocated objects are handled
        call arh_deconstructor()
        if (allocated(arh_object)) then
            write (stderr, *) "test_arh_deconstructor failed: Already deallocated "// &
                "ARH object not handled."
            test_arh_deconstructor = .false.
        end if
        if (allocated(oao_object)) then
            write (stderr, *) "test_arh_deconstructor failed: Already deallocated "// &
                "OAO object not handled."
            test_arh_deconstructor = .false.
        end if

    end function test_arh_deconstructor

    logical(c_bool) function test_cache_history_projections() bind(C)
        !
        ! this function tests the routine which caches the packed history-projection
        ! directions the low-rank part of the approximate Hessian is built from, by
        ! projecting and packing every history entry
        !
        use otr_arh, only: cache_history_projections
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao_unit_tests, only: generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_list = 2, n_electrons = 1

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    v_diff(n_ao, n_ao, n_particle, n_list), &
                    empty_v_diff(n_ao, n_ao, n_particle, 0)
        real(rp), allocatable :: dirs(:, :), expected(:, :)
        integer(ip) :: j, k

        ! assume tests pass
        test_cache_history_projections = .true.

        ! generate a random density matrix and random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)
            do k = 1, n_list
                v_diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! expected projections
        expected = ref_cache_dirs(v_diff, dm_oao, n_param)

        ! call routine and determine if dimensions and values of the resulting
        ! projections match
        call cache_history_projections(v_diff, dm_oao, n_list, n_param, dirs)
        if (size(dirs, 1) /= n_param .or. size(dirs, 2) /= n_list) then
            write (stderr, *) "test_cache_history_projections failed: Incorrect "// &
                "dimensions of directions."
            test_cache_history_projections = .false.
        else if (norm2(dirs - expected) > tol) then
            write (stderr, *) "test_cache_history_projections failed: Incorrect "// &
                "directions."
            test_cache_history_projections = .false.
        end if
        deallocate(dirs, expected)

        ! call routine for an empty history and determine if no directions are
        ! returned
        call cache_history_projections(empty_v_diff, dm_oao, 0_ip, n_param, dirs)
        if (size(dirs, 1) /= n_param .or. size(dirs, 2) /= 0) then
            write (stderr, *) "test_cache_history_projections failed: Incorrect "// &
                "dimensions of directions for empty history."
            test_cache_history_projections = .false.
        end if
        deallocate(dirs)

    end function test_cache_history_projections

    logical(c_bool) function test_cache_history_projections_channel() bind(C)
        !
        ! this function tests the open-shell routine which caches the packed
        ! history-projection directions of a single spin channel, embedding every
        ! history entry into that channel alone before projecting and packing it
        !
        use otr_arh, only: cache_history_projections_channel
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao_unit_tests, only: ref_project_asymm, ref_pack_asymm, &
                                      generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_list = 2, n_electrons = 1

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    v_diff(n_ao, n_ao, n_particle, n_list)
        real(rp), allocatable :: dirs(:, :), expected(:, :)
        integer(ip) :: channel, j, k

        ! assume tests pass
        test_cache_history_projections_channel = .true.

        ! generate a random density matrix and random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)
            do k = 1, n_list
                v_diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! every channel has to reproduce the history entry of that channel alone
        allocate(expected(n_param, n_list), dirs(n_param, n_list))
        do channel = 1, n_particle
            do k = 1, n_list
                expected(:, k) = ref_pack_asymm(ref_project_asymm( &
                    embed_channel(v_diff(:, :, channel, k), channel, n_ao, &
                                  n_particle), dm_oao), n_param)
            end do
            deallocate(dirs)
            call cache_history_projections_channel(v_diff, channel, dm_oao, n_list, &
                                                   n_param, n_particle, dirs)
            if (size(dirs, 1) /= n_param .or. size(dirs, 2) /= n_list) then
                write (stderr, *) "test_cache_history_projections_channel failed: "// &
                    "Incorrect dimensions of directions."
                test_cache_history_projections_channel = .false.
            else if (norm2(dirs - expected) > tol) then
                write (stderr, *) "test_cache_history_projections_channel failed: "// &
                    "Incorrect directions."
                test_cache_history_projections_channel = .false.
            end if
        end do
        deallocate(dirs, expected)

    end function test_cache_history_projections_channel

    logical(c_bool) function test_cache_channel_split_projections() bind(C)
        !
        ! this function tests the open-shell history-projection caching routine which
        ! embeds a history entry into one spin channel before projecting and packing
        ! it, keeping the two channels as separate columns
        !
        use otr_arh, only: cache_channel_split_projections
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao_unit_tests, only: ref_project_asymm, ref_pack_asymm, &
                                      generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_list = 2, n_electrons = 1

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    v_same(n_ao, n_ao, n_particle, n_list)
        real(rp), allocatable :: u(:, :), expected(:, :)
        integer(ip) :: j, k

        ! assume test passes
        test_cache_channel_split_projections = .true.

        ! generate a random density matrix and random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)
            do k = 1, n_list
                v_same(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! expected projections
        allocate(expected(n_param, n_particle * n_list))
        do k = 1, n_list
            do j = 1, n_particle
                expected(:, (j - 1)*n_list + k) = ref_pack_asymm(ref_project_asymm( &
                    embed_channel(v_same(:, :, j, k), j, n_ao, n_particle), dm_oao), &
                    n_param)
            end do
        end do

        ! call routine and determine if dimensions and values of the resulting
        ! projections match
        call cache_channel_split_projections(v_same, dm_oao, n_list, n_param, &
                                             n_particle, u)
        if (size(u, 1) /= n_param .or. size(u, 2) /= n_particle * n_list) then
            write (stderr, *) "test_cache_channel_split_projections failed: "// &
                "Incorrect dimensions of split projections."
            test_cache_channel_split_projections = .false.
        else if (norm2(u - expected) > tol) then
            write (stderr, *) "test_cache_channel_split_projections failed: "// &
                "Incorrect split projections."
            test_cache_channel_split_projections = .false.
        end if
        deallocate(u, expected)

    end function test_cache_channel_split_projections

    logical(c_bool) function test_cache_combined_channel_projections() bind(C)
        !
        ! this function tests the open-shell history-projection caching routine which
        ! embeds a history entry into one spin channel before projecting and packing
        ! it, summing a same-spin channel with the opposite-spin channel of the other
        ! spin
        !
        use otr_arh, only: cache_combined_channel_projections
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao_unit_tests, only: ref_project_asymm, ref_pack_asymm, &
                                      generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_list = 2, n_electrons = 1

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    v_same(n_ao, n_ao, n_particle, n_list), &
                    v_opp(n_ao, n_ao, n_particle, n_list)
        real(rp), allocatable :: u(:, :), expected(:, :)
        integer(ip) :: j, k

        ! assume test passes
        test_cache_combined_channel_projections = .true.

        ! generate a random density matrix and random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)
            do k = 1, n_list
                v_same(:, :, j, k) = generate_random_symm_matrix(n_ao)
                v_opp(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! expected projections
        allocate(expected(n_param, n_particle * n_list))
        do k = 1, n_list
            expected(:, k) = ref_pack_asymm(ref_project_asymm( &
                embed_channel(v_same(:, :, 1, k), 1_ip, n_ao, n_particle) + &
                embed_channel(v_opp(:, :, 2, k), 2_ip, n_ao, n_particle), dm_oao), &
                n_param)
            expected(:, n_list + k) = ref_pack_asymm(ref_project_asymm( &
                embed_channel(v_same(:, :, 2, k), 2_ip, n_ao, n_particle) + &
                embed_channel(v_opp(:, :, 1, k), 1_ip, n_ao, n_particle), dm_oao), &
                n_param)
        end do

        ! call routine and determine if dimensions and values of the resulting
        ! projections match
        call cache_combined_channel_projections(v_same, v_opp, dm_oao, n_list, &
                                                n_param, n_particle, u)
        if (size(u, 1) /= n_param .or. size(u, 2) /= n_particle * n_list) then
            write (stderr, *) "test_cache_combined_channel_projections failed: "// &
                "Incorrect dimensions of combined projections."
            test_cache_combined_channel_projections = .false.
        else if (norm2(u - expected) > tol) then
            write (stderr, *) "test_cache_combined_channel_projections failed: "// &
                "Incorrect combined projections."
            test_cache_combined_channel_projections = .false.
        end if
        deallocate(u, expected)

    end function test_cache_combined_channel_projections

    logical(c_bool) function test_get_low_rank_hess_factors() bind(C)
        !
        ! this function tests the subroutine which assembles the low-rank part of the
        ! approximate Hessian
        !
        use otr_arh, only: get_low_rank_hess_factors, arh_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_unit_tests, only: generate_random_symm_matrix

        integer(ip), parameter :: n_param = 4, n_diff = 2

        integer(ip), target :: n_param_target, n_particle_target
        real(rp) :: density(n_param, n_diff), potential(n_param, n_diff), &
                    linear(n_param, n_diff), nonlinear(n_param, n_diff), &
                    metric_inv(n_diff, n_diff), a_sym(n_diff, n_diff), &
                    a_inv(n_diff, n_diff), a_inv_comb(n_diff, n_diff), &
                    metric_weighted_a_sym(n_diff, n_diff)

        ! assume tests pass
        test_get_low_rank_hess_factors = .true.

        ! generate the packed history directions and the small dense matrices the
        ! coupling matrices are assembled from
        call random_number(density)
        call random_number(potential)
        call random_number(linear)
        call random_number(nonlinear)
        metric_inv = generate_random_symm_matrix(n_diff)
        a_sym = generate_random_symm_matrix(n_diff)
        a_inv = generate_random_symm_matrix(n_diff)
        a_inv_comb = generate_random_symm_matrix(n_diff)

        ! independently reconstruct the metric-weighted curvature
        metric_weighted_a_sym = matmul(metric_inv, matmul(a_sym, metric_inv))

        ! set up the ARH object with the cached quantities the assembly requires
        allocate(arh_object)
        call setup_settings(arh_object%settings)
        n_param_target = n_param
        arh_object%n_param => n_param_target
        n_particle_target = 1
        arh_object%n_particle => n_particle_target
        arh_object%dm_dirs = density
        arh_object%potential_dirs = potential
        arh_object%linear_potential_dirs = linear
        arh_object%nonlinear_potential_dirs = nonlinear
        arh_object%metric_inv = metric_inv
        arh_object%a_sym = a_sym
        arh_object%a_inv = a_inv
        arh_object%a_inv_comb = a_inv_comb

        ! multisecant SR1 stacks the linear and non-linear directions and couples
        ! each block through its own separately regularized system
        arh_object%settings%arh_type = "ms_sr1"
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%expansion_dirs(:, :n_diff) - linear) > tol) .or. &
            any(abs(arh_object%expansion_dirs(:, n_diff + 1:) - nonlinear) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "expansion directions for multisecant SR1."
            test_get_low_rank_hess_factors = .false.
        end if
        if (any(abs(arh_object%coupling_matrix(:n_diff, :n_diff) - 8.0_rp * a_inv) > &
                tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, n_diff + 1:) - 8.0_rp * &
                    a_inv_comb) > tol) .or. &
            any(abs(arh_object%coupling_matrix(:n_diff, n_diff + 1:)) > tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, :n_diff)) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "coupling matrix for multisecant SR1."
            test_get_low_rank_hess_factors = .false.
        end if
        if (any(abs(arh_object%projection_dirs - arh_object%expansion_dirs) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Expansion "// &
                "and projection directions do not coincide for multisecant SR1."
            test_get_low_rank_hess_factors = .false.
        end if

        ! subspace-projected multisecant expands in and contracts against the density
        ! difference history alone
        arh_object%settings%arh_type = "ms_sp"
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%expansion_dirs - density) > tol) .or. &
            any(abs(arh_object%projection_dirs - density) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "directions for subspace-projected multisecant."
            test_get_low_rank_hess_factors = .false.
        end if
        if (any(abs(arh_object%coupling_matrix - 8.0_rp * metric_weighted_a_sym) > &
                tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "coupling matrix for subspace-projected multisecant."
            test_get_low_rank_hess_factors = .false.
        end if

        ! symmetrized ARH couples the density and potential difference histories in
        ! both directions, leaving the diagonal blocks empty
        arh_object%settings%arh_type = "symm_arh"
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%expansion_dirs(:, :n_diff) - density) > tol) .or. &
            any(abs(arh_object%expansion_dirs(:, n_diff + 1:) - potential) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "expansion directions for symmetrized ARH."
            test_get_low_rank_hess_factors = .false.
        end if
        if (any(abs(arh_object%coupling_matrix(:n_diff, n_diff + 1:) - 4.0_rp * &
                    metric_inv) > tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, :n_diff) - 4.0_rp * &
                    metric_inv) > tol) .or. &
            any(abs(arh_object%coupling_matrix(:n_diff, :n_diff)) > tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, n_diff + 1:)) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "coupling matrix for symmetrized ARH."
            test_get_low_rank_hess_factors = .false.
        end if

        ! multisecant PSB adds a density-density block subtracting the doubly counted
        ! curvature to the symmetrized ARH coupling
        arh_object%settings%arh_type = "ms_psb"
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%coupling_matrix(:n_diff, :n_diff) + 8.0_rp * &
                    metric_weighted_a_sym) > tol) .or. &
            any(abs(arh_object%coupling_matrix(:n_diff, n_diff + 1:) - 8.0_rp * &
                    metric_inv) > tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, :n_diff) - 8.0_rp * &
                    metric_inv) > tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, n_diff + 1:)) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "coupling matrix for multisecant PSB."
            test_get_low_rank_hess_factors = .false.
        end if

        ! standard ARH is the only type whose expansion and projection directions
        ! differ, expanding in the potential and contracting against the density
        ! difference history
        arh_object%settings%arh_type = "arh"
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%expansion_dirs - potential) > tol) .or. &
            any(abs(arh_object%projection_dirs - density) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "directions for standard ARH."
            test_get_low_rank_hess_factors = .false.
        end if
        if (any(abs(arh_object%coupling_matrix - 8.0_rp * metric_inv) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "coupling matrix for standard ARH."
            test_get_low_rank_hess_factors = .false.
        end if

        ! test the open-shell coupling matrix which is exactly half its closed-shell 
        ! counterpart,
        n_particle_target = 2
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%coupling_matrix - 4.0_rp * metric_inv) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Open-shell "// &
                "coupling matrix is not half the closed-shell one."
            test_get_low_rank_hess_factors = .false.
        end if

        ! an empty history has to leave every factor unallocated
        n_particle_target = 1
        deallocate(arh_object%dm_dirs)
        allocate(arh_object%dm_dirs(n_param, 0))
        call get_low_rank_hess_factors()
        if (allocated(arh_object%expansion_dirs) .or. &
            allocated(arh_object%projection_dirs) .or. &
            allocated(arh_object%coupling_matrix)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Factors were "// &
                "assembled for an empty history."
            test_get_low_rank_hess_factors = .false.
        end if

        ! deallocate ARH object
        deallocate(arh_object)

    end function test_get_low_rank_hess_factors

    logical(c_bool) function test_build_a_sym_cs() bind(C)
        !
        ! this function tests the function which builds the weighted-symmetrized 
        ! A = S^T Y matrix
        !
        use otr_arh, only: build_a_sym_cs

        integer(ip), parameter :: n_ao = 2, n_particle = 1, n_diff = 2

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_diff(n_ao, n_ao, n_particle, n_diff), &
                    expected_a_sym(n_diff, n_diff)
        real(rp), allocatable :: a_sym(:, :)

        ! assume tests pass
        test_build_a_sym_cs = .true.

        ! initialize density matrix and potential differences and calculate expected A 
        ! matrix
        dm_diff(:, :, 1, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 1, 2) = reshape([0.0_rp, 0.0_rp, &
                                       0.0_rp, 2.0_rp], [n_ao, n_ao])
        v_diff(:, :, 1, 1) = reshape([10.0_rp, 0.0_rp, &
                                      0.0_rp, -0.3_rp], [n_ao, n_ao])
        v_diff(:, :, 1, 2) = reshape([0.6_rp, 0.0_rp, &
                                      0.0_rp, 0.1_rp], [n_ao, n_ao])
        expected_a_sym = reshape([10.0_rp, -0.2_rp, &
                                  -0.2_rp, 0.2_rp], [n_diff, n_diff])

        ! generate A matrix and verify
        a_sym = build_a_sym_cs(dm_diff, v_diff, n_ao)
        if (norm2(a_sym - expected_a_sym) > tol) then
            write (stderr, *) "test_build_a_sym_cs failed: Incorrect symmetrized "// &
                "matrix."
            test_build_a_sym_cs = .false.
        end if

    end function test_build_a_sym_cs

    logical(c_bool) function test_build_a_block_sym_os() bind(C)
        !
        ! this function tests the function which builds the dense,
        ! cross-channel-symmetrized open-shell A matrix
        !
        use otr_arh, only: build_a_block_sym_os
        use otr_oao_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 2, n_diff = 1

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_same_eff(n_ao, n_ao, n_particle, n_diff), &
                    v_opp(n_ao, n_ao, n_particle, n_diff), &
                    expected_a_block(2 * n_diff, 2 * n_diff)
        real(rp), allocatable :: a_block(:, :)

        ! assume tests pass
        test_build_a_block_sym_os = .true.

        ! initialize density matrix and potential differences and calculate expected A 
        ! matrix
        dm_diff(:, :, 1, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 1) = reshape([0.0_rp, 0.0_rp, &
                                       0.0_rp, 1.0_rp], [n_ao, n_ao])
        v_same_eff(:, :, 1, 1) = reshape([2.0_rp, 0.0_rp, &
                                          0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_same_eff(:, :, 2, 1) = reshape([0.0_rp, 0.0_rp, &
                                          0.0_rp, 4.0_rp], [n_ao, n_ao])
        v_opp(:, :, 1, 1) = reshape([3.0_rp, 0.0_rp, &
                                     0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_opp(:, :, 2, 1) = reshape([0.0_rp, 0.0_rp, &
                                     0.0_rp, 5.0_rp], [n_ao, n_ao])
        expected_a_block = reshape([2.0_rp, 4.0_rp, &
                                    4.0_rp, 4.0_rp], [2 * n_diff, 2 * n_diff])

        ! generate A matrix and verify
        a_block = build_a_block_sym_os(dm_diff, v_same_eff, v_opp, n_ao)
        if (norm2(a_block - expected_a_block) > tol) then
            write (stderr, *) "test_build_a_block_sym_os failed: Incorrect "// &
                "cross-channel-symmetrized A matrix."
            test_build_a_block_sym_os = .false.
        end if

    end function test_build_a_block_sym_os

    logical(c_bool) function test_symmetrize_weighted() bind(C)
        !
        ! this function tests the subroutine which performs a weighted symmetrization
        ! of a square matrix
        !
        use otr_arh, only: symmetrize_weighted

        real(rp) :: a(2, 2), expected(2, 2), step_norms(2)

        ! assume tests pass
        test_symmetrize_weighted = .true.

        ! initialize matrix with an antisymmetric contribution and step norms which
        ! bias the weight towards the element associated with the larger step norm
        a = reshape([1.0_rp, 3.0_rp, &
                     2.0_rp, 4.0_rp], [2, 2])
        step_norms = [1.0_rp, 3.0_rp]

        ! initialize expected matrix, where the off-diagonal elements are blended with
        ! a weight of 3 / (1 + 3) on the element of the larger step norm
        expected = reshape([1.0_rp, 2.75_rp, &
                            2.75_rp, 4.0_rp], [2, 2])

        ! call routine and determine if values of resulting matrix match
        call symmetrize_weighted(a, step_norms)
        if (norm2(a - expected) > tol) then
            write (stderr, *) "test_symmetrize_weighted failed: Incorrect matrix "// &
                "values after weighted symmetrization."
            test_symmetrize_weighted = .false.
        end if

        ! initialize matrix and vanishing step norms
        a = reshape([1.0_rp, 3.0_rp, &
                     2.0_rp, 4.0_rp], [2, 2])
        step_norms = [0.0_rp, 0.0_rp]

        ! initialize expected matrix, where the off-diagonal elements are averaged
        expected = reshape([1.0_rp, 2.5_rp, &
                            2.5_rp, 4.0_rp], [2, 2])

        ! call routine and determine if values of resulting matrix match
        call symmetrize_weighted(a, step_norms)
        if (norm2(a - expected) > tol) then
            write (stderr, *) "test_symmetrize_weighted failed: Incorrect matrix "// &
                "values after weighted symmetrization for vanishing step norms."
            test_symmetrize_weighted = .false.
        end if

        ! initialize symmetric matrix and step norms
        a = reshape([1.0_rp, 2.0_rp, &
                     2.0_rp, 4.0_rp], [2, 2])
        expected = a
        step_norms = [1.0_rp, 3.0_rp]

        ! call routine and determine if an already symmetric matrix is left unchanged
        call symmetrize_weighted(a, step_norms)
        if (norm2(a - expected) > tol) then
            write (stderr, *) "test_symmetrize_weighted failed: Symmetric matrix "// &
                "not left unchanged by weighted symmetrization."
            test_symmetrize_weighted = .false.
        end if

    end function test_symmetrize_weighted

    logical(c_bool) function test_cross_symmetrize_weighted() bind(C)
        !
        ! this function tests the subroutine which performs a weighted symmetrization
        ! between two related off-diagonal blocks of a larger matrix
        !
        use otr_arh, only: cross_symmetrize_weighted

        real(rp) :: a12(2, 2), a21(2, 2), expected12(2, 2), expected21(2, 2), &
                    step_norms1(2), step_norms2(2)

        ! assume tests pass
        test_cross_symmetrize_weighted = .true.

        ! initialize blocks which are not transposes of each other and step norms which
        ! bias the weight towards the element associated with the larger step norm
        a12 = reshape([1.0_rp, 3.0_rp, &
                       2.0_rp, 4.0_rp], [2, 2])
        a21 = reshape([5.0_rp, 7.0_rp, &
                       6.0_rp, 8.0_rp], [2, 2])
        step_norms1 = [1.0_rp, 1.0_rp]
        step_norms2 = [1.0_rp, 3.0_rp]

        ! initialize expected blocks, where the elements of the second column are
        ! blended with a weight of 3 / (1 + 3) on the element of the larger step norm
        ! while the elements of the first column are averaged
        expected12 = reshape([3.0_rp, 4.5_rp, &
                              5.75_rp, 7.0_rp], [2, 2])
        expected21 = transpose(expected12)

        ! call routine and determine if values of resulting blocks match and if these
        ! are transposes of each other
        call cross_symmetrize_weighted(a12, a21, step_norms1, step_norms2)
        if (norm2(a12 - expected12) > tol) then
            write (stderr, *) "test_cross_symmetrize_weighted failed: Incorrect "// &
                "first block values after cross-symmetrization."
            test_cross_symmetrize_weighted = .false.
        end if
        if (norm2(a21 - expected21) > tol) then
            write (stderr, *) "test_cross_symmetrize_weighted failed: Incorrect "// &
                "second block values after cross-symmetrization."
            test_cross_symmetrize_weighted = .false.
        end if
        if (norm2(a21 - transpose(a12)) > tol) then
            write (stderr, *) "test_cross_symmetrize_weighted failed: Blocks are "// &
                "not transposes of each other after cross-symmetrization."
            test_cross_symmetrize_weighted = .false.
        end if

        ! initialize blocks and vanishing step norms
        a12 = reshape([1.0_rp, 3.0_rp, &
                       2.0_rp, 4.0_rp], [2, 2])
        a21 = reshape([5.0_rp, 7.0_rp, &
                       6.0_rp, 8.0_rp], [2, 2])
        step_norms1 = [0.0_rp, 0.0_rp]
        step_norms2 = [0.0_rp, 0.0_rp]

        ! initialize expected blocks, where all elements are averaged
        expected12 = reshape([3.0_rp, 4.5_rp, &
                              4.5_rp, 6.0_rp], [2, 2])
        expected21 = transpose(expected12)

        ! call routine and determine if values of resulting blocks match
        call cross_symmetrize_weighted(a12, a21, step_norms1, step_norms2)
        if (norm2(a12 - expected12) > tol) then
            write (stderr, *) "test_cross_symmetrize_weighted failed: Incorrect "// &
                "first block values after cross-symmetrization for vanishing step "// &
                "norms."
            test_cross_symmetrize_weighted = .false.
        end if
        if (norm2(a21 - expected21) > tol) then
            write (stderr, *) "test_cross_symmetrize_weighted failed: Incorrect "// &
                "second block values after cross-symmetrization for vanishing step "// &
                "norms."
            test_cross_symmetrize_weighted = .false.
        end if

    end function test_cross_symmetrize_weighted

    logical(c_bool) function test_symmetrize_exact() bind(C)
        !
        ! this function tests the subroutine which performs a plain, unweighted
        ! symmetrization of a square matrix
        !
        use otr_arh, only: symmetrize_exact

        real(rp) :: a(3, 3), expected(3, 3)

        ! assume tests pass
        test_symmetrize_exact = .true.

        ! initialize matrix with an antisymmetric contribution
        a = reshape([1.0_rp, 4.0_rp, 5.0_rp, &
                     2.0_rp, 6.0_rp, 9.0_rp, &
                     3.0_rp, 7.0_rp, 8.0_rp], [3, 3])

        ! initialize expected matrix, which averages each pair of off-diagonal elements
        expected = reshape([1.0_rp, 3.0_rp, 4.0_rp, &
                            3.0_rp, 6.0_rp, 8.0_rp, &
                            4.0_rp, 8.0_rp, 8.0_rp], [3, 3])

        ! call routine and determine if values of resulting matrix match
        call symmetrize_exact(a)
        if (norm2(a - expected) > tol) then
            write (stderr, *) "test_symmetrize_exact failed: Incorrect matrix "// &
                "values after symmetrization."
            test_symmetrize_exact = .false.
        end if

    end function test_symmetrize_exact

    logical(c_bool) function test_noise_threshold() bind(C)
        !
        ! this function tests the function which estimates the numerical noise floor of
        ! a matrix that is exactly symmetric in exact arithmetic
        !
        use otr_arh, only: noise_threshold, arh_settings_type, symm_noise_safety_margin
        use opentrustregion, only: numerical_zero
        use opentrustregion_unit_tests, only: setup_settings

        real(rp) :: a(3, 3), thresh
        integer(ip) :: error
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_noise_threshold = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize matrix whose antisymmetric part has a single non-vanishing pair of
        ! elements of magnitude 0.3 and therefore a largest singular value of 0.3
        a = reshape([1.0_rp, 1.7_rp, 3.0_rp, &
                     2.3_rp, 4.0_rp, 5.0_rp, &
                     3.0_rp, 5.0_rp, 6.0_rp], [3, 3])

        ! call routine and determine if the threshold is the safety margin times the
        ! measured asymmetry
        thresh = noise_threshold(a, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_noise_threshold failed: Produced error."
            test_noise_threshold = .false.
        end if
        if (abs(thresh - symm_noise_safety_margin * 0.3_rp) > tol) then
            write (stderr, *) "test_noise_threshold failed: Incorrect threshold "// &
                "for asymmetric matrix."
            test_noise_threshold = .false.
        end if

        ! initialize exactly symmetric matrix
        a = reshape([1.0_rp, 2.0_rp, 3.0_rp, &
                     2.0_rp, 4.0_rp, 5.0_rp, &
                     3.0_rp, 5.0_rp, 6.0_rp], [3, 3])

        ! call routine and determine if the threshold falls back to the numerical zero
        ! floor instead of vanishing
        thresh = noise_threshold(a, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_noise_threshold failed: Produced error for "// &
                "symmetric matrix."
            test_noise_threshold = .false.
        end if
        if (abs(thresh - symm_noise_safety_margin * numerical_zero) > &
            tol * numerical_zero) then
            write (stderr, *) "test_noise_threshold failed: Incorrect threshold "// &
                "for symmetric matrix."
            test_noise_threshold = .false.
        end if

    end function test_noise_threshold

    logical(c_bool) function test_regularized_eigval_inv() bind(C)
        !
        ! this function tests the function which returns Tikhonov-regularized
        ! pseudoinverse eigenvalues
        !
        use otr_arh, only: regularized_eigval_inv

        real(rp) :: eig_vals(4), eig_vals_inv(4), expected(4)

        ! assume tests pass
        test_regularized_eigval_inv = .true.

        ! initialize eigenvalues spanning eigenvalues far above, at the order of and
        ! below the threshold
        eig_vals = [2.0_rp, -0.5_rp, 0.1_rp, 0.0_rp]

        ! initialize expected inverted eigenvalues, where eigenvalues far above the
        ! threshold are essentially inverted exactly, eigenvalues at the order of the
        ! threshold are damped and vanishing eigenvalues do not diverge
        expected = [0.4987531172069825_rp, -1.9230769230769231_rp, &
                    4.9999999999999999_rp, 0.0_rp]

        ! call routine and determine if values of resulting eigenvalues match
        eig_vals_inv = regularized_eigval_inv(eig_vals, 0.1_rp)
        if (norm2(eig_vals_inv - expected) > tol) then
            write (stderr, *) "test_regularized_eigval_inv failed: Incorrect "// &
                "regularized inverse eigenvalues."
            test_regularized_eigval_inv = .false.
        end if

    end function test_regularized_eigval_inv

    logical(c_bool) function test_truncated_eigval_inv() bind(C)
        !
        ! this function tests the function which returns hard-truncated pseudoinverse
        ! eigenvalues
        !
        use otr_arh, only: truncated_eigval_inv

        real(rp) :: eig_vals(5), eig_vals_inv(5), expected(5)

        ! assume tests pass
        test_truncated_eigval_inv = .true.

        ! initialize eigenvalues spanning eigenvalues above, at and below the threshold
        eig_vals = [2.0_rp, -0.5_rp, 0.1_rp, 0.05_rp, 0.0_rp]

        ! initialize expected inverted eigenvalues, where eigenvalues above the
        ! threshold are inverted exactly while eigenvalues at or below the threshold
        ! are discarded
        expected = [0.5_rp, -2.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]

        ! call routine and determine if values of resulting eigenvalues match
        eig_vals_inv = truncated_eigval_inv(eig_vals, 0.1_rp)
        if (norm2(eig_vals_inv - expected) > tol) then
            write (stderr, *) "test_truncated_eigval_inv failed: Incorrect "// &
                "truncated inverse eigenvalues."
            test_truncated_eigval_inv = .false.
        end if

    end function test_truncated_eigval_inv

    logical(c_bool) function test_get_arh_metric_inv() bind(C)
        !
        ! this function tests the subroutine which calculates the block-diagonal
        ! pseudoinverse of the ARH metric
        !
        use otr_arh, only: get_arh_metric_inv
        use otr_oao_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 2, n_diff = 3

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    dm_diff_empty(n_ao, n_ao, n_particle, 0), &
                    expected(n_particle * n_diff, n_particle * n_diff)
        real(rp), allocatable :: metric_inv(:, :)

        ! assume tests pass
        test_get_arh_metric_inv = .true.

        ! initialize density matrix differences, where the second difference of the
        ! first particle repeats the first one and is therefore linearly dependent,
        ! while the second particle is linearly independent but not orthogonal
        dm_diff(:, :, 1, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 1, 2) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 1, 3) = reshape([0.0_rp, 1.0_rp, &
                                       1.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 2) = reshape([1.0_rp, 1.0_rp, &
                                       1.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 3) = reshape([0.0_rp, 0.0_rp, &
                                       0.0_rp, 1.0_rp], [n_ao, n_ao])

        ! initialize the expected block-diagonal pseudoinverse: the first particle
        ! has the metric [[1, 0], [0, 2]] over its accepted differences 1 and 3, so
        ! its inverse lands on the diagonal entries 1 and 3 while the rejected
        ! difference 2 carries a vanishing row and column, and the second particle
        ! has the metric [[1, 1, 0], [1, 3, 0], [0, 0, 1]] of full rank
        expected = 0.0_rp
        expected(1, 1) = 1.0_rp
        expected(3, 3) = 0.5_rp
        expected(4:5, 4:5) = reshape([1.5_rp, -0.5_rp, &
                                      -0.5_rp, 0.5_rp], [2, 2])
        expected(6, 6) = 1.0_rp

        ! call routine and determine if dimensions and values of the resulting
        ! pseudoinverse match
        call get_arh_metric_inv(dm_diff, metric_inv)
        if (size(metric_inv, 1) /= n_particle * n_diff .or. &
            size(metric_inv, 2) /= n_particle * n_diff) then
            write (stderr, *) "test_get_arh_metric_inv failed: Incorrect "// &
                "dimensions of metric pseudoinverse."
            test_get_arh_metric_inv = .false.
            return
        end if
        if (norm2(metric_inv - expected) > tol) then
            write (stderr, *) "test_get_arh_metric_inv failed: Incorrect metric "// &
                "pseudoinverse."
            test_get_arh_metric_inv = .false.
        end if
        deallocate(metric_inv)

        ! call routine for an empty history and determine if an empty pseudoinverse
        ! is returned
        call get_arh_metric_inv(dm_diff_empty, metric_inv)
        if (size(metric_inv, 1) /= 0 .or. size(metric_inv, 2) /= 0) then
            write (stderr, *) "test_get_arh_metric_inv failed: Incorrect "// &
                "dimensions of metric pseudoinverse for empty history."
            test_get_arh_metric_inv = .false.
        end if
        deallocate(metric_inv)

    end function test_get_arh_metric_inv

    logical(c_bool) function test_get_ms_a_inv_cs() bind(C)
        !
        ! this function tests the subroutine which computes the pseudoinverse
        ! multisecant SR1 matrix for the closed-shell case
        !
        use otr_arh, only: get_ms_a_inv_cs, arh_settings_type
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_ao = 2, n_particle = 1, n_diff = 2

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    fock_diff(n_ao, n_ao, n_particle, n_diff), &
                    empty_dm_diff(n_ao, n_ao, n_particle, 0), &
                    empty_fock_diff(n_ao, n_ao, n_particle, 0), &
                    expected_a_inv(n_diff, n_diff)
        real(rp), allocatable :: a_inv(:, :)
        integer(ip) :: error
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_cs = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize orthogonal density matrix differences whose step norms differ, so
        ! that the weighted symmetrization of the non-linear part is biased rather than
        ! a plain average, and Fock matrix differences which produce a multisecant SR1
        ! matrix with the symmetric part 10 and 0.2 on the diagonal and an
        ! antisymmetric part of magnitude 0.6
        dm_diff(:, :, 1, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 1, 2) = reshape([0.0_rp, 0.0_rp, &
                                       0.0_rp, 2.0_rp], [n_ao, n_ao])
        fock_diff(:, :, 1, 1) = reshape([10.0_rp, 0.0_rp, &
                                         0.0_rp, -0.3_rp], [n_ao, n_ao])
        fock_diff(:, :, 1, 2) = reshape([0.6_rp, 0.0_rp, &
                                         0.0_rp, 0.1_rp], [n_ao, n_ao])

        ! initialize expected pseudoinverse for the linear part, where the exact
        ! symmetrization removes the off-diagonal elements and the measured asymmetry
        ! produces a threshold of 3.0 which discards the eigenvalue 0.2
        expected_a_inv = reshape([0.1_rp, 0.0_rp, &
                                  0.0_rp, 0.0_rp], [n_diff, n_diff])

        ! call routine for the linear part and determine if values of resulting
        ! pseudoinverse match
        call get_ms_a_inv_cs(dm_diff, fock_diff, .true., a_inv, n_ao, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Produced error for "// &
                "linear part."
            test_get_ms_a_inv_cs = .false.
            return
        end if
        if (norm2(a_inv - expected_a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Incorrect "// &
                "pseudoinverse for linear part."
            test_get_ms_a_inv_cs = .false.
        end if

        ! initialize expected pseudoinverse for the non-linear part, where the step
        ! norms 1 and 2 blend the off-diagonal elements to -0.2 with a weight of
        ! 2 / (1 + 2) on the element of the larger step norm, and both eigenvalues are
        ! retained and only marginally damped, so that the pseudoinverse is the exact
        ! inverse of the symmetrized matrix
        expected_a_inv = reshape([0.10204081632653061_rp, 0.10204081632653061_rp, &
                                  0.10204081632653061_rp, 5.1020408163265306_rp], &
                                 [n_diff, n_diff])

        ! call routine for the non-linear part and determine if values of resulting
        ! pseudoinverse match
        call get_ms_a_inv_cs(dm_diff, fock_diff, .false., a_inv, n_ao, settings, &
                             error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Produced error for "// &
                "non-linear part."
            test_get_ms_a_inv_cs = .false.
            return
        end if
        if (norm2(a_inv - expected_a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Incorrect "// &
                "pseudoinverse for non-linear part."
            test_get_ms_a_inv_cs = .false.
        end if

        ! call routine for an empty history and determine if dimensions of the
        ! resulting pseudoinverse vanish
        deallocate(a_inv)
        call get_ms_a_inv_cs(empty_dm_diff, empty_fock_diff, .true., a_inv, n_ao, &
                             settings, error)
        if (size(a_inv, 1) /= 0 .or. size(a_inv, 2) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Incorrect "// &
                "pseudoinverse dimensions for empty history."
            test_get_ms_a_inv_cs = .false.
        end if

        ! deallocate multisecant SR1 quantities
        deallocate(a_inv)

    end function test_get_ms_a_inv_cs

    logical(c_bool) function test_get_ms_a_inv_os_linear() bind(C)
        !
        ! this function tests the subroutine which computes the pseudoinverse
        ! multisecant SR1 matrix in a spin-separated manner for the linear part in the 
        ! open-shell case
        !
        use otr_arh, only: get_ms_a_inv_os_linear, arh_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 2, n_diff = 1

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_same_spin_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_opposite_spin_diff(n_ao, n_ao, n_particle, n_diff), &
                    empty_dm_diff(n_ao, n_ao, n_particle, 0), &
                    empty_v_diff(n_ao, n_ao, n_particle, 0), &
                    expected_a_inv(2 * n_diff, 2 * n_diff)
        real(rp), allocatable :: a_inv(:, :)
        integer(ip) :: error
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_os_linear = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize normalized density matrix differences and spin-resolved potential
        ! matrix differences which produce a spin-combined multisecant SR1 matrix with
        ! the symmetric part 10 and 0.2 on the diagonal and an antisymmetric part of
        ! magnitude 0.1
        dm_diff(:, :, 1, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_same_spin_diff(:, :, 1, 1) = reshape([10.0_rp, 0.0_rp, &
                                                0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_same_spin_diff(:, :, 2, 1) = reshape([0.2_rp, 0.0_rp, &
                                                0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_opposite_spin_diff(:, :, 1, 1) = reshape([0.1_rp, 0.0_rp, &
                                                    0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_opposite_spin_diff(:, :, 2, 1) = reshape([-0.1_rp, 0.0_rp, &
                                                    0.0_rp, 0.0_rp], [n_ao, n_ao])

        ! initialize expected pseudoinverse, where the measured asymmetry produces a
        ! threshold of 0.5 which discards the eigenvalue 0.2
        expected_a_inv = reshape([0.1_rp, 0.0_rp, &
                                  0.0_rp, 0.0_rp], [2 * n_diff, 2 * n_diff])

        ! call routine and determine if values of resulting pseudoinverse match
        call get_ms_a_inv_os_linear(dm_diff, v_same_spin_diff, v_opposite_spin_diff, &
                                    a_inv, n_ao, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Produced error."
            test_get_ms_a_inv_os_linear = .false.
            return
        end if
        if (norm2(a_inv - expected_a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "pseudoinverse."
            test_get_ms_a_inv_os_linear = .false.
        end if

        ! call routine for an empty history and determine if dimensions of the
        ! resulting pseudoinverse vanish
        deallocate(a_inv)
        call get_ms_a_inv_os_linear(empty_dm_diff, empty_v_diff, empty_v_diff, &
                                    a_inv, n_ao, settings, error)
        if (size(a_inv, 1) /= 0 .or. size(a_inv, 2) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "pseudoinverse dimensions for empty history."
            test_get_ms_a_inv_os_linear = .false.
        end if

        ! deallocate multisecant SR1 quantities
        deallocate(a_inv)

    end function test_get_ms_a_inv_os_linear

    logical(c_bool) function test_get_ms_a_inv_os_nonlinear() bind(C)
        !
        ! this function tests the subroutine which computes the pseudoinverse
        ! multisecant SR1 matrix in a spin-combined manner for the non-linear part in 
        ! the open-shell case
        !
        use otr_arh, only: get_ms_a_inv_os_nonlinear, arh_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 2, n_diff = 2

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_diff(n_ao, n_ao, n_particle, n_diff), &
                    empty_dm_diff(n_ao, n_ao, n_particle, 0), &
                    empty_v_diff(n_ao, n_ao, n_particle, 0), &
                    expected_a_inv(n_diff, n_diff)
        real(rp), allocatable :: a_inv(:, :)
        integer(ip) :: error
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_os_nonlinear = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize density matrix differences which are confined to a single spin
        ! channel each and whose step norms differ, so that the weighted symmetrization
        ! is biased rather than a plain average, and potential matrix differences which
        ! produce a spin-combined multisecant SR1 matrix with the symmetric part 10 and
        ! 0.2 on the diagonal and an antisymmetric part of magnitude 0.6
        dm_diff = 0.0_rp
        dm_diff(:, :, 1, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 2) = reshape([2.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_diff = 0.0_rp
        v_diff(:, :, 1, 1) = reshape([10.0_rp, 0.0_rp, &
                                      0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_diff(:, :, 2, 1) = reshape([-0.3_rp, 0.0_rp, &
                                      0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_diff(:, :, 1, 2) = reshape([0.6_rp, 0.0_rp, &
                                      0.0_rp, 0.0_rp], [n_ao, n_ao])
        v_diff(:, :, 2, 2) = reshape([0.1_rp, 0.0_rp, &
                                      0.0_rp, 0.0_rp], [n_ao, n_ao])

        ! initialize expected pseudoinverse, where the step norms 1 and 2 blend the
        ! off-diagonal elements to -0.2 with a weight of 2 / (1 + 2) on the element of
        ! the larger step norm, and both eigenvalues are retained and only marginally
        ! damped, so that the pseudoinverse is the exact inverse of the symmetrized
        ! matrix
        expected_a_inv = reshape([0.10204081632653061_rp, 0.10204081632653061_rp, &
                                  0.10204081632653061_rp, 5.1020408163265306_rp], &
                                 [n_diff, n_diff])

        ! call routine and determine if values of resulting pseudoinverse match
        call get_ms_a_inv_os_nonlinear(dm_diff, v_diff, a_inv, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Produced error."
            test_get_ms_a_inv_os_nonlinear = .false.
            return
        end if
        if (norm2(a_inv - expected_a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Incorrect "// &
                "pseudoinverse."
            test_get_ms_a_inv_os_nonlinear = .false.
        end if

        ! call routine for an empty history and determine if dimensions of the
        ! resulting pseudoinverse vanish
        deallocate(a_inv)
        call get_ms_a_inv_os_nonlinear(empty_dm_diff, empty_v_diff, a_inv, settings, &
                                       error)
        if (size(a_inv, 1) /= 0 .or. size(a_inv, 2) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Incorrect "// &
                "pseudoinverse dimensions for empty history."
            test_get_ms_a_inv_os_nonlinear = .false.
        end if

        ! deallocate multisecant SR1 quantities
        deallocate(a_inv)

    end function test_get_ms_a_inv_os_nonlinear

    logical(c_bool) function test_spectral_to_dense() bind(C)
        !
        ! this function tests the function which reconstructs a dense symmetric
        ! matrix from its eigenvectors and inverted eigenvalues
        !
        use otr_arh, only: spectral_to_dense

        integer(ip), parameter :: n = 2

        real(rp) :: eigvecs(n, n), inv_eigvals(n), expected(n, n)
        real(rp), allocatable :: mat(:, :)

        ! assume tests pass
        test_spectral_to_dense = .true.

        ! initialize an orthonormal but non-symmetric eigenvector matrix, so that a
        ! transposed reconstruction would be caught, and distinct inverted eigenvalues
        eigvecs = reshape([0.6_rp, 0.8_rp, &
                           -0.8_rp, 0.6_rp], [n, n])
        inv_eigvals = [0.25_rp, 4.0_rp]

        ! initialize the expected matrix 0.25 * v1 v1^T + 4 * v2 v2^T
        expected = reshape([2.65_rp, -1.8_rp, &
                            -1.8_rp, 1.6_rp], [n, n])

        ! call routine and determine if dimensions and values of the reconstructed
        ! matrix match
        mat = spectral_to_dense(eigvecs, inv_eigvals)
        if (size(mat, 1) /= n .or. size(mat, 2) /= n) then
            write (stderr, *) "test_spectral_to_dense failed: Incorrect dimensions "// &
                "of reconstructed matrix."
            test_spectral_to_dense = .false.
        else if (norm2(mat - expected) > tol) then
            write (stderr, *) "test_spectral_to_dense failed: Incorrect "// &
                "reconstructed matrix."
            test_spectral_to_dense = .false.
        end if

    end function test_spectral_to_dense

    logical(c_bool) function test_prepend() bind(C)
        !
        ! this function tests the subroutine which prepends an array to a list of arrays
        !
        use otr_arh, only: prepend

        real(rp), allocatable :: list(:, :, :, :)
        real(rp) :: new_array(2, 1, 1), expected(2, 1, 1, 3)

        ! assume tests pass
        test_prepend = .true.

        ! allocate empty list and initialize array to be prepended
        allocate(list(2, 1, 1, 0))
        new_array = reshape([1.0_rp, 2.0_rp], [2, 1, 1])

        ! prepend array to empty list and determine if dimensions and values of
        ! resulting list match
        call prepend(list, new_array)
        if (size(list, 4) /= 1) then
            write (stderr, *) "test_prepend failed: Incorrect list dimensions "// &
                "after prepending to empty list."
            test_prepend = .false.
        end if
        if (norm2(list(:, :, :, 1) - new_array) > tol) then
            write (stderr, *) "test_prepend failed: Incorrect list values after "// &
                "prepending to empty list."
            test_prepend = .false.
        end if

        ! initialize expected list after prepending two further arrays
        expected = reshape([5.0_rp, 6.0_rp, &
                            3.0_rp, 4.0_rp, &
                            1.0_rp, 2.0_rp], [2, 1, 1, 3])

        ! prepend two further arrays and determine if dimensions and values of
        ! resulting list match, so that new arrays are added at the front while the
        ! order of the existing arrays is retained
        call prepend(list, reshape([3.0_rp, 4.0_rp], [2, 1, 1]))
        call prepend(list, reshape([5.0_rp, 6.0_rp], [2, 1, 1]))
        if (size(list, 4) /= 3) then
            write (stderr, *) "test_prepend failed: Incorrect list dimensions "// &
                "after prepending to non-empty list."
            test_prepend = .false.
        end if
        if (norm2(list - expected) > tol) then
            write (stderr, *) "test_prepend failed: Incorrect list values after "// &
                "prepending to non-empty list."
            test_prepend = .false.
        end if

        ! deallocate list
        deallocate(list)

    end function test_prepend

end module otr_arh_unit_tests
