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

    function generate_random_upper_triangular(n) result(matrix)
        !
        ! this function generates a random upper triangular matrix whose diagonal is
        ! bounded away from zero, so that it is invertible and well enough conditioned
        ! to stand in for the Cholesky factor a history factorization would produce
        !
        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        integer(ip) :: i, j

        call random_number(matrix)
        do j = 1, n
            do i = j + 1, n
                matrix(i, j) = 0.0_rp
            end do
            matrix(j, j) = matrix(j, j) + 1.0_rp
        end do

    end function generate_random_upper_triangular

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

    subroutine ref_pivoted_cholesky(diff1, diff2, dm_oao, n_param, pivot_order, &
                                    full_dirs_check, chol_ref_check, info)
        !
        ! this subroutine independently reproduces, for a two-entry raw 
        ! history-difference pair, the magnitude-based pivot order the history 
        ! factorization would select, the resulting packed and projected columns in 
        ! that order, and the Cholesky factor of the pivoted pair's raw Gram matrix
        !
        use otr_oao_unit_tests, only: ref_project_asymm, ref_pack_asymm

        real(rp), intent(in) :: diff1(:, :, :), diff2(:, :, :), dm_oao(:, :, :)
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: pivot_order(2), info
        real(rp), intent(out) :: full_dirs_check(n_param, 2), chol_ref_check(2, 2)

        real(rp) :: flats(size(diff1), 2)
        integer(ip) :: j, k
        external :: dpotrf

        flats(:, 1) = reshape(diff1, [size(diff1)])
        flats(:, 2) = reshape(diff2, [size(diff2)])

        ! magnitude decides which entry is pivoted first
        if (dot_product(flats(:, 1), flats(:, 1)) >= &
            dot_product(flats(:, 2), flats(:, 2))) then
            pivot_order = [1, 2]
        else
            pivot_order = [2, 1]
        end if

        do k = 1, 2
            if (pivot_order(k) == 1) then
                full_dirs_check(:, k) = &
                    ref_pack_asymm(ref_project_asymm(diff1, dm_oao), n_param)
            else
                full_dirs_check(:, k) = &
                    ref_pack_asymm(ref_project_asymm(diff2, dm_oao), n_param)
            end if
        end do

        do j = 1, 2
            do k = 1, 2
                chol_ref_check(j, k) = dot_product(flats(:, pivot_order(j)), &
                                                   flats(:, pivot_order(k)))
            end do
        end do
        call dpotrf("U", 2_ip, chol_ref_check, 2_ip, info)

        ! clear lower triangle
        chol_ref_check(2, 1) = 0.0_rp

    end subroutine ref_pivoted_cholesky

    function ref_build_a_part(dm_diff, v_diff, linear) result(a)
        !
        ! this function independently reproduces the raw product A = S^T Y over the 
        ! flattened AO and particle dimensions, together with its symmetrization: 
        ! unweighted for the linear part, and for the non-linear part blended toward 
        ! whichever entry of each off-diagonal pair draws its response from the shorter 
        ! step
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :), v_diff(:, :, :, :)
        logical, intent(in) :: linear
        real(rp), allocatable :: a(:, :)

        integer(ip) :: n_diff, flat_len, i, k
        real(rp) :: weight
        real(rp), allocatable :: raw(:, :), step_norms(:)

        n_diff = size(dm_diff, 4)
        flat_len = size(dm_diff, 1) * size(dm_diff, 2) * size(dm_diff, 3)
        raw = matmul(transpose(reshape(dm_diff, [flat_len, n_diff])), &
                     reshape(v_diff, [flat_len, n_diff]))

        if (linear) then
            a = 0.5_rp * (raw + transpose(raw))
        else
            allocate(step_norms(n_diff))
            do i = 1, n_diff
                step_norms(i) = norm2(dm_diff(:, :, :, i))
            end do
            a = raw
            do k = 2, n_diff
                do i = 1, k - 1
                    weight = step_norms(k) / (step_norms(i) + step_norms(k))
                    a(i, k) = weight * raw(k, i) + (1.0_rp - weight) * raw(i, k)
                    a(k, i) = a(i, k)
                end do
            end do
            deallocate(step_norms)
        end if
        deallocate(raw)

    end function ref_build_a_part

    function ref_ms_a_inv(a_tilde, y_gram) result(a_inv)
        !
        ! this function independently reproduces the screened pseudoinverse the 
        ! multisecant routines build from a congruence-transformed A: eigenvalues at 
        ! the level of numerical noise are discarded, as are directions whose 
        ! eigenvalue is small against the norm of the response they divide
        !
        use otr_arh, only: eig_val_noise_factor, ms_sr1_skip_thresh

        real(rp), intent(in) :: a_tilde(:, :), y_gram(:, :)
        real(rp), allocatable :: a_inv(:, :)

        integer(ip) :: n, i, info, lwork
        real(rp) :: thresh, y_norm
        real(rp), allocatable :: vecs(:, :), vals(:), diag(:, :), work(:)
        external :: dsyev

        n = size(a_tilde, 1)
        allocate(vecs(n, n), vals(n), diag(n, n))
        vecs = a_tilde
        lwork = 3_ip * n
        allocate(work(lwork))
        call dsyev("V", "U", n, vecs, n, vals, work, lwork, info)
        deallocate(work)

        ! keep the exact inverse only where the eigenvalue is above the noise floor and
        ! the direction is not screened out by the skipping criterion
        thresh = eig_val_noise_factor * maxval(abs(vals)) * epsilon(1.0_rp)
        diag = 0.0_rp
        do i = 1, n
            y_norm = sqrt(max(dot_product(vecs(:, i), matmul(y_gram, vecs(:, i))), &
                              0.0_rp))
            if (abs(vals(i)) > thresh .and. &
                abs(vals(i)) >= ms_sr1_skip_thresh * y_norm) &
                diag(i, i) = 1.0_rp / vals(i)
        end do

        a_inv = matmul(vecs, matmul(diag, transpose(vecs)))
        deallocate(vecs, vals, diag)

    end function ref_ms_a_inv

    subroutine ref_stack_os_linear(dm_diff, v_same_spin_diff, v_opposite_spin_diff, &
                                   n_ao, s_full, y_full)
        !
        ! this subroutine independently reproduces the stacked history and response of 
        ! the open-shell linear multisecant system: each history column touches only 
        ! its own spin block, while each response column pairs the same-spin potential 
        ! of one channel with the opposite-spin potential of the other, so that S^T Y
        ! reproduces the four blocks of A and Y^T Y the response Gram matrix
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :), v_same_spin_diff(:, :, :, :), &
                                v_opposite_spin_diff(:, :, :, :)
        integer(ip), intent(in) :: n_ao
        real(rp), intent(out) :: s_full(:, :), y_full(:, :)

        integer(ip) :: n_diff, n_ao2, k

        n_diff = size(dm_diff, 4)
        n_ao2 = n_ao * n_ao
        s_full = 0.0_rp
        do k = 1, n_diff
            s_full(:n_ao2, k) = reshape(dm_diff(:, :, 1, k), [n_ao2])
            s_full(n_ao2 + 1:, n_diff + k) = reshape(dm_diff(:, :, 2, k), [n_ao2])
            y_full(:n_ao2, k) = reshape(v_same_spin_diff(:, :, 1, k), [n_ao2])
            y_full(n_ao2 + 1:, k) = reshape(v_opposite_spin_diff(:, :, 2, k), [n_ao2])
            y_full(:n_ao2, n_diff + k) = &
                reshape(v_opposite_spin_diff(:, :, 1, k), [n_ao2])
            y_full(n_ao2 + 1:, n_diff + k) = &
                reshape(v_same_spin_diff(:, :, 2, k), [n_ao2])
        end do

    end subroutine ref_stack_os_linear

    function ref_congruence_transform(a, map, chol) result(a_tilde)
        !
        ! this function independently reproduces the congruence transformation used to 
        ! rebase a matrix into the orthonormalized S-basis, A -> R^-T (P A P^T) R^-1, 
        ! with the two triangular solves written out as explicit substitutions so that 
        ! no production routine is involved
        !
        real(rp), intent(in) :: a(:, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: a_tilde(:, :)

        integer(ip) :: n, i, j, k

        ! select and reorder rows/columns according to map
        n = size(map)
        allocate(a_tilde(n, n))
        do j = 1, n
            do i = 1, n
                a_tilde(i, j) = a(map(i), map(j))
            end do
        end do

        ! right solve X R = P A P^T by forward substitution along each row
        do i = 1, n
            do j = 1, n
                do k = 1, j - 1
                    a_tilde(i, j) = a_tilde(i, j) - a_tilde(i, k) * chol(k, j)
                end do
                a_tilde(i, j) = a_tilde(i, j) / chol(j, j)
            end do
        end do

        ! left solve R^T Y = X by forward substitution along each column
        do j = 1, n
            do i = 1, n
                do k = 1, i - 1
                    a_tilde(i, j) = a_tilde(i, j) - chol(k, i) * a_tilde(k, j)
                end do
                a_tilde(i, j) = a_tilde(i, j) / chol(i, i)
            end do
        end do

    end function ref_congruence_transform

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
                                      generate_random_density_matrix

        integer(ip), parameter :: n_particle = 1, n_electrons = 2, &
                                  n_param = n_ao * (n_ao - 1) / 2

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: dm_saved(n_ao, n_ao, n_particle), &
                    fock_saved(n_ao, n_ao, n_particle), &
                    v_nonlinear_saved(n_ao, n_ao, n_particle), &
                    dm_saved_2(n_ao, n_ao, n_particle), &
                    fock_saved_2(n_ao, n_ao, n_particle), &
                    v_nonlinear_saved_2(n_ao, n_ao, n_particle), &
                    kappa(n_param), grad(n_param), h_diag(n_param), func, &
                    dm_diff_check(n_ao, n_ao, n_particle, 2), &
                    full_dirs_check(n_param, 2), chol_ref_check(2, 2)
        integer(ip) :: n_diff, info, error, pivot_order(2)
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
        n_diff = size(arh_object%dm_dirs, 2)
        if (n_diff /= 2) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Two independent "// &
                "history entries were not both accepted."
            test_update_orbs_arh_cs = .false.
            return
        end if
        dm_diff_check(:, :, :, 1) = dm_saved_2 - arh_object%dm_oao
        dm_diff_check(:, :, :, 2) = dm_saved - arh_object%dm_oao
        call ref_pivoted_cholesky(dm_diff_check(:, :, :, 1), &
                                  dm_diff_check(:, :, :, 2), arh_object%dm_oao, &
                                  n_param, pivot_order, full_dirs_check, &
                                  chol_ref_check, info)
        if (info /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Reference Cholesky "// &
                "factorization of the raw history Gram matrix failed."
            test_update_orbs_arh_cs = .false.
        else
            if (norm2(matmul(arh_object%dm_dirs, chol_ref_check) - full_dirs_check) > &
                tol) then
                write (stderr, *) "test_update_orbs_arh_cs failed: Density matrix "// &
                    "difference directions are not rebased by the inverse of the "// &
                    "Cholesky factor of the raw history Gram matrix."
                test_update_orbs_arh_cs = .false.
            end if
        end if
        if (size(arh_object%potential_dirs, 2) /= n_diff) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Potential "// &
                "difference directions are not rebased onto the same basis as the "// &
                "density matrix difference directions."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%a_sym - transpose(arh_object%a_sym)) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Symmetrized A "// &
                "matrix is not symmetric."
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
        if (size(arh_object%linear_potential_dirs, 2) /= n_diff .or. &
            size(arh_object%nonlinear_potential_dirs, 2) /= n_diff) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Two independent "// &
                "history entries were not both accepted for multisecant SR1."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (size(arh_object%a_inv, 1) /= n_diff .or. &
            norm2(arh_object%a_inv - transpose(arh_object%a_inv)) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Multisecant SR1 "// &
                "pseudoinverse is not a symmetric matrix matching the linear "// &
                "potential difference directions."
            test_update_orbs_arh_cs = .false.
        end if
        if (size(arh_object%a_inv_comb, 1) /= n_diff .or. &
            norm2(arh_object%a_inv_comb - transpose(arh_object%a_inv_comb)) > tol) &
        then
            write (stderr, *) "test_update_orbs_arh_cs failed: Spin-combined "// &
                "multisecant SR1 pseudoinverse is not a symmetric matrix matching "// &
                "the non-linear potential difference directions."
            test_update_orbs_arh_cs = .false.
        end if

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
                                      generate_random_density_matrix

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
                    grad(n_param), h_diag(n_param), func, &
                    dm_diff_check(n_ao, n_ao, n_particle, 2), &
                    embedded_check(n_ao, n_ao, n_particle, 2), &
                    full_dirs_check(n_param, 2 * n_particle), &
                    full_dirs_block(n_param, 2), chol_ref_check(2, 2)
        integer(ip) :: i, j, col, n_diff, info, error, pivot_order(2)
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
        col = size(arh_object%dm_dirs, 2)
        if (col /= n_particle * n_diff) then
            write (stderr, *) "test_update_orbs_arh_os failed: Two independent "// &
                "history entries were not both accepted in both channels."
            test_update_orbs_arh_os = .false.
            return
        end if
        dm_diff_check(:, :, :, 1) = dm_saved_2 - arh_object%dm_oao
        dm_diff_check(:, :, :, 2) = dm_saved - arh_object%dm_oao
        do j = 1, n_particle
            embedded_check = 0.0_rp
            embedded_check(:, :, j, :) = dm_diff_check(:, :, j, :)
            call ref_pivoted_cholesky(embedded_check(:, :, :, 1), &
                                      embedded_check(:, :, :, 2), arh_object%dm_oao, &
                                      n_param, pivot_order, full_dirs_block, &
                                      chol_ref_check, info)
            full_dirs_check(:, (j - 1) * n_diff + 1:j * n_diff) = full_dirs_block
            if (info /= 0) then
                write (stderr, *) "test_update_orbs_arh_os failed: Reference "// &
                    "Cholesky factorization of the raw per-channel history Gram "// &
                    "matrix failed."
                test_update_orbs_arh_os = .false.
            else
                if (norm2(matmul( &
                    arh_object%dm_dirs(:, (j - 1) * n_diff + 1:j * n_diff), &
                    chol_ref_check) - full_dirs_block) > tol) then
                    write (stderr, *) "test_update_orbs_arh_os failed: Density "// &
                        "matrix difference directions are not rebased by the "// &
                        "inverse of the per-channel Cholesky factor."
                    test_update_orbs_arh_os = .false.
                end if
            end if
        end do
        if (size(arh_object%potential_dirs, 2) /= col) then
            write (stderr, *) "test_update_orbs_arh_os failed: Potential "// &
                "difference directions are not rebased onto the same basis as the "// &
                "density matrix difference directions."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%a_sym - transpose(arh_object%a_sym)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Symmetrized A "// &
                "matrix is not symmetric."
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

        ! determine that both channels' history entries were accepted as independent 
        ! by their respective factorizations, and that the inverse A matrices are 
        ! symmetric matrices matching the direction sets they are meant to couple
        n_diff = size(arh_object%dm_list, 4)
        if (size(arh_object%nonlinear_potential_dirs, 2) /= n_diff) then
            write (stderr, *) "test_update_orbs_arh_os failed: Two independent "// &
                "history entries were not both accepted for the spin-combined "// &
                "non-linear multisecant SR1 system."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (size(arh_object%linear_potential_dirs, 2) /= n_particle * n_diff) then
            write (stderr, *) "test_update_orbs_arh_os failed: Two independent "// &
                "history entries were not both accepted in both channels for the "// &
                "spin-separated linear multisecant SR1 system."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (size(arh_object%a_inv, 1) /= n_particle * n_diff .or. &
            norm2(arh_object%a_inv - transpose(arh_object%a_inv)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Spin-separated "// &
                "multisecant SR1 pseudoinverse is not a symmetric matrix matching "// &
                "the linear potential difference directions."
            test_update_orbs_arh_os = .false.
        end if
        if (size(arh_object%a_inv_comb, 1) /= n_diff .or. &
            norm2(arh_object%a_inv_comb - transpose(arh_object%a_inv_comb)) > tol) &
        then
            write (stderr, *) "test_update_orbs_arh_os failed: Spin-combined "// &
                "multisecant SR1 pseudoinverse is not a symmetric matrix matching "// &
                "the non-linear potential difference directions."
            test_update_orbs_arh_os = .false.
        end if

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

    logical(c_bool) function test_cache_history_dirs() bind(C)
        !
        ! this function tests the routine which caches the packed history-projection
        ! directions the low-rank part of the approximate Hessian is built from, by
        ! projecting and packing every history entry and rebasing the result via a
        ! given (map, chol) pair; undoing the rebasing by right-multiplying by chol
        ! again must reproduce the gathered, reordered raw projections exactly
        !
        use otr_arh, only: cache_history_dirs
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao_unit_tests, only: generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_list = 2, n_electrons = 1

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    v_diff(n_ao, n_ao, n_particle, n_list), &
                    empty_v_diff(n_ao, n_ao, n_particle, 0), chol(n_list, n_list), &
                    empty_chol(0, 0)
        real(rp), allocatable :: dirs(:, :), raw(:, :)
        integer(ip) :: j, k, map(n_list), empty_map(0)

        ! assume tests pass
        test_cache_history_dirs = .true.

        ! generate a random density matrix and random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)
            do k = 1, n_list
                v_diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! independently reproduce the raw (un-rebased) projections
        raw = ref_cache_dirs(v_diff, dm_oao, n_param)

        ! reorder and rebase into an arbitrary orthonormalized basis
        map = [2, 1]
        chol = generate_random_upper_triangular(n_list)

        ! call routine and determine if dimensions are correct and undoing the 
        ! reproduces the gathered, reordered raw projections
        call cache_history_dirs(v_diff, dm_oao, n_list, n_param, map, chol, dirs)
        if (size(dirs, 1) /= n_param .or. size(dirs, 2) /= n_list) then
            write (stderr, *) "test_cache_history_dirs failed: Incorrect "// &
                "dimensions of directions."
            test_cache_history_dirs = .false.
        else if (norm2(matmul(dirs, chol) - raw(:, map)) > tol) then
            write (stderr, *) "test_cache_history_dirs failed: Directions are not "// &
                "rebased directions of the gathered, reordered raw projections."
            test_cache_history_dirs = .false.
        end if
        deallocate(dirs, raw)

        ! call routine for an empty history and determine if no directions are
        ! returned
        call cache_history_dirs(empty_v_diff, dm_oao, 0_ip, n_param, empty_map, &
                                empty_chol, dirs)
        if (size(dirs, 1) /= n_param .or. size(dirs, 2) /= 0) then
            write (stderr, *) "test_cache_history_dirs failed: Incorrect "// &
                "dimensions of directions for empty history."
            test_cache_history_dirs = .false.
        end if
        deallocate(dirs)

    end function test_cache_history_dirs

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

    logical(c_bool) function test_cache_channel_split_dirs() bind(C)
        !
        ! this function tests the open-shell history-projection caching routine which
        ! embeds a history entry into one spin channel before projecting and packing
        ! it, keeping the two channels as separate columns, then rebases the result
        ! via a given (map, chol) pair; undoing the rebasing by right-multiplying by
        ! chol again must reproduce the gathered, reordered raw projections exactly
        !
        use otr_arh, only: cache_channel_split_dirs
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao_unit_tests, only: ref_project_asymm, ref_pack_asymm, &
                                      generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_list = 2, n_electrons = 1, &
                                  n_col = n_particle * n_list

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    v_same(n_ao, n_ao, n_particle, n_list), chol(n_col, n_col)
        real(rp), allocatable :: u(:, :), raw(:, :)
        integer(ip) :: j, k, map(n_col)

        ! assume test passes
        test_cache_channel_split_dirs = .true.

        ! generate a random density matrix and random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)
            do k = 1, n_list
                v_same(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! independently reproduce the raw (un-rebased) split projections
        allocate(raw(n_param, n_col))
        do k = 1, n_list
            do j = 1, n_particle
                raw(:, (j - 1)*n_list + k) = ref_pack_asymm(ref_project_asymm( &
                    embed_channel(v_same(:, :, j, k), j, n_ao, n_particle), dm_oao), &
                    n_param)
            end do
        end do

        ! reorder and rebase into an arbitrary orthonormalized basis
        map = [3, 1, 4, 2]
        chol = generate_random_upper_triangular(n_col)

        ! call routine and determine if dimensions are correct and undoing the
        ! rebasing reproduces the gathered, reordered raw split projections
        call cache_channel_split_dirs(v_same, dm_oao, n_list, n_param, n_particle, &
                                      map, chol, u)
        if (size(u, 1) /= n_param .or. size(u, 2) /= n_col) then
            write (stderr, *) "test_cache_channel_split_dirs failed: "// &
                "Incorrect dimensions of split directions."
            test_cache_channel_split_dirs = .false.
        else if (norm2(matmul(u, chol) - raw(:, map)) > tol) then
            write (stderr, *) "test_cache_channel_split_dirs failed: Split "// &
                "directions are not rebased directions of the gathered, reordered "// &
                "raw split projections."
            test_cache_channel_split_dirs = .false.
        end if
        deallocate(u, raw)

    end function test_cache_channel_split_dirs

    logical(c_bool) function test_cache_combined_channel_dirs() bind(C)
        !
        ! this function tests the open-shell history-projection caching routine which
        ! embeds a history entry into one spin channel before projecting and packing
        ! it, summing a same-spin channel with the opposite-spin channel of the other
        ! spin, then rebases the result via a given (map, chol) pair; undoing the
        ! rebasing by right-multiplying by chol again must reproduce the gathered,
        ! reordered raw projections exactly
        !
        use otr_arh, only: cache_combined_channel_dirs
        use otr_oao_test_reference, only: n_ao, n_particle, n_param
        use otr_oao_unit_tests, only: ref_project_asymm, ref_pack_asymm, &
                                      generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_list = 2, n_electrons = 1, &
                                  n_col = n_particle * n_list

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    v_same(n_ao, n_ao, n_particle, n_list), &
                    v_opp(n_ao, n_ao, n_particle, n_list), chol(n_col, n_col)
        real(rp), allocatable :: u(:, :), raw(:, :)
        integer(ip) :: j, k, map(n_col)

        ! assume test passes
        test_cache_combined_channel_dirs = .true.

        ! generate a random density matrix and random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)
            do k = 1, n_list
                v_same(:, :, j, k) = generate_random_symm_matrix(n_ao)
                v_opp(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do

        ! independently reproduce the raw (un-rebased) combined projections
        allocate(raw(n_param, n_col))
        do k = 1, n_list
            raw(:, k) = ref_pack_asymm(ref_project_asymm( &
                embed_channel(v_same(:, :, 1, k), 1_ip, n_ao, n_particle) + &
                embed_channel(v_opp(:, :, 2, k), 2_ip, n_ao, n_particle), dm_oao), &
                n_param)
            raw(:, n_list + k) = ref_pack_asymm(ref_project_asymm( &
                embed_channel(v_same(:, :, 2, k), 2_ip, n_ao, n_particle) + &
                embed_channel(v_opp(:, :, 1, k), 1_ip, n_ao, n_particle), dm_oao), &
                n_param)
        end do

        ! reorder and rebase into an arbitrary orthonormalized basis
        map = [3, 1, 4, 2]
        chol = generate_random_upper_triangular(n_col)

        ! call routine and determine if dimensions are correct and undoing the
        ! rebasing reproduces the gathered, reordered raw combined projections
        call cache_combined_channel_dirs(v_same, v_opp, dm_oao, n_list, n_param, &
                                         n_particle, map, chol, u)
        if (size(u, 1) /= n_param .or. size(u, 2) /= n_col) then
            write (stderr, *) "test_cache_combined_channel_dirs failed: "// &
                "Incorrect dimensions of combined directions."
            test_cache_combined_channel_dirs = .false.
        else if (norm2(matmul(u, chol) - raw(:, map)) > tol) then
            write (stderr, *) "test_cache_combined_channel_dirs failed: Combined "// &
                "directions are not rebased directions of the gathered, reordered "// &
                "raw combined projections."
            test_cache_combined_channel_dirs = .false.
        end if
        deallocate(u, raw)

    end function test_cache_combined_channel_dirs

    logical(c_bool) function test_get_low_rank_hess_factors() bind(C)
        !
        ! this function tests the subroutine which assembles the low-rank part of the
        ! approximate Hessian
        !
        use otr_arh, only: get_low_rank_hess_factors, arh_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_unit_tests, only: generate_random_symm_matrix, identity_matrix

        integer(ip), parameter :: n_param = 4, n_diff = 2

        integer(ip), target :: n_param_target, n_particle_target
        real(rp) :: density(n_param, n_diff), potential(n_param, n_diff), &
                    linear(n_param, n_diff), nonlinear(n_param, n_diff), &
                    a_sym(n_diff, n_diff), a_inv(n_diff, n_diff), &
                    a_inv_comb(n_diff, n_diff)
        logical :: passed_eq

        ! assume tests pass
        test_get_low_rank_hess_factors = .true.

        ! generate the packed history directions and the small dense matrices the
        ! coupling matrices are assembled from
        call random_number(density)
        call random_number(potential)
        call random_number(linear)
        call random_number(nonlinear)
        a_sym = generate_random_symm_matrix(n_diff)
        a_inv = generate_random_symm_matrix(n_diff)
        a_inv_comb = generate_random_symm_matrix(n_diff)

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
        arh_object%a_sym = a_sym
        arh_object%a_inv = a_inv
        arh_object%a_inv_comb = a_inv_comb

        ! multisecant SR1 stacks the linear and non-linear directions and couples each 
        ! block through its own separately regularized system
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
        ! difference history alone, coupled directly through the (already
        ! congruence-transformed) symmetrized A matrix
        arh_object%settings%arh_type = "ms_sp"
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%expansion_dirs - density) > tol) .or. &
            any(abs(arh_object%projection_dirs - density) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "directions for subspace-projected multisecant."
            test_get_low_rank_hess_factors = .false.
        end if
        if (any(abs(arh_object%coupling_matrix - 8.0_rp * a_sym) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "coupling matrix for subspace-projected multisecant."
            test_get_low_rank_hess_factors = .false.
        end if

        ! symmetrized ARH couples the density and potential difference histories in
        ! both directions with a plain identity block, leaving the diagonal blocks empty
        arh_object%settings%arh_type = "symm_arh"
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%expansion_dirs(:, :n_diff) - density) > tol) .or. &
            any(abs(arh_object%expansion_dirs(:, n_diff + 1:) - potential) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "expansion directions for symmetrized ARH."
            test_get_low_rank_hess_factors = .false.
        end if
        if (any(abs(arh_object%coupling_matrix(:n_diff, n_diff + 1:) - 4.0_rp * &
                    identity_matrix(n_diff)) > tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, :n_diff) - 4.0_rp * &
                    identity_matrix(n_diff)) > tol) .or. &
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
        if (any(abs(arh_object%coupling_matrix(:n_diff, :n_diff) + 8.0_rp * a_sym) > &
                tol) .or. &
            any(abs(arh_object%coupling_matrix(:n_diff, n_diff + 1:) - 8.0_rp * &
                    identity_matrix(n_diff)) > tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, :n_diff) - 8.0_rp * &
                    identity_matrix(n_diff)) > tol) .or. &
            any(abs(arh_object%coupling_matrix(n_diff + 1:, n_diff + 1:)) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "coupling matrix for multisecant PSB."
            test_get_low_rank_hess_factors = .false.
        end if

        ! standard ARH is the only type whose expansion and projection directions
        ! differ, expanding in the potential and contracting against the density
        ! difference history, coupled by a plain identity
        arh_object%settings%arh_type = "arh"
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%expansion_dirs - potential) > tol) .or. &
            any(abs(arh_object%projection_dirs - density) > tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "directions for standard ARH."
            test_get_low_rank_hess_factors = .false.
        end if
        if (any(abs(arh_object%coupling_matrix - 8.0_rp * identity_matrix(n_diff)) > &
                tol)) then
            write (stderr, *) "test_get_low_rank_hess_factors failed: Incorrect "// &
                "coupling matrix for standard ARH."
            test_get_low_rank_hess_factors = .false.
        end if

        ! test the open-shell coupling matrix which is exactly half its closed-shell
        ! counterpart,
        n_particle_target = 2
        call get_low_rank_hess_factors()
        if (any(abs(arh_object%coupling_matrix - 4.0_rp * identity_matrix(n_diff)) > &
                tol)) then
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

    logical(c_bool) function test_build_a_part() bind(C)
        !
        ! this function tests the subroutine which constructs A = S^T Y for a single
        ! part of a coupled potential-difference response and symmetrizes it, for both
        ! the linear and the non-linear case
        !
        use otr_arh, only: build_a_part
        use otr_oao_test_reference, only: n_ao, n_particle

        integer(ip), parameter :: n_diff = 3
        character(10), parameter :: case_names(2) = &
            [character(10) :: "linear", "non-linear"]
        integer(ip) :: i_case
        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_diff(n_ao, n_ao, n_particle, n_diff), a(n_diff, n_diff), &
                    a_linear(n_diff, n_diff), expected(n_diff, n_diff)
        logical :: linear

        ! assume tests pass
        test_build_a_part = .true.

        ! random history and response so that the raw product is generically asymmetric
        ! and the two symmetrizations are distinguishable
        call random_number(dm_diff)
        call random_number(v_diff)

        do i_case = 1, 2
            linear = i_case == 1

            ! get expected matrix
            expected = ref_build_a_part(dm_diff, v_diff, linear)

            ! call routine and determine if the symmetrized matrix matches
            call build_a_part(dm_diff, v_diff, linear, a)
            if (maxval(abs(a - expected)) > tol) then
                write (stderr, *) "test_build_a_part failed for the "// &
                    trim(case_names(i_case))//" part: Incorrect symmetrized A matrix."
                test_build_a_part = .false.
            end if
            if (linear) a_linear = a
        end do

        ! the weighting has to change the result, or the non-linear case has silently
        ! reduced to the linear one and asserts nothing
        if (maxval(abs(a - a_linear)) < tol) then
            write (stderr, *) "test_build_a_part failed: Weighted symmetrization "// &
                "reproduced the exact symmetrization."
            test_build_a_part = .false.
        end if

    end function test_build_a_part

    logical(c_bool) function test_build_a_sym_cs() bind(C)
        !
        ! this function tests the function which builds the A = S^T Y matrix,
        ! symmetrizing its linear and non-linear contributions separately and then
        ! congruence-transforming the result: undoing that transform by left- and
        ! right-multiplying by chol again must reproduce the gathered, reordered
        ! raw symmetrized matrix exactly
        !
        use otr_arh, only: build_a_sym_cs
        use otr_oao_test_reference, only: n_ao

        integer(ip), parameter :: n_particle = 1, n_diff = 3

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_linear_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_nonlinear_diff(n_ao, n_ao, n_particle, n_diff), &
                    expected_a_sym(n_diff, n_diff), chol(n_diff, n_diff)
        real(rp), allocatable :: a_sym(:, :)
        integer(ip) :: map(n_diff)

        ! assume tests pass
        test_build_a_sym_cs = .true.

        ! random density matrix differences and linear and non-linear potential
        ! differences, so that both symmetrizations act on generically asymmetric
        ! contributions
        call random_number(dm_diff)
        call random_number(v_linear_diff)
        call random_number(v_nonlinear_diff)

        ! the two contributions are symmetrized separately and summed, reproduced here
        ! independently of the routine
        expected_a_sym = ref_build_a_part(dm_diff, v_linear_diff, .true.) + &
                         ref_build_a_part(dm_diff, v_nonlinear_diff, .false.)

        ! reorder and congruence-transform into an arbitrary orthonormalized basis
        map = [2, 3, 1]
        chol = generate_random_upper_triangular(n_diff)

        ! generate A matrix and verify
        a_sym = build_a_sym_cs(dm_diff, v_linear_diff, v_nonlinear_diff, map, chol)
        if (size(a_sym, 1) /= n_diff .or. size(a_sym, 2) /= n_diff) then
            write (stderr, *) "test_build_a_sym_cs failed: Incorrect dimensions."
            test_build_a_sym_cs = .false.
            return
        end if
        if (norm2(matmul(transpose(chol), matmul(a_sym, chol)) - &
                  expected_a_sym(map, map)) > tol) then
            write (stderr, *) "test_build_a_sym_cs failed: Result does not "// &
                "invert back to the gathered, reordered raw symmetrized matrix."
            test_build_a_sym_cs = .false.
        end if

    end function test_build_a_sym_cs

    logical(c_bool) function test_build_a_block_sym_os() bind(C)
        !
        ! this function tests the function which builds the cross-channel-symmetrized 
        ! open-shell A = S^T Y matrix and congruence-transforms it: undoing that
        ! transform by left- and right-multiplying by chol again must reproduce the
        ! gathered, reordered raw block matrix exactly
        !
        use otr_arh, only: build_a_block_sym_os
        use otr_oao_test_reference, only: n_ao, n_particle

        integer(ip), parameter :: n_diff = 3, n_col = n_particle * n_diff

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_same_linear(n_ao, n_ao, n_particle, n_diff), &
                    v_nonlinear(n_ao, n_ao, n_particle, n_diff), &
                    v_opp(n_ao, n_ao, n_particle, n_diff), &
                    expected_a_block(n_col, n_col), chol(n_col, n_col), &
                    a_opp(n_diff, n_diff, n_particle), averaged(n_diff, n_diff), &
                    dm_j(n_ao, n_ao, 1, n_diff), v_j(n_ao, n_ao, 1, n_diff)
        real(rp), allocatable :: a_block(:, :)
        integer(ip) :: map(n_col), j, lo, hi

        ! assume tests pass
        test_build_a_block_sym_os = .true.

        ! random per-channel histories and same-spin, non-linear and opposite-spin
        ! potentials, so that every block is generically asymmetric
        call random_number(dm_diff)
        call random_number(v_same_linear)
        call random_number(v_nonlinear)
        call random_number(v_opp)

        ! each diagonal block sums the separately symmetrized linear and non-linear
        ! same-spin contributions of one channel, while the raw opposite-spin products
        ! form the off-diagonal blocks
        do j = 1, n_particle
            lo = (j - 1) * n_diff + 1
            hi = j * n_diff
            dm_j = reshape(dm_diff(:, :, j, :), [n_ao, n_ao, 1_ip, n_diff])
            v_j = reshape(v_same_linear(:, :, j, :), [n_ao, n_ao, 1_ip, n_diff])
            expected_a_block(lo:hi, lo:hi) = ref_build_a_part(dm_j, v_j, .true.)
            v_j = reshape(v_nonlinear(:, :, j, :), [n_ao, n_ao, 1_ip, n_diff])
            expected_a_block(lo:hi, lo:hi) = expected_a_block(lo:hi, lo:hi) + &
                                             ref_build_a_part(dm_j, v_j, .false.)
            a_opp(:, :, j) = matmul( &
                transpose(reshape(dm_diff(:, :, j, :), [n_ao * n_ao, n_diff])), &
                reshape(v_opp(:, :, j, :), [n_ao * n_ao, n_diff]))
        end do

        ! the two opposite-spin blocks are exact transposes of one another, so any
        ! mismatch between them is noise and is averaged away rather than blended
        averaged = 0.5_rp * (a_opp(:, :, 1) + transpose(a_opp(:, :, 2)))
        expected_a_block(:n_diff, n_diff + 1:) = averaged
        expected_a_block(n_diff + 1:, :n_diff) = transpose(averaged)

        ! reorder and congruence-transform into an arbitrary orthonormalized basis
        map = [4, 1, 6, 2, 5, 3]
        chol = generate_random_upper_triangular(n_col)

        ! generate A matrix and verify
        a_block = build_a_block_sym_os(dm_diff, v_same_linear, v_nonlinear, v_opp, &
                                       n_ao, map, chol)
        if (size(a_block, 1) /= n_col .or. size(a_block, 2) /= n_col) then
            write (stderr, *) "test_build_a_block_sym_os failed: Incorrect dimensions."
            test_build_a_block_sym_os = .false.
            return
        end if
        if (norm2(matmul(transpose(chol), matmul(a_block, chol)) - &
                  expected_a_block(map, map)) > tol) then
            write (stderr, *) "test_build_a_block_sym_os failed: Result does not "// &
                "invert back to the gathered, reordered raw block matrix."
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

    logical(c_bool) function test_cross_symmetrize_exact() bind(C)
        !
        ! this function tests the subroutine which performs a plain, unweighted
        ! cross-symmetrization between two related off-diagonal blocks of a larger 
        ! matrix
        !
        use otr_arh, only: cross_symmetrize_exact

        real(rp) :: a12(2, 2), a21(2, 2), expected12(2, 2), expected21(2, 2)

        ! assume tests pass
        test_cross_symmetrize_exact = .true.

        ! initialize blocks which are not transposes of each other
        a12 = reshape([1.0_rp, 3.0_rp, &
                       2.0_rp, 4.0_rp], [2, 2])
        a21 = reshape([5.0_rp, 6.0_rp, &
                       7.0_rp, 9.0_rp], [2, 2])

        ! initialize expected blocks, averaging a12(i, k) with a21(k, i) directly
        expected12 = reshape([3.0_rp, 5.0_rp, &
                              4.0_rp, 6.5_rp], [2, 2])
        expected21 = transpose(expected12)

        ! call routine and determine if values of resulting blocks match and if
        ! these are transposes of each other
        call cross_symmetrize_exact(a12, a21)
        if (norm2(a12 - expected12) > tol) then
            write (stderr, *) "test_cross_symmetrize_exact failed: Incorrect first "// &
                "block values after cross-symmetrization."
            test_cross_symmetrize_exact = .false.
        end if
        if (norm2(a21 - expected21) > tol) then
            write (stderr, *) "test_cross_symmetrize_exact failed: Incorrect "// &
                "second block values after cross-symmetrization."
            test_cross_symmetrize_exact = .false.
        end if
        if (norm2(a21 - transpose(a12)) > tol) then
            write (stderr, *) "test_cross_symmetrize_exact failed: Blocks are not "// &
                "transposes of each other after cross-symmetrization."
            test_cross_symmetrize_exact = .false.
        end if

    end function test_cross_symmetrize_exact

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

    logical(c_bool) function test_factorize_history() bind(C)
        !
        ! this function tests the subroutine which performs a pivoted, rank-revealing 
        ! Cholesky factorization of the Gram matrix of a set of flattened history 
        ! vectors
        !
        use otr_arh, only: factorize_history

        integer(ip), parameter :: n_ao = 2, n_diff = 3

        real(rp) :: dm_diff(n_ao, n_ao, n_diff), dm_diff_empty(n_ao, n_ao, 0), &
                    expected_chol(2, 2)
        real(rp), allocatable :: chol(:, :)
        integer(ip), allocatable :: map(:)
        integer(ip) :: n_accepted

        ! assume tests pass
        test_factorize_history = .true.

        ! initialize density matrix differences: the first column is a smaller multiple 
        ! of the same direction as the second, larger column, and the  third column is 
        ! orthogonal to both
        dm_diff(:, :, 1) = reshape([0.5_rp, 0.0_rp, &
                                    0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2) = reshape([2.0_rp, 0.0_rp, &
                                    0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 3) = reshape([0.0_rp, 1.0_rp, &
                                    1.0_rp, 0.0_rp], [n_ao, n_ao])

        ! the Gram matrix is [[0.25, 1, 0], [1, 4, 0], [0, 0, 2]]: columns 1 and 2
        ! are parallel, so pivoting accepts the larger-magnitude column 2 first,
        ! then the orthogonal column 3 (residual 2, exceeding column 1's now-zero
        ! residual), permanently rejecting the smaller column 1
        expected_chol = reshape([2.0_rp, 0.0_rp, &
                                 0.0_rp, sqrt(2.0_rp)], [2, 2])

        ! call routine and determine if the accepted count, map and Cholesky factor
        ! match
        call factorize_history(reshape(dm_diff, [n_ao * n_ao, n_diff]), chol, map, &
                               n_accepted)
        if (n_accepted /= 2) then
            write (stderr, *) "test_factorize_history failed: Incorrect number of "// &
                "accepted columns."
            test_factorize_history = .false.
            return
        end if
        if (any(map /= [2, 3])) then
            write (stderr, *) "test_factorize_history failed: Incorrect map back "// &
                "to original history indices."
            test_factorize_history = .false.
            return
        end if
        if (norm2(chol - expected_chol) > tol) then
            write (stderr, *) "test_factorize_history failed: Incorrect Cholesky "// &
                "factor."
            test_factorize_history = .false.
        end if
        if (norm2(matmul(transpose(chol), chol) - &
                  reshape([4.0_rp, 0.0_rp, 0.0_rp, 2.0_rp], [2, 2])) > tol) then
            write (stderr, *) "test_factorize_history failed: Cholesky factor "// &
                "does not reproduce the Gram matrix of the accepted columns."
            test_factorize_history = .false.
        end if
        deallocate(chol, map)

        ! call routine for an empty history and determine if an empty factorization is 
        ! returned
        call factorize_history(reshape(dm_diff_empty, [n_ao * n_ao, 0_ip]), chol, map, &
                               n_accepted)
        if (n_accepted /= 0 .or. size(chol, 1) /= 0 .or. size(map) /= 0) then
            write (stderr, *) "test_factorize_history failed: Incorrect "// &
                "factorization of an empty history."
            test_factorize_history = .false.
        end if
        deallocate(chol, map)

    end function test_factorize_history

    logical(c_bool) function test_rebase_dirs() bind(C)
        !
        ! this function tests the function which re-expresses a set of packed
        ! history-direction columns in the orthonormalized basis defined by a Cholesky 
        ! factor: selecting, reordering and right-dividing by the factor should be 
        ! exactly undone by right-multiplying by the factor again
        !
        use otr_arh, only: rebase_dirs

        integer(ip), parameter :: n_param = 3, n_dm = 3, n_accepted = 2

        real(rp) :: dirs(n_param, n_dm), chol(n_accepted, n_accepted)
        integer(ip) :: map(n_accepted)
        real(rp), allocatable :: rebased(:, :)

        ! assume tests pass
        test_rebase_dirs = .true.

        call random_number(dirs)
        chol = generate_random_upper_triangular(n_accepted)
        map = [3, 1]

        rebased = rebase_dirs(dirs, map, chol)
        if (size(rebased, 1) /= n_param .or. size(rebased, 2) /= n_accepted) then
            write (stderr, *) "test_rebase_dirs failed: Incorrect dimensions."
            test_rebase_dirs = .false.
            return
        end if

        ! undoing the right-division by right-multiplying by the same factor must
        ! reproduce the gathered, reordered columns exactly
        if (norm2(matmul(rebased, chol) - dirs(:, map)) > tol) then
            write (stderr, *) "test_rebase_dirs failed: Rebased directions do not "// &
                "invert back to the gathered, reordered original columns."
            test_rebase_dirs = .false.
        end if

    end function test_rebase_dirs

    logical(c_bool) function test_congruence_transform() bind(C)
        !
        ! this function tests the function which applies the congruence transformation 
        ! A -> R^-T (P A P^T) R^-1: undoing it by left- and right-multiplying by the 
        ! same factor must reproduce the gathered, reordered original matrix exactly
        !
        use otr_arh, only: congruence_transform
        use otr_oao_unit_tests, only: generate_random_symm_matrix

        integer(ip), parameter :: n_dm = 3, n_accepted = 2

        real(rp) :: a(n_dm, n_dm), chol(n_accepted, n_accepted), &
                    a_gathered(n_accepted, n_accepted)
        integer(ip) :: map(n_accepted)
        real(rp), allocatable :: a_tilde(:, :)

        ! assume tests pass
        test_congruence_transform = .true.

        a = generate_random_symm_matrix(n_dm)
        chol = generate_random_upper_triangular(n_accepted)
        map = [3, 1]
        a_gathered = a(map, map)

        a_tilde = congruence_transform(a, map, chol)
        if (size(a_tilde, 1) /= n_accepted .or. size(a_tilde, 2) /= n_accepted) then
            write (stderr, *) "test_congruence_transform failed: Incorrect "// &
                "dimensions."
            test_congruence_transform = .false.
            return
        end if
        if (norm2(a_tilde - transpose(a_tilde)) > tol) then
            write (stderr, *) "test_congruence_transform failed: Result is not "// &
                "symmetric even though the input and the transform preserve symmetry."
            test_congruence_transform = .false.
        end if
        if (norm2(matmul(transpose(chol), matmul(a_tilde, chol)) - a_gathered) > &
            tol) then
            write (stderr, *) "test_congruence_transform failed: Result does not "// &
                "invert back to the gathered, reordered original matrix."
            test_congruence_transform = .false.
        end if

    end function test_congruence_transform

    logical(c_bool) function test_combine_channels() bind(C)
        !
        ! this function tests the subroutine which assembles a block-diagonal Cholesky 
        ! factor and concatenated, offset index map from two independent per-channel 
        ! history factorizations
        !
        use otr_arh, only: combine_channels

        integer(ip), parameter :: n1 = 2, n2 = 1, n_offset = 5

        real(rp) :: chol1(n1, n1), chol2(n2, n2)
        integer(ip) :: map1(n1), map2(n2)
        real(rp), allocatable :: chol_comb(:, :)
        integer(ip), allocatable :: map_comb(:)

        ! assume tests pass
        test_combine_channels = .true.

        chol1 = reshape([1.0_rp, 0.0_rp, 2.0_rp, 3.0_rp], [n1, n1])
        chol2 = reshape([4.0_rp], [n2, n2])
        map1 = [2, 4]
        map2 = [1]

        call combine_channels(chol1, map1, chol2, map2, n_offset, chol_comb, map_comb)
        if (size(chol_comb, 1) /= n1 + n2 .or. size(map_comb) /= n1 + n2) then
            write (stderr, *) "test_combine_channels failed: Incorrect dimensions."
            test_combine_channels = .false.
            return
        end if
        if (any(abs(chol_comb(1:n1, 1:n1) - chol1) > tol) .or. &
            any(abs(chol_comb(n1 + 1:, n1 + 1:) - chol2) > tol) .or. &
            any(abs(chol_comb(1:n1, n1 + 1:)) > tol) .or. &
            any(abs(chol_comb(n1 + 1:, 1:n1)) > tol)) then
            write (stderr, *) "test_combine_channels failed: Incorrect "// &
                "block-diagonal Cholesky factor."
            test_combine_channels = .false.
        end if
        if (any(map_comb(1:n1) /= map1) .or. &
            any(map_comb(n1 + 1:) /= n_offset + map2)) then
            write (stderr, *) "test_combine_channels failed: Incorrect combined map."
            test_combine_channels = .false.
        end if

    end function test_combine_channels

    logical(c_bool) function test_get_ms_a_inv_cs() bind(C)
        !
        ! this function tests the subroutine which computes the pseudoinverse
        ! multisecant SR1 matrix for the closed-shell case, for both the linear and the
        ! non-linear part
        !
        use otr_arh, only: get_ms_a_inv_cs, arh_settings_type
        use otr_oao_test_reference, only: n_ao
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_particle = 1, n_diff = 3, n_accepted = 2, &
                                  flat_len = n_ao * n_ao * n_particle
        character(10), parameter :: case_names(2) = &
            [character(10) :: "linear", "non-linear"]

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    fock_diff(n_ao, n_ao, n_particle, n_diff), &
                    empty_dm_diff(n_ao, n_ao, n_particle, 0), &
                    empty_fock_diff(n_ao, n_ao, n_particle, 0), &
                    chol(n_accepted, n_accepted), flat(flat_len, n_diff), &
                    flat_v(flat_len, n_diff)
        real(rp), allocatable :: a_inv(:, :), expected(:, :), a_tilde(:, :), &
                                 y_gram(:, :)
        integer(ip) :: map(n_accepted), error, i_case
        logical :: linear
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_cs = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize random history and response, reordered and rebased into an 
        ! arbitrary orthonormalized basis that also drops one history entry; the 
        ! entries are centered on zero so that a history direction can end up poorly 
        ! aligned with its own response and be screened out
        call random_number(dm_diff)
        call random_number(fock_diff)
        dm_diff = dm_diff - 0.5_rp
        fock_diff = fock_diff - 0.5_rp
        map = [3, 1]
        chol = generate_random_upper_triangular(n_accepted)
        flat = reshape(fock_diff, [flat_len, n_diff])

        do i_case = 1, 2
            linear = i_case == 1

            ! the expected pseudoinverse is assembled independently of the routine
            a_tilde = ref_congruence_transform( &
                ref_build_a_part(dm_diff, fock_diff, linear), map, chol)
            y_gram = ref_congruence_transform(matmul(transpose(flat), flat), map, chol)
            expected = ref_ms_a_inv(a_tilde, y_gram)

            ! call routine and determine if the pseudoinverse matches
            call get_ms_a_inv_cs(dm_diff, fock_diff, linear, map, chol, a_inv, &
                                 settings, error)
            if (error /= 0) then
                write (stderr, *) "test_get_ms_a_inv_cs failed for the "// &
                    trim(case_names(i_case))//" part: Produced error."
                test_get_ms_a_inv_cs = .false.
            else if (norm2(a_inv - expected) > tol) then
                write (stderr, *) "test_get_ms_a_inv_cs failed for the "// &
                    trim(case_names(i_case))//" part: Incorrect pseudoinverse."
                test_get_ms_a_inv_cs = .false.
            end if
            deallocate(a_inv, expected, a_tilde, y_gram)
        end do

        ! a history whose responses live almost entirely outside the span of the steps 
        ! leaves every direction badly aligned with its own response, so the screening
        ! criterion has to discard all of them and the pseudoinverse has to vanish
        call random_number(flat)
        call random_number(flat_v)
        flat(flat_len / 2 + 1:, :) = 0.0_rp
        flat_v(:flat_len / 2, :) = 1e-3_rp * flat_v(:flat_len / 2, :)
        dm_diff = reshape(flat, shape(dm_diff))
        fock_diff = reshape(flat_v, shape(fock_diff))

        ! call routine and determine if the pseudoinverse vanishes
        call get_ms_a_inv_cs(dm_diff, fock_diff, .false., map, chol, a_inv, settings, &
                             error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Produced error for a "// &
                "badly aligned history."
            test_get_ms_a_inv_cs = .false.
        else if (norm2(a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Screening did not "// &
                "discard every direction of a badly aligned history."
            test_get_ms_a_inv_cs = .false.
        end if
        deallocate(a_inv)

        ! call routine for an empty history and determine if dimensions of the
        ! resulting pseudoinverse vanish
        call get_ms_a_inv_cs(empty_dm_diff, empty_fock_diff, .true., &
                             [integer(ip) ::], reshape([real(rp) ::], [0, 0]), a_inv, &
                             settings, error)
        if (size(a_inv, 1) /= 0 .or. size(a_inv, 2) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Incorrect "// &
                "pseudoinverse dimensions for empty history."
            test_get_ms_a_inv_cs = .false.
        end if
        deallocate(a_inv)

    end function test_get_ms_a_inv_cs

    logical(c_bool) function test_get_ms_a_inv_os_linear() bind(C)
        !
        ! this function tests the subroutine which computes the pseudoinverse 
        ! multisecant SR1 matrix in a spin-separated manner for the linear part in the 
        ! open-shell case
        !
        use otr_arh, only: get_ms_a_inv_os_linear, arh_settings_type
        use otr_oao_test_reference, only: n_ao, n_particle
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_diff = 2, n_col = n_particle * n_diff, &
                                  n_accepted = 3, n_ao2 = n_ao * n_ao

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_same_spin_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_opposite_spin_diff(n_ao, n_ao, n_particle, n_diff), &
                    empty_dm_diff(n_ao, n_ao, n_particle, 0), &
                    empty_v_diff(n_ao, n_ao, n_particle, 0), &
                    chol(n_accepted, n_accepted), s_full(n_particle * n_ao2, n_col), &
                    y_full(n_particle * n_ao2, n_col)
        real(rp), allocatable :: a_inv(:, :), expected(:, :), a_tilde(:, :), &
                                 y_gram(:, :)
        integer(ip) :: map(n_accepted), error
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_os_linear = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize random history and spin-resolved potentials centered on zero, 
        ! reordered and rebased into an arbitrary orthonormalized basis that also drops 
        ! one column
        call random_number(dm_diff)
        call random_number(v_same_spin_diff)
        call random_number(v_opposite_spin_diff)
        dm_diff = dm_diff - 0.5_rp
        v_same_spin_diff = v_same_spin_diff - 0.5_rp
        v_opposite_spin_diff = v_opposite_spin_diff - 0.5_rp
        map = [4, 1, 3]
        chol = generate_random_upper_triangular(n_accepted)

        ! the expected pseudoinverse is assembled independently of the routine from the
        ! stacked history and response, whose products give both the blocks of A and
        ! the response Gram matrix; A is symmetrized only after the transform, since
        ! its two halves are exactly transposes of one another in exact arithmetic
        call ref_stack_os_linear(dm_diff, v_same_spin_diff, v_opposite_spin_diff, &
                                 n_ao, s_full, y_full)
        a_tilde = ref_congruence_transform(matmul(transpose(s_full), y_full), map, chol)
        a_tilde = 0.5_rp * (a_tilde + transpose(a_tilde))
        y_gram = ref_congruence_transform(matmul(transpose(y_full), y_full), map, chol)
        expected = ref_ms_a_inv(a_tilde, y_gram)

        ! call routine and determine if the pseudoinverse matches
        call get_ms_a_inv_os_linear(dm_diff, v_same_spin_diff, v_opposite_spin_diff, &
                                    map, chol, a_inv, n_ao, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Produced error."
            test_get_ms_a_inv_os_linear = .false.
        else if (norm2(a_inv - expected) > tol) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "pseudoinverse."
            test_get_ms_a_inv_os_linear = .false.
        end if
        deallocate(a_inv, expected, a_tilde, y_gram)

        ! confining the history to the first atomic orbital row while suppressing the
        ! potentials there leaves every direction badly aligned with its own response,
        ! so the screening criterion has to discard all of them
        dm_diff(2:, :, :, :) = 0.0_rp
        v_same_spin_diff(:1, :, :, :) = 1e-3_rp * v_same_spin_diff(:1, :, :, :)
        v_opposite_spin_diff(:1, :, :, :) = 1e-3_rp * v_opposite_spin_diff(:1, :, :, :)

        ! call routine and determine if the pseudoinverse vanishes
        call get_ms_a_inv_os_linear(dm_diff, v_same_spin_diff, v_opposite_spin_diff, &
                                    map, chol, a_inv, n_ao, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Produced error "// &
                "for a badly aligned history."
            test_get_ms_a_inv_os_linear = .false.
        else if (norm2(a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Screening did "// &
                "not discard every direction of a badly aligned history."
            test_get_ms_a_inv_os_linear = .false.
        end if
        deallocate(a_inv)

        ! call routine for an empty history and determine if dimensions of the
        ! resulting pseudoinverse vanish
        call get_ms_a_inv_os_linear(empty_dm_diff, empty_v_diff, empty_v_diff, &
                                    [integer(ip) ::], reshape([real(rp) ::], [0, 0]), &
                                    a_inv, n_ao, settings, error)
        if (size(a_inv, 1) /= 0 .or. size(a_inv, 2) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "pseudoinverse dimensions for empty history."
            test_get_ms_a_inv_os_linear = .false.
        end if
        deallocate(a_inv)

    end function test_get_ms_a_inv_os_linear

    logical(c_bool) function test_get_ms_a_inv_os_nonlinear() bind(C)
        !
        ! this function tests the subroutine which computes the pseudoinverse
        ! multisecant SR1 matrix in a spin-combined manner for the non-linear part in 
        ! the open-shell case
        !
        use otr_arh, only: get_ms_a_inv_os_nonlinear, arh_settings_type
        use otr_oao_test_reference, only: n_ao, n_particle
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_diff = 3, n_accepted = 2, &
                                  flat_len = n_ao * n_ao * n_particle

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    v_diff(n_ao, n_ao, n_particle, n_diff), &
                    empty_dm_diff(n_ao, n_ao, n_particle, 0), &
                    empty_v_diff(n_ao, n_ao, n_particle, 0), &
                    chol(n_accepted, n_accepted), flat(flat_len, n_diff), &
                    flat_v(flat_len, n_diff)
        real(rp), allocatable :: a_inv(:, :), expected(:, :), a_tilde(:, :), &
                                 y_gram(:, :)
        integer(ip) :: map(n_accepted), error
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_os_nonlinear = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize random history and response centered on zero, reordered and 
        ! rebased into an arbitrary orthonormalized basis that also drops one history 
        ! entry
        call random_number(dm_diff)
        call random_number(v_diff)
        dm_diff = dm_diff - 0.5_rp
        v_diff = v_diff - 0.5_rp
        map = [3, 1]
        chol = generate_random_upper_triangular(n_accepted)
        flat_v = reshape(v_diff, [flat_len, n_diff])

        ! the expected pseudoinverse is assembled independently of the routine
        a_tilde = ref_congruence_transform( &
            ref_build_a_part(dm_diff, v_diff, .false.), map, chol)
        y_gram = ref_congruence_transform(matmul(transpose(flat_v), flat_v), map, chol)
        expected = ref_ms_a_inv(a_tilde, y_gram)

        ! call routine and determine if the pseudoinverse matches
        call get_ms_a_inv_os_nonlinear(dm_diff, v_diff, map, chol, a_inv, settings, &
                                       error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Produced error."
            test_get_ms_a_inv_os_nonlinear = .false.
        else if (norm2(a_inv - expected) > tol) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Incorrect "// &
                "pseudoinverse."
            test_get_ms_a_inv_os_nonlinear = .false.
        end if
        deallocate(a_inv, expected, a_tilde, y_gram)

        ! a history whose responses live almost entirely outside the span of the steps
        ! leaves every direction badly aligned with its own response, so the screening
        ! criterion has to discard all of them and the pseudoinverse has to vanish
        call random_number(flat)
        call random_number(flat_v)
        flat(flat_len / 2 + 1:, :) = 0.0_rp
        flat_v(:flat_len / 2, :) = 1e-3_rp * flat_v(:flat_len / 2, :)
        dm_diff = reshape(flat, shape(dm_diff))
        v_diff = reshape(flat_v, shape(v_diff))

        ! call routine and determine if the pseudoinverse vanishes
        call get_ms_a_inv_os_nonlinear(dm_diff, v_diff, map, chol, a_inv, settings, &
                                       error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Produced "// &
                "error for a badly aligned history."
            test_get_ms_a_inv_os_nonlinear = .false.
        else if (norm2(a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Screening "// &
                "did not discard every direction of a badly aligned history."
            test_get_ms_a_inv_os_nonlinear = .false.
        end if
        deallocate(a_inv)

        ! call routine for an empty history and determine if dimensions of the
        ! resulting pseudoinverse vanish
        call get_ms_a_inv_os_nonlinear(empty_dm_diff, empty_v_diff, [integer(ip) ::], &
                                       reshape([real(rp) ::], [0, 0]), a_inv, &
                                       settings, error)
        if (size(a_inv, 1) /= 0 .or. size(a_inv, 2) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Incorrect "// &
                "pseudoinverse dimensions for empty history."
            test_get_ms_a_inv_os_nonlinear = .false.
        end if
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

    logical(c_bool) function test_apply_ms_sr1_skip() bind(C)
        !
        ! this function tests the subroutine which discards history directions failing
        ! the multisecant SR1 skipping criterion
        !
        use otr_arh, only: apply_ms_sr1_skip, ms_sr1_skip_thresh
        use otr_oao_unit_tests, only: identity_matrix

        integer(ip), parameter :: n = 4
        integer(ip) :: i
        real(rp) :: eig_vals(n), eig_vecs(n, n), y_gram(n, n), eig_vals_inv(n), &
                    expected(n), y_norm(n), factor(n), b(n, n), theta, c, s

        ! assume tests pass
        test_apply_ms_sr1_skip = .true.

        ! each eigenvalue is placed at a fixed multiple of the criterion applied to its
        ! own response norm, so which directions are discarded does not depend on the
        ! value of the skipping threshold; the multiples sit just either side of one so
        ! that a response norm computed even slightly wrongly flips a decision
        factor = [-0.95_rp, 0.95_rp, 1.05_rp, 1.05_rp]

        ! a symmetric positive semi-definite response Gram matrix, built as B^T B so
        ! that v^T y_gram v is a genuine squared response norm, scaled so that the norm
        ! is clearly different from its own square
        call random_number(b)
        y_gram = 4.0_rp * matmul(transpose(b), b)

        ! a Givens rotation in the (1, 3) plane, leaving directions 2 and 4 along the
        ! axes, so that the quadratic form is a plain diagonal entry for some
        ! directions and a genuine mixture for the others
        eig_vecs = identity_matrix(n)
        theta = 0.7_rp
        c = cos(theta)
        s = sin(theta)
        eig_vecs(1, 1) = c
        eig_vecs(3, 1) = s
        eig_vecs(1, 3) = -s
        eig_vecs(3, 3) = c

        ! scale each eigenvalue against its own response norm
        do i = 1, n
            y_norm(i) = sqrt(max(dot_product(eig_vecs(:, i), &
                                             matmul(y_gram, eig_vecs(:, i))), 0.0_rp))
            eig_vals(i) = factor(i) * ms_sr1_skip_thresh * y_norm(i)
        end do

        ! the exact inverse before skipping, and the expected result afterwards
        eig_vals_inv = 1.0_rp / eig_vals
        expected = eig_vals_inv
        expected(1) = 0.0_rp
        expected(2) = 0.0_rp

        ! call routine and determine if the surviving inverse eigenvalues match
        call apply_ms_sr1_skip(eig_vals, eig_vecs, y_gram, eig_vals_inv)
        if (norm2(eig_vals_inv - expected) > tol) then
            write (stderr, *) "test_apply_ms_sr1_skip failed: Incorrect inverse "// &
                "eigenvalues after skipping."
            test_apply_ms_sr1_skip = .false.
        end if

    end function test_apply_ms_sr1_skip

    logical(c_bool) function test_response_gram() bind(C)
        !
        ! this function tests the function which returns the Gram matrix of the
        ! response history rebased into the orthonormalized S-basis
        !
        use otr_arh, only: response_gram
        use otr_oao_test_reference, only: n_ao, n_particle

        integer(ip), parameter :: n_dm = 3, n_accepted = 2
        integer(ip) :: map(n_accepted), i, j
        real(rp) :: v_diff(n_ao, n_ao, n_particle, n_dm), &
                    chol(n_accepted, n_accepted), &
                    flat(n_ao * n_ao * n_particle, n_dm), gram(n_dm, n_dm)
        real(rp), allocatable :: y_gram(:, :), expected(:, :)

        ! assume tests pass
        test_response_gram = .true.

        ! random response history and an upper-triangular Cholesky factor with a
        ! non-trivial map that both selects a subset and reorders it
        call random_number(v_diff)
        map = [3_ip, 1_ip]
        chol = generate_random_upper_triangular(n_accepted)

        ! build the expected Gram matrix independently from the flattened history
        flat = reshape(v_diff, [n_ao * n_ao * n_particle, n_dm])
        gram = matmul(transpose(flat), flat)
        expected = ref_congruence_transform(gram, map, chol)

        ! call routine and determine if the rebased Gram matrix matches
        y_gram = response_gram(v_diff, map, chol)
        if (size(y_gram, 1) /= n_accepted .or. size(y_gram, 2) /= n_accepted) then
            write (stderr, *) "test_response_gram failed: Incorrect shape."
            test_response_gram = .false.
        else if (maxval(abs(y_gram - expected)) > tol) then
            write (stderr, *) "test_response_gram failed: Incorrect rebased "// &
                "response Gram matrix."
            test_response_gram = .false.
        end if

    end function test_response_gram

    logical(c_bool) function test_response_gram_os_linear() bind(C)
        !
        ! this function tests the function which returns the Gram matrix of the
        ! open-shell linear response history, whose same-spin and opposite-spin
        ! potentials are interleaved exactly as the rows of A pair them
        !
        use otr_arh, only: response_gram_os_linear
        use otr_oao_test_reference, only: n_ao, n_particle

        integer(ip), parameter :: n_dm = 2, n_accepted = 3, n_ao2 = n_ao * n_ao
        integer(ip) :: map(n_accepted), k
        real(rp) :: v_same(n_ao, n_ao, n_particle, n_dm), &
                    v_opp(n_ao, n_ao, n_particle, n_dm), chol(n_accepted, n_accepted), &
                    y_full(2 * n_ao2, 2 * n_dm), gram(2 * n_dm, 2 * n_dm)
        real(rp), allocatable :: y_gram(:, :), expected(:, :)

        ! assume tests pass
        test_response_gram_os_linear = .true.

        ! random same-spin and opposite-spin response histories
        call random_number(v_same)
        call random_number(v_opp)
        map = [4_ip, 1_ip, 3_ip]
        chol = generate_random_upper_triangular(n_accepted)

        ! stack the alpha and beta blocks independently, in the interleaved column
        ! order the routine documents
        do k = 1, n_dm
            y_full(:n_ao2, k) = reshape(v_same(:, :, 1, k), [n_ao2])
            y_full(n_ao2 + 1:, k) = reshape(v_opp(:, :, 2, k), [n_ao2])
            y_full(:n_ao2, n_dm + k) = reshape(v_opp(:, :, 1, k), [n_ao2])
            y_full(n_ao2 + 1:, n_dm + k) = reshape(v_same(:, :, 2, k), [n_ao2])
        end do
        gram = matmul(transpose(y_full), y_full)
        expected = ref_congruence_transform(gram, map, chol)

        ! call routine and determine if the rebased Gram matrix matches
        y_gram = response_gram_os_linear(v_same, v_opp, n_ao, map, chol)
        if (size(y_gram, 1) /= n_accepted .or. size(y_gram, 2) /= n_accepted) then
            write (stderr, *) "test_response_gram_os_linear failed: Incorrect shape."
            test_response_gram_os_linear = .false.
        else if (maxval(abs(y_gram - expected)) > tol) then
            write (stderr, *) "test_response_gram_os_linear failed: Incorrect "// &
                "rebased open-shell linear response Gram matrix."
            test_response_gram_os_linear = .false.
        end if

    end function test_response_gram_os_linear

end module otr_arh_unit_tests
