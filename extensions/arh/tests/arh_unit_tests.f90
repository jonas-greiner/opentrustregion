! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_density_qn_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use otr_common_unit_tests, only: mock_orbital_basis_type
    use otr_density_qn, only: density_qn_mo_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! multipliers of the density matrix returned by the mock density matrix evaluating
    ! functions, which differ between the first and any subsequent call, following the
    ! same convention as the shared multiplier for the Fock matrix; the Coulomb and
    ! exact-exchange multipliers are positive, so that both mock kernels are positive
    ! definite
    real(rp), parameter :: mock_v_coulomb_factor(2) = [3.0_rp, 7.0_rp], &
                           mock_v_exchange_factor(2) = [1.0_rp, 2.0_rp]

    ! the non-linear potential keeps a single multiplier across calls, unlike the
    ! potentials above, so that it stays one consistent function of the density, which
    ! is what the history screening measures and would otherwise reject as noise
    real(rp), parameter :: mock_v_nonlinear_factor = 6.0_rp

    ! multipliers the mock density matrix evaluating functions apply on their first
    ! call to the potentials the density QN object keeps, in the order of
    ! stacked_potentials
    real(rp), parameter :: mock_potential_factors(3) = &
        [mock_v_coulomb_factor(1), mock_v_exchange_factor(1), mock_v_nonlinear_factor]

    ! orbital basis which records what the roles a density QN object hands on to its
    ! orbitals receive and returns three times its input with the error code it is set
    ! to, for the tests of these roles
    type, extends(mock_orbital_basis_type) :: recording_basis_type
        real(rp), allocatable :: received(:)
        class(*), pointer :: received_context => null()
        integer(ip) :: error = 0
    contains
        procedure :: precond_pd => recording_basis_precond_pd
        procedure :: project => recording_basis_project
        procedure :: get_extra_trial_vectors => recording_basis_extra_trial_vectors
    end type

    ! density QN object in the MO basis which records the calls of the inverse of its
    ! approximate Hessian and of the builders of its model and returns three times the
    ! vector it receives with the error code it is set to, for the tests of the
    ! routines which hand on to them
    type, extends(density_qn_mo_type) :: recording_density_qn_type
        real(rp), allocatable :: received(:)
        logical :: level_shift_received = .false.
        real(rp) :: received_level_shift = 0.0_rp
        character(len=2) :: built_shell = ""
        integer(ip) :: error = 0
    contains
        procedure :: inv_hess_x => recording_inv_hess_x
        procedure :: build_hess_model_cs => recording_build_hess_model_cs
        procedure :: build_hess_model_os => recording_build_hess_model_os
    end type

    ! dispatch the mock potential by density matrix rank
    interface mock_potential
        module procedure mock_potential_cs, mock_potential_os
    end interface

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

    subroutine mock_evaluate_dm_cs(dm, energy, fock, v_coulomb, v_exchange, &
                                   v_nonlinear, error, context)
        !
        ! this subroutine is a mock density matrix evaluating function with Coulomb,
        ! exact-exchange and non-linear potential contributions for the closed-shell
        ! case, which returns multiples of the density matrix that change between calls
        ! so that non-vanishing differences are produced
        !
        use otr_common_unit_tests, only: record_mock_call, mock_factor, mock_fock_factor
        use otr_common_mock, only: record_host_call

        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target, contiguous :: &
            fock(:, :), v_coulomb(:, :), v_exchange(:, :), v_nonlinear(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! check host context and fail with the error code set in it
        error = record_host_call("evaluate_dm", context)
        if (error /= 0) then
            energy = sum(dm)
            if (present(fock)) fock = 0.0_rp
            if (present(v_coulomb)) v_coulomb = 0.0_rp
            if (present(v_exchange)) v_exchange = 0.0_rp
            if (present(v_nonlinear)) v_nonlinear = 0.0_rp
            return
        end if

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + merge( &
            2_ip, 0_ip, present(v_coulomb)) + merge(4_ip, 0_ip, present(v_exchange)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))

        energy = sum(dm)
        if (present(fock)) fock = mock_potential(mock_factor(mock_fock_factor), dm)
        if (present(v_coulomb)) &
            v_coulomb = mock_potential(mock_factor(mock_v_coulomb_factor), dm)
        if (present(v_exchange)) &
            v_exchange = mock_potential(mock_factor(mock_v_exchange_factor), dm)
        if (present(v_nonlinear)) &
            v_nonlinear = mock_potential(mock_v_nonlinear_factor, dm)

    end subroutine mock_evaluate_dm_cs

    subroutine mock_evaluate_dm_os(dm, energy, fock, v_coulomb, v_exchange, &
                                   v_nonlinear, error, context)
        !
        ! this subroutine is a mock density matrix evaluating function with Coulomb,
        ! exact-exchange and non-linear potential contributions for the open-shell case,
        ! which returns multiples of the density matrix that change between calls so
        ! that non-vanishing differences are produced
        !
        use otr_common_unit_tests, only: record_mock_call, mock_factor, mock_fock_factor
        use otr_common_mock, only: record_host_call

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target :: &
            fock(:, :, :), v_coulomb(:, :, :), v_exchange(:, :, :), v_nonlinear(:, :, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: i

        ! check host context and fail with the error code set in it
        error = record_host_call("evaluate_dm", context)
        if (error /= 0) then
            energy = sum(dm)
            if (present(fock)) fock = 0.0_rp
            if (present(v_coulomb)) v_coulomb = 0.0_rp
            if (present(v_exchange)) v_exchange = 0.0_rp
            if (present(v_nonlinear)) v_nonlinear = 0.0_rp
            return
        end if

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + merge( &
            2_ip, 0_ip, present(v_coulomb)) + merge(4_ip, 0_ip, present(v_exchange)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))

        energy = sum(dm)
        do i = 1, size(dm, 3)
            if (present(fock)) fock(:, :, i) = &
                mock_potential(mock_factor(mock_fock_factor), dm(:, :, i))
            if (present(v_coulomb)) v_coulomb(:, :, i) = &
                mock_potential(mock_factor(mock_v_coulomb_factor), dm(:, :, i))
            if (present(v_exchange)) v_exchange(:, :, i) = &
                mock_potential(mock_factor(mock_v_exchange_factor), dm(:, :, i))
            if (present(v_nonlinear)) v_nonlinear(:, :, i) = &
                mock_potential(mock_v_nonlinear_factor, dm(:, :, i))
        end do

    end subroutine mock_evaluate_dm_os

    subroutine recording_basis_precond_pd(self, residual, precond_residual, error, &
                                          context)
        !
        ! this subroutine records the residual and the context the positive-definite
        ! preconditioner of the recording orbital basis receives
        !
        class(recording_basis_type), intent(inout) :: self
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        self%received = residual
        self%received_context => context
        precond_residual = 3.0_rp * residual
        error = self%error

    end subroutine recording_basis_precond_pd

    subroutine recording_basis_project(self, vector, error, context)
        !
        ! this subroutine records the vector and the context the projection of the
        ! recording orbital basis receives
        !
        class(recording_basis_type), intent(inout) :: self
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        self%received = vector
        self%received_context => context
        vector = 3.0_rp * vector
        error = self%error

    end subroutine recording_basis_project

    subroutine recording_basis_extra_trial_vectors(self, trial_vectors, error, context)
        !
        ! this subroutine records the context the extra trial vector function of the
        ! recording orbital basis receives and returns trial vectors of threes
        !
        class(recording_basis_type), intent(inout) :: self
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        self%received_context => context
        trial_vectors = 3.0_rp
        error = self%error

    end subroutine recording_basis_extra_trial_vectors

    subroutine recording_inv_hess_x(self, x, inv_hess_x, error, level_shift)
        !
        ! this subroutine records the vector and the level shift the inverse of the
        ! approximate Hessian of the recording density QN object receives
        !
        class(recording_density_qn_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: inv_hess_x(:)
        integer(ip), intent(out) :: error
        real(rp), intent(in), optional :: level_shift

        self%received = x
        self%level_shift_received = present(level_shift)
        if (present(level_shift)) self%received_level_shift = level_shift
        inv_hess_x = 3.0_rp * x
        error = self%error

    end subroutine recording_inv_hess_x

    subroutine recording_build_hess_model_cs(self, error)
        !
        ! this subroutine records that the closed-shell model builder of the recording
        ! density QN object was called
        !
        class(recording_density_qn_type), intent(inout) :: self
        integer(ip), intent(out) :: error

        self%built_shell = "cs"
        error = self%error

    end subroutine recording_build_hess_model_cs

    subroutine recording_build_hess_model_os(self, error)
        !
        ! this subroutine records that the open-shell model builder of the recording
        ! density QN object was called
        !
        class(recording_density_qn_type), intent(inout) :: self
        integer(ip), intent(out) :: error

        self%built_shell = "os"
        error = self%error

    end subroutine recording_build_hess_model_os

    subroutine generate_random_fock_partition(n_ao, n_occ, dm_oao, fock_oo, fock_vv)
        !
        ! this subroutine generates a random density matrix with the given number of
        ! occupied orbitals per particle/spin channel and the corresponding
        ! occupied-occupied and virtual-virtual Fock matrix blocks
        !
        use opentrustregion_unit_tests, only: identity_matrix, &
                                              generate_random_symm_matrix
        use otr_common_unit_tests, only: generate_random_density_matrix

        integer(ip), intent(in) :: n_ao, n_occ(:)
        real(rp), intent(out) :: dm_oao(:, :, :), fock_oo(:, :, :), fock_vv(:, :, :)

        real(rp) :: full_fock(n_ao, n_ao), complement(n_ao, n_ao)
        integer(ip) :: j

        do j = 1, size(n_occ, kind=ip)
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_occ(j))
            full_fock = generate_random_symm_matrix(n_ao)
            complement = identity_matrix(n_ao) - dm_oao(:, :, j)
            fock_oo(:, :, j) = matmul(dm_oao(:, :, j), &
                                      matmul(full_fock, dm_oao(:, :, j)))
            fock_vv(:, :, j) = matmul(complement, matmul(full_fock, complement))
        end do

    end subroutine generate_random_fock_partition

    function generate_random_cholesky_factor(n) result(matrix)
        !
        ! this function generates a random upper triangular matrix whose diagonal is
        ! bounded away from zero, so that it is invertible and well enough conditioned
        ! to stand in for the Cholesky factor a history factorization would produce
        !
        use opentrustregion_unit_tests, only: identity_matrix, &
                                              generate_random_upper_triangular_matrix

        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        matrix = generate_random_upper_triangular_matrix(n) + identity_matrix(n)

    end function generate_random_cholesky_factor

    function generate_random_psd_matrix(n, rank) result(matrix)
        !
        ! this function generates a random symmetric positive semidefinite matrix of
        ! the given rank as B B^T
        !
        integer(ip), intent(in) :: n, rank
        real(rp) :: matrix(n, n)

        real(rp) :: b(n, rank)

        call random_number(b)
        b = b - 0.5_rp
        matrix = matmul(b, transpose(b))

    end function generate_random_psd_matrix

    subroutine setup_density_qn_and_mo_objects(density_qn, mo, mo_coeff, ao_overlap, &
                                               n_occ, context)
        !
        ! this subroutine sets up the given MO and density QN objects for orbitals
        ! parameterized in the MO basis the way the MO and density QN factories would,
        ! pointing to the test-local MO coefficients
        !
        use otr_common, only: error_density_qn
        use otr_density_qn, only: density_qn_type, density_qn_mo_type
        use otr_mo, only: mo_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_mo_unit_tests, only: setup_mo_object

        class(density_qn_type), allocatable, intent(out) :: density_qn
        type(mo_type), allocatable, intent(out), target :: mo
        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:)
        type(test_context_type), intent(inout), target :: context

        ! set up the MO object
        call setup_mo_object(mo, mo_coeff, ao_overlap, n_occ)

        ! set up the density QN object in the MO basis
        allocate(density_qn, source=density_qn_mo_type(mo))
        density_qn%error_origin = error_density_qn
        call setup_settings(density_qn%settings, context)

    end subroutine setup_density_qn_and_mo_objects

    subroutine setup_density_qn_and_scaled_oao_objects(density_qn, oao, basis_scale, &
                                                       dm_ao, dm_oao, context)
        !
        ! this subroutine sets up the given OAO object for random density matrices in
        ! an AO basis whose overlap has the inverse square root c I, so that the
        ! transformation to the OAO basis scales by c^2, and the given density QN
        ! object pointing to it the way the density QN factory would
        !
        use otr_common, only: error_oao, error_density_qn
        use otr_density_qn, only: density_qn_type, density_qn_oao_type
        use otr_oao, only: oao_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix
        use otr_common_unit_tests, only: generate_random_density_matrix
        use otr_common_test_reference, only: n_occ

        class(density_qn_type), allocatable, intent(out) :: density_qn
        type(oao_type), allocatable, intent(out), target :: oao
        real(rp), intent(in) :: basis_scale
        real(rp), intent(out), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(out) :: dm_oao(:, :, :)
        type(test_context_type), intent(inout), target :: context

        integer(ip) :: n_ao, n_particle, n_param, i

        ! set up the OAO object
        n_ao = size(dm_ao, 1, kind=ip)
        n_particle = size(dm_ao, 3, kind=ip)
        n_param = n_particle * n_ao * (n_ao - 1) / 2
        allocate(oao)
        call setup_settings(oao%settings, context)
        oao%error_origin = error_oao
        oao%n_ao = n_ao
        oao%n_particle = n_particle
        oao%n_param = n_param
        oao%s_inv_sqrt = basis_scale * identity_matrix(n_ao)
        do i = 1, n_particle
            dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do
        dm_ao = basis_scale**2 * dm_oao
        oao%dm_ao => dm_ao
        oao%dm_oao = dm_oao
        allocate(oao%fock_oo(n_ao, n_ao, n_particle), &
                 oao%fock_vv(n_ao, n_ao, n_particle), oao%grad(n_param), &
                 oao%h_diag(n_param))

        ! set up the density QN object
        allocate(density_qn, source=density_qn_oao_type(oao))
        density_qn%error_origin = error_density_qn
        call setup_settings(density_qn%settings, context)

    end subroutine setup_density_qn_and_scaled_oao_objects

    subroutine setup_empty_history(density_qn, n_particle)
        !
        ! this subroutine sets up the given density QN object with an empty history and
        ! vanishing potentials at the current point for the closed- or open-shell case
        !
        use otr_density_qn, only: density_qn_type

        class(density_qn_type), intent(inout) :: density_qn
        integer(ip), intent(in) :: n_particle

        integer(ip) :: n_ao

        n_ao = density_qn%orbitals%n_ao
        allocate(density_qn%dm_list(n_ao, n_ao, n_particle, 0), &
                 density_qn%v_coulomb_list(n_ao, n_ao, n_particle, 0), &
                 density_qn%v_exchange_list(n_ao, n_ao, n_particle, 0), &
                 density_qn%v_nonlinear_list(n_ao, n_ao, n_particle, 0))
        allocate(density_qn%v_coulomb(n_ao, n_ao, n_particle), &
                 density_qn%v_exchange(n_ao, n_ao, n_particle), &
                 density_qn%v_nonlinear(n_ao, n_ao, n_particle))
        density_qn%v_coulomb = 0.0_rp
        density_qn%v_exchange = 0.0_rp
        density_qn%v_nonlinear = 0.0_rp
        if (n_particle == 1) then
            density_qn%evaluate_dm_cs => mock_evaluate_dm_cs
        else
            density_qn%evaluate_dm_os => mock_evaluate_dm_os
        end if

    end subroutine setup_empty_history

    function stacked_potentials(density_qn, k) result(potentials)
        !
        ! this function stacks the potentials the given density QN object keeps along
        ! the last dimension, the current ones for k = 0 and those of the k-th history
        ! entry otherwise
        !
        use otr_density_qn, only: density_qn_type

        class(density_qn_type), intent(in) :: density_qn
        integer(ip), intent(in) :: k
        real(rp), allocatable :: potentials(:, :, :, :)

        if (k == 0) then
            potentials = reshape( &
                [density_qn%v_coulomb, density_qn%v_exchange, density_qn%v_nonlinear], &
                [shape(density_qn%v_nonlinear), 3])
        else
            potentials = reshape([density_qn%v_coulomb_list(:, :, :, k), &
                                  density_qn%v_exchange_list(:, :, :, k), &
                                  density_qn%v_nonlinear_list(:, :, :, k)], &
                                 [shape(density_qn%v_nonlinear_list(:, :, :, k)), 3])
        end if

    end function stacked_potentials

    function ref_build_a_part(dm_cols, v_cols) result(a)
        !
        ! this function independently reproduces the raw product A = S^T Y of the
        ! history columns, together with its symmetrization
        !
        real(rp), intent(in) :: dm_cols(:, :), v_cols(:, :)
        real(rp), allocatable :: a(:, :)

        real(rp), allocatable :: raw(:, :)

        raw = matmul(transpose(dm_cols), v_cols)
        a = 0.5_rp * (raw + transpose(raw))
        deallocate(raw)

    end function ref_build_a_part

    function ref_ms_a_inv(a_tilde, y_gram) result(a_inv)
        !
        ! this function independently reproduces the screened pseudoinverse the
        ! multisecant routines build from a congruence-transformed A: eigenvalues at
        ! the level of numerical noise are discarded, as are directions whose
        ! eigenvalue is small against the norm of the response they divide
        !
        use otr_density_qn, only: eig_val_noise_factor, ms_sr1_skip_thresh

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

    function ref_psd_pinv(a_tilde, scale) result(a_inv)
        !
        ! this function independently reproduces the pseudoinverse of a matrix positive
        ! semidefinite in exact arithmetic: eigenvalues above the noise floor relative
        ! to the given spectral scale are inverted, all others discarded
        !
        use otr_density_qn, only: eig_val_noise_factor

        real(rp), intent(in) :: a_tilde(:, :), scale
        real(rp), allocatable :: a_inv(:, :)

        integer(ip) :: n, i, info, lwork
        real(rp), allocatable :: vecs(:, :), vals(:), diag(:, :), work(:)
        external :: dsyev

        n = size(a_tilde, 1)
        allocate(vecs(n, n), vals(n), diag(n, n))
        vecs = 0.5_rp * (a_tilde + transpose(a_tilde))
        lwork = 3_ip * n
        allocate(work(lwork))
        call dsyev("V", "U", n, vecs, n, vals, work, lwork, info)
        deallocate(work)

        diag = 0.0_rp
        do i = 1, n
            if (vals(i) > eig_val_noise_factor * epsilon(1.0_rp) * scale) &
                diag(i, i) = 1.0_rp / vals(i)
        end do
        a_inv = matmul(vecs, matmul(diag, transpose(vecs)))

    end function ref_psd_pinv

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

    function check_build_hess_model(n_particle, test_name) result(passed)
        !
        ! this function checks the assembly of the approximate Hessian model from the
        ! history for the closed- or open-shell case: standard ARH has to reproduce any
        ! history in the MO basis exactly and every other QN method one that symmetric
        ! linear Coulomb and exact-exchange kernels produce, subspace-projected
        ! multisecant and symmetrized ARH only in the closed-shell case, where its
        ! potential differences lie in the span of the density matrix differences,
        ! while every QN method in the OAO basis has to keep the whole history in the
        ! linear system, drop the entries beyond the window of the trust radius from
        ! the non-linear system, keep them without a state, leave the non-linear system
        ! out without a non-linear response and cache every quantity in the basis of
        ! its own system, without evaluating the density matrix or changing the history
        !
        use opentrustregion, only: solver_state_type
        use otr_density_qn, only: density_qn_type, qn_methods, density_qn_oao_type
        use otr_mo, only: mo_type
        use otr_oao, only: oao_type
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix
        use otr_common_unit_tests, only: mock_requests, &
                                         generate_random_density_matrix, &
                                         generate_random_orthogonal_matrix
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, ref_mo_transform, &
                                     ref_pack_ov, ref_unpack_ov

        type(test_context_type), target :: context
        type(oao_type), allocatable, target :: oao
        type(mo_type), allocatable, target :: mo
        class(density_qn_type), allocatable, target :: density_qn
        type(solver_state_type), target :: state
        integer(ip), intent(in) :: n_particle
        character(len=*), intent(in) :: test_name
        logical :: passed

        ! trust radius of every case of the OAO basis, of which the third has no state
        integer(ip), parameter :: n_hist = 3, n_cases = 4
        real(rp), parameter :: case_trust_radii(n_cases) = &
            [1.0_rp, 1e-4_rp, 0.0_rp, 1.0_rp]

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), step_mo(n_mo, n_mo, n_particle), &
                    v_diff_mo(n_mo, n_mo, n_particle), perturbation(n_ao, n_ao), &
                    dm_step(n_ao, n_ao), dm_hist(n_ao, n_ao, n_particle)
        real(rp), allocatable :: steps(:, :), orth(:, :), response(:), &
                                 dm_list(:, :, :, :), current(:, :, :, :), &
                                 lists(:, :, :, :, :)
        integer(ip) :: n_pot, n_list, n_linear, n_nonlinear, n_systems, n_terms, i, j, &
                       k, i_case, i_method, error
        logical :: built
        character(len=:), allocatable :: case_qn_method, case_name

        ! assume test passes
        passed = .true.

        ! the Coulomb and exact-exchange potentials are kept besides the non-linear
        ! potential; the linear potential of a channel is the Coulomb potential of all
        ! channels minus its exact-exchange potential
        n_pot = 3

        ! set up the density QN object in the MO basis at random orbitals, whose MOs
        ! span the same space for both channels as they do in practice, so that the
        ! history inner product of the two channels agrees and a Coulomb kernel acting
        ! on both channels is symmetric
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        if (n_particle == 2) mo_coeff(:, :, 2) = &
            matmul(mo_coeff(:, :, 1), generate_random_orthogonal_matrix(n_mo))
        call setup_density_qn_and_mo_objects(density_qn, mo, mo_coeff, ao_overlap, &
                                             n_occ(:n_particle), context)

        ! current density matrix, multiplied by the AO overlap matrix from both sides
        do j = 1, n_particle
            dm_hist(:, :, j) = matmul(ao_overlap, matmul(mo%dm_ao(:, :, j), ao_overlap))
        end do

        ! steps of every channel which are not orthogonal, so that the pivoted
        ! factorization of the history reorders them and its factor is not diagonal
        allocate(steps(mo%n_param, n_hist))
        k = 0
        do j = 1, n_particle
            i = n_occ(j) * (n_mo - n_occ(j))
            orth = generate_random_orthogonal_matrix(i)
            steps(k + 1:k + i, 1) = orth(:, 1)
            steps(k + 1:k + i, 2) = orth(:, 1) + 0.3_rp * orth(:, 2)
            steps(k + 1:k + i, 3) = orth(:, 3) + 0.5_rp * orth(:, 1)
            k = k + i
        end do

        ! history of the density matrices of these steps, stored multiplied by the AO
        ! overlap matrix from both sides, with random current and history potentials of
        ! the linear part and a non-linear potential which does not change
        allocate(dm_list(n_ao, n_ao, n_particle, n_hist), &
                 current(n_ao, n_ao, n_particle, n_pot), &
                 lists(n_ao, n_ao, n_particle, n_hist, n_pot))
        do j = 1, n_particle
            do i = 1, n_pot
                current(:, :, j, i) = generate_random_symm_matrix(n_ao)
            end do
        end do
        do k = 1, n_hist
            step_mo = ref_unpack_ov(steps(:, k), mo%mo_channels)
            do j = 1, n_particle
                dm_list(:, :, j, k) = matmul(ao_overlap, matmul( &
                    mo%dm_ao(:, :, j) + matmul(mo_coeff(:, :, j), matmul( &
                        step_mo(:, :, j) + transpose(step_mo(:, :, j)), &
                        transpose(mo_coeff(:, :, j)))), ao_overlap))
                do i = 1, n_pot - 1
                    lists(:, :, j, k, i) = generate_random_symm_matrix(n_ao)
                end do
            end do
            lists(:, :, :, k, n_pot) = current(:, :, :, n_pot)
        end do

        ! build the model for every QN method and determine if the low-rank part
        ! reproduces the linear potential difference of every step, which pins the
        ! scaling of the low-rank part, the matrices it couples the directions with and
        ! the (map, chol) pair its directions are rebased with
        do i_method = 1, size(qn_methods)
            case_qn_method = trim(qn_methods(i_method))
            density_qn%settings%qn_method = case_qn_method

            ! subspace-projected multisecant and symmetrized ARH only reproduce
            ! potential differences in the span of the density matrix differences,
            ! which the Coulomb potential of the other channel leaves in the open-shell
            ! case
            if ((case_qn_method == "ms_sp" .or. case_qn_method == "symm_arh") .and. &
                n_particle == 2) cycle

            ! every method but standard ARH only reproduces a history that symmetric
            ! linear kernels produce: Coulomb and exact-exchange potentials of every
            ! channel that are positive multiples of its density matrix difference, so
            ! that the Coulomb potentials of both channels act on every channel
            if (case_qn_method /= "arh") then
                do k = 1, n_hist
                    do j = 1, n_particle
                        dm_step = dm_list(:, :, j, k) - dm_hist(:, :, j)
                        lists(:, :, j, k, 1) = current(:, :, j, 1) + 3.0_rp * dm_step
                        lists(:, :, j, k, 2) = current(:, :, j, 2) + dm_step
                    end do
                end do
            end if

            call build_from_history()
            if (error /= 0) then
                write(stderr, *) test_name//" failed: Produced error for "// &
                    case_qn_method//" in the MO basis."
                passed = .false.
            end if
            if (density_qn%model_stale) then
                write(stderr, *) test_name//" failed: Hessian model still marked "// &
                    "stale for "//case_qn_method//" in the MO basis."
                passed = .false.
            end if
            if (.not. allocated(density_qn%coupling_matrix)) then
                write(stderr, *) test_name//" failed: No low-rank part for "// &
                    case_qn_method//" in the MO basis."
                passed = .false.
                cycle
            end if
            do k = 1, n_hist
                response = matmul(density_qn%expansion_dirs, matmul( &
                    density_qn%coupling_matrix, &
                    matmul(transpose(density_qn%projection_dirs), steps(:, k))))
                do j = 1, n_particle
                    v_diff_mo(:, :, j) = ref_mo_transform(mo_coeff(:, :, j), sum( &
                        lists(:, :, :, k, 1) - current(:, :, :, 1), &
                        dim=3) - (lists(:, :, j, k, 2) - current(:, :, j, 2)))
                end do
                if (norm2(response - merge(4.0_rp, 2.0_rp, n_particle == 1) * &
                          ref_pack_ov(v_diff_mo, mo%mo_channels)) > &
                    tol * (1.0_rp + norm2(response))) then
                    write(stderr, *) test_name//" failed: History not reproduced "// &
                        "exactly for "//case_qn_method//" in the MO basis."
                    passed = .false.
                end if
            end do
        end do
        deallocate(density_qn, mo, dm_list, current, lists)

        ! set up the OAO object with an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        allocate(oao)
        call setup_settings(oao%settings, context)
        oao%n_ao = n_ao
        oao%n_particle = n_particle
        oao%n_param = n_particle * n_ao * (n_ao - 1) / 2
        oao%s_inv_sqrt = identity_matrix(n_ao)
        allocate(oao%dm_oao(n_ao, n_ao, n_particle))
        do j = 1, n_particle
            oao%dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_occ(j))
        end do
        dm_hist = oao%dm_oao

        ! build the model from a history of independent points within the window of
        ! the trust radius, from one spanning very different step lengths for a trust
        ! radius whose window only holds the short step and without a state, and from
        ! a history of independent points without a non-linear response
        n_list = 2
        do i_case = 1, n_cases
            allocate(dm_list(n_ao, n_ao, n_particle, n_list), &
                     current(n_ao, n_ao, n_particle, n_pot), &
                     lists(n_ao, n_ao, n_particle, n_list, n_pot))
            select case (i_case)
            case (1)
                case_name = "a history of independent points"
            case (2)
                case_name = "a history spanning very different step lengths"
            case (3)
                case_name = &
                    "a history spanning very different step lengths without a state"
            case default
                case_name = "a history without a non-linear response"
            end select
            do j = 1, n_particle
                if (i_case == 2 .or. i_case == 3) then
                    call random_number(perturbation)
                    dm_list(:, :, j, 1) = oao%dm_oao(:, :, j) + 1e-6_rp * &
                                          (perturbation + transpose(perturbation))
                    dm_list(:, :, j, 2) = generate_random_density_matrix(n_ao, n_occ(j))
                else
                    do k = 1, n_list
                        call random_number(perturbation)
                        dm_list(:, :, j, k) = oao%dm_oao(:, :, j) + 1e-2_rp * &
                                              (perturbation + transpose(perturbation))
                    end do
                end if
            end do

            ! potentials which are functions of the density, except for a non-linear
            ! potential which does not change with the density
            do i = 1, n_pot
                current(:, :, :, i) = mock_potential(mock_potential_factors(i), &
                                                     oao%dm_oao)
                do k = 1, n_list
                    lists(:, :, :, k, i) = mock_potential(mock_potential_factors(i), &
                                                          dm_list(:, :, :, k))
                end do
            end do
            if (i_case == 4) &
                lists(:, :, :, :, n_pot) = spread(current(:, :, :, n_pot), 4, n_list)

            ! build the model for every QN method on a newly set up density QN object
            do i_method = 1, size(qn_methods)
                case_qn_method = trim(qn_methods(i_method))
                allocate(density_qn, source=density_qn_oao_type(oao))
                call setup_settings(density_qn%settings, context)
                density_qn%settings%qn_method = case_qn_method
                if (i_case /= 3) then
                    density_qn%state => state
                    state%trust_radius = case_trust_radii(i_case)
                end if
                call build_from_history()

                ! determine if the model is assembled from the history without
                ! evaluating the density matrix or changing the history
                if (error /= 0) then
                    write(stderr, *) test_name//" failed: Produced error for "// &
                        case_qn_method//" and "//case_name//"."
                    passed = .false.
                end if
                if (density_qn%model_stale) then
                    write(stderr, *) test_name//" failed: Hessian model still "// &
                        "marked stale for "//case_qn_method//" and "//case_name//"."
                    passed = .false.
                end if
                if (size(mock_requests) /= 0) then
                    write(stderr, *) test_name// &
                        " failed: Density matrix evaluating function called for "// &
                        case_qn_method//" and "//case_name//"."
                    passed = .false.
                end if
                do k = 1, n_list
                    if (norm2(density_qn%dm_list(:, :, :, k) - dm_list(:, :, :, k)) > &
                        tol .or. norm2(stacked_potentials(density_qn, k) - &
                                       lists(:, :, :, k, :)) > tol) then
                        write(stderr, *) test_name//" failed: History changed for "// &
                            case_qn_method//" and "//case_name//"."
                        passed = .false.
                    end if
                end do
                if (.not. (allocated(density_qn%expansion_dirs) .and. &
                           allocated(density_qn%projection_dirs) .and. &
                           allocated(density_qn%coupling_matrix))) then
                    write(stderr, *) test_name//" failed: Low-rank Hessian factors "// &
                        "not assembled for "//case_qn_method//" and "//case_name//"."
                    passed = .false.
                end if

                ! determine if the linear and non-linear systems of the QN method are
                ! constructed
                if (case_qn_method == "ms_sr1") then
                    built = allocated(density_qn%a_inv) .and. &
                            allocated(density_qn%a_inv_comb) .and. &
                            allocated(density_qn%linear_potential_dirs) .and. &
                            allocated(density_qn%nonlinear_potential_dirs)
                else
                    built = allocated(density_qn%dm_dirs) .and. &
                            allocated(density_qn%dm_dirs_nonlinear)
                    if (case_qn_method /= "ms_sp") built = &
                        built .and. allocated(density_qn%linear_potential_dirs) .and. &
                        allocated(density_qn%nonlinear_potential_dirs)
                    if (case_qn_method == "ms_sp" .or. case_qn_method == "ms_psb") &
                        built = built .and. allocated(density_qn%a_sym) .and. &
                                allocated(density_qn%a_sym_nonlinear)
                end if
                if (.not. built) then
                    write(stderr, *) test_name//" failed: Linear and non-linear "// &
                        "systems not constructed for "//case_qn_method//" and "// &
                        case_name//"."
                    passed = .false.
                    deallocate(density_qn)
                    cycle
                end if
                if (case_qn_method == "ms_sr1") then
                    n_linear = size(density_qn%linear_potential_dirs, 2)
                    n_nonlinear = size(density_qn%nonlinear_potential_dirs, 2)
                else
                    n_linear = size(density_qn%dm_dirs, 2)
                    n_nonlinear = size(density_qn%dm_dirs_nonlinear, 2)
                end if

                ! determine if the linear system keeps the entire history, for
                ! multisecant SR1 in both its Coulomb and its exact-exchange term,
                ! while the non-linear system, which multisecant SR1 factorizes for
                ! both channels together and the ARH family per channel, drops the
                ! entries beyond the window of the trust radius and is left out without
                ! a non-linear response
                n_systems = merge(1_ip, n_particle, case_qn_method == "ms_sr1")
                n_terms = merge(2_ip, 1_ip, case_qn_method == "ms_sr1")
                if (n_linear /= n_terms * n_particle * n_list) then
                    write(stderr, *) test_name//" failed: The linear system does "// &
                        "not keep the entire history for "//case_qn_method//" and "// &
                        case_name//"."
                    passed = .false.
                end if
                if ((i_case == 1 .or. i_case == 3) .and. &
                    n_nonlinear /= n_systems * n_list) then
                    write(stderr, *) test_name//" failed: The non-linear system "// &
                        "does not keep every independent entry for "//case_qn_method// &
                        " and "//case_name//"."
                    passed = .false.
                else if (i_case == 2 .and. n_nonlinear /= n_systems) then
                    write(stderr, *) test_name//" failed: The non-linear system "// &
                        "does not drop the history entry beyond the window of the "// &
                        "trust radius for "//case_qn_method//" and "//case_name//"."
                    passed = .false.
                else if (i_case == 4 .and. n_nonlinear /= 0) then
                    write(stderr, *) test_name//" failed: The non-linear system is "// &
                        "not left out for "//case_qn_method//" and "//case_name//"."
                    passed = .false.
                end if

                ! determine if every quantity follows the basis of its own system
                if (case_qn_method == "ms_sr1") then
                    if (size(density_qn%a_inv, 1) /= n_linear .or. &
                        size(density_qn%a_inv_comb, 1) /= n_nonlinear) then
                        write(stderr, *) test_name//" failed: Multisecant SR1 "// &
                            "pseudoinverses do not match the potential difference "// &
                            "directions of their system for "//case_name//"."
                        passed = .false.
                    end if
                else if (case_qn_method /= "ms_sp") then
                    if (size(density_qn%linear_potential_dirs, 2) /= n_linear .or. &
                        size(density_qn%nonlinear_potential_dirs, 2) /= n_nonlinear) &
                        then
                        write(stderr, *) test_name//" failed: Potential difference "// &
                            "directions are not rebased onto the basis of the "// &
                            "density matrix difference directions of their system "// &
                            "for "//case_qn_method//" and "//case_name//"."
                        passed = .false.
                    end if
                end if
                if (case_qn_method == "ms_sp" .or. case_qn_method == "ms_psb") then
                    if (size(density_qn%a_sym, 1) /= n_linear .or. &
                        size(density_qn%a_sym_nonlinear, 1) /= n_nonlinear) then
                        write(stderr, *) test_name//" failed: Symmetrized A "// &
                            "matrices do not match the density matrix difference "// &
                            "directions of their system for "//case_qn_method// &
                            " and "//case_name//"."
                        passed = .false.
                    end if
                end if

                ! deallocate density QN object
                deallocate(density_qn)
            end do
            deallocate(dm_list, current, lists)
        end do

        ! deallocate OAO object
        deallocate(oao)

    contains

        subroutine build_from_history()
            !
            ! this subroutine hands the history, the current density matrix and the
            ! current potentials to the density QN object, marks its model stale and
            ! builds the model for the shell
            !
            ! history, current density matrix and current potentials
            density_qn%dm_list = dm_list
            density_qn%dm_current = dm_hist
            if (n_particle == 1) then
                density_qn%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                density_qn%evaluate_dm_os => mock_evaluate_dm_os
            end if
            density_qn%v_coulomb = current(:, :, :, 1)
            density_qn%v_coulomb_list = lists(:, :, :, :, 1)
            density_qn%v_exchange = current(:, :, :, 2)
            density_qn%v_exchange_list = lists(:, :, :, :, 2)
            density_qn%v_nonlinear = current(:, :, :, n_pot)
            density_qn%v_nonlinear_list = lists(:, :, :, :, n_pot)
            density_qn%model_stale = .true.

            ! build the model
            mock_requests = [integer(ip) :: ]
            if (n_particle == 1) then
                call density_qn%build_hess_model_cs(error)
            else
                call density_qn%build_hess_model_os(error)
            end if

        end subroutine build_from_history

    end function check_build_hess_model

    logical(c_bool) function test_density_qn_factory_cs() bind(C)
        !
        ! this function tests the subroutine which returns the density QN orbital
        ! updating function and objective function for the closed-shell case, which
        ! hands the level of the solver settings' context and the closed-shell density
        ! matrix evaluating function to the common setup
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_state_type
        use otr_density_qn, only: density_qn_factory_cs, density_qn_type, &
                                  density_qn_settings_type, &
                                  density_qn_evaluate_dm_cs_type
        use otr_mo, only: mo_type
        use otr_common, only: extension_type, extension_solver_settings_type, &
                              update_orbs_callback, obj_func_callback, error_density_qn
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_common_unit_tests, only: setup_owning_solver_settings
        use otr_density_qn_test_reference, only: operator(==)

        type(test_context_type), target :: context
        type(extension_solver_settings_type), target :: solver_settings
        type(solver_state_type), target :: state
        class(extension_type), allocatable, target :: basis
        type(density_qn_settings_type) :: settings
        procedure(density_qn_evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), pointer :: update_orbs_density_qn_funptr
        procedure(obj_func_type), pointer :: obj_func_density_qn_funptr
        integer(ip) :: error

        ! assume tests pass
        test_density_qn_factory_cs = .true.

        ! set up solver settings wired with a closed-shell MO object
        allocate(mo_type :: basis)
        basis%n_param = 1
        select type (basis)
        type is (mo_type)
            basis%n_particle = 1
        end select
        call setup_owning_solver_settings(solver_settings, basis)
        call setup_settings(settings, context)
        evaluate_dm_funptr => mock_evaluate_dm_cs

        ! call routine and determine if the density QN object built on the MO object
        ! holds only the closed-shell density matrix evaluating function and if the
        ! role callback functions are returned, the remaining setup is covered by the
        ! test of the common setup
        call density_qn_factory_cs(evaluate_dm_funptr, update_orbs_density_qn_funptr, &
                                   obj_func_density_qn_funptr, solver_settings, state, &
                                   error, settings)
        if (error /= 0) then
            write(stderr, *) "test_density_qn_factory_cs failed: Produced error."
            test_density_qn_factory_cs = .false.
        end if
        if (.not. associated(update_orbs_density_qn_funptr, update_orbs_callback) .or. &
            .not. associated(obj_func_density_qn_funptr, obj_func_callback)) then
            write(stderr, *) "test_density_qn_factory_cs failed: Role callback "// &
                "functions not returned."
            test_density_qn_factory_cs = .false.
        end if
        if (size(solver_settings%extensions) /= 2) then
            write(stderr, *) "test_density_qn_factory_cs failed: Solver settings "// &
                "not wired with an object."
            test_density_qn_factory_cs = .false.
        else
            select type (density_qn => solver_settings%extensions(2)%object)
            class is (density_qn_type)
                select type (mo => solver_settings%extensions(1)%object)
                type is (mo_type)
                    if (.not. associated(density_qn%orbitals, mo)) then
                        write(stderr, *) "test_density_qn_factory_cs failed: "// &
                            "density QN object not built on the orbital basis "// &
                            "object of the solver settings."
                        test_density_qn_factory_cs = .false.
                    end if
                end select
                if (.not. associated(density_qn%state, state)) then
                    write(stderr, *) "test_density_qn_factory_cs failed: Record of "// &
                        "the state not passed on."
                    test_density_qn_factory_cs = .false.
                end if
                if (.not. (density_qn%settings == settings)) then
                    write(stderr, *) "test_density_qn_factory_cs failed: Settings "// &
                        "not passed on."
                    test_density_qn_factory_cs = .false.
                end if
                if (.not. associated(density_qn%evaluate_dm_cs, &
                                     mock_evaluate_dm_cs) .or. &
                    associated(density_qn%evaluate_dm_os)) then
                    write(stderr, *) "test_density_qn_factory_cs failed: "// &
                        "Closed-shell density matrix evaluating function not "// &
                        "handed over."
                    test_density_qn_factory_cs = .false.
                end if
            class default
                write(stderr, *) "test_density_qn_factory_cs failed: Solver "// &
                    "settings not wired with a density QN object."
                test_density_qn_factory_cs = .false.
            end select
        end if

        ! call routine for solver settings no factory wired and determine if the error
        ! of the common setup is passed on
        solver_settings%context => null()
        call density_qn_factory_cs(evaluate_dm_funptr, update_orbs_density_qn_funptr, &
                                   obj_func_density_qn_funptr, solver_settings, state, &
                                   error, settings)
        if (error /= error_density_qn + 1) then
            write(stderr, *) "test_density_qn_factory_cs failed: Error of the "// &
                "common setup not passed on."
            test_density_qn_factory_cs = .false.
        end if

    end function test_density_qn_factory_cs

    logical(c_bool) function test_density_qn_factory_os() bind(C)
        !
        ! this function tests the subroutine which returns the density QN orbital
        ! updating function and objective function for the open-shell case, which hands
        ! the level of the solver settings' context and the open-shell density matrix
        ! evaluating function to the common setup
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_state_type
        use otr_density_qn, only: density_qn_factory_os, density_qn_type, &
                                  density_qn_settings_type, &
                                  density_qn_evaluate_dm_os_type
        use otr_oao, only: oao_type
        use otr_common, only: extension_type, extension_solver_settings_type, &
                              update_orbs_callback, obj_func_callback, error_density_qn
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_common_unit_tests, only: setup_owning_solver_settings
        use otr_density_qn_test_reference, only: operator(==)

        type(test_context_type), target :: context
        type(extension_solver_settings_type), target :: solver_settings
        type(solver_state_type), target :: state
        class(extension_type), allocatable, target :: basis
        type(density_qn_settings_type) :: settings
        procedure(density_qn_evaluate_dm_os_type), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), pointer :: update_orbs_density_qn_funptr
        procedure(obj_func_type), pointer :: obj_func_density_qn_funptr
        integer(ip) :: error

        ! assume tests pass
        test_density_qn_factory_os = .true.

        ! set up solver settings wired with an open-shell OAO object
        allocate(oao_type :: basis)
        basis%n_param = 1
        select type (basis)
        type is (oao_type)
            basis%n_particle = 2
        end select
        call setup_owning_solver_settings(solver_settings, basis)
        call setup_settings(settings, context)
        evaluate_dm_funptr => mock_evaluate_dm_os

        ! call routine and determine if the density QN object built on the OAO object
        ! holds only the open-shell density matrix evaluating function and if the role
        ! callback functions are returned, the remaining setup is covered by the test
        ! of the common setup
        call density_qn_factory_os(evaluate_dm_funptr, update_orbs_density_qn_funptr, &
                                   obj_func_density_qn_funptr, solver_settings, state, &
                                   error, settings)
        if (error /= 0) then
            write(stderr, *) "test_density_qn_factory_os failed: Produced error."
            test_density_qn_factory_os = .false.
        end if
        if (.not. associated(update_orbs_density_qn_funptr, update_orbs_callback) .or. &
            .not. associated(obj_func_density_qn_funptr, obj_func_callback)) then
            write(stderr, *) "test_density_qn_factory_os failed: Role callback "// &
                "functions not returned."
            test_density_qn_factory_os = .false.
        end if
        if (size(solver_settings%extensions) /= 2) then
            write(stderr, *) "test_density_qn_factory_os failed: Solver settings "// &
                "not wired with an object."
            test_density_qn_factory_os = .false.
        else
            select type (density_qn => solver_settings%extensions(2)%object)
            class is (density_qn_type)
                select type (oao => solver_settings%extensions(1)%object)
                type is (oao_type)
                    if (.not. associated(density_qn%orbitals, oao)) then
                        write(stderr, *) "test_density_qn_factory_os failed: "// &
                            "density QN object not built on the orbital basis "// &
                            "object of the solver settings."
                        test_density_qn_factory_os = .false.
                    end if
                end select
                if (.not. associated(density_qn%state, state)) then
                    write(stderr, *) "test_density_qn_factory_os failed: Record of "// &
                        "the state not passed on."
                    test_density_qn_factory_os = .false.
                end if
                if (.not. (density_qn%settings == settings)) then
                    write(stderr, *) "test_density_qn_factory_os failed: Settings "// &
                        "not passed on."
                    test_density_qn_factory_os = .false.
                end if
                if (.not. associated(density_qn%evaluate_dm_os, &
                                     mock_evaluate_dm_os) .or. &
                    associated(density_qn%evaluate_dm_cs)) then
                    write(stderr, *) "test_density_qn_factory_os failed: "// &
                        "Open-shell density matrix evaluating function not handed over."
                    test_density_qn_factory_os = .false.
                end if
            class default
                write(stderr, *) "test_density_qn_factory_os failed: Solver "// &
                    "settings not wired with a density QN object."
                test_density_qn_factory_os = .false.
            end select
        end if

        ! call routine for solver settings no factory wired and determine if the error
        ! of the common setup is passed on
        solver_settings%context => null()
        call density_qn_factory_os(evaluate_dm_funptr, update_orbs_density_qn_funptr, &
                                   obj_func_density_qn_funptr, solver_settings, state, &
                                   error, settings)
        if (error /= error_density_qn + 1) then
            write(stderr, *) "test_density_qn_factory_os failed: Error of the "// &
                "common setup not passed on."
            test_density_qn_factory_os = .false.
        end if

    end function test_density_qn_factory_os

    logical(c_bool) function test_density_qn_factory_common() bind(C)
        !
        ! this function tests the subroutine which sets up a new density QN object on
        ! the orbital basis object of a level with the density matrix evaluating
        ! function of the shell of its orbitals, after the sanity check of its
        ! settings, wires the solver settings with it and the record of the state,
        ! adjusts them to the approximate Hessian and returns its orbital updating
        ! function and objective function
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_state_type
        use otr_density_qn, only: &
            density_qn_factory_common, density_qn_type, density_qn_mo_type, &
            density_qn_oao_type, density_qn_settings_type, &
            default_density_qn_settings, density_qn_n_micro, &
            density_qn_evaluate_dm_cs_type, density_qn_evaluate_dm_os_type
        use otr_mo, only: mo_type
        use otr_oao, only: oao_type
        use otr_common, only: extension_type, extension_level_type, &
                              extension_solver_settings_type, orbital_basis_type, &
                              update_orbs_callback, obj_func_callback, error_density_qn
        use otr_common_unit_tests, only: mock_extension_type, &
                                         setup_owning_solver_settings
        use otr_common_test_reference, only: same_logger_and_context
        use opentrustregion_unit_tests, only: setup_settings, setup_error_logging, &
                                              test_context_type
        use otr_density_qn_test_reference, only: operator(==)

        type(extension_solver_settings_type), allocatable, target :: solver_settings, &
                                                                     other_settings
        type(solver_state_type), target :: state, other_state
        class(extension_type), allocatable, target :: basis
        type(extension_level_type), pointer :: level
        type(density_qn_settings_type) :: settings, expected_settings
        type(test_context_type), target :: context
        procedure(density_qn_evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(density_qn_evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), pointer :: update_orbs_density_qn_funptr
        procedure(obj_func_type), pointer :: obj_func_density_qn_funptr
        integer(ip) :: i_case, error
        logical :: mo_basis, right_type, built, stored
        character(len=41), parameter :: case_names(3) = &
            [character(len=41) :: "for the closed-shell case in the MO basis", &
             "for the open-shell case in the OAO basis", "without settings"]
        character(len=74), parameter :: refusal_names(7) = [ &
            character(len=74) :: "an unknown QN method", "solver settings no "// &
            "factory wired", "solver settings another extension wired", "an "// &
            "open-shell density matrix evaluating function for closed-shell orbitals", &
            "a closed-shell density matrix evaluating function for open-shell "// &
            "orbitals", "no density matrix evaluating function for closed-shell "// &
            "orbitals", "no density matrix evaluating function for open-shell orbitals"]

        ! assume tests pass
        test_density_qn_factory_common = .true.

        ! create objects on a closed-shell MO object, on an open-shell OAO object and
        ! for no settings, each of which has to be the density QN object of the basis
        ! of the orbital basis object, built on it, with the settings and the density
        ! matrix evaluating function, and to wire the solver settings, which own it
        ! after the orbital basis object, with the record of the state and adjust them
        ! to the approximate Hessian of its QN method, which settings in uppercase
        ! select, while the orbital basis object keeps its own settings
        do i_case = 1, size(case_names)
            mo_basis = i_case /= 2
            call setup_basis(mo_basis, merge(2_ip, 1_ip, i_case == 2))
            call setup_settings(settings, context)
            settings%qn_method = "ARH"
            expected_settings = settings
            expected_settings%qn_method = "arh"
            if (i_case == 3) expected_settings = default_density_qn_settings
            evaluate_dm_cs_funptr => null()
            evaluate_dm_os_funptr => null()
            if (i_case == 2) then
                evaluate_dm_os_funptr => mock_evaluate_dm_os
            else
                evaluate_dm_cs_funptr => mock_evaluate_dm_cs
            end if
            level => solver_settings%extensions(1)%object%solver_level
            if (i_case == 3) then
                call density_qn_factory_common( &
                    level, evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                    update_orbs_density_qn_funptr, obj_func_density_qn_funptr, &
                    solver_settings, state, error)
            else
                call density_qn_factory_common( &
                    level, evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                    update_orbs_density_qn_funptr, obj_func_density_qn_funptr, &
                    solver_settings, state, error, settings)
            end if
            if (error /= 0) then
                write(stderr, *) "test_density_qn_factory_common failed: Produced "// &
                    "error "//trim(case_names(i_case))//"."
                test_density_qn_factory_common = .false.
                deallocate(solver_settings)
                cycle
            end if
            if (.not. associated(update_orbs_density_qn_funptr, &
                                 update_orbs_callback) .or. &
                .not. associated(obj_func_density_qn_funptr, obj_func_callback)) then
                write(stderr, *) "test_density_qn_factory_common failed: Role "// &
                    "callback functions not returned "//trim(case_names(i_case))//"."
                test_density_qn_factory_common = .false.
            end if
            if (size(solver_settings%extensions) /= 2) then
                write(stderr, *) "test_density_qn_factory_common failed: Object "// &
                    "not owned by the solver settings "//trim(case_names(i_case))//"."
                test_density_qn_factory_common = .false.
                deallocate(solver_settings)
                cycle
            end if
            select type (density_qn => solver_settings%extensions(2)%object)
            class is (density_qn_type)
                select type (density_qn)
                type is (density_qn_mo_type)
                    right_type = mo_basis
                type is (density_qn_oao_type)
                    right_type = .not. mo_basis
                class default
                    right_type = .false.
                end select
                built = .false.
                select type (orbitals => solver_settings%extensions(1)%object)
                class is (orbital_basis_type)
                    built = associated(density_qn%orbitals, orbitals)
                    if (orbitals%settings%verbose /= 7) then
                        write(stderr, *) "test_density_qn_factory_common failed: "// &
                            "Settings of the orbital basis object changed "// &
                            trim(case_names(i_case))//"."
                        test_density_qn_factory_common = .false.
                    end if
                end select
                if (.not. (right_type .and. built)) then
                    write(stderr, *) "test_density_qn_factory_common failed: "// &
                        "density QN object of the basis of the orbital basis "// &
                        "object not built on it "//trim(case_names(i_case))//"."
                    test_density_qn_factory_common = .false.
                end if
                if (density_qn%error_origin /= error_density_qn) then
                    write(stderr, *) "test_density_qn_factory_common failed: Error "// &
                        "origin of the object not set "//trim(case_names(i_case))//"."
                    test_density_qn_factory_common = .false.
                end if
                if (i_case == 2) then
                    stored = associated(density_qn%evaluate_dm_os, &
                                        mock_evaluate_dm_os) .and. &
                             .not. associated(density_qn%evaluate_dm_cs)
                else
                    stored = associated(density_qn%evaluate_dm_cs, &
                                        mock_evaluate_dm_cs) .and. &
                             .not. associated(density_qn%evaluate_dm_os)
                end if
                if (.not. (stored .and. density_qn%settings == expected_settings .and. &
                           same_logger_and_context(density_qn%settings, &
                                                   expected_settings))) then
                    write(stderr, *) "test_density_qn_factory_common failed: "// &
                        "Settings or density matrix evaluating function not stored "// &
                        trim(case_names(i_case))//"."
                    test_density_qn_factory_common = .false.
                end if
                if (.not. associated(density_qn%state, state) .or. .not. associated( &
                    solver_settings%context, density_qn%solver_level)) then
                    write(stderr, *) "test_density_qn_factory_common failed: "// &
                        "Solver settings not wired with the object reading the "// &
                        "record of the state "//trim(case_names(i_case))//"."
                    test_density_qn_factory_common = .false.
                end if
            class default
                write(stderr, *) "test_density_qn_factory_common failed: Solver "// &
                    "settings not wired with a density QN object "// &
                    trim(case_names(i_case))//"."
                test_density_qn_factory_common = .false.
            end select
            if (.not. solver_settings%refresh_hess .or. &
                solver_settings%n_micro /= density_qn_n_micro .or. &
                (solver_settings%hess_symm .neqv. &
                 expected_settings%qn_method /= "arh")) then
                write(stderr, *) "test_density_qn_factory_common failed: Solver "// &
                    "settings not adjusted to the approximate Hessian "// &
                    trim(case_names(i_case))//"."
                test_density_qn_factory_common = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create objects for an unknown QN method, for solver settings no factory
        ! wired, for solver settings another extension than an orbital basis wired and
        ! with the density matrix evaluating function of the other shell or without any
        ! for either shell of the orbitals, which have to return an error with an error
        ! message and leave the solver settings untouched
        do i_case = 1, size(refusal_names)
            call setup_basis(.true., merge(2_ip, 1_ip, i_case == 5 .or. i_case == 7))
            call setup_settings(settings, context)
            call setup_error_logging(settings, context)
            level => solver_settings%extensions(1)%object%solver_level
            evaluate_dm_cs_funptr => mock_evaluate_dm_cs
            evaluate_dm_os_funptr => null()
            select case (i_case)
            case (1)
                settings%qn_method = "unknown"
            case (2)
                level => null()
            case (3)
                allocate(mock_extension_type :: basis)
                solver_settings = extension_solver_settings_type()
                call setup_owning_solver_settings(solver_settings, basis)
                level => solver_settings%extensions(1)%object%solver_level
            case (4)
                evaluate_dm_cs_funptr => null()
                evaluate_dm_os_funptr => mock_evaluate_dm_os
            case (6, 7)
                evaluate_dm_cs_funptr => null()
            end select
            call density_qn_factory_common( &
                level, evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                update_orbs_density_qn_funptr, obj_func_density_qn_funptr, &
                solver_settings, state, error, settings)
            if (error /= error_density_qn + 1) then
                write(stderr, *) "test_density_qn_factory_common failed: No error "// &
                    "returned for "//trim(refusal_names(i_case))//"."
                test_density_qn_factory_common = .false.
            end if
            if (context%log_message == "") then
                write(stderr, *) "test_density_qn_factory_common failed: No error "// &
                    "message printed for "//trim(refusal_names(i_case))//"."
                test_density_qn_factory_common = .false.
            end if
            if (size(solver_settings%extensions) /= 1 .or. .not. associated( &
                solver_settings%context, &
                solver_settings%extensions(1)%object%solver_level)) then
                write(stderr, *) "test_density_qn_factory_common failed: Solver "// &
                    "settings changed for "//trim(refusal_names(i_case))//"."
                test_density_qn_factory_common = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create objects for two separate solver settings, which have to be distinct
        ! objects reading their own record of the state, so that calculations with
        ! separate settings can be used alternately
        call setup_basis(.true., 1_ip)
        call move_alloc(solver_settings, other_settings)
        call setup_basis(.true., 1_ip)
        call setup_settings(settings, context)
        evaluate_dm_cs_funptr => mock_evaluate_dm_cs
        evaluate_dm_os_funptr => null()
        call density_qn_factory_common( &
            solver_settings%extensions(1)%object%solver_level, evaluate_dm_cs_funptr, &
            evaluate_dm_os_funptr, update_orbs_density_qn_funptr, &
            obj_func_density_qn_funptr, solver_settings, state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_density_qn_factory_common failed: Produced "// &
                "error for the first of two separate solver settings."
            test_density_qn_factory_common = .false.
        end if
        call density_qn_factory_common( &
            other_settings%extensions(1)%object%solver_level, evaluate_dm_cs_funptr, &
            evaluate_dm_os_funptr, update_orbs_density_qn_funptr, &
            obj_func_density_qn_funptr, other_settings, other_state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_density_qn_factory_common failed: Produced "// &
                "error for the second of two separate solver settings."
            test_density_qn_factory_common = .false.
        end if
        if (size(solver_settings%extensions) /= 2 .or. &
            size(other_settings%extensions) /= 2) then
            write(stderr, *) "test_density_qn_factory_common failed: Separate "// &
                "solver settings not wired."
            test_density_qn_factory_common = .false.
        else if (.not. associated( &
            solver_settings%extensions(2)%object%state, state) .or. .not. associated( &
                other_settings%extensions(2)%object%state, other_state)) then
            write(stderr, *) "test_density_qn_factory_common failed: Separate "// &
                "solver settings not wired with distinct objects reading their own "// &
                "state."
            test_density_qn_factory_common = .false.
        end if
        deallocate(solver_settings, other_settings)

    contains

        subroutine setup_basis(mo, n_particle)
            !
            ! this subroutine sets up solver settings wired with an MO or OAO object
            ! for the given number of particles, whose settings carry a verbosity of
            ! their own
            !
            logical, intent(in) :: mo
            integer(ip), intent(in) :: n_particle

            if (mo) then
                allocate(mo_type :: basis)
            else
                allocate(oao_type :: basis)
            end if
            basis%n_param = 1
            select type (basis)
            class is (orbital_basis_type)
                basis%n_particle = n_particle
                call setup_settings(basis%settings, context)
                basis%settings%verbose = 7
            end select
            allocate(solver_settings)
            call setup_owning_solver_settings(solver_settings, basis)

        end subroutine setup_basis

    end function test_density_qn_factory_common

    logical(c_bool) function test_density_qn_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! density QN input parameters
        !
        use otr_density_qn, only: density_qn_settings_type, density_qn_sanity_check, &
                                  qn_methods
        use opentrustregion_unit_tests, only: test_context_type, setup_error_logging

        type(test_context_type), target :: context
        type(density_qn_settings_type) :: settings
        integer(ip) :: i, error

        ! assume tests pass
        test_density_qn_sanity_check = .true.

        ! setup settings object which passes only error messages on
        call setup_error_logging(settings, context)

        ! check if all available QN methods are accepted without a message
        do i = 1, size(qn_methods)
            settings%qn_method = qn_methods(i)
            call density_qn_sanity_check(settings, error)
            if (error /= 0 .or. context%log_message /= "") then
                write(stderr, *) "test_density_qn_sanity_check failed: Error "// &
                    "thrown for "//trim(qn_methods(i))//" QN method."
                test_density_qn_sanity_check = .false.
            end if
        end do

        ! check if QN method is converted to lowercase
        settings%qn_method = "MS_PSB"
        call density_qn_sanity_check(settings, error)
        if (settings%qn_method /= "ms_psb") then
            write(stderr, *) "test_density_qn_sanity_check failed: QN method not "// &
                "converted to lowercase."
            test_density_qn_sanity_check = .false.
        end if

        ! check if unknown QN method is rejected with an error message
        call setup_error_logging(settings, context)
        settings%qn_method = "unknown"
        call density_qn_sanity_check(settings, error)
        if (error /= 1) then
            write(stderr, *) "test_density_qn_sanity_check failed: Error not "// &
                "thrown for unknown QN method."
            test_density_qn_sanity_check = .false.
        end if
        if (context%log_message == "") then
            write(stderr, *) "test_density_qn_sanity_check failed: No error "// &
                "message printed for unknown QN method."
            test_density_qn_sanity_check = .false.
        end if

    end function test_density_qn_sanity_check

    logical(c_bool) function test_density_qn_set_solver_settings() bind(C)
        !
        ! this function tests the subroutine which requests the Hessian refresh, raises
        ! the micro iteration limit and reports the symmetry of the approximate Hessian
        ! in the solver settings
        !
        use otr_density_qn, only: density_qn_set_solver_settings, &
                                  density_qn_settings_type, density_qn_n_micro
        use opentrustregion, only: solver_settings_type
        use test_reference, only: ref_solver_settings, operator(/=)

        type(solver_settings_type) :: solver_settings
        type(density_qn_settings_type) :: settings
        integer(ip) :: i_method
        character(len=6), parameter :: qn_methods(2) = &
            [character(len=6) :: "arh", "ms_sr1"]

        ! assume tests pass
        test_density_qn_set_solver_settings = .true.

        ! adjust settings with the reference values, without the Hessian refresh and
        ! with the opposite symmetry, so that every adjustment shows, for density QN
        ! settings of the standard ARH method, whose approximate Hessian is not
        ! symmetric, and of one whose is, which apart from the Hessian refresh, the
        ! micro iteration limit and the symmetry of the approximate Hessian have to be
        ! kept
        do i_method = 1, size(qn_methods)
            settings%qn_method = qn_methods(i_method)
            solver_settings = ref_solver_settings
            solver_settings%refresh_hess = .false.
            solver_settings%hess_symm = i_method == 1
            call density_qn_set_solver_settings(solver_settings, settings)
            if (.not. solver_settings%refresh_hess) then
                write(stderr, *) "test_density_qn_set_solver_settings failed: "// &
                    "Hessian refresh not requested for QN method "// &
                    trim(qn_methods(i_method))//"."
                test_density_qn_set_solver_settings = .false.
            end if
            if (solver_settings%n_micro /= density_qn_n_micro) then
                write(stderr, *) "test_density_qn_set_solver_settings failed: "// &
                    "Micro iteration limit not raised for QN method "// &
                    trim(qn_methods(i_method))//"."
                test_density_qn_set_solver_settings = .false.
            end if
            if (solver_settings%hess_symm .neqv. (i_method /= 1)) then
                write(stderr, *) "test_density_qn_set_solver_settings failed: "// &
                    "Symmetry of the approximate Hessian reported wrongly for QN "// &
                    "method "//trim(qn_methods(i_method))//"."
                test_density_qn_set_solver_settings = .false.
            end if
            solver_settings%refresh_hess = ref_solver_settings%refresh_hess
            solver_settings%n_micro = ref_solver_settings%n_micro
            solver_settings%hess_symm = ref_solver_settings%hess_symm
            if (solver_settings /= ref_solver_settings) then
                write(stderr, *) "test_density_qn_set_solver_settings failed: "// &
                    "Settings not kept for QN method "//trim(qn_methods(i_method))//"."
                test_density_qn_set_solver_settings = .false.
            end if
        end do

    end function test_density_qn_set_solver_settings

    logical(c_bool) function test_obj_func_density_qn() bind(C)
        !
        ! this function tests the function which defines the energy evaluation, which
        ! also adds the evaluated point to the history, for the closed- and the
        ! open-shell case
        !
        use otr_common_test_reference, only: n_particle_ref => n_particle

        test_obj_func_density_qn = &
            check_obj_func_density_qn(1_ip, "test_obj_func_density_qn")
        if (.not. check_obj_func_density_qn(n_particle_ref, &
                                            "test_obj_func_density_qn")) &
            test_obj_func_density_qn = .false.

    contains

        function check_obj_func_density_qn(n_particle, test_name) result(passed)
            !
            ! this function checks the energy evaluation for the closed- or open-shell
            ! case for orbitals parameterized in the OAO basis, which also adds the
            ! evaluated point together with the potentials the shell keeps to the
            ! history
            !
            use test_reference, only: arm_host_context, host_context_reached
            use otr_density_qn, only: density_qn_type
            use otr_oao, only: oao_type
            use otr_common_test_reference, only: n_ao
            use otr_common_unit_tests, only: level_context, mock_requests
            use opentrustregion_unit_tests, only: test_context_type

            type(test_context_type), target :: context
            type(oao_type), allocatable, target :: oao
            class(density_qn_type), allocatable, target :: density_qn
            integer(ip), intent(in) :: n_particle
            character(len=*), intent(in) :: test_name
            logical :: passed

            real(rp), parameter :: basis_scale = 0.5_rp

            real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
            real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                        kappa(n_particle * n_ao * (n_ao - 1) / 2), energy
            real(rp), allocatable :: potentials(:, :, :, :)
            integer(ip) :: k, error

            ! assume test passes
            passed = .true.

            ! set up the OAO and density QN objects with an empty history, which is
            ! removed for the first call
            call setup_density_qn_and_scaled_oao_objects(density_qn, oao, basis_scale, &
                                                         dm_ao, dm_oao, context)
            call setup_empty_history(density_qn, n_particle)
            call arm_host_context(density_qn%settings, context)
            deallocate(density_qn%dm_list)

            ! call routine with an orbital rotation before the history exists and
            ! determine if the potentials the shell keeps are requested and nothing is
            ! added to the history
            call random_number(kappa)
            kappa = 0.1_rp * kappa
            mock_requests = [integer(ip) :: ]
            energy = density_qn%obj_func(kappa, error, level_context(density_qn, 1_ip))
            if (error /= 0) then
                write(stderr, *) test_name//" failed: Produced error without history."
                passed = .false.
            end if
            if (size(mock_requests) /= 1) then
                write(stderr, *) test_name//" failed: Density matrix evaluating "// &
                    "function not called exactly once."
                passed = .false.
            end if
            if (any(mock_requests /= 14_ip)) then
                write(stderr, *) test_name//" failed: Incorrect outputs requested "// &
                    "from density matrix evaluating function."
                passed = .false.
            end if
            if (allocated(density_qn%dm_list)) then
                write(stderr, *) test_name// &
                    " failed: History created by energy evaluation."
                passed = .false.
            end if
            if (density_qn%model_stale) then
                write(stderr, *) test_name// &
                    " failed: Hessian model marked stale without history."
                passed = .false.
            end if

            ! call routine again with an empty history and determine if the rotated
            ! density matrix is added together with the potentials evaluated at it in
            ! the history basis, while the current density matrix is left untouched
            allocate(density_qn%dm_list(n_ao, n_ao, n_particle, 0))
            mock_requests = [integer(ip) :: ]
            energy = density_qn%obj_func(kappa, error, level_context(density_qn, 1_ip))
            if (error /= 0) then
                write(stderr, *) test_name//" failed: Produced error with history."
                passed = .false.
            else if (size(density_qn%dm_list, 4) /= 1) then
                write(stderr, *) test_name// &
                    " failed: Evaluated point not added to history."
                passed = .false.
            else
                if (abs(energy - basis_scale**2 * &
                        sum(density_qn%dm_list(:, :, :, 1))) > tol) then
                    write(stderr, *) test_name//" failed: Incorrect energy."
                    passed = .false.
                end if
                if (norm2(density_qn%dm_list(:, :, :, 1) - dm_oao) < tol) then
                    write(stderr, *) test_name// &
                        " failed: Incorrect density matrix history."
                    passed = .false.
                end if
                potentials = stacked_potentials(density_qn, 1_ip)
                do k = 1, size(mock_potential_factors)
                    if (norm2(potentials(:, :, :, k) - basis_scale**4 * &
                              mock_potential(mock_potential_factors(k), &
                                             density_qn%dm_list(:, :, :, 1))) > tol) &
                        then
                        write(stderr, *) test_name// &
                            " failed: Incorrect potential history."
                        passed = .false.
                    end if
                end do
            end if
            if (norm2(oao%dm_oao - dm_oao) > tol .or. &
                norm2(density_qn%orbitals%dm_ao - basis_scale**2 * dm_oao) > tol) then
                write(stderr, *) test_name// &
                    " failed: Current density matrix changed by energy evaluation."
                passed = .false.
            end if
            if (.not. density_qn%model_stale) then
                write(stderr, *) test_name//" failed: Hessian model not marked "// &
                    "stale after adding the evaluated point to the history."
                passed = .false.
            end if

            ! call routine at the same point and determine if it is not added twice and
            ! leaves the Hessian model untouched
            density_qn%model_stale = .false.
            energy = density_qn%obj_func(kappa, error, level_context(density_qn, 1_ip))
            if (error /= 0) then
                write(stderr, *) test_name// &
                    " failed: Produced error for a point the history already holds."
                passed = .false.
            end if
            if (size(density_qn%dm_list, 4) /= 1) then
                write(stderr, *) test_name// &
                    " failed: History extended for a point it already holds."
                passed = .false.
            end if
            if (density_qn%model_stale) then
                write(stderr, *) test_name//" failed: Hessian model marked stale "// &
                    "for a point the history already holds."
                passed = .false.
            end if

            ! call routine at another point whose evaluation fails and determine if the
            ! error is passed on without adding the point to the history
            context%mock_error = 1
            energy = density_qn%obj_func(2.0_rp * kappa, error, &
                                         level_context(density_qn, 1_ip))
            if (error /= 1) then
                write(stderr, *) test_name//" failed: Error of the density matrix "// &
                    "evaluating function not passed on."
                passed = .false.
            end if
            if (size(density_qn%dm_list, 4) /= 1 .or. density_qn%model_stale) then
                write(stderr, *) test_name// &
                    " failed: Point of a failed evaluation added to the history."
                passed = .false.
            end if
            context%mock_error = 0

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("obj_func_density_qn", context)) &
                passed = .false.

            ! deallocate density QN and OAO objects
            deallocate(density_qn, oao)

        end function check_obj_func_density_qn

    end function test_obj_func_density_qn

    logical(c_bool) function test_update_orbs_density_qn() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation for the closed- and the open-shell case
        !
        use otr_common_test_reference, only: n_particle_ref => n_particle

        test_update_orbs_density_qn = &
            check_update_orbs_density_qn(1_ip, "test_update_orbs_density_qn")
        if (.not. check_update_orbs_density_qn(n_particle_ref, &
                                               "test_update_orbs_density_qn")) &
            test_update_orbs_density_qn = .false.

    contains

        function check_update_orbs_density_qn(n_particle, test_name) result(passed)
            !
            ! this function checks the energy, gradient and Hessian diagonal evaluation
            ! for the closed- or open-shell case for orbitals parameterized in the OAO
            ! basis, which passes the Fock matrix in the history basis on to the
            ! orbital basis, adds the point the orbitals are rotated away from to the
            ! history and relates the step to its size in the trust-region metric the
            ! state holds
            !
            use test_reference, only: arm_host_context, host_context_reached
            use opentrustregion, only: hess_x_type, solver_state_type
            use otr_common, only: hess_x_callback
            use otr_density_qn, only: density_qn_type
            use otr_oao, only: oao_type
            use otr_common_test_reference, only: n_ao
            use otr_common_unit_tests, only: level_context, mock_requests, &
                                             mock_fock_factor
            use opentrustregion_unit_tests, only: test_context_type, identity_matrix

            type(test_context_type), target :: context
            type(oao_type), allocatable, target :: oao
            class(density_qn_type), allocatable, target :: density_qn
            type(solver_state_type), target :: state
            integer(ip), intent(in) :: n_particle
            character(len=*), intent(in) :: test_name
            logical :: passed

            real(rp), parameter :: basis_scale = 0.5_rp, metric_ratio = 4.0_rp

            real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
            real(rp) :: dm_oao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle), &
                        dm_saved(n_ao, n_ao, n_particle, 2), complement(n_ao, n_ao), &
                        func
            real(rp), allocatable :: kappa(:), grad(:), h_diag(:), &
                                     potentials(:, :, :, :), &
                                     potentials_saved(:, :, :, :, :)
            integer(ip) :: n_param, n_list, k, error
            procedure(hess_x_type), pointer :: hess_x_funptr

            ! assume test passes
            passed = .true.

            ! set up the OAO and density QN objects
            call setup_density_qn_and_scaled_oao_objects(density_qn, oao, basis_scale, &
                                                         dm_ao, dm_oao, context)
            density_qn%settings%qn_method = "ms_psb"
            call arm_host_context(density_qn%settings, context)
            if (n_particle == 1) then
                density_qn%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                density_qn%evaluate_dm_os => mock_evaluate_dm_os
            end if
            n_param = density_qn%orbitals%n_param
            allocate(kappa(n_param), grad(n_param), h_diag(n_param), potentials_saved( &
                n_ao, n_ao, n_particle, size(mock_potential_factors), 2))

            ! the checks after a failed call would read quantities it did not set
            checks: block
                ! call routine without an orbital rotation for an uninitialized object
                kappa = 0.0_rp
                mock_requests = [integer(ip) :: ]
                call density_qn%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, &
                                            error, level_context(density_qn, 1_ip))
                if (error /= 0) then
                    write(stderr, *) test_name//" failed: Produced error."
                    passed = .false.
                    exit checks
                end if

                ! determine if the energy and the potentials of the density matrix
                ! evaluating function are stored in the history basis and if the Fock
                ! matrix in the history basis is passed on to the orbital basis, which
                ! keeps its occupied-occupied and virtual-virtual parts
                if (abs(func - sum(dm_ao)) > tol) then
                    write(stderr, *) test_name//" failed: Incorrect energy."
                    passed = .false.
                end if
                potentials = stacked_potentials(density_qn, 0_ip)
                do k = 1, size(mock_potential_factors)
                    if (norm2(potentials(:, :, :, k) - basis_scale**2 * &
                              mock_potential(mock_potential_factors(k), dm_ao)) > tol) &
                        then
                        write(stderr, *) test_name//" failed: Incorrect potential."
                        passed = .false.
                    end if
                end do
                fock = basis_scale**2 * mock_potential(mock_fock_factor(1), dm_ao)
                do k = 1, n_particle
                    complement = identity_matrix(n_ao) - oao%dm_oao(:, :, k)
                    if (norm2(oao%fock_oo(:, :, k) - &
                              matmul(oao%dm_oao(:, :, k), &
                                     matmul(fock(:, :, k), oao%dm_oao(:, :, k)))) + &
                        norm2(oao%fock_vv(:, :, k) - matmul( &
                            complement, matmul(fock(:, :, k), complement))) > tol) then
                        write(stderr, *) test_name//" failed: Fock matrix in the "// &
                            "history basis not passed on to the orbital basis."
                        passed = .false.
                    end if
                end do

                ! determine if the gradient, Hessian diagonal and Hessian linear
                ! transformation are returned
                if (norm2(grad - density_qn%orbitals%grad) > tol) then
                    write(stderr, *) test_name//" failed: Gradient not returned."
                    passed = .false.
                end if
                if (norm2(h_diag - density_qn%orbitals%h_diag) > tol) then
                    write(stderr, *) test_name// &
                        " failed: Hessian diagonal not returned."
                    passed = .false.
                end if
                if (.not. associated(hess_x_funptr, hess_x_callback)) then
                    write(stderr, *) test_name// &
                        " failed: Returned Hessian linear transformation is wrong."
                    passed = .false.
                end if

                ! determine if the history is initialized empty and the approximate
                ! Hessian model is assembled from it
                if (size(density_qn%dm_list, 4) /= 0) then
                    write(stderr, *) test_name// &
                        " failed: Density matrix history not initialized empty."
                    passed = .false.
                end if
                if (.not. allocated(density_qn%dm_dirs)) then
                    write(stderr, *) test_name//" failed: Approximate Hessian "// &
                        "model not assembled for the initial point."
                    passed = .false.
                end if

                ! determine if a call without an orbital rotation, which takes no step,
                ! leaves the ratio of the step to its size in the trust-region metric
                ! unchanged
                if (abs(density_qn%step_norm_ratio - 1.0_rp) > tol) then
                    write(stderr, *) test_name//" failed: Ratio of the step to its "// &
                        "size in the trust-region metric changed without a step."
                    passed = .false.
                end if

                ! call routine again without an orbital rotation and determine if the
                ! quantities of the already initialized object are reused, even for a
                ! vanishing energy, without rebuilding the static part of the Hessian
                oao%energy = 0.0_rp
                oao%hess_eigen_stale = .false.
                call density_qn%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, &
                                            error, level_context(density_qn, 1_ip))
                if (error /= 0) then
                    write(stderr, *) test_name// &
                        " failed: Produced error for an initialized object."
                    passed = .false.
                end if
                if (size(mock_requests) /= 1) then
                    write(stderr, *) test_name// &
                        " failed: Quantities recomputed without an orbital rotation."
                    passed = .false.
                end if
                if (size(density_qn%dm_list, 4) /= 0) then
                    write(stderr, *) test_name// &
                        " failed: History extended without an orbital rotation."
                    passed = .false.
                end if
                if (oao%hess_eigen_stale) then
                    write(stderr, *) test_name//" failed: Static part of the "// &
                        "Hessian rebuilt without an orbital rotation."
                    passed = .false.
                end if

                ! rotate twice in a row, the second time for a model marked stale, and
                ! save the quantities at the points the orbitals are rotated away from,
                ! which the history has to retain with the latest first, by steps whose
                ! size in the trust-region metric the state holds
                kappa = 0.1_rp
                density_qn%state => state
                state%step_size = norm2(kappa) / metric_ratio
                do k = 2, 1, -1
                    dm_saved(:, :, :, k) = oao%dm_oao
                    potentials_saved(:, :, :, :, k) = &
                        stacked_potentials(density_qn, 0_ip)
                    density_qn%model_stale = k == 1
                    call density_qn%update_orbs(kappa, func, grad, h_diag, &
                                                hess_x_funptr, error, &
                                                level_context(density_qn, 1_ip))
                    if (error /= 0) then
                        write(stderr, *) test_name// &
                            " failed: Produced error for an orbital rotation."
                        passed = .false.
                        exit checks
                    end if
                end do
                if (size(density_qn%dm_list, 4) /= 2) then
                    write(stderr, *) test_name// &
                        " failed: History not extended to two entries."
                    passed = .false.
                    exit checks
                end if
                do k = 1, 2
                    if (norm2(density_qn%dm_list(:, :, :, k) - dm_saved(:, :, :, k)) > &
                        tol) then
                        write(stderr, *) test_name// &
                            " failed: Incorrect density matrix history."
                        passed = .false.
                    end if
                    if (norm2(stacked_potentials(density_qn, k) - &
                              potentials_saved(:, :, :, :, k)) > tol) then
                        write(stderr, *) test_name// &
                            " failed: Incorrect potential history."
                        passed = .false.
                    end if
                end do

                ! determine if the ratio of the step to its size in the trust-region
                ! metric is recorded
                if (abs(density_qn%step_norm_ratio - metric_ratio) > tol) then
                    write(stderr, *) test_name//" failed: Ratio of the step to its "// &
                        "size in the trust-region metric not recorded."
                    passed = .false.
                end if

                ! determine if the approximate Hessian model is assembled from the
                ! history extended by the point the orbitals were rotated away from
                if (density_qn%model_stale .or. .not. allocated(density_qn%dm_dirs)) &
                    then
                    write(stderr, *) test_name// &
                        " failed: Approximate Hessian model not assembled."
                    passed = .false.
                else if (size(density_qn%dm_dirs, 2) /= 2 * n_particle) then
                    write(stderr, *) test_name//" failed: Approximate Hessian "// &
                        "model not assembled from the extended history."
                    passed = .false.
                end if

                ! determine if a rotation starting from a density matrix which the
                ! history already holds, as after the energy evaluation of an accepted
                ! trial point, does not add that density matrix a second time
                density_qn%dm_list(:, :, :, 1) = oao%dm_oao
                call density_qn%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, &
                                            error, level_context(density_qn, 1_ip))
                if (error /= 0) then
                    write(stderr, *) test_name//" failed: Produced error for a "// &
                        "density matrix already held in the history."
                    passed = .false.
                else if (size(density_qn%dm_list, 4) /= 2) then
                    write(stderr, *) test_name//" failed: History extended for a "// &
                        "density matrix it already holds."
                    passed = .false.
                end if

                ! call routine with an orbital rotation whose evaluation fails and
                ! determine if the error is passed on and the rotated density is marked
                ! as not evaluated, so that a following call without an orbital
                ! rotation evaluates it without adding it to the history
                context%mock_error = 1
                call density_qn%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, &
                                            error, level_context(density_qn, 1_ip))
                if (error /= 1) then
                    write(stderr, *) test_name//" failed: Error of the density "// &
                        "matrix evaluating function not passed on."
                    passed = .false.
                end if
                if (.not. density_qn%evaluation_stale) then
                    write(stderr, *) test_name//" failed: Evaluation not marked "// &
                        "stale after a failed evaluation."
                    passed = .false.
                end if
                context%mock_error = 0
                n_list = size(density_qn%dm_list, 4, kind=ip)
                kappa = 0.0_rp
                call density_qn%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, &
                                            error, level_context(density_qn, 1_ip))
                if (error /= 0) then
                    write(stderr, *) test_name// &
                        " failed: Produced error after a failed evaluation."
                    passed = .false.
                end if
                if (density_qn%evaluation_stale) then
                    write(stderr, *) test_name//" failed: Evaluation still marked "// &
                        "stale after a failed evaluation was recovered."
                    passed = .false.
                end if
                if (size(density_qn%dm_list, 4) /= n_list) then
                    write(stderr, *) test_name//" failed: History extended by a "// &
                        "density matrix which has not been evaluated."
                    passed = .false.
                end if

                ! determine if every recompute asked for all potentials the shell needs
                if (any(mock_requests /= 15_ip)) then
                    write(stderr, *) test_name//" failed: Incorrect outputs "// &
                        "requested from density matrix evaluating function."
                    passed = .false.
                end if
            end block checks

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("update_orbs_density_qn", context)) &
                passed = .false.

            ! deallocate density QN and OAO objects
            deallocate(density_qn, oao)

        end function check_update_orbs_density_qn

    end function test_update_orbs_density_qn

    logical(c_bool) function test_build_hess_model_cs() bind(C)
        !
        ! this function tests the subroutine which assembles the closed-shell
        ! approximate Hessian model from the history relative to the current point
        !
        test_build_hess_model_cs = &
            check_build_hess_model(1_ip, "test_build_hess_model_cs")

    end function test_build_hess_model_cs

    logical(c_bool) function test_build_hess_model_os() bind(C)
        !
        ! this function tests the subroutine which assembles the open-shell approximate
        ! Hessian model from the history relative to the current point
        !
        use otr_common_test_reference, only: n_particle

        test_build_hess_model_os = &
            check_build_hess_model(n_particle, "test_build_hess_model_os")

    end function test_build_hess_model_os

    logical(c_bool) function test_rebuild_stale_hess_model() bind(C)
        !
        ! this function tests the subroutine which rebuilds the approximate Hessian
        ! model if the objective function has added points to the history since it was
        ! last assembled, which calls the model builder of the shell of the orbitals and
        ! passes its error on
        !
        use otr_mo, only: mo_type

        type(recording_density_qn_type) :: density_qn
        type(mo_type), target :: mo
        integer(ip) :: i_case, error
        character(len=2), parameter :: expected_shells(3) = ["  ", "cs", "os"]
        character(len=36), parameter :: case_names(3) = &
            [character(len=36) :: "for a model which is not stale", &
             "for a stale closed-shell model", "for a stale open-shell model"]

        ! assume tests pass
        test_rebuild_stale_hess_model = .true.

        ! call routine for a model which is not stale and for stale models of orbitals
        ! of either shell, whose builders fail, and determine if only the builder of the
        ! shell is called and its error passed on
        density_qn%orbitals => mo
        do i_case = 1, size(case_names)
            density_qn%built_shell = ""
            density_qn%error = 1
            density_qn%model_stale = i_case /= 1
            mo%n_particle = merge(1_ip, 2_ip, i_case == 2)
            call density_qn%rebuild_stale_hess_model(error)
            if (density_qn%built_shell /= expected_shells(i_case)) then
                write(stderr, *) "test_rebuild_stale_hess_model failed: Wrong "// &
                    "model builder called "//trim(case_names(i_case))//"."
                test_rebuild_stale_hess_model = .false.
            end if
            if (error /= merge(0_ip, 1_ip, i_case == 1)) then
                write(stderr, *) "test_rebuild_stale_hess_model failed: Error of "// &
                    "the model builder not passed on "//trim(case_names(i_case))//"."
                test_rebuild_stale_hess_model = .false.
            end if
        end do

    end function test_rebuild_stale_hess_model

    logical(c_bool) function test_hess_x_density_qn() bind(C)
        !
        ! this function tests the subroutine which defines the Hessian linear
        ! transformation on the basis of augmented Roothaan-Hall and related methods,
        ! which adds the low-rank part to the static part of the orbital basis
        !
        use otr_common_unit_tests, only: level_context
        use otr_density_qn, only: density_qn_type
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, &
                                     setup_random_mo_channels, ref_hess_x_static_mo
        use opentrustregion_unit_tests, only: test_context_type

        type(test_context_type), target :: context
        integer(ip), parameter :: n_dirs = 2

        type(mo_type), allocatable, target :: mo
        class(density_qn_type), allocatable, target :: density_qn
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        real(rp), allocatable :: x(:), hess_x(:), expected_hess_x(:)
        integer(ip) :: n_param, i_case, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_hess_x_density_qn = .true.

        ! set up the density QN object in the MO basis with random Fock matrix blocks
        ! of the differently occupied spin channels
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_density_qn_and_mo_objects(density_qn, mo, mo_coeff, ao_overlap, &
                                             n_occ, context)
        call setup_random_mo_channels(mo, n_occ, n_mo)
        n_param = density_qn%orbitals%n_param
        allocate(x(n_param), hess_x(n_param))
        call random_number(x)

        ! call routine for a low-rank part with random, and hence non-symmetric,
        ! coupling matrix, for a stale model on an empty history, which has to be
        ! rebuilt before it is applied and thereby loses the low-rank part, and for the
        ! model without low-rank part which this leaves behind
        do i_case = 1, 3
            expected_hess_x = ref_hess_x_static_mo(x, mo%mo_channels)
            if (i_case == 1) then
                case_name = "with low-rank part"
                allocate(density_qn%expansion_dirs(n_param, n_dirs), &
                         density_qn%projection_dirs(n_param, n_dirs), &
                         density_qn%coupling_matrix(n_dirs, n_dirs))
                call random_number(density_qn%expansion_dirs)
                call random_number(density_qn%projection_dirs)
                call random_number(density_qn%coupling_matrix)
                expected_hess_x = expected_hess_x + matmul( &
                    density_qn%expansion_dirs, &
                    matmul(density_qn%coupling_matrix, &
                           matmul(transpose(density_qn%projection_dirs), x)))
            else if (i_case == 2) then
                case_name = "for stale model"
                call setup_empty_history(density_qn, n_particle)
                density_qn%model_stale = .true.
            else
                case_name = "without low-rank part"
            end if
            call density_qn%hess_x(x, hess_x, error, level_context(density_qn, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_hess_x_density_qn failed: Produced error "// &
                    case_name//"."
                test_hess_x_density_qn = .false.
            end if
            if (norm2(hess_x - expected_hess_x) > tol) then
                write(stderr, *) "test_hess_x_density_qn failed: Incorrect Hessian "// &
                    "linear transformation "//case_name//"."
                test_hess_x_density_qn = .false.
            end if
            if (density_qn%model_stale) then
                write(stderr, *) "test_hess_x_density_qn failed: Hessian model not "// &
                    "rebuilt "//case_name//"."
                test_hess_x_density_qn = .false.
            end if
        end do

        ! deallocate density QN and MO objects
        deallocate(density_qn, mo)

    end function test_hess_x_density_qn

    logical(c_bool) function test_inv_hess_x_density_qn() bind(C)
        !
        ! this function tests the exact, optionally level-shifted inverse of the
        ! approximate Hessian in both the closed- and the open-shell case and for both
        ! orbital bases
        !
        use otr_common, only: error_density_qn
        use otr_common_test_reference, only: n_ao, n_occ, n_particle_ref => n_particle
        use otr_mo_test_reference, only: n_mo
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, &
                                     setup_random_mo_channels, ref_hess_x_static_mo
        use otr_density_qn, only: density_qn_type
        use otr_mo, only: mo_type
        use otr_oao, only: oao_type
        use opentrustregion_unit_tests, only: test_context_type, setup_error_logging, &
                                              identity_matrix

        type(test_context_type), target :: context
        class(density_qn_type), allocatable, target :: density_qn
        type(mo_type), allocatable, target :: mo
        type(oao_type), allocatable, target :: oao
        real(rp), target :: bracket_mo_coeff(n_ao, n_mo, n_particle_ref)
        real(rp) :: bracket_overlap(n_ao, n_ao)
        real(rp), allocatable :: x(:), inv_x(:)
        integer(ip) :: i_basis, n_particle, inv_error, i_channel

        ! assume tests pass
        test_inv_hess_x_density_qn = .true.

        ! the inverse only depends on the low-rank factors and the orbital basis, not
        ! on the QN method which assembled them, so every basis and shell is checked
        ! once
        do i_basis = 1, 2
            do n_particle = 1, n_particle_ref
                if (.not. check_inv_hess_x_density_qn_case(i_basis == 1)) &
                    test_inv_hess_x_density_qn = .false.
            end do
        end do

        ! set up the density QN object in the MO basis with an up-to-date
        ! eigendecomposition in the MO basis itself, whose eigenvalue pairs, the
        ! open-shell differences 2 (e_v - e_o), are all 2, and a low-rank part along
        ! the first parameter with a coupling of -2, for which the matrix I + C P^T
        ! D^-1 E of the Woodbury identity vanishes, and determine if the inverse
        ! returns an error with an error message
        bracket_overlap = generate_random_ao_overlap(n_ao)
        bracket_mo_coeff = &
            generate_random_mo_coeff(bracket_overlap, n_mo, n_particle_ref)
        call setup_density_qn_and_mo_objects(density_qn, mo, bracket_mo_coeff, &
                                             bracket_overlap, n_occ, context)
        do i_channel = 1, mo%n_particle
            associate (channel => mo%mo_channels(i_channel))
                channel%occ_eigvecs = identity_matrix(channel%n_occ)
                channel%occ_eigvals = spread(0.0_rp, 1, channel%n_occ)
                channel%virt_eigvecs = identity_matrix(channel%n_virt)
                channel%virt_eigvals = spread(1.0_rp, 1, channel%n_virt)
            end associate
        end do
        mo%hess_eigen_stale = .false.
        call setup_error_logging(density_qn%settings, context)
        allocate(density_qn%expansion_dirs(mo%n_param, 1), &
                 density_qn%coupling_matrix(1, 1))
        density_qn%expansion_dirs = 0.0_rp
        density_qn%expansion_dirs(1, 1) = 1.0_rp
        density_qn%projection_dirs = density_qn%expansion_dirs
        density_qn%coupling_matrix = -2.0_rp
        density_qn%error_origin = error_density_qn
        allocate(x(mo%n_param), inv_x(mo%n_param))
        call random_number(x)
        call density_qn%inv_hess_x(x, inv_x, inv_error)
        if (inv_error /= error_density_qn + 1) then
            write(stderr, *) "test_inv_hess_x_density_qn failed: No error returned "// &
                "for a singular approximate Hessian."
            test_inv_hess_x_density_qn = .false.
        end if
        if (context%log_message == "") then
            write(stderr, *) "test_inv_hess_x_density_qn failed: No error message "// &
                "printed for a singular approximate Hessian."
            test_inv_hess_x_density_qn = .false.
        end if
        deallocate(density_qn, mo)

    contains

        subroutine check_inv_hess_x_density_qn(hess, vector, mu, case_name, passed)
            !
            ! this subroutine checks the exact-inverse of the density QN Woodbury
            ! preconditioner against a reference Hessian by applying the reference (B -
            ! mu*I) to the preconditioned vector which has to recover the vector, both
            ! with and without the level shift
            !
            real(rp), intent(in) :: hess(:, :), vector(:), mu
            character(len=*), intent(in) :: case_name
            logical, intent(inout) :: passed

            real(rp), allocatable :: actual(:)
            integer(ip) :: error

            allocate(actual(size(vector)))

            ! shifted entry point: applying the reference (B - mu*I) to the
            ! preconditioned vector has to recover the vector
            call density_qn%inv_hess_x(vector, actual, error, mu)
            if (error /= 0) then
                write(stderr, *) "test_inv_hess_x_density_qn failed: Shifted "// &
                    "inverse produced error for the "//case_name//" case."
                passed = .false.
            end if
            if (norm2(matmul(hess, actual) - mu * actual - vector) > &
                tol * (1.0_rp + norm2(hess) * norm2(actual))) then
                write(stderr, *) "test_inv_hess_x_density_qn failed: Shifted "// &
                    "inverse is not inverse of shifted reference Hessian for the "// &
                    case_name//" case."
                passed = .false.
            end if

            ! the unshifted entry point has to reproduce the same round trip against the
            ! reference Hessian itself, without a level shift
            call density_qn%inv_hess_x(vector, actual, error)
            if (error /= 0) then
                write(stderr, *) "test_inv_hess_x_density_qn failed: Unshifted "// &
                    "inverse produced error for the "//case_name//" case."
                passed = .false.
            end if
            if (norm2(matmul(hess, actual) - vector) > &
                tol * (1.0_rp + norm2(hess) * norm2(actual))) then
                write(stderr, *) "test_inv_hess_x_density_qn failed: Inverse is "// &
                    "not inverse of reference Hessian for the "//case_name//" case."
                passed = .false.
            end if

        end subroutine check_inv_hess_x_density_qn

        function check_inv_hess_x_density_qn_case(mo_basis) result(passed)
            !
            ! this function performs the exact-inverse check of the density QN Woodbury
            ! preconditioner for one orbital basis and the current shell with random
            ! expansion and projection directions and a random coupling matrix, and for
            ! a stale model on an empty history, which has to be rebuilt without
            ! low-rank part before it is inverted
            !
            use otr_mo_test_reference, only: n_mo
            use otr_common_test_reference, only: n_ao, n_occ, shell_names
            use opentrustregion_unit_tests, only: identity_matrix, setup_settings
            use otr_oao_unit_tests, only: ref_unpack_asymm, ref_hess_x_oao, &
                                          generate_random_nonredundant_vector
            use otr_density_qn, only: density_qn_oao_type

            logical, intent(in) :: mo_basis
            logical :: passed

            ! shifts of the occupied and virtual orbital energies keep the eigenvalue
            ! pairs of the static part away from zero, while the low-rank part is scaled
            ! down so that it cannot make the approximate Hessian singular
            integer(ip), parameter :: n_dirs = 2
            real(rp), parameter :: mu = 0.1_rp, mo_shift = 1.5_rp, &
                                   oao_shift = 10.0_rp, low_rank_scale = 0.25_rp

            real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
            real(rp) :: ao_overlap(n_ao, n_ao), dm_oao(n_ao, n_ao, n_particle), &
                        fock_oo(n_ao, n_ao, n_particle), &
                        fock_vv(n_ao, n_ao, n_particle), &
                        zero_response(n_ao, n_ao, n_particle), coupling(n_dirs, n_dirs)
            real(rp), allocatable :: static_hess(:, :), hess(:, :), e_i(:), dirs(:, :)
            integer(ip) :: n_param, i, k
            character(len=:), allocatable :: case_name

            ! assume test passes
            passed = .true.
            case_name = trim(merge("MO ", "OAO", mo_basis))//"-basis "// &
                        trim(shell_names(n_particle))

            ! set up the orbital basis with a static part whose eigenvalue pairs are
            ! bounded away from zero
            if (mo_basis) then
                ao_overlap = generate_random_ao_overlap(n_ao)
                mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
                call setup_density_qn_and_mo_objects( &
                    density_qn, mo, mo_coeff, ao_overlap, n_occ(:n_particle), context)
                call setup_random_mo_channels(mo, n_occ(:n_particle), n_mo)
                do k = 1, n_particle
                    associate (channel => mo%mo_channels(k))
                        channel%fock_oo = channel%fock_oo - &
                                          mo_shift * identity_matrix(channel%n_occ)
                        channel%occ_eigvals = channel%occ_eigvals - mo_shift
                        channel%fock_vv = channel%fock_vv + &
                                          mo_shift * identity_matrix(channel%n_virt)
                        channel%virt_eigvals = channel%virt_eigvals + mo_shift
                    end associate
                end do
            else
                call generate_random_fock_partition(n_ao, n_occ(:n_particle), dm_oao, &
                                                    fock_oo, fock_vv)
                do k = 1, n_particle
                    fock_oo(:, :, k) = fock_oo(:, :, k) - oao_shift * dm_oao(:, :, k)
                    fock_vv(:, :, k) = fock_vv(:, :, k) + oao_shift * &
                                       (identity_matrix(n_ao) - dm_oao(:, :, k))
                end do

                ! set up the OAO object, whose eigendecomposition of the static part
                ! the inverse refreshes, and the density QN object in the OAO basis
                allocate(oao)
                call setup_settings(oao%settings, context)
                oao%n_ao = n_ao
                oao%n_particle = n_particle
                oao%n_param = n_particle * n_ao * (n_ao - 1) / 2
                oao%dm_oao = dm_oao
                oao%fock_oo = fock_oo
                oao%fock_vv = fock_vv
                oao%hess_eigen_stale = .true.
                allocate(density_qn, source=density_qn_oao_type(oao))
                call setup_settings(density_qn%settings, context)

                ! cache a wrong eigendecomposition and wrong directions rotated into
                ! its eigenbasis, which have to be refreshed since the static part is
                ! marked as changed
                allocate(oao%hess_eigvecs(n_ao, n_ao, n_particle), &
                         oao%hess_eigvals(n_ao, n_particle), &
                         density_qn%rotated_expansion_dirs(oao%n_param, n_dirs), &
                         density_qn%rotated_projection_dirs(oao%n_param, n_dirs))
                oao%hess_eigvecs = 0.0_rp
                oao%hess_eigvals = 1.0_rp
                density_qn%rotated_expansion_dirs = 0.0_rp
                density_qn%rotated_projection_dirs = 0.0_rp
            end if
            n_param = density_qn%orbitals%n_param

            ! build the dense static part
            allocate(static_hess(n_param, n_param), e_i(n_param))
            zero_response = 0.0_rp
            do i = 1, n_param
                e_i = 0.0_rp
                e_i(i) = 1.0_rp
                if (mo_basis) then
                    static_hess(:, i) = ref_hess_x_static_mo(e_i, mo%mo_channels)
                else
                    static_hess(:, i) = ref_hess_x_oao( &
                        ref_unpack_asymm(e_i, n_particle, n_ao), zero_response, &
                        dm_oao, fock_oo, fock_vv, n_param)
                end if
            end do

            ! random expansion and projection directions and a random vector, which are
            ! confined to the non-redundant subspace in the OAO basis
            allocate(dirs(n_param, 2 * n_dirs + 1))
            do i = 1, 2 * n_dirs + 1
                if (mo_basis) then
                    call random_number(dirs(:, i))
                    dirs(:, i) = dirs(:, i) - 0.5_rp
                else
                    dirs(:, i) = generate_random_nonredundant_vector( &
                        n_param, n_particle, n_ao, dm_oao)
                end if
            end do
            density_qn%expansion_dirs = dirs(:, :n_dirs)
            density_qn%projection_dirs = dirs(:, n_dirs + 1:2 * n_dirs)

            ! random, and hence non-symmetric, coupling matrix, so that a transposed
            ! coupling matrix would not go unnoticed
            call random_number(coupling)
            density_qn%coupling_matrix = low_rank_scale * (coupling - 0.5_rp)

            ! check inverse against the dense approximate Hessian
            hess = static_hess + matmul(density_qn%expansion_dirs, matmul( &
                density_qn%coupling_matrix, transpose(density_qn%projection_dirs)))
            call check_inv_hess_x_density_qn(hess, dirs(:, 2 * n_dirs + 1), mu, &
                                             case_name, passed)
            if (.not. (allocated(density_qn%rotated_expansion_dirs) .and. &
                       allocated(density_qn%rotated_projection_dirs))) then
                write(stderr, *) "test_inv_hess_x_density_qn failed: Directions "// &
                    "rotated into the eigenbasis not kept for the "//case_name//" case."
                passed = .false.
            end if

            ! mark the model stale on an empty history and determine if it is rebuilt
            ! before it is inverted, which discards the low-rank part and leaves the
            ! inverse of the static part
            call setup_empty_history(density_qn, n_particle)
            density_qn%model_stale = .true.
            call check_inv_hess_x_density_qn(static_hess, dirs(:, 2 * n_dirs + 1), mu, &
                                             "stale-model "//case_name, passed)
            if (density_qn%model_stale .or. allocated(density_qn%coupling_matrix)) then
                write(stderr, *) "test_inv_hess_x_density_qn failed: Stale Hessian "// &
                    "model not rebuilt before it is inverted for the "//case_name// &
                    " case."
                passed = .false.
            end if

            ! deallocate density QN and orbital objects
            deallocate(density_qn)
            if (allocated(mo)) deallocate(mo)
            if (allocated(oao)) deallocate(oao)

        end function check_inv_hess_x_density_qn_case

    end function test_inv_hess_x_density_qn

    logical(c_bool) function test_precond_density_qn() bind(C)
        !
        ! this function tests the subroutine which defines the level-shifted
        ! preconditioner of the density QN approximate Hessian, which hands the
        ! residual and the level shift to the inverse of the approximate Hessian and
        ! passes on what it returns
        !
        use otr_common_unit_tests, only: level_context
        use otr_mo_test_reference, only: n_param

        real(rp), parameter :: mu = 0.3_rp

        type(recording_density_qn_type), target :: density_qn
        real(rp) :: residual(n_param), precond_residual(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_precond_density_qn = .true.

        ! call routine for an inverse which fails and determine if the residual and the
        ! level shift are handed to it and its result and error are passed on
        call random_number(residual)
        density_qn%error = 1
        call density_qn%precond(residual, mu, precond_residual, error, &
                                level_context(density_qn, 1_ip))
        if (.not. allocated(density_qn%received)) then
            write(stderr, *) "test_precond_density_qn failed: Inverse not called."
            test_precond_density_qn = .false.
        else if (any(shape(density_qn%received) /= [n_param])) then
            write(stderr, *) "test_precond_density_qn failed: Residual not handed "// &
                "to the inverse."
            test_precond_density_qn = .false.
        else if (norm2(density_qn%received - residual) > tol) then
            write(stderr, *) "test_precond_density_qn failed: Residual not handed "// &
                "to the inverse."
            test_precond_density_qn = .false.
        end if
        if (.not. density_qn%level_shift_received .or. &
            abs(density_qn%received_level_shift - mu) > tol) then
            write(stderr, *) "test_precond_density_qn failed: Level shift not "// &
                "handed to the inverse."
            test_precond_density_qn = .false.
        end if
        if (norm2(precond_residual - 3.0_rp * residual) > tol) then
            write(stderr, *) "test_precond_density_qn failed: Result of the "// &
                "inverse not passed on."
            test_precond_density_qn = .false.
        end if
        if (error /= 1) then
            write(stderr, *) "test_precond_density_qn failed: Error of the inverse "// &
                "not passed on."
            test_precond_density_qn = .false.
        end if

    end function test_precond_density_qn

    logical(c_bool) function test_precond_pd_density_qn() bind(C)
        !
        ! this function tests the subroutine which defines the positive-definite
        ! preconditioner of the density QN approximate Hessian, which hands the
        ! residual and the context to the one of the orbital basis and passes on what
        ! it returns
        !
        use otr_common_unit_tests, only: level_context
        use otr_mo_test_reference, only: n_param

        type(recording_density_qn_type), target :: density_qn
        type(recording_basis_type), target :: basis
        class(*), pointer :: context
        real(rp) :: residual(n_param), precond_residual(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_precond_pd_density_qn = .true.

        ! point the density QN object at a recording orbital basis whose preconditioner
        ! fails
        density_qn%orbitals => basis
        basis%error = 1
        context => level_context(density_qn, 1_ip)
        call random_number(residual)

        ! call routine and determine if the residual and the context are handed to the
        ! orbital basis and its result and error are passed on
        call density_qn%precond_pd(residual, precond_residual, error, context)
        if (.not. allocated(basis%received)) then
            write(stderr, *) "test_precond_pd_density_qn failed: Preconditioner of "// &
                "the orbital basis not called."
            test_precond_pd_density_qn = .false.
        else if (any(shape(basis%received) /= [n_param])) then
            write(stderr, *) "test_precond_pd_density_qn failed: Residual not "// &
                "handed to the orbital basis."
            test_precond_pd_density_qn = .false.
        else if (norm2(basis%received - residual) > tol) then
            write(stderr, *) "test_precond_pd_density_qn failed: Residual not "// &
                "handed to the orbital basis."
            test_precond_pd_density_qn = .false.
        end if
        if (.not. associated(basis%received_context, context)) then
            write(stderr, *) "test_precond_pd_density_qn failed: Context not "// &
                "handed to the orbital basis."
            test_precond_pd_density_qn = .false.
        end if
        if (norm2(precond_residual - 3.0_rp * residual) > tol) then
            write(stderr, *) "test_precond_pd_density_qn failed: Result of the "// &
                "orbital basis not passed on."
            test_precond_pd_density_qn = .false.
        end if
        if (error /= 1) then
            write(stderr, *) "test_precond_pd_density_qn failed: Error of the "// &
                "orbital basis not passed on."
            test_precond_pd_density_qn = .false.
        end if

    end function test_precond_pd_density_qn

    logical(c_bool) function test_get_extra_trial_vectors_density_qn() bind(C)
        !
        ! this function tests the subroutine which returns the extra trial vectors of
        ! the orbital basis for the solver's initial trial space, which hands the
        ! context to the orbital basis and passes on what it returns
        !
        use otr_common_unit_tests, only: level_context
        use otr_mo_test_reference, only: n_param

        integer(ip), parameter :: n_extra = 2

        type(recording_density_qn_type), target :: density_qn
        type(recording_basis_type), target :: basis
        class(*), pointer :: context
        real(rp) :: trial_vectors(n_param, n_extra)
        integer(ip) :: error

        ! assume tests pass
        test_get_extra_trial_vectors_density_qn = .true.

        ! point the density QN object at a recording orbital basis whose extra trial
        ! vector function fails
        density_qn%orbitals => basis
        basis%error = 1
        context => level_context(density_qn, 1_ip)

        ! call routine and determine if the context is handed to the orbital basis and
        ! its trial vectors and error are passed on
        trial_vectors = 0.0_rp
        call density_qn%get_extra_trial_vectors(trial_vectors, error, context)
        if (.not. associated(basis%received_context, context)) then
            write(stderr, *) "test_get_extra_trial_vectors_density_qn failed: "// &
                "Context not handed to the orbital basis."
            test_get_extra_trial_vectors_density_qn = .false.
        end if
        if (any(abs(trial_vectors - 3.0_rp) > tol)) then
            write(stderr, *) "test_get_extra_trial_vectors_density_qn failed: "// &
                "Trial vectors of the orbital basis not passed on."
            test_get_extra_trial_vectors_density_qn = .false.
        end if
        if (error /= 1) then
            write(stderr, *) "test_get_extra_trial_vectors_density_qn failed: "// &
                "Error of the orbital basis not passed on."
            test_get_extra_trial_vectors_density_qn = .false.
        end if

    end function test_get_extra_trial_vectors_density_qn

    logical(c_bool) function test_project_density_qn_oao() bind(C)
        !
        ! this function tests the subroutine which applies the projection of the
        ! orbitals in the OAO basis, which hands the vector and the context to the
        ! projection of the orbital basis and passes on what it returns
        !
        use otr_density_qn, only: density_qn_oao_type
        use otr_common_unit_tests, only: level_context
        use otr_oao_test_reference, only: n_param

        type(density_qn_oao_type), target :: density_qn
        type(recording_basis_type), target :: basis
        class(*), pointer :: context
        real(rp) :: vector(n_param), original(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_project_density_qn_oao = .true.

        ! point the density QN object at a recording orbital basis whose projection
        ! fails
        density_qn%orbitals => basis
        basis%error = 1
        context => level_context(density_qn, 1_ip)
        call random_number(vector)
        original = vector

        ! call routine and determine if the vector and the context are handed to the
        ! orbital basis and its projected vector and error are passed on
        call density_qn%project(vector, error, context)
        if (.not. allocated(basis%received)) then
            write(stderr, *) "test_project_density_qn_oao failed: Projection of "// &
                "the orbital basis not called."
            test_project_density_qn_oao = .false.
        else if (any(shape(basis%received) /= [n_param])) then
            write(stderr, *) "test_project_density_qn_oao failed: Vector not "// &
                "handed to the orbital basis."
            test_project_density_qn_oao = .false.
        else if (norm2(basis%received - original) > tol) then
            write(stderr, *) "test_project_density_qn_oao failed: Vector not "// &
                "handed to the orbital basis."
            test_project_density_qn_oao = .false.
        end if
        if (.not. associated(basis%received_context, context)) then
            write(stderr, *) "test_project_density_qn_oao failed: Context not "// &
                "handed to the orbital basis."
            test_project_density_qn_oao = .false.
        end if
        if (norm2(vector - 3.0_rp * original) > tol) then
            write(stderr, *) "test_project_density_qn_oao failed: Projected vector "// &
                "of the orbital basis not passed on."
            test_project_density_qn_oao = .false.
        end if
        if (error /= 1) then
            write(stderr, *) "test_project_density_qn_oao failed: Error of the "// &
                "orbital basis not passed on."
            test_project_density_qn_oao = .false.
        end if

    end function test_project_density_qn_oao

    logical(c_bool) function test_construct_density_qn_mo() bind(C)
        !
        ! this function tests the function which returns the density QN object for the
        ! MO basis pointing to an MO object and its quantities
        !
        use otr_density_qn, only: density_qn_mo_type
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object

        type(mo_type), allocatable, target :: mo
        type(density_qn_mo_type) :: density_qn
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)

        ! assume tests pass
        test_construct_density_qn_mo = .true.

        ! set up the MO object
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_mo_object(mo, mo_coeff, ao_overlap, n_occ)

        ! call routine and determine if the quantities of the MO object are associated
        density_qn = density_qn_mo_type(mo)
        if (.not. associated(density_qn%orbitals, mo)) then
            write(stderr, *) "test_construct_density_qn_mo failed: Orbitals not "// &
                "associated."
            test_construct_density_qn_mo = .false.
        end if
        if (.not. associated(density_qn%mo_coeff, mo_coeff)) then
            write(stderr, *) "test_construct_density_qn_mo failed: MO coefficients "// &
                "not associated."
            test_construct_density_qn_mo = .false.
        end if
        if (.not. associated(density_qn%ao_overlap, mo%ao_overlap)) then
            write(stderr, *) "test_construct_density_qn_mo failed: AO overlap "// &
                "matrix not associated."
            test_construct_density_qn_mo = .false.
        end if
        if (.not. associated(density_qn%mo_channels, mo%mo_channels)) then
            write(stderr, *) "test_construct_density_qn_mo failed: Particle "// &
                "channels not associated."
            test_construct_density_qn_mo = .false.
        end if

        ! determine if the number of parameters is taken over and the roles the object
        ! performs itself are set, without a projection
        if (density_qn%n_param /= mo%n_param) then
            write(stderr, *) "test_construct_density_qn_mo failed: Number of "// &
                "parameters not taken over."
            test_construct_density_qn_mo = .false.
        end if
        if (.not. all([density_qn%provides%update_orbs, density_qn%provides%obj_func, &
                       density_qn%provides%hess_x, density_qn%provides%precond, &
                       density_qn%provides%precond_pd, &
                       density_qn%provides%get_extra_trial_vectors]) .or. &
            density_qn%provides%project) then
            write(stderr, *) "test_construct_density_qn_mo failed: Roles the "// &
                "object performs itself set wrongly."
            test_construct_density_qn_mo = .false.
        end if

        ! call routine for an MO object without MO coefficients, overlap matrix and
        ! particle channels and determine if only the orbitals are associated
        deallocate(mo)
        allocate(mo)
        density_qn = density_qn_mo_type(mo)
        if (.not. associated(density_qn%orbitals, mo)) then
            write(stderr, *) "test_construct_density_qn_mo failed: Orbitals not "// &
                "associated for an empty MO object."
            test_construct_density_qn_mo = .false.
        end if
        if (associated(density_qn%mo_coeff) .or. &
            associated(density_qn%ao_overlap) .or. associated(density_qn%mo_channels)) &
            then
            write(stderr, *) "test_construct_density_qn_mo failed: Missing "// &
                "quantities associated for an empty MO object."
            test_construct_density_qn_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo)

    end function test_construct_density_qn_mo

    logical(c_bool) function test_construct_density_qn_oao() bind(C)
        !
        ! this function tests the function which returns the density QN object for the
        ! OAO basis pointing to an OAO object and its quantities
        !
        use otr_density_qn, only: density_qn_oao_type
        use otr_oao, only: oao_type
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle

        type(oao_type), allocatable, target :: oao
        type(density_qn_oao_type) :: density_qn

        ! assume tests pass
        test_construct_density_qn_oao = .true.

        ! set up an OAO object whose basis-specific quantities are allocated except for
        ! the Fock matrix blocks
        allocate(oao)
        oao%n_param = n_param
        allocate(oao%s_inv_sqrt(n_ao, n_ao), oao%dm_oao(n_ao, n_ao, n_particle))

        ! call routine and determine if the allocated quantities are associated and the
        ! others are not
        density_qn = density_qn_oao_type(oao)
        if (.not. associated(density_qn%orbitals, oao)) then
            write(stderr, *) "test_construct_density_qn_oao failed: Orbitals not "// &
                "associated."
            test_construct_density_qn_oao = .false.
        end if
        if (.not. associated(density_qn%s_inv_sqrt, oao%s_inv_sqrt)) then
            write(stderr, *) "test_construct_density_qn_oao failed: Inverse square "// &
                "root of the overlap matrix not associated."
            test_construct_density_qn_oao = .false.
        end if
        if (.not. associated(density_qn%dm_oao, oao%dm_oao)) then
            write(stderr, *) "test_construct_density_qn_oao failed: Density matrix "// &
                "in the OAO basis not associated."
            test_construct_density_qn_oao = .false.
        end if
        if (associated(density_qn%fock_oo) .or. associated(density_qn%fock_vv)) then
            write(stderr, *) "test_construct_density_qn_oao failed: Unallocated "// &
                "Fock matrix blocks associated."
            test_construct_density_qn_oao = .false.
        end if

        ! determine if the number of parameters is taken over and the roles the object
        ! performs itself are set, including the projection
        if (density_qn%n_param /= n_param) then
            write(stderr, *) "test_construct_density_qn_oao failed: Number of "// &
                "parameters not taken over."
            test_construct_density_qn_oao = .false.
        end if
        if (.not. all([density_qn%provides%update_orbs, density_qn%provides%obj_func, &
                       density_qn%provides%hess_x, density_qn%provides%precond, &
                       density_qn%provides%precond_pd, density_qn%provides%project, &
                       density_qn%provides%get_extra_trial_vectors])) then
            write(stderr, *) "test_construct_density_qn_oao failed: Roles the "// &
                "object performs itself not set."
            test_construct_density_qn_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao)

    end function test_construct_density_qn_oao

    logical(c_bool) function test_rotate_trial_density_qn_mo() bind(C)
        !
        ! this function tests the subroutine which returns the density matrix rotated
        ! by an orbital rotation in the AO basis and in the history basis, S D S, for
        ! the MO basis without moving the current orbitals
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use otr_common, only: error_density_qn
        use otr_density_qn, only: density_qn_mo_type
        use otr_mo, only: mo_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_mo_test_reference, only: n_mo, n_param
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object, &
                                     ref_rotate_mo_coeff

        type(test_context_type), target :: context
        type(mo_type), allocatable, target :: mo
        type(density_qn_mo_type) :: density_qn
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), mo_coeff_before(n_ao, n_mo, n_particle), &
                    expected(n_ao, n_mo, n_particle), kappa(n_param), &
                    rot_dm_ao(n_ao, n_ao, n_particle), &
                    rot_dm_hist(n_ao, n_ao, n_particle), expected_dm(n_ao, n_ao)
        integer(ip) :: error, k

        ! assume tests pass
        test_rotate_trial_density_qn_mo = .true.

        ! set up the MO object with random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        mo_coeff_before = mo_coeff
        call setup_mo_object(mo, mo_coeff, ao_overlap, n_occ)
        density_qn = density_qn_mo_type(mo)
        density_qn%error_origin = error_density_qn
        call setup_settings(density_qn%settings, context)

        ! call routine and determine if the density matrix of the independently rotated
        ! orbitals is returned in the AO basis and multiplied by the AO overlap matrix
        ! from both sides, while the current orbitals are left untouched
        call random_number(kappa)
        kappa = 0.2_rp * (kappa - 0.5_rp)
        expected = ref_rotate_mo_coeff(kappa, mo_coeff, mo%mo_channels)
        call density_qn%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_trial_density_qn_mo failed: Produced error."
            test_rotate_trial_density_qn_mo = .false.
        end if
        do k = 1, n_particle
            expected_dm = matmul(expected(:, :n_occ(k), k), &
                                 transpose(expected(:, :n_occ(k), k)))
            if (norm2(rot_dm_ao(:, :, k) - expected_dm) > tol) then
                write(stderr, *) "test_rotate_trial_density_qn_mo failed: "// &
                    "Incorrect rotated density matrix in the AO basis."
                test_rotate_trial_density_qn_mo = .false.
            end if
            if (norm2(rot_dm_hist(:, :, k) - &
                      matmul(ao_overlap, matmul(expected_dm, ao_overlap))) > tol) then
                write(stderr, *) "test_rotate_trial_density_qn_mo failed: "// &
                    "Incorrect rotated density matrix in the history basis."
                test_rotate_trial_density_qn_mo = .false.
            end if
        end do
        if (norm2(mo_coeff - mo_coeff_before) > tol) then
            write(stderr, *) "test_rotate_trial_density_qn_mo failed: Current MO "// &
                "coefficients changed."
            test_rotate_trial_density_qn_mo = .false.
        end if

        ! call routine with a rotation that is not finite, which the matrix exponential
        ! refuses, and determine if its error is returned with the origin of the object
        kappa(1) = ieee_value(1.0_rp, ieee_quiet_nan)
        call density_qn%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= error_density_qn + 1) then
            write(stderr, *) "test_rotate_trial_density_qn_mo failed: Error for a "// &
                "rotation that is not finite not returned with the origin of the "// &
                "object."
            test_rotate_trial_density_qn_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo)

    end function test_rotate_trial_density_qn_mo

    logical(c_bool) function test_rotate_trial_density_qn_oao() bind(C)
        !
        ! this function tests the subroutine which returns the density matrix rotated
        ! by an orbital rotation in the AO and OAO basis without moving the current
        ! orbitals
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use otr_common, only: error_density_qn
        use otr_density_qn, only: density_qn_oao_type
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_occ
        use otr_common_unit_tests, only: generate_random_density_matrix
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix

        type(test_context_type), target :: context
        integer(ip), parameter :: n_param_oao = n_ao * (n_ao - 1) / 2
        real(rp), parameter :: basis_scale = 1.5_rp, angle = 0.3_rp

        type(oao_type), allocatable, target :: oao
        type(density_qn_oao_type) :: density_qn
        real(rp), target :: dm_ao(n_ao, n_ao, 1)
        real(rp) :: dm_oao(n_ao, n_ao, 1), kappa(n_param_oao), &
                    rot_dm_ao(n_ao, n_ao, 1), rot_dm_hist(n_ao, n_ao, 1), &
                    rotation(n_ao, n_ao)
        integer(ip) :: error

        ! assume tests pass
        test_rotate_trial_density_qn_oao = .true.

        ! set up the OAO object with an AO basis whose overlap has the inverse square
        ! root c I, so that the density matrix in the AO basis is c^2 times the one in
        ! the OAO basis
        allocate(oao)
        call setup_settings(oao%settings, context)
        oao%n_ao = n_ao
        oao%n_particle = 1
        oao%n_param = n_param_oao
        oao%s_inv_sqrt = basis_scale * identity_matrix(n_ao)
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_occ(1))
        dm_ao = basis_scale**2 * dm_oao
        oao%dm_ao => dm_ao
        oao%dm_oao = dm_oao
        density_qn = density_qn_oao_type(oao)
        density_qn%error_origin = error_density_qn
        call setup_settings(density_qn%settings, context)

        ! call routine and determine if the density matrix is rotated consistently in
        ! both bases while the current density matrix is left untouched
        call random_number(kappa)
        kappa = 0.2_rp * kappa
        call density_qn%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Produced error."
            test_rotate_trial_density_qn_oao = .false.
        end if
        if (norm2(rot_dm_hist - dm_oao) < tol) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Density "// &
                "matrix not rotated."
            test_rotate_trial_density_qn_oao = .false.
        end if
        if (norm2(rot_dm_ao - basis_scale**2 * rot_dm_hist) > tol) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Density "// &
                "matrix not rotated consistently in the AO and OAO basis."
            test_rotate_trial_density_qn_oao = .false.
        end if
        if (norm2(matmul(rot_dm_hist(:, :, 1), rot_dm_hist(:, :, 1)) - &
                  rot_dm_hist(:, :, 1)) > tol) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Rotated "// &
                "density matrix not idempotent."
            test_rotate_trial_density_qn_oao = .false.
        end if
        if (norm2(oao%dm_oao - dm_oao) > tol) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Current "// &
                "density matrix in the OAO basis moved."
            test_rotate_trial_density_qn_oao = .false.
        end if
        if (norm2(dm_ao - basis_scale**2 * dm_oao) > tol) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Current "// &
                "density matrix in the AO basis moved."
            test_rotate_trial_density_qn_oao = .false.
        end if

        ! call routine with a rotation of the first pair of orbitals, whose rotation
        ! matrix is a plane rotation, and determine if the density matrix is rotated as
        ! expected; the closed form of the rotation catches a sign-flipped orbital
        ! rotation, which the consistency and idempotency checks above do not
        kappa = 0.0_rp
        kappa(1) = angle
        rotation = identity_matrix(n_ao)
        rotation(1:2, 1:2) = reshape([cos(angle), -sin(angle), &
                                      sin(angle), cos(angle)], [2, 2])
        call density_qn%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Produced "// &
                "error for a plane rotation."
            test_rotate_trial_density_qn_oao = .false.
        end if
        if (norm2(rot_dm_hist(:, :, 1) - matmul( &
            transpose(rotation), matmul(dm_oao(:, :, 1), rotation))) > tol) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Incorrect "// &
                "density matrix for a plane rotation."
            test_rotate_trial_density_qn_oao = .false.
        end if

        ! call routine with a rotation that is not finite, which the matrix exponential
        ! refuses, and determine if its error is returned with the origin of the object
        kappa(1) = ieee_value(1.0_rp, ieee_quiet_nan)
        call density_qn%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= error_density_qn + 1) then
            write(stderr, *) "test_rotate_trial_density_qn_oao failed: Error for a "// &
                "rotation that is not finite not returned with the origin of the "// &
                "object."
            test_rotate_trial_density_qn_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao)

    end function test_rotate_trial_density_qn_oao

    logical(c_bool) function test_history_dm_density_qn_mo() bind(C)
        !
        ! this function tests the function which returns the current density matrix in
        ! the history basis for the MO basis
        !
        use otr_density_qn, only: density_qn_mo_type
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object

        type(mo_type), allocatable, target :: mo
        type(density_qn_mo_type) :: density_qn
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), dm(n_ao, n_ao, n_particle)
        integer(ip) :: k

        ! assume tests pass
        test_history_dm_density_qn_mo = .true.

        ! set up the MO object with random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_mo_object(mo, mo_coeff, ao_overlap, n_occ)
        density_qn = density_qn_mo_type(mo)

        ! call routine and determine if the density matrix multiplied by the AO overlap
        ! matrix from both sides is returned
        dm = density_qn%history_dm()
        do k = 1, n_particle
            if (norm2(dm(:, :, k) - matmul(ao_overlap, matmul(mo%dm_ao(:, :, k), &
                                                              ao_overlap))) > tol) then
                write(stderr, *) "test_history_dm_density_qn_mo failed: Incorrect "// &
                    "density matrix."
                test_history_dm_density_qn_mo = .false.
            end if
        end do

        ! deallocate MO object
        deallocate(mo)

    end function test_history_dm_density_qn_mo

    logical(c_bool) function test_history_dm_density_qn_oao() bind(C)
        !
        ! this function tests the function which returns the current density matrix in
        ! the history basis for the OAO basis
        !
        use otr_density_qn, only: density_qn_oao_type
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_particle

        type(oao_type), allocatable, target :: oao
        type(density_qn_oao_type) :: density_qn
        real(rp) :: dm_oao(n_ao, n_ao, n_particle)

        ! assume tests pass
        test_history_dm_density_qn_oao = .true.

        ! set up the OAO object with a random density matrix
        call random_number(dm_oao)
        allocate(oao)
        oao%dm_oao = dm_oao
        density_qn = density_qn_oao_type(oao)

        ! call routine and determine if the density matrix in the OAO basis is returned
        if (norm2(density_qn%history_dm() - dm_oao) > tol) then
            write(stderr, *) "test_history_dm_density_qn_oao failed: Incorrect "// &
                "density matrix."
            test_history_dm_density_qn_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao)

    end function test_history_dm_density_qn_oao

    logical(c_bool) function test_to_history_basis_density_qn_mo() bind(C)
        !
        ! this function tests the function which transforms a potential from the AO
        ! basis to the history basis for the MO basis
        !
        use otr_density_qn, only: density_qn_mo_type
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use opentrustregion_unit_tests, only: generate_random_symm_matrix
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object

        type(mo_type), allocatable, target :: mo
        type(density_qn_mo_type) :: density_qn
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), matrix(n_ao, n_ao, n_particle)
        integer(ip) :: k

        ! assume tests pass
        test_to_history_basis_density_qn_mo = .true.

        ! set up the MO object with random orthonormal MO coefficients and generate
        ! random potentials
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_mo_object(mo, mo_coeff, ao_overlap, n_occ)
        density_qn = density_qn_mo_type(mo)
        do k = 1, n_particle
            matrix(:, :, k) = generate_random_symm_matrix(n_ao)
        end do

        ! call routine and determine if the potentials are kept in the AO basis
        if (norm2(density_qn%to_history_basis(matrix) - matrix) > tol) then
            write(stderr, *) "test_to_history_basis_density_qn_mo failed: "// &
                "Incorrect potential."
            test_to_history_basis_density_qn_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo)

    end function test_to_history_basis_density_qn_mo

    logical(c_bool) function test_to_history_basis_density_qn_oao() bind(C)
        !
        ! this function tests the function which transforms a potential from the AO
        ! basis to the history basis for the OAO basis
        !
        use otr_density_qn, only: density_qn_oao_type
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_particle
        use opentrustregion_unit_tests, only: generate_random_symm_matrix

        type(oao_type), allocatable, target :: oao
        type(density_qn_oao_type) :: density_qn
        real(rp) :: s_inv_sqrt(n_ao, n_ao), matrix(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle)
        integer(ip) :: k

        ! assume tests pass
        test_to_history_basis_density_qn_oao = .true.

        ! set up the OAO object with a random inverse square root of the overlap matrix
        ! and independently transform random potentials
        s_inv_sqrt = generate_random_symm_matrix(n_ao)
        do k = 1, n_particle
            matrix(:, :, k) = generate_random_symm_matrix(n_ao)
            expected(:, :, k) = matmul(s_inv_sqrt, matmul(matrix(:, :, k), s_inv_sqrt))
        end do
        allocate(oao)
        oao%s_inv_sqrt = s_inv_sqrt
        density_qn = density_qn_oao_type(oao)

        ! call routine and determine if the potentials are transformed with the inverse
        ! square root of the overlap matrix
        if (norm2(density_qn%to_history_basis(matrix) - expected) > tol) then
            write(stderr, *) "test_to_history_basis_density_qn_oao failed: "// &
                "Incorrect potential."
            test_to_history_basis_density_qn_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao)

    end function test_to_history_basis_density_qn_oao

    logical(c_bool) function test_history_columns_density_qn_mo() bind(C)
        !
        ! this function tests the subroutine which expresses a history of difference
        ! matrices as history columns and as packed columns in the non-redundant
        ! parameter space for the MO basis
        !
        use otr_density_qn, only: density_qn_mo_type
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_particle
        use opentrustregion_unit_tests, only: generate_random_symm_matrix
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object, &
                                     ref_mo_transform, ref_pack_ov

        integer(ip), parameter :: n_list = 2

        type(mo_type), allocatable, target :: mo
        type(density_qn_mo_type) :: density_qn
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        real(rp), allocatable :: diff(:, :, :, :), diff_mo(:, :, :), cols(:, :), &
                                 packed(:, :), expected_cols(:, :), &
                                 expected_packed(:, :)
        integer(ip), allocatable :: n_occ(:)
        integer(ip) :: i_case, n_part, j, k
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_history_columns_density_qn_mo = .true.

        ! loop over every occupation case: the history columns are the full difference
        ! matrices in the current MO basis and the packed columns the entries of their
        ! occupied-virtual blocks between orbitals of the same irrep
        do i_case = 1, n_cases
            case_name = trim(case_names(i_case))
            n_part = case_n_particle(i_case)
            n_occ = case_n_occ(:n_part, i_case)
            ao_overlap = generate_random_ao_overlap(n_ao)
            mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
            call setup_mo_object(mo, mo_coeff(:, :, :n_part), ao_overlap, n_occ, &
                                 case_irreps(:, :n_part, i_case))
            density_qn = density_qn_mo_type(mo)
            allocate(diff(n_ao, n_ao, n_part, n_list), diff_mo(n_mo, n_mo, n_part), &
                     expected_cols(n_part * n_mo**2, n_list), &
                     expected_packed(mo%n_param, n_list))
            do j = 1, n_part
                do k = 1, n_list
                    diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
                end do
            end do

            ! independently transform the difference matrices
            do k = 1, n_list
                do j = 1, n_part
                    diff_mo(:, :, j) = ref_mo_transform(mo_coeff(:, :, j), &
                                                        diff(:, :, j, k))
                    expected_cols((j - 1) * n_mo**2 + 1:j * n_mo**2, k) = &
                        reshape(diff_mo(:, :, j), [n_mo**2])
                end do
                expected_packed(:, k) = ref_pack_ov(diff_mo, mo%mo_channels)
            end do

            ! call routine and determine if both sets of columns are correct
            call density_qn%history_columns(diff, cols, packed)
            if (any(shape(cols) /= shape(expected_cols))) then
                write(stderr, *) "test_history_columns_density_qn_mo failed: "// &
                    "Incorrect dimensions of history columns in the "//case_name// &
                    " case."
                test_history_columns_density_qn_mo = .false.
            else if (norm2(cols - expected_cols) > tol) then
                write(stderr, *) "test_history_columns_density_qn_mo failed: "// &
                    "Incorrect history columns in the "//case_name//" case."
                test_history_columns_density_qn_mo = .false.
            end if
            if (any(shape(packed) /= shape(expected_packed))) then
                write(stderr, *) "test_history_columns_density_qn_mo failed: "// &
                    "Incorrect dimensions of packed columns in the "//case_name// &
                    " case."
                test_history_columns_density_qn_mo = .false.
            else if (norm2(packed - expected_packed) > tol) then
                write(stderr, *) "test_history_columns_density_qn_mo failed: "// &
                    "Incorrect packed columns in the "//case_name//" case."
                test_history_columns_density_qn_mo = .false.
            end if
            deallocate(diff, diff_mo, expected_cols, expected_packed, mo)
        end do

    end function test_history_columns_density_qn_mo

    logical(c_bool) function test_history_columns_density_qn_oao() bind(C)
        !
        ! this function tests the subroutine which expresses a history of difference
        ! matrices as history columns and as packed columns in the non-redundant
        ! parameter space for the OAO basis
        !
        use otr_density_qn, only: density_qn_oao_type
        use otr_oao, only: oao_type
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle, n_occ
        use opentrustregion_unit_tests, only: generate_random_symm_matrix
        use otr_common_unit_tests, only: generate_random_density_matrix

        integer(ip), parameter :: n_list = 2

        type(oao_type), allocatable, target :: oao
        type(density_qn_oao_type) :: density_qn
        real(rp) :: dm_oao(n_ao, n_ao, n_particle), diff(n_ao, n_ao, n_particle, n_list)
        real(rp), allocatable :: cols(:, :), packed(:, :)
        integer(ip) :: j, k

        ! assume tests pass
        test_history_columns_density_qn_oao = .true.

        ! set up the OAO object with a random density matrix, which the packed columns
        ! are projected with, and generate random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_occ(j))
            do k = 1, n_list
                diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do
        allocate(oao)
        oao%n_ao = n_ao
        oao%n_particle = n_particle
        oao%n_param = n_param
        oao%dm_oao = dm_oao
        density_qn = density_qn_oao_type(oao)

        ! call routine and determine if the history columns are the flattened
        ! difference matrices and the packed columns their independently reproduced
        ! projections
        call density_qn%history_columns(diff, cols, packed)
        if (size(cols, 1) /= n_ao * n_ao * n_particle .or. size(cols, 2) /= n_list) then
            write(stderr, *) "test_history_columns_density_qn_oao failed: "// &
                "Incorrect dimensions of history columns."
            test_history_columns_density_qn_oao = .false.
        else if ( &
            norm2(cols - reshape(diff, [n_ao * n_ao * n_particle, n_list])) > tol) then
            write(stderr, *) "test_history_columns_density_qn_oao failed: History "// &
                "columns are not the flattened difference matrices."
            test_history_columns_density_qn_oao = .false.
        end if
        if (size(packed, 1) /= n_param .or. size(packed, 2) /= n_list) then
            write(stderr, *) "test_history_columns_density_qn_oao failed: "// &
                "Incorrect dimensions of packed columns."
            test_history_columns_density_qn_oao = .false.
        else if (norm2(packed - ref_packed_columns(diff)) > tol) then
            write(stderr, *) "test_history_columns_density_qn_oao failed: Packed "// &
                "columns are not the projected and packed difference matrices."
            test_history_columns_density_qn_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao)

    contains

        function ref_packed_columns(v_diff) result(dirs)
            !
            ! this function independently projects and packs every column of a
            ! history-difference array into the non-redundant subspace of the current
            ! density matrix
            !
            use otr_oao_unit_tests, only: ref_project_asymm, ref_pack_asymm

            real(rp), intent(in) :: v_diff(:, :, :, :)
            real(rp) :: dirs(n_param, size(v_diff, 4))

            integer(ip) :: col

            do col = 1, size(v_diff, 4, kind=ip)
                dirs(:, col) = ref_pack_asymm( &
                    ref_project_asymm(v_diff(:, :, :, col), dm_oao), n_param)
            end do

        end function ref_packed_columns

    end function test_history_columns_density_qn_oao

    logical(c_bool) function test_history_channel_rows_density_qn_mo() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in a history column for the MO basis
        !
        use otr_density_qn, only: density_qn_mo_type
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_particle

        type(mo_type), allocatable, target :: mo
        type(density_qn_mo_type) :: density_qn
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        integer(ip), allocatable :: rows(:, :)
        integer(ip) :: k

        ! assume tests pass
        test_history_channel_rows_density_qn_mo = .true.

        ! set up the MO object with the number of particle channels and the MO
        ! coefficients, whose shape is all the routine reads of them
        allocate(mo)
        mo%n_particle = n_particle
        mo%mo_coeff => mo_coeff
        density_qn = density_qn_mo_type(mo)

        ! call routine and determine if every channel holds the full difference matrix
        ! in the MO basis
        rows = density_qn%history_channel_rows()
        if (any(shape(rows) /= [2_ip, n_particle])) then
            write(stderr, *) "test_history_channel_rows_density_qn_mo failed: "// &
                "Incorrect dimensions."
            test_history_channel_rows_density_qn_mo = .false.
        else if (any( &
            rows /= reshape([([(k - 1) * n_mo**2 + 1, k * n_mo**2], k=1, n_particle)], &
                            [2_ip, n_particle]))) then
            write(stderr, *) "test_history_channel_rows_density_qn_mo failed: "// &
                "Incorrect rows."
            test_history_channel_rows_density_qn_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo)

    end function test_history_channel_rows_density_qn_mo

    logical(c_bool) function test_history_channel_rows_density_qn_oao() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in a history column for the OAO basis
        !
        use otr_density_qn, only: density_qn_oao_type
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_particle

        type(oao_type), allocatable, target :: oao
        type(density_qn_oao_type) :: density_qn
        integer(ip), allocatable :: rows(:, :)
        integer(ip) :: j

        ! assume tests pass
        test_history_channel_rows_density_qn_oao = .true.

        ! set up the OAO object with the dimensions the routine reads
        allocate(oao)
        oao%n_ao = n_ao
        oao%n_particle = n_particle
        density_qn = density_qn_oao_type(oao)

        ! call routine and determine if every channel holds a full difference matrix
        rows = density_qn%history_channel_rows()
        if (any(shape(rows) /= [2_ip, n_particle])) then
            write(stderr, *) "test_history_channel_rows_density_qn_oao failed: "// &
                "Incorrect dimensions."
            test_history_channel_rows_density_qn_oao = .false.
        else if (any(rows /= &
                     reshape([((j - 1) * n_ao**2 + 1, j * n_ao**2, j=1, n_particle)], &
                             [2_ip, n_particle]))) then
            write(stderr, *) "test_history_channel_rows_density_qn_oao failed: "// &
                "Incorrect rows."
            test_history_channel_rows_density_qn_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao)

    end function test_history_channel_rows_density_qn_oao

    logical(c_bool) function test_packed_channel_rows_density_qn_mo() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in a packed column for the MO basis
        !
        use otr_density_qn, only: density_qn_mo_type
        use otr_mo, only: mo_type
        use otr_mo_test_reference, only: case_irreps, case_names
        use otr_common_test_reference, only: n_occ
        use otr_mo_unit_tests, only: setup_minimal_mo_object

        type(mo_type), allocatable, target :: mo
        type(density_qn_mo_type) :: density_qn
        integer(ip), allocatable :: rows(:, :)
        integer(ip) :: n1, n2, i_sym

        ! assume tests pass
        test_packed_channel_rows_density_qn_mo = .true.

        ! set up the MO object with the differently occupied spin channels of the
        ! open-shell occupation case with symmetry, whose parameters are all the
        ! routine reads
        i_sym = findloc(case_names, "open-shell with symmetry", dim=1, kind=ip)
        call setup_minimal_mo_object(mo, n_occ, case_irreps(:, :, i_sym))
        density_qn = density_qn_mo_type(mo)

        ! call routine and determine if every channel holds its occupied-virtual
        ! rotations between orbitals of the same irrep, whose number differs between
        ! the spin channels
        n1 = count(mo%mo_channels(1)%param_mask, kind=ip)
        n2 = count(mo%mo_channels(2)%param_mask, kind=ip)
        rows = density_qn%packed_channel_rows()
        if (any(shape(rows) /= [2, 2])) then
            write(stderr, *) "test_packed_channel_rows_density_qn_mo failed: "// &
                "Incorrect dimensions."
            test_packed_channel_rows_density_qn_mo = .false.
        else if (any(rows /= reshape([1_ip, n1, &
                                      n1 + 1_ip, n1 + n2], [2, 2]))) then
            write(stderr, *) "test_packed_channel_rows_density_qn_mo failed: "// &
                "Incorrect rows."
            test_packed_channel_rows_density_qn_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo)

    end function test_packed_channel_rows_density_qn_mo

    logical(c_bool) function test_packed_channel_rows_density_qn_oao() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in a packed column for the OAO basis
        !
        use otr_density_qn, only: density_qn_oao_type
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_particle

        integer(ip), parameter :: n_channel = n_ao * (n_ao - 1) / 2

        type(oao_type), allocatable, target :: oao
        type(density_qn_oao_type) :: density_qn
        integer(ip), allocatable :: rows(:, :)
        integer(ip) :: j

        ! assume tests pass
        test_packed_channel_rows_density_qn_oao = .true.

        ! set up the OAO object with the dimensions the routine reads
        allocate(oao)
        oao%n_ao = n_ao
        oao%n_particle = n_particle
        density_qn = density_qn_oao_type(oao)

        ! call routine and determine if every channel holds the antisymmetric
        ! parameters of that channel
        rows = density_qn%packed_channel_rows()
        if (any(shape(rows) /= [2_ip, n_particle])) then
            write(stderr, *) "test_packed_channel_rows_density_qn_oao failed: "// &
                "Incorrect dimensions."
            test_packed_channel_rows_density_qn_oao = .false.
        else if (any(rows /= reshape( &
            [((j - 1) * n_channel + 1, j * n_channel, j=1, n_particle)], &
            [2_ip, n_particle]))) then
            write(stderr, *) "test_packed_channel_rows_density_qn_oao failed: "// &
                "Incorrect rows."
            test_packed_channel_rows_density_qn_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao)

    end function test_packed_channel_rows_density_qn_oao

    logical(c_bool) function test_hess_x_static_density_qn_mo() bind(C)
        !
        ! this function tests the function which applies the static part of the Hessian
        ! to a trial vector for the MO basis, which applies the one of the MO object
        !
        use otr_density_qn, only: density_qn_mo_type
        use otr_mo, only: mo_type
        use otr_common_test_reference, only: n_occ
        use otr_mo_unit_tests, only: setup_minimal_mo_object
        use opentrustregion_unit_tests, only: identity_matrix

        type(mo_type), allocatable, target :: mo
        type(density_qn_mo_type) :: density_qn
        real(rp), allocatable :: x(:)
        integer(ip) :: k

        ! assume tests pass
        test_hess_x_static_density_qn_mo = .true.

        ! set up the MO object with vanishing occupied-occupied and unit
        ! virtual-virtual blocks of the Fock matrix, so that the static part of the
        ! Hessian, the open-shell 2 (X F_vv - F_oo X), doubles the trial vector
        call setup_minimal_mo_object(mo, n_occ)
        do k = 1, size(n_occ, kind=ip)
            associate (channel => mo%mo_channels(k))
                channel%fock_oo = 0.0_rp
                channel%fock_vv = identity_matrix(channel%n_virt)
            end associate
        end do
        density_qn = density_qn_mo_type(mo)
        allocate(x(mo%n_param))
        call random_number(x)

        ! call routine and determine if the trial vector is doubled
        if (norm2(density_qn%hess_x_static(x) - 2.0_rp * x) > tol) then
            write(stderr, *) "test_hess_x_static_density_qn_mo failed: Incorrect "// &
                "static part."
            test_hess_x_static_density_qn_mo = .false.
        end if
        deallocate(mo)

    end function test_hess_x_static_density_qn_mo

    logical(c_bool) function test_hess_x_static_density_qn_oao() bind(C)
        !
        ! this function tests the function which applies the static part of the Hessian
        ! to a trial vector for the OAO basis, for the closed- and the open-shell case
        !
        use otr_density_qn, only: density_qn_oao_type
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, shell_names
        use otr_oao_unit_tests, only: ref_unpack_asymm, ref_hess_x_oao, &
                                      generate_random_nonredundant_vector

        type(oao_type), allocatable, target :: oao
        type(density_qn_oao_type) :: density_qn
        real(rp), allocatable :: dm_oao(:, :, :), fock_oo(:, :, :), fock_vv(:, :, :), &
                                 x(:), response(:, :, :)
        integer(ip) :: n_particle, n_param

        ! assume tests pass
        test_hess_x_static_density_qn_oao = .true.

        ! loop over the closed-shell and the open-shell case, whose static parts are
        ! scaled differently
        do n_particle = 1, n_particle_ref
            n_param = n_particle * n_ao * (n_ao - 1) / 2

            ! set up the OAO object with a random partition of the Fock matrix
            allocate(dm_oao(n_ao, n_ao, n_particle), fock_oo(n_ao, n_ao, n_particle), &
                     fock_vv(n_ao, n_ao, n_particle), response(n_ao, n_ao, n_particle))
            call generate_random_fock_partition(n_ao, n_occ(:n_particle), dm_oao, &
                                                fock_oo, fock_vv)
            allocate(oao)
            oao%n_ao = n_ao
            oao%n_particle = n_particle
            oao%dm_oao = dm_oao
            oao%fock_oo = fock_oo
            oao%fock_vv = fock_vv
            density_qn = density_qn_oao_type(oao)

            ! call routine and determine if the static part matches the independently
            ! assembled one
            x = generate_random_nonredundant_vector(n_param, n_particle, n_ao, dm_oao)
            response = 0.0_rp
            if (norm2(density_qn%hess_x_static(x) - &
                      ref_hess_x_oao(ref_unpack_asymm(x, n_particle, n_ao), response, &
                                     dm_oao, fock_oo, fock_vv, n_param)) > tol) then
                write(stderr, *) "test_hess_x_static_density_qn_oao failed: "// &
                    "Incorrect static part for the "//trim(shell_names(n_particle))// &
                    " case."
                test_hess_x_static_density_qn_oao = .false.
            end if

            ! deallocate OAO object and arrays of the shell
            deallocate(oao, dm_oao, fock_oo, fock_vv, response)
        end do

    end function test_hess_x_static_density_qn_oao

    logical(c_bool) function test_cache_channel_split_dirs() bind(C)
        !
        ! this function tests the open-shell history-projection caching routine which
        ! keeps the two spin channels of every packed column as separate columns,
        ! zeroing the rows of the respective other channel, then rebases the result via
        ! a given (map, chol) pair; undoing the rebasing by right-multiplying by chol
        ! again must reproduce the gathered, reordered raw split columns exactly
        !
        use otr_density_qn, only: cache_channel_split_dirs
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_list = 2, n_col = n_particle * n_list, &
                                  n_rows(n_particle) = [4, 3]

        real(rp) :: packed(sum(n_rows), n_list), raw(sum(n_rows), n_col), &
                    chol(n_col, n_col)
        real(rp), allocatable :: u(:, :)
        integer(ip) :: rows(2, n_particle), map(n_col), j, k

        ! assume test passes
        test_cache_channel_split_dirs = .true.

        ! generate random packed columns and the rows of both channels
        call random_number(packed)
        rows = reshape([1_ip, n_rows(1), &
                        n_rows(1) + 1_ip, sum(n_rows)], [2, 2])

        ! independently reproduce the raw (un-rebased) split columns
        raw = 0.0_rp
        do k = 1, n_list
            do j = 1, 2
                raw(rows(1, j):rows(2, j), (j - 1) * n_list + k) = &
                    packed(rows(1, j):rows(2, j), k)
            end do
        end do

        ! reorder and rebase into an arbitrary orthonormalized basis
        map = [3, 1, 4, 2]
        chol = generate_random_cholesky_factor(n_col)

        ! call routine and determine if dimensions are correct and undoing the rebasing
        ! reproduces the gathered, reordered raw split columns
        call cache_channel_split_dirs(packed, rows, map, chol, u)
        if (size(u, 1) /= sum(n_rows) .or. size(u, 2) /= n_col) then
            write(stderr, *) "test_cache_channel_split_dirs failed: Incorrect "// &
                "dimensions of split directions."
            test_cache_channel_split_dirs = .false.
        else if (norm2(matmul(u, chol) - raw(:, map)) > tol) then
            write(stderr, *) "test_cache_channel_split_dirs failed: Split "// &
                "directions are not rebased directions of the gathered, reordered "// &
                "raw split columns."
            test_cache_channel_split_dirs = .false.
        end if
        deallocate(u)

    end function test_cache_channel_split_dirs

    logical(c_bool) function test_cache_combined_channel_dirs() bind(C)
        !
        ! this function tests the open-shell history-projection caching routine which
        ! combines a same-spin channel of every packed column with the opposite-spin
        ! channel of the other spin, then rebases the result via a given (map, chol)
        ! pair; undoing the rebasing by right-multiplying by chol again must reproduce
        ! the gathered, reordered raw combined columns exactly
        !
        use otr_density_qn, only: cache_combined_channel_dirs
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_list = 2, n_col = n_particle * n_list, &
                                  n_rows(n_particle) = [4, 3]

        real(rp) :: same(sum(n_rows), n_list), opp(sum(n_rows), n_list), &
                    raw(sum(n_rows), n_col), chol(n_col, n_col)
        real(rp), allocatable :: u(:, :)
        integer(ip) :: rows(2, n_particle), map(n_col), k

        ! assume test passes
        test_cache_combined_channel_dirs = .true.

        ! generate random packed columns and the rows of both channels
        call random_number(same)
        call random_number(opp)
        rows = reshape([1_ip, n_rows(1), &
                        n_rows(1) + 1_ip, sum(n_rows)], [2, 2])

        ! independently reproduce the raw (un-rebased) combined columns
        do k = 1, n_list
            raw(:n_rows(1), k) = same(:n_rows(1), k)
            raw(n_rows(1) + 1:, k) = opp(n_rows(1) + 1:, k)
            raw(:n_rows(1), n_list + k) = opp(:n_rows(1), k)
            raw(n_rows(1) + 1:, n_list + k) = same(n_rows(1) + 1:, k)
        end do

        ! reorder and rebase into an arbitrary orthonormalized basis
        map = [3, 1, 4, 2]
        chol = generate_random_cholesky_factor(n_col)

        ! call routine and determine if dimensions are correct and undoing the rebasing
        ! reproduces the gathered, reordered raw combined columns
        call cache_combined_channel_dirs(same, opp, rows, map, chol, u)
        if (size(u, 1) /= sum(n_rows) .or. size(u, 2) /= n_col) then
            write(stderr, *) "test_cache_combined_channel_dirs failed: Incorrect "// &
                "dimensions of combined directions."
            test_cache_combined_channel_dirs = .false.
        else if (norm2(matmul(u, chol) - raw(:, map)) > tol) then
            write(stderr, *) "test_cache_combined_channel_dirs failed: Combined "// &
                "directions are not rebased directions of the gathered, reordered "// &
                "raw combined columns."
            test_cache_combined_channel_dirs = .false.
        end if
        deallocate(u)

    end function test_cache_combined_channel_dirs

    logical(c_bool) function test_get_low_rank_hess_factors() bind(C)
        !
        ! this function tests the subroutine which assembles the low-rank part of the
        ! approximate Hessian
        !
        use otr_density_qn, only: density_qn_type, density_qn_oao_type, qn_methods
        use otr_oao, only: oao_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix
        use otr_common_test_reference, only: n_particle_ref => n_particle, shell_names

        type(test_context_type), target :: context
        integer(ip), parameter :: n_param = 5, n_diff = 3, n_diff_nl = 2, &
                                  n_both = n_diff + n_diff_nl, dm_nl = 2 * n_diff, &
                                  pot_nl = 2 * n_diff + n_diff_nl

        type(oao_type), allocatable, target :: oao
        class(density_qn_type), allocatable, target :: density_qn
        real(rp) :: density(n_param, n_diff), density_nl(n_param, n_diff_nl), &
                    linear(n_param, n_diff), nonlinear(n_param, n_diff_nl), &
                    a_sym(n_diff, n_diff), a_sym_nl(n_diff_nl, n_diff_nl), &
                    a_inv(n_diff, n_diff), a_inv_comb(n_diff_nl, n_diff_nl), scale
        real(rp), allocatable :: expansion(:, :), projection(:, :), coupling(:, :)
        integer(ip) :: i_method, n_particle
        character(len=:), allocatable :: case_qn_method, case_name

        ! assume tests pass
        test_get_low_rank_hess_factors = .true.

        ! generate the packed history directions and the small dense matrices the
        ! coupling matrices are assembled from
        call random_number(density)
        call random_number(density_nl)
        call random_number(linear)
        call random_number(nonlinear)
        a_sym = generate_random_symm_matrix(n_diff)
        a_sym_nl = generate_random_symm_matrix(n_diff_nl)
        a_inv = generate_random_symm_matrix(n_diff)
        a_inv_comb = generate_random_symm_matrix(n_diff_nl)

        ! set up the density QN object with the cached quantities the assembly requires
        allocate(oao)
        oao%n_param = n_param
        allocate(density_qn, source=density_qn_oao_type(oao))
        call setup_settings(density_qn%settings, context)
        density_qn%dm_dirs = density
        density_qn%dm_dirs_nonlinear = density_nl
        density_qn%linear_potential_dirs = linear
        density_qn%nonlinear_potential_dirs = nonlinear
        density_qn%a_sym = a_sym
        density_qn%a_sym_nonlinear = a_sym_nl
        density_qn%a_inv = a_inv
        density_qn%a_inv_comb = a_inv_comb

        ! call routine for every QN method and both shells, whose coupling matrices
        ! differ by a factor of two, and for an empty history, and determine if the
        ! factors match
        do i_method = 1, size(qn_methods)
            case_qn_method = trim(qn_methods(i_method))
            do n_particle = 1, n_particle_ref
                case_name = case_qn_method//" in the "// &
                            trim(shell_names(n_particle))//" case"
                scale = merge(8.0_rp, 4.0_rp, n_particle == 1)

                ! the expected factors: every system (linear and non-linear) has its
                ! own block of directions, which no coupling mixes across systems
                select case (case_qn_method)
                ! multisecant SR1 couples the potential difference histories of each
                ! system through its own pseudoinverse
                case ("ms_sr1")
                    expansion = reshape([linear, nonlinear], [n_param, n_both])
                    projection = expansion
                    coupling = block_diagonal_matrix(scale * a_inv, scale * a_inv_comb)
                ! subspace-projected multisecant couples the density difference
                ! histories of each system through its symmetrized A matrix
                case ("ms_sp")
                    expansion = reshape([density, density_nl], [n_param, n_both])
                    projection = expansion
                    coupling = block_diagonal_matrix(scale * a_sym, scale * a_sym_nl)
                ! symmetrized ARH couples the density and potential difference
                ! histories of each system in both directions with half the identity,
                ! and multisecant PSB with the identity, adding a density-density block
                ! subtracting the doubly counted curvature
                case ("symm_arh", "ms_psb")
                    expansion = reshape([density, linear, density_nl, nonlinear], &
                                        [n_param, 2 * n_both])
                    projection = expansion
                    coupling = merge(0.5_rp, 1.0_rp, case_qn_method == "symm_arh") * &
                               scale * block_diagonal_matrix( &
                                   cshift(identity_matrix(2 * n_diff), n_diff), &
                                   cshift(identity_matrix(2 * n_diff_nl), n_diff_nl))
                    if (case_qn_method == "ms_psb") then
                        coupling(:n_diff, :n_diff) = -scale * a_sym
                        coupling(dm_nl + 1:pot_nl, dm_nl + 1:pot_nl) = -scale * a_sym_nl
                    end if
                ! standard ARH expands in the potential and contracts against the
                ! density difference history of each system with the identity
                case ("arh")
                    expansion = reshape([linear, nonlinear], [n_param, n_both])
                    projection = reshape([density, density_nl], [n_param, n_both])
                    coupling = scale * identity_matrix(n_both)
                end select

                ! call routine with directions rotated into the eigenbasis of the
                ! static part, and determine if the factors match and these are
                ! discarded
                oao%n_particle = n_particle
                density_qn%settings%qn_method = case_qn_method
                density_qn%rotated_expansion_dirs = density
                density_qn%rotated_projection_dirs = density
                call density_qn%get_low_rank_hess_factors()
                if (allocated(density_qn%rotated_expansion_dirs) .or. &
                    allocated(density_qn%rotated_projection_dirs)) then
                    write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                        "Rotated directions of the previous factors not discarded "// &
                        "for "//case_name//"."
                    test_get_low_rank_hess_factors = .false.
                end if
                if (.not. (allocated(density_qn%expansion_dirs) .and. &
                           allocated(density_qn%projection_dirs) .and. &
                           allocated(density_qn%coupling_matrix))) then
                    write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                        "Factors not assembled for "//case_name//"."
                    test_get_low_rank_hess_factors = .false.
                else if ( &
                    any(shape(density_qn%expansion_dirs) /= shape(expansion)) .or. &
                    any(shape(density_qn%projection_dirs) /= shape(projection)) .or. &
                    any(shape(density_qn%coupling_matrix) /= shape(coupling))) then
                    write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                        "Incorrect dimensions of the factors for "//case_name//"."
                    test_get_low_rank_hess_factors = .false.
                else
                    if (norm2(density_qn%expansion_dirs - expansion) > tol) then
                        write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                            "Incorrect expansion directions for "//case_name//"."
                        test_get_low_rank_hess_factors = .false.
                    end if
                    if (norm2(density_qn%projection_dirs - projection) > tol) then
                        write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                            "Incorrect projection directions for "//case_name//"."
                        test_get_low_rank_hess_factors = .false.
                    end if
                    if (norm2(density_qn%coupling_matrix - coupling) > tol) then
                        write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                            "Incorrect coupling matrix for "//case_name//"."
                        test_get_low_rank_hess_factors = .false.
                    end if
                end if
            end do

            ! an empty history has to leave every factor of the method unallocated
            density_qn%dm_dirs = density(:, :0)
            density_qn%dm_dirs_nonlinear = density_nl(:, :0)
            density_qn%linear_potential_dirs = linear(:, :0)
            density_qn%nonlinear_potential_dirs = nonlinear(:, :0)
            call density_qn%get_low_rank_hess_factors()
            if (allocated(density_qn%expansion_dirs) .or. &
                allocated(density_qn%projection_dirs) .or. &
                allocated(density_qn%coupling_matrix)) then
                write(stderr, *) "test_get_low_rank_hess_factors failed: Factors "// &
                    "were assembled for an empty history for "//case_qn_method//"."
                test_get_low_rank_hess_factors = .false.
            end if
            density_qn%dm_dirs = density
            density_qn%dm_dirs_nonlinear = density_nl
            density_qn%linear_potential_dirs = linear
            density_qn%nonlinear_potential_dirs = nonlinear
        end do

        ! deallocate density QN and OAO objects
        deallocate(density_qn, oao)

    contains

        function block_diagonal_matrix(a, b) result(matrix)
            !
            ! this function returns the block-diagonal matrix of the given blocks
            !
            real(rp), intent(in) :: a(:, :), b(:, :)
            real(rp) :: matrix(size(a, 1) + size(b, 1), size(a, 2) + size(b, 2))

            matrix = 0.0_rp
            matrix(:size(a, 1), :size(a, 2)) = a
            matrix(size(a, 1) + 1:, size(a, 2) + 1:) = b

        end function block_diagonal_matrix

    end function test_get_low_rank_hess_factors

    logical(c_bool) function test_build_a_part() bind(C)
        !
        ! this function tests the subroutine which constructs A = S^T Y for a single
        ! part of a coupled potential-difference response and symmetrizes it
        !
        use otr_density_qn, only: build_a_part

        integer(ip), parameter :: n_diff = 3, n_rows = 7
        real(rp) :: dm_cols(n_rows, n_diff), v_cols(n_rows, n_diff), &
                    a(n_diff, n_diff), expected(n_diff, n_diff), a_empty(0, 0)

        ! assume tests pass
        test_build_a_part = .true.

        ! random history and response so that the raw product is generically asymmetric
        ! and the symmetrization is exercised
        call random_number(dm_cols)
        call random_number(v_cols)

        ! get expected matrix
        expected = ref_build_a_part(dm_cols, v_cols)

        ! call routine and determine if the symmetrized matrix matches
        call build_a_part(dm_cols, v_cols, a)
        if (maxval(abs(a - expected)) > tol) then
            write(stderr, *) "test_build_a_part failed: Incorrect symmetrized A matrix."
            test_build_a_part = .false.
        end if

        ! call routine for an empty history, which has to return without handing BLAS
        ! the vanishing leading dimension of its result, which reference BLAS rejects
        call build_a_part(dm_cols(:, :0), v_cols(:, :0), a_empty)

    end function test_build_a_part

    logical(c_bool) function test_build_a_transformed() bind(C)
        !
        ! this function tests the function which builds the A = S^T Y matrix of a
        ! single part of the response and congruence-transforms it into the
        ! orthonormalized basis of that part: undoing the transform by left- and
        ! right-multiplying by chol again must reproduce the gathered, reordered raw
        ! symmetrized matrix exactly
        !
        use otr_density_qn, only: build_a_transformed

        integer(ip), parameter :: n_diff = 3, n_rows = 7

        real(rp) :: dm_cols(n_rows, n_diff), v_cols(n_rows, n_diff), &
                    expected(n_diff, n_diff), chol(n_diff, n_diff)
        real(rp), allocatable :: a_t(:, :)
        integer(ip) :: map(n_diff)

        ! assume tests pass
        test_build_a_transformed = .true.

        ! random density matrix and potential differences, so that the symmetrization
        ! acts on a generically asymmetric contribution
        call random_number(dm_cols)
        call random_number(v_cols)
        expected = ref_build_a_part(dm_cols, v_cols)

        ! reorder and congruence-transform into an arbitrary orthonormalized basis
        map = [2, 3, 1]
        chol = generate_random_cholesky_factor(n_diff)

        a_t = build_a_transformed(dm_cols, v_cols, map, chol)
        if (size(a_t, 1) /= n_diff .or. size(a_t, 2) /= n_diff) then
            write(stderr, *) "test_build_a_transformed failed: Incorrect dimensions."
            test_build_a_transformed = .false.
            return
        end if
        if (norm2(matmul(transpose(chol), matmul(a_t, chol)) - expected(map, map)) > &
            tol) then
            write(stderr, *) "test_build_a_transformed failed: Result does not "// &
                "invert back to the gathered, reordered raw symmetrized matrix."
            test_build_a_transformed = .false.
        end if
        deallocate(a_t)

    end function test_build_a_transformed

    logical(c_bool) function test_build_a_block_linear_os() bind(C)
        !
        ! this function tests the function which builds the cross-channel-symmetrized
        ! open-shell A = S^T Y matrix of the linear response and congruence-transforms
        ! it: undoing that transform must reproduce the gathered, reordered raw block
        ! matrix, whose diagonal blocks hold the symmetrized same-spin contributions
        ! and whose off-diagonal blocks hold the cross-symmetrized opposite-spin ones
        !
        use otr_density_qn, only: build_a_block_linear_os
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_diff = 3, n_col = n_particle * n_diff, &
                                  n_rows(n_particle) = [4, 3]

        real(rp) :: dm_cols(sum(n_rows), n_diff), v_same_cols(sum(n_rows), n_diff), &
                    v_opp_cols(sum(n_rows), n_diff), expected_a_block(n_col, n_col), &
                    chol(n_col, n_col), a_opp(n_diff, n_diff, 2), &
                    averaged(n_diff, n_diff)
        real(rp), allocatable :: a_block(:, :)
        integer(ip) :: rows(2, n_particle), map(n_col), j, lo, hi

        ! assume tests pass
        test_build_a_block_linear_os = .true.

        ! random per-channel histories and same-spin and opposite-spin potentials, so
        ! that every block is generically asymmetric
        call random_number(dm_cols)
        call random_number(v_same_cols)
        call random_number(v_opp_cols)
        rows = reshape([1_ip, n_rows(1), &
                        n_rows(1) + 1_ip, sum(n_rows)], [2, 2])

        ! each diagonal block holds the symmetrized same-spin contribution of one
        ! channel, while the raw opposite-spin products form the off-diagonal blocks
        do j = 1, 2
            lo = (j - 1) * n_diff + 1
            hi = j * n_diff
            expected_a_block(lo:hi, lo:hi) = ref_build_a_part( &
                dm_cols(rows(1, j):rows(2, j), :), &
                v_same_cols(rows(1, j):rows(2, j), :))
            a_opp(:, :, j) = matmul(transpose(dm_cols(rows(1, j):rows(2, j), :)), &
                                    v_opp_cols(rows(1, j):rows(2, j), :))
        end do

        ! the two opposite-spin blocks are exact transposes of one another, so any
        ! mismatch between them is noise and is averaged away rather than blended
        averaged = 0.5_rp * (a_opp(:, :, 1) + transpose(a_opp(:, :, 2)))
        expected_a_block(:n_diff, n_diff + 1:) = averaged
        expected_a_block(n_diff + 1:, :n_diff) = transpose(averaged)

        ! reorder and congruence-transform into an arbitrary orthonormalized basis
        map = [4, 1, 6, 2, 5, 3]
        chol = generate_random_cholesky_factor(n_col)

        ! generate A matrix and verify
        a_block = &
            build_a_block_linear_os(dm_cols, v_same_cols, v_opp_cols, rows, map, chol)
        if (size(a_block, 1) /= n_col .or. size(a_block, 2) /= n_col) then
            write(stderr, *) "test_build_a_block_linear_os failed: Incorrect "// &
                "dimensions."
            test_build_a_block_linear_os = .false.
        else if (norm2(matmul(transpose(chol), matmul(a_block, chol)) - &
                       expected_a_block(map, map)) > tol) then
            write(stderr, *) "test_build_a_block_linear_os failed: Result does not "// &
                "invert back to the gathered, reordered raw block matrix."
            test_build_a_block_linear_os = .false.
        end if

        ! generate A matrix for an empty history and determine if it is empty
        a_block = build_a_block_linear_os(dm_cols(:, :0), v_same_cols(:, :0), &
                                          v_opp_cols(:, :0), rows, [integer(ip) :: ], &
                                          reshape([real(rp) :: ], [0, 0]))
        if (size(a_block) /= 0) then
            write(stderr, *) "test_build_a_block_linear_os failed: A matrix not "// &
                "empty for empty history."
            test_build_a_block_linear_os = .false.
        end if

    end function test_build_a_block_linear_os

    logical(c_bool) function test_build_a_block_nonlinear_os() bind(C)
        !
        ! this function tests the function which builds the open-shell A = S^T Y matrix
        ! of the non-linear response and congruence-transforms it: the non-linear
        ! response has no opposite-spin counterpart, so the off-diagonal blocks must
        ! vanish and only the symmetrized per-channel diagonal blocks survive
        !
        use otr_density_qn, only: build_a_block_nonlinear_os
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_diff = 3, n_col = n_particle * n_diff, &
                                  n_rows(n_particle) = [4, 3]

        real(rp) :: dm_cols(sum(n_rows), n_diff), &
                    v_nonlinear_cols(sum(n_rows), n_diff), &
                    expected_a_block(n_col, n_col), chol(n_col, n_col)
        real(rp), allocatable :: a_block(:, :)
        integer(ip) :: rows(2, n_particle), map(n_col), j, lo, hi

        ! assume tests pass
        test_build_a_block_nonlinear_os = .true.

        ! random per-channel histories and non-linear potentials, so that every
        ! diagonal block is generically asymmetric
        call random_number(dm_cols)
        call random_number(v_nonlinear_cols)
        rows = reshape([1_ip, n_rows(1), &
                        n_rows(1) + 1_ip, sum(n_rows)], [2, 2])

        ! only the diagonal blocks are populated, each holding the symmetrized
        ! non-linear contribution of one channel
        expected_a_block = 0.0_rp
        do j = 1, 2
            lo = (j - 1) * n_diff + 1
            hi = j * n_diff
            expected_a_block(lo:hi, lo:hi) = ref_build_a_part( &
                dm_cols(rows(1, j):rows(2, j), :), &
                v_nonlinear_cols(rows(1, j):rows(2, j), :))
        end do

        ! reorder and congruence-transform into an arbitrary orthonormalized basis
        map = [4, 1, 6, 2, 5, 3]
        chol = generate_random_cholesky_factor(n_col)

        ! generate A matrix and verify
        a_block = build_a_block_nonlinear_os(dm_cols, v_nonlinear_cols, rows, map, chol)
        if (size(a_block, 1) /= n_col .or. size(a_block, 2) /= n_col) then
            write(stderr, *) "test_build_a_block_nonlinear_os failed: Incorrect "// &
                "dimensions."
            test_build_a_block_nonlinear_os = .false.
            return
        end if
        if (norm2(matmul(transpose(chol), matmul(a_block, chol)) - &
                  expected_a_block(map, map)) > tol) then
            write(stderr, *) "test_build_a_block_nonlinear_os failed: Result does "// &
                "not invert back to the gathered, reordered raw block matrix."
            test_build_a_block_nonlinear_os = .false.
        end if

    end function test_build_a_block_nonlinear_os

    logical(c_bool) function test_cross_symmetrize() bind(C)
        !
        ! this function tests the subroutine which cross-symmetrizes two related
        ! off-diagonal blocks of a larger matrix
        !
        use otr_density_qn, only: cross_symmetrize

        integer(ip), parameter :: n_rows = 2, n_cols = 3

        real(rp) :: a12(n_rows, n_cols), a21(n_cols, n_rows), &
                    expected12(n_rows, n_cols), expected21(n_cols, n_rows)

        ! assume tests pass
        test_cross_symmetrize = .true.

        ! initialize random rectangular blocks which are not transposes of each other
        call random_number(a12)
        call random_number(a21)

        ! initialize expected blocks, averaging a12(i, k) with a21(k, i) directly
        expected12 = 0.5_rp * (a12 + transpose(a21))
        expected21 = transpose(expected12)

        ! call routine and determine if values of resulting blocks match and if these
        ! are transposes of each other
        call cross_symmetrize(a12, a21)
        if (norm2(a12 - expected12) > tol) then
            write(stderr, *) "test_cross_symmetrize failed: Incorrect first block "// &
                "values after cross-symmetrization."
            test_cross_symmetrize = .false.
        end if
        if (norm2(a21 - expected21) > tol) then
            write(stderr, *) "test_cross_symmetrize failed: Incorrect second block "// &
                "values after cross-symmetrization."
            test_cross_symmetrize = .false.
        end if
        if (norm2(a21 - transpose(a12)) > tol) then
            write(stderr, *) "test_cross_symmetrize failed: Blocks are not "// &
                "transposes of each other after cross-symmetrization."
            test_cross_symmetrize = .false.
        end if

    end function test_cross_symmetrize

    logical(c_bool) function test_truncated_eigval_inv() bind(C)
        !
        ! this function tests the function which returns hard-truncated pseudoinverse
        ! eigenvalues
        !
        use otr_density_qn, only: truncated_eigval_inv

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
            write(stderr, *) "test_truncated_eigval_inv failed: Incorrect "// &
                "truncated inverse eigenvalues."
            test_truncated_eigval_inv = .false.
        end if

    end function test_truncated_eigval_inv

    logical(c_bool) function test_history_window() bind(C)
        !
        ! this function tests the function which returns the longest step the
        ! non-linear system of a density QN object may use, from the trust radius the
        ! record of the solver state holds
        !
        use opentrustregion, only: solver_state_type
        use otr_density_qn, only: density_qn_type, history_window_factor
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao
        use opentrustregion_unit_tests, only: test_context_type

        real(rp), parameter :: trust_radius = 0.3_rp, metric_ratio = 0.5_rp

        type(test_context_type), target :: context
        type(oao_type), allocatable, target :: oao
        class(density_qn_type), allocatable, target :: density_qn
        type(solver_state_type), target :: state
        real(rp), target :: dm_ao(n_ao, n_ao, 1)
        real(rp) :: dm_oao(n_ao, n_ao, 1)

        ! assume tests pass
        test_history_window = .true.

        ! set up the OAO and density QN objects
        call setup_density_qn_and_scaled_oao_objects(density_qn, oao, 1.0_rp, dm_ao, &
                                                     dm_oao, context)

        ! without a state, the whole history is used
        if (density_qn%history_window() < huge(1.0_rp)) then
            write(stderr, *) "test_history_window failed: The history is "// &
                "restricted without a state."
            test_history_window = .false.
        end if

        ! before the solver sets a trust radius, the whole history is used
        density_qn%state => state
        if (density_qn%history_window() < huge(1.0_rp)) then
            write(stderr, *) "test_history_window failed: The history is "// &
                "restricted before the solver set a trust radius."
            test_history_window = .false.
        end if

        ! the trust radius is converted into a distance between density matrices
        state%trust_radius = trust_radius
        if (abs(density_qn%history_window() - &
                history_window_factor * sqrt(2.0_rp) * trust_radius) > tol) then
            write(stderr, *) "test_history_window failed: Incorrect window for a "// &
                "spherical trust region."
            test_history_window = .false.
        end if

        ! and into a step length with the ratio of the last step to its size in the
        ! trust-region metric
        density_qn%step_norm_ratio = metric_ratio
        if (abs(density_qn%history_window() - history_window_factor * sqrt(2.0_rp) * &
                metric_ratio * trust_radius) > tol) then
            write(stderr, *) "test_history_window failed: Incorrect window for a "// &
                "trust region in another metric."
            test_history_window = .false.
        end if

        ! deallocate density QN and OAO objects
        deallocate(density_qn, oao)

    end function test_history_window

    logical(c_bool) function test_history_step_mask() bind(C)
        !
        ! this function tests the function which decides, from step length, which
        ! history entries the non-linear multisecant system is allowed to use
        !
        use otr_density_qn, only: history_step_mask

        integer(ip), parameter :: n_rows = 4, n_diff = 3, n_rows_os = 8, n_diff_os = 2
        real(rp), parameter :: max_step = 5.0_rp

        real(rp) :: dm_cols(n_rows, n_diff), v_cols(n_rows, n_diff), &
                    dm_cols_os(n_rows_os, n_diff_os), v_cols_os(n_rows_os, n_diff_os)
        logical, allocatable :: keep(:)

        ! assume tests pass
        test_history_step_mask = .true.

        ! initialize a history whose second step is exactly as long as the longest step
        ! the window admits and whose third is far longer, with random non-linear
        ! responses
        dm_cols = 0.0_rp
        dm_cols(1, 1) = 1.0_rp
        dm_cols(3, 2) = max_step
        dm_cols(2, 3) = 1e4_rp
        call random_number(v_cols)

        ! call routine and determine that only the entry beyond the window is dropped
        keep = history_step_mask(dm_cols, v_cols, max_step)
        if (size(keep) /= n_diff) then
            write(stderr, *) "test_history_step_mask failed: Incorrect number of "// &
                "entries."
            test_history_step_mask = .false.
            return
        end if
        if (.not. (keep(1) .and. keep(2))) then
            write(stderr, *) "test_history_step_mask failed: An entry within the "// &
                "window was not kept."
            test_history_step_mask = .false.
        end if
        if (keep(3)) then
            write(stderr, *) "test_history_step_mask failed: An entry beyond the "// &
                "window was not dropped."
            test_history_step_mask = .false.
        end if

        ! a window shorter than every step drops the entire history, which an
        ! unbounded window keeps
        keep = history_step_mask(dm_cols, v_cols, 0.5_rp)
        if (any(keep)) then
            write(stderr, *) "test_history_step_mask failed: An entry beyond a "// &
                "window shorter than every step was kept."
            test_history_step_mask = .false.
        end if
        keep = history_step_mask(dm_cols, v_cols, huge(1.0_rp))
        if (.not. all(keep)) then
            write(stderr, *) "test_history_step_mask failed: The history was not "// &
                "kept in full for an unbounded window."
            test_history_step_mask = .false.
        end if

        ! the open-shell step length spans both spin channels, so an entry that is long
        ! only in the rows of the second one still has to be dropped
        dm_cols_os = 0.0_rp
        dm_cols_os(1, 1) = 1.0_rp
        dm_cols_os(1, 2) = 1.0_rp
        dm_cols_os(6, 2) = 1e4_rp
        call random_number(v_cols_os)
        keep = history_step_mask(dm_cols_os, v_cols_os, max_step)
        if (.not. keep(1) .or. keep(2)) then
            write(stderr, *) "test_history_step_mask failed: An entry beyond the "// &
                "window only through its second spin channel was not dropped."
            test_history_step_mask = .false.
        end if

        ! without non-linear responses, every entry is dropped, also within the window,
        ! while a single non-vanishing response keeps the decision to the window
        v_cols = 0.0_rp
        keep = history_step_mask(dm_cols, v_cols, huge(1.0_rp))
        if (any(keep)) then
            write(stderr, *) "test_history_step_mask failed: An entry was kept "// &
                "without non-linear responses."
            test_history_step_mask = .false.
        end if
        v_cols(4, 3) = 1.0_rp
        keep = history_step_mask(dm_cols, v_cols, max_step)
        if (.not. (keep(1) .and. keep(2)) .or. keep(3)) then
            write(stderr, *) "test_history_step_mask failed: The window does not "// &
                "decide for a single non-vanishing non-linear response."
            test_history_step_mask = .false.
        end if

    end function test_history_step_mask

    logical(c_bool) function test_factorize_history() bind(C)
        !
        ! this function tests the subroutine which performs a pivoted, rank-revealing
        ! Cholesky factorization of the Gram matrix of a set of flattened history
        ! vectors
        !
        integer(ip), parameter :: n_rows = 4, n_diff = 3
        real(rp), parameter :: small_norm = 1e-8_rp, large_norm = 1e6_rp

        real(rp) :: cols(n_rows, n_diff)

        ! assume tests pass
        test_factorize_history = .true.

        ! the Gram matrix of a smaller and a larger multiple of the same direction and
        ! an orthogonal column normalizes to a unit diagonal with the first two columns
        ! fully overlapping: the magnitude difference no longer breaks the tie, so
        ! pivoting accepts the lower-indexed column first, then the orthogonal one,
        ! rejecting the parallel column, while the factor is reported for the unscaled
        ! columns
        cols(:, 1) = [0.5_rp, 0.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 2) = [2.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 3) = [0.0_rp, 1.0_rp, 1.0_rp, 0.0_rp]
        if (.not. check_factorization("a parallel pair of different length", cols, &
                                      [1_ip, 3_ip])) test_factorize_history = .false.

        ! an independent column many orders of magnitude shorter than the others is
        ! retained
        cols(:, 1) = [small_norm, 0.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 2) = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 3) = [0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp]
        if (.not. check_factorization("an independent much shorter column", cols, &
                                      [1_ip, 2_ip, 3_ip])) &
            test_factorize_history = .false.

        ! a vanishing column is rejected rather than admitted with a zero-norm factor
        ! column
        cols(:, 2) = 0.0_rp
        if (.not. check_factorization("a vanishing column", cols, [1_ip, 3_ip])) &
            test_factorize_history = .false.

        ! an excluded column is passed over in favour of its parallel twin, which an
        ! exact tie would otherwise resolve towards the lower index
        cols(:, 1) = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 2) = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 3) = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
        if (.not. check_factorization("an excluded column with a parallel twin", cols, &
                                      [2_ip, 3_ip], [.false., .true., .true.])) &
            test_factorize_history = .false.

        ! an excluded column is rejected even though it is linearly independent of the
        ! others
        cols(:, 1) = [3.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 2) = [0.0_rp, 0.0_rp, 0.0_rp, 7.0_rp]
        cols(:, 3) = [0.0_rp, 4.0_rp, 0.0_rp, 0.0_rp]
        if (.not. check_factorization("an excluded but independent column", cols, &
                                      [1_ip, 3_ip], [.true., .false., .true.])) &
            test_factorize_history = .false.

        ! whether a column far shorter than the others, but still far above any
        ! absolute round-off constant, has vanished has to be judged against the
        ! largest column rather than against an absolute constant
        cols(:, 1) = [large_norm, 0.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 2) = [0.0_rp, small_norm, 0.0_rp, 0.0_rp]
        cols(:, 3) = [0.0_rp, 0.0_rp, large_norm, 0.0_rp]
        if (.not. check_factorization( &
            "a column vanishing relative to the largest one", cols, [1_ip, 3_ip])) &
            test_factorize_history = .false.

        ! omitting the mask has to mean using the whole history, so that it agrees with
        ! a mask keeping every column
        cols(:, 1) = [2.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 2) = [0.0_rp, 3.0_rp, 0.0_rp, 0.0_rp]
        cols(:, 3) = [0.0_rp, 0.0_rp, 5.0_rp, 0.0_rp]
        if (.not. check_factorization("an omitted mask", cols, [1_ip, 2_ip, 3_ip])) &
            test_factorize_history = .false.
        if (.not. check_factorization("a mask keeping every column", cols, &
                                      [1_ip, 2_ip, 3_ip], spread(.true., 1, n_diff))) &
            test_factorize_history = .false.

        ! an empty history returns an empty factorization
        if (.not. check_factorization("an empty history", cols(:, :0), &
                                      [integer(ip) :: ])) &
            test_factorize_history = .false.

    contains

        function check_factorization(case_name, history_cols, expected_map, keep) &
            result(passed)
            !
            ! this function checks the pivoted Cholesky factorization of the Gram matrix
            ! of a set of history columns for one case: the accepted columns have to
            ! follow the expected map back to the original history indices and the
            ! factor has to reproduce the Gram matrix of the accepted columns, which
            ! pins it since it is upper triangular with a positive diagonal
            !
            use otr_density_qn, only: factorize_history

            character(len=*), intent(in) :: case_name
            real(rp), intent(in) :: history_cols(:, :)
            integer(ip), intent(in) :: expected_map(:)
            logical, intent(in), optional :: keep(:)
            logical :: passed

            real(rp), allocatable :: chol(:, :), gram(:, :)
            integer(ip), allocatable :: map(:)
            integer(ip) :: n_accepted

            ! assume test passes
            passed = .true.

            ! call routine and determine if the accepted columns and the factor match
            call factorize_history(history_cols, chol, map, n_accepted, keep)
            gram = matmul(transpose(history_cols), history_cols)
            if (n_accepted /= size(expected_map) .or. &
                size(map) /= size(expected_map) .or. &
                any(shape(chol) /= size(expected_map))) then
                write(stderr, *) "test_factorize_history failed: Incorrect number "// &
                    "of accepted columns for "//case_name//"."
                passed = .false.
            else if (any(map /= expected_map)) then
                write(stderr, *) "test_factorize_history failed: Incorrect map "// &
                    "back to original history indices for "//case_name//"."
                passed = .false.
            else if (norm2(matmul(transpose(chol), chol) - gram(map, map)) > &
                     tol * (1.0_rp + norm2(gram))) then
                write(stderr, *) "test_factorize_history failed: Cholesky factor "// &
                    "does not reproduce the Gram matrix of the accepted columns "// &
                    "for "//case_name//"."
                passed = .false.
            end if

        end function check_factorization

    end function test_factorize_history

    logical(c_bool) function test_rebase_dirs() bind(C)
        !
        ! this function tests the function which re-expresses a set of packed
        ! history-direction columns in the orthonormalized basis defined by a Cholesky
        ! factor: selecting, reordering and right-dividing by the factor should be
        ! exactly undone by right-multiplying by the factor again
        !
        use otr_density_qn, only: rebase_dirs

        integer(ip), parameter :: n_param = 3, n_dm = 3, n_accepted = 2

        real(rp) :: dirs(n_param, n_dm), chol(n_accepted, n_accepted)
        integer(ip) :: map(n_accepted)
        real(rp), allocatable :: rebased(:, :)

        ! assume tests pass
        test_rebase_dirs = .true.

        call random_number(dirs)
        chol = generate_random_cholesky_factor(n_accepted)
        map = [3, 1]

        rebased = rebase_dirs(dirs, map, chol)
        if (size(rebased, 1) /= n_param .or. size(rebased, 2) /= n_accepted) then
            write(stderr, *) "test_rebase_dirs failed: Incorrect dimensions."
            test_rebase_dirs = .false.
            return
        end if

        ! undoing the right-division by right-multiplying by the same factor must
        ! reproduce the gathered, reordered columns exactly
        if (norm2(matmul(rebased, chol) - dirs(:, map)) > tol) then
            write(stderr, *) "test_rebase_dirs failed: Rebased directions do not "// &
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
        use otr_density_qn, only: congruence_transform
        use opentrustregion_unit_tests, only: generate_random_symm_matrix

        integer(ip), parameter :: n_dm = 3, n_accepted = 2

        real(rp) :: a(n_dm, n_dm), chol(n_accepted, n_accepted), &
                    a_gathered(n_accepted, n_accepted)
        integer(ip) :: map(n_accepted)
        real(rp), allocatable :: a_tilde(:, :)

        ! assume tests pass
        test_congruence_transform = .true.

        a = generate_random_symm_matrix(n_dm)
        chol = generate_random_cholesky_factor(n_accepted)
        map = [3, 1]
        a_gathered = a(map, map)

        a_tilde = congruence_transform(a, map, chol)
        if (size(a_tilde, 1) /= n_accepted .or. size(a_tilde, 2) /= n_accepted) then
            write(stderr, *) "test_congruence_transform failed: Incorrect dimensions."
            test_congruence_transform = .false.
            return
        end if
        if (norm2(a_tilde - transpose(a_tilde)) > tol) then
            write(stderr, *) "test_congruence_transform failed: Result is not "// &
                "symmetric even though the input and the transform preserve symmetry."
            test_congruence_transform = .false.
        end if
        if (norm2(matmul(transpose(chol), matmul(a_tilde, chol)) - a_gathered) > tol) &
            then
            write(stderr, *) "test_congruence_transform failed: Result does not "// &
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
        use otr_density_qn, only: combine_channels

        integer(ip), parameter :: n1 = 2, n2 = 1, n_offset = 5

        real(rp) :: chol1(n1, n1), chol2(n2, n2)
        integer(ip) :: map1(n1), map2(n2)
        real(rp), allocatable :: chol_comb(:, :)
        integer(ip), allocatable :: map_comb(:)

        ! assume tests pass
        test_combine_channels = .true.

        ! initialize random upper triangular factors of both channels and maps into
        ! their histories
        chol1 = generate_random_cholesky_factor(n1)
        chol2 = generate_random_cholesky_factor(n2)
        map1 = [2, 4]
        map2 = [1]

        call combine_channels(chol1, map1, chol2, map2, n_offset, chol_comb, map_comb)
        if (size(chol_comb, 1) /= n1 + n2 .or. size(map_comb) /= n1 + n2) then
            write(stderr, *) "test_combine_channels failed: Incorrect dimensions."
            test_combine_channels = .false.
            return
        end if
        if (any(abs(chol_comb(1:n1, 1:n1) - chol1) > tol) .or. &
            any(abs(chol_comb(n1 + 1:, n1 + 1:) - chol2) > tol) .or. &
            any(abs(chol_comb(1:n1, n1 + 1:)) > tol) .or. &
            any(abs(chol_comb(n1 + 1:, 1:n1)) > tol)) then
            write(stderr, *) "test_combine_channels failed: Incorrect "// &
                "block-diagonal Cholesky factor."
            test_combine_channels = .false.
        end if
        if (any(map_comb(1:n1) /= map1) .or. &
            any(map_comb(n1 + 1:) /= n_offset + map2)) then
            write(stderr, *) "test_combine_channels failed: Incorrect combined map."
            test_combine_channels = .false.
        end if

    end function test_combine_channels

    logical(c_bool) function test_get_ms_a_inv() bind(C)
        !
        ! this function tests the subroutine which computes the pseudoinverse
        ! multisecant SR1 matrix of one part of the response
        !
        use otr_density_qn, only: get_ms_a_inv, density_qn_settings_type
        use otr_common_test_reference, only: n_ao
        use opentrustregion_unit_tests, only: test_context_type, setup_settings

        type(test_context_type), target :: context
        integer(ip), parameter :: n_diff = 3, n_accepted = 2, flat_len = n_ao * n_ao

        real(rp) :: chol(n_accepted, n_accepted), dm_cols(flat_len, n_diff), &
                    v_cols(flat_len, n_diff)
        real(rp), allocatable :: a_inv(:, :), expected(:, :), a_tilde(:, :), &
                                 y_gram(:, :)
        integer(ip) :: map(n_accepted), error
        type(density_qn_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! reorder and rebase into an arbitrary orthonormalized basis that also drops
        ! one history entry
        map = [3, 1]
        chol = generate_random_cholesky_factor(n_accepted)

        ! initialize random history and response; the entries are centered on zero so
        ! that a history direction can end up poorly aligned with its own response and
        ! be screened out
        call random_number(dm_cols)
        call random_number(v_cols)
        dm_cols = dm_cols - 0.5_rp
        v_cols = v_cols - 0.5_rp

        ! the expected pseudoinverse is assembled independently of the routine
        a_tilde = ref_congruence_transform(ref_build_a_part(dm_cols, v_cols), map, chol)
        y_gram = ref_congruence_transform(matmul(transpose(v_cols), v_cols), map, chol)
        expected = ref_ms_a_inv(a_tilde, y_gram)

        ! call routine and determine if the pseudoinverse matches
        call get_ms_a_inv(dm_cols, v_cols, map, chol, a_inv, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_a_inv failed: Produced error."
            test_get_ms_a_inv = .false.
        else if (norm2(a_inv - expected) > tol) then
            write(stderr, *) "test_get_ms_a_inv failed: Incorrect pseudoinverse."
            test_get_ms_a_inv = .false.
        end if

        ! a history whose responses live almost entirely outside the span of the steps
        ! leaves every direction badly aligned with its own response, so the screening
        ! criterion has to discard all of them and the pseudoinverse has to vanish
        call random_number(dm_cols)
        dm_cols(n_ao + 1:, :) = 0.0_rp
        call random_number(v_cols)
        v_cols(:n_ao, :) = 1e-3_rp * v_cols(:n_ao, :)

        ! call routine and determine if the pseudoinverse vanishes
        call get_ms_a_inv(dm_cols, v_cols, map, chol, a_inv, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_a_inv failed: Produced error for a badly "// &
                "aligned history."
            test_get_ms_a_inv = .false.
        else if (norm2(a_inv) > tol) then
            write(stderr, *) "test_get_ms_a_inv failed: Screening did not discard "// &
                "every direction of a badly aligned history."
            test_get_ms_a_inv = .false.
        end if

        ! call routine for an empty history and determine if dimensions of the
        ! resulting pseudoinverse vanish
        call get_ms_a_inv(dm_cols(:, :0), v_cols(:, :0), [integer(ip) :: ], &
                          reshape([real(rp) :: ], [0, 0]), a_inv, settings, error)
        if (size(a_inv, 1) /= 0 .or. size(a_inv, 2) /= 0) then
            write(stderr, *) "test_get_ms_a_inv failed: Incorrect pseudoinverse "// &
                "dimensions for empty history."
            test_get_ms_a_inv = .false.
        end if

    end function test_get_ms_a_inv

    logical(c_bool) function test_get_ms_a_inv_jk_cs() bind(C)
        !
        ! this function tests the subroutine which computes the closed-shell
        ! multisecant SR1 term of the linear potential as the difference of a Coulomb
        ! and an exact-exchange term: for a history in an orthonormalized basis that is
        ! known exactly (S = Q R) and positive semidefinite kernels, the model has to
        ! fulfill the secant condition of the linear kernel on the history, which only
        ! holds if both terms are assembled from the responses of their own kernel
        !
        use otr_density_qn, only: get_ms_a_inv_jk_cs, density_qn_settings_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix

        type(test_context_type), target :: context
        integer(ip), parameter :: n_rows = 6, n_diff = 3

        real(rp) :: coulomb(n_rows, n_rows), exchange(n_rows, n_rows), &
                    dm_cols(n_rows, n_diff), chol(n_diff, n_diff), &
                    v_coulomb_cols(n_rows, n_diff), v_exchange_cols(n_rows, n_diff)
        real(rp), allocatable :: a_inv(:, :), dirs(:, :), orth(:, :)
        integer(ip) :: map(n_diff), i, error
        type(density_qn_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_jk_cs = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! history S = Q R in identity order, so that R is the Cholesky factor of its
        ! metric, positive semidefinite kernels of a lower rank than the history space
        ! and their responses, which also serve as packed columns
        orth = generate_random_orthogonal_matrix(n_rows)
        chol = generate_random_cholesky_factor(n_diff)
        dm_cols = matmul(orth(:, :n_diff), chol)
        map = [(i, i=1, n_diff)]
        coulomb = generate_random_psd_matrix(n_rows, 4_ip)
        exchange = generate_random_psd_matrix(n_rows, 2_ip)
        v_coulomb_cols = matmul(coulomb, dm_cols)
        v_exchange_cols = matmul(exchange, dm_cols)

        ! call routine and determine if it fulfills the secant condition and the bounds
        ! of both terms
        call get_ms_a_inv_jk_cs(dm_cols, v_coulomb_cols, v_exchange_cols, &
                                v_coulomb_cols, v_exchange_cols, map, chol, a_inv, &
                                dirs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_a_inv_jk_cs failed: Produced error."
            test_get_ms_a_inv_jk_cs = .false.
        else if (size(a_inv, 1) /= 2 * n_diff .or. size(dirs, 2) /= 2 * n_diff) then
            write(stderr, *) "test_get_ms_a_inv_jk_cs failed: Incorrect dimensions."
            test_get_ms_a_inv_jk_cs = .false.
        else if (norm2(matmul(matmul(dirs, matmul(a_inv, transpose(dirs))), dm_cols) - &
                       (v_coulomb_cols - v_exchange_cols)) > tol) then
            write(stderr, *) "test_get_ms_a_inv_jk_cs failed: Secant condition not "// &
                "fulfilled."
            test_get_ms_a_inv_jk_cs = .false.
        end if

        ! call routine without exact exchange and determine if the exact-exchange term
        ! is left out and the Coulomb term alone fulfills the secant condition
        v_exchange_cols = 0.0_rp
        call get_ms_a_inv_jk_cs(dm_cols, v_coulomb_cols, v_exchange_cols, &
                                v_coulomb_cols, v_exchange_cols, map, chol, a_inv, &
                                dirs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_a_inv_jk_cs failed: Produced error "// &
                "without exact exchange."
            test_get_ms_a_inv_jk_cs = .false.
        else if (size(a_inv, 1) /= n_diff .or. size(dirs, 2) /= n_diff) then
            write(stderr, *) "test_get_ms_a_inv_jk_cs failed: Exact-exchange term "// &
                "not left out without exact exchange."
            test_get_ms_a_inv_jk_cs = .false.
        else if (norm2(matmul(matmul(dirs, matmul(a_inv, transpose(dirs))), dm_cols) - &
                       v_coulomb_cols) > tol) then
            write(stderr, *) "test_get_ms_a_inv_jk_cs failed: Secant condition not "// &
                "fulfilled without exact exchange."
            test_get_ms_a_inv_jk_cs = .false.
        end if

    end function test_get_ms_a_inv_jk_cs

    logical(c_bool) function test_get_ms_a_inv_jk_os() bind(C)
        !
        ! this function tests the subroutine which computes the open-shell multisecant
        ! SR1 term of the linear potential as the difference of a Coulomb and an
        ! exact-exchange term on the spin-separated density difference directions: for
        ! per-channel histories in orthonormalized bases that are known exactly, a
        ! Coulomb kernel acting through the total density on both channels and an
        ! exact-exchange kernel acting within every channel, the model has to fulfill
        ! the secant condition of the linear kernel on every spin-separated direction,
        ! which only holds if both terms are assembled from the responses of their own
        ! kernel and channels
        !
        use otr_density_qn, only: get_ms_a_inv_jk_os, density_qn_settings_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix

        type(test_context_type), target :: context
        integer(ip), parameter :: n_rows(2) = [4, 3], n_diff = 2, n_common = 3, &
                                  n_flat = sum(n_rows)

        real(rp) :: transform(n_common, n_flat), coulomb(n_flat, n_flat), &
                    exchange(n_flat, n_flat), dm_cols(n_flat, n_diff), &
                    directions(n_flat, 2 * n_diff), chol(2 * n_diff, 2 * n_diff), &
                    v_coulomb_alpha_cols(n_flat, n_diff), &
                    v_coulomb_beta_cols(n_flat, n_diff), v_exchange_cols(n_flat, n_diff)
        real(rp), allocatable :: a_inv(:, :), dirs(:, :), orth(:, :)
        integer(ip) :: rows(2, 2), map(2 * n_diff), i, error
        type(density_qn_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_jk_os = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! rows of the two channels, which the packed columns share
        rows = reshape([1_ip, n_rows(1), &
                        n_rows(1) + 1_ip, n_flat], [2, 2])

        ! per-channel histories S = Q R in identity order, whose factors combine
        ! block-diagonally as the per-channel factorizations do
        chol = 0.0_rp
        dm_cols = 0.0_rp
        orth = generate_random_orthogonal_matrix(n_rows(1))
        chol(:n_diff, :n_diff) = generate_random_cholesky_factor(n_diff)
        dm_cols(rows(1, 1):rows(2, 1), :) = matmul(orth(:, :n_diff), &
                                                   chol(:n_diff, :n_diff))
        orth = generate_random_orthogonal_matrix(n_rows(2))
        chol(n_diff + 1:, n_diff + 1:) = generate_random_cholesky_factor(n_diff)
        dm_cols(rows(1, 2):rows(2, 2), :) = matmul(orth(:, :n_diff), &
                                                   chol(n_diff + 1:, n_diff + 1:))
        map = [(i, i=1, 2 * n_diff)]

        ! spin-separated directions, the alpha directions first
        directions = 0.0_rp
        directions(rows(1, 1):rows(2, 1), :n_diff) = dm_cols(rows(1, 1):rows(2, 1), :)
        directions(rows(1, 2):rows(2, 2), n_diff + 1:) = &
            dm_cols(rows(1, 2):rows(2, 2), :)

        ! positive semidefinite kernels: the Coulomb kernel maps the density of both
        ! channels to a total density, on which a positive semidefinite kernel acts, and
        ! back to both channels, the exact-exchange kernel is block-diagonal
        call random_number(transform)
        coulomb = matmul(transpose(transform), matmul( &
            generate_random_psd_matrix(n_common, n_common), transform))
        exchange = 0.0_rp
        exchange(rows(1, 1):rows(2, 1), rows(1, 1):rows(2, 1)) = &
            generate_random_psd_matrix(n_rows(1), 2_ip)
        exchange(rows(1, 2):rows(2, 2), rows(1, 2):rows(2, 2)) = &
            generate_random_psd_matrix(n_rows(2), 2_ip)

        ! Coulomb responses of the alpha and beta directions in both channels and
        ! exact-exchange responses of every channel's direction in that channel, which
        ! also serve as packed columns
        v_coulomb_alpha_cols = matmul(coulomb, directions(:, :n_diff))
        v_coulomb_beta_cols = matmul(coulomb, directions(:, n_diff + 1:))
        v_exchange_cols = matmul(exchange, directions(:, :n_diff)) + &
                          matmul(exchange, directions(:, n_diff + 1:))

        ! call routine and determine if it fulfills the secant condition and the bounds
        ! of both terms
        call get_ms_a_inv_jk_os(dm_cols, v_coulomb_alpha_cols, v_coulomb_beta_cols, &
                                v_exchange_cols, v_coulomb_alpha_cols, &
                                v_coulomb_beta_cols, v_exchange_cols, rows, rows, map, &
                                chol, a_inv, dirs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_a_inv_jk_os failed: Produced error."
            test_get_ms_a_inv_jk_os = .false.
        else if (size(a_inv, 1) /= 4 * n_diff .or. size(dirs, 2) /= 4 * n_diff) then
            write(stderr, *) "test_get_ms_a_inv_jk_os failed: Incorrect dimensions."
            test_get_ms_a_inv_jk_os = .false.
        else if ( &
            norm2(matmul(matmul(dirs, matmul(a_inv, transpose(dirs))), directions) - &
                  matmul(coulomb - exchange, directions)) > tol) then
            write(stderr, *) "test_get_ms_a_inv_jk_os failed: Secant condition not "// &
                "fulfilled."
            test_get_ms_a_inv_jk_os = .false.
        end if

        ! call routine without exact exchange and determine if the exact-exchange term
        ! is left out and the Coulomb term alone fulfills the secant condition
        v_exchange_cols = 0.0_rp
        call get_ms_a_inv_jk_os(dm_cols, v_coulomb_alpha_cols, v_coulomb_beta_cols, &
                                v_exchange_cols, v_coulomb_alpha_cols, &
                                v_coulomb_beta_cols, v_exchange_cols, rows, rows, map, &
                                chol, a_inv, dirs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_a_inv_jk_os failed: Produced error "// &
                "without exact exchange."
            test_get_ms_a_inv_jk_os = .false.
        else if (size(a_inv, 1) /= 2 * n_diff .or. size(dirs, 2) /= 2 * n_diff) then
            write(stderr, *) "test_get_ms_a_inv_jk_os failed: Exact-exchange term "// &
                "not left out without exact exchange."
            test_get_ms_a_inv_jk_os = .false.
        else if (norm2(matmul(matmul(dirs, matmul(a_inv, transpose(dirs))), &
                              directions) - matmul(coulomb, directions)) > tol) then
            write(stderr, *) "test_get_ms_a_inv_jk_os failed: Secant condition not "// &
                "fulfilled without exact exchange."
            test_get_ms_a_inv_jk_os = .false.
        end if

    end function test_get_ms_a_inv_jk_os

    logical(c_bool) function test_get_ms_jk_inv() bind(C)
        !
        ! this function tests the subroutine which assembles the coupling and the
        ! directions of the multisecant SR1 term of the linear potential from a Coulomb
        ! and an exact-exchange A: both pseudoinverses with a joint noise floor, the
        ! exact-exchange term left out without exact exchange, an exact-exchange A of
        ! noise against the Coulomb A and negative eigenvalues not inverted, and an
        ! empty history
        !
        use otr_density_qn, only: get_ms_jk_inv, density_qn_settings_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings

        type(test_context_type), target :: context
        integer(ip), parameter :: n_diff = 3, n_accepted = 2, n_param = 5
        real(rp), parameter :: noise = 1e-20_rp

        real(rp) :: a_coulomb(n_diff, n_diff), a_exchange(n_diff, n_diff), &
                    chol(n_accepted, n_accepted), dirs_coulomb(n_param, n_accepted), &
                    dirs_exchange(n_param, n_accepted), scale
        real(rp), allocatable :: a_inv(:, :), dirs(:, :), expected(:, :), &
                                 a_tilde_coulomb(:, :), a_tilde_exchange(:, :)
        integer(ip) :: map(n_accepted), error
        type(density_qn_settings_type) :: settings

        ! assume tests pass
        test_get_ms_jk_inv = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! reorder and rebase into an arbitrary orthonormalized basis that also drops
        ! one history entry, and random directions of both terms
        map = [3, 1]
        chol = generate_random_cholesky_factor(n_accepted)
        a_coulomb = generate_random_psd_matrix(n_diff, n_diff)
        a_exchange = generate_random_psd_matrix(n_diff, n_diff)
        call random_number(dirs_coulomb)
        call random_number(dirs_exchange)

        ! the expected coupling is assembled independently of the routine, with the
        ! noise floor relative to the larger of both spectra
        a_tilde_coulomb = ref_congruence_transform(a_coulomb, map, chol)
        a_tilde_exchange = ref_congruence_transform(a_exchange, map, chol)
        scale = max(spectral_scale(a_tilde_coulomb), spectral_scale(a_tilde_exchange))
        allocate(expected(2 * n_accepted, 2 * n_accepted))
        expected = 0.0_rp
        expected(:n_accepted, :n_accepted) = ref_psd_pinv(a_tilde_coulomb, scale)
        expected(n_accepted + 1:, n_accepted + 1:) = &
            -ref_psd_pinv(a_tilde_exchange, scale)

        ! call routine and determine if the coupling and directions of both terms match
        call get_ms_jk_inv(a_coulomb, a_exchange, dirs_coulomb, dirs_exchange, &
                           .false., map, chol, a_inv, dirs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_jk_inv failed: Produced error."
            test_get_ms_jk_inv = .false.
        else if (norm2(a_inv - expected) > tol * (1.0_rp + norm2(expected)) .or. &
                 norm2(dirs(:, :n_accepted) - dirs_coulomb) > tol .or. &
                 norm2(dirs(:, n_accepted + 1:) - dirs_exchange) > tol) then
            write(stderr, *) "test_get_ms_jk_inv failed: Incorrect coupling or "// &
                "directions."
            test_get_ms_jk_inv = .false.
        end if

        ! call routine without exact exchange and determine if only the Coulomb term is
        ! returned
        call get_ms_jk_inv(a_coulomb, 0.0_rp * a_exchange, dirs_coulomb, &
                           dirs_exchange, .true., map, chol, a_inv, dirs, settings, &
                           error)
        expected = ref_psd_pinv(a_tilde_coulomb, spectral_scale(a_tilde_coulomb))
        if (error /= 0) then
            write(stderr, *) "test_get_ms_jk_inv failed: Produced error without "// &
                "exact exchange."
            test_get_ms_jk_inv = .false.
        else if (size(a_inv, 1) /= n_accepted .or. size(dirs, 2) /= n_accepted) then
            write(stderr, *) "test_get_ms_jk_inv failed: Exact-exchange term not "// &
                "left out without exact exchange."
            test_get_ms_jk_inv = .false.
        else if (norm2(a_inv - expected) > tol * (1.0_rp + norm2(expected)) .or. &
                 norm2(dirs - dirs_coulomb) > tol) then
            write(stderr, *) "test_get_ms_jk_inv failed: Incorrect coupling or "// &
                "directions without exact exchange."
            test_get_ms_jk_inv = .false.
        end if

        ! call routine for an exact-exchange A at the level of rounding noise against
        ! the Coulomb A and determine if it is not inverted, which a noise floor
        ! relative to its own spectrum would do
        call get_ms_jk_inv(a_coulomb, noise * a_exchange, dirs_coulomb, dirs_exchange, &
                           .false., map, chol, a_inv, dirs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_jk_inv failed: Produced error for an "// &
                "exact-exchange A of noise."
            test_get_ms_jk_inv = .false.
        else if (norm2(a_inv(n_accepted + 1:, n_accepted + 1:)) > tol) then
            write(stderr, *) "test_get_ms_jk_inv failed: Exact-exchange A of noise "// &
                "inverted."
            test_get_ms_jk_inv = .false.
        end if

        ! call routine for a negative definite exact-exchange A and determine if its
        ! eigenvalues, which a positive semidefinite matrix only has below the noise
        ! floor, are not inverted however large they are
        call get_ms_jk_inv(a_coulomb, -a_exchange, dirs_coulomb, dirs_exchange, &
                           .false., map, chol, a_inv, dirs, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_jk_inv failed: Produced error for a "// &
                "negative definite exact-exchange A."
            test_get_ms_jk_inv = .false.
        else if (norm2(a_inv(n_accepted + 1:, n_accepted + 1:)) > tol) then
            write(stderr, *) "test_get_ms_jk_inv failed: Negative eigenvalues of "// &
                "exact-exchange A inverted."
            test_get_ms_jk_inv = .false.
        end if

        ! call routine for an empty history and determine if the coupling and the
        ! directions are empty
        call get_ms_jk_inv(a_coulomb(:0, :0), a_exchange(:0, :0), dirs_coulomb(:, :0), &
                           dirs_exchange(:, :0), .false., [integer(ip) :: ], reshape( &
                               [real(rp) :: ], [0, 0]), a_inv, dirs, settings, error)
        if (error /= 0 .or. size(a_inv) /= 0 .or. size(dirs, 1) /= n_param .or. &
            size(dirs, 2) /= 0) then
            write(stderr, *) "test_get_ms_jk_inv failed: Incorrect coupling or "// &
                "directions for empty history."
            test_get_ms_jk_inv = .false.
        end if

    contains

        function spectral_scale(a) result(scale_out)
            !
            ! this function returns the largest absolute eigenvalue of a symmetric
            ! matrix, the scale the noise floor is relative to
            !
            real(rp), intent(in) :: a(:, :)
            real(rp) :: scale_out

            integer(ip) :: n, info
            real(rp), allocatable :: vecs(:, :), vals(:), work(:)
            external :: dsyev

            n = size(a, 1)
            allocate(vecs(n, n), vals(n), work(3 * n))
            vecs = 0.5_rp * (a + transpose(a))
            call dsyev("N", "U", n, vecs, n, vals, work, 3_ip * n, info)
            scale_out = maxval(abs(vals))

        end function spectral_scale

    end function test_get_ms_jk_inv

    logical(c_bool) function test_response_gram() bind(C)
        !
        ! this function tests the function which returns the Gram matrix of the
        ! response history rebased into the orthonormalized S-basis
        !
        use otr_density_qn, only: response_gram

        integer(ip), parameter :: n_dm = 3, n_accepted = 2, n_rows = 7
        integer(ip) :: map(n_accepted)
        real(rp) :: v_cols(n_rows, n_dm), chol(n_accepted, n_accepted), gram(n_dm, n_dm)
        real(rp), allocatable :: y_gram(:, :), expected(:, :)

        ! assume tests pass
        test_response_gram = .true.

        ! random response history and an upper-triangular Cholesky factor with a
        ! non-trivial map that both selects a subset and reorders it
        call random_number(v_cols)
        map = [3_ip, 1_ip]
        chol = generate_random_cholesky_factor(n_accepted)

        ! build the expected Gram matrix independently from the history columns
        gram = matmul(transpose(v_cols), v_cols)
        expected = ref_congruence_transform(gram, map, chol)

        ! call routine and determine if the rebased Gram matrix matches
        y_gram = response_gram(v_cols, map, chol)
        if (size(y_gram, 1) /= n_accepted .or. size(y_gram, 2) /= n_accepted) then
            write(stderr, *) "test_response_gram failed: Incorrect shape."
            test_response_gram = .false.
        else if (maxval(abs(y_gram - expected)) > tol) then
            write(stderr, *) "test_response_gram failed: Incorrect rebased "// &
                "response Gram matrix."
            test_response_gram = .false.
        end if

    end function test_response_gram

    logical(c_bool) function test_apply_ms_sr1_skip() bind(C)
        !
        ! this function tests the subroutine which discards history directions failing
        ! the multisecant SR1 skipping criterion
        !
        use otr_density_qn, only: apply_ms_sr1_skip, ms_sr1_skip_thresh
        use opentrustregion_unit_tests, only: identity_matrix

        integer(ip), parameter :: n = 4
        integer(ip) :: i
        real(rp) :: eig_vals(n), eig_vecs(n, n), y_gram(n, n), eig_vals_inv(n), &
                    expected(n), y_norm(n), factor(n), b(n, n), theta, c, s

        ! assume tests pass
        test_apply_ms_sr1_skip = .true.

        ! each eigenvalue is placed at a fixed multiple of the criterion applied to its
        ! own response norm, so which directions are discarded does not depend on the
        ! value of the skipping threshold; the multiples sit just either side of one in
        ! magnitude, with either sign, so that a response norm computed even slightly
        ! wrongly flips a decision and only the magnitude of a negative eigenvalue
        ! decides whether it is kept
        factor = [-0.95_rp, 0.95_rp, 1.05_rp, -1.05_rp]

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
            y_norm(i) = sqrt(max( &
                dot_product(eig_vecs(:, i), matmul(y_gram, eig_vecs(:, i))), 0.0_rp))
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
            write(stderr, *) "test_apply_ms_sr1_skip failed: Incorrect inverse "// &
                "eigenvalues after skipping."
            test_apply_ms_sr1_skip = .false.
        end if

    end function test_apply_ms_sr1_skip

    logical(c_bool) function test_spectral_to_dense() bind(C)
        !
        ! this function tests the function which reconstructs a dense symmetric matrix
        ! from its eigenvectors and inverted eigenvalues
        !
        use otr_density_qn, only: spectral_to_dense
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix

        integer(ip), parameter :: n = 3

        real(rp) :: eigvecs(n, n), inv_eigvals(n), expected(n, n)
        real(rp), allocatable :: mat(:, :)
        integer(ip) :: i

        ! assume tests pass
        test_spectral_to_dense = .true.

        ! initialize a random orthonormal eigenvector matrix, which is not symmetric, so
        ! that a transposed reconstruction would be caught, and distinct random
        ! inverted eigenvalues
        eigvecs = generate_random_orthogonal_matrix(n)
        call random_number(inv_eigvals)

        ! initialize the expected matrix as the sum of the outer products of the
        ! eigenvectors weighted by the inverted eigenvalues
        expected = 0.0_rp
        do i = 1, n
            expected = expected + inv_eigvals(i) * spread(eigvecs(:, i), 2, n) * &
                       spread(eigvecs(:, i), 1, n)
        end do

        ! call routine and determine if dimensions and values of the reconstructed
        ! matrix match
        mat = spectral_to_dense(eigvecs, inv_eigvals)
        if (size(mat, 1) /= n .or. size(mat, 2) /= n) then
            write(stderr, *) "test_spectral_to_dense failed: Incorrect dimensions "// &
                "of reconstructed matrix."
            test_spectral_to_dense = .false.
        else if (norm2(mat - expected) > tol) then
            write(stderr, *) "test_spectral_to_dense failed: Incorrect "// &
                "reconstructed matrix."
            test_spectral_to_dense = .false.
        end if

    end function test_spectral_to_dense

    logical(c_bool) function test_density_in_history() bind(C)
        !
        ! this function tests the function which reports whether a density is already
        ! held in the history
        !
        use opentrustregion, only: numerical_zero
        use otr_density_qn, only: density_qn_type, density_qn_oao_type
        use otr_common_test_reference, only: n_ao, n_particle

        real(rp), parameter :: scale = 1e6_rp

        class(density_qn_type), allocatable :: density_qn
        real(rp) :: dm(n_ao, n_ao, n_particle), other(n_ao, n_ao, n_particle)

        ! assume tests pass
        test_density_in_history = .true.

        ! generate two different densities
        call random_number(dm)
        call random_number(other)
        other = other + 1.0_rp

        ! determine if no density is reported before the history exists
        allocate(density_qn_oao_type :: density_qn)
        if (density_qn%density_in_history(dm)) then
            write(stderr, *) "test_density_in_history failed: Density matrix "// &
                "reported without a history."
            test_density_in_history = .false.
        end if

        ! determine if a density is found in a history holding it besides another one,
        ! also when it differs by much less than numerical zero relative to its size,
        ! but not when it differs by much more
        allocate(density_qn%dm_list(n_ao, n_ao, n_particle, 2))
        density_qn%dm_list(:, :, :, 1) = other
        density_qn%dm_list(:, :, :, 2) = dm
        if (.not. density_qn%density_in_history(dm)) then
            write(stderr, *) "test_density_in_history failed: Density matrix held "// &
                "in history not found."
            test_density_in_history = .false.
        end if
        if (.not. density_qn%density_in_history(dm + 0.1_rp * numerical_zero)) then
            write(stderr, *) "test_density_in_history failed: Density matrix "// &
                "differing by much less than numerical zero not found."
            test_density_in_history = .false.
        end if
        if (density_qn%density_in_history(dm + 10.0_rp * numerical_zero)) then
            write(stderr, *) "test_density_in_history failed: Density matrix "// &
                "differing by much more than numerical zero found."
            test_density_in_history = .false.
        end if

        ! determine if the comparison is relative to the size of the density, where the
        ! densities are scaled after perturbing them so that the differences scale with
        ! them instead of vanishing in their precision
        density_qn%dm_list(:, :, :, 2) = scale * dm
        if (.not. density_qn%density_in_history(scale * &
                                                (dm + 0.1_rp * numerical_zero))) then
            write(stderr, *) "test_density_in_history failed: Large density matrix "// &
                "differing by much less than numerical zero relative to its size "// &
                "not found."
            test_density_in_history = .false.
        end if
        if (density_qn%density_in_history(scale * (dm + 10.0_rp * numerical_zero))) then
            write(stderr, *) "test_density_in_history failed: Large density matrix "// &
                "differing by much more than numerical zero relative to its size found."
            test_density_in_history = .false.
        end if

        ! deallocate density QN object
        deallocate(density_qn)

    end function test_density_in_history

    logical(c_bool) function test_prepend() bind(C)
        !
        ! this function tests the subroutine which prepends an array to a list of arrays
        !
        use otr_density_qn, only: prepend

        real(rp), allocatable :: list(:, :, :, :)
        real(rp) :: new_array(2, 1, 1), new_arrays(2, 1, 1, 2), expected(2, 1, 1, 3)

        ! assume tests pass
        test_prepend = .true.

        ! allocate empty list and initialize random arrays to be prepended
        allocate(list(2, 1, 1, 0))
        call random_number(new_array)
        call random_number(new_arrays)

        ! prepend array to empty list and determine if dimensions and values of
        ! resulting list match
        call prepend(list, new_array)
        if (size(list, 4) /= 1) then
            write(stderr, *) "test_prepend failed: Incorrect list dimensions after "// &
                "prepending to empty list."
            test_prepend = .false.
        end if
        if (norm2(list(:, :, :, 1) - new_array) > tol) then
            write(stderr, *) "test_prepend failed: Incorrect list values after "// &
                "prepending to empty list."
            test_prepend = .false.
        end if

        ! initialize expected list after prepending two further arrays, the latest
        ! first
        expected(:, :, :, 1) = new_arrays(:, :, :, 2)
        expected(:, :, :, 2) = new_arrays(:, :, :, 1)
        expected(:, :, :, 3) = new_array

        ! prepend two further arrays and determine if dimensions and values of
        ! resulting list match, so that new arrays are added at the front while the
        ! order of the existing arrays is retained
        call prepend(list, new_arrays(:, :, :, 1))
        call prepend(list, new_arrays(:, :, :, 2))
        if (size(list, 4) /= 3) then
            write(stderr, *) "test_prepend failed: Incorrect list dimensions after "// &
                "prepending to non-empty list."
            test_prepend = .false.
        end if
        if (norm2(list - expected) > tol) then
            write(stderr, *) "test_prepend failed: Incorrect list values after "// &
                "prepending to non-empty list."
            test_prepend = .false.
        end if

        ! deallocate list
        deallocate(list)

    end function test_prepend

end module otr_density_qn_unit_tests
