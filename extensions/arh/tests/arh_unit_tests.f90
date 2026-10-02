! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! multipliers of the density matrix returned by the mock density matrix evaluating
    ! functions, which differ between the first and any subsequent call, following the
    ! same convention as the shared multiplier for the Fock matrix
    real(rp), parameter :: mock_v_same_spin_factor(2) = [3.0_rp, 7.0_rp], &
                           mock_v_opposite_spin_factor(2) = [4.0_rp, 9.0_rp]

    ! the non-linear potential keeps a single multiplier across calls, unlike the
    ! potentials above, so that it stays one consistent function of the density, which
    ! is what the history screening measures and would otherwise reject as noise
    real(rp), parameter :: mock_v_nonlinear_factor = 6.0_rp

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

    subroutine mock_evaluate_dm_cs(dm, energy, fock, v_nonlinear, error)
        !
        ! this subroutine is a mock density matrix evaluating function with a separate
        ! non-linear potential contribution for the closed-shell case, which returns
        ! multiples of the density matrix that change between calls so that
        ! non-vanishing differences are produced
        !
        use otr_common_unit_tests, only: record_mock_call, mock_factor, mock_fock_factor

        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target, contiguous :: fock(:, :), &
                                                               v_nonlinear(:, :)
        integer(ip), intent(out) :: error

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(v_nonlinear)))

        error = 0
        energy = sum(dm)
        if (present(fock)) fock = mock_potential(mock_factor(mock_fock_factor), dm)
        if (present(v_nonlinear)) &
            v_nonlinear = mock_potential(mock_v_nonlinear_factor, dm)

    end subroutine mock_evaluate_dm_cs

    subroutine mock_evaluate_dm_os(dm, energy, fock, v_same_spin, v_opposite_spin, &
                                   v_nonlinear, error)
        !
        ! this subroutine is a mock density matrix evaluating function with
        ! spin-resolved and non-linear potential contributions for the open-shell case,
        ! which returns multiples of the density matrix that change between calls so
        ! that non-vanishing differences are produced
        !
        use otr_common_unit_tests, only: record_mock_call, mock_factor, mock_fock_factor

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target :: &
            fock(:, :, :), v_same_spin(:, :, :), v_opposite_spin(:, :, :), &
            v_nonlinear(:, :, :)
        integer(ip), intent(out) :: error

        integer(ip) :: i

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(v_same_spin)) + &
                              merge(4_ip, 0_ip, present(v_opposite_spin)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))

        error = 0
        energy = sum(dm)
        do i = 1, size(dm, 3)
            if (present(fock)) fock(:, :, i) = &
                mock_potential(mock_factor(mock_fock_factor), dm(:, :, i))
            if (present(v_same_spin)) v_same_spin(:, :, i) = &
                mock_potential(mock_factor(mock_v_same_spin_factor), dm(:, :, i))
            if (present(v_opposite_spin)) v_opposite_spin(:, :, i) = &
                mock_potential(mock_factor(mock_v_opposite_spin_factor), dm(:, :, i))
            if (present(v_nonlinear)) v_nonlinear(:, :, i) = &
                mock_potential(mock_v_nonlinear_factor, dm(:, :, i))
        end do

    end subroutine mock_evaluate_dm_os

    subroutine generate_random_fock_partition(n_ao, n_occ, dm_oao, fock_oo, fock_vv)
        !
        ! this subroutine generates a random density matrix with the given number of
        ! occupied orbitals per particle/spin channel and the corresponding
        ! occupied-occupied and virtual-virtual Fock matrix blocks
        !
        use otr_common_unit_tests, only: &
            identity_matrix, generate_random_density_matrix, generate_random_symm_matrix

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

    function generate_random_nonredundant_vector(n_param, n_particle, n_ao, dm_oao) &
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
        vector = ref_pack_asymm(ref_project_asymm( &
            ref_unpack_asymm(vector, n_particle, n_ao), dm_oao), n_param)

    end function generate_random_nonredundant_vector

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

    subroutine setup_arh_and_mo_objects(mo_coeff, ao_overlap, n_occ)
        !
        ! this subroutine sets up the module-global MO and ARH objects for orbitals
        ! parameterized in the MO basis the way the MO factory would, pointing to the
        ! test-local MO coefficients
        !
        use otr_arh, only: arh_object, arh_mo_type
        use otr_mo, only: mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_unit_tests, only: setup_mo_object

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:)

        ! set up the MO object
        call setup_mo_object(mo_coeff, ao_overlap, n_occ)

        ! set up the ARH object in the MO basis
        allocate(arh_object, source=arh_mo_type(mo_object))
        call setup_settings(arh_object%settings)

    end subroutine setup_arh_and_mo_objects

    subroutine setup_arh_and_oao_objects(arh_type, dm_oao, fock_oo, fock_vv, n_ao, &
                                         n_particle, n_param)
        !
        ! this subroutine sets up the module-global OAO and ARH objects for the
        ! preconditioner tests; the allocatable caches the individual preconditioner
        ! branches need are filled in by the caller
        !
        use otr_arh, only: arh_object, arh_oao_type
        use otr_oao, only: oao_object
        use opentrustregion_unit_tests, only: setup_settings

        character(len=*), intent(in) :: arh_type
        real(rp), intent(in) :: dm_oao(:, :, :), fock_oo(:, :, :), fock_vv(:, :, :)
        integer(ip), intent(in) :: n_ao, n_particle, n_param

        ! set up the OAO object so that the static-Hessian eigendecomposition used by
        ! the preconditioner can be refreshed
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%dm_oao = dm_oao
        oao_object%fock_oo = fock_oo
        oao_object%fock_vv = fock_vv
        oao_object%hess_eigen_stale = .true.

        ! set up the ARH object in the OAO basis
        allocate(arh_object, source=arh_oao_type(oao_object))
        call setup_settings(arh_object%settings)
        arh_object%settings%arh_type = arh_type

    end subroutine setup_arh_and_oao_objects

    subroutine setup_arh_and_scaled_oao_objects(basis_scale, dm_ao, dm_oao)
        !
        ! this subroutine sets up the module-global OAO object for random density
        ! matrices in an AO basis whose overlap has the inverse square root c I, so
        ! that the transformation to the OAO basis scales by c^2, and the ARH object
        ! pointing to it the way the ARH factory would
        !
        use otr_arh, only: arh_object, arh_oao_type
        use otr_oao, only: oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_unit_tests, only: identity_matrix, generate_random_density_matrix
        use otr_common_test_reference, only: n_occ

        real(rp), intent(in) :: basis_scale
        real(rp), intent(out), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(out) :: dm_oao(:, :, :)

        integer(ip) :: n_ao, n_particle, n_param, i

        ! set up the OAO object
        n_ao = size(dm_ao, 1, kind=ip)
        n_particle = size(dm_ao, 3, kind=ip)
        n_param = n_particle * n_ao * (n_ao - 1) / 2
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%s_inv_sqrt = basis_scale * identity_matrix(n_ao)
        do i = 1, n_particle
            dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do
        dm_ao = basis_scale**2 * dm_oao
        oao_object%dm_ao => dm_ao
        oao_object%dm_oao = dm_oao
        allocate(oao_object%fock_oo(n_ao, n_ao, n_particle), &
                 oao_object%fock_vv(n_ao, n_ao, n_particle), oao_object%grad(n_param), &
                 oao_object%h_diag(n_param))

        ! set up the ARH object
        allocate(arh_object, source=arh_oao_type(oao_object))
        call setup_settings(arh_object%settings)

    end subroutine setup_arh_and_scaled_oao_objects

    subroutine setup_empty_history(n_particle)
        !
        ! this subroutine sets up the module-global ARH object with an empty history
        ! and vanishing potentials at the current point for the closed- or open-shell
        ! case
        !
        use otr_arh, only: arh_object

        integer(ip), intent(in) :: n_particle

        integer(ip) :: n_ao

        n_ao = arh_object%orbitals%n_ao
        allocate(arh_object%dm_list(n_ao, n_ao, n_particle, 0), &
                 arh_object%v_nonlinear_list(n_ao, n_ao, n_particle, 0))
        allocate(arh_object%v_nonlinear(n_ao, n_ao, n_particle))
        arh_object%v_nonlinear = 0.0_rp
        if (n_particle == 1) then
            arh_object%evaluate_dm_cs => mock_evaluate_dm_cs
            allocate(arh_object%fock_list(n_ao, n_ao, n_particle, 0), &
                     arh_object%fock(n_ao, n_ao, n_particle))
            arh_object%fock = 0.0_rp
        else
            arh_object%evaluate_dm_os => mock_evaluate_dm_os
            allocate(arh_object%v_same_spin_list(n_ao, n_ao, n_particle, 0), &
                     arh_object%v_opposite_spin_list(n_ao, n_ao, n_particle, 0), &
                     arh_object%v_same_spin(n_ao, n_ao, n_particle), &
                     arh_object%v_opposite_spin(n_ao, n_ao, n_particle))
            arh_object%v_same_spin = 0.0_rp
            arh_object%v_opposite_spin = 0.0_rp
        end if

    end subroutine setup_empty_history

    function stacked_potentials(n_particle, k) result(potentials)
        !
        ! this function stacks the potentials the ARH object keeps for the closed- or
        ! open-shell case along the last dimension, the current ones for k = 0 and
        ! those of the k-th history entry otherwise
        !
        use otr_arh, only: arh_object

        integer(ip), intent(in) :: n_particle, k
        real(rp), allocatable :: potentials(:, :, :, :)

        if (n_particle == 1 .and. k == 0) then
            potentials = reshape([arh_object%fock, arh_object%v_nonlinear], &
                                 [shape(arh_object%v_nonlinear), 2])
        else if (n_particle == 1) then
            potentials = reshape([arh_object%fock_list(:, :, :, k), &
                                  arh_object%v_nonlinear_list(:, :, :, k)], &
                                 [shape(arh_object%v_nonlinear_list(:, :, :, k)), 2])
        else if (k == 0) then
            potentials = reshape([arh_object%v_same_spin, arh_object%v_opposite_spin, &
                                  arh_object%v_nonlinear], &
                                 [shape(arh_object%v_nonlinear), 3])
        else
            potentials = reshape([arh_object%v_same_spin_list(:, :, :, k), &
                                  arh_object%v_opposite_spin_list(:, :, :, k), &
                                  arh_object%v_nonlinear_list(:, :, :, k)], &
                                 [shape(arh_object%v_nonlinear_list(:, :, :, k)), 3])
        end if

    end function stacked_potentials

    function potential_factors(n_particle) result(factors)
        !
        ! this function returns the multipliers the mock density matrix evaluating
        ! functions apply on their first call to the potentials the ARH object keeps
        ! for the closed- or open-shell case, in the order of stacked_potentials
        !
        use otr_common_unit_tests, only: mock_fock_factor

        integer(ip), intent(in) :: n_particle
        real(rp), allocatable :: factors(:)

        if (n_particle == 1) then
            factors = [mock_fock_factor(1), mock_v_nonlinear_factor]
        else
            factors = [mock_v_same_spin_factor(1), mock_v_opposite_spin_factor(1), &
                       mock_v_nonlinear_factor]
        end if

    end function potential_factors

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

    function check_arh_factory(stage, test_name, mo_basis, n_particle, error, &
                               settings, solver_settings, obj_func_funptr, &
                               update_orbs_funptr) result(passed)
        !
        ! this function checks the state an ARH factory leaves behind: for a new
        ! calculation (stage 1) the ARH object for the orbital basis pointing to the
        ! orbital object, the stored settings and density matrix evaluating function,
        ! the returned functions of the shell and the ARH routines wired into the
        ! solver settings, for a new calculation after a previous one (stage 2) the
        ! discarded history and evaluation state, and for an unknown ARH type (stage 3)
        ! the rejection before the orbital object is changed
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_object, arh_mo_type, arh_oao_type, arh_settings_type, &
                           obj_func_arh_cs_callback_ptr, obj_func_arh_os_callback_ptr, &
                           update_orbs_arh_cs_callback_ptr, &
                           update_orbs_arh_os_callback_ptr, precond_arh_callback_ptr
        use otr_mo, only: mo_object
        use otr_oao, only: oao_object, project_oao_callback_ptr
        use otr_arh_test_reference, only: operator(==)

        integer(ip), intent(in) :: stage
        character(len=*), intent(in) :: test_name
        logical, intent(in) :: mo_basis
        integer(ip), intent(in) :: n_particle, error
        type(arh_settings_type), intent(in) :: settings
        type(solver_settings_type), intent(in) :: solver_settings
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        logical :: passed

        logical :: correct_basis, correct_orbitals, stored, returned

        ! assume test passes
        passed = .true.

        ! an unknown ARH type has to be rejected before the orbital object is changed
        if (stage == 3) then
            if (error == 0) then
                write(stderr, *) test_name// &
                    " failed: Error not thrown for unknown ARH type."
                passed = .false.
            end if
            if (arh_object%orbitals%hess_eigen_stale) then
                write(stderr, *) test_name// &
                    " failed: Orbital object changed for an unknown ARH type."
                passed = .false.
            end if
            return
        end if

        ! the checks below read the ARH object
        if (error /= 0) then
            write(stderr, *) test_name//" failed: Produced error."
            passed = .false.
            return
        end if
        if (.not. allocated(arh_object)) then
            write(stderr, *) test_name//" failed: ARH object not allocated."
            passed = .false.
            return
        end if

        ! a new calculation after a previous one has to discard its history and
        ! evaluation state
        if (stage == 2) then
            if (allocated(arh_object%dm_list)) then
                write(stderr, *) test_name// &
                    " failed: History of the previous calculation kept."
                passed = .false.
            end if
            if (.not. arh_object%evaluation_stale) then
                write(stderr, *) test_name// &
                    " failed: Evaluation state of the previous calculation kept."
                passed = .false.
            end if
            return
        end if

        ! determine if the orbitals are parameterized in the right basis and point to
        ! the object holding them
        correct_basis = .false.
        select type (arh => arh_object)
        type is (arh_mo_type)
            correct_basis = mo_basis
        type is (arh_oao_type)
            correct_basis = .not. mo_basis
        end select
        if (.not. correct_basis) then
            write(stderr, *) test_name// &
                " failed: Orbitals not parameterized in the right basis."
            passed = .false.
        end if
        correct_orbitals = .false.
        if (mo_basis) then
            if (allocated(mo_object)) &
                correct_orbitals = associated(arh_object%orbitals, mo_object)
        else
            if (allocated(oao_object)) &
                correct_orbitals = associated(arh_object%orbitals, oao_object)
        end if
        if (.not. correct_orbitals) then
            write(stderr, *) test_name// &
                " failed: Orbitals not associated with the orbital object."
            passed = .false.
        end if

        ! determine if the settings and the density matrix evaluating function are
        ! stored, the functions of the shell are returned and the ARH routines are
        ! wired into the solver settings
        if (.not. (arh_object%settings == settings)) then
            write(stderr, *) test_name//" failed: Settings not stored correctly."
            passed = .false.
        end if
        if (n_particle == 1) then
            stored = associated(arh_object%evaluate_dm_cs, mock_evaluate_dm_cs)
            returned = associated(obj_func_funptr, obj_func_arh_cs_callback_ptr) .and. &
                       associated(update_orbs_funptr, update_orbs_arh_cs_callback_ptr)
        else
            stored = associated(arh_object%evaluate_dm_os, mock_evaluate_dm_os)
            returned = associated(obj_func_funptr, obj_func_arh_os_callback_ptr) .and. &
                       associated(update_orbs_funptr, update_orbs_arh_os_callback_ptr)
        end if
        if (.not. stored) then
            write(stderr, *) test_name// &
                " failed: Density matrix evaluating function not stored correctly."
            passed = .false.
        end if
        if (.not. returned) then
            write(stderr, *) test_name//" failed: Returned function pointers are wrong."
            passed = .false.
        end if
        if (.not. associated(solver_settings%precond, precond_arh_callback_ptr)) then
            write(stderr, *) test_name// &
                " failed: ARH routines not wired into solver settings."
            passed = .false.
        end if
        if ((associated(solver_settings%project, project_oao_callback_ptr) .and. &
             associated(solver_settings%stability_settings%project, &
                        project_oao_callback_ptr)) .eqv. mo_basis) then
            write(stderr, *) test_name//" failed: Projection of the OAO basis not "// &
                "wired in the OAO basis or wired in the MO basis."
            passed = .false.
        end if

    end function check_arh_factory

    function check_obj_func_arh(n_particle, test_name) result(passed)
        !
        ! this function checks the energy evaluation for the closed- or open-shell case
        ! for orbitals parameterized in the OAO basis, which also adds the evaluated
        ! point together with the potentials the shell keeps to the history
        !
        use opentrustregion, only: obj_func_type
        use otr_arh, only: obj_func_arh_cs_callback, obj_func_arh_os_callback, &
                           arh_object
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao
        use otr_common_unit_tests, only: mock_requests

        integer(ip), intent(in) :: n_particle
        character(len=*), intent(in) :: test_name
        logical :: passed

        real(rp), parameter :: basis_scale = 0.5_rp

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    kappa(n_particle * n_ao * (n_ao - 1) / 2), energy
        real(rp), allocatable :: factors(:), potentials(:, :, :, :)
        integer(ip) :: k, error
        procedure(obj_func_type), pointer :: obj_func

        ! assume test passes
        passed = .true.

        ! set up the OAO and ARH objects with an empty history, which is removed for
        ! the first call, and select the routine of the shell
        call setup_arh_and_scaled_oao_objects(basis_scale, dm_ao, dm_oao)
        call setup_empty_history(n_particle)
        deallocate(arh_object%dm_list)
        if (n_particle == 1) then
            obj_func => obj_func_arh_cs_callback
        else
            obj_func => obj_func_arh_os_callback
        end if
        factors = potential_factors(n_particle)

        ! call routine with an orbital rotation before the history exists and determine
        ! if the potentials the shell keeps are requested and nothing is added to the
        ! history
        call random_number(kappa)
        kappa = 0.1_rp * kappa
        mock_requests = [integer(ip) :: ]
        energy = obj_func(kappa, error)
        if (error /= 0) then
            write(stderr, *) test_name//" failed: Produced error without history."
            passed = .false.
        end if
        if (size(mock_requests) /= 1) then
            write(stderr, *) test_name//" failed: Density matrix evaluating "// &
                "function not called exactly once."
            passed = .false.
        end if
        if (any(mock_requests /= merge(3_ip, 14_ip, n_particle == 1))) then
            write(stderr, *) test_name//" failed: Incorrect outputs requested from "// &
                "density matrix evaluating function."
            passed = .false.
        end if
        if (allocated(arh_object%dm_list)) then
            write(stderr, *) test_name//" failed: History created by energy evaluation."
            passed = .false.
        end if
        if (arh_object%model_stale) then
            write(stderr, *) test_name// &
                " failed: Hessian model marked stale without history."
            passed = .false.
        end if

        ! call routine again with an empty history and determine if the rotated density
        ! matrix is added together with the potentials evaluated at it in the history
        ! basis, while the current density matrix is left untouched
        allocate(arh_object%dm_list(n_ao, n_ao, n_particle, 0))
        mock_requests = [integer(ip) :: ]
        energy = obj_func(kappa, error)
        if (error /= 0) then
            write(stderr, *) test_name//" failed: Produced error with history."
            passed = .false.
        else if (size(arh_object%dm_list, 4) /= 1) then
            write(stderr, *) test_name//" failed: Evaluated point not added to history."
            passed = .false.
        else
            if (abs(energy - basis_scale**2 * sum(arh_object%dm_list(:, :, :, 1))) > &
                tol) then
                write(stderr, *) test_name//" failed: Incorrect energy."
                passed = .false.
            end if
            if (norm2(arh_object%dm_list(:, :, :, 1) - dm_oao) < tol) then
                write(stderr, *) test_name//" failed: Incorrect density matrix history."
                passed = .false.
            end if
            potentials = stacked_potentials(n_particle, 1_ip)
            do k = 1, size(factors)
                if (norm2(potentials(:, :, :, k) - basis_scale**4 * mock_potential( &
                    factors(k), arh_object%dm_list(:, :, :, 1))) > tol) then
                    write(stderr, *) test_name//" failed: Incorrect potential history."
                    passed = .false.
                end if
            end do
        end if
        if (norm2(oao_object%dm_oao - dm_oao) > tol .or. &
            norm2(arh_object%orbitals%dm_ao - basis_scale**2 * dm_oao) > tol) then
            write(stderr, *) test_name// &
                " failed: Current density matrix changed by energy evaluation."
            passed = .false.
        end if
        if (.not. arh_object%model_stale) then
            write(stderr, *) test_name//" failed: Hessian model not marked stale "// &
                "after adding the evaluated point to the history."
            passed = .false.
        end if

        ! call routine at the same point and determine if it is not added twice and
        ! leaves the Hessian model untouched
        arh_object%model_stale = .false.
        energy = obj_func(kappa, error)
        if (error /= 0) then
            write(stderr, *) test_name// &
                " failed: Produced error for a point the history already holds."
            passed = .false.
        end if
        if (size(arh_object%dm_list, 4) /= 1) then
            write(stderr, *) test_name// &
                " failed: History extended for a point it already holds."
            passed = .false.
        end if
        if (arh_object%model_stale) then
            write(stderr, *) test_name//" failed: Hessian model marked stale for a "// &
                "point the history already holds."
            passed = .false.
        end if

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    end function check_obj_func_arh

    function check_update_orbs_arh(n_particle, test_name) result(passed)
        !
        ! this function checks the energy, gradient and Hessian diagonal evaluation for
        ! the closed- or open-shell case for orbitals parameterized in the OAO basis,
        ! which passes the Fock matrix in the history basis on to the orbital basis and
        ! adds the point the orbitals are rotated away from to the history
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_arh, only: update_orbs_arh_cs_callback, update_orbs_arh_os_callback, &
                           arh_object, hess_x_arh_callback_ptr
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao
        use otr_common_unit_tests, only: identity_matrix, mock_requests, &
                                         mock_fock_factor

        integer(ip), intent(in) :: n_particle
        character(len=*), intent(in) :: test_name
        logical :: passed

        real(rp), parameter :: basis_scale = 0.5_rp

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: dm_oao(n_ao, n_ao, n_particle), fock(n_ao, n_ao, n_particle), &
                    dm_saved(n_ao, n_ao, n_particle, 2), complement(n_ao, n_ao), func
        real(rp), allocatable :: kappa(:), grad(:), h_diag(:), factors(:), &
                                 potentials(:, :, :, :), potentials_saved(:, :, :, :, :)
        integer(ip) :: n_param, n_list, k, error
        procedure(update_orbs_type), pointer :: update_orbs
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume test passes
        passed = .true.

        ! set up the OAO and ARH objects and select the routine of the shell
        call setup_arh_and_scaled_oao_objects(basis_scale, dm_ao, dm_oao)
        arh_object%settings%arh_type = "ms_psb"
        if (n_particle == 1) then
            update_orbs => update_orbs_arh_cs_callback
            arh_object%evaluate_dm_cs => mock_evaluate_dm_cs
        else
            update_orbs => update_orbs_arh_os_callback
            arh_object%evaluate_dm_os => mock_evaluate_dm_os
        end if
        factors = potential_factors(n_particle)
        n_param = arh_object%orbitals%n_param
        allocate(kappa(n_param), grad(n_param), h_diag(n_param), &
                 potentials_saved(n_ao, n_ao, n_particle, size(factors), 2))

        ! the checks after a failed call would read quantities it did not set
        checks: block
            ! call routine without an orbital rotation for an uninitialized object
            kappa = 0.0_rp
            mock_requests = [integer(ip) :: ]
            call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
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
            potentials = stacked_potentials(n_particle, 0_ip)
            do k = 1, size(factors)
                if (norm2(potentials(:, :, :, k) - basis_scale**2 * &
                          mock_potential(factors(k), dm_ao)) > tol) then
                    write(stderr, *) test_name//" failed: Incorrect potential."
                    passed = .false.
                end if
            end do
            fock = basis_scale**2 * mock_potential(mock_fock_factor(1), dm_ao)
            do k = 1, n_particle
                complement = identity_matrix(n_ao) - oao_object%dm_oao(:, :, k)
                if (norm2(oao_object%fock_oo(:, :, k) - &
                          matmul(oao_object%dm_oao(:, :, k), &
                                 matmul(fock(:, :, k), oao_object%dm_oao(:, :, k)))) + &
                    norm2(oao_object%fock_vv(:, :, k) - matmul( &
                        complement, matmul(fock(:, :, k), complement))) > tol) then
                    write(stderr, *) test_name//" failed: Fock matrix in the "// &
                        "history basis not passed on to the orbital basis."
                    passed = .false.
                end if
            end do

            ! determine if the gradient, Hessian diagonal and Hessian linear
            ! transformation are returned
            if (norm2(grad - arh_object%orbitals%grad) > tol) then
                write(stderr, *) test_name//" failed: Gradient not returned."
                passed = .false.
            end if
            if (norm2(h_diag - arh_object%orbitals%h_diag) > tol) then
                write(stderr, *) test_name//" failed: Hessian diagonal not returned."
                passed = .false.
            end if
            if (.not. associated(hess_x_funptr, hess_x_arh_callback_ptr)) then
                write(stderr, *) test_name// &
                    " failed: Returned Hessian linear transformation is wrong."
                passed = .false.
            end if

            ! determine if the history is initialized empty and the approximate Hessian
            ! model is assembled from it
            if (size(arh_object%dm_list, 4) /= 0) then
                write(stderr, *) test_name// &
                    " failed: Density matrix history not initialized empty."
                passed = .false.
            end if
            if (.not. allocated(arh_object%dm_dirs)) then
                write(stderr, *) test_name//" failed: Approximate Hessian model "// &
                    "not assembled for the initial point."
                passed = .false.
            end if

            ! call routine again without an orbital rotation and determine if the
            ! quantities of the already initialized object are reused, even for a
            ! vanishing energy, without rebuilding the static part of the Hessian
            oao_object%energy = 0.0_rp
            oao_object%hess_eigen_stale = .false.
            call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
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
            if (size(arh_object%dm_list, 4) /= 0) then
                write(stderr, *) test_name// &
                    " failed: History extended without an orbital rotation."
                passed = .false.
            end if
            if (oao_object%hess_eigen_stale) then
                write(stderr, *) test_name//" failed: Static part of the Hessian "// &
                    "rebuilt without an orbital rotation."
                passed = .false.
            end if

            ! rotate twice in a row, the second time for a model marked stale, and save
            ! the quantities at the points the orbitals are rotated away from, which
            ! the history has to retain with the latest first
            kappa = 0.1_rp
            do k = 2, 1, -1
                dm_saved(:, :, :, k) = oao_object%dm_oao
                potentials_saved(:, :, :, :, k) = stacked_potentials(n_particle, 0_ip)
                arh_object%model_stale = k == 1
                call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
                if (error /= 0) then
                    write(stderr, *) test_name// &
                        " failed: Produced error for an orbital rotation."
                    passed = .false.
                    exit checks
                end if
            end do
            if (size(arh_object%dm_list, 4) /= 2) then
                write(stderr, *) test_name// &
                    " failed: History not extended to two entries."
                passed = .false.
                exit checks
            end if
            do k = 1, 2
                if (norm2(arh_object%dm_list(:, :, :, k) - dm_saved(:, :, :, k)) > &
                    tol) then
                    write(stderr, *) test_name// &
                        " failed: Incorrect density matrix history."
                    passed = .false.
                end if
                if (norm2(stacked_potentials(n_particle, k) - &
                          potentials_saved(:, :, :, :, k)) > tol) then
                    write(stderr, *) test_name//" failed: Incorrect potential history."
                    passed = .false.
                end if
            end do

            ! determine if the approximate Hessian model is assembled from the history
            ! extended by the point the orbitals were rotated away from
            if (arh_object%model_stale .or. .not. allocated(arh_object%dm_dirs)) then
                write(stderr, *) test_name// &
                    " failed: Approximate Hessian model not assembled."
                passed = .false.
            else if (size(arh_object%dm_dirs, 2) /= 2 * n_particle) then
                write(stderr, *) test_name//" failed: Approximate Hessian model "// &
                    "not assembled from the extended history."
                passed = .false.
            end if

            ! determine if a rotation starting from a density matrix which the history
            ! already holds, as after the energy evaluation of an accepted trial point,
            ! does not add that density matrix a second time
            arh_object%dm_list(:, :, :, 1) = oao_object%dm_oao
            call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
            if (error /= 0) then
                write(stderr, *) test_name//" failed: Produced error for a density "// &
                    "matrix already held in the history."
                passed = .false.
            else if (size(arh_object%dm_list, 4) /= 2) then
                write(stderr, *) test_name//" failed: History extended for a "// &
                    "density matrix it already holds."
                passed = .false.
            end if

            ! call routine with an orbital rotation whose evaluation fails and
            ! determine if the error is passed on and the rotated density is marked as
            ! not evaluated, so that a following call without an orbital rotation
            ! evaluates it without adding it to the history
            if (n_particle == 1) then
                arh_object%evaluate_dm_cs => mock_evaluate_dm_failing_cs
            else
                arh_object%evaluate_dm_os => mock_evaluate_dm_failing_os
            end if
            call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
            if (error == 0) then
                write(stderr, *) test_name//" failed: Error of the density matrix "// &
                    "evaluating function not passed on."
                passed = .false.
            end if
            if (.not. arh_object%evaluation_stale) then
                write(stderr, *) test_name// &
                    " failed: Evaluation not marked stale after a failed evaluation."
                passed = .false.
            end if
            if (n_particle == 1) then
                arh_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                arh_object%evaluate_dm_os => mock_evaluate_dm_os
            end if
            n_list = size(arh_object%dm_list, 4, kind=ip)
            kappa = 0.0_rp
            call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
            if (error /= 0) then
                write(stderr, *) test_name// &
                    " failed: Produced error after a failed evaluation."
                passed = .false.
            end if
            if (arh_object%evaluation_stale) then
                write(stderr, *) test_name//" failed: Evaluation still marked "// &
                    "stale after a failed evaluation was recovered."
                passed = .false.
            end if
            if (size(arh_object%dm_list, 4) /= n_list) then
                write(stderr, *) test_name//" failed: History extended by a "// &
                    "density matrix which has not been evaluated."
                passed = .false.
            end if

            ! determine if every recompute asked for all potentials the shell needs
            if (any(mock_requests /= merge(3_ip, 15_ip, n_particle == 1))) then
                write(stderr, *) test_name//" failed: Incorrect outputs requested "// &
                    "from density matrix evaluating function."
                passed = .false.
            end if
        end block checks

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

    contains

        subroutine mock_evaluate_dm_failing_cs(dm, energy_out, fock_out, &
                                               v_nonlinear_out, error_out)
            !
            ! this subroutine is a mock density matrix evaluating function for the
            ! closed-shell case, which fails
            !
            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out) :: energy_out
            real(rp), intent(out), optional, target, contiguous :: fock_out(:, :), &
                                                                   v_nonlinear_out(:, :)
            integer(ip), intent(out) :: error_out

            error_out = 1
            energy_out = sum(dm)
            if (present(fock_out)) fock_out = 0.0_rp
            if (present(v_nonlinear_out)) v_nonlinear_out = 0.0_rp

        end subroutine mock_evaluate_dm_failing_cs

        subroutine mock_evaluate_dm_failing_os(dm, energy_out, fock_out, &
                                               v_same_spin_out, v_opposite_spin_out, &
                                               v_nonlinear_out, error_out)
            !
            ! this subroutine is a mock density matrix evaluating function for the
            ! open-shell case, which fails
            !
            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out) :: energy_out
            real(rp), intent(out), optional, target :: &
                fock_out(:, :, :), v_same_spin_out(:, :, :), &
                v_opposite_spin_out(:, :, :), v_nonlinear_out(:, :, :)
            integer(ip), intent(out) :: error_out

            error_out = 1
            energy_out = sum(dm)
            if (present(fock_out)) fock_out = 0.0_rp
            if (present(v_same_spin_out)) v_same_spin_out = 0.0_rp
            if (present(v_opposite_spin_out)) v_opposite_spin_out = 0.0_rp
            if (present(v_nonlinear_out)) v_nonlinear_out = 0.0_rp

        end subroutine mock_evaluate_dm_failing_os

    end function check_update_orbs_arh

    function check_build_hess_model(n_particle, test_name) result(passed)
        !
        ! this function checks the assembly of the approximate Hessian model from the
        ! history for the closed- or open-shell case: standard ARH has to reproduce a
        ! history in the MO basis exactly, while every ARH type in the OAO basis has to
        ! keep the whole history in the linear system, drop the entries beyond the
        ! step-length cutoff and those whose response is noise from the non-linear
        ! system and cache every quantity in the basis of its own system, without
        ! evaluating the density matrix or changing the history
        !
        use otr_arh, only: build_hess_model_cs, build_hess_model_os, arh_object, &
                           arh_types, arh_oao_type
        use otr_mo, only: mo_object
        use otr_oao, only: oao_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_unit_tests, only: &
            identity_matrix, mock_requests, generate_random_density_matrix, &
            generate_random_symm_matrix, generate_random_orthogonal_matrix
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, ref_mo_transform, &
                                     ref_pack_ov, ref_unpack_ov

        integer(ip), intent(in) :: n_particle
        character(len=*), intent(in) :: test_name
        logical :: passed

        integer(ip), parameter :: n_hist = 3, n_noisy = 3
        real(rp), parameter :: noisy_scales(n_noisy) = [1e-3_rp, 3e-3_rp, 1e-2_rp], &
                               duplicate_offset = 1e-3_rp

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), step_mo(n_mo, n_mo, n_particle), &
                    v_diff_mo(n_mo, n_mo, n_particle), perturbation(n_ao, n_ao)
        real(rp), allocatable :: steps(:, :), orth(:, :), response(:), factors(:), &
                                 dm_list(:, :, :, :), current(:, :, :, :), &
                                 lists(:, :, :, :, :)
        integer(ip) :: n_pot, n_list, n_linear, n_nonlinear, n_systems, i, j, k, &
                       i_case, i_type, i_noisy, error
        logical :: built
        character(len=:), allocatable :: arh_type, case_name

        ! assume test passes
        passed = .true.

        ! the closed-shell case keeps the Fock matrix and the open-shell case the
        ! same-spin and opposite-spin potentials besides the non-linear potential
        n_pot = n_particle + 1

        ! set up the ARH object in the MO basis at random orbitals
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_arh_and_mo_objects(mo_coeff, ao_overlap, n_occ(:n_particle))
        arh_object%settings%arh_type = "arh"

        ! steps of every channel which are not orthogonal, so that the pivoted
        ! factorization of the history reorders them and its factor is not diagonal
        allocate(steps(mo_object%n_param, n_hist))
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
            step_mo = ref_unpack_ov(steps(:, k), n_occ(:n_particle), n_mo)
            do j = 1, n_particle
                dm_list(:, :, j, k) = matmul(ao_overlap, matmul( &
                    mo_object%dm_ao(:, :, j) + matmul(mo_coeff(:, :, j), matmul( &
                        step_mo(:, :, j) + transpose(step_mo(:, :, j)), &
                        transpose(mo_coeff(:, :, j)))), ao_overlap))
                do i = 1, n_pot - 1
                    lists(:, :, j, k, i) = generate_random_symm_matrix(n_ao)
                end do
            end do
            lists(:, :, :, k, n_pot) = current(:, :, :, n_pot)
        end do

        ! build the model and determine if the low-rank part reproduces the linear
        ! potential difference of every step, which pins the scaling of the low-rank
        ! part and the (map, chol) pair its directions are rebased with
        call build_from_history()
        if (error /= 0) then
            write(stderr, *) test_name//" failed: Produced error for the MO basis."
            passed = .false.
        end if
        if (arh_object%model_stale) then
            write(stderr, *) test_name// &
                " failed: Hessian model still marked stale for the MO basis."
            passed = .false.
        end if
        if (.not. allocated(arh_object%coupling_matrix)) then
            write(stderr, *) test_name//" failed: No low-rank part for the MO basis."
            passed = .false.
        else
            do k = 1, n_hist
                response = matmul(arh_object%expansion_dirs, matmul( &
                    arh_object%coupling_matrix, &
                    matmul(transpose(arh_object%projection_dirs), steps(:, k))))
                do j = 1, n_particle
                    v_diff_mo(:, :, j) = ref_mo_transform( &
                        mo_coeff(:, :, j), sum(lists(:, :, j, k, :n_pot - 1) - &
                                               current(:, :, j, :n_pot - 1), dim=3))
                end do
                if (norm2(response - merge(4.0_rp, 2.0_rp, n_particle == 1) * &
                          ref_pack_ov(v_diff_mo, n_occ(:n_particle))) > &
                    tol * (1.0_rp + norm2(response))) then
                    write(stderr, *) test_name// &
                        " failed: History not reproduced exactly for the MO basis."
                    passed = .false.
                end if
            end do
        end if
        deallocate(arh_object, mo_object, dm_list, current, lists)

        ! set up the OAO object with an orthonormal AO basis, so that the AO and the
        ! OAO basis coincide
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_particle * n_ao * (n_ao - 1) / 2
        oao_object%s_inv_sqrt = identity_matrix(n_ao)
        allocate(oao_object%dm_oao(n_ao, n_ao, n_particle))
        do j = 1, n_particle
            oao_object%dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_occ(j))
        end do

        ! build the model from a history of independent points, one spanning very
        ! different step lengths and one whose non-linear response is noise in one
        ! channel at a time
        factors = potential_factors(n_particle)
        do i_case = 1, 2 + n_particle
            i_noisy = i_case - 2
            n_list = merge(n_noisy, 2_ip, i_noisy > 0)
            allocate(dm_list(n_ao, n_ao, n_particle, n_list), &
                     current(n_ao, n_ao, n_particle, n_pot), &
                     lists(n_ao, n_ao, n_particle, n_list, n_pot))
            do j = 1, n_particle
                if (i_case == 1) then
                    case_name = "a history of independent points"
                    do k = 1, n_list
                        call random_number(perturbation)
                        dm_list(:, :, j, k) = oao_object%dm_oao(:, :, j) + 1e-2_rp * &
                                              (perturbation + transpose(perturbation))
                    end do
                else if (i_case == 2) then
                    case_name = "a history spanning very different step lengths"
                    call random_number(perturbation)
                    dm_list(:, :, j, 1) = oao_object%dm_oao(:, :, j) + 1e-6_rp * &
                                          (perturbation + transpose(perturbation))
                    dm_list(:, :, j, 2) = generate_random_density_matrix(n_ao, n_occ(j))
                else
                    case_name = "a history whose non-linear response is noise in "// &
                                "channel "//achar(iachar("0") + i_noisy)
                    do k = 1, n_list
                        call random_number(perturbation)
                        dm_list(:, :, j, k) = &
                            oao_object%dm_oao(:, :, j) + &
                            noisy_scales(k) * (perturbation + transpose(perturbation))
                    end do

                    ! the last entry barely departs from the direction of the one
                    ! before it, so that it stays linearly independent but contributes
                    ! a residual far below the noise of the response, while still being
                    ! admitted by the step-length screen
                    call random_number(perturbation)
                    dm_list(:, :, j, n_list) = &
                        dm_list(:, :, j, n_list - 1) + &
                        duplicate_offset * noisy_scales(n_list - 1) * &
                        (perturbation + transpose(perturbation))
                end if
            end do

            ! potentials which are functions of the density, except for the non-linear
            ! response of the channel under test, which is noise
            do i = 1, n_pot
                current(:, :, :, i) = mock_potential(factors(i), oao_object%dm_oao)
                do k = 1, n_list
                    lists(:, :, :, k, i) = mock_potential(factors(i), &
                                                          dm_list(:, :, :, k))
                end do
            end do
            if (i_noisy > 0) call random_number(lists(:, :, i_noisy, :, n_pot))

            ! build the model for every ARH type on a newly set up ARH object
            do i_type = 1, size(arh_types)
                arh_type = trim(arh_types(i_type))
                allocate(arh_object, source=arh_oao_type(oao_object))
                call setup_settings(arh_object%settings)
                arh_object%settings%arh_type = arh_type
                call build_from_history()

                ! determine if the model is assembled from the history without
                ! evaluating the density matrix or changing the history
                if (error /= 0) then
                    write(stderr, *) test_name//" failed: Produced error for "// &
                        arh_type//" and "//case_name//"."
                    passed = .false.
                end if
                if (arh_object%model_stale) then
                    write(stderr, *) test_name//" failed: Hessian model still "// &
                        "marked stale for "//arh_type//" and "//case_name//"."
                    passed = .false.
                end if
                if (size(mock_requests) /= 0) then
                    write(stderr, *) test_name// &
                        " failed: Density matrix evaluating function called for "// &
                        arh_type//" and "//case_name//"."
                    passed = .false.
                end if
                do k = 1, n_list
                    if (norm2(arh_object%dm_list(:, :, :, k) - dm_list(:, :, :, k)) > &
                        tol .or. norm2(stacked_potentials(n_particle, k) - &
                                       lists(:, :, :, k, :)) > tol) then
                        write(stderr, *) test_name//" failed: History changed for "// &
                            arh_type//" and "//case_name//"."
                        passed = .false.
                    end if
                end do
                if (.not. (allocated(arh_object%expansion_dirs) .and. &
                           allocated(arh_object%projection_dirs) .and. &
                           allocated(arh_object%coupling_matrix))) then
                    write(stderr, *) test_name//" failed: Low-rank Hessian factors "// &
                        "not assembled for "//arh_type//" and "//case_name//"."
                    passed = .false.
                end if

                ! determine if the linear and non-linear systems of the ARH type are
                ! constructed
                if (arh_type == "ms_sr1") then
                    built = allocated(arh_object%a_inv) .and. &
                            allocated(arh_object%a_inv_comb) .and. &
                            allocated(arh_object%linear_potential_dirs) .and. &
                            allocated(arh_object%nonlinear_potential_dirs)
                else
                    built = allocated(arh_object%dm_dirs) .and. &
                            allocated(arh_object%dm_dirs_nonlinear)
                    if (arh_type /= "ms_sp") built = &
                        built .and. allocated(arh_object%linear_potential_dirs) .and. &
                        allocated(arh_object%nonlinear_potential_dirs)
                    if (arh_type == "ms_sp" .or. arh_type == "ms_psb") &
                        built = built .and. allocated(arh_object%a_sym) .and. &
                                allocated(arh_object%a_sym_nonlinear)
                end if
                if (.not. built) then
                    write(stderr, *) test_name//" failed: Linear and non-linear "// &
                        "systems not constructed for "//arh_type//" and "//case_name// &
                        "."
                    passed = .false.
                    deallocate(arh_object)
                    cycle
                end if
                if (arh_type == "ms_sr1") then
                    n_linear = size(arh_object%linear_potential_dirs, 2)
                    n_nonlinear = size(arh_object%nonlinear_potential_dirs, 2)
                else
                    n_linear = size(arh_object%dm_dirs, 2)
                    n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
                end if

                ! determine if the linear system keeps the entire history, while the
                ! non-linear system, which multisecant SR1 factorizes for both channels
                ! together and the ARH family per channel, drops the far and the noisy
                ! entries
                n_systems = merge(1_ip, n_particle, arh_type == "ms_sr1")
                if (n_linear /= n_particle * n_list) then
                    write(stderr, *) test_name//" failed: The linear system does "// &
                        "not keep the entire history for "//arh_type//" and "// &
                        case_name//"."
                    passed = .false.
                end if
                if (i_case == 1 .and. n_nonlinear /= n_systems * n_list) then
                    write(stderr, *) test_name//" failed: The non-linear system "// &
                        "does not keep every independent entry for "//arh_type// &
                        " and "//case_name//"."
                    passed = .false.
                else if (i_case == 2 .and. n_nonlinear /= n_systems) then
                    write(stderr, *) test_name//" failed: The non-linear system "// &
                        "does not drop the history entry beyond the step-length "// &
                        "cutoff for "//arh_type//" and "//case_name//"."
                    passed = .false.
                else if (i_noisy > 0 .and. n_nonlinear >= n_systems * n_list) then
                    write(stderr, *) test_name//" failed: The non-linear system "// &
                        "does not drop a history entry whose response is noise for "// &
                        arh_type//" and "//case_name//"."
                    passed = .false.
                end if

                ! determine if every quantity follows the basis of its own system
                if (arh_type == "ms_sr1") then
                    if (size(arh_object%a_inv, 1) /= n_linear .or. &
                        size(arh_object%a_inv_comb, 1) /= n_nonlinear) then
                        write(stderr, *) test_name//" failed: Multisecant SR1 "// &
                            "pseudoinverses do not match the potential difference "// &
                            "directions of their system for "//case_name//"."
                        passed = .false.
                    end if
                else if (arh_type /= "ms_sp") then
                    if (size(arh_object%linear_potential_dirs, 2) /= n_linear .or. &
                        size(arh_object%nonlinear_potential_dirs, 2) /= n_nonlinear) &
                        then
                        write(stderr, *) test_name//" failed: Potential difference "// &
                            "directions are not rebased onto the basis of the "// &
                            "density matrix difference directions of their system "// &
                            "for "//arh_type//" and "//case_name//"."
                        passed = .false.
                    end if
                end if
                if (arh_type == "ms_sp" .or. arh_type == "ms_psb") then
                    if (size(arh_object%a_sym, 1) /= n_linear .or. &
                        size(arh_object%a_sym_nonlinear, 1) /= n_nonlinear) then
                        write(stderr, *) test_name//" failed: Symmetrized A "// &
                            "matrices do not match the density matrix difference "// &
                            "directions of their system for "//arh_type//" and "// &
                            case_name//"."
                        passed = .false.
                    end if
                end if

                ! deallocate ARH object
                deallocate(arh_object)
            end do
            deallocate(dm_list, current, lists)
        end do

        ! deallocate OAO object
        deallocate(oao_object)

    contains

        subroutine build_from_history()
            !
            ! this subroutine hands the history and the current potentials to the ARH
            ! object, marks its model stale and builds the model for the shell
            !
            ! history and current potentials the shell keeps
            arh_object%dm_list = dm_list
            if (n_particle == 1) then
                arh_object%evaluate_dm_cs => mock_evaluate_dm_cs
                arh_object%fock = current(:, :, :, 1)
                arh_object%fock_list = lists(:, :, :, :, 1)
            else
                arh_object%evaluate_dm_os => mock_evaluate_dm_os
                arh_object%v_same_spin = current(:, :, :, 1)
                arh_object%v_same_spin_list = lists(:, :, :, :, 1)
                arh_object%v_opposite_spin = current(:, :, :, 2)
                arh_object%v_opposite_spin_list = lists(:, :, :, :, 2)
            end if
            arh_object%v_nonlinear = current(:, :, :, n_pot)
            arh_object%v_nonlinear_list = lists(:, :, :, :, n_pot)
            arh_object%model_stale = .true.

            ! build the model
            mock_requests = [integer(ip) :: ]
            if (n_particle == 1) then
                call build_hess_model_cs(error)
            else
                call build_hess_model_os(error)
            end if

        end subroutine build_from_history

    end function check_build_hess_model

    logical(c_bool) function test_arh_factory_mo_cs() bind(C)
        !
        ! this function tests the subroutine which returns the modified ARH orbital
        ! updating function for the closed-shell case with the orbitals parameterized
        ! in the MO basis
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_factory, arh_object, arh_settings_type, &
                           evaluate_dm_cs_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff

        real(rp), target :: mo_coeff(n_ao, n_mo)
        real(rp) :: ao_overlap(n_ao, n_ao), mo_coeff_3d(n_ao, n_mo, 1)
        integer(ip) :: error
        type(arh_settings_type) :: settings
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        type(solver_settings_type) :: solver_settings

        ! setup settings object
        call setup_settings(settings)
        settings%arh_type = "ms_psb"

        ! initialize random orthonormal MO coefficients and callback function pointers
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff_3d = generate_random_mo_coeff(ao_overlap, n_mo, 1_ip)
        mo_coeff = mo_coeff_3d(:, :, 1)
        evaluate_dm_funptr => mock_evaluate_dm_cs

        ! call routine for a new calculation and determine if the ARH object is set up
        call arh_factory(mo_coeff, ao_overlap, n_occ(1), 1_ip, n_ao, n_mo, &
                         evaluate_dm_funptr, obj_func_arh_funptr, &
                         update_orbs_arh_funptr, solver_settings, error, settings)
        test_arh_factory_mo_cs = check_arh_factory( &
            1_ip, "test_arh_factory_mo_cs", .true., 1_ip, error, settings, &
            solver_settings, obj_func_arh_funptr, update_orbs_arh_funptr)

        ! the calls below read the ARH and MO objects
        if (test_arh_factory_mo_cs) then
            ! determine if the stored MO coefficients are those of the caller, which
            ! are rotated in place, although they are passed without the dimension of
            ! the particle channels
            mo_object%mo_coeff(1, 1, 1) = 42.0_rp
            if (abs(mo_coeff(1, 1) - 42.0_rp) > tol) then
                write(stderr, *) "test_arh_factory_mo_cs failed: MO coefficients "// &
                    "not associated with those of the caller."
                test_arh_factory_mo_cs = .false.
            end if

            ! leave the history and evaluation state of a previous calculation behind
            ! and call routine again for a new calculation
            arh_object%evaluation_stale = .false.
            allocate(arh_object%dm_list(1, 1, 1, 1))
            call arh_factory(mo_coeff, ao_overlap, n_occ(1), 1_ip, n_ao, n_mo, &
                             evaluate_dm_funptr, obj_func_arh_funptr, &
                             update_orbs_arh_funptr, solver_settings, error, settings)
            if (.not. check_arh_factory(2_ip, "test_arh_factory_mo_cs", .true., 1_ip, &
                                        error, settings, solver_settings, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr)) &
                test_arh_factory_mo_cs = .false.

            ! call routine with an unknown ARH type
            arh_object%orbitals%hess_eigen_stale = .false.
            settings%arh_type = "unknown"
            call arh_factory(mo_coeff, ao_overlap, n_occ(1), 1_ip, n_ao, n_mo, &
                             evaluate_dm_funptr, obj_func_arh_funptr, &
                             update_orbs_arh_funptr, solver_settings, error, settings)
            if (.not. check_arh_factory(3_ip, "test_arh_factory_mo_cs", .true., 1_ip, &
                                        error, settings, solver_settings, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr)) &
                test_arh_factory_mo_cs = .false.
        end if

        ! deallocate ARH and MO objects
        if (allocated(arh_object)) deallocate(arh_object)
        if (allocated(mo_object)) deallocate(mo_object)

    end function test_arh_factory_mo_cs

    logical(c_bool) function test_arh_factory_mo_os() bind(C)
        !
        ! this function tests the subroutine which returns the modified ARH orbital
        ! updating function for the open-shell case with the orbitals parameterized in
        ! the MO basis, whose spin channels are differently occupied
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_factory, arh_object, arh_settings_type, &
                           evaluate_dm_os_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: error
        type(arh_settings_type) :: settings
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        type(solver_settings_type) :: solver_settings

        ! setup settings object
        call setup_settings(settings)
        settings%arh_type = "ms_psb"

        ! initialize random orthonormal MO coefficients and callback function pointers
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        evaluate_dm_funptr => mock_evaluate_dm_os

        ! call routine for a new calculation and determine if the ARH object is set up
        call arh_factory(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                         evaluate_dm_funptr, obj_func_arh_funptr, &
                         update_orbs_arh_funptr, solver_settings, error, settings)
        test_arh_factory_mo_os = check_arh_factory( &
            1_ip, "test_arh_factory_mo_os", .true., n_particle, error, settings, &
            solver_settings, obj_func_arh_funptr, update_orbs_arh_funptr)

        ! the calls below read the ARH object
        if (test_arh_factory_mo_os) then
            ! leave the history and evaluation state of a previous calculation behind
            ! and call routine again for a new calculation
            arh_object%evaluation_stale = .false.
            allocate(arh_object%dm_list(1, 1, 1, 1))
            call arh_factory(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                             evaluate_dm_funptr, obj_func_arh_funptr, &
                             update_orbs_arh_funptr, solver_settings, error, settings)
            if (.not. check_arh_factory(2_ip, "test_arh_factory_mo_os", .true., &
                                        n_particle, error, settings, solver_settings, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr)) &
                test_arh_factory_mo_os = .false.

            ! call routine with an unknown ARH type
            arh_object%orbitals%hess_eigen_stale = .false.
            settings%arh_type = "unknown"
            call arh_factory(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                             evaluate_dm_funptr, obj_func_arh_funptr, &
                             update_orbs_arh_funptr, solver_settings, error, settings)
            if (.not. check_arh_factory(3_ip, "test_arh_factory_mo_os", .true., &
                                        n_particle, error, settings, solver_settings, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr)) &
                test_arh_factory_mo_os = .false.
        end if

        ! deallocate ARH and MO objects
        if (allocated(arh_object)) deallocate(arh_object)
        if (allocated(mo_object)) deallocate(mo_object)

    end function test_arh_factory_mo_os

    logical(c_bool) function test_arh_factory_oao_cs() bind(C)
        !
        ! this function tests the subroutine which returns the modified ARH orbital
        ! updating function for the closed-shell case with the orbitals parameterized
        ! in the OAO basis
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_factory, arh_object, arh_settings_type, &
                           evaluate_dm_cs_type
        use otr_common_test_reference, only: n_ao, n_occ
        use otr_oao, only: oao_object
        use otr_common_unit_tests, only: identity_matrix, generate_random_density_matrix
        use opentrustregion_unit_tests, only: setup_settings

        real(rp), target :: dm_ao(n_ao, n_ao)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: error
        type(arh_settings_type) :: settings
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        type(solver_settings_type) :: solver_settings

        ! setup settings object
        call setup_settings(settings)
        settings%arh_type = "arh"

        ! initialize density matrix, orthonormal AO basis and callback function pointers
        dm_ao = generate_random_density_matrix(n_ao, n_occ(1))
        ao_overlap = identity_matrix(n_ao)
        evaluate_dm_funptr => mock_evaluate_dm_cs

        ! call routine for a new calculation and determine if the ARH object is set up
        call arh_factory(dm_ao, ao_overlap, 1_ip, n_ao, evaluate_dm_funptr, &
                         obj_func_arh_funptr, update_orbs_arh_funptr, solver_settings, &
                         error, settings)
        test_arh_factory_oao_cs = check_arh_factory( &
            1_ip, "test_arh_factory_oao_cs", .false., 1_ip, error, settings, &
            solver_settings, obj_func_arh_funptr, update_orbs_arh_funptr)

        ! the calls below read the ARH object
        if (test_arh_factory_oao_cs) then
            ! leave the history and evaluation state of a previous calculation behind
            ! and call routine again for a new calculation
            arh_object%evaluation_stale = .false.
            allocate(arh_object%dm_list(1, 1, 1, 1))
            call arh_factory(dm_ao, ao_overlap, 1_ip, n_ao, evaluate_dm_funptr, &
                             obj_func_arh_funptr, update_orbs_arh_funptr, &
                             solver_settings, error, settings)
            if (.not. check_arh_factory(2_ip, "test_arh_factory_oao_cs", .false., &
                                        1_ip, error, settings, solver_settings, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr)) &
                test_arh_factory_oao_cs = .false.

            ! call routine with an unknown ARH type
            arh_object%orbitals%hess_eigen_stale = .false.
            settings%arh_type = "unknown"
            call arh_factory(dm_ao, ao_overlap, 1_ip, n_ao, evaluate_dm_funptr, &
                             obj_func_arh_funptr, update_orbs_arh_funptr, &
                             solver_settings, error, settings)
            if (.not. check_arh_factory(3_ip, "test_arh_factory_oao_cs", .false., &
                                        1_ip, error, settings, solver_settings, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr)) &
                test_arh_factory_oao_cs = .false.
        end if

        ! deallocate ARH and OAO objects
        if (allocated(arh_object)) deallocate(arh_object)
        if (allocated(oao_object)) deallocate(oao_object)

    end function test_arh_factory_oao_cs

    logical(c_bool) function test_arh_factory_oao_os() bind(C)
        !
        ! this function tests the subroutine which returns the modified ARH orbital
        ! updating function for the open-shell case with the orbitals parameterized in
        ! the OAO basis
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use otr_arh, only: arh_factory, arh_object, arh_settings_type, &
                           evaluate_dm_os_type
        use otr_common_test_reference, only: n_ao, n_particle, n_occ
        use otr_oao, only: oao_object
        use otr_common_unit_tests, only: identity_matrix, generate_random_density_matrix
        use opentrustregion_unit_tests, only: setup_settings

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: i, error
        type(arh_settings_type) :: settings
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        type(solver_settings_type) :: solver_settings

        ! setup settings object
        call setup_settings(settings)
        settings%arh_type = "arh"

        ! initialize density matrices, orthonormal AO basis and callback function
        ! pointers
        do i = 1, n_particle
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do
        ao_overlap = identity_matrix(n_ao)
        evaluate_dm_funptr => mock_evaluate_dm_os

        ! call routine for a new calculation and determine if the ARH object is set up
        call arh_factory(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                         obj_func_arh_funptr, update_orbs_arh_funptr, solver_settings, &
                         error, settings)
        test_arh_factory_oao_os = check_arh_factory( &
            1_ip, "test_arh_factory_oao_os", .false., n_particle, error, settings, &
            solver_settings, obj_func_arh_funptr, update_orbs_arh_funptr)

        ! the calls below read the ARH object
        if (test_arh_factory_oao_os) then
            ! leave the history and evaluation state of a previous calculation behind
            ! and call routine again for a new calculation
            arh_object%evaluation_stale = .false.
            allocate(arh_object%dm_list(1, 1, 1, 1))
            call arh_factory(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                             obj_func_arh_funptr, update_orbs_arh_funptr, &
                             solver_settings, error, settings)
            if (.not. check_arh_factory(2_ip, "test_arh_factory_oao_os", .false., &
                                        n_particle, error, settings, solver_settings, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr)) &
                test_arh_factory_oao_os = .false.

            ! call routine with an unknown ARH type
            arh_object%orbitals%hess_eigen_stale = .false.
            settings%arh_type = "unknown"
            call arh_factory(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                             obj_func_arh_funptr, update_orbs_arh_funptr, &
                             solver_settings, error, settings)
            if (.not. check_arh_factory(3_ip, "test_arh_factory_oao_os", .false., &
                                        n_particle, error, settings, solver_settings, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr)) &
                test_arh_factory_oao_os = .false.
        end if

        ! deallocate ARH and OAO objects
        if (allocated(arh_object)) deallocate(arh_object)
        if (allocated(oao_object)) deallocate(oao_object)

    end function test_arh_factory_oao_os

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
                write(stderr, *) "test_arh_sanity_check failed: Error thrown for "// &
                    trim(arh_types(i))//" ARH type."
                test_arh_sanity_check = .false.
            end if
        end do

        ! check if ARH type is converted to lowercase
        settings%arh_type = "MS_PSB"
        call arh_sanity_check(settings, error)
        if (settings%arh_type /= "ms_psb") then
            write(stderr, *) "test_arh_sanity_check failed: ARH type not converted "// &
                "to lowercase."
            test_arh_sanity_check = .false.
        end if

        ! check if unknown ARH type is rejected
        settings%arh_type = "unknown"
        call arh_sanity_check(settings, error)
        if (error == 0) then
            write(stderr, *) "test_arh_sanity_check failed: Error not thrown for "// &
                "unknown ARH type."
            test_arh_sanity_check = .false.
        end if

    end function test_arh_sanity_check

    logical(c_bool) function test_arh_set_solver_settings() bind(C)
        !
        ! this function tests the subroutine which wires the ARH preconditioners, the
        ! extra trial vectors and a given projection of the orbital basis into the
        ! solver settings, requests the Hessian refresh and reports the symmetry of the
        ! approximate Hessian
        !
        use otr_arh, only: arh_set_solver_settings, arh_type, arh_oao_type, &
                           precond_arh_callback_ptr, precond_pd_arh_callback_ptr, &
                           get_extra_trial_vectors_arh_callback_ptr, arh_n_micro
        use otr_oao, only: oao_object
        use opentrustregion, only: solver_settings_type, default_solver_settings
        use opentrustregion_unit_tests, only: mock_project
        use test_reference, only: ref_settings, assignment(=), operator(/=)

        type(solver_settings_type) :: solver_settings
        class(arh_type), allocatable :: arh
        integer(ip) :: i_case, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_arh_set_solver_settings = .true.

        ! wire the ARH routines into uninitialized settings, which have to be
        ! initialized to their defaults, together with a given projection, as the OAO
        ! basis does, and into settings initialized to the reference values with a
        ! projection supplied by the caller, which have to be kept, without a given
        ! projection, as the MO basis does
        allocate(oao_object)
        allocate(arh, source=arh_oao_type(oao_object))
        do i_case = 1, 2
            if (i_case == 1) then
                case_name = "for uninitialized settings"
                arh%settings%arh_type = "arh"
                solver_settings%initialized = .false.
                call arh_set_solver_settings(solver_settings, arh, error, mock_project)
            else
                case_name = "for initialized settings"
                arh%settings%arh_type = "ms_sr1"
                solver_settings = ref_settings
                solver_settings%project => mock_project
                solver_settings%stability_settings%project => mock_project
                call arh_set_solver_settings(solver_settings, arh, error)
            end if
            if (error /= 0) then
                write(stderr, *) "test_arh_set_solver_settings failed: Produced "// &
                    "error "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (.not. solver_settings%initialized) then
                write(stderr, *) "test_arh_set_solver_settings failed: Settings "// &
                    "not initialized "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (.not. associated(solver_settings%precond, precond_arh_callback_ptr)) &
                then
                write(stderr, *) "test_arh_set_solver_settings failed: ARH "// &
                    "preconditioner not wired "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (.not. associated(solver_settings%stability_settings%precond, &
                                 precond_arh_callback_ptr)) then
                write(stderr, *) "test_arh_set_solver_settings failed: ARH "// &
                    "preconditioner not wired into the stability check "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (.not. associated(solver_settings%precond_pd, &
                                 precond_pd_arh_callback_ptr)) then
                write(stderr, *) "test_arh_set_solver_settings failed: ARH "// &
                    "positive-definite preconditioner not wired "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (.not. associated(solver_settings%get_extra_trial_vectors, &
                                 get_extra_trial_vectors_arh_callback_ptr)) then
                write(stderr, *) "test_arh_set_solver_settings failed: ARH extra "// &
                    "trial vectors not wired "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (.not. associated( &
                solver_settings%stability_settings%get_extra_trial_vectors, &
                get_extra_trial_vectors_arh_callback_ptr)) then
                write(stderr, *) "test_arh_set_solver_settings failed: ARH extra "// &
                    "trial vectors not wired into the stability check "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (.not. (associated(solver_settings%project, mock_project) .and. &
                       associated(solver_settings%stability_settings%project, &
                                  mock_project))) then
                write(stderr, *) "test_arh_set_solver_settings failed: Given "// &
                    "projection not wired or projection supplied by the caller not "// &
                    "kept "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (.not. solver_settings%refresh_hess) then
                write(stderr, *) "test_arh_set_solver_settings failed: Hessian "// &
                    "refresh not requested "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (solver_settings%n_micro /= arh_n_micro) then
                write(stderr, *) "test_arh_set_solver_settings failed: Micro "// &
                    "iteration limit not raised "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if
            if (solver_settings%hess_symm .neqv. &
                (trim(arh%settings%arh_type) /= "arh")) then
                write(stderr, *) "test_arh_set_solver_settings failed: Symmetry of "// &
                    "the approximate Hessian of ARH type "// &
                    trim(arh%settings%arh_type)//" reported wrongly "//case_name//"."
                test_arh_set_solver_settings = .false.
            end if

            ! apart from the Hessian refresh, the micro iteration limit and the symmetry
            ! of the approximate Hessian, uninitialized settings have to be set to their
            ! defaults and initialized settings have to be kept
            if (i_case == 1) then
                solver_settings%refresh_hess = default_solver_settings%refresh_hess
                solver_settings%n_micro = default_solver_settings%n_micro
                solver_settings%hess_symm = default_solver_settings%hess_symm
                if (solver_settings /= default_solver_settings) then
                    write(stderr, *) "test_arh_set_solver_settings failed: "// &
                        "Settings not set to their defaults "//case_name//"."
                    test_arh_set_solver_settings = .false.
                end if
            else
                solver_settings%refresh_hess = ref_settings%refresh_hess
                solver_settings%n_micro = ref_settings%n_micro
                solver_settings%hess_symm = ref_settings%hess_symm
                if (solver_settings /= ref_settings) then
                    write(stderr, *) "test_arh_set_solver_settings failed: "// &
                        "Settings not kept "//case_name//"."
                    test_arh_set_solver_settings = .false.
                end if
            end if
        end do
        deallocate(arh, oao_object)

    end function test_arh_set_solver_settings

    logical(c_bool) function test_obj_func_arh_cs_callback() bind(C)
        !
        ! this function tests the function which defines the energy evaluation for the
        ! closed-shell case, which also adds the evaluated point to the history
        !
        test_obj_func_arh_cs_callback = &
            check_obj_func_arh(1_ip, "test_obj_func_arh_cs_callback")

    end function test_obj_func_arh_cs_callback

    logical(c_bool) function test_obj_func_arh_os_callback() bind(C)
        !
        ! this function tests the function which defines the energy evaluation for the
        ! open-shell case, which also adds the evaluated point to the history
        !
        use otr_common_test_reference, only: n_particle

        test_obj_func_arh_os_callback = &
            check_obj_func_arh(n_particle, "test_obj_func_arh_os_callback")

    end function test_obj_func_arh_os_callback

    logical(c_bool) function test_update_orbs_arh_cs_callback() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation for the closed-shell case
        !
        test_update_orbs_arh_cs_callback = &
            check_update_orbs_arh(1_ip, "test_update_orbs_arh_cs_callback")

    end function test_update_orbs_arh_cs_callback

    logical(c_bool) function test_update_orbs_arh_os_callback() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation for the open-shell case
        !
        use otr_common_test_reference, only: n_particle

        test_update_orbs_arh_os_callback = &
            check_update_orbs_arh(n_particle, "test_update_orbs_arh_os_callback")

    end function test_update_orbs_arh_os_callback

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
        ! last assembled
        !
        use otr_arh, only: rebuild_stale_hess_model, arh_object, arh_oao_type
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_unit_tests, only: mock_requests, generate_random_density_matrix

        integer(ip) :: n_particle_case, n_param_case
        real(rp) :: dm_oao(n_ao, n_ao, n_particle)
        integer(ip) :: i, i_case, error
        logical :: stale
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_rebuild_stale_hess_model = .true.

        ! a stale model on an empty history has to be rebuilt for the closed- and the
        ! open-shell case, which discards its low-rank part, while a model that is not
        ! stale has to be left untouched
        do i_case = 1, 3
            if (i_case == 1) then
                case_name = "for stale closed-shell model"
                n_particle_case = 1
            else if (i_case == 2) then
                case_name = "for stale open-shell model"
                n_particle_case = 2
            else
                case_name = "for model which is not stale"
                n_particle_case = 2
            end if
            stale = i_case /= 3
            n_param_case = n_particle_case * n_ao * (n_ao - 1) / 2
            do i = 1, n_particle_case
                dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
            end do
            allocate(oao_object)
            oao_object%n_ao = n_ao
            oao_object%n_particle = n_particle_case
            oao_object%n_param = n_param_case
            oao_object%dm_oao = dm_oao(:, :, :n_particle_case)
            allocate(arh_object, source=arh_oao_type(oao_object))
            call setup_settings(arh_object%settings)
            arh_object%settings%arh_type = "ms_sr1"
            call setup_empty_history(n_particle_case)
            allocate(arh_object%coupling_matrix(1, 1))
            arh_object%model_stale = stale
            mock_requests = [integer(ip) :: ]
            call rebuild_stale_hess_model(error)
            if (error /= 0) then
                write(stderr, *) "test_rebuild_stale_hess_model failed: Produced "// &
                    "error "//case_name//"."
                test_rebuild_stale_hess_model = .false.
            end if
            if (arh_object%model_stale) then
                write(stderr, *) "test_rebuild_stale_hess_model failed: Hessian "// &
                    "model still marked stale "//case_name//"."
                test_rebuild_stale_hess_model = .false.
            end if
            if (stale .eqv. allocated(arh_object%coupling_matrix)) then
                write(stderr, *) "test_rebuild_stale_hess_model failed: Hessian "// &
                    "model rebuilt wrongly "//case_name//"."
                test_rebuild_stale_hess_model = .false.
            end if
            if (size(mock_requests) /= 0) then
                write(stderr, *) "test_rebuild_stale_hess_model failed: Density "// &
                    "matrix evaluating function called "//case_name//"."
                test_rebuild_stale_hess_model = .false.
            end if
            deallocate(arh_object, oao_object)
        end do

    end function test_rebuild_stale_hess_model

    logical(c_bool) function test_hess_x_arh_callback() bind(C)
        !
        ! this function tests the subroutine which defines the Hessian linear
        ! transformation on the basis of augmented Roothaan-Hall and related methods,
        ! which adds the low-rank part to the static part of the orbital basis
        !
        use otr_arh, only: hess_x_arh_callback, arh_object
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, &
                                     setup_random_mo_channels, ref_hess_x_static_mo

        integer(ip), parameter :: n_dirs = 2

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        real(rp), allocatable :: x(:), hess_x(:), expected_hess_x(:)
        integer(ip) :: n_param, i_case, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_hess_x_arh_callback = .true.

        ! set up the ARH object in the MO basis with random Fock matrix blocks of the
        ! differently occupied spin channels
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_arh_and_mo_objects(mo_coeff, ao_overlap, n_occ)
        call setup_random_mo_channels(n_occ, n_mo)
        n_param = arh_object%orbitals%n_param
        allocate(x(n_param), hess_x(n_param))
        call random_number(x)

        ! call routine for a low-rank part with random, and hence non-symmetric,
        ! coupling matrix, for a stale model on an empty history, which has to be
        ! rebuilt before it is applied and thereby loses the low-rank part, and for the
        ! model without low-rank part which this leaves behind
        do i_case = 1, 3
            expected_hess_x = ref_hess_x_static_mo(x, mo_object%mo_channels)
            if (i_case == 1) then
                case_name = "with low-rank part"
                allocate(arh_object%expansion_dirs(n_param, n_dirs), &
                         arh_object%projection_dirs(n_param, n_dirs), &
                         arh_object%coupling_matrix(n_dirs, n_dirs))
                call random_number(arh_object%expansion_dirs)
                call random_number(arh_object%projection_dirs)
                call random_number(arh_object%coupling_matrix)
                expected_hess_x = expected_hess_x + matmul( &
                    arh_object%expansion_dirs, &
                    matmul(arh_object%coupling_matrix, &
                           matmul(transpose(arh_object%projection_dirs), x)))
            else if (i_case == 2) then
                case_name = "for stale model"
                call setup_empty_history(n_particle)
                arh_object%model_stale = .true.
            else
                case_name = "without low-rank part"
            end if
            call hess_x_arh_callback(x, hess_x, error)
            if (error /= 0) then
                write(stderr, *) "test_hess_x_arh_callback failed: Produced error "// &
                    case_name//"."
                test_hess_x_arh_callback = .false.
            end if
            if (norm2(hess_x - expected_hess_x) > tol) then
                write(stderr, *) "test_hess_x_arh_callback failed: Incorrect "// &
                    "Hessian linear transformation "//case_name//"."
                test_hess_x_arh_callback = .false.
            end if
            if (arh_object%model_stale) then
                write(stderr, *) "test_hess_x_arh_callback failed: Hessian model "// &
                    "not rebuilt "//case_name//"."
                test_hess_x_arh_callback = .false.
            end if
        end do

        ! deallocate ARH and MO objects
        deallocate(arh_object, mo_object)

    end function test_hess_x_arh_callback

    logical(c_bool) function test_inv_hess_x_arh() bind(C)
        !
        ! this function tests the exact, optionally level-shifted inverse of the
        ! approximate Hessian in both the closed- and the open-shell case and for both
        ! orbital bases
        !
        use otr_common_test_reference, only: n_particle_ref => n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, &
                                     setup_random_mo_channels, ref_hess_x_static_mo

        integer(ip) :: i_basis, n_particle

        ! assume tests pass
        test_inv_hess_x_arh = .true.

        ! the inverse only depends on the low-rank factors and the orbital basis, not
        ! on the ARH type which assembled them, so every basis and shell is checked
        ! once
        do i_basis = 1, 2
            do n_particle = 1, n_particle_ref
                if (.not. check_inv_hess_x_arh_case(i_basis == 1)) &
                    test_inv_hess_x_arh = .false.
            end do
        end do

    contains

        subroutine check_inv_hess_x_arh(hess, vector, mu, case_name, passed)
            !
            ! this subroutine checks the exact-inverse of the ARH Woodbury
            ! preconditioner against a reference Hessian by applying the reference
            ! (B - mu*I) to the preconditioned vector which has to recover the vector,
            ! both with and without the level shift
            !
            use otr_arh, only: inv_hess_x_arh

            real(rp), intent(in) :: hess(:, :), vector(:), mu
            character(len=*), intent(in) :: case_name
            logical, intent(inout) :: passed

            real(rp), allocatable :: actual(:)
            integer(ip) :: error

            allocate(actual(size(vector)))

            ! shifted entry point: applying the reference (B - mu*I) to the
            ! preconditioned vector has to recover the vector
            call inv_hess_x_arh(vector, actual, error, mu)
            if (error /= 0) then
                write(stderr, *) "test_inv_hess_x_arh failed: Shifted inverse "// &
                    "produced error for the "//case_name//" case."
                passed = .false.
            end if
            if (norm2(matmul(hess, actual) - mu * actual - vector) > &
                tol * (1.0_rp + norm2(hess) * norm2(actual))) then
                write(stderr, *) "test_inv_hess_x_arh failed: Shifted inverse is "// &
                    "not inverse of shifted reference Hessian for the "//case_name// &
                    " case."
                passed = .false.
            end if

            ! the unshifted entry point has to reproduce the same round trip against the
            ! reference Hessian itself, without a level shift
            call inv_hess_x_arh(vector, actual, error)
            if (error /= 0) then
                write(stderr, *) "test_inv_hess_x_arh failed: Unshifted inverse "// &
                    "produced error for the "//case_name//" case."
                passed = .false.
            end if
            if (norm2(matmul(hess, actual) - vector) > &
                tol * (1.0_rp + norm2(hess) * norm2(actual))) then
                write(stderr, *) "test_inv_hess_x_arh failed: Inverse is not "// &
                    "inverse of reference Hessian for the "//case_name//" case."
                passed = .false.
            end if

        end subroutine check_inv_hess_x_arh

        function check_inv_hess_x_arh_case(mo_basis) result(passed)
            !
            ! this function performs the exact-inverse check of the ARH Woodbury
            ! preconditioner for one orbital basis and the current shell with random
            ! expansion and projection directions and a random coupling matrix, and for
            ! a stale model on an empty history, which has to be rebuilt without
            ! low-rank part before it is inverted
            !
            use otr_arh, only: arh_object
            use otr_mo, only: mo_object
            use otr_oao, only: oao_object
            use otr_mo_test_reference, only: n_mo
            use otr_common_test_reference, only: n_ao, n_occ
            use otr_common_unit_tests, only: identity_matrix, shell_names
            use otr_oao_unit_tests, only: ref_unpack_asymm, ref_hess_x_oao

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
                call setup_arh_and_mo_objects(mo_coeff, ao_overlap, n_occ(:n_particle))
                call setup_random_mo_channels(n_occ(:n_particle), n_mo)
                do k = 1, n_particle
                    associate (channel => mo_object%mo_channels(k))
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
                call setup_arh_and_oao_objects("ms_sr1", dm_oao, fock_oo, fock_vv, &
                                               n_ao, n_particle, &
                                               n_particle * n_ao * (n_ao - 1) / 2)

                ! cache a wrong eigendecomposition, which has to be refreshed since the
                ! static part is marked as changed
                allocate(oao_object%hess_eigvecs(n_ao, n_ao, n_particle), &
                         oao_object%hess_eigvals(n_ao, n_particle))
                oao_object%hess_eigvecs = 0.0_rp
                oao_object%hess_eigvals = 1.0_rp
            end if
            n_param = arh_object%orbitals%n_param

            ! build the dense static part
            allocate(static_hess(n_param, n_param), e_i(n_param))
            zero_response = 0.0_rp
            do i = 1, n_param
                e_i = 0.0_rp
                e_i(i) = 1.0_rp
                if (mo_basis) then
                    static_hess(:, i) = ref_hess_x_static_mo(e_i, mo_object%mo_channels)
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
            arh_object%expansion_dirs = dirs(:, :n_dirs)
            arh_object%projection_dirs = dirs(:, n_dirs + 1:2 * n_dirs)

            ! random, and hence non-symmetric, coupling matrix, so that a transposed
            ! coupling matrix would not go unnoticed
            call random_number(coupling)
            arh_object%coupling_matrix = low_rank_scale * (coupling - 0.5_rp)

            ! check inverse against the dense approximate Hessian
            hess = static_hess + matmul(arh_object%expansion_dirs, matmul( &
                arh_object%coupling_matrix, transpose(arh_object%projection_dirs)))
            call check_inv_hess_x_arh(hess, dirs(:, 2 * n_dirs + 1), mu, case_name, &
                                      passed)

            ! mark the model stale on an empty history and determine if it is rebuilt
            ! before it is inverted, which discards the low-rank part and leaves the
            ! inverse of the static part
            call setup_empty_history(n_particle)
            arh_object%model_stale = .true.
            call check_inv_hess_x_arh(static_hess, dirs(:, 2 * n_dirs + 1), mu, &
                                      "stale-model "//case_name, passed)
            if (arh_object%model_stale .or. allocated(arh_object%coupling_matrix)) then
                write(stderr, *) "test_inv_hess_x_arh failed: Stale Hessian model "// &
                    "not rebuilt before it is inverted for the "//case_name//" case."
                passed = .false.
            end if

            ! deallocate ARH and orbital objects
            deallocate(arh_object)
            if (allocated(mo_object)) deallocate(mo_object)
            if (allocated(oao_object)) deallocate(oao_object)

        end function check_inv_hess_x_arh_case

    end function test_inv_hess_x_arh

    logical(c_bool) function test_precond_arh_callback() bind(C)
        !
        ! this function tests the subroutine which defines the level-shifted
        ! preconditioner of the ARH approximate Hessian, which applies the
        ! level-shifted inverse of the approximate Hessian
        !
        use otr_arh, only: precond_arh_callback, arh_object
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, &
                                     setup_identity_mo_eigenbasis

        real(rp), parameter :: mu = 0.3_rp

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        real(rp), allocatable :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! assume tests pass
        test_precond_arh_callback = .true.

        ! set up the ARH object in the MO basis without a low-rank part and with an
        ! eigendecomposition whose eigenvalue pairs, the open-shell differences 2 (e_v
        ! - e_o), shifted by the level shift are all 2, so that the level-shifted
        ! inverse of the approximate Hessian halves the residual
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_arh_and_mo_objects(mo_coeff, ao_overlap, n_occ)
        call setup_identity_mo_eigenbasis(0.0_rp, (2.0_rp + mu) / 2.0_rp)
        allocate(residual(mo_object%n_param), precond_residual(mo_object%n_param))
        call random_number(residual)

        ! call routine and determine if the residual is halved
        call precond_arh_callback(residual, mu, precond_residual, error)
        if (error /= 0) then
            write(stderr, *) "test_precond_arh_callback failed: Produced error."
            test_precond_arh_callback = .false.
        end if
        if (norm2(precond_residual - 0.5_rp * residual) > tol) then
            write(stderr, *) "test_precond_arh_callback failed: Incorrect "// &
                "preconditioned residual."
            test_precond_arh_callback = .false.
        end if

        ! deallocate ARH and MO objects
        deallocate(arh_object, mo_object)

    end function test_precond_arh_callback

    logical(c_bool) function test_precond_pd_arh_callback() bind(C)
        !
        ! this function tests the subroutine which defines the positive-definite
        ! preconditioner of the ARH approximate Hessian, which applies the one of the
        ! orbital basis to the orbital object
        !
        use otr_arh, only: precond_pd_arh_callback, arh_object
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, &
                                     setup_identity_mo_eigenbasis

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        real(rp), allocatable :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! assume tests pass
        test_precond_pd_arh_callback = .true.

        ! set up the ARH object in the MO basis with an eigendecomposition whose
        ! eigenvalue pairs, the open-shell differences 2 (e_v - e_o), are all -2, so
        ! that the positive-definite preconditioner of the orbital basis, which divides
        ! by their magnitudes, halves the residual, unlike the level-shifted one
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_arh_and_mo_objects(mo_coeff, ao_overlap, n_occ)
        call setup_identity_mo_eigenbasis(1.0_rp, 0.0_rp)
        allocate(residual(mo_object%n_param), precond_residual(mo_object%n_param))
        call random_number(residual)

        ! call routine and determine if the residual is halved
        call precond_pd_arh_callback(residual, precond_residual, error)
        if (error /= 0) then
            write(stderr, *) "test_precond_pd_arh_callback failed: Produced error."
            test_precond_pd_arh_callback = .false.
        end if
        if (norm2(precond_residual - 0.5_rp * residual) > tol) then
            write(stderr, *) "test_precond_pd_arh_callback failed: Incorrect "// &
                "preconditioned residual."
            test_precond_pd_arh_callback = .false.
        end if

        ! deallocate ARH and MO objects
        deallocate(arh_object, mo_object)

    end function test_precond_pd_arh_callback

    logical(c_bool) function test_get_extra_trial_vectors_arh_callback() bind(C)
        !
        ! this function tests the subroutine which returns the extra trial vectors of
        ! the orbital basis for the solver's initial trial space
        !
        use otr_arh, only: get_extra_trial_vectors_arh_callback, arh_object
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, &
                                     setup_identity_mo_eigenbasis

        integer(ip), parameter :: n_extra = 2

        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        real(rp), allocatable :: trial_vectors(:, :), expected(:, :)
        integer(ip) :: error

        ! assume tests pass
        test_get_extra_trial_vectors_arh_callback = .true.

        ! set up the ARH object in the MO basis with an eigendecomposition in the MO
        ! basis itself whose only negative eigenvalue pair belongs to the first
        ! occupied and the first virtual orbital of the first particle channel, at
        ! packed index 1, so that the first extra trial vector is the unit vector along
        ! it and the second vanishes
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_arh_and_mo_objects(mo_coeff, ao_overlap, n_occ)
        call setup_identity_mo_eigenbasis(0.0_rp, 1.0_rp)
        mo_object%mo_channels(1)%occ_eigvals(1) = 1.5_rp
        mo_object%mo_channels(1)%virt_eigvals(2) = 2.0_rp
        allocate(trial_vectors(mo_object%n_param, n_extra), &
                 expected(mo_object%n_param, n_extra))
        expected = 0.0_rp
        expected(1, 1) = 1.0_rp

        ! call routine and determine if the extra trial vectors of the orbital basis
        ! are returned
        call get_extra_trial_vectors_arh_callback(trial_vectors, error)
        if (error /= 0) then
            write(stderr, *) "test_get_extra_trial_vectors_arh_callback failed: "// &
                "Produced error."
            test_get_extra_trial_vectors_arh_callback = .false.
        end if
        if (norm2(trial_vectors - expected) > tol) then
            write(stderr, *) "test_get_extra_trial_vectors_arh_callback failed: "// &
                "Incorrect trial vectors."
            test_get_extra_trial_vectors_arh_callback = .false.
        end if
        deallocate(arh_object, mo_object)

    end function test_get_extra_trial_vectors_arh_callback

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
            write(stderr, *) "test_init_arh_settings failed: Function raised error."
            test_init_arh_settings = .false.
        end if

        ! check settings
        if (.not. (settings == default_settings)) then
            write(stderr, *) "test_init_arh_settings failed: Settings not "// &
                "initialized correctly."
            test_init_arh_settings = .false.
        end if

    end function test_init_arh_settings

    logical(c_bool) function test_arh_deconstructor() bind(C)
        !
        ! this function tests the subroutine which deallocates the ARH objects
        !
        use otr_arh, only: arh_deconstructor, arh_object, arh_oao_type
        use otr_mo, only: mo_object
        use otr_oao, only: oao_object

        ! assume tests pass
        test_arh_deconstructor = .true.

        ! allocate ARH, MO and OAO objects
        if (.not. allocated(arh_object)) allocate(arh_oao_type :: arh_object)
        if (.not. allocated(mo_object)) allocate(mo_object)
        if (.not. allocated(oao_object)) allocate(oao_object)

        ! call routine and determine if all objects are deallocated
        call arh_deconstructor()
        if (allocated(arh_object)) then
            write(stderr, *) "test_arh_deconstructor failed: ARH object not "// &
                "deallocated."
            test_arh_deconstructor = .false.
        end if
        if (allocated(mo_object)) then
            write(stderr, *) "test_arh_deconstructor failed: MO object not deallocated."
            test_arh_deconstructor = .false.
        end if
        if (allocated(oao_object)) then
            write(stderr, *) "test_arh_deconstructor failed: OAO object not "// &
                "deallocated."
            test_arh_deconstructor = .false.
        end if

        ! call routine again, which has to handle the already deallocated objects
        call arh_deconstructor()

    end function test_arh_deconstructor

    logical(c_bool) function test_construct_arh_mo() bind(C)
        !
        ! this function tests the function which returns the ARH object for the MO
        ! basis pointing to an MO object and its quantities
        !
        use otr_arh, only: arh_mo_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object

        type(arh_mo_type) :: arh
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)

        ! assume tests pass
        test_construct_arh_mo = .true.

        ! set up the MO object
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_mo_object(mo_coeff, ao_overlap, n_occ)

        ! call routine and determine if the quantities of the MO object are associated
        arh = arh_mo_type(mo_object)
        if (.not. associated(arh%orbitals, mo_object)) then
            write(stderr, *) "test_construct_arh_mo failed: Orbitals not associated."
            test_construct_arh_mo = .false.
        end if
        if (.not. associated(arh%mo_coeff, mo_coeff)) then
            write(stderr, *) "test_construct_arh_mo failed: MO coefficients not "// &
                "associated."
            test_construct_arh_mo = .false.
        end if
        if (.not. associated(arh%ao_overlap, mo_object%ao_overlap)) then
            write(stderr, *) "test_construct_arh_mo failed: AO overlap matrix not "// &
                "associated."
            test_construct_arh_mo = .false.
        end if
        if (.not. associated(arh%mo_channels, mo_object%mo_channels)) then
            write(stderr, *) "test_construct_arh_mo failed: Particle channels not "// &
                "associated."
            test_construct_arh_mo = .false.
        end if

        ! call routine for an MO object without MO coefficients, overlap matrix and
        ! particle channels and determine if only the orbitals are associated
        deallocate(mo_object)
        allocate(mo_object)
        arh = arh_mo_type(mo_object)
        if (.not. associated(arh%orbitals, mo_object)) then
            write(stderr, *) "test_construct_arh_mo failed: Orbitals not "// &
                "associated for an empty MO object."
            test_construct_arh_mo = .false.
        end if
        if (associated(arh%mo_coeff) .or. associated(arh%ao_overlap) .or. &
            associated(arh%mo_channels)) then
            write(stderr, *) "test_construct_arh_mo failed: Missing quantities "// &
                "associated for an empty MO object."
            test_construct_arh_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object)

    end function test_construct_arh_mo

    logical(c_bool) function test_construct_arh_oao() bind(C)
        !
        ! this function tests the function which returns the ARH object for the OAO
        ! basis pointing to an OAO object and its quantities
        !
        use otr_arh, only: arh_oao_type
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao, n_particle

        type(arh_oao_type) :: arh

        ! assume tests pass
        test_construct_arh_oao = .true.

        ! set up an OAO object whose basis-specific quantities are allocated except for
        ! the Fock matrix blocks
        allocate(oao_object)
        allocate(oao_object%s_inv_sqrt(n_ao, n_ao), &
                 oao_object%dm_oao(n_ao, n_ao, n_particle))

        ! call routine and determine if the allocated quantities are associated and the
        ! others are not
        arh = arh_oao_type(oao_object)
        if (.not. associated(arh%orbitals, oao_object)) then
            write(stderr, *) "test_construct_arh_oao failed: Orbitals not associated."
            test_construct_arh_oao = .false.
        end if
        if (.not. associated(arh%s_inv_sqrt, oao_object%s_inv_sqrt)) then
            write(stderr, *) "test_construct_arh_oao failed: Inverse square root "// &
                "of the overlap matrix not associated."
            test_construct_arh_oao = .false.
        end if
        if (.not. associated(arh%dm_oao, oao_object%dm_oao)) then
            write(stderr, *) "test_construct_arh_oao failed: Density matrix in the "// &
                "OAO basis not associated."
            test_construct_arh_oao = .false.
        end if
        if (associated(arh%fock_oo) .or. associated(arh%fock_vv)) then
            write(stderr, *) "test_construct_arh_oao failed: Unallocated Fock "// &
                "matrix blocks associated."
            test_construct_arh_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_construct_arh_oao

    logical(c_bool) function test_rotate_trial_arh_mo() bind(C)
        !
        ! this function tests the subroutine which returns the density matrix rotated
        ! by an orbital rotation in the AO basis and in the history basis, S D S, for
        ! the MO basis without moving the current orbitals
        !
        use otr_arh, only: arh_mo_type
        use otr_mo, only: mo_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_mo_test_reference, only: n_mo, n_param_os
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object, &
                                     ref_rotate_mo_coeff

        type(arh_mo_type) :: arh
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), mo_coeff_before(n_ao, n_mo, n_particle), &
                    expected(n_ao, n_mo, n_particle), kappa(n_param_os), &
                    rot_dm_ao(n_ao, n_ao, n_particle), &
                    rot_dm_hist(n_ao, n_ao, n_particle), expected_dm(n_ao, n_ao)
        integer(ip) :: error, k

        ! assume tests pass
        test_rotate_trial_arh_mo = .true.

        ! set up the MO object with random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        mo_coeff_before = mo_coeff
        call setup_mo_object(mo_coeff, ao_overlap, n_occ)
        arh = arh_mo_type(mo_object)
        call setup_settings(arh%settings)

        ! call routine and determine if the density matrix of the independently rotated
        ! orbitals is returned in the AO basis and multiplied by the AO overlap matrix
        ! from both sides, while the current orbitals are left untouched
        call random_number(kappa)
        kappa = 0.2_rp * (kappa - 0.5_rp)
        expected = ref_rotate_mo_coeff(kappa, mo_coeff, n_occ)
        call arh%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_trial_arh_mo failed: Produced error."
            test_rotate_trial_arh_mo = .false.
        end if
        do k = 1, 2
            expected_dm = matmul(expected(:, :n_occ(k), k), &
                                 transpose(expected(:, :n_occ(k), k)))
            if (norm2(rot_dm_ao(:, :, k) - expected_dm) > tol) then
                write(stderr, *) "test_rotate_trial_arh_mo failed: Incorrect "// &
                    "rotated density matrix in the AO basis."
                test_rotate_trial_arh_mo = .false.
            end if
            if (norm2(rot_dm_hist(:, :, k) - &
                      matmul(ao_overlap, matmul(expected_dm, ao_overlap))) > tol) then
                write(stderr, *) "test_rotate_trial_arh_mo failed: Incorrect "// &
                    "rotated density matrix in the history basis."
                test_rotate_trial_arh_mo = .false.
            end if
        end do
        if (norm2(mo_coeff - mo_coeff_before) > tol) then
            write(stderr, *) "test_rotate_trial_arh_mo failed: Current MO "// &
                "coefficients changed."
            test_rotate_trial_arh_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object)

    end function test_rotate_trial_arh_mo

    logical(c_bool) function test_rotate_trial_arh_oao() bind(C)
        !
        ! this function tests the subroutine which returns the density matrix rotated
        ! by an orbital rotation in the AO and OAO basis without moving the current
        ! orbitals
        !
        use otr_arh, only: arh_oao_type
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao, n_occ
        use otr_common_unit_tests, only: identity_matrix, generate_random_density_matrix
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_param_oao = n_ao * (n_ao - 1) / 2
        real(rp), parameter :: basis_scale = 1.5_rp, angle = 0.3_rp

        type(arh_oao_type) :: arh
        real(rp), target :: dm_ao(n_ao, n_ao, 1)
        real(rp) :: dm_oao(n_ao, n_ao, 1), kappa(n_param_oao), &
                    rot_dm_ao(n_ao, n_ao, 1), rot_dm_hist(n_ao, n_ao, 1), &
                    rotation(n_ao, n_ao)
        integer(ip) :: error

        ! assume tests pass
        test_rotate_trial_arh_oao = .true.

        ! set up the OAO object with an AO basis whose overlap has the inverse square
        ! root c I, so that the density matrix in the AO basis is c^2 times the one in
        ! the OAO basis
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = 1
        oao_object%n_param = n_param_oao
        oao_object%s_inv_sqrt = basis_scale * identity_matrix(n_ao)
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_occ(1))
        dm_ao = basis_scale**2 * dm_oao
        oao_object%dm_ao => dm_ao
        oao_object%dm_oao = dm_oao
        arh = arh_oao_type(oao_object)
        call setup_settings(arh%settings)

        ! call routine and determine if the density matrix is rotated consistently in
        ! both bases while the current density matrix is left untouched
        call random_number(kappa)
        kappa = 0.2_rp * kappa
        call arh%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_trial_arh_oao failed: Produced error."
            test_rotate_trial_arh_oao = .false.
        end if
        if (norm2(rot_dm_hist - dm_oao) < tol) then
            write(stderr, *) "test_rotate_trial_arh_oao failed: Density matrix not "// &
                "rotated."
            test_rotate_trial_arh_oao = .false.
        end if
        if (norm2(rot_dm_ao - basis_scale**2 * rot_dm_hist) > tol) then
            write(stderr, *) "test_rotate_trial_arh_oao failed: Density matrix not "// &
                "rotated consistently in the AO and OAO basis."
            test_rotate_trial_arh_oao = .false.
        end if
        if (norm2(matmul(rot_dm_hist(:, :, 1), rot_dm_hist(:, :, 1)) - &
                  rot_dm_hist(:, :, 1)) > tol) then
            write(stderr, *) "test_rotate_trial_arh_oao failed: Rotated density "// &
                "matrix not idempotent."
            test_rotate_trial_arh_oao = .false.
        end if
        if (norm2(oao_object%dm_oao - dm_oao) > tol) then
            write(stderr, *) "test_rotate_trial_arh_oao failed: Current density "// &
                "matrix in the OAO basis moved."
            test_rotate_trial_arh_oao = .false.
        end if
        if (norm2(dm_ao - basis_scale**2 * dm_oao) > tol) then
            write(stderr, *) "test_rotate_trial_arh_oao failed: Current density "// &
                "matrix in the AO basis moved."
            test_rotate_trial_arh_oao = .false.
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
        call arh%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_trial_arh_oao failed: Produced error for "// &
                "a plane rotation."
            test_rotate_trial_arh_oao = .false.
        end if
        if (norm2(rot_dm_hist(:, :, 1) - matmul( &
            transpose(rotation), matmul(dm_oao(:, :, 1), rotation))) > tol) then
            write(stderr, *) "test_rotate_trial_arh_oao failed: Incorrect density "// &
                "matrix for a plane rotation."
            test_rotate_trial_arh_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_rotate_trial_arh_oao

    logical(c_bool) function test_history_dm_arh_mo() bind(C)
        !
        ! this function tests the function which returns the current density matrix in
        ! the history basis for the MO basis
        !
        use otr_arh, only: arh_mo_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object

        type(arh_mo_type) :: arh
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), dm(n_ao, n_ao, n_particle)
        integer(ip) :: k

        ! assume tests pass
        test_history_dm_arh_mo = .true.

        ! set up the MO object with random orthonormal MO coefficients
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_mo_object(mo_coeff, ao_overlap, n_occ)
        arh = arh_mo_type(mo_object)

        ! call routine and determine if the density matrix multiplied by the AO overlap
        ! matrix from both sides is returned
        dm = arh%history_dm()
        do k = 1, 2
            if (norm2(dm(:, :, k) - matmul( &
                ao_overlap, matmul(mo_object%dm_ao(:, :, k), ao_overlap))) > tol) then
                write(stderr, *) "test_history_dm_arh_mo failed: Incorrect density "// &
                    "matrix."
                test_history_dm_arh_mo = .false.
            end if
        end do

        ! deallocate MO object
        deallocate(mo_object)

    end function test_history_dm_arh_mo

    logical(c_bool) function test_history_dm_arh_oao() bind(C)
        !
        ! this function tests the function which returns the current density matrix in
        ! the history basis for the OAO basis
        !
        use otr_arh, only: arh_oao_type
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao, n_particle

        type(arh_oao_type) :: arh
        real(rp) :: dm_oao(n_ao, n_ao, n_particle)

        ! assume tests pass
        test_history_dm_arh_oao = .true.

        ! set up the OAO object with a random density matrix
        call random_number(dm_oao)
        allocate(oao_object)
        oao_object%dm_oao = dm_oao
        arh = arh_oao_type(oao_object)

        ! call routine and determine if the density matrix in the OAO basis is returned
        if (norm2(arh%history_dm() - dm_oao) > tol) then
            write(stderr, *) "test_history_dm_arh_oao failed: Incorrect density matrix."
            test_history_dm_arh_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_history_dm_arh_oao

    logical(c_bool) function test_to_history_basis_arh_mo() bind(C)
        !
        ! this function tests the function which transforms a potential from the AO
        ! basis to the history basis for the MO basis
        !
        use otr_arh, only: arh_mo_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_occ, n_particle
        use otr_common_unit_tests, only: generate_random_symm_matrix
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object

        type(arh_mo_type) :: arh
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao), matrix(n_ao, n_ao, n_particle)
        integer(ip) :: k

        ! assume tests pass
        test_to_history_basis_arh_mo = .true.

        ! set up the MO object with random orthonormal MO coefficients and generate
        ! random potentials
        ao_overlap = generate_random_ao_overlap(n_ao)
        mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
        call setup_mo_object(mo_coeff, ao_overlap, n_occ)
        arh = arh_mo_type(mo_object)
        do k = 1, 2
            matrix(:, :, k) = generate_random_symm_matrix(n_ao)
        end do

        ! call routine and determine if the potentials are kept in the AO basis
        if (norm2(arh%to_history_basis(matrix) - matrix) > tol) then
            write(stderr, *) "test_to_history_basis_arh_mo failed: Incorrect potential."
            test_to_history_basis_arh_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object)

    end function test_to_history_basis_arh_mo

    logical(c_bool) function test_to_history_basis_arh_oao() bind(C)
        !
        ! this function tests the function which transforms a potential from the AO
        ! basis to the history basis for the OAO basis
        !
        use otr_arh, only: arh_oao_type
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_common_unit_tests, only: generate_random_symm_matrix

        type(arh_oao_type) :: arh
        real(rp) :: s_inv_sqrt(n_ao, n_ao), matrix(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle)
        integer(ip) :: k

        ! assume tests pass
        test_to_history_basis_arh_oao = .true.

        ! set up the OAO object with a random inverse square root of the overlap matrix
        ! and independently transform random potentials
        s_inv_sqrt = generate_random_symm_matrix(n_ao)
        do k = 1, n_particle
            matrix(:, :, k) = generate_random_symm_matrix(n_ao)
            expected(:, :, k) = matmul(s_inv_sqrt, matmul(matrix(:, :, k), s_inv_sqrt))
        end do
        allocate(oao_object)
        oao_object%s_inv_sqrt = s_inv_sqrt
        arh = arh_oao_type(oao_object)

        ! call routine and determine if the potentials are transformed with the inverse
        ! square root of the overlap matrix
        if (norm2(arh%to_history_basis(matrix) - expected) > tol) then
            write(stderr, *) "test_to_history_basis_arh_oao failed: Incorrect "// &
                "potential."
            test_to_history_basis_arh_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_to_history_basis_arh_oao

    logical(c_bool) function test_history_columns_arh_mo() bind(C)
        !
        ! this function tests the subroutine which expresses a history of difference
        ! matrices as history columns and as packed columns in the non-redundant
        ! parameter space for the MO basis
        !
        use otr_arh, only: arh_mo_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo, n_cases, case_n_particle, case_n_occ, &
                                         case_names
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_common_unit_tests, only: generate_random_symm_matrix
        use otr_mo_unit_tests, only: generate_random_ao_overlap, &
                                     generate_random_mo_coeff, setup_mo_object, &
                                     ref_mo_transform, ref_pack_ov

        integer(ip), parameter :: n_list = 2

        type(arh_mo_type) :: arh
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        real(rp), allocatable :: diff(:, :, :, :), diff_mo(:, :, :), cols(:, :), &
                                 packed(:, :), expected_cols(:, :), &
                                 expected_packed(:, :)
        integer(ip), allocatable :: n_occ(:)
        integer(ip) :: i_case, n_part, j, k
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_history_columns_arh_mo = .true.

        ! loop over every occupation case: the history columns are the full difference
        ! matrices in the current MO basis and the packed columns their
        ! occupied-virtual blocks
        do i_case = 1, n_cases
            case_name = trim(case_names(i_case))
            n_part = case_n_particle(i_case)
            n_occ = case_n_occ(:n_part, i_case)
            ao_overlap = generate_random_ao_overlap(n_ao)
            mo_coeff = generate_random_mo_coeff(ao_overlap, n_mo, n_particle)
            call setup_mo_object(mo_coeff(:, :, :n_part), ao_overlap, n_occ)
            arh = arh_mo_type(mo_object)
            allocate(diff(n_ao, n_ao, n_part, n_list), diff_mo(n_mo, n_mo, n_part), &
                     expected_cols(n_part * n_mo**2, n_list), &
                     expected_packed(mo_object%n_param, n_list))
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
                expected_packed(:, k) = ref_pack_ov(diff_mo, n_occ)
            end do

            ! call routine and determine if both sets of columns are correct
            call arh%history_columns(diff, cols, packed)
            if (any(shape(cols) /= shape(expected_cols))) then
                write(stderr, *) "test_history_columns_arh_mo failed: Incorrect "// &
                    "dimensions of history columns in the "//case_name//" case."
                test_history_columns_arh_mo = .false.
            else if (norm2(cols - expected_cols) > tol) then
                write(stderr, *) "test_history_columns_arh_mo failed: Incorrect "// &
                    "history columns in the "//case_name//" case."
                test_history_columns_arh_mo = .false.
            end if
            if (any(shape(packed) /= shape(expected_packed))) then
                write(stderr, *) "test_history_columns_arh_mo failed: Incorrect "// &
                    "dimensions of packed columns in the "//case_name//" case."
                test_history_columns_arh_mo = .false.
            else if (norm2(packed - expected_packed) > tol) then
                write(stderr, *) "test_history_columns_arh_mo failed: Incorrect "// &
                    "packed columns in the "//case_name//" case."
                test_history_columns_arh_mo = .false.
            end if
            deallocate(diff, diff_mo, expected_cols, expected_packed, mo_object)
        end do

    end function test_history_columns_arh_mo

    logical(c_bool) function test_history_columns_arh_oao() bind(C)
        !
        ! this function tests the subroutine which expresses a history of difference
        ! matrices as history columns and as packed columns in the non-redundant
        ! parameter space for the OAO basis
        !
        use otr_arh, only: arh_oao_type
        use otr_oao, only: oao_object
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle, n_occ
        use otr_common_unit_tests, only: generate_random_density_matrix, &
                                         generate_random_symm_matrix

        integer(ip), parameter :: n_list = 2

        type(arh_oao_type) :: arh
        real(rp) :: dm_oao(n_ao, n_ao, n_particle), diff(n_ao, n_ao, n_particle, n_list)
        real(rp), allocatable :: cols(:, :), packed(:, :)
        integer(ip) :: j, k

        ! assume tests pass
        test_history_columns_arh_oao = .true.

        ! set up the OAO object with a random density matrix, which the packed columns
        ! are projected with, and generate random history entries
        do j = 1, n_particle
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_occ(j))
            do k = 1, n_list
                diff(:, :, j, k) = generate_random_symm_matrix(n_ao)
            end do
        end do
        allocate(oao_object)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%dm_oao = dm_oao
        arh = arh_oao_type(oao_object)

        ! call routine and determine if the history columns are the flattened
        ! difference matrices and the packed columns their independently reproduced
        ! projections
        call arh%history_columns(diff, cols, packed)
        if (size(cols, 1) /= n_ao * n_ao * n_particle .or. size(cols, 2) /= n_list) then
            write(stderr, *) "test_history_columns_arh_oao failed: Incorrect "// &
                "dimensions of history columns."
            test_history_columns_arh_oao = .false.
        else if ( &
            norm2(cols - reshape(diff, [n_ao * n_ao * n_particle, n_list])) > tol) then
            write(stderr, *) "test_history_columns_arh_oao failed: History columns "// &
                "are not the flattened difference matrices."
            test_history_columns_arh_oao = .false.
        end if
        if (size(packed, 1) /= n_param .or. size(packed, 2) /= n_list) then
            write(stderr, *) "test_history_columns_arh_oao failed: Incorrect "// &
                "dimensions of packed columns."
            test_history_columns_arh_oao = .false.
        else if (norm2(packed - ref_cache_dirs(diff)) > tol) then
            write(stderr, *) "test_history_columns_arh_oao failed: Packed columns "// &
                "are not the projected and packed difference matrices."
            test_history_columns_arh_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    contains

        function ref_cache_dirs(v_diff) result(dirs)
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

        end function ref_cache_dirs

    end function test_history_columns_arh_oao

    logical(c_bool) function test_history_channel_rows_arh_mo() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in a history column for the MO basis
        !
        use otr_arh, only: arh_mo_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_ao, n_particle

        type(arh_mo_type) :: arh
        real(rp), target :: mo_coeff(n_ao, n_mo, n_particle)
        integer(ip), allocatable :: rows(:, :)

        ! assume tests pass
        test_history_channel_rows_arh_mo = .true.

        ! set up the MO object with the number of particle channels and the MO
        ! coefficients, whose shape is all the routine reads of them
        allocate(mo_object)
        mo_object%n_particle = 2
        mo_object%mo_coeff => mo_coeff
        arh = arh_mo_type(mo_object)

        ! call routine and determine if every channel holds the full difference matrix
        ! in the MO basis
        rows = arh%history_channel_rows()
        if (any(shape(rows) /= [2, 2])) then
            write(stderr, *) "test_history_channel_rows_arh_mo failed: Incorrect "// &
                "dimensions."
            test_history_channel_rows_arh_mo = .false.
        else if (any(rows /= reshape([1_ip, n_mo**2, &
                                      n_mo**2 + 1_ip, 2_ip * n_mo**2], [2, 2]))) then
            write(stderr, *) "test_history_channel_rows_arh_mo failed: Incorrect rows."
            test_history_channel_rows_arh_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object)

    end function test_history_channel_rows_arh_mo

    logical(c_bool) function test_history_channel_rows_arh_oao() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in a history column for the OAO basis
        !
        use otr_arh, only: arh_oao_type
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao, n_particle

        type(arh_oao_type) :: arh
        integer(ip), allocatable :: rows(:, :)
        integer(ip) :: j

        ! assume tests pass
        test_history_channel_rows_arh_oao = .true.

        ! set up the OAO object with the dimensions the routine reads
        allocate(oao_object)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        arh = arh_oao_type(oao_object)

        ! call routine and determine if every channel holds a full difference matrix
        rows = arh%history_channel_rows()
        if (any(shape(rows) /= [2_ip, n_particle])) then
            write(stderr, *) "test_history_channel_rows_arh_oao failed: Incorrect "// &
                "dimensions."
            test_history_channel_rows_arh_oao = .false.
        else if (any(rows /= &
                     reshape([((j - 1) * n_ao**2 + 1, j * n_ao**2, j=1, n_particle)], &
                             [2_ip, n_particle]))) then
            write(stderr, *) "test_history_channel_rows_arh_oao failed: Incorrect rows."
            test_history_channel_rows_arh_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_history_channel_rows_arh_oao

    logical(c_bool) function test_packed_channel_rows_arh_mo() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in a packed column for the MO basis
        !
        use otr_arh, only: arh_mo_type
        use otr_mo, only: mo_object
        use otr_mo_test_reference, only: n_mo
        use otr_common_test_reference, only: n_occ
        use otr_mo_unit_tests, only: setup_minimal_mo_object

        type(arh_mo_type) :: arh
        integer(ip), allocatable :: rows(:, :)
        integer(ip) :: n1, n2

        ! assume tests pass
        test_packed_channel_rows_arh_mo = .true.

        ! set up the MO object with the differently occupied spin channels, whose
        ! occupations are all the routine reads
        call setup_minimal_mo_object(n_occ)
        arh = arh_mo_type(mo_object)

        ! call routine and determine if every channel holds its occupied-virtual
        ! rotations, whose number differs between the spin channels
        n1 = n_occ(1) * (n_mo - n_occ(1))
        n2 = n_occ(2) * (n_mo - n_occ(2))
        rows = arh%packed_channel_rows()
        if (any(shape(rows) /= [2, 2])) then
            write(stderr, *) "test_packed_channel_rows_arh_mo failed: Incorrect "// &
                "dimensions."
            test_packed_channel_rows_arh_mo = .false.
        else if (any(rows /= reshape([1_ip, n1, &
                                      n1 + 1_ip, n1 + n2], [2, 2]))) then
            write(stderr, *) "test_packed_channel_rows_arh_mo failed: Incorrect rows."
            test_packed_channel_rows_arh_mo = .false.
        end if

        ! deallocate MO object
        deallocate(mo_object)

    end function test_packed_channel_rows_arh_mo

    logical(c_bool) function test_packed_channel_rows_arh_oao() bind(C)
        !
        ! this function tests the function which returns the rows every particle
        ! channel occupies in a packed column for the OAO basis
        !
        use otr_arh, only: arh_oao_type
        use otr_oao, only: oao_object
        use otr_common_test_reference, only: n_ao, n_particle

        integer(ip), parameter :: n_channel = n_ao * (n_ao - 1) / 2

        type(arh_oao_type) :: arh
        integer(ip), allocatable :: rows(:, :)
        integer(ip) :: j

        ! assume tests pass
        test_packed_channel_rows_arh_oao = .true.

        ! set up the OAO object with the dimensions the routine reads
        allocate(oao_object)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        arh = arh_oao_type(oao_object)

        ! call routine and determine if every channel holds the antisymmetric
        ! parameters of that channel
        rows = arh%packed_channel_rows()
        if (any(shape(rows) /= [2_ip, n_particle])) then
            write(stderr, *) "test_packed_channel_rows_arh_oao failed: Incorrect "// &
                "dimensions."
            test_packed_channel_rows_arh_oao = .false.
        else if (any(rows /= reshape( &
            [((j - 1) * n_channel + 1, j * n_channel, j=1, n_particle)], &
            [2_ip, n_particle]))) then
            write(stderr, *) "test_packed_channel_rows_arh_oao failed: Incorrect rows."
            test_packed_channel_rows_arh_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_packed_channel_rows_arh_oao

    logical(c_bool) function test_hess_x_static_arh_mo() bind(C)
        !
        ! this function tests the function which applies the static part of the Hessian
        ! to a trial vector for the MO basis, which applies the one of the MO object
        !
        use otr_arh, only: arh_mo_type
        use otr_mo, only: mo_object
        use otr_common_test_reference, only: n_occ
        use otr_mo_unit_tests, only: setup_minimal_mo_object
        use otr_common_unit_tests, only: identity_matrix

        type(arh_mo_type) :: arh
        real(rp), allocatable :: x(:)
        integer(ip) :: k

        ! assume tests pass
        test_hess_x_static_arh_mo = .true.

        ! set up the MO object with vanishing occupied-occupied and unit
        ! virtual-virtual blocks of the Fock matrix, so that the static part of the
        ! Hessian, the open-shell 2 (X F_vv - F_oo X), doubles the trial vector
        call setup_minimal_mo_object(n_occ)
        do k = 1, size(n_occ, kind=ip)
            associate (channel => mo_object%mo_channels(k))
                channel%fock_oo = 0.0_rp
                channel%fock_vv = identity_matrix(channel%n_virt)
            end associate
        end do
        arh = arh_mo_type(mo_object)
        allocate(x(mo_object%n_param))
        call random_number(x)

        ! call routine and determine if the trial vector is doubled
        if (norm2(arh%hess_x_static(x) - 2.0_rp * x) > tol) then
            write(stderr, *) "test_hess_x_static_arh_mo failed: Incorrect static part."
            test_hess_x_static_arh_mo = .false.
        end if
        deallocate(mo_object)

    end function test_hess_x_static_arh_mo

    logical(c_bool) function test_hess_x_static_arh_oao() bind(C)
        !
        ! this function tests the function which applies the static part of the Hessian
        ! to a trial vector for the OAO basis
        !
        use otr_arh, only: arh_oao_type
        use otr_oao, only: oao_object
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle, n_occ
        use otr_oao_unit_tests, only: ref_unpack_asymm, ref_hess_x_oao

        type(arh_oao_type) :: arh
        real(rp) :: dm_oao(n_ao, n_ao, n_particle), fock_oo(n_ao, n_ao, n_particle), &
                    fock_vv(n_ao, n_ao, n_particle), x(n_param), &
                    response(n_ao, n_ao, n_particle)

        ! assume tests pass
        test_hess_x_static_arh_oao = .true.

        ! set up the OAO object with a random partition of the Fock matrix
        call generate_random_fock_partition(n_ao, n_occ(:n_particle), dm_oao, fock_oo, &
                                            fock_vv)
        allocate(oao_object)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%dm_oao = dm_oao
        oao_object%fock_oo = fock_oo
        oao_object%fock_vv = fock_vv
        arh = arh_oao_type(oao_object)

        ! call routine and determine if the static part matches the independently
        ! assembled one
        x = generate_random_nonredundant_vector(n_param, n_particle, n_ao, dm_oao)
        response = 0.0_rp
        if (norm2(arh%hess_x_static(x) - &
                  ref_hess_x_oao(ref_unpack_asymm(x, n_particle, n_ao), response, &
                                 dm_oao, fock_oo, fock_vv, n_param)) > tol) then
            write(stderr, *) "test_hess_x_static_arh_oao failed: Incorrect static part."
            test_hess_x_static_arh_oao = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_hess_x_static_arh_oao

    logical(c_bool) function test_cache_channel_split_dirs() bind(C)
        !
        ! this function tests the open-shell history-projection caching routine which
        ! keeps the two spin channels of every packed column as separate columns,
        ! zeroing the rows of the respective other channel, then rebases the result via
        ! a given (map, chol) pair; undoing the rebasing by right-multiplying by chol
        ! again must reproduce the gathered, reordered raw split columns exactly
        !
        use otr_arh, only: cache_channel_split_dirs
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
        chol = generate_random_upper_triangular(n_col)

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
        use otr_arh, only: cache_combined_channel_dirs
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
        chol = generate_random_upper_triangular(n_col)

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
        use otr_arh, only: get_low_rank_hess_factors, arh_object, arh_oao_type, &
                           arh_types
        use otr_oao, only: oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_unit_tests, only: identity_matrix, generate_random_symm_matrix, &
                                         shell_names
        use otr_common_test_reference, only: n_particle_ref => n_particle

        integer(ip), parameter :: n_param = 5, n_diff = 3, n_diff_nl = 2, &
                                  n_both = n_diff + n_diff_nl, dm_nl = 2 * n_diff, &
                                  pot_nl = 2 * n_diff + n_diff_nl

        real(rp) :: density(n_param, n_diff), density_nl(n_param, n_diff_nl), &
                    linear(n_param, n_diff), nonlinear(n_param, n_diff_nl), &
                    a_sym(n_diff, n_diff), a_sym_nl(n_diff_nl, n_diff_nl), &
                    a_inv(n_diff, n_diff), a_inv_comb(n_diff_nl, n_diff_nl), scale
        real(rp), allocatable :: expansion(:, :), projection(:, :), coupling(:, :)
        integer(ip) :: i_type, n_particle
        character(len=:), allocatable :: arh_type, case_name

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

        ! set up the ARH object with the cached quantities the assembly requires
        allocate(oao_object)
        oao_object%n_param = n_param
        allocate(arh_object, source=arh_oao_type(oao_object))
        call setup_settings(arh_object%settings)
        arh_object%dm_dirs = density
        arh_object%dm_dirs_nonlinear = density_nl
        arh_object%linear_potential_dirs = linear
        arh_object%nonlinear_potential_dirs = nonlinear
        arh_object%a_sym = a_sym
        arh_object%a_sym_nonlinear = a_sym_nl
        arh_object%a_inv = a_inv
        arh_object%a_inv_comb = a_inv_comb

        ! call routine for every ARH type and both shells, whose coupling matrices
        ! differ by a factor of two, and determine if the factors match
        do i_type = 1, size(arh_types)
            arh_type = trim(arh_types(i_type))
            do n_particle = 1, n_particle_ref
                case_name = arh_type//" in the "//trim(shell_names(n_particle))//" case"
                scale = merge(8.0_rp, 4.0_rp, n_particle == 1)

                ! the expected factors: every system (linear and non-linear) has its
                ! own block of directions, which no coupling mixes across systems
                select case (arh_type)
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
                    coupling = merge(0.5_rp, 1.0_rp, arh_type == "symm_arh") * scale * &
                               block_diagonal_matrix( &
                                   cshift(identity_matrix(2 * n_diff), n_diff), &
                                   cshift(identity_matrix(2 * n_diff_nl), n_diff_nl))
                    if (arh_type == "ms_psb") then
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

                ! call routine and determine if the factors match
                oao_object%n_particle = n_particle
                arh_object%settings%arh_type = arh_type
                call get_low_rank_hess_factors()
                if (.not. (allocated(arh_object%expansion_dirs) .and. &
                           allocated(arh_object%projection_dirs) .and. &
                           allocated(arh_object%coupling_matrix))) then
                    write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                        "Factors not assembled for "//case_name//"."
                    test_get_low_rank_hess_factors = .false.
                else if ( &
                    any(shape(arh_object%expansion_dirs) /= shape(expansion)) .or. &
                    any(shape(arh_object%projection_dirs) /= shape(projection)) .or. &
                    any(shape(arh_object%coupling_matrix) /= shape(coupling))) then
                    write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                        "Incorrect dimensions of the factors for "//case_name//"."
                    test_get_low_rank_hess_factors = .false.
                else
                    if (norm2(arh_object%expansion_dirs - expansion) > tol) then
                        write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                            "Incorrect expansion directions for "//case_name//"."
                        test_get_low_rank_hess_factors = .false.
                    end if
                    if (norm2(arh_object%projection_dirs - projection) > tol) then
                        write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                            "Incorrect projection directions for "//case_name//"."
                        test_get_low_rank_hess_factors = .false.
                    end if
                    if (norm2(arh_object%coupling_matrix - coupling) > tol) then
                        write(stderr, *) "test_get_low_rank_hess_factors failed: "// &
                            "Incorrect coupling matrix for "//case_name//"."
                        test_get_low_rank_hess_factors = .false.
                    end if
                end if
            end do
        end do

        ! an empty history has to leave every factor unallocated
        deallocate(arh_object%dm_dirs, arh_object%dm_dirs_nonlinear, &
                   arh_object%linear_potential_dirs, &
                   arh_object%nonlinear_potential_dirs)
        allocate(arh_object%dm_dirs(n_param, 0), &
                 arh_object%dm_dirs_nonlinear(n_param, 0), &
                 arh_object%linear_potential_dirs(n_param, 0), &
                 arh_object%nonlinear_potential_dirs(n_param, 0))
        call get_low_rank_hess_factors()
        if (allocated(arh_object%expansion_dirs) .or. &
            allocated(arh_object%projection_dirs) .or. &
            allocated(arh_object%coupling_matrix)) then
            write(stderr, *) "test_get_low_rank_hess_factors failed: Factors were "// &
                "assembled for an empty history."
            test_get_low_rank_hess_factors = .false.
        end if

        ! deallocate ARH and OAO objects
        deallocate(arh_object, oao_object)

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
        use otr_arh, only: build_a_part

        integer(ip), parameter :: n_diff = 3, n_rows = 7
        real(rp) :: dm_cols(n_rows, n_diff), v_cols(n_rows, n_diff), &
                    a(n_diff, n_diff), expected(n_diff, n_diff)

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

    end function test_build_a_part

    logical(c_bool) function test_build_a_transformed() bind(C)
        !
        ! this function tests the function which builds the A = S^T Y matrix of a
        ! single part of the response and congruence-transforms it into the
        ! orthonormalized basis of that part: undoing the transform by left- and
        ! right-multiplying by chol again must reproduce the gathered, reordered raw
        ! symmetrized matrix exactly
        !
        use otr_arh, only: build_a_transformed

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
        chol = generate_random_upper_triangular(n_diff)

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
        use otr_arh, only: build_a_block_linear_os
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
        chol = generate_random_upper_triangular(n_col)

        ! generate A matrix and verify
        a_block = &
            build_a_block_linear_os(dm_cols, v_same_cols, v_opp_cols, rows, map, chol)
        if (size(a_block, 1) /= n_col .or. size(a_block, 2) /= n_col) then
            write(stderr, *) "test_build_a_block_linear_os failed: Incorrect "// &
                "dimensions."
            test_build_a_block_linear_os = .false.
            return
        end if
        if (norm2(matmul(transpose(chol), matmul(a_block, chol)) - &
                  expected_a_block(map, map)) > tol) then
            write(stderr, *) "test_build_a_block_linear_os failed: Result does not "// &
                "invert back to the gathered, reordered raw block matrix."
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
        use otr_arh, only: build_a_block_nonlinear_os
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
        chol = generate_random_upper_triangular(n_col)

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
        use otr_arh, only: cross_symmetrize

        real(rp) :: a12(2, 2), a21(2, 2), expected12(2, 2), expected21(2, 2)

        ! assume tests pass
        test_cross_symmetrize = .true.

        ! initialize blocks which are not transposes of each other
        a12 = reshape([1.0_rp, 3.0_rp, &
                       2.0_rp, 4.0_rp], [2, 2])
        a21 = reshape([5.0_rp, 6.0_rp, &
                       7.0_rp, 9.0_rp], [2, 2])

        ! initialize expected blocks, averaging a12(i, k) with a21(k, i) directly
        expected12 = reshape([3.0_rp, 5.0_rp, &
                              4.0_rp, 6.5_rp], [2, 2])
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
            write(stderr, *) "test_truncated_eigval_inv failed: Incorrect "// &
                "truncated inverse eigenvalues."
            test_truncated_eigval_inv = .false.
        end if

    end function test_truncated_eigval_inv

    logical(c_bool) function test_history_step_mask() bind(C)
        !
        ! this function tests the function which decides, from step length alone, which
        ! history entries the non-linear multisecant system is allowed to use
        !
        use otr_arh, only: history_step_mask

        integer(ip), parameter :: n_rows = 4, n_diff = 3, n_rows_os = 8, n_diff_os = 2

        real(rp) :: dm_cols(n_rows, n_diff), dm_cols_os(n_rows_os, n_diff_os), &
                    dm_cols_single(n_rows, 1)
        logical, allocatable :: keep(:)

        ! assume tests pass
        test_history_step_mask = .true.

        ! initialize a history whose third step is far longer than the shortest one,
        ! and whose second is comfortably within reach of it
        dm_cols = 0.0_rp
        dm_cols(1, 1) = 1.0_rp
        dm_cols(3, 2) = 5.0_rp
        dm_cols(2, 3) = 1e4_rp

        ! call routine and determine that only the entry beyond the cutoff is dropped
        keep = history_step_mask(dm_cols)
        if (size(keep) /= n_diff) then
            write(stderr, *) "test_history_step_mask failed: Incorrect number of "// &
                "entries."
            test_history_step_mask = .false.
            return
        end if
        if (.not. (keep(1) .and. keep(2))) then
            write(stderr, *) "test_history_step_mask failed: An entry within the "// &
                "cutoff was not kept."
            test_history_step_mask = .false.
        end if
        if (keep(3)) then
            write(stderr, *) "test_history_step_mask failed: An entry beyond the "// &
                "cutoff was not dropped."
            test_history_step_mask = .false.
        end if
        deallocate(keep)

        ! scaling the whole history leaves the decision unchanged, since the cutoff is
        ! a ratio to the shortest step rather than an absolute length
        dm_cols = 1e-6_rp * dm_cols
        keep = history_step_mask(dm_cols)
        if (.not. (keep(1) .and. keep(2)) .or. keep(3)) then
            write(stderr, *) "test_history_step_mask failed: The decision is not "// &
                "invariant under an overall scaling of the history."
            test_history_step_mask = .false.
        end if
        deallocate(keep)

        ! bringing the long step within the cutoff keeps the entire history
        dm_cols(2, 3) = dm_cols(1, 1)
        keep = history_step_mask(dm_cols)
        if (.not. all(keep)) then
            write(stderr, *) "test_history_step_mask failed: A history whose steps "// &
                "are all comparable was not kept in full."
            test_history_step_mask = .false.
        end if
        deallocate(keep)

        ! the open-shell step length spans both spin channels, so an entry that is long
        ! only in the rows of the second one still has to be dropped
        dm_cols_os = 0.0_rp
        dm_cols_os(1, 1) = 1.0_rp
        dm_cols_os(1, 2) = 1.0_rp
        dm_cols_os(6, 2) = 1e4_rp
        keep = history_step_mask(dm_cols_os)
        if (.not. keep(1) .or. keep(2)) then
            write(stderr, *) "test_history_step_mask failed: An entry beyond the "// &
                "cutoff only through its second spin channel was not dropped."
            test_history_step_mask = .false.
        end if
        deallocate(keep)

        ! a single entry has nothing to be judged against and is always kept
        dm_cols_single = 0.0_rp
        dm_cols_single(1, 1) = 1.0_rp
        keep = history_step_mask(dm_cols_single)
        if (size(keep) /= 1 .or. .not. keep(1)) then
            write(stderr, *) "test_history_step_mask failed: A single history "// &
                "entry was not kept."
            test_history_step_mask = .false.
        end if
        deallocate(keep)

    end function test_history_step_mask

    logical(c_bool) function test_median() bind(C)
        !
        ! this function tests the function which returns the median of an array, for an
        ! odd and an even number of entries and for an empty array
        !
        use otr_arh, only: median

        integer(ip), parameter :: n_odd = 5, n_even = 6
        real(rp) :: odd(n_odd), even(n_even)

        ! assume tests pass
        test_median = .true.

        ! an unordered odd number of entries has a single central one, which is not the
        ! mean of the array
        odd = [3.0_rp, -1.0_rp, 10.0_rp, 0.5_rp, 2.0_rp]
        if (abs(median(odd) - 2.0_rp) > tol) then
            write(stderr, *) "test_median failed: Incorrect median for an odd "// &
                "number of entries."
            test_median = .false.
        end if

        ! an even number of entries averages the two central ones, which a single large
        ! outlier must not move
        even = [3.0_rp, -1.0_rp, 10.0_rp, 0.5_rp, 2.0_rp, 100.0_rp]
        if (abs(median(even) - 2.5_rp) > tol) then
            write(stderr, *) "test_median failed: Incorrect median for an even "// &
                "number of entries."
            test_median = .false.
        end if

        ! an empty array has no entries to take a median of
        if (abs(median([real(rp) :: ])) > tol) then
            write(stderr, *) "test_median failed: Median of an empty array does "// &
                "not vanish."
            test_median = .false.
        end if

    end function test_median

    logical(c_bool) function test_resolvable_residual() bind(C)
        !
        ! this function tests the function which returns the shortest residual norm a
        ! history direction has to contribute for its response to describe curvature
        ! rather than the error in it
        !
        use otr_arh, only: resolvable_residual

        integer(ip), parameter :: flat_len = 4, n_diff = 4
        real(rp), parameter :: lengths(n_diff) = [1.0_rp, 2.0_rp, 0.5_rp, 4.0_rp]

        real(rp) :: basis(flat_len, n_diff), steps(flat_len, n_diff), &
                    responses(flat_len, n_diff), q(flat_len, n_diff), &
                    z(flat_len, n_diff), a(n_diff, n_diff), asymmetries(n_diff), &
                    curvatures(n_diff), expected
        logical :: keep(n_diff)
        integer(ip) :: i, k

        ! assume tests pass
        test_resolvable_residual = .true.

        ! orthogonal steps of different lengths, so that the orthonormalization the
        ! routine performs internally is a scaling by those lengths regardless of the
        ! order it accepts them in, while every quantity below is a median over
        ! directions and therefore independent of that order
        basis = 0.5_rp * reshape([1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp, 1.0_rp, -1.0_rp, &
                                  1.0_rp, -1.0_rp, 1.0_rp, 1.0_rp, -1.0_rp, -1.0_rp, &
                                  1.0_rp, -1.0_rp, -1.0_rp, 1.0_rp], [flat_len, n_diff])
        do i = 1, n_diff
            steps(:, i) = lengths(i) * basis(:, i)
        end do

        ! random responses, so that the transformed response matrix is generically
        ! asymmetric and its asymmetry is a genuine error estimate
        call random_number(responses)
        responses = responses - 0.5_rp

        ! the transformed response matrix follows from the same scaling
        q = basis
        do i = 1, n_diff
            z(:, i) = responses(:, i) / lengths(i)
        end do
        a = matmul(transpose(q), z)

        ! the error of a direction is the typical asymmetry of its own column, and
        ! undoing the amplification of every column recovers the error of the response
        ! itself, which the typical curvature turns into a length
        do k = 1, n_diff
            asymmetries(k) = ref_median(pack(abs(a(:, k) - a(k, :)), &
                                             [(i /= k, i=1, n_diff)]))
            curvatures(k) = abs(a(k, k))
        end do
        expected = ref_median(asymmetries * lengths) / ref_median(curvatures)

        ! call routine and determine if the threshold matches
        keep = .true.
        if (abs(resolvable_residual(steps, responses, keep) - expected) > tol) then
            write(stderr, *) "test_resolvable_residual failed: Incorrect shortest "// &
                "resolvable residual norm."
            test_resolvable_residual = .false.
        end if

        ! a direction needs at least two partners for the asymmetry of its column to be
        ! a typical value, so a shorter history yields no threshold at all
        keep = [.true., .true., .false., .false.]
        if (abs(resolvable_residual(steps, responses, keep)) > tol) then
            write(stderr, *) "test_resolvable_residual failed: Threshold does not "// &
                "vanish for a history too short to estimate an error scale from."
            test_resolvable_residual = .false.
        end if

    contains

        function ref_order_statistic(x, rank) result(val)
            !
            ! this function independently reproduces the entry of a given rank, counted
            ! from the smallest, of an array by counting, for every entry, how many
            ! entries lie below and alongside it
            !
            real(rp), intent(in) :: x(:)
            integer(ip), intent(in) :: rank
            real(rp) :: val

            integer(ip) :: j, n_below, n_at_most

            val = 0.0_rp
            do j = 1, size(x)
                n_below = count(x < x(j))
                n_at_most = count(x <= x(j))
                if (n_below < rank .and. rank <= n_at_most) then
                    val = x(j)
                    return
                end if
            end do

        end function ref_order_statistic

        function ref_median(x) result(med)
            !
            ! this function independently reproduces the median of an array from its
            ! order statistics
            !
            real(rp), intent(in) :: x(:)
            real(rp) :: med

            integer(ip) :: n

            n = size(x)
            if (n == 0) then
                med = 0.0_rp
            else if (mod(n, 2_ip) == 1) then
                med = ref_order_statistic(x, (n + 1) / 2)
            else
                med = 0.5_rp * (ref_order_statistic(x, n / 2) + &
                                ref_order_statistic(x, n / 2 + 1))
            end if

        end function ref_median

    end function test_resolvable_residual

    logical(c_bool) function test_factorize_history() bind(C)
        !
        ! this function tests the subroutine which performs a pivoted, rank-revealing
        ! Cholesky factorization of the Gram matrix of a set of flattened history
        ! vectors
        !
        integer(ip), parameter :: n_rows = 4, n_diff = 3
        real(rp), parameter :: small_norm = 1e-8_rp, large_norm = 1e6_rp, &
                               resolvable_norm = 1e-3_rp

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

        ! a column repeating the direction of the first one but contributing a new
        ! component of resolvable norm is rejected when more than that norm is demanded
        ! and accepted when less is demanded, while the pivoting then ranks the columns
        ! by length; without a demand, it is accepted as independent, so that it is the
        ! demand rather than linear dependence that rejects it
        cols(:, 3) = [0.5_rp, 0.0_rp, resolvable_norm, 0.0_rp]
        if (.not. check_factorization( &
            "a column contributing less than the demanded residual", cols, &
            [2_ip, 1_ip], min_residual=1e2_rp * resolvable_norm)) &
            test_factorize_history = .false.
        if (.not. check_factorization( &
            "a column contributing more than the demanded residual", cols, &
            [2_ip, 1_ip, 3_ip], min_residual=1e-2_rp * resolvable_norm)) &
            test_factorize_history = .false.
        if (.not. check_factorization("an omitted demanded residual", cols, &
                                      [1_ip, 2_ip, 3_ip])) &
            test_factorize_history = .false.

        ! an empty history returns an empty factorization
        if (.not. check_factorization("an empty history", cols(:, :0), &
                                      [integer(ip) :: ])) &
            test_factorize_history = .false.

    contains

        function check_factorization(case_name, history_cols, expected_map, keep, &
                                     min_residual) result(passed)
            !
            ! this function checks the pivoted Cholesky factorization of the Gram matrix
            ! of a set of history columns for one case: the accepted columns have to
            ! follow the expected map back to the original history indices and the
            ! factor has to reproduce the Gram matrix of the accepted columns, which
            ! pins it since it is upper triangular with a positive diagonal
            !
            use otr_arh, only: factorize_history

            character(len=*), intent(in) :: case_name
            real(rp), intent(in) :: history_cols(:, :)
            integer(ip), intent(in) :: expected_map(:)
            logical, intent(in), optional :: keep(:)
            real(rp), intent(in), optional :: min_residual
            logical :: passed

            real(rp), allocatable :: chol(:, :), gram(:, :)
            integer(ip), allocatable :: map(:)
            integer(ip) :: n_accepted

            ! assume test passes
            passed = .true.

            ! call routine and determine if the accepted columns and the factor match
            call factorize_history(history_cols, chol, map, n_accepted, keep, &
                                   min_residual)
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
        use otr_arh, only: congruence_transform
        use otr_common_unit_tests, only: generate_random_symm_matrix

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
        use otr_arh, only: combine_channels

        integer(ip), parameter :: n1 = 2, n2 = 1, n_offset = 5

        real(rp) :: chol1(n1, n1), chol2(n2, n2)
        integer(ip) :: map1(n1), map2(n2)
        real(rp), allocatable :: chol_comb(:, :)
        integer(ip), allocatable :: map_comb(:)

        ! assume tests pass
        test_combine_channels = .true.

        chol1 = reshape([1.0_rp, 0.0_rp, &
                         2.0_rp, 3.0_rp], [n1, n1])
        chol2 = reshape([4.0_rp], [n2, n2])
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
        use otr_arh, only: get_ms_a_inv, arh_settings_type
        use otr_common_test_reference, only: n_ao
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_diff = 3, n_accepted = 2, flat_len = n_ao * n_ao

        real(rp) :: chol(n_accepted, n_accepted), dm_cols(flat_len, n_diff), &
                    v_cols(flat_len, n_diff)
        real(rp), allocatable :: a_inv(:, :), expected(:, :), a_tilde(:, :), &
                                 y_gram(:, :)
        integer(ip) :: map(n_accepted), error
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv = .true.

        ! setup settings object
        call setup_settings(settings)

        ! reorder and rebase into an arbitrary orthonormalized basis that also drops
        ! one history entry
        map = [3, 1]
        chol = generate_random_upper_triangular(n_accepted)

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

    logical(c_bool) function test_get_ms_a_inv_os_linear() bind(C)
        !
        ! this function tests the subroutine which computes the pseudoinverse
        ! multisecant SR1 matrix in a spin-separated manner for the linear part in the
        ! open-shell case
        !
        use otr_arh, only: get_ms_a_inv_os_linear, arh_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_diff = 2, n_col = n_particle * n_diff, &
                                  n_accepted = 3, n_rows(n_particle) = [4, 3]

        real(rp) :: dm_cols(sum(n_rows), n_diff), &
                    v_same_spin_cols(sum(n_rows), n_diff), &
                    v_opposite_spin_cols(sum(n_rows), n_diff), &
                    empty_cols(sum(n_rows), 0), chol(n_accepted, n_accepted), &
                    s_full(sum(n_rows), n_col), y_full(sum(n_rows), n_col)
        real(rp), allocatable :: a_inv(:, :), expected(:, :), a_tilde(:, :), &
                                 y_gram(:, :)
        integer(ip) :: rows(2, n_particle), map(n_accepted), error, i, j, k
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_ms_a_inv_os_linear = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize random history and spin-resolved potentials centered on zero,
        ! reordered and rebased into an arbitrary orthonormalized basis that also drops
        ! one column
        call random_number(dm_cols)
        call random_number(v_same_spin_cols)
        call random_number(v_opposite_spin_cols)
        dm_cols = dm_cols - 0.5_rp
        v_same_spin_cols = v_same_spin_cols - 0.5_rp
        v_opposite_spin_cols = v_opposite_spin_cols - 0.5_rp
        rows = reshape([1_ip, n_rows(1), &
                        n_rows(1) + 1_ip, sum(n_rows)], [2, 2])
        map = [4, 1, 3]
        chol = generate_random_upper_triangular(n_accepted)

        ! the expected pseudoinverse is assembled independently of the routine from the
        ! stacked history and response, whose products give both the blocks of A and
        ! the response Gram matrix; A is symmetrized only after the transform, since
        ! its two halves are exactly transposes of one another in exact arithmetic;
        ! each history column touches only its own spin block, while each response
        ! column pairs the same-spin potential of one channel with the opposite-spin
        ! potential of the other, so that S^T Y reproduces the four blocks of A and Y^T
        ! Y the response Gram matrix
        s_full = 0.0_rp
        do k = 1, n_diff
            do i = rows(1, 1), rows(2, 1)
                s_full(i, k) = dm_cols(i, k)
                y_full(i, k) = v_same_spin_cols(i, k)
                y_full(i, n_diff + k) = v_opposite_spin_cols(i, k)
            end do
            do i = rows(1, 2), rows(2, n_particle)
                s_full(i, n_diff + k) = dm_cols(i, k)
                y_full(i, k) = v_opposite_spin_cols(i, k)
                y_full(i, n_diff + k) = v_same_spin_cols(i, k)
            end do
        end do
        a_tilde = ref_congruence_transform(matmul(transpose(s_full), y_full), map, chol)
        a_tilde = 0.5_rp * (a_tilde + transpose(a_tilde))
        y_gram = ref_congruence_transform(matmul(transpose(y_full), y_full), map, chol)
        expected = ref_ms_a_inv(a_tilde, y_gram)

        ! call routine and determine if the pseudoinverse matches
        call get_ms_a_inv_os_linear(dm_cols, v_same_spin_cols, v_opposite_spin_cols, &
                                    rows, map, chol, a_inv, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_a_inv_os_linear failed: Produced error."
            test_get_ms_a_inv_os_linear = .false.
        else if (norm2(a_inv - expected) > tol) then
            write(stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "pseudoinverse."
            test_get_ms_a_inv_os_linear = .false.
        end if
        deallocate(a_inv, expected, a_tilde, y_gram)

        ! confining the history to the first row of every channel while suppressing the
        ! potentials there leaves every direction badly aligned with its own response,
        ! so the screening criterion has to discard all of them
        do j = 1, 2
            dm_cols(rows(1, j) + 1:rows(2, j), :) = 0.0_rp
            v_same_spin_cols(rows(1, j), :) = 1e-3_rp * v_same_spin_cols(rows(1, j), :)
            v_opposite_spin_cols(rows(1, j), :) = 1e-3_rp * &
                                                  v_opposite_spin_cols(rows(1, j), :)
        end do

        ! call routine and determine if the pseudoinverse vanishes
        call get_ms_a_inv_os_linear(dm_cols, v_same_spin_cols, v_opposite_spin_cols, &
                                    rows, map, chol, a_inv, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_get_ms_a_inv_os_linear failed: Produced error "// &
                "for a badly aligned history."
            test_get_ms_a_inv_os_linear = .false.
        else if (norm2(a_inv) > tol) then
            write(stderr, *) "test_get_ms_a_inv_os_linear failed: Screening did "// &
                "not discard every direction of a badly aligned history."
            test_get_ms_a_inv_os_linear = .false.
        end if
        deallocate(a_inv)

        ! call routine for an empty history and determine if dimensions of the
        ! resulting pseudoinverse vanish
        call get_ms_a_inv_os_linear( &
            empty_cols, empty_cols, empty_cols, rows, [integer(ip) :: ], &
            reshape([real(rp) :: ], [0, 0]), a_inv, settings, error)
        if (size(a_inv, 1) /= 0 .or. size(a_inv, 2) /= 0) then
            write(stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "pseudoinverse dimensions for empty history."
            test_get_ms_a_inv_os_linear = .false.
        end if
        deallocate(a_inv)

    end function test_get_ms_a_inv_os_linear

    logical(c_bool) function test_response_gram() bind(C)
        !
        ! this function tests the function which returns the Gram matrix of the
        ! response history rebased into the orthonormalized S-basis
        !
        use otr_arh, only: response_gram

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
        chol = generate_random_upper_triangular(n_accepted)

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

    logical(c_bool) function test_response_gram_os_linear() bind(C)
        !
        ! this function tests the function which returns the Gram matrix of the
        ! open-shell linear response history, whose same-spin and opposite-spin
        ! potentials are interleaved exactly as the rows of A pair them
        !
        use otr_arh, only: response_gram_os_linear
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_dm = 2, n_accepted = 3, n_rows(n_particle) = [4, 3]
        integer(ip) :: rows(2, n_particle), map(n_accepted), k
        real(rp) :: v_same(sum(n_rows), n_dm), v_opp(sum(n_rows), n_dm), &
                    chol(n_accepted, n_accepted), y_full(sum(n_rows), 2 * n_dm), &
                    gram(2 * n_dm, 2 * n_dm)
        real(rp), allocatable :: y_gram(:, :), expected(:, :)

        ! assume tests pass
        test_response_gram_os_linear = .true.

        ! random same-spin and opposite-spin response histories
        call random_number(v_same)
        call random_number(v_opp)
        rows = reshape([1_ip, n_rows(1), &
                        n_rows(1) + 1_ip, sum(n_rows)], [2, 2])
        map = [4_ip, 1_ip, 3_ip]
        chol = generate_random_upper_triangular(n_accepted)

        ! stack the alpha and beta blocks independently, in the interleaved column
        ! order the routine documents
        do k = 1, n_dm
            y_full(:n_rows(1), k) = v_same(:n_rows(1), k)
            y_full(n_rows(1) + 1:, k) = v_opp(n_rows(1) + 1:, k)
            y_full(:n_rows(1), n_dm + k) = v_opp(:n_rows(1), k)
            y_full(n_rows(1) + 1:, n_dm + k) = v_same(n_rows(1) + 1:, k)
        end do
        gram = matmul(transpose(y_full), y_full)
        expected = ref_congruence_transform(gram, map, chol)

        ! call routine and determine if the rebased Gram matrix matches
        y_gram = response_gram_os_linear(v_same, v_opp, rows, map, chol)
        if (size(y_gram, 1) /= n_accepted .or. size(y_gram, 2) /= n_accepted) then
            write(stderr, *) "test_response_gram_os_linear failed: Incorrect shape."
            test_response_gram_os_linear = .false.
        else if (maxval(abs(y_gram - expected)) > tol) then
            write(stderr, *) "test_response_gram_os_linear failed: Incorrect "// &
                "rebased open-shell linear response Gram matrix."
            test_response_gram_os_linear = .false.
        end if

    end function test_response_gram_os_linear

    logical(c_bool) function test_apply_ms_sr1_skip() bind(C)
        !
        ! this function tests the subroutine which discards history directions failing
        ! the multisecant SR1 skipping criterion
        !
        use otr_arh, only: apply_ms_sr1_skip, ms_sr1_skip_thresh
        use otr_common_unit_tests, only: identity_matrix

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
        use otr_arh, only: spectral_to_dense

        integer(ip), parameter :: n = 2

        real(rp) :: eigvecs(n, n), inv_eigvals(n), expected(n, n)
        real(rp), allocatable :: mat(:, :)

        ! assume tests pass
        test_spectral_to_dense = .true.

        ! initialize an orthonormal but non-symmetric eigenvector matrix, so that a
        ! transposed reconstruction would be caught, and distinct inverted eigenvalues
        eigvecs = reshape([0.6_rp, 0.8_rp, -0.8_rp, 0.6_rp], [n, n])
        inv_eigvals = [0.25_rp, 4.0_rp]

        ! initialize the expected matrix 0.25 * v1 v1^T + 4 * v2 v2^T
        expected = reshape([2.65_rp, -1.8_rp, -1.8_rp, 1.6_rp], [n, n])

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
        use otr_arh, only: density_in_history, arh_object, arh_oao_type
        use otr_common_test_reference, only: n_ao, n_particle

        real(rp), parameter :: scale = 1e6_rp

        real(rp) :: dm(n_ao, n_ao, n_particle), other(n_ao, n_ao, n_particle)

        ! assume tests pass
        test_density_in_history = .true.

        ! generate two different densities
        call random_number(dm)
        call random_number(other)
        other = other + 1.0_rp

        ! determine if no density is reported before the history exists
        allocate(arh_oao_type :: arh_object)
        if (density_in_history(dm)) then
            write(stderr, *) "test_density_in_history failed: Density matrix "// &
                "reported without a history."
            test_density_in_history = .false.
        end if

        ! determine if a density is found in a history holding it besides another one,
        ! also when it differs by much less than numerical zero relative to its size,
        ! but not when it differs by much more
        allocate(arh_object%dm_list(n_ao, n_ao, n_particle, 2))
        arh_object%dm_list(:, :, :, 1) = other
        arh_object%dm_list(:, :, :, 2) = dm
        if (.not. density_in_history(dm)) then
            write(stderr, *) "test_density_in_history failed: Density matrix held "// &
                "in history not found."
            test_density_in_history = .false.
        end if
        if (.not. density_in_history(dm + 0.1_rp * numerical_zero)) then
            write(stderr, *) "test_density_in_history failed: Density matrix "// &
                "differing by much less than numerical zero not found."
            test_density_in_history = .false.
        end if
        if (density_in_history(dm + 10.0_rp * numerical_zero)) then
            write(stderr, *) "test_density_in_history failed: Density matrix "// &
                "differing by much more than numerical zero found."
            test_density_in_history = .false.
        end if

        ! determine if the comparison is relative to the size of the density, where the
        ! densities are scaled after perturbing them so that the differences scale with
        ! them instead of vanishing in their precision
        arh_object%dm_list(:, :, :, 2) = scale * dm
        if (.not. density_in_history(scale * (dm + 0.1_rp * numerical_zero))) then
            write(stderr, *) "test_density_in_history failed: Large density matrix "// &
                "differing by much less than numerical zero relative to its size "// &
                "not found."
            test_density_in_history = .false.
        end if
        if (density_in_history(scale * (dm + 10.0_rp * numerical_zero))) then
            write(stderr, *) "test_density_in_history failed: Large density matrix "// &
                "differing by much more than numerical zero relative to its size found."
            test_density_in_history = .false.
        end if

        ! deallocate ARH object
        deallocate(arh_object)

    end function test_density_in_history

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
            write(stderr, *) "test_prepend failed: Incorrect list dimensions after "// &
                "prepending to empty list."
            test_prepend = .false.
        end if
        if (norm2(list(:, :, :, 1) - new_array) > tol) then
            write(stderr, *) "test_prepend failed: Incorrect list values after "// &
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

end module otr_arh_unit_tests
