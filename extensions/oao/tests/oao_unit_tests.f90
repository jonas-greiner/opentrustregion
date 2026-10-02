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

contains

    subroutine setup_identity_oao_eigenbasis(eigval)
        !
        ! this subroutine sets up the module-global OAO object with the shared
        ! dimensions and an up-to-date eigendecomposition of the static part of the
        ! Hessian in the OAO basis itself with the given eigenvalue for every
        ! eigenvector, so that all eigenvalue pairs are the same and routines using the
        ! eigendecomposition reduce to a scaling
        !
        use otr_oao, only: oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_common_unit_tests, only: identity_matrix

        real(rp), intent(in) :: eigval

        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        oao_object%hess_eigvecs = spread(identity_matrix(n_ao), 3, n_particle)
        allocate(oao_object%hess_eigvals(n_ao, n_particle), source=eigval)
        oao_object%hess_eigen_stale = .false.

    end subroutine setup_identity_oao_eigenbasis

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
        use otr_common_unit_tests, only: identity_matrix

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

    function ref_hess_x_static_oao(x_full, fock_oo, fock_vv) result(hess_x_full)
        !
        ! this function independently reproduces the static part of the Hessian in the
        ! OAO basis, (F_vv - F_oo) X - h.c. for the unpacked trial vector X of every
        ! particle channel, before its projection, scaling and packing
        !
        real(rp), intent(in) :: x_full(:, :, :), fock_oo(:, :, :), fock_vv(:, :, :)
        real(rp) :: hess_x_full(size(x_full, 1), size(x_full, 2), size(x_full, 3))

        integer(ip) :: i

        do i = 1, size(x_full, 3, kind=ip)
            hess_x_full(:, :, i) = matmul(fock_vv(:, :, i) - fock_oo(:, :, i), &
                                          x_full(:, :, i))
            hess_x_full(:, :, i) = hess_x_full(:, :, i) - &
                                   transpose(hess_x_full(:, :, i))
        end do

    end function ref_hess_x_static_oao

    function ref_hess_x_oao(x_full, response, dm_oao, fock_oo, fock_vv, n_param) &
        result(hess_x)
        !
        ! this function assembles a Hessian linear transformation in the OAO basis from
        ! a given response contribution
        !
        real(rp), intent(in) :: x_full(:, :, :), response(:, :, :), dm_oao(:, :, :), &
                                fock_oo(:, :, :), fock_vv(:, :, :)
        integer(ip), intent(in) :: n_param
        real(rp) :: hess_x(n_param)

        real(rp) :: hess_x_full(size(x_full, 1), size(x_full, 2), size(x_full, 3))

        ! project the combined static and response contributions onto the
        ! occupied-virtual and virtual-occupied subspace, matching the production code
        hess_x_full = ref_project_asymm( &
            ref_hess_x_static_oao(x_full, fock_oo, fock_vv) + response, dm_oao)

        ! pack Hessian linear transformation
        hess_x = merge(4.0_rp, 2.0_rp, size(x_full, 3) == 1) * &
                 ref_pack_asymm(hess_x_full, n_param)

    end function ref_hess_x_oao

    subroutine ref_diagonalize_static_part(fock_oo, fock_vv, eigvecs, eigvals)
        !
        ! this subroutine independently diagonalizes the static part of the Hessian for
        ! each particle channel
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

    function ref_rotate_eigenbasis_oao(vector, eigvecs, n_particle, n_ao, &
                                       to_eigenbasis) result(rotated)
        !
        ! this function independently rotates a packed antisymmetric vector into or out
        ! of a given eigenbasis, by rotating the unpacked antisymmetric matrix X of
        ! every particle channel as U^T X U or as U X U^T, reproducing the
        ! corresponding OAO routines so that tests of routines which rotate into or out
        ! of the eigenbasis internally do not depend on them
        !
        real(rp), intent(in) :: vector(:), eigvecs(:, :, :)
        integer(ip), intent(in) :: n_particle, n_ao
        logical, intent(in) :: to_eigenbasis
        real(rp) :: rotated(size(vector))

        real(rp) :: full(n_ao, n_ao, n_particle), rotated_full(n_ao, n_ao, n_particle)
        integer(ip) :: i

        full = ref_unpack_asymm(vector, n_particle, n_ao)
        do i = 1, n_particle
            if (to_eigenbasis) then
                rotated_full(:, :, i) = matmul(transpose(eigvecs(:, :, i)), &
                                               matmul(full(:, :, i), eigvecs(:, :, i)))
            else
                rotated_full(:, :, i) = matmul(eigvecs(:, :, i), matmul( &
                    full(:, :, i), transpose(eigvecs(:, :, i))))
            end if
        end do
        rotated = ref_pack_asymm(rotated_full, size(vector, kind=ip))

    end function ref_rotate_eigenbasis_oao

    function ref_hess_eigval_pairs_oao(eigvals, n_particle, n_ao, n_param) &
        result(eigval_pairs)
        !
        ! this function independently constructs the pairwise sums of the static
        ! Hessian part eigenvalues, reproducing the corresponding OAO routine so that
        ! tests of routines which build the eigenvalue pairs internally do not depend
        ! on it
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

    end function ref_hess_eigval_pairs_oao

    function check_rotate_hess_eigenbasis_oao(to_eigenbasis) result(passed)
        !
        ! this function checks the rotation of a random packed antisymmetric vector
        ! into or out of the eigenbasis of the static part of the Hessian for the OAO
        ! basis
        !
        use otr_oao, only: oao_type
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix

        logical, intent(in) :: to_eigenbasis
        logical :: passed

        type(oao_type) :: oao
        real(rp) :: eigvecs(n_ao, n_ao, n_particle), vector(n_param), rotated(n_param)
        integer(ip) :: j
        character(len=:), allocatable :: test_name

        ! assume test passes
        passed = .true.
        test_name = trim(merge("test_rotate_to_hess_eigenbasis_oao  ", &
                               "test_rotate_from_hess_eigenbasis_oao", to_eigenbasis))

        ! generate random orthogonal eigenvectors for every particle channel
        do j = 1, n_particle
            eigvecs(:, :, j) = generate_random_orthogonal_matrix(n_ao)
        end do
        call random_number(vector)

        ! set up the OAO object
        oao%n_ao = n_ao
        oao%n_particle = n_particle
        oao%n_param = n_param
        oao%hess_eigvecs = eigvecs

        ! call routine and determine if the rotation matches
        if (to_eigenbasis) then
            rotated = oao%rotate_to_hess_eigenbasis(vector)
        else
            rotated = oao%rotate_from_hess_eigenbasis(vector)
        end if
        if (norm2(rotated - ref_rotate_eigenbasis_oao(vector, eigvecs, n_particle, &
                                                      n_ao, to_eigenbasis)) > tol) then
            write(stderr, *) test_name//" failed: Incorrect rotation."
            passed = .false.
        end if

    end function check_rotate_hess_eigenbasis_oao

    logical(c_bool) function test_oao_factory_cs() bind(C)
        !
        ! this function tests the subroutine which returns the modified OAO orbital
        ! updating function for the closed-shell case
        !
        use otr_oao, only: oao_factory_cs, oao_object, oao_settings_type, &
                           obj_func_oao_callback_ptr, update_orbs_oao_callback_ptr, &
                           precond_oao_callback_ptr
        use otr_common, only: evaluate_dm_cs_type
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_ao, n_occ
        use otr_common_unit_tests, only: identity_matrix, &
                                         generate_random_density_matrix, &
                                         mock_evaluate_dm_cs, mock_evaluate_dm_os

        integer(ip), parameter :: n_particle = 1

        real(rp), target :: dm_ao(n_ao, n_ao)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: error
        type(oao_settings_type) :: settings
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        type(solver_settings_type) :: solver_settings

        ! assume tests pass
        test_oao_factory_cs = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize density matrix and an orthonormal AO basis
        dm_ao = generate_random_density_matrix(n_ao, n_occ(1))
        ao_overlap = identity_matrix(n_ao)

        ! initialize callback function pointers
        evaluate_dm_funptr => mock_evaluate_dm_cs

        ! leave an OAO object from an open-shell calculation behind, whose density
        ! matrix evaluating function has to be cleared
        allocate(oao_object)
        oao_object%evaluate_dm_os => mock_evaluate_dm_os

        ! call routine and determine if an error is produced
        call oao_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                            obj_func_oao_funptr, update_orbs_oao_funptr, &
                            solver_settings, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_oao_factory_cs failed: Produced error."
            test_oao_factory_cs = .false.
            if (allocated(oao_object)) deallocate(oao_object)
            return
        end if

        ! determine if the OAO object points to the density matrix of the caller as its
        ! only particle channel, the remaining common setup is covered by the test of
        ! the common setup
        if (.not. associated(oao_object%dm_ao)) then
            write(stderr, *) "test_oao_factory_cs failed: Density matrix not "// &
                "associated."
            test_oao_factory_cs = .false.
            deallocate(oao_object)
            return
        end if
        if (any(shape(oao_object%dm_ao) /= [n_ao, n_ao, n_particle])) then
            write(stderr, *) "test_oao_factory_cs failed: Density matrix not "// &
                "associated as a single particle channel."
            test_oao_factory_cs = .false.
            deallocate(oao_object)
            return
        end if
        dm_ao(1, 1) = dm_ao(1, 1) + 1.0_rp
        if (norm2(oao_object%dm_ao(:, :, 1) - dm_ao) > tol) then
            write(stderr, *) "test_oao_factory_cs failed: Density matrix of the "// &
                "caller not taken over."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(oao_object%evaluate_dm_cs, mock_evaluate_dm_cs)) then
            write(stderr, *) "test_oao_factory_cs failed: Density matrix "// &
                "evaluating function not stored correctly."
            test_oao_factory_cs = .false.
        end if
        if (associated(oao_object%evaluate_dm_os)) then
            write(stderr, *) "test_oao_factory_cs failed: Open-shell density "// &
                "matrix evaluating function of a previous calculation kept."
            test_oao_factory_cs = .false.
        end if

        ! determine if returned function pointers point to the correct routines
        if (.not. associated(obj_func_oao_funptr, obj_func_oao_callback_ptr)) then
            write(stderr, *) "test_oao_factory_cs failed: Returned objective "// &
                "function is wrong."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(update_orbs_oao_funptr, update_orbs_oao_callback_ptr)) then
            write(stderr, *) "test_oao_factory_cs failed: Returned orbital "// &
                "updating function is wrong."
            test_oao_factory_cs = .false.
        end if

        ! determine if the OAO routines are wired into the solver settings
        if (.not. associated(solver_settings%precond, precond_oao_callback_ptr)) then
            write(stderr, *) "test_oao_factory_cs failed: OAO routines not wired "// &
                "into solver settings."
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
                           obj_func_oao_callback_ptr, update_orbs_oao_callback_ptr, &
                           precond_oao_callback_ptr
        use otr_common, only: evaluate_dm_os_type
        use opentrustregion, only: obj_func_type, update_orbs_type, solver_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_ao, n_particle, n_occ
        use otr_common_unit_tests, only: identity_matrix, &
                                         generate_random_density_matrix, &
                                         mock_evaluate_dm_cs, mock_evaluate_dm_os

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: i, error
        type(oao_settings_type) :: settings
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        type(solver_settings_type) :: solver_settings

        ! assume tests pass
        test_oao_factory_os = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize density matrices and an orthonormal AO basis
        do i = 1, n_particle
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do
        ao_overlap = identity_matrix(n_ao)

        ! initialize callback function pointers
        evaluate_dm_funptr => mock_evaluate_dm_os

        ! leave an OAO object from a closed-shell calculation behind, whose density
        ! matrix evaluating function has to be cleared
        allocate(oao_object)
        oao_object%evaluate_dm_cs => mock_evaluate_dm_cs

        ! call routine and determine if an error is produced
        call oao_factory_os(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                            obj_func_oao_funptr, update_orbs_oao_funptr, &
                            solver_settings, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_oao_factory_os failed: Produced error."
            test_oao_factory_os = .false.
            if (allocated(oao_object)) deallocate(oao_object)
            return
        end if

        ! determine if the OAO object points to the density matrix of the caller, the
        ! remaining common setup is covered by the test of the common setup
        if (.not. associated(oao_object%dm_ao, dm_ao)) then
            write(stderr, *) "test_oao_factory_os failed: Density matrix of the "// &
                "caller not taken over."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(oao_object%evaluate_dm_os, mock_evaluate_dm_os)) then
            write(stderr, *) "test_oao_factory_os failed: Density matrix "// &
                "evaluating function not stored correctly."
            test_oao_factory_os = .false.
        end if
        if (associated(oao_object%evaluate_dm_cs)) then
            write(stderr, *) "test_oao_factory_os failed: Closed-shell density "// &
                "matrix evaluating function of a previous calculation kept."
            test_oao_factory_os = .false.
        end if

        ! determine if returned function pointers point to the correct routines
        if (.not. associated(obj_func_oao_funptr, obj_func_oao_callback_ptr)) then
            write(stderr, *) "test_oao_factory_os failed: Returned objective "// &
                "function is wrong."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(update_orbs_oao_funptr, update_orbs_oao_callback_ptr)) then
            write(stderr, *) "test_oao_factory_os failed: Returned orbital "// &
                "updating function is wrong."
            test_oao_factory_os = .false.
        end if

        ! determine if the OAO routines are wired into the solver settings
        if (.not. associated(solver_settings%precond, precond_oao_callback_ptr)) then
            write(stderr, *) "test_oao_factory_os failed: OAO routines not wired "// &
                "into solver settings."
            test_oao_factory_os = .false.
        end if
        deallocate(oao_object)

    end function test_oao_factory_os

    logical(c_bool) function test_oao_factory_common() bind(C)
        !
        ! this function tests the subroutine which performs the common OAO
        ! initialization operations, which sets the OAO object up anew unless it was
        ! already set up for the same dimensions, for the closed- and the open-shell
        ! case
        !
        use otr_common, only: orbital_settings_type
        use otr_oao, only: oao_factory_common, oao_object
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, operator(==)
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_unit_tests, only: generate_random_density_matrix, shell_names, &
                                         mock_get_response_cs, mock_get_response_os, &
                                         mock_evaluate_dm_cs, mock_evaluate_dm_os

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle_ref), &
                            dm_ao_new(n_ao, n_ao, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), overlap_diag(n_ao), s_sqrt(n_ao, n_ao), &
                    s_inv_sqrt(n_ao, n_ao)
        integer(ip) :: n_particle, leftover, n_ao_old, n_particle_old, i, error
        type(orbital_settings_type) :: settings
        character(len=:), allocatable :: shell, leftover_case
        character(len=39), parameter :: leftover_names(2) = &
            [character(len=39) :: "a different number of AOs", &
             "a different number of particle channels"]

        ! assume tests pass
        test_oao_factory_common = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize a random diagonal AO overlap matrix, whose square roots are
        ! diagonal with the square roots of its elements, so that the density matrix in
        ! the OAO basis differs from the one in the AO basis
        call random_number(overlap_diag)
        overlap_diag = 0.5_rp + overlap_diag
        ao_overlap = 0.0_rp
        s_sqrt = 0.0_rp
        s_inv_sqrt = 0.0_rp
        do i = 1, n_ao
            ao_overlap(i, i) = overlap_diag(i)
            s_sqrt(i, i) = sqrt(overlap_diag(i))
            s_inv_sqrt(i, i) = 1.0_rp / sqrt(overlap_diag(i))
        end do

        do n_particle = 1, n_particle_ref
            shell = trim(shell_names(n_particle))

            ! initialize random density matrices
            do i = 1, n_particle_ref
                dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
                dm_ao_new(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
            end do

            ! leave OAO objects from calculations which differ only in the number of
            ! AOs or of particle channels behind and determine if each is set up anew
            do leftover = 1, size(leftover_names, kind=ip)
                n_ao_old = merge(n_ao + 1, n_ao, leftover == 1)
                n_particle_old = merge(3 - n_particle, n_particle, leftover == 2)
                leftover_case = shell//" object left from a calculation with "// &
                                trim(leftover_names(leftover))
                if (allocated(oao_object)) deallocate(oao_object)
                allocate(oao_object)
                oao_object%n_ao = n_ao_old
                oao_object%n_particle = n_particle_old
                allocate(oao_object%s_sqrt(n_ao_old, n_ao_old), &
                         oao_object%s_inv_sqrt(n_ao_old, n_ao_old), &
                         oao_object%fock_oo(n_ao_old, n_ao_old, n_particle_old), &
                         oao_object%grad(1), oao_object%h_diag(1))
                oao_object%s_sqrt = 0.0_rp
                oao_object%s_inv_sqrt = 0.0_rp
                call oao_factory_common(dm_ao(:, :, :n_particle), ao_overlap, &
                                        n_particle, n_ao, error, settings)
                if (error /= 0) then
                    write(stderr, *) "test_oao_factory_common failed: Produced "// &
                        "error for the "//leftover_case//"."
                    test_oao_factory_common = .false.
                end if
                if (.not. check_oao_object(dm_ao(:, :, :n_particle), s_sqrt, &
                                           leftover_case)) &
                    test_oao_factory_common = .false.
                if (norm2(oao_object%s_inv_sqrt - s_inv_sqrt) > tol) then
                    write(stderr, *) "test_oao_factory_common failed: Inverse "// &
                        "square root of the overlap matrix not computed for the "// &
                        leftover_case//"."
                    test_oao_factory_common = .false.
                end if
            end do

            ! leave evaluated quantities with their response functions, the density
            ! matrix evaluating function, marked square roots of the overlap matrix and
            ! only the gradient allocated behind and call routine again with new
            ! settings and a new starting density, and determine if the object is kept
            ! while the new settings and density are taken over, its evaluated state is
            ! discarded, the density matrix evaluating function, which only the
            ! factories set, is kept and the Hessian diagonal is allocated
            oao_object%evaluation_stale = .false.
            oao_object%response_stale = .false.
            oao_object%hess_eigen_stale = .false.
            oao_object%get_response_cs => mock_get_response_cs
            oao_object%get_response_os => mock_get_response_os
            if (n_particle == 1) then
                oao_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                oao_object%evaluate_dm_os => mock_evaluate_dm_os
            end if
            oao_object%s_sqrt = 2.0_rp * s_sqrt
            oao_object%s_inv_sqrt = 0.5_rp * s_inv_sqrt
            deallocate(oao_object%h_diag)
            settings%verbose = settings%verbose + 1
            call oao_factory_common(dm_ao_new(:, :, :n_particle), ao_overlap, &
                                    n_particle, n_ao, error, settings)
            if (error /= 0) then
                write(stderr, *) "test_oao_factory_common failed: Produced error "// &
                    "for the reused "//shell//" object."
                test_oao_factory_common = .false.
            end if
            if (.not. check_oao_object(dm_ao_new(:, :, :n_particle), 2.0_rp * s_sqrt, &
                                       "reused "//shell//" object")) &
                test_oao_factory_common = .false.
            if (norm2(oao_object%s_inv_sqrt - 0.5_rp * s_inv_sqrt) > tol) then
                write(stderr, *) "test_oao_factory_common failed: Square roots of "// &
                    "the overlap matrix recomputed for the reused "//shell//" object."
                test_oao_factory_common = .false.
            end if
            if (.not. (associated(oao_object%evaluate_dm_cs, mock_evaluate_dm_cs) .or. &
                       associated(oao_object%evaluate_dm_os, mock_evaluate_dm_os))) then
                write(stderr, *) "test_oao_factory_common failed: Density matrix "// &
                    "evaluating function not kept for the reused "//shell//" object."
                test_oao_factory_common = .false.
            end if

            ! leave an OAO object behind whose previous setup failed before the square
            ! roots of the overlap matrix were computed and determine if it is set up
            ! anew although the dimensions match
            deallocate(oao_object%s_inv_sqrt)
            call oao_factory_common(dm_ao(:, :, :n_particle), ao_overlap, n_particle, &
                                    n_ao, error, settings)
            if (error /= 0) then
                write(stderr, *) "test_oao_factory_common failed: Produced error "// &
                    "for the "//shell//" object whose previous setup failed."
                test_oao_factory_common = .false.
            end if
            if (.not. check_oao_object(dm_ao(:, :, :n_particle), s_sqrt, &
                                       shell//" object whose previous setup failed")) &
                test_oao_factory_common = .false.
            if (.not. allocated(oao_object%s_inv_sqrt)) then
                write(stderr, *) "test_oao_factory_common failed: Inverse square "// &
                    "root of the overlap matrix not computed for the "//shell// &
                    " object whose previous setup failed."
                test_oao_factory_common = .false.
            else if (norm2(oao_object%s_inv_sqrt - s_inv_sqrt) > tol) then
                write(stderr, *) "test_oao_factory_common failed: Incorrect "// &
                    "inverse square root of the overlap matrix for the "//shell// &
                    " object whose previous setup failed."
                test_oao_factory_common = .false.
            end if

            ! call routine for an invalid number of particle channels and determine if
            ! the sanity check rejects it before the object is changed
            oao_object%hess_eigen_stale = .false.
            call oao_factory_common(dm_ao(:, :, :n_particle), ao_overlap, 3_ip, n_ao, &
                                    error, settings)
            if (error == 0) then
                write(stderr, *) "test_oao_factory_common failed: Error not thrown "// &
                    "for an invalid number of particle channels for the "//shell// &
                    " object."
                test_oao_factory_common = .false.
            end if
            if (oao_object%hess_eigen_stale) then
                write(stderr, *) "test_oao_factory_common failed: OAO object "// &
                    "changed for an invalid number of particle channels for the "// &
                    shell//" object."
                test_oao_factory_common = .false.
            end if

            ! deallocate OAO object
            deallocate(oao_object)
        end do

    contains

        function check_oao_object(caller_dm_ao, expected_s_sqrt, case_name) &
            result(passed)
            !
            ! this function checks if the OAO object is set up for the given density
            ! matrix of the caller and the number of AOs and particle channels of the
            ! current shell, with the density matrix transformed to the OAO basis with
            ! the given square root of the overlap matrix, the evaluated state with its
            ! response functions discarded and the settings stored
            !
            real(rp), intent(in), target :: caller_dm_ao(:, :, :)
            real(rp), intent(in) :: expected_s_sqrt(:, :)
            character(len=*), intent(in) :: case_name
            logical :: passed

            integer(ip) :: n_param, k

            ! assume test passes
            passed = .true.

            ! the checks below compare arrays of these dimensions
            n_param = n_particle * n_ao * (n_ao - 1) / 2
            if (any([oao_object%n_ao, oao_object%n_particle, oao_object%n_param] /= &
                    [n_ao, n_particle, n_param])) then
                write(stderr, *) "test_oao_factory_common failed: Incorrect "// &
                    "dimensions for the "//case_name//"."
                passed = .false.
                return
            end if

            ! determine if the density matrix of the caller, its transformation to the
            ! OAO basis, the Fock matrix blocks, the gradient and Hessian diagonal, the
            ! discarded evaluated state and the settings are set up
            if (.not. associated(oao_object%dm_ao, caller_dm_ao)) then
                write(stderr, *) "test_oao_factory_common failed: Density matrix "// &
                    "of the caller not taken over for the "//case_name//"."
                passed = .false.
            end if
            do k = 1, n_particle
                if (norm2(oao_object%dm_oao(:, :, k) - matmul(expected_s_sqrt, matmul( &
                    caller_dm_ao(:, :, k), expected_s_sqrt))) > tol) then
                    write(stderr, *) "test_oao_factory_common failed: Density "// &
                        "matrix not transformed to the OAO basis for the "// &
                        case_name//"."
                    passed = .false.
                end if
            end do
            if (any([size(oao_object%fock_oo, 3), size(oao_object%fock_vv, 3)] /= &
                    n_particle)) then
                write(stderr, *) "test_oao_factory_common failed: Fock matrix "// &
                    "blocks not allocated for every particle channel for the "// &
                    case_name//"."
                passed = .false.
            end if
            if (.not. (allocated(oao_object%grad) .and. allocated(oao_object%h_diag))) &
                then
                write(stderr, *) "test_oao_factory_common failed: Gradient or "// &
                    "Hessian diagonal not allocated for the "//case_name//"."
                passed = .false.
            else if (any([size(oao_object%grad), size(oao_object%h_diag)] /= n_param)) &
                then
                write(stderr, *) "test_oao_factory_common failed: Gradient or "// &
                    "Hessian diagonal not allocated with the number of parameters "// &
                    "for the "//case_name//"."
                passed = .false.
            end if
            if (.not. (oao_object%evaluation_stale .and. &
                       oao_object%response_stale .and. oao_object%hess_eigen_stale)) &
                then
                write(stderr, *) "test_oao_factory_common failed: Evaluated state "// &
                    "not discarded for the "//case_name//"."
                passed = .false.
            end if
            if (associated(oao_object%get_response_cs) .or. &
                associated(oao_object%get_response_os)) then
                write(stderr, *) "test_oao_factory_common failed: Response "// &
                    "function kept for the "//case_name//"."
                passed = .false.
            end if
            if (.not. (oao_object%settings == settings)) then
                write(stderr, *) "test_oao_factory_common failed: Settings not "// &
                    "stored for the "//case_name//"."
                passed = .false.
            end if

        end function check_oao_object

    end function test_oao_factory_common

    logical(c_bool) function test_oao_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the OAO
        ! input parameters
        !
        use otr_oao, only: oao_settings_type, oao_sanity_check
        use opentrustregion_unit_tests, only: setup_settings

        type(oao_settings_type) :: settings
        integer(ip) :: n_particle, error

        ! assume tests pass
        test_oao_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if a positive number of AOs is accepted for one and for two particle
        ! channels
        do n_particle = 1, 2
            call oao_sanity_check(settings, n_particle, 1_ip, error)
            if (error /= 0) then
                write(stderr, *) "test_oao_sanity_check failed: Error thrown for "// &
                    "valid dimensions."
                test_oao_sanity_check = .false.
            end if
        end do

        ! check if a vanishing number of AOs is rejected
        call oao_sanity_check(settings, 1_ip, 0_ip, error)
        if (error == 0) then
            write(stderr, *) "test_oao_sanity_check failed: Error not thrown for "// &
                "vanishing number of AOs."
            test_oao_sanity_check = .false.
        end if

        ! check if a number of particle channels other than one or two is rejected
        do n_particle = 0, 3, 3
            call oao_sanity_check(settings, n_particle, 1_ip, error)
            if (error == 0) then
                write(stderr, *) "test_oao_sanity_check failed: Error not thrown "// &
                    "for invalid number of particle channels."
                test_oao_sanity_check = .false.
            end if
        end do

    end function test_oao_sanity_check

    logical(c_bool) function test_oao_set_solver_settings() bind(C)
        !
        ! this function tests the subroutine which wires the OAO preconditioners,
        ! projection and extra trial vectors into the solver settings
        !
        use otr_oao, only: oao_set_solver_settings, precond_oao_callback_ptr, &
                           precond_pd_oao_callback_ptr, project_oao_callback_ptr, &
                           get_extra_trial_vectors_oao_callback_ptr
        use opentrustregion, only: solver_settings_type, default_solver_settings
        use test_reference, only: ref_settings, assignment(=), operator(/=)

        type(solver_settings_type) :: solver_settings
        integer(ip) :: i_case, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_oao_set_solver_settings = .true.

        ! wire the OAO routines into uninitialized settings, which have to be
        ! initialized to their defaults, and into settings initialized to the reference
        ! values, which have to be kept
        do i_case = 1, 2
            if (i_case == 1) then
                case_name = "for uninitialized settings"
            else
                case_name = "for initialized settings"
                solver_settings = ref_settings
            end if
            call oao_set_solver_settings(solver_settings, error)
            if (error /= 0) then
                write(stderr, *) "test_oao_set_solver_settings failed: Produced "// &
                    "error "//case_name//"."
                test_oao_set_solver_settings = .false.
            end if
            if (.not. solver_settings%initialized) then
                write(stderr, *) "test_oao_set_solver_settings failed: Settings "// &
                    "not initialized "//case_name//"."
                test_oao_set_solver_settings = .false.
            end if
            if (i_case == 1 .and. solver_settings /= default_solver_settings) then
                write(stderr, *) "test_oao_set_solver_settings failed: Settings "// &
                    "not set to their defaults "//case_name//"."
                test_oao_set_solver_settings = .false.
            end if
            if (i_case == 2 .and. solver_settings /= ref_settings) then
                write(stderr, *) "test_oao_set_solver_settings failed: Settings "// &
                    "not kept "//case_name//"."
                test_oao_set_solver_settings = .false.
            end if
            if (.not. (associated( &
                solver_settings%precond, precond_oao_callback_ptr) .and. associated( &
                    solver_settings%precond_pd, precond_pd_oao_callback_ptr) .and. &
                associated(solver_settings%project, project_oao_callback_ptr) .and. &
                associated(solver_settings%get_extra_trial_vectors, &
                           get_extra_trial_vectors_oao_callback_ptr))) then
                write(stderr, *) "test_oao_set_solver_settings failed: OAO "// &
                    "routines not wired into solver settings "//case_name//"."
                test_oao_set_solver_settings = .false.
            end if
            if (.not. ( &
                associated(solver_settings%stability_settings%precond, &
                           precond_oao_callback_ptr) .and. &
                associated(solver_settings%stability_settings%project, &
                           project_oao_callback_ptr) .and. &
                associated(solver_settings%stability_settings%get_extra_trial_vectors, &
                           get_extra_trial_vectors_oao_callback_ptr))) then
                write(stderr, *) "test_oao_set_solver_settings failed: OAO "// &
                    "routines not wired into stability check settings "//case_name//"."
                test_oao_set_solver_settings = .false.
            end if
        end do

    end function test_oao_set_solver_settings

    logical(c_bool) function test_obj_func_oao_callback() bind(C)
        !
        ! this function tests the function which defines the energy evaluation in the
        ! OAO basis
        !
        use otr_common_unit_tests, only: &
            identity_matrix, generate_random_density_matrix, &
            generate_random_symm_matrix, mock_requests, shell_names, &
            mock_evaluate_dm_cs, mock_evaluate_dm_os
        use otr_oao, only: obj_func_oao_callback, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle, n_occ

        real(rp), parameter :: angle = 0.3_rp

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), s_inv_sqrt(n_ao, n_ao), &
                    kappa(n_param), rotation(n_ao, n_ao), energy, expected_energy
        integer(ip) :: i, i_shell, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_obj_func_oao_callback = .true.

        ! set up the OAO object with a non-orthonormal AO basis, so that the density
        ! matrix has to be transformed to the AO basis before it is evaluated, and a
        ! different square root of the overlap matrix, which must not be used for it
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%n_param = n_param
        do i = 1, n_particle
            dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do
        s_inv_sqrt = identity_matrix(n_ao) + 0.1_rp * generate_random_symm_matrix(n_ao)
        oao_object%dm_oao = dm_oao
        oao_object%s_inv_sqrt = s_inv_sqrt
        oao_object%s_sqrt = identity_matrix(n_ao)

        ! loop over the closed-shell and the open-shell case, where the closed-shell
        ! case only passes on the first spin channel of the density matrix, so that
        ! only that channel enters the energy
        kappa = 0.0_rp
        do i_shell = 1, 2
            case_name = trim(shell_names(i_shell))
            if (i_shell == 1) then
                oao_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                oao_object%evaluate_dm_cs => null()
                oao_object%evaluate_dm_os => mock_evaluate_dm_os
            end if

            ! call routine without an orbital rotation and determine if only the energy
            ! of the density matrix evaluating function is requested and returned for
            ! the density matrix in the AO basis
            mock_requests = [integer(ip) :: ]
            energy = obj_func_oao_callback(kappa, error)
            if (error /= 0) then
                write(stderr, *) "test_obj_func_oao_callback failed: Produced "// &
                    "error for the "//case_name//" case."
                test_obj_func_oao_callback = .false.
            end if
            expected_energy = 0.0_rp
            do i = 1, i_shell
                expected_energy = expected_energy + sum(matmul(s_inv_sqrt, matmul( &
                    dm_oao(:, :, i), s_inv_sqrt)))
            end do
            if (abs(energy - expected_energy) > tol) then
                write(stderr, *) "test_obj_func_oao_callback failed: Incorrect "// &
                    "energy for the "//case_name//" case."
                test_obj_func_oao_callback = .false.
            end if
            if (size(mock_requests) /= 1) then
                write(stderr, *) "test_obj_func_oao_callback failed: Density "// &
                    "matrix evaluating function not called exactly once for the "// &
                    case_name//" case."
                test_obj_func_oao_callback = .false.
            end if
            if (any(mock_requests /= 0)) then
                write(stderr, *) "test_obj_func_oao_callback failed: Incorrect "// &
                    "outputs requested from density matrix evaluating function for "// &
                    "the "//case_name//" case."
                test_obj_func_oao_callback = .false.
            end if
        end do

        ! call routine with an orbital rotation and determine if the current density
        ! matrix is left untouched
        call random_number(kappa)
        kappa = 0.1_rp * kappa
        energy = obj_func_oao_callback(kappa, error)
        if (error /= 0) then
            write(stderr, *) "test_obj_func_oao_callback failed: Produced error "// &
                "for an orbital rotation."
            test_obj_func_oao_callback = .false.
        end if
        if (norm2(oao_object%dm_oao - dm_oao) > tol) then
            write(stderr, *) "test_obj_func_oao_callback failed: Current density "// &
                "matrix changed by an orbital rotation."
            test_obj_func_oao_callback = .false.
        end if

        ! call routine with a rotation of the first pair of orbitals of the first spin
        ! channel, whose rotation matrix is a plane rotation, and determine if the
        ! energy of the rotated density matrix is returned; the closed form of the
        ! rotation catches an ignored or sign-flipped orbital rotation, which the
        ! checks above do not
        kappa = 0.0_rp
        kappa(1) = angle
        rotation = identity_matrix(n_ao)
        rotation(1:2, 1:2) = reshape([cos(angle), -sin(angle), &
                                      sin(angle), cos(angle)], [2, 2])
        energy = obj_func_oao_callback(kappa, error)
        if (error /= 0) then
            write(stderr, *) "test_obj_func_oao_callback failed: Produced error "// &
                "for a plane rotation."
            test_obj_func_oao_callback = .false.
        end if
        expected_energy = sum(matmul(s_inv_sqrt, matmul( &
            matmul(transpose(rotation), matmul(dm_oao(:, :, 1), rotation)), &
            s_inv_sqrt))) + sum(matmul(s_inv_sqrt, matmul(dm_oao(:, :, 2), s_inv_sqrt)))
        if (abs(energy - expected_energy) > tol) then
            write(stderr, *) "test_obj_func_oao_callback failed: Incorrect energy "// &
                "for a plane rotation."
            test_obj_func_oao_callback = .false.
        end if
        deallocate(oao_object)

    end function test_obj_func_oao_callback

    logical(c_bool) function test_update_orbs_oao_callback() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation in the OAO basis, for the closed-shell and the
        ! open-shell case
        !
        use otr_common_unit_tests, only: &
            identity_matrix, generate_random_density_matrix, &
            generate_random_symm_matrix, mock_fock_factor, mock_requests, shell_names, &
            mock_get_response_cs, mock_get_response_os, mock_evaluate_dm_cs, &
            mock_evaluate_dm_os, mock_evaluate_dm_failing_cs, &
            mock_evaluate_dm_failing_os
        use otr_oao, only: update_orbs_oao_callback, oao_object, hess_x_oao_callback_ptr
        use opentrustregion, only: hess_x_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_ao, n_occ, n_particle_ref => n_particle

        real(rp), allocatable, target :: dm_ao(:, :, :)
        real(rp), allocatable :: dm_oao(:, :, :), kappa(:), grad(:), h_diag(:)
        real(rp) :: s_inv_sqrt(n_ao, n_ao), fock_oao(n_ao, n_ao), func
        integer(ip) :: n_particle, n_param, i, error
        logical :: response_set
        character(len=:), allocatable :: case_name
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_update_orbs_oao_callback = .true.

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            n_param = n_particle * n_ao * (n_ao - 1) / 2

            ! set up the OAO object with a non-orthonormal AO basis, whose density
            ! matrix in the AO basis follows from the one in the OAO basis, so that the
            ! Fock matrix has to be transformed to the OAO basis
            s_inv_sqrt = identity_matrix(n_ao) + &
                         0.1_rp * generate_random_symm_matrix(n_ao)
            allocate(dm_ao(n_ao, n_ao, n_particle), dm_oao(n_ao, n_ao, n_particle), &
                     kappa(n_param), grad(n_param), h_diag(n_param))
            do i = 1, n_particle
                dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
                dm_ao(:, :, i) = matmul(s_inv_sqrt, matmul(dm_oao(:, :, i), s_inv_sqrt))
            end do
            allocate(oao_object)
            call setup_settings(oao_object%settings)
            oao_object%n_ao = n_ao
            oao_object%n_particle = n_particle
            oao_object%n_param = n_param
            oao_object%dm_ao => dm_ao
            oao_object%dm_oao = dm_oao
            oao_object%s_inv_sqrt = s_inv_sqrt
            allocate(oao_object%fock_oo(n_ao, n_ao, n_particle), &
                     oao_object%fock_vv(n_ao, n_ao, n_particle), &
                     oao_object%grad(n_param), oao_object%h_diag(n_param))
            if (n_particle == 1) then
                oao_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                oao_object%evaluate_dm_os => mock_evaluate_dm_os
            end if
            mock_requests = [integer(ip) :: ]

            ! call routine without an orbital rotation for a starting density which has
            ! not been evaluated yet, with the response marked current so that only the
            ! evaluation flag forces the evaluation, and determine if the energy of the
            ! density matrix in the AO basis, the static part of the Hessian from the
            ! Fock matrix transformed to the OAO basis, the gradient, the Hessian
            ! diagonal, the response function and the Hessian linear transformation are
            ! obtained
            kappa = 0.0_rp
            oao_object%response_stale = .false.
            if (.not. check_update_orbs_oao_stage(1_ip, "a starting density")) &
                test_update_orbs_oao_callback = .false.
            if (abs(func - sum(dm_ao)) > tol) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Incorrect "// &
                    "energy for the "//case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if
            do i = 1, n_particle
                fock_oao = matmul(s_inv_sqrt, matmul(mock_fock_factor(1) * &
                                                     dm_ao(:, :, i), s_inv_sqrt))
                if (norm2(oao_object%fock_oo(:, :, i) - matmul( &
                    dm_oao(:, :, i), matmul(fock_oao, dm_oao(:, :, i)))) > tol) then
                    write(stderr, *) "test_update_orbs_oao_callback failed: Static "// &
                        "Hessian part not built from the Fock matrix in the OAO "// &
                        "basis for the "//case_name//" case."
                    test_update_orbs_oao_callback = .false.
                end if
            end do
            if (norm2(grad - oao_object%grad) > tol) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Gradient "// &
                    "not returned for the "//case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if
            if (norm2(h_diag - oao_object%h_diag) > tol) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Hessian "// &
                    "diagonal not returned for the "//case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if
            if (n_particle == 1) then
                response_set = &
                    associated(oao_object%get_response_cs, mock_get_response_cs)
            else
                response_set = &
                    associated(oao_object%get_response_os, mock_get_response_os)
            end if
            if (.not. response_set) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Response "// &
                    "function not stored for the "//case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if
            if (.not. associated(hess_x_funptr, hess_x_oao_callback_ptr)) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Returned "// &
                    "Hessian linear transformation is wrong for the "//case_name// &
                    " case."
                test_update_orbs_oao_callback = .false.
            end if
            if (oao_object%evaluation_stale) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Evaluation "// &
                    "still marked stale after the quantities were computed for the "// &
                    case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if

            ! call routine again without an orbital rotation after clearing the outputs
            ! and setting a vanishing stored energy, and determine if the quantities of
            ! the evaluated density are returned without an evaluation, even for the
            ! vanishing energy
            oao_object%energy = 0.0_rp
            grad = 0.0_rp
            h_diag = 0.0_rp
            if (.not. check_update_orbs_oao_stage(1_ip, "an evaluated density")) &
                test_update_orbs_oao_callback = .false.
            if (abs(func) > tol) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Stored "// &
                    "energy not returned for an evaluated density for the "// &
                    case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if
            if (norm2(grad - oao_object%grad) > tol) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Stored "// &
                    "gradient not returned for an evaluated density for the "// &
                    case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if
            if (norm2(h_diag - oao_object%h_diag) > tol) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Stored "// &
                    "Hessian diagonal not returned for an evaluated density for "// &
                    "the "//case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if

            ! mark the response as stale, as an approximate-Hessian extension such as
            ! ARH would after moving the density without going through this routine,
            ! and call again without an orbital rotation, and determine if this forces
            ! an evaluation and clears the flag
            oao_object%response_stale = .true.
            if (.not. check_update_orbs_oao_stage(2_ip, "a stale response")) &
                test_update_orbs_oao_callback = .false.
            if (oao_object%response_stale) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Response "// &
                    "still marked stale after being recomputed for the "//case_name// &
                    " case."
                test_update_orbs_oao_callback = .false.
            end if

            ! mark the evaluation as stale, as the factory does for a new starting
            ! density, and call again without an orbital rotation, and determine if
            ! this forces an evaluation and clears the flag
            oao_object%evaluation_stale = .true.
            if (.not. check_update_orbs_oao_stage(3_ip, "a stale evaluation")) &
                test_update_orbs_oao_callback = .false.
            if (oao_object%evaluation_stale) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Evaluation "// &
                    "still marked stale after being recomputed for the "//case_name// &
                    " case."
                test_update_orbs_oao_callback = .false.
            end if

            ! call routine with an orbital rotation whose evaluation fails and
            ! determine if the error is passed on and the rotated density is marked as
            ! not evaluated
            if (n_particle == 1) then
                oao_object%evaluate_dm_cs => mock_evaluate_dm_failing_cs
            else
                oao_object%evaluate_dm_os => mock_evaluate_dm_failing_os
            end if
            kappa = 0.1_rp
            call update_orbs_oao_callback(kappa, func, grad, h_diag, hess_x_funptr, &
                                          error)
            if (error == 0) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Error of "// &
                    "the density matrix evaluating function not passed on for the "// &
                    case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if
            if (.not. oao_object%evaluation_stale) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Evaluation "// &
                    "not marked stale after a failed evaluation for the "//case_name// &
                    " case."
                test_update_orbs_oao_callback = .false.
            end if
            if (n_particle == 1) then
                oao_object%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                oao_object%evaluate_dm_os => mock_evaluate_dm_os
            end if

            ! call routine with an orbital rotation, keeping the density matrix in the
            ! OAO basis it starts from, and determine if the density is moved,
            ! consistently in the AO and the OAO basis, and evaluated, and if every
            ! evaluation requested the Fock matrix and the response function
            dm_oao = oao_object%dm_oao
            if (.not. check_update_orbs_oao_stage(4_ip, "an orbital rotation")) &
                test_update_orbs_oao_callback = .false.
            if (norm2(oao_object%dm_oao - dm_oao) < tol) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Density "// &
                    "matrix not moved by an orbital rotation for the "//case_name// &
                    " case."
                test_update_orbs_oao_callback = .false.
            end if
            do i = 1, n_particle
                if (norm2(dm_ao(:, :, i) - matmul(s_inv_sqrt, matmul( &
                    oao_object%dm_oao(:, :, i), s_inv_sqrt))) > tol) then
                    write(stderr, *) "test_update_orbs_oao_callback failed: "// &
                        "Density matrix not moved consistently in the AO and the "// &
                        "OAO basis for the "//case_name//" case."
                    test_update_orbs_oao_callback = .false.
                end if
            end do
            if (any(mock_requests /= 3)) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Incorrect "// &
                    "outputs requested from density matrix evaluating function for "// &
                    "the "//case_name//" case."
                test_update_orbs_oao_callback = .false.
            end if

            ! deallocate OAO object and case arrays
            deallocate(oao_object, dm_ao, dm_oao, kappa, grad, h_diag)
        end do

    contains

        function check_update_orbs_oao_stage(n_calls, stage) result(passed)
            !
            ! this function calls the OAO orbital updating function at the current
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
            call update_orbs_oao_callback(kappa, func, grad, h_diag, hess_x_funptr, &
                                          error)
            if (error /= 0) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Produced "// &
                    "error for "//stage//" for the "//case_name//" case."
                passed = .false.
            end if
            if (size(mock_requests) /= n_calls) then
                write(stderr, *) "test_update_orbs_oao_callback failed: Density "// &
                    "matrix evaluating function not called the expected number of "// &
                    "times for "//stage//" for the "//case_name//" case."
                passed = .false.
            end if

        end function check_update_orbs_oao_stage

    end function test_update_orbs_oao_callback

    logical(c_bool) function test_hess_x_oao_callback() bind(C)
        !
        ! this function tests the subroutine which defines the Hessian linear
        ! transformation in the OAO basis
        !
        use otr_common_unit_tests, only: &
            identity_matrix, generate_random_density_matrix, &
            generate_random_symm_matrix, mock_requests, mock_response_factor, &
            shell_names, mock_get_response_cs, mock_get_response_os, mock_evaluate_dm_os
        use otr_oao, only: hess_x_oao_callback, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_ao, n_occ, n_particle_ref => n_particle

        integer(ip) :: n_particle, n_param
        real(rp), target :: dm_oao(n_ao, n_ao, n_particle_ref)
        real(rp) :: fock_oo(n_ao, n_ao, n_particle_ref), &
                    fock_vv(n_ao, n_ao, n_particle_ref), s_inv_sqrt(n_ao, n_ao)
        real(rp), allocatable :: x(:), x_full(:, :, :), delta_dm(:, :, :), &
                                 response(:, :, :), hess_x(:), expected_hess_x(:)
        integer(ip) :: i, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_hess_x_oao_callback = .true.

        ! generate random density matrices, Fock matrix contributions and a
        ! non-orthonormal AO basis, so that the density matrix response has to be
        ! transformed to the AO basis and the Fock matrix response back to the OAO
        ! basis
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_occ(1))
        dm_oao(:, :, 2) = generate_random_density_matrix(n_ao, n_occ(2))
        do i = 1, 2
            fock_oo(:, :, i) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, i) = generate_random_symm_matrix(n_ao)
        end do
        s_inv_sqrt = identity_matrix(n_ao) + 0.1_rp * generate_random_symm_matrix(n_ao)

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%dm_oao = dm_oao
        oao_object%fock_oo = fock_oo
        oao_object%fock_vv = fock_vv
        oao_object%s_inv_sqrt = s_inv_sqrt
        oao_object%response_stale = .false.

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            n_param = n_particle * n_ao * (n_ao - 1) / 2
            oao_object%n_particle = n_particle
            oao_object%n_param = n_param
            if (n_particle == 1) then
                oao_object%get_response_cs => mock_get_response_cs
            else
                oao_object%get_response_cs => null()
                oao_object%get_response_os => mock_get_response_os
            end if

            ! the response is the mock response of the density matrix displacement of a
            ! random trial vector in the AO basis, transformed back to the OAO basis
            if (allocated(x)) deallocate(x, hess_x)
            allocate(x(n_param), hess_x(n_param))
            call random_number(x)
            x_full = ref_unpack_asymm(x, n_particle, n_ao)
            delta_dm = ref_project_symm(x_full, dm_oao(:, :, :n_particle))
            response = delta_dm
            do i = 1, n_particle
                response(:, :, i) = matmul( &
                    s_inv_sqrt, matmul(mock_response_factor * matmul( &
                        s_inv_sqrt, matmul(delta_dm(:, :, i), s_inv_sqrt)), s_inv_sqrt))
            end do
            expected_hess_x = ref_hess_x_oao( &
                x_full, response, dm_oao(:, :, :n_particle), &
                fock_oo(:, :, :n_particle), fock_vv(:, :, :n_particle), n_param)

            ! call routine and determine if values of resulting Hessian linear
            ! transformation match
            call hess_x_oao_callback(x, hess_x, error)
            if (error /= 0) then
                write(stderr, *) "test_hess_x_oao_callback failed: Produced error "// &
                    "for the "//case_name//" case."
                test_hess_x_oao_callback = .false.
            end if
            if (norm2(hess_x - expected_hess_x) > tol) then
                write(stderr, *) "test_hess_x_oao_callback failed: Incorrect "// &
                    "Hessian linear transformation for the "//case_name//" case."
                test_hess_x_oao_callback = .false.
            end if
        end do

        ! mark the response as stale, as ARH would after moving the density without
        ! going through the orbital updating routine, replace the response function by
        ! one of an outdated density and call again for the open-shell trial vector,
        ! and determine if this triggers the density matrix evaluating function to
        ! refresh the response, which is then used
        oao_object%dm_ao => dm_oao
        oao_object%evaluate_dm_os => mock_evaluate_dm_os
        oao_object%get_response_os => mock_get_stale_response_os
        oao_object%response_stale = .true.
        mock_requests = [integer(ip) :: ]
        call hess_x_oao_callback(x, hess_x, error)
        if (error /= 0) then
            write(stderr, *) "test_hess_x_oao_callback failed: Produced error "// &
                "while refreshing a stale response."
            test_hess_x_oao_callback = .false.
        end if
        if (size(mock_requests) /= 1) then
            write(stderr, *) "test_hess_x_oao_callback failed: Stale response was "// &
                "not refreshed."
            test_hess_x_oao_callback = .false.
        end if
        if (any(mock_requests /= 2)) then
            write(stderr, *) "test_hess_x_oao_callback failed: Incorrect outputs "// &
                "requested from density matrix evaluating function while "// &
                "refreshing a stale response."
            test_hess_x_oao_callback = .false.
        end if
        if (oao_object%response_stale) then
            write(stderr, *) "test_hess_x_oao_callback failed: Response still "// &
                "marked stale after being refreshed."
            test_hess_x_oao_callback = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write(stderr, *) "test_hess_x_oao_callback failed: Incorrect Hessian "// &
                "linear transformation with a refreshed response."
            test_hess_x_oao_callback = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    contains

        subroutine mock_get_stale_response_os(dm, response_out, error_out)
            !
            ! this subroutine is a mock response function for the open-shell case which
            ! belongs to an outdated density and returns a vanishing response
            !
            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out), target :: response_out(:, :, :)
            integer(ip), intent(out) :: error_out

            error_out = 0
            response_out = 0.0_rp * dm

        end subroutine mock_get_stale_response_os

    end function test_hess_x_oao_callback

    logical(c_bool) function test_project_oao_callback() bind(C)
        !
        ! this function tests the subroutine which discards the redundant
        ! occupied-occupied and virtual-virtual rotations from a vector, retaining only
        ! its occupied-virtual and virtual-occupied contributions
        !
        use otr_oao, only: project_oao_callback, oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle, n_occ
        use otr_common_unit_tests, only: generate_random_density_matrix

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), vector(n_param), expected(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_project_oao_callback = .true.

        ! generate random density matrices and vector
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_occ(1))
        dm_oao(:, :, 2) = generate_random_density_matrix(n_ao, n_occ(2))
        call random_number(vector)

        ! set up the OAO object
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = n_particle
        oao_object%dm_oao = dm_oao

        ! initialize expected vector
        expected = ref_pack_asymm(ref_project_asymm( &
            ref_unpack_asymm(vector, n_particle, n_ao), dm_oao), n_param)

        ! call routine and determine if values of resulting vector match
        call project_oao_callback(vector, error)
        if (error /= 0) then
            write(stderr, *) "test_project_oao_callback failed: Produced error."
            test_project_oao_callback = .false.
        end if
        if (norm2(vector - expected) > tol) then
            write(stderr, *) "test_project_oao_callback failed: Incorrect vector."
            test_project_oao_callback = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_project_oao_callback

    logical(c_bool) function test_precond_oao_callback() bind(C)
        !
        ! this function tests the subroutine which defines a level-shifted
        ! preconditioner of the OAO basis, which applies the one of the orbital basis
        ! to the OAO object
        !
        use otr_oao, only: precond_oao_callback, oao_object

        real(rp), parameter :: mu = 0.3_rp

        real(rp), allocatable :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! assume tests pass
        test_precond_oao_callback = .true.

        ! set up the OAO object with an eigendecomposition whose eigenvalue pairs, the
        ! open-shell sums 2 (e_i + e_j), shifted by the level shift are all 2, so that
        ! the preconditioner of the orbital basis halves the residual
        call setup_identity_oao_eigenbasis((2.0_rp + mu) / 4.0_rp)
        allocate(residual(oao_object%n_param), precond_residual(oao_object%n_param))
        call random_number(residual)

        ! call routine and determine if the residual is halved
        call precond_oao_callback(residual, mu, precond_residual, error)
        if (error /= 0) then
            write(stderr, *) "test_precond_oao_callback failed: Produced error."
            test_precond_oao_callback = .false.
        end if
        if (norm2(precond_residual - 0.5_rp * residual) > tol) then
            write(stderr, *) "test_precond_oao_callback failed: Incorrect "// &
                "preconditioned residual."
            test_precond_oao_callback = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_precond_oao_callback

    logical(c_bool) function test_precond_pd_oao_callback() bind(C)
        !
        ! this function tests the subroutine which defines the positive-definite
        ! preconditioner of the OAO basis, which applies the one of the orbital basis
        ! to the OAO object
        !
        use otr_oao, only: precond_pd_oao_callback, oao_object

        real(rp), allocatable :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! assume tests pass
        test_precond_pd_oao_callback = .true.

        ! set up the OAO object with an eigendecomposition whose eigenvalue pairs, the
        ! open-shell sums 2 (e_i + e_j), are all -2, so that the positive-definite
        ! preconditioner of the orbital basis, which divides by their magnitudes,
        ! halves the residual, unlike the level-shifted one
        call setup_identity_oao_eigenbasis(-0.5_rp)
        allocate(residual(oao_object%n_param), precond_residual(oao_object%n_param))
        call random_number(residual)

        ! call routine and determine if the residual is halved
        call precond_pd_oao_callback(residual, precond_residual, error)
        if (error /= 0) then
            write(stderr, *) "test_precond_pd_oao_callback failed: Produced error."
            test_precond_pd_oao_callback = .false.
        end if
        if (norm2(precond_residual - 0.5_rp * residual) > tol) then
            write(stderr, *) "test_precond_pd_oao_callback failed: Incorrect "// &
                "preconditioned residual."
            test_precond_pd_oao_callback = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_precond_pd_oao_callback

    logical(c_bool) function test_get_extra_trial_vectors_oao_callback() bind(C)
        !
        ! this function tests the subroutine which returns the extra trial vectors of
        ! the OAO basis for the solver's initial trial space
        !
        use otr_common_unit_tests, only: identity_matrix
        use otr_oao, only: get_extra_trial_vectors_oao_callback, oao_object
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_ao = 3, n_param = n_ao * (n_ao - 1) / 2, n_extra = 2

        real(rp) :: trial_vectors(n_param, n_extra), unit_vector(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_get_extra_trial_vectors_oao_callback = .true.

        ! set up the OAO object with an up-to-date eigendecomposition in the
        ! orthogonalized AO basis itself whose first two eigenvectors span the occupied
        ! space, so that the only negative non-redundant eigenvalue sum belongs to the
        ! pair at packed index 2 and its rotation is the unit vector
        allocate(oao_object)
        call setup_settings(oao_object%settings)
        oao_object%n_ao = n_ao
        oao_object%n_particle = 1
        oao_object%n_param = n_param
        oao_object%hess_eigvecs = reshape(identity_matrix(n_ao), [n_ao, n_ao, 1_ip])
        oao_object%hess_eigvals = reshape([-2.0_rp, -1.0_rp, 1.5_rp], [n_ao, 1_ip])
        oao_object%dm_oao = oao_object%hess_eigvecs
        oao_object%dm_oao(3, 3, 1) = 0.0_rp
        oao_object%hess_eigen_stale = .false.

        ! call routine and determine if the extra trial vectors of the OAO basis are
        ! returned
        call get_extra_trial_vectors_oao_callback(trial_vectors, error)
        if (error /= 0) then
            write(stderr, *) "test_get_extra_trial_vectors_oao_callback failed: "// &
                "Produced error."
            test_get_extra_trial_vectors_oao_callback = .false.
        end if
        unit_vector = 0.0_rp
        unit_vector(2) = 1.0_rp
        if (norm2(trial_vectors(:, 1) - unit_vector) > tol) then
            write(stderr, *) "test_get_extra_trial_vectors_oao_callback failed: "// &
                "Incorrect extra trial vector."
            test_get_extra_trial_vectors_oao_callback = .false.
        end if
        if (norm2(trial_vectors(:, 2)) > tol) then
            write(stderr, *) "test_get_extra_trial_vectors_oao_callback failed: "// &
                "Slot without a negative eigenvalue sum does not vanish."
            test_get_extra_trial_vectors_oao_callback = .false.
        end if

        ! deallocate OAO object
        deallocate(oao_object)

    end function test_get_extra_trial_vectors_oao_callback

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
            write(stderr, *) "test_init_oao_settings failed: Function raised error."
            test_init_oao_settings = .false.
        end if

        ! check settings
        if (.not. (settings == default_settings)) then
            write(stderr, *) "test_init_oao_settings failed: Settings not "// &
                "initialized correctly."
            test_init_oao_settings = .false.
        end if

    end function test_init_oao_settings

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
            write(stderr, *) "test_oao_deconstructor failed: Object not deallocated."
            test_oao_deconstructor = .false.
        end if

        ! call routine again and determine if an already deallocated object is handled
        call oao_deconstructor()
        if (allocated(oao_object)) then
            write(stderr, *) "test_oao_deconstructor failed: Already deallocated "// &
                "object not handled."
            test_oao_deconstructor = .false.
        end if

    end function test_oao_deconstructor

    logical(c_bool) function test_rotate_orbitals_oao() bind(C)
        !
        ! this function tests the subroutine which moves the current orbitals by an
        ! orbital rotation, updating the density matrix in the AO and in the OAO basis
        !
        use otr_common_unit_tests, only: &
            identity_matrix, generate_random_density_matrix, generate_random_symm_matrix
        use otr_oao, only: oao_type, oao_settings_type
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_ao = 2, n_particle = 1, n_occ = 1
        real(rp), parameter :: angle = 0.3_rp

        type(oao_type) :: oao
        type(oao_settings_type) :: settings
        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: dm_oao(n_ao, n_ao, n_particle), s_inv_sqrt(n_ao, n_ao), &
                    rotation(n_ao, n_ao), expected_dm_oao(n_ao, n_ao, n_particle), &
                    expected_dm_ao(n_ao, n_ao, n_particle)
        integer(ip) :: error

        ! assume tests pass
        test_rotate_orbitals_oao = .true.

        ! setup settings object
        call setup_settings(settings)

        ! set up the OAO object with a random transformation to the AO basis
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_occ)
        s_inv_sqrt = identity_matrix(n_ao) + 0.1_rp * generate_random_symm_matrix(n_ao)
        dm_ao(:, :, 1) = matmul(s_inv_sqrt, matmul(dm_oao(:, :, 1), s_inv_sqrt))
        oao%dm_ao => dm_ao
        oao%dm_oao = dm_oao
        oao%s_inv_sqrt = s_inv_sqrt
        oao%response_stale = .false.

        ! initialize expected density matrices, where the rotation matrix is the
        ! exponential of the antisymmetric matrix the rotation is unpacked into
        rotation = reshape([cos(angle), -sin(angle), sin(angle), cos(angle)], &
                           [n_ao, n_ao])
        expected_dm_oao(:, :, 1) = matmul(transpose(rotation), &
                                          matmul(dm_oao(:, :, 1), rotation))
        expected_dm_ao(:, :, 1) = matmul(s_inv_sqrt, &
                                         matmul(expected_dm_oao(:, :, 1), s_inv_sqrt))

        ! call routine and determine if the density matrices are moved and the response
        ! is marked stale
        call oao%rotate_orbitals([angle], settings, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_orbitals_oao failed: Produced error."
            test_rotate_orbitals_oao = .false.
        end if
        if (norm2(oao%dm_oao - expected_dm_oao) > tol) then
            write(stderr, *) "test_rotate_orbitals_oao failed: Incorrect density "// &
                "matrix in the OAO basis."
            test_rotate_orbitals_oao = .false.
        end if
        if (norm2(dm_ao - expected_dm_ao) > tol) then
            write(stderr, *) "test_rotate_orbitals_oao failed: Incorrect density "// &
                "matrix in the AO basis."
            test_rotate_orbitals_oao = .false.
        end if
        if (.not. oao%response_stale) then
            write(stderr, *) "test_rotate_orbitals_oao failed: Response not marked "// &
                "stale."
            test_rotate_orbitals_oao = .false.
        end if

    end function test_rotate_orbitals_oao

    logical(c_bool) function test_calculate_grad_h_diag_oao() bind(C)
        !
        ! this function tests the subroutine which calculates the gradient, Hessian
        ! diagonal and the occupied-occupied and virtual-virtual parts of the Fock
        ! matrix in the OAO basis
        !
        use otr_common_unit_tests, only: identity_matrix, &
                                         generate_random_density_matrix, &
                                         generate_random_symm_matrix, shell_names
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_occ, n_particle_ref => n_particle

        type(oao_type) :: oao
        integer(ip) :: n_particle, n_param
        real(rp) :: dm_oao(n_ao, n_ao, n_particle_ref), &
                    fock_oao(n_ao, n_ao, n_particle_ref), proj_v(n_ao, n_ao), &
                    expected_fock_oo(n_ao, n_ao, n_particle_ref), &
                    expected_fock_vv(n_ao, n_ao, n_particle_ref), &
                    fock_ov(n_ao, n_ao, n_particle_ref), &
                    grad_full(n_ao, n_ao, n_particle_ref), &
                    static_diag(n_ao, n_particle_ref), shell_scale
        real(rp), allocatable :: expected_grad(:), expected_h_diag(:)
        integer(ip) :: i, j
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_calculate_grad_h_diag_oao = .true.

        ! initialize density and Fock matrices for both particle slots, since the
        ! closed- and open-shell cases below both read from these shared arrays
        do i = 1, size(dm_oao, 3)
            dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
            fock_oao(:, :, i) = generate_random_symm_matrix(n_ao)
        end do

        ! initialize expected occupancy-resolved parts of the Fock matrix and the
        ! diagonal of the static Hessian part
        do i = 1, size(dm_oao, 3)
            proj_v = identity_matrix(n_ao) - dm_oao(:, :, i)
            expected_fock_oo(:, :, i) = &
                matmul(dm_oao(:, :, i), matmul(fock_oao(:, :, i), dm_oao(:, :, i)))
            expected_fock_vv(:, :, i) = matmul(proj_v, &
                                               matmul(fock_oao(:, :, i), proj_v))
            fock_ov(:, :, i) = matmul(dm_oao(:, :, i), &
                                      matmul(fock_oao(:, :, i), proj_v))
            static_diag(:, i) = &
                [(expected_fock_vv(j, j, i) - expected_fock_oo(j, j, i), j=1, n_ao)]
        end do

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            shell_scale = merge(4.0_rp, 2.0_rp, n_particle == 1)
            n_param = n_particle * n_ao * (n_ao - 1) / 2

            ! initialize expected gradient and Hessian diagonal, whose elements are the
            ! scaled pairwise sums of the diagonal of the static Hessian part
            do i = 1, n_particle
                grad_full(:, :, i) = shell_scale * &
                                     (fock_ov(:, :, i) - transpose(fock_ov(:, :, i)))
            end do
            expected_grad = ref_pack_asymm(grad_full(:, :, :n_particle), n_param)
            expected_h_diag = ref_hess_eigval_pairs_oao(static_diag(:, :n_particle), &
                                                        n_particle, n_ao, n_param)

            ! set up the OAO object
            oao%n_ao = n_ao
            oao%n_particle = n_particle
            oao%n_param = n_param
            oao%dm_oao = dm_oao(:, :, :n_particle)
            allocate(oao%grad(n_param), oao%h_diag(n_param), &
                     oao%fock_oo(n_ao, n_ao, n_particle), &
                     oao%fock_vv(n_ao, n_ao, n_particle))
            oao%hess_eigen_stale = .false.

            ! call routine and determine if values of resulting quantities match and
            ! the eigendecomposition of the static Hessian part is marked stale
            call oao%calculate_grad_h_diag(fock_oao(:, :, :n_particle))
            if (norm2(oao%fock_oo - expected_fock_oo(:, :, :n_particle)) > tol) then
                write(stderr, *) "test_calculate_grad_h_diag_oao failed: Incorrect "// &
                    "occupied-occupied part of the Fock matrix for the "//case_name// &
                    " case."
                test_calculate_grad_h_diag_oao = .false.
            end if
            if (norm2(oao%fock_vv - expected_fock_vv(:, :, :n_particle)) > tol) then
                write(stderr, *) "test_calculate_grad_h_diag_oao failed: Incorrect "// &
                    "virtual-virtual part of the Fock matrix for the "//case_name// &
                    " case."
                test_calculate_grad_h_diag_oao = .false.
            end if
            if (norm2(oao%grad - expected_grad) > tol) then
                write(stderr, *) "test_calculate_grad_h_diag_oao failed: Incorrect "// &
                    "gradient for the "//case_name//" case."
                test_calculate_grad_h_diag_oao = .false.
            end if
            if (norm2(oao%h_diag - expected_h_diag) > tol) then
                write(stderr, *) "test_calculate_grad_h_diag_oao failed: Incorrect "// &
                    "Hessian diagonal for the "//case_name//" case."
                test_calculate_grad_h_diag_oao = .false.
            end if
            if (.not. oao%hess_eigen_stale) then
                write(stderr, *) "test_calculate_grad_h_diag_oao failed: "// &
                    "Eigendecomposition not marked stale for the "//case_name//" case."
                test_calculate_grad_h_diag_oao = .false.
            end if
            deallocate(oao%grad, oao%h_diag, oao%fock_oo, oao%fock_vv)
        end do

    end function test_calculate_grad_h_diag_oao

    logical(c_bool) function test_refresh_hess_eigen_oao() bind(C)
        !
        ! this function tests the subroutine that diagonalizes the static part of the
        ! Hessian and caches the result if the static part has changed
        !
        use otr_oao, only: oao_type, oao_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_common_unit_tests, only: generate_random_symm_matrix

        type(oao_type) :: oao
        type(oao_settings_type) :: settings
        real(rp) :: fock_oo(n_ao, n_ao, n_particle), fock_vv(n_ao, n_ao, n_particle), &
                    expected_eigvecs(n_ao, n_ao, n_particle), &
                    expected_eigvals(n_ao, n_particle)
        integer(ip) :: j, error

        ! assume tests pass
        test_refresh_hess_eigen_oao = .true.

        ! setup settings object
        call setup_settings(settings)

        ! generate random Fock matrix contributions
        do j = 1, n_particle
            fock_oo(:, :, j) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, j) = generate_random_symm_matrix(n_ao)
        end do

        ! set up the OAO object with a stale eigendecomposition
        oao%n_ao = n_ao
        oao%n_particle = n_particle
        oao%fock_oo = fock_oo
        oao%fock_vv = fock_vv
        oao%hess_eigen_stale = .true.

        ! independently diagonalize the static part
        call ref_diagonalize_static_part(fock_oo, fock_vv, expected_eigvecs, &
                                         expected_eigvals)

        ! call routine and determine if the cached eigendecomposition matches
        call oao%refresh_hess_eigen(settings, error)
        if (error /= 0) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: Produced error."
            test_refresh_hess_eigen_oao = .false.
        end if
        if (norm2(oao%hess_eigvecs - expected_eigvecs) > tol) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: Incorrect "// &
                "eigenvectors."
            test_refresh_hess_eigen_oao = .false.
        end if
        if (norm2(oao%hess_eigvals - expected_eigvals) > tol) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: Incorrect "// &
                "eigenvalues."
            test_refresh_hess_eigen_oao = .false.
        end if
        if (oao%hess_eigen_stale) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: "// &
                "Eigendecomposition still marked stale after being refreshed."
            test_refresh_hess_eigen_oao = .false.
        end if

        ! change the static part without marking the eigendecomposition stale and
        ! determine if the cached eigendecomposition is kept
        oao%fock_vv = 2.0_rp * fock_vv
        call oao%refresh_hess_eigen(settings, error)
        if (error /= 0) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: Produced error "// &
                "for an up-to-date eigendecomposition."
            test_refresh_hess_eigen_oao = .false.
        end if
        if (norm2(oao%hess_eigvals - expected_eigvals) > tol) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: Up-to-date "// &
                "eigendecomposition recomputed."
            test_refresh_hess_eigen_oao = .false.
        end if

    end function test_refresh_hess_eigen_oao

    logical(c_bool) function test_rotate_to_hess_eigenbasis_oao() bind(C)
        !
        ! this function tests the function which rotates a packed antisymmetric vector
        ! into the eigenbasis of the static part of the Hessian for the OAO basis
        !
        test_rotate_to_hess_eigenbasis_oao = check_rotate_hess_eigenbasis_oao(.true.)

    end function test_rotate_to_hess_eigenbasis_oao

    logical(c_bool) function test_rotate_from_hess_eigenbasis_oao() bind(C)
        !
        ! this function tests the function which rotates a packed antisymmetric vector
        ! out of the eigenbasis of the static part of the Hessian for the OAO basis
        !
        test_rotate_from_hess_eigenbasis_oao = check_rotate_hess_eigenbasis_oao(.false.)

    end function test_rotate_from_hess_eigenbasis_oao

    logical(c_bool) function test_get_hess_eigval_pairs_oao() bind(C)
        !
        ! this function tests the function that returns the pairwise sums of the cached
        ! static Hessian part eigenvalues; the closed-shell and open-shell cases apply
        ! different scaling factors
        !
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_particle_ref => n_particle
        use otr_common_unit_tests, only: shell_names

        type(oao_type) :: oao
        real(rp), allocatable :: eigvals(:, :), expected(:), eigval_pairs(:)
        integer(ip) :: n_particle, n_param

        ! assume tests pass
        test_get_hess_eigval_pairs_oao = .true.

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            n_param = n_particle * n_ao * (n_ao - 1) / 2

            ! set up the OAO object with random eigenvalues
            oao%n_ao = n_ao
            oao%n_particle = n_particle
            oao%n_param = n_param
            allocate(eigvals(n_ao, n_particle))
            call random_number(eigvals)
            oao%hess_eigvals = eigvals

            ! independently construct the expected pairwise sums
            expected = ref_hess_eigval_pairs_oao(eigvals, n_particle, n_ao, n_param)

            ! call routine and determine if values match
            eigval_pairs = oao%get_hess_eigval_pairs()
            if (norm2(eigval_pairs - expected) > tol) then
                write(stderr, *) "test_get_hess_eigval_pairs_oao failed: Incorrect "// &
                    "pairwise eigenvalue sums for the "// &
                    trim(shell_names(n_particle))//" case."
                test_get_hess_eigval_pairs_oao = .false.
            end if
            deallocate(eigvals)
        end do

    end function test_get_hess_eigval_pairs_oao

    logical(c_bool) function test_get_extra_trial_vectors_oao() bind(C)
        !
        ! this function tests the subroutine which returns curvature-informed extra
        ! trial vectors, for the closed-shell and the open-shell case
        !
        use otr_oao, only: oao_type, oao_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix, &
                                         generate_random_symm_matrix, shell_names
        use otr_common_test_reference, only: n_particle_ref => n_particle

        integer(ip), parameter :: n_ao = 3
        type(oao_type) :: oao
        type(oao_settings_type) :: settings
        integer(ip) :: n_particle, n_param, n_extra, i, k, error
        integer(ip), allocatable :: ref_idx(:), n_occ(:)
        real(rp), allocatable :: eigvals(:, :), eigvecs(:, :, :), dm_oao(:, :, :), &
                                 trial_vectors(:, :), unit_vector(:), expected(:)
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_get_extra_trial_vectors_oao = .true.

        ! setup settings object
        call setup_settings(settings)

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            n_param = n_particle * n_ao * (n_ao - 1) / 2
            allocate(eigvals(n_ao, n_particle))
            if (n_particle == 1) then
                ! the first two eigenvectors span the occupied space, so the pairwise
                ! sums are -3.0 for the redundant occupied-occupied pair (1, 2) at
                ! packed index 1, -0.5 for (1, 3) at index 2 and +0.5 for (2, 3) at
                ! index 3; only one non-redundant sum is negative while the other is
                ! positive, so the trailing slot can only be left vanishing by
                ! rejecting a positive sum rather than by running out of candidates
                eigvals(:, 1) = [-2.0_rp, -1.0_rp, 1.5_rp]
                n_occ = [2_ip]
                ref_idx = [2_ip]
            else
                ! the two channels are given a different number of occupied
                ! eigenvectors, so that the redundant pair of the first channel is
                ! virtual-virtual while that of the second is occupied-occupied; the
                ! four remaining negative sums are all distinct, so their expected
                ! order is unambiguous and spans both channels, and the trailing slot
                ! stays empty because only the two redundant pairs are left
                eigvals(:, 1) = [0.5_rp, -1.0_rp, -2.0_rp]
                eigvals(:, 2) = [-1.1_rp, -2.3_rp, 0.7_rp]
                n_occ = [1_ip, 2_ip]
                ref_idx = [6_ip, 2_ip, 1_ip, 5_ip]
            end if
            n_extra = size(ref_idx, kind=ip) + 1

            ! random eigenbasis per channel, with the leading eigenvectors spanning the
            ! occupied space so that the density matrix is their projector and the
            ! occupied-virtual split of the eigenvectors is known
            allocate(eigvecs(n_ao, n_ao, n_particle), dm_oao(n_ao, n_ao, n_particle))
            do k = 1, n_particle
                eigvecs(:, :, k) = generate_random_orthogonal_matrix(n_ao)
                dm_oao(:, :, k) = matmul(eigvecs(:, :n_occ(k), k), &
                                         transpose(eigvecs(:, :n_occ(k), k)))
            end do

            ! set up the OAO object with the eigendecomposition injected directly
            oao%n_ao = n_ao
            oao%n_particle = n_particle
            oao%n_param = n_param
            oao%hess_eigvecs = eigvecs
            oao%hess_eigvals = eigvals
            oao%dm_oao = dm_oao
            oao%hess_eigen_stale = .false.

            ! call routine and determine if the expected vectors, the rotations of the
            ! unit vectors at the expected packed indices ordered by increasing
            ! eigenvalue sum, are returned and the trailing slot vanishes
            allocate(trial_vectors(n_param, n_extra), unit_vector(n_param))
            call oao%get_extra_trial_vectors(trial_vectors, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_get_extra_trial_vectors_oao failed: "// &
                    "Produced error for the "//case_name//" case."
                test_get_extra_trial_vectors_oao = .false.
            end if
            do i = 1, size(ref_idx, kind=ip)
                unit_vector = 0.0_rp
                unit_vector(ref_idx(i)) = 1.0_rp
                expected = ref_rotate_eigenbasis_oao(unit_vector, eigvecs, n_particle, &
                                                     n_ao, .false.)
                if (norm2(trial_vectors(:, i) - expected) > tol) then
                    write(stderr, *) "test_get_extra_trial_vectors_oao failed: "// &
                        "Incorrect extra trial vector for the "//case_name//" case."
                    test_get_extra_trial_vectors_oao = .false.
                end if
            end do
            if (norm2(trial_vectors(:, n_extra)) > tol) then
                write(stderr, *) "test_get_extra_trial_vectors_oao failed: Slot "// &
                    "without a negative eigenvalue sum does not vanish for the "// &
                    case_name//" case."
                test_get_extra_trial_vectors_oao = .false.
            end if

            ! mark the eigendecomposition stale for random Fock matrix blocks and
            ! determine if it is refreshed
            oao%fock_oo = dm_oao
            oao%fock_vv = dm_oao
            do k = 1, n_particle
                oao%fock_oo(:, :, k) = generate_random_symm_matrix(n_ao)
                oao%fock_vv(:, :, k) = generate_random_symm_matrix(n_ao)
            end do
            oao%hess_eigen_stale = .true.
            call oao%get_extra_trial_vectors(trial_vectors, settings, error)
            if (error /= 0) then
                write(stderr, *) "test_get_extra_trial_vectors_oao failed: "// &
                    "Produced error for a stale eigendecomposition for the "// &
                    case_name//" case."
                test_get_extra_trial_vectors_oao = .false.
            end if
            if (oao%hess_eigen_stale) then
                write(stderr, *) "test_get_extra_trial_vectors_oao failed: "// &
                    "Eigendecomposition not marked refreshed for the "//case_name// &
                    " case."
                test_get_extra_trial_vectors_oao = .false.
            end if
            deallocate(eigvals, eigvecs, dm_oao, trial_vectors, unit_vector, expected)
        end do

    end function test_get_extra_trial_vectors_oao

    logical(c_bool) function test_rotate_dm_ao() bind(C)
        !
        ! this function tests the subroutine which returns the rotated density matrix
        ! in the AO basis
        !
        use otr_common_unit_tests, only: &
            identity_matrix, generate_random_density_matrix, generate_random_symm_matrix
        use otr_oao, only: rotate_dm_ao, oao_settings_type
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_ao = 2, n_particle = 1, n_occ = 1
        real(rp), parameter :: angle = 0.3_rp

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), s_inv_sqrt(n_ao, n_ao), &
                    rotation(n_ao, n_ao), expected_dm_oao(n_ao, n_ao, n_particle), &
                    expected_dm_ao(n_ao, n_ao, n_particle), &
                    rot_dm_ao(n_ao, n_ao, n_particle), &
                    rot_dm_oao(n_ao, n_ao, n_particle)
        integer(ip) :: error
        type(oao_settings_type) :: settings

        ! assume tests pass
        test_rotate_dm_ao = .true.

        ! setup settings object
        call setup_settings(settings)

        ! set up a density matrix and a random transformation to the AO basis
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_occ)
        s_inv_sqrt = identity_matrix(n_ao) + 0.1_rp * generate_random_symm_matrix(n_ao)

        ! initialize expected density matrices, where the rotation matrix is the
        ! exponential of the antisymmetric matrix the rotation is unpacked into and the
        ! rotated density matrix stays idempotent so that purification leaves it
        ! unchanged
        rotation = reshape([cos(angle), -sin(angle), sin(angle), cos(angle)], &
                           [n_ao, n_ao])
        expected_dm_oao(:, :, 1) = matmul(transpose(rotation), &
                                          matmul(dm_oao(:, :, 1), rotation))
        expected_dm_ao(:, :, 1) = matmul(s_inv_sqrt, &
                                         matmul(expected_dm_oao(:, :, 1), s_inv_sqrt))

        ! call routine and determine if values of the resulting density matrices match
        call rotate_dm_ao([angle], dm_oao, s_inv_sqrt, rot_dm_ao, settings, error, &
                          rot_dm_oao)
        if (error /= 0) then
            write(stderr, *) "test_rotate_dm_ao failed: Produced error."
            test_rotate_dm_ao = .false.
        end if
        if (norm2(rot_dm_oao - expected_dm_oao) > tol) then
            write(stderr, *) "test_rotate_dm_ao failed: Incorrect rotated density "// &
                "matrix in the OAO basis."
            test_rotate_dm_ao = .false.
        end if
        if (norm2(rot_dm_ao - expected_dm_ao) > tol) then
            write(stderr, *) "test_rotate_dm_ao failed: Incorrect rotated density "// &
                "matrix in the AO basis."
            test_rotate_dm_ao = .false.
        end if

        ! call routine without the optional rotated density matrix in the OAO basis and
        ! determine if values of the resulting density matrix match
        rot_dm_ao = 0.0_rp
        call rotate_dm_ao([angle], dm_oao, s_inv_sqrt, rot_dm_ao, settings, error)
        if (error /= 0) then
            write(stderr, *) "test_rotate_dm_ao failed: Produced error without the "// &
                "optional argument."
            test_rotate_dm_ao = .false.
        end if
        if (norm2(rot_dm_ao - expected_dm_ao) > tol) then
            write(stderr, *) "test_rotate_dm_ao failed: Incorrect rotated density "// &
                "matrix in the AO basis without the optional argument."
            test_rotate_dm_ao = .false.
        end if

    end function test_rotate_dm_ao

    logical(c_bool) function test_hess_x_static_oao() bind(C)
        !
        ! this function tests the function which applies the static part of the Hessian
        ! in the OAO basis to unpacked trial vectors
        !
        use otr_common_unit_tests, only: generate_random_symm_matrix, shell_names
        use otr_oao, only: hess_x_static_oao
        use otr_common_test_reference, only: n_ao, n_particle_ref => n_particle

        real(rp) :: fock_oo(n_ao, n_ao, n_particle_ref), &
                    fock_vv(n_ao, n_ao, n_particle_ref)
        real(rp), allocatable :: x(:), x_full(:, :, :)
        integer(ip) :: n_particle, i

        ! assume tests pass
        test_hess_x_static_oao = .true.

        ! generate random Fock matrix contributions, which differ between the particle
        ! channels
        do i = 1, n_particle_ref
            fock_oo(:, :, i) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, i) = generate_random_symm_matrix(n_ao)
        end do

        ! loop over the closed-shell and the open-shell case and determine if the
        ! static part of random unpacked trial vectors matches
        do n_particle = 1, n_particle_ref
            allocate(x(n_particle * n_ao * (n_ao - 1) / 2))
            call random_number(x)
            x_full = ref_unpack_asymm(x, n_particle, n_ao)
            if (norm2(hess_x_static_oao(x_full, fock_oo(:, :, :n_particle), &
                                        fock_vv(:, :, :n_particle)) - &
                      ref_hess_x_static_oao(x_full, fock_oo(:, :, :n_particle), &
                                            fock_vv(:, :, :n_particle))) > tol) then
                write(stderr, *) "test_hess_x_static_oao failed: Incorrect static "// &
                    "part for the "//trim(shell_names(n_particle))//" case."
                test_hess_x_static_oao = .false.
            end if
            deallocate(x)
        end do

    end function test_hess_x_static_oao

    logical(c_bool) function test_project_asymm() bind(C)
        !
        ! this function tests the function which retains only the occupied-virtual and
        ! virtual-occupied contributions to a matrix in antisymmetric form
        !
        use otr_oao, only: project_asymm
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 3
        real(rp) :: matrix(n_ao, n_ao, n_particle), dm_oao(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle)
        real(rp), allocatable :: projected_matrix(:, :, :)

        ! assume tests pass
        test_project_asymm = .true.

        ! initialize matrices and density matrices, each occupying a single, distinct
        ! orbital only
        matrix(:, :, 1) = reshape( &
            [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp, 7.0_rp, 8.0_rp, 9.0_rp], &
            [n_ao, n_ao])
        matrix(:, :, 2) = reshape( &
            [9.0_rp, 8.0_rp, 7.0_rp, 6.0_rp, 5.0_rp, 4.0_rp, 3.0_rp, 2.0_rp, 1.0_rp], &
            [n_ao, n_ao])
        dm_oao = 0.0_rp
        dm_oao(1, 1, 1) = 1.0_rp
        dm_oao(2, 2, 2) = 1.0_rp

        ! initialize expected matrices, where only the antisymmetrized occupied-virtual
        ! elements survive
        expected(:, :, 1) = reshape([0.0_rp, -4.0_rp, -7.0_rp, 4.0_rp, 0.0_rp, 0.0_rp, &
                                     7.0_rp, 0.0_rp, 0.0_rp], [n_ao, n_ao])
        expected(:, :, 2) = reshape([0.0_rp, 8.0_rp, 0.0_rp, -8.0_rp, 0.0_rp, -2.0_rp, &
                                     0.0_rp, 2.0_rp, 0.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrix match
        projected_matrix = project_asymm(matrix, dm_oao)
        if (size(projected_matrix, 1) /= n_ao .or. &
            size(projected_matrix, 3) /= n_particle) then
            write(stderr, *) "test_project_asymm failed: Incorrect matrix dimensions."
            test_project_asymm = .false.
            return
        end if
        if (norm2(projected_matrix - expected) > tol) then
            write(stderr, *) "test_project_asymm failed: Incorrect matrix values."
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
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 3
        real(rp) :: x_full(n_ao, n_ao, n_particle), dm_oao(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle)
        real(rp), allocatable :: projected_matrix(:, :, :)

        ! assume tests pass
        test_project_symm = .true.

        ! initialize antisymmetric trial vectors and density matrices, each occupying a
        ! single, distinct orbital only
        x_full(:, :, 1) = reshape([0.0_rp, -1.0_rp, -2.0_rp, 1.0_rp, 0.0_rp, -3.0_rp, &
                                   2.0_rp, 3.0_rp, 0.0_rp], [n_ao, n_ao])
        x_full(:, :, 2) = reshape([0.0_rp, -4.0_rp, -5.0_rp, 4.0_rp, 0.0_rp, -6.0_rp, &
                                   5.0_rp, 6.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_oao = 0.0_rp
        dm_oao(1, 1, 1) = 1.0_rp
        dm_oao(2, 2, 2) = 1.0_rp

        ! initialize expected matrices, where only the symmetrized occupied-virtual
        ! elements survive
        expected(:, :, 1) = reshape( &
            [0.0_rp, 1.0_rp, 2.0_rp, 1.0_rp, 0.0_rp, 0.0_rp, 2.0_rp, 0.0_rp, 0.0_rp], &
            [n_ao, n_ao])
        expected(:, :, 2) = reshape([0.0_rp, -4.0_rp, 0.0_rp, -4.0_rp, 0.0_rp, 6.0_rp, &
                                     0.0_rp, 6.0_rp, 0.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrix match
        projected_matrix = project_symm(x_full, dm_oao)
        if (size(projected_matrix, 1) /= n_ao .or. &
            size(projected_matrix, 3) /= n_particle) then
            write(stderr, *) "test_project_symm failed: Incorrect matrix dimensions."
            test_project_symm = .false.
            return
        end if
        if (norm2(projected_matrix - expected) > tol) then
            write(stderr, *) "test_project_symm failed: Incorrect matrix values."
            test_project_symm = .false.
        end if
        deallocate(projected_matrix)

    end function test_project_symm

    logical(c_bool) function test_purify() bind(C)
        !
        ! this function tests the subroutine which purifies a density matrix
        !
        use otr_oao, only: purify
        use otr_common_unit_tests, only: generate_random_density_matrix
        use otr_common_test_reference, only: n_particle, n_occ

        integer(ip), parameter :: n_ao = 3
        real(rp) :: dm(n_ao, n_ao, n_particle), expected(n_ao, n_ao, n_particle), &
                    idempotent_dm(n_ao, n_ao, n_particle), &
                    purified_dm(n_ao, n_ao, n_particle)
        integer(ip) :: i

        ! assume tests pass
        test_purify = .true.

        ! initialize density matrices with fractional occupations
        dm(:, :, 1) = reshape( &
            [0.6_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.2_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.3_rp], &
            [n_ao, n_ao])
        dm(:, :, 2) = reshape( &
            [0.7_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.4_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.5_rp], &
            [n_ao, n_ao])

        ! initialize expected density matrices, where the occupations are driven
        ! towards zero and one
        expected(:, :, 1) = reshape([0.648_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.104_rp, &
                                     0.0_rp, 0.0_rp, 0.0_rp, 0.216_rp], [n_ao, n_ao])
        expected(:, :, 2) = reshape([0.784_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.352_rp, &
                                     0.0_rp, 0.0_rp, 0.0_rp, 0.5_rp], [n_ao, n_ao])

        ! call routine and determine if values of resulting density matrices match
        call purify(dm)
        if (norm2(dm - expected) > tol) then
            write(stderr, *) "test_purify failed: Incorrect density matrix values."
            test_purify = .false.
        end if

        ! initialize idempotent density matrices
        do i = 1, n_particle
            idempotent_dm(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do
        purified_dm = idempotent_dm

        ! call routine and determine if idempotent density matrices are left unchanged
        call purify(purified_dm)
        if (norm2(purified_dm - idempotent_dm) > tol) then
            write(stderr, *) "test_purify failed: Idempotent density matrices not "// &
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
        trans_matrix = reshape([1.0_rp, 0.0_rp, 2.0_rp, 1.0_rp], [n_ao, n_ao])
        matrix(:, :, 1) = reshape([1.0_rp, 0.0_rp, 0.0_rp, 2.0_rp], [n_ao, n_ao])

        ! initialize expected matrix
        expected(:, :, 1) = reshape([1.0_rp, 0.0_rp, 6.0_rp, 2.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrix match
        matrix_transformed = symmetric_transformation(trans_matrix, matrix)
        if (size(matrix_transformed, 1) /= n_ao .or. &
            size(matrix_transformed, 3) /= n_particle) then
            write(stderr, *) "test_symmetric_transformation failed: Incorrect "// &
                "matrix dimensions."
            test_symmetric_transformation = .false.
            return
        end if
        if (norm2(matrix_transformed - expected) > tol) then
            write(stderr, *) "test_symmetric_transformation failed: Incorrect "// &
                "matrix values."
            test_symmetric_transformation = .false.
        end if
        deallocate(matrix_transformed)

    end function test_symmetric_transformation

    logical(c_bool) function test_unpack_asymm() bind(C)
        !
        ! this function tests the function which unpacks an antisymmetric matrix
        !
        use otr_oao, only: unpack_asymm
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 3
        real(rp) :: expected(n_ao, n_ao, n_particle)
        real(rp), allocatable :: matrix(:, :, :)

        ! assume tests pass
        test_unpack_asymm = .true.

        ! initialize expected matrices
        expected(:, :, 1) = reshape([0.0_rp, -1.0_rp, -2.0_rp, 1.0_rp, 0.0_rp, &
                                     -3.0_rp, 2.0_rp, 3.0_rp, 0.0_rp], [n_ao, n_ao])
        expected(:, :, 2) = reshape([0.0_rp, -4.0_rp, -5.0_rp, 4.0_rp, 0.0_rp, &
                                     -6.0_rp, 5.0_rp, 6.0_rp, 0.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting matrices
        ! match
        matrix = unpack_asymm([1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp], &
                              n_particle, n_ao)
        if (size(matrix, 1) /= n_ao .or. size(matrix, 2) /= n_ao .or. &
            size(matrix, 3) /= n_particle) then
            write(stderr, *) "test_unpack_asymm failed: Incorrect matrix dimensions."
            test_unpack_asymm = .false.
            return
        end if
        if (norm2(matrix - expected) > tol) then
            write(stderr, *) "test_unpack_asymm failed: Incorrect matrix values."
            test_unpack_asymm = .false.
        end if
        deallocate(matrix)

    end function test_unpack_asymm

    logical(c_bool) function test_pack_asymm() bind(C)
        !
        ! this function tests the function which packs an antisymmetric matrix
        !
        use otr_oao, only: pack_asymm
        use otr_common_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 3, n_param = n_particle * n_ao * (n_ao - 1) / 2
        real(rp) :: matrix(n_ao, n_ao, n_particle)
        real(rp), allocatable :: matrix_nonred(:)

        ! assume tests pass
        test_pack_asymm = .true.

        ! initialize antisymmetric matrices
        matrix(:, :, 1) = reshape([0.0_rp, -1.0_rp, -2.0_rp, 1.0_rp, 0.0_rp, -3.0_rp, &
                                   2.0_rp, 3.0_rp, 0.0_rp], [n_ao, n_ao])
        matrix(:, :, 2) = reshape([0.0_rp, -3.0_rp, -4.0_rp, 3.0_rp, 0.0_rp, -5.0_rp, &
                                   4.0_rp, 5.0_rp, 0.0_rp], [n_ao, n_ao])

        ! call routine and determine if dimensions and values of resulting vector
        ! match, where each spin is packed separately
        matrix_nonred = pack_asymm(matrix, n_param)
        if (size(matrix_nonred) /= n_param) then
            write(stderr, *) "test_pack_asymm failed: Incorrect vector dimension."
            test_pack_asymm = .false.
            return
        end if
        if (norm2(matrix_nonred - [1.0_rp, 2.0_rp, 3.0_rp, 3.0_rp, 4.0_rp, 5.0_rp]) > &
            tol) then
            write(stderr, *) "test_pack_asymm failed: Incorrect vector values."
            test_pack_asymm = .false.
        end if
        deallocate(matrix_nonred)

    end function test_pack_asymm

end module otr_oao_unit_tests
