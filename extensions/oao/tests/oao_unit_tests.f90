! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_unit_tests

    use opentrustregion_core, only: rp, ip, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

contains

    subroutine mock_host_project(vector, error, context)
        !
        ! this subroutine is a mock host projection, which records its call, checks the
        ! host context it receives, returns the error code set in it and discards the
        ! first element of the vector, so that it does not commute with the OAO
        ! projection
        !
        use otr_common_mock, only: record_host_call

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        error = record_host_call("project", context)
        vector(1) = 0.0_rp

    end subroutine mock_host_project

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
        use opentrustregion_unit_tests, only: identity_matrix

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

    function generate_random_nonredundant_vector(n_param, n_particle, n_ao, dm_oao) &
        result(vector)
        !
        ! this function generates a random vector confined to the non-redundant
        ! subspace
        !
        integer(ip), intent(in) :: n_param, n_particle, n_ao
        real(rp), intent(in) :: dm_oao(:, :, :)
        real(rp) :: vector(n_param)

        call random_number(vector)
        vector = ref_pack_asymm(ref_project_asymm( &
            ref_unpack_asymm(vector, n_particle, n_ao), dm_oao), n_param)

    end function generate_random_nonredundant_vector

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
        ! this function tests the subroutine which returns the OAO orbital updating
        ! function and objective function for the closed-shell case, which hands the
        ! density matrix of the caller as the only particle channel and the closed-shell
        ! density matrix evaluating function to the common setup
        !
        use opentrustregion_core, only: obj_func_type, update_orbs_type, &
                                        solver_state_type
        use otr_oao, only: oao_factory_cs, oao_type, oao_settings_type, error_oao
        use otr_common, only: evaluate_dm_cs_type, extension_solver_settings_type, &
                              extension_level_type, update_orbs_callback, &
                              obj_func_callback
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix
        use otr_common_test_reference, only: n_ao, n_occ, operator(==)
        use otr_common_unit_tests, only: generate_random_density_matrix, &
                                         mock_evaluate_dm_cs, wired_level

        type(test_context_type), target :: context
        integer(ip), parameter :: n_particle = 1

        real(rp), target :: dm_ao(n_ao, n_ao), dm_ao_fewer(n_ao - 1, n_ao - 1)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: error
        type(oao_settings_type) :: settings
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        type(extension_solver_settings_type), target :: solver_settings
        type(solver_state_type), target :: state
        type(extension_level_type), pointer :: level

        ! assume tests pass
        test_oao_factory_cs = .true.

        ! setup settings objects
        call setup_settings(settings, context)

        ! initialize density matrix and an orthonormal AO basis
        dm_ao = generate_random_density_matrix(n_ao, n_occ(1))
        ao_overlap = identity_matrix(n_ao)

        ! initialize callback function pointers
        evaluate_dm_funptr => mock_evaluate_dm_cs

        ! call routine and determine if the role callback functions are returned and
        ! the solver settings wired with an OAO object, the remaining setup is covered
        ! by the test of the common setup
        call oao_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                            update_orbs_oao_funptr, obj_func_oao_funptr, &
                            solver_settings, state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_oao_factory_cs failed: Produced error."
            test_oao_factory_cs = .false.
        end if
        if (.not. associated(update_orbs_oao_funptr, update_orbs_callback) .or. &
            .not. associated(obj_func_oao_funptr, obj_func_callback)) then
            write(stderr, *) "test_oao_factory_cs failed: Role callback functions "// &
                "not returned."
            test_oao_factory_cs = .false.
        end if
        level => wired_level(solver_settings)
        if (.not. associated(level)) then
            write(stderr, *) "test_oao_factory_cs failed: Solver settings not "// &
                "wired with an object."
            test_oao_factory_cs = .false.
        else
            ! determine if the OAO object points to the density matrix of the caller as
            ! its only particle channel, reads the record of the state and holds the
            ! given settings and only the closed-shell density matrix evaluating
            ! function
            select type (oao => level%extension)
            type is (oao_type)
                if (.not. associated(oao%dm_ao)) then
                    write(stderr, *) "test_oao_factory_cs failed: Density matrix "// &
                        "not associated."
                    test_oao_factory_cs = .false.
                else if (any(shape(oao%dm_ao) /= [n_ao, n_ao, n_particle])) then
                    write(stderr, *) "test_oao_factory_cs failed: Density matrix "// &
                        "not associated as a single particle channel."
                    test_oao_factory_cs = .false.
                else
                    dm_ao(1, 1) = dm_ao(1, 1) + 1.0_rp
                    if (norm2(oao%dm_ao(:, :, 1) - dm_ao) > tol) then
                        write(stderr, *) "test_oao_factory_cs failed: Density "// &
                            "matrix of the caller not taken over."
                        test_oao_factory_cs = .false.
                    end if
                end if
                if (.not. associated(oao%state, state)) then
                    write(stderr, *) "test_oao_factory_cs failed: Record of the "// &
                        "state not passed on."
                    test_oao_factory_cs = .false.
                end if
                if (.not. (oao%settings == settings%orbital_settings_type)) then
                    write(stderr, *) "test_oao_factory_cs failed: Settings not "// &
                        "passed on."
                    test_oao_factory_cs = .false.
                end if
                if (.not. associated(oao%evaluate_dm_cs, mock_evaluate_dm_cs) .or. &
                    associated(oao%evaluate_dm_os)) then
                    write(stderr, *) "test_oao_factory_cs failed: Closed-shell "// &
                        "density matrix evaluating function not handed over."
                    test_oao_factory_cs = .false.
                end if
            class default
                write(stderr, *) "test_oao_factory_cs failed: Solver settings not "// &
                    "wired with an OAO object."
                test_oao_factory_cs = .false.
            end select
        end if

        ! call routine for a density matrix with fewer AOs than given, whose shape the
        ! common setup has to be handed, and determine if its error is passed on
        dm_ao_fewer = 0.0_rp
        call oao_factory_cs(dm_ao_fewer, ao_overlap, n_particle, n_ao, &
                            evaluate_dm_funptr, update_orbs_oao_funptr, &
                            obj_func_oao_funptr, solver_settings, state, error, &
                            settings)
        if (error /= error_oao + 1) then
            write(stderr, *) "test_oao_factory_cs failed: Error of the common "// &
                "setup not passed on for a density matrix with fewer AOs than given."
            test_oao_factory_cs = .false.
        end if

    end function test_oao_factory_cs

    logical(c_bool) function test_oao_factory_os() bind(C)
        !
        ! this function tests the subroutine which returns the OAO orbital updating
        ! function and objective function for the open-shell case, which hands the
        ! density matrix of the caller and the open-shell density matrix evaluating
        ! function to the common setup
        !
        use opentrustregion_core, only: obj_func_type, update_orbs_type, &
                                        solver_state_type
        use otr_oao, only: oao_factory_os, oao_type, oao_settings_type, error_oao
        use otr_common, only: evaluate_dm_os_type, extension_solver_settings_type, &
                              extension_level_type, update_orbs_callback, &
                              obj_func_callback
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix
        use otr_common_test_reference, only: n_ao, n_particle, n_occ, operator(==)
        use otr_common_unit_tests, only: generate_random_density_matrix, &
                                         mock_evaluate_dm_os, wired_level

        type(test_context_type), target :: context
        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: ao_overlap(n_ao, n_ao)
        integer(ip) :: i, error
        type(oao_settings_type) :: settings
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        type(extension_solver_settings_type), target :: solver_settings
        type(solver_state_type), target :: state
        type(extension_level_type), pointer :: level

        ! assume tests pass
        test_oao_factory_os = .true.

        ! setup settings objects
        call setup_settings(settings, context)

        ! initialize density matrices and an orthonormal AO basis
        do i = 1, n_particle
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do
        ao_overlap = identity_matrix(n_ao)

        ! initialize callback function pointers
        evaluate_dm_funptr => mock_evaluate_dm_os

        ! call routine and determine if the role callback functions are returned and
        ! the solver settings wired with an OAO object, the remaining setup is covered
        ! by the test of the common setup
        call oao_factory_os(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_funptr, &
                            update_orbs_oao_funptr, obj_func_oao_funptr, &
                            solver_settings, state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_oao_factory_os failed: Produced error."
            test_oao_factory_os = .false.
        end if
        if (.not. associated(update_orbs_oao_funptr, update_orbs_callback) .or. &
            .not. associated(obj_func_oao_funptr, obj_func_callback)) then
            write(stderr, *) "test_oao_factory_os failed: Role callback functions "// &
                "not returned."
            test_oao_factory_os = .false.
        end if
        level => wired_level(solver_settings)
        if (.not. associated(level)) then
            write(stderr, *) "test_oao_factory_os failed: Solver settings not "// &
                "wired with an object."
            test_oao_factory_os = .false.
        else
            ! determine if the OAO object points to the density matrix of the caller,
            ! reads the record of the state and holds the given settings and only the
            ! open-shell density matrix evaluating function
            select type (oao => level%extension)
            type is (oao_type)
                if (.not. associated(oao%dm_ao, dm_ao)) then
                    write(stderr, *) "test_oao_factory_os failed: Density matrix "// &
                        "of the caller not taken over."
                    test_oao_factory_os = .false.
                end if
                if (.not. associated(oao%state, state)) then
                    write(stderr, *) "test_oao_factory_os failed: Record of the "// &
                        "state not passed on."
                    test_oao_factory_os = .false.
                end if
                if (.not. (oao%settings == settings%orbital_settings_type)) then
                    write(stderr, *) "test_oao_factory_os failed: Settings not "// &
                        "passed on."
                    test_oao_factory_os = .false.
                end if
                if (.not. associated(oao%evaluate_dm_os, mock_evaluate_dm_os) .or. &
                    associated(oao%evaluate_dm_cs)) then
                    write(stderr, *) "test_oao_factory_os failed: Open-shell "// &
                        "density matrix evaluating function not handed over."
                    test_oao_factory_os = .false.
                end if
            class default
                write(stderr, *) "test_oao_factory_os failed: Solver settings not "// &
                    "wired with an OAO object."
                test_oao_factory_os = .false.
            end select
        end if

        ! call routine for an invalid number of particle channels and determine if the
        ! error of the common setup is passed on
        call oao_factory_os(dm_ao, ao_overlap, 3_ip, n_ao, evaluate_dm_funptr, &
                            update_orbs_oao_funptr, obj_func_oao_funptr, &
                            solver_settings, state, error, settings)
        if (error /= error_oao + 1) then
            write(stderr, *) "test_oao_factory_os failed: Error of the common "// &
                "setup not passed on."
            test_oao_factory_os = .false.
        end if

    end function test_oao_factory_os

    logical(c_bool) function test_oao_factory_common() bind(C)
        !
        ! this function tests the subroutine which sets up a new OAO object for the
        ! density matrix of the caller with the density matrix evaluating function of
        ! either shell, if one is given, after the sanity check of its settings, wires
        ! the solver settings with it and the record of the state and returns its
        ! orbital updating function and objective function
        !
        use opentrustregion_core, only: obj_func_type, update_orbs_type, &
                                        solver_state_type, verbosity_warning
        use otr_common, only: evaluate_dm_cs_type, evaluate_dm_os_type, &
                              extension_solver_settings_type, extension_level_type, &
                              update_orbs_callback, obj_func_callback, &
                              evaluate_dm_missing_warning_msg
        use otr_oao, only: oao_factory_common, oao_type, oao_settings_type, &
                           default_oao_settings, error_oao
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, &
                                             operator(==), same_logger_and_context
        use opentrustregion_unit_tests, only: setup_settings, test_context_type
        use otr_common_unit_tests, only: generate_random_density_matrix, &
                                         mock_evaluate_dm_cs, mock_evaluate_dm_os, &
                                         wired_level

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle_ref)
        real(rp) :: ao_overlap(n_ao, n_ao), overlap_diag(n_ao), s_sqrt(n_ao, n_ao), &
                    s_inv_sqrt(n_ao, n_ao)
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        type(extension_solver_settings_type), allocatable, target :: solver_settings, &
                                                                     other_settings
        type(solver_state_type), target :: state, other_state
        type(oao_settings_type) :: settings, expected_settings
        type(test_context_type), target :: context
        type(extension_level_type), pointer :: level, other_level
        integer(ip) :: n_particle, n_param, i, i_case, error
        logical :: owned, stored
        character(len=45), parameter :: case_names(4) = &
            [character(len=45) :: "for the closed-shell case", "for the open-shell "// &
             "case", "without settings", "without a density matrix evaluating function"]

        ! assume tests pass
        test_oao_factory_common = .true.

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

        ! initialize random density matrices
        do i = 1, n_particle_ref
            dm_ao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
        end do

        ! create objects for either shell, for no settings and without a density matrix
        ! evaluating function, which prints a warning at the verbosity of warnings,
        ! each of which has to be an OAO object
        ! set up for the density matrix of the caller with the settings and the density
        ! matrix evaluating function and to wire the solver settings, which own it,
        ! with the record of the state
        do i_case = 1, size(case_names)
            allocate(solver_settings)
            n_particle = merge(2_ip, 1_ip, i_case == 2)
            call setup_settings(settings, context)
            if (i_case == 4) settings%verbose = verbosity_warning
            expected_settings = settings
            if (i_case == 3) expected_settings = default_oao_settings
            evaluate_dm_cs_funptr => null()
            evaluate_dm_os_funptr => null()
            if (i_case == 2) then
                evaluate_dm_os_funptr => mock_evaluate_dm_os
            else if (i_case /= 4) then
                evaluate_dm_cs_funptr => mock_evaluate_dm_cs
            end if
            if (i_case == 3) then
                call oao_factory_common(dm_ao(:, :, :n_particle), ao_overlap, &
                                        n_particle, n_ao, evaluate_dm_cs_funptr, &
                                        evaluate_dm_os_funptr, update_orbs_oao_funptr, &
                                        obj_func_oao_funptr, solver_settings, state, &
                                        error)
            else
                call oao_factory_common(dm_ao(:, :, :n_particle), ao_overlap, &
                                        n_particle, n_ao, evaluate_dm_cs_funptr, &
                                        evaluate_dm_os_funptr, update_orbs_oao_funptr, &
                                        obj_func_oao_funptr, solver_settings, state, &
                                        error, settings)
            end if
            if (error /= 0) then
                write(stderr, *) "test_oao_factory_common failed: Produced error "// &
                    trim(case_names(i_case))//"."
                test_oao_factory_common = .false.
                deallocate(solver_settings)
                cycle
            end if
            if (.not. associated(update_orbs_oao_funptr, update_orbs_callback) .or. &
                .not. associated(obj_func_oao_funptr, obj_func_callback)) then
                write(stderr, *) "test_oao_factory_common failed: Role callback "// &
                    "functions not returned "//trim(case_names(i_case))//"."
                test_oao_factory_common = .false.
            end if
            level => wired_level(solver_settings)
            if (.not. associated(level)) then
                write(stderr, *) "test_oao_factory_common failed: Solver settings "// &
                    "not wired with an object "//trim(case_names(i_case))//"."
                test_oao_factory_common = .false.
                deallocate(solver_settings)
                cycle
            end if
            owned = .false.
            if (allocated(solver_settings%extensions)) then
                if (size(solver_settings%extensions) == 1) owned = &
                    associated(level%extension, solver_settings%extensions(1)%object)
            end if
            if (.not. owned .or. .not. associated(level%extension%state, state)) then
                write(stderr, *) "test_oao_factory_common failed: Object not owned "// &
                    "by the solver settings or not reading the record of the state "// &
                    trim(case_names(i_case))//"."
                test_oao_factory_common = .false.
            end if
            select type (oao => level%extension)
            type is (oao_type)
                ! determine if the object is set up for the density matrix of the
                ! caller and the numbers of AOs and particle channels of the case,
                ! whose dimensions the remaining checks compare arrays of, with the
                ! square roots of the overlap matrix, the density matrix transformed to
                ! the OAO basis, the Fock matrix blocks, the gradient and Hessian
                ! diagonal, the roles it performs itself, nothing evaluated yet and the
                ! expected settings stored
                n_param = n_particle * n_ao * (n_ao - 1) / 2
                if (any([oao%n_ao, oao%n_particle, oao%n_param] /= &
                        [n_ao, n_particle, n_param])) then
                    write(stderr, *) "test_oao_factory_common failed: Incorrect "// &
                        "dimensions "//trim(case_names(i_case))//"."
                    test_oao_factory_common = .false.
                else
                    if (.not. associated(oao%dm_ao, dm_ao(:, :, :n_particle))) then
                        write(stderr, *) "test_oao_factory_common failed: Density "// &
                            "matrix of the caller not taken over "// &
                            trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    end if
                    if (norm2(oao%s_sqrt - s_sqrt) > tol .or. &
                        norm2(oao%s_inv_sqrt - s_inv_sqrt) > tol) then
                        write(stderr, *) "test_oao_factory_common failed: Square "// &
                            "roots of the overlap matrix not computed "// &
                            trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    end if
                    do i = 1, n_particle
                        if (norm2(oao%dm_oao(:, :, i) - matmul( &
                            s_sqrt, matmul(dm_ao(:, :, i), s_sqrt))) > tol) then
                            write(stderr, *) "test_oao_factory_common failed: "// &
                                "Density matrix not transformed to the OAO basis "// &
                                trim(case_names(i_case))//"."
                            test_oao_factory_common = .false.
                        end if
                    end do
                    if (any( &
                        [size(oao%fock_oo, 3), size(oao%fock_vv, 3)] /= n_particle)) &
                        then
                        write(stderr, *) "test_oao_factory_common failed: Fock "// &
                            "matrix blocks not allocated for every particle "// &
                            "channel "//trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    end if
                    if (.not. all([oao%provides%update_orbs, oao%provides%obj_func, &
                                   oao%provides%hess_x, oao%provides%precond, &
                                   oao%provides%precond_pd, oao%provides%project, &
                                   oao%provides%get_extra_trial_vectors])) then
                        write(stderr, *) "test_oao_factory_common failed: Roles "// &
                            "the object performs itself not set "// &
                            trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    end if
                    if (oao%error_origin /= error_oao) then
                        write(stderr, *) "test_oao_factory_common failed: Error "// &
                            "origin of the object not set "// &
                            trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    end if
                    if (.not. (allocated(oao%grad) .and. allocated(oao%h_diag))) then
                        write(stderr, *) "test_oao_factory_common failed: Gradient "// &
                            "or Hessian diagonal not allocated "// &
                            trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    else if (any([size(oao%grad), size(oao%h_diag)] /= n_param)) then
                        write(stderr, *) "test_oao_factory_common failed: Gradient "// &
                            "or Hessian diagonal not allocated with the number of "// &
                            "parameters "//trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    end if
                    if (.not. (oao%evaluation_stale .and. oao%response_stale .and. &
                               oao%hess_eigen_stale)) then
                        write(stderr, *) "test_oao_factory_common failed: Object "// &
                            "not marked unevaluated "//trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    end if
                    if (.not. ( &
                        oao%settings == expected_settings%orbital_settings_type .and. &
                        same_logger_and_context(oao%settings, expected_settings))) then
                        write(stderr, *) "test_oao_factory_common failed: Settings "// &
                            "not stored "//trim(case_names(i_case))//"."
                        test_oao_factory_common = .false.
                    end if
                end if
                if (i_case == 2) then
                    stored = associated(oao%evaluate_dm_os, mock_evaluate_dm_os) .and. &
                             .not. associated(oao%evaluate_dm_cs)
                else if (i_case == 4) then
                    stored = .not. (associated(oao%evaluate_dm_cs) .or. &
                                    associated(oao%evaluate_dm_os))
                else
                    stored = associated(oao%evaluate_dm_cs, mock_evaluate_dm_cs) .and. &
                             .not. associated(oao%evaluate_dm_os)
                end if
                if (.not. stored) then
                    write(stderr, *) "test_oao_factory_common failed: Density "// &
                        "matrix evaluating function not stored "// &
                        trim(case_names(i_case))//"."
                    test_oao_factory_common = .false.
                end if
            class default
                write(stderr, *) "test_oao_factory_common failed: Solver settings "// &
                    "not wired with an OAO object "//trim(case_names(i_case))//"."
                test_oao_factory_common = .false.
            end select
            if ((i_case == 4) .neqv. (index( &
                context%log_message, " "//evaluate_dm_missing_warning_msg) /= 0)) then
                write(stderr, *) "test_oao_factory_common failed: Warning for a "// &
                    "missing density matrix evaluating function printed wrongly "// &
                    trim(case_names(i_case))//"."
                test_oao_factory_common = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create an object for an invalid number of particle channels, which the sanity
        ! check rejects, which has to return its error, return no functions in place of
        ! those returned before and leave the solver settings untouched
        allocate(solver_settings)
        call setup_settings(settings, context)
        evaluate_dm_cs_funptr => mock_evaluate_dm_cs
        evaluate_dm_os_funptr => null()
        update_orbs_oao_funptr => update_orbs_callback
        obj_func_oao_funptr => obj_func_callback
        call oao_factory_common(dm_ao, ao_overlap, 3_ip, n_ao, evaluate_dm_cs_funptr, &
                                evaluate_dm_os_funptr, update_orbs_oao_funptr, &
                                obj_func_oao_funptr, solver_settings, state, error, &
                                settings)
        if (error /= error_oao + 1) then
            write(stderr, *) "test_oao_factory_common failed: Error not returned "// &
                "for an invalid number of particle channels."
            test_oao_factory_common = .false.
        end if
        if (associated(update_orbs_oao_funptr) .or. associated(obj_func_oao_funptr)) &
            then
            write(stderr, *) "test_oao_factory_common failed: Functions returned "// &
                "for an invalid number of particle channels."
            test_oao_factory_common = .false.
        end if
        if (associated(solver_settings%context) .or. &
            allocated(solver_settings%extensions)) then
            write(stderr, *) "test_oao_factory_common failed: Solver settings "// &
                "changed for an invalid number of particle channels."
            test_oao_factory_common = .false.
        end if
        deallocate(solver_settings)

        ! create objects for two separate solver settings, which have to be distinct
        ! objects reading their own record of the state, so that calculations with
        ! separate settings can be used alternately
        allocate(solver_settings, other_settings)
        call oao_factory_common(dm_ao(:, :, :1), ao_overlap, 1_ip, n_ao, &
                                evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                                update_orbs_oao_funptr, obj_func_oao_funptr, &
                                solver_settings, state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_oao_factory_common failed: Produced error for "// &
                "the first of two separate solver settings."
            test_oao_factory_common = .false.
        end if
        call oao_factory_common(dm_ao(:, :, :1), ao_overlap, 1_ip, n_ao, &
                                evaluate_dm_cs_funptr, evaluate_dm_os_funptr, &
                                update_orbs_oao_funptr, obj_func_oao_funptr, &
                                other_settings, other_state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_oao_factory_common failed: Produced error for "// &
                "the second of two separate solver settings."
            test_oao_factory_common = .false.
        end if
        level => wired_level(solver_settings)
        other_level => wired_level(other_settings)
        if (.not. (associated(level) .and. associated(other_level))) then
            write(stderr, *) "test_oao_factory_common failed: Separate solver "// &
                "settings not wired."
            test_oao_factory_common = .false.
        else if (associated(level, other_level) .or. &
                 .not. associated(level%extension%state, state) .or. &
                 .not. associated(other_level%extension%state, other_state)) then
            write(stderr, *) "test_oao_factory_common failed: Separate solver "// &
                "settings not wired with distinct objects reading their own state."
            test_oao_factory_common = .false.
        end if
        deallocate(solver_settings, other_settings)

    end function test_oao_factory_common

    logical(c_bool) function test_oao_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! dimensions of orbitals parameterized in the OAO basis
        !
        use otr_oao, only: oao_settings_type, oao_sanity_check
        use opentrustregion_unit_tests, only: test_context_type, setup_error_logging

        type :: oao_sanity_case_type
            integer(ip) :: dm_ao_shape(3), ao_overlap_shape(2), n_particle, n_ao
            logical :: valid
            character(len=53) :: name
        end type

        ! every invalid case fails a single check
        type(oao_sanity_case_type), parameter :: cases(*) = [ &
            ! dm_ao_shape, ao_overlap_shape, n_particle, n_ao, valid, then name
            oao_sanity_case_type( &
                [3, 3, 1], [3, 3], 1, 3, .true., "valid closed-shell dimensions"), &
            oao_sanity_case_type( &
                [3, 3, 2], [3, 3], 2, 3, .true., "valid open-shell dimensions"), &
            oao_sanity_case_type( &
                [0, 0, 1], [0, 0], 1, 0, .false., "vanishing number of AOs"), &
            oao_sanity_case_type( &
                [1, 1, 1], [1, 1], 1, 1, .false., "single AO"), &
            oao_sanity_case_type( &
                [3, 3, 0], [3, 3], 0, 3, .false., &
                "vanishing number of particle channels"), &
            oao_sanity_case_type( &
                [3, 3, 3], [3, 3], 3, 3, .false., "three particle channels"), &
            oao_sanity_case_type( &
                [2, 2, 1], [3, 3], 1, 3, .false., &
                "density matrix with wrong number of AOs"), &
            oao_sanity_case_type( &
                [3, 3, 1], [3, 3], 2, 3, .false., &
                "density matrix with wrong number of particle channels"), &
            oao_sanity_case_type( &
                [3, 3, 1], [3, 2], 1, 3, .false., "wrongly shaped AO overlap matrix")]

        type(oao_sanity_case_type) :: c
        type(test_context_type), target :: context
        type(oao_settings_type) :: settings
        integer(ip) :: i_case, error

        ! assume tests pass
        test_oao_sanity_check = .true.

        ! call routine for every case and determine if valid dimensions are accepted
        ! without a message and invalid ones rejected with an error message
        do i_case = 1, size(cases)
            c = cases(i_case)
            call setup_error_logging(settings, context)
            call oao_sanity_check(settings, c%dm_ao_shape, c%ao_overlap_shape, &
                                  c%n_particle, c%n_ao, error)
            if (c%valid) then
                if (error /= 0 .or. context%log_message /= "") then
                    write(stderr, *) "test_oao_sanity_check failed: Valid input "// &
                        "rejected for "//trim(c%name)//"."
                    test_oao_sanity_check = .false.
                end if
            else
                if (error /= 1) then
                    write(stderr, *) "test_oao_sanity_check failed: No error "// &
                        "returned for "//trim(c%name)//"."
                    test_oao_sanity_check = .false.
                end if
                if (context%log_message == "") then
                    write(stderr, *) "test_oao_sanity_check failed: No error "// &
                        "message printed for "//trim(c%name)//"."
                    test_oao_sanity_check = .false.
                end if
            end if
        end do

    end function test_oao_sanity_check

    logical(c_bool) function test_obj_func_oao() bind(C)
        !
        ! this function tests the function which defines the energy evaluation in the
        ! OAO basis
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use test_reference, only: arm_host_context, host_context_reached
        use otr_common_unit_tests, only: &
            level_context, generate_random_density_matrix, mock_requests, &
            mock_evaluate_dm_cs, mock_evaluate_dm_os
        use otr_oao, only: oao_type, error_oao
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, shell_names

        type(test_context_type), target :: context
        real(rp), parameter :: angle = 0.3_rp

        type(oao_type), allocatable, target :: oao
        real(rp) :: s_inv_sqrt(n_ao, n_ao), rotation(n_ao, n_ao), energy, &
                    expected_energy
        real(rp), allocatable :: dm_oao(:, :, :), kappa(:)
        integer(ip) :: n_particle, n_param, i, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_obj_func_oao = .true.

        ! loop over the closed-shell and the open-shell case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            n_param = n_particle * n_ao * (n_ao - 1) / 2

            ! set up the OAO object of the shell with a non-orthonormal AO basis, so
            ! that the density matrix has to be transformed to the AO basis before it
            ! is evaluated, and a different square root of the overlap matrix, which
            ! must not be used for it
            allocate(oao, dm_oao(n_ao, n_ao, n_particle), kappa(n_param))
            call setup_settings(oao%settings, context)
            oao%error_origin = error_oao
            call arm_host_context(oao%settings, context)
            oao%n_ao = n_ao
            oao%n_particle = n_particle
            oao%n_param = n_param
            do i = 1, n_particle
                dm_oao(:, :, i) = generate_random_density_matrix(n_ao, n_occ(i))
            end do
            s_inv_sqrt = identity_matrix(n_ao) + &
                         0.1_rp * generate_random_symm_matrix(n_ao)
            oao%dm_oao = dm_oao
            oao%s_inv_sqrt = s_inv_sqrt
            oao%s_sqrt = identity_matrix(n_ao)
            if (n_particle == 1) then
                oao%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                oao%evaluate_dm_os => mock_evaluate_dm_os
            end if

            ! call routine without an orbital rotation and determine if only the energy
            ! of the density matrix evaluating function is requested and returned for
            ! the density matrix in the AO basis
            kappa = 0.0_rp
            mock_requests = [integer(ip) :: ]
            energy = oao%obj_func(kappa, error, level_context(oao, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_obj_func_oao failed: Produced error for the "// &
                    case_name//" case."
                test_obj_func_oao = .false.
            end if
            expected_energy = 0.0_rp
            do i = 1, n_particle
                expected_energy = expected_energy + sum(matmul(s_inv_sqrt, matmul( &
                    dm_oao(:, :, i), s_inv_sqrt)))
            end do
            if (abs(energy - expected_energy) > tol) then
                write(stderr, *) "test_obj_func_oao failed: Incorrect energy for "// &
                    "the "//case_name//" case."
                test_obj_func_oao = .false.
            end if
            if (size(mock_requests) /= 1) then
                write(stderr, *) "test_obj_func_oao failed: Density matrix "// &
                    "evaluating function not called exactly once for the "// &
                    case_name//" case."
                test_obj_func_oao = .false.
            end if
            if (any(mock_requests /= 0)) then
                write(stderr, *) "test_obj_func_oao failed: Incorrect outputs "// &
                    "requested from density matrix evaluating function for the "// &
                    case_name//" case."
                test_obj_func_oao = .false.
            end if

            ! call routine with an orbital rotation and determine if the current
            ! density matrix is left untouched
            call random_number(kappa)
            kappa = 0.1_rp * kappa
            energy = oao%obj_func(kappa, error, level_context(oao, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_obj_func_oao failed: Produced error for an "// &
                    "orbital rotation for the "//case_name//" case."
                test_obj_func_oao = .false.
            end if
            if (norm2(oao%dm_oao - dm_oao) > tol) then
                write(stderr, *) "test_obj_func_oao failed: Current density matrix "// &
                    "changed by an orbital rotation for the "//case_name//" case."
                test_obj_func_oao = .false.
            end if

            ! call routine with a rotation of the first pair of orbitals of the first
            ! particle channel, whose rotation matrix is a plane rotation, and
            ! determine if the energy of the rotated density matrix is returned; the
            ! closed form of the rotation catches an ignored or sign-flipped orbital
            ! rotation, which the checks above do not
            kappa = 0.0_rp
            kappa(1) = angle
            rotation = identity_matrix(n_ao)
            rotation(1:2, 1:2) = reshape([cos(angle), -sin(angle), &
                                          sin(angle), cos(angle)], [2, 2])
            energy = oao%obj_func(kappa, error, level_context(oao, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_obj_func_oao failed: Produced error for a "// &
                    "plane rotation for the "//case_name//" case."
                test_obj_func_oao = .false.
            end if
            expected_energy = sum(matmul(s_inv_sqrt, matmul(matmul( &
                transpose(rotation), matmul(dm_oao(:, :, 1), rotation)), s_inv_sqrt)))
            do i = 2, n_particle
                expected_energy = expected_energy + sum(matmul(s_inv_sqrt, matmul( &
                    dm_oao(:, :, i), s_inv_sqrt)))
            end do
            if (abs(energy - expected_energy) > tol) then
                write(stderr, *) "test_obj_func_oao failed: Incorrect energy for a "// &
                    "plane rotation for the "//case_name//" case."
                test_obj_func_oao = .false.
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("obj_func_oao", context)) &
                test_obj_func_oao = .false.

            ! call routine with a density matrix evaluating function that fails and
            ! determine if its error is returned unchanged
            context%mock_error = 3
            energy = oao%obj_func(kappa, error, level_context(oao, 1_ip))
            if (error /= 3) then
                write(stderr, *) "test_obj_func_oao failed: Error of the density "// &
                    "matrix evaluating function not returned for the "//case_name// &
                    " case."
                test_obj_func_oao = .false.
            end if
            context%mock_error = 0

            ! call routine with a rotation that is not finite, which the matrix
            ! exponential refuses, and determine if its error is returned with the
            ! origin of the object
            kappa(1) = ieee_value(1.0_rp, ieee_quiet_nan)
            energy = oao%obj_func(kappa, error, level_context(oao, 1_ip))
            if (error /= error_oao + 1) then
                write(stderr, *) "test_obj_func_oao failed: Error for a rotation "// &
                    "that is not finite not returned with the origin of the object "// &
                    "for the "//case_name//" case."
                test_obj_func_oao = .false.
            end if

            ! call routine without a density matrix evaluating function and determine
            ! if the error of its check is returned
            oao%evaluate_dm_cs => null()
            oao%evaluate_dm_os => null()
            energy = oao%obj_func(kappa, error, level_context(oao, 1_ip))
            if (error /= error_oao + 1) then
                write(stderr, *) "test_obj_func_oao failed: No error returned "// &
                    "without a density matrix evaluating function for the "// &
                    case_name//" case."
                test_obj_func_oao = .false.
            end if

            deallocate(oao, dm_oao, kappa)
        end do

    end function test_obj_func_oao

    logical(c_bool) function test_update_orbs_oao() bind(C)
        !
        ! this function tests the subroutine which defines the energy, gradient and
        ! Hessian diagonal evaluation in the OAO basis, for the closed-shell and the
        ! open-shell case
        !
        use test_reference, only: arm_host_context, host_context_reached
        use otr_common_unit_tests, only: &
            level_context, generate_random_density_matrix, mock_fock_factor, &
            mock_requests, mock_get_response_cs, mock_get_response_os, &
            mock_evaluate_dm_cs, mock_evaluate_dm_os
        use otr_oao, only: oao_type, error_oao
        use otr_common, only: hess_x_callback
        use opentrustregion_core, only: hess_x_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, shell_names

        type(test_context_type), target :: context
        type(oao_type), allocatable, target :: oao
        real(rp), allocatable, target :: dm_ao(:, :, :)
        real(rp), allocatable :: dm_oao(:, :, :), kappa(:), grad(:), h_diag(:)
        real(rp) :: s_inv_sqrt(n_ao, n_ao), fock_oao(n_ao, n_ao), &
                    rotation(n_ao, n_ao), angle, func
        integer(ip) :: n_particle, n_param, i, error
        logical :: response_set
        character(len=:), allocatable :: case_name
        procedure(hess_x_type), pointer :: hess_x_funptr

        ! assume tests pass
        test_update_orbs_oao = .true.

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
            allocate(oao)
            call setup_settings(oao%settings, context)
            oao%error_origin = error_oao
            call arm_host_context(oao%settings, context)
            oao%n_ao = n_ao
            oao%n_particle = n_particle
            oao%n_param = n_param
            oao%dm_ao => dm_ao
            oao%dm_oao = dm_oao
            oao%s_inv_sqrt = s_inv_sqrt
            allocate(oao%fock_oo(n_ao, n_ao, n_particle), &
                     oao%fock_vv(n_ao, n_ao, n_particle), oao%grad(n_param), &
                     oao%h_diag(n_param))
            if (n_particle == 1) then
                oao%evaluate_dm_cs => mock_evaluate_dm_cs
            else
                oao%evaluate_dm_os => mock_evaluate_dm_os
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
            oao%response_stale = .false.
            if (.not. check_update_orbs_oao_stage(1_ip, "a starting density")) &
                test_update_orbs_oao = .false.
            if (abs(func - sum(dm_ao)) > tol) then
                write(stderr, *) "test_update_orbs_oao failed: Incorrect energy "// &
                    "for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if
            do i = 1, n_particle
                fock_oao = matmul(s_inv_sqrt, matmul(mock_fock_factor(1) * &
                                                     dm_ao(:, :, i), s_inv_sqrt))
                if (norm2(oao%fock_oo(:, :, i) - matmul( &
                    dm_oao(:, :, i), matmul(fock_oao, dm_oao(:, :, i)))) > tol) then
                    write(stderr, *) "test_update_orbs_oao failed: Static Hessian "// &
                        "part not built from the Fock matrix in the OAO basis for "// &
                        "the "//case_name//" case."
                    test_update_orbs_oao = .false.
                end if
            end do
            if (norm2(grad - oao%grad) > tol) then
                write(stderr, *) "test_update_orbs_oao failed: Gradient not "// &
                    "returned for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if
            if (norm2(h_diag - oao%h_diag) > tol) then
                write(stderr, *) "test_update_orbs_oao failed: Hessian diagonal "// &
                    "not returned for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if
            if (n_particle == 1) then
                response_set = associated(oao%get_response_cs, mock_get_response_cs)
            else
                response_set = associated(oao%get_response_os, mock_get_response_os)
            end if
            if (.not. response_set) then
                write(stderr, *) "test_update_orbs_oao failed: Response function "// &
                    "not stored for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if
            if (.not. associated(hess_x_funptr, hess_x_callback)) then
                write(stderr, *) "test_update_orbs_oao failed: Returned Hessian "// &
                    "linear transformation is wrong for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if
            if (oao%evaluation_stale) then
                write(stderr, *) "test_update_orbs_oao failed: Evaluation still "// &
                    "marked stale after the quantities were computed for the "// &
                    case_name//" case."
                test_update_orbs_oao = .false.
            end if

            ! call routine again without an orbital rotation after clearing the outputs
            ! and setting a vanishing stored energy, and determine if the quantities of
            ! the evaluated density are returned without an evaluation, even for the
            ! vanishing energy
            oao%energy = 0.0_rp
            grad = 0.0_rp
            h_diag = 0.0_rp
            if (.not. check_update_orbs_oao_stage(1_ip, "an evaluated density")) &
                test_update_orbs_oao = .false.
            if (abs(func) > tol) then
                write(stderr, *) "test_update_orbs_oao failed: Stored energy not "// &
                    "returned for an evaluated density for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if
            if (norm2(grad - oao%grad) > tol) then
                write(stderr, *) "test_update_orbs_oao failed: Stored gradient not "// &
                    "returned for an evaluated density for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if
            if (norm2(h_diag - oao%h_diag) > tol) then
                write(stderr, *) "test_update_orbs_oao failed: Stored Hessian "// &
                    "diagonal not returned for an evaluated density for the "// &
                    case_name//" case."
                test_update_orbs_oao = .false.
            end if

            ! mark the response as stale, as an approximate-Hessian extension such as
            ! density QN would after moving the density without going through this
            ! routine, and call again without an orbital rotation, and determine if
            ! this forces an evaluation and clears the flag
            oao%response_stale = .true.
            if (.not. check_update_orbs_oao_stage(2_ip, "a stale response")) &
                test_update_orbs_oao = .false.
            if (oao%response_stale) then
                write(stderr, *) "test_update_orbs_oao failed: Response still "// &
                    "marked stale after being recomputed for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if

            ! mark the evaluation as stale, as the factory does for a new starting
            ! density, and call again without an orbital rotation, and determine if
            ! this forces an evaluation and clears the flag
            oao%evaluation_stale = .true.
            if (.not. check_update_orbs_oao_stage(3_ip, "a stale evaluation")) &
                test_update_orbs_oao = .false.
            if (oao%evaluation_stale) then
                write(stderr, *) "test_update_orbs_oao failed: Evaluation still "// &
                    "marked stale after being recomputed for the "//case_name//" case."
                test_update_orbs_oao = .false.
            end if

            ! call routine with plane rotations of the first pair of orbitals of every
            ! particle channel by negative angles at evaluated orbitals, so that only
            ! the rotation itself, whose elements sum to a negative value, calls for an
            ! evaluation, keeping the density matrix in the OAO basis it starts from,
            ! and determine if the density is rotated by the closed form of the plane
            ! rotations, consistently in the AO and the OAO basis, and evaluated, and
            ! if every evaluation requested the Fock matrix and the response function
            kappa = 0.0_rp
            do i = 1, n_particle
                kappa((i - 1) * (n_param / n_particle) + 1) = -0.1_rp * &
                                                              real(i + 2, kind=rp)
            end do
            dm_oao = oao%dm_oao
            if (.not. check_update_orbs_oao_stage(4_ip, "an orbital rotation")) &
                test_update_orbs_oao = .false.
            do i = 1, n_particle
                angle = kappa((i - 1) * (n_param / n_particle) + 1)
                rotation = identity_matrix(n_ao)
                rotation(1:2, 1:2) = reshape([cos(angle), -sin(angle), &
                                              sin(angle), cos(angle)], [2, 2])
                if (norm2(oao%dm_oao(:, :, i) - matmul( &
                    transpose(rotation), matmul(dm_oao(:, :, i), rotation))) > tol) then
                    write(stderr, *) "test_update_orbs_oao failed: Density matrix "// &
                        "not rotated correctly by an orbital rotation for the "// &
                        case_name//" case."
                    test_update_orbs_oao = .false.
                end if
            end do
            do i = 1, n_particle
                if (norm2(dm_ao(:, :, i) - matmul( &
                    s_inv_sqrt, matmul(oao%dm_oao(:, :, i), s_inv_sqrt))) > tol) then
                    write(stderr, *) "test_update_orbs_oao failed: Density matrix "// &
                        "not moved consistently in the AO and the OAO basis for "// &
                        "the "//case_name//" case."
                    test_update_orbs_oao = .false.
                end if
            end do
            if (any(mock_requests /= 3)) then
                write(stderr, *) "test_update_orbs_oao failed: Incorrect outputs "// &
                    "requested from density matrix evaluating function for the "// &
                    case_name//" case."
                test_update_orbs_oao = .false.
            end if

            ! call routine with the orbital rotation again, whose evaluation fails, and
            ! determine if the error is passed on and the rotated density is marked as
            ! not evaluated
            context%mock_error = 1
            call oao%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                 level_context(oao, 1_ip))
            if (error /= 1) then
                write(stderr, *) "test_update_orbs_oao failed: Error of the "// &
                    "density matrix evaluating function not passed on for the "// &
                    case_name//" case."
                test_update_orbs_oao = .false.
            end if
            if (.not. oao%evaluation_stale) then
                write(stderr, *) "test_update_orbs_oao failed: Evaluation not "// &
                    "marked stale after a failed evaluation for the "//case_name// &
                    " case."
                test_update_orbs_oao = .false.
            end if
            context%mock_error = 0

            ! call routine without a density matrix evaluating function and determine
            ! if the error of its check is returned
            oao%evaluate_dm_cs => null()
            oao%evaluate_dm_os => null()
            call oao%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                 level_context(oao, 1_ip))
            if (error /= error_oao + 1) then
                write(stderr, *) "test_update_orbs_oao failed: No error returned "// &
                    "without a density matrix evaluating function for the "// &
                    case_name//" case."
                test_update_orbs_oao = .false.
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("update_orbs_oao", context)) &
                test_update_orbs_oao = .false.

            ! deallocate OAO object and case arrays
            deallocate(oao, dm_ao, dm_oao, kappa, grad, h_diag)
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
            call oao%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                 level_context(oao, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_update_orbs_oao failed: Produced error for "// &
                    stage//" for the "//case_name//" case."
                passed = .false.
            end if
            if (size(mock_requests) /= n_calls) then
                write(stderr, *) "test_update_orbs_oao failed: Density matrix "// &
                    "evaluating function not called the expected number of times "// &
                    "for "//stage//" for the "//case_name//" case."
                passed = .false.
            end if

        end function check_update_orbs_oao_stage

    end function test_update_orbs_oao

    logical(c_bool) function test_hess_x_oao() bind(C)
        !
        ! this function tests the subroutine which defines the Hessian linear
        ! transformation in the OAO basis
        !
        use test_reference, only: arm_host_context, host_context_reached
        use otr_common_unit_tests, only: &
            level_context, generate_random_density_matrix, mock_requests, &
            mock_response_factor, mock_response_coupling, model_response_scale, &
            mock_get_response_cs, mock_get_coupled_response_os, &
            mock_get_stale_response_os, model_response_cs, model_response_os, &
            mock_evaluate_dm_cs, mock_evaluate_dm_os
        use otr_oao, only: oao_type, error_oao
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, shell_names

        type(test_context_type), target :: context
        type(oao_type), target :: oao
        ! step of the finite difference and its tolerance relative to the second
        ! derivative it approximates
        real(rp), parameter :: step = 1e-3_rp, fd_tol = 1e-5_rp

        integer(ip) :: n_particle, n_param
        real(rp), target :: dm_oao(n_ao, n_ao, n_particle_ref)
        real(rp) :: fock_oo(n_ao, n_ao, n_particle_ref), &
                    fock_vv(n_ao, n_ao, n_particle_ref), s_inv_sqrt(n_ao, n_ao), &
                    model_h(n_ao, n_ao), fock(n_ao, n_ao), complement(n_ao, n_ao), &
                    model_scale, second_diff
        real(rp), allocatable :: x(:), hess_x(:), expected_hess_x(:)
        integer(ip) :: i, error
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_hess_x_oao = .true.

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
        call setup_settings(oao%settings, context)
        oao%error_origin = error_oao
        call arm_host_context(oao%settings, context)
        oao%n_ao = n_ao
        oao%s_inv_sqrt = s_inv_sqrt
        oao%response_stale = .false.

        ! loop over the closed-shell and the open-shell case, for which the object
        ! holds the density matrices and Fock matrix contributions of the shell and the
        ! mock response function, which couples the particle channels in the open-shell
        ! case
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            n_param = n_particle * n_ao * (n_ao - 1) / 2
            oao%n_particle = n_particle
            oao%n_param = n_param
            oao%dm_oao = dm_oao(:, :, :n_particle)
            oao%fock_oo = fock_oo(:, :, :n_particle)
            oao%fock_vv = fock_vv(:, :, :n_particle)
            if (n_particle == 1) then
                oao%evaluate_dm_cs => mock_evaluate_dm_cs
                oao%get_response_cs => mock_get_response_cs
            else
                oao%evaluate_dm_cs => null()
                oao%get_response_cs => null()
                oao%evaluate_dm_os => mock_evaluate_dm_os
                oao%get_response_os => mock_get_coupled_response_os
            end if

            ! a random trial vector
            if (allocated(x)) deallocate(x, hess_x)
            allocate(x(n_param), hess_x(n_param))
            call random_number(x)
            expected_hess_x = ref_coupled_hess_x_oao(mock_response_coupling)

            ! call routine and determine if values of resulting Hessian linear
            ! transformation match
            call oao%hess_x(x, hess_x, error, level_context(oao, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_hess_x_oao failed: Produced error for the "// &
                    case_name//" case."
                test_hess_x_oao = .false.
            end if
            if (norm2(hess_x - expected_hess_x) > tol) then
                write(stderr, *) "test_hess_x_oao failed: Incorrect Hessian linear "// &
                    "transformation for the "//case_name//" case."
                test_hess_x_oao = .false.
            end if

            ! call routine with a response function that fails and determine if its
            ! error is returned unchanged
            context%mock_error = 3
            call oao%hess_x(x, hess_x, error, level_context(oao, 1_ip))
            if (error /= 3) then
                write(stderr, *) "test_hess_x_oao failed: Error of the response "// &
                    "function not returned for the "//case_name//" case."
                test_hess_x_oao = .false.
            end if
            context%mock_error = 0
        end do

        ! mark the response as stale, as density QN would after moving the density
        ! without going through the orbital updating routine, replace the response
        ! function by one of an outdated density and call again for the open-shell
        ! trial vector, and determine if this triggers the density matrix evaluating
        ! function to refresh the response, which is then used
        oao%dm_ao => dm_oao
        oao%evaluate_dm_os => mock_evaluate_dm_os
        oao%get_response_os => mock_get_stale_response_os
        oao%response_stale = .true.
        n_particle = n_particle_ref
        expected_hess_x = ref_coupled_hess_x_oao(0.0_rp)
        mock_requests = [integer(ip) :: ]
        call oao%hess_x(x, hess_x, error, level_context(oao, 1_ip))
        if (error /= 0) then
            write(stderr, *) "test_hess_x_oao failed: Produced error while "// &
                "refreshing a stale response."
            test_hess_x_oao = .false.
        end if
        if (size(mock_requests) /= 1) then
            write(stderr, *) "test_hess_x_oao failed: Stale response was not refreshed."
            test_hess_x_oao = .false.
        end if
        if (any(mock_requests /= 2)) then
            write(stderr, *) "test_hess_x_oao failed: Incorrect outputs requested "// &
                "from density matrix evaluating function while refreshing a stale "// &
                "response."
            test_hess_x_oao = .false.
        end if
        if (oao%response_stale) then
            write(stderr, *) "test_hess_x_oao failed: Response still marked stale "// &
                "after being refreshed."
            test_hess_x_oao = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write(stderr, *) "test_hess_x_oao failed: Incorrect Hessian linear "// &
                "transformation with a refreshed response."
            test_hess_x_oao = .false.
        end if

        ! check that the host context reached the callback functions unchanged
        if (.not. host_context_reached("hess_x_oao", context)) test_hess_x_oao = .false.

        ! call routine without a density matrix evaluating function and determine if
        ! the error of its check is returned
        oao%evaluate_dm_cs => null()
        oao%evaluate_dm_os => null()
        call oao%hess_x(x, hess_x, error, level_context(oao, 1_ip))
        if (error /= error_oao + 1) then
            write(stderr, *) "test_hess_x_oao failed: No error returned without a "// &
                "density matrix evaluating function."
            test_hess_x_oao = .false.
        end if

        ! set up the OAO object with the Fock matrix and the response of a quadratic
        ! model energy at its current density, and determine if the quadratic form of
        ! the Hessian linear transformation along a random non-redundant trial vector
        ! matches the second finite difference of the model energy along it, which
        ! checks the static and the response part together with their scaling
        ! independently of how the routine assembles them
        model_h = generate_random_symm_matrix(n_ao)
        do n_particle = 1, n_particle_ref
            case_name = trim(shell_names(n_particle))
            model_scale = model_response_scale(n_particle)
            n_param = n_particle * n_ao * (n_ao - 1) / 2
            do i = 1, n_particle
                fock = matmul(s_inv_sqrt, matmul(model_h + model_scale * matmul( &
                    s_inv_sqrt, matmul(dm_oao(:, :, i), s_inv_sqrt)), s_inv_sqrt))
                complement = identity_matrix(n_ao) - dm_oao(:, :, i)
                fock_oo(:, :, i) = matmul(dm_oao(:, :, i), &
                                          matmul(fock, dm_oao(:, :, i)))
                fock_vv(:, :, i) = matmul(complement, matmul(fock, complement))
            end do
            oao%n_particle = n_particle
            oao%n_param = n_param
            oao%dm_oao = dm_oao(:, :, :n_particle)
            oao%fock_oo = fock_oo(:, :, :n_particle)
            oao%fock_vv = fock_vv(:, :, :n_particle)
            oao%evaluate_dm_cs => null()
            oao%get_response_cs => null()
            oao%evaluate_dm_os => null()
            oao%get_response_os => null()
            if (n_particle == 1) then
                oao%evaluate_dm_cs => mock_evaluate_dm_cs
                oao%get_response_cs => model_response_cs
            else
                oao%evaluate_dm_os => mock_evaluate_dm_os
                oao%get_response_os => model_response_os
            end if
            oao%response_stale = .false.
            deallocate(x, hess_x)
            allocate(hess_x(n_param))
            x = generate_random_nonredundant_vector(n_param, n_particle, n_ao, &
                                                    dm_oao(:, :, :n_particle))
            call oao%hess_x(x, hess_x, error, level_context(oao, 1_ip))
            if (error /= 0) then
                write(stderr, *) "test_hess_x_oao failed: Produced error for the "// &
                    "quadratic model energy for the "//case_name//" case."
                test_hess_x_oao = .false.
            end if
            second_diff = (model_energy(step * x) - 2.0_rp * &
                           model_energy(0.0_rp * x) + model_energy(-step * x)) / step**2
            if (abs(dot_product(x, hess_x) - second_diff) > &
                fd_tol * max(1.0_rp, abs(second_diff))) then
                write(stderr, *) "test_hess_x_oao failed: Hessian linear "// &
                    "transformation does not match the finite difference of the "// &
                    "quadratic model energy for the "//case_name//" case."
                test_hess_x_oao = .false.
            end if
        end do

    contains

        function ref_coupled_hess_x_oao(coupling) result(expected)
            !
            ! this function returns the Hessian linear transformation of the trial
            ! vector of the current shell for the mock response of the density matrix
            ! displacement in the AO basis, which adds the given multiple of the
            ! displacements of the other particle channels to the response of every
            ! channel, transformed back to the OAO basis
            !
            real(rp), intent(in) :: coupling
            real(rp), allocatable :: expected(:)

            real(rp), allocatable :: x_full(:, :, :), delta_dm(:, :, :), &
                                     delta_dm_ao(:, :, :), response(:, :, :)
            integer(ip) :: k

            x_full = ref_unpack_asymm(x, n_particle, n_ao)
            delta_dm = ref_project_symm(x_full, dm_oao(:, :, :n_particle))
            delta_dm_ao = delta_dm
            do k = 1, n_particle
                delta_dm_ao(:, :, k) = matmul(s_inv_sqrt, &
                                              matmul(delta_dm(:, :, k), s_inv_sqrt))
            end do
            response = delta_dm
            do k = 1, n_particle
                response(:, :, k) = matmul(s_inv_sqrt, matmul( &
                    mock_response_factor * delta_dm_ao(:, :, k) + coupling * &
                    (sum(delta_dm_ao, dim=3) - delta_dm_ao(:, :, k)), s_inv_sqrt))
            end do
            expected = ref_hess_x_oao(x_full, response, dm_oao(:, :, :n_particle), &
                                      fock_oo(:, :, :n_particle), &
                                      fock_vv(:, :, :n_particle), n_param)

        end function ref_coupled_hess_x_oao

        function model_energy(kappa) result(energy)
            !
            ! this function evaluates the quadratic model energy at the current density
            ! matrix of the shell rotated by kappa, as in test_hess_x_mo, for the
            ! density matrix in the AO basis; the rotation exp(-K) D exp(K) by the
            ! antisymmetric matrix K kappa unpacks to is expanded to second order, whose
            ! error is odd in kappa and drops out of the second central difference
            !
            real(rp), intent(in) :: kappa(:)
            real(rp) :: energy

            real(rp) :: kappa_full(n_ao, n_ao, n_particle), dm(n_ao, n_ao), &
                        commutator(n_ao, n_ao)
            integer(ip) :: k

            kappa_full = ref_unpack_asymm(kappa, n_particle, n_ao)
            energy = 0.0_rp
            do k = 1, n_particle
                commutator = matmul(dm_oao(:, :, k), kappa_full(:, :, k)) - &
                             matmul(kappa_full(:, :, k), dm_oao(:, :, k))
                dm = dm_oao(:, :, k) + commutator + &
                     0.5_rp * (matmul(commutator, kappa_full(:, :, k)) - &
                               matmul(kappa_full(:, :, k), commutator))
                dm = matmul(s_inv_sqrt, matmul(dm, s_inv_sqrt))
                energy = energy + merge(2.0_rp, 1.0_rp, n_particle == 1) * &
                         (sum(model_h * dm) + 0.5_rp * model_scale * sum(dm * dm))
            end do

        end function model_energy

    end function test_hess_x_oao

    logical(c_bool) function test_project_oao() bind(C)
        !
        ! this function tests the subroutine which discards the redundant
        ! occupied-occupied and virtual-virtual rotations from a vector, retaining only
        ! its occupied-virtual and virtual-occupied contributions, and then applies the
        ! host's projection saved in the level it receives as context
        !
        use otr_oao, only: oao_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings
        use otr_oao_test_reference, only: n_param
        use otr_common_test_reference, only: n_ao, n_particle, n_occ
        use otr_common_unit_tests, only: level_context, generate_random_density_matrix
        use test_reference, only: host_context, reset_host_context, host_context_reached

        type(test_context_type), target :: context
        type(oao_type), target :: oao
        real(rp) :: dm_oao(n_ao, n_ao, n_particle), vector(n_param), &
                    projected(n_param), expected(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_project_oao = .true.

        ! generate random density matrices and vector
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_occ(1))
        dm_oao(:, :, 2) = generate_random_density_matrix(n_ao, n_occ(2))
        call random_number(vector)

        ! set up the OAO object, whose levels point back at it, as the factory sets
        ! them up, and whose level of the solver holds a host projection with the host
        ! context, which does not commute with the OAO projection, so that the order of
        ! the two shows, while the level of the stability check holds none
        call setup_settings(oao%settings, context)
        oao%n_ao = n_ao
        oao%n_particle = n_particle
        oao%dm_oao = dm_oao
        oao%solver_level%extension => oao
        oao%stability_level%extension => oao
        oao%solver_level%host%project => mock_host_project
        oao%solver_level%host%context => host_context

        ! initialize expected vector
        expected = ref_pack_asymm(ref_project_asymm( &
            ref_unpack_asymm(vector, n_particle, n_ao), dm_oao), n_param)

        ! call routine with the level of the stability check and determine if values of
        ! resulting vector match without a projection of the host
        projected = vector
        call oao%project(projected, error, level_context(oao, 2_ip))
        if (error /= 0) then
            write(stderr, *) "test_project_oao failed: Produced error without a "// &
                "projection of the host."
            test_project_oao = .false.
        end if
        if (norm2(projected - expected) > tol) then
            write(stderr, *) "test_project_oao failed: Incorrect vector without a "// &
                "projection of the host."
            test_project_oao = .false.
        end if

        ! call routine with the level of the solver and determine if the host's
        ! projection is applied after the OAO projection with the host context
        call reset_host_context()
        projected = vector
        call oao%project(projected, error, level_context(oao, 1_ip))
        expected(1) = 0.0_rp
        if (error /= 0) then
            write(stderr, *) "test_project_oao failed: Produced error with a "// &
                "projection of the host."
            test_project_oao = .false.
        end if
        if (norm2(projected - expected) > tol) then
            write(stderr, *) "test_project_oao failed: Projection of the host not "// &
                "applied after the OAO projection."
            test_project_oao = .false.
        end if
        if (.not. host_context_reached("project_oao")) test_project_oao = .false.

        ! call routine with a failing projection of the host and determine if its error
        ! is passed on
        call reset_host_context()
        host_context%mock_error = 3
        call oao%project(projected, error, level_context(oao, 1_ip))
        if (error /= 3) then
            write(stderr, *) "test_project_oao failed: Error of the projection of "// &
                "the host not passed on."
            test_project_oao = .false.
        end if

    end function test_project_oao

    logical(c_bool) function test_rotate_orbitals_oao() bind(C)
        !
        ! this function tests the subroutine which moves the current orbitals by an
        ! orbital rotation, updating the density matrix in the AO and in the OAO basis
        !
        use otr_common_unit_tests, only: generate_random_density_matrix
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use otr_oao, only: oao_type, error_oao
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix

        type(test_context_type), target :: context
        integer(ip), parameter :: n_ao = 2, n_particle = 1, n_occ = 1
        real(rp), parameter :: angle = 0.3_rp

        type(oao_type) :: oao
        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: dm_oao(n_ao, n_ao, n_particle), s_inv_sqrt(n_ao, n_ao), &
                    rotation(n_ao, n_ao), expected_dm_oao(n_ao, n_ao, n_particle), &
                    expected_dm_ao(n_ao, n_ao, n_particle)
        integer(ip) :: error

        ! assume tests pass
        test_rotate_orbitals_oao = .true.

        ! setup settings object
        call setup_settings(oao%settings, context)

        ! set up the OAO object with a random transformation to the AO basis
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_occ)
        s_inv_sqrt = identity_matrix(n_ao) + 0.1_rp * generate_random_symm_matrix(n_ao)
        dm_ao(:, :, 1) = matmul(s_inv_sqrt, matmul(dm_oao(:, :, 1), s_inv_sqrt))
        oao%dm_ao => dm_ao
        oao%dm_oao = dm_oao
        oao%s_inv_sqrt = s_inv_sqrt
        oao%response_stale = .false.
        oao%error_origin = error_oao

        ! initialize expected density matrices, where the rotation matrix is the
        ! exponential of the antisymmetric matrix the rotation is unpacked into
        rotation = reshape([cos(angle), -sin(angle), &
                            sin(angle), cos(angle)], [n_ao, n_ao])
        expected_dm_oao(:, :, 1) = matmul(transpose(rotation), &
                                          matmul(dm_oao(:, :, 1), rotation))
        expected_dm_ao(:, :, 1) = matmul(s_inv_sqrt, &
                                         matmul(expected_dm_oao(:, :, 1), s_inv_sqrt))

        ! call routine and determine if the density matrices are moved and the response
        ! is marked stale
        call oao%rotate_orbitals([angle], error)
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

        ! call routine with a rotation that is not finite, which the matrix exponential
        ! refuses, and determine if its error is returned with the origin of the object
        ! while the density matrices stay in place
        call oao%rotate_orbitals([ieee_value(1.0_rp, ieee_quiet_nan)], error)
        if (error /= error_oao + 1) then
            write(stderr, *) "test_rotate_orbitals_oao failed: Error for a "// &
                "rotation that is not finite not returned with the origin of the "// &
                "object."
            test_rotate_orbitals_oao = .false.
        end if
        if (norm2(oao%dm_oao - expected_dm_oao) > tol .or. &
            norm2(dm_ao - expected_dm_ao) > tol) then
            write(stderr, *) "test_rotate_orbitals_oao failed: Density matrices "// &
                "moved by a rotation that is not finite."
            test_rotate_orbitals_oao = .false.
        end if

    end function test_rotate_orbitals_oao

    logical(c_bool) function test_calculate_grad_h_diag_oao() bind(C)
        !
        ! this function tests the subroutine which calculates the gradient, Hessian
        ! diagonal and the occupied-occupied and virtual-virtual parts of the Fock
        ! matrix in the OAO basis
        !
        use opentrustregion_unit_tests, only: identity_matrix, &
                                              generate_random_symm_matrix
        use otr_common_unit_tests, only: generate_random_density_matrix
        use otr_oao, only: oao_type
        use otr_common_test_reference, only: n_ao, n_occ, &
                                             n_particle_ref => n_particle, shell_names

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
        use otr_oao, only: oao_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              generate_random_symm_matrix
        use otr_common_test_reference, only: n_ao, n_particle

        type(test_context_type), target :: context
        type(oao_type) :: oao
        real(rp) :: fock_oo(n_ao, n_ao, n_particle), fock_vv(n_ao, n_ao, n_particle), &
                    expected_eigvecs(n_ao, n_ao, n_particle), &
                    expected_eigvals(n_ao, n_particle)
        integer(ip) :: j, error

        ! assume tests pass
        test_refresh_hess_eigen_oao = .true.

        ! setup settings object
        call setup_settings(oao%settings, context)

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
        call oao%refresh_hess_eigen(error)
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
        call oao%refresh_hess_eigen(error)
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

        ! mark the eigendecomposition stale and determine if the one computed before is
        ! replaced by that of the changed static part
        call ref_diagonalize_static_part(fock_oo, 2.0_rp * fock_vv, expected_eigvecs, &
                                         expected_eigvals)
        oao%hess_eigen_stale = .true.
        call oao%refresh_hess_eigen(error)
        if (error /= 0) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: Produced error "// &
                "when replacing an eigendecomposition."
            test_refresh_hess_eigen_oao = .false.
        end if
        if (norm2(oao%hess_eigvals - expected_eigvals) > tol) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: "// &
                "Eigendecomposition not replaced although stale."
            test_refresh_hess_eigen_oao = .false.
        end if
        if (oao%hess_eigen_stale) then
            write(stderr, *) "test_refresh_hess_eigen_oao failed: "// &
                "Eigendecomposition still marked stale after being replaced."
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
        use otr_common_test_reference, only: n_ao, n_particle_ref => n_particle, &
                                             shell_names

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
        use otr_oao, only: oao_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              generate_random_symm_matrix
        use otr_common_unit_tests, only: level_context, &
                                         generate_random_orthogonal_matrix
        use otr_common_test_reference, only: n_particle_ref => n_particle, shell_names

        type(test_context_type), target :: context
        integer(ip), parameter :: n_ao = 3
        type(oao_type), target :: oao
        integer(ip) :: n_particle, n_param, n_extra, i, k, error
        integer(ip), allocatable :: ref_idx(:), n_occ(:)
        real(rp), allocatable :: eigvals(:, :), eigvecs(:, :, :), dm_oao(:, :, :), &
                                 trial_vectors(:, :), unit_vector(:), expected(:)
        character(len=:), allocatable :: case_name

        ! assume tests pass
        test_get_extra_trial_vectors_oao = .true.

        ! setup settings object
        call setup_settings(oao%settings, context)

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
            call oao%get_extra_trial_vectors(trial_vectors, error, &
                                             level_context(oao, 1_ip))
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
            call oao%get_extra_trial_vectors(trial_vectors, error, &
                                             level_context(oao, 1_ip))
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
        use otr_common_unit_tests, only: generate_random_density_matrix
        use otr_oao, only: rotate_dm_ao, oao_settings_type
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              identity_matrix, &
                                              generate_random_symm_matrix

        type(test_context_type), target :: context
        integer(ip), parameter :: n_ao = 2, n_particle = 2
        real(rp), parameter :: angles(n_particle) = [0.3_rp, -0.7_rp]

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), s_inv_sqrt(n_ao, n_ao), &
                    rotation(n_ao, n_ao), rotated(n_ao, n_ao), &
                    expected_dm_oao(n_ao, n_ao, n_particle), &
                    expected_dm_ao(n_ao, n_ao, n_particle), &
                    rot_dm_ao(n_ao, n_ao, n_particle), &
                    rot_dm_oao(n_ao, n_ao, n_particle)
        integer(ip) :: i, error
        type(oao_settings_type) :: settings

        ! assume tests pass
        test_rotate_dm_ao = .true.

        ! setup settings object
        call setup_settings(settings, context)

        ! set up density matrices of both particle channels whose occupations are not
        ! zero and one, so that the purification of the rotated density matrices
        ! changes them, and a random transformation to the AO basis
        do i = 1, n_particle
            dm_oao(:, :, i) = 0.8_rp * generate_random_density_matrix(n_ao, 1_ip) + &
                              0.1_rp * identity_matrix(n_ao)
        end do
        s_inv_sqrt = identity_matrix(n_ao) + 0.1_rp * generate_random_symm_matrix(n_ao)

        ! initialize expected density matrices, the density matrix of every particle
        ! channel rotated by the exponential of the antisymmetric matrix its angle is
        ! unpacked into and purified as 3 P^2 - 2 P^3
        do i = 1, n_particle
            rotation = reshape([cos(angles(i)), -sin(angles(i)), &
                                sin(angles(i)), cos(angles(i))], [n_ao, n_ao])
            rotated = matmul(transpose(rotation), matmul(dm_oao(:, :, i), rotation))
            expected_dm_oao(:, :, i) = 3.0_rp * matmul(rotated, rotated) - 2.0_rp * &
                                       matmul(rotated, matmul(rotated, rotated))
            expected_dm_ao(:, :, i) = &
                matmul(s_inv_sqrt, matmul(expected_dm_oao(:, :, i), s_inv_sqrt))
        end do

        ! call routine and determine if values of the resulting density matrices match
        call rotate_dm_ao(angles, dm_oao, s_inv_sqrt, rot_dm_ao, settings, error, &
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
        call rotate_dm_ao(angles, dm_oao, s_inv_sqrt, rot_dm_ao, settings, error)
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
        use opentrustregion_unit_tests, only: generate_random_symm_matrix
        use otr_oao, only: hess_x_static_oao
        use otr_common_test_reference, only: n_ao, n_particle_ref => n_particle, &
                                             shell_names

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
        use otr_common_test_reference, only: n_particle, n_occ, n_ao_random => n_ao
        use otr_common_unit_tests, only: generate_random_density_matrix

        integer(ip), parameter :: n_ao = 3
        real(rp) :: matrix(n_ao, n_ao, n_particle), dm_oao(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle), &
                    random_matrix(n_ao_random, n_ao_random, n_particle), &
                    random_dm(n_ao_random, n_ao_random, n_particle)
        real(rp), allocatable :: projected_matrix(:, :, :)
        integer(ip) :: i

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
            write(stderr, *) "test_project_asymm failed: Incorrect matrix dimensions."
            test_project_asymm = .false.
        else if (norm2(projected_matrix - expected) > tol) then
            write(stderr, *) "test_project_asymm failed: Incorrect matrix values."
            test_project_asymm = .false.
        end if

        ! call routine for random matrices which are not symmetric and random density
        ! matrices with the shared occupations, of which some occupy several orbitals,
        ! so that the occupied-occupied part, which the projection onto the virtual
        ! space removes, does not vanish after the antisymmetrization, and determine if
        ! the result matches the reference
        call random_number(random_matrix)
        do i = 1, n_particle
            random_dm(:, :, i) = generate_random_density_matrix(n_ao_random, n_occ(i))
        end do
        projected_matrix = project_asymm(random_matrix, random_dm)
        if (any(shape(projected_matrix) /= shape(random_matrix))) then
            write(stderr, *) "test_project_asymm failed: Incorrect matrix "// &
                "dimensions for random density matrices."
            test_project_asymm = .false.
        else if (norm2(projected_matrix - &
                       ref_project_asymm(random_matrix, random_dm)) > tol) then
            write(stderr, *) "test_project_asymm failed: Incorrect matrix values "// &
                "for random density matrices."
            test_project_asymm = .false.
        end if

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

        integer(ip), parameter :: n_ao = 2, n_particle = 2
        real(rp) :: trans_matrix(n_ao, n_ao), matrix(n_ao, n_ao, n_particle), &
                    expected(n_ao, n_ao, n_particle)
        real(rp), allocatable :: matrix_transformed(:, :, :)

        ! assume tests pass
        test_symmetric_transformation = .true.

        ! initialize transformation matrix and different matrices of both particle
        ! channels to be transformed
        trans_matrix = reshape([1.0_rp, 0.0_rp, &
                                2.0_rp, 1.0_rp], [n_ao, n_ao])
        matrix(:, :, 1) = reshape([1.0_rp, 0.0_rp, &
                                   0.0_rp, 2.0_rp], [n_ao, n_ao])
        matrix(:, :, 2) = reshape([0.0_rp, 1.0_rp, &
                                   1.0_rp, 3.0_rp], [n_ao, n_ao])

        ! initialize expected matrices
        expected(:, :, 1) = reshape([1.0_rp, 0.0_rp, &
                                     6.0_rp, 2.0_rp], [n_ao, n_ao])
        expected(:, :, 2) = reshape([2.0_rp, 1.0_rp, &
                                     11.0_rp, 5.0_rp], [n_ao, n_ao])

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
        matrix(:, :, 1) = reshape([0.0_rp, -1.0_rp, -2.0_rp, &
                                   1.0_rp, 0.0_rp, -3.0_rp, &
                                   2.0_rp, 3.0_rp, 0.0_rp], [n_ao, n_ao])
        matrix(:, :, 2) = reshape([0.0_rp, -3.0_rp, -4.0_rp, &
                                   3.0_rp, 0.0_rp, -5.0_rp, &
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
