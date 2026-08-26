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

contains

    subroutine build_metric_chol(dm_diff, chol, rank, map)
        !
        ! this subroutine constructs the ARH metric and its upper triangular Cholesky
        ! factor for linearly independent density matrix differences, so that tests
        ! which require a metric do not depend on the corresponding ARH routine
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :)
        real(rp), intent(out) :: chol(:, :, :)
        integer(ip), intent(out) :: rank(:), map(:, :)

        integer(ip) :: n_diff, i, j, k, info
        external :: dpotrf

        n_diff = size(dm_diff, 4)
        do k = 1, size(dm_diff, 3)
            ! generate metric
            do j = 1, n_diff
                do i = 1, n_diff
                    chol(i, j, k) = sum(dm_diff(:, :, k, i) * dm_diff(:, :, k, j))
                end do
            end do

            ! factorize metric in place
            call dpotrf("U", n_diff, chol(:, :, k), n_diff, info)

            ! all differences are linearly independent and keep their order
            rank(k) = n_diff
            map(:, k) = [(i, i = 1, n_diff)]
        end do

    end subroutine build_metric_chol

    subroutine build_a_inv(a, eig_vecs, eig_vals_inv)
        !
        ! this subroutine constructs the spectral decomposition of the inverse of a
        ! symmetric matrix, so that tests which require a multisecant SR1 matrix do not
        ! depend on the corresponding ARH routines
        !
        real(rp), intent(in) :: a(:, :)
        real(rp), intent(out) :: eig_vecs(:, :), eig_vals_inv(:)

        integer(ip) :: n, lwork, info
        real(rp), allocatable :: eig_vals(:), work(:)
        external :: dsyev

        n = size(a, 1)
        lwork = 10 * n
        allocate(eig_vals(n), work(lwork))

        ! symmetrize to remove numerical noise and diagonalize
        eig_vecs = 0.5_rp * (a + transpose(a))
        call dsyev("V", "U", n, eig_vecs, n, eig_vals, work, lwork, info)
        eig_vals_inv = 1.0_rp / eig_vals

        deallocate(eig_vals, work)

    end subroutine build_a_inv

    function get_pseudoinverse(eig_vecs, eig_vals_inv) result(a_inv)
        !
        ! this function assembles the pseudoinverse from the eigenvectors and inverted
        ! eigenvalues returned by the multisecant SR1 routines
        !
        real(rp), intent(in) :: eig_vecs(:, :), eig_vals_inv(:)
        real(rp) :: a_inv(size(eig_vecs, 1), size(eig_vecs, 1))

        a_inv = matmul(eig_vecs * spread(eig_vals_inv, 1, size(eig_vecs, 1)), &
                       transpose(eig_vecs))

    end function get_pseudoinverse

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
        fock = mock_factor(mock_fock_factor) * dm
        v_nonlinear = mock_factor(mock_v_nonlinear_factor) * dm

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

        n_mock_calls = n_mock_calls + 1

        error = 0
        energy = sum(dm)
        fock = mock_factor(mock_fock_factor) * dm
        v_same_spin = mock_factor(mock_v_same_spin_factor) * dm
        v_opposite_spin = mock_factor(mock_v_opposite_spin_factor) * dm
        v_nonlinear = mock_factor(mock_v_nonlinear_factor) * dm

    end subroutine mock_update_dm_os

    logical(c_bool) function test_arh_factory_cs() bind(C)
        !
        ! this function tests the subroutine which returns the modified ARH orbital
        ! updating function for the closed-shell case
        !
        use opentrustregion, only: obj_func_type, update_orbs_type, precond_type, &
                                   precond_pd_type, project_type
        use otr_arh, only: arh_factory, arh_object, arh_settings_type, &
                           update_dm_cs_type, update_orbs_arh_cs_ptr
        use otr_oao_test_reference, only: n_ao
        use otr_arh_test_reference, only: operator(==)
        use otr_oao, only: oao_object, get_energy_cs_type, obj_func_oao_ptr, &
                           precond_oao_ptr, precond_pd_oao_ptr, project_oao_ptr
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
        if (.not. associated(arh_object%get_energy_cs, mock_get_energy_cs)) then
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
        if (.not. associated(precond_arh_funptr, precond_oao_ptr)) then
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
                           update_dm_os_type, update_orbs_arh_os_ptr
        use otr_oao_test_reference, only: n_ao, n_particle
        use otr_arh_test_reference, only: operator(==)
        use otr_oao, only: oao_object, get_energy_os_type, obj_func_oao_ptr, &
                           precond_oao_ptr, precond_pd_oao_ptr, project_oao_ptr
        use otr_oao_unit_tests, only: mock_get_energy_os, identity_matrix, &
                                      generate_random_density_matrix
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_electrons = 2, &
                                  n_param = n_particle * n_ao * (n_ao - 1) / 2

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
        if (.not. associated(arh_object%get_energy_os, mock_get_energy_os)) then
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
        if (.not. associated(precond_arh_funptr, precond_oao_ptr)) then
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

    logical(c_bool) function test_get_arh_metric() bind(C)
        !
        ! this function tests the subroutine which calculates the ARH metric and its
        ! Cholesky factorization
        !
        use otr_arh, only: get_arh_metric
        use otr_oao_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 2, n_diff = 3

        real(rp) :: dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    expected_chol(n_diff, n_diff, n_particle)
        real(rp), allocatable :: chol(:, :, :)
        integer(ip) :: expected_rank(n_particle), expected_map(n_diff, n_particle)
        integer(ip), allocatable :: rank(:), map(:, :)

        ! assume tests pass
        test_get_arh_metric = .true.

        ! initialize density matrix differences, where the second difference of the
        ! first particle repeats the first one and is therefore linearly dependent
        dm_diff(:, :, 1, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 1, 2) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 1, 3) = reshape([0.0_rp, 1.0_rp, &
                                       1.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 1) = reshape([1.0_rp, 0.0_rp, &
                                       0.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 2) = reshape([0.0_rp, 1.0_rp, &
                                       1.0_rp, 0.0_rp], [n_ao, n_ao])
        dm_diff(:, :, 2, 3) = reshape([0.0_rp, 0.0_rp, &
                                       0.0_rp, 1.0_rp], [n_ao, n_ao])

        ! initialize expected numerical ranks and expected maps, where the linearly
        ! dependent difference is stored at the back while the accepted differences
        ! retain their chronological order
        expected_rank = [2, 3]
        expected_map(:, 1) = [1, 3, 2]
        expected_map(:, 2) = [1, 2, 3]

        ! initialize expected Cholesky factors, where only the accepted differences
        ! contribute
        expected_chol = 0.0_rp
        expected_chol(1, 1, 1) = 1.0_rp
        expected_chol(2, 2, 1) = sqrt(2.0_rp)
        expected_chol(1, 1, 2) = 1.0_rp
        expected_chol(2, 2, 2) = sqrt(2.0_rp)
        expected_chol(3, 3, 2) = 1.0_rp

        ! call routine and determine if dimensions, numerical ranks, maps and values of
        ! resulting Cholesky factors match
        call get_arh_metric(dm_diff, chol, rank, map)
        if (size(chol, 1) /= n_diff .or. size(chol, 2) /= n_diff .or. &
            size(chol, 3) /= n_particle .or. size(rank) /= n_particle .or. &
            size(map, 1) /= n_diff .or. size(map, 2) /= n_particle) then
            write (stderr, *) "test_get_arh_metric failed: Incorrect dimensions of "// &
                "returned metric quantities."
            test_get_arh_metric = .false.
            return
        end if
        if (any(rank /= expected_rank)) then
            write (stderr, *) "test_get_arh_metric failed: Incorrect numerical "// &
                "rank of metric."
            test_get_arh_metric = .false.
        end if
        if (any(map /= expected_map)) then
            write (stderr, *) "test_get_arh_metric failed: Incorrect map of metric."
            test_get_arh_metric = .false.
        end if
        if (norm2(chol - expected_chol) > tol) then
            write (stderr, *) "test_get_arh_metric failed: Incorrect Cholesky "// &
                "factor of metric."
            test_get_arh_metric = .false.
        end if

        ! deallocate metric quantities
        deallocate(chol, rank, map)

    end function test_get_arh_metric

    logical(c_bool) function test_multiply_with_inverse_metric() bind(C)
        !
        ! this function tests the function which multiplies a vector with the
        ! pseudoinverse of the ARH metric
        !
        use otr_arh, only: multiply_with_inverse_metric

        integer(ip), parameter :: n_diff = 2

        real(rp) :: chol(n_diff, n_diff), vec(n_diff), result_vec(n_diff), &
                    expected(n_diff)
        integer(ip) :: map(n_diff)

        ! assume tests pass
        test_multiply_with_inverse_metric = .true.

        ! initialize upper triangular Cholesky factor, which corresponds to a metric
        ! with the elements 4, 2 and 10, and vector
        chol = reshape([2.0_rp, 0.0_rp, &
                        1.0_rp, 3.0_rp], [n_diff, n_diff])
        vec = [1.0_rp, 2.0_rp]

        ! call routine for a metric of full rank and determine if values of resulting
        ! vector match
        map = [1, 2]
        expected = [1.0_rp / 6.0_rp, 1.0_rp / 6.0_rp]
        result_vec = multiply_with_inverse_metric(vec, chol, 2_ip, map)
        if (norm2(result_vec - expected) > tol) then
            write (stderr, *) "test_multiply_with_inverse_metric failed: Incorrect "// &
                "vector for metric of full rank."
            test_multiply_with_inverse_metric = .false.
        end if

        ! call routine for a rank-deficient metric and determine if the linearly
        ! dependent direction is filtered out
        expected = [0.25_rp, 0.0_rp]
        result_vec = multiply_with_inverse_metric(vec, chol, 1_ip, map)
        if (norm2(result_vec - expected) > tol) then
            write (stderr, *) "test_multiply_with_inverse_metric failed: Incorrect "// &
                "vector for rank-deficient metric."
            test_multiply_with_inverse_metric = .false.
        end if

        ! call routine for a permuted metric and determine if the permutation is
        ! applied to both the incoming and the outgoing vector
        map = [2, 1]
        expected = [0.0_rp, 0.5_rp]
        result_vec = multiply_with_inverse_metric(vec, chol, 2_ip, map)
        if (norm2(result_vec - expected) > tol) then
            write (stderr, *) "test_multiply_with_inverse_metric failed: Incorrect "// &
                "vector for permuted metric."
            test_multiply_with_inverse_metric = .false.
        end if

        ! call routine for a vanishing rank and determine if the resulting vector
        ! vanishes
        expected = [0.0_rp, 0.0_rp]
        result_vec = multiply_with_inverse_metric(vec, chol, 0_ip, map)
        if (norm2(result_vec - expected) > tol) then
            write (stderr, *) "test_multiply_with_inverse_metric failed: Incorrect "// &
                "vector for vanishing rank."
            test_multiply_with_inverse_metric = .false.
        end if

    end function test_multiply_with_inverse_metric

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
                    a_inv(n_diff, n_diff), expected_a_inv(n_diff, n_diff)
        real(rp), allocatable :: eig_vecs(:, :), eig_vals_inv(:)
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
        call get_ms_a_inv_cs(dm_diff, fock_diff, .true., eig_vecs, eig_vals_inv, n_ao, &
                             settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Produced error for "// &
                "linear part."
            test_get_ms_a_inv_cs = .false.
            return
        end if
        a_inv = get_pseudoinverse(eig_vecs, eig_vals_inv)
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
        call get_ms_a_inv_cs(dm_diff, fock_diff, .false., eig_vecs, eig_vals_inv, &
                             n_ao, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Produced error for "// &
                "non-linear part."
            test_get_ms_a_inv_cs = .false.
            return
        end if
        a_inv = get_pseudoinverse(eig_vecs, eig_vals_inv)
        if (norm2(a_inv - expected_a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Incorrect "// &
                "pseudoinverse for non-linear part."
            test_get_ms_a_inv_cs = .false.
        end if

        ! call routine for an empty history and determine if dimensions of resulting
        ! quantities vanish
        call get_ms_a_inv_cs(empty_dm_diff, empty_fock_diff, .true., eig_vecs, &
                             eig_vals_inv, n_ao, settings, error)
        if (size(eig_vecs, 1) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Incorrect eigenvector "// &
                "dimensions for empty history."
            test_get_ms_a_inv_cs = .false.
        end if
        if (size(eig_vals_inv) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_cs failed: Incorrect inverted "// &
                "eigenvalue dimensions for empty history."
            test_get_ms_a_inv_cs = .false.
        end if

        ! deallocate multisecant SR1 quantities
        deallocate(eig_vecs, eig_vals_inv)

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
                    a_inv(2 * n_diff, 2 * n_diff), &
                    expected_a_inv(2 * n_diff, 2 * n_diff)
        real(rp), allocatable :: eig_vecs(:, :), eig_vals_inv(:)
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
                                    eig_vecs, eig_vals_inv, n_ao, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Produced error."
            test_get_ms_a_inv_os_linear = .false.
            return
        end if
        a_inv = get_pseudoinverse(eig_vecs, eig_vals_inv)
        if (norm2(a_inv - expected_a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "pseudoinverse."
            test_get_ms_a_inv_os_linear = .false.
        end if

        ! call routine for an empty history and determine if dimensions of resulting
        ! quantities vanish
        call get_ms_a_inv_os_linear(empty_dm_diff, empty_v_diff, empty_v_diff, &
                                    eig_vecs, eig_vals_inv, n_ao, settings, error)
        if (size(eig_vecs, 1) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "eigenvector dimensions for empty history."
            test_get_ms_a_inv_os_linear = .false.
        end if
        if (size(eig_vals_inv) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_linear failed: Incorrect "// &
                "inverted eigenvalue dimensions for empty history."
            test_get_ms_a_inv_os_linear = .false.
        end if

        ! deallocate multisecant SR1 quantities
        deallocate(eig_vecs, eig_vals_inv)

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
                    empty_v_diff(n_ao, n_ao, n_particle, 0), a_inv(n_diff, n_diff), &
                    expected_a_inv(n_diff, n_diff)
        real(rp), allocatable :: eig_vecs(:, :), eig_vals_inv(:)
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
        call get_ms_a_inv_os_nonlinear(dm_diff, v_diff, eig_vecs, &
                                               eig_vals_inv, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Produced error."
            test_get_ms_a_inv_os_nonlinear = .false.
            return
        end if
        a_inv = get_pseudoinverse(eig_vecs, eig_vals_inv)
        if (norm2(a_inv - expected_a_inv) > tol) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Incorrect "// &
                "pseudoinverse."
            test_get_ms_a_inv_os_nonlinear = .false.
        end if

        ! call routine for an empty history and determine if dimensions of resulting
        ! quantities vanish
        call get_ms_a_inv_os_nonlinear(empty_dm_diff, empty_v_diff, eig_vecs, &
                                       eig_vals_inv, settings, error)
        if (size(eig_vecs, 1) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Incorrect "// &
                "eigenvector dimensions for empty history."
            test_get_ms_a_inv_os_nonlinear = .false.
        end if
        if (size(eig_vals_inv) /= 0) then
            write (stderr, *) "test_get_ms_a_inv_os_nonlinear failed: Incorrect "// &
                "inverted eigenvalue dimensions for empty history."
            test_get_ms_a_inv_os_nonlinear = .false.
        end if

        ! deallocate multisecant SR1 quantities
        deallocate(eig_vecs, eig_vals_inv)

    end function test_get_ms_a_inv_os_nonlinear

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
                    v_nonlinear_saved(n_ao, n_ao, n_particle), kappa(n_param), &
                    grad(n_param), h_diag(n_param), func
        integer(ip) :: error
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
        oao_object%s_sqrt = identity_matrix(n_ao)
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
        if (norm2(arh_object%fock_oao - mock_fock_factor(1) * dm_ao) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect Fock matrix."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%v_nonlinear_oao - mock_v_nonlinear_factor(1) * dm_ao) > &
            tol) then
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
        if (size(arh_object%dm_diff, 4) /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Difference history "// &
                "not initialized empty."
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

        ! save the current quantities, which the history has to retain
        dm_saved = arh_object%dm_oao
        fock_saved = arh_object%fock_oao
        v_nonlinear_saved = arh_object%v_nonlinear_oao

        ! call routine with an orbital rotation and determine if the history is
        ! extended
        kappa = 0.1_rp
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_cs(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Produced error "// &
                "after orbital rotation."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after an orbital rotation."
            test_update_orbs_arh_cs = .false.
        end if
        if (size(arh_object%dm_list, 4) /= 1) then
            write (stderr, *) "test_update_orbs_arh_cs failed: History not extended."
            test_update_orbs_arh_cs = .false.
            return
        end if
        if (norm2(arh_object%dm_list(:, :, :, 1) - dm_saved) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect density "// &
                "matrix history."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%fock_list(:, :, :, 1) - fock_saved) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect Fock "// &
                "matrix history."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%v_nonlinear_list(:, :, :, 1) - v_nonlinear_saved) > tol) &
        then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect "// &
                "non-linear potential history."
            test_update_orbs_arh_cs = .false.
        end if

        ! determine if the differences are constructed from the history and the current 
        ! quantities
        if (norm2(arh_object%dm_diff(:, :, :, 1) - (dm_saved - arh_object%dm_oao)) > &
            tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect density "// &
                "matrix difference."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%fock_diff(:, :, :, 1) - &
                  (fock_saved - arh_object%fock_oao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect Fock "// &
                "matrix difference."
            test_update_orbs_arh_cs = .false.
        end if
        if (norm2(arh_object%v_nonlinear_diff(:, :, :, 1) - &
                  (v_nonlinear_saved - arh_object%v_nonlinear_oao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect "// &
                "non-linear potential difference."
            test_update_orbs_arh_cs = .false.
        end if

        ! determine if the linear part of the potential difference is obtained as the
        ! remainder of the total Fock matrix difference
        if (norm2(arh_object%v_linear_diff - (arh_object%fock_diff - &
                                              arh_object%v_nonlinear_diff)) > tol) &
        then
            write (stderr, *) "test_update_orbs_arh_cs failed: Incorrect linear "// &
                "part of potential difference."
            test_update_orbs_arh_cs = .false.
        end if

        ! determine if the ARH metric is constructed
        if (.not. allocated(arh_object%metric_chol)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: ARH metric Cholesky "// &
                "factor not constructed."
            test_update_orbs_arh_cs = .false.
        end if
        if (.not. allocated(arh_object%metric_rank)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: ARH metric rank not "// &
                "constructed."
            test_update_orbs_arh_cs = .false.
        end if
        if (.not. allocated(arh_object%metric_map)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: ARH metric map not "// &
                "constructed."
            test_update_orbs_arh_cs = .false.
        end if

        ! call routine for multisecant SR1 and determine if the separately regularized
        ! multisecant SR1 systems are constructed instead of the ARH metric
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
        if (.not. allocated(arh_object%a_eigvecs)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Multisecant SR1 "// &
                "system not constructed."
            test_update_orbs_arh_cs = .false.
        end if
        if (.not. allocated(arh_object%a_eigvecs_comb)) then
            write (stderr, *) "test_update_orbs_arh_cs failed: Spin-combined "// &
                "multisecant SR1 system not constructed."
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
        use otr_oao_test_reference, only: n_ao, n_particle
        use otr_oao, only: oao_object
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_unit_tests, only: n_mock_calls, identity_matrix, &
                                      generate_random_density_matrix

        integer(ip), parameter :: n_electrons = 2, &
                                  n_param = n_particle * n_ao * (n_ao - 1) / 2

        real(rp), target :: dm_ao(n_ao, n_ao, n_particle)
        real(rp) :: dm_saved(n_ao, n_ao, n_particle), &
                    v_same_spin_saved(n_ao, n_ao, n_particle), &
                    v_opposite_spin_saved(n_ao, n_ao, n_particle), &
                    v_nonlinear_saved(n_ao, n_ao, n_particle), kappa(n_param), &
                    grad(n_param), h_diag(n_param), func
        integer(ip) :: i, error
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
        oao_object%s_sqrt = identity_matrix(n_ao)
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
        if (norm2(arh_object%v_same_spin_oao - mock_v_same_spin_factor(1) * dm_ao) > &
            tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect same-spin "// &
                "potential."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_opposite_spin_oao - &
                  mock_v_opposite_spin_factor(1) * dm_ao) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "opposite-spin potential."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_nonlinear_oao - mock_v_nonlinear_factor(1) * dm_ao) > &
            tol) then
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
        if (size(arh_object%dm_diff, 4) /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Difference history "// &
                "not initialized empty."
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

        ! save the current quantities, which the history has to retain
        dm_saved = arh_object%dm_oao
        v_same_spin_saved = arh_object%v_same_spin_oao
        v_opposite_spin_saved = arh_object%v_opposite_spin_oao
        v_nonlinear_saved = arh_object%v_nonlinear_oao

        ! call routine with an orbital rotation and determine if the history is
        ! extended
        kappa = 0.1_rp
        oao_object%hess_eigen_stale = .false.
        call update_orbs_arh_os(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Produced error "// &
                "after orbital rotation."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (.not. oao_object%hess_eigen_stale) then
            write (stderr, *) "test_update_orbs_arh_os failed: Cached "// &
                "eigendecomposition of the static Hessian part not marked stale "// &
                "after an orbital rotation."
            test_update_orbs_arh_os = .false.
        end if
        if (size(arh_object%dm_list, 4) /= 1) then
            write (stderr, *) "test_update_orbs_arh_os failed: History not extended."
            test_update_orbs_arh_os = .false.
            return
        end if
        if (norm2(arh_object%dm_list(:, :, :, 1) - dm_saved) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect density "// &
                "matrix history."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_same_spin_list(:, :, :, 1) - v_same_spin_saved) > tol) &
        then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect same-spin "// &
                "potential history."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_opposite_spin_list(:, :, :, 1) - v_opposite_spin_saved) &
            > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "opposite-spin potential history."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_nonlinear_list(:, :, :, 1) - v_nonlinear_saved) > tol) &
        then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "non-linear potential history."
            test_update_orbs_arh_os = .false.
        end if

        ! determine if the differences are constructed from the history and the
        ! current quantities
        if (norm2(arh_object%dm_diff(:, :, :, 1) - (dm_saved - arh_object%dm_oao)) > &
            tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect density "// &
                "matrix difference."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_same_spin_diff(:, :, :, 1) - &
                  (v_same_spin_saved - arh_object%v_same_spin_oao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect same-spin "// &
                "potential difference."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_opposite_spin_diff(:, :, :, 1) - &
                  (v_opposite_spin_saved - arh_object%v_opposite_spin_oao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "opposite-spin potential difference."
            test_update_orbs_arh_os = .false.
        end if
        if (norm2(arh_object%v_nonlinear_diff(:, :, :, 1) - &
                  (v_nonlinear_saved - arh_object%v_nonlinear_oao)) > tol) then
            write (stderr, *) "test_update_orbs_arh_os failed: Incorrect "// &
                "non-linear potential difference."
            test_update_orbs_arh_os = .false.
        end if

        ! determine if the ARH metric is constructed
        if (.not. allocated(arh_object%metric_chol)) then
            write (stderr, *) "test_update_orbs_arh_os failed: ARH metric Cholesky "// &
                "factor not constructed."
            test_update_orbs_arh_os = .false.
        end if
        if (.not. allocated(arh_object%metric_rank)) then
            write (stderr, *) "test_update_orbs_arh_os failed: ARH metric rank not "// &
                "constructed."
            test_update_orbs_arh_os = .false.
        end if
        if (.not. allocated(arh_object%metric_map)) then
            write (stderr, *) "test_update_orbs_arh_os failed: ARH metric map not "// &
                "constructed."
            test_update_orbs_arh_os = .false.
        end if

        ! call routine for multisecant SR1 and determine if the spin-separated and
        ! spin-combined multisecant SR1 systems are constructed instead of the ARH
        ! metric
        arh_object%settings%arh_type = "ms_sr1"
        call update_orbs_arh_os(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) then
            write (stderr, *) "test_update_orbs_arh_os failed: Produced error for "// &
                "multisecant SR1."
            test_update_orbs_arh_os = .false.
        end if
        if (.not. allocated(arh_object%a_eigvecs)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Spin-separated "// &
                "multisecant SR1 system not constructed."
            test_update_orbs_arh_os = .false.
        end if
        if (.not. allocated(arh_object%a_eigvecs_comb)) then
            write (stderr, *) "test_update_orbs_arh_os failed: Spin-combined "// &
                "multisecant SR1 system not constructed."
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
        use otr_oao_test_reference, only: n_ao
        use otr_oao_unit_tests, only: ref_unpack_asymm, ref_project_symm, ref_hess_x, &
                                      identity_matrix, generate_random_density_matrix, &
                                      generate_random_symm_matrix
        use opentrustregion_unit_tests, only: setup_settings

        integer(ip), parameter :: n_diff = 1

        ! inverted multisecant SR1 eigenvalues, which are chosen independently of the
        ! density matrix history since only the product of the eigenvectors and the
        ! inverted eigenvalues enters the response
        real(rp), parameter :: a_inv_eigval = 0.5_rp, a_inv_eigval_comb = 0.25_rp

        integer(ip), target :: n_particle, n_param, n_ao_target
        real(rp), target :: dm_oao(n_ao, n_ao, 2), fock_oo(n_ao, n_ao, 2), &
                            fock_vv(n_ao, n_ao, 2)
        real(rp) :: dm_diff(n_ao, n_ao, 2, n_diff), fock_diff(n_ao, n_ao, 2, n_diff), &
                    v_linear_diff(n_ao, n_ao, 2, n_diff), &
                    v_same_spin_diff(n_ao, n_ao, 2, n_diff), &
                    v_opposite_spin_diff(n_ao, n_ao, 2, n_diff), &
                    v_nonlinear_diff(n_ao, n_ao, 2, n_diff), &
                    metric_chol(n_diff, n_diff, 2), alpha(2)
        integer(ip) :: metric_rank(2), metric_map(n_diff, 2), j, error
        real(rp), allocatable :: x(:), x_full(:, :, :), delta_dm(:, :, :), &
                                 response(:, :, :), hess_x(:), expected_hess_x(:)

        ! assume tests pass
        test_hess_x_arh = .true.

        ! generate random density matrices, Fock matrix contributions, density matrix
        ! and potential matrix differences, metrics and trial vector
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, 2_ip)
        dm_oao(:, :, 2) = generate_random_density_matrix(n_ao, 1_ip)
        do j = 1, 2
            fock_oo(:, :, j) = generate_random_symm_matrix(n_ao)
            fock_vv(:, :, j) = generate_random_symm_matrix(n_ao)
            dm_diff(:, :, j, 1) = generate_random_dm_diff(dm_oao(:, :, j), n_ao)
            fock_diff(:, :, j, 1) = generate_random_symm_matrix(n_ao)
            v_linear_diff(:, :, j, 1) = generate_random_symm_matrix(n_ao)
            v_same_spin_diff(:, :, j, 1) = generate_random_symm_matrix(n_ao)
            v_opposite_spin_diff(:, :, j, 1) = generate_random_symm_matrix(n_ao)
            v_nonlinear_diff(:, :, j, 1) = generate_random_symm_matrix(n_ao)
            metric_chol(1, 1, j) = sqrt(sum(dm_diff(:, :, j, 1)**2))
            metric_rank(j) = 1
            metric_map(1, j) = 1
        end do

        ! set up the ARH object with the quantities the Hessian linear transformation
        ! requires
        allocate(arh_object)
        call setup_settings(arh_object%settings)
        n_ao_target = n_ao
        arh_object%n_ao => n_ao_target
        arh_object%n_particle => n_particle
        arh_object%n_param => n_param
        arh_object%dm_oao => dm_oao
        arh_object%fock_oo => fock_oo
        arh_object%fock_vv => fock_vv
        arh_object%a_eigvecs = identity_matrix(n_diff)
        arh_object%a_inv_eigvals = [a_inv_eigval]
        arh_object%a_eigvecs_comb = identity_matrix(n_diff)
        arh_object%a_inv_eigvals_comb = [a_inv_eigval_comb]

        ! set up the closed-shell case with standard ARH
        n_particle = 1
        n_param = n_ao * (n_ao - 1) / 2
        arh_object%settings%arh_type = "arh"
        arh_object%dm_diff = dm_diff(:, :, 1:1, :)
        arh_object%fock_diff = fock_diff(:, :, 1:1, :)
        arh_object%v_linear_diff = v_linear_diff(:, :, 1:1, :)
        arh_object%v_nonlinear_diff = v_nonlinear_diff(:, :, 1:1, :)
        arh_object%metric_chol = metric_chol(:, :, 1:1)
        arh_object%metric_rank = metric_rank(1:1)
        arh_object%metric_map = metric_map(:, 1:1)

        ! the response reproduces the Fock matrix difference scaled by the projection
        ! of the displacement onto the density matrix difference
        allocate(x(n_param))
        call random_number(x)
        x_full = ref_unpack_asymm(x, n_particle, n_ao)
        delta_dm = ref_project_symm(x_full, dm_oao(:, :, 1:1))
        alpha(1) = sum(dm_diff(:, :, 1, 1) * delta_dm(:, :, 1)) / &
                   sum(dm_diff(:, :, 1, 1)**2)
        allocate(response(n_ao, n_ao, 1))
        response(:, :, 1) = alpha(1) * fock_diff(:, :, 1, 1)
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

        ! set up the closed-shell case with multisecant SR1, where the linear part of
        ! the potential difference is treated separately from the non-linear part
        arh_object%settings%arh_type = "ms_sr1"

        ! the response reproduces the linear and non-linear parts of the potential
        ! difference scaled by the projections of the displacement onto these
        response(:, :, 1) = &
            a_inv_eigval * sum(v_linear_diff(:, :, 1, 1) * delta_dm(:, :, 1)) * &
            v_linear_diff(:, :, 1, 1) + a_inv_eigval_comb * &
            sum(v_nonlinear_diff(:, :, 1, 1) * delta_dm(:, :, 1)) * &
            v_nonlinear_diff(:, :, 1, 1)
        expected_hess_x = ref_hess_x(x_full, response, dm_oao(:, :, 1:1), &
                                     fock_oo(:, :, 1:1), fock_vv(:, :, 1:1), &
                                     n_param)

        ! call routine and determine if values of resulting Hessian linear
        ! transformation match
        call hess_x_arh(x, hess_x, error)
        if (error /= 0) then
            write (stderr, *) "test_hess_x_arh failed: Produced error for "// &
                "closed-shell case with multisecant SR1."
            test_hess_x_arh = .false.
        end if
        if (norm2(hess_x - expected_hess_x) > tol) then
            write (stderr, *) "test_hess_x_arh failed: Incorrect Hessian linear "// &
                "transformation for closed-shell case with multisecant SR1."
            test_hess_x_arh = .false.
        end if
        deallocate(x, hess_x, response)

        ! set up the open-shell case with standard ARH
        n_particle = 2
        n_param = n_particle * n_param
        arh_object%settings%arh_type = "arh"
        arh_object%dm_diff = dm_diff
        arh_object%v_same_spin_diff = v_same_spin_diff
        arh_object%v_opposite_spin_diff = v_opposite_spin_diff
        arh_object%v_nonlinear_diff = v_nonlinear_diff
        arh_object%metric_chol = metric_chol
        arh_object%metric_rank = metric_rank
        arh_object%metric_map = metric_map

        ! the response reproduces the spin-resolved potential differences scaled by the
        ! projections of the displacements onto the density matrix differences, where
        ! the non-linear part is folded into the same-spin part without any
        ! same-/opposite-spin cross-mixing
        allocate(x(n_param))
        call random_number(x)
        x_full = ref_unpack_asymm(x, n_particle, n_ao)
        delta_dm = ref_project_symm(x_full, dm_oao)
        allocate(response(n_ao, n_ao, n_particle))
        do j = 1, n_particle
            alpha(j) = sum(dm_diff(:, :, j, 1) * delta_dm(:, :, j)) / &
                       sum(dm_diff(:, :, j, 1)**2)
        end do
        do j = 1, n_particle
            response(:, :, j) = alpha(j) * (v_same_spin_diff(:, :, j, 1) + &
                                            v_nonlinear_diff(:, :, j, 1)) + &
                                alpha(3 - j) * v_opposite_spin_diff(:, :, j, 1)
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

    logical(c_bool) function test_get_response_contribution_cs() bind(C)
        !
        ! this function tests the function that computes the response contribution to 
        ! the ARH Hessian for the closed-shell case
        !
        use otr_arh, only: get_response_contribution_cs, arh_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_unit_tests, only: generate_random_density_matrix, &
                                      generate_random_symm_matrix

        integer(ip), parameter :: n_ao = 6, n_electrons = 3, n_particle = 1, n_diff = 2

        ! ARH types for which the response operator is symmetric
        character(kw_len), parameter :: symmetric_types(4) = &
            [character(kw_len) :: "symm_arh", "ms_sp", "ms_psb", "ms_sr1"]

        ! ARH types which reproduce the multisecant conditions exactly
        character(kw_len), parameter :: multisecant_types(3) = &
            [character(kw_len) :: "arh", "ms_psb", "ms_sr1"]

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    fock_diff(n_ao, n_ao, n_particle, n_diff), &
                    hess(n_ao, n_ao, n_ao, n_ao), &
                    metric_chol(n_diff, n_diff, n_particle), a(n_diff, n_diff), &
                    a_eigvecs(n_diff, n_diff), a_inv_eigvals(n_diff), &
                    delta_x(n_ao, n_ao, n_particle), delta_y(n_ao, n_ao, n_particle), &
                    g_x(n_ao, n_ao, n_particle), g_y(n_ao, n_ao, n_particle), &
                    coeff(n_diff), delta_dm(n_ao, n_ao, n_particle), &
                    expected_response(n_ao, n_ao, n_particle), &
                    response(n_ao, n_ao, n_particle)
        integer(ip) :: metric_rank(n_particle), metric_map(n_diff, n_particle)
        integer(ip) :: i, j
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_response_contribution_cs = .true.

        ! setup settings object
        call setup_settings(settings)

        ! generate random density matrix
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_electrons)

        ! generate valid density matrix differences
        do i = 1, n_diff
            dm_diff(:, :, 1, i) = generate_random_dm_diff(dm_oao(:, :, 1), n_ao)
        end do

        ! create a mock Hessian operator and generate the corresponding Fock matrix
        ! differences
        hess = generate_random_symm_hessian(n_ao)
        do i = 1, n_diff
            fock_diff(:, :, 1, i) = contract_symm_hessian(hess, dm_diff(:, :, 1, i))
        end do

        ! construct metric and its Cholesky factorization
        call build_metric_chol(dm_diff, metric_chol, metric_rank, metric_map)

        ! construct multisecant SR1 matrix and the spectral decomposition of its inverse
        do j = 1, n_diff
            do i = 1, n_diff
                a(i, j) = sum(dm_diff(:, :, 1, i) * fock_diff(:, :, 1, j))
            end do
        end do
        call build_a_inv(a, a_eigvecs, a_inv_eigvals)

        ! generate two random density matrix displacements
        delta_x(:, :, 1) = generate_random_symm_matrix(n_ao)
        delta_y(:, :, 1) = generate_random_symm_matrix(n_ao)

        ! compute the Hessian actions G(x) and G(y) and test symmetry of Hessian
        ! operator for all symmetric ARH methods
        do i = 1, size(symmetric_types)
            settings%arh_type = symmetric_types(i)
            g_x = get_response_contribution_cs(delta_x, dm_diff, fock_diff, &
                                               metric_chol, metric_rank, metric_map, &
                                               a_eigvecs, a_inv_eigvals, n_ao, &
                                               n_particle, settings)
            g_y = get_response_contribution_cs(delta_y, dm_diff, fock_diff, &
                                               metric_chol, metric_rank, metric_map, &
                                               a_eigvecs, a_inv_eigvals, n_ao, &
                                               n_particle, settings)
            if (abs(sum(delta_y * g_x) - sum(delta_x * g_y)) > tol) then
                write (stderr, *) "test_get_response_contribution_cs failed: "// &
                    "Hessian operator is not symmetric for "// &
                    trim(symmetric_types(i))//" ARH method."
                test_get_response_contribution_cs = .false.
            end if
        end do

        ! create a displacement as a linear combination of previous density matrix
        ! differences, for which the response has to reproduce the corresponding linear
        ! combination of Fock matrix differences
        call random_number(coeff)
        delta_dm = 0.0_rp
        expected_response = 0.0_rp
        do i = 1, n_diff
            delta_dm(:, :, 1) = delta_dm(:, :, 1) + coeff(i) * dm_diff(:, :, 1, i)
            expected_response(:, :, 1) = expected_response(:, :, 1) + coeff(i) * &
                                         fock_diff(:, :, 1, i)
        end do

        ! check if multisecant conditions are fulfilled
        do i = 1, size(multisecant_types)
            settings%arh_type = multisecant_types(i)
            response = get_response_contribution_cs(delta_dm, dm_diff, fock_diff, &
                                                    metric_chol, metric_rank, &
                                                    metric_map, a_eigvecs, &
                                                    a_inv_eigvals, n_ao, n_particle, &
                                                    settings)
            if (norm2(response - expected_response) > tol) then
                write (stderr, *) "test_get_response_contribution_cs failed: "// &
                    "Hessian operator does not fulfill multisecant conditions for "// &
                    trim(multisecant_types(i))//" ARH method."
                test_get_response_contribution_cs = .false.
            end if
        end do

    end function test_get_response_contribution_cs

    logical(c_bool) function test_get_response_contribution_os_separated() &
        bind(C)
        !
        ! this function tests the function that computes the spin-separated response
        ! contribution to the ARH Hessian for the open-shell case
        !
        use otr_arh, only: get_response_contribution_os_separated, &
                           arh_settings_type
        use opentrustregion_unit_tests, only: setup_settings
        use otr_oao_unit_tests, only: generate_random_density_matrix, &
                                      generate_random_symm_matrix
        use otr_oao_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 6, n_electrons = 3, n_diff = 2

        ! ARH types for which the response operator is symmetric
        character(kw_len), parameter :: symmetric_types(4) = &
            [character(kw_len) :: "symm_arh", "ms_sp", "ms_psb", "ms_sr1"]

        ! ARH types which reproduce the multisecant conditions exactly
        character(kw_len), parameter :: multisecant_types(3) = &
            [character(kw_len) :: "arh", "ms_psb", "ms_sr1"]

        ! ARH types which reproduce the multisecant conditions exactly when the
        ! non-linear potential contribution is folded in, which multisecant SR1
        ! deliberately ignores here since it requires a dedicated spin-combined system
        character(kw_len), parameter :: nonlinear_multisecant_types(2) = &
            [character(kw_len) :: "arh", "ms_psb"]

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    same_hess(n_ao, n_ao, n_ao, n_ao), &
                    opposite_hess(n_ao, n_ao, n_ao, n_ao), &
                    nonlinear_hess(n_ao, n_ao, n_ao, n_ao), &
                    same_v_diff(n_ao, n_ao, n_particle, n_diff), &
                    opposite_v_diff(n_ao, n_ao, n_particle, n_diff), &
                    nonlinear_v_diff(n_ao, n_ao, n_particle, n_diff), &
                    metric_chol(n_diff, n_diff, n_particle), &
                    a(2 * n_diff, 2 * n_diff), a_eigvecs(2 * n_diff, 2 * n_diff), &
                    a_inv_eigvals(2 * n_diff), delta_x(n_ao, n_ao, n_particle), &
                    delta_y(n_ao, n_ao, n_particle), g_x(n_ao, n_ao, n_particle), &
                    g_y(n_ao, n_ao, n_particle), coeff(n_diff), &
                    delta_dm(n_ao, n_ao, n_particle), &
                    expected_response(n_ao, n_ao, n_particle), &
                    expected_nonlinear_response(n_ao, n_ao, n_particle), &
                    response(n_ao, n_ao, n_particle)
        integer(ip) :: metric_rank(n_particle), metric_map(n_diff, n_particle)
        integer(ip) :: i, j
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_response_contribution_os_separated = .true.

        ! setup settings object
        call setup_settings(settings)

        do j = 1, n_particle
            ! generate random density matrix
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)

            ! generate valid density matrix differences
            do i = 1, n_diff
                dm_diff(:, :, j, i) = generate_random_dm_diff(dm_oao(:, :, j), n_ao)
            end do
        end do

        ! create mock Hessian operators
        same_hess = generate_random_symm_hessian(n_ao)
        opposite_hess = generate_random_symm_hessian(n_ao)
        nonlinear_hess = generate_random_symm_hessian(n_ao)

        ! generate corresponding spin-resolved potential matrix differences, where the
        ! opposite-spin potential of one spin channel is driven by the density matrix
        ! difference of the other spin channel while the non-linear potential does not
        ! mix spin channels
        do i = 1, n_diff
            do j = 1, n_particle
                same_v_diff(:, :, j, i) = contract_symm_hessian(same_hess, &
                                                                dm_diff(:, :, j, i))
                opposite_v_diff(:, :, j, i) = &
                    contract_symm_hessian(opposite_hess, dm_diff(:, :, 3 - j, i))
                nonlinear_v_diff(:, :, j, i) = &
                    contract_symm_hessian(nonlinear_hess, dm_diff(:, :, j, i))
            end do
        end do

        ! construct metrics and their Cholesky factorizations
        call build_metric_chol(dm_diff, metric_chol, metric_rank, metric_map)

        ! construct spin-combined multisecant SR1 matrix and the spectral decomposition
        ! of its inverse
        do j = 1, n_diff
            do i = 1, n_diff
                a(i, j) = sum(dm_diff(:, :, 1, i) * same_v_diff(:, :, 1, j))
                a(i, n_diff + j) = sum(dm_diff(:, :, 1, i) * &
                                       opposite_v_diff(:, :, 1, j))
                a(n_diff + i, j) = sum(dm_diff(:, :, 2, i) * &
                                       opposite_v_diff(:, :, 2, j))
                a(n_diff + i, n_diff + j) = sum(dm_diff(:, :, 2, i) * &
                                                same_v_diff(:, :, 2, j))
            end do
        end do
        call build_a_inv(a, a_eigvecs, a_inv_eigvals)

        ! generate two random density matrix displacements
        do j = 1, n_particle
            delta_x(:, :, j) = generate_random_symm_matrix(n_ao)
            delta_y(:, :, j) = generate_random_symm_matrix(n_ao)
        end do

        ! compute the Hessian actions G(x) and G(y) and test symmetry of Hessian
        ! operator for all symmetric ARH methods
        do i = 1, size(symmetric_types)
            settings%arh_type = symmetric_types(i)
            g_x = get_response_contribution_os_separated(delta_x, dm_diff, &
                                                         same_v_diff, opposite_v_diff, &
                                                         metric_chol, metric_rank, &
                                                         metric_map, a_eigvecs, &
                                                         a_inv_eigvals, n_ao, &
                                                         n_particle, settings)
            g_y = get_response_contribution_os_separated(delta_y, dm_diff, &
                                                         same_v_diff, opposite_v_diff, &
                                                         metric_chol, metric_rank, &
                                                         metric_map, a_eigvecs, &
                                                         a_inv_eigvals, n_ao, &
                                                         n_particle, settings)
            if (abs(sum(delta_y * g_x) - sum(delta_x * g_y)) > tol) then
                write (stderr, *) "test_get_response_contribution_os_separated "// &
                    "failed: Hessian operator is not symmetric for "// &
                    trim(symmetric_types(i))//" ARH method."
                test_get_response_contribution_os_separated = .false.
            end if
        end do

        ! create a displacement as a linear combination of previous density matrix
        ! differences, using the same expansion coefficients in both spin channels, for
        ! which the response has to reproduce the corresponding linear combination of
        ! potential matrix differences
        call random_number(coeff)
        delta_dm = 0.0_rp
        expected_response = 0.0_rp
        expected_nonlinear_response = 0.0_rp
        do j = 1, n_particle
            do i = 1, n_diff
                delta_dm(:, :, j) = delta_dm(:, :, j) + coeff(i) * dm_diff(:, :, j, i)
                expected_response(:, :, j) = expected_response(:, :, j) + coeff(i) * &
                                             (same_v_diff(:, :, j, i) + &
                                              opposite_v_diff(:, :, j, i))
                expected_nonlinear_response(:, :, j) = &
                    expected_nonlinear_response(:, :, j) + coeff(i) * &
                    nonlinear_v_diff(:, :, j, i)
            end do
        end do

        ! check if multisecant conditions are fulfilled
        do i = 1, size(multisecant_types)
            settings%arh_type = multisecant_types(i)
            response = &
                get_response_contribution_os_separated(delta_dm, dm_diff, same_v_diff, &
                                                       opposite_v_diff, metric_chol, &
                                                       metric_rank, metric_map, &
                                                       a_eigvecs, a_inv_eigvals, n_ao, &
                                                       n_particle, settings)
            if (norm2(response - expected_response) > tol) then
                write (stderr, *) "test_get_response_contribution_os_separated "// &
                    "failed: Hessian operator does not fulfill multisecant "// &
                    "conditions for "//trim(multisecant_types(i))//" ARH method."
                test_get_response_contribution_os_separated = .false.
            end if
        end do

        ! check if multisecant conditions are also fulfilled when the non-linear
        ! potential contribution is folded into the same-spin potential
        do i = 1, size(nonlinear_multisecant_types)
            settings%arh_type = nonlinear_multisecant_types(i)
            response = &
                get_response_contribution_os_separated(delta_dm, dm_diff, same_v_diff, &
                                                       opposite_v_diff, metric_chol, &
                                                       metric_rank, metric_map, &
                                                       a_eigvecs, a_inv_eigvals, n_ao, &
                                                       n_particle, settings, &
                                                       nonlinear_v_diff)
            if (norm2(response - expected_response - expected_nonlinear_response) > &
                tol) then
                write (stderr, *) "test_get_response_contribution_os_separated "// &
                    "failed: Hessian operator does not fulfill multisecant "// &
                    "conditions for non-linear potential contribution for "// &
                    trim(nonlinear_multisecant_types(i))//" ARH method."
                test_get_response_contribution_os_separated = .false.
            end if
        end do

    end function test_get_response_contribution_os_separated

    logical(c_bool) function test_get_response_contribution_ms_sr1_nonlinear() bind(C)
        !
        ! this function tests the function that computes the spin-combined multisecant
        ! SR1 response contribution to the ARH Hessian arising from the non-linear part
        ! of the potential
        !
        use otr_arh, only: get_response_contribution_ms_sr1_nonlinear, arh_settings_type
        use otr_oao_unit_tests, only: generate_random_density_matrix, &
                                      generate_random_symm_matrix
        use otr_oao_test_reference, only: n_particle

        integer(ip), parameter :: n_ao = 6, n_electrons = 3, n_diff = 2

        real(rp) :: dm_oao(n_ao, n_ao, n_particle), &
                    dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    hess(n_ao, n_ao, n_ao, n_ao), &
                    v_diff(n_ao, n_ao, n_particle, n_diff), a(n_diff, n_diff), &
                    a_eigvecs(n_diff, n_diff), a_inv_eigvals(n_diff), &
                    delta_x(n_ao, n_ao, n_particle), delta_y(n_ao, n_ao, n_particle), &
                    g_x(n_ao, n_ao, n_particle), g_y(n_ao, n_ao, n_particle), &
                    coeff(n_diff), delta_dm(n_ao, n_ao, n_particle), &
                    expected_response(n_ao, n_ao, n_particle), &
                    response(n_ao, n_ao, n_particle)
        integer(ip) :: i, j

        ! assume tests pass
        test_get_response_contribution_ms_sr1_nonlinear = .true.

        do j = 1, n_particle
            ! generate random density matrix
            dm_oao(:, :, j) = generate_random_density_matrix(n_ao, n_electrons)

            ! generate valid density matrix differences
            do i = 1, n_diff
                dm_diff(:, :, j, i) = generate_random_dm_diff(dm_oao(:, :, j), n_ao)
            end do
        end do

        ! create a mock Hessian operator and generate the corresponding non-linear
        ! potential matrix differences, which do not mix spin channels
        hess = generate_random_symm_hessian(n_ao)
        do i = 1, n_diff
            do j = 1, n_particle
                v_diff(:, :, j, i) = contract_symm_hessian(hess, dm_diff(:, :, j, i))
            end do
        end do

        ! construct spin-combined multisecant SR1 matrix and the spectral decomposition
        ! of its inverse
        do j = 1, n_diff
            do i = 1, n_diff
                a(i, j) = sum(dm_diff(:, :, :, i) * v_diff(:, :, :, j))
            end do
        end do
        call build_a_inv(a, a_eigvecs, a_inv_eigvals)

        ! generate two random density matrix displacements
        do j = 1, n_particle
            delta_x(:, :, j) = generate_random_symm_matrix(n_ao)
            delta_y(:, :, j) = generate_random_symm_matrix(n_ao)
        end do

        ! compute the Hessian actions G(x) and G(y) and test symmetry of Hessian
        ! operator
        g_x = get_response_contribution_ms_sr1_nonlinear(delta_x, v_diff, a_eigvecs, &
                                                         a_inv_eigvals, n_ao, &
                                                         n_particle)
        g_y = get_response_contribution_ms_sr1_nonlinear(delta_y, v_diff, a_eigvecs, &
                                                         a_inv_eigvals, n_ao, &
                                                         n_particle)
        if (abs(sum(delta_y * g_x) - sum(delta_x * g_y)) > tol) then
            write (stderr, *) "test_get_response_contribution_ms_sr1_nonlinear "// &
                "failed: Hessian operator is not symmetric."
            test_get_response_contribution_ms_sr1_nonlinear = .false.
        end if

        ! create a displacement as a linear combination of previous density matrix
        ! differences, using the same expansion coefficients in both spin channels, for
        ! which the response has to reproduce the corresponding linear combination of
        ! potential matrix differences
        call random_number(coeff)
        delta_dm = 0.0_rp
        expected_response = 0.0_rp
        do j = 1, n_particle
            do i = 1, n_diff
                delta_dm(:, :, j) = delta_dm(:, :, j) + coeff(i) * dm_diff(:, :, j, i)
                expected_response(:, :, j) = expected_response(:, :, j) + coeff(i) * &
                                             v_diff(:, :, j, i)
            end do
        end do

        ! check if multisecant conditions are fulfilled
        response = get_response_contribution_ms_sr1_nonlinear(delta_dm, v_diff, &
                                                              a_eigvecs, &
                                                              a_inv_eigvals, n_ao, &
                                                              n_particle)
        if (norm2(response - expected_response) > tol) then
            write (stderr, *) "test_get_response_contribution_ms_sr1_nonlinear "// &
                "failed: Hessian operator does not fulfill multisecant conditions."
            test_get_response_contribution_ms_sr1_nonlinear = .false.
        end if

    end function test_get_response_contribution_ms_sr1_nonlinear

end module otr_arh_unit_tests
