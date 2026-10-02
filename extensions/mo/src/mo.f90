! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo

    use opentrustregion, only: rp, ip, settings_type
    use otr_common, only: orbital_settings_type, orbital_basis_type, channel_rows

    implicit none

    ! occupations of a particle channel in the MO basis together with the
    ! occupied-occupied and virtual-virtual blocks of its Fock matrix, from which the
    ! static part of the Hessian is built, and their eigendecompositions
    type :: mo_channel_type
        integer(ip) :: n_occ = 0, n_virt = 0
        real(rp), allocatable :: fock_oo(:, :), fock_vv(:, :), occ_eigvecs(:, :), &
                                 occ_eigvals(:), virt_eigvecs(:, :), virt_eigvals(:)
    end type mo_channel_type

    ! orbitals parameterized in the MO basis, the counterpart of the OAO object for the
    ! OAO basis; unlike the OAO object, which points to the density matrix of the
    ! caller, the MO object allocates its density matrix, which its finalizer frees
    type, extends(orbital_basis_type) :: mo_type
        integer(ip) :: n_mo
        real(rp), pointer, contiguous :: mo_coeff(:, :, :) => null()
        real(rp), allocatable :: ao_overlap(:, :)
        type(mo_channel_type), allocatable :: mo_channels(:)
    contains
        procedure :: rotate_orbitals => rotate_orbitals_mo
        procedure :: calculate_grad_h_diag => calculate_grad_h_diag_mo
        procedure :: refresh_hess_eigen => refresh_hess_eigen_mo
        procedure :: rotate_to_hess_eigenbasis => rotate_to_hess_eigenbasis_mo
        procedure :: rotate_from_hess_eigenbasis => rotate_from_hess_eigenbasis_mo
        procedure :: get_hess_eigval_pairs => get_hess_eigval_pairs_mo
        procedure :: get_extra_trial_vectors => get_extra_trial_vectors_mo
        final :: finalize_mo
    end type mo_type

    ! global variables
    type(mo_type), allocatable, target :: mo_object

contains

    subroutine mo_factory_common(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                                 error, settings)
        !
        ! this subroutine performs common MO initialization operations
        !
        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        integer(ip), intent(out) :: error
        class(orbital_settings_type), intent(in) :: settings

        logical :: reuse
        integer(ip) :: i
        external :: dgemm

        ! perform sanity check
        call mo_sanity_check(settings, shape(mo_coeff, kind=ip), shape( &
            ao_overlap, kind=ip), n_occ, n_particle, n_ao, n_mo, error)
        if (error /= 0) return

        ! allocate object
        if (.not. allocated(mo_object)) allocate(mo_object)

        ! determine whether the object was already set up for the same dimensions and
        ! occupations, reading them only once it was set up at all
        reuse = allocated(mo_object%mo_channels)
        if (reuse) reuse = mo_object%n_ao == n_ao .and. mo_object%n_mo == n_mo .and. &
                           size(mo_object%mo_channels, kind=ip) == n_particle
        if (reuse) reuse = all(mo_object%mo_channels%n_occ == n_occ)

        ! set up the object anew otherwise
        if (.not. reuse) then
            ! deallocate arrays if they are already allocated
            if (allocated(mo_object%mo_channels)) deallocate(mo_object%mo_channels)
            if (associated(mo_object%dm_ao)) deallocate(mo_object%dm_ao)
            if (allocated(mo_object%grad)) deallocate(mo_object%grad)
            if (allocated(mo_object%h_diag)) deallocate(mo_object%h_diag)

            ! number of atomic orbitals, MOs and particles
            mo_object%n_ao = n_ao
            mo_object%n_mo = n_mo
            mo_object%n_particle = n_particle

            ! number of parameters is the number of occupied-virtual pairs of every
            ! particle channel
            mo_object%n_param = sum(n_occ * (n_mo - n_occ))

            ! set occupations of every particle channel
            allocate(mo_object%mo_channels(n_particle))
            mo_object%mo_channels%n_occ = n_occ
            mo_object%mo_channels%n_virt = n_mo - n_occ

            ! allocate density matrix
            allocate(mo_object%dm_ao(n_ao, n_ao, n_particle))
        end if

        ! starting MO coefficients, which are rotated in place, and AO overlap matrix
        mo_object%mo_coeff => mo_coeff
        mo_object%ao_overlap = ao_overlap

        ! construct the starting density matrix in the AO basis from the occupied
        ! orbitals
        do i = 1, n_particle
            call dgemm("N", "T", n_ao, n_ao, n_occ(i), 1.0_rp, mo_coeff(:, :, i), &
                       n_ao, mo_coeff(:, :, i), n_ao, 0.0_rp, &
                       mo_object%dm_ao(:, :, i), n_ao)
        end do

        ! nothing has been evaluated at the starting orbitals yet, so drop the static
        ! Hessian part left from a previous calculation
        mo_object%hess_eigen_stale = .true.

        ! allocate gradient and Hessian diagonal
        if (.not. allocated(mo_object%grad)) allocate(mo_object%grad(mo_object%n_param))
        if (.not. allocated(mo_object%h_diag)) &
            allocate(mo_object%h_diag(mo_object%n_param))

    end subroutine mo_factory_common

    subroutine mo_sanity_check(settings, mo_coeff_shape, ao_overlap_shape, n_occ, &
                               n_particle, n_ao, n_mo, error)
        !
        ! this subroutine performs a sanity check for the dimensions of orbitals
        ! parameterized in the MO basis, given the shapes of the MO coefficients (with
        ! the particle channels along the last dimension) and of the AO overlap matrix
        ! and the number of particle channels of the factory
        !
        use opentrustregion, only: verbosity_error

        class(orbital_settings_type), intent(in) :: settings
        integer(ip), intent(in) :: mo_coeff_shape(3), ao_overlap_shape(2), n_occ(:), &
                                   n_particle, n_ao, n_mo
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! check that number of AOs is positive
        if (n_ao < 1) then
            call settings%log("Number of AOs should be larger than 0.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that number of MOs is positive and does not exceed the number of AOs
        if (n_mo < 1 .or. n_mo > n_ao) then
            call settings%log("Number of MOs should be larger than 0 and not "// &
                              "larger than the number of AOs.", verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that there is one particle channel for the closed-shell and two for
        ! the open-shell case
        if (n_particle < 1 .or. n_particle > 2) then
            call settings%log("Number of particles should be 1 or 2.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that there is one occupation per particle channel
        if (size(n_occ) /= n_particle) then
            call settings%log("Number of occupations should match the number of "// &
                              "particle channels.", verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that the arrays have the given dimensions
        if (any(mo_coeff_shape /= [n_ao, n_mo, n_particle]) .or. &
            any(ao_overlap_shape /= [n_ao, n_ao])) then
            call settings%log("Shapes of MO coefficients and AO overlap matrix "// &
                              "should match the given dimensions.", verbosity_error, &
                              .true.)
            error = 1
            return
        end if

        ! check that the occupations fit into the MOs
        if (any(n_occ < 0) .or. any(n_occ > n_mo)) then
            call settings%log("Number of occupied orbitals should not be negative "// &
                              "and not larger than the number of MOs.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that there is at least one occupied-virtual rotation
        if (sum(n_occ * (n_mo - n_occ)) < 1) then
            call settings%log("There should be at least one occupied-virtual "// &
                              "rotation.", verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine mo_sanity_check

    subroutine mo_deconstructor()
        !
        ! this subroutine deallocates the MO objects
        !
        if (allocated(mo_object)) deallocate(mo_object)

    end subroutine mo_deconstructor

    subroutine rotate_orbitals_mo(self, kappa, settings, error)
        !
        ! this subroutine moves the current orbitals by the orbital rotation kappa,
        ! updating the MO coefficients in place together with the density matrix in the
        ! AO basis
        !
        class(mo_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: rot_mo_coeff(:, :, :), rot_dm_ao(:, :, :)

        ! rotate MO coefficients and update the density matrix accordingly
        allocate(rot_mo_coeff, mold=self%mo_coeff)
        allocate(rot_dm_ao, mold=self%dm_ao)
        call rotate_mo_coeff(kappa, self%mo_coeff, self%ao_overlap, self%mo_channels, &
                             rot_mo_coeff, rot_dm_ao, settings, error)
        if (error /= 0) return
        self%mo_coeff = rot_mo_coeff
        self%dm_ao = rot_dm_ao

    end subroutine rotate_orbitals_mo

    subroutine calculate_grad_h_diag_mo(self, fock)
        !
        ! this subroutine calculates the gradient, Hessian diagonal and static part of
        ! the Hessian at the current density from its Fock matrix in the AO basis,
        ! storing the static part as the occupied-occupied and virtual-virtual blocks
        ! of the Fock matrix in the MO basis for every particle channel
        !
        class(mo_type), intent(inout) :: self
        real(rp), intent(in) :: fock(:, :, :)

        integer(ip) :: n_occ, n_virt, offset, i, j, k
        real(rp) :: shell_scale
        real(rp), allocatable :: fock_mo(:, :)

        ! set scaling factor for closed- and open-shell systems
        shell_scale = merge(4.0_rp, 2.0_rp, self%n_particle == 1)

        offset = 0
        do k = 1, self%n_particle
            n_occ = self%mo_channels(k)%n_occ
            n_virt = self%mo_channels(k)%n_virt

            ! transform Fock matrix to the MO basis
            fock_mo = mo_transform(self%mo_coeff(:, :, k), fock(:, :, k))

            ! occupied-occupied and virtual-virtual blocks of the Fock matrix
            self%mo_channels(k)%fock_oo = fock_mo(:n_occ, :n_occ)
            self%mo_channels(k)%fock_vv = fock_mo(n_occ + 1:, n_occ + 1:)

            ! construct gradient from the occupied-virtual block
            self%grad(offset + 1:offset + n_occ * n_virt) = &
                shell_scale * reshape(fock_mo(:n_occ, n_occ + 1:), [n_occ * n_virt])

            ! construct Hessian diagonal
            do j = 1, n_virt
                do i = 1, n_occ
                    self%h_diag(offset + (j - 1) * n_occ + i) = &
                        shell_scale * (fock_mo(n_occ + j, n_occ + j) - fock_mo(i, i))
                end do
            end do
            offset = offset + n_occ * n_virt
        end do

        ! the static Hessian part was just rebuilt, so any cached eigendecomposition of
        ! it is now stale
        self%hess_eigen_stale = .true.

    end subroutine calculate_grad_h_diag_mo

    subroutine refresh_hess_eigen_mo(self, settings, error)
        !
        ! this subroutine refreshes the eigendecomposition of the static part of the
        ! Hessian in the MO basis, if the static part has changed since it was last
        ! computed, by diagonalizing the occupied-occupied and virtual-virtual blocks
        ! of the Fock matrix of every particle channel, whose eigenvectors diagonalize
        ! the static part
        !
        use opentrustregion, only: symm_mat_diag

        class(mo_type), intent(inout) :: self
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_occ, n_virt, i

        ! initialize error flag
        error = 0

        ! nothing to do if the cached eigendecomposition is up to date
        if (.not. self%hess_eigen_stale) return

        do i = 1, self%n_particle
            associate(channel => self%mo_channels(i))
                n_occ = channel%n_occ
                n_virt = channel%n_virt

                ! allocate cache arrays
                if (allocated(channel%occ_eigvecs)) &
                    deallocate(channel%occ_eigvecs, channel%occ_eigvals, &
                               channel%virt_eigvecs, channel%virt_eigvals)
                allocate( &
                    channel%occ_eigvecs(n_occ, n_occ), channel%occ_eigvals(n_occ), &
                    channel%virt_eigvecs(n_virt, n_virt), channel%virt_eigvals(n_virt))

                ! diagonalize occupied-occupied and virtual-virtual blocks
                if (n_occ > 0) then
                    call symm_mat_diag(channel%fock_oo, channel%occ_eigvals, &
                                       channel%occ_eigvecs, settings, error)
                    if (error /= 0) return
                end if
                if (n_virt > 0) then
                    call symm_mat_diag(channel%fock_vv, channel%virt_eigvals, &
                                       channel%virt_eigvecs, settings, error)
                    if (error /= 0) return
                end if
            end associate
        end do

        ! the cached eigendecomposition now matches the current static Hessian part
        self%hess_eigen_stale = .false.

    end subroutine refresh_hess_eigen_mo

    function rotate_to_hess_eigenbasis_mo(self, vector) result(rotated)
        !
        ! this function rotates a parameter vector in the MO basis into the eigenbasis
        ! of the cached static Hessian part by rotating the occupied-virtual block of
        ! every particle channel as U_o^T X U_v
        !
        class(mo_type), intent(in) :: self
        real(rp), intent(in) :: vector(:)
        real(rp), allocatable :: rotated(:)

        integer(ip) :: n_occ, n_virt, i, rows(2, self%n_particle)
        real(rp), allocatable :: ov_block(:, :), temp(:, :)
        external :: dgemm

        ! rows of every particle channel in the parameter vector
        rows = channel_rows(self%mo_channels%n_occ * self%mo_channels%n_virt)

        ! rotate each particle channel by the cached eigenvectors
        allocate(rotated(self%n_param))
        do i = 1, self%n_particle
            n_occ = self%mo_channels(i)%n_occ
            n_virt = self%mo_channels(i)%n_virt
            if (n_occ == 0 .or. n_virt == 0) cycle
            ov_block = reshape(vector(rows(1, i):rows(2, i)), [n_occ, n_virt])
            allocate(temp(n_occ, n_virt))
            call dgemm("T", "N", n_occ, n_virt, n_occ, 1.0_rp, &
                       self%mo_channels(i)%occ_eigvecs, n_occ, ov_block, n_occ, &
                       0.0_rp, temp, n_occ)
            call dgemm("N", "N", n_occ, n_virt, n_virt, 1.0_rp, temp, n_occ, &
                       self%mo_channels(i)%virt_eigvecs, n_virt, 0.0_rp, ov_block, &
                       n_occ)
            rotated(rows(1, i):rows(2, i)) = reshape(ov_block, [n_occ * n_virt])
            deallocate(ov_block, temp)
        end do

    end function rotate_to_hess_eigenbasis_mo

    function rotate_from_hess_eigenbasis_mo(self, vector) result(rotated)
        !
        ! this function rotates a parameter vector in the MO basis out of the
        ! eigenbasis of the cached static Hessian part by rotating the occupied-virtual
        ! block of every particle channel as U_o X U_v^T
        !
        class(mo_type), intent(in) :: self
        real(rp), intent(in) :: vector(:)
        real(rp), allocatable :: rotated(:)

        integer(ip) :: n_occ, n_virt, i, rows(2, self%n_particle)
        real(rp), allocatable :: ov_block(:, :), temp(:, :)
        external :: dgemm

        ! rows of every particle channel in the parameter vector
        rows = channel_rows(self%mo_channels%n_occ * self%mo_channels%n_virt)

        ! rotate each particle channel by the cached eigenvectors
        allocate(rotated(self%n_param))
        do i = 1, self%n_particle
            n_occ = self%mo_channels(i)%n_occ
            n_virt = self%mo_channels(i)%n_virt
            if (n_occ == 0 .or. n_virt == 0) cycle
            ov_block = reshape(vector(rows(1, i):rows(2, i)), [n_occ, n_virt])
            allocate(temp(n_occ, n_virt))
            call dgemm("N", "N", n_occ, n_virt, n_occ, 1.0_rp, &
                       self%mo_channels(i)%occ_eigvecs, n_occ, ov_block, n_occ, &
                       0.0_rp, temp, n_occ)
            call dgemm("N", "T", n_occ, n_virt, n_virt, 1.0_rp, temp, n_occ, &
                       self%mo_channels(i)%virt_eigvecs, n_virt, 0.0_rp, ov_block, &
                       n_occ)
            rotated(rows(1, i):rows(2, i)) = reshape(ov_block, [n_occ * n_virt])
            deallocate(ov_block, temp)
        end do

    end function rotate_from_hess_eigenbasis_mo

    function get_hess_eigval_pairs_mo(self) result(eigval_pairs)
        !
        ! this function returns the eigenvalues of the cached static Hessian part in
        ! the MO basis, the scaled differences of virtual and occupied Fock matrix
        ! block eigenvalues, packed in the same order as h_diag
        !
        class(mo_type), intent(in) :: self
        real(rp), allocatable :: eigval_pairs(:)

        integer(ip) :: offset, i, j, k
        real(rp) :: shell_scale

        ! set scaling factor for closed- and open-shell systems
        shell_scale = merge(4.0_rp, 2.0_rp, self%n_particle == 1)

        ! construct differences of eigenvalues
        allocate(eigval_pairs(self%n_param))
        offset = 0
        do k = 1, self%n_particle
            associate(channel => self%mo_channels(k))
                do j = 1, channel%n_virt
                    do i = 1, channel%n_occ
                        eigval_pairs(offset + (j - 1) * channel%n_occ + i) = &
                            shell_scale * &
                            (channel%virt_eigvals(j) - channel%occ_eigvals(i))
                    end do
                end do
                offset = offset + channel%n_occ * channel%n_virt
            end associate
        end do

    end function get_hess_eigval_pairs_mo

    subroutine get_extra_trial_vectors_mo(self, trial_vectors, settings, error)
        !
        ! this subroutine returns extra trial vectors for the MO basis along the
        ! eigenvectors of the static Hessian part with the most negative eigenvalues,
        ! the rotations between the pseudo-canonical occupied and virtual orbitals with
        ! the most negative orbital energy differences
        !
        class(mo_type), intent(inout) :: self
        real(rp), intent(out) :: trial_vectors(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: ivec, min_idx
        real(rp), allocatable :: eigval_pairs(:), unit_vector(:)

        ! initialize error flag
        error = 0

        ! vanishing vectors tell the solver that the slots left unfilled below
        ! contribute no direction
        trial_vectors = 0.0_rp

        ! refresh the eigendecomposition if the static Hessian part has changed
        call self%refresh_hess_eigen(settings, error)
        if (error /= 0) return

        ! fill the requested slots with the rotations belonging to the most negative
        ! eigenvalue pairs, stopping as soon as none are negative any more
        eigval_pairs = self%get_hess_eigval_pairs()
        allocate(unit_vector(self%n_param))
        do ivec = 1, size(trial_vectors, 2, kind=ip)
            min_idx = minloc(eigval_pairs, dim=1)
            if (eigval_pairs(min_idx) >= 0.0_rp) exit
            unit_vector = 0.0_rp
            unit_vector(min_idx) = 1.0_rp
            trial_vectors(:, ivec) = self%rotate_from_hess_eigenbasis(unit_vector)
            eigval_pairs(min_idx) = huge(1.0_rp)
        end do
        deallocate(eigval_pairs, unit_vector)

    end subroutine get_extra_trial_vectors_mo

    subroutine finalize_mo(self)
        !
        ! this subroutine frees the density matrix the MO object allocated
        !
        type(mo_type), intent(inout) :: self

        if (associated(self%dm_ao)) deallocate(self%dm_ao)

    end subroutine finalize_mo

    subroutine rotate_mo_coeff(kappa, mo_coeff, ao_overlap, mo_channels, rot_mo_coeff, &
                               rot_dm_ao, settings, error)
        !
        ! this subroutine rotates the given MO coefficients of every particle channel
        ! by the exponential of the antisymmetric matrix whose occupied-virtual blocks
        ! are given by kappa, symmetrically reorthonormalizes the rotated orbitals so
        ! that rounding errors do not accumulate over iterations, and returns them
        ! together with the resulting density matrix in the AO basis
        !
        use otr_common, only: matrix_exponential, compute_sqrt_and_inv_sqrt

        real(rp), intent(in) :: kappa(:), mo_coeff(:, :, :), ao_overlap(:, :)
        type(mo_channel_type), intent(in) :: mo_channels(:)
        real(rp), intent(out) :: rot_mo_coeff(:, :, :), rot_dm_ao(:, :, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_mo, n_occ, n_virt, i
        integer(ip), allocatable :: rows(:, :)
        real(rp), allocatable :: kappa_full(:, :), u(:, :), rotated(:, :), &
                                 overlap_rotated(:, :), metric(:, :), &
                                 metric_sqrt(:, :), metric_inv_sqrt(:, :)
        external :: dgemm

        ! initialize error flag
        error = 0

        ! number of AOs and MOs
        n_ao = size(mo_coeff, 1, kind=ip)
        n_mo = size(mo_coeff, 2, kind=ip)

        ! rows of every particle channel in the parameter vector
        rows = channel_rows(mo_channels%n_occ * mo_channels%n_virt)

        allocate(kappa_full(n_mo, n_mo), rotated(n_ao, n_mo), &
                 overlap_rotated(n_ao, n_mo), metric(n_mo, n_mo))
        do i = 1, size(mo_channels, kind=ip)
            n_occ = mo_channels(i)%n_occ
            n_virt = mo_channels(i)%n_virt

            ! construct the antisymmetric rotation generator from its occupied-virtual
            ! block
            kappa_full = 0.0_rp
            kappa_full(:n_occ, n_occ + 1:) = reshape(kappa(rows(1, i):rows(2, i)), &
                                                     [n_occ, n_virt])
            kappa_full(n_occ + 1:, :n_occ) = -transpose(kappa_full(:n_occ, n_occ + 1:))

            ! get rotation matrix
            u = matrix_exponential(kappa_full, settings, error)
            if (error /= 0) return

            ! rotate orbitals as C U^T, which rotates the density matrix in the MO
            ! basis as U^T D U
            call dgemm("N", "T", n_ao, n_mo, n_mo, 1.0_rp, mo_coeff(:, :, i), n_ao, u, &
                       n_mo, 0.0_rp, rotated, n_ao)

            ! reorthonormalize the rotated orbitals with the inverse square root of
            ! their overlap matrix
            call dgemm("N", "N", n_ao, n_mo, n_ao, 1.0_rp, ao_overlap, n_ao, rotated, &
                       n_ao, 0.0_rp, overlap_rotated, n_ao)
            call dgemm("T", "N", n_mo, n_mo, n_ao, 1.0_rp, rotated, n_ao, &
                       overlap_rotated, n_ao, 0.0_rp, metric, n_mo)
            call compute_sqrt_and_inv_sqrt(metric, metric_sqrt, metric_inv_sqrt, &
                                           settings, error)
            if (error /= 0) return
            call dgemm("N", "N", n_ao, n_mo, n_mo, 1.0_rp, rotated, n_ao, &
                       metric_inv_sqrt, n_mo, 0.0_rp, rot_mo_coeff(:, :, i), n_ao)

            ! construct the density matrix from the occupied orbitals
            call dgemm("N", "T", n_ao, n_ao, n_occ, 1.0_rp, rot_mo_coeff(:, :, i), &
                       n_ao, rot_mo_coeff(:, :, i), n_ao, 0.0_rp, rot_dm_ao(:, :, i), &
                       n_ao)
        end do

    end subroutine rotate_mo_coeff

    function mo_transform(coeff, matrix) result(transformed)
        !
        ! this function transforms a matrix given in the AO basis with the given
        ! orbital coefficients as C^T M C
        !
        real(rp), intent(in) :: coeff(:, :), matrix(:, :)
        real(rp), allocatable :: transformed(:, :)

        integer(ip) :: n_ao, n_mo
        real(rp), allocatable :: temp(:, :)
        external :: dgemm

        ! dimensions
        n_ao = size(coeff, 1, kind=ip)
        n_mo = size(coeff, 2, kind=ip)

        ! transform as C^T (M C)
        allocate(temp(n_ao, n_mo), transformed(n_mo, n_mo))
        call dgemm("N", "N", n_ao, n_mo, n_ao, 1.0_rp, matrix, n_ao, coeff, n_ao, &
                   0.0_rp, temp, n_ao)
        call dgemm("T", "N", n_mo, n_mo, n_ao, 1.0_rp, coeff, n_ao, temp, n_ao, &
                   0.0_rp, transformed, n_mo)
        deallocate(temp)

    end function mo_transform

end module otr_mo
