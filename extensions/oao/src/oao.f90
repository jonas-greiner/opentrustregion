! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao

    use opentrustregion, only: rp, ip, settings_type, obj_func_type, update_orbs_type, &
                               hess_x_type, precond_type, precond_pd_type, &
                               project_type, get_extra_trial_vectors_type, &
                               solver_settings_type
    use otr_common, only: orbital_settings_type, default_orbital_settings, &
                          orbital_basis_type, evaluate_dm_cs_type, evaluate_dm_os_type

    implicit none

    type, extends(orbital_settings_type) :: oao_settings_type
    contains
        procedure :: init => init_oao_settings
    end type oao_settings_type

    type(oao_settings_type), parameter :: default_oao_settings = &
        oao_settings_type(orbital_settings_type = default_orbital_settings)

    ! orbitals parameterized in the OAO basis, which points to the density matrix of
    ! the caller, which it updates in place
    type, extends(orbital_basis_type) :: oao_type
        real(rp), allocatable :: s_sqrt(:, :), s_inv_sqrt(:, :), dm_oao(:, :, :), &
                                 fock_oo(:, :, :), fock_vv(:, :, :), &
                                 hess_eigvecs(:, :, :), hess_eigvals(:, :)
    contains
        procedure :: rotate_orbitals => rotate_orbitals_oao
        procedure :: calculate_grad_h_diag => calculate_grad_h_diag_oao
        procedure :: refresh_hess_eigen => refresh_hess_eigen_oao
        procedure :: rotate_to_hess_eigenbasis => rotate_to_hess_eigenbasis_oao
        procedure :: rotate_from_hess_eigenbasis => rotate_from_hess_eigenbasis_oao
        procedure :: get_hess_eigval_pairs => get_hess_eigval_pairs_oao
        procedure :: get_extra_trial_vectors => get_extra_trial_vectors_oao
    end type oao_type

    ! global variables
    type(oao_type), allocatable, target :: oao_object

    ! create function pointers to ensure that routines comply with interface
    procedure(obj_func_type), pointer :: obj_func_oao_callback_ptr => &
        obj_func_oao_callback
    procedure(update_orbs_type), pointer :: update_orbs_oao_callback_ptr => &
        update_orbs_oao_callback
    procedure(hess_x_type), pointer :: hess_x_oao_callback_ptr => hess_x_oao_callback
    procedure(precond_type), pointer :: precond_oao_callback_ptr => precond_oao_callback
    procedure(precond_pd_type), pointer :: precond_pd_oao_callback_ptr => &
        precond_pd_oao_callback
    procedure(project_type), pointer :: project_oao_callback_ptr => project_oao_callback
    procedure(get_extra_trial_vectors_type), pointer :: &
        get_extra_trial_vectors_oao_callback_ptr => get_extra_trial_vectors_oao_callback

    ! define module procedures for different spin cases
    interface oao_factory
        module procedure oao_factory_cs, oao_factory_os
    end interface oao_factory

contains

    subroutine oao_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_cs, &
                              obj_func_oao_funptr, update_orbs_oao_funptr, &
                              solver_settings, error, settings)
        !
        ! this function returns a modified OAO orbital updating function for the
        ! closed-shell case and wires the OAO preconditioners, projection and extra
        ! trial vectors into the solver settings; the density matrix is updated in
        ! place, so it has to outlive the calculation
        !
        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_cs
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(inout) :: settings

        real(rp), pointer, contiguous :: dm_ao_3d(:, :, :)

        ! initialize error flag
        error = 0

        ! call common setup
        dm_ao_3d(1:n_ao, 1:n_ao, 1:1) => dm_ao
        call oao_factory_common(dm_ao_3d, ao_overlap, n_particle, n_ao, error, settings)
        if (error /= 0) return
        nullify(dm_ao_3d)

        ! set pointers to functions
        oao_object%evaluate_dm_cs => evaluate_dm_cs
        oao_object%evaluate_dm_os => null()

        ! get pointers to modified function
        obj_func_oao_funptr => obj_func_oao_callback
        update_orbs_oao_funptr => update_orbs_oao_callback

        ! wire the remaining OAO routines into the solver settings
        call oao_set_solver_settings(solver_settings, error)

    end subroutine oao_factory_cs

    subroutine oao_factory_os(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_os, &
                              obj_func_oao_funptr, update_orbs_oao_funptr, &
                              solver_settings, error, settings)
        !
        ! this function returns a modified OAO orbital updating function for the
        ! open-shell case and wires the OAO preconditioners, projection and extra trial
        ! vectors into the solver settings; the density matrix is updated in place, so
        ! it has to outlive the calculation
        !
        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(inout) :: settings

        ! initialize error flag
        error = 0

        ! call common setup
        call oao_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        if (error /= 0) return

        ! set pointers to functions
        oao_object%evaluate_dm_os => evaluate_dm_os
        oao_object%evaluate_dm_cs => null()

        ! get pointers to modified function
        obj_func_oao_funptr => obj_func_oao_callback
        update_orbs_oao_funptr => update_orbs_oao_callback

        ! wire the remaining OAO routines into the solver settings
        call oao_set_solver_settings(solver_settings, error)

    end subroutine oao_factory_os

    subroutine oao_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        !
        ! this subroutine performs common OAO initialization operations
        !
        use otr_common, only: compute_sqrt_and_inv_sqrt

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        integer(ip), intent(out) :: error
        class(orbital_settings_type), intent(in) :: settings

        logical :: reuse

        ! perform sanity check
        call oao_sanity_check(settings, n_particle, n_ao, error)
        if (error /= 0) return

        ! allocate object
        if (.not. allocated(oao_object)) allocate(oao_object)

        ! set (potentially new) settings
        oao_object%settings = settings

        ! determine whether the object was already set up for the same dimensions,
        ! which the square roots of the overlap matrix, allocated only once computed,
        ! indicate, reading the dimensions only once it was set up at all
        reuse = allocated(oao_object%s_inv_sqrt)
        if (reuse) &
            reuse = oao_object%n_particle == n_particle .and. oao_object%n_ao == n_ao

        ! set up the object anew otherwise
        if (.not. reuse) then
            ! deallocate arrays if they are already allocated
            if (allocated(oao_object%s_sqrt)) deallocate(oao_object%s_sqrt)
            if (allocated(oao_object%s_inv_sqrt)) deallocate(oao_object%s_inv_sqrt)
            if (allocated(oao_object%dm_oao)) deallocate(oao_object%dm_oao)
            if (allocated(oao_object%fock_oo)) deallocate(oao_object%fock_oo)
            if (allocated(oao_object%fock_vv)) deallocate(oao_object%fock_vv)
            if (allocated(oao_object%hess_eigvecs)) &
                deallocate(oao_object%hess_eigvecs, oao_object%hess_eigvals)
            if (allocated(oao_object%grad)) deallocate(oao_object%grad)
            if (allocated(oao_object%h_diag)) deallocate(oao_object%h_diag)

            ! number of particles
            oao_object%n_particle = n_particle

            ! get number of atomic orbitals
            oao_object%n_ao = n_ao

            ! get number of non-redundant parameters
            oao_object%n_param = n_particle * n_ao * (n_ao - 1) / 2

            ! get square root and inverse square root of AO overlap matrix
            call compute_sqrt_and_inv_sqrt(ao_overlap, oao_object%s_sqrt, &
                                           oao_object%s_inv_sqrt, oao_object%settings, &
                                           error)
            if (error /= 0) return

            ! allocate matrices
            allocate(oao_object%fock_oo(n_ao, n_ao, oao_object%n_particle), &
                     oao_object%fock_vv(n_ao, n_ao, oao_object%n_particle))
        end if

        ! starting density matrix, also in the orthogonalized AO basis
        oao_object%dm_ao => dm_ao
        oao_object%dm_oao = symmetric_transformation(oao_object%s_sqrt, dm_ao)

        ! nothing has been evaluated at the starting density yet, so drop any
        ! quantities left from a previous calculation
        oao_object%evaluation_stale = .true.
        oao_object%response_stale = .true.
        oao_object%hess_eigen_stale = .true.
        oao_object%get_response_cs => null()
        oao_object%get_response_os => null()

        ! allocate gradient and Hessian diagonal
        if (.not. allocated(oao_object%grad)) &
            allocate(oao_object%grad(oao_object%n_param))
        if (.not. allocated(oao_object%h_diag)) &
            allocate(oao_object%h_diag(oao_object%n_param))

    end subroutine oao_factory_common

    subroutine oao_sanity_check(settings, n_particle, n_ao, error)
        !
        ! this subroutine performs a sanity check for OAO input parameters
        !
        use opentrustregion, only: verbosity_error, string_to_lowercase

        class(orbital_settings_type), intent(in) :: settings
        integer(ip), intent(in) :: n_particle, n_ao
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

        ! check that there is one particle channel for the closed-shell and two for
        ! the open-shell case
        if (n_particle < 1 .or. n_particle > 2) then
            call settings%log("Number of particles should be 1 or 2.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine oao_sanity_check

    subroutine oao_set_solver_settings(solver_settings, error)
        !
        ! this subroutine wires the OAO preconditioners, projection and extra trial
        ! vectors into the solver settings and those of its stability check
        !
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! initialize settings
        if (.not. solver_settings%initialized) then
            call solver_settings%init(error)
            if (error /= 0) return
        end if

        ! set callback functions of the solver and its stability check
        solver_settings%precond => precond_oao_callback
        solver_settings%precond_pd => precond_pd_oao_callback
        solver_settings%project => project_oao_callback
        solver_settings%get_extra_trial_vectors => get_extra_trial_vectors_oao_callback
        solver_settings%stability_settings%precond => precond_oao_callback
        solver_settings%stability_settings%project => project_oao_callback
        solver_settings%stability_settings%get_extra_trial_vectors => &
            get_extra_trial_vectors_oao_callback

    end subroutine oao_set_solver_settings

    function obj_func_oao_callback(kappa, error) result(energy)
        !
        ! this function defines the energy evaluation in the OAO basis
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        real(rp), allocatable :: rot_dm_ao(:, :, :)

        ! initialize energy in case of error
        energy = 0.0_rp

        ! get rotated density matrix in AO basis
        allocate(rot_dm_ao(oao_object%n_ao, oao_object%n_ao, oao_object%n_particle))
        call rotate_dm_ao(kappa, oao_object%dm_oao, oao_object%s_inv_sqrt, rot_dm_ao, &
                          oao_object%settings, error)
        if (error /= 0) return

        ! calculate the mean-field energy
        if (associated(oao_object%evaluate_dm_os)) then
            call oao_object%evaluate_dm_os(rot_dm_ao, energy, error=error)
        else
            call oao_object%evaluate_dm_cs(rot_dm_ao(:, :, 1), energy, error=error)
        end if
        if (error /= 0) return

    end function obj_func_oao_callback

    subroutine update_orbs_oao_callback(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation
        ! and the Hessian linear transformation in the OAO basis
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle
        real(rp), allocatable :: fock_ao(:, :, :)

        ! initialize error flag
        error = 0

        ! evaluate at the rotated orbitals, or at the current ones if these have not
        ! been evaluated yet or their response is stale
        if ((sum(abs(kappa)) > 0.0_rp) .or. oao_object%evaluation_stale .or. &
            oao_object%response_stale) then
            ! number of AOs
            n_ao = oao_object%n_ao

            ! number of particles
            n_particle = oao_object%n_particle

            ! rotate orbitals, which have not been evaluated until this succeeds
            oao_object%evaluation_stale = .true.
            call oao_object%rotate_orbitals(kappa, oao_object%settings, error)
            if (error /= 0) return

            ! get energy, Fock matrix, and response function
            allocate(fock_ao(n_ao, n_ao, n_particle))
            if (associated(oao_object%evaluate_dm_os)) then
                call oao_object%evaluate_dm_os(oao_object%dm_ao, oao_object%energy, &
                                               fock_ao, oao_object%get_response_os, &
                                               error)
            else
                call oao_object%evaluate_dm_cs(oao_object%dm_ao(:, :, 1), &
                                               oao_object%energy, fock_ao(:, :, 1), &
                                               oao_object%get_response_cs, error)
            end if
            if (error /= 0) then
                deallocate(fock_ao)
                return
            end if

            ! calculate gradient, Hessian diagonal and static part of the Hessian from
            ! the Fock matrix in the OAO basis
            call oao_object%calculate_grad_h_diag( &
                symmetric_transformation(oao_object%s_inv_sqrt, fock_ao))
            deallocate(fock_ao)

            ! the rotated density matrix and its response have been evaluated
            oao_object%evaluation_stale = .false.
            oao_object%response_stale = .false.
        end if

        ! set outputs
        func = oao_object%energy
        grad = oao_object%grad
        h_diag = oao_object%h_diag
        hess_x_funptr => hess_x_oao_callback
        
    end subroutine update_orbs_oao_callback

    subroutine hess_x_oao_callback(x, hess_x, error)
        !
        ! this function defines the Hessian linear transformation in the OAO basis, the
        ! static part from the occupied-occupied and virtual-virtual parts of the Fock
        ! matrix together with the response of the Fock matrix to the density matrix
        ! response to the trial vector, projected onto the occupied-virtual and
        ! virtual-occupied subspace
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, n_param
        real(rp), allocatable :: x_full(:, :, :), dm_response(:, :, :), &
                                 fock_response(:, :, :), hess_x_full(:, :, :)

        ! initialize error flag
        error = 0

        ! number of AOs
        n_ao = oao_object%n_ao

        ! number of particles
        n_particle = oao_object%n_particle

        ! number of parameters
        n_param = oao_object%n_param

        ! rebuild the response if the density was moved without it being updated, since
        ! the static and response parts would otherwise refer to different points
        if (oao_object%response_stale) then
            call oao_object%refresh_response(error)
            if (error /= 0) return
        end if

        ! unpack trial vector
        x_full = unpack_asymm(x, n_particle, n_ao)

        ! get static part
        hess_x_full = hess_x_static_oao(x_full, oao_object%fock_oo, oao_object%fock_vv)

        ! get density matrix response to trial vector
        dm_response = project_symm(x_full, oao_object%dm_oao)
        deallocate(x_full)

        ! transform density matrix from OAO basis to AO basis
        dm_response = symmetric_transformation(oao_object%s_inv_sqrt, dm_response)

        ! get response of Fock matrix to density matrix response
        allocate(fock_response(n_ao, n_ao, n_particle))
        if (associated(oao_object%get_response_os)) then
            call oao_object%get_response_os(dm_response, fock_response, error)
        else
            call oao_object%get_response_cs(dm_response(:, :, 1), &
                                            fock_response(:, :, 1), error)
        end if
        if (error /= 0) return

        ! transform Fock response to OAO basis
        fock_response = symmetric_transformation(oao_object%s_inv_sqrt, fock_response)

        ! project the combined static and response contributions onto the
        ! occupied-virtual and virtual-occupied subspace; the static part is already
        ! confined to that subspace in exact arithmetic, but the projection on the full
        ! Hessian linear transformation is free since the Fock response has to be
        ! projected anyways and the full projection can prevent some numerical leakage
        ! into the redundant subspace
        hess_x_full = project_asymm(hess_x_full + fock_response, oao_object%dm_oao)
        deallocate(fock_response)

        ! pack Hessian linear transformation
        if (oao_object%n_particle == 1) then
            hess_x = 4.0_rp * pack_asymm(hess_x_full, oao_object%n_param)
        else
            hess_x = 2.0_rp * pack_asymm(hess_x_full, oao_object%n_param)
        end if
        deallocate(hess_x_full)

    end subroutine hess_x_oao_callback

    subroutine project_oao_callback(vector, error)
        !
        ! this subroutine discards the redundant occupied-occupied and virtual-virtual
        ! rotations from a vector describing an orbital rotation in-place, retaining
        ! only its occupied-virtual and virtual-occupied contributions
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: vector_full(:, :, :), projected_vector_full(:, :, :)

        ! initialize error flag
        error = 0

        ! unpack vector
        vector_full = unpack_asymm(vector, oao_object%n_particle, oao_object%n_ao)

        ! project vector onto occupied-virtual and virtual-occupied subspace
        projected_vector_full = project_asymm(vector_full, oao_object%dm_oao)
        deallocate(vector_full)

        ! pack vector
        vector = pack_asymm(projected_vector_full, size(vector, kind=ip))
        deallocate(projected_vector_full)

    end subroutine project_oao_callback

    subroutine precond_oao_callback(residual, mu, precond_residual, error)
        !
        ! this subroutine defines a level-shifted preconditioner based on the exact
        ! eigendecomposition of the static part of the Hessian
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        call oao_object%precond(residual, mu, precond_residual, oao_object%settings, &
                                error)

    end subroutine precond_oao_callback

    subroutine precond_pd_oao_callback(residual, precond_residual, error)
        !
        ! this subroutine defines the positive-definite preconditioner based on the
        ! exact eigendecomposition of the static part of the Hessian
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        call oao_object%precond_pd(residual, precond_residual, oao_object%settings, &
                                   error)

    end subroutine precond_pd_oao_callback

    subroutine get_extra_trial_vectors_oao_callback(trial_vectors, error)
        !
        ! this subroutine returns the extra trial vectors of the OAO basis for the
        ! solver's initial trial space
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error

        call oao_object%get_extra_trial_vectors(trial_vectors, oao_object%settings, &
                                                error)

    end subroutine get_extra_trial_vectors_oao_callback

    subroutine init_oao_settings(self, error)
        !
        ! this subroutine initializes the OAO settings
        !
        use opentrustregion, only: verbosity_error

        class(oao_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        select type(settings => self)
        type is (oao_settings_type)
            settings = default_oao_settings
        class default
            call settings%log("Orthogonal atomic orbital settings could not be "// &
                              "initialized because initialization routine received "// &
                              "the wrong type. The type oao_settings_type was "// &
                              "likely subclassed without providing an "// &
                              "initialization routine.", verbosity_error, .true.)
            error = 1
        end select

    end subroutine init_oao_settings

    subroutine oao_deconstructor()
        !
        ! this subroutine deallocates the OAO objects
        !
        if (allocated(oao_object)) deallocate(oao_object)

    end subroutine oao_deconstructor

    subroutine rotate_orbitals_oao(self, kappa, settings, error)
        !
        ! this subroutine moves the current orbitals by the orbital rotation kappa,
        ! updating the density matrix in the AO and in the OAO basis
        !
        class(oao_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: rot_dm_ao(:, :, :), rot_dm_oao(:, :, :)

        ! rotate density matrix
        allocate(rot_dm_ao, mold=self%dm_ao)
        allocate(rot_dm_oao, mold=self%dm_oao)
        call rotate_dm_ao(kappa, self%dm_oao, self%s_inv_sqrt, rot_dm_ao, settings, &
                          error, rot_dm_oao)
        if (error /= 0) return
        self%dm_ao = rot_dm_ao
        self%dm_oao = rot_dm_oao

        ! the density was moved but the response was not rebuilt
        self%response_stale = .true.

    end subroutine rotate_orbitals_oao

    subroutine calculate_grad_h_diag_oao(self, fock)
        !
        ! this subroutine calculates the gradient, Hessian diagonal and static part of
        ! the Hessian at the current density from its Fock matrix in the OAO basis,
        ! storing the static part as the occupied-occupied and virtual-virtual parts of
        ! the Fock matrix
        !
        class(oao_type), intent(inout) :: self
        real(rp), intent(in) :: fock(:, :, :)

        integer(ip) :: n_ao, n_particle, i, j, k, idx
        real(rp), allocatable :: dm_fock_oao(:, :), fock_dm_oao(:, :), &
                                 fock_ov(:, :, :), fock_vo(:, :, :), grad_full(:, :, :)
        external :: dgemm

        ! number of AOs and particles
        n_ao = self%n_ao
        n_particle = self%n_particle

        ! get contributions to Fock matrix based on occupancies
        allocate(dm_fock_oao(n_ao, n_ao), fock_dm_oao(n_ao, n_ao), &
                 fock_ov(n_ao, n_ao, n_particle), fock_vo(n_ao, n_ao, n_particle))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, self%dm_oao(:, :, i), n_ao, &
                       fock(:, :, i), n_ao, 0.0_rp, dm_fock_oao, n_ao)
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_fock_oao, n_ao, &
                       self%dm_oao(:, :, i), n_ao, 0.0_rp, self%fock_oo(:, :, i), &
                       n_ao) ! DFD
            fock_ov(:, :, i) = dm_fock_oao - self%fock_oo(:, :, i) ! DF(I-D)
            fock_vo(:, :, i) = transpose(fock_ov(:, :, i)) ! (I_D)FD
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, fock(:, :, i), n_ao, &
                       self%dm_oao(:, :, i), n_ao, 0.0_rp, fock_dm_oao, n_ao)
            self%fock_vv(:, :, i) = fock(:, :, i) - dm_fock_oao - fock_dm_oao + &
                                    self%fock_oo(:, :, i) ! (I_D)F(I_D)
        end do
        deallocate(dm_fock_oao, fock_dm_oao)

        ! construct gradient
        if (n_particle == 1) then
            grad_full = 4.0_rp * (fock_ov - fock_vo)
        else
            grad_full = 2.0_rp * (fock_ov - fock_vo)
        end if
        deallocate(fock_ov, fock_vo)

        ! pack gradient
        self%grad = pack_asymm(grad_full, self%n_param)

        ! construct Hessian diagonal
        idx = 1
        do k = 1, n_particle
            do j = 1, n_ao
                do i = 1, j - 1
                    self%h_diag(idx) = self%fock_vv(i, i, k) + self%fock_vv(j, j, k) - &
                                       self%fock_oo(i, i, k) - self%fock_oo(j, j, k)
                    idx = idx + 1
                end do
            end do
        end do
        if (n_particle == 1) then
            self%h_diag = 4.0_rp * self%h_diag
        else
            self%h_diag = 2.0_rp * self%h_diag
        end if

        ! the static Hessian part was just rebuilt, so any cached eigendecomposition of
        ! it is now stale
        self%hess_eigen_stale = .true.

    end subroutine calculate_grad_h_diag_oao

    subroutine refresh_hess_eigen_oao(self, settings, error)
        !
        ! this subroutine diagonalizes the static part of the Hessian for each particle
        ! channel and caches the result if the static part has changed since it was
        ! last diagonalized; this operator is symmetric and, since the
        ! occupied-occupied/virtual-virtual blocks of the Fock matrix are sandwiched
        ! between the (idempotent) density matrix and its complement, its eigenspaces
        ! coincide with the occupied and virtual subspaces even though
        ! occupied-occupied and virtual-virtual parts of the Fock matrix are not
        ! individually diagonal in the OAO basis
        !
        use opentrustregion, only: verbosity_error

        class(oao_type), intent(inout) :: self
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, lwork, info
        real(rp), allocatable :: a(:, :), work(:)
        character(300) :: msg
        external :: dsyev

        ! initialize error flag
        error = 0

        ! keep the cached eigendecomposition if the static part has not changed
        if (.not. self%hess_eigen_stale) return

        ! number of AOs and particles
        n_ao = self%n_ao
        n_particle = self%n_particle

        ! allocate cache arrays if necessary
        if (.not. allocated(self%hess_eigvecs)) allocate(self%hess_eigvecs( &
            n_ao, n_ao, n_particle), self%hess_eigvals(n_ao, n_particle))

        allocate(a(n_ao, n_ao))
        do i = 1, n_particle
            ! form the static part for this particle channel, dsyev overwrites it
            a = self%fock_vv(:, :, i) - self%fock_oo(:, :, i)

            ! query optimal workspace size
            lwork = -1
            allocate(work(1))
            call dsyev("V", "U", n_ao, a, n_ao, self%hess_eigvals(:, i), work, lwork, &
                       info)
            lwork = int(work(1))
            deallocate(work)
            allocate(work(lwork))

            ! perform eigendecomposition
            call dsyev("V", "U", n_ao, a, n_ao, self%hess_eigvals(:, i), work, lwork, &
                       info)
            deallocate(work)

            ! check for successful execution
            if (info /= 0) then
                write (msg, '(A, I0)') "Eigendecomposition of static Hessian part "// &
                    "failed: Error in DSYEV, info = ", info
                call settings%log(msg, verbosity_error, .true.)
                error = 1
                deallocate(a)
                return
            end if

            self%hess_eigvecs(:, :, i) = a
        end do
        deallocate(a)

        ! the cached eigendecomposition now matches the current static Hessian part
        self%hess_eigen_stale = .false.

    end subroutine refresh_hess_eigen_oao

    function rotate_to_hess_eigenbasis_oao(self, vector) result(rotated)
        !
        ! this function rotates a packed antisymmetric orbital-rotation vector into the
        ! eigenbasis of the cached static Hessian part, one particle channel at a time
        !
        class(oao_type), intent(in) :: self
        real(rp), intent(in) :: vector(:)
        real(rp), allocatable :: rotated(:)

        integer(ip) :: n_ao, n_particle, i
        real(rp), allocatable :: full(:, :, :), rotated_full(:, :, :), temp(:, :)
        external :: dgemm

        ! number of AOs and particles
        n_ao = self%n_ao
        n_particle = self%n_particle

        ! unpack vector
        full = unpack_asymm(vector, n_particle, n_ao)

        ! rotate each particle channel by the cached eigenvectors
        allocate(rotated_full(n_ao, n_ao, n_particle), temp(n_ao, n_ao))
        do i = 1, n_particle
            call dgemm("T", "N", n_ao, n_ao, n_ao, 1.0_rp, self%hess_eigvecs(:, :, i), &
                       n_ao, full(:, :, i), n_ao, 0.0_rp, temp, n_ao)
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, temp, n_ao, &
                       self%hess_eigvecs(:, :, i), n_ao, 0.0_rp, &
                       rotated_full(:, :, i), n_ao)
        end do
        deallocate(full, temp)

        ! pack vector
        rotated = pack_asymm(rotated_full, self%n_param)
        deallocate(rotated_full)

    end function rotate_to_hess_eigenbasis_oao

    function rotate_from_hess_eigenbasis_oao(self, vector) result(rotated)
        !
        ! this function rotates a packed antisymmetric orbital-rotation vector out of
        ! the eigenbasis of the cached static Hessian part, one particle channel at a
        ! time
        !
        class(oao_type), intent(in) :: self
        real(rp), intent(in) :: vector(:)
        real(rp), allocatable :: rotated(:)

        integer(ip) :: n_ao, n_particle, i
        real(rp), allocatable :: full(:, :, :), rotated_full(:, :, :), temp(:, :)
        external :: dgemm

        ! number of AOs and particles
        n_ao = self%n_ao
        n_particle = self%n_particle

        ! unpack vector
        full = unpack_asymm(vector, n_particle, n_ao)

        ! rotate each particle channel by the cached eigenvectors
        allocate(rotated_full(n_ao, n_ao, n_particle), temp(n_ao, n_ao))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, self%hess_eigvecs(:, :, i), &
                       n_ao, full(:, :, i), n_ao, 0.0_rp, temp, n_ao)
            call dgemm("N", "T", n_ao, n_ao, n_ao, 1.0_rp, temp, n_ao, &
                       self%hess_eigvecs(:, :, i), n_ao, 0.0_rp, &
                       rotated_full(:, :, i), n_ao)
        end do
        deallocate(full, temp)

        ! pack vector
        rotated = pack_asymm(rotated_full, self%n_param)
        deallocate(rotated_full)

    end function rotate_from_hess_eigenbasis_oao

    function get_hess_eigval_pairs_oao(self) result(eigval_pairs)
        !
        ! this function returns the pairwise sums of the cached static-Hessian-part
        ! eigenvalues, packed in the same order as h_diag
        !
        class(oao_type), intent(in) :: self
        real(rp), allocatable :: eigval_pairs(:)

        integer(ip) :: n_ao, n_particle, i, j, k, idx

        ! number of AOs and particles
        n_ao = self%n_ao
        n_particle = self%n_particle

        ! construct pairwise sums of eigenvalues
        allocate(eigval_pairs(self%n_param))
        idx = 1
        do k = 1, n_particle
            do j = 1, n_ao
                do i = 1, j - 1
                    eigval_pairs(idx) = self%hess_eigvals(i, k) + &
                                        self%hess_eigvals(j, k)
                    idx = idx + 1
                end do
            end do
        end do
        if (n_particle == 1) then
            eigval_pairs = 4.0_rp * eigval_pairs
        else
            eigval_pairs = 2.0_rp * eigval_pairs
        end if

    end function get_hess_eigval_pairs_oao

    subroutine get_extra_trial_vectors_oao(self, trial_vectors, settings, error)
        !
        ! this subroutine returns curvature-informed extra trial vectors for the
        ! solver's initial trial space: the orbital rotations between those
        ! occupied-virtual eigenvector pairs of the static Hessian part whose
        ! eigenvalue sums are most negative
        !
        class(oao_type), intent(inout) :: self
        real(rp), intent(out) :: trial_vectors(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, j, k, idx
        real(rp), allocatable :: eigval_pairs(:), dm_eigvecs(:, :)
        logical, allocatable :: is_occupied(:, :)
        real(rp), external :: ddot
        external :: dgemm

        ! initialize error flag
        error = 0

        ! number of AOs and particles
        n_ao = self%n_ao
        n_particle = self%n_particle

        ! refresh the eigendecomposition if the static Hessian part has changed
        call self%refresh_hess_eigen(settings, error)
        if (error /= 0) return

        ! determine which eigenvectors span the occupied space; the eigenspaces of the
        ! static Hessian part coincide with the occupied and virtual subspaces, so the
        ! density matrix expectation value of each eigenvector is either one or zero
        allocate(is_occupied(n_ao, n_particle), dm_eigvecs(n_ao, n_ao))
        do k = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, self%dm_oao(:, :, k), n_ao, &
                       self%hess_eigvecs(:, :, k), n_ao, 0.0_rp, dm_eigvecs, n_ao)
            do i = 1, n_ao
                is_occupied(i, k) = ddot(n_ao, self%hess_eigvecs(:, i, k), 1_ip, &
                                         dm_eigvecs(:, i), 1_ip) > 0.5_rp
            end do
        end do
        deallocate(dm_eigvecs)

        ! get pairwise sums of eigenvalues and exclude the redundant pairs, which are
        ! the ones whose two eigenvectors lie in the same subspace
        eigval_pairs = self%get_hess_eigval_pairs()
        idx = 1
        do k = 1, n_particle
            do j = 1, n_ao
                do i = 1, j - 1
                    if (is_occupied(i, k) .eqv. is_occupied(j, k)) &
                        eigval_pairs(idx) = huge(1.0_rp)
                    idx = idx + 1
                end do
            end do
        end do
        deallocate(is_occupied)

        ! fill the requested slots with the rotations belonging to the most negative
        ! remaining eigenvalue pairs
        call self%fill_extra_trial_vectors(eigval_pairs, trial_vectors)

    end subroutine get_extra_trial_vectors_oao

    subroutine rotate_dm_ao(kappa, dm_oao, s_inv_sqrt, rot_dm_ao, settings, error, &
                            rot_dm_oao)
        !
        ! this subroutine returns the density matrix in the OAO basis rotated by the
        ! orbital rotation kappa, transformed to the AO basis by the inverse square
        ! root of the AO overlap matrix
        !
        use otr_common, only: matrix_exponential

        real(rp), intent(in) :: kappa(:), dm_oao(:, :, :), s_inv_sqrt(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        real(rp), intent(out) :: rot_dm_ao(:, :, :)
        real(rp), intent(out), target, optional :: rot_dm_oao(:, :, :)

        real(rp), allocatable :: kappa_full(:, :, :), temp(:, :)
        real(rp), pointer :: rot_dm_oao_ptr(:, :, :)
        real(rp), allocatable, target :: u(:, :, :), rot_dm_oao_local(:, :, :)
        integer(ip) :: n_ao, n_particle, i

        ! initialize error flag
        error = 0

        ! number of AOs and particles
        n_ao = size(dm_oao, 1, kind=ip)
        n_particle = size(dm_oao, 3, kind=ip)

        ! get rotation matrix
        allocate(u(n_ao, n_ao, n_particle))
        kappa_full = unpack_asymm(kappa, n_particle, n_ao)
        do i = 1, n_particle
            u(:, :, i) = matrix_exponential(kappa_full(:, :, i), settings, error)
            if (error /= 0) return
        end do

        ! prepare rotated density matrix array
        if (present(rot_dm_oao)) then
            rot_dm_oao_ptr => rot_dm_oao
        else
            allocate(rot_dm_oao_local(n_ao, n_ao, n_particle))
            rot_dm_oao_ptr => rot_dm_oao_local
        end if

        ! rotate density matrix
        allocate(temp(n_ao, n_ao))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, i), n_ao, &
                       u(:, :, i), n_ao, 0.0_rp, temp, n_ao)
            call dgemm("T", "N", n_ao, n_ao, n_ao, 1.0_rp, u(:, :, i), n_ao, temp, &
                       n_ao, 0.0_rp, rot_dm_oao_ptr(:, :, i), n_ao)
        end do

        ! purify density matrix
        call purify(rot_dm_oao_ptr)

        ! transform density matrix from OAO basis to AO basis
        rot_dm_ao = symmetric_transformation(s_inv_sqrt, rot_dm_oao_ptr)

        ! deallocate local memory if needed
        if (.not. present(rot_dm_oao)) deallocate(rot_dm_oao_local)

    end subroutine rotate_dm_ao

    function hess_x_static_oao(x_full, fock_oo, fock_vv) result(hess_x_full)
        !
        ! this function applies the static part of the Hessian in the OAO basis (which
        ! is built from the occupied-occupied and virtual-virtual parts of the Fock
        ! matrix) to the unpacked trial vector of every particle channel, leaving the
        ! projection onto the occupied-virtual and virtual-occupied subspace, the
        ! scaling and the packing to the caller, so that the exact Hessian linear
        ! transformation can project the static and response parts together
        !
        real(rp), intent(in) :: x_full(:, :, :), fock_oo(:, :, :), fock_vv(:, :, :)
        real(rp), allocatable :: hess_x_full(:, :, :)

        integer(ip) :: n_ao, i
        external :: dgemm

        ! number of AOs
        n_ao = size(x_full, 1, kind=ip)

        ! apply the static part (F_vv - F_oo) X - h.c. to every particle channel
        allocate(hess_x_full, mold=x_full)
        do i = 1, size(x_full, 3, kind=ip)
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, &
                       fock_vv(:, :, i) - fock_oo(:, :, i), n_ao, x_full(:, :, i), &
                       n_ao, 0.0_rp, hess_x_full(:, :, i), n_ao)
            hess_x_full(:, :, i) = hess_x_full(:, :, i) - &
                                   transpose(hess_x_full(:, :, i))
        end do

    end function hess_x_static_oao

    function project_asymm(matrix, dm_oao) result(projected_matrix)
        !
        ! this function projects a matrix onto the occupied-virtual and
        ! virtual-occupied subspace of the provided density matrix and returns the
        ! result in antisymmetric form, discarding the redundant occupied-occupied and
        ! virtual-virtual contributions; for a symmetric matrix (e.g. a Fock response)
        ! this antisymmetrizes the retained occupied-virtual block, while for an
        ! already antisymmetric matrix (e.g. an orbital rotation) it reproduces its
        ! occupied-virtual and virtual-occupied contributions unchanged
        !
        real(rp), intent(in) :: matrix(:, :, :), dm_oao(:, :, :)
        real(rp), allocatable :: projected_matrix(:, :, :)

        integer(ip) :: n_ao, i, j
        real(rp), allocatable :: proj_v(:, :), temp(:, :)
        external :: dgemm

        ! number of AOs
        n_ao = size(matrix, 1)

        allocate(projected_matrix(n_ao, n_ao, size(matrix, 3)), proj_v(n_ao, n_ao), &
                 temp(n_ao, n_ao))
        do i = 1, size(matrix, 3)
            ! construct projection matrix on virtual space (I-D)
            proj_v = 0.0_rp
            do j = 1, n_ao
                proj_v(j, j) = 1.0_rp
            end do
            proj_v = proj_v - dm_oao(:, :, i)

            ! construct occupied-virtual contributions DM(I-D)
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, matrix(:, :, i), n_ao, &
                       proj_v, n_ao, 0.0_rp, temp, n_ao)
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, i), n_ao, &
                       temp, n_ao, 0.0_rp, projected_matrix(:, :, i), n_ao)

            ! antisymmetrize to add the virtual-occupied contributions
            projected_matrix(:, :, i) = projected_matrix(:, :, i) - &
                                        transpose(projected_matrix(:, :, i))
        end do
        deallocate(proj_v, temp)

    end function project_asymm

    function project_symm(x_full, dm_oao) result(projected_matrix)
        !
        ! this function projects an antisymmetric orbital-rotation matrix onto the
        ! occupied-virtual and virtual-occupied subspace of the provided density matrix
        ! and returns the result in symmetric form, i.e. the density matrix response to
        ! the rotation, discarding the redundant occupied-occupied and virtual-virtual
        ! contributions
        !
        real(rp), intent(in) :: x_full(:, :, :), dm_oao(:, :, :)
        real(rp), allocatable :: projected_matrix(:, :, :)

        integer(ip) :: n_ao, i
        real(rp), allocatable :: temp(:, :)
        external :: dgemm

        ! number of AOs
        n_ao = size(x_full, 1)

        allocate(projected_matrix(n_ao, n_ao, size(x_full, 3)), temp(n_ao, n_ao))
        do i = 1, size(x_full, 3)
            ! construct product of density matrix and trial vector D*X
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, i), n_ao, &
                       x_full(:, :, i), n_ao, 0.0_rp, temp, n_ao)

            ! symmetrize to get density matrix response
            projected_matrix(:, :, i) = temp + transpose(temp)
        end do
        deallocate(temp)

    end function project_symm

    subroutine purify(dm)
        !
        ! this function purifies a density matrix (dm_purified = 3*dm^2 - 2*dm^3)
        !
        real(rp), intent(inout) :: dm(:, :, :)

        real(rp), allocatable :: dm_squared(:, :), dm_cubed(:, :)
        integer(ip) :: n, i
        external :: dgemm

        ! size
        n = size(dm, 1)

        ! allocate arrays
        allocate(dm_squared(n, n), dm_cubed(n, n))

        do i = 1, size(dm, 3)
            ! square density matrix
            call dgemm("N", "N", n, n, n, 1.0_rp, dm(:, :, i), n, dm(:, :, i), n, &
                       0.0_rp, dm_squared, n)

            ! cube density matrix
            call dgemm("N", "N", n, n, n, 1.0_rp, dm_squared, n, dm(:, :, i), n, &
                       0.0_rp, dm_cubed, n)

            ! purify density matrix
            dm(:, :, i) = 3.0_rp * dm_squared - 2.0_rp * dm_cubed
        end do
        deallocate(dm_squared, dm_cubed)

    end subroutine purify

    function symmetric_transformation(trans_matrix, matrix) result(matrix_transformed)
        !
        ! this function performs a symmetric transformation (X' = U * X * U)
        !
        real(rp), intent(in) :: trans_matrix(:, :), matrix(:, :, :)
        real(rp), allocatable :: matrix_transformed(:, :, :)

        real(rp), allocatable :: temp(:, :)
        integer(ip) :: n, i
        external :: dgemm

        n = size(matrix, 1)

        allocate(temp(n, n), matrix_transformed(n, n, size(matrix, 3)))
        do i = 1, size(matrix, 3)
            call dgemm("N", "N", n, n, n, 1.0_rp, trans_matrix, n, matrix(:, :, i), n, &
                       0.0_rp, temp, n)
            call dgemm("N", "N", n, n, n, 1.0_rp, temp, n, trans_matrix, n, 0.0_rp, &
                       matrix_transformed(:, :, i), n)
        end do
        deallocate(temp)

    end function symmetric_transformation

    function unpack_asymm(matrix_nonred, n_particle, n_ao) result(matrix)
        !
        ! this function unpacks an antisymmetric matrix and returns the resulting
        ! unpacked matrix
        !
        real(rp), intent(in) :: matrix_nonred(:)
        integer(ip), intent(in) :: n_particle, n_ao
        real(rp), allocatable :: matrix(:, :, :)

        integer(ip) :: i, j, k, idx

        ! initialize full matrix
        allocate(matrix(n_ao, n_ao, n_particle))
        matrix = 0.0_rp

        ! unpack asymmetric matrix
        idx = 1
        do k = 1, n_particle
            do j = 1, n_ao
                do i = 1, j - 1
                    matrix(i, j, k) = matrix_nonred(idx)
                    matrix(j, i, k) = -matrix_nonred(idx)
                    idx = idx + 1
                end do
            end do
        end do

    end function unpack_asymm

    function pack_asymm(matrix, n_param) result(matrix_nonred)
        !
        ! this function packs an antisymmetric matrix for RHF while deallocating the
        ! original unpacked matrix
        !
        real(rp), intent(in) :: matrix(:, :, :)
        integer(ip), intent(in) :: n_param
        real(rp), allocatable :: matrix_nonred(:)

        integer(ip) :: idx, i, j, k

        ! allocate redundant matrix
        allocate(matrix_nonred(n_param))

        ! pack asymmetric matrix
        idx = 1
        do k = 1, size(matrix, 3)
            do j = 1, size(matrix, 2)
                do i = 1, j - 1
                    matrix_nonred(idx) = matrix(i, j, k)
                    idx = idx + 1
                end do
            end do
        end do

    end function pack_asymm

end module otr_oao
