! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo

    use opentrustregion_core, only: rp, ip, settings_type, obj_func_type, &
                                    update_orbs_type, hess_x_type, solver_state_type, &
                                    add_error_origin
    use otr_common, only: orbital_settings_type, orbital_basis_type, channel_rows, &
                          evaluate_dm_cs_type, evaluate_dm_os_type, &
                          extension_solver_settings_type

    implicit none

    type, extends(orbital_settings_type) :: mo_settings_type
    end type

    type(mo_settings_type), parameter :: default_mo_settings = mo_settings_type()

    ! warning for settings that were not initialized
    character(len=*), parameter :: mo_settings_uninitialized_warning_msg = &
        "MO settings were not initialized. All MO settings are set to default values"

    ! error origin of the MO extension, with which its factory and object tag the
    ! errors they detect themselves
    integer(ip), parameter :: error_mo = 5100

    ! occupations and orbital irreps of a particle channel in the MO basis, the mask of
    ! its occupied-virtual pairs of the same irrep, which are its parameters, and the
    ! occupied-occupied and virtual-virtual blocks of its Fock matrix, from which the
    ! static part of the Hessian is built, together with their eigendecompositions
    type :: mo_channel_type
        integer(ip) :: n_occ = 0, n_virt = 0
        integer(ip), allocatable :: irreps(:)
        logical, allocatable :: param_mask(:, :)
        real(rp), allocatable :: fock_oo(:, :), fock_vv(:, :), occ_eigvecs(:, :), &
                                 occ_eigvals(:), virt_eigvecs(:, :), virt_eigvals(:)
    end type

    ! orbitals parameterized in the MO basis, which holds its density matrix itself, at
    ! which the pointer of the orbital basis points
    type, extends(orbital_basis_type) :: mo_type
        integer(ip) :: n_mo
        real(rp), pointer, contiguous :: mo_coeff(:, :, :) => null()
        real(rp), allocatable :: ao_overlap(:, :), own_dm_ao(:, :, :)
        type(mo_channel_type), allocatable :: mo_channels(:)
    contains
        procedure :: update_orbs => update_orbs_mo
        procedure :: obj_func => obj_func_mo
        procedure :: hess_x => hess_x_mo
        procedure :: get_extra_trial_vectors => get_extra_trial_vectors_mo
        procedure :: rotate_orbitals => rotate_orbitals_mo
        procedure :: calculate_grad_h_diag => calculate_grad_h_diag_mo
        procedure :: refresh_hess_eigen => refresh_hess_eigen_mo
        procedure :: rotate_to_hess_eigenbasis => rotate_to_hess_eigenbasis_mo
        procedure :: rotate_from_hess_eigenbasis => rotate_from_hess_eigenbasis_mo
        procedure :: get_hess_eigval_pairs => get_hess_eigval_pairs_mo
    end type

    ! define module procedures for different spin cases
    interface mo_factory
        module procedure mo_factory_cs, mo_factory_os
    end interface

contains

    subroutine mo_factory_cs(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                             evaluate_dm_cs, update_orbs_mo_funptr, &
                             obj_func_mo_funptr, solver_settings, solver_state, error, &
                             settings, orbsym)
        !
        ! this function returns the MO orbital updating function and objective function
        ! for the closed-shell case and wires the solver settings with a new MO object,
        ! which they own, whose preconditioners and extra trial vectors the solver then
        ! uses, while the projection is left to the caller, since the parameters are
        ! non-redundant, and which reads the state of the solver from the record the
        ! host also passes to the solver; the MO coefficients are rotated in place, so
        ! they have to outlive the calculation; if the irreps of the MOs are given,
        ! only the occupied-virtual rotations between orbitals of the same irrep are
        ! parameters; without a density matrix evaluating function, the orbital
        ! updating function, objective function and Hessian linear transformation of
        ! the object return an error when called, so that only the density QN factory
        ! can build on it
        !
        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ, n_particle, n_ao, n_mo
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_cs
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_mo_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(mo_settings_type), intent(in), optional :: settings
        integer(ip), intent(in), optional :: orbsym(:)

        real(rp), pointer, contiguous :: mo_coeff_3d(:, :, :)
        integer(ip), allocatable :: orbsym_2d(:, :)
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os

        ! call common setup with the closed-shell density matrix evaluating function,
        ! passing the irreps of the MOs only if they are given
        mo_coeff_3d(1:size(mo_coeff, 1), 1:size(mo_coeff, 2), 1:1) => mo_coeff
        if (present(orbsym)) orbsym_2d = reshape(orbsym, [size(orbsym, kind=ip), 1_ip])
        evaluate_dm_os => null()
        call mo_factory_common( &
            mo_coeff_3d, ao_overlap, [n_occ], n_particle, n_ao, n_mo, evaluate_dm_cs, &
            evaluate_dm_os, update_orbs_mo_funptr, obj_func_mo_funptr, &
            solver_settings, solver_state, error, settings, orbsym_2d)

    end subroutine mo_factory_cs

    subroutine mo_factory_os(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                             evaluate_dm_os, update_orbs_mo_funptr, &
                             obj_func_mo_funptr, solver_settings, solver_state, error, &
                             settings, orbsym)
        !
        ! this function returns the MO orbital updating function and objective function
        ! for the open-shell case and wires the solver settings with a new MO object,
        ! which they own, whose preconditioners and extra trial vectors the solver then
        ! uses, while the projection is left to the caller, since the parameters are
        ! non-redundant, and which reads the state of the solver from the record the
        ! host also passes to the solver; the MO coefficients are rotated in place, so
        ! they have to outlive the calculation; if the irreps of the MOs of every
        ! particle channel are given, only the occupied-virtual rotations between
        ! orbitals of the same irrep are parameters; without a density matrix
        ! evaluating function, the orbital updating function, objective function and
        ! Hessian linear transformation of the object return an error when called, so
        ! that only the density QN factory can build on it
        !
        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_mo_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(mo_settings_type), intent(in), optional :: settings
        integer(ip), intent(in), optional :: orbsym(:, :)

        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs

        ! call common setup with the open-shell density matrix evaluating function
        evaluate_dm_cs => null()
        call mo_factory_common(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                               evaluate_dm_cs, evaluate_dm_os, update_orbs_mo_funptr, &
                               obj_func_mo_funptr, solver_settings, solver_state, &
                               error, settings, orbsym)

    end subroutine mo_factory_os

    subroutine mo_factory_common(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                                 evaluate_dm_cs, evaluate_dm_os, &
                                 update_orbs_mo_funptr, obj_func_mo_funptr, &
                                 solver_settings, solver_state, error, settings, orbsym)
        !
        ! this subroutine sets up a new MO object with the density matrix evaluating
        ! function of either shell, if one is given, wires the solver settings with it
        ! and returns its orbital updating function and objective function
        !
        use opentrustregion_core, only: verbosity_warning
        use otr_common, only: extension_type, provided_roles_type, &
                              update_orbs_callback, obj_func_callback, &
                              evaluate_dm_missing_warning_msg

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_cs
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_mo_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(mo_settings_type), intent(in), optional :: settings
        integer(ip), intent(in), optional :: orbsym(:, :)

        type(mo_settings_type) :: local_settings
        class(extension_type), allocatable, target :: object
        integer(ip) :: i, j
        integer(ip), allocatable :: irreps(:, :)
        external :: dgemm

        ! initialize error flag and the returned functions
        error = 0
        update_orbs_mo_funptr => null()
        obj_func_mo_funptr => null()

        ! work on a copy of the given settings or, without any, of the default settings
        if (present(settings)) then
            local_settings = settings
        else
            local_settings = default_mo_settings
        end if

        ! perform sanity check
        call mo_sanity_check(local_settings, shape(mo_coeff, kind=ip), shape( &
            ao_overlap, kind=ip), n_occ, n_particle, n_ao, n_mo, error, orbsym)
        call add_error_origin(error, error_mo, local_settings)
        if (error /= 0) return

        ! without a density matrix evaluating function, the roles which need one are
        ! not available
        if (.not. (associated(evaluate_dm_cs) .or. associated(evaluate_dm_os))) &
            call local_settings%log(evaluate_dm_missing_warning_msg, verbosity_warning)

        ! irreps of the MOs of every particle channel
        if (present(orbsym)) then
            irreps = orbsym
        else
            allocate(irreps(n_mo, n_particle), source=0_ip)
        end if

        ! set up the object, which performs every role but the projection itself, in
        ! place, since its density matrix points at the one it holds, for the MO
        ! coefficients of the caller, which it points to
        allocate(mo_type :: object)
        select type (object)
        type is (mo_type)
            object%provides = provided_roles_type( &
                update_orbs=.true., obj_func=.true., hess_x=.true., precond=.true., &
                precond_pd=.true., get_extra_trial_vectors=.true.)
            object%error_origin = error_mo
            object%settings = local_settings%orbital_settings_type
            object%evaluate_dm_cs => evaluate_dm_cs
            object%evaluate_dm_os => evaluate_dm_os

            ! number of atomic orbitals, MOs and particles
            object%n_ao = n_ao
            object%n_mo = n_mo
            object%n_particle = n_particle

            ! set occupations and irreps of every particle channel, whose parameters
            ! are the occupied-virtual pairs of the same irrep
            allocate(object%mo_channels(n_particle))
            object%mo_channels%n_occ = n_occ
            object%mo_channels%n_virt = n_mo - n_occ
            do i = 1, n_particle
                associate (channel => object%mo_channels(i))
                    channel%irreps = irreps(:, i)
                    allocate(channel%param_mask(n_occ(i), n_mo - n_occ(i)))
                    do j = 1, n_mo - n_occ(i)
                        channel%param_mask(:, j) = &
                            irreps(:n_occ(i), i) == irreps(n_occ(i) + j, i)
                    end do
                end associate
            end do

            ! number of parameters is the number of these pairs of every channel
            object%n_param = count_mo_params(n_occ, n_mo, orbsym)

            ! starting MO coefficients, which are rotated in place, and AO overlap
            ! matrix
            object%mo_coeff => mo_coeff
            object%ao_overlap = ao_overlap

            ! construct the starting density matrix in the AO basis from the occupied
            ! orbitals in the one the object holds
            allocate(object%own_dm_ao(n_ao, n_ao, n_particle))
            object%dm_ao => object%own_dm_ao
            do i = 1, n_particle
                call dgemm("N", "T", n_ao, n_ao, n_occ(i), 1.0_rp, mo_coeff(:, :, i), &
                           n_ao, mo_coeff(:, :, i), n_ao, 0.0_rp, &
                           object%dm_ao(:, :, i), n_ao)
            end do

            ! allocate gradient and Hessian diagonal
            allocate(object%grad(object%n_param), object%h_diag(object%n_param))
        end select

        ! wire the solver settings with the object, which reads the record of the state
        ! and saves the host's callback functions and host contexts, and let the
        ! settings own it
        call solver_settings%wire(object, solver_state, local_settings, error)
        call add_error_origin(error, error_mo, local_settings)
        if (error /= 0) return

        ! get the orbital updating function and the objective function of the object
        update_orbs_mo_funptr => update_orbs_callback
        obj_func_mo_funptr => obj_func_callback

    end subroutine mo_factory_common

    subroutine mo_sanity_check(settings, mo_coeff_shape, ao_overlap_shape, n_occ, &
                               n_particle, n_ao, n_mo, error, orbsym)
        !
        ! this subroutine performs a sanity check for the dimensions of orbitals
        ! parameterized in the MO basis, given the shapes of the MO coefficients (with
        ! the particle channels along the last dimension) and of the AO overlap matrix,
        ! the number of particle channels of the factory and, if given, the irreps of
        ! the MOs of every particle channel
        !
        use opentrustregion_core, only: verbosity_error

        class(orbital_settings_type), intent(in) :: settings
        integer(ip), intent(in) :: mo_coeff_shape(3), ao_overlap_shape(2), n_occ(:), &
                                   n_particle, n_ao, n_mo
        integer(ip), intent(out) :: error
        integer(ip), intent(in), optional :: orbsym(:, :)

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

        ! check that there is one particle channel for the closed-shell and two for the
        ! open-shell case
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

        ! check that the irreps are given for every MO of every particle channel
        if (present(orbsym)) then
            if (any(shape(orbsym, kind=ip) /= [n_mo, n_particle])) then
                call settings%log("Shape of the MO irreps should match the number "// &
                                  "of MOs and particle channels.", verbosity_error, &
                                  .true.)
                error = 1
                return
            end if
        end if

        ! check that the occupations fit into the MOs
        if (any(n_occ < 0) .or. any(n_occ > n_mo)) then
            call settings%log("Number of occupied orbitals should not be negative "// &
                              "and not larger than the number of MOs.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that there is at least one occupied-virtual rotation, between orbitals
        ! of the same irrep if the irreps are given
        if (count_mo_params(n_occ, n_mo, orbsym) < 1) then
            if (present(orbsym)) then
                call settings%log("There should be at least one occupied-virtual "// &
                                  "rotation between orbitals of the same irrep.", &
                                  verbosity_error, .true.)
            else
                call settings%log("There should be at least one occupied-virtual "// &
                                  "rotation.", verbosity_error, .true.)
            end if
            error = 1
            return
        end if

    end subroutine mo_sanity_check

    function obj_func_mo(self, kappa, error, context) result(energy)
        !
        ! this function defines the energy evaluation in the MO basis
        !
        class(mo_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: energy

        real(rp), allocatable :: rot_mo_coeff(:, :, :), rot_dm_ao(:, :, :)

        ! initialize energy in case of error
        energy = 0.0_rp

        ! the energy needs the density matrix evaluating function
        call self%check_evaluate_dm(error)
        if (error /= 0) return

        ! get rotated density matrix in AO basis without moving the current orbitals
        allocate(rot_mo_coeff, mold=self%mo_coeff)
        allocate(rot_dm_ao, mold=self%dm_ao)
        call rotate_mo_coeff(kappa, self%mo_coeff, self%ao_overlap, self%mo_channels, &
                             rot_mo_coeff, rot_dm_ao, self%settings, error)
        call add_error_origin(error, self%error_origin, self%settings)
        if (error /= 0) return

        ! calculate the mean-field energy
        if (associated(self%evaluate_dm_cs)) then
            call self%evaluate_dm_cs(rot_dm_ao(:, :, 1), energy, error=error, &
                                     context=self%settings%context)
        else
            call self%evaluate_dm_os(rot_dm_ao, energy, error=error, &
                                     context=self%settings%context)
        end if
        if (error /= 0) return

    end function obj_func_mo

    subroutine update_orbs_mo(self, kappa, func, grad, h_diag, hess_x_funptr, error, &
                              context)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation
        ! and the Hessian linear transformation in the MO basis
        !
        use otr_common, only: hess_x_callback

        class(mo_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: n_ao, n_particle
        real(rp), allocatable :: fock_ao(:, :, :)

        ! initialize error flag
        error = 0

        ! the evaluation needs the density matrix evaluating function
        call self%check_evaluate_dm(error)
        if (error /= 0) return

        ! evaluate at the rotated orbitals, or at the current ones if these have not
        ! been evaluated yet or their response is stale
        if ((sum(abs(kappa)) > 0.0_rp) .or. self%evaluation_stale .or. &
            self%response_stale) then
            ! number of AOs
            n_ao = self%n_ao

            ! number of particles
            n_particle = self%n_particle

            ! rotate orbitals, which have not been evaluated until this succeeds
            self%evaluation_stale = .true.
            call self%rotate_orbitals(kappa, error)
            if (error /= 0) return

            ! get energy, Fock matrix, and response function
            allocate(fock_ao(n_ao, n_ao, n_particle))
            if (associated(self%evaluate_dm_cs)) then
                call self%evaluate_dm_cs(self%dm_ao(:, :, 1), self%energy, &
                                         fock_ao(:, :, 1), self%get_response_cs, &
                                         error, self%settings%context)
            else
                call self%evaluate_dm_os(self%dm_ao, self%energy, fock_ao, &
                                         self%get_response_os, error, &
                                         self%settings%context)
            end if
            if (error /= 0) then
                deallocate(fock_ao)
                return
            end if

            ! calculate gradient, Hessian diagonal and static part of the Hessian from
            ! the Fock matrix in the AO basis
            call self%calculate_grad_h_diag(fock_ao)
            deallocate(fock_ao)

            ! the rotated orbitals and their response have been evaluated
            self%evaluation_stale = .false.
            self%response_stale = .false.
        end if

        ! set outputs
        func = self%energy
        grad = self%grad
        h_diag = self%h_diag
        hess_x_funptr => hess_x_callback

    end subroutine update_orbs_mo

    subroutine hess_x_mo(self, x, hess_x, error, context)
        !
        ! this function defines the Hessian linear transformation in the MO basis, the
        ! static part from the occupied-occupied and virtual-virtual blocks of the Fock
        ! matrix together with the response of the Fock matrix to the density matrix
        ! displacement C_o X C_v^T + C_v X^T C_o^T of the occupied-virtual block X of
        ! every particle channel, transformed back to the occupied-virtual block
        !
        class(mo_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: n_ao, n_occ, n_virt, i, rows(2, self%n_particle)
        real(rp) :: shell_scale
        real(rp), allocatable :: x_block(:, :), temp(:, :), dm_response(:, :, :), &
                                 fock_response(:, :, :), hess_x_block(:, :)
        external :: dgemm

        ! initialize error flag
        error = 0

        ! the response part needs the density matrix evaluating function
        call self%check_evaluate_dm(error)
        if (error /= 0) return

        ! number of AOs
        n_ao = self%n_ao

        ! set scaling factor for closed- and open-shell systems
        shell_scale = merge(4.0_rp, 2.0_rp, self%n_particle == 1)

        ! rows of every particle channel in the parameter vector
        rows = mo_param_rows(self%mo_channels)

        ! rebuild the response if the orbitals were moved without it being updated,
        ! since the static and response parts would otherwise refer to different points
        if (self%response_stale) then
            call self%refresh_response(error)
            if (error /= 0) return
        end if

        ! get static part
        hess_x = hess_x_static_mo(x, self%mo_channels)

        ! get density matrix response to trial vector in the AO basis for every
        ! particle channel, before the response of the Fock matrix couples them
        allocate(dm_response(n_ao, n_ao, self%n_particle))
        dm_response = 0.0_rp
        do i = 1, self%n_particle
            n_occ = self%mo_channels(i)%n_occ
            n_virt = self%mo_channels(i)%n_virt
            if (n_occ == 0 .or. n_virt == 0) cycle
            x_block = unpack_ov(x(rows(1, i):rows(2, i)), self%mo_channels(i))
            allocate(temp(n_occ, n_ao))
            call dgemm("N", "T", n_occ, n_ao, n_virt, 1.0_rp, x_block, n_occ, &
                       self%mo_coeff(:, n_occ + 1:, i), n_ao, 0.0_rp, temp, n_occ)
            call dgemm("N", "N", n_ao, n_ao, n_occ, 1.0_rp, &
                       self%mo_coeff(:, :n_occ, i), n_ao, temp, n_occ, 0.0_rp, &
                       dm_response(:, :, i), n_ao)
            dm_response(:, :, i) = dm_response(:, :, i) + &
                                   transpose(dm_response(:, :, i))
            deallocate(x_block, temp)
        end do

        ! get response of Fock matrix to density matrix response
        allocate(fock_response(n_ao, n_ao, self%n_particle))
        if (associated(self%get_response_os)) then
            call self%get_response_os(dm_response, fock_response, error, &
                                      self%settings%context)
        else
            call self%get_response_cs(dm_response(:, :, 1), fock_response(:, :, 1), &
                                      error, self%settings%context)
        end if
        deallocate(dm_response)
        if (error /= 0) return

        ! add the occupied-virtual block of the Fock matrix response in the MO basis
        do i = 1, self%n_particle
            n_occ = self%mo_channels(i)%n_occ
            n_virt = self%mo_channels(i)%n_virt
            if (n_occ == 0 .or. n_virt == 0) cycle
            allocate(temp(n_occ, n_ao), hess_x_block(n_occ, n_virt))
            call dgemm("T", "N", n_occ, n_ao, n_ao, 1.0_rp, &
                       self%mo_coeff(:, :n_occ, i), n_ao, fock_response(:, :, i), &
                       n_ao, 0.0_rp, temp, n_occ)
            call dgemm("N", "N", n_occ, n_virt, n_ao, 1.0_rp, temp, n_occ, &
                       self%mo_coeff(:, n_occ + 1:, i), n_ao, 0.0_rp, hess_x_block, &
                       n_occ)
            hess_x(rows(1, i):rows(2, i)) = &
                hess_x(rows(1, i):rows(2, i)) + &
                shell_scale * pack_ov(hess_x_block, self%mo_channels(i))
            deallocate(temp, hess_x_block)
        end do
        deallocate(fock_response)

    end subroutine hess_x_mo

    subroutine rotate_orbitals_mo(self, kappa, error)
        !
        ! this subroutine moves the current orbitals by the orbital rotation kappa,
        ! updating the MO coefficients in place together with the density matrix in the
        ! AO basis
        !
        class(mo_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: rot_mo_coeff(:, :, :), rot_dm_ao(:, :, :)

        ! rotate MO coefficients and update the density matrix accordingly
        allocate(rot_mo_coeff, mold=self%mo_coeff)
        allocate(rot_dm_ao, mold=self%dm_ao)
        call rotate_mo_coeff(kappa, self%mo_coeff, self%ao_overlap, self%mo_channels, &
                             rot_mo_coeff, rot_dm_ao, self%settings, error)
        call add_error_origin(error, self%error_origin, self%settings)
        if (error /= 0) return
        self%mo_coeff = rot_mo_coeff
        self%dm_ao = rot_dm_ao

        ! the orbitals were moved but the response was not rebuilt
        self%response_stale = .true.

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

        integer(ip) :: n_occ, n_virt, i, j, k, rows(2, self%n_particle)
        real(rp) :: shell_scale
        real(rp), allocatable :: fock_mo(:, :), diag_diff(:, :)

        ! set scaling factor for closed- and open-shell systems
        shell_scale = merge(4.0_rp, 2.0_rp, self%n_particle == 1)

        ! rows of every particle channel in the parameter vector
        rows = mo_param_rows(self%mo_channels)

        do k = 1, self%n_particle
            n_occ = self%mo_channels(k)%n_occ
            n_virt = self%mo_channels(k)%n_virt

            ! transform Fock matrix to the MO basis
            fock_mo = mo_transform(self%mo_coeff(:, :, k), fock(:, :, k))

            ! occupied-occupied and virtual-virtual blocks of the Fock matrix
            self%mo_channels(k)%fock_oo = fock_mo(:n_occ, :n_occ)
            self%mo_channels(k)%fock_vv = fock_mo(n_occ + 1:, n_occ + 1:)

            ! construct gradient from the occupied-virtual block
            self%grad(rows(1, k):rows(2, k)) = &
                shell_scale * pack_ov(fock_mo(:n_occ, n_occ + 1:), self%mo_channels(k))

            ! construct Hessian diagonal from the differences of the diagonal elements
            ! of the virtual-virtual and occupied-occupied blocks
            allocate(diag_diff(n_occ, n_virt))
            do j = 1, n_virt
                do i = 1, n_occ
                    diag_diff(i, j) = fock_mo(n_occ + j, n_occ + j) - fock_mo(i, i)
                end do
            end do
            self%h_diag(rows(1, k):rows(2, k)) = shell_scale * &
                                                 pack_ov(diag_diff, self%mo_channels(k))
            deallocate(diag_diff)
        end do

        ! the static Hessian part was just rebuilt, so any cached eigendecomposition of
        ! it is now stale
        self%hess_eigen_stale = .true.

    end subroutine calculate_grad_h_diag_mo

    subroutine refresh_hess_eigen_mo(self, error)
        !
        ! this subroutine refreshes the eigendecomposition of the static part of the
        ! Hessian in the MO basis, if the static part has changed since it was last
        ! computed, by diagonalizing, per irrep, the occupied-occupied and
        ! virtual-virtual blocks of the Fock matrix of every particle channel, whose
        ! eigenvectors diagonalize the static part
        !
        class(mo_type), intent(inout) :: self
        integer(ip), intent(out) :: error

        integer(ip) :: n_occ, n_virt, i

        ! initialize error flag
        error = 0

        ! nothing to do if the cached eigendecomposition is up to date
        if (.not. self%hess_eigen_stale) return

        do i = 1, self%n_particle
            associate (channel => self%mo_channels(i))
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
                call diagonalize_per_irrep(channel%fock_oo, channel%irreps(:n_occ), &
                                           channel%occ_eigvals, channel%occ_eigvecs, &
                                           self%settings, error)
                call add_error_origin(error, self%error_origin, self%settings)
                if (error /= 0) return
                call diagonalize_per_irrep( &
                    channel%fock_vv, channel%irreps(n_occ + 1:), channel%virt_eigvals, &
                    channel%virt_eigvecs, self%settings, error)
                call add_error_origin(error, self%error_origin, self%settings)
                if (error /= 0) return
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
        rows = mo_param_rows(self%mo_channels)

        ! rotate each particle channel by the cached eigenvectors
        allocate(rotated(self%n_param))
        do i = 1, self%n_particle
            n_occ = self%mo_channels(i)%n_occ
            n_virt = self%mo_channels(i)%n_virt
            if (n_occ == 0 .or. n_virt == 0) cycle
            ov_block = unpack_ov(vector(rows(1, i):rows(2, i)), self%mo_channels(i))
            allocate(temp(n_occ, n_virt))
            call dgemm("T", "N", n_occ, n_virt, n_occ, 1.0_rp, &
                       self%mo_channels(i)%occ_eigvecs, n_occ, ov_block, n_occ, &
                       0.0_rp, temp, n_occ)
            call dgemm("N", "N", n_occ, n_virt, n_virt, 1.0_rp, temp, n_occ, &
                       self%mo_channels(i)%virt_eigvecs, n_virt, 0.0_rp, ov_block, &
                       n_occ)
            rotated(rows(1, i):rows(2, i)) = pack_ov(ov_block, self%mo_channels(i))
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
        rows = mo_param_rows(self%mo_channels)

        ! rotate each particle channel by the cached eigenvectors
        allocate(rotated(self%n_param))
        do i = 1, self%n_particle
            n_occ = self%mo_channels(i)%n_occ
            n_virt = self%mo_channels(i)%n_virt
            if (n_occ == 0 .or. n_virt == 0) cycle
            ov_block = unpack_ov(vector(rows(1, i):rows(2, i)), self%mo_channels(i))
            allocate(temp(n_occ, n_virt))
            call dgemm("N", "N", n_occ, n_virt, n_occ, 1.0_rp, &
                       self%mo_channels(i)%occ_eigvecs, n_occ, ov_block, n_occ, &
                       0.0_rp, temp, n_occ)
            call dgemm("N", "T", n_occ, n_virt, n_virt, 1.0_rp, temp, n_occ, &
                       self%mo_channels(i)%virt_eigvecs, n_virt, 0.0_rp, ov_block, &
                       n_occ)
            rotated(rows(1, i):rows(2, i)) = pack_ov(ov_block, self%mo_channels(i))
            deallocate(ov_block, temp)
        end do

    end function rotate_from_hess_eigenbasis_mo

    function get_hess_eigval_pairs_mo(self) result(eigval_pairs)
        !
        ! this function returns the eigenvalues of the cached static Hessian part in
        ! the MO basis, the scaled differences of virtual and occupied Fock matrix
        ! block eigenvalues of the same irrep, packed in the same order as h_diag
        !
        class(mo_type), intent(in) :: self
        real(rp), allocatable :: eigval_pairs(:)

        integer(ip) :: i, j, k, rows(2, self%n_particle)
        real(rp) :: shell_scale
        real(rp), allocatable :: eigval_diff(:, :)

        ! set scaling factor for closed- and open-shell systems
        shell_scale = merge(4.0_rp, 2.0_rp, self%n_particle == 1)

        ! rows of every particle channel in the parameter vector
        rows = mo_param_rows(self%mo_channels)

        ! construct differences of eigenvalues, each stored at the orbitals of the
        ! irrep of its eigenvector, so that packing keeps the pairs of the same irrep
        allocate(eigval_pairs(self%n_param))
        do k = 1, self%n_particle
            associate (channel => self%mo_channels(k))
                allocate(eigval_diff(channel%n_occ, channel%n_virt))
                do j = 1, channel%n_virt
                    do i = 1, channel%n_occ
                        eigval_diff(i, j) = channel%virt_eigvals(j) - &
                                            channel%occ_eigvals(i)
                    end do
                end do
                eigval_pairs(rows(1, k):rows(2, k)) = shell_scale * &
                                                      pack_ov(eigval_diff, channel)
                deallocate(eigval_diff)
            end associate
        end do

    end function get_hess_eigval_pairs_mo

    subroutine get_extra_trial_vectors_mo(self, trial_vectors, error, context)
        !
        ! this subroutine returns extra trial vectors for the MO basis along the
        ! eigenvectors of the static Hessian part with the most negative eigenvalues,
        ! the rotations between the pseudo-canonical occupied and virtual orbitals with
        ! the most negative orbital energy differences
        !
        class(mo_type), intent(inout) :: self
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! initialize error flag
        error = 0

        ! refresh the eigendecomposition if the static Hessian part has changed
        call self%refresh_hess_eigen(error)
        if (error /= 0) return

        ! fill the requested slots with the rotations belonging to the most negative
        ! eigenvalue pairs
        call self%fill_extra_trial_vectors(self%get_hess_eigval_pairs(), trial_vectors)

    end subroutine get_extra_trial_vectors_mo

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
        rows = mo_param_rows(mo_channels)

        allocate(kappa_full(n_mo, n_mo), rotated(n_ao, n_mo), &
                 overlap_rotated(n_ao, n_mo), metric(n_mo, n_mo))
        do i = 1, size(mo_channels, kind=ip)
            n_occ = mo_channels(i)%n_occ
            n_virt = mo_channels(i)%n_virt

            ! construct the antisymmetric rotation generator from its occupied-virtual
            ! block
            kappa_full = 0.0_rp
            kappa_full(:n_occ, n_occ + 1:) = unpack_ov(kappa(rows(1, i):rows(2, i)), &
                                                       mo_channels(i))
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

    function hess_x_static_mo(x, mo_channels) result(hess_x)
        !
        ! this function applies the static part of the Hessian in the MO basis (which
        ! is built from the occupied-occupied and virtual-virtual blocks of the Fock
        ! matrix of every particle channel) to a trial vector
        !
        real(rp), intent(in) :: x(:)
        type(mo_channel_type), intent(in) :: mo_channels(:)
        real(rp), allocatable :: hess_x(:)

        integer(ip) :: n_occ, n_virt, i, rows(2, size(mo_channels))
        real(rp) :: shell_scale
        real(rp), allocatable :: x_block(:, :), hess_x_block(:, :)
        external :: dgemm

        ! set scaling factor for closed- and open-shell systems
        shell_scale = merge(4.0_rp, 2.0_rp, size(mo_channels) == 1)

        ! rows of every particle channel in the parameter vector
        rows = mo_param_rows(mo_channels)

        ! apply the static part X F_vv - F_oo X to the occupied-virtual block of every
        ! particle channel
        allocate(hess_x(size(x)))
        do i = 1, size(mo_channels, kind=ip)
            n_occ = mo_channels(i)%n_occ
            n_virt = mo_channels(i)%n_virt
            if (n_occ == 0 .or. n_virt == 0) cycle
            x_block = unpack_ov(x(rows(1, i):rows(2, i)), mo_channels(i))
            allocate(hess_x_block(n_occ, n_virt))
            call dgemm("N", "N", n_occ, n_virt, n_virt, 1.0_rp, x_block, n_occ, &
                       mo_channels(i)%fock_vv, n_virt, 0.0_rp, hess_x_block, n_occ)
            call dgemm("N", "N", n_occ, n_virt, n_occ, -1.0_rp, &
                       mo_channels(i)%fock_oo, n_occ, x_block, n_occ, 1.0_rp, &
                       hess_x_block, n_occ)
            hess_x(rows(1, i):rows(2, i)) = shell_scale * &
                                            pack_ov(hess_x_block, mo_channels(i))
            deallocate(x_block, hess_x_block)
        end do

    end function hess_x_static_mo

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

    function count_mo_params(n_occ, n_mo, orbsym) result(n_param)
        !
        ! this function returns the number of parameters of orbitals parameterized in
        ! the MO basis, the occupied-virtual pairs of every particle channel, only
        ! those between orbitals of the same irrep if the irreps of the MOs are given
        !
        integer(ip), intent(in) :: n_occ(:), n_mo
        integer(ip), intent(in), optional :: orbsym(:, :)
        integer(ip) :: n_param

        integer(ip) :: i, j

        if (present(orbsym)) then
            n_param = 0
            do i = 1, size(n_occ, kind=ip)
                do j = n_occ(i) + 1, n_mo
                    n_param = n_param + &
                              count(orbsym(:n_occ(i), i) == orbsym(j, i), kind=ip)
                end do
            end do
        else
            n_param = sum(n_occ * (n_mo - n_occ))
        end if

    end function count_mo_params

    function mo_param_rows(mo_channels) result(rows)
        !
        ! this function returns the first and last row of every particle channel in the
        ! parameter vector, which holds the parameters of every channel, its
        ! occupied-virtual pairs of the same irrep, one after another
        !
        type(mo_channel_type), intent(in) :: mo_channels(:)
        integer(ip) :: rows(2, size(mo_channels))

        integer(ip) :: i

        rows = channel_rows([(count(mo_channels(i)%param_mask, kind=ip), i=1, &
                              size(mo_channels, kind=ip))])

    end function mo_param_rows

    function pack_ov(ov_block, mo_channel) result(packed)
        !
        ! this function packs the parameters of a particle channel from its
        ! occupied-virtual block, the entries between orbitals of the same irrep, with
        ! the occupied index running fastest
        !
        real(rp), intent(in) :: ov_block(:, :)
        type(mo_channel_type), intent(in) :: mo_channel
        real(rp), allocatable :: packed(:)

        packed = pack(ov_block, mo_channel%param_mask)

    end function pack_ov

    function unpack_ov(packed, mo_channel) result(ov_block)
        !
        ! this function unpacks the parameters of a particle channel into its
        ! occupied-virtual block, which vanishes between orbitals of different irreps
        !
        real(rp), intent(in) :: packed(:)
        type(mo_channel_type), intent(in) :: mo_channel
        real(rp), allocatable :: ov_block(:, :)

        ov_block = unpack(packed, mo_channel%param_mask, 0.0_rp)

    end function unpack_ov

    subroutine diagonalize_per_irrep(matrix, irreps, eigvals, eigvecs, settings, error)
        !
        ! this subroutine diagonalizes a symmetric matrix between orbitals of the given
        ! irreps separately within every irrep, ignoring any coupling between different
        ! irreps, which the static part of the Hessian never sees, and stores the
        ! eigenvectors of every irrep at the orbitals of that irrep, so that they never
        ! mix irreps, even where eigenvalues of different irreps are degenerate
        !
        use opentrustregion_core, only: symm_mat_diag

        real(rp), intent(in) :: matrix(:, :)
        integer(ip), intent(in) :: irreps(:)
        real(rp), intent(out) :: eigvals(:), eigvecs(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n, i, j
        integer(ip), allocatable :: idx(:)
        logical, allocatable :: done(:)
        real(rp), allocatable :: irrep_eigvals(:), irrep_eigvecs(:, :)

        ! initialize error flag
        error = 0

        ! dimension of the matrix
        n = size(irreps, kind=ip)

        ! diagonalize the block of every irrep
        eigvecs = 0.0_rp
        allocate(done(n), source=.false.)
        do i = 1, n
            if (done(i)) cycle
            idx = pack([(j, j=1, n)], irreps == irreps(i))
            allocate(irrep_eigvals(size(idx)), irrep_eigvecs(size(idx), size(idx)))
            call symm_mat_diag(matrix(idx, idx), irrep_eigvals, irrep_eigvecs, &
                               settings, error)
            if (error /= 0) return
            eigvals(idx) = irrep_eigvals
            eigvecs(idx, idx) = irrep_eigvecs
            done(idx) = .true.
            deallocate(irrep_eigvals, irrep_eigvecs)
        end do

    end subroutine diagonalize_per_irrep

end module otr_mo
