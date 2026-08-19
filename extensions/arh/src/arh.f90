! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh

    use opentrustregion, only: rp, ip, kw_len, obj_func_type, update_orbs_type, &
                               hess_x_type, project_type
    use otr_oao, only: oao_settings_type, default_oao_settings, get_energy_os_type, &
                       get_energy_cs_type

    implicit none

    ! define useful parameters
    real(rp), parameter :: rel_regularization_thresh = 1e-10_rp, &
                           symm_noise_safety_margin = 5.0_rp

    type, extends(oao_settings_type) :: arh_settings_type
        character(kw_len) :: arh_type
    contains
        procedure :: init => init_arh_settings
    end type arh_settings_type

    type(arh_settings_type), parameter :: default_arh_settings = &
        arh_settings_type(oao_settings_type = default_oao_settings, &
                          arh_type = "ms_sr1")

    ! define setting options
    character(kw_len), parameter :: arh_types(5) = &
            [character(len=kw_len) :: "arh", "symm_arh", "ms_psb", "ms_sp", "ms_sr1"]

    abstract interface
        subroutine update_dm_cs_type(dm, energy, fock, v_nonlinear, error)
            import :: rp, ip

            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), target, contiguous :: fock(:, :), v_nonlinear(:, :)
            integer(ip), intent(out) :: error
        end subroutine update_dm_cs_type

        subroutine update_dm_os_type(dm, energy, fock, v_same_spin, v_opposite_spin, &
                                       v_nonlinear, error)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), target :: fock(:, :, :), v_same_spin(:, :, :), &
                                             v_opposite_spin(:, :, :), &
                                             v_nonlinear(:, :, :)
            integer(ip), intent(out) :: error
        end subroutine update_dm_os_type
    end interface

    type :: arh_type
        type(arh_settings_type) :: settings
        integer(ip), pointer :: n_ao => null(), n_param => null(), n_particle => null()
        real(rp), pointer, contiguous :: dm_ao(:, :, :) => null()
        real(rp), pointer :: s_inv_sqrt(:, :) => null(), dm_oao(:, :, :) => null(), &
                             fock_oo(:, :, :) => null(), fock_vv(:, :, :) => null(), &
                             energy => null(), grad(:) => null(), h_diag(:) => null()
        real(rp), allocatable :: fock_oao(:, :, :), v_same_spin_oao(:, :, :), &
                                 v_opposite_spin_oao(:, :, :), &
                                 v_nonlinear_oao(:, :, :), metric_chol(:, :, :), &
                                 a_eigvecs(:, :), a_inv_eigvals(:), &
                                 a_eigvecs_comb(:, :), a_inv_eigvals_comb(:), &
                                 dm_list(:, :, :, :), fock_list(:, :, :, :), &
                                 v_same_spin_list(:, :, :, :), &
                                 v_opposite_spin_list(:, :, :, :), &
                                 v_nonlinear_list(:, :, :, :), &
                                 dm_diff(:, :, :, :), fock_diff(:, :, :, :), &
                                 v_linear_diff(:, :, :, :), &
                                 v_same_spin_diff(:, :, :, :), &
                                 v_opposite_spin_diff(:, :, :, :), &
                                 v_nonlinear_diff(:, :, :, :)
        integer(ip), allocatable :: metric_rank(:), metric_map(:, :)
        procedure(get_energy_os_type), pointer, nopass :: get_energy_os => null()
        procedure(update_dm_os_type), pointer, nopass :: update_dm_os => null()
        procedure(get_energy_cs_type), pointer, nopass :: get_energy_cs => null()
        procedure(update_dm_cs_type), pointer, nopass :: update_dm_cs => null()
    end type arh_type

    ! global variables
    type(arh_type), allocatable :: arh_object

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_arh_cs_ptr => update_orbs_arh_cs
    procedure(update_orbs_type), pointer :: update_orbs_arh_os_ptr => update_orbs_arh_os
    procedure(hess_x_type), pointer :: hess_x_arh_ptr => hess_x_arh

    ! define module procedures for different spin cases
    interface arh_factory
        module procedure arh_factory_cs, arh_factory_os
    end interface arh_factory

    contains

    subroutine arh_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, get_energy_cs, &
                              update_dm_cs, obj_func_arh_funptr, &
                              update_orbs_arh_funptr, project_arh_funptr, error, &
                              settings)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! closed-shell case
        !
        use otr_oao, only: get_energy_cs_type, obj_func_oao, project_oao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_cs_type), intent(in), pointer :: get_energy_cs
        procedure(update_dm_cs_type), intent(in), pointer :: update_dm_cs
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        real(rp), pointer, contiguous :: dm_ao_3d(:, :, :)

        ! initialize error flag
        error = 0

        ! call common setup
        dm_ao_3d(1:n_ao, 1:n_ao, 1:1) => dm_ao
        call arh_factory_common(dm_ao_3d, ao_overlap, n_particle, n_ao, error, settings)
        if (error /= 0) return
        nullify(dm_ao_3d)

        ! set pointers to functions
        arh_object%get_energy_cs => get_energy_cs
        arh_object%update_dm_cs => update_dm_cs

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_oao
        update_orbs_arh_funptr => update_orbs_arh_cs
        project_arh_funptr => project_oao

    end subroutine arh_factory_cs

    subroutine arh_factory_os(dm_ao, ao_overlap, n_particle, n_ao, get_energy_os, &
                              update_dm_os, obj_func_arh_funptr, &
                              update_orbs_arh_funptr, project_arh_funptr, error, &
                              settings)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! open-shell case
        !
        use otr_oao, only: get_energy_os_type, obj_func_oao, project_oao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_os_type), intent(in), pointer :: get_energy_os
        procedure(update_dm_os_type), intent(in), pointer :: update_dm_os
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize error flag
        error = 0

        ! call common setup
        call arh_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        if (error /= 0) return

        ! set pointers to functions
        arh_object%get_energy_os => get_energy_os
        arh_object%update_dm_os => update_dm_os

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_oao
        update_orbs_arh_funptr => update_orbs_arh_os
        project_arh_funptr => project_oao

    end subroutine arh_factory_os

    subroutine arh_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        !
        ! this function performs common ARH initialization operations
        !
        use otr_oao, only: oao_factory_common, oao_object

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! call common OAO setup
        call oao_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        if (error /= 0) return

        ! perform sanity check
        call arh_sanity_check(settings, error)
        if (error /= 0) return

        ! discard any state (in particular history and derived quantities) from a
        ! previous calculation: arh_factory is only ever called to start a new
        ! calculation, never to resume one, so state from an unrelated trajectory must
        ! never be reused, regardless of whether the new dimensions happen to match
        ! the old ones
        if (allocated(arh_object)) deallocate(arh_object)
        allocate(arh_object)

        ! set (potentially new) settings
        arh_object%settings = settings

        ! associate OAO values
        arh_object%n_ao => oao_object%n_ao
        arh_object%n_param => oao_object%n_param
        arh_object%n_particle => oao_object%n_particle
        arh_object%dm_ao => oao_object%dm_ao
        arh_object%s_inv_sqrt => oao_object%s_inv_sqrt
        arh_object%dm_oao => oao_object%dm_oao
        arh_object%fock_oo => oao_object%fock_oo
        arh_object%fock_vv => oao_object%fock_vv
        arh_object%energy => oao_object%energy
        if (allocated(oao_object%grad)) arh_object%grad => oao_object%grad
        if (allocated(oao_object%h_diag)) arh_object%h_diag => oao_object%h_diag

    end subroutine arh_factory_common

    subroutine arh_sanity_check(settings, error)
        !
        ! this subroutine performs a sanity check for ARH input parameters
        !
        use opentrustregion, only: verbosity_error, string_to_lowercase

        type(arh_settings_type), intent(inout) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! convert strings to lowercase
        settings%arh_type = string_to_lowercase(settings%arh_type)

        ! check for character options
        if (.not. any(settings%arh_type == arh_types)) then
            call settings%log("ARH type option unknown. Possible values are "// &
                              """arh"" (standard ARH), ""symm_arh"" (symmetrized "// &
                              "ARH), ""ms_sp"" (subspace-projected multisecant "// &
                              "method), ""ms_psb"" (multisecant PSB), and "// &
                              """ms_sr1"" (multisecant SR1 version).", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine arh_sanity_check

    subroutine update_orbs_arh_cs(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation 
        ! in the OAO basis and the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall for the closed-shell case
        !
        use opentrustregion, only: hess_x_type, numerical_zero
        use otr_oao, only: rotate_dm_ao, symmetric_transformation, oao_object, &
                           calculate_grad_h_diag

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list
        real(rp), allocatable :: fock_ao(:, :, :), v_nonlinear_ao(:, :, :)

        external :: dgemm

        ! initialize error flag
        error = 0

        ! check if orbitals are actually rotated
        if ((sum(abs(kappa)) > 0.0_rp) .or. &
            (abs(arh_object%energy) <= numerical_zero) .or. &
            (.not. (allocated(oao_object%grad) .and. allocated(oao_object%h_diag) &
                    .and. (allocated(arh_object%dm_list))))) then
            ! number of AOs
            n_ao = arh_object%n_ao

            ! number of particles
            n_particle = arh_object%n_particle

            ! update list of density, Fock and non-linear potential matrices
            if (allocated(arh_object%dm_list)) then
                call prepend(arh_object%dm_list, arh_object%dm_oao)
                call prepend(arh_object%fock_list, arh_object%fock_oao)
                call prepend(arh_object%v_nonlinear_list, arh_object%v_nonlinear_oao)
            else
                allocate(arh_object%dm_list(n_ao, n_ao, n_particle, 0), &
                         arh_object%fock_list(n_ao, n_ao, n_particle, 0), &
                         arh_object%v_nonlinear_list(n_ao, n_ao, n_particle, 0))
            end if

            ! rotate density matrix
            call rotate_dm_ao(kappa, n_particle, n_ao, arh_object%dm_ao, error, &
                              arh_object%dm_oao)
            if (error /= 0) return

            ! get energy, Fock matrix and non-linear potential
            allocate(fock_ao(n_ao, n_ao, n_particle), &
                     v_nonlinear_ao(n_ao, n_ao, n_particle))
            call arh_object%update_dm_cs(arh_object%dm_ao(:, :, 1), arh_object%energy, &
                                         fock_ao(:, :, 1), v_nonlinear_ao(:, :, 1), &
                                         error)
            if (error /= 0) then
                deallocate(fock_ao, v_nonlinear_ao)
                return
            end if

            ! transform Fock matrix and non-linear potential to OAO basis
            arh_object%fock_oao = &
                symmetric_transformation(arh_object%s_inv_sqrt, fock_ao)
            arh_object%v_nonlinear_oao = &
                symmetric_transformation(arh_object%s_inv_sqrt, v_nonlinear_ao)
            deallocate(fock_ao, v_nonlinear_ao)

            ! calculate gradient and Hessian diagonal
            if (.not. associated(arh_object%grad)) then
                if (.not. allocated(oao_object%grad)) &
                    allocate(oao_object%grad(arh_object%n_param))
                arh_object%grad => oao_object%grad
            end if
            if (.not. associated(arh_object%h_diag)) then
                if (.not. allocated(oao_object%h_diag)) &
                    allocate(oao_object%h_diag(arh_object%n_param))
                arh_object%h_diag => oao_object%h_diag
            end if
            call calculate_grad_h_diag(arh_object%dm_oao, arh_object%fock_oao, &
                                       n_particle, n_ao, arh_object%grad, &
                                       arh_object%h_diag, arh_object%fock_oo, &
                                       arh_object%fock_vv)

            ! prepare differences for response part of ARH Hessian
            n_list = size(arh_object%dm_list, 4)
            if (allocated(arh_object%dm_diff)) &
                deallocate(arh_object%dm_diff, arh_object%fock_diff, &
                           arh_object%v_linear_diff, arh_object%v_nonlinear_diff)
            allocate(arh_object%dm_diff(n_ao, n_ao, n_particle, n_list), &
                     arh_object%fock_diff(n_ao, n_ao, n_particle, n_list), &
                     arh_object%v_linear_diff(n_ao, n_ao, n_particle, n_list), &
                     arh_object%v_nonlinear_diff(n_ao, n_ao, n_particle, n_list))
            do i = 1, n_list
                arh_object%dm_diff(:, :, :, i) = arh_object%dm_list(:, :, :, i) - &
                                                 arh_object%dm_oao
                arh_object%fock_diff(:, :, :, i) = arh_object%fock_list(:, :, :, i) - &
                                                   arh_object%fock_oao
                arh_object%v_nonlinear_diff(:, :, :, i) = &
                    arh_object%v_nonlinear_list(:, :, :, i) - arh_object%v_nonlinear_oao
            end do

            ! linear (Coulomb and exact exchange) part of the potential difference,
            ! obtained as the remainder of the total Fock matrix difference, since the
            ! one-electron contribution cancels in the differences
            arh_object%v_linear_diff = arh_object%fock_diff - &
                                       arh_object%v_nonlinear_diff

            ! construct and diagonalize ARH metric
            if (arh_object%settings%arh_type == "ms_sr1") then
                ! linear part: this is exact since Coulomb and exact exchange are
                ! linear in the density matrix
                call get_ms_a_inv_cs(arh_object%dm_diff, arh_object%v_linear_diff, &
                                     .true., arh_object%a_eigvecs, &
                                     arh_object%a_inv_eigvals, n_ao, &
                                     arh_object%settings, error)
                if (error /= 0) return
                ! non-linear part: kept on a dedicated, separately regularized
                ! multisecant SR1 system, since the exact linear secant relationship
                ! and the non-linear one cannot be represented well by a single
                ! regularized eigenbasis
                call get_ms_a_inv_cs(arh_object%dm_diff, arh_object%v_nonlinear_diff, &
                                     .false., arh_object%a_eigvecs_comb, &
                                     arh_object%a_inv_eigvals_comb, n_ao, &
                                     arh_object%settings, error)
                if (error /= 0) return
            else
                call get_arh_metric(arh_object%dm_diff, arh_object%metric_chol, &
                                    arh_object%metric_rank, arh_object%metric_map)
            end if
        end if

        ! set outputs
        func = arh_object%energy
        grad = arh_object%grad
        h_diag = arh_object%h_diag
        hess_x_funptr => hess_x_arh
        
    end subroutine update_orbs_arh_cs

    subroutine update_orbs_arh_os(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation 
        ! in the OAO basis and the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall for the open-shell case
        !
        use opentrustregion, only: hess_x_type, numerical_zero
        use otr_oao, only: rotate_dm_ao, symmetric_transformation, oao_object, &
                           calculate_grad_h_diag

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list
        real(rp), allocatable :: fock_ao(:, :, :), fock_oao(:, :, :), &
                                 v_same_spin_ao(:, :, :), v_opposite_spin_ao(:, :, :), &
                                 v_nonlinear_ao(:, :, :)

        external :: dgemm

        ! initialize error flag
        error = 0

        ! check if orbitals are actually rotated
        if ((sum(abs(kappa)) > 0.0_rp) .or. &
            (abs(arh_object%energy) <= numerical_zero) .or. &
            (.not. (allocated(oao_object%grad) .and. allocated(oao_object%h_diag) &
                    .and. (allocated(arh_object%dm_list))))) then
            ! number of AOs
            n_ao = arh_object%n_ao

            ! number of particles
            n_particle = arh_object%n_particle

            ! update list of density and potential matrices
            if (allocated(arh_object%dm_list)) then
                call prepend(arh_object%dm_list, arh_object%dm_oao)
                call prepend(arh_object%v_same_spin_list, arh_object%v_same_spin_oao)
                call prepend(arh_object%v_opposite_spin_list, &
                             arh_object%v_opposite_spin_oao)
                call prepend(arh_object%v_nonlinear_list, arh_object%v_nonlinear_oao)
            else
                allocate(arh_object%dm_list(n_ao, n_ao, n_particle, 0), &
                         arh_object%v_same_spin_list(n_ao, n_ao, n_particle, 0), &
                         arh_object%v_opposite_spin_list(n_ao, n_ao, n_particle, 0), &
                         arh_object%v_nonlinear_list(n_ao, n_ao, n_particle, 0))
            end if

            ! rotate density matrix
            call rotate_dm_ao(kappa, n_particle, n_ao, arh_object%dm_ao, error, &
                              arh_object%dm_oao)
            if (error /= 0) return

            ! get energy, Fock matrix, same and opposite spin potentials, and 
            ! non-linear potential
            allocate(fock_ao(n_ao, n_ao, n_particle), &
                     v_same_spin_ao(n_ao, n_ao, n_particle), &
                     v_opposite_spin_ao(n_ao, n_ao, n_particle), &
                     v_nonlinear_ao(n_ao, n_ao, n_particle))
            call arh_object%update_dm_os(arh_object%dm_ao, arh_object%energy, fock_ao, &
                                         v_same_spin_ao, v_opposite_spin_ao, &
                                         v_nonlinear_ao, error)
            if (error /= 0) then
                deallocate(fock_ao, v_same_spin_ao, v_opposite_spin_ao, v_nonlinear_ao)
                return
            end if

            ! transform Fock matrix to OAO basis
            fock_oao = symmetric_transformation(arh_object%s_inv_sqrt, fock_ao)
            deallocate(fock_ao)

            ! transform same and opposite spin and non-linear potentials to OAO basis
            arh_object%v_same_spin_oao = &
                symmetric_transformation(arh_object%s_inv_sqrt, v_same_spin_ao)
            arh_object%v_opposite_spin_oao = &
                symmetric_transformation(arh_object%s_inv_sqrt, v_opposite_spin_ao)
            arh_object%v_nonlinear_oao = &
                symmetric_transformation(arh_object%s_inv_sqrt, v_nonlinear_ao)
            deallocate(v_same_spin_ao, v_opposite_spin_ao, v_nonlinear_ao)

            ! calculate gradient and Hessian diagonal
            if (.not. associated(arh_object%grad)) then
                if (.not. allocated(oao_object%grad)) &
                    allocate(oao_object%grad(arh_object%n_param))
                arh_object%grad => oao_object%grad
            end if
            if (.not. associated(arh_object%h_diag)) then
                if (.not. allocated(oao_object%h_diag)) &
                    allocate(oao_object%h_diag(arh_object%n_param))
                arh_object%h_diag => oao_object%h_diag
            end if
            call calculate_grad_h_diag(arh_object%dm_oao, fock_oao, n_particle, n_ao, &
                                       arh_object%grad, arh_object%h_diag, &
                                       arh_object%fock_oo, arh_object%fock_vv)
            deallocate(fock_oao)

            ! prepare differences for response part of ARH Hessian
            n_list = size(arh_object%dm_list, 4)
            if (allocated(arh_object%dm_diff)) &
                deallocate(arh_object%dm_diff, arh_object%v_same_spin_diff, &
                           arh_object%v_opposite_spin_diff, arh_object%v_nonlinear_diff)
            allocate(arh_object%dm_diff(n_ao, n_ao, n_particle, n_list), &
                     arh_object%v_same_spin_diff(n_ao, n_ao, n_particle, n_list), &
                     arh_object%v_opposite_spin_diff(n_ao, n_ao, n_particle, n_list), &
                     arh_object%v_nonlinear_diff(n_ao, n_ao, n_particle, n_list))
            do i = 1, n_list
                arh_object%dm_diff(:, :, :, i) = arh_object%dm_list(:, :, :, i) - &
                                                 arh_object%dm_oao
                arh_object%v_same_spin_diff(:, :, :, i) = &
                    arh_object%v_same_spin_list(:, :, :, i) - &
                    arh_object%v_same_spin_oao
                arh_object%v_opposite_spin_diff(:, :, :, i) = &
                    arh_object%v_opposite_spin_list(:, :, :, i) - &
                    arh_object%v_opposite_spin_oao
                arh_object%v_nonlinear_diff(:, :, :, i) = &
                    arh_object%v_nonlinear_list(:, :, :, i) - arh_object%v_nonlinear_oao
            end do

            ! get inverted A matrix for MS-SR1
            if (arh_object%settings%arh_type == "ms_sr1") then
                ! linear part: get spin-separated multisecant SR1 matrix for which 
                ! separation is exact since Coulomb and exact exchange are linear in 
                ! the density matrix
                call get_ms_a_inv_os_linear(arh_object%dm_diff, &
                                            arh_object%v_same_spin_diff, &
                                            arh_object%v_opposite_spin_diff, &
                                            arh_object%a_eigvecs, &
                                            arh_object%a_inv_eigvals, n_ao, &
                                            arh_object%settings, error)
                if (error /= 0) return
                ! non-linear part: get spin-combined multisecant SR1 matrix
                call get_ms_a_inv_os_nonlinear(arh_object%dm_diff, &
                                               arh_object%v_nonlinear_diff, &
                                               arh_object%a_eigvecs_comb, &
                                               arh_object%a_inv_eigvals_comb, &
                                               arh_object%settings, error)
                if (error /= 0) return
            ! get inverse of metric for ARH and related methods
            else
                call get_arh_metric(arh_object%dm_diff, arh_object%metric_chol, &
                                    arh_object%metric_rank, arh_object%metric_map)
            end if
        end if

        ! set outputs
        func = arh_object%energy
        grad = arh_object%grad
        h_diag = arh_object%h_diag
        hess_x_funptr => hess_x_arh
        
    end subroutine update_orbs_arh_os

    subroutine hess_x_arh(x, hess_x, error)
        !
        ! this function defines the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall and related methods
        !
        use otr_oao, only: unpack_asymm, project_asymm, project_symm, pack_asymm, &
                           symmetric_transformation

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i
        real(rp), allocatable :: x_full(:, :, :), hess_x_full(:, :, :), &
                                 delta_dm(:, :, :), fock_response(:, :, :)

        external :: dgemm

        ! initialize error flag
        error = 0

        ! number of AOs
        n_ao = arh_object%n_ao

        ! number of particles
        n_particle = arh_object%n_particle

        ! unpack trial vector
        x_full = unpack_asymm(x, n_particle, n_ao)

        ! get static part
        allocate(hess_x_full(n_ao, n_ao, n_particle))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, arh_object%fock_vv(:, :, i) &
                       - arh_object%fock_oo(:, :, i), n_ao, x_full(:, :, i), n_ao, &
                       0.0_rp, hess_x_full(:, :, i), n_ao)
            hess_x_full(:, :, i) = hess_x_full(:, :, i) - &
                                   transpose(hess_x_full(:, :, i))
        end do

        ! get density matrix response to trial vector
        delta_dm = project_symm(x_full, arh_object%dm_oao)

        ! approximate response contributions for linear part
        if (n_particle == 1) then
            ! MS-SR1: exact treatment of the linear part plus a dedicated,
            ! separately regularized system for the non-linear part, since a single
            ! regularized eigenbasis cannot represent both the exact linear secant
            ! relationship and the non-linear one well
            if (arh_object%settings%arh_type == "ms_sr1") then
                fock_response = &
                    get_response_contribution_cs(delta_dm, arh_object%dm_diff, &
                                                 arh_object%v_linear_diff, &
                                                 arh_object%metric_chol, &
                                                 arh_object%metric_rank, &
                                                 arh_object%metric_map, &
                                                 arh_object%a_eigvecs, &
                                                 arh_object%a_inv_eigvals, n_ao, &
                                                 n_particle, arh_object%settings) &
                    + get_response_contribution_ms_sr1_nonlinear( &
                        delta_dm, arh_object%v_nonlinear_diff, &
                        arh_object%a_eigvecs_comb, arh_object%a_inv_eigvals_comb, &
                        n_ao, n_particle)
            ! ARH and symmetric variants: the response is linear in the potential
            ! difference history and the coefficients only depend on the density
            ! matrix history, so splitting the potential into a linear and a
            ! non-linear part would leave the response unchanged and the total Fock
            ! matrix difference is used directly
            else
                fock_response = &
                    get_response_contribution_cs(delta_dm, &
                                                 arh_object%dm_diff, &
                                                 arh_object%fock_diff, &
                                                 arh_object%metric_chol, &
                                                 arh_object%metric_rank, &
                                                 arh_object%metric_map, &
                                                 arh_object%a_eigvecs, &
                                                 arh_object%a_inv_eigvals, n_ao, &
                                                 n_particle, arh_object%settings)
            end if
        else
            ! MS-SR1: exact spin-separated treatment of the linear part plus a 
            ! spin-combined treatment of the non-linear part, since the joint 
            ! spin-separated (S,Y) system of multisecant SR1 requires a 
            ! same-/opposite-spin cross relation that does not exist for non-linear 
            ! (e.g. XC) contributions
            if (arh_object%settings%arh_type == "ms_sr1") then
                fock_response = &
                    get_response_contribution_os_separated(delta_dm, &
                        arh_object%dm_diff, arh_object%v_same_spin_diff, &
                        arh_object%v_opposite_spin_diff, arh_object%metric_chol, &
                        arh_object%metric_rank, arh_object%metric_map, &
                        arh_object%a_eigvecs, arh_object%a_inv_eigvals, n_ao, &
                        n_particle, arh_object%settings) &
                    + get_response_contribution_ms_sr1_nonlinear(delta_dm, &
                        arh_object%v_nonlinear_diff, arh_object%a_eigvecs_comb, &
                        arh_object%a_inv_eigvals_comb, n_ao, n_particle)
            ! ARH and symmetric variants: coefficients only depend on the 
            ! (spin-separated) density matrix history, so the non-linear part is folded 
            ! into the same per-channel treatment as the linear part, without requiring 
            ! a dedicated spin-combined metric
            else
                fock_response = &
                    get_response_contribution_os_separated(delta_dm, &
                        arh_object%dm_diff, arh_object%v_same_spin_diff, &
                        arh_object%v_opposite_spin_diff, arh_object%metric_chol, &
                        arh_object%metric_rank, arh_object%metric_map, &
                        arh_object%a_eigvecs, arh_object%a_inv_eigvals, n_ao, &
                        n_particle, arh_object%settings, &
                        arh_object%v_nonlinear_diff)
            end if
        end if
        deallocate(x_full, delta_dm)

        ! project Fock response onto occupied-virtual and virtual-occupied subspace
        hess_x_full = hess_x_full + project_asymm(fock_response, arh_object%dm_oao)

        ! pack Hessian linear transformation
        if (arh_object%n_particle == 1) then
            hess_x = 4.0_rp * pack_asymm(hess_x_full, size(hess_x))
        else
            hess_x = 2.0_rp * pack_asymm(hess_x_full, size(hess_x))
        end if
        deallocate(hess_x_full)

    end subroutine hess_x_arh

    subroutine init_arh_settings(self, error)
        !
        ! this subroutine initializes the ARH settings
        !
        use opentrustregion, only: verbosity_error

        class(arh_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        select type(settings => self)
        type is (arh_settings_type)
            settings = default_arh_settings
        class default
            call settings%log("Augmented Roothaan-Hall settings could not be "// &
                              "initialized because initialization routine received "// &
                              "the wrong type. The type arh_settings_type was "// &
                              "likely subclassed without providing an "// &
                              "initialization routine.", verbosity_error, .true.)
            error = 1
        end select

    end subroutine init_arh_settings

    subroutine arh_deconstructor()
        !
        ! this subroutine deallocates the ARH objects
        !
        use otr_oao, only: oao_deconstructor

        if (allocated(arh_object)) deallocate(arh_object)
        call oao_deconstructor()

    end subroutine arh_deconstructor

    function get_response_contribution_cs(delta_dm, dm_diff, v_diff, metric_chol, &
                                          metric_rank, metric_map, a_eigvecs, &
                                          a_inv_eigvals, n_ao, n_particle, settings) &
        result(response)
        !
        ! this function computes the response contribution to the ARH Hessian for the 
        ! closed-shell case
        !
        real(rp), intent(in) :: delta_dm(:, :, :), dm_diff(:, :, :, :), &
                                v_diff(:, :, :, :), metric_chol(:, :, :), &
                                a_eigvecs(:, :), a_inv_eigvals(:)
        integer(ip), intent(in) :: metric_rank(:), metric_map(:, :)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp) :: response(n_ao, n_ao, n_particle)
        type(arh_settings_type), intent(in) :: settings

        integer(ip) :: n_diff, i
        real(rp), allocatable :: s_proj(:), y_proj(:), alpha_s(:), alpha_y(:), &
                                 alpha_sy(:), Y_s(:, :), sy(:), y_tilde(:), &
                                 step_norms(:), a(:, :)
        external :: dgemm, dgemv
        real(rp), external :: ddot, dnrm2

        ! number of density matrix differences
        n_diff = size(dm_diff, 4)

        ! check if there are any density matrix differences to process
        response = 0.0_rp
        if (n_diff == 0) return

        ! multisecant SR1
        if (settings%arh_type == "ms_sr1") then
            if (n_diff > 0) then
                allocate(y_proj(n_diff), alpha_y(n_diff), y_tilde(n_diff))
                ! y_proj = Y^T * delta_dm
                do i = 1, n_diff
                    y_proj(i) = ddot(n_ao * n_ao, v_diff(:, :, 1, i), 1_ip, &
                                     delta_dm(:, :, 1), 1_ip)
                end do

                ! alpha_y = M^-1 * Y^T * delta_dm
                call dgemv("T", n_diff, n_diff, 1.0_rp, a_eigvecs, n_diff, y_proj, &
                           1_ip, 0.0_rp, y_tilde, 1_ip)
                y_tilde = a_inv_eigvals * y_tilde
                call dgemv("N", n_diff, n_diff, 1.0_rp, a_eigvecs, n_diff, y_tilde, &
                           1_ip, 0.0_rp, alpha_y, 1_ip)

                ! response = Y * M^-1 * Y^T * delta_dm
                do i = 1, n_diff
                    response(:, :, 1) = response(:, :, 1) + alpha_y(i) * &
                                        v_diff(:, :, 1, i)
                end do
                deallocate(y_proj, alpha_y, y_tilde)
            end if
        ! ARH and symmetric variants
        else
            ! alpha_s = M^-1 * S^T * delta_dm
            allocate(s_proj(n_diff), alpha_s(n_diff))
            call dgemv("T", n_ao * n_ao, n_diff, 1.0_rp, dm_diff, n_ao * n_ao, &
                       delta_dm, 1_ip, 0.0_rp, s_proj, 1_ip)
            alpha_s = multiply_with_inverse_metric(s_proj, metric_chol(:, :, 1), &
                                                   metric_rank(1), metric_map(:, 1))
            deallocate(s_proj)

            ! compute components based on ARH type
            select case (settings%arh_type)
            ! symmetrized ARH
            case ("symm_arh")
                allocate(Y_s(n_ao, n_ao), y_proj(n_diff), alpha_y(n_diff))
                ! Y_s = Y * M^-1 * S^T * delta_dm
                call dgemv("N", n_ao * n_ao, n_diff, 1.0_rp, v_diff, n_ao * n_ao, &
                           alpha_s, 1_ip, 0.0_rp, Y_s, 1_ip)
                ! y_proj = Y^T * delta_dm
                call dgemv("T", n_ao * n_ao, n_diff, 1.0_rp, v_diff, n_ao * n_ao, &
                           delta_dm, 1_ip, 0.0_rp, y_proj, 1_ip)
                ! alpha_y = M^-1 * Y^T * delta_dm
                alpha_y = multiply_with_inverse_metric(y_proj, metric_chol(:, :, 1), &
                                                       metric_rank(1), metric_map(:, 1))
                ! response = 0.5 * (Y * M^-1 * S^T * delta_dm + 
                !                   S * M^-1 * Y^T * delta_dm)
                response(:, :, 1) = 0.5_rp * Y_s
                call dgemv("N", n_ao * n_ao, n_diff, 0.5_rp, dm_diff, n_ao * n_ao, &
                           alpha_y, 1_ip, 1.0_rp, response(:, :, 1), 1_ip)
                deallocate(Y_s, y_proj, alpha_y)

            ! subspace-projected multisecant and multisecant PSB
            case ("ms_sp", "ms_psb")
                allocate(Y_s(n_ao, n_ao), sy(n_diff), alpha_sy(n_diff), &
                         a(n_diff, n_diff), step_norms(n_diff))

                ! Y_s = Y * M^-1 * S^T * delta_dm
                call dgemv("N", n_ao * n_ao, n_diff, 1.0_rp, v_diff, n_ao * n_ao, &
                           alpha_s, 1_ip, 0.0_rp, Y_s, 1_ip)

                ! calculate step norms for weighting
                do i = 1, n_diff
                    step_norms(i) = dnrm2(n_ao * n_ao, dm_diff(:, :, 1, i), 1_ip)
                end do

                ! A = S^T Y
                call dgemm("T", "N", n_diff, n_diff, n_ao * n_ao, 1.0_rp, dm_diff, &
                           n_ao * n_ao, v_diff, n_ao * n_ao, 0.0_rp, a, n_diff)

                ! apply weighted symmetrization to A
                call symmetrize_weighted(a, step_norms)
                deallocate(step_norms)

                ! sy = A_sym * M^-1 * S^T * delta_dm
                call dgemv("N", n_diff, n_diff, 1.0_rp, a, n_diff, alpha_s, 1_ip, &
                           0.0_rp, sy, 1_ip)
                deallocate(a)

                ! alpha_sy = M^-1 * A_sym * M^-1 * S^T * delta_dm
                alpha_sy = multiply_with_inverse_metric(sy, metric_chol(:, :, 1), &
                                                        metric_rank(1), &
                                                        metric_map(:, 1))
                
                if (settings%arh_type == "ms_sp") then
                    ! response = S * M^-1 * A_sym * M^-1 * S^T * delta_dm
                    call dgemv("N", n_ao * n_ao, n_diff, 1.0_rp, dm_diff, n_ao * n_ao, &
                               alpha_sy, 1_ip, 0.0_rp, response(:, :, 1), 1_ip)
                else
                    allocate(y_proj(n_diff), alpha_y(n_diff))
                    ! y_proj = Y^T * delta_dm
                    call dgemv("T", n_ao * n_ao, n_diff, 1.0_rp, v_diff, n_ao * n_ao, &
                               delta_dm, 1_ip, 0.0_rp, y_proj, 1_ip)
                    ! alpha_y = M^-1 * Y^T * delta_dm
                    alpha_y = multiply_with_inverse_metric(y_proj, &
                                                           metric_chol(:, :, 1), &
                                                           metric_rank(1), &
                                                           metric_map(:, 1))
                    ! alpha_y = M^-1 * Y^T * delta_dm - M^-1 * A_sym * M^-1 * S^T * 
                    !           delta_dm
                    alpha_y = alpha_y - alpha_sy
                    ! response = Y * M^-1 * S^T * delta_dm + S * M^-1 * Y^T * delta_dm -
                    !            S * M^-1 * A_sym * M^-1 * S^T * delta_dm
                    response(:, :, 1) = Y_s
                    call dgemv("N", n_ao * n_ao, n_diff, 1.0_rp, dm_diff, n_ao * n_ao, &
                               alpha_y, 1_ip, 1.0_rp, response(:, :, 1), 1_ip)
                    deallocate(y_proj, alpha_y)
                end if
                deallocate(Y_s, sy, alpha_sy)

            ! standard ARH
            case default
                ! response = Y * M^-1 * S^T * delta_dm
                call dgemv("N", n_ao * n_ao, n_diff, 1.0_rp, v_diff, n_ao * n_ao, &
                           alpha_s, 1_ip, 0.0_rp, response(:, :, 1), 1_ip)
            end select
            deallocate(alpha_s)
        end if

    end function get_response_contribution_cs

    function get_response_contribution_ms_sr1_nonlinear(delta_dm, v_diff, a_eigvecs, &
                                                        a_inv_eigvals, n_ao, &
                                                        n_particle) result(response)
        !
        ! this function computes the MS-SR1 response contribution to the ARH Hessian
        ! arising from the non-linear part of the potential, using a single set of
        ! expansion coefficients obtained from a dedicated, separately regularized
        ! multisecant SR1 system. This is used for both the closed- and the open-shell
        ! case, since the non-linear part is never separated by spin channel
        !
        real(rp), intent(in) :: delta_dm(:, :, :), v_diff(:, :, :, :), &
                                a_eigvecs(:, :), a_inv_eigvals(:)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp) :: response(n_ao, n_ao, n_particle)

        integer(ip) :: n_diff, i, j
        real(rp), allocatable :: y_proj(:), alpha_y(:), y_tilde(:)

        external :: dgemv
        real(rp), external :: ddot

        ! number of density matrix differences
        n_diff = size(v_diff, 4)

        ! check if there are any density matrix differences to process
        response = 0.0_rp
        if (n_diff == 0) return

        allocate(y_proj(n_diff), alpha_y(n_diff), y_tilde(n_diff))
        ! y_proj = Y^T * delta_dm
        do i = 1, n_diff
            y_proj(i) = ddot(n_ao * n_ao * n_particle, v_diff(:, :, :, i), 1_ip, &
                             delta_dm, 1_ip)
        end do
        ! alpha_y = M^-1 * Y^T * delta_dm
        call dgemv("T", n_diff, n_diff, 1.0_rp, a_eigvecs, n_diff, y_proj, 1_ip, &
                   0.0_rp, y_tilde, 1_ip)
        y_tilde = a_inv_eigvals * y_tilde
        call dgemv("N", n_diff, n_diff, 1.0_rp, a_eigvecs, n_diff, y_tilde, 1_ip, &
                   0.0_rp, alpha_y, 1_ip)
        ! response = Y * M^-1 * Y^T * delta_dm
        do i = 1, n_diff
            do j = 1, n_particle
                response(:, :, j) = response(:, :, j) + alpha_y(i) * v_diff(:, :, j, i)
            end do
        end do
        deallocate(y_proj, alpha_y, y_tilde)

    end function get_response_contribution_ms_sr1_nonlinear

    function get_response_contribution_os_separated(delta_dm, dm_diff, &
                                                    v_same_spin_diff, &
                                                    v_opposite_spin_diff, &
                                                    metric_chol, metric_rank, &
                                                    metric_map, a_eigvecs, &
                                                    a_inv_eigvals, n_ao, n_particle, &
                                                    settings, v_nonlinear_diff) &
        result(response)
        !
        ! this function computes the spin-separated response contributions to the ARH 
        ! Hessian, for ARH and related symmetric methods; if the optional non-linear 
        ! potential difference is present, it is also included, using the same 
        ! per-channel expansion coefficients as the linear part, without any 
        ! same-/opposite-spin cross-mixing, since no such decomposition exists for the 
        ! non-linear part; for MS-SR1, the non-linear contribution is ignored here, 
        ! since its non-linear contribution instead requires a dedicated spin-combined 
        ! MS-SR1 system
        !
        real(rp), intent(in) :: delta_dm(:, :, :), dm_diff(:, :, :, :), &
                                v_same_spin_diff(:, :, :, :), &
                                v_opposite_spin_diff(:, :, :, :), &
                                metric_chol(:, :, :), a_eigvecs(:, :), a_inv_eigvals(:)
        integer(ip), intent(in) :: metric_rank(:), metric_map(:, :)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp), intent(in), optional :: v_nonlinear_diff(:, :, :, :)
        real(rp) :: response(n_ao, n_ao, n_particle)
        type(arh_settings_type), intent(in) :: settings

        integer(ip) :: n_diff, i, j, k
        real(rp), allocatable :: y_proj(:), alpha_y(:), y_tilde(:), s_proj(:, :), &
                                 alpha_s(:, :), Y_s(:, :), sy(:), alpha_sy(:), &
                                 step_norms(:, :), a_same(:, :, :), a_opp(:, :, :), &
                                 v_same_spin_diff_eff(:, :, :, :)

        real(rp), external :: ddot, dnrm2

        n_diff = size(dm_diff, 4)

        response = 0.0_rp
        if (n_diff == 0) return

        ! multisecant SR1
        if (settings%arh_type == "ms_sr1") then
            allocate(y_proj(2 * n_diff), alpha_y(2 * n_diff), y_tilde(2 * n_diff))
            ! y_proj = Y^T * delta_dm
            do i = 1, n_diff
                y_proj(i) = ddot(n_ao * n_ao, v_same_spin_diff(:, :, 1, i), 1_ip, &
                                 delta_dm(:, :, 1), 1_ip) + &
                            ddot(n_ao * n_ao, v_opposite_spin_diff(:, :, 2, i), 1_ip, &
                                 delta_dm(:, :, 2), 1_ip)
                y_proj(n_diff + i) = ddot(n_ao * n_ao, &
                                          v_opposite_spin_diff(:, :, 1, i), 1_ip, &
                                          delta_dm(:, :, 1), 1_ip) + &
                                     ddot(n_ao * n_ao, v_same_spin_diff(:, :, 2, i), &
                                          1_ip, delta_dm(:, :, 2), 1_ip)
            end do
            ! alpha_y = M^-1 * Y^T * delta_dm
            call dgemv("T", 2 * n_diff, 2 * n_diff, 1.0_rp, a_eigvecs, 2 * n_diff, &
                       y_proj, 1_ip, 0.0_rp, y_tilde, 1_ip)
            y_tilde = a_inv_eigvals * y_tilde
            call dgemv("N", 2 * n_diff, 2 * n_diff, 1.0_rp, a_eigvecs, 2 * n_diff, &
                       y_tilde, 1_ip, 0.0_rp, alpha_y, 1_ip)
            ! response = Y * alpha
            do i = 1, n_diff
                response(:, :, 1) = response(:, :, 1) + alpha_y(i) * &
                                    v_same_spin_diff(:, :, 1, i) + &
                                    alpha_y(n_diff + i) * &
                                    v_opposite_spin_diff(:, :, 1, i)
                response(:, :, 2) = response(:, :, 2) + alpha_y(i) * &
                                    v_opposite_spin_diff(:, :, 2, i) + &
                                    alpha_y(n_diff + i) * v_same_spin_diff(:, :, 2, i)
            end do
            deallocate(y_proj, alpha_y, y_tilde)

        ! ARH and symmetric variants
        else
            ! fold the non-linear potential difference, if present, into the same-spin
            ! potential difference without any same-/opposite-spin cross-mixing, since
            ! no such decomposition exists for it
            allocate(v_same_spin_diff_eff, source=v_same_spin_diff)
            if (present(v_nonlinear_diff)) &
                v_same_spin_diff_eff = v_same_spin_diff_eff + v_nonlinear_diff

            ! alpha_s = M^-1 * S^T * delta_dm
            allocate(s_proj(n_diff, n_particle), alpha_s(n_diff, n_particle))
            do i = 1, n_diff
                do j = 1, n_particle
                    s_proj(i, j) = ddot(n_ao * n_ao, dm_diff(:, :, j, i), 1_ip, &
                                        delta_dm(:, :, j), 1_ip)
                end do
            end do
            do j = 1, n_particle
                alpha_s(:, j) = &
                    multiply_with_inverse_metric(s_proj(:, j), metric_chol(:, :, j), &
                                                 metric_rank(j), metric_map(:, j))
            end do
            deallocate(s_proj)

            ! compute components based on ARH type
            select case (settings%arh_type)
            ! symmetrized ARH
            case ("symm_arh")
                allocate(Y_s(n_ao, n_ao), y_proj(n_diff), alpha_y(n_diff))
                do j = 1, n_particle
                    ! Y_s = Y * M^-1 * S^T * delta_dm
                    Y_s = 0.0_rp
                    do i = 1, n_diff
                        Y_s = Y_s + alpha_s(i, j) * v_same_spin_diff_eff(:, :, j, i) + &
                              alpha_s(i, 3 - j) * v_opposite_spin_diff(:, :, j, i)
                    end do
                    ! y_proj = Y^T * delta_dm
                    do i = 1, n_diff
                        y_proj(i) = ddot(n_ao * n_ao, &
                                         v_same_spin_diff_eff(:, :, j, i), 1_ip, &
                                         delta_dm(:, :, j), 1_ip) + &
                                    ddot(n_ao * n_ao, &
                                         v_opposite_spin_diff(:, :, 3 - j, i), 1_ip, &
                                         delta_dm(:, :, 3 - j), 1_ip)
                    end do
                    ! alpha_y = M^-1 * Y^T * delta_dm
                    alpha_y = &
                        multiply_with_inverse_metric(y_proj, metric_chol(:, :, j), &
                                                     metric_rank(j), metric_map(:, j))
                    ! response = 0.5 * (Y * M^-1 * S^T * delta_dm + 
                    !                   S * M^-1 * Y^T * delta_dm)
                    response(:, :, j) = 0.5_rp * Y_s
                    do i = 1, n_diff
                        response(:, :, j) = response(:, :, j) + 0.5_rp * alpha_y(i) * &
                                            dm_diff(:, :, j, i)
                    end do  
                end do
                deallocate(Y_s, y_proj, alpha_y)

            ! subspace-projected multisecant and multisecant PSB
            case ("ms_sp", "ms_psb")
                allocate(Y_s(n_ao, n_ao), sy(n_diff), alpha_sy(n_diff), &
                         a_same(n_diff, n_diff, n_particle), &
                         a_opp(n_diff, n_diff, n_particle), &
                         step_norms(n_diff, n_particle))

                ! compute step norms for weighing
                do j = 1, n_particle
                    do i = 1, n_diff
                        step_norms(i, j) = dnrm2(n_ao * n_ao, dm_diff(:, :, j, i), 1_ip)
                    end do
                end do

                ! A = S^T Y
                do j = 1, n_particle
                    do k = 1, n_diff
                        do i = 1, n_diff
                            a_same(i, k, j) = &
                                ddot(n_ao * n_ao, dm_diff(:, :, j, i), 1_ip, &
                                     v_same_spin_diff_eff(:, :, j, k), 1_ip)
                            a_opp(i, k, j) = &
                                ddot(n_ao * n_ao, dm_diff(:, :, j, i), 1_ip, &
                                     v_opposite_spin_diff(:, :, j, k), 1_ip)
                        end do
                    end do
                end do

                ! symmetrize A_same within each spin channel
                do j = 1, n_particle
                    call symmetrize_weighted(a_same(:, :, j), step_norms(:, j))
                end do

                ! cross-symmetrize A_opp between spin channels
                call cross_symmetrize_weighted(a_opp(:, :, 1), a_opp(:, :, 2), &
                                               step_norms(:, 1), step_norms(:, 2))
                deallocate(step_norms)

                do j = 1, n_particle
                    ! Y_s = Y * M^-1 * S^T * delta_dm
                    Y_s = 0.0_rp
                    do i = 1, n_diff
                        Y_s = Y_s + alpha_s(i, j) * v_same_spin_diff_eff(:, :, j, i) + &
                                    alpha_s(i, 3 - j) * v_opposite_spin_diff(:, :, j, i)
                    end do

                    ! sy = A_sym * M^-1 * S^T * delta_dm
                    do i = 1, n_diff
                        sy(i) = ddot(n_diff, a_same(i, :, j), 1_ip, alpha_s(:, j), &
                                     1_ip) + &
                                ddot(n_diff, a_opp(i, :, j), 1_ip, alpha_s(:, 3 - j), &
                                     1_ip)
                    end do

                    ! alpha_sy = M^-1 * A_sym * M^-1 * S^T * delta_dm
                    alpha_sy = &
                        multiply_with_inverse_metric(sy, metric_chol(:, :, j), &
                                                     metric_rank(j), metric_map(:, j))

                    if (settings%arh_type == "ms_sp") then
                        ! response = S * M^-1 * A_sym * M^-1 * S^T * delta_dm
                        do i = 1, n_diff
                            response(:, :, j) = response(:, :, j) + alpha_sy(i) * &
                                                dm_diff(:, :, j, i)
                        end do
                    else
                        allocate(y_proj(n_diff), alpha_y(n_diff))
                        ! y_proj = Y^T * delta_dm
                        do i = 1, n_diff
                            y_proj(i) = &
                                ddot(n_ao * n_ao, v_same_spin_diff_eff(:, :, j, i), &
                                     1_ip, delta_dm(:, :, j), 1_ip) + &
                                ddot(n_ao * n_ao, &
                                     v_opposite_spin_diff(:, :, 3 - j, i), 1_ip, &
                                     delta_dm(:, :, 3 - j), 1_ip)
                        end do
                        ! alpha_y = M^-1 * Y^T * delta_dm
                        alpha_y = &
                            multiply_with_inverse_metric(y_proj, metric_chol(:, :, j), &
                                                         metric_rank(j), &
                                                         metric_map(:, j))
                        ! alpha_y = M^-1 * Y^T * delta_dm - 
                        !           M^-1 * A_sym * M^-1 * S^T * delta_dm
                        alpha_y = alpha_y - alpha_sy
                        ! response = Y * M^-1 * S^T * delta_dm + 
                        !            S * M^-1 * Y^T * delta_dm -
                        !            S * M^-1 * A_sym * M^-1 * S^T * delta_dm
                        response(:, :, j) = Y_s
                        do i = 1, n_diff
                            response(:, :, j) = response(:, :, j) + alpha_y(i) * &
                                                dm_diff(:, :, j, i)
                        end do  
                        deallocate(y_proj, alpha_y)
                    end if
                end do
                deallocate(Y_s, sy, alpha_sy, a_same, a_opp)

            ! standard ARH
            case default
                ! response = Y * M^-1 * S^T * delta_dm
                do j = 1, n_particle
                    do i = 1, n_diff
                        response(:, :, j) = response(:, :, j) + alpha_s(i, j) * &
                                            v_same_spin_diff_eff(:, :, j, i) + &
                                            alpha_s(i, 3 - j) * &
                                            v_opposite_spin_diff(:, :, j, i)
                    end do
                end do
            end select
            deallocate(alpha_s)
        end if

    end function get_response_contribution_os_separated

    function multiply_with_inverse_metric(vec, chol, rank, map) result(result_vec)
        !
        ! this function multiplies a vector with the pseudoinverse of the ARH metric
        ! using forward and backward substitution over the numerical rank
        !
        real(rp), intent(in) :: vec(:), chol(:, :)
        integer(ip), intent(in) :: rank, map(:)
        real(rp) :: result_vec(size(vec))

        integer(ip) :: n_dm, i
        real(rp), allocatable :: perm_vec(:)
        external :: dtrsv

        n_dm = size(vec, 1)
        result_vec = 0.0_rp
        if (rank == 0) return

        ! forward permutation: perm_vec = P^T * vec
        allocate(perm_vec(n_dm))
        do i = 1, n_dm
            perm_vec(i) = vec(map(i))
        end do

        ! forward substitution: solve R^T * y = perm_vec(1:rank)
        call dtrsv("U", "T", "N", rank, chol, n_dm, perm_vec, 1_ip)

        ! filter out linear dependencies
        if (rank < n_dm) perm_vec(rank+1:n_dm) = 0.0_rp

        ! backward substitution: solve R * c = y
        call dtrsv("U", "N", "N", rank, chol, n_dm, perm_vec, 1_ip)

        ! backward permutation to original basis: result_vec = P * perm_vec
        do i = 1, n_dm
            result_vec(map(i)) = perm_vec(i)
        end do
        deallocate(perm_vec)

    end function multiply_with_inverse_metric

    subroutine symmetrize_weighted(a, step_norms)
        !
        ! this subroutine performs a weighted symmetrization of a square matrix 
        ! a(i, k), blending a(i, k) with a(k, i) using a weight biased toward the entry 
        ! associated with the larger step norm
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(inout) :: a(:, :)
        real(rp), intent(in) :: step_norms(:)

        integer(ip) :: i, k, n
        real(rp) :: weight

        n = size(a, 1)
        do k = 2, n
            do i = 1, k - 1
                if (step_norms(i) + step_norms(k) > numerical_zero) then
                    weight = step_norms(k) / (step_norms(i) + step_norms(k))
                else
                    weight = 0.5_rp
                end if
                a(i, k) = weight * a(k, i) + (1.0_rp - weight) * a(i, k)
                a(k, i) = a(i, k)
            end do
        end do

    end subroutine symmetrize_weighted

    subroutine cross_symmetrize_weighted(a12, a21, step_norms1, step_norms2)
        !
        ! this subroutine performs a weighted symmetrization between two related 
        ! off-diagonal blocks a12(i, k) and a21(k, i) of a larger matrix, analogous to
        ! symmetrize_weighted but for blocks that are not self-transposed
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(inout) :: a12(:, :), a21(:, :)
        real(rp), intent(in) :: step_norms1(:), step_norms2(:)

        integer(ip) :: i, k
        real(rp) :: weight

        do k = 1, size(a12, 2)
            do i = 1, size(a12, 1)
                if (step_norms1(i) + step_norms2(k) > numerical_zero) then
                    weight = step_norms2(k) / (step_norms1(i) + step_norms2(k))
                else
                    weight = 0.5_rp
                end if
                a12(i, k) = weight * a21(k, i) + (1.0_rp - weight) * a12(i, k)
                a21(k, i) = a12(i, k)
            end do
        end do

    end subroutine cross_symmetrize_weighted

    subroutine symmetrize_exact(a)
        !
        ! this subroutine performs a plain, unweighted symmetrization of a square
        ! matrix, appropriate when a is exactly symmetric in exact arithmetic so that 
        ! any observed asymmetry is pure numerical noise rather than a genuine
        ! inconsistency to average over
        !
        real(rp), intent(inout) :: a(:, :)

        a = 0.5_rp * (a + transpose(a))

    end subroutine symmetrize_exact

    function noise_threshold(a, settings, error) result(thresh)
        !
        ! this function estimates the numerical noise floor of a matrix a that is
        ! exactly symmetric in exact arithmetic, using the spectral norm (largest
        ! singular value) of its observed antisymmetric part as a direct measurement of 
        ! that noise
        !
        use opentrustregion, only: numerical_zero, verbosity_error

        real(rp), intent(in) :: a(:, :)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        real(rp) :: thresh

        real(rp), allocatable :: anti(:, :), s(:), work(:)
        real(rp) :: dummy(1, 1)
        integer(ip) :: n, lwork, info
        character(300) :: msg
        external :: dgesvd

        ! initialize error flag and threshold in case of error
        error = 0
        thresh = 0.0_rp

        n = size(a, 1)

        ! get antisymmetric part
        allocate (anti(n, n), s(n))
        anti = 0.5_rp * (a - transpose(a))

        ! workspace query
        allocate (work(1))
        call dgesvd('N', 'N', n, n, anti, n, s, dummy, 1_ip, dummy, 1_ip, work, &
                    -1_ip, info)
        lwork = int(work(1), kind=ip)
        deallocate (work)
        allocate (work(lwork))

        ! compute singular values only
        call dgesvd('N', 'N', n, n, anti, n, s, dummy, 1_ip, dummy, 1_ip, work, &
                    lwork, info)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Singular value decomposition failed: Error in "// &
                "DGESVD, info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! set threshold to largest singular value times safety margin
        thresh = symm_noise_safety_margin * max(s(1), numerical_zero)
        deallocate (anti, s, work)

    end function noise_threshold

    function regularized_eigval_inv(eig_vals, eig_val_thresh) result(eig_vals_inv)
        !
        ! this function returns Tikhonov-regularized pseudoinverse eigenvalues, damping
        ! directions close to numerical noise relative to the given threshold
        !
        real(rp), intent(in) :: eig_vals(:), eig_val_thresh
        real(rp) :: eig_vals_inv(size(eig_vals))

        eig_vals_inv = eig_vals / (eig_vals**2 + eig_val_thresh**2)

    end function regularized_eigval_inv

    function truncated_eigval_inv(eig_vals, eig_val_thresh) result(eig_vals_inv)
        !
        ! this function returns a hard-truncated pseudoinverse: the exact inverse is
        ! kept for eigenvalues above the given threshold, while eigenvalues at or
        ! below it are discarded entirely rather than smoothly damped
        !
        real(rp), intent(in) :: eig_vals(:), eig_val_thresh
        real(rp) :: eig_vals_inv(size(eig_vals))

        where (abs(eig_vals) > eig_val_thresh)
            eig_vals_inv = 1.0_rp / eig_vals
        elsewhere
            eig_vals_inv = 0.0_rp
        end where

    end function truncated_eigval_inv

    subroutine get_arh_metric(dm_diff, chol, rank, map)
        !
        ! this subroutine calculates the augmented Roothaan-Hall metric factorized via 
        ! an unpivoted Cholesky decomposition, since the density matrix differences are 
        ! saved in reverse order, linearly dependent older vectors are skipped while 
        ! preserving the chronological sequence for the active basis
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: dm_diff(:, :, :, :)
        real(rp), intent(out), allocatable :: chol(:, :, :)
        integer(ip), intent(out), allocatable :: rank(:), map(:, :)

        integer(ip) :: n_ao, n_particle, n_dm, i, j, k
        integer(ip) :: n_accepted, n_rejected
        real(rp) :: raw_diagonal
        real(rp), allocatable :: metric(:, :), work(:)
        real(rp) :: tol
        external :: dtrsv
        real(rp), external :: ddot

        n_ao = size(dm_diff, 1)
        n_particle = size(dm_diff, 3)
        n_dm = size(dm_diff, 4)

        ! allocate metric arrays
        allocate(chol(n_dm, n_dm, n_particle), rank(n_particle), &
                 map(n_dm, n_particle), metric(n_dm, n_dm), work(n_dm))

        chol = 0.0_rp
        rank = 0
        do k = 1, n_particle
            ! initialize tolerance with maximum diagonal element
            tol = 0.0_rp

            ! generate full Gram matrix
            do j = 1, n_dm
                do i = 1, j
                    ! compute Tr(dm_diff_i * dm_diff_j)
                    metric(i, j) = ddot(n_ao * n_ao, dm_diff(:, :, k, j), 1_ip, &
                                        dm_diff(:, :, k, i), 1_ip)
                    metric(j, i) = metric(i, j)
                end do
                tol = max(tol, metric(j, j))
            end do

            if (n_dm > 0) then
                ! tolerance for linear dependencies
                tol = n_dm * numerical_zero * tol

                ! loop chronologically through history
                n_accepted = 0
                n_rejected = 0
                do i = 1, n_dm
                    raw_diagonal = metric(i, i)

                    ! project current column onto the accepted columns
                    if (n_accepted > 0) then
                        do j = 1, n_accepted
                            work(j) = metric(map(j, k), i)
                        end do

                        ! solve R_accepted^T * work = T_accepted_vs_current
                        call dtrsv("U", "T", "N", n_accepted, chol(:, :, k), n_dm, &
                                   work, 1_ip)
                        
                        ! subtract projections to find remaining orthogonal magnitude
                        raw_diagonal = raw_diagonal - &
                                      ddot(n_accepted, work(1:n_accepted), 1_ip, &
                                           work(1:n_accepted), 1_ip)
                    end if
                    
                    ! check for linear dependency
                    if (raw_diagonal < tol) then
                        ! dependency found: store index at the back of map and skip
                        ! column
                        n_rejected = n_rejected + 1
                        map(n_dm - n_rejected + 1, k) = i
                    else
                        ! independent: accept column and append to the active factors
                        n_accepted = n_accepted + 1
                        map(n_accepted, k) = i

                        if (n_accepted > 1) then
                            chol(1:n_accepted-1, n_accepted, k) = work(1:n_accepted-1)
                        end if
                        chol(n_accepted, n_accepted, k) = sqrt(max(0.0_rp, &
                                                                   raw_diagonal))
                    end if
                end do
                rank(k) = n_accepted
            end if
        end do
        deallocate(metric, work)

    end subroutine get_arh_metric

    subroutine get_ms_a_inv_cs(dm_diff, v_diff, linear, eig_vecs, eig_vals_inv, n_ao, &
                               settings, error)
        !
        ! this subroutine computes the pseudoinverse multisecant SR1 matrix for the 
        ! closed-shell case; for the linear (Coulomb and exact exchange) part; for the 
        ! linear part, A is exactly symmetric in exact arithmetic, so its observed 
        ! asymmetry is used directly as a calibrated numerical noise floor to form the 
        ! pseudoinverse; for the non-linear part, no such exact relationship holds, so 
        ! a weighted symmetrization is used with a relative threshold for the 
        ! regularization
        !
        use opentrustregion, only: symm_mat_diag, numerical_zero

        real(rp), intent(in) :: dm_diff(:, :, :, :), v_diff(:, :, :, :)
        logical, intent(in) :: linear
        real(rp), intent(out), allocatable :: eig_vecs(:, :), eig_vals_inv(:)
        integer(ip), intent(in) :: n_ao
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, i
        real(rp), allocatable :: a(:, :), eig_vals(:), step_norms(:)
        real(rp) :: eig_val_thresh
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! handle empty history
        n_dm = size(dm_diff, 4)
        if (n_dm == 0) then
            allocate(eig_vecs(0, 0), eig_vals_inv(0))
            return
        end if

        ! A = S^T Y
        allocate(eig_vecs(n_dm, n_dm), eig_vals_inv(n_dm), a(n_dm, n_dm), &
                 eig_vals(n_dm))
        call dgemm("T", "N", n_dm, n_dm, n_ao * n_ao, 1.0_rp, dm_diff, &
                   n_ao * n_ao, v_diff, n_ao * n_ao, 0.0_rp, a, n_dm)

        ! symmetrize
        if (linear) then
            ! measure the numerical noise floor from the observed asymmetry before
            ! enforcing symmetry
            eig_val_thresh = noise_threshold(a, settings, error)
            if (error /= 0) return
            call symmetrize_exact(a)
        else
            ! calculate step norms for weighing
            allocate(step_norms(n_dm))
            do i = 1, n_dm
                step_norms(i) = dnrm2(n_ao * n_ao, dm_diff(:, :, 1, i), 1_ip)
            end do

            ! enforce symmetry
            call symmetrize_weighted(a, step_norms)
            deallocate(step_norms)
        end if

        ! perform spectral decomposition
        call symm_mat_diag(a, eig_vals, eig_vecs, settings, error)
        if (error /= 0) return

        ! construct inverse
        if (linear) then
            eig_vals_inv = truncated_eigval_inv(eig_vals, eig_val_thresh)
        else
            eig_val_thresh = max(rel_regularization_thresh * maxval(abs(eig_vals)), &
                                 numerical_zero)
            eig_vals_inv = regularized_eigval_inv(eig_vals, eig_val_thresh)
        end if
        deallocate(a, eig_vals)

    end subroutine get_ms_a_inv_cs

    subroutine get_ms_a_inv_os_linear(dm_diff, v_same_spin_diff, v_opposite_spin_diff, &
                                      eig_vecs, eig_vals_inv, n_ao, settings, error)
            !
            ! this subroutine computes the pseudoinverse multisecant SR1 matrix in a
            ! spin-separated manner for the linear (Coulomb and exact exchange) part in the 
            ! open-shell case; A is exactly symmetric in exact arithmetic, so its observed 
            ! asymmetry is used directly as a calibrated numerical noise floor to form the 
            ! pseudoinverse
            !
            use opentrustregion, only: symm_mat_diag

            real(rp), intent(in) :: dm_diff(:, :, :, :), v_same_spin_diff(:, :, :, :), &
                                    v_opposite_spin_diff(:, :, :, :)
            real(rp), intent(out), allocatable :: eig_vecs(:, :), eig_vals_inv(:)
            integer(ip), intent(in) :: n_ao
            type(arh_settings_type), intent(in) :: settings
            integer(ip), intent(out) :: error

            integer(ip) :: n_dm, i, k
            real(rp), allocatable :: a(:, :), eig_vals(:)
            real(rp) :: eig_val_thresh

            real(rp), external :: ddot

            ! initialize error flag
            error = 0

            ! handle empty history
            n_dm = size(dm_diff, 4)
            if (n_dm == 0) then
                allocate(eig_vecs(0, 0), eig_vals_inv(0))
                return
            end if

            ! A = S^T Y
            allocate(a(2 * n_dm, 2 * n_dm))
            do k = 1, n_dm
                do i = 1, n_dm
                    a(i, k) = ddot(n_ao * n_ao, dm_diff(:, :, 1, i), 1_ip, &
                                   v_same_spin_diff(:, :, 1, k), 1_ip)
                    a(i, n_dm + k) = ddot(n_ao * n_ao, dm_diff(:, :, 1, i), 1_ip, &
                                          v_opposite_spin_diff(:, :, 1, k), 1_ip)
                    a(n_dm + i, k) = ddot(n_ao * n_ao, dm_diff(:, :, 2, i), 1_ip, &
                                          v_opposite_spin_diff(:, :, 2, k), 1_ip)
                    a(n_dm + i, n_dm + k) = ddot(n_ao * n_ao, dm_diff(:, :, 2, i), &
                                                 1_ip, v_same_spin_diff(:, :, 2, k), &
                                                 1_ip)
                end do
            end do

            ! measure the numerical noise floor from the observed asymmetry before
            ! enforcing symmetry
            eig_val_thresh = noise_threshold(a, settings, error)
            if (error /= 0) return
            call symmetrize_exact(a)

            ! perform spectral decomposition
            allocate(eig_vecs(2 * n_dm, 2 * n_dm), eig_vals(2 * n_dm))
            call symm_mat_diag(a, eig_vals, eig_vecs, settings, error)
            if (error /= 0) return

            ! construct regularized inverse
            allocate(eig_vals_inv(2 * n_dm))
            eig_vals_inv = truncated_eigval_inv(eig_vals, eig_val_thresh)

            deallocate(a, eig_vals)

    end subroutine get_ms_a_inv_os_linear

    subroutine get_ms_a_inv_os_nonlinear(dm_diff, v_diff, eig_vecs, eig_vals_inv, &
                                         settings, error)
        !
        ! this subroutine computes the pseudoinverse multisecant SR1 matrix in a
        ! spin-combined manner for the non-linear part in the open-shell case; a 
        ! weighted symmetrization is used with a relative threshold for the 
        ! regularization
        !
        use opentrustregion, only: symm_mat_diag, numerical_zero

        real(rp), intent(in) :: dm_diff(:, :, :, :), v_diff(:, :, :, :)
        real(rp), intent(out), allocatable :: eig_vecs(:, :), eig_vals_inv(:)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, flat_len, i, j
        real(rp), allocatable :: a(:, :), eig_vals(:), step_norms(:)
        real(rp) :: eig_val_thresh
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        ! handle empty history
        n_dm = size(dm_diff, 4)
        if (n_dm == 0) then
            allocate(eig_vecs(0, 0), eig_vals_inv(0))
            return
        end if

        ! A = S^T Y
        allocate(a(n_dm, n_dm))
        flat_len = size(dm_diff, 1) * size(dm_diff, 2) * size(dm_diff, 3)
        do j = 1, n_dm
            do i = 1, n_dm
                a(i, j) = ddot(flat_len, dm_diff(:, :, :, i), 1_ip, &
                               v_diff(:, :, :, j), 1_ip)
            end do
        end do

        ! calculate step norms for weighing
        allocate(step_norms(n_dm))
        do i = 1, n_dm
            step_norms(i) = dnrm2(flat_len, dm_diff(:, :, :, i), 1_ip)
        end do

        ! enforce symmetry
        call symmetrize_weighted(a, step_norms)
        deallocate(step_norms)

        ! perform spectral decomposition
        allocate(eig_vecs(n_dm, n_dm), eig_vals(n_dm))
        call symm_mat_diag(a, eig_vals, eig_vecs, settings, error)
        if (error /= 0) return

        ! construct regularized inverse
        allocate(eig_vals_inv(n_dm))
        eig_val_thresh = max(rel_regularization_thresh * maxval(abs(eig_vals)), &
                             numerical_zero)
        eig_vals_inv = regularized_eigval_inv(eig_vals, eig_val_thresh)
        deallocate(a, eig_vals)

    end subroutine get_ms_a_inv_os_nonlinear

    subroutine prepend(list, new_array)
        !
        ! this subroutine prepends an array to a list of arrays of equal dimension
        !
        real(rp), intent(inout), allocatable :: list(:, :, :, :)
        real(rp), intent(in) :: new_array(:, :, :)

        integer(ip) :: n1, n2, n3
        real(rp), allocatable :: temp(:, :, :, :)

        n1 = size(new_array, 1)
        n2 = size(new_array, 2)
        n3 = size(new_array, 3)

        allocate(temp(n1, n2, n3, size(list, 4) + 1))
        temp(:, :, :, 1) = new_array
        temp(:, :, :, 2:size(list, 4) + 1) = list
        deallocate(list)
        list = temp

    end subroutine prepend

end module otr_arh
