! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh

    use opentrustregion, only: rp, ip, kw_len, obj_func_type, update_orbs_type, &
                               hess_x_type, precond_type, precond_pd_type, project_type
    use otr_oao, only: oao_settings_type, default_oao_settings

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
                                 v_nonlinear_oao(:, :, :), metric_inv(:, :), &
                                 a_sym(:, :), a_inv(:, :), a_inv_comb(:, :), &
                                 dm_list(:, :, :, :), fock_list(:, :, :, :), &
                                 v_same_spin_list(:, :, :, :), &
                                 v_opposite_spin_list(:, :, :, :), &
                                 v_nonlinear_list(:, :, :, :), &
                                 linear_potential_dirs(:, :), &
                                 nonlinear_potential_dirs(:, :), dm_dirs(:, :), &
                                 potential_dirs(:, :), expansion_dirs(:, :), &
                                 projection_dirs(:, :), coupling_matrix(:, :)
        procedure(update_dm_os_type), pointer, nopass :: update_dm_os => null()
        procedure(update_dm_cs_type), pointer, nopass :: update_dm_cs => null()
    end type arh_type

    ! global variables
    type(arh_type), allocatable :: arh_object

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_arh_cs_ptr => update_orbs_arh_cs
    procedure(update_orbs_type), pointer :: update_orbs_arh_os_ptr => update_orbs_arh_os
    procedure(hess_x_type), pointer :: hess_x_arh_ptr => hess_x_arh
    procedure(precond_type), pointer :: precond_arh_ptr => precond_arh

    ! define module procedures for different spin cases
    interface arh_factory
        module procedure arh_factory_cs, arh_factory_os
    end interface arh_factory

    contains

    subroutine arh_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, get_energy_cs, &
                              update_dm_cs, obj_func_arh_funptr, &
                              update_orbs_arh_funptr, precond_arh_funptr, &
                              precond_pd_arh_funptr, project_arh_funptr, error, &
                              settings)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! closed-shell case
        !
        use otr_oao, only: get_energy_cs_type, obj_func_oao, precond_pd_oao, &
                           project_oao, oao_object

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_cs_type), intent(in), pointer :: get_energy_cs
        procedure(update_dm_cs_type), intent(in), pointer :: update_dm_cs
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(precond_type), intent(out), pointer :: precond_arh_funptr
        procedure(precond_pd_type), intent(out), pointer :: precond_pd_arh_funptr
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
        oao_object%get_energy_cs => get_energy_cs
        arh_object%update_dm_cs => update_dm_cs

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_oao
        update_orbs_arh_funptr => update_orbs_arh_cs
        precond_arh_funptr => precond_arh
        precond_pd_arh_funptr => precond_pd_oao
        project_arh_funptr => project_oao

    end subroutine arh_factory_cs

    subroutine arh_factory_os(dm_ao, ao_overlap, n_particle, n_ao, get_energy_os, &
                              update_dm_os, obj_func_arh_funptr, &
                              update_orbs_arh_funptr, precond_arh_funptr, &
                              precond_pd_arh_funptr, project_arh_funptr, error, &
                              settings)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! open-shell case
        !
        use otr_oao, only: get_energy_os_type, obj_func_oao, precond_pd_oao, &
                           project_oao, oao_object

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_os_type), intent(in), pointer :: get_energy_os
        procedure(update_dm_os_type), intent(in), pointer :: update_dm_os
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(precond_type), intent(out), pointer :: precond_arh_funptr
        procedure(precond_pd_type), intent(out), pointer :: precond_pd_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize error flag
        error = 0

        ! call common setup
        call arh_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        if (error /= 0) return

        ! set pointers to functions
        oao_object%get_energy_os => get_energy_os
        arh_object%update_dm_os => update_dm_os

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_oao
        update_orbs_arh_funptr => update_orbs_arh_os
        precond_arh_funptr => precond_arh
        precond_pd_arh_funptr => precond_pd_oao
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
        real(rp), allocatable :: fock_ao(:, :, :), v_nonlinear_ao(:, :, :), &
                                 dm_diff(:, :, :, :), fock_diff(:, :, :, :), &
                                 v_linear_diff(:, :, :, :), v_nonlinear_diff(:, :, :, :)

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

            ! the density was updated but the response was not rebuilt
            oao_object%response_stale = .true.

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

            ! the static Hessian part was just rebuilt, so any cached
            ! eigendecomposition of it is now stale
            oao_object%hess_eigen_stale = .true.

            ! prepare the density matrix difference
            n_list = size(arh_object%dm_list, 4)
            allocate(dm_diff(n_ao, n_ao, n_particle, n_list))
            do i = 1, n_list
                dm_diff(:, :, :, i) = arh_object%dm_list(:, :, :, i) - arh_object%dm_oao
            end do

            ! for MS-SR1 and cache the direction vectors
            if (arh_object%settings%arh_type == "ms_sr1") then
                ! prepare the linear and non-linear potential differences
                allocate(v_linear_diff(n_ao, n_ao, n_particle, n_list), &
                         v_nonlinear_diff(n_ao, n_ao, n_particle, n_list))
                do i = 1, n_list
                    v_nonlinear_diff(:, :, :, i) = &
                        arh_object%v_nonlinear_list(:, :, :, i) - &
                        arh_object%v_nonlinear_oao
                    v_linear_diff(:, :, :, i) = &
                        arh_object%fock_list(:, :, :, i) - arh_object%fock_oao - &
                        v_nonlinear_diff(:, :, :, i)
                end do

                ! get inverted A matrix
                ! linear part: this is exact since Coulomb and exact exchange are
                ! linear in the density matrix
                call get_ms_a_inv_cs(dm_diff, v_linear_diff, .true., arh_object%a_inv, &
                                     n_ao, arh_object%settings, error)
                if (error /= 0) return
                ! non-linear part: kept on a dedicated, separately regularized
                ! multisecant SR1 system, since the exact linear secant relationship
                ! and the non-linear one cannot be represented well by a single
                ! regularized eigenbasis
                call get_ms_a_inv_cs(dm_diff, v_nonlinear_diff, .false., &
                                     arh_object%a_inv_comb, n_ao, arh_object%settings, &
                                     error)
                if (error /= 0) return

                ! cache the packed history-projection directions the low-rank Hessian 
                ! factors are assembled from; kept in the original packed basis
                call cache_history_projections(v_linear_diff, arh_object%dm_oao, &
                                               n_list, arh_object%n_param, &
                                               arh_object%linear_potential_dirs)
                call cache_history_projections(v_nonlinear_diff, arh_object%dm_oao, &
                                               n_list, arh_object%n_param, &
                                               arh_object%nonlinear_potential_dirs)
            ! ARH and related methods
            else
                ! prepare the Fock matrix differences
                allocate(fock_diff(n_ao, n_ao, n_particle, n_list))
                do i = 1, n_list
                    fock_diff(:, :, :, i) = arh_object%fock_list(:, :, :, i) - &
                                            arh_object%fock_oao
                end do

                ! get inverse of metric
                call get_arh_metric_inv(dm_diff, arh_object%metric_inv)

                ! cache the packed history-projection directions the low-rank Hessian 
                ! factors are assembled from; kept in the original packed basis
                call cache_history_projections(dm_diff, arh_object%dm_oao, n_list, &
                                               arh_object%n_param, arh_object%dm_dirs)
                if (arh_object%settings%arh_type /= "ms_sp") &
                    call cache_history_projections(fock_diff, arh_object%dm_oao, &
                                                   n_list, arh_object%n_param, &
                                                   arh_object%potential_dirs)

                ! construct A = S^T Y
                if (arh_object%settings%arh_type == "ms_sp" .or. &
                    arh_object%settings%arh_type == "ms_psb") then
                    if (n_list > 0) then
                        arh_object%a_sym = build_a_sym_cs(dm_diff, fock_diff, n_ao)
                    else
                        if (allocated(arh_object%a_sym)) deallocate(arh_object%a_sym)
                        allocate(arh_object%a_sym(0, 0))
                    end if
                end if
            end if

            ! assemble the low-rank (response) part of the approximate Hessian
            call get_low_rank_hess_factors()
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
                                 v_nonlinear_ao(:, :, :), dm_diff(:, :, :, :), &
                                 v_same_spin_diff(:, :, :, :), &
                                 v_opposite_spin_diff(:, :, :, :), &
                                 v_nonlinear_diff(:, :, :, :), &
                                 v_same_spin_diff_eff(:, :, :, :)

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

            ! the density was updated but the response was not rebuilt
            oao_object%response_stale = .true.

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

            ! the static Hessian part was just rebuilt, so any cached
            ! eigendecomposition of it is now stale
            oao_object%hess_eigen_stale = .true.

            ! prepare the density matrix and opposite-spin potential differences
            n_list = size(arh_object%dm_list, 4)
            allocate(dm_diff(n_ao, n_ao, n_particle, n_list), &
                     v_opposite_spin_diff(n_ao, n_ao, n_particle, n_list))
            do i = 1, n_list
                dm_diff(:, :, :, i) = arh_object%dm_list(:, :, :, i) - arh_object%dm_oao
                v_opposite_spin_diff(:, :, :, i) = &
                    arh_object%v_opposite_spin_list(:, :, :, i) - &
                    arh_object%v_opposite_spin_oao
            end do

            ! MS-SR1
            if (arh_object%settings%arh_type == "ms_sr1") then
                ! prepare the same-spin and non-linear potential differences
                allocate(v_same_spin_diff(n_ao, n_ao, n_particle, n_list), &
                         v_nonlinear_diff(n_ao, n_ao, n_particle, n_list))
                do i = 1, n_list
                    v_same_spin_diff(:, :, :, i) = &
                        arh_object%v_same_spin_list(:, :, :, i) - &
                        arh_object%v_same_spin_oao
                    v_nonlinear_diff(:, :, :, i) = &
                        arh_object%v_nonlinear_list(:, :, :, i) - &
                        arh_object%v_nonlinear_oao
                end do

                ! get inverted A matrix
                ! linear part: get spin-separated multisecant SR1 matrix for which 
                ! separation is exact since Coulomb and exact exchange are linear in 
                ! the density matrix
                call get_ms_a_inv_os_linear(dm_diff, v_same_spin_diff, &
                                            v_opposite_spin_diff, arh_object%a_inv, &
                                            n_ao, arh_object%settings, error)
                if (error /= 0) return
                ! non-linear part: get spin-combined multisecant SR1 matrix
                call get_ms_a_inv_os_nonlinear(dm_diff, v_nonlinear_diff, &
                                               arh_object%a_inv_comb, &
                                               arh_object%settings, error)
                if (error /= 0) return

                ! cache the packed history-projection directions the low-rank Hessian 
                ! factors are assembled from; the linear potential directions combine 
                ! the same-/opposite-spin channels, while the non-linear potential 
                ! directions need no channel-splitting
                call cache_combined_channel_projections( &
                    v_same_spin_diff, v_opposite_spin_diff, arh_object%dm_oao, &
                    n_list, arh_object%n_param, n_particle, &
                    arh_object%linear_potential_dirs)
                call cache_history_projections(v_nonlinear_diff, arh_object%dm_oao, &
                                               n_list, arh_object%n_param, &
                                               arh_object%nonlinear_potential_dirs)
            ! ARH and related methods
            else
                ! fold in non-linear potantial into same-spin potential
                allocate(v_same_spin_diff_eff(n_ao, n_ao, n_particle, n_list))
                do i = 1, n_list
                    v_same_spin_diff_eff(:, :, :, i) = &
                        arh_object%v_same_spin_list(:, :, :, i) - &
                        arh_object%v_same_spin_oao + &
                        arh_object%v_nonlinear_list(:, :, :, i) - &
                        arh_object%v_nonlinear_oao
                end do

                ! get inverse of metric which is block-diagonal accross channels
                call get_arh_metric_inv(dm_diff, arh_object%metric_inv)

                ! cache the packed history-projection directions the low-rank Hessian 
                ! factors are assembled from; the density matrix directions isolate
                ! each channel separately, the effective potential directions combine 
                ! channels
                call cache_channel_split_projections( &
                    dm_diff, arh_object%dm_oao, n_list, &
                    arh_object%n_param, n_particle, arh_object%dm_dirs)
                if (arh_object%settings%arh_type /= "ms_sp") &
                    call cache_combined_channel_projections( &
                        v_same_spin_diff_eff, v_opposite_spin_diff, &
                        arh_object%dm_oao, n_list, arh_object%n_param, n_particle, &
                        arh_object%potential_dirs)

                ! construct A = S^T Y which is cross-channel symmetrized
                if (arh_object%settings%arh_type == "ms_sp" .or. &
                    arh_object%settings%arh_type == "ms_psb") then
                    if (n_list > 0) then
                        arh_object%a_sym = build_a_block_sym_os( &
                            dm_diff, v_same_spin_diff_eff, v_opposite_spin_diff, n_ao)
                    else
                        if (allocated(arh_object%a_sym)) deallocate(arh_object%a_sym)
                        allocate(arh_object%a_sym(0, 0))
                    end if
                end if
            end if

            ! assemble the low-rank (response) part of the approximate Hessian
            call get_low_rank_hess_factors()
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
        use otr_oao, only: unpack_asymm, project_asymm, pack_asymm, &
                           symmetric_transformation

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, n_dirs, i
        real(rp), allocatable :: x_full(:, :, :), hess_x_full(:, :, :), &
                                 projected_x(:), coupled_x(:)

        external :: dgemm, dgemv

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

        ! project, scale and pack the static part; the static part is already
        ! confined to the occupied-virtual and virtual-occupied subspace in exact
        ! arithmetic, but projecting anyway prevents numerical leakage into the
        ! redundant subspace
        hess_x_full = project_asymm(hess_x_full, arh_object%dm_oao)
        if (n_particle == 1) then
            hess_x = 4.0_rp * pack_asymm(hess_x_full, size(hess_x, kind=ip))
        else
            hess_x = 2.0_rp * pack_asymm(hess_x_full, size(hess_x, kind=ip))
        end if
        deallocate(x_full, hess_x_full)

        ! add the response part, which the shared low-rank factors express directly
        ! in the packed parameter space as expansion * coupling * projection^T,
        ! already carrying both the projection onto the non-redundant subspace and
        ! the scaling applied to the static part above
        if (allocated(arh_object%coupling_matrix)) then
            n_dirs = size(arh_object%expansion_dirs, 2)
            allocate(projected_x(n_dirs), coupled_x(n_dirs))
            call dgemv("T", size(x, kind=ip), n_dirs, 1.0_rp, &
                       arh_object%projection_dirs, size(x, kind=ip), x, 1_ip, 0.0_rp, &
                       projected_x, 1_ip)
            call dgemv("N", n_dirs, n_dirs, 1.0_rp, &
                       arh_object%coupling_matrix, n_dirs, projected_x, 1_ip, 0.0_rp, &
                       coupled_x, 1_ip)
            call dgemv("N", size(hess_x, kind=ip), n_dirs, 1.0_rp, &
                       arh_object%expansion_dirs, size(hess_x, kind=ip), coupled_x, &
                       1_ip, 1.0_rp, hess_x, 1_ip)
            deallocate(projected_x, coupled_x)
        end if

    end subroutine hess_x_arh

    subroutine inv_hess_x_arh(x, inv_hess_x, error, level_shift)
        !
        ! this subroutine applies the exact inverse of the approximate Hessian to a
        ! vector, optionally level-shifted as (G - level_shift*I)^-1; G is the static 
        ! part D plus the low-rank correction which is inverted using the 
        ! Sherman-Morrison-Woodbury identity written as
        !
        !     (D + E C P^T)^-1 = D^-1 - D^-1 E (I + C P^T D^-1 E)^-1 C P^T D^-1
        !
        ! with E the expansion directions, P the projection directions and C the
        ! coupling matrix
        !
        use opentrustregion, only: precond_floor, verbosity_error
        use otr_oao, only: precond_oao, rotate_to_hess_eigenbasis, &
                           rotate_from_hess_eigenbasis, get_hess_eigval_pairs, &
                           refresh_hess_eigen, oao_object

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: inv_hess_x(:)
        integer(ip), intent(out) :: error
        real(rp), intent(in), optional :: level_shift

        integer(ip) :: n_param, n_dirs, i, info
        real(rp) :: mu
        character(300) :: msg
        real(rp), allocatable :: rotated_x(:), eigval_pairs(:), scaled_x(:), &
                                 rotated_expansion(:, :), rotated_projection(:, :), &
                                 weighted_projection(:, :), dirs_overlap(:, :), &
                                 bracket_matrix(:, :), projected_x(:), &
                                 bracket_rhs(:), bracket_solution(:), correction(:)
        integer(ip), allocatable :: ipiv(:)
        external :: dgemv, dgemm, dgesv

        ! initialize error flag
        error = 0

        ! an absent level shift gives the plain inverse of the approximate Hessian
        mu = 0.0_rp
        if (present(level_shift)) mu = level_shift

        ! fall back to the static part alone if there is no low-rank correction
        if (.not. allocated(arh_object%coupling_matrix)) then
            call precond_oao(x, mu, inv_hess_x, error)
            return
        end if

        ! refresh the eigendecomposition if the static Hessian part has changed
        if (oao_object%hess_eigen_stale) then
            call refresh_hess_eigen(error)
            if (error /= 0) return
        end if

        ! get parameters
        n_param = size(x)
        n_dirs = size(arh_object%expansion_dirs, 2)

        ! rotate x into the static-Hessian eigenbasis and apply the level-shifted 
        ! diagonal D^-1
        rotated_x = rotate_to_hess_eigenbasis(x)
        eigval_pairs = get_hess_eigval_pairs() - mu
        where (abs(eigval_pairs) < precond_floor) eigval_pairs = precond_floor
        scaled_x = rotated_x / eigval_pairs

        ! rotate both sets of directions into the same eigenbasis
        allocate(rotated_expansion(n_param, n_dirs), &
                 rotated_projection(n_param, n_dirs))
        do i = 1, n_dirs
            rotated_expansion(:, i) = &
                rotate_to_hess_eigenbasis(arh_object%expansion_dirs(:, i))
            rotated_projection(:, i) = &
                rotate_to_hess_eigenbasis(arh_object%projection_dirs(:, i))
        end do

        ! projected_x = P^T D^-1 x
        allocate(projected_x(n_dirs))
        call dgemv("T", n_param, n_dirs, 1.0_rp, rotated_projection, n_param, &
                   scaled_x, 1_ip, 0.0_rp, projected_x, 1_ip)

        ! dirs_overlap = P^T D^-1 E
        allocate(weighted_projection(n_param, n_dirs))
        do i = 1, n_param
            weighted_projection(i, :) = rotated_projection(i, :) / eigval_pairs(i)
        end do
        allocate(dirs_overlap(n_dirs, n_dirs))
        call dgemm("T", "N", n_dirs, n_dirs, n_param, 1.0_rp, weighted_projection, &
                   n_param, rotated_expansion, n_param, 0.0_rp, dirs_overlap, n_dirs)
        deallocate(weighted_projection)

        ! solve (I + C * P^T D^-1 E) bracket_solution = C * projected_x
        allocate(bracket_matrix(n_dirs, n_dirs))
        call dgemm("N", "N", n_dirs, n_dirs, n_dirs, 1.0_rp, &
                   arh_object%coupling_matrix, n_dirs, dirs_overlap, n_dirs, 0.0_rp, &
                   bracket_matrix, n_dirs)
        do i = 1, n_dirs
            bracket_matrix(i, i) = bracket_matrix(i, i) + 1.0_rp
        end do
        allocate(bracket_rhs(n_dirs))
        call dgemv("N", n_dirs, n_dirs, 1.0_rp, arh_object%coupling_matrix, n_dirs, &
                   projected_x, 1_ip, 0.0_rp, bracket_rhs, 1_ip)
        allocate(bracket_solution(n_dirs), ipiv(n_dirs))
        bracket_solution = bracket_rhs
        call dgesv(n_dirs, 1_ip, bracket_matrix, n_dirs, ipiv, bracket_solution, &
                   n_dirs, info)
        if (info /= 0) then
            write (msg, '(A, I0)') "Level-shifted approximate Hessian is singular: "// &
                "Error in DGESV, info = ", info
            call arh_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! correction = D^-1 E bracket_solution, result = D^-1 x - correction
        allocate(correction(n_param))
        call dgemv("N", n_param, n_dirs, 1.0_rp, rotated_expansion, n_param, &
                   bracket_solution, 1_ip, 0.0_rp, correction, 1_ip)
        correction = correction / eigval_pairs
        inv_hess_x = rotate_from_hess_eigenbasis(scaled_x - correction)

    end subroutine inv_hess_x_arh

    subroutine precond_arh(residual, mu, precond_residual, error)
        !
        ! this subroutine defines the preconditioner of the ARH approximate Hessian
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        call inv_hess_x_arh(residual, precond_residual, error, mu)

    end subroutine precond_arh

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

    subroutine cache_history_projections(v_diff, dm_oao, n_list, n_param, dirs)
        !
        ! this subroutine caches the packed history-projection directions the
        ! low-rank part of the approximate Hessian is built from: each history entry
        ! is projected onto the occupied-virtual/virtual-occupied subspace and packed
        ! into the same antisymmetric parameter space as a trial vector
        !
        use otr_oao, only: project_asymm, pack_asymm

        real(rp), intent(in) :: v_diff(:, :, :, :), dm_oao(:, :, :)
        integer(ip), intent(in) :: n_list, n_param
        real(rp), intent(out), allocatable :: dirs(:, :)

        integer(ip) :: k
        real(rp), allocatable :: projected(:, :, :)

        allocate(dirs(n_param, n_list))
        do k = 1, n_list
            projected = project_asymm(v_diff(:, :, :, k), dm_oao)
            dirs(:, k) = pack_asymm(projected, n_param)
        end do

    end subroutine cache_history_projections

    subroutine cache_history_projections_channel(v_diff, channel, dm_oao, n_list, &
                                                 n_param, n_particle, dirs)
        !
        ! this subroutine caches the packed history-projection directions the low-rank 
        ! part of the approximate Hessian is built from for open-shell systems; each 
        ! history column is embedded into only the given channel (the other channels 
        ! zeroed) before projecting and packing, since the open-shell response 
        ! contracts a channel's own difference against only that channel's density 
        ! response
        !
        use otr_oao, only: project_asymm, pack_asymm

        real(rp), intent(in) :: v_diff(:, :, :, :), dm_oao(:, :, :)
        integer(ip), intent(in) :: channel, n_list, n_param, n_particle
        real(rp), intent(out), allocatable :: dirs(:, :)

        integer(ip) :: k, n_ao
        real(rp), allocatable :: embedded(:, :, :), projected(:, :, :)

        n_ao = size(dm_oao, 1)
        allocate(dirs(n_param, n_list), embedded(n_ao, n_ao, n_particle))
        embedded = 0.0_rp
        do k = 1, n_list
            embedded(:, :, channel) = v_diff(:, :, channel, k)
            projected = project_asymm(embedded, dm_oao)
            dirs(:, k) = pack_asymm(projected, n_param)
        end do
        deallocate(embedded)

    end subroutine cache_history_projections_channel

    subroutine cache_channel_split_projections(v_diff, dm_oao, n_list, n_param, &
                                               n_particle, dirs)
        !
        ! this subroutine caches the open-shell history projections one spin channel
        ! at a time: column i holds channel 1 of entry i, column n_list+i holds
        ! channel 2; the two are never summed, since the per-channel metric
        ! contraction does not mix channels
        !
        real(rp), intent(in) :: v_diff(:, :, :, :), dm_oao(:, :, :)
        integer(ip), intent(in) :: n_list, n_param, n_particle
        real(rp), intent(out), allocatable :: dirs(:, :)

        real(rp), allocatable :: channel1(:, :), channel2(:, :)

        call cache_history_projections_channel(v_diff, 1_ip, dm_oao, n_list, n_param, &
                                               n_particle, channel1)
        call cache_history_projections_channel(v_diff, 2_ip, dm_oao, n_list, n_param, &
                                               n_particle, channel2)

        allocate(dirs(n_param, 2 * n_list))
        dirs(:, 1:n_list) = channel1
        dirs(:, n_list + 1:2 * n_list) = channel2
        deallocate(channel1, channel2)

    end subroutine cache_channel_split_projections

    subroutine cache_combined_channel_projections(v_same, v_opp, dm_oao, n_list, &
                                                  n_param, n_particle, dirs)
        !
        ! this subroutine caches the open-shell history projections two channels at a
        ! time: column i sums channel 1 of the same-spin potential with channel 2 of 
        ! the opposite-spin potential, column n_list+i the mirror image; every caller 
        ! pairs a same-spin difference with the opposite-spin one, differing only in 
        ! what is passed as the same-spin potential
        !
        real(rp), intent(in) :: v_same(:, :, :, :), v_opp(:, :, :, :), dm_oao(:, :, :)
        integer(ip), intent(in) :: n_list, n_param, n_particle
        real(rp), intent(out), allocatable :: dirs(:, :)

        real(rp), allocatable :: same1(:, :), same2(:, :), opp1(:, :), opp2(:, :)

        call cache_history_projections_channel(v_same, 1_ip, dm_oao, n_list, n_param, &
                                               n_particle, same1)
        call cache_history_projections_channel(v_same, 2_ip, dm_oao, n_list, n_param, &
                                               n_particle, same2)
        call cache_history_projections_channel(v_opp, 1_ip, dm_oao, n_list, n_param, &
                                               n_particle, opp1)
        call cache_history_projections_channel(v_opp, 2_ip, dm_oao, n_list, n_param, &
                                               n_particle, opp2)

        allocate(dirs(n_param, 2 * n_list))
        dirs(:, 1:n_list) = same1 + opp2
        dirs(:, n_list + 1:2 * n_list) = same2 + opp1
        deallocate(same1, same2, opp1, opp2)

    end subroutine cache_combined_channel_projections

    subroutine get_low_rank_hess_factors()
        !
        ! this subroutine assembles the low-rank part of the approximate Hessian in
        ! the packed parameter space, as
        !
        !     G_low_rank = expansion_dirs * coupling_matrix * transpose(projection_dirs)
        !
        ! by constructing expansion_dirs, coupling_matrix, and projection_dirs
        
        ! the response is defined through Frobenius inner products
        ! <history_k, delta_dm(x)> of a history matrix with the density response, yet
        ! the whole correction can be expressed on packed parameter vectors alone:
        ! project_asymm (on symmetric input) and project_symm (on antisymmetric
        ! input) are adjoint with respect to that inner product, so projecting and
        ! packing a history matrix into dirs(:, k) turns the contraction into twice
        ! the plain dot product <dirs(:, k), x>; packing and projecting are linear,
        ! so the response's output expansion collapses the same way
        !
        ! every closed-shell coupling matrix therefore carries an overall factor of
        ! 8: the factor of 2 above, times a factor of 4 (or 2) for the closed- (or
        ! open-)shell scaling of the Hessian linear transformation
        !
        integer(ip) :: n_linear, n_nonlinear, n_total, n_diff
        real(rp) :: shell_scale
        real(rp), allocatable :: metric_weighted_a_sym(:, :)

        ! discard the factors assembled for the previous history
        if (allocated(arh_object%expansion_dirs)) deallocate(arh_object%expansion_dirs)
        if (allocated(arh_object%projection_dirs)) &
            deallocate(arh_object%projection_dirs)
        if (allocated(arh_object%coupling_matrix)) &
            deallocate(arh_object%coupling_matrix)

        ! set scaling factor for closed- and open-shell systems
        shell_scale = merge(1.0_rp, 0.5_rp, arh_object%n_particle == 1)

        select case (arh_object%settings%arh_type)
        ! multisecant SR1: the linear and non-linear parts of the potential enter as
        ! independent blocks, each with its own separately regularized system
        case ("ms_sr1")
            if (.not. allocated(arh_object%linear_potential_dirs) .or. &
                .not. allocated(arh_object%nonlinear_potential_dirs)) return
            n_linear = size(arh_object%linear_potential_dirs, 2)
            n_nonlinear = size(arh_object%nonlinear_potential_dirs, 2)
            n_total = n_linear + n_nonlinear
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, n_linear + 1:) = &
                arh_object%nonlinear_potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_linear, :n_linear) = 8.0_rp * shell_scale * &
                                                               arh_object%a_inv
            arh_object%coupling_matrix(n_linear + 1:, n_linear + 1:) = &
                8.0_rp * shell_scale * arh_object%a_inv_comb

        ! subspace-projected multisecant: the response both expands in and contracts
        ! against the density difference history alone
        case ("ms_sp")
            if (.not. allocated(arh_object%dm_dirs)) return
            n_diff = size(arh_object%dm_dirs, 2)
            if (n_diff == 0) return

            metric_weighted_a_sym = matmul( &
                arh_object%metric_inv, matmul(arh_object%a_sym, arh_object%metric_inv))
            arh_object%expansion_dirs = arh_object%dm_dirs
            arh_object%coupling_matrix = 8.0_rp * shell_scale * metric_weighted_a_sym

        ! symmetrized ARH: the density and potential difference histories couple in
        ! both directions, so the coupling matrix is purely off-diagonal; also 
        ! introduces an additional factor of 1/2
        case ("symm_arh")
            if (.not. allocated(arh_object%dm_dirs)) return
            n_diff = size(arh_object%dm_dirs, 2)
            if (n_diff == 0) return

            allocate(arh_object%expansion_dirs(arh_object%n_param, 2 * n_diff), &
                     arh_object%coupling_matrix(2 * n_diff, 2 * n_diff))
            arh_object%expansion_dirs(:, :n_diff) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_diff + 1:) = arh_object%potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_diff, n_diff + 1:) = &
                4.0_rp * shell_scale * arh_object%metric_inv
            arh_object%coupling_matrix(n_diff + 1:, :n_diff) = &
                4.0_rp * shell_scale * arh_object%metric_inv

        ! multisecant PSB: the symmetrized ARH coupling with an additional
        ! density-density block subtracting the doubly counted curvature
        case ("ms_psb")
            if (.not. allocated(arh_object%dm_dirs)) return
            n_diff = size(arh_object%dm_dirs, 2)
            if (n_diff == 0) return

            metric_weighted_a_sym = matmul( &
                arh_object%metric_inv, matmul(arh_object%a_sym, arh_object%metric_inv))
            allocate(arh_object%expansion_dirs(arh_object%n_param, 2 * n_diff), &
                     arh_object%coupling_matrix(2 * n_diff, 2 * n_diff))
            arh_object%expansion_dirs(:, :n_diff) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_diff + 1:) = arh_object%potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_diff, :n_diff) = &
                -8.0_rp * shell_scale * metric_weighted_a_sym
            arh_object%coupling_matrix(:n_diff, n_diff + 1:) = &
                8.0_rp * shell_scale * arh_object%metric_inv
            arh_object%coupling_matrix(n_diff + 1:, :n_diff) = &
                8.0_rp * shell_scale * arh_object%metric_inv

        ! standard ARH: the response expands in the potential difference history
        ! while contracting against the density difference history, so unlike every
        ! other type the two sets of directions differ
        case ("arh")
            if (.not. allocated(arh_object%dm_dirs)) return
            n_diff = size(arh_object%dm_dirs, 2)
            if (n_diff == 0) return

            arh_object%expansion_dirs = arh_object%potential_dirs
            arh_object%projection_dirs = arh_object%dm_dirs
            arh_object%coupling_matrix = 8.0_rp * shell_scale * arh_object%metric_inv

        case default
            return
        end select

        ! every type except standard ARH expands in and contracts against the same
        ! set of directions
        if (allocated(arh_object%expansion_dirs) .and. &
            .not. allocated(arh_object%projection_dirs)) &
            arh_object%projection_dirs = arh_object%expansion_dirs

    end subroutine get_low_rank_hess_factors

    function build_a_sym_cs(dm_diff, v_diff, n_ao) result(a_sym)
        !
        ! this function builds the dense, weighted-symmetrized A = S^T Y matrix shared 
        ! by the MS-SP and MS_PSB response contributions for the closed-shell case
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :), v_diff(:, :, :, :)
        integer(ip), intent(in) :: n_ao
        real(rp), allocatable :: a_sym(:, :)

        integer(ip) :: n_diff, i
        real(rp), allocatable :: step_norms(:)
        external :: dgemm
        real(rp), external :: dnrm2

        n_diff = size(dm_diff, 4)
        allocate(a_sym(n_diff, n_diff), step_norms(n_diff))

        ! build A = S^T Y
        call dgemm("T", "N", n_diff, n_diff, n_ao * n_ao, 1.0_rp, dm_diff, &
                   n_ao * n_ao, v_diff, n_ao * n_ao, 0.0_rp, a_sym, n_diff)

        ! symmetrize
        do i = 1, n_diff
            step_norms(i) = dnrm2(n_ao * n_ao, dm_diff(:, :, 1, i), 1_ip)
        end do
        call symmetrize_weighted(a_sym, step_norms)
        deallocate(step_norms)

    end function build_a_sym_cs

    function build_a_block_sym_os(dm_diff, v_same_eff, v_opp, n_ao) result(a_block)
        !
        ! this function builds the dense, cross-channel- and weighted-symmetrized 
        ! open-shell A = S^T Y matrix entering the MS-SP and MS-PSB response 
        ! contributions for the open-shell case: each diagonal block is a 
        ! within-channel weighted symmetrization, while the two off-diagonal blocks are 
        ! cross-symmetrized together
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :), v_same_eff(:, :, :, :), &
                                v_opp(:, :, :, :)
        integer(ip), intent(in) :: n_ao
        real(rp), allocatable :: a_block(:, :)

        integer(ip) :: n_diff, i, j
        real(rp), allocatable :: step_norms(:, :), a_same(:, :, :), a_opp(:, :, :), &
                                 dm_diff_j(:, :, :), v_same_eff_j(:, :, :), &
                                 v_opp_j(:, :, :)
        external :: dgemm
        real(rp), external :: dnrm2

        n_diff = size(dm_diff, 4)
        allocate(step_norms(n_diff, 2), a_same(n_diff, n_diff, 2), &
                 a_opp(n_diff, n_diff, 2))

        do j = 1, 2
            do i = 1, n_diff
                step_norms(i, j) = dnrm2(n_ao * n_ao, dm_diff(:, :, j, i), 1_ip)
            end do
        end do

        ! build A = S^T Y and symmetrize diagonal blocks
        do j = 1, 2
            dm_diff_j = dm_diff(:, :, j, :)
            v_same_eff_j = v_same_eff(:, :, j, :)
            v_opp_j = v_opp(:, :, j, :)
            call dgemm("T", "N", n_diff, n_diff, n_ao * n_ao, 1.0_rp, dm_diff_j, &
                       n_ao * n_ao, v_same_eff_j, n_ao * n_ao, 0.0_rp, &
                       a_same(:, :, j), n_diff)
            call symmetrize_weighted(a_same(:, :, j), step_norms(:, j))
            call dgemm("T", "N", n_diff, n_diff, n_ao * n_ao, 1.0_rp, dm_diff_j, &
                       n_ao * n_ao, v_opp_j, n_ao * n_ao, 0.0_rp, a_opp(:, :, j), &
                       n_diff)
        end do
        deallocate(dm_diff_j, v_same_eff_j, v_opp_j)

        ! cross-symmetrize off-diagonal blocks
        call cross_symmetrize_weighted(a_opp(:, :, 1), a_opp(:, :, 2), &
                                       step_norms(:, 1), step_norms(:, 2))
        deallocate(step_norms)

        allocate(a_block(2 * n_diff, 2 * n_diff))
        a_block(1:n_diff, 1:n_diff) = a_same(:, :, 1)
        a_block(1:n_diff, n_diff + 1:2 * n_diff) = a_opp(:, :, 1)
        a_block(n_diff + 1:2 * n_diff, 1:n_diff) = a_opp(:, :, 2)
        a_block(n_diff + 1:2 * n_diff, n_diff + 1:2 * n_diff) = a_same(:, :, 2)
        deallocate(a_same, a_opp)

    end function build_a_block_sym_os

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

    subroutine get_arh_metric_inv(dm_diff, metric_inv)
        !
        ! this subroutine calculates the pseudoinverse of the augmented Roothaan-Hall
        ! metric, block-diagonal across particle channels since the channels never mix
        ! when the metric itself is inverted (only the response's input and output
        ! directions do); the metric is factorized via an unpivoted Cholesky
        ! decomposition and, since the density matrix differences are saved in reverse
        ! order, linearly dependent older vectors are skipped while preserving the
        ! chronological sequence for the active basis; the accepted block is then
        ! inverted directly from its Cholesky factor and scattered back to the original
        ! history ordering, so that rejected directions carry zero rows and columns
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: dm_diff(:, :, :, :)
        real(rp), intent(out), allocatable :: metric_inv(:, :)

        integer(ip) :: n_ao, n_particle, n_dm, i, j, k, offset, info
        integer(ip) :: n_accepted, n_rejected
        integer(ip), allocatable :: map(:)
        real(rp) :: raw_diagonal, tol
        real(rp), allocatable :: metric(:, :), chol(:, :), work(:)
        external :: dtrsv, dpotri
        real(rp), external :: ddot

        n_ao = size(dm_diff, 1)
        n_particle = size(dm_diff, 3)
        n_dm = size(dm_diff, 4)

        ! allocate the block-diagonal pseudoinverse and handle an empty history
        allocate(metric_inv(n_particle * n_dm, n_particle * n_dm))
        metric_inv = 0.0_rp
        if (n_dm == 0) return

        allocate(metric(n_dm, n_dm), chol(n_dm, n_dm), map(n_dm), work(n_dm))

        do k = 1, n_particle
            ! initialize tolerance with maximum diagonal element
            chol = 0.0_rp
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
                        work(j) = metric(map(j), i)
                    end do

                    ! solve R_accepted^T * work = T_accepted_vs_current
                    call dtrsv("U", "T", "N", n_accepted, chol, n_dm, work, 1_ip)

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
                    map(n_dm - n_rejected + 1) = i
                else
                    ! independent: accept column and append to the active factors
                    n_accepted = n_accepted + 1
                    map(n_accepted) = i

                    if (n_accepted > 1) then
                        chol(1:n_accepted-1, n_accepted) = work(1:n_accepted-1)
                    end if
                    chol(n_accepted, n_accepted) = sqrt(max(0.0_rp, raw_diagonal))
                end if
            end do
            if (n_accepted == 0) cycle

            ! invert the accepted block in place from its Cholesky factor; a failed
            ! inversion means the block is singular despite the dependency screening,
            ! in which case the channel simply contributes nothing
            call dpotri("U", n_accepted, chol, n_dm, info)
            if (info /= 0) cycle

            ! scatter the inverted upper triangle back to the original history
            ! ordering, leaving the rejected directions as zero rows and columns
            offset = (k - 1) * n_dm
            do j = 1, n_accepted
                do i = 1, j
                    metric_inv(offset + map(i), offset + map(j)) = chol(i, j)
                    metric_inv(offset + map(j), offset + map(i)) = chol(i, j)
                end do
            end do
        end do
        deallocate(metric, chol, map, work)

    end subroutine get_arh_metric_inv

    subroutine get_ms_a_inv_cs(dm_diff, v_diff, linear, a_inv, n_ao, settings, error)
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
        real(rp), intent(out), allocatable :: a_inv(:, :)
        integer(ip), intent(in) :: n_ao
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, i
        real(rp), allocatable :: a(:, :), eig_vecs(:, :), eig_vals(:), &
                                 eig_vals_inv(:), step_norms(:)
        real(rp) :: eig_val_thresh
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! handle empty history
        n_dm = size(dm_diff, 4)
        if (n_dm == 0) then
            allocate(a_inv(0, 0))
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

        ! reassemble the pseudoinverse
        a_inv = spectral_to_dense(eig_vecs, eig_vals_inv)
        deallocate(a, eig_vecs, eig_vals, eig_vals_inv)

    end subroutine get_ms_a_inv_cs

    subroutine get_ms_a_inv_os_linear(dm_diff, v_same_spin_diff, v_opposite_spin_diff, &
                                      a_inv, n_ao, settings, error)
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
        real(rp), intent(out), allocatable :: a_inv(:, :)
        integer(ip), intent(in) :: n_ao
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, i, k
        real(rp), allocatable :: a(:, :), eig_vecs(:, :), eig_vals(:), &
                                 eig_vals_inv(:)
        real(rp) :: eig_val_thresh

        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! handle empty history
        n_dm = size(dm_diff, 4)
        if (n_dm == 0) then
            allocate(a_inv(0, 0))
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

        ! reassemble the pseudoinverse
        a_inv = spectral_to_dense(eig_vecs, eig_vals_inv)
        deallocate(a, eig_vecs, eig_vals, eig_vals_inv)

    end subroutine get_ms_a_inv_os_linear

    subroutine get_ms_a_inv_os_nonlinear(dm_diff, v_diff, a_inv, settings, error)
        !
        ! this subroutine computes the pseudoinverse multisecant SR1 matrix in a
        ! spin-combined manner for the non-linear part in the open-shell case; a 
        ! weighted symmetrization is used with a relative threshold for the 
        ! regularization
        !
        use opentrustregion, only: symm_mat_diag, numerical_zero

        real(rp), intent(in) :: dm_diff(:, :, :, :), v_diff(:, :, :, :)
        real(rp), intent(out), allocatable :: a_inv(:, :)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, flat_len, i, j
        real(rp), allocatable :: a(:, :), eig_vecs(:, :), eig_vals(:), &
                                 eig_vals_inv(:), step_norms(:)
        real(rp) :: eig_val_thresh
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        ! handle empty history
        n_dm = size(dm_diff, 4)
        if (n_dm == 0) then
            allocate(a_inv(0, 0))
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

        ! reassemble the pseudoinverse
        a_inv = spectral_to_dense(eig_vecs, eig_vals_inv)
        deallocate(a, eig_vecs, eig_vals, eig_vals_inv)

    end subroutine get_ms_a_inv_os_nonlinear

    function spectral_to_dense(eigvecs, inv_eigvals) result(mat)
        !
        ! this function reconstructs the dense symmetric matrix V * diag(lambda) * V^T
        ! from its eigenvectors and eigenvalues
        !
        real(rp), intent(in) :: eigvecs(:, :), inv_eigvals(:)
        real(rp), allocatable :: mat(:, :)

        integer(ip) :: n, i
        real(rp), allocatable :: scaled(:, :)
        external :: dgemm

        n = size(eigvecs, 1)
        allocate(scaled(n, n), mat(n, n))
        scaled = eigvecs
        do i = 1, n
            scaled(:, i) = scaled(:, i) * inv_eigvals(i)
        end do
        call dgemm("N", "T", n, n, n, 1.0_rp, scaled, n, eigvecs, n, 0.0_rp, mat, n)
        deallocate(scaled)

    end function spectral_to_dense

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
