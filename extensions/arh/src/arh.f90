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
    real(rp), parameter :: history_step_ratio = 1e2_rp, ms_sr1_skip_thresh = 0.05_rp, &
                           eig_val_noise_factor = 10.0_rp

    type, extends(oao_settings_type) :: arh_settings_type
        character(kw_len) :: arh_type
    contains
        procedure :: init => init_arh_settings
    end type arh_settings_type

    type(arh_settings_type), parameter :: default_arh_settings = &
        arh_settings_type(oao_settings_type = default_oao_settings, arh_type = "ms_sr1")

    ! define setting options
    character(kw_len), parameter :: arh_types(5) = &
        [character(len=kw_len) :: "arh", "symm_arh", "ms_psb", "ms_sp", "ms_sr1"]

    abstract interface
        subroutine evaluate_dm_cs_type(dm, energy, fock, v_nonlinear, error)
            import :: rp, ip

            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), optional, target, contiguous :: fock(:, :), &
                                                                   v_nonlinear(:, :)
            integer(ip), intent(out) :: error
        end subroutine evaluate_dm_cs_type

        subroutine evaluate_dm_os_type(dm, energy, fock, v_same_spin, v_opposite_spin, &
                                       v_nonlinear, error)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), optional, target :: &
                fock(:, :, :), v_same_spin(:, :, :), v_opposite_spin(:, :, :), &
                v_nonlinear(:, :, :)
            integer(ip), intent(out) :: error
        end subroutine evaluate_dm_os_type
    end interface

    type :: arh_type
        type(arh_settings_type) :: settings
        integer(ip), pointer :: n_ao => null(), n_param => null(), n_particle => null()
        real(rp), pointer, contiguous :: dm_ao(:, :, :) => null()
        real(rp), pointer :: s_inv_sqrt(:, :) => null(), dm_oao(:, :, :) => null(), &
                             fock_oo(:, :, :) => null(), fock_vv(:, :, :) => null(), &
                             energy => null(), grad(:) => null(), h_diag(:) => null()
        real(rp), allocatable :: &
            fock_oao(:, :, :), v_same_spin_oao(:, :, :), v_opposite_spin_oao(:, :, :), &
            v_nonlinear_oao(:, :, :), a_sym(:, :), a_sym_nonlinear(:, :), a_inv(:, :), &
            a_inv_comb(:, :), dm_list(:, :, :, :), fock_list(:, :, :, :), &
            v_same_spin_list(:, :, :, :), v_opposite_spin_list(:, :, :, :), &
            v_nonlinear_list(:, :, :, :), linear_potential_dirs(:, :), &
            nonlinear_potential_dirs(:, :), dm_dirs(:, :), dm_dirs_nonlinear(:, :), &
            expansion_dirs(:, :), projection_dirs(:, :), coupling_matrix(:, :)
        procedure(evaluate_dm_os_type), pointer, nopass :: evaluate_dm_os => null()
        procedure(evaluate_dm_cs_type), pointer, nopass :: evaluate_dm_cs => null()
    end type arh_type

    ! global variables
    type(arh_type), allocatable :: arh_object

    ! create function pointers to ensure that routines comply with interface
    procedure(obj_func_type), pointer :: obj_func_arh_cs_ptr => obj_func_arh_cs
    procedure(obj_func_type), pointer :: obj_func_arh_os_ptr => obj_func_arh_os
    procedure(update_orbs_type), pointer :: update_orbs_arh_cs_ptr => update_orbs_arh_cs
    procedure(update_orbs_type), pointer :: update_orbs_arh_os_ptr => update_orbs_arh_os
    procedure(hess_x_type), pointer :: hess_x_arh_ptr => hess_x_arh
    procedure(precond_type), pointer :: precond_arh_ptr => precond_arh

    ! define module procedures for different spin cases
    interface arh_factory
        module procedure arh_factory_cs, arh_factory_os
    end interface arh_factory

contains

    subroutine arh_factory_cs(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_cs, &
                              obj_func_arh_funptr, update_orbs_arh_funptr, &
                              precond_arh_funptr, precond_pd_arh_funptr, &
                              project_arh_funptr, error, settings)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! closed-shell case
        !
        use otr_oao, only: precond_pd_oao, project_oao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_cs
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
        arh_object%evaluate_dm_cs => evaluate_dm_cs

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh_cs
        update_orbs_arh_funptr => update_orbs_arh_cs
        precond_arh_funptr => precond_arh
        precond_pd_arh_funptr => precond_pd_oao
        project_arh_funptr => project_oao

    end subroutine arh_factory_cs

    subroutine arh_factory_os(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_os, &
                              obj_func_arh_funptr, update_orbs_arh_funptr, &
                              precond_arh_funptr, precond_pd_arh_funptr, &
                              project_arh_funptr, error, settings)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! open-shell case
        !
        use otr_oao, only: precond_pd_oao, project_oao

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os
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
        arh_object%evaluate_dm_os => evaluate_dm_os

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh_os
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

    function obj_func_arh_cs(kappa, error) result(energy)
        !
        ! this function defines the energy evaluation in the OAO basis for the
        ! closed-shell case, which also adds the evaluated point with its Fock matrix
        ! and non-linear potential to the history, unless it is already there
        !
        use otr_oao, only: rotate_dm_ao, symmetric_transformation

        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        integer(ip) :: n_ao
        real(rp), allocatable :: rot_dm_ao(:, :, :), rot_dm_oao(:, :, :), &
                                 fock_ao(:, :, :), v_nonlinear_ao(:, :, :)

        ! initialize energy in case of error
        energy = 0.0_rp

        ! number of AOs
        n_ao = arh_object%n_ao

        ! get rotated density matrix in AO and OAO basis
        allocate(rot_dm_ao(n_ao, n_ao, 1), rot_dm_oao(n_ao, n_ao, 1), &
                 fock_ao(n_ao, n_ao, 1), v_nonlinear_ao(n_ao, n_ao, 1))
        call rotate_dm_ao(kappa, arh_object%n_particle, n_ao, rot_dm_ao, error, &
                          rot_dm_oao)
        if (error /= 0) return

        ! calculate mean-field energy
        call arh_object%evaluate_dm_cs(rot_dm_ao(:, :, 1), energy, fock_ao(:, :, 1), &
                                       v_nonlinear_ao(:, :, 1), error)
        if (error /= 0) return

        ! update list of density, Fock and non-linear potential matrices
        if (allocated(arh_object%dm_list)) then
            if (.not. density_in_history(rot_dm_oao)) then
                call prepend(arh_object%dm_list, rot_dm_oao)
                call prepend(arh_object%fock_list, &
                             symmetric_transformation(arh_object%s_inv_sqrt, fock_ao))
                call prepend(arh_object%v_nonlinear_list, symmetric_transformation( &
                    arh_object%s_inv_sqrt, v_nonlinear_ao))
            end if
        end if

    end function obj_func_arh_cs

    function obj_func_arh_os(kappa, error) result(energy)
        !
        ! this function defines the energy evaluation in the OAO basis for the
        ! open-shell case, which also adds the evaluated point with its same-spin,
        ! opposite-spin and non-linear potentials to the history, unless it is already
        ! there
        !
        use otr_oao, only: rotate_dm_ao, symmetric_transformation

        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        integer(ip) :: n_ao, n_particle
        real(rp), allocatable :: rot_dm_ao(:, :, :), rot_dm_oao(:, :, :), &
                                 v_same_spin_ao(:, :, :), v_opposite_spin_ao(:, :, :), &
                                 v_nonlinear_ao(:, :, :)

        ! initialize energy in case of error
        energy = 0.0_rp

        ! number of AOs and number of particles
        n_ao = arh_object%n_ao
        n_particle = arh_object%n_particle

        ! get rotated density matrix in AO and OAO basis
        allocate(rot_dm_ao(n_ao, n_ao, n_particle), &
                 rot_dm_oao(n_ao, n_ao, n_particle), &
                 v_same_spin_ao(n_ao, n_ao, n_particle), &
                 v_opposite_spin_ao(n_ao, n_ao, n_particle), &
                 v_nonlinear_ao(n_ao, n_ao, n_particle))
        call rotate_dm_ao(kappa, n_particle, n_ao, rot_dm_ao, error, rot_dm_oao)
        if (error /= 0) return

        ! calculate mean-field energy
        call arh_object%evaluate_dm_os(rot_dm_ao, energy, v_same_spin=v_same_spin_ao, &
                                       v_opposite_spin=v_opposite_spin_ao, &
                                       v_nonlinear=v_nonlinear_ao, error=error)
        if (error /= 0) return

        ! update list of density and potential matrices
        if (allocated(arh_object%dm_list)) then
            if (.not. density_in_history(rot_dm_oao)) then
                call prepend(arh_object%dm_list, rot_dm_oao)
                call prepend(arh_object%v_same_spin_list, symmetric_transformation( &
                    arh_object%s_inv_sqrt, v_same_spin_ao))
                call prepend( &
                    arh_object%v_opposite_spin_list, &
                    symmetric_transformation(arh_object%s_inv_sqrt, v_opposite_spin_ao))
                call prepend(arh_object%v_nonlinear_list, symmetric_transformation( &
                    arh_object%s_inv_sqrt, v_nonlinear_ao))
            end if
        end if

    end function obj_func_arh_os

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

        integer(ip) :: n_ao, n_particle, i, n_list, n_acc, n_acc_nonlinear
        real(rp) :: min_residual
        real(rp), allocatable :: fock_ao(:, :, :), v_nonlinear_ao(:, :, :), &
                                 dm_diff(:, :, :, :), v_linear_diff(:, :, :, :), &
                                 v_nonlinear_diff(:, :, :, :), chol(:, :), &
                                 chol_nonlinear(:, :)
        integer(ip), allocatable :: map(:), map_nonlinear(:)
        logical, allocatable :: keep_nonlinear(:)

        external :: dgemm

        ! initialize error flag
        error = 0

        ! check if orbitals are actually rotated
        if ((sum(abs(kappa)) > 0.0_rp) .or. &
            (abs(arh_object%energy) <= numerical_zero) .or. (.not. ( &
                allocated(oao_object%grad) .and. allocated(oao_object%h_diag) .and. &
                (allocated(arh_object%dm_list))))) then
            ! number of AOs
            n_ao = arh_object%n_ao

            ! number of particles
            n_particle = arh_object%n_particle

            ! update list of density, Fock and non-linear potential matrices
            if (allocated(arh_object%dm_list)) then
                if (.not. density_in_history(arh_object%dm_oao)) then
                    call prepend(arh_object%dm_list, arh_object%dm_oao)
                    call prepend(arh_object%fock_list, arh_object%fock_oao)
                    call prepend(arh_object%v_nonlinear_list, &
                                 arh_object%v_nonlinear_oao)
                end if
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
            call arh_object%evaluate_dm_cs(arh_object%dm_ao(:, :, 1), &
                                           arh_object%energy, fock_ao(:, :, 1), &
                                           v_nonlinear_ao(:, :, 1), error)
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

            ! factorize the density-matrix-difference history for the linear part which 
            ! resolves linear dependencies in the history; the coupling formulas below 
            ! never need an explicit metric inverse since everything is expressed in 
            ! the resulting orthonormalized basis via a single triangular solve
            call factorize_history(reshape(dm_diff, [n_ao * n_ao, n_list]), chol, map, &
                                   n_acc)

            ! get potential differences for linear (Coulomb and exact exchange) and 
            ! non-linear (XC) parts
            allocate(v_linear_diff(n_ao, n_ao, n_particle, n_list), &
                     v_nonlinear_diff(n_ao, n_ao, n_particle, n_list))
            do i = 1, n_list
                v_nonlinear_diff(:, :, :, i) = &
                    arh_object%v_nonlinear_list(:, :, :, i) - arh_object%v_nonlinear_oao
                v_linear_diff(:, :, :, i) = arh_object%fock_list(:, :, :, i) - &
                                            arh_object%fock_oao - &
                                            v_nonlinear_diff(:, :, :, i)
            end do

            ! factorize the same history for the non-linear part, which can only be
            ! done once its response is known since the error in that response sets
            ! the shortest residual a direction has to contribute
            keep_nonlinear = history_step_mask(dm_diff)
            min_residual = resolvable_residual( &
                reshape(dm_diff, [n_ao * n_ao, n_list]), &
                reshape(v_nonlinear_diff, [n_ao * n_ao, n_list]), keep_nonlinear)
            call factorize_history(reshape(dm_diff, [n_ao * n_ao, n_list]), &
                                   chol_nonlinear, map_nonlinear, n_acc_nonlinear, &
                                   keep_nonlinear, min_residual)

            ! for MS-SR1 and cache the direction vectors
            if (arh_object%settings%arh_type == "ms_sr1") then
                ! get inverted A matrix
                ! linear part: this is exact since Coulomb and exact exchange are
                ! linear in the density matrix
                call get_ms_a_inv(dm_diff, v_linear_diff, map, chol, arh_object%a_inv, &
                                  arh_object%settings, error)
                if (error /= 0) return
                ! non-linear part: kept on its own system so that it contracts against
                ! its own inverse and the exact linear secant relationship is not
                ! averaged with the approximate one
                call get_ms_a_inv(dm_diff, v_nonlinear_diff, map_nonlinear, &
                                  chol_nonlinear, arh_object%a_inv_comb, &
                                  arh_object%settings, error)
                if (error /= 0) return

                ! cache the packed history-projection directions the low-rank Hessian 
                ! factors are assembled from, rebased into the orthonormalized S-basis
                call cache_history_dirs(v_linear_diff, arh_object%dm_oao, n_list, &
                                        arh_object%n_param, map, chol, &
                                        arh_object%linear_potential_dirs)
                call cache_history_dirs( &
                    v_nonlinear_diff, arh_object%dm_oao, n_list, arh_object%n_param, &
                    map_nonlinear, chol_nonlinear, arh_object%nonlinear_potential_dirs)
            ! ARH and related methods
            else
                ! cache the packed history-projection directions the low-rank Hessian 
                ! factors are assembled from, rebased into the orthonormalized S-basis
                ! belonging to the linear and non-linear systems
                call cache_history_dirs(dm_diff, arh_object%dm_oao, n_list, &
                                        arh_object%n_param, map, chol, &
                                        arh_object%dm_dirs)
                call cache_history_dirs(dm_diff, arh_object%dm_oao, n_list, &
                                        arh_object%n_param, map_nonlinear, &
                                        chol_nonlinear, arh_object%dm_dirs_nonlinear)
                if (arh_object%settings%arh_type /= "ms_sp") then
                    call cache_history_dirs(v_linear_diff, arh_object%dm_oao, n_list, &
                                            arh_object%n_param, map, chol, &
                                            arh_object%linear_potential_dirs)
                    call cache_history_dirs(v_nonlinear_diff, arh_object%dm_oao, &
                                            n_list, arh_object%n_param, map_nonlinear, &
                                            chol_nonlinear, &
                                            arh_object%nonlinear_potential_dirs)
                end if

                ! construct A = S^T Y for the linear and non-linear system, 
                ! congruence-transformed into its own orthonormalized S-basis
                if (arh_object%settings%arh_type == "ms_sp" .or. &
                    arh_object%settings%arh_type == "ms_psb") then
                    arh_object%a_sym = &
                        build_a_transformed(dm_diff, v_linear_diff, map, chol)
                    arh_object%a_sym_nonlinear = build_a_transformed( &
                        dm_diff, v_nonlinear_diff, map_nonlinear, chol_nonlinear)
                end if
                deallocate(v_nonlinear_diff, v_linear_diff)
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

        integer(ip) :: n_ao, n_particle, i, n_list, n_acc1, n_acc2, n_acc_nl, &
                       n_acc1_nl, n_acc2_nl
        real(rp) :: min_residual
        integer(ip), allocatable :: map1_nl(:), map2_nl(:), map_comb_nl(:)
        logical, allocatable :: keep_nonlinear(:)
        real(rp), allocatable :: &
            fock_ao(:, :, :), fock_oao(:, :, :), v_same_spin_ao(:, :, :), &
            v_opposite_spin_ao(:, :, :), v_nonlinear_ao(:, :, :), dm_diff(:, :, :, :), &
            v_same_spin_diff(:, :, :, :), v_opposite_spin_diff(:, :, :, :), &
            v_nonlinear_diff(:, :, :, :), v_zero(:, :, :, :), chol1(:, :), &
            chol2(:, :), chol1_nl(:, :), chol2_nl(:, :), chol_comb(:, :), &
            chol_comb_nl(:, :), chol_nl(:, :)
        integer(ip), allocatable :: map1(:), map2(:), map_comb(:), map_nl(:)

        external :: dgemm

        ! initialize error flag
        error = 0

        ! check if orbitals are actually rotated
        if ((sum(abs(kappa)) > 0.0_rp) .or. &
            (abs(arh_object%energy) <= numerical_zero) .or. (.not. ( &
                allocated(oao_object%grad) .and. allocated(oao_object%h_diag) .and. &
                (allocated(arh_object%dm_list))))) then
            ! number of AOs
            n_ao = arh_object%n_ao

            ! number of particles
            n_particle = arh_object%n_particle

            ! update list of density and potential matrices
            if (allocated(arh_object%dm_list)) then
                if (.not. density_in_history(arh_object%dm_oao)) then
                    call prepend(arh_object%dm_list, arh_object%dm_oao)
                    call prepend(arh_object%v_same_spin_list, &
                                 arh_object%v_same_spin_oao)
                    call prepend(arh_object%v_opposite_spin_list, &
                                 arh_object%v_opposite_spin_oao)
                    call prepend(arh_object%v_nonlinear_list, &
                                 arh_object%v_nonlinear_oao)
                end if
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
            call arh_object%evaluate_dm_os(arh_object%dm_ao, arh_object%energy, &
                                           fock_ao, v_same_spin_ao, &
                                           v_opposite_spin_ao, v_nonlinear_ao, error)
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

            ! factorize the density-matrix-difference history once per channel which 
            ! resolves linear dependencies in the history; the coupling formulas below 
            ! never need an explicit metric inverse since everything is expressed in 
            ! the resulting orthonormalized basis via a single triangular solve; the 
            ! two channels never mix in this factorization (the metric itself is
            ! block-diagonal across channels), so they are factorized independently
            ! and then combined into one block-diagonal (chol, map) pair wherever the 
            ! coupled types or the linear multisecant system need to act on both 
            ! channels together
            call factorize_history(reshape(dm_diff(:, :, 1, :), &
                                           [n_ao * n_ao, n_list]), chol1, map1, n_acc1)
            call factorize_history(reshape(dm_diff(:, :, 2, :), &
                                           [n_ao * n_ao, n_list]), chol2, map2, n_acc2)
            call combine_channels(chol1, map1, chol2, map2, n_list, chol_comb, map_comb)

            ! keep the same-spin and opposite-spin potential differences separate 
            ! from the non-linear one since the former are exact at any distance from 
            ! the current density and therefore take the full history, while the 
            ! non-linear one describes a drifting Hessian and is fitted only to the 
            ! history entries the screening below admits
            allocate(v_same_spin_diff(n_ao, n_ao, n_particle, n_list), &
                     v_nonlinear_diff(n_ao, n_ao, n_particle, n_list))
            do i = 1, n_list
                v_same_spin_diff(:, :, :, i) = &
                    arh_object%v_same_spin_list(:, :, :, i) - arh_object%v_same_spin_oao
                v_nonlinear_diff(:, :, :, i) = &
                    arh_object%v_nonlinear_list(:, :, :, i) - arh_object%v_nonlinear_oao
            end do

            ! MS-SR1
            if (arh_object%settings%arh_type == "ms_sr1") then
                ! get inverted A matrix
                ! linear part: get spin-separated multisecant SR1 matrix for which 
                ! separation is exact since Coulomb and exact exchange are linear in 
                ! the density matrix
                call get_ms_a_inv_os_linear( &
                    dm_diff, v_same_spin_diff, v_opposite_spin_diff, map_comb, &
                    chol_comb, arh_object%a_inv, n_ao, arh_object%settings, error)
                if (error /= 0) return
                ! non-linear part: get spin-combined multisecant SR1 matrix; the
                ! non-linear response mixes both channels at once, so this needs its
                ! own, separate combined-flat factorization of the same history
                keep_nonlinear = history_step_mask(dm_diff)
                min_residual = resolvable_residual( &
                    reshape(dm_diff, [n_ao * n_ao * n_particle, n_list]), &
                    reshape(v_nonlinear_diff, [n_ao * n_ao * n_particle, n_list]), &
                    keep_nonlinear)
                call factorize_history( &
                    reshape(dm_diff, [n_ao * n_ao * n_particle, n_list]), chol_nl, &
                    map_nl, n_acc_nl, keep_nonlinear, min_residual)
                call get_ms_a_inv(dm_diff, v_nonlinear_diff, map_nl, chol_nl, &
                                  arh_object%a_inv_comb, arh_object%settings, error)
                if (error /= 0) return

                ! cache the packed history-projection directions the low-rank Hessian 
                ! factors are assembled from; the linear potential directions combine 
                ! the same-/opposite-spin channels, while the non-linear potential 
                ! directions need no channel-splitting; rebased into the
                ! orthonormalized S-basis
                call cache_combined_channel_dirs( &
                    v_same_spin_diff, v_opposite_spin_diff, arh_object%dm_oao, n_list, &
                    arh_object%n_param, n_particle, map_comb, chol_comb, &
                    arh_object%linear_potential_dirs)
                call cache_history_dirs(v_nonlinear_diff, arh_object%dm_oao, n_list, &
                                        arh_object%n_param, map_nl, chol_nl, &
                                        arh_object%nonlinear_potential_dirs)
            ! ARH and related methods
            else
                ! the non-linear response has no opposite-spin counterpart, so the
                ! combined-channel routine is handed a vanishing one
                allocate(v_zero(n_ao, n_ao, n_particle, n_list))
                v_zero = 0.0_rp

                ! screened per-channel factorization for the non-linear system; both
                ! channels are screened on the step length of the whole density matrix
                ! rather than on their own channel's, since the non-linear potential of
                ! either channel is a functional of both spin densities and its
                ! staleness is therefore set by the total step
                keep_nonlinear = history_step_mask(dm_diff)
                min_residual = resolvable_residual( &
                    reshape(dm_diff(:, :, 1, :), [n_ao * n_ao, n_list]), &
                    reshape(v_nonlinear_diff(:, :, 1, :), [n_ao * n_ao, n_list]), &
                    keep_nonlinear)
                call factorize_history( &
                    reshape(dm_diff(:, :, 1, :), [n_ao * n_ao, n_list]), chol1_nl, &
                    map1_nl, n_acc1_nl, keep_nonlinear, min_residual)
                min_residual = resolvable_residual( &
                    reshape(dm_diff(:, :, 2, :), [n_ao * n_ao, n_list]), &
                    reshape(v_nonlinear_diff(:, :, 2, :), [n_ao * n_ao, n_list]), &
                    keep_nonlinear)
                call factorize_history( &
                    reshape(dm_diff(:, :, 2, :), [n_ao * n_ao, n_list]), chol2_nl, &
                    map2_nl, n_acc2_nl, keep_nonlinear, min_residual)
                call combine_channels(chol1_nl, map1_nl, chol2_nl, map2_nl, n_list, &
                                      chol_comb_nl, map_comb_nl)

                ! cache the packed history-projection directions the low-rank Hessian 
                ! factors are assembled from; the density matrix directions isolate
                ! each channel separately, the potential directions combine channels,
                ! each rebased into the S-basis belonging to its own system
                call cache_channel_split_dirs(dm_diff, arh_object%dm_oao, n_list, &
                                              arh_object%n_param, n_particle, &
                                              map_comb, chol_comb, arh_object%dm_dirs)
                call cache_channel_split_dirs( &
                    dm_diff, arh_object%dm_oao, n_list, arh_object%n_param, &
                    n_particle, map_comb_nl, chol_comb_nl, arh_object%dm_dirs_nonlinear)
                if (arh_object%settings%arh_type /= "ms_sp") then
                    call cache_combined_channel_dirs( &
                        v_same_spin_diff, v_opposite_spin_diff, arh_object%dm_oao, &
                        n_list, arh_object%n_param, n_particle, map_comb, chol_comb, &
                        arh_object%linear_potential_dirs)
                    call cache_combined_channel_dirs( &
                        v_nonlinear_diff, v_zero, arh_object%dm_oao, n_list, &
                        arh_object%n_param, n_particle, map_comb_nl, chol_comb_nl, &
                        arh_object%nonlinear_potential_dirs)
                end if

                ! construct A = S^T Y for the linear and non-linear system, 
                ! congruence-transformed into its own combined orthonormalized S-basis
                if (arh_object%settings%arh_type == "ms_sp" .or. &
                    arh_object%settings%arh_type == "ms_psb") then
                    arh_object%a_sym = build_a_block_linear_os( &
                        dm_diff, v_same_spin_diff, v_opposite_spin_diff, n_ao, &
                        map_comb, chol_comb)
                    arh_object%a_sym_nonlinear = build_a_block_nonlinear_os( &
                        dm_diff, v_nonlinear_diff, n_ao, map_comb_nl, chol_comb_nl)
                end if
                deallocate(v_same_spin_diff, v_nonlinear_diff, v_zero)
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
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, &
                       arh_object%fock_vv(:, :, i) - arh_object%fock_oo(:, :, i), &
                       n_ao, x_full(:, :, i), n_ao, 0.0_rp, hess_x_full(:, :, i), n_ao)
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
            call dgemv("N", n_dirs, n_dirs, 1.0_rp, arh_object%coupling_matrix, &
                       n_dirs, projected_x, 1_ip, 0.0_rp, coupled_x, 1_ip)
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
                                 bracket_matrix(:, :), projected_x(:), bracket_rhs(:), &
                                 bracket_solution(:), correction(:)
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

    subroutine cache_history_dirs(v_diff, dm_oao, n_list, n_param, map, chol, dirs)
        !
        ! this subroutine caches the packed history-projection directions the
        ! low-rank part of the approximate Hessian is built from: each history entry
        ! is projected onto the occupied-virtual/virtual-occupied subspace and packed
        ! into the same antisymmetric parameter space as a trial vector, then rebased
        ! into the orthonormalized basis defined by map and chol
        !
        use otr_oao, only: project_asymm, pack_asymm

        real(rp), intent(in) :: v_diff(:, :, :, :), dm_oao(:, :, :), chol(:, :)
        integer(ip), intent(in) :: n_list, n_param, map(:)
        real(rp), intent(out), allocatable :: dirs(:, :)

        integer(ip) :: k
        real(rp), allocatable :: projected(:, :, :), raw(:, :)

        allocate(raw(n_param, n_list))
        do k = 1, n_list
            projected = project_asymm(v_diff(:, :, :, k), dm_oao)
            raw(:, k) = pack_asymm(projected, n_param)
        end do
        dirs = rebase_dirs(raw, map, chol)

    end subroutine cache_history_dirs

    subroutine cache_history_projections_channel(v_diff, channel, dm_oao, n_list, &
                                                 n_param, n_particle, projections)
        !
        ! this subroutine caches the packed, projected history columns of a single spin 
        ! channel for open-shell systems; each history column is embedded into only the 
        ! given channel (the other channels zeroed) before projecting and packing, 
        ! since the open-shell response contracts a channel's own difference against 
        ! only that channel's density response
        !
        use otr_oao, only: project_asymm, pack_asymm

        real(rp), intent(in) :: v_diff(:, :, :, :), dm_oao(:, :, :)
        integer(ip), intent(in) :: channel, n_list, n_param, n_particle
        real(rp), intent(out), allocatable :: projections(:, :)

        integer(ip) :: k, n_ao
        real(rp), allocatable :: embedded(:, :, :), projected(:, :, :)

        n_ao = size(dm_oao, 1)
        allocate(projections(n_param, n_list), embedded(n_ao, n_ao, n_particle))
        embedded = 0.0_rp
        do k = 1, n_list
            embedded(:, :, channel) = v_diff(:, :, channel, k)
            projected = project_asymm(embedded, dm_oao)
            projections(:, k) = pack_asymm(projected, n_param)
        end do
        deallocate(embedded)

    end subroutine cache_history_projections_channel

    subroutine cache_channel_split_dirs(v_diff, dm_oao, n_list, n_param, n_particle, &
                                        map, chol, dirs)
        !
        ! this subroutine caches the open-shell history projections one spin channel
        ! at a time: column i holds channel 1 of entry i, column n_list+i holds
        ! channel 2; the two are never summed, since the per-channel metric
        ! contraction does not mix channels; the concatenated columns are then rebased 
        ! into the orthonormalized basis defined by map and chol
        !
        real(rp), intent(in) :: v_diff(:, :, :, :), dm_oao(:, :, :), chol(:, :)
        integer(ip), intent(in) :: n_list, n_param, n_particle, map(:)
        real(rp), intent(out), allocatable :: dirs(:, :)

        real(rp), allocatable :: channel1(:, :), channel2(:, :), raw(:, :)

        call cache_history_projections_channel(v_diff, 1_ip, dm_oao, n_list, n_param, &
                                               n_particle, channel1)
        call cache_history_projections_channel(v_diff, 2_ip, dm_oao, n_list, n_param, &
                                               n_particle, channel2)

        allocate(raw(n_param, 2 * n_list))
        raw(:, 1:n_list) = channel1
        raw(:, n_list + 1:2 * n_list) = channel2
        deallocate(channel1, channel2)
        dirs = rebase_dirs(raw, map, chol)

    end subroutine cache_channel_split_dirs

    subroutine cache_combined_channel_dirs(v_same, v_opp, dm_oao, n_list, n_param, &
                                           n_particle, map, chol, dirs)
        !
        ! this subroutine caches the open-shell history projections two channels at a
        ! time: column i sums channel 1 of the same-spin potential with channel 2 of
        ! the opposite-spin potential, column n_list+i the mirror image; every caller
        ! pairs a same-spin difference with the opposite-spin one, differing only in
        ! what is passed as the same-spin potential; the combined columns are then
        ! rebased into the orthonormalized basis defined by map and chol
        !
        real(rp), intent(in) :: v_same(:, :, :, :), v_opp(:, :, :, :), &
                                dm_oao(:, :, :), chol(:, :)
        integer(ip), intent(in) :: n_list, n_param, n_particle, map(:)
        real(rp), intent(out), allocatable :: dirs(:, :)

        real(rp), allocatable :: same1(:, :), same2(:, :), opp1(:, :), opp2(:, :), &
                                 raw(:, :)

        call cache_history_projections_channel(v_same, 1_ip, dm_oao, n_list, n_param, &
                                               n_particle, same1)
        call cache_history_projections_channel(v_same, 2_ip, dm_oao, n_list, n_param, &
                                               n_particle, same2)
        call cache_history_projections_channel(v_opp, 1_ip, dm_oao, n_list, n_param, &
                                               n_particle, opp1)
        call cache_history_projections_channel(v_opp, 2_ip, dm_oao, n_list, n_param, &
                                               n_particle, opp2)

        allocate(raw(n_param, 2 * n_list))
        raw(:, 1:n_list) = same1 + opp2
        raw(:, n_list + 1:2 * n_list) = same2 + opp1
        deallocate(same1, same2, opp1, opp2)
        dirs = rebase_dirs(raw, map, chol)

    end subroutine cache_combined_channel_dirs

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
        integer(ip) :: n_linear, n_nonlinear, n_total, i
        real(rp) :: shell_scale

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
        ! independent blocks, each contracted against its own inverse rather than
        ! against a single inverse of the summed potential
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
            if (.not. allocated(arh_object%dm_dirs) .or. &
                .not. allocated(arh_object%dm_dirs_nonlinear)) return
            n_linear = size(arh_object%dm_dirs, 2)
            n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
            n_total = n_linear + n_nonlinear
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_linear + 1:) = arh_object%dm_dirs_nonlinear
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_linear, :n_linear) = 8.0_rp * shell_scale * &
                                                               arh_object%a_sym
            arh_object%coupling_matrix(n_linear + 1:, n_linear + 1:) = &
                8.0_rp * shell_scale * arh_object%a_sym_nonlinear

        ! symmetrized ARH: the density and potential difference histories couple in
        ! both directions, so the coupling matrix is purely off-diagonal; also 
        ! introduces an additional factor of 1/2
        case ("symm_arh")
            if (.not. allocated(arh_object%dm_dirs) .or. &
                .not. allocated(arh_object%linear_potential_dirs)) return
            n_linear = size(arh_object%dm_dirs, 2)
            n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
            n_total = 2 * (n_linear + n_nonlinear)
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_linear + 1:2 * n_linear) = &
                arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, 2 * n_linear + 1:2 * n_linear + &
                                      n_nonlinear) = arh_object%dm_dirs_nonlinear
            arh_object%expansion_dirs(:, 2 * n_linear + n_nonlinear + 1:) = &
                arh_object%nonlinear_potential_dirs
            ! the metric scaling collapses in the orthonormalized S-basis
            arh_object%coupling_matrix = 0.0_rp
            do i = 1, n_linear
                arh_object%coupling_matrix(i, n_linear + i) = 4.0_rp * shell_scale
                arh_object%coupling_matrix(n_linear + i, i) = 4.0_rp * shell_scale
            end do
            do i = 1, n_nonlinear
                arh_object%coupling_matrix(2 * n_linear + i, 2 * n_linear + &
                                           n_nonlinear + i) = 4.0_rp * shell_scale
                arh_object%coupling_matrix(2 * n_linear + n_nonlinear + i, &
                                           2 * n_linear + i) = 4.0_rp * shell_scale
            end do

        ! multisecant PSB: the symmetrized ARH coupling with an additional
        ! density-density block subtracting the doubly counted curvature
        case ("ms_psb")
            if (.not. allocated(arh_object%dm_dirs) .or. &
                .not. allocated(arh_object%linear_potential_dirs)) return
            n_linear = size(arh_object%dm_dirs, 2)
            n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
            n_total = 2 * (n_linear + n_nonlinear)
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_linear + 1:2 * n_linear) = &
                arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, 2 * n_linear + 1:2 * n_linear + &
                                      n_nonlinear) = arh_object%dm_dirs_nonlinear
            arh_object%expansion_dirs(:, 2 * n_linear + n_nonlinear + 1:) = &
                arh_object%nonlinear_potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_linear, :n_linear) = -8.0_rp * shell_scale * &
                                                               arh_object%a_sym
            arh_object%coupling_matrix(2 * n_linear + 1:2 * n_linear + n_nonlinear, &
                                       2 * n_linear + 1:2 * n_linear + n_nonlinear) = &
                -8.0_rp * shell_scale * arh_object%a_sym_nonlinear
            ! the metric scaling in the off-diagonal blocks collapses in the 
            ! orthonormalized S-basis
            do i = 1, n_linear
                arh_object%coupling_matrix(i, n_linear + i) = 8.0_rp * shell_scale
                arh_object%coupling_matrix(n_linear + i, i) = 8.0_rp * shell_scale
            end do
            do i = 1, n_nonlinear
                arh_object%coupling_matrix(2 * n_linear + i, 2 * n_linear + &
                                           n_nonlinear + i) = 8.0_rp * shell_scale
                arh_object%coupling_matrix(2 * n_linear + n_nonlinear + i, &
                                           2 * n_linear + i) = 8.0_rp * shell_scale
            end do

        ! standard ARH: the response expands in the potential difference history
        ! while contracting against the density difference history, so unlike every
        ! other type the two sets of directions differ
        case ("arh")
            if (.not. allocated(arh_object%dm_dirs) .or. &
                .not. allocated(arh_object%linear_potential_dirs)) return
            n_linear = size(arh_object%dm_dirs, 2)
            n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
            n_total = n_linear + n_nonlinear
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%n_param, n_total), &
                     arh_object%projection_dirs(arh_object%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, n_linear + 1:) = &
                arh_object%nonlinear_potential_dirs
            arh_object%projection_dirs(:, :n_linear) = arh_object%dm_dirs
            arh_object%projection_dirs(:, n_linear + 1:) = arh_object%dm_dirs_nonlinear
            arh_object%coupling_matrix = 0.0_rp
            ! the metric scaling collapses in the orthonormalized S-basis
            do i = 1, n_total
                arh_object%coupling_matrix(i, i) = 8.0_rp * shell_scale
            end do

        case default
            return
        end select

        ! every type except standard ARH expands in and contracts against the same
        ! set of directions
        if (allocated(arh_object%expansion_dirs) .and. &
            .not. allocated(arh_object%projection_dirs)) &
            arh_object%projection_dirs = arh_object%expansion_dirs

    end subroutine get_low_rank_hess_factors

    subroutine build_a_part(dm_diff, v_diff, a)
        !
        ! this subroutine constructs A = S^T Y for a single part of a coupled
        ! potential-difference response, flattened over the AO and particle
        ! dimensions, and symmetrizes it
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :), v_diff(:, :, :, :)
        real(rp), intent(out) :: a(:, :)

        integer(ip) :: n_dm, flat_len
        external :: dgemm

        ! A = S^T Y on the full history, flattened over the AO and particle
        ! dimensions
        n_dm = size(dm_diff, 4)
        flat_len = size(dm_diff, 1) * size(dm_diff, 2) * size(dm_diff, 3)
        call dgemm("T", "N", n_dm, n_dm, flat_len, 1.0_rp, dm_diff, flat_len, v_diff, &
                   flat_len, 0.0_rp, a, n_dm)

        ! symmetrize A
        a = 0.5_rp * (a + transpose(a))

    end subroutine build_a_part

    function build_a_transformed(dm_diff, v_diff, map, chol) result(a_t)
        !
        ! this function builds A = S^T Y for a single (linear or non-linear) part of
        ! the coupled potential-difference response and congruence-transforms it into
        ! the orthonormalized S-basis belonging to that part
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :), v_diff(:, :, :, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: a_t(:, :)

        integer(ip) :: n_diff
        real(rp), allocatable :: a(:, :)

        n_diff = size(dm_diff, 4)
        allocate(a(n_diff, n_diff))
        call build_a_part(dm_diff, v_diff, a)
        a_t = congruence_transform(a, map, chol)
        deallocate(a)

    end function build_a_transformed

    function build_a_block_linear_os(dm_diff, v_same_linear, v_opp, n_ao, map, chol) &
        result(a_block)
        !
        ! this function builds the dense, cross-channel-symmetrized open-shell
        ! A = S^T Y matrix of the linear response entering the MS-SP and MS-PSB 
        ! contributions; the same-spin response fills the per-channel diagonal blocks 
        ! and the opposite-spin response the off-diagonal ones, and the result is 
        ! congruence-transformed to the orthonormalized, rank-independent S-basis
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :), v_same_linear(:, :, :, :), &
                                v_opp(:, :, :, :), chol(:, :)
        integer(ip), intent(in) :: n_ao, map(:)
        real(rp), allocatable :: a_block(:, :)

        integer(ip) :: n_diff, j
        real(rp), allocatable :: a_same(:, :, :), a_opp(:, :, :), a_full(:, :)
        external :: dgemm

        n_diff = size(dm_diff, 4)
        allocate(a_same(n_diff, n_diff, 2), a_opp(n_diff, n_diff, 2))

        do j = 1, 2
            call build_a_part( &
                reshape(dm_diff(:, :, j, :), [n_ao, n_ao, 1_ip, n_diff]), &
                reshape(v_same_linear(:, :, j, :), [n_ao, n_ao, 1_ip, n_diff]), &
                a_same(:, :, j))
            call dgemm("T", "N", n_diff, n_diff, n_ao * n_ao, 1.0_rp, &
                       reshape(dm_diff(:, :, j, :), [n_ao * n_ao, n_diff]), &
                       n_ao * n_ao, reshape(v_opp(:, :, j, :), [n_ao * n_ao, n_diff]), &
                       n_ao * n_ao, 0.0_rp, a_opp(:, :, j), n_diff)
        end do

        ! cross-symmetrize off-diagonal blocks exactly, since the opposite-spin
        ! potential difference is purely linear
        call cross_symmetrize(a_opp(:, :, 1), a_opp(:, :, 2))

        allocate(a_full(2 * n_diff, 2 * n_diff))
        a_full(1:n_diff, 1:n_diff) = a_same(:, :, 1)
        a_full(1:n_diff, n_diff + 1:2 * n_diff) = a_opp(:, :, 1)
        a_full(n_diff + 1:2 * n_diff, 1:n_diff) = a_opp(:, :, 2)
        a_full(n_diff + 1:2 * n_diff, n_diff + 1:2 * n_diff) = a_same(:, :, 2)
        deallocate(a_same, a_opp)
        a_block = congruence_transform(a_full, map, chol)
        deallocate(a_full)

    end function build_a_block_linear_os

    function build_a_block_nonlinear_os(dm_diff, v_nonlinear, n_ao, map, chol) &
        result(a_block)
        !
        ! this function builds the open-shell A = S^T Y matrix of the non-linear 
        ! response entering the MS-SP and MS-PSB contributions; the non-linear response
        ! has no opposite-spin counterpart to cross-symmetrize against, so it fills the
        ! per-channel diagonal blocks only, and the result is congruence-transformed to
        ! the orthonormalized, rank-independent S-basis
        !
        real(rp), intent(in) :: dm_diff(:, :, :, :), v_nonlinear(:, :, :, :), chol(:, :)
        integer(ip), intent(in) :: n_ao, map(:)
        real(rp), allocatable :: a_block(:, :)

        integer(ip) :: n_diff, j
        real(rp), allocatable :: a_same(:, :, :), a_full(:, :)

        n_diff = size(dm_diff, 4)
        allocate(a_same(n_diff, n_diff, 2))
        do j = 1, 2
            call build_a_part( &
                reshape(dm_diff(:, :, j, :), [n_ao, n_ao, 1_ip, n_diff]), &
                reshape(v_nonlinear(:, :, j, :), [n_ao, n_ao, 1_ip, n_diff]), &
                a_same(:, :, j))
        end do

        allocate(a_full(2 * n_diff, 2 * n_diff))
        a_full = 0.0_rp
        a_full(1:n_diff, 1:n_diff) = a_same(:, :, 1)
        a_full(n_diff + 1:2 * n_diff, n_diff + 1:2 * n_diff) = a_same(:, :, 2)
        deallocate(a_same)
        a_block = congruence_transform(a_full, map, chol)
        deallocate(a_full)

    end function build_a_block_nonlinear_os

    subroutine cross_symmetrize(a12, a21)
        !
        ! this subroutine cross-symmetrizes two related off-diagonal blocks a12(i, k)
        ! and a21(k, i) of a larger matrix by averaging each with the transpose of the
        ! other, which is appropriate here since a12 and transpose(a21) are equal in
        ! exact arithmetic and any observed mismatch is therefore numerical noise
        !
        real(rp), intent(inout) :: a12(:, :), a21(:, :)

        real(rp), allocatable :: averaged(:, :)

        averaged = 0.5_rp * (a12 + transpose(a21))
        a12 = averaged
        a21 = transpose(averaged)
        deallocate(averaged)

    end subroutine cross_symmetrize

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

    function history_step_mask(dm_diff) result(keep)
        !
        ! this function returns, for every history entry, whether the non-linear
        ! multisecant system may use it, depending on step length
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: dm_diff(:, :, :, :)
        logical, allocatable :: keep(:)

        integer(ip) :: n_list, flat_len, i
        real(rp) :: shortest, cutoff
        real(rp), allocatable :: steps(:)
        real(rp), external :: dnrm2

        n_list = size(dm_diff, 4)
        allocate(keep(n_list), steps(n_list))
        keep = .true.
        if (n_list < 2) return

        ! get step lengths
        flat_len = size(dm_diff, 1) * size(dm_diff, 2) * size(dm_diff, 3)
        do i = 1, n_list
            steps(i) = dnrm2(flat_len, dm_diff(:, :, :, i), 1_ip)
        end do

        ! get shortest numerically relevant step
        shortest = huge(1.0_rp)
        do i = 1, n_list
            if (steps(i) > numerical_zero * maxval(steps)) &
                shortest = min(shortest, steps(i))
        end do
        if (shortest >= huge(1.0_rp)) return

        ! define cut off and decide on steps to keep
        cutoff = history_step_ratio * shortest
        keep = steps <= cutoff
        deallocate(steps)

    end function history_step_mask

    function median(x) result(med)
        !
        ! this function returns the median of an array
        !
        real(rp), intent(in) :: x(:)
        real(rp) :: med

        integer(ip) :: i, j, n
        real(rp) :: val
        real(rp), allocatable :: sorted(:)

        n = size(x)
        if (n == 0) then
            med = 0.0_rp
            return
        end if

        ! sort a copy by insertion, which is short and fast enough for the few hundred
        ! history directions this is called with at most
        sorted = x
        do i = 2, n
            val = sorted(i)
            j = i - 1
            do while (j >= 1)
                if (sorted(j) <= val) exit
                sorted(j + 1) = sorted(j)
                j = j - 1
            end do
            sorted(j + 1) = val
        end do

        ! average the two central entries if there is no single one
        if (mod(n, 2_ip) == 1) then
            med = sorted((n + 1) / 2)
        else
            med = 0.5_rp * (sorted(n / 2) + sorted(n / 2 + 1))
        end if
        deallocate(sorted)

    end function median

    function resolvable_residual(steps, responses, keep) result(min_residual)
        !
        ! this function returns the shortest residual norm a history direction has to
        ! contribute for its response to describe curvature rather than the error in
        ! it: the response is an exactly symmetric function of the step only for a
        ! vanishing step, so the antisymmetric part of the transformed response matrix
        ! Q^T Z collects both the numerical noise of the response and the finite-step
        ! error of a curvature that changes along the step, neither of which a
        ! symmetric fit can use; every column of Z is divided by the residual norm of
        ! its own direction, so multiplying the typical size of column k of that part
        ! by that residual norm undoes the division, and the median over the
        ! directions makes a single error scale of it; a direction shorter than that
        ! scale divided by the typical curvature contributes more amplified error than
        ! curvature
        !
        real(rp), intent(in) :: steps(:, :), responses(:, :)
        logical, intent(in) :: keep(:)
        real(rp) :: min_residual

        integer(ip) :: n_accepted, i, k, m
        integer(ip), allocatable :: map(:)
        real(rp), allocatable :: chol(:, :), q(:, :), z(:, :), a(:, :), &
                                 pair_asymmetries(:), asymmetries(:), residuals(:), &
                                 curvatures(:)
        real(rp) :: response_error, curvature
        external :: dgemm

        ! a direction needs at least two partners for a typical asymmetry of its own
        ! column to mean anything
        min_residual = 0.0_rp
        call factorize_history(steps, chol, map, n_accepted, keep)
        if (n_accepted < 3) return

        ! transformed response matrix Q^T Z in the orthonormalized step basis
        q = rebase_dirs(steps, map, chol)
        z = rebase_dirs(responses, map, chol)
        allocate(a(n_accepted, n_accepted))
        call dgemm("T", "N", n_accepted, n_accepted, size(q, 1, kind=ip), 1.0_rp, q, &
                   size(q, 1, kind=ip), z, size(z, 1, kind=ip), 0.0_rp, a, n_accepted)

        ! per-direction asymmetry, residual norm and curvature
        allocate(pair_asymmetries(n_accepted - 1), asymmetries(n_accepted), &
                 residuals(n_accepted), curvatures(n_accepted))
        do k = 1, n_accepted
            m = 0
            do i = 1, n_accepted
                if (i == k) cycle
                m = m + 1
                pair_asymmetries(m) = abs(a(i, k) - a(k, i))
            end do
            asymmetries(k) = median(pair_asymmetries)
            residuals(k) = abs(chol(k, k))
            curvatures(k) = abs(a(k, k))
        end do

        ! undo the amplification of every column to get the error of the response
        ! itself, and take the curvature the matrix describes at its typical size
        response_error = median(asymmetries * residuals)
        curvature = median(curvatures)

        ! residual norm at which the amplified error reaches the curvature scale
        if (curvature > 0.0_rp) min_residual = response_error / curvature
        deallocate(chol, map, q, z, a, pair_asymmetries, asymmetries, residuals, &
                   curvatures)

    end function resolvable_residual

    subroutine factorize_history(vecs, chol, map, n_accepted, keep, min_residual)
        !
        ! this subroutine performs a pivoted, rank-revealing Cholesky factorization of
        ! the Gram matrix of a set of history vectors, taken on the Gram normalized to
        ! a unit diagonal so that the pivot order and the rank tolerance test genuine
        ! linear dependence rather than step length: the secant conditions are
        ! homogeneous, so rescaling a column with its response leaves the approximate
        ! Hessian unchanged, while a length-driven pivot would discard the newest,
        ! shortest steps as dependent purely for being small; the normalization is
        ! undone in the returned factor, so together with the map back to original
        ! indices any other history can be re-expressed in the same orthonormalized
        ! basis by a single triangular solve, without forming a metric inverse
        !
        ! keep optionally restricts which columns a caller trusts; omitting it uses the 
        ! whole history, which is right for a response that is exact at any distance
        !
        ! min_residual optionally demands a shortest residual norm a column has to
        ! contribute to be accepted, which is right for a response whose own error
        ! limits how short a direction can be and still be resolved; the Gram is then 
        ! left unnormalized so that its pivots are comparable to that absolute length, 
        ! which also makes the pivoting rank the columns by the length they contribute 
        ! rather than by independence, since that length is what the tolerance bounds, 
        ! and leaves the returned factor already carrying the norms; a non-positive
        ! value demands nothing
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: vecs(:, :)
        real(rp), intent(out), allocatable :: chol(:, :)
        integer(ip), intent(out), allocatable :: map(:)
        integer(ip), intent(out) :: n_accepted
        logical, intent(in), optional :: keep(:)
        real(rp), intent(in), optional :: min_residual

        integer(ip) :: vec_len, n_vecs, i, j, info
        logical :: normalize
        real(rp) :: tol
        real(rp), allocatable :: metric(:, :), work(:), col_norm(:)
        logical, allocatable :: usable(:)
        integer(ip), allocatable :: piv(:)
        external :: dsyrk, dpstrf

        vec_len = size(vecs, 1)
        n_vecs = size(vecs, 2)

        n_accepted = 0
        if (n_vecs == 0) then
            allocate(chol(0, 0), map(0))
            return
        end if

        ! decide whether the pivoting ranks the columns by independence or by the
        ! length they contribute
        normalize = .true.
        if (present(min_residual)) normalize = min_residual <= 0.0_rp

        ! generate the Gram matrix
        allocate(metric(n_vecs, n_vecs))
        metric = 0.0_rp
        call dsyrk("U", "T", n_vecs, vec_len, 1.0_rp, vecs, vec_len, 0.0_rp, metric, &
                   n_vecs)

        ! determine the norm of every step
        allocate(col_norm(n_vecs), usable(n_vecs))
        do j = 1, n_vecs
            col_norm(j) = sqrt(metric(j, j))
        end do

        ! check whether any steps are numerically vanishing
        usable = col_norm > numerical_zero * maxval(col_norm)
        if (present(keep)) usable = usable .and. keep

        ! normalize the usable steps to a unit diagonal,when requested, so that the 
        ! pivoting ranks them by independence rather than by magnitude, and zero out 
        ! the rest so that they are rejected rather than admitted with a zero pivot
        do j = 1, n_vecs
            do i = 1, j - 1
                if (.not. (usable(i) .and. usable(j))) then
                    metric(i, j) = 0.0_rp
                else if (normalize) then
                    metric(i, j) = metric(i, j) / (col_norm(i) * col_norm(j))
                end if
            end do
            if (.not. usable(j)) then
                metric(j, j) = 0.0_rp
            else if (normalize) then
                metric(j, j) = 1.0_rp
            end if
        end do

        ! tolerance for linear dependencies, relative to the unit diagonal, or the
        ! squared shortest resolvable residual on the unnormalized Gram
        if (normalize) then
            tol = n_vecs * numerical_zero
        else
            tol = max(min_residual**2, n_vecs * numerical_zero * maxval(col_norm)**2)
        end if

        ! perform pivoted rank-revealing Cholesky
        allocate(piv(n_vecs), work(2 * n_vecs))
        call dpstrf("U", n_vecs, metric, n_vecs, piv, n_accepted, tol, work, info)

        ! shrink to the accepted block only and undo the normalization by scaling,
        ! where it was applied
        allocate(chol(n_accepted, n_accepted), map(n_accepted))
        chol = metric(1:n_accepted, 1:n_accepted)
        map = piv(1:n_accepted)
        if (normalize) then
            do j = 1, n_accepted
                chol(:, j) = chol(:, j) * col_norm(map(j))
            end do
        end if
        deallocate(metric, piv, work, col_norm, usable)

    end subroutine factorize_history

    function rebase_dirs(dirs, map, chol) result(rebased)
        !
        ! this function re-expresses a set of packed history-direction columns in the 
        ! orthonormalized basis by selecting and reordering the accepted columns via 
        ! map, then applying a single triangular solve by the Cholesky factor, 
        ! equivalent to right-multiplying by the inverse of the (implicit) 
        ! upper-triangular factor without ever forming that inverse explicitly
        !
        real(rp), intent(in) :: dirs(:, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: rebased(:, :)

        integer(ip) :: n_param, n_accepted, j
        external :: dtrsm

        n_param = size(dirs, 1)
        n_accepted = size(map)
        allocate(rebased(n_param, n_accepted))
        do j = 1, n_accepted
            rebased(:, j) = dirs(:, map(j))
        end do

        if (n_accepted > 0) call dtrsm("R", "U", "N", "N", n_param, n_accepted, &
                                       1.0_rp, chol, n_accepted, rebased, n_param)

    end function rebase_dirs

    function congruence_transform(a, map, chol) result(a_tilde)
        !
        ! this function applies the congruence transformation A -> R^-T (P A P^T) R^-1, 
        ! where P selects and reorders the rows/columns of A according to map and R is 
        ! the Cholesky factor
        !
        real(rp), intent(in) :: a(:, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: a_tilde(:, :)

        integer(ip) :: n_accepted, i, j
        external :: dtrsm

        n_accepted = size(map)
        allocate(a_tilde(n_accepted, n_accepted))
        do j = 1, n_accepted
            do i = 1, n_accepted
                a_tilde(i, j) = a(map(i), map(j))
            end do
        end do

        if (n_accepted > 0) then
            call dtrsm("R", "U", "N", "N", n_accepted, n_accepted, 1.0_rp, chol, &
                       n_accepted, a_tilde, n_accepted)
            call dtrsm("L", "U", "T", "N", n_accepted, n_accepted, 1.0_rp, chol, &
                       n_accepted, a_tilde, n_accepted)
        end if

    end function congruence_transform

    subroutine combine_channels(chol1, map1, chol2, map2, n_offset, chol_comb, map_comb)
        !
        ! this subroutine assembles a block-diagonal Cholesky factor and concatenated, 
        ! offset index map from two independent per-channel history factorizations, 
        ! so that a single rebase or congruence-transform call can act on data that 
        ! combines both channels
        !
        real(rp), intent(in) :: chol1(:, :), chol2(:, :)
        integer(ip), intent(in) :: map1(:), map2(:), n_offset
        real(rp), intent(out), allocatable :: chol_comb(:, :)
        integer(ip), intent(out), allocatable :: map_comb(:)

        integer(ip) :: n1, n2, n_total

        n1 = size(map1)
        n2 = size(map2)
        n_total = n1 + n2
        allocate(chol_comb(n_total, n_total), map_comb(n_total))
        chol_comb = 0.0_rp
        chol_comb(1:n1, 1:n1) = chol1
        chol_comb(n1 + 1:n_total, n1 + 1:n_total) = chol2
        map_comb(1:n1) = map1
        map_comb(n1 + 1:n_total) = n_offset + map2

    end subroutine combine_channels

    subroutine get_ms_a_inv(dm_diff, v_diff, map, chol, a_inv, settings, error)
        !
        ! this subroutine computes the pseudoinverse multisecant SR1 matrix for a
        ! single part of a coupled potential-difference response, from the symmetrized
        ! A of that part congruence-transformed into its own orthonormalized S-basis
        !
        use opentrustregion, only: symm_mat_diag

        real(rp), intent(in) :: dm_diff(:, :, :, :), v_diff(:, :, :, :)
        integer(ip), intent(in) :: map(:)
        real(rp), intent(in) :: chol(:, :)
        real(rp), intent(out), allocatable :: a_inv(:, :)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, n_accepted
        real(rp), allocatable :: a(:, :), a_tilde(:, :), eig_vecs(:, :), eig_vals(:), &
                                 eig_vals_inv(:)
        real(rp), allocatable :: y_gram(:, :)
        real(rp) :: eig_val_thresh

        ! initialize error flag
        error = 0

        ! handle empty history
        n_dm = size(dm_diff, 4)
        n_accepted = size(map)
        if (n_dm == 0 .or. n_accepted == 0) then
            allocate(a_inv(n_accepted, n_accepted))
            a_inv = 0.0_rp
            return
        end if

        ! build and symmetrize A = S^T Y
        allocate(a(n_dm, n_dm))
        call build_a_part(dm_diff, v_diff, a)

        ! congruence-transform to the orthonormalized, rank-independent S-basis
        a_tilde = congruence_transform(a, map, chol)
        deallocate(a)

        ! perform spectral decomposition
        allocate(eig_vecs(n_accepted, n_accepted), eig_vals(n_accepted))
        call symm_mat_diag(a_tilde, eig_vals, eig_vecs, settings, error)
        if (error /= 0) return

        ! construct inverse, discarding only eigenvalues at the level of numerical 
        ! noise relative to the largest one
        allocate(eig_vals_inv(n_accepted))
        eig_val_thresh = eig_val_noise_factor * maxval(abs(eig_vals)) * epsilon(1.0_rp)
        eig_vals_inv = truncated_eigval_inv(eig_vals, eig_val_thresh)

        ! discard directions failing the MS-SR1 skipping criterion
        y_gram = response_gram(v_diff, map, chol)
        call apply_ms_sr1_skip(eig_vals, eig_vecs, y_gram, eig_vals_inv)
        deallocate(y_gram)

        ! reassemble the pseudoinverse
        a_inv = spectral_to_dense(eig_vecs, eig_vals_inv)
        deallocate(a_tilde, eig_vecs, eig_vals, eig_vals_inv)

    end subroutine get_ms_a_inv

    subroutine get_ms_a_inv_os_linear(dm_diff, v_same_spin_diff, v_opposite_spin_diff, &
                                      map, chol, a_inv, n_ao, settings, error)
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
        integer(ip), intent(in) :: map(:)
        real(rp), intent(in) :: chol(:, :)
        real(rp), intent(out), allocatable :: a_inv(:, :)
        integer(ip), intent(in) :: n_ao
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, n_accepted, i, k
        real(rp), allocatable :: a(:, :), a_tilde(:, :), eig_vecs(:, :), eig_vals(:), &
                                 eig_vals_inv(:)
        real(rp), allocatable :: y_gram(:, :)
        real(rp) :: eig_val_thresh

        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! handle empty history
        n_dm = size(dm_diff, 4)
        n_accepted = size(map)
        if (n_dm == 0 .or. n_accepted == 0) then
            allocate(a_inv(n_accepted, n_accepted))
            a_inv = 0.0_rp
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
                a(n_dm + i, n_dm + k) = ddot(n_ao * n_ao, dm_diff(:, :, 2, i), 1_ip, &
                                             v_same_spin_diff(:, :, 2, k), 1_ip)
            end do
        end do

        ! congruence-transform to the orthonormalized, rank-independent S-basis
        a_tilde = congruence_transform(a, map, chol)
        deallocate(a)

        ! enforce the symmetry A has in exact arithmetic, which this routine does
        ! after the transformation rather than on A itself as the routines building A
        ! from its parts do
        a_tilde = 0.5_rp * (a_tilde + transpose(a_tilde))

        ! perform spectral decomposition
        allocate(eig_vecs(n_accepted, n_accepted), eig_vals(n_accepted))
        call symm_mat_diag(a_tilde, eig_vals, eig_vecs, settings, error)
        if (error /= 0) return

        ! construct inverse, discarding only eigenvalues at the level of numerical 
        ! noise relative to the largest one
        allocate(eig_vals_inv(n_accepted))
        eig_val_thresh = eig_val_noise_factor * maxval(abs(eig_vals)) * epsilon(1.0_rp)
        eig_vals_inv = truncated_eigval_inv(eig_vals, eig_val_thresh)

        ! discard directions failing the multisecant SR1 skipping criterion
        y_gram = response_gram_os_linear(v_same_spin_diff, v_opposite_spin_diff, n_ao, &
                                         map, chol)
        call apply_ms_sr1_skip(eig_vals, eig_vecs, y_gram, eig_vals_inv)
        deallocate(y_gram)

        ! reassemble the pseudoinverse
        a_inv = spectral_to_dense(eig_vecs, eig_vals_inv)
        deallocate(a_tilde, eig_vecs, eig_vals, eig_vals_inv)

    end subroutine get_ms_a_inv_os_linear

    function response_gram(v_diff, map, chol) result(y_gram)
        !
        ! this function returns the Gram matrix of the response history rebased into
        ! the orthonormalized S-basis, so that the squared norm of the response along
        ! an eigenvector v of the congruence-transformed A is v^T y_gram v
        !
        real(rp), intent(in) :: v_diff(:, :, :, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: y_gram(:, :)

        integer(ip) :: n_dm, flat_len
        real(rp), allocatable :: gram(:, :)
        external :: dgemm

        n_dm = size(v_diff, 4)
        flat_len = size(v_diff, 1) * size(v_diff, 2) * size(v_diff, 3)
        allocate(gram(n_dm, n_dm))
        call dgemm("T", "N", n_dm, n_dm, flat_len, 1.0_rp, v_diff, flat_len, v_diff, &
                   flat_len, 0.0_rp, gram, n_dm)
        y_gram = congruence_transform(gram, map, chol)
        deallocate(gram)

    end function response_gram

    function response_gram_os_linear(v_same_spin_diff, v_opposite_spin_diff, n_ao, &
                                     map, chol) result(y_gram)
        !
        ! this function returns the Gram matrix of the open-shell linear response
        ! history rebased into the orthonormalized S-basis; the same-spin and
        ! opposite-spin potentials are interleaved exactly as the rows of A pair them, 
        ! so that column k of the implicit response matrix stacks the alpha and beta 
        ! blocks that A(:, k) contracts against
        !
        real(rp), intent(in) :: v_same_spin_diff(:, :, :, :), &
                                v_opposite_spin_diff(:, :, :, :), chol(:, :)
        integer(ip), intent(in) :: n_ao, map(:)
        real(rp), allocatable :: y_gram(:, :)

        integer(ip) :: n_dm, n_ao2, k
        real(rp), allocatable :: y_full(:, :), gram(:, :)
        external :: dgemm

        n_dm = size(v_same_spin_diff, 4)
        n_ao2 = n_ao * n_ao
        allocate(y_full(2 * n_ao2, 2 * n_dm), gram(2 * n_dm, 2 * n_dm))
        do k = 1, n_dm
            y_full(:n_ao2, k) = reshape(v_same_spin_diff(:, :, 1, k), [n_ao2])
            y_full(n_ao2 + 1:, k) = reshape(v_opposite_spin_diff(:, :, 2, k), [n_ao2])
            y_full(:n_ao2, n_dm + k) = reshape(v_opposite_spin_diff(:, :, 1, k), &
                                               [n_ao2])
            y_full(n_ao2 + 1:, n_dm + k) = reshape(v_same_spin_diff(:, :, 2, k), &
                                                   [n_ao2])
        end do
        call dgemm("T", "N", 2 * n_dm, 2 * n_dm, 2 * n_ao2, 1.0_rp, y_full, 2 * n_ao2, &
                   y_full, 2 * n_ao2, 0.0_rp, gram, 2 * n_dm)
        y_gram = congruence_transform(gram, map, chol)
        deallocate(y_full, gram)

    end function response_gram_os_linear

    subroutine apply_ms_sr1_skip(eig_vals, eig_vecs, y_gram, eig_vals_inv)
        !
        ! this subroutine applies the multisecant analogue of the SR1 skipping
        ! criterion: an eigenvalue of the congruence-transformed A plays the role of 
        ! the classical denominator s^T (y - H s), and a direction is discarded when it 
        ! is small against the norm of the response it divides; the direction has unit 
        ! norm in the orthonormalized S-basis, so the ratio is a dimensionless angle
        !
        real(rp), intent(in) :: eig_vals(:), eig_vecs(:, :), y_gram(:, :)
        real(rp), intent(inout) :: eig_vals_inv(:)

        integer(ip) :: i, n
        real(rp) :: y_norm
        real(rp), allocatable :: temp(:)
        real(rp), external :: ddot
        external :: dgemv

        n = size(eig_vals)
        allocate(temp(n))
        do i = 1, n
            call dgemv("N", n, n, 1.0_rp, y_gram, n, eig_vecs(:, i), 1_ip, 0.0_rp, &
                       temp, 1_ip)
            y_norm = sqrt(max(ddot(n, eig_vecs(:, i), 1_ip, temp, 1_ip), 0.0_rp))
            if (abs(eig_vals(i)) < ms_sr1_skip_thresh * y_norm) eig_vals_inv(i) = 0.0_rp
        end do
        deallocate(temp)

    end subroutine apply_ms_sr1_skip

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

    function density_in_history(dm_oao) result(in_history)
        !
        ! this function reports whether a density is already held in the history
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: dm_oao(:, :, :)
        logical :: in_history

        integer(ip) :: i
        real(rp) :: dm_scale

        in_history = .false.
        if (.not. allocated(arh_object%dm_list)) return

        ! judge the difference against the size of the density, so that the test is a
        ! relative one
        dm_scale = max(maxval(abs(dm_oao)), numerical_zero)
        do i = 1, size(arh_object%dm_list, 4)
            if (maxval(abs(arh_object%dm_list(:, :, :, i) - dm_oao)) <= &
                numerical_zero * dm_scale) then
                in_history = .true.
                return
            end if
        end do

    end function density_in_history

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
