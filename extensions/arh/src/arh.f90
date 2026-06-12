! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh

    use opentrustregion, only: rp, ip, kw_len, obj_func_type, update_orbs_type, &
                               hess_x_type, project_type
    use otr_oao, only: oao_settings_type, default_oao_settings, oao_object, &
                       get_energy_3d_type, get_energy_2d_type, update_dm_2d_type

    implicit none

    type, extends(oao_settings_type) :: arh_settings_type
        character(kw_len) :: arh_type
    contains
        procedure :: init => init_arh_settings
    end type arh_settings_type

    type(arh_settings_type), parameter :: default_arh_settings = &
        arh_settings_type(oao_settings_type = default_oao_settings, &
                          arh_type = "symmetric")

    ! define setting options
    character(kw_len), parameter :: arh_types(3) = &
            [character(len=kw_len) :: "standard", "symmetric", "multisecant_psb"]

    abstract interface
        subroutine update_dm_jk_type(dm, energy, fock, coulomb, exchange, &
                                     get_response_funptr, error)
            use otr_oao, only: get_response_3d_type
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), target :: fock(:, :, :), coulomb(:, :, :), &
                                             exchange(:, :, :)
            procedure(get_response_3d_type), intent(out), pointer :: get_response_funptr
            integer(ip), intent(out) :: error
        end subroutine update_dm_jk_type
    end interface

    type :: arh_type
        type(arh_settings_type) :: settings
        integer(ip), pointer :: n_ao => null(), n_param => null(), n_particle => null()
        real(rp), pointer :: dm_ao(:, :, :) => null(), s_inv_sqrt(:, :) => null(), &
                             dm_oao(:, :, :) => null(), fock_oo(:, :, :) => null(), &
                             fock_vv(:, :, :) => null()
        real(rp), allocatable :: fock_oao(:, :, :), same_v_oao(:, :, :), &
                                 opposite_v_oao(:, :, :), metric_eigvals(:, :), &
                                 metric_eigvecs(:, :, :), dm_list(:, :, :, :), &
                                 fock_list(:, :, :, :), same_v_list(:, :, :, :), &
                                 opposite_v_list(:, :, :, :), dm_diff(:, :, :, :), &
                                 fock_diff(:, :, :, :), same_v_diff(:, :, :, :), &
                                 opposite_v_diff(:, :, :, :)
        procedure(get_energy_3d_type), pointer, nopass :: get_energy => null()
        procedure(update_dm_jk_type), pointer, nopass :: update_dm_jk => null()
        procedure(get_energy_2d_type), pointer, nopass :: get_energy_2d => null()
        procedure(update_dm_2d_type), pointer, nopass :: update_dm_2d => null()
    end type arh_type

    ! global variables
    type(arh_type), allocatable :: arh_object

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_arh_closed_shell_ptr => &
        update_orbs_arh_closed_shell
    procedure(update_orbs_type), pointer :: update_orbs_arh_open_shell_ptr => &
        update_orbs_arh_open_shell
    procedure(hess_x_type), pointer :: hess_x_arh_ptr => hess_x_arh

    ! define module procedures for different spin cases
    interface arh_factory
        module procedure arh_factory_closed_shell, arh_factory_open_shell
    end interface arh_factory

    contains

    subroutine arh_factory_closed_shell(dm_ao, ao_overlap, n_particle, n_ao, &
                                        get_energy, update_dm, obj_func_arh_funptr, &
                                        update_orbs_arh_funptr, project_arh_funptr, &
                                        error, settings)
        !
        ! this function returns a modified ARH orbital updating function for the 
        ! closed-shell case
        !
        use otr_oao, only: get_energy_2d_type, update_dm_2d_type, obj_func_oao, &
                           project_oao

        real(rp), intent(in) :: dm_ao(:, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_2d_type), intent(in), pointer :: get_energy
        procedure(update_dm_2d_type), intent(in), pointer :: update_dm
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize error flag
        error = 0

        ! call common setup
        call arh_factory_common(reshape(dm_ao, [n_ao, n_ao, 1]), ao_overlap, &
                                n_particle, n_ao, error, settings)
        if (error /= 0) return

        ! set pointers to functions
        arh_object%get_energy_2d => get_energy
        arh_object%update_dm_2d => update_dm

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_oao
        update_orbs_arh_funptr => update_orbs_arh_closed_shell
        project_arh_funptr => project_oao

    end subroutine arh_factory_closed_shell

    subroutine arh_factory_open_shell(dm_ao, ao_overlap, n_particle, n_ao, get_energy, &
                                      update_dm_jk, obj_func_arh_funptr, &
                                      update_orbs_arh_funptr, project_arh_funptr, &
                                      error, settings)
        !
        ! this function returns a modified ARH orbital updating function for the 
        ! open-shell case
        !
        use otr_oao, only: get_energy_3d_type, obj_func_oao, project_oao

        real(rp), intent(in) :: dm_ao(:, :, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_3d_type), intent(in), pointer :: get_energy
        procedure(update_dm_jk_type), intent(in), pointer :: update_dm_jk
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
        arh_object%get_energy => get_energy
        arh_object%update_dm_jk => update_dm_jk

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_oao
        update_orbs_arh_funptr => update_orbs_arh_open_shell
        project_arh_funptr => project_oao

    end subroutine arh_factory_open_shell

    subroutine arh_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        !
        ! this function performs common ARH initialization operations
        !
        use otr_oao, only: oao_factory_common, oao_object

        real(rp), intent(in) :: dm_ao(:, :, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! call common OAO setup
        call oao_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)

        ! perform sanity check
        call arh_sanity_check(settings, error)
        if (error /= 0) return

        ! allocate objects
        if (.not. allocated(arh_object)) allocate(arh_object)

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
                              """standard"" (standard ARH), ""symmetric"" (simple "// &
                              "symmetrized ARH), and ""multisecant_psb"" "// &
                              "(multisecant PSB version of ARH).", verbosity_error, &
                              .true.)
            error = 1
            return
        end if

    end subroutine arh_sanity_check

    subroutine update_orbs_arh_closed_shell(kappa, func, grad, h_diag, hess_x_funptr, &
                                            error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation 
        ! in the OAO basis and the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall for the closed-shell case
        !
        use opentrustregion, only: hess_x_type
        use otr_oao, only: get_response_2d_type, rotate_dm_ao, &
                           symmetric_transformation, calculate_grad_h_diag

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list
        real(rp), allocatable :: fock_ao(:, :, :)
        procedure(get_response_2d_type), pointer :: get_response_2d_funptr

        external :: dgemm

        ! number of AOs
        n_ao = arh_object%n_ao

        ! number of particles
        n_particle = arh_object%n_particle

        ! update list of density and Fock matrices
        if (allocated(arh_object%dm_list)) then
            call append(arh_object%dm_list, arh_object%dm_oao)
            call append(arh_object%fock_list, arh_object%fock_oao)
        else
            allocate(arh_object%dm_list(n_ao, n_ao, n_particle, 0), &
                     arh_object%fock_list(n_ao, n_ao, n_particle, 0))
        end if

        ! rotate density matrix
        call rotate_dm_ao(kappa, n_particle, n_ao, &
                          arh_object%settings%restricted, arh_object%dm_ao, error, &
                          arh_object%dm_oao)
        if (error /= 0) return

        ! get energyand Fock matrix
        allocate(fock_ao(n_ao, n_ao, n_particle))
        call arh_object%update_dm_2d(arh_object%dm_ao(:, :, 1), func, &
                                     fock_ao(:, :, 1), get_response_2d_funptr, error)
        if (error /= 0) then
            deallocate(fock_ao)
            return
        end if

        ! transform Fock matrix to OAO basis
        arh_object%fock_oao = &
            symmetric_transformation(arh_object%s_inv_sqrt, fock_ao)
        deallocate(fock_ao)

        ! calculate gradient and Hessian diagonal
        call calculate_grad_h_diag(arh_object%dm_oao, arh_object%fock_oao, n_particle, &
                                   n_ao, arh_object%settings%restricted, grad, h_diag, &
                                   arh_object%fock_oo, arh_object%fock_vv)

        ! construct and diagonalize ARH metric
        call get_arh_metric(arh_object%dm_list, arh_object%dm_oao, &
                            arh_object%metric_eigvals, arh_object%metric_eigvecs, &
                            arh_object%settings, error)
        if (error /= 0) return

        ! prepare differences for response part of ARH Hessian
        n_list = size(arh_object%dm_list, 4)
        if (allocated(arh_object%dm_diff)) deallocate(arh_object%dm_diff, &
                                                      arh_object%fock_diff)
        allocate(arh_object%dm_diff(n_ao, n_ao, n_particle, n_list), &
                 arh_object%fock_diff(n_ao, n_ao, n_particle, n_list))
        do i = 1, n_list
            arh_object%dm_diff(:, :, :, i) = arh_object%dm_list(:, :, :, i) - &
                                             arh_object%dm_oao
            arh_object%fock_diff(:, :, :, i) = arh_object%fock_list(:, :, :, i) - &
                                               arh_object%fock_oao
        end do

        ! define pointer to ARH Hessian linear transformation function
        hess_x_funptr => hess_x_arh
        
    end subroutine update_orbs_arh_closed_shell

    subroutine update_orbs_arh_open_shell(kappa, func, grad, h_diag, hess_x_funptr, &
                                          error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation 
        ! in the OAO basis and the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall for the open-shell case
        !
        use opentrustregion, only: hess_x_type
        use otr_oao, only: get_response_3d_type, rotate_dm_ao, &
                           symmetric_transformation, calculate_grad_h_diag

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list
        real(rp), allocatable :: fock_ao(:, :, :), fock_oao(:, :, :), &
                                 coulomb_ao(:, :, :), exchange_ao(:, :, :), &
                                 same_v_ao(:, :, :), opposite_v_ao(:, :, :)
        procedure(get_response_3d_type), pointer :: get_response_3d_funptr
                                 
        external :: dgemm

        ! number of AOs
        n_ao = arh_object%n_ao

        ! number of particles
        n_particle = arh_object%n_particle

        ! update list of density and potential matrices
        if (allocated(arh_object%dm_list)) then
            call append(arh_object%dm_list, arh_object%dm_oao)
            call append(arh_object%same_v_list, arh_object%same_v_oao)
            call append(arh_object%opposite_v_list, arh_object%opposite_v_oao)
        else
            allocate(arh_object%dm_list(n_ao, n_ao, n_particle, 0), &
                     arh_object%same_v_list(n_ao, n_ao, n_particle, 0), &
                     arh_object%opposite_v_list(n_ao, n_ao, n_particle, 0))
        end if

        ! rotate density matrix
        call rotate_dm_ao(kappa, n_particle, n_ao, arh_object%settings%restricted, &
                          arh_object%dm_ao, error, arh_object%dm_oao)
        if (error /= 0) return

        ! get energy, Fock matrix, and same and opposite spin potentials
        allocate(fock_ao(n_ao, n_ao, n_particle), coulomb_ao(n_ao, n_ao, n_particle), &
                 exchange_ao(n_ao, n_ao, n_particle))
        call arh_object%update_dm_jk(arh_object%dm_ao, func, fock_ao, coulomb_ao, &
                                     exchange_ao, get_response_3d_funptr, error)
        if (error /= 0) then
            deallocate(fock_ao, coulomb_ao, exchange_ao)
            return
        end if

        ! transform Fock matrix to OAO basis
        fock_oao = symmetric_transformation(arh_object%s_inv_sqrt, fock_ao)
        deallocate(fock_ao)

        ! get same and opposite spin potentials
        allocate(same_v_ao(n_ao, n_ao, n_particle), &
                 opposite_v_ao(n_ao, n_ao, n_particle))
        same_v_ao = coulomb_ao - exchange_ao
        opposite_v_ao = coulomb_ao
        deallocate(coulomb_ao, exchange_ao)

        ! transform same and opposite spin potentials to OAO basis
        arh_object%same_v_oao = symmetric_transformation(arh_object%s_inv_sqrt, &
                                                         same_v_ao)
        arh_object%opposite_v_oao = symmetric_transformation(arh_object%s_inv_sqrt, &
                                                             opposite_v_ao)
        deallocate(same_v_ao, opposite_v_ao)

        ! calculate gradient and Hessian diagonal
        call calculate_grad_h_diag(arh_object%dm_oao, fock_oao, n_particle, n_ao, &
                                   arh_object%settings%restricted, grad, h_diag, &
                                   arh_object%fock_oo, arh_object%fock_vv)
        deallocate(fock_oao)

        ! construct and diagonalize ARH metric
        call get_arh_metric(arh_object%dm_list, arh_object%dm_oao, &
                            arh_object%metric_eigvals, arh_object%metric_eigvecs, &
                            arh_object%settings, error)
        if (error /= 0) return

        ! prepare differences for response part of ARH Hessian
        n_list = size(arh_object%dm_list, 4)
        if (allocated(arh_object%dm_diff)) deallocate(arh_object%dm_diff, &
                                                      arh_object%same_v_diff, &
                                                      arh_object%opposite_v_diff)
        allocate(arh_object%dm_diff(n_ao, n_ao, n_particle, n_list), &
                 arh_object%same_v_diff(n_ao, n_ao, n_particle, n_list), &
                 arh_object%opposite_v_diff(n_ao, n_ao, n_particle, n_list))
        do i = 1, n_list
            arh_object%dm_diff(:, :, :, i) = arh_object%dm_list(:, :, :, i) - &
                                             arh_object%dm_oao
            arh_object%same_v_diff(:, :, :, i) = arh_object%same_v_list(:, :, :, i) - &
                                                 arh_object%same_v_oao
            arh_object%opposite_v_diff(:, :, :, i) = &
                arh_object%opposite_v_list(:, :, :, i) - arh_object%opposite_v_oao
        end do

        ! define pointer to ARH Hessian linear transformation function
        hess_x_funptr => hess_x_arh
        
    end subroutine update_orbs_arh_open_shell

    subroutine hess_x_arh(x, hess_x, error)
        !
        ! this function defines the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall
        !
        use otr_oao, only: unpack_asymm, project, pack_asymm

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, n_param, i
        real(rp), allocatable :: x_full(:, :, :), hess_x_full(:, :, :)

        external :: dgemm

        ! initialize error flag
        error = 0

        ! number of AOs
        n_ao = arh_object%n_ao

        ! number of particles
        n_particle = arh_object%n_particle

        ! number of parameters
        n_param = arh_object%n_param

        ! unpack trial vector
        x_full = unpack_asymm(x, n_particle, n_ao, arh_object%settings%restricted)

        ! for ROHF, we must explicitly project the trial vector to the valid [c-v], 
        ! [o-v], and [c-o] rotation blocks. Because we reuse UHF machinery where 
        ! alpha-occ = [c+o] and beta-occ = [c], a unified spatial trial vector that is 
        ! unpacked into alpha and beta spin channels contains rotations that are 
        ! internal and thus redundant in the ROHF context. Without this projection, 
        ! the alpha Hessian vector product would leak non-zero values into 
        ! the [c-c] block and the beta Hessian vector product into the [v-v] blocks. 
        ! This projection thus prevents singular dimensions in the solver.
        if (arh_object%settings%restricted .and. arh_object%n_particle > 1) &
            x_full = project(x_full, arh_object%dm_oao)

        ! get static part
        allocate(hess_x_full(n_ao, n_ao, n_particle))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, arh_object%fock_vv(:, :, i) &
                       - arh_object%fock_oo(:, :, i), n_ao, x_full(:, :, i), n_ao, &
                       0.0_rp, hess_x_full(:, :, i), n_ao)
            hess_x_full(:, :, i) = hess_x_full(:, :, i) - &
                                   transpose(hess_x_full(:, :, i))
        end do

        ! get response contributions
        if (n_particle == 1) then
            hess_x_full = hess_x_full + &
                get_response_contribution_closed_shell(arh_object%dm_oao, x_full, &
                                                       arh_object%dm_diff, &
                                                       arh_object%fock_diff, &
                                                       arh_object%metric_eigvals, &
                                                       arh_object%metric_eigvecs, &
                                                       n_ao, n_particle, &
                                                       arh_object%settings)
        else
            hess_x_full = hess_x_full + &
                get_response_contribution_open_shell(arh_object%dm_oao, x_full, &
                                                     arh_object%dm_diff, &
                                                     arh_object%same_v_diff, &
                                                     arh_object%opposite_v_diff, &
                                                     arh_object%metric_eigvals, &
                                                     arh_object%metric_eigvecs, n_ao, &
                                                     n_particle, arh_object%settings)
        end if
        if (error /= 0) return
        deallocate(x_full)

        ! pack Hessian linear transformation
        hess_x = pack_asymm(hess_x_full, size(hess_x), arh_object%settings%restricted)
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

    function get_response_contribution_closed_shell(dm_oao, x, dm_diff, fock_diff, &
                                                    metric_eigvals, metric_eigvecs, &
                                                    n_ao, n_particle, settings) &
                                                    result(response)
        !
        ! this function computes the response contribution to the ARH Hessian for the 
        ! closed-shell case
        !
        use opentrustregion, only: numerical_zero
        use otr_oao, only: project

        real(rp), intent(in) :: dm_oao(:, :, :), x(:, :, :), dm_diff(:, :, :, :), &
                                fock_diff(:, :, :, :), metric_eigvals(:, :), &
                                metric_eigvecs(:, :, :)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp) :: response(n_ao, n_ao, n_particle)
        type(arh_settings_type), intent(in) :: settings

        integer(ip) :: n_diff
        real(rp), allocatable :: delta_dm(:, :), s_delta_dm(:), t_s_delta_dm(:), &
                                 y_delta_dm(:), t_y_delta_dm(:), y_t_s_delta_dm(:, :), &
                                 sy_t_s_delta_dm(:), t_sy_t_s_delta_dm(:), val(:)
        real(rp) :: factor
        external :: dgemm, dgemv

        ! number of density matrix differences
        n_diff = size(dm_diff, 4)

        ! get factor for symmetrization
        if (settings%arh_type == "symmetric") then
            factor = 0.5_rp
        else
            factor = 1.0_rp
        end if

        ! calculate response contributions
        response = 0.0_rp
        if (n_diff > 0) then
            ! get current displacement as commutator of trial vector and current 
            ! density matrix
            allocate(delta_dm(n_ao, n_ao))
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, 1), n_ao, &
                       x(:, :, 1), n_ao, 0.0_rp, delta_dm, n_ao)
            delta_dm = delta_dm + transpose(delta_dm)

            ! get traces of density and Fock matrix differences with current 
            ! displacement
            allocate(s_delta_dm(n_diff))
            call dgemv("T", n_ao * n_ao, n_diff, 1.0_rp, dm_diff, n_ao * n_ao, &
                       delta_dm, 1, 0.0_rp, s_delta_dm, 1)
            if (settings%arh_type == "symmetric" .or. &
                settings%arh_type == "multisecant_psb") then
                allocate(y_delta_dm(n_diff))
                call dgemv("T", n_ao * n_ao, n_diff, 1.0_rp, fock_diff, n_ao * n_ao, &
                           delta_dm, 1, 0.0_rp, y_delta_dm, 1)
            end if
            deallocate(delta_dm)

            ! multiply pseudoinverse metric with traces of density matrix and Fock 
            ! matrix differences and current displacement
            allocate(t_s_delta_dm(n_diff))
            t_s_delta_dm = &
                multiply_with_inverse_metric(s_delta_dm, metric_eigvals(:, 1), &
                                             metric_eigvecs(:, :, 1))
            if (settings%arh_type == "symmetric" .or. &
                settings%arh_type == "multisecant_psb") then
                allocate(t_y_delta_dm(n_diff))
                t_y_delta_dm = &
                    multiply_with_inverse_metric(y_delta_dm, metric_eigvals(:, 1), &
                                                 metric_eigvecs(:, :, 1))
            end if

            ! contract with Fock matrix differences to get ARH contribution
            allocate(y_t_s_delta_dm(n_ao, n_ao))
            call dgemv("N", n_ao * n_ao, n_diff, 1.0_rp, fock_diff, n_ao * n_ao, &
                       t_s_delta_dm, 1, 0.0_rp, y_t_s_delta_dm, 1)
            deallocate(t_s_delta_dm)
            response(:, :, 1) = factor * y_t_s_delta_dm

            ! calculate multisecant PSB contribution
            if (settings%arh_type == "multisecant_psb") then
                ! get trace of density matrix differences with matrix
                allocate(sy_t_s_delta_dm(n_diff))
                call dgemv("T", n_ao * n_ao, n_diff, 1.0_rp, dm_diff, n_ao * n_ao, &
                           y_t_s_delta_dm, 1, 0.0_rp, sy_t_s_delta_dm, 1)
                deallocate(y_t_s_delta_dm)

                ! multiply pseudoinverse metric with vector
                allocate(t_sy_t_s_delta_dm(n_diff))
                t_sy_t_s_delta_dm = &
                    multiply_with_inverse_metric(sy_t_s_delta_dm, &
                                                 metric_eigvals(:, 1), &
                                                 metric_eigvecs(:, :, 1))
            else
                deallocate(y_t_s_delta_dm)
            end if

            ! contract with density matrix differences to get symmetric and multisecant 
            ! PSB contributions
            if (settings%arh_type == "symmetric" .or. &
                settings%arh_type == "multisecant_psb") then
                val = t_y_delta_dm
                deallocate(t_y_delta_dm)
                if (settings%arh_type == "multisecant_psb") then
                    val = val - t_sy_t_s_delta_dm
                    deallocate(t_sy_t_s_delta_dm)
                end if
                call dgemv("N", n_ao * n_ao, n_diff, factor, dm_diff, n_ao * n_ao, &
                           val, 1, 1.0_rp, response(:, :, 1), 1)
                deallocate(val)
            end if

            ! add only v-o and o-v contributions of ARH response part
            response = project(response, dm_oao)
        end if
        
    end function get_response_contribution_closed_shell

    function get_response_contribution_open_shell(dm_oao, x, dm_diff, same_v_diff, &
                                                  opposite_v_diff, metric_eigvals, &
                                                  metric_eigvecs, n_ao,  n_particle, &
                                                  settings) result(response)
        !
        ! this function computes the response contribution to the ARH Hessian for the 
        ! open-shell case
        !
        use opentrustregion, only: numerical_zero
        use otr_oao, only: project

        real(rp), intent(in) :: dm_oao(:, :, :), x(:, :, :), dm_diff(:, :, :, :), &
                                same_v_diff(:, :, :, :), opposite_v_diff(:, :, :, :), &
                                metric_eigvals(:, :), metric_eigvecs(:, :, :)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp) :: response(n_ao, n_ao, n_particle)
        type(arh_settings_type), intent(in) :: settings

        integer(ip) :: n_diff, i, j, k
        real(rp), allocatable :: delta_dm(:, :, :), s_delta_dm(:, :), &
                                 t_s_delta_dm(:, :), same_y_delta_dm(:), &
                                 t_same_y_delta_dm(:), opposite_y_delta_dm(:, :), &
                                 t_opposite_y_delta_dm(:, :), y_t_s_delta_dm(:, :), &
                                 sy_t_s_delta_dm(:), t_sy_t_s_delta_dm(:)
        real(rp) :: factor, val

        real(rp), external :: ddot

        ! number of density matrix differences
        n_diff = size(dm_diff, 4)

        ! get factor for symmetrization
        if (settings%arh_type == "symmetric") then
            factor = 0.5_rp
        else
            factor = 1.0_rp
        end if

        ! calculate response contributions
        if (n_diff > 0) then
            ! get current displacement as commutator of trial vector and current 
            ! density matrix
            allocate(delta_dm(n_ao, n_ao, n_particle))
            do j = 1, n_particle
                call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, j), n_ao, &
                           x(:, :, j), n_ao, 0.0_rp, delta_dm(:, :, j), n_ao)
                delta_dm(:, :, j) = delta_dm(:, :, j) + transpose(delta_dm(:, :, j))
            end do
            
            ! allocate necessary arrays
            allocate(s_delta_dm(n_diff, n_particle), t_s_delta_dm(n_diff, n_particle), &
                     y_t_s_delta_dm(n_ao, n_ao))
            if (settings%arh_type == "symmetric" .or. &
                settings%arh_type == "multisecant_psb") &
                allocate(same_y_delta_dm(n_diff), t_same_y_delta_dm(n_diff), &
                         opposite_y_delta_dm(n_diff, n_particle), &
                         t_opposite_y_delta_dm(n_diff, n_particle))
            if (settings%arh_type == "multisecant_psb") then
                allocate(sy_t_s_delta_dm(n_diff), t_sy_t_s_delta_dm(n_diff))
            end if

            ! get traces of density matrix differences with current displacement
            do i = 1, n_diff
                do j = 1, n_particle
                    s_delta_dm(i, j) = ddot(n_ao * n_ao, dm_diff(:, :, j, i), 1, &
                                            delta_dm(:, :, j), 1)
                end do
            end do

            ! multiply pseudoinverse metric with traces of density matrix differences 
            ! and current displacement
            do j = 1, n_particle
                t_s_delta_dm(:, j) = &
                    multiply_with_inverse_metric(s_delta_dm(:, j), &
                                                 metric_eigvals(:, j), &
                                                 metric_eigvecs(:, :, j))
            end do

            ! loop over particles
            do j = 1, n_particle
                if (settings%arh_type == "symmetric" .or. &
                    settings%arh_type == "multisecant_psb") then
                    ! get traces of potential differences with current displacement
                    do i = 1, n_diff
                        same_y_delta_dm(i) = ddot(n_ao * n_ao, &
                                                  same_v_diff(:, :, j, i), 1, &
                                                  delta_dm(:, :, j), 1)
                        do k = 1, n_particle
                            if (k == j) cycle
                            opposite_y_delta_dm(i, k) = &
                                ddot(n_ao * n_ao, opposite_v_diff(:, :, j, i), 1, &
                                     delta_dm(:, :, k), 1)
                        end do
                    end do
                
                    ! multiply pseudoinverse metric with traces of potential 
                    ! differences and current displacement
                    t_opposite_y_delta_dm = 0.0_rp
                    t_same_y_delta_dm = &
                        multiply_with_inverse_metric(same_y_delta_dm, &
                                                     metric_eigvals(:, j), &
                                                     metric_eigvecs(:, :, j))
                    do k = 1, n_particle
                        if (k == j) cycle
                        t_opposite_y_delta_dm(:, k) = &
                            multiply_with_inverse_metric(opposite_y_delta_dm(:, k), &
                                                         metric_eigvals(:, j), &
                                                         metric_eigvecs(:, :, j))
                    end do
                end if

                ! contract with same and opposite spin potential to get ARH contribution
                y_t_s_delta_dm = 0.0_rp
                do i = 1, n_diff
                    y_t_s_delta_dm = y_t_s_delta_dm + t_s_delta_dm(i, j) * &
                                     same_v_diff(:, :, j, i)
                    do k = 1, n_particle
                        if (k == j) cycle
                        y_t_s_delta_dm = y_t_s_delta_dm + t_s_delta_dm(i, k) * &
                                         opposite_v_diff(:, :, k, i)
                    end do
                end do
                response(:, :, j) = factor * y_t_s_delta_dm

                ! calculate multisecant PSB contribution
                if (settings%arh_type == "multisecant_psb") then
                    ! get trace of density matrix differences with matrix
                    do i = 1, n_diff
                        sy_t_s_delta_dm(i) = ddot(n_ao * n_ao, dm_diff(:, :, j, i), 1, &
                                                  y_t_s_delta_dm, 1)
                    end do

                    ! multiply pseudoinverse metric with vector
                    t_sy_t_s_delta_dm = &
                        multiply_with_inverse_metric(sy_t_s_delta_dm, &
                                                     metric_eigvals(:, j), &
                                                     metric_eigvecs(:, :, j))
                end if

                ! contract with density matrix differences to get symmetric and 
                ! multisecant PSB contributions
                if (settings%arh_type == "symmetric" .or. &
                    settings%arh_type == "multisecant_psb") then
                    do i = 1, n_diff
                        val = factor * t_same_y_delta_dm(i)
                        do k = 1, n_particle
                            if (k == j) cycle
                            val = val + factor * t_opposite_y_delta_dm(i, k)
                        end do
                        if (settings%arh_type == "multisecant_psb") &
                            val = val - factor * t_sy_t_s_delta_dm(i)
                        response(:, :, j) = response(:, :, j) + val * &
                                            dm_diff(:, :, j, i)
                    end do
                end if
            end do
            deallocate(delta_dm, s_delta_dm, t_s_delta_dm, y_t_s_delta_dm)
            if (settings%arh_type == "symmetric" .or. &
                settings%arh_type == "multisecant_psb") &
                deallocate(same_y_delta_dm, opposite_y_delta_dm, t_same_y_delta_dm, &
                           t_opposite_y_delta_dm)
            if (settings%arh_type == "multisecant_psb") then
                deallocate(sy_t_s_delta_dm, t_sy_t_s_delta_dm)
            end if

            ! add only v-o and o-v contributions of ARH response part
            response = project(response, dm_oao)
        end if

    end function get_response_contribution_open_shell

    function multiply_with_inverse_metric(vec, metric_eigvals, metric_eigvecs) &
        result(result_vec)
        !
        ! this function multiplies a vector with the pseudoinverse of the ARH metric by 
        ! transforming to the metric eigenbasis, applying the pseudoinverse and then 
        ! transforming back to the original basis
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: vec(:), metric_eigvals(:), metric_eigvecs(:, :)
        real(rp) :: result_vec(size(vec))

        integer(ip) :: n_diff, i
        real(rp), allocatable :: temp(:)

        n_diff = size(metric_eigvals, 1)

        allocate(temp(n_diff))
        call dgemv("T", n_diff, n_diff, 1.0_rp, metric_eigvecs, n_diff, vec, 1, &
                   0.0_rp, temp, 1)
        do i = 1, n_diff
            if (metric_eigvals(i) > numerical_zero) then
                temp(i) = temp(i) / metric_eigvals(i)
            else
                temp(i) = 0.0_rp
            end if
        end do
        call dgemv("N", n_diff, n_diff, 1.0_rp, metric_eigvecs, n_diff, temp, 1, &
                   0.0_rp, result_vec, 1)
        deallocate(temp)

    end function multiply_with_inverse_metric

    subroutine get_arh_metric(dm_list, dm_oao, eigvals, eigvecs, settings, error)
        !
        ! this function calculates the augmented Roothaan-Hall metric
        !
        use opentrustregion, only: symm_mat_diag

        real(rp), intent(in) :: dm_list(:, :, :, :), dm_oao(:, :, :)
        real(rp), intent(out), allocatable :: eigvals(:, :), eigvecs(:, :, :)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_particle, n_dm, i, j, k
        real(rp), allocatable :: metric(:, :), delta_i(:, :), delta_j(:, :)

        n_particle = size(dm_list, 3)
        n_dm = size(dm_list, 4)

        if (allocated(arh_object%metric_eigvals)) deallocate(arh_object%metric_eigvals)
        if (allocated(arh_object%metric_eigvecs)) deallocate(arh_object%metric_eigvecs)
        allocate(arh_object%metric_eigvals(n_dm, n_particle), &
                 arh_object%metric_eigvecs(n_dm, n_dm, n_particle), metric(n_dm, n_dm))

        do k = 1, n_particle
            ! generate ARH metric
            do j = 1, n_dm
                delta_j = dm_list(:, :, k, j) - dm_oao(:, :, k)
                do i = 1, j
                    delta_i = dm_list(:, :, k, i) - dm_oao(:, :, k)
                    ! compute Tr(delta_i * delta_j)
                    metric(i, j) = sum(delta_j * transpose(delta_i))
                    metric(j, i) = metric(i, j)
                end do
            end do
            ! diagonalize metric
            if (n_dm > 0) call symm_mat_diag(metric, eigvals(:, k), eigvecs(:, :, k), &
                                             settings, error)
            if (error /= 0) exit
        end do

        if (allocated(delta_i)) deallocate(delta_i)
        if (allocated(delta_j)) deallocate(delta_j)
        deallocate(metric)

    end subroutine get_arh_metric

    subroutine append(list, new_array)
        !
        ! this subroutine appends an array to a list of arrays of equal dimension
        !
        real(rp), intent(inout), allocatable :: list(:, :, :, :)
        real(rp), intent(in) :: new_array(:, :, :)

        integer(ip) :: n1, n2, n3
        real(rp), allocatable :: temp(:, :, :, :)

        n1 = size(new_array, 1)
        n2 = size(new_array, 2)
        n3 = size(new_array, 3)

        allocate(temp(n1, n2, n3, size(list, 4) + 1))
        temp(:, :, :, :size(list, 4)) = list
        temp(:, :, :, size(list, 4) + 1) = new_array
        deallocate(list)
        list = temp

    end subroutine append

end module otr_arh
