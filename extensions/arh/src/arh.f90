! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh

    use opentrustregion, only: rp, ip, settings_type, obj_func_type, update_orbs_type, &
                               hess_x_type, project_type
    use otr_oao, only: oao_settings_type, default_oao_settings, oao_type

    implicit none

    type, extends(oao_settings_type) :: arh_settings_type
        logical :: symm_arh
    contains
        procedure :: init => init_arh_settings
    end type arh_settings_type

    type(arh_settings_type), parameter :: default_arh_settings = &
        arh_settings_type(oao_settings_type = default_oao_settings, symm_arh = .true.)

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

    type, extends(oao_type) :: arh_type
        real(rp), allocatable :: fock_oao(:, :, :), same_v_oao(:, :, :), &
                                 opposite_v_oao(:, :, :), metric_eigvals(:, :), &
                                 metric_eigvecs(:, :, :), dm_list(:, :, :, :), &
                                 fock_list(:, :, :, :), same_v_list(:, :, :, :), &
                                 opposite_v_list(:, :, :, :), dm_diff(:, :, :, :), &
                                 fock_diff(:, :, :, :), same_v_diff(:, :, :, :), &
                                 opposite_v_diff(:, :, :, :)
        procedure(update_dm_jk_type), pointer, nopass :: update_dm_jk => null()
    end type arh_type

    ! global variables
    type(arh_type), pointer :: arh_object

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
        use otr_oao, only: get_energy_2d_type, update_dm_2d_type, oao_object, &
                           oao_factory_common, obj_func_oao, project_oao

        real(rp), intent(in) :: dm_ao(:, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_2d_type), intent(in), pointer :: get_energy
        procedure(update_dm_2d_type), intent(in), pointer :: update_dm
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        type(oao_type) :: temp_object
        logical :: copy_object

        ! initialize error flag
        error = 0

        ! allocate objects
        copy_object = .false.
        if (allocated(oao_object)) then
            select type (object => oao_object)
                type is (oao_type)
                    temp_object = oao_object
                    deallocate(oao_object)
                    copy_object = .true.
            end select
        end if
        if (.not. allocated(oao_object)) then
            allocate(arh_type :: oao_object)
            allocate(arh_settings_type :: oao_object%settings)
        end if
        select type (object => oao_object)
            type is (arh_type)
                if (copy_object) object%oao_type = temp_object
                arh_object => object
        end select

        ! call common setup
        call oao_factory_common(reshape(dm_ao, [n_ao, n_ao, 1]), ao_overlap, &
                                n_particle, n_ao, error, settings)

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
        use otr_oao, only: get_energy_3d_type, oao_factory_common, oao_object, &
                           obj_func_oao, project_oao

        real(rp), intent(in) :: dm_ao(:, :, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_3d_type), intent(in), pointer :: get_energy
        procedure(update_dm_jk_type), intent(in), pointer :: update_dm_jk
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(project_type), intent(out), pointer :: project_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        type(oao_type) :: temp_object
        logical :: copy_object

        ! initialize error flag
        error = 0

        ! allocate objects
        copy_object = .false.
        if (allocated(oao_object)) then
            select type (object => oao_object)
                type is (oao_type)
                    temp_object = oao_object
                    deallocate(oao_object)
                    copy_object = .true.
            end select
        end if
        if (.not. allocated(oao_object)) then
            allocate(arh_type :: oao_object)
            allocate(arh_settings_type :: oao_object%settings)
        end if
        select type (object => oao_object)
            type is (arh_type)
                if (copy_object) object%oao_type = temp_object
                arh_object => object
        end select

        ! call common setup
        call oao_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)

        ! set pointers to functions
        arh_object%get_energy => get_energy
        arh_object%update_dm_jk => update_dm_jk

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_oao
        update_orbs_arh_funptr => update_orbs_arh_open_shell
        project_arh_funptr => project_oao

    end subroutine arh_factory_open_shell

    subroutine update_orbs_arh_closed_shell(kappa, func, grad, h_diag, hess_x_funptr, &
                                            error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation 
        ! in the OAO basis and the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall for the closed-shell case
        !
        use opentrustregion, only: hess_x_type
        use otr_oao, only: rotate_dm_ao, symmetric_transformation, calculate_grad_h_diag

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list
        real(rp), allocatable :: fock_ao(:, :, :)
        type(arh_settings_type), pointer :: arh_settings

        external :: dgemm

        ! type narrowing
        select type (settings => arh_object%settings)
            type is (arh_settings_type)
                arh_settings => settings
        end select

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

        ! get energy, Fock matrix, and response function
        allocate(fock_ao(n_ao, n_ao, n_particle))
        call arh_object%update_dm_2d(arh_object%dm_ao(:, :, 1), func, &
                                     fock_ao(:, :, 1), arh_object%get_response_2d, &
                                     error)
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
                            arh_settings, error)
        if (error /= 0) return

        ! prepare differences for two-electron part of Hessian
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
        use otr_oao, only: rotate_dm_ao, symmetric_transformation, calculate_grad_h_diag

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list
        real(rp), allocatable :: fock_ao(:, :, :), fock_oao(:, :, :), &
                                 coulomb_ao(:, :, :), exchange_ao(:, :, :), &
                                 same_v_ao(:, :, :), opposite_v_ao(:, :, :)
        type(arh_settings_type), pointer :: arh_settings
                                 
        external :: dgemm

        ! type narrowing
        select type (settings => arh_object%settings)
            type is (arh_settings_type)
                arh_settings => settings
        end select

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

        ! get energy, Fock matrix, same and opposite spin potentials, and response 
        ! function
        allocate(fock_ao(n_ao, n_ao, n_particle), coulomb_ao(n_ao, n_ao, n_particle), &
                 exchange_ao(n_ao, n_ao, n_particle))
        call arh_object%update_dm_jk(arh_object%dm_ao, func, fock_ao, coulomb_ao, &
                                     exchange_ao, arh_object%get_response, error)
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
                            arh_settings, error)
        if (error /= 0) return

        ! prepare differences for two-electron part of Hessian
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
        type(arh_settings_type), pointer :: arh_settings

        external :: dgemm

        ! type narrowing
        select type (settings => arh_object%settings)
            type is (arh_settings_type)
                arh_settings => settings
        end select

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

        ! get one electron part
        allocate(hess_x_full(n_ao, n_ao, n_particle))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, arh_object%fock_vv(:, :, i) &
                       - arh_object%fock_oo(:, :, i), n_ao, x_full(:, :, i), n_ao, &
                       0.0_rp, hess_x_full(:, :, i), n_ao)
            hess_x_full(:, :, i) = hess_x_full(:, :, i) - &
                                   transpose(hess_x_full(:, :, i))
        end do

        ! get two electron contributions
        if (n_particle == 1) then
            hess_x_full = hess_x_full + &
                get_two_el_contribution_closed_shell(arh_object%dm_oao, x_full, &
                                                     arh_object%dm_diff, &
                                                     arh_object%fock_diff, &
                                                     arh_object%metric_eigvals, &
                                                     arh_object%metric_eigvecs, n_ao, &
                                                     n_particle, arh_settings)
        else
            hess_x_full = hess_x_full + &
                get_two_el_contribution_open_shell(arh_object%dm_oao, x_full, &
                                                   arh_object%dm_diff, &
                                                   arh_object%same_v_diff, &
                                                   arh_object%opposite_v_diff, &
                                                   arh_object%metric_eigvals, &
                                                   arh_object%metric_eigvecs, n_ao, &
                                                   n_particle, arh_settings)
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

    function get_two_el_contribution_closed_shell(dm_oao, x, dm_diff, fock_diff, &
                                                  metric_eigvals, metric_eigvecs, &
                                                  n_ao, n_particle, settings) &
                                                  result(two_el)
        !
        ! this function computes the two-electron contribution to the ARH Hessian for
        ! the closed-shell case
        !
        use opentrustregion, only: numerical_zero
        use otr_oao, only: project

        real(rp), intent(in) :: dm_oao(:, :, :), x(:, :, :), dm_diff(:, :, :, :), &
                                fock_diff(:, :, :, :), metric_eigvals(:, :), &
                                metric_eigvecs(:, :, :)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp) :: two_el(n_ao, n_ao, n_particle)
        type(arh_settings_type), intent(in) :: settings

        integer(ip) :: n_diff, i
        real(rp), allocatable :: dm_oao_x(:, :), vec1(:, :), vec2(:, :), tmp_vec1(:), &
                                 tmp_vec2(:)
        real(rp) :: factor

        ! number of density matrix differences
        n_diff = size(dm_diff, 4)

        ! get factor for symmetrization
        if (settings%symm_arh) then
            factor = 0.5_rp
        else
            factor = 1.0_rp
        end if

        ! calculate two-electron contributions
        two_el = 0.0_rp
        if (n_diff > 0) then
            ! get commutator of trial vector and current density matrix
            allocate(dm_oao_x(n_ao, n_ao))
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, 1), n_ao, &
                       x(:, :, 1), n_ao, 0.0_rp, dm_oao_x, n_ao)
            dm_oao_x = dm_oao_x + transpose(dm_oao_x)

            ! get traces of density and Fock matrix differences with commutator of 
            ! trial vector and current density matrix
            allocate(vec1(n_diff, 1))
            if (settings%symm_arh) allocate(vec2(n_diff, 1))
            do i = 1, n_diff
                vec1(i, 1) = sum(dm_diff(:, :, 1, i) * dm_oao_x)
                if (settings%symm_arh) vec2(i, 1) = sum(fock_diff(:, :, 1, i) * &
                                                        dm_oao_x)
            end do
            deallocate(dm_oao_x)

            ! multiply pseudoinverse metric with vectors
            allocate(tmp_vec1(n_diff))
            tmp_vec1 = 0.0_rp
            if (settings%symm_arh) then
                allocate(tmp_vec2(n_diff))
                tmp_vec2 = 0.0_rp
            end if
            do i = 1, n_diff
                if (metric_eigvals(i, 1) > numerical_zero) then
                    tmp_vec1 = tmp_vec1 + (dot_product(metric_eigvecs(:, i, 1), &
                                                       vec1(:, 1)) / &
                                           metric_eigvals(i, 1)) * &
                               metric_eigvecs(:, i, 1)
                    if (settings%symm_arh) then
                        tmp_vec2 = tmp_vec2 + (dot_product(metric_eigvecs(:, i, 1), &
                                                           vec2(:, 1)) / &
                                               metric_eigvals(i, 1)) * &
                                   metric_eigvecs(:, i, 1)
                    end if
                end if
            end do
            vec1(:, 1) = tmp_vec1
            deallocate(tmp_vec1)
            if (settings%symm_arh) then
                vec2(:, 1) = tmp_vec2
                deallocate(tmp_vec2)
            end if

            ! contract with Fock and density matrix differences
            do i = 1, n_diff
                two_el(:, :, 1) = two_el(:, :, 1) + factor * vec1(i, 1) * &
                                  fock_diff(:, :, 1, i)
                if (settings%symm_arh) then
                    two_el(:, :, 1) = two_el(:, :, 1) + factor * vec2(i, 1) * &
                                      dm_diff(:, :, 1, i)
                end if
            end do
            deallocate(vec1)
            if (settings%symm_arh) deallocate(vec2)

            ! add only v-o and o-v contributions of ARH two-electron part
            two_el = project(two_el, dm_oao)
        end if
        
    end function get_two_el_contribution_closed_shell

    function get_two_el_contribution_open_shell(dm_oao, x, dm_diff, same_v_diff, &
                                                opposite_v_diff, metric_eigvals, &
                                                metric_eigvecs, n_ao,  n_particle, &
                                                settings) result(two_el)
        !
        ! this function computes the two-electron contribution to the ARH Hessian for
        ! the open-shell case
        !
        use opentrustregion, only: numerical_zero
        use otr_oao, only: project

        real(rp), intent(in) :: dm_oao(:, :, :), x(:, :, :), dm_diff(:, :, :, :), &
                                same_v_diff(:, :, :, :), opposite_v_diff(:, :, :, :), &
                                metric_eigvals(:, :), metric_eigvecs(:, :, :)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp) :: two_el(n_ao, n_ao, n_particle)
        type(arh_settings_type), intent(in) :: settings

        integer(ip) :: n_diff, i, j, k
        real(rp), allocatable :: dm_oao_x(:, :, :), vec1(:, :), vec2(:), vec3(:, :), &
                                 tmp_vec1(:, :), tmp_vec2(:), tmp_vec3(:, :)
        real(rp) :: factor

        real(rp), external :: ddot

        ! number of density matrix differences
        n_diff = size(dm_diff, 4)

        ! get factor for symmetrization
        if (settings%symm_arh) then
            factor = 0.5_rp
        else
            factor = 1.0_rp
        end if

        ! calculate two-electron contributions
        two_el = 0.0_rp
        if (n_diff > 0) then
            ! get commutator of trial vector and current density matrix
            allocate(dm_oao_x(n_ao, n_ao, n_particle))
            do j = 1, n_particle
                call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, j), n_ao, &
                           x(:, :, j), n_ao, 0.0_rp, dm_oao_x(:, :, j), n_ao)
                dm_oao_x(:, :, j) = dm_oao_x(:, :, j) + transpose(dm_oao_x(:, :, j))
            end do
            
            ! loop over particles
            allocate(vec1(n_diff, n_particle))
            if (settings%symm_arh) allocate(vec2(n_diff), vec3(n_diff, n_particle))
            do j = 1, n_particle
                ! get traces of density matrix and integral differences with commutator 
                ! of trial vector and current density matrix
                do i = 1, n_diff
                    if (settings%symm_arh) vec2(i) = &
                        sum(same_v_diff(:, :, j, i) * dm_oao_x(:, :, j))
                    do k = 1, n_particle
                        vec1(i, k) = sum(dm_diff(:, :, k, i) * dm_oao_x(:, :, k))
                        if (settings%symm_arh) vec3(i, k) = &
                            sum(opposite_v_diff(:, :, j, i) * dm_oao_x(:, :, k))
                    end do
                end do

                ! multiply pseudoinverse metric with vectors
                allocate(tmp_vec1(n_diff, n_particle))
                tmp_vec1 = 0.0_rp
                if (settings%symm_arh) then
                    allocate(tmp_vec2(n_diff), tmp_vec3(n_diff, n_particle))
                    tmp_vec2 = 0.0_rp
                    tmp_vec3 = 0.0_rp
                end if
                do i = 1, n_diff
                    if (metric_eigvals(i, j) > numerical_zero) then
                        if (settings%symm_arh) then
                            tmp_vec2 = &
                                tmp_vec2 + &
                                ddot(n_diff, metric_eigvecs(:, i, j), 1, vec2, 1) &
                                / metric_eigvals(i, j) * metric_eigvecs(:, i, j)
                        end if
                    end if
                    do k = 1, n_particle
                        if (metric_eigvals(i, k) > numerical_zero) then
                            tmp_vec1(:, k) = &
                                tmp_vec1(:, k) + &
                                ddot(n_diff, metric_eigvecs(:, i, k), 1, vec1(:, k), &
                                     1) / metric_eigvals(i, k) * metric_eigvecs(:, i, k)
                        end if
                        if (settings%symm_arh) then
                            if (metric_eigvals(i, j) > numerical_zero) then
                                tmp_vec3(:, k) = &
                                    tmp_vec3(:, k) + &
                                    ddot(n_diff, metric_eigvecs(:, i, j), 1, &
                                         vec3(:, k), 1) / metric_eigvals(i, j) * &
                                    metric_eigvecs(:, i, j)
                            end if
                        end if
                    end do
                end do
                vec1 = tmp_vec1
                deallocate(tmp_vec1)
                if (settings%symm_arh) then
                    vec2 = tmp_vec2
                    vec3 = tmp_vec3
                    deallocate(tmp_vec2, tmp_vec3)
                end if

                ! contract with same and opposite spin potential and density matrix 
                ! differences
                do i = 1, n_diff
                    two_el(:, :, j) = two_el(:, :, j) + factor * vec1(i, j) * &
                                      same_v_diff(:, :, j, i)
                    if (settings%symm_arh) then
                        two_el(:, :, j) = two_el(:, :, j) + factor * vec2(i) * &
                                          dm_diff(:, :, j, i)
                    end if
                    do k = 1, n_particle
                        if (k == j) cycle
                        two_el(:, :, j) = two_el(:, :, j) + factor * vec1(i, k) * &
                                          opposite_v_diff(:, :, k, i)
                        if (settings%symm_arh) then
                            two_el(:, :, j) = two_el(:, :, j) + factor * vec3(i, k) * &
                                              dm_diff(:, :, j, i)
                        end if
                    end do
                end do
            end do
            deallocate(dm_oao_x, vec1)
            if (settings%symm_arh) deallocate(vec2, vec3)

            ! add only v-o and o-v contributions of ARH two-electron part
            two_el = project(two_el, dm_oao)
        end if

    end function get_two_el_contribution_open_shell

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
