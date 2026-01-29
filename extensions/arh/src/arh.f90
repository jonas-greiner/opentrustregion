! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh

    ! Things to do:
    ! Add project function to OTR
    ! Replace preconditioner with project function
    ! Enable guess vectors (h_diag and random) with project function

    use opentrustregion, only: rp, ip, settings_type, obj_func_type, update_orbs_type, &
                               hess_x_type, precond_type

    implicit none

    type, extends(settings_type) :: arh_settings_type
        logical :: restricted
    contains
        procedure :: init => init_arh_settings
    end type arh_settings_type

    type(arh_settings_type), parameter :: default_arh_settings = &
        arh_settings_type(logger = null(), initialized = .true., restricted = .false., &
                          verbose = 0)

    abstract interface
        function get_energy_closed_shell_type(dm, error) result(energy)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :)
            integer(ip), intent(out) :: error

            real(rp) :: energy
        end function get_energy_closed_shell_type
    end interface

    abstract interface
        function get_energy_open_shell_type(dm, error) result(energy)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            integer(ip), intent(out) :: error

            real(rp) :: energy
        end function get_energy_open_shell_type
    end interface

    abstract interface
        subroutine get_fock_type(dm, energy, fock, error)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :)
            integer(ip), intent(out) :: error
            real(rp), intent(out) :: energy
            real(rp), intent(out), target :: fock(:, :)
        end subroutine get_fock_type
    end interface

    abstract interface
        subroutine get_fock_jk_type(dm, energy, fock, coulomb, exchange, error)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            integer(ip), intent(out) :: error
            real(rp), intent(out) :: energy
            real(rp), intent(out), target :: fock(:, :, :), coulomb(:, :, :), &
                                             exchange(:, :, :)
        end subroutine get_fock_jk_type
    end interface

    type :: arh_type
        type(arh_settings_type) :: settings
        integer(ip) :: n_ao, n_param, n_particle
        real(rp), allocatable :: dm_ao(:, :, :), s_sqrt(:, :), s_inv_sqrt(:, :), &
                                 dm_oao(:, :, :), fock_oao(:, :, :), &
                                 same_v_oao(:, :, :), opposite_v_oao(:, :, :), &
                                 fock_oo(:, :, :), fock_vv(:, :, :), metric(:, :, :), &
                                 dm_list(:, :, :, :), fock_list(:, :, :, :), &
                                 same_v_list(:, :, :, :), opposite_v_list(:, :, :, :), &
                                 dm_diff(:, :, :, :), fock_diff(:, :, :, :), &
                                 same_v_diff(:, :, :, :), opposite_v_diff(:, :, :, :), &
                                 h_diag_test(:)
        procedure(get_energy_closed_shell_type), pointer, nopass :: &
            get_energy_closed_shell
        procedure(get_energy_open_shell_type), pointer, nopass :: get_energy_open_shell
        procedure(get_fock_type), pointer, nopass :: get_fock
        procedure(get_fock_jk_type), pointer, nopass :: get_fock_jk
    end type arh_type

    ! global variables
    type(arh_type) :: arh_object

    ! create function pointers to ensure that routines comply with interface
    procedure(obj_func_type), pointer :: obj_func_arh_ptr => obj_func_arh
    procedure(update_orbs_type), pointer :: update_orbs_arh_closed_shell_ptr => &
        update_orbs_arh_closed_shell
    procedure(update_orbs_type), pointer :: update_orbs_arh_open_shell_ptr => &
        update_orbs_arh_open_shell
    procedure(hess_x_type), pointer :: hess_x_arh_ptr => hess_x_arh
    procedure(precond_type), pointer :: precond_arh_ptr => &
        level_shifted_diag_precond_arh

    ! define module procedures for different spin cases
    interface arh_factory
        module procedure arh_factory_closed_shell, arh_factory_open_shell
    end interface arh_factory

    interface arh_deconstructor
        module procedure arh_deconstructor_closed_shell, arh_deconstructor_open_shell
    end interface arh_deconstructor

    contains

    subroutine arh_factory_closed_shell(dm_ao, ao_overlap, n_particle, n_ao, &
                                        get_energy, get_fock, obj_func_arh_funptr, &
                                        update_orbs_arh_funptr, precond_arh_funptr, &
                                        error, settings)
        !
        ! this function returns a modified ARH orbital updating function for the 
        ! closed-shell case
        !
        real(rp), intent(in) :: dm_ao(:, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_closed_shell_type), intent(in), pointer :: get_energy
        procedure(get_fock_type), intent(in), pointer :: get_fock
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(precond_type), intent(out), pointer :: precond_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

        ! call common setup
        call arh_factory_common(reshape(dm_ao, [n_ao, n_ao, 1]), ao_overlap, &
                                n_particle, n_ao, error, settings)

        ! set pointers to functions
        arh_object%get_energy_closed_shell => get_energy
        arh_object%get_fock => get_fock

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh
        update_orbs_arh_funptr => update_orbs_arh_closed_shell
        precond_arh_funptr => level_shifted_diag_precond_arh

    end subroutine arh_factory_closed_shell

    subroutine arh_factory_open_shell(dm_ao, ao_overlap, n_particle, n_ao, get_energy, &
                                      get_fock_jk, obj_func_arh_funptr, &
                                      update_orbs_arh_funptr, precond_arh_funptr, &
                                      error, settings)
        !
        ! this function returns a modified ARH orbital updating function for the 
        ! open-shell case
        !
        real(rp), intent(in) :: dm_ao(:, :, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_open_shell_type), intent(in), pointer :: get_energy
        procedure(get_fock_jk_type), intent(in), pointer :: get_fock_jk
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(precond_type), intent(out), pointer :: precond_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

        ! call common setup
        call arh_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)

        ! set pointers to functions
        arh_object%get_energy_open_shell => get_energy
        arh_object%get_fock_jk => get_fock_jk

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh
        update_orbs_arh_funptr => update_orbs_arh_open_shell
        precond_arh_funptr => level_shifted_diag_precond_arh

    end subroutine arh_factory_open_shell

    subroutine arh_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        !
        ! this function returns a modified ARH orbital updating function
        !
        real(rp), intent(in) :: dm_ao(:, :, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

        ! set settings
        arh_object%settings = settings

        ! number of particles
        arh_object%n_particle = n_particle

        ! get number of atomic orbitals
        arh_object%n_ao = n_ao

        ! get number of non-redundant parameters
        arh_object%n_param = n_ao * (n_ao - 1) / 2
        if (.not. arh_object%settings%restricted) &
            arh_object%n_param = n_particle * arh_object%n_param

        ! starting density matrix
        arh_object%dm_ao = dm_ao

        ! get square root and inverse square root of AO overlap matrix
        call compute_sqrt_and_inv_sqrt(ao_overlap, arh_object%s_sqrt, &
                                       arh_object%s_inv_sqrt, error)
        if (error /= 0) return

        ! get per spin contribution to density matrix in orthogonalized AO basis
        arh_object%dm_oao = symmetric_transformation(arh_object%s_sqrt, dm_ao)

        ! allocate matrices
        allocate(arh_object%fock_oo(n_ao, n_ao, arh_object%n_particle), &
                 arh_object%fock_vv(n_ao, n_ao, arh_object%n_particle), &
                 arh_object%h_diag_test(arh_object%n_param))

    end subroutine arh_factory_common

    function obj_func_arh(kappa, error) result(energy)
        !
        ! this function defines the energy evaluation in OAO basis
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        real(rp), allocatable :: rot_dm_ao(:, :, :)

        ! initialize energy in case of error
        energy = 0.0_rp

        ! get rotated density matrix in AO basis
        allocate(rot_dm_ao(arh_object%n_ao, arh_object%n_ao, arh_object%n_particle))
        call rotate_dm_ao(kappa, arh_object%n_particle, arh_object%n_ao, &
                          arh_object%settings%restricted, rot_dm_ao, error)
        if (error /= 0) return

        ! calculate mean-field energy
        if (arh_object%n_particle == 1) then
            energy = arh_object%get_energy_closed_shell(rot_dm_ao(:, :, 1), error)
        else
            energy = arh_object%get_energy_open_shell(rot_dm_ao, error)
        end if
        if (error /= 0) return

    end function obj_func_arh

    subroutine update_orbs_arh_closed_shell(kappa, func, grad, h_diag, hess_x_funptr, &
                                          error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation 
        ! in the OAO basis and the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall for the closed-shell case
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list
        real(rp), allocatable :: fock_ao(:, :)
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
        call rotate_dm_ao(kappa, n_particle, n_ao, arh_object%settings%restricted, &
                          arh_object%dm_ao, error, arh_object%dm_oao)
        if (error /= 0) return

        ! get mean-field energy and Fock matrix
        allocate(fock_ao(n_ao, n_ao))
        call arh_object%get_fock(reshape(arh_object%dm_ao, [n_ao, n_ao]), func, &
                                 fock_ao, error)
        if (error /= 0) then
            deallocate(fock_ao)
            return
        end if

        ! transform Fock matrix to OAO basis
        arh_object%fock_oao = &
            symmetric_transformation(arh_object%s_inv_sqrt, &
                                     reshape(fock_ao, [n_ao, n_ao, n_particle]))
        deallocate(fock_ao)

        ! calculate gradient and Hessian diagonal
        call calculate_grad_h_diag(arh_object%dm_oao, arh_object%fock_oao, n_particle, &
                                   n_ao, arh_object%settings%restricted, grad, h_diag, &
                                   arh_object%fock_oo, arh_object%fock_vv)
        arh_object%h_diag_test = h_diag

        ! construct ARH metric
        arh_object%metric = get_arh_metric(arh_object%dm_list, arh_object%dm_oao)

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

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list
        real(rp), allocatable :: fock_ao(:, :, :), fock_oao(:, :, :), &
                                 coulomb_ao(:, :, :), exchange_ao(:, :, :), &
                                 same_v_ao(:, :, :), opposite_v_ao(:, :, :)
                                 
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

        ! get mean-field energy, Fock matrix and same and opposite spin potentials
        allocate(fock_ao(n_ao, n_ao, n_particle), coulomb_ao(n_ao, n_ao, n_particle), &
                 exchange_ao(n_ao, n_ao, n_particle))
        call arh_object%get_fock_jk(arh_object%dm_ao, func, fock_ao, coulomb_ao, &
                                    exchange_ao, error)
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
        opposite_v_ao(:, :, 1) = coulomb_ao(:, :, 2)
        opposite_v_ao(:, :, 2) = coulomb_ao(:, :, 1)

        ! transform same and opposite spin potentials to OAO basis
        arh_object%same_v_oao = symmetric_transformation(arh_object%s_inv_sqrt, &
                                                         same_v_ao)
        arh_object%opposite_v_oao = symmetric_transformation(arh_object%s_inv_sqrt, &
                                                             opposite_v_ao)
        deallocate(coulomb_ao, exchange_ao)

        ! calculate gradient and Hessian diagonal
        call calculate_grad_h_diag(arh_object%dm_oao, fock_oao, n_particle, n_ao, &
                                   arh_object%settings%restricted, grad, h_diag, &
                                   arh_object%fock_oo, arh_object%fock_vv)
        arh_object%h_diag_test = h_diag
        deallocate(fock_oao)

        ! construct ARH metric
        arh_object%metric = get_arh_metric(arh_object%dm_list, arh_object%dm_oao)

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
        ! internal and this redundant in the UHF context. Without this projection, 
        ! the alpha Hessian vector product would leak non-zero values into 
        ! the [c-c] block and the beta Hessian vector product into the [v-v] blocks. 
        ! This projection thus prevents singular dimensions in the ARH solver.
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
                                                     arh_object%metric, n_ao, &
                                                     n_particle, arh_object%settings, &
                                                     error)
        else
            hess_x_full = hess_x_full + &
                get_two_el_contribution_open_shell(arh_object%dm_oao, x_full, &
                                                   arh_object%dm_diff, &
                                                   arh_object%same_v_diff, &
                                                   arh_object%opposite_v_diff, &
                                                   arh_object%metric, n_ao, &
                                                   n_particle, arh_object%settings, &
                                                   error)
        end if
        if (error /= 0) return
        deallocate(x_full)

        ! pack Hessian linear transformation
        call pack_asymm(hess_x_full, hess_x, arh_object%settings%restricted)

    end subroutine hess_x_arh

    subroutine level_shifted_diag_precond_arh(vector, mu, precond_vector, error)
        !
        ! this function defines the default level-shifted diagonal preconditioner
        ! but only extracting v-o and o-v contributions
        !
        real(rp), intent(in), target :: vector(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_vector(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: precond_arr(:),  precond_vector_full(:, :, :)

        ! initialize error flag
        error = 0
        
        ! construct level-shifted preconditioner
        precond_arr = arh_object%h_diag_test - mu
        where (abs(precond_arr) < 1d-10)
            precond_arr = 1d-10
        end where

        ! precondition vector
        precond_vector = vector / precond_arr
        deallocate(precond_arr)

        ! unpack vector
        precond_vector_full = unpack_asymm(precond_vector, arh_object%n_particle, &
                                           arh_object%n_ao, &
                                           arh_object%settings%restricted)

        ! extract only v-o and o-v contributions
        precond_vector_full = project(precond_vector_full, arh_object%dm_oao)

        ! pack vector
        call pack_asymm(precond_vector_full, precond_vector, &
                        arh_object%settings%restricted)
        
    end subroutine level_shifted_diag_precond_arh

    function project(matrix, dm_oao) result(projected_matrix)
        !
        ! this function only retains occupied-virtual and virtual-occupied 
        ! contributions to a matrix
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

            ! construct virtual-occupied contributions DM(I-D)
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, matrix(:, :, i), n_ao, &
                       proj_v, n_ao, 0.0_rp, temp, n_ao)
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, i), n_ao, &
                       temp, n_ao, 0.0_rp, projected_matrix(:, :, i), n_ao)

            ! add occupied-virtual contributions (I-D)MD
            projected_matrix(:, :, i) = projected_matrix(:, :, i) - &
                                        transpose(projected_matrix(:, :, i))
        end do
        deallocate(proj_v, temp)

    end function project

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

    subroutine arh_deconstructor_closed_shell(dm_ao, error)
        !
        ! this subroutine deallocates the ARH objects for the closed-shell case
        !
        use opentrustregion, only: verbosity_error

        real(rp), intent(out) :: dm_ao(:, :)
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! get final density matrix in AO basis
        if (.not. allocated(arh_object%dm_ao)) then
            call arh_object%settings%log("AO density matrix is not allocated. The "// &
                                         "ARH deconstructor should only be run "// &
                                         "after both the ARH factory and the OTR "// &
                                         "solver have been called.", verbosity_error, &
                                         .true.)
            error = 1
            return
        end if
        dm_ao = arh_object%dm_ao(:, :, 1)

        ! deallocate all allocated arrays
        deallocate(arh_object%dm_ao, arh_object%s_sqrt, arh_object%s_inv_sqrt, &
                   arh_object%dm_oao, arh_object%fock_oao, arh_object%fock_oo, &
                   arh_object%fock_vv, arh_object%dm_list, arh_object%fock_list, &
                   arh_object%dm_diff, arh_object%fock_diff)

    end subroutine arh_deconstructor_closed_shell

    subroutine arh_deconstructor_open_shell(dm_ao, error)
        !
        ! this subroutine deallocates the ARH objects for the open-shell case
        !
        use opentrustregion, only: verbosity_error

        real(rp), intent(out) :: dm_ao(:, :, :)
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! get final density matrix in AO basis
        if (.not. allocated(arh_object%dm_ao)) then
            call arh_object%settings%log("AO density matrix is not allocated. The "// &
                                         "ARH deconstructor should only be run "// &
                                         "after both the ARH factory and the OTR "// &
                                         "solver have been called.", verbosity_error, &
                                         .true.)
            error = 1
            return
        end if
        dm_ao = arh_object%dm_ao

        ! deallocate all allocated arrays
        deallocate(arh_object%dm_ao, arh_object%s_sqrt, arh_object%s_inv_sqrt, &
                   arh_object%dm_oao, arh_object%fock_oo, arh_object%fock_vv, &
                   arh_object%same_v_oao, arh_object%opposite_v_oao, &
                   arh_object%dm_list, arh_object%same_v_list, &
                   arh_object%opposite_v_list, arh_object%dm_diff, &
                   arh_object%same_v_diff, arh_object%opposite_v_diff)

    end subroutine arh_deconstructor_open_shell

    subroutine rotate_dm_ao(kappa, n_particle, n_ao, restricted, rot_dm_ao, error, &
                            rot_dm_oao)
        !
        ! this subroutine returns the rotated density matrix in AO basis
        !
        real(rp), intent(in) :: kappa(:)
        integer(ip), intent(in) :: n_particle, n_ao
        logical, intent(in) :: restricted
        integer(ip), intent(out) :: error
        real(rp), intent(out) :: rot_dm_ao(:, :, :)
        real(rp), intent(out), target, optional :: rot_dm_oao(:, :, :)

        real(rp), allocatable :: kappa_full(:, :, :), temp(:, :)
        real(rp), pointer :: rot_dm_oao_ptr(:, :, :), u_ptr(:, :)
        real(rp), allocatable, target :: u(:, :, :), rot_dm_oao_local(:, :, :)
        integer(ip) :: n_rot_mat, i

        ! initialize error flag
        error = 0

        ! get rotation matrix
        n_rot_mat = n_particle
        if (restricted) n_rot_mat = 1
        allocate(u(n_ao, n_ao, n_rot_mat))
        kappa_full = unpack_asymm(kappa, n_rot_mat, n_ao, restricted)
        do i = 1, n_rot_mat
            u(:, :, i) = matrix_exponential(kappa_full(:, :, i), error)
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
        ! spin-restricted case: use same rotation matrix for both spins
        if (restricted) u_ptr => u(:, :, n_rot_mat)
        do i = 1, n_particle
            ! spin-unrestricted case: use corresponding rotation matrix
            if (.not. restricted) u_ptr => u(:, :, i)
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, arh_object%dm_oao(:, :, i), &
                       n_ao, u_ptr, n_ao, 0.0_rp, temp, n_ao)
            call dgemm("T", "N", n_ao, n_ao, n_ao, 1.0_rp, u_ptr, n_ao, temp, n_ao, &
                       0.0_rp, rot_dm_oao_ptr(:, :, i), n_ao)
        end do

        ! purify density matrix
        call purify(rot_dm_oao_ptr)

        ! transform density matrix from OAO basis to AO basis
        rot_dm_ao = symmetric_transformation(arh_object%s_inv_sqrt, rot_dm_oao_ptr)

        ! deallocate local memory if needed
        if (.not. present(rot_dm_oao)) deallocate(rot_dm_oao_local)

    end subroutine rotate_dm_ao

    subroutine calculate_grad_h_diag(dm_oao, fock_oao, n_particle, n_ao, restricted, &
                                     grad, h_diag, fock_oo, fock_vv)
        !
        ! this function calculates the gradient and Hessian diagonal in OAO basis while 
        ! also returning the occupied-occupied and virtual-virtual parts of the Fock 
        ! matrix
        !
        real(rp), intent(in) :: dm_oao(:, :, :), fock_oao(:, :, :)
        integer(ip), intent(in) :: n_particle, n_ao
        logical, intent(in) :: restricted
        real(rp), intent(out) :: grad(:), h_diag(:)
        real(rp), intent(out) :: fock_oo(:, :, :), fock_vv(:, :, :)

        integer(ip) :: i, j, k, idx
        real(rp), allocatable :: dm_fock_oao(:, :), fock_dm_oao(:, :), &
                                 fock_ov(:, :, :), fock_vo(:, :, :), grad_full(:, :, :)
        external :: dgemm

        ! get contributions to Fock matrix based on occupancies
        allocate(dm_fock_oao(n_ao, n_ao), fock_dm_oao(n_ao, n_ao), &
                 fock_ov(n_ao, n_ao, n_particle), fock_vo(n_ao, n_ao, n_particle))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, i), n_ao, &
                       fock_oao(:, :, i), n_ao, 0.0_rp, dm_fock_oao, n_ao)   
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_fock_oao, n_ao, &
                       dm_oao(:, :, i), n_ao, 0.0_rp, fock_oo(:, :, i), &
                       n_ao) ! DFD
            fock_ov(:, :, i) = dm_fock_oao - fock_oo(:, :, i) ! DF(I-D)
            fock_vo(:, :, i) = transpose(fock_ov(:, :, i)) ! (I_D)FD
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, fock_oao(:, :, i), n_ao, &
                       dm_oao(:, :, i), n_ao, 0.0_rp, fock_dm_oao, n_ao)
            fock_vv(:, :, i) = fock_oao(:, :, i) - dm_fock_oao - fock_dm_oao + &
                               fock_oo(:, :, i) ! (I_D)F(I_D)
        end do
        deallocate(dm_fock_oao, fock_dm_oao)

        ! construct gradient
        grad_full = fock_ov - fock_vo
        deallocate(fock_ov, fock_vo)

        ! pack gradient
        call pack_asymm(grad_full, grad, restricted)

        ! construct Hessian diagonal
        idx = 1
        if (restricted) then
            do j = 1, n_ao
                do i = 1, j - 1
                    h_diag(idx) = sum(fock_vv(i, i, :) + fock_vv(j, j, :) - &
                                      fock_oo(i, i, :) - fock_oo(j, j, :))
                    idx = idx + 1
                end do
            end do
        else
            do k = 1, n_particle
                do j = 1, n_ao
                    do i = 1, j - 1
                        h_diag(idx) = fock_vv(i, i, k) + fock_vv(j, j, k) - &
                                      fock_oo(i, i, k) - fock_oo(j, j, k)
                        idx = idx + 1
                    end do
                end do
            end do
        end if

    end subroutine calculate_grad_h_diag

    function get_two_el_contribution_closed_shell(dm_oao, x, dm_diff, fock_diff, &
                                                  metric, n_ao, n_particle, settings, &
                                                  error) result(two_el)
        !
        ! this function computes the two-electron contribution to the ARH Hessian for
        ! the closed-shell case
        !
        real(rp), intent(in) :: dm_oao(:, :, :), x(:, :, :), dm_diff(:, :, :, :), &
                                fock_diff(:, :, :, :), metric(:, :, :)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp) :: two_el(n_ao, n_ao, n_particle)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_diff, i
        real(rp), allocatable :: dm_oao_x(:, :), vec(:, :)

        ! number of density matrix differences
        n_diff = size(dm_diff, 4)

        ! calculate two-electron contributions
        two_el = 0.0_rp
        if (n_diff > 0) then
            ! contract density matrix with trial vector
            allocate(dm_oao_x(n_ao, n_ao))
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, 1), n_ao, &
                       x(:, :, 1), n_ao, 0.0_rp, dm_oao_x, n_ao)
            dm_oao_x = dm_oao_x + transpose(dm_oao_x)

            ! get vector containing density matrix difference
            allocate(vec(n_diff, 1))
            do i = 1, n_diff
                vec(i, 1) = sum(dm_diff(:, :, 1, i) * dm_oao_x)
            end do
            deallocate(dm_oao_x)

            ! multiply inverse metric with vector through linear solver
            call solve_symm_linear_system(metric(:, :, 1), vec, settings, error)
            if (error /= 0) return

            ! contract with Fock matrix differences
            do i = 1, n_diff
                two_el(:, :, 1) = two_el(:, :, 1) + vec(i, 1) * fock_diff(:, :, 1, i)
            end do
            deallocate(vec)

            ! add only v-o and o-v contributions of ARH two-electron part
            two_el = project(two_el, dm_oao)
        end if
        
    end function get_two_el_contribution_closed_shell

    function get_two_el_contribution_open_shell(dm_oao, x, dm_diff, same_v_diff, &
                                                opposite_v_diff, metric, n_ao, &
                                                n_particle, settings, error) &
            result(two_el)
        !
        ! this function computes the two-electron contribution to the ARH Hessian for
        ! the open-shell case
        !
        real(rp), intent(in) :: dm_oao(:, :, :), x(:, :, :), dm_diff(:, :, :, :), &
                                same_v_diff(:, :, :, :), opposite_v_diff(:, :, :, :), &
                                metric(:, :, :)
        integer(ip), intent(in) :: n_ao, n_particle
        real(rp) :: two_el(n_ao, n_ao, n_particle)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_diff, i, j, k
        real(rp), allocatable :: dm_oao_x(:, :, :), vec(:, :)

        ! number of density matrix differences
        n_diff = size(dm_diff, 4)

        ! calculate two-electron contributions
        two_el = 0.0_rp
        if (n_diff > 0) then
            ! contract density matrix with trial vector
            allocate(dm_oao_x(n_ao, n_ao, n_particle))
            do j = 1, n_particle
                call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_oao(:, :, j), n_ao, &
                           x(:, :, j), n_ao, 0.0_rp, dm_oao_x(:, :, j), n_ao)
                dm_oao_x(:, :, j) = dm_oao_x(:, :, j) + transpose(dm_oao_x(:, :, j))
            end do
            
            ! loop over particles
            allocate(vec(n_diff, n_particle))
            do j = 1, n_particle
                ! get vector containing density matrix difference
                do i = 1, n_diff
                    do k = 1, n_particle
                        vec(i, k) = sum(dm_diff(:, :, j, i) * dm_oao_x(:, :, k))
                    end do
                end do

                ! multiply inverse metric with vector through linear solver
                call solve_symm_linear_system(metric(:, :, j), vec, settings, error)
                if (error /= 0) return

                ! contract with same and opposite spin potential differences
                do i = 1, n_diff
                    two_el(:, :, j) = two_el(:, :, j) + vec(i, j) * &
                                      (same_v_diff(:, :, j, i) - &
                                       opposite_v_diff(:, :, j, i))
                    do k = 1, n_particle
                        two_el(:, :, j) = two_el(:, :, j) + vec(i, k) * &
                                          opposite_v_diff(:, :, k, i)
                    end do
                end do
            end do
            deallocate(dm_oao_x, vec)

            ! add only v-o and o-v contributions of ARH two-electron part
            two_el = project(two_el, dm_oao)
        end if

    end function get_two_el_contribution_open_shell

    function get_arh_metric(dm_list, dm_oao) result(metric)
        !
        ! this function calculates the augmented Roothaan-Hall metric
        !
        real(rp), intent(in) :: dm_list(:, :, :, :)
        real(rp), intent(in) :: dm_oao(:, :, :)
        real(rp), allocatable :: metric(:, :, :)

        integer(ip) :: n_particle, n_dm, i, j, k
        real(rp), allocatable :: delta_i(:, :), delta_j(:, :)

        n_particle = size(dm_list, 3)
        n_dm = size(dm_list, 4)

        if (allocated(metric)) deallocate(metric)
        allocate(metric(n_dm, n_dm, n_particle))

        ! generate ARH metric
        do k = 1, n_particle
            do j = 1, n_dm
                delta_j = dm_list(:, :, k, j) - dm_oao(:, :, k)
                do i = 1, j
                    delta_i = dm_list(:, :, k, i) - dm_oao(:, :, k)
                    ! compute Tr(delta_i * delta_j)
                    metric(i, j, k) = sum(delta_j * transpose(delta_i))
                    metric(j, i, k) = metric(i, j, k)
                end do
            end do
        end do

    end function get_arh_metric

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

    function unpack_asymm(matrix_nonred, n_particle, n_ao, spin_sum) result(matrix)
        !
        ! this function unpacks an antisymmetric matrix and returns the resulting 
        ! unpacked matrix
        !
        real(rp), intent(in) :: matrix_nonred(:)
        integer(ip), intent(in) :: n_particle, n_ao
        logical, intent(in) :: spin_sum
        
        real(rp), allocatable :: matrix(:, :, :)

        integer(ip) :: i, j, k, idx

        ! allocate and initialize full matrix
        allocate(matrix(n_ao, n_ao, n_particle))
        matrix = 0.0_rp

        ! initialize index
        idx = 1

        ! spin-restricted case: same matrix for both spins
        if (spin_sum) then
            do j = 1, n_ao
                do i = 1, j - 1
                    matrix(i, j, 1) = matrix_nonred(idx)
                    matrix(j, i, 1) = -matrix_nonred(idx)
                    idx = idx + 1
                end do
            end do
            do i = 2, n_particle
                matrix(:, :, i) = matrix(:, :, 1)
            end do
        ! spin-unrestricted case: unpack each spin separately
        else
            do k = 1, n_particle
                do j = 1, n_ao
                    do i = 1, j - 1
                        matrix(i, j, k) = matrix_nonred(idx)
                        matrix(j, i, k) = -matrix_nonred(idx)
                        idx = idx + 1
                    end do
                end do
            end do
        end if

    end function unpack_asymm

    subroutine pack_asymm(matrix, matrix_nonred, spin_sum)
        !
        ! this subroutine packs an antisymmetric matrix for RHF while deallocating the 
        ! original unpacked matrix
        !
        real(rp), intent(inout), allocatable :: matrix(:, :, :)
        real(rp), intent(out) :: matrix_nonred(:)
        logical, intent(in) :: spin_sum

        integer(ip) :: idx, i, j, k

        ! initialize index
        idx = 1

        ! spin-restricted case: sum over both spins
        if (spin_sum) then
            do i = 1, size(matrix, 2)
                do j = 1, i - 1
                    matrix_nonred(idx) = sum(matrix(j, i, :))
                    idx = idx + 1
                end do
            end do
        ! spin-unrestricted case: pack each spin separately
        else
            do i = 1, size(matrix, 3)
                do j = 1, size(matrix, 2)
                    do k = 1, j - 1
                        matrix_nonred(idx) = matrix(k, j, i)
                        idx = idx + 1
                    end do
                end do
            end do
        end if
        deallocate(matrix)

    end subroutine pack_asymm

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

    subroutine compute_sqrt_and_inv_sqrt(A, sqrtA, inv_sqrtA, error)
        ! 
        ! this subroutine calculates the square root and inverse square root of a 
        ! matrix
        !
        use opentrustregion, only: solver_settings_type, verbosity_error

        real(rp), intent(in)  :: A(:, :)
        real(rp), allocatable, intent(out) :: sqrtA(:, :), inv_sqrtA(:, :)
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, lwork, info, i
        real(rp), allocatable :: eigvecs(:, :), eigvals(:), work(:)
        character(300) :: msg
        external :: dsyev, dgemm

        ! initialize error flag
        error = 0

        ! get number of AOs
        n_ao = size(A, 1)

        ! allocate eigenvector and eigenvalue arrays
        allocate(eigvecs(n_ao, n_ao), eigvals(n_ao))

        ! copy input because dsyev overwrites it
        eigvecs = A

        ! query optimal workspace size
        lwork = -1
        allocate(work(1))
        call dsyev("V", "U", n_ao, eigvecs, n_ao, eigvals, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call dsyev("V", "U", n_ao, eigvecs, n_ao, eigvals, work, lwork, info)

        ! deallocatework array
        deallocate(work)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call arh_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! get square roots of eigenvalues
        eigvals = sqrt(eigvals)

        ! allocateand initialize output matrices
        allocate(sqrtA(n_ao, n_ao), inv_sqrtA(n_ao, n_ao))
        sqrtA = 0.0_rp
        inv_sqrtA = 0.0_rp

        ! construct the square root and inverse square root of A
        do i = 1, n_ao
            call dgemm("N","T", n_ao, n_ao, 1_ip, eigvals(i), eigvecs(:, i), n_ao, &
                       eigvecs(:, i), n_ao, 1.0_rp, sqrtA, n_ao)
            call dgemm("N","T", n_ao, n_ao, 1_ip, 1.0_rp / eigvals(i), eigvecs(:, i), &
                       n_ao, eigvecs(:, i), n_ao, 1.0_rp, inv_sqrtA, n_ao)
        end do

        deallocate(eigvecs, eigvals)

    end subroutine compute_sqrt_and_inv_sqrt

    subroutine solve_symm_linear_system(matrix, rhs, settings, error)
        !
        ! this subroutine solves a symmetric linear system A x = B using LAPACK's DSYSV
        !
        use opentrustregion, only: verbosity_error

        real(rp), intent(in) :: matrix(:, :)
        real(rp), intent(inout) :: rhs(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: temp(:, :), work(:)
        integer(ip), allocatable :: ipiv(:)
        integer(ip) :: n, n_rhs, lwork, info
        character(300) :: msg

        external :: dsysv

        ! initialize error flag
        error = 0

        ! linear system size
        n = size(matrix, 1)

        ! number of right-hand sides
        n_rhs = size(rhs, 2)

        ! make a local copy because DSYSV overwrites the matrix
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate(ipiv(n), work(1))
        call dsysv("U", n, n_rhs, temp, n, ipiv, rhs, n, work, lwork, info)
        if (info /= 0) then
            write (msg, '(A, I0)') "Linear solver failed: Error in DSYSV, info = ", info
            call settings%log(msg, verbosity_error, .true.)
            deallocate(temp, ipiv, work)
            error = 1
            return
        end if
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! actual solve
        call dsysv("U", n, n_rhs, temp, n, ipiv, rhs, n, work, lwork, info)

        ! clean up
        deallocate(temp, ipiv, work)

        if (info /= 0) then
            write (msg, '(A, I0)') "Linear solver failed: Error in DSYSV, info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
        end if

    end subroutine solve_symm_linear_system

    function matrix_exponential(A, error) result(expA)
        !
        ! this function calculates the matrix exponential of a real antisymmetric 
        ! matrix using the scaling and squaring method applied to the Taylor expansion
        ! of the exponential, the scale factor is derived from the Frobenius norm which 
        ! is an upper bound for the spectral norm, convergence is tested against the
        ! last term of the expansion which works because the sum of the Frobenius norms
        ! of two matrices is larger than the Frobenius norm of the sum of both matrices
        !
        use opentrustregion, only: solver_settings_type, verbosity_error

        real(rp), intent(in) :: A(:, :)
        integer(ip), intent(out) :: error
        real(rp), allocatable :: expA(:, :)

        integer(ip) :: n, i, power
        real(rp) :: scale, fac, A_norm
        real(rp), allocatable :: An(:, :), tmp(:, :)
        external :: dgemm

        ! initialize error flag
        error = 0

        ! matrix size
        n = size(A, 1)

        ! workspace allocation
        allocate(An(n, n), tmp(n, n), expA(n, n))

        ! compute Frobenius norm of A
        A_norm = sqrt(sum(A ** 2))

        ! determine scale factor
        power = 3
        if (A_norm > 1.0_rp) then
            power = power + int(ceiling(log(A_norm) / log(2.0_rp)))
        end if
        scale = 2.0_rp ** (-power)

        ! initialize exponential and product of matrices
        expA = 0.0_rp
        An = 0.0_rp
        do i = 1, n
            expA(i, i) = 1.0_rp
            An(i, i) = 1.0_rp
        end do

        ! perform Taylor expansion
        i = 1
        fac = 1.0_rp
        do while (A_norm > 1e-12_rp)
            ! get factorial
            fac = fac / real(i, rp)

            ! multiply another matrix and change scale factor accordingly
            call dgemm('N','N', n, n, n, scale, A, n, An, n, 0.0_rp, tmp, n)
            An = tmp

            ! add next expansion order
            expA = expA + fac * An

            ! convergence check for last expansion order
            A_norm = fac * sqrt(sum(An ** 2))

            ! check if maximum number of iterations is reached
            i = i + 1
            ! check for errors
            if (i > 100) then
                call arh_object%settings%log("Maximum number of iterations for "// &
                                             "Taylor expansion of matrix "// &
                                             "exponential reached.", verbosity_error, &
                                             .true.)
                error = 1
                return
            end if
        end do
        deallocate(An)

        ! squaring step
        do i = 1, power
            call dgemm('N', 'N', n, n, n, 1.0_rp, expA, n, expA, n, 0.0_rp, tmp, n)
            expA = tmp
        end do
        deallocate(tmp)

    end function matrix_exponential

end module otr_arh
