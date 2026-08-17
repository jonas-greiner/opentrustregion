! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao

    use opentrustregion, only: rp, ip, settings_type, obj_func_type, update_orbs_type, &
                               hess_x_type, project_type

    implicit none

    type, extends(settings_type) :: oao_settings_type
        logical :: restricted
    contains
        procedure :: init => init_oao_settings
    end type oao_settings_type

    type(oao_settings_type), parameter :: default_oao_settings = &
        oao_settings_type(logger = null(), initialized = .true., restricted = .false., &
                          verbose = 0)

    abstract interface
        function get_energy_2d_type(dm, error) result(energy)
            import :: rp, ip

            real(rp), intent(in), target, contiguous :: dm(:, :)
            integer(ip), intent(out) :: error

            real(rp) :: energy
        end function get_energy_2d_type
    end interface

    abstract interface
        function get_energy_3d_type(dm, error) result(energy)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            integer(ip), intent(out) :: error

            real(rp) :: energy
        end function get_energy_3d_type
    end interface

    abstract interface
        subroutine get_response_2d_type(dm, response, error)
            import :: rp, ip

            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out), target, contiguous :: response(:, :)
            integer(ip), intent(out) :: error
        end subroutine get_response_2d_type
    end interface

    abstract interface
        subroutine get_response_3d_type(dm, response, error)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out), target :: response(:, :, :)
            integer(ip), intent(out) :: error
        end subroutine get_response_3d_type
    end interface

    abstract interface
        subroutine update_dm_2d_type(dm, energy, fock, get_response_funptr, error)
            import :: rp, ip

            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), target, contiguous :: fock(:, :)
            procedure(get_response_2d_type), intent(out), pointer :: get_response_funptr
            integer(ip), intent(out) :: error
        end subroutine update_dm_2d_type
    end interface

    abstract interface
        subroutine update_dm_3d_type(dm, energy, fock, get_response_funptr, error)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), target :: fock(:, :, :)
            procedure(get_response_3d_type), intent(out), pointer :: get_response_funptr
            integer(ip), intent(out) :: error
        end subroutine update_dm_3d_type
    end interface

    type :: oao_type
        type(oao_settings_type) :: settings
        integer(ip) :: n_ao, n_param, n_particle
        real(rp) :: energy = 0.0_rp
        real(rp), pointer, contiguous :: dm_ao(:, :, :)
        real(rp), allocatable :: s_sqrt(:, :), s_inv_sqrt(:, :), dm_oao(:, :, :), &
                                 fock_oo(:, :, :), fock_vv(:, :, :), grad(:), &
                                 h_diag(:)
        procedure(get_energy_3d_type), pointer, nopass :: get_energy => null()
        procedure(update_dm_3d_type), pointer, nopass :: update_dm => null()
        procedure(get_response_3d_type), pointer, nopass :: get_response => null()
        procedure(get_energy_2d_type), pointer, nopass :: get_energy_2d => null()
        procedure(update_dm_2d_type), pointer, nopass :: update_dm_2d => null()
        procedure(get_response_2d_type), pointer, nopass :: get_response_2d => null()
    end type oao_type

    ! global variables
    type(oao_type), allocatable, target :: oao_object

    ! create function pointers to ensure that routines comply with interface
    procedure(obj_func_type), pointer :: obj_func_oao_ptr => obj_func_oao
    procedure(update_orbs_type), pointer :: update_orbs_oao_ptr => update_orbs_oao
    procedure(hess_x_type), pointer :: hess_x_oao_ptr => hess_x_oao
    procedure(project_type), pointer :: project_oao_ptr => project_oao

    contains

    subroutine oao_factory_closed_shell(dm_ao, ao_overlap, n_particle, n_ao, &
                                        get_energy, update_dm, obj_func_oao_funptr, &
                                        update_orbs_oao_funptr, project_oao_funptr, &
                                        error, settings)
        !
        ! this function returns a modified OAO orbital updating function for the 
        ! closed-shell case
        !
        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_2d_type), intent(in), pointer :: get_energy
        procedure(update_dm_2d_type), intent(in), pointer :: update_dm
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(project_type), intent(out), pointer :: project_oao_funptr
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
        oao_object%get_energy_2d => get_energy
        oao_object%update_dm_2d => update_dm

        ! get pointers to modified function
        obj_func_oao_funptr => obj_func_oao
        update_orbs_oao_funptr => update_orbs_oao
        project_oao_funptr => project_oao

    end subroutine oao_factory_closed_shell

    subroutine oao_factory_open_shell(dm_ao, ao_overlap, n_particle, n_ao, get_energy, &
                                      update_dm, obj_func_oao_funptr, &
                                      update_orbs_oao_funptr, project_oao_funptr, &
                                      error, settings)
        !
        ! this function returns a modified OAO orbital updating function for the 
        ! open-shell case
        !
        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(get_energy_3d_type), intent(in), pointer :: get_energy
        procedure(update_dm_3d_type), intent(in), pointer :: update_dm
        procedure(obj_func_type), intent(out), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_oao_funptr
        procedure(project_type), intent(out), pointer :: project_oao_funptr
        integer(ip), intent(out) :: error
        type(oao_settings_type), intent(inout) :: settings

        ! initialize error flag
        error = 0

        ! call common setup
        call oao_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        if (error /= 0) return

        ! set pointers to functions
        oao_object%get_energy => get_energy
        oao_object%update_dm => update_dm

        ! get pointers to modified function
        obj_func_oao_funptr => obj_func_oao
        update_orbs_oao_funptr => update_orbs_oao
        project_oao_funptr => project_oao

    end subroutine oao_factory_open_shell

    subroutine oao_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        !
        ! this function performs common OAO initialization operations
        !
        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        integer(ip), intent(out) :: error
        class(oao_settings_type), intent(inout) :: settings

        logical :: initialized

        ! perform sanity check
        call oao_sanity_check(settings, n_ao, error)
        if (error /= 0) return

        ! allocate objects
        if (.not. allocated(oao_object)) allocate(oao_object)

        ! determine whether object has been initialized
        initialized = oao_object%settings%initialized

        ! set (potentially new) settings
        oao_object%settings = settings

        ! check if object has not been initialized or has been initialize with wrong 
        ! settings
        if (.not. initialized .or. (initialized .and. &
                                    (oao_object%n_particle /= n_particle .or. &
                                     oao_object%n_ao /= n_ao))) then
            ! deallocate arrays if they are already allocated
            if (allocated(oao_object%s_sqrt)) deallocate(oao_object%s_sqrt)
            if (allocated(oao_object%s_inv_sqrt)) deallocate(oao_object%s_inv_sqrt)
            if (allocated(oao_object%dm_oao)) deallocate(oao_object%dm_oao)
            if (allocated(oao_object%fock_oo)) deallocate(oao_object%fock_oo)
            if (allocated(oao_object%fock_vv)) deallocate(oao_object%fock_vv)

            ! number of particles
            oao_object%n_particle = n_particle

            ! get number of atomic orbitals
            oao_object%n_ao = n_ao
            
            ! starting density matrix
            oao_object%dm_ao => dm_ao

            ! get square root and inverse square root of AO overlap matrix
            call compute_sqrt_and_inv_sqrt(ao_overlap, oao_object%s_sqrt, &
                                           oao_object%s_inv_sqrt, error)
            if (error /= 0) return

            ! get per spin contribution to density matrix in orthogonalized AO basis
            oao_object%dm_oao = symmetric_transformation(oao_object%s_sqrt, dm_ao)

            ! allocate matrices
            allocate(oao_object%fock_oo(n_ao, n_ao, oao_object%n_particle), &
                     oao_object%fock_vv(n_ao, n_ao, oao_object%n_particle))
        end if

        ! get number of non-redundant parameters
        oao_object%n_param = n_ao * (n_ao - 1) / 2
        if (.not. oao_object%settings%restricted) &
            oao_object%n_param = n_particle * oao_object%n_param

    end subroutine oao_factory_common

    subroutine oao_sanity_check(settings, n_ao, error)
        !
        ! this subroutine performs a sanity check for OAO input parameters
        !
        use opentrustregion, only: verbosity_error, string_to_lowercase

        class(oao_settings_type), intent(inout) :: settings
        integer(ip), intent(in) :: n_ao
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

    end subroutine oao_sanity_check

    function obj_func_oao(kappa, error) result(energy)
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
        allocate(rot_dm_ao(oao_object%n_ao, oao_object%n_ao, oao_object%n_particle))
        call rotate_dm_ao(kappa, oao_object%n_particle, oao_object%n_ao, &
                          oao_object%settings%restricted, rot_dm_ao, error)
        if (error /= 0) return

        ! calculate mean-field energy
        if (associated(oao_object%get_energy)) then
            energy = oao_object%get_energy(rot_dm_ao, error)
        else
            energy = oao_object%get_energy_2d(rot_dm_ao(:, :, 1), error)
        end if
        if (error /= 0) return

    end function obj_func_oao

    subroutine update_orbs_oao(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation 
        ! and the Hessian linear transformation in the OAO basis
        !
        use opentrustregion, only: hess_x_type, numerical_zero

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle
        real(rp), allocatable :: fock_ao(:, :, :), fock_oao(:, :, :)
        external :: dgemm

        ! initialize error flag
        error = 0

        ! check if orbitals are actually rotated
        if ((sum(abs(kappa)) > 0.0_rp) .or. &
            (abs(oao_object%energy) <= numerical_zero) .or. &
            (.not. (allocated(oao_object%grad) .and. allocated(oao_object%h_diag) &
                    .and. (associated(oao_object%get_response_2d) .or. &
                           associated(oao_object%get_response))))) then
            ! number of AOs
            n_ao = oao_object%n_ao

            ! number of particles
            n_particle = oao_object%n_particle

            ! rotate density matrix
            call rotate_dm_ao(kappa, n_particle, n_ao, oao_object%settings%restricted, &
                              oao_object%dm_ao, error, oao_object%dm_oao)
            if (error /= 0) return

            ! get energy, Fock matrix, and response function
            allocate(fock_ao(n_ao, n_ao, n_particle))
            if (associated(oao_object%update_dm)) then
                call oao_object%update_dm(oao_object%dm_ao, oao_object%energy, &
                                          fock_ao, oao_object%get_response, error)
            else
                call oao_object%update_dm_2d(oao_object%dm_ao(:, :, 1), &
                                             oao_object%energy, fock_ao(:, :, 1), &
                                             oao_object%get_response_2d, error)
            end if
            if (error /= 0) then
                deallocate(fock_ao)
                return
            end if

            ! transform Fock matrix to OAO basis
            fock_oao = symmetric_transformation(oao_object%s_inv_sqrt, fock_ao)
            deallocate(fock_ao)

            ! calculate gradient and Hessian diagonal
            if (.not. allocated(oao_object%grad)) &
                allocate(oao_object%grad(oao_object%n_param))
            if (.not. allocated(oao_object%h_diag)) &
                allocate(oao_object%h_diag(oao_object%n_param))
            call calculate_grad_h_diag(oao_object%dm_oao, fock_oao, n_particle, n_ao, &
                                       oao_object%settings%restricted, &
                                       oao_object%grad, oao_object%h_diag, &
                                       oao_object%fock_oo, oao_object%fock_vv)
            deallocate(fock_oao)
        end if

        ! set outputs
        func = oao_object%energy
        grad = oao_object%grad
        h_diag = oao_object%h_diag
        hess_x_funptr => hess_x_oao
        
    end subroutine update_orbs_oao

    subroutine hess_x_oao(x, hess_x, error)
        !
        ! this function defines the Hessian linear transformation in the orthogonal AO 
        ! basis
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, n_param, i
        real(rp), allocatable :: x_full(:, :, :), dm_response(:, :, :), &
                                 fock_response(:, :, :), hess_x_full(:, :, :)

        external :: dgemm

        ! initialize error flag
        error = 0

        ! number of AOs
        n_ao = oao_object%n_ao

        ! number of particles
        n_particle = oao_object%n_particle

        ! number of parameters
        n_param = oao_object%n_param

        ! unpack trial vector
        x_full = unpack_asymm(x, n_particle, n_ao, oao_object%settings%restricted)

        ! for ROHF, we must explicitly project the trial vector to the valid [c-v], 
        ! [o-v], and [c-o] rotation blocks. Because we reuse UHF machinery where 
        ! alpha-occ = [c+o] and beta-occ = [c], a unified spatial trial vector that is 
        ! unpacked into alpha and beta spin channels contains rotations that are 
        ! internal and thus redundant in the UHF context. Without this projection, 
        ! the alpha Hessian vector product would leak non-zero values into 
        ! the [c-c] block and the beta Hessian vector product into the [v-v] blocks. 
        ! This projection thus prevents singular dimensions in the solver.
        if (oao_object%settings%restricted .and. oao_object%n_particle > 1) &
            x_full = project(x_full, oao_object%dm_oao)

        ! get one electron part
        allocate(hess_x_full(n_ao, n_ao, n_particle))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, oao_object%fock_vv(:, :, i) &
                       - oao_object%fock_oo(:, :, i), n_ao, x_full(:, :, i), n_ao, &
                       0.0_rp, hess_x_full(:, :, i), n_ao)
            hess_x_full(:, :, i) = hess_x_full(:, :, i) - &
                                   transpose(hess_x_full(:, :, i))
        end do

        ! get commutator of trial vector and current density matrix
        allocate(dm_response(n_ao, n_ao, n_particle))
        do i = 1, n_particle
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, oao_object%dm_oao(:, :, i), &
                       n_ao, x_full(:, :, i), n_ao, 0.0_rp, dm_response(:, :, i), n_ao)
            dm_response(:, :, i) = dm_response(:, :, i) + &
                                   transpose(dm_response(:, :, i))
        end do
        deallocate(x_full)

        ! transform density matrix from OAO basis to AO basis
        dm_response = symmetric_transformation(oao_object%s_inv_sqrt, dm_response)

        ! get response of Fock matrix to density matrix response
        allocate(fock_response(n_ao, n_ao, n_particle))
        if (associated(oao_object%get_response)) then
            call oao_object%get_response(dm_response, fock_response, error)
        else
            call oao_object%get_response_2d(dm_response(:, :, 1), &
                                            fock_response(:, :, 1), error)
        end if
        if (error /= 0) return

        ! transform Fock response to OAO basis
        fock_response = symmetric_transformation(oao_object%s_inv_sqrt, fock_response)

        ! project function (is this necessary?)
        hess_x_full = hess_x_full + project(fock_response, oao_object%dm_oao)
        deallocate(fock_response)
        if (error /= 0) return

        ! pack Hessian linear transformation
        hess_x = pack_asymm(hess_x_full, oao_object%n_param, &
                            oao_object%settings%restricted)
        deallocate(hess_x_full)

    end subroutine hess_x_oao

    subroutine project_oao(vector, error)
        !
        ! this subroutine projects out the occupied-occupied and virtual-virtual 
        ! contributions from a vector in-place
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: vector_full(:, :, :), projected_vector_full(:, :, :)

        ! initialize error flag
        error = 0

        ! unpack vector
        vector_full = unpack_asymm(vector, oao_object%n_particle, &
                                   oao_object%n_ao, oao_object%settings%restricted)

        ! project out o-o and v-v contributions
        projected_vector_full = project(vector_full, oao_object%dm_oao)
        deallocate(vector_full)

        ! pack vector
        vector = pack_asymm(projected_vector_full, size(vector, kind=ip), &
                            oao_object%settings%restricted)
        deallocate(projected_vector_full)

    end subroutine project_oao

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
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, oao_object%dm_oao(:, :, i), &
                       n_ao, u_ptr, n_ao, 0.0_rp, temp, n_ao)
            call dgemm("T", "N", n_ao, n_ao, n_ao, 1.0_rp, u_ptr, n_ao, temp, n_ao, &
                       0.0_rp, rot_dm_oao_ptr(:, :, i), n_ao)
        end do

        ! purify density matrix
        call purify(rot_dm_oao_ptr)

        ! transform density matrix from OAO basis to AO basis
        rot_dm_ao = symmetric_transformation(oao_object%s_inv_sqrt, rot_dm_oao_ptr)

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
        grad = pack_asymm(grad_full, size(grad, kind=ip), restricted)

        ! construct Hessian diagonal
        idx = 1
        if (restricted) then
            do j = 1, n_ao
                do i = 1, j - 1
                    h_diag(idx) = sum(fock_vv(i, i, :) + fock_vv(j, j, :) - &
                                      fock_oo(i, i, :) - fock_oo(j, j, :)) / n_particle
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

        ! initialize full matrix
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

    function pack_asymm(matrix, n_param, spin_sum) result(matrix_nonred)
        !
        ! this function packs an antisymmetric matrix for RHF while deallocating the 
        ! original unpacked matrix
        !
        real(rp), intent(in) :: matrix(:, :, :)
        integer(ip), intent(in) :: n_param
        logical, intent(in) :: spin_sum
        real(rp), allocatable :: matrix_nonred(:)

        integer(ip) :: idx, i, j, k

        ! allocate redundant matrix
        allocate(matrix_nonred(n_param))

        ! initialize index
        idx = 1

        ! spin-restricted case: sum over both spins
        if (spin_sum) then
            do j = 1, size(matrix, 2)
                do i = 1, j - 1
                    matrix_nonred(idx) = sum(matrix(i, j, :)) / size(matrix, 3)
                    idx = idx + 1
                end do
            end do
        ! spin-unrestricted case: pack each spin separately
        else
            do k = 1, size(matrix, 3)
                do j = 1, size(matrix, 2)
                    do i = 1, j - 1
                        matrix_nonred(idx) = matrix(i, j, k)
                        idx = idx + 1
                    end do
                end do
            end do
        end if

    end function pack_asymm

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
            call oao_object%settings%log(msg, verbosity_error, .true.)
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
                call oao_object%settings%log("Maximum number of iterations for "// &
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

end module otr_oao
