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
    ! Implement for ROHF and UHF

    use opentrustregion, only: rp, ip, settings_type, obj_func_type, update_orbs_type, &
                               hess_x_type, precond_type

    implicit none

    type, extends(settings_type) :: arh_settings_type
    contains
        procedure :: init => init_arh_settings
    end type arh_settings_type

    type(arh_settings_type), parameter :: default_arh_settings = &
        arh_settings_type(logger = null(), initialized = .true., verbose = 0)

    abstract interface
        function get_energy_type(dm, error) result(energy)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :)
            integer(ip), intent(out) :: error

            real(rp) :: energy
        end function get_energy_type
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

    type :: arh_type
        type(arh_settings_type) :: settings
        integer(ip) :: n_ao, n_param
        real(rp), allocatable :: dm_ao(:, :), s_sqrt(:, :), s_inv_sqrt(:, :), &
                                 dm_per_spin_oao(:, :), fock_oao(:, :), fock_oo(:, :), &
                                 fock_vv(:, :), metric(:, :), dm_list(:, :, :), &
                                 fock_list(:, :, :), dm_diff(:, :, :), &
                                 fock_diff(:, :, :), h_diag_test(:)
        procedure(get_energy_type), pointer, nopass :: get_energy
        procedure(get_fock_type), pointer, nopass :: get_fock
    end type arh_type

    ! global variables
    type(arh_type) :: arh_object

    ! create function pointers to ensure that routines comply with interface
    procedure(obj_func_type), pointer :: obj_func_arh_ptr => obj_func_arh
    procedure(update_orbs_type), pointer :: update_orbs_arh_ptr => update_orbs_arh
    procedure(hess_x_type), pointer :: hess_x_arh_ptr => hess_x_arh
    procedure(precond_type), pointer :: precond_arh_ptr => &
        level_shifted_diag_precond_arh

    contains

    subroutine arh_factory(dm_ao, ao_overlap, n_ao, get_energy, get_fock, &
                           obj_func_arh_funptr, update_orbs_arh_funptr, &
                           precond_arh_funptr, error, settings)
        !
        ! this function returns a modified ARH orbital updating function
        !
        real(rp), intent(in) :: dm_ao(:, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_ao
        procedure(get_energy_type), intent(in), pointer :: get_energy
        procedure(get_fock_type), intent(in), pointer :: get_fock
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(precond_type), intent(out), pointer :: precond_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

        ! set settings
        arh_object%settings = settings

        ! get square root and inverse square root of AO overlap matrix
        call compute_sqrt_and_inv_sqrt(ao_overlap, arh_object%s_sqrt, &
                                       arh_object%s_inv_sqrt, error)
        if (error /= 0) return

        ! get per spin contribution to density matrix in orthogonalized AO basis
        arh_object%dm_per_spin_oao = get_dm_per_spin_oao(dm_ao, arh_object%s_sqrt)

        ! get number of non-redundant parameters
        arh_object%n_ao = n_ao
        arh_object%n_param = n_ao * (n_ao - 1) / 2

        ! starting density matrix
        arh_object%dm_ao = dm_ao

        ! allocate matrices
        allocate(arh_object%fock_oo(n_ao, n_ao), &
                 arh_object%h_diag_test(arh_object%n_param))

        ! set pointers to functions
        arh_object%get_energy => get_energy
        arh_object%get_fock => get_fock

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh
        update_orbs_arh_funptr => update_orbs_arh
        precond_arh_funptr => level_shifted_diag_precond_arh

    end subroutine arh_factory

    function obj_func_arh(kappa, error) result(energy)
        !
        ! this function defines the energy evaluation in OAO basis
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        real(rp), allocatable :: u(:, :), rot_dm_per_spin_oao(:, :), rot_dm_ao(:, :)

        ! initialize error flag
        error = 0

        ! initialize energy in case of error
        energy = 0.0_rp

        ! get rotation matrix
        u = matrix_exponential(unpack_rhf(kappa), error)
        if (error /= 0) return

        ! rotate density matrix
        rot_dm_per_spin_oao = rotate_dm(arh_object%dm_per_spin_oao, u)

        ! purify density matrix
        rot_dm_per_spin_oao = purify(rot_dm_per_spin_oao)

        ! transform density matrix from OAO basis to AO basis
        rot_dm_ao = get_dm_ao(rot_dm_per_spin_oao, arh_object%s_inv_sqrt)

        ! calculate mean-field energy
        energy = arh_object%get_energy(rot_dm_ao, error)
        if (error /= 0) return

    end function obj_func_arh

    subroutine update_orbs_arh(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal 
        ! evaluation in the OAO basis and the Hessian linear transformation on the 
        ! basis of augmented Roothaan-Hall
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, i, j, idx, n_list
        real(rp), allocatable :: u(:, :), fock_ao(:, :), &
                                 dm_per_spin_oao_fock_oao(:, :), &
                                 fock_oao_dm_per_spin_oao(:, :), fock_ov(:, :), &
                                 fock_vo(:, :), grad_full(:, :)
        external :: dgemm

        ! number of AOs
        n_ao = arh_object%n_ao

        ! get rotation matrix
        u = matrix_exponential(unpack_rhf(kappa), error)
        if (error /= 0) return

        ! update list of density and Fock matrices
        if (allocated(arh_object%dm_list)) then
            call append(arh_object%dm_list, arh_object%dm_per_spin_oao)
            call append(arh_object%fock_list, arh_object%fock_oao)
        else
            allocate(arh_object%dm_list(n_ao, n_ao, 0), &
                     arh_object%fock_list(n_ao, n_ao, 0))
        end if

        ! rotate density matrix
        arh_object%dm_per_spin_oao = rotate_dm(arh_object%dm_per_spin_oao, u)

        ! purify density matrix
        arh_object%dm_per_spin_oao = purify(arh_object%dm_per_spin_oao)

        ! transform density matrix from OAO basis to AO basis
        arh_object%dm_ao = get_dm_ao(arh_object%dm_per_spin_oao, arh_object%s_inv_sqrt)

        ! get mean-field energy and Fock matrix
        allocate(fock_ao(n_ao, n_ao))
        call arh_object%get_fock(arh_object%dm_ao, func, fock_ao, error)
        if (error /= 0) then
            deallocate(fock_ao)
            return
        end if

        ! transform Fock matrix to OAO basis
        arh_object%fock_oao = symmetric_transformation(arh_object%s_inv_sqrt, fock_ao)
        deallocate(fock_ao)

        ! get contributions to Fock matrix based on occupancies
        allocate(dm_per_spin_oao_fock_oao(n_ao, n_ao), &
                 fock_oao_dm_per_spin_oao(n_ao, n_ao))
        call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, arh_object%dm_per_spin_oao, &
                   n_ao, arh_object%fock_oao, n_ao, 0.0_rp, dm_per_spin_oao_fock_oao, &
                   n_ao)
        call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_per_spin_oao_fock_oao, n_ao, &
                   arh_object%dm_per_spin_oao, n_ao, 0.0_rp, arh_object%fock_oo, &
                   n_ao) ! DFD
        fock_ov = dm_per_spin_oao_fock_oao - arh_object%fock_oo ! DF(I-D)
        fock_vo = transpose(fock_ov) ! (I_D)FD
        call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, arh_object%fock_oao, n_ao, &
                   arh_object%dm_per_spin_oao, n_ao, 0.0_rp, fock_oao_dm_per_spin_oao, &
                   n_ao)
        arh_object%fock_vv = arh_object%fock_oao - dm_per_spin_oao_fock_oao - &
                  fock_oao_dm_per_spin_oao + arh_object%fock_oo ! (I_D)F(I_D)
        deallocate(dm_per_spin_oao_fock_oao, fock_oao_dm_per_spin_oao)

        ! construct gradient
        grad_full = fock_ov - fock_vo
        deallocate(fock_ov, fock_vo)

        ! pack gradient
        call pack_rhf(grad_full, grad)

        ! construct Hessian diagonal
        idx = 1
        do j = 1, n_ao
            do i = 1, j - 1
                h_diag(idx) = arh_object%fock_vv(i, i) + arh_object%fock_vv(j, j) &
                              - arh_object%fock_oo(i, i) - arh_object%fock_oo(j, j)
                idx = idx + 1
            end do
        end do
        arh_object%h_diag_test = h_diag

        ! construct ARH metric
        arh_object%metric = get_arh_metric(arh_object%dm_list, &
                                           arh_object%dm_per_spin_oao)

        ! prepare differences for two-electron part of Hessian
        n_list = size(arh_object%dm_list, 3)
        if (allocated(arh_object%dm_diff)) deallocate(arh_object%dm_diff, &
                                                      arh_object%fock_diff)
        allocate(arh_object%dm_diff(n_ao, n_ao, n_list), &
                 arh_object%fock_diff(n_ao, n_ao, n_list))
        do i = 1, n_list
            arh_object%dm_diff(:, :, i) = arh_object%dm_list(:, :, i) - &
                                          arh_object%dm_per_spin_oao
            arh_object%fock_diff(:, :, i) = arh_object%fock_list(:, :, i) - &
                                            arh_object%fock_oao
        end do

        ! define pointer to ARH Hessian linear transformation function
        hess_x_funptr => hess_x_arh
        
    end subroutine update_orbs_arh

    subroutine hess_x_arh(x, hess_x, error)
        !
        ! this function defines the Hessian linear transformation on the basis of 
        ! augmented Roothaan-Hall
        !
        use opentrustregion, only: verbosity_error

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_diff, i, lwork, info
        real(rp), allocatable :: x_full(:, :), dm_per_spin_oao_x(:, :), vec(:), &
                                 two_el(:, :), hess_x_full(:, :), temp(:, :), work(:)

        integer(ip), allocatable :: ipiv(:)
        character(300) :: msg
        external :: dgemm, dsysv

        ! initialize error flag
        error = 0

        ! number of AOs
        n_ao = arh_object%n_ao

        ! unpack trial vector
        x_full = unpack_rhf(x)

        ! get one electron part
        allocate(hess_x_full(n_ao, n_ao))
        call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, arh_object%fock_vv - &
                   arh_object%fock_oo, n_ao, x_full, n_ao, 0.0_rp, hess_x_full, n_ao)
        hess_x_full = hess_x_full - transpose(hess_x_full)

        ! number of density matrix differences
        n_diff = size(arh_object%dm_diff, 3)

        ! get two electron part
        if (n_diff > 0) then
            ! contract density matrix with trial vector
            allocate(dm_per_spin_oao_x(n_ao, n_ao))
            call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, arh_object%dm_per_spin_oao, &
                       n_ao, x_full, n_ao, 0.0_rp, dm_per_spin_oao_x, n_ao)
            dm_per_spin_oao_x = dm_per_spin_oao_x + transpose(dm_per_spin_oao_x)

            ! get vector containing density matrix difference
            allocate(vec(n_diff))
            do i = 1, n_diff
                vec(i) = sum(arh_object%dm_diff(:, :, i) * dm_per_spin_oao_x)
            end do
            deallocate(dm_per_spin_oao_x)

            ! multiply inverse metric with vector through linear solver, query 
            ! optimal workspace size
            temp = arh_object%metric
            lwork = -1
            allocate(ipiv(n_diff), work(1))
            call dsysv("U", n_diff, 1_ip, temp, n_diff, ipiv, vec, n_diff, work, &
                       lwork, info)
            lwork = int(work(1))
            deallocate(work)
            allocate(work(lwork))

            ! solve linear system
            call dsysv("U", n_diff, 1_ip, temp, n_diff, ipiv, vec, n_diff, work, &
                       lwork, info)

            ! deallocatework array
            deallocate(temp, ipiv, work)

            ! check for errors
            if (info /= 0) then
                write (msg, '(A, I0)') "Linear solver failed: Error in DSYSV, "// &
                                        "info = ", info
                call arh_object%settings%log(msg, verbosity_error, .true.)
                error = 0
                return
            end if

            ! contract with fock matrix differences
            allocate(two_el(n_ao, n_ao))
            two_el = 0.0_rp
            do i = 1, n_diff
                two_el = two_el + vec(i) * arh_object%fock_diff(:, :, i)
            end do
            deallocate(vec)

            ! add only v-o and o-v contributions of ARH two-electron part
            hess_x_full = hess_x_full + project(two_el, arh_object%dm_per_spin_oao)
            deallocate(two_el)
        end if
        deallocate(x_full)

        ! pack Hessian linear transformation
        call pack_rhf(hess_x_full, hess_x)

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

        real(rp), allocatable :: precond_arr(:),  precond_vector_full(:, :)

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
        precond_vector_full = unpack_rhf(precond_vector)

        ! extract only v-o and o-v contributions
        precond_vector_full = project(precond_vector_full, arh_object%dm_per_spin_oao)

        ! pack vector
        call pack_rhf(precond_vector_full, precond_vector)
        
    end subroutine level_shifted_diag_precond_arh

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

    subroutine arh_deconstructor(dm_ao, error)
        !
        ! this subroutine deallocates the ARH objects
        !
        use opentrustregion, only: verbosity_error

        real(rp), intent(out) :: dm_ao(:, :)
        integer(ip), intent(out) :: error

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

        deallocate(arh_object%dm_ao, arh_object%s_sqrt, arh_object%s_inv_sqrt, &
                   arh_object%dm_per_spin_oao, arh_object%fock_oao, &
                   arh_object%fock_oo, arh_object%fock_vv, arh_object%dm_list, &
                   arh_object%fock_list, arh_object%dm_diff, arh_object%fock_diff)

    end subroutine arh_deconstructor

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

    function unpack_rhf(matrix_nonred) result(matrix)
        !
        ! this function unpacks an antisymmetric matrix for RHF and returns the 
        ! resulting unpacked matrix
        !
        real(rp), intent(in) :: matrix_nonred(:)
        
        real(rp), allocatable :: matrix(:, :)

        integer(ip) :: n_ao, i, j, idx

        ! get number of AOs
        n_ao = int((1.0_rp + sqrt(1.0_rp + 8.0_rp * size(matrix_nonred))) / 2.0_rp)

        ! allocateand initialize full matrix
        allocate(matrix(n_ao, n_ao))
        matrix = 0.0_rp

        ! construct antisymmetric matrix
        idx = 1
        do j = 1, n_ao
            do i = 1, j - 1
                matrix(i, j) = matrix_nonred(idx)
                matrix(j, i) = -matrix_nonred(idx)
                idx = idx + 1
            end do
        end do

    end function unpack_rhf

    subroutine pack_rhf(matrix, matrix_nonred)
        !
        ! this subroutine packs an antisymmetric matrix for RHF while deallocating the 
        ! original unpacked matrix
        !
        real(rp), intent(inout), allocatable :: matrix(:, :)
        real(rp), intent(inout) :: matrix_nonred(:)

        integer(ip) :: n_ao, i, j, idx

        ! get number of AOs
        n_ao = size(matrix, 1)

        ! pack vector
        idx = 1
        do j = 1, n_ao
            do i = 1, j - 1
                matrix_nonred(idx) = matrix(i, j)
                idx = idx + 1
            end do
        end do
        deallocate(matrix)

    end subroutine pack_rhf

    function rotate_dm(dm, u) result(rot_dm)
        !
        ! this function rotates a density matrix (dm_rotated = U^T * dm * U)
        !
        real(rp), intent(in) :: dm(:, :), u(:, :)
        real(rp), allocatable :: rot_dm(:, :)

        real(rp), allocatable :: temp(:, :)
        integer(ip) :: n
        external :: dgemm

        ! size
        n = size(dm, 1)

        ! allocatearrays
        allocate(rot_dm(n, n), temp(n, n))

        ! rotate density matrix
        call dgemm("N","N", n, n, n, 1.0_rp, dm, n, u, n, 0.0_rp, temp, n)
        call dgemm("T","N", n, n, n, 1.0_rp, u, n, temp, n, 0.0_rp, rot_dm, n)

        deallocate(temp)

    end function rotate_dm

    function purify(dm) result(dm_purified)
        !
        ! this function purifies a density matrix (dm_purified = 3*dm^2 - 2*dm^3)
        !
        real(rp), intent(in) :: dm(:, :)
        real(rp), allocatable :: dm_purified(:, :)

        real(rp), allocatable :: dm_squared(:, :)
        integer(ip) :: n
        external :: dgemm

        ! size
        n = size(dm, 1)

        ! allocatearrays
        allocate(dm_squared(n, n), dm_purified(n, n))

        ! square density matrix
        call dgemm("N", "N", n, n, n, 1.0_rp, dm, n, dm, n, 0.0_rp, dm_squared, n)

        ! purify density matrix
        call dgemm("N", "N", n, n, n, 2.0_rp, dm_squared, n, dm, n, 0.0_rp, &
                   dm_purified, n)
        dm_purified = 3.0_rp*dm_squared - dm_purified

        deallocate(dm_squared)

    end function purify

    function symmetric_transformation(trans_matrix, matrix) result(matrix_transformed)
        !
        ! this function performs a symmetric transformation (X' = U * X * U)
        !
        real(rp), intent(in) :: trans_matrix(:, :), matrix(:, :)
        real(rp), allocatable :: matrix_transformed(:, :)

        real(rp), allocatable :: temp(:, :)
        integer(ip) :: n
        external :: dgemm

        n = size(matrix, 1)
        allocate(temp(n, n), matrix_transformed(n, n))

        call dgemm("N", "N", n, n, n, 1.0_rp, trans_matrix, n, matrix, n, 0.0_rp, &
                   temp, n)
        call dgemm("N", "N", n, n, n, 1.0_rp, temp, n, trans_matrix, n, 0.0_rp, &
                   matrix_transformed, n)

        deallocate(temp)

    end function symmetric_transformation

    function get_dm_per_spin_oao(dm_ao, s_sqrt) result(dm_per_spin_oao)
        !
        ! this function constructs the density matrix per spin eigenvalue for RHF and 
        ! transforms it from the AO basis to the orthogonal AO basis
        !
        real(rp), intent(in) :: dm_ao(:, :), s_sqrt(:, :)
        real(rp), allocatable :: dm_per_spin_oao(:, :)

        allocate(dm_per_spin_oao(size(dm_ao, 1), size(dm_ao, 2)))

        dm_per_spin_oao = symmetric_transformation(s_sqrt, dm_ao / 2)

    end function get_dm_per_spin_oao

    function get_dm_ao(dm_per_spin_oao, s_inv_sqrt) result(dm_ao)
        !
        ! this function constructs the density matrix for RHF and transforms it from 
        ! the orthogonal AO basis to the AO basis
        !
        real(rp), intent(in) :: dm_per_spin_oao(:, :), s_inv_sqrt(:, :)
        real(rp), allocatable :: dm_ao(:, :)

        allocate(dm_ao(size(dm_per_spin_oao, 1), size(dm_per_spin_oao, 2)))

        dm_ao = 2.0_rp * symmetric_transformation(s_inv_sqrt, dm_per_spin_oao)

    end function get_dm_ao

    subroutine append(list, new_matrix)
        !
        ! this subroutine appends a matrix to a list of matrices of equal dimension
        !
        real(rp), intent(inout), allocatable :: list(:, :, :)
        real(rp), intent(in) :: new_matrix(:, :)

        integer(ip) :: nrows, ncols
        real(rp), allocatable :: temp(:, :, :)

        nrows = size(new_matrix, 1)
        ncols = size(new_matrix, 2)

        allocate(temp(nrows, ncols, size(list, 3) + 1))
        temp(:, :, :size(list, 3)) = list
        temp(:, :, size(list, 3) + 1) = new_matrix
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

        ! allocateeigenvector and eigenvalue arrays
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

    function get_arh_metric(dm_list, dm_per_spin_oao) result(metric)
        !
        ! this function calculates the augmented Roothaan-Hall metric
        !
        real(rp), intent(in) :: dm_list(:, :, :)
        real(rp), intent(in) :: dm_per_spin_oao(:, :)
        real(rp), allocatable :: metric(:, :)

        integer(ip) :: n_dm, i, j
        real(rp), allocatable :: delta_i(:, :), delta_j(:, :)
        real(rp) :: trace_val

        n_dm = size(dm_list, 3)

        if (allocated(metric)) deallocate(metric)
        allocate(metric(n_dm, n_dm))

        ! generate ARH metric
        do i = 1, n_dm
            delta_i = dm_list(:, :, i) - dm_per_spin_oao
            do j = 1, i
                delta_j = dm_list(:, :, j) - dm_per_spin_oao
                ! compute Tr(delta_i * delta_j)
                trace_val = sum(delta_i * transpose(delta_j))
                metric(i, j) = trace_val
                metric(j, i) = trace_val
            end do
        end do

    end function get_arh_metric

    function project(matrix, dm_per_spin_oao) result(projected_matrix)
        !
        ! this function only retains occupied-virtual and virtual-occupied 
        ! contributions to a matrix
        !
        real(rp), intent(in) :: matrix(:, :), dm_per_spin_oao(:, :)
        real(rp), allocatable :: projected_matrix(:, :)

        integer(ip) :: n_ao, i
        real(rp), allocatable :: proj_v(:, :), temp(:, :)
        external :: dgemm

        n_ao = size(matrix, 1)

        allocate(proj_v(n_ao, n_ao))

        ! construct projection matrix on virtual space (I-D)
        proj_v = 0.0_rp
        do i = 1, n_ao
            proj_v(i, i) = 1.0_rp
        end do
        proj_v = proj_v - dm_per_spin_oao

        ! construct virtual-occupied contributions DM(I-D)
        allocate(temp(n_ao, n_ao))
        call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, matrix, n_ao, proj_v, n_ao, &
                   0.0_rp, temp, n_ao)
        deallocate(proj_v)
        allocate(projected_matrix(n_ao, n_ao))
        call dgemm("N", "N", n_ao, n_ao, n_ao, 1.0_rp, dm_per_spin_oao, n_ao, temp, &
                   n_ao, 0.0_rp, projected_matrix, n_ao)
        deallocate(temp)

        ! add occupied-virtual contributions (I_D)MD
        projected_matrix = projected_matrix - transpose(projected_matrix)

    end function project

end module otr_arh
