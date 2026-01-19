! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn

    ! Things to do:
    ! Implement for other manifolds

    use opentrustregion, only: rp, ip, kw_len, settings_type, update_orbs_type, &
                               hess_x_type
    use otr_common, only: change_reference_type

    implicit none

    type, extends(settings_type) :: qn_settings_type
        character(kw_len) :: hess_update_scheme
    contains
        procedure :: init => init_qn_settings
    end type qn_settings_type

    type(qn_settings_type), parameter :: default_qn_settings = &
        qn_settings_type(logger = null(), initialized = .true., &
                         hess_update_scheme = "sr1", verbose = 0)

    type, abstract :: updating_type
        real(rp), allocatable :: kappa_list(:, :), local_grad_list(:, :), &
                                 grad_list(:, :), h_diag(:)
        integer(ip) :: n_points = 0
        type(qn_settings_type) :: settings
    contains
        procedure(init_interface), deferred :: init
        procedure :: add
        procedure(clear_interface), deferred :: clear
    end type updating_type

    type, extends(updating_type) :: sr1_updating_type
    contains
        procedure :: init => sr1_init
        procedure :: clear => sr1_clear
    end type sr1_updating_type

    type, extends(updating_type) :: bfgs_updating_type
        real(rp), allocatable :: y_list(:, :)
    contains
        procedure :: init => bfgs_init
        procedure :: clear => bfgs_clear
    end type bfgs_updating_type

    abstract interface
        subroutine init_interface(self, n_param)
            import updating_type, ip
            class(updating_type), intent(inout) :: self
            integer(ip), intent(in) :: n_param
        end subroutine init_interface

        subroutine clear_interface(self)
            import updating_type
            class(updating_type), intent(inout) :: self
        end subroutine clear_interface
    end interface

    ! global variables
    class(updating_type), pointer :: update_object
    type(sr1_updating_type), target :: sr1_object
    type(bfgs_updating_type), target :: bfgs_object
    procedure(update_orbs_type), pointer :: update_orbs_orig_funptr
    procedure(hess_x_type), pointer :: hess_x_qn_funptr
    procedure(change_reference_type), pointer :: change_reference_funptr

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_qn_ptr => update_orbs_qn

    contains

    subroutine update_orbs_qn_factory(update_orbs_funptr_in, &
                                      change_reference_funptr_in, n_param, settings, &
                                      error, update_orbs_qn_funptr)
        !
        ! this subroutine returns a modified quasi-Newton orbital updating function
        !
        use opentrustregion, only: verbosity_error

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr_in
        procedure(change_reference_type), intent(in), pointer :: &
            change_reference_funptr_in
        integer(ip), intent(in) :: n_param
        type(qn_settings_type), intent(inout) :: settings
        integer(ip), intent(out) :: error
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_qn_funptr

        ! initialize error flag
        error = 0

        ! get object for quasi-Newton and corresponding Hessian linear transformation 
        ! function
        if (settings%hess_update_scheme == "sr1") then
            update_object => sr1_object
            hess_x_qn_funptr => sr1_hess_x_fun
        else if (settings%hess_update_scheme == "bfgs") then
            update_object => bfgs_object
            hess_x_qn_funptr => bfgs_hess_x_fun
        else
            call settings%log("Quasi-Newton updating scheme not implemented.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! initialize updating object
        call update_object%init(n_param)

        ! set settings
        update_object%settings = settings

        ! initialize settings
        if (.not. update_object%settings%initialized) then
            call update_object%settings%init(error)
            if (error /= 0) return
        end if

        ! set pointer to original orbital updating function
        update_orbs_orig_funptr => update_orbs_funptr_in

        ! set pointer to change of reference function
        change_reference_funptr => change_reference_funptr_in

        ! get pointer to modified orbital updating function
        update_orbs_qn_funptr => update_orbs_qn

    end subroutine update_orbs_qn_factory

    subroutine update_orbs_qn(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this subroutine is a modified quasi-Newton orbital updating function
        !
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! update orbitals
        call update_orbs_orig_funptr(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) return

        ! get pointer to modified Hessian linear transformation function
        hess_x_funptr => hess_x_qn_funptr

        ! add new step
        call update_object%add(kappa, grad, error)
        if (error /= 0) return

        ! set Hessian diagonal
        update_object%h_diag = h_diag
        
    end subroutine update_orbs_qn

    subroutine init_qn_settings(self, error)
        !
        ! this subroutine initializes the quasi-Newton settings
        !
        use opentrustregion, only: verbosity_error

        class(qn_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        select type(settings => self)
        type is (qn_settings_type)
            settings = default_qn_settings
        class default
            call settings%log("Quasi-Newton settings could not be initialized "// &
                              "because initialization routine received the wrong "// &
                              "type. The type qn_settings_type was likely "// &
                              "subclassed without providing an initialization "// &
                              "routine.", verbosity_error, .true.)
            error = 1
        end select

    end subroutine init_qn_settings

    subroutine update_orbs_qn_deconstructor()
        !
        ! this subroutine deallocates the quasi-Newton objects
        !
        if (associated(update_object)) then
            call update_object%clear()
            nullify(update_object)
        end if

    end subroutine update_orbs_qn_deconstructor

    subroutine add(self, kappa, grad, error)
        !
        ! this subroutine adds a new rotation and gradient to the corresponding lists 
        ! and updates the reference accordingly
        !
        class(updating_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), grad(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: list(:, :)

        ! initialize error flag
        error = 0

        ! extend rotations
        allocate(list(size(kappa), self%n_points + 1))
        if (self%n_points > 0) list(:, 1:self%n_points) = &
            self%kappa_list(:, 1:self%n_points)
        list(:, self%n_points + 1) = kappa
        call move_alloc(list, self%kappa_list)

        ! extend local gradients
        allocate(list(size(kappa), self%n_points + 1))
        if (self%n_points > 0) list(:, 1:self%n_points) = &
            self%local_grad_list(:, 1:self%n_points)
        list(:, self%n_points + 1) = grad
        call move_alloc(list, self%local_grad_list)

        ! allocate gradients
        deallocate(self%grad_list)
        allocate(self%grad_list(size(kappa), self%n_points + 1))

        self%n_points = self%n_points + 1

        ! move displacement and gradient to current reference
        call change_reference_funptr(kappa, self%n_points, self%kappa_list, &
                                     self%local_grad_list, self%grad_list, error)
        if (error /= 0) return

    end subroutine add

    subroutine sr1_init(self, n_param)
        ! 
        ! this subroutine initializes the SR1 object
        !
        class(sr1_updating_type), intent(inout) :: self
        integer(ip), intent(in) :: n_param

        ! allocate arrays
        allocate(self%kappa_list(n_param, 0))
        allocate(self%local_grad_list(n_param, 0))
        allocate(self%grad_list(n_param, 0))

    end subroutine sr1_init

    subroutine sr1_clear(self)
        ! 
        ! this subroutine deallocates the SR1 object
        !
        class(sr1_updating_type), intent(inout) :: self

        deallocate(self%kappa_list)
        deallocate(self%local_grad_list)
        deallocate(self%grad_list)
        self%n_points = 0

    end subroutine sr1_clear

    subroutine sr1_hess_x_fun(x, hess_x, error)
        !
        ! this subroutines computes the product of a symmetric rank-1 updated Hessian
        ! with a trial vector, this cannot use the two-loop recursion like BFGS, so it
        ! uses L-SR1 instead
        !
        use opentrustregion, only: numerical_zero, verbosity_error

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, problem_dim, i, info, lwork
        real(rp), allocatable :: kappa_diff(:, :), psi(:, :), U(:, :), q(:), z(:), &
                                 work(:)
        integer(ip), allocatable :: ipiv(:)
        character(300) :: msg
        external :: dgemv, dgemm, dgetrf, dgetri
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! get number of parameters
        n_param = size(x)

        ! check if trial vector vanishes
        if (abs(ddot(n_param, x, 1_ip, x, 1_ip)) < numerical_zero) then
            hess_x = 0.0_rp
            return
        end if

        ! initialize with Hessian diagonal
        hess_x = sr1_object%h_diag * x

        ! get problem dimension
        problem_dim = sr1_object%n_points - 1

        ! only current point is available
        if (problem_dim == 0) return

        ! build psi = grad_diff - h_diag * kappa_diff
        allocate(kappa_diff(n_param, problem_dim), psi(n_param, problem_dim))
        do i = 1, problem_dim
            kappa_diff(:, i) = sr1_object%kappa_list(:, i + 1) - &
                sr1_object%kappa_list(:, i)
            psi(:, i) = sr1_object%grad_list(:, i + 1) - sr1_object%grad_list(:, i) - &
                sr1_object%h_diag * kappa_diff(:, i)
        end do

        ! build U = psi^T * kappa_diff
        allocate(U(problem_dim, problem_dim))
        call dgemm('T','N', problem_dim, problem_dim, n_param, 1.0_rp, psi, n_param, &
                   kappa_diff, n_param, 0.0_rp, U, problem_dim)
        deallocate(kappa_diff)

        ! symmetrize U
        U = 0.5_rp * (U + transpose(U))

        ! invert U
        allocate(ipiv(problem_dim))
        call dgetrf(problem_dim, problem_dim, U, problem_dim, ipiv, info)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Matrix inversion failed: Error in DGETRF, "// &
                "info = ", info
            call sr1_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        lwork = -1
        allocate(work(1))
        call dgetri(problem_dim, U, problem_dim, ipiv, work, lwork, info)

        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        call dgetri(problem_dim, U, problem_dim, ipiv, work, lwork, info)
        deallocate(ipiv)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Matrix inversion failed: Error in DGETRI, "// &
                "info = ", info
            call sr1_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! q = psi^T * x
        allocate(q(problem_dim))
        call dgemv('T', n_param, problem_dim, 1.0_rp, psi, n_param, x, 1_ip, 0.0_rp, &
                   q, 1_ip)

        ! z = U^-1 * q
        allocate(z(problem_dim))
        call dgemv('N', problem_dim, problem_dim, 1.0_rp, U, problem_dim, q, 1_ip, &
                   0.0_rp, z, 1_ip)
        deallocate(U, q)

        ! w += psi * z
        call dgemv('N', n_param, problem_dim, 1.0_rp, psi, n_param, z, 1_ip, 1.0_rp, &
                   hess_x, 1_ip)

        ! deallocate arrays
        deallocate(psi, z)

    end subroutine sr1_hess_x_fun

    subroutine bfgs_init(self, n_param)
        !
        ! this subroutine initializes the BFGS object
        !
        class(bfgs_updating_type), intent(inout) :: self
        integer(ip), intent(in) :: n_param

        ! allocate arrays
        allocate(self%kappa_list(n_param, 0))
        allocate(self%local_grad_list(n_param, 0))
        allocate(self%grad_list(n_param, 0))
        allocate(self%y_list(n_param, 0))

    end subroutine bfgs_init

    subroutine bfgs_clear(self)
        !
        ! this subroutine deallocates the BFGS object
        !
        class(bfgs_updating_type), intent(inout) :: self

        deallocate(self%kappa_list)
        deallocate(self%local_grad_list)
        deallocate(self%grad_list)
        deallocate(self%y_list)
        self%n_points = 0

    end subroutine bfgs_clear

    subroutine bfgs_hess_x_fun(x, hess_x, error)
        !
        ! this implements the two-loop recursion for the BFGS update of the Hessian 
        ! applied to a vector x according to T.H. Fischer and J. Almloef, JPC 96, 9768 
        ! (1992), doi:10.1021/j100203a036
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, it
        real(rp) :: S(6), T(4)
        logical :: update_y
        real(rp), allocatable :: list(:, :), kappa_diff(:), grad_diff(:), &
                                 curr_kappa_diff(:), curr_grad_diff(:)
        real(rp), parameter :: vanish_denom_thr = 1.0e-9_rp
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(x)

        ! check if trial vector vanishes
        if (abs(ddot(n_param, x, 1_ip, x, 1_ip)) < numerical_zero) then
            hess_x = 0.0_rp
            return
        end if

        ! initialize dot products
        S = 0.0_rp
        T = 0.0_rp

        ! initialize with Hessian diagonal
        hess_x = bfgs_object%h_diag * x

        ! only current point is available
        if (bfgs_object%n_points == 1) return

        ! get current displacement and gradient difference
        curr_kappa_diff = bfgs_object%kappa_list(:, bfgs_object%n_points) - &
                          bfgs_object%kappa_list(:, bfgs_object%n_points - 1)
        curr_grad_diff = bfgs_object%grad_list(:, bfgs_object%n_points) - &
                         bfgs_object%grad_list(:, bfgs_object%n_points - 1)

        ! initialize initial Hessian approximation multiplied with displacement
        update_y = .false.
        if (size(bfgs_object%y_list, 2) < bfgs_object%n_points - 1) then
            update_y = .true.
            allocate(list(size(x), bfgs_object%n_points - 1))
            if (bfgs_object%n_points > 2) list(:, :bfgs_object%n_points - 2) = &
                bfgs_object%y_list
            list(:, bfgs_object%n_points - 1) = bfgs_object%h_diag * curr_kappa_diff
            call move_alloc(list, bfgs_object%y_list)
        end if

        ! loop over up to (n-2)th iteration
        do it = 1, bfgs_object%n_points - 2
            ! get displacement and gradient difference for iteration
            kappa_diff = bfgs_object%kappa_list(:, it + 1) - &
                         bfgs_object%kappa_list(:, it)
            grad_diff = bfgs_object%grad_list(:, it + 1) - bfgs_object%grad_list(:, it)

            ! calculate dot products (S(2) is the inverse of the paper)
            S(1) = ddot(n_param, kappa_diff, 1_ip, grad_diff, 1_ip)
            if (abs(S(1)) < vanish_denom_thr) then
                S(1) = 0.0_rp
            else
                S(1) = 1.0_rp / S(1)
            end if
            S(2) = ddot(n_param, kappa_diff, 1_ip, bfgs_object%y_list(:, it), 1_ip)
            if (abs(S(2)) < vanish_denom_thr) then
                S(2) = 1.0_rp / vanish_denom_thr
            else
                S(2) = 1.0_rp / S(2)
            end if
            S(3) = ddot(n_param, grad_diff, 1_ip, x, 1_ip)
            S(4) = ddot(n_param, bfgs_object%y_list(:, it), 1_ip, x, 1_ip)

            ! get dot products for current Hessian approximation multiplied with 
            ! displacement
            if (update_y) then
                S(5) = ddot(n_param, grad_diff, 1_ip, curr_kappa_diff, 1_ip)
                S(6) = ddot(n_param, bfgs_object%y_list(:, it), 1_ip, curr_kappa_diff, &
                            1_ip)
            end if

            ! calculate more dot products
            T(1) = S(1) * S(3)
            T(2) = S(2) * S(4)
            if (update_y) then
                T(3) = S(1) * S(5)
                T(4) = S(2) * S(6)
            end if

            ! calculate current Hessian approximation multiplied with trial vector
            hess_x = hess_x + T(1) * grad_diff - T(2) * bfgs_object%y_list(:, it)

            ! get current Hessian approximation multiplied with displacement
            if (update_y) then
                bfgs_object%y_list(:, bfgs_object%n_points - 1) = &
                    bfgs_object%y_list(:, bfgs_object%n_points - 1) + T(3) * &
                    grad_diff(:) - T(4) * bfgs_object%y_list(:, it)
            end if
        end do

        ! deallocate arrays
        if (allocated(kappa_diff)) deallocate(kappa_diff)
        if (allocated(grad_diff)) deallocate(grad_diff)

        ! calculate dot products
        S(1) = ddot(n_param, curr_kappa_diff, 1_ip, curr_grad_diff, 1_ip)
        if (abs(S(1)) < vanish_denom_thr) then
            S(1) = 0.0_rp
        else
            S(1) = 1.0_rp / S(1)
        end if
        S(2) = ddot(n_param, curr_kappa_diff, 1_ip, &
                    bfgs_object%y_list(:, bfgs_object%n_points - 1), 1_ip)
        if (abs(S(2)) < vanish_denom_thr) then
            S(2) = 1.0_rp / vanish_denom_thr
        else
            S(2) = 1.0_rp / S(2)
        end if
        S(3) = ddot(n_param, curr_grad_diff, 1_ip, x, 1_ip)
        S(4) = ddot(n_param, bfgs_object%y_list(:, bfgs_object%n_points - 1), 1_ip, x, &
                    1_ip)
        T(1) = S(1) * S(3)
        T(2) = S(2) * S(4)

        ! calculate final Hessian approximation multiplied with trial vector
        hess_x = hess_x + T(1) * curr_grad_diff - T(2) * &
            bfgs_object%y_list(:, bfgs_object%n_points - 1)

        ! deallocate arrays
        deallocate(curr_kappa_diff, curr_grad_diff)
    end subroutine bfgs_hess_x_fun

end module otr_qn
