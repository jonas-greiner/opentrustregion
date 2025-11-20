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
        real(rp), allocatable :: start(:), prev_grad(:)
        integer(ip) :: n = 0
    contains
        procedure(init_interface), deferred :: init
        procedure(add_interface), deferred :: add
        procedure(clear_interface), deferred :: clear
    end type updating_type

    type, extends(updating_type) :: sr1_updating_type
        real(rp), allocatable :: u_list(:, :), d_list(:)
    contains
        procedure :: init => sr1_init
        procedure :: add => sr1_add
        procedure :: clear => sr1_clear
    end type sr1_updating_type

    type, extends(updating_type) :: bfgs_updating_type
        real(rp), allocatable :: p_list(:, :), y_list(:, :), d_p_list(:), d_y_list(:)
    contains
        procedure :: init => bfgs_init
        procedure :: add => bfgs_add
        procedure :: clear => bfgs_clear
    end type bfgs_updating_type

    abstract interface
        subroutine init_interface(self, n_param)
            import updating_type, ip
            class(updating_type), intent(inout) :: self
            integer(ip), intent(in) :: n_param
        end subroutine init_interface

        subroutine add_interface(self, step, grad_diff, hess_x_funptr, error)
            import updating_type, rp, hess_x_type, ip
            class(updating_type), intent(inout) :: self
            real(rp), intent(in) :: step(:), grad_diff(:)
            procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
            integer(ip), intent(out) :: error
        end subroutine add_interface

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

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_qn_ptr => update_orbs_qn

    contains

    subroutine update_orbs_qn_factory(update_orbs_funptr, n_param, settings, error, &
                                      update_orbs_qn_funptr)
        !
        ! this function returns a modified quasi-Newton orbital updating function
        !
        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        integer(ip), intent(in) :: n_param
        type(qn_settings_type), intent(inout) :: settings
        integer(ip), intent(out) :: error
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_qn_funptr

        ! initialize error flag
        error = 0

        ! set pointer to original orbital updating function
        update_orbs_orig_funptr => update_orbs_funptr

        ! get pointer to modified orbital updating function
        update_orbs_qn_funptr => update_orbs_qn

        ! initialize settings
        if (.not. settings%initialized) then
            call settings%init(error)
            if (error /= 0) return
        end if

        ! get object for quasi-Newton and coresponding Hessian linear transformation 
        ! function
        if (settings%hess_update_scheme == "sr1") then
            update_object => sr1_object
            hess_x_qn_funptr => sr1_hess_x_fun
        else if (settings%hess_update_scheme == "bfgs") then
            update_object => bfgs_object
            hess_x_qn_funptr => bfgs_hess_x_fun
        else
            call settings%log("Quasi-Newton updating scheme not implemented.", 1, &
                              .true.)
            error = 1
            return
        end if

        ! initialize updating object
        call update_object%init(n_param)

    end subroutine update_orbs_qn_factory

    subroutine update_orbs_qn(kappa, func, grad, h_diag, hess_x_funptr, error)

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        real(rp), allocatable :: step(:), grad_diff(:)

        ! initialize error flag
        error = 0

        ! update orbitals
        call update_orbs_orig_funptr(kappa, func, grad, h_diag, hess_x_funptr, error)
        if (error /= 0) return

        ! get pointer to modified Hessian linear transformation function
        hess_x_funptr => hess_x_qn_funptr

        ! check if updating has been started
        if (allocated(update_object%prev_grad)) then
            ! the step technically needs to be transported
            step = kappa

            ! the previous step technically needs to be transported
            grad_diff = grad - update_object%prev_grad

            ! add new step
            call update_object%add(step, grad_diff, hess_x_funptr, error)
            if (error /= 0) return
        else
            ! initialize start
            update_object%start = max(h_diag, 1d-3)
        end if

        ! update gradient from previous step
        update_object%prev_grad = grad
        
    end subroutine update_orbs_qn

    subroutine init_qn_settings(self, error)
        !
        ! this subroutine initializes the quasi-Newton settings
        !
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
                              "routine.", 1, .true.)
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

    subroutine sr1_init(self, n_param)

        class(sr1_updating_type), intent(inout) :: self
        integer(ip), intent(in) :: n_param

        ! allocate arrays
        allocate(self%u_list(n_param, 0))
        allocate(self%d_list(0))

    end subroutine sr1_init

    subroutine sr1_add(self, step, grad_diff, hess_x_funptr, error)

        class(sr1_updating_type), intent(inout) :: self
        real(rp), intent(in) :: step(:), grad_diff(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        real(rp), allocatable :: u_list(:, :)
        real(rp) :: denom

        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! append new values
        allocate(u_list(size(step), self%n + 1))
        if (self%n > 0) u_list(:, 1:self%n) = self%u_list(:, 1:self%n)
        call hess_x_funptr(step, u_list(:, self%n + 1), error)
        if (error /= 0) return
        u_list(:, self%n + 1) = grad_diff - u_list(:, self%n + 1)
        call move_alloc(u_list, self%u_list)

        denom = ddot(size(step), self%u_list(:, self%n + 1), 1, step, 1)
        if (abs(denom) <= 1d-12) denom = 1d-12
        self%d_list = [self%d_list, denom]
        
        self%n = self%n + 1

    end subroutine sr1_add

    subroutine sr1_clear(self)

        class(sr1_updating_type), intent(inout) :: self

        if (allocated(self%u_list)) deallocate(self%u_list)
        if (allocated(self%d_list)) deallocate(self%d_list)
        self%n = 0

    end subroutine sr1_clear

    subroutine sr1_hess_x_fun(x, hess_x, error)

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: i
        real(rp) :: coeff
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        hess_x = sr1_object%start * x
        do i = 1, sr1_object%n
            coeff = ddot(size(x), sr1_object%u_list(:, i), 1, x, 1) / &
                    sr1_object%d_list(i)
            hess_x = hess_x + coeff * sr1_object%u_list(:, i)
        end do

    end subroutine sr1_hess_x_fun

    subroutine bfgs_init(self, n_param)

        class(bfgs_updating_type), intent(inout) :: self
        integer(ip), intent(in) :: n_param

        ! allocate arrays
        allocate(self%p_list(n_param, 0))
        allocate(self%y_list(n_param, 0))
        allocate(self%d_p_list(0))
        allocate(self%d_y_list(0))

    end subroutine bfgs_init

    subroutine bfgs_add(self, step, grad_diff, hess_x_funptr, error)

        class(bfgs_updating_type), intent(inout) :: self
        real(rp), intent(in) :: step(:), grad_diff(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        real(rp), allocatable :: hess_x_step(:), list(:, :)
        real(rp) :: denom_p, denom_y

        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! check if denominators vanish, if yes skip update
        allocate(hess_x_step(size(step)))
        call hess_x_funptr(step, hess_x_step, error)
        if (error /= 0) return
        denom_p = ddot(size(step), hess_x_step, 1, step, 1)
        denom_y = ddot(size(step), grad_diff, 1, step, 1)
        if (abs(denom_p) <= 1d-12 .or. abs(denom_y) <= 1d-12) then
            deallocate(hess_x_step)
            return
        end if

        ! append new values
        allocate(list(size(step), self%n + 1))
        if (self%n > 0) list(:, 1:self%n) = self%p_list(:, 1:self%n)
        list(:, self%n + 1) = hess_x_step
        deallocate(hess_x_step)
        call move_alloc(list, self%p_list)
        allocate(list(size(step), self%n + 1))
        if (self%n > 0) list(:, 1:self%n) = self%y_list(:, 1:self%n)
        list(:, self%n + 1) = grad_diff
        call move_alloc(list, self%y_list)
        self%d_p_list = [self%d_p_list, denom_p]
        self%d_y_list = [self%d_y_list, denom_y]
        
        self%n = self%n + 1

    end subroutine bfgs_add

    subroutine bfgs_clear(self)

        class(bfgs_updating_type), intent(inout) :: self

        if (allocated(self%p_list)) deallocate(self%p_list)
        if (allocated(self%y_list)) deallocate(self%y_list)
        if (allocated(self%d_p_list)) deallocate(self%d_p_list)
        if (allocated(self%d_y_list)) deallocate(self%d_y_list)
        self%n = 0

    end subroutine bfgs_clear

    subroutine bfgs_hess_x_fun(x, hess_x, error)

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: i
        real(rp) :: coeff_p, coeff_y
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! equation (7.33) in Numerical Optimization
        hess_x = bfgs_object%start * x
        do i = 1, bfgs_object%n
            coeff_p = ddot(size(x), bfgs_object%p_list(:, i), 1, x, 1) / &
                      bfgs_object%d_p_list(i)
            coeff_y = ddot(size(x), bfgs_object%y_list(:, i), 1, x, 1) / &
                      bfgs_object%d_y_list(i)
            hess_x = hess_x - coeff_p * bfgs_object%p_list(:, i) + coeff_y * &
                     bfgs_object%y_list(:, i)
        end do

    end subroutine bfgs_hess_x_fun

end module otr_qn
