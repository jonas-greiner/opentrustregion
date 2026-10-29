! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek

    use opentrustregion_core, only: rp, ip, settings_type, update_orbs_type, &
                                    hess_x_type, add_error_origin
    use otr_common, only: extension_type, provided_roles_type

    implicit none

    abstract interface
        subroutine change_reference_type(new_ref, n_points, kappa_list, &
                                         local_grad_list, grad_list, error, context)
            import :: rp, ip

            real(rp), intent(in), target :: new_ref(:)
            integer(ip), intent(in) :: n_points
            real(rp), intent(inout), target :: kappa_list(:, :), &
                                               local_grad_list(:, :), grad_list(:, :)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine change_reference_type
    end interface

    type, extends(settings_type) :: s_gek_settings_type
        logical :: use_subspace = .true.
        integer(ip) :: max_points = 20
    end type

    type(s_gek_settings_type), parameter :: default_s_gek_settings = &
        s_gek_settings_type()

    ! warning for settings that were not initialized
    character(len=*), parameter :: s_gek_settings_uninitialized_warning_msg = &
        "S-GEK settings were not initialized. All S-GEK settings are set to "// &
        "default values"

    ! error origin of the S-GEK extension, with which its factory and object tag the
    ! errors they detect themselves
    integer(ip), parameter :: error_s_gek = 5500

    ! S-GEK object, holding the sample points, the surrogate model, the settings and
    ! the host's change of reference function
    type, extends(extension_type) :: s_gek_type
        real(rp), allocatable :: kappa_list(:, :), disp_list(:, :), obj_func_list(:), &
                                 grad_list(:, :), local_grad_list(:, :), y_list(:, :), &
                                 h_diag(:), char_length(:), transform_matrix(:, :), &
                                 kriging_weights(:), hess(:, :), subspace(:, :), &
                                 actual_disp_list(:, :), actual_grad_list(:, :), &
                                 actual_hess_approx(:, :), actual_trans_kappa_list(:, :)
        integer(ip) :: n_points, n_covariance, n_space
        real(rp) :: loosen_factor
        type(s_gek_settings_type) :: settings
        procedure(change_reference_type), pointer, nopass :: change_reference => null()
    contains
        procedure :: update_orbs => update_orbs_s_gek
        procedure :: hess_x => hess_x_s_gek
        procedure :: init => init_s_gek
        procedure :: add => add_s_gek
        procedure :: setup_kriging
        procedure :: kriging_model
        procedure :: predict_hess
        procedure :: bfgs_inv_hess_x
    end type

    ! mathematical constants
    real(rp), parameter :: pi = 4.0_rp * atan(1.0_rp)

    ! bias
    real(rp), parameter :: bias = 10.0_rp

    ! controls nu parameter of Matern covariance function (1 => 3/2, 2 => 5/2, 3 => 7/2)
    integer(ip), parameter :: p_matern = 2

    ! undershoot avoidance parameters
    real(rp), parameter :: angle_lower_thres = cos(5.0_rp * pi / 180.0_rp), &
                           angle_upper_thres = cos(20.0_rp * pi / 180.0_rp), &
                           loosen_step = 0.5_rp * (1.0_rp + sqrt(5.0_rp))

contains

    subroutine s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                             change_reference_funptr, update_orbs_s_gek_funptr, &
                             obj_func_s_gek_funptr, solver_settings, solver_state, &
                             error, settings)
        !
        ! this subroutine returns the S-GEK orbital updating function, which wraps the
        ! host's, and the objective function, and wires the solver settings with a new
        ! S-GEK object, which they own and through which the host's callback functions
        ! keep receiving the host context, and which reads the state of the solver from
        ! the record the host also passes to the solver
        !
        use opentrustregion_core, only: obj_func_type, solver_state_type
        use otr_common, only: extension_solver_settings_type, update_orbs_callback, &
                              obj_func_callback

        integer(ip), intent(in) :: n_param
        procedure(update_orbs_type), intent(in), pointer :: update_orbs_funptr
        procedure(obj_func_type), intent(in), pointer :: obj_func_funptr
        procedure(change_reference_type), intent(in), pointer :: change_reference_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_s_gek_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_s_gek_funptr
        type(extension_solver_settings_type), intent(inout), target :: solver_settings
        type(solver_state_type), intent(inout), target :: solver_state
        integer(ip), intent(out) :: error
        type(s_gek_settings_type), intent(in), optional :: settings

        type(s_gek_settings_type) :: local_settings
        class(extension_type), allocatable, target :: object

        ! initialize error flag and the returned functions
        error = 0
        update_orbs_s_gek_funptr => null()
        obj_func_s_gek_funptr => null()

        ! work on a copy of the given settings or, without any, of the default settings
        if (present(settings)) then
            local_settings = settings
        else
            local_settings = default_s_gek_settings
        end if

        ! perform sanity check
        call s_gek_sanity_check(local_settings, n_param, error)
        call add_error_origin(error, error_s_gek, local_settings)
        if (error /= 0) return

        ! create the object and initialize it with the host's change of reference
        ! function
        allocate(s_gek_type :: object)
        select type (object)
        type is (s_gek_type)
            call object%init(n_param, local_settings)
            object%change_reference => change_reference_funptr
        end select

        ! wire the solver settings with the object, which reads the record of the state
        ! and saves the host's orbital update, objective function, callback functions
        ! and host contexts, and let the settings own it
        call solver_settings%wire(object, solver_state, local_settings, error, &
                                  update_orbs_funptr, obj_func_funptr)
        call add_error_origin(error, error_s_gek, local_settings)
        if (error /= 0) return

        ! get the orbital updating function and the objective function of the object
        update_orbs_s_gek_funptr => update_orbs_callback
        obj_func_s_gek_funptr => obj_func_callback

    end subroutine s_gek_factory

    subroutine s_gek_sanity_check(settings, n_param, error)
        !
        ! this subroutine performs a sanity check for S-GEK input parameters
        !
        use opentrustregion_core, only: verbosity_error

        type(s_gek_settings_type), intent(in) :: settings
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! check that the number of parameters is positive
        if (n_param < 1) then
            call settings%log("Number of parameters should be larger than 0.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check that at least the current point is stored
        if (settings%max_points < 1) then
            call settings%log("Maximum number of points should be larger than 0.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine s_gek_sanity_check

    subroutine invoke_change_reference(change_reference, new_ref, n_points, &
                                       kappa_list, local_grad_list, grad_list, &
                                       settings, error_origin, error)
        !
        ! this subroutine calls the host's change of reference function with the host
        ! context of the settings and passes on its error, and returns an error of its
        ! own with the origin of the object if the function returns values that are not
        ! finite for the stored points
        !
        use opentrustregion_core, only: verbosity_error
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        procedure(change_reference_type), intent(in), pointer :: change_reference
        real(rp), intent(in), target :: new_ref(:)
        integer(ip), intent(in) :: n_points
        real(rp), intent(inout), target :: kappa_list(:, :), local_grad_list(:, :), &
                                           grad_list(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(in) :: error_origin
        integer(ip), intent(out) :: error

        ! call the change of reference function
        call change_reference(new_ref, n_points, kappa_list, local_grad_list, &
                              grad_list, error, settings%context)
        if (error /= 0) return

        ! check the result
        if (.not. (all(ieee_is_finite(kappa_list(:, :n_points))) .and. &
                   all(ieee_is_finite(local_grad_list(:, :n_points))) .and. &
                   all(ieee_is_finite(grad_list(:, :n_points))))) then
            call settings%log("Change of reference function returned values that "// &
                              "are not finite.", verbosity_error, .true.)
            error = 1
            call add_error_origin(error, error_origin, settings)
        end if

    end subroutine invoke_change_reference

    recursive subroutine update_orbs_s_gek(self, kappa, func, grad, h_diag, &
                                           hess_x_funptr, error, context)
        !
        ! this subroutine is the S-GEK orbital updating function, which updates the
        ! orbitals with the host's orbital update and host context, replaces the
        ! Hessian linear transformation the host returns by the S-GEK one and adds the
        ! new point to the history
        !
        use otr_common, only: hess_x_callback

        class(s_gek_type), intent(inout) :: self
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! initialize error flag
        error = 0

        ! update orbitals
        call self%solver_level%host%update_orbs(kappa, func, grad, h_diag, &
                                                hess_x_funptr, error, &
                                                self%solver_level%host%context)
        if (error /= 0) return

        ! get pointer to the S-GEK Hessian linear transformation of the object
        hess_x_funptr => hess_x_callback

        ! add the new point, also the starting point the first call reaches without a
        ! step
        if (sum(abs(kappa)) > 0.0_rp .or. self%n_points == 0) then
            call self%add(kappa, func, grad, h_diag, error)
            if (error /= 0) return
        end if

    end subroutine update_orbs_s_gek

    subroutine init_s_gek(self, n_param, settings)
        !
        ! this subroutine initializes the S-GEK object with the given settings
        !
        class(s_gek_type), intent(inout) :: self
        integer(ip), intent(in) :: n_param
        type(s_gek_settings_type), intent(in) :: settings

        ! set the roles the object performs itself, the orbital update and the Hessian
        ! linear transformation, and the settings
        self%provides = provided_roles_type(update_orbs=.true., hess_x=.true.)
        self%error_origin = error_s_gek
        self%settings = settings

        ! initialize dimensions
        self%n_points = 0
        self%n_param = n_param
        self%n_covariance = 0

        ! initialize loosen factor for undershoot avoidance
        self%loosen_factor = 1.0_rp

        ! allocate arrays
        allocate(self%h_diag(self%n_param))
        allocate(self%kappa_list(self%n_param, 0))
        allocate(self%disp_list(self%n_param, 0))
        allocate(self%obj_func_list(0))
        allocate(self%grad_list(self%n_param, 0))
        allocate(self%local_grad_list(self%n_param, 0))
        allocate(self%y_list(self%n_param, 0))

    end subroutine init_s_gek

    subroutine add_s_gek(self, kappa, obj_func, grad, h_diag, error)
        !
        ! this subroutine adds a new data point to the S-GEK object and updates the
        ! surrogate model
        !
        class(s_gek_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), obj_func, grad(:), h_diag(:)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: kappa_list(:, :), obj_func_list(:), &
                                 local_grad_list(:, :), scale_matrix(:, :), tmp(:, :), &
                                 bfgs_disp(:), actual_bfgs_disp(:)
        real(rp) :: angle, norm
        integer(ip) :: i, j, k, l
        real(rp), external :: ddot, dnrm2
        external :: dgemm, dgemv

        ! initialize error flag
        error = 0

        ! append new values
        if (self%n_points < self%settings%max_points) then
            allocate(kappa_list(self%n_param, self%n_points + 1))
            if (self%n_points > 0) &
                kappa_list(:, 1:self%n_points) = self%kappa_list(:, 1:self%n_points)
            kappa_list(:, self%n_points + 1) = kappa
            call move_alloc(kappa_list, self%kappa_list)

            deallocate(self%disp_list)
            allocate(self%disp_list(self%n_param, self%n_points))

            allocate(obj_func_list(self%n_points + 1))
            if (self%n_points > 0) &
                obj_func_list(1:self%n_points) = self%obj_func_list(1:self%n_points)
            obj_func_list(self%n_points + 1) = obj_func
            call move_alloc(obj_func_list, self%obj_func_list)

            deallocate(self%grad_list)
            allocate(self%grad_list(self%n_param, self%n_points + 1))

            allocate(local_grad_list(self%n_param, self%n_points + 1))
            if (self%n_points > 0) local_grad_list(:, 1:self%n_points) = &
                self%local_grad_list(:, 1:self%n_points)
            local_grad_list(:, self%n_points + 1) = grad
            call move_alloc(local_grad_list, self%local_grad_list)

            self%n_points = self%n_points + 1

        ! shift values
        else
            self%kappa_list(:, :self%n_points - 1) = self%kappa_list(:, 2:self%n_points)
            self%kappa_list(:, self%n_points) = kappa

            self%obj_func_list(:self%n_points - 1) = self%obj_func_list(2:self%n_points)
            self%obj_func_list(self%n_points) = obj_func

            self%local_grad_list(:, :self%n_points - 1) = &
                self%local_grad_list(:, 2:self%n_points)
            self%local_grad_list(:, self%n_points) = grad

        end if

        ! move displacement and gradient to current reference
        call invoke_change_reference(self%change_reference, kappa, self%n_points, &
                                     self%kappa_list, self%local_grad_list, &
                                     self%grad_list, self%settings, self%error_origin, &
                                     error)
        if (error /= 0) return

        ! set Hessian diagonal
        self%h_diag = h_diag

        ! calculate displacements
        do i = 1, self%n_points - 1
            self%disp_list(:, i) = self%kappa_list(:, i + 1) - self%kappa_list(:, i)
        end do

        ! get BFGS guess by reducing the update depth until the step is reasonable
        allocate(bfgs_disp(self%n_param))
        do i = 1, self%n_points
            call self%bfgs_inv_hess_x(self%grad_list(:, self%n_points), bfgs_disp, i, &
                                      error)
            if (error /= 0) return
            if (dnrm2(self%n_param, bfgs_disp, 1_ip) <= pi) exit
        end do

        ! check for undershoot avoidance
        if (self%n_points > 2) then
            ! calculate angle between last two displacements
            angle = ddot(self%n_param, self%disp_list(:, self%n_points - 2), 1_ip, &
                         self%disp_list(:, self%n_points - 1), 1_ip) / &
                    (dnrm2(self%n_param, self%disp_list(:, self%n_points - 1), 1_ip) * &
                     dnrm2(self%n_param, self%disp_list(:, self%n_points - 2), 1_ip))
            if (angle < angle_upper_thres) then
                self%loosen_factor = 1.0_rp
            else if (angle > angle_lower_thres) then
                self%loosen_factor = loosen_step * self%loosen_factor
            end if
        end if

        ! transform to subspace
        if (self%settings%use_subspace) then
            ! all (n-1) displacements, n gradients (or n-1 gradient differences and nth
            ! gradient) and BFGS direction
            self%n_space = 2 * (self%n_points - 1) + 2
            if (allocated(self%subspace)) deallocate(self%subspace)
            allocate(self%subspace(self%n_param, self%n_space))

            ! add displacements and gradient differences to the subspace and normalize
            ! these
            j = 1
            do k = 1, self%n_points - 1
                self%subspace(:, j) = self%grad_list(:, k + 1) - self%grad_list(:, k)
                self%subspace(:, j) = self%subspace(:, j) / &
                                      dnrm2(self%n_param, self%subspace(:, j), 1_ip)
                j = j + 1

                self%subspace(:, j) = self%disp_list(:, k) / &
                                      dnrm2(self%n_param, self%disp_list(:, k), 1_ip)
                j = j + 1
            end do

            ! add gradient to the subspace
            self%subspace(:, self%n_space - 1) = &
                self%grad_list(:, self%n_points) / &
                dnrm2(self%n_param, self%grad_list(:, self%n_points), 1_ip)

            ! add BFGS displacement to the subspace
            self%subspace(:, self%n_space) = bfgs_disp / &
                                             dnrm2(self%n_param, bfgs_disp, 1_ip)

            ! orthogonalize all unit vectors
            do l = 1, 2
                j = 1
                do i = 2, self%n_space
                    do k = 1, j
                        self%subspace(:, i) = self%subspace(:, i) - ddot( &
                            self%n_param, self%subspace(:, i), 1_ip, &
                            self%subspace(:, k), 1_ip) * self%subspace(:, k)
                    end do
                    ! renormalize and only add vector if it is not linearly dependent
                    norm = dnrm2(self%n_param, self%subspace(:, i), 1_ip)
                    if (norm**2 > 1e-17_rp) then
                        j = j + 1
                        self%subspace(:, j) = self%subspace(:, i) / norm
                    end if
                end do
            end do

            ! subspace can be smaller after orthogonalization whenever subspace vectors
            ! are linearly dependent
            self%n_space = j
            allocate(tmp(self%n_param, self%n_space))
            tmp = self%subspace(:, 1:self%n_space)
            call move_alloc(tmp, self%subspace)

        else
            ! consider full space
            self%n_space = self%n_param
        end if

        ! set size of covariance matrix
        self%n_covariance = self%n_points + self%n_space * self%n_points

        ! deallocate projected arrays
        if (allocated(self%actual_disp_list)) deallocate(self%actual_disp_list)
        if (allocated(self%actual_grad_list)) deallocate(self%actual_grad_list)
        if (allocated(self%actual_hess_approx)) deallocate(self%actual_hess_approx)
        if (allocated(self%actual_trans_kappa_list)) &
            deallocate(self%actual_trans_kappa_list)

        ! allocate projected arrays
        allocate(self%actual_disp_list(self%n_space, self%n_points), &
                 self%actual_grad_list(self%n_space, self%n_points), &
                 self%actual_hess_approx(self%n_space, self%n_space), &
                 self%actual_trans_kappa_list(self%n_space, self%n_points))

        ! transform to subspace
        if (self%settings%use_subspace) then
            ! project the displacement coordinates, gradients, and the approximate
            ! Hessian into the subspace
            do i = 1, self%n_points
                call dgemv("T", self%n_param, self%n_space, 1.0_rp, self%subspace, &
                           self%n_param, &
                           self%kappa_list(:, i) - self%kappa_list(:, self%n_points), &
                           1_ip, 0.0_rp, self%actual_disp_list(:, i), 1_ip)
                call dgemv("T", self%n_param, self%n_space, 1.0_rp, self%subspace, &
                           self%n_param, self%grad_list(:, i), 1_ip, 0.0_rp, &
                           self%actual_grad_list(:, i), 1_ip)
            end do
            allocate(tmp(self%n_param, self%n_space))
            do j = 1, self%n_space
                tmp(:, j) = self%h_diag * self%subspace(:, j)
            end do
            call dgemm("T", "N", self%n_space, self%n_space, self%n_param, 1.0_rp, &
                       tmp, self%n_param, self%subspace, self%n_param, 0.0_rp, &
                       self%actual_hess_approx, self%n_space)
            deallocate(tmp)

            ! construct BFGS displacement into the subspace
            if (self%loosen_factor > 1.0_rp) then
                allocate(actual_bfgs_disp(self%n_space))
                call dgemv("T", self%n_param, self%n_space, 1.0_rp, self%subspace, &
                           self%n_param, bfgs_disp, 1_ip, 0.0_rp, actual_bfgs_disp, &
                           1_ip)
            end if

        else
            ! set actual displacements, gradients, and Hessian approximation
            do i = 1, self%n_points
                self%actual_disp_list(:, i) = self%kappa_list(:, i) - &
                                              self%kappa_list(:, self%n_points)
            end do
            self%actual_grad_list = self%grad_list
            self%actual_hess_approx = 0.0_rp
            do i = 1, self%n_space
                self%actual_hess_approx(i, i) = self%h_diag(i)
            end do

            ! set BFGS displacement
            if (self%loosen_factor > 1.0_rp) then
                actual_bfgs_disp = bfgs_disp
            end if

        end if
        deallocate(bfgs_disp)

        ! if undershoot avoidance triggers scale along BFGS displacement
        if (self%loosen_factor > 1.0_rp) then
            ! construct scale matrix
            allocate(scale_matrix(self%n_space, self%n_space), &
                     tmp(self%n_space, self%n_space))
            do i = 1, self%n_space
                scale_matrix(:, i) = &
                    (1.0_rp / self%loosen_factor - 1.0_rp) * actual_bfgs_disp(i) * &
                    actual_bfgs_disp / &
                    ddot(self%n_space, actual_bfgs_disp, 1_ip, actual_bfgs_disp, 1_ip)
                scale_matrix(i, i) = scale_matrix(i, i) + 1.0_rp
            end do
            deallocate(actual_bfgs_disp)

            ! scale Hessian approximation
            call dgemm("N", "N", self%n_space, self%n_space, self%n_space, 1.0_rp, &
                       self%actual_hess_approx, self%n_space, scale_matrix, &
                       self%n_space, 0.0_rp, tmp, self%n_space)
            call dgemm("N", "N", self%n_space, self%n_space, self%n_space, 1.0_rp, &
                       scale_matrix, self%n_space, tmp, self%n_space, 0.0_rp, &
                       self%actual_hess_approx, self%n_space)
            deallocate(scale_matrix, tmp)
        end if

        ! allocate kriging weights
        if (allocated(self%kriging_weights)) deallocate(self%kriging_weights)
        allocate(self%kriging_weights(self%n_covariance))

        ! setup new kriging model
        call self%setup_kriging(self%actual_disp_list, self%actual_grad_list, &
                                self%obj_func_list, self%actual_hess_approx, error)
        if (error /= 0) return

        ! predict Hessian
        if (allocated(self%hess)) deallocate(self%hess)
        allocate(self%hess(self%n_space, self%n_space))
        self%hess = self%predict_hess(self%actual_disp_list(:, self%n_points))

    end subroutine add_s_gek

    subroutine hess_x_s_gek(self, x, hess_x, error, context)
        !
        ! this subroutine is the S-GEK Hessian linear transformation of the object
        !
        class(s_gek_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: i
        real(rp), allocatable :: actual_x(:), actual_hess_x(:)
        real(rp), external :: ddot
        external :: dgemv

        ! initialize error flag
        error = 0

        ! compute the trial vector in the subspace
        allocate(actual_x(self%n_space), actual_hess_x(self%n_space))
        if (self%settings%use_subspace) then
            do i = 1, self%n_space
                actual_x(i) = ddot(size(x, kind=ip), x, 1_ip, self%subspace(:, i), 1_ip)
            end do
        else
            actual_x = x
        end if

        ! compute the Hessian vector product in the subspace
        call dgemv("N", self%n_space, self%n_space, 1.0_rp, self%hess, self%n_space, &
                   actual_x, 1_ip, 0.0_rp, actual_hess_x, 1_ip)

        ! compute the Hessian linear transformation in the full space
        if (self%settings%use_subspace) then
            hess_x = 0.0_rp
            do i = 1, self%n_space
                hess_x = hess_x + actual_hess_x(i) * self%subspace(:, i)
            end do
        else
            hess_x = actual_hess_x
        end if

        deallocate(actual_x, actual_hess_x)

    end subroutine hess_x_s_gek

    subroutine setup_kriging(self, kappa_list, grad_list, obj_func_list, &
                             approx_hessian, error)
        !
        ! this subroutine sets up the kriging model for the S-GEK method
        !
        use opentrustregion_core, only: symm_mat_diag

        class(s_gek_type), intent(inout) :: self
        real(rp), intent(in) :: kappa_list(:, :), grad_list(:, :), obj_func_list(:)
        real(rp), intent(inout) :: approx_hessian(:, :)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: trans_grad_list(:, :), hess_eigvals(:)
        external :: dgemm

        ! initialize error flag
        error = 0

        ! diagonalize approximate Hessian and set transformation matrix to transform to
        ! approximate Hessian eigenvector basis
        if (allocated(self%transform_matrix)) deallocate(self%transform_matrix)
        allocate(hess_eigvals(self%n_space), &
                 self%transform_matrix(self%n_space, self%n_space))
        call symm_mat_diag(approx_hessian, hess_eigvals, self%transform_matrix, &
                           self%settings, error)
        call add_error_origin(error, self%error_origin, self%settings)
        if (error /= 0) return

        ! set the characteristic length such that the kriging hessian reproduces the
        ! diagonal value of the approximate Hessian
        if (allocated(self%char_length)) deallocate(self%char_length)
        allocate(self%char_length(self%n_space))
        call set_char_length(self%char_length, bias, hess_eigvals, self%n_space)
        deallocate(hess_eigvals)

        ! allocate arrays for transformed gradients
        allocate(trans_grad_list(self%n_space, self%n_points))

        ! transform to the basis which diagonalizes the approximate Hessian and set the
        ! coordinates of the sample points and gradients of the function at the sample
        ! points
        call dgemm("T", "N", self%n_space, self%n_points, self%n_space, 1.0_rp, &
                   self%transform_matrix, self%n_space, kappa_list, self%n_space, &
                   0.0_rp, self%actual_trans_kappa_list, self%n_space)
        call dgemm("T", "N", self%n_space, self%n_points, self%n_space, 1.0_rp, &
                   self%transform_matrix, self%n_space, grad_list, self%n_space, &
                   0.0_rp, trans_grad_list, self%n_space)

        ! form the inverse of the covariance matrix times the generalized value vector
        ! to get kriging weights
        call self%kriging_model(obj_func_list, trans_grad_list, error)
        if (error /= 0) return

        deallocate(trans_grad_list)

    end subroutine setup_kriging

    subroutine set_char_length(char_length, baseline, hess_eigvals, n_space)
        !
        ! this subroutine sets the characteristic length for each parameter such that
        ! the kriging Hessian for a single point of Kriging is identical to the
        ! diagonal values of the approximate Hessian
        !
        real(rp), intent(out) :: char_length(:)
        real(rp), intent(in) :: baseline, hess_eigvals(:)
        integer(ip), intent(in) :: n_space

        integer(ip) :: i
        real(rp) :: hess_bounded
        real(rp), parameter :: hess_min = 0.025_rp

        do i = 1, n_space
            ! make sure that the characteristic length is not too long
            hess_bounded = max(abs(hess_eigvals(i)), hess_min)
            char_length(i) = sqrt(5.0_rp / 3.0_rp * baseline / hess_bounded)
        end do

    end subroutine set_char_length

    subroutine kriging_model(self, obj_func_list, trans_grad_list, error)
        !
        ! this function generates the covariance matrix and then constructs kriging
        ! weights by solving x = R^{-1} (y - mu f) where R is the covariance matrix, y
        ! is a vector including values and gradients and mu is the trend function
        !
        use opentrustregion_core, only: verbosity_error

        class(s_gek_type), intent(inout) :: self
        real(rp), intent(in) :: obj_func_list(:), trans_grad_list(:, :)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: covariance_matrix(:, :)
        integer(ip) :: is, ie, i, info
        real(rp) :: mu
        character(len=300) :: msg
        external :: dposv

        ! initialize error flag
        error = 0

        ! generate the covariance matrix
        allocate(covariance_matrix(self%n_covariance, self%n_covariance))
        covariance_matrix = get_covariance_matrix( &
            self%char_length, self%actual_trans_kappa_list, self%n_points, &
            self%n_space, self%n_covariance)
        ! establish the trend function (baseline) mu, make sure the base line is above
        ! any data point, this to make sure the surrogate model is bound
        mu = -huge(mu)
        do i = 1, self%n_points
            mu = max(mu, obj_func_list(i) + bias)
        end do

        ! add the biased value vector
        self%kriging_weights(:self%n_points) = obj_func_list - mu

        ! add the gradients
        do i = 1, self%n_space
            is = i * self%n_points + 1
            ie = is + self%n_points - 1
            self%kriging_weights(is:ie) = trans_grad_list(i, :)
        end do

        ! solve R x = (y - mu f), i.e. x = R^{-1} (y - mu f)
        call dposv("U", self%n_covariance, 1_ip, covariance_matrix, self%n_covariance, &
                   self%kriging_weights, self%n_covariance, info)

        ! check for successful execution
        if (info /= 0) then
            write(msg, '(A, I0)') "Linear solve failed: Error in DPOSV, info = ", info
            call self%settings%log(msg, verbosity_error, .true.)
            error = 1
            call add_error_origin(error, self%error_origin, self%settings)
            return
        end if

        ! deallocate covariance matrix
        deallocate(covariance_matrix)

    end subroutine kriging_model

    function get_covariance_matrix(char_length, kappa_list, n_points, n_param, &
                                   n_covariance) result(covariance_matrix)
        !
        ! this function computes the covariance function for all distances between the
        ! sample points, including gradients and second derivatives of the covariance
        ! function with respect to all distances for all distances compare the
        ! resulting covariance matrix with eq 2 of doi:10.1007/s00366-015-0397
        !
        real(rp), intent(in) :: char_length(:), kappa_list(:, :)
        integer(ip), intent(in) :: n_points, n_param, n_covariance
        real(rp) :: covariance_matrix(n_covariance, n_covariance)

        integer(ip) :: i, j, i0, i1, j0, j1, k
        real(rp), allocatable :: diffx_j(:, :), diffx_i(:, :), matern_1st_deriv(:, :), &
                                 matern_2nd_deriv(:, :), distances_per_param(:, :, :), &
                                 distances(:, :)
        ! shifts to avoid becoming singular in first and second derivatives
        real(rp), parameter :: eps = 1e-13_rp, eps2 = 1e-10_rp

        ! allocate temporary memory
        allocate(diffx_j(n_points, n_points), diffx_i(n_points, n_points), &
                 matern_1st_deriv(n_points, n_points), &
                 matern_2nd_deriv(n_points, n_points), &
                 distances_per_param(n_points, n_points, n_param), &
                 distances(n_points, n_points))

        ! compute the distances between the sample points using the characteristic
        ! length
        covariance_matrix = 0.0_rp
        distances = 0.0_rp
        do i = 1, n_param
            do k = 1, n_points
                do j = 1, n_points
                    distances_per_param(k, j, i) = &
                        (kappa_list(i, k) - kappa_list(i, j)) / char_length(i)
                end do
            end do
            ! accumulate contributions to the square of the individual distances
            distances = distances + distances_per_param(:, :, i)**2
        end do
        distances = sqrt(distances)

        ! evaluate the covariance function for all the distances
        call get_matern(distances, covariance_matrix(:n_points, :n_points), n_points, &
                        n_points)

        ! evaluate first derivatives of the covariance function with respect to d at
        ! all distances, note that for the full derivative this has to be complemented
        ! by the derivative of d with respect to the individual components of the
        ! coordinates
        call get_matern_deriv(1_ip, distances, matern_1st_deriv, n_points, n_points)

        ! loop over component of the coordinate to differentiate
        do i = 1, n_param
            ! compute the range of the block in the covariance matrix
            i0 = i * n_points + 1
            i1 = i0 + n_points - 1

            ! do an on-the-fly evaluation of the derivative with respect to x_i
            diffx_i = -2.0_rp * distances_per_param(:, :, i) / char_length(i)

            ! writing the 1st row of 1st derivatives with respect the coordinates
            covariance_matrix(:n_points, i0:i1) = matern_1st_deriv * diffx_i
        end do

        ! fill in the transpose block
        covariance_matrix(n_points + 1:n_covariance, :n_points) = &
            transpose(covariance_matrix(:n_points, n_points + 1:n_covariance))

        ! evaluate the second derivatives of the covariance function with respect to d
        ! at all distances
        call get_matern_deriv(2_ip, distances, matern_2nd_deriv, n_points, n_points)

        ! loop over component of the first coordinate to differentiate
        do i = 1, n_param
            i0 = i * n_points + 1
            i1 = i0 + n_points - 1

            diffx_i = -2.0_rp * distances_per_param(:, :, i) / char_length(i)

            ! loop over component of the second coordinate to differentiate
            do j = i, n_param
                j0 = j * n_points + 1
                j1 = j0 + n_points - 1

                diffx_j = 2.0_rp * distances_per_param(:, :, j) / char_length(j)

                ! if differentiating twice on the same dimension
                covariance_matrix(i0:i1, j0:j1) = matern_2nd_deriv * diffx_j * diffx_i
                if (i == j) covariance_matrix(i0:i1, j0:j1) = &
                    covariance_matrix(i0:i1, j0:j1) - matern_1st_deriv * &
                    (2.0_rp / (char_length(i) * char_length(j)))

                ! fill in the transpose block
                if (i /= j) covariance_matrix(j0:j1, i0:i1) = &
                    transpose(covariance_matrix(i0:i1, j0:j1))
            end do

        end do

        ! add constants to reflect the error in the cost function and the gradient,
        ! respectively
        do j = 1, n_covariance
            if (j <= n_points) then
                covariance_matrix(j, j) = covariance_matrix(j, j) + eps
            else
                covariance_matrix(j, j) = covariance_matrix(j, j) + eps2
            end if
        end do

        deallocate(diffx_j, diffx_i, matern_1st_deriv, matern_2nd_deriv, &
                   distances_per_param, distances)

    end function get_covariance_matrix

    function predict_hess(self, kappa) result(hess)
        !
        ! this function predicts the Hessian at a given coordinate using the S-GEK
        ! model
        !
        class(s_gek_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:)
        real(rp) :: hess(self%n_space, self%n_space)

        real(rp), allocatable :: kappa_s(:), hess_s(:, :), covariance_vector(:, :, :)
        integer(ip) :: i, k
        real(rp), external :: ddot
        external :: dgemm

        allocate(kappa_s(self%n_space), hess_s(self%n_space, self%n_space), &
                 covariance_vector(self%n_covariance, self%n_space, self%n_space))

        ! transform coordinates to eigenvectors of approximate Hessian
        call dgemm("T", "N", self%n_space, 1_ip, self%n_space, 1.0_rp, &
                   self%transform_matrix, self%n_space, kappa, self%n_space, 0.0_rp, &
                   kappa_s, self%n_space)

        ! generate covariance vector for Hessian at current point which contains the
        ! correlation function values and gradients between the sample data and the
        ! prediction
        covariance_vector = covariance_vector_hess( &
            kappa_s, self%char_length, self%actual_trans_kappa_list, self%n_points, &
            self%n_space, self%n_covariance)

        ! predict Hessian for surrogate model
        do k = 1, self%n_space
            do i = k, self%n_space
                hess(k, i) = ddot(self%n_covariance, covariance_vector(:, i, k), 1_ip, &
                                  self%kriging_weights, 1_ip)
                if (i /= k) hess(i, k) = hess(k, i)
            end do
        end do

        ! backtransform Hessian to original Hessian
        call dgemm("N", "N", self%n_space, self%n_space, self%n_space, 1.0_rp, &
                   self%transform_matrix, self%n_space, hess, self%n_space, 0.0_rp, &
                   hess_s, self%n_space)
        call dgemm("N", "T", self%n_space, self%n_space, self%n_space, 1.0_rp, hess_s, &
                   self%n_space, self%transform_matrix, self%n_space, 0.0_rp, hess, &
                   self%n_space)

        deallocate(kappa_s, hess_s, covariance_vector)

    end function predict_hess

    function covariance_vector_hess(curr_kappa, char_length, kappa_list, n_points, &
                                    n_param, n_covariance) result(covariance_vector)
        !
        ! this function calculates the covariance vector for the Hessian prediction
        !
        real(rp), intent(in) :: curr_kappa(:), char_length(:), kappa_list(:, :)
        integer(ip), intent(in) :: n_points, n_param, n_covariance
        real(rp) :: covariance_vector(n_covariance, n_param, n_param)

        integer(ip) :: i, j, k, k0, k1
        real(rp) :: sdiffxi, sdiffxj, sdiffxk
        real(rp), allocatable :: distances_per_param(:, :), distances(:), &
                                 matern_1st_deriv(:), matern_2nd_deriv(:), &
                                 matern_3rd_deriv(:), diffxi(:), diffxj(:), diffxk(:)

        ! allocate temporary arrays
        allocate(distances_per_param(n_points, n_param), distances(n_points), &
                 matern_1st_deriv(n_points), matern_2nd_deriv(n_points), &
                 matern_3rd_deriv(n_points), diffxi(n_points), diffxj(n_points), &
                 diffxk(n_points))

        ! calculate distances_per_param and distances which are temporary matrices for
        ! the construction of Psi which is inside of Grad-Psi
        distances = 0.0_rp
        do i = 1, n_param
            do j = 1, n_points
                distances_per_param(j, i) = (kappa_list(i, j) - curr_kappa(i)) / &
                                            char_length(i)
            end do
            distances = distances + distances_per_param(:, i)**2
        end do
        distances = sqrt(distances)

        ! calculate derivatives of covariance matrix
        call get_matern_deriv(1_ip, distances, matern_1st_deriv, n_points, 1_ip)
        call get_matern_deriv(2_ip, distances, matern_2nd_deriv, n_points, 1_ip)
        call get_matern_deriv(3_ip, distances, matern_3rd_deriv, n_points, 1_ip)

        do i = 1, n_param
            diffxi = 2.0_rp * distances_per_param(:, i) / char_length(i)
            sdiffxi = 2.0_rp / char_length(i)**2
            do j = 1, n_param
                diffxj = 2.0_rp * distances_per_param(:, j) / char_length(j)
                sdiffxj = 2.0_rp / char_length(j)**2

                ! get value part
                if (i == j) then
                    covariance_vector(:n_points, i, j) = &
                        matern_2nd_deriv * diffxi * diffxj + &
                        matern_1st_deriv * 2.0_rp / (char_length(i) * char_length(j))
                else
                    covariance_vector(:n_points, i, j) = matern_2nd_deriv * diffxi * &
                                                         diffxj
                end if

                ! get gradient part
                do k = 1, n_param
                    diffxk(:) = 2.0_rp * distances_per_param(:, k) / char_length(k)
                    sdiffxk = 2.0_rp / char_length(k)**2
                    k0 = k * n_points + 1
                    k1 = k0 + n_points - 1
                    if (i == j .and. j == k) then
                        covariance_vector(k0:k1, i, j) = &
                            matern_3rd_deriv * diffxi * diffxj * diffxk + &
                            3.0_rp * matern_2nd_deriv * diffxi * sdiffxj
                    else if (i == j) then
                        covariance_vector(k0:k1, i, j) = &
                            matern_3rd_deriv * diffxi * diffxj * diffxk + &
                            matern_2nd_deriv * diffxk * sdiffxi
                    else if (i == k) then
                        covariance_vector(k0:k1, i, j) = &
                            matern_3rd_deriv * diffxi * diffxj * diffxk + &
                            matern_2nd_deriv * diffxj * sdiffxi
                    else if (j == k) then
                        covariance_vector(k0:k1, i, j) = &
                            matern_3rd_deriv * diffxi * diffxj * diffxk + &
                            matern_2nd_deriv * diffxi * sdiffxk
                    else
                        covariance_vector(k0:k1, i, j) = matern_3rd_deriv * diffxi * &
                                                         diffxj * diffxk
                    end if
                end do
            end do
        end do

        ! deallocate temporary arrays
        deallocate(distances_per_param, distances, matern_1st_deriv, matern_2nd_deriv, &
                   matern_3rd_deriv, diffxi, diffxj, diffxk)

    end function covariance_vector_hess

    subroutine get_matern(distance, matern, dim_1, dim_2)
        !
        ! this subroutine computes the Matern covariance function
        !
        integer(ip), intent(in) :: dim_1, dim_2
        real(rp), intent(in) :: distance(dim_1, dim_2)
        real(rp), intent(out) :: matern(dim_1, dim_2)
        integer(ip) :: i
        real(rp) :: prefactor, real_i

        select case (p_matern)
        ! v = 1/2
        case (0)
            matern = exp(-distance)
        ! v = 3/2
        case (1)
            matern = exp(-sqrt(3.0_rp) * distance) * (sqrt(3.0_rp) * distance + 1.0_rp)
        ! v = 5/2
        case (2)
            matern = exp(-sqrt(5.0_rp) * distance) * &
                     (5.0_rp / 3.0_rp * distance**2 + sqrt(5.0_rp) * distance + 1.0_rp)
        ! v = 7/2
        case (3)
            matern = exp(-sqrt(7.0_rp) * distance) * &
                     (7.0_rp / 15.0_rp * sqrt(7.0_rp) * distance**3 + &
                      14.0_rp / 5.0_rp * distance**2 + sqrt(7.0_rp) * distance + 1.0_rp)
        ! general expression
        ! https://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
        case default
            prefactor = gamma(p_matern + 1.0_rp) / gamma(2.0_rp * p_matern + 1.0_rp)
            matern = 0.0_rp
            do i = 0, p_matern
                real_i = real(i, kind=rp)
                matern = matern + ( &
                    gamma(p_matern + 1.0_rp + real_i) / &
                    (gamma(real_i + 1.0_rp) * gamma(p_matern + 1.0_rp - real_i))) * &
                         (2.0_rp * sqrt(2.0_rp * p_matern + 1.0_rp) * distance)** &
                         (p_matern - i)
            end do
            matern = prefactor * matern * &
                     exp(-sqrt(2.0_rp * p_matern + 1.0_rp) * distance)
        end select

    end subroutine get_matern

    subroutine get_matern_deriv(diff_order, distance, matern_deriv, dim_1, dim_2)
        !
        ! this subroutine computes the derivatives of the Matern covariance function
        !
        integer(ip), intent(in) :: diff_order, dim_1, dim_2
        real(rp), intent(in) :: distance(dim_1, dim_2)
        real(rp), intent(out) :: matern_deriv(dim_1, dim_2)

        integer(ip) :: i, j
        real(rp) :: matern_const
        real(rp), allocatable :: matern_weighted(:, :)
        real(rp), parameter :: Thr = epsilon(Thr)

        allocate(matern_weighted(dim_1, dim_2))

        matern_const = sqrt(2.0_rp * p_matern + 1.0_rp)
        matern_weighted = (2.0_rp * p_matern + 1.0_rp) / &
                          (2.0_rp * p_matern - 1.0_rp) * exp(-matern_const * distance)
        select case (p_matern)
        case (1)
            select case (diff_order)
            case (1)
                matern_deriv = -matern_weighted / 2.0_rp
            case (2)
                do j = 1, dim_2
                    do i = 1, dim_1
                        if (abs(distance(i, j)) > Thr) then
                            matern_deriv(i, j) = matern_const / distance(i, j) * &
                                                 matern_weighted(i, j) / 4.0_rp
                        else
                            matern_deriv(i, j) = 0.0_rp
                        end if
                    end do
                end do
            case (3)
                matern_deriv = -(2.0_rp * matern_const - 3.0_rp * distance) * &
                               matern_weighted
            end select
        case (2)
            select case (diff_order)
            case (1)
                matern_deriv = -matern_weighted * (1.0_rp + matern_const * distance) / &
                               2.0_rp
            case (2)
                matern_deriv = matern_weighted * 5.0_rp / 4.0_rp
            case (3)
                do j = 1, dim_2
                    do i = 1, dim_1
                        if (abs(distance(i, j)) > Thr) then
                            matern_deriv(i, j) = -5.0_rp / 8.0_rp * matern_const / &
                                                 distance(i, j) * matern_weighted(i, j)
                        else
                            matern_deriv(i, j) = 0.0_rp
                        end if
                    end do
                end do
            end select
        case (3)
            select case (diff_order)
            case (1)
                matern_deriv = -matern_weighted * &
                               (1.0_rp + matern_const * distance + distance**2) / 2.0_rp
            case (2)
                matern_deriv = matern_weighted * 7.0_rp * &
                               (1.0_rp + matern_const * distance) / 12.0_rp
            case (3)
                matern_deriv = -matern_weighted * 49.0_rp / 24.0_rp
            end select
        end select

        deallocate(matern_weighted)

    end subroutine get_matern_deriv

    subroutine bfgs_inv_hess_x(self, x, inv_hess_x, start_iter, error)
        !
        ! this implements the two-loop recursion for the BFGS update of the inverse
        ! Hessian of the object applied to a vector x according to T. H. Fischer and J.
        ! Almloef, JPC 96, 9768 (1992), doi:10.1021/j100203a036
        !
        use opentrustregion_core, only: numerical_zero

        class(s_gek_type), intent(inout) :: self
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: inv_hess_x(:)
        integer(ip), intent(in) :: start_iter
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, i, it
        real(rp) :: s_mat(6), t_mat(4)
        logical :: update_y
        real(rp), allocatable :: list(:, :), kappa_diff(:), grad_diff(:), &
                                 curr_kappa_diff(:), curr_grad_diff(:)
        real(rp), parameter :: vanish_denom_thr = 1e-9_rp
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(x)

        ! check if trial vector vanishes
        if (abs(ddot(n_param, x, 1_ip, x, 1_ip)) < numerical_zero) then
            inv_hess_x = 0.0_rp
            return
        end if

        ! initialize dot products
        s_mat = 0.0_rp
        t_mat = 0.0_rp

        ! initialize with Hessian diagonal
        do i = 1, n_param
            if (abs(self%h_diag(i)) < vanish_denom_thr) then
                inv_hess_x(i) = 1e2_rp * x(i)
            else
                inv_hess_x(i) = x(i) / self%h_diag(i)
            end if
        end do

        ! only current point is available
        if (self%n_points == start_iter) return

        ! get current displacement and gradient difference
        curr_kappa_diff = self%kappa_list(:, self%n_points) - &
                          self%kappa_list(:, self%n_points - 1)
        curr_grad_diff = self%grad_list(:, self%n_points) - &
                         self%grad_list(:, self%n_points - 1)

        ! initialize initial Hessian approximation multiplied with gradient
        update_y = .false.
        if (size(self%y_list, 2) < self%n_points - 1) then
            update_y = .true.
            allocate(list(size(x), self%n_points - 1))
            if (self%n_points > 2) list(:, :self%n_points - 2) = self%y_list
            do i = 1, n_param
                if (abs(self%h_diag(i)) < vanish_denom_thr) then
                    list(i, self%n_points - 1) = 1e2_rp * curr_grad_diff(i)
                else
                    list(i, self%n_points - 1) = curr_grad_diff(i) / self%h_diag(i)
                end if
            end do
            call move_alloc(list, self%y_list)
        end if

        ! loop until (n-2)th iteraion
        do it = start_iter, self%n_points - 2
            ! get displacement and gradient difference for iteration
            kappa_diff = self%kappa_list(:, it + 1) - self%kappa_list(:, it)
            grad_diff = self%grad_list(:, it + 1) - self%grad_list(:, it)

            ! calculate dot products (s_mat(2) is the inverse of the paper)
            s_mat(1) = ddot(n_param, kappa_diff, 1_ip, grad_diff, 1_ip)
            if (abs(s_mat(1)) < vanish_denom_thr) then
                s_mat(1) = 0.0_rp
            else
                s_mat(1) = 1.0_rp / s_mat(1)
            end if
            s_mat(2) = ddot(n_param, grad_diff, 1_ip, self%y_list(:, it), 1_ip)
            if (abs(s_mat(2)) < vanish_denom_thr) then
                s_mat(2) = 1.0_rp
            else
                s_mat(2) = 1.0_rp + s_mat(1) * s_mat(2)
            end if
            s_mat(3) = ddot(n_param, kappa_diff, 1_ip, x, 1_ip)
            s_mat(4) = ddot(n_param, self%y_list(:, it), 1_ip, x, 1_ip)

            ! get dot products for current Hessian approximation multiplied with
            ! gradient
            if (update_y) then
                s_mat(5) = ddot(n_param, kappa_diff, 1_ip, curr_grad_diff, 1_ip)
                s_mat(6) = ddot(n_param, self%y_list(:, it), 1_ip, curr_grad_diff, 1_ip)
            end if

            ! calculate more dot products
            t_mat(1) = s_mat(2) * t_mat(2) - s_mat(1) * s_mat(4)
            t_mat(2) = s_mat(1) * s_mat(3)
            if (update_y) then
                t_mat(3) = s_mat(2) * t_mat(4) - s_mat(1) * s_mat(6)
                t_mat(4) = s_mat(1) * s_mat(5)
            end if

            ! calculate current Hessian approximation multiplied with trial vector
            inv_hess_x = inv_hess_x + t_mat(1) * kappa_diff - &
                         t_mat(2) * self%y_list(:, it)

            ! get current Hessian approximation multiplied with gradient
            if (update_y) then
                self%y_list(:, self%n_points - 1) = &
                    self%y_list(:, self%n_points - 1) + t_mat(3) * kappa_diff(:) - &
                    t_mat(4) * self%y_list(:, it)
            end if
        end do

        ! deallocate arrays
        if (allocated(kappa_diff)) deallocate(kappa_diff)
        if (allocated(grad_diff)) deallocate(grad_diff)

        ! calculate dot products
        s_mat(1) = ddot(n_param, curr_kappa_diff, 1_ip, curr_grad_diff, 1_ip)
        if (abs(s_mat(1)) < vanish_denom_thr) then
            s_mat(1) = 0.0_rp
        else
            s_mat(1) = 1.0_rp / s_mat(1)
        end if
        s_mat(2) = ddot(n_param, curr_grad_diff, 1_ip, &
                        self%y_list(:, self%n_points - 1), 1_ip)
        if (abs(s_mat(2)) < vanish_denom_thr) then
            s_mat(2) = 1.0_rp
        else
            s_mat(2) = 1.0_rp + s_mat(1) * s_mat(2)
        end if
        s_mat(3) = ddot(n_param, curr_kappa_diff, 1_ip, x, 1_ip)
        s_mat(4) = ddot(n_param, self%y_list(:, self%n_points - 1), 1_ip, x, 1_ip)
        t_mat(2) = s_mat(1) * s_mat(3)
        t_mat(1) = s_mat(2) * t_mat(2) - s_mat(1) * s_mat(4)

        ! calculate final Hessian approximation multiplied with trial vector
        inv_hess_x = inv_hess_x + t_mat(1) * curr_kappa_diff - &
                     t_mat(2) * self%y_list(:, self%n_points - 1)

        ! deallocate arrays
        deallocate(curr_kappa_diff, curr_grad_diff)

    end subroutine bfgs_inv_hess_x

end module otr_s_gek
