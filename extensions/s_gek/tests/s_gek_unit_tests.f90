! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_unit_tests

    use opentrustregion_core, only: rp, ip, stderr, update_orbs_type, obj_func_type, &
                                    hess_x_type, conv_check_type, project_type
    use test_reference, only: is_close, is_below, is_at_most
    use otr_common, only: obj_func_grad_type
    use otr_s_gek, only: change_reference_type
    use opentrustregion_unit_tests, only: test_context_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: quartic_update_orbs_ptr => &
        quartic_update_orbs
    procedure(obj_func_type), pointer :: quartic_obj_func_ptr => quartic_obj_func
    procedure(hess_x_type), pointer :: quartic_hess_x_ptr => quartic_hess_x
    procedure(obj_func_grad_type), pointer :: quartic_obj_func_grad_ptr => &
        quartic_obj_func_grad
    procedure(hess_x_type), pointer :: quartic_init_hess_ptr => quartic_init_hess
    procedure(change_reference_type), pointer :: quartic_change_reference_ptr => &
        quartic_change_reference
    procedure(conv_check_type), pointer :: mock_conv_check_ptr => mock_conv_check
    procedure(project_type), pointer :: quartic_project_ptr => quartic_project

    ! dimension and number of data points of the GEK the tests of its routines set up,
    ! the eigenvalues of the prior of the test of its setup, one negative and one below
    ! the floor of the length scales, and of the other tests, whose length scales keep
    ! data points spread over a few units apart, and the dimension of the full space of
    ! the other tests, of which the tests of the S-GEK step project a window of
    ! n_points points
    integer(ip), parameter :: n_dim = 3, n_points = 3, n_full = 8
    real(rp), parameter :: prior_eigenvalues(n_dim) = [-2.0_rp, 0.01_rp, 1.5_rp], &
                           random_gek_eigenvalues(n_dim) = [-2.0_rp, 0.5_rp, 1.5_rp]

    ! quartic test problem 1/2 y^T A y + d sum(y^4) + offset of y = x - minimum in a
    ! flat chart, whose callback functions move and evaluate the current point x
    real(rp), allocatable :: quartic_x(:), quartic_minimum(:), quartic_a(:, :)
    real(rp) :: quartic_d = 0.0_rp, quartic_offset = 0.0_rp

    ! callback functions of the quartic problem, whose calls are counted, the one the
    ! fault-injection sweep lets fail at its call number fault_call with an error of
    ! code fault_error_code or a value that is not a number, for the chart
    ! transformation in its list number nan_list, the number of calls of any of them
    ! after that, whether one received another context than that of the settings it
    ! belongs to, solver_context or s_gek_context, the number of new trial points the
    ! objective function still raises above the current point, or to infinity, to
    ! reject them, these points, at which it keeps raising it, whether the orbital
    ! update moved to one of them, the last trial point and whether the orbital update
    ! moved elsewhere after one, the largest number of points the chart transformation
    ! received, the moves it saw from the current point, last in the history, and the
    ! others, the trial points it received alone and the raised ones among them, and
    ! the first step from the starting point, which the energy and gradient at a trial
    ! point or the orbital update received, and whether it was evaluated as a trial
    ! point
    integer(ip), parameter :: update_orbs_callback = 1, hess_x_callback = 2, &
                              change_reference_callback = 3, &
                              obj_func_grad_callback = 4, init_hess_callback = 5, &
                              obj_func_callback = 6, n_quartic_callbacks = 6
    integer(ip), parameter :: error_fault = 1, nan_fault = 2, fault_error_code = 7
    real(rp), allocatable :: rejected_points(:, :), first_step(:), last_trial(:)
    integer(ip) :: n_callback_calls(n_quartic_callbacks) = 0, fault_callback = 0, &
                   fault_call = 0, fault_kind = error_fault, nan_list = 1, &
                   n_calls_after_fault = 0, n_rejected_trials = 0, max_history = 0, &
                   moves_from_current = 0, moves_wrong = 0, kept_trials = 0, &
                   kept_rejected = 0
    logical :: fault_injected = .false., context_wrong = .false., &
               infinite_rejections = .false., moved_to_rejected = .false., &
               moved_off_trial = .false., first_step_trial = .false.
    type(test_context_type), target :: solver_context, s_gek_context

    ! the error code the mock convergence check returns, the number of its calls and
    ! the number of the call from which it reports convergence
    integer(ip) :: conv_check_error = 0, n_conv_check_calls = 0, &
                   conv_check_converged_from = 1

    ! the unit direction the projection of the quartic problem removes, if allocated,
    ! along which the callback functions then return no gradient, the largest component
    ! along it of the steps the callback functions received, the number of calls of the
    ! projection and the fault it returns, none, an error or a value that is not a
    ! number
    real(rp), allocatable :: projected_direction(:)
    real(rp) :: max_off_projection = 0.0_rp
    integer(ip) :: n_project_calls = 0, project_fault_kind = 0

contains

    function quartic_func(x) result(func)
        !
        ! this function evaluates the quartic problem at a point
        !
        real(rp), intent(in) :: x(:)
        real(rp) :: func

        real(rp) :: y(size(x))

        y = x - quartic_minimum
        func = 0.5_rp * dot_product(y, matmul(quartic_a, y)) + quartic_d * sum(y**4) + &
               quartic_offset

    end function quartic_func

    function quartic_grad(x) result(grad)
        !
        ! this function evaluates the gradient of the quartic problem at a point
        !
        real(rp), intent(in) :: x(:)
        real(rp) :: grad(size(x))

        real(rp) :: y(size(x))

        y = x - quartic_minimum
        grad = matmul(quartic_a, y) + 4.0_rp * quartic_d * y**3

    end function quartic_grad

    function quartic_hess(x) result(hess)
        !
        ! this function evaluates the Hessian of the quartic problem at a point
        !
        real(rp), intent(in) :: x(:)
        real(rp) :: hess(size(x), size(x))

        integer(ip) :: i

        hess = quartic_a
        do i = 1, size(x)
            hess(i, i) = hess(i, i) + &
                         12.0_rp * quartic_d * (x(i) - quartic_minimum(i))**2
        end do

    end function quartic_hess

    subroutine record_quartic_call(callback, context, expected_context, error, output)
        !
        ! this subroutine records the call of a callback function of the quartic
        ! problem, checks that it received the context of the settings it belongs to
        ! and lets the call fail with an error or a value that is not a number in
        ! output if the fault-injection sweep chose it
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        integer(ip), intent(in) :: callback
        class(*), intent(in), pointer :: context
        type(test_context_type), intent(in), target :: expected_context
        integer(ip), intent(out) :: error
        real(rp), intent(inout) :: output

        if (.not. associated(context, expected_context)) context_wrong = .true.
        if (fault_injected) n_calls_after_fault = n_calls_after_fault + 1
        n_callback_calls(callback) = n_callback_calls(callback) + 1
        error = 0
        if (callback /= fault_callback .or. n_callback_calls(callback) /= fault_call) &
            return
        fault_injected = .true.
        if (fault_kind == error_fault) then
            error = fault_error_code
        else
            output = ieee_value(output, ieee_quiet_nan)
        end if

    end subroutine record_quartic_call

    logical function is_rejected(point)
        !
        ! this function returns whether the objective function of the quartic problem
        ! raised its value at a point to reject it
        !
        real(rp), intent(in) :: point(:)

        integer(ip) :: k

        is_rejected = .false.
        if (.not. allocated(rejected_points)) return
        do k = 1, size(rejected_points, 2)
            if (all(is_close(point, rejected_points(:, k), 0.0_rp))) &
                is_rejected = .true.
        end do

    end function is_rejected

    subroutine quartic_update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                   context)
        !
        ! this subroutine is the orbital update of the quartic problem, which moves the
        ! current point by kappa and records whether it moved elsewhere than to the
        ! trial point evaluated since the last move, the first step from the starting
        ! point and its call
        !
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp) :: hess(size(kappa), size(kappa))
        integer(ip) :: i

        if (is_rejected(quartic_x + kappa)) moved_to_rejected = .true.
        if (allocated(last_trial)) then
            if (.not. all(is_close(quartic_x + kappa, last_trial, 0.0_rp))) &
                moved_off_trial = .true.
            deallocate(last_trial)
        end if
        if (.not. allocated(first_step) .and. &
            n_callback_calls(update_orbs_callback) > 0) first_step = kappa
        call record_off_projection(kappa)
        quartic_x = quartic_x + kappa
        func = quartic_func(quartic_x)
        grad = projected(quartic_grad(quartic_x))
        hess = quartic_hess(quartic_x)
        do i = 1, size(kappa)
            h_diag(i) = hess(i, i)
        end do
        hess_x_funptr => quartic_hess_x
        call record_quartic_call(update_orbs_callback, context, solver_context, error, &
                                 grad(1))

    end subroutine quartic_update_orbs

    subroutine quartic_hess_x(x, hess_x, error, context)
        !
        ! this subroutine is the Hessian linear transformation of the quartic problem
        ! at the current point, which records its call
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        hess_x = matmul(quartic_hess(quartic_x), x)
        call record_quartic_call(hess_x_callback, context, solver_context, error, &
                                 hess_x(1))

    end subroutine quartic_hess_x

    function trial_point_func(kappa) result(func)
        !
        ! this function returns the energy of the quartic problem at a trial point,
        ! raised above that of the current point, or to infinity, at the points to be
        ! rejected, and records the trial point and the first step from the starting
        ! point
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf

        real(rp), intent(in) :: kappa(:)
        real(rp) :: func

        if (.not. allocated(first_step)) then
            first_step = kappa
            first_step_trial = .true.
        end if
        call record_off_projection(kappa)
        last_trial = quartic_x + kappa
        func = quartic_func(quartic_x + kappa)
        if (n_rejected_trials > 0 .and. .not. is_rejected(quartic_x + kappa)) then
            n_rejected_trials = n_rejected_trials - 1
            if (.not. allocated(rejected_points)) &
                allocate(rejected_points(size(kappa), 0))
            rejected_points = reshape([rejected_points, quartic_x + kappa], &
                                      [size(kappa), size(rejected_points, 2) + 1])
        end if
        if (is_rejected(quartic_x + kappa)) then
            func = max(func, quartic_func(quartic_x)) + 1.0_rp
            if (infinite_rejections) func = ieee_value(func, ieee_positive_inf)
        end if

    end function trial_point_func

    function quartic_obj_func(kappa, error, context) result(func)
        !
        ! this function is the objective function of the quartic problem at a trial
        ! point, raised at the points to be rejected, which records the first step from
        ! the starting point and its call
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        func = trial_point_func(kappa)
        call record_quartic_call(obj_func_callback, context, solver_context, error, &
                                 func)

    end function quartic_obj_func

    function quartic_obj_func_grad(kappa, grad, error, context) result(func)
        !
        ! this function is the energy and gradient of the quartic problem at a trial
        ! point, the energy raised at the points to be rejected, which records the
        ! first step from the starting point and its call
        !
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out), target :: grad(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        func = trial_point_func(kappa)
        grad = projected(quartic_grad(quartic_x + kappa))
        call record_quartic_call(obj_func_grad_callback, context, s_gek_context, &
                                 error, grad(1))

    end function quartic_obj_func_grad

    subroutine quartic_init_hess(x, hess_x, error, context)
        !
        ! this subroutine is the static Hessian of the quartic problem, the Hessian of
        ! its quadratic term, which records its call
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        hess_x = matmul(quartic_a, x)
        call record_quartic_call(init_hess_callback, context, s_gek_context, error, &
                                 hess_x(1))

    end subroutine quartic_init_hess

    subroutine quartic_change_reference(new_ref, kappa, local_grad, grad, error, &
                                        context)
        !
        ! this subroutine is the chart transformation of the flat chart of the quartic
        ! problem, which shifts the positions by the new reference and returns the
        ! local gradients as the gradients, records the number of points it received,
        ! whether a move with more than one point in the history started from the
        ! current point, last in the history, and the trial points it received alone
        ! with a vanishing reference, the raised ones among them, and records its call
        !
        real(rp), intent(in), target :: new_ref(:)
        real(rp), intent(inout), target :: kappa(:, :), local_grad(:, :), grad(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(ip) :: k, n

        ! moves with more than one point in the history and trial points received alone
        ! with a vanishing reference
        n = size(kappa, 2)
        max_history = max(max_history, n)
        if (n > 1) then
            if (all(is_at_most(abs(kappa(:, n)), 0.0_rp))) then
                moves_from_current = moves_from_current + 1
            else
                moves_wrong = moves_wrong + 1
            end if
        else if (all(is_at_most(abs(new_ref), 0.0_rp))) then
            kept_trials = kept_trials + 1
            if (is_rejected(quartic_x + kappa(:, 1))) kept_rejected = kept_rejected + 1
        end if

        ! shift the positions
        do k = 1, n
            kappa(:, k) = kappa(:, k) - new_ref
        end do
        grad = local_grad
        select case (nan_list)
        case (2)
            call record_quartic_call(change_reference_callback, context, &
                                     s_gek_context, error, local_grad(1, 1))
        case (3)
            call record_quartic_call(change_reference_callback, context, &
                                     s_gek_context, error, grad(1, 1))
        case default
            call record_quartic_call(change_reference_callback, context, &
                                     s_gek_context, error, kappa(1, 1))
        end select

    end subroutine quartic_change_reference

    function mock_conv_check(error, context) result(converged)
        !
        ! this function is a mock of the convergence check, which reports convergence
        ! from its call number conv_check_converged_from with the error code
        ! conv_check_error, checks the context and counts its calls
        !
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        if (.not. associated(context, solver_context)) context_wrong = .true.
        n_conv_check_calls = n_conv_check_calls + 1
        converged = n_conv_check_calls >= conv_check_converged_from
        error = conv_check_error

    end function mock_conv_check

    subroutine quartic_project(vector, error, context)
        !
        ! this subroutine is the projection of the quartic problem, which removes the
        ! projected direction, checks the context, counts its calls and returns the
        ! fault project_fault_kind
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        if (.not. associated(context, solver_context)) context_wrong = .true.
        n_project_calls = n_project_calls + 1
        vector = projected(vector)
        error = 0
        if (project_fault_kind == error_fault) error = fault_error_code
        if (project_fault_kind == nan_fault) &
            vector(1) = ieee_value(1.0_rp, ieee_quiet_nan)

    end subroutine quartic_project

    function projected(vector) result(projected_vector)
        !
        ! this function removes the projected direction of the quartic problem from a
        ! vector, if one is set
        !
        real(rp), intent(in) :: vector(:)
        real(rp) :: projected_vector(size(vector))

        projected_vector = vector
        if (allocated(projected_direction)) projected_vector = &
            vector - dot_product(projected_direction, vector) * projected_direction

    end function projected

    subroutine record_off_projection(step)
        !
        ! this subroutine records the component of a step along the projected direction
        ! of the quartic problem, if one is set
        !
        real(rp), intent(in) :: step(:)

        if (allocated(projected_direction)) max_off_projection = &
            max(max_off_projection, abs(dot_product(projected_direction, step)))

    end subroutine record_off_projection

    function generate_random_spd_matrix(n) result(matrix)
        !
        ! this function generates a random symmetric positive definite matrix with
        ! eigenvalues of at least one
        !
        use opentrustregion_unit_tests, only: identity_matrix

        integer(ip), intent(in) :: n
        real(rp) :: matrix(n, n)

        real(rp) :: r(n, n)

        call random_number(r)
        matrix = matmul(r, transpose(r)) + identity_matrix(n)

    end function generate_random_spd_matrix

    function generate_random_prior(eigenvalues) result(prior)
        !
        ! this function generates a random prior Hessian of the GEK with the given
        ! eigenvalues
        !
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix, &
                                         diagonal_matrix

        real(rp), intent(in) :: eigenvalues(n_dim)
        real(rp) :: prior(n_dim, n_dim)

        real(rp) :: eigvecs(n_dim, n_dim)

        eigvecs = generate_random_orthogonal_matrix(n_dim)
        prior = matmul(eigvecs, &
                       matmul(diagonal_matrix(eigenvalues), transpose(eigvecs)))

    end function generate_random_prior

    subroutine ref_setup_gek(gek, kappa_list, grad_list, func_list, hess)
        !
        ! this subroutine sets up the GEK of the positions kappa_list, gradients
        ! grad_list and energies func_list of data points with a prior hess: eigenbasis
        ! and length scales of the prior, the Cholesky factor of the covariance matrix
        ! and the solution vector of the energies relative to the trend and the
        ! gradients
        !
        use otr_s_gek, only: gek_type, trend_bias, prior_eigval_floor

        type(gek_type), intent(out) :: gek
        real(rp), intent(in) :: kappa_list(:, :), grad_list(:, :), func_list(:), &
                                hess(:, :)

        real(rp), allocatable :: eigvals(:), work(:), b(:), gs(:, :)
        integer(ip) :: n, p, m, i, lwork, info
        external :: dsyev, dpotrf, dpotrs

        n = size(kappa_list, 1)
        p = size(kappa_list, 2)
        m = p * (n + 1)
        gek%n_dim = n
        gek%n_points = p

        ! eigenbasis and length scales
        gek%eigvecs = 0.5_rp * (hess + transpose(hess))
        allocate(eigvals(n), work(1))
        call dsyev("V", "U", n, gek%eigvecs, n, eigvals, work, -1_ip, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))
        call dsyev("V", "U", n, gek%eigvecs, n, eigvals, work, lwork, info)
        gek%length_scales = sqrt(5.0_rp / 3.0_rp * trend_bias / &
                                 max(abs(eigvals), prior_eigval_floor))
        gek%kappa_list = matmul(transpose(gek%eigvecs), kappa_list)

        ! Cholesky factor of the covariance matrix and solution vector
        gek%cov_chol = ref_covariance_gek(gek)
        call dpotrf("U", m, gek%cov_chol, m, info)
        gek%trend = maxval(func_list) + trend_bias
        gs = matmul(transpose(gek%eigvecs), grad_list)
        allocate(b(m))
        b(:p) = func_list - gek%trend
        do i = 1, n
            b(p + (i - 1) * p + 1:p + i * p) = gs(i, :)
        end do
        gek%weights = b
        call dpotrs("U", m, 1_ip, gek%cov_chol, m, gek%weights, m, info)
        gek%process_variance = dot_product(b, gek%weights) / real(m, kind=rp)
        gek%excess_dispersion = .false.
        gek%current_dispersion = 0.0_rp

    end subroutine ref_setup_gek

    function ref_covariance_gek(gek) result(cov)
        !
        ! this function returns the regularized covariance matrix of the energies and
        ! gradients of the data points of a GEK from the derivatives of the kernel: the
        ! covariance of two energies is the kernel, that of an energy and a gradient
        ! its derivative with respect to the position of the gradient and that of two
        ! gradients its mixed second derivative
        !
        use otr_s_gek, only: gek_type, func_regularization, grad_regularization

        type(gek_type), intent(in) :: gek
        real(rp) :: cov(gek%n_points * (gek%n_dim + 1), gek%n_points * (gek%n_dim + 1))

        real(rp) :: kernel, d1(gek%n_dim), d2(gek%n_dim, gek%n_dim), &
                    d3(gek%n_dim, gek%n_dim, gek%n_dim)
        integer(ip) :: p, i, j, k, mm

        p = gek%n_points
        do j = 1, p
            do k = 1, p
                call ref_matern(gek%kappa_list(:, k), gek%kappa_list(:, j), &
                                gek%length_scales, kernel, d1, d2, d3)
                cov(k, j) = kernel
                do i = 1, gek%n_dim
                    cov(k, p + (i - 1) * p + j) = -d1(i)
                    cov(p + (i - 1) * p + k, j) = d1(i)
                    do mm = 1, gek%n_dim
                        cov(p + (i - 1) * p + k, p + (mm - 1) * p + j) = -d2(i, mm)
                    end do
                end do
            end do
        end do
        do k = 1, p
            cov(k, k) = cov(k, k) + func_regularization
        end do
        do k = p + 1, size(cov, 1)
            cov(k, k) = cov(k, k) + grad_regularization
        end do

    end function ref_covariance_gek

    function ref_covariance_vector_gek(gek, eig_kappa) result(c)
        !
        ! this function returns the covariance of the energy at a point eig_kappa of
        ! the eigenbasis of the prior with the energies and gradients of the data
        ! points of a GEK
        !
        use otr_s_gek, only: gek_type

        type(gek_type), intent(in) :: gek
        real(rp), intent(in) :: eig_kappa(:)
        real(rp) :: c(gek%n_points * (gek%n_dim + 1))

        real(rp) :: kernel, d1(gek%n_dim), d2(gek%n_dim, gek%n_dim), &
                    d3(gek%n_dim, gek%n_dim, gek%n_dim)
        integer(ip) :: p, i, k

        p = gek%n_points
        do k = 1, p
            call ref_matern(eig_kappa, gek%kappa_list(:, k), gek%length_scales, &
                            kernel, d1, d2, d3)
            c(k) = kernel
            do i = 1, gek%n_dim
                c(p + (i - 1) * p + k) = -d1(i)
            end do
        end do

    end function ref_covariance_vector_gek

    function ref_energy_gek(gek, kappa) result(energy)
        !
        ! this function returns the energy a GEK predicts at a point kappa, the trend
        ! plus the covariance vector of the point times the solution vector
        !
        use otr_s_gek, only: gek_type

        type(gek_type), intent(in) :: gek
        real(rp), intent(in) :: kappa(:)
        real(rp) :: energy

        energy = gek%trend + dot_product(ref_covariance_vector_gek( &
            gek, matmul(transpose(gek%eigvecs), kappa)), gek%weights)

    end function ref_energy_gek

    function ref_gradient_gek(gek, kappa) result(gradient)
        !
        ! this function returns the gradient a GEK predicts at a point kappa, the
        ! derivative of the covariance vector times the solution vector
        !
        use otr_s_gek, only: gek_type

        type(gek_type), intent(in) :: gek
        real(rp), intent(in) :: kappa(:)
        real(rp) :: gradient(gek%n_dim)

        real(rp) :: eig_kappa(gek%n_dim), eig_gradient(gek%n_dim), kernel, &
                    d1(gek%n_dim), d2(gek%n_dim, gek%n_dim), &
                    d3(gek%n_dim, gek%n_dim, gek%n_dim)
        integer(ip) :: p, i, k

        p = gek%n_points
        eig_kappa = matmul(transpose(gek%eigvecs), kappa)
        eig_gradient = 0.0_rp
        do k = 1, p
            call ref_matern(eig_kappa, gek%kappa_list(:, k), gek%length_scales, &
                            kernel, d1, d2, d3)
            eig_gradient = eig_gradient + gek%weights(k) * d1
            do i = 1, gek%n_dim
                eig_gradient = eig_gradient - &
                               gek%weights(p + (i - 1) * p + k) * d2(:, i)
            end do
        end do
        gradient = matmul(gek%eigvecs, eig_gradient)

    end function ref_gradient_gek

    function ref_hessian_gek(gek, kappa) result(hessian)
        !
        ! this function returns the Hessian a GEK predicts at a point kappa, the second
        ! derivative of the covariance vector times the solution vector
        !
        use otr_s_gek, only: gek_type

        type(gek_type), intent(in) :: gek
        real(rp), intent(in) :: kappa(:)
        real(rp) :: hessian(gek%n_dim, gek%n_dim)

        real(rp) :: eig_kappa(gek%n_dim), eig_hessian(gek%n_dim, gek%n_dim), kernel, &
                    d1(gek%n_dim), d2(gek%n_dim, gek%n_dim), &
                    d3(gek%n_dim, gek%n_dim, gek%n_dim)
        integer(ip) :: p, i, k

        p = gek%n_points
        eig_kappa = matmul(transpose(gek%eigvecs), kappa)
        eig_hessian = 0.0_rp
        do k = 1, p
            call ref_matern(eig_kappa, gek%kappa_list(:, k), gek%length_scales, &
                            kernel, d1, d2, d3)
            eig_hessian = eig_hessian + gek%weights(k) * d2
            do i = 1, gek%n_dim
                eig_hessian = eig_hessian - &
                              gek%weights(p + (i - 1) * p + k) * d3(:, :, i)
            end do
        end do
        hessian = matmul(gek%eigvecs, matmul(eig_hessian, transpose(gek%eigvecs)))

    end function ref_hessian_gek

    function ref_dispersion_gek(gek, kappa) result(dispersion)
        !
        ! this function returns the half width of the 95% confidence interval of the
        ! energy a GEK predicts at a point kappa from the posterior variance of the
        ! kernel, or its excess over the dispersion at the current point
        !
        use otr_s_gek, only: gek_type, confidence_factor

        type(gek_type), intent(in) :: gek
        real(rp), intent(in) :: kappa(:)
        real(rp) :: dispersion

        real(rp) :: c(gek%n_points * (gek%n_dim + 1)), &
                    solved(gek%n_points * (gek%n_dim + 1)), &
                    cov(gek%n_points * (gek%n_dim + 1), gek%n_points * (gek%n_dim + 1))
        integer(ip) :: m, info
        external :: dpotrf, dpotrs

        m = gek%n_points * (gek%n_dim + 1)
        c = ref_covariance_vector_gek(gek, matmul(transpose(gek%eigvecs), kappa))
        cov = ref_covariance_gek(gek)
        solved = c
        call dpotrf("U", m, cov, m, info)
        call dpotrs("U", m, 1_ip, cov, m, solved, m, info)
        dispersion = confidence_factor * sqrt(max(1.0_rp - dot_product(c, solved), &
                                                  0.0_rp) * gek%process_variance)
        if (gek%excess_dispersion) &
            dispersion = max(0.0_rp, dispersion - gek%current_dispersion)

    end function ref_dispersion_gek

    subroutine ref_matern(kappa, data_kappa, length_scales, kernel, d1, d2, d3)
        !
        ! this subroutine returns the Matern 5/2 kernel of the scaled distance r of a
        ! point kappa from a data point data_kappa and its first three derivatives with
        ! respect to kappa, from the closed-form derivatives of the kernel with respect
        ! to r
        !
        real(rp), intent(in) :: kappa(:), data_kappa(:), length_scales(:)
        real(rp), intent(out) :: kernel, d1(:), d2(:, :), d3(:, :, :)

        real(rp), parameter :: a = sqrt(5.0_rp)
        real(rp) :: scaled_diff(size(kappa)), r, decay, f1, f2, f3
        integer(ip) :: i, j, k

        scaled_diff = (kappa - data_kappa) / length_scales**2
        r = sqrt(sum(((kappa - data_kappa) / length_scales)**2))
        decay = exp(-a * r)
        kernel = (1.0_rp + a * r + a**2 * r**2 / 3.0_rp) * decay

        ! derivative of the kernel divided by r, the derivative of that divided by r
        ! and the derivative of that divided by r, which vanishes at r = 0 together
        ! with the terms it enters
        f1 = -a**2 / 3.0_rp * (1.0_rp + a * r) * decay
        f2 = a**4 / 3.0_rp * decay
        f3 = 0.0_rp
        if (r > 0.0_rp) f3 = -a**5 / 3.0_rp * decay / r

        ! derivatives with respect to the point
        d1 = f1 * scaled_diff
        do j = 1, size(kappa)
            do i = 1, size(kappa)
                d2(i, j) = f2 * scaled_diff(i) * scaled_diff(j)
                do k = 1, size(kappa)
                    d3(i, j, k) = f3 * scaled_diff(i) * scaled_diff(j) * scaled_diff(k)
                    if (i == k) d3(i, j, k) = d3(i, j, k) + f2 * scaled_diff(j) / &
                                              length_scales(i)**2
                    if (j == k) d3(i, j, k) = d3(i, j, k) + f2 * scaled_diff(i) / &
                                              length_scales(j)**2
                    if (i == j) d3(i, j, k) = d3(i, j, k) + f2 * scaled_diff(k) / &
                                              length_scales(i)**2
                end do
            end do
            d2(j, j) = d2(j, j) + f1 / length_scales(j)**2
        end do

    end subroutine ref_matern

    subroutine setup_quartic(a, d, offset, minimum)
        !
        ! this subroutine sets up the quartic problem with the current point at the
        ! origin and clears the records of the calls of its callback functions
        !
        real(rp), intent(in) :: a(:, :), d, offset, minimum(:)

        quartic_a = a
        quartic_d = d
        quartic_offset = offset
        quartic_minimum = minimum
        quartic_x = 0.0_rp * minimum
        n_callback_calls = 0
        fault_callback = 0
        fault_call = 0
        fault_kind = error_fault
        nan_list = 1
        n_calls_after_fault = 0
        fault_injected = .false.
        context_wrong = .false.
        n_rejected_trials = 0
        infinite_rejections = .false.
        if (allocated(rejected_points)) deallocate(rejected_points)
        moved_to_rejected = .false.
        if (allocated(last_trial)) deallocate(last_trial)
        moved_off_trial = .false.
        max_history = 0
        moves_from_current = 0
        moves_wrong = 0
        kept_trials = 0
        kept_rejected = 0
        if (allocated(first_step)) deallocate(first_step)
        first_step_trial = .false.
        if (allocated(projected_direction)) deallocate(projected_direction)
        max_off_projection = 0.0_rp
        n_project_calls = 0
        project_fault_kind = 0

    end subroutine setup_quartic

    subroutine setup_random_gek(gek, kappa_list, grad_list, func_list)
        !
        ! this subroutine sets up a GEK of random data points with random gradients and
        ! energies and a random prior with the eigenvalues random_gek_eigenvalues
        !
        use otr_s_gek, only: gek_type

        type(gek_type), intent(out) :: gek
        real(rp), intent(out) :: kappa_list(n_dim, n_points), &
                                 grad_list(n_dim, n_points), func_list(n_points)

        call random_number(kappa_list)
        kappa_list = 4.0_rp * kappa_list - 2.0_rp
        call random_number(grad_list)
        grad_list = 2.0_rp * grad_list - 1.0_rp
        call random_number(func_list)
        call ref_setup_gek(gek, kappa_list, grad_list, func_list, &
                           generate_random_prior(random_gek_eigenvalues))

    end subroutine setup_random_gek

    subroutine setup_quadratic_gek(gek, current_grad, minimum)
        !
        ! this subroutine sets up a GEK of data points of a random convex quadratic,
        ! whose Hessian is the prior, at random points around the current point, the
        ! origin and last data point, and returns its gradient and the minimum of the
        ! quadratic
        !
        use otr_s_gek, only: gek_type

        type(gek_type), intent(out) :: gek
        real(rp), intent(out) :: current_grad(n_dim), minimum(n_dim)

        real(rp) :: hess(n_dim, n_dim), kappa_list(n_dim, n_points), &
                    grad_list(n_dim, n_points), func_list(n_points)
        integer(ip) :: k

        hess = generate_random_spd_matrix(n_dim)
        call random_number(minimum)
        minimum = 0.5_rp * (minimum - 0.5_rp) / norm2(minimum - 0.5_rp)
        call random_number(kappa_list)
        kappa_list = 0.4_rp * kappa_list - 0.2_rp
        kappa_list(:, n_points) = 0.0_rp
        do k = 1, n_points
            grad_list(:, k) = matmul(hess, kappa_list(:, k) - minimum)
            func_list(k) = 0.5_rp * &
                           dot_product(kappa_list(:, k) - minimum, grad_list(:, k))
        end do
        call ref_setup_gek(gek, kappa_list, grad_list, func_list, hess)
        current_grad = grad_list(:, n_points)

    end subroutine setup_quadratic_gek

    function roundoff_tolerance(gek) result(tolerance)
        !
        ! this function returns the tolerance of the values a GEK predicts, which carry
        ! the rounding errors of its solution vector
        !
        use otr_s_gek, only: gek_type

        type(gek_type), intent(in) :: gek
        real(rp) :: tolerance

        tolerance = 1e-12_rp * max(1.0_rp, sum(abs(gek%weights)))

    end function roundoff_tolerance

    logical(c_bool) function test_s_gek_solver() bind(C)
        !
        ! this function tests the S-GEK solver on the quartic problem in a flat chart:
        ! its convergence with the default and with OpenMolcas's options, with a static
        ! prior and with a rejected trial point, the first step, the contexts its
        ! callback functions receive, the order of its history, the state, the
        ! convergence check, the escape step of the stability check, the gradient
        ! noise, maximum precision, the maximum number of iterations, its sanity checks
        ! and the origins of the errors of every callback call
        !
        use otr_s_gek, only: s_gek_solver, s_gek_settings_type, error_s_gek, &
                             error_change_reference, max_history_points
        use otr_common, only: error_obj_func_grad, error_init_hess
        use opentrustregion_core, only: &
            solver_settings_type, solver_state_type, error_update_orbs, error_hess_x, &
            error_conv_check, error_stability_check, error_solver_max_iter, &
            error_obj_func, error_project, function_unchanged_warning_msg
        use opentrustregion_unit_tests, only: setup_settings, setup_error_logging

        integer(ip), parameter :: n_quartic = 4
        character(len=37), parameter :: rejection_names(3) = &
            [character(len=37) :: "rejected trial point", &
             "trial point of infinite energy", "rejected trial point without gradient"]
        character(len=26), parameter :: acceptance_names(2) = &
            [character(len=26) :: "with trial acceptance", "without trial acceptance"]
        character(len=31), parameter :: noise_names(4) = &
            [character(len=31) :: "without gradient noise", "with gradient noise", &
             "with gradient noise only global", "with gradient noise only local"]
        character(len=16), parameter :: callback_names(n_quartic_callbacks) = &
            [character(len=16) :: "update_orbs", "hess_x", "change_reference", &
             "obj_func_grad", "init_hess", "obj_func"]
        type(solver_settings_type) :: settings
        type(s_gek_settings_type) :: gek_settings
        type(solver_state_type) :: state
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(change_reference_type), pointer :: change_reference_funptr
        real(rp) :: convex_a(n_quartic, n_quartic), saddle_a(n_quartic, n_quartic), &
                    far_a(n_quartic, n_quartic), diagonal_a(n_quartic, n_quartic), &
                    minimum(n_quartic), proposal(n_quartic), &
                    start_hess(n_quartic, n_quartic)
        integer(ip) :: error, callback, i_call, kind, expected, i_case, j, n_evals, &
                       stop_check, n_calls(n_quartic_callbacks)
        character(len=80) :: fault_name
        character(len=37) :: case_name

        ! assume tests pass
        test_s_gek_solver = .true.

        update_orbs_funptr => quartic_update_orbs
        obj_func_funptr => quartic_obj_func
        change_reference_funptr => quartic_change_reference
        convex_a = generate_random_spd_matrix(n_quartic)
        far_a = 0.0_rp
        diagonal_a = 0.0_rp
        do i_case = 1, n_quartic
            far_a(i_case, i_case) = 10.0_rp**(i_case - 3)
            diagonal_a(i_case, i_case) = merge(0.5_rp * real(i_case, kind=rp), &
                                               0.02_rp, i_case < n_quartic)
        end do
        call random_number(minimum)
        minimum = minimum - 0.5_rp

        ! default options: S-GEK reaches the minimum, every callback function receives
        ! the context of its settings, every step is taken to an accepted trial point
        ! from the current point, last in the history, the state counts the calls and
        ! the rows are logged
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        call run()
        call check_converged("default options")
        if (moves_from_current == 0 .or. moves_wrong > 0 .or. moved_off_trial .or. &
            max_history /= 1 + state%n_obj_func) then
            write(stderr, *) "test_s_gek_solver failed: Step not taken to an "// &
                "accepted trial point from the current point or the trial point "// &
                "kept in the history."
            test_s_gek_solver = .false.
        end if
        if (state%n_update_orbs /= n_callback_calls(update_orbs_callback) .or. &
            state%n_obj_func /= n_callback_calls(obj_func_grad_callback)) then
            write(stderr, *) "test_s_gek_solver failed: State does not count the calls."
            test_s_gek_solver = .false.
        end if
        if (index(solver_context%log_message, "Iteration") == 0 .or. &
            index(solver_context%log_message, &
                  "Total number of objective function evaluations") == 0) then
            write(stderr, *) "test_s_gek_solver failed: Header or totals not logged."
            test_s_gek_solver = .false.
        end if

        ! OpenMolcas's options, which need no energy and gradient at trial points,
        ! without the gradient noise OpenMolcas lacks
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        settings%grad_noise = .false.
        gek_settings%obj_func_grad => null()
        gek_settings%inner_solver = "rs-rfo"
        gek_settings%prior_update = "none"
        gek_settings%acceptance = "none"
        gek_settings%excess_dispersion = .false.
        call run()
        call check_converged("OpenMolcas's options")
        if (kept_trials > 0 .or. moves_wrong > 0) then
            write(stderr, *) "test_s_gek_solver failed: Trial point kept without "// &
                "trial acceptance."
            test_s_gek_solver = .false.
        end if

        ! the first step, without the gradient noise, is the proposal of the floored
        ! Hessian diagonal, here the Hessian diagonal of the convex problem, which
        ! S-GEK moves along without evaluating it as a trial point
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        settings%grad_noise = .false.
        start_hess = quartic_hess(quartic_x)
        proposal = -quartic_grad(quartic_x) / [(start_hess(j, j), j=1, n_quartic)]
        call run()
        if (error /= 0 .or. context_wrong .or. .not. allocated(first_step)) then
            write(stderr, *) "test_s_gek_solver failed: First step not taken "// &
                "without error."
            test_s_gek_solver = .false.
        else if (.not. all(is_close(first_step, proposal)) .or. first_step_trial) then
            write(stderr, *) "test_s_gek_solver failed: Wrong first step."
            test_s_gek_solver = .false.
        end if

        ! static prior
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        gek_settings%init_hess => quartic_init_hess
        call run()
        call check_converged("static prior")
        if (n_callback_calls(init_hess_callback) == 0) then
            write(stderr, *) "test_s_gek_solver failed: Static prior not used."
            test_s_gek_solver = .false.
        end if

        ! a rejected trial point, whose row the log marks, stays in the history before
        ! the current point, while one with an infinite energy leaves the history
        ! unchanged, after which the halved variance threshold moves the next trial
        ! point, which without the gradient noise would otherwise repeat the rejected
        ! one; the energy raised at a point without its gradient can leave the
        ! surrogate without a step, at which S-GEK stops at maximum precision, so only
        ! the infinite energy, which hosts return for steps they rule out, has to reach
        ! the minimum; without the energy and gradient at trial points, the objective
        ! function of the solver evaluates them, whose rejected point leaves the
        ! history unchanged and whose accepted one S-GEK moves to from the current point
        do i_case = 1, 3
            call setup_case(convex_a, 0.1_rp, 0.0_rp)
            settings%grad_noise = .false.
            n_rejected_trials = 1
            infinite_rejections = i_case == 2
            if (i_case == 3) gek_settings%obj_func_grad => null()
            case_name = rejection_names(i_case)
            call run()
            if (i_case == 1) then
                if (error /= 0 .or. context_wrong) then
                    write(stderr, *) "test_s_gek_solver failed: Produced error "// &
                        "with "//trim(case_name)//"."
                    test_s_gek_solver = .false.
                end if
            else
                call check_converged(trim(case_name))
            end if
            if (n_rejected_trials > 0 .or. moves_wrong > 0 .or. &
                kept_rejected /= merge(1, 0, i_case == 1) .or. &
                (i_case == 3 .and. kept_trials > 0)) then
                write(stderr, *) "test_s_gek_solver failed: Rejected trial point "// &
                    "not kept, or kept with an infinite energy or without its "// &
                    "gradient, with "//trim(case_name)//"."
                test_s_gek_solver = .false.
            end if
            if (index(solver_context%log_message, "|  x") == 0) then
                write(stderr, *) "test_s_gek_solver failed: Row of the rejected "// &
                    "trial point not marked with "//trim(case_name)//"."
                test_s_gek_solver = .false.
            end if
            if (moved_to_rejected .or. moved_off_trial) then
                write(stderr, *) "test_s_gek_solver failed: Moved to a rejected "// &
                    "point or elsewhere than to the trial point with "// &
                    trim(case_name)//"."
                test_s_gek_solver = .false.
            end if
            if (i_case == 3 .and. &
                (state%n_obj_func /= n_callback_calls(obj_func_callback) .or. &
                 n_callback_calls(obj_func_grad_callback) > 0)) then
                write(stderr, *) "test_s_gek_solver failed: Trial points not "// &
                    "evaluated by the objective function of the solver or its "// &
                    "calls not counted with "//trim(case_name)//"."
                test_s_gek_solver = .false.
            end if
        end do

        ! a run of rejected trial points kept before the current point that fills the
        ! history and drops every point S-GEK moved to before the current point, whose
        ! energy the variance threshold still takes as that of the last move
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        settings%grad_noise = .false.
        n_rejected_trials = max_history_points - 1
        call run()
        if (error /= 0 .or. n_rejected_trials > 0 .or. &
            kept_rejected /= max_history_points - 1) then
            write(stderr, *) "test_s_gek_solver failed: History of rejected trial "// &
                "points without a previous point moved to not handled."
            test_s_gek_solver = .false.
        end if

        ! an error of the objective function of the solver at a trial point is reported
        ! with its origin
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        gek_settings%obj_func_grad => null()
        fault_callback = obj_func_callback
        fault_call = 1
        call run()
        fault_callback = 0
        if (error /= error_obj_func + fault_error_code .or. n_calls_after_fault > 0) &
            then
            write(stderr, *) "test_s_gek_solver failed: Error of the objective "// &
                "function at a trial point not reported with its origin at once."
            test_s_gek_solver = .false.
        end if

        ! with a projection that removes a direction mixing two coordinates, along
        ! which the callback functions return no gradient, from a start whose offset
        ! from the minimum lies in the projected space, every step S-GEK takes, the
        ! proposal of the floored Hessian diagonal among them, lies in that space, and
        ! an error of the projection is reported with its origin at once and values
        ! that are not finite with the origin of S-GEK
        do i_case = 1, 3
            call setup_case(convex_a, 0.1_rp, 0.0_rp)
            projected_direction = [1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp] / sqrt(2.0_rp)
            quartic_x = quartic_minimum + projected(quartic_x - quartic_minimum)
            settings%project => quartic_project
            project_fault_kind = &
                merge(0_ip, merge(error_fault, nan_fault, i_case == 2), i_case == 1)
            call run()
            select case (i_case)
            case (1)
                call check_converged("projection")
                if (.not. is_at_most(max_off_projection, 1e-10_rp) .or. &
                    n_project_calls == 0) then
                    write(stderr, *) "test_s_gek_solver failed: Step taken out of "// &
                        "the space of the projection."
                    test_s_gek_solver = .false.
                end if
            case (2)
                if (error /= error_project + fault_error_code .or. &
                    state%n_update_orbs /= 1) then
                    write(stderr, *) "test_s_gek_solver failed: Error of the "// &
                        "projection not reported with its origin at once."
                    test_s_gek_solver = .false.
                end if
            case (3)
                if (error /= error_s_gek + 1 .or. state%n_update_orbs /= 1) then
                    write(stderr, *) "test_s_gek_solver failed: Projection "// &
                        "returning values that are not finite not reported with "// &
                        "the origin of S-GEK at once."
                    test_s_gek_solver = .false.
                end if
            end select
        end do

        ! far from the minimum of an ill-conditioned problem, S-GEK reaches a
        ! convergence threshold the surrogate resolves, which cannot resolve the energy
        ! changes near the minimum far from the energies of the start, with a history
        ! of the starting point, the first step and every trial point, all of them
        ! accepted or kept, of which the move to the last point sees all but that
        ! point, up to its maximum number of points, beyond which, as almost every run
        ! that far from the minimum takes, it drops the oldest
        call setup_case(far_a, 0.1_rp, 0.0_rp)
        quartic_minimum = 30.0_rp * minimum
        settings%n_macro = 300
        settings%conv_tol = 1e-5_rp
        call run()
        if (error /= 0 .or. .not. is_below(state%grad_rms, settings%conv_tol) .or. &
            context_wrong .or. moves_wrong > 0) then
            write(stderr, *) "test_s_gek_solver failed: Minimum not reached far "// &
                "from it."
            test_s_gek_solver = .false.
        end if
        if (max_history /= min(max_history_points, 1 + state%n_obj_func)) then
            write(stderr, *) "test_s_gek_solver failed: History not filled up to "// &
                "its maximum number of points far from the minimum."
            test_s_gek_solver = .false.
        end if

        ! the convergence check stops at the starting point and its error is reported
        ! with its origin
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        settings%conv_check => mock_conv_check
        conv_check_error = 0
        n_conv_check_calls = 0
        call run()
        if (error /= 0 .or. state%n_update_orbs /= 1 .or. n_conv_check_calls /= 1 .or. &
            context_wrong) then
            write(stderr, *) "test_s_gek_solver failed: Convergence check not used."
            test_s_gek_solver = .false.
        end if
        conv_check_error = 3
        call run()
        conv_check_error = 0
        if (error /= error_conv_check + 3) then
            write(stderr, *) "test_s_gek_solver failed: Error of the convergence "// &
                "check not reported with its origin."
            test_s_gek_solver = .false.
        end if

        ! at the saddle point of a double well, without the gradient noise, the
        ! stability check finds the unstable mode, along which S-GEK escapes to one of
        ! the minima
        saddle_a = 0.0_rp
        saddle_a(1, 1) = -2.0_rp
        saddle_a(2, 2) = 3.0_rp
        saddle_a(3, 3) = 1.0_rp
        saddle_a(4, 4) = 4.0_rp
        call setup_case(saddle_a, 1.0_rp, 0.0_rp)
        quartic_x(1) = quartic_minimum(1)
        settings%stability = .true.
        settings%grad_noise = .false.
        call run()
        if (error /= 0 .or. .not. is_close(abs(quartic_x(1) - quartic_minimum(1)), &
                                           sqrt(0.5_rp), 1e-4_rp) .or. &
            .not. all(is_close(quartic_x(2:), quartic_minimum(2:), 1e-4_rp)) .or. &
            n_callback_calls(hess_x_callback) == 0) then
            write(stderr, *) "test_s_gek_solver failed: No escape from the saddle "// &
                "point to a minimum."
            test_s_gek_solver = .false.
        end if

        ! exactly at the saddle point, whose gradient vanishes, and displaced from it
        ! along a stable mode so little that the gradient is converged, S-GEK performs
        ! the stability check also when the settings do not ask for it, as at any
        ! converged starting point, and escapes to one of the minima
        do i_case = 1, 2
            call setup_case(saddle_a, 1.0_rp, 0.0_rp)
            quartic_x = quartic_minimum
            if (i_case == 2) quartic_x(2) = quartic_x(2) + 1e-7_rp
            settings%grad_noise = .false.
            call run()
            if (error /= 0 .or. .not. is_close(abs(quartic_x(1) - quartic_minimum(1)), &
                                               sqrt(0.5_rp), 1e-4_rp) .or. &
                .not. all(is_close(quartic_x(2:), quartic_minimum(2:), 1e-4_rp)) .or. &
                n_callback_calls(hess_x_callback) == 0) then
                write(stderr, *) "test_s_gek_solver failed: No escape from the "// &
                    "saddle point it starts "// &
                    trim(merge("at  ", "near", i_case == 1))// &
                    " without the stability check of the settings."
                test_s_gek_solver = .false.
            end if
        end do

        ! with an energy so large that the energy changes of the steps onto the saddle
        ! point can not be resolved, S-GEK stops at maximum precision near it, where
        ! the stability check finds the unstable mode, and continues after the escape
        ! step, which is no step whose energy change can not be resolved, to a minimum
        call setup_case(saddle_a, 1.0_rp, 1e8_rp)
        quartic_x(1) = quartic_minimum(1)
        settings%stability = .true.
        settings%grad_noise = .false.
        settings%conv_tol = 0.0_rp
        call run()
        if (error /= 0 .or. .not. is_close(abs(quartic_x(1) - quartic_minimum(1)), &
                                           sqrt(0.5_rp), 1e-2_rp) .or. &
            n_callback_calls(hess_x_callback) == 0) then
            write(stderr, *) "test_s_gek_solver failed: No escape to a minimum "// &
                "after maximum precision near the saddle point."
            test_s_gek_solver = .false.
        end if

        ! from a point on the saddle point's symmetry subspace, which the gradient
        ! never leaves, S-GEK converges without the gradient noise and the stability
        ! check onto the saddle point and with the noise, here a large one, so that it
        ! leaves the subspace before converging, to one of the minima, also with noise
        ! only in the global region, after the first step and a step the variance
        ! threshold restricts, or only in the local region, after a step within it
        do i_case = 1, 4
            call setup_case(saddle_a, 1.0_rp, 0.0_rp)
            quartic_x(1) = quartic_minimum(1)
            settings%grad_noise = i_case > 1
            settings%global_red_factor = merge(0.0_rp, 0.1_rp, i_case == 4)
            settings%local_red_factor = merge(0.0_rp, 0.1_rp, i_case == 3)
            case_name = noise_names(i_case)
            call run()
            if (error /= 0 .or. &
                .not. is_close(abs(quartic_x(1) - quartic_minimum(1)), &
                               merge(0.0_rp, sqrt(0.5_rp), i_case == 1), 1e-4_rp) .or. &
                .not. all(is_close(quartic_x(2:), quartic_minimum(2:), 1e-4_rp))) then
                write(stderr, *) "test_s_gek_solver failed: Wrong stationary point "// &
                    "reached from the symmetry subspace "//trim(case_name)//"."
                test_s_gek_solver = .false.
            end if
        end do

        ! displaced from the minimum of a quadratic problem with a diagonal Hessian
        ! along a coordinate whose curvature lies below the floor of the Hessian
        ! diagonal, so that the first step, without the gradient noise, only shortens
        ! the displacement, with an energy so large that the energies of the points
        ! round to the same value and a convergence threshold that can not be met, the
        ! next step, whose energy change and its prediction can not be resolved, stops
        ! at maximum precision with a warning, with the trial acceptance without taking
        ! the step and without it after taking it
        do i_case = 1, 2
            call setup_case(diagonal_a, 0.0_rp, 1e8_rp)
            quartic_x = quartic_minimum
            quartic_x(n_quartic) = quartic_x(n_quartic) + 1e-5_rp
            settings%grad_noise = .false.
            settings%conv_tol = 0.0_rp
            if (i_case == 2) gek_settings%acceptance = "none"
            case_name = acceptance_names(i_case)
            call run()
            if (error /= 0 .or. .not. state%max_precision_reached .or. index( &
                solver_context%log_message, function_unchanged_warning_msg) == 0) then
                write(stderr, *) "test_s_gek_solver failed: Maximum precision not "// &
                    "reached with a warning "//trim(case_name)//"."
                test_s_gek_solver = .false.
            end if
            if (state%n_update_orbs /= i_case + 1) then
                write(stderr, *) "test_s_gek_solver failed: Step whose energy "// &
                    "change can not be resolved taken wrong "//trim(case_name)//"."
                test_s_gek_solver = .false.
            end if
        end do

        ! convergence takes precedence over maximum precision: a convergence check that
        ! passes first at the point at which the same run stops at maximum precision,
        ! the number of whose convergence check the first run finds, ends the second
        ! converged, without the warning
        do i_case = 1, 2
            call setup_case(diagonal_a, 0.0_rp, 1e8_rp)
            quartic_x = quartic_minimum
            quartic_x(n_quartic) = quartic_x(n_quartic) + 1e-5_rp
            settings%grad_noise = .false.
            settings%conv_tol = 0.0_rp
            if (i_case == 2) then
                settings%conv_check => mock_conv_check
                conv_check_converged_from = stop_check
                n_conv_check_calls = 0
            end if
            call run()
            stop_check = state%macro_iteration + 1
        end do
        conv_check_converged_from = 1
        if (error /= 0 .or. state%max_precision_reached .or. &
            index(solver_context%log_message, function_unchanged_warning_msg) > 0) then
            write(stderr, *) "test_s_gek_solver failed: Convergence at a point of "// &
                "maximum precision not reported as convergence."
            test_s_gek_solver = .false.
        end if

        ! near the minimum of the same problem, with a convergence threshold that can
        ! not be met, a trial point whose energy is infinite is rejected while the
        ! energy change the step predicts can not be resolved, which stops at maximum
        ! precision without moving to such a point
        call setup_case(diagonal_a, 0.0_rp, 1e8_rp)
        quartic_x = quartic_minimum
        quartic_x(n_quartic) = quartic_x(n_quartic) + 1e-5_rp
        settings%grad_noise = .false.
        settings%conv_tol = 0.0_rp
        n_rejected_trials = huge(1_ip)
        infinite_rejections = .true.
        call run()
        if (error /= 0 .or. .not. state%max_precision_reached .or. moved_to_rejected) &
            then
            write(stderr, *) "test_s_gek_solver failed: Rejected steps not stopped "// &
                "at maximum precision once their prediction can not be resolved."
            test_s_gek_solver = .false.
        end if

        ! a run without the gradient noise, which repeats itself, converges within as
        ! many macro iterations as it evaluated points, the row of its last point plus
        ! one, so that the point the last evaluation reaches is checked, and ends with
        ! one fewer with the error of the solver and an error message
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        settings%grad_noise = .false.
        call run()
        n_evals = state%macro_iteration + 1
        call check_converged("as many macro iterations as evaluated points")
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        settings%grad_noise = .false.
        settings%n_macro = n_evals
        call run()
        if (error /= 0 .or. state%macro_iteration /= n_evals - 1) then
            write(stderr, *) "test_s_gek_solver failed: Point of the last macro "// &
                "iteration not checked or not counted as its row."
            test_s_gek_solver = .false.
        end if
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        settings%grad_noise = .false.
        settings%n_macro = n_evals - 1
        call setup_error_logging(settings, solver_context)
        call run()
        if (error /= error_solver_max_iter .or. &
            len_trim(solver_context%log_message) == 0) then
            write(stderr, *) "test_s_gek_solver failed: Maximum number of "// &
                "iterations not reported with an error message."
            test_s_gek_solver = .false.
        end if

        ! settings the sanity checks refuse are reported with the origin of S-GEK and
        ! an error message
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        call setup_error_logging(settings, solver_context)
        call s_gek_solver(update_orbs_funptr, obj_func_funptr, &
                          change_reference_funptr, 0_ip, error, settings, state, &
                          gek_settings)
        if (error /= error_s_gek + 1 .or. len_trim(solver_context%log_message) == 0) &
            then
            write(stderr, *) "test_s_gek_solver failed: Solver settings refused by "// &
                "the sanity check not reported with the origin of S-GEK and an "// &
                "error message."
            test_s_gek_solver = .false.
        end if
        call setup_case(convex_a, 0.1_rp, 0.0_rp)
        gek_settings%inner_solver = "unknown"
        call setup_error_logging(gek_settings, s_gek_context)
        call run()
        if (error /= error_s_gek + 1 .or. len_trim(s_gek_context%log_message) == 0) then
            write(stderr, *) "test_s_gek_solver failed: S-GEK settings refused by "// &
                "the sanity check not reported with the origin of S-GEK and an "// &
                "error message."
            test_s_gek_solver = .false.
        end if

        ! fault-injection sweep: an unfaulted run from the saddle point of the double
        ! well with the stability check, whose escape calls the objective function, and
        ! the static prior calls every callback function, then every call fails in
        ! turn, with an error and with a value that is not a number, which only rejects
        ! the step for the objective function
        call run_with_fault(0_ip, 0_ip, error_fault)
        if (error /= 0) then
            write(stderr, *) "test_s_gek_solver failed: Produced error without fault."
            test_s_gek_solver = .false.
            return
        end if
        n_calls = n_callback_calls
        do callback = 1, n_quartic_callbacks
            if (n_calls(callback) == 0) then
                write(stderr, *) "test_s_gek_solver failed: Unfaulted run did not "// &
                    "call "//trim(callback_names(callback))//"."
                test_s_gek_solver = .false.
            end if
            do i_call = 1, n_calls(callback)
                do kind = error_fault, nan_fault
                    call run_with_fault(callback, i_call, kind)
                    write(fault_name, '(A, A, A, I0)') trim(callback_names(callback)), &
                        trim(merge(" error at call", " NaN at call  ", &
                                   kind == error_fault)), " ", i_call
                    if (kind == error_fault) then
                        select case (callback)
                        case (update_orbs_callback)
                            expected = error_update_orbs + fault_error_code
                        case (hess_x_callback)
                            expected = error_hess_x + fault_error_code
                        case (change_reference_callback)
                            expected = error_change_reference + fault_error_code
                        case (obj_func_grad_callback)
                            expected = error_obj_func_grad + fault_error_code
                        case (init_hess_callback)
                            expected = error_init_hess + fault_error_code
                        case default
                            expected = error_obj_func + fault_error_code
                        end select
                    else if (callback == obj_func_callback) then
                        expected = 0
                    else if (callback == hess_x_callback) then
                        expected = error_stability_check + 1
                    else
                        expected = error_s_gek + 1
                    end if
                    if (error /= expected) then
                        write(stderr, *) "test_s_gek_solver failed: Failing "// &
                            trim(fault_name)//" not reported with its origin."
                        test_s_gek_solver = .false.
                    end if
                    if (error /= 0 .and. n_calls_after_fault > 0) then
                        write(stderr, *) "test_s_gek_solver failed: Callback "// &
                            "function called after failing "//trim(fault_name)//"."
                        test_s_gek_solver = .false.
                    end if
                    if (state%n_update_orbs /= &
                        n_callback_calls(update_orbs_callback) .or. &
                        state%n_obj_func /= n_callback_calls(obj_func_grad_callback) + &
                        n_callback_calls(obj_func_callback)) then
                        write(stderr, *) "test_s_gek_solver failed: State does not "// &
                            "count the calls for failing "//trim(fault_name)//"."
                        test_s_gek_solver = .false.
                    end if
                end do
            end do
        end do

    contains

        subroutine setup_case(a, d, offset)
            !
            ! this subroutine sets up the quartic problem and settings with the default
            ! options, the chart transformation and the energy and gradient at trial
            ! points, whose callback functions receive distinct contexts and log every
            ! message
            !
            real(rp), intent(in) :: a(:, :), d, offset

            call setup_quartic(a, d, offset, minimum)
            call setup_settings(settings, solver_context)
            settings%conv_tol = 1e-6_rp
            settings%n_macro = 100
            settings%stability = .false.
            call setup_settings(gek_settings, s_gek_context)
            gek_settings%obj_func_grad => quartic_obj_func_grad

        end subroutine setup_case

        subroutine run()
            !
            ! this subroutine runs the S-GEK solver on the quartic problem
            !
            call s_gek_solver(update_orbs_funptr, obj_func_funptr, &
                              change_reference_funptr, n_quartic, error, settings, &
                              state, gek_settings)

        end subroutine run

        subroutine check_converged(name)
            !
            ! this subroutine checks that the solver converged to the minimum of the
            ! convex quartic problem without error and with every callback function
            ! receiving the context of its settings
            !
            character(len=*), intent(in) :: name

            if (error /= 0 .or. &
                .not. all(is_close(quartic_x, quartic_minimum, 1e-5_rp)) .or. &
                .not. is_below(state%grad_rms, settings%conv_tol)) then
                write(stderr, *) "test_s_gek_solver failed: Minimum not reached "// &
                    "with "//name//"."
                test_s_gek_solver = .false.
            end if
            if (context_wrong) then
                write(stderr, *) "test_s_gek_solver failed: Callback function "// &
                    "received another context than that of its settings with "//name// &
                    "."
                test_s_gek_solver = .false.
            end if

        end subroutine check_converged

        subroutine run_with_fault(fault_callback_in, fault_call_in, kind_in)
            !
            ! this subroutine runs the solver from the saddle point of the double well
            ! without the gradient noise, with the stability check and the static prior
            ! to a looser convergence threshold and lets the given call of the given
            ! callback function fail with a fault of the given kind
            !
            integer(ip), intent(in) :: fault_callback_in, fault_call_in, kind_in

            call setup_case(saddle_a, 1.0_rp, 0.0_rp)
            quartic_x(1) = quartic_minimum(1)
            settings%grad_noise = .false.
            settings%conv_tol = 1e-4_rp
            settings%stability = .true.
            settings%verbose = 0
            gek_settings%init_hess => quartic_init_hess
            fault_callback = fault_callback_in
            fault_call = fault_call_in
            fault_kind = kind_in
            call run()

        end subroutine run_with_fault

    end function test_s_gek_solver

    logical(c_bool) function test_s_gek_sanity_check() bind(C)
        !
        ! this function tests the subroutine which checks the S-GEK settings
        !
        use otr_s_gek, only: s_gek_settings_type, s_gek_sanity_check
        use opentrustregion_unit_tests, only: test_context_type, setup_error_logging

        character(len=*), parameter :: option_names(3) = &
            [character(len=12) :: "inner_solver", "prior_update", "acceptance"]
        type(s_gek_settings_type) :: settings
        type(test_context_type), target :: context
        integer(ip) :: error, option

        ! assume tests pass
        test_s_gek_sanity_check = .true.

        ! options in upper case are accepted and converted to lower case
        call setup_valid_settings()
        settings%inner_solver = "RS-RFO"
        settings%prior_update = "NONE"
        settings%acceptance = "Trial"
        call s_gek_sanity_check(settings, error)
        if (error /= 0 .or. settings%inner_solver /= "rs-rfo" .or. &
            settings%prior_update /= "none" .or. settings%acceptance /= "trial") then
            write(stderr, *) "test_s_gek_sanity_check failed: Valid options not "// &
                "accepted in lower case."
            test_s_gek_sanity_check = .false.
        end if

        ! settings with an unknown value of any option are refused with an error message
        do option = 1, 3
            call setup_valid_settings()
            select case (option)
            case (1)
                settings%inner_solver = "unknown"
            case (2)
                settings%prior_update = "unknown"
            case (3)
                settings%acceptance = "unknown"
            end select
            call s_gek_sanity_check(settings, error)
            if (error /= 1) then
                write(stderr, *) "test_s_gek_sanity_check failed: Settings with an "// &
                    "unknown value of "//trim(option_names(option))//" not refused."
                test_s_gek_sanity_check = .false.
            end if
            if (len_trim(context%log_message) == 0) then
                write(stderr, *) "test_s_gek_sanity_check failed: No error message "// &
                    "printed for settings with an unknown value of "// &
                    trim(option_names(option))//"."
                test_s_gek_sanity_check = .false.
            end if
        end do

    contains

        subroutine setup_valid_settings()
            !
            ! this subroutine sets up valid settings, the default settings without
            ! callback functions, which log only error messages
            !
            settings = s_gek_settings_type()
            call setup_error_logging(settings, context)

        end subroutine setup_valid_settings

    end function test_s_gek_sanity_check

    logical(c_bool) function test_invoke_change_reference() bind(C)
        !
        ! this function tests the subroutine which calls the host's chart transformation
        !
        use otr_s_gek, only: s_gek_settings_type, invoke_change_reference, &
                             error_s_gek, error_change_reference
        use opentrustregion_unit_tests, only: setup_error_logging

        integer(ip), parameter :: n_hist = 3
        type(s_gek_settings_type) :: settings
        procedure(change_reference_type), pointer :: change_reference_funptr
        real(rp) :: new_ref(n_full), kappa(n_full, n_hist), &
                    local_grad(n_full, n_hist), grad(n_full, n_hist), &
                    expected(n_full, n_hist), minimum(n_full)
        integer(ip) :: error, k

        ! assume tests pass
        test_invoke_change_reference = .true.

        ! transformation of a random history in the flat chart of the quartic problem
        call random_number(minimum)
        call setup_quartic(generate_random_spd_matrix(n_full), 0.0_rp, 0.0_rp, minimum)
        change_reference_funptr => quartic_change_reference
        settings%context => s_gek_context
        call random_number(new_ref)
        call random_number(kappa)
        call random_number(local_grad)
        call random_number(grad)
        do k = 1, n_hist
            expected(:, k) = kappa(:, k) - new_ref
        end do
        call invoke_change_reference(change_reference_funptr, new_ref, kappa, &
                                     local_grad, grad, settings, error)
        if (error /= 0 .or. .not. all(is_close(kappa, expected)) .or. &
            .not. all(is_close(grad, local_grad))) then
            write(stderr, *) "test_invoke_change_reference failed: Transformed "// &
                "history not returned."
            test_invoke_change_reference = .false.
        end if
        if (context_wrong .or. n_callback_calls(change_reference_callback) /= 1) then
            write(stderr, *) "test_invoke_change_reference failed: Chart "// &
                "transformation not called once with the context of the settings."
            test_invoke_change_reference = .false.
        end if

        ! an error of the transformation is passed on with its origin
        fault_callback = change_reference_callback
        fault_call = 2
        call invoke_change_reference(change_reference_funptr, new_ref, kappa, &
                                     local_grad, grad, settings, error)
        if (error /= error_change_reference + fault_error_code) then
            write(stderr, *) "test_invoke_change_reference failed: Error of the "// &
                "chart transformation not passed on with its origin."
            test_invoke_change_reference = .false.
        end if

        ! values that are not finite in any of the lists are reported with the origin
        ! of S-GEK and an error message
        fault_kind = nan_fault
        do k = 1, 3
            call setup_error_logging(settings, s_gek_context)
            fault_call = n_callback_calls(change_reference_callback) + 1
            nan_list = k
            call invoke_change_reference(change_reference_funptr, new_ref, kappa, &
                                         local_grad, grad, settings, error)
            if (error /= error_s_gek + 1) then
                write(stderr, '(A, I0, A)') " test_invoke_change_reference failed: "// &
                    "Values that are not finite in list ", k, &
                    " not reported with the origin of S-GEK."
                test_invoke_change_reference = .false.
            end if
            if (len_trim(s_gek_context%log_message) == 0) then
                write(stderr, '(A, I0, A)') " test_invoke_change_reference failed: "// &
                    "No error message printed for values that are not finite in "// &
                    "list ", k, "."
                test_invoke_change_reference = .false.
            end if
            kappa = 0.0_rp
            local_grad = 0.0_rp
        end do
        nan_list = 1

    end function test_invoke_change_reference

    logical(c_bool) function test_floored_h_diag() bind(C)
        !
        ! this function tests the function which floors a diagonal Hessian
        !
        use otr_s_gek, only: floored_h_diag

        ! assume tests pass
        test_floored_h_diag = .true.

        ! negative elements of magnitude above and below one, small positive elements
        ! below the floor, at the floor and above it
        if (.not. all(is_close( &
            floored_h_diag([-3.0_rp, -0.5_rp, 0.01_rp, 0.05_rp, 0.2_rp]), &
            [3.0_rp, 1.0_rp, 0.05_rp, 0.05_rp, 0.2_rp]))) then
            write(stderr, *) "test_floored_h_diag failed: Wrong floored diagonal."
            test_floored_h_diag = .false.
        end if

    end function test_floored_h_diag

    logical(c_bool) function test_variance_threshold_of() bind(C)
        !
        ! this function tests the function which returns the variance threshold of an
        ! energy change: its magnitude between the floor and the cap, the floor below
        ! it and the cap, scaled by the bound factor, above it
        !
        use otr_s_gek, only: variance_threshold_of, variance_threshold_floor, &
                             variance_threshold_cap

        ! assume tests pass
        test_variance_threshold_of = .true.

        ! the magnitude of an energy change, rising or falling, between floor and cap
        if (.not. is_close(variance_threshold_of(-0.02_rp, 1.0_rp), 0.02_rp) .or. &
            .not. is_close(variance_threshold_of(0.02_rp, 1.0_rp), 0.02_rp)) then
            write(stderr, *) "test_variance_threshold_of failed: Energy change "// &
                "between floor and cap not returned."
            test_variance_threshold_of = .false.
        end if

        ! the floor below it
        if (.not. is_close(variance_threshold_of(1e-4_rp, 1.0_rp), &
                           variance_threshold_floor)) then
            write(stderr, *) "test_variance_threshold_of failed: Floor not "// &
                "returned for a small energy change."
            test_variance_threshold_of = .false.
        end if

        ! the cap above it, which the bound factor scales
        if (.not. is_close(variance_threshold_of(10.0_rp, 1.0_rp), &
                           variance_threshold_cap) .or. &
            .not. is_close(variance_threshold_of(10.0_rp, 20.0_rp), &
                           20.0_rp * variance_threshold_cap) .or. &
            .not. is_close(variance_threshold_of(0.02_rp, 1e4_rp), 0.02_rp)) then
            write(stderr, *) "test_variance_threshold_of failed: Cap scaled by the "// &
                "bound factor not returned for a large energy change."
            test_variance_threshold_of = .false.
        end if

    end function test_variance_threshold_of

    logical(c_bool) function test_updated_threshold_factor() bind(C)
        !
        ! this function tests the function which returns the factor of the variance
        ! threshold after a step: one once the step is taken, whether the trial point
        ! was kept or not, half the factor after a rejected trial point that leaves the
        ! history unchanged and the factor unchanged after one that is kept
        !
        use otr_s_gek, only: updated_threshold_factor
        use opentrustregion_core, only: trust_radius_shrink_factor

        ! assume tests pass
        test_updated_threshold_factor = .true.

        ! a step taken resets the factor
        if (.not. is_close(updated_threshold_factor(0.25_rp, .true., .true.), &
                           1.0_rp) .or. .not. is_close(updated_threshold_factor( &
                               0.25_rp, .true., .false.), 1.0_rp)) then
            write(stderr, *) "test_updated_threshold_factor failed: Factor not "// &
                "reset once a step is taken."
            test_updated_threshold_factor = .false.
        end if

        ! a rejected trial point that leaves the history unchanged shrinks the factor
        if (.not. is_close(updated_threshold_factor(0.25_rp, .false., .false.), &
                           0.25_rp * trust_radius_shrink_factor)) then
            write(stderr, *) "test_updated_threshold_factor failed: Factor not "// &
                "shrunk after a rejected trial point that leaves the history unchanged."
            test_updated_threshold_factor = .false.
        end if

        ! a rejected trial point kept as data leaves it unchanged
        if (.not. is_close(updated_threshold_factor(0.25_rp, .false., .true.), &
                           0.25_rp)) then
            write(stderr, *) "test_updated_threshold_factor failed: Factor changed "// &
                "after a rejected trial point kept as data."
            test_updated_threshold_factor = .false.
        end if

    end function test_updated_threshold_factor

    logical(c_bool) function test_inverse_bfgs_step() bind(C)
        !
        ! this function tests the function which returns the step of the inverse BFGS
        ! update of a diagonal Hessian
        !
        use otr_s_gek, only: inverse_bfgs_step
        use otr_common, only: bfgs_skip_thresh
        use opentrustregion_unit_tests, only: identity_matrix

        integer(ip), parameter :: n_hist = 5
        real(rp) :: kappa_list(n_dim, n_hist), grad_list(n_dim, n_hist), &
                    h_diag(n_dim), s(n_dim), t(n_dim), expected(n_dim), &
                    hinv(n_dim, n_dim), y(n_dim), rho, left(n_dim, n_dim)
        integer(ip) :: i, k

        ! assume tests pass
        test_inverse_bfgs_step = .true.

        ! history of four pairs, the first and the last with positive curvature, the
        ! second orthogonal and the third with negative curvature, both of which are
        ! skipped, with a random positive diagonal, compared relative to the step,
        ! whose rounding the inverse of a small product s^T y amplifies
        call random_number(kappa_list)
        call random_number(grad_list)
        call random_number(t)
        s = kappa_list(:, 2) - kappa_list(:, 1)
        grad_list(:, 2) = grad_list(:, 1) + 2.0_rp * s + 0.1_rp * t
        s = kappa_list(:, 3) - kappa_list(:, 2)
        grad_list(:, 3) = grad_list(:, 2) + t - &
                          dot_product(t, s) / dot_product(s, s) * s
        s = kappa_list(:, 4) - kappa_list(:, 3)
        grad_list(:, 4) = grad_list(:, 3) - s
        s = kappa_list(:, 5) - kappa_list(:, 4)
        grad_list(:, 5) = grad_list(:, 4) + 3.0_rp * s + 0.1_rp * t
        call random_number(h_diag)
        h_diag = h_diag + 0.05_rp

        ! step of the dense inverse BFGS updates of the inverse diagonal over the pairs
        ! with sufficiently positive curvature
        hinv = 0.0_rp
        do i = 1, n_dim
            hinv(i, i) = 1.0_rp / h_diag(i)
        end do
        do k = 1, n_hist - 1
            s = kappa_list(:, k + 1) - kappa_list(:, k)
            y = grad_list(:, k + 1) - grad_list(:, k)
            if (dot_product(s, y) <= bfgs_skip_thresh * norm2(s) * norm2(y)) cycle
            rho = 1.0_rp / dot_product(s, y)
            left = identity_matrix(n_dim) - &
                   rho * spread(s, 2, n_dim) * spread(y, 1, n_dim)
            hinv = matmul(left, matmul(hinv, transpose(left))) + &
                   rho * spread(s, 2, n_dim) * spread(s, 1, n_dim)
        end do
        expected = -matmul(hinv, grad_list(:, n_hist))
        if (.not. all(is_close(inverse_bfgs_step(kappa_list, grad_list, h_diag), &
                               expected, &
                               1e-8_rp * max(1.0_rp, maxval(abs(expected)))))) then
            write(stderr, *) "test_inverse_bfgs_step failed: Wrong step."
            test_inverse_bfgs_step = .false.
        end if

        ! a single point gives the step of the inverse diagonal
        if (.not. all(is_close( &
            inverse_bfgs_step(kappa_list(:, 1:1), grad_list(:, 1:1), h_diag), &
            -grad_list(:, 1) / h_diag))) then
            write(stderr, *) "test_inverse_bfgs_step failed: Wrong step without pairs."
            test_inverse_bfgs_step = .false.
        end if

    end function test_inverse_bfgs_step

    logical(c_bool) function test_s_gek_step() bind(C)
        !
        ! this function tests the subroutine which computes the S-GEK step of a window
        ! of points: the step is the minimum of the GEK of the data projected onto the
        ! subspace of the differences, the current gradient, the proposal and the
        ! perturbation with the projected prior, the floored diagonal or the host's
        ! static Hessian, updated by SR1 or BFGS if asked for and stretched along the
        ! proposal by the stretch factor, within the given variance threshold, found by
        ! the chosen inner solver
        !
        use otr_s_gek, only: s_gek_settings_type, s_gek_step, gek_type, error_s_gek
        use otr_common, only: error_init_hess
        use otr_common_unit_tests, only: ref_sr1_update, ref_bfgs_update
        use opentrustregion_unit_tests, only: setup_error_logging
        use opentrustregion_core, only: verbosity_warning, verbosity_error
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        type(s_gek_settings_type) :: settings
        type(gek_type) :: gek
        real(rp) :: kappa_list(n_full, n_points), grad_list(n_full, n_points), &
                    func_list(n_points), h_diag(n_full), proposal(n_full), &
                    perturbation(n_full), step(n_full), minimum(n_full), &
                    pred_func_diff, stretch_factor, dispersion
        real(rp), allocatable :: basis(:, :)
        integer(ip) :: error, i_prior, i_update, i_stretch, i_inner, k
        logical :: restricted, converged
        character(len=60) :: case_name
        character(len=4), parameter :: updates(3) = &
            [character(len=4) :: "none", "sr1", "bfgs"]

        ! assume tests pass
        test_s_gek_step = .true.

        ! window of three points of a convex quadratic, the current point last at the
        ! origin, its diagonal and a proposal step
        call random_number(minimum)
        call setup_quartic(generate_random_spd_matrix(n_full), 0.0_rp, 0.0_rp, &
                           minimum - 0.5_rp)
        call random_number(kappa_list)
        kappa_list = kappa_list - 0.5_rp
        kappa_list(:, n_points) = 0.0_rp
        do k = 1, n_points
            grad_list(:, k) = quartic_grad(kappa_list(:, k))
            func_list(k) = quartic_func(kappa_list(:, k))
        end do
        do k = 1, n_full
            h_diag(k) = quartic_a(k, k)
        end do
        proposal = -grad_list(:, n_points) / h_diag
        perturbation = 0.0_rp
        basis = ref_basis()
        settings%context => s_gek_context

        ! every prior without, with SR1 and with BFGS update and with and without
        ! stretch factor: the unrestricted step is a stationary point of the GEK of the
        ! projected data, lies in the subspace and comes with the energy the GEK
        ! predicts for it
        do i_prior = 1, 2
            do i_update = 1, 3
                do i_stretch = 1, 2
                    settings%init_hess => null()
                    if (i_prior == 2) settings%init_hess => quartic_init_hess
                    settings%prior_update = updates(i_update)
                    stretch_factor = merge(1.0_rp, 1.6_rp, i_stretch == 1)
                    write(case_name, '(A, A, A, A, F3.1)') &
                        trim(merge("diagonal prior", "static prior  ", i_prior == 1)), &
                        ", update ", trim(settings%prior_update), ", stretch factor ", &
                        stretch_factor
                    n_callback_calls = 0
                    call s_gek_step(settings, kappa_list, grad_list, func_list, &
                                    h_diag, proposal, perturbation, stretch_factor, &
                                    1.0_rp, 1e10_rp, step, restricted, converged, &
                                    pred_func_diff, error)
                    if (error /= 0 .or. .not. converged) then
                        write(stderr, *) "test_s_gek_step failed: Unrestricted "// &
                            "step produced error or was restricted for "// &
                            trim(case_name)//"."
                        test_s_gek_step = .false.
                        cycle
                    end if
                    call ref_window_gek(i_prior == 2, updates(i_update), stretch_factor)
                    call check_step(trim(case_name))
                    if (i_prior == 2 .and. &
                        (context_wrong .or. &
                         n_callback_calls(init_hess_callback) /= size(basis, 2))) then
                        write(stderr, *) "test_s_gek_step failed: Static Hessian "// &
                            "not applied once to every basis vector with the "// &
                            "context of the settings for "//trim(case_name)//"."
                        test_s_gek_step = .false.
                    end if
                end do
            end do
        end do

        ! a perturbation of the gradient joins the subspace and perturbs the gradient
        ! of the surrogate, with the static Hessian as prior, which keeps it convex,
        ! for either inner solver: the step is a stationary point of the reference GEK
        ! for the gradient perturbed by the projected perturbation in the subspace that
        ! includes it and comes with the energy change of the unperturbed GEK
        settings%init_hess => quartic_init_hess
        settings%prior_update = "sr1"
        call random_number(perturbation)
        perturbation = 0.1_rp * norm2(grad_list(:, n_points)) * (perturbation - 0.5_rp)
        basis = ref_basis()
        call ref_window_gek(.true., "sr1", 1.0_rp)
        do i_inner = 1, 2
            settings%inner_solver = trim(merge("otr   ", "rs-rfo", i_inner == 1))
            call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, &
                            proposal, perturbation, 1.0_rp, 1.0_rp, 1e10_rp, step, &
                            restricted, converged, pred_func_diff, error)
            if (error /= 0 .or. .not. converged) then
                write(stderr, *) "test_s_gek_step failed: Unrestricted step with a "// &
                    "perturbation produced error or was restricted with "// &
                    trim(settings%inner_solver)//"."
                test_s_gek_step = .false.
            else
                call check_step("perturbation with "//trim(settings%inner_solver))
            end if
        end do
        perturbation = 0.0_rp
        basis = ref_basis()

        ! with the static Hessian as prior, which keeps the surrogate convex, the micro
        ! iterations of restricted-step RFO end at the same stationary point, and with
        ! a step bound too small to reach it without convergence, which the inner
        ! OpenTrustRegion solve does not depend on
        settings%init_hess => quartic_init_hess
        settings%prior_update = "sr1"
        call ref_window_gek(.true., "sr1", 1.0_rp)
        settings%inner_solver = "rs-rfo"
        call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, proposal, &
                        perturbation, 1.0_rp, 1.0_rp, 1e10_rp, step, restricted, &
                        converged, pred_func_diff, error)
        if (error /= 0 .or. .not. converged) then
            write(stderr, *) "test_s_gek_step failed: Unrestricted step of the "// &
                "micro iterations produced error or was restricted."
            test_s_gek_step = .false.
        else
            call check_step("micro iterations")
        end if
        call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, proposal, &
                        perturbation, 1.0_rp, 1e-3_rp, 1e10_rp, step, restricted, &
                        converged, pred_func_diff, error)
        if (error /= 0 .or. converged .or. restricted) then
            write(stderr, *) "test_s_gek_step failed: Micro iterations not selected."
            test_s_gek_step = .false.
        end if
        settings%inner_solver = "otr"
        call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, proposal, &
                        perturbation, 1.0_rp, 1e-3_rp, 1e10_rp, step, restricted, &
                        converged, pred_func_diff, error)
        if (error /= 0 .or. .not. converged) then
            write(stderr, *) "test_s_gek_step failed: Inner OpenTrustRegion solve "// &
                "not selected."
            test_s_gek_step = .false.
        end if

        ! a variance threshold below the excess dispersion of the unrestricted step
        ! restricts the step to within it, where it still lowers the energy of the GEK,
        ! with either inner solver
        gek%current_dispersion = ref_dispersion_gek(gek, &
                                                    [(0.0_rp, k=1, size(basis, 2))])
        gek%excess_dispersion = .true.
        dispersion = ref_dispersion_gek(gek, matmul(transpose(basis), step))
        do i_inner = 1, 2
            settings%inner_solver = trim(merge("otr   ", "rs-rfo", i_inner == 1))
            call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, &
                            proposal, perturbation, 1.0_rp, 1.0_rp, &
                            0.5_rp * dispersion, step, restricted, converged, &
                            pred_func_diff, error)
            if (error /= 0 .or. .not. restricted .or. converged .or. .not. is_below( &
                ref_dispersion_gek(gek, matmul(transpose(basis), step)), &
                0.5_rp * dispersion) .or. .not. is_below(pred_func_diff, 0.0_rp)) then
                write(stderr, *) "test_s_gek_step failed: Step not restricted to "// &
                    "within the variance threshold with the "// &
                    trim(settings%inner_solver)//" inner solver."
                test_s_gek_step = .false.
            end if
        end do
        settings%inner_solver = "otr"

        ! an error of the static Hessian is passed on with its origin and values that
        ! are not finite are reported with the origin of S-GEK and an error message
        settings%init_hess => quartic_init_hess
        n_callback_calls = 0
        fault_callback = init_hess_callback
        fault_call = 2
        fault_kind = error_fault
        call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, proposal, &
                        perturbation, 1.0_rp, 1.0_rp, 1e10_rp, step, restricted, &
                        converged, pred_func_diff, error)
        if (error /= error_init_hess + fault_error_code .or. n_calls_after_fault > 0) &
            then
            write(stderr, *) "test_s_gek_step failed: Error of the static Hessian "// &
                "not passed on with its origin at once."
            test_s_gek_step = .false.
        end if
        call setup_error_logging(settings, s_gek_context)
        n_callback_calls = 0
        fault_kind = nan_fault
        call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, proposal, &
                        perturbation, 1.0_rp, 1.0_rp, 1e10_rp, step, restricted, &
                        converged, pred_func_diff, error)
        if (error /= error_s_gek + 1 .or. len_trim(s_gek_context%log_message) == 0) then
            write(stderr, *) "test_s_gek_step failed: Static Hessian with values "// &
                "that are not finite not reported with the origin of S-GEK and an "// &
                "error message."
            test_s_gek_step = .false.
        end if
        fault_callback = 0
        settings%init_hess => null()

        ! a proposal along the current gradient is dropped from the subspace without a
        ! warning, since its vectors are routinely linearly dependent
        s_gek_context%log_message = ""
        settings%verbose = verbosity_warning
        call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, &
                        grad_list(:, n_points), perturbation, 1.0_rp, 1.0_rp, 1e10_rp, &
                        step, restricted, converged, pred_func_diff, error)
        if (error /= 0 .or. len(s_gek_context%log_message) /= 0) then
            write(stderr, *) "test_s_gek_step failed: Linearly dependent proposal "// &
                "produced an error or a message."
            test_s_gek_step = .false.
        end if
        settings%verbose = verbosity_error

        ! a GEK that cannot be set up, from a diagonal that is not finite and enters
        ! only the prior, is reported with the origin of S-GEK and an error message
        s_gek_context%log_message = ""
        h_diag(1) = ieee_value(1.0_rp, ieee_quiet_nan)
        call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, proposal, &
                        perturbation, 1.0_rp, 1.0_rp, 1e10_rp, step, restricted, &
                        converged, pred_func_diff, error)
        if (error /= error_s_gek + 1 .or. len_trim(s_gek_context%log_message) == 0) then
            write(stderr, *) "test_s_gek_step failed: GEK that cannot be set up "// &
                "not reported with the origin of S-GEK and an error message."
            test_s_gek_step = .false.
        end if
        h_diag(1) = quartic_a(1, 1)

        ! a subspace that cannot be orthonormalized, from a position that is not
        ! finite, is reported with the origin of S-GEK and an error message
        s_gek_context%log_message = ""
        kappa_list(1, 1) = ieee_value(1.0_rp, ieee_quiet_nan)
        call s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, proposal, &
                        perturbation, 1.0_rp, 1.0_rp, 1e10_rp, step, restricted, &
                        converged, pred_func_diff, error)
        if (error /= error_s_gek + 1 .or. len_trim(s_gek_context%log_message) == 0) then
            write(stderr, *) "test_s_gek_step failed: Subspace that cannot be "// &
                "orthonormalized not reported with the origin of S-GEK and an "// &
                "error message."
            test_s_gek_step = .false.
        end if

    contains

        function ref_basis() result(columns)
            !
            ! this function returns an orthonormal basis of the space of the gradient
            ! and position differences, the current gradient, the proposal and the
            ! perturbation, the eigenvectors of their projector with a non-vanishing
            ! eigenvalue
            !
            real(rp), allocatable :: columns(:, :)

            real(rp) :: vectors(n_full, 2 * n_points + 1), projector(n_full, n_full), &
                        eigvals(n_full), work(10 * n_full)
            integer(ip) :: info, j

            do j = 1, n_points - 1
                vectors(:, 2 * j - 1) = grad_list(:, j + 1) - grad_list(:, j)
                vectors(:, 2 * j) = kappa_list(:, j + 1) - kappa_list(:, j)
            end do
            vectors(:, 2 * n_points - 1) = grad_list(:, n_points)
            vectors(:, 2 * n_points) = proposal
            vectors(:, 2 * n_points + 1) = perturbation
            projector = matmul(vectors, transpose(vectors))
            call dsyev("V", "U", n_full, projector, n_full, eigvals, work, &
                       size(work, kind=ip), info)
            columns = projector(:, count(eigvals <= 1e-10_rp * maxval(eigvals)) + 1:)

        end function ref_basis

        subroutine ref_window_gek(static, update, factor)
            !
            ! this subroutine sets up the GEK of the window projected onto the basis
            ! with the projected prior, the diagonal or the static Hessian, updated by
            ! SR1 or BFGS if asked for and stretched along the proposal by the stretch
            ! factor
            !
            logical, intent(in) :: static
            character(len=*), intent(in) :: update
            real(rp), intent(in) :: factor

            real(rp), allocatable :: q_s(:, :), g_s(:, :), prior(:, :), w(:), d(:, :), &
                                     steps(:, :), grad_diffs(:, :)
            integer(ip) :: m, j

            m = size(basis, 2)
            q_s = matmul(transpose(basis), &
                         kappa_list - spread(kappa_list(:, n_points), 2, n_points))
            g_s = matmul(transpose(basis), grad_list)
            if (static) then
                prior = matmul(transpose(basis), matmul(quartic_a, basis))
            else
                prior = matmul(transpose(basis), spread(h_diag, 2, m) * basis)
            end if
            steps = q_s(:, 2:) - q_s(:, :n_points - 1)
            grad_diffs = g_s(:, 2:) - g_s(:, :n_points - 1)
            if (update == "sr1") then
                prior = ref_sr1_update(prior, steps, grad_diffs)
            else if (update == "bfgs") then
                prior = ref_bfgs_update(prior, steps, grad_diffs)
            end if
            w = matmul(transpose(basis), proposal) / norm2(proposal)
            d = (1.0_rp / factor - 1.0_rp) * spread(w, 2, m) * spread(w, 1, m)
            do j = 1, m
                d(j, j) = d(j, j) + 1.0_rp
            end do
            call ref_setup_gek(gek, q_s, g_s, func_list, matmul(d, matmul(prior, d)))

        end subroutine ref_window_gek

        subroutine check_step(name)
            !
            ! this subroutine checks that the step lies in the subspace, is a
            ! stationary point of the reference GEK for the gradient perturbed by the
            ! projected perturbation and comes with the energy change the unperturbed
            ! GEK predicts
            !
            character(len=*), intent(in) :: name

            real(rp), allocatable :: sub_kappa(:)
            real(rp) :: expected

            sub_kappa = matmul(transpose(basis), step)
            if (.not. all(is_close(matmul(basis, sub_kappa), step, 1e-10_rp))) then
                write(stderr, *) "test_s_gek_step failed: Step does not lie in the "// &
                    "subspace for "//name//"."
                test_s_gek_step = .false.
            end if
            if (.not. is_below(norm2(ref_gradient_gek(gek, sub_kappa) + matmul( &
                transpose(basis), perturbation)), 1e-6_rp)) then
                write(stderr, *) "test_s_gek_step failed: Step is not a stationary "// &
                    "point of the GEK for "//name//"."
                test_s_gek_step = .false.
            end if
            ! the predicted change, relative to which the prior, updated in the basis
            ! of the production and of the reference, rounds differently
            expected = ref_energy_gek(gek, sub_kappa) - &
                       ref_energy_gek(gek, 0.0_rp * sub_kappa)
            if (.not. is_close(pred_func_diff, expected, max( &
                roundoff_tolerance(gek), 1e-8_rp * abs(expected)))) then
                write(stderr, *) "test_s_gek_step failed: Wrong predicted energy "// &
                    "change for "//name//"."
                test_s_gek_step = .false.
            end if

        end subroutine check_step

    end function test_s_gek_step

    logical(c_bool) function test_updated_prior() bind(C)
        !
        ! this function tests the function which returns the prior updated by the SR1
        ! or BFGS updates over consecutive pairs of positions and gradients
        !
        use otr_s_gek, only: updated_prior
        use otr_common_unit_tests, only: ref_sr1_update, ref_bfgs_update
        use opentrustregion_unit_tests, only: generate_random_symm_matrix

        real(rp) :: hess(n_dim, n_dim), asymm(n_dim, n_dim), &
                    kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    updated(n_dim, n_dim), expected(n_dim, n_dim), s(n_dim), t(n_dim)

        ! assume tests pass
        test_updated_prior = .true.

        ! SR1 update over random pairs, which satisfies the secant condition of the
        ! last pair, of a prior given with an antisymmetric part, compared relative to
        ! the update, whose rounding a small denominator amplifies
        hess = generate_random_symm_matrix(n_dim)
        call random_number(asymm)
        asymm = asymm - transpose(asymm)
        call random_number(kappa_list)
        call random_number(grad_list)
        updated = updated_prior(hess + asymm, kappa_list, grad_list, "sr1")
        expected = ref_sr1_update(hess, &
                                  kappa_list(:, 2:) - kappa_list(:, :n_points - 1), &
                                  grad_list(:, 2:) - grad_list(:, :n_points - 1))
        if (.not. all(is_close(updated, expected, &
                               1e-10_rp * max(1.0_rp, maxval(abs(expected)))))) then
            write(stderr, *) "test_updated_prior failed: Wrong SR1 update."
            test_updated_prior = .false.
        end if
        if (.not. all(is_close(matmul(updated, kappa_list(:, 3) - kappa_list(:, 2)), &
                               grad_list(:, 3) - grad_list(:, 2), &
                               1e-10_rp * max(1.0_rp, maxval(abs(updated)))))) then
            write(stderr, *) "test_updated_prior failed: SR1 update does not "// &
                "satisfy the secant condition of the last pair."
            test_updated_prior = .false.
        end if

        ! BFGS update of a positive definite prior given with an antisymmetric part
        ! over a pair of negative curvature, which is skipped, and a pair of positive
        ! curvature, whose secant condition the update satisfies
        hess = generate_random_spd_matrix(n_dim)
        call random_number(kappa_list)
        call random_number(grad_list)
        call random_number(t)
        s = kappa_list(:, 2) - kappa_list(:, 1)
        grad_list(:, 2) = grad_list(:, 1) - s
        s = kappa_list(:, 3) - kappa_list(:, 2)
        grad_list(:, 3) = grad_list(:, 2) + matmul(hess, s) + 0.1_rp * t
        updated = updated_prior(hess + asymm, kappa_list, grad_list, "bfgs")
        if (.not. all(is_close(updated, ref_bfgs_update( &
            hess, kappa_list(:, 2:) - kappa_list(:, :n_points - 1), &
            grad_list(:, 2:) - grad_list(:, :n_points - 1)), 1e-10_rp))) then
            write(stderr, *) "test_updated_prior failed: Wrong BFGS update."
            test_updated_prior = .false.
        end if
        if (.not. all(is_close(matmul(updated, s), grad_list(:, 3) - grad_list(:, 2), &
                               1e-10_rp))) then
            write(stderr, *) "test_updated_prior failed: BFGS update does not "// &
                "satisfy the secant condition of the last pair."
            test_updated_prior = .false.
        end if

    end function test_updated_prior

    logical(c_bool) function test_setup_gek() bind(C)
        !
        ! this function tests the subroutine which sets up the GEK
        !
        use otr_s_gek, only: gek_type, trend_bias, prior_eigval_floor, confidence_factor
        use opentrustregion_unit_tests, only: identity_matrix
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        integer(ip), parameter :: m = n_points * (n_dim + 1)
        type(gek_type) :: gek
        real(rp) :: kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    func_list(n_points), hess(n_dim, n_dim), asymm(n_dim, n_dim), &
                    diag(n_dim, n_dim), upper(m, m), cov(m, m), b(m), &
                    gs(n_dim, n_points)
        integer(ip) :: error, i

        ! assume tests pass
        test_setup_gek = .true.

        ! set up the GEK of random data with a prior with a negative eigenvalue and one
        ! below the floor of the length scales, given with an antisymmetric part
        call random_number(kappa_list)
        kappa_list = 2.0_rp * kappa_list - 1.0_rp
        call random_number(grad_list)
        grad_list = 2.0_rp * grad_list - 1.0_rp
        call random_number(func_list)
        hess = generate_random_prior(prior_eigenvalues)
        call random_number(asymm)
        asymm = asymm - transpose(asymm)
        call gek%setup(kappa_list, grad_list, func_list, hess + asymm, .false., error)
        if (error /= 0) then
            write(stderr, *) "test_setup_gek failed: Produced error."
            test_setup_gek = .false.
            return
        end if
        if (gek%n_dim /= n_dim .or. gek%n_points /= n_points) then
            write(stderr, *) "test_setup_gek failed: Wrong dimensions."
            test_setup_gek = .false.
        end if

        ! eigenbasis of the symmetric part of the prior, length scales from its
        ! eigenvalues, floored, and positions in the eigenbasis
        diag = matmul(transpose(gek%eigvecs), matmul(hess, gek%eigvecs))
        if (.not. all(is_close(matmul(transpose(gek%eigvecs), gek%eigvecs), &
                               identity_matrix(n_dim))) .or. &
            .not. all(is_close(diag - diag * identity_matrix(n_dim), 0.0_rp))) then
            write(stderr, *) "test_setup_gek failed: Wrong eigenbasis of the prior."
            test_setup_gek = .false.
        end if
        do i = 1, n_dim
            if (.not. is_close(gek%length_scales(i), &
                               sqrt(5.0_rp / 3.0_rp * trend_bias / &
                                    max(abs(diag(i, i)), prior_eigval_floor)))) then
                write(stderr, *) "test_setup_gek failed: Wrong length scale."
                test_setup_gek = .false.
            end if
        end do
        if (.not. all(is_close(gek%kappa_list, &
                               matmul(transpose(gek%eigvecs), kappa_list)))) then
            write(stderr, *) "test_setup_gek failed: Wrong positions in the eigenbasis."
            test_setup_gek = .false.
        end if

        ! Cholesky factor of the regularized covariance matrix, trend, solution vector
        ! and process variance
        upper = 0.0_rp
        do i = 1, m
            upper(:i, i) = gek%cov_chol(:i, i)
        end do
        cov = ref_covariance_gek(gek)
        if (.not. all(is_close(matmul(transpose(upper), upper), cov, 1e-12_rp))) then
            write(stderr, *) "test_setup_gek failed: Wrong Cholesky factor of the "// &
                "covariance matrix."
            test_setup_gek = .false.
        end if
        if (.not. is_close(gek%trend, maxval(func_list) + trend_bias)) then
            write(stderr, *) "test_setup_gek failed: Wrong trend."
            test_setup_gek = .false.
        end if
        gs = matmul(transpose(gek%eigvecs), grad_list)
        b(:n_points) = func_list - maxval(func_list) - trend_bias
        do i = 1, n_dim
            b(n_points + (i - 1) * n_points + 1:n_points + i * n_points) = gs(i, :)
        end do
        if (.not. all(is_close(matmul(cov, gek%weights), b, &
                               1e-12_rp * sum(abs(gek%weights))))) then
            write(stderr, *) "test_setup_gek failed: Wrong solution vector."
            test_setup_gek = .false.
        end if
        if (.not. is_close(gek%process_variance, dot_product(b, gek%weights) / &
                           real(m, kind=rp), 1e-12_rp * sum(abs(b * gek%weights)))) then
            write(stderr, *) "test_setup_gek failed: Wrong process variance."
            test_setup_gek = .false.
        end if
        if (gek%excess_dispersion .or. &
            .not. is_at_most(abs(gek%current_dispersion), 0.0_rp)) then
            write(stderr, *) "test_setup_gek failed: Excess dispersion set up "// &
                "without being asked for."
            test_setup_gek = .false.
        end if

        ! the excess dispersion is measured from the dispersion at the current point,
        ! which the regularization keeps small but above zero, that of the GEK without
        ! the excess dispersion
        call gek%setup(kappa_list, grad_list, func_list, hess, .true., error)
        if (error /= 0 .or. .not. gek%excess_dispersion) then
            write(stderr, *) "test_setup_gek failed: Excess dispersion not set up."
            test_setup_gek = .false.
        end if
        gek%excess_dispersion = .false.
        if (.not. is_close( &
            gek%current_dispersion, ref_dispersion_gek(gek, kappa_list(:, n_points)), &
            1e-8_rp * confidence_factor * sqrt(gek%process_variance))) then
            write(stderr, *) "test_setup_gek failed: Excess dispersion not "// &
                "measured from the current point."
            test_setup_gek = .false.
        end if

        ! a position that is not a number leaves a covariance matrix that is not
        ! positive definite
        kappa_list(1, 1) = ieee_value(1.0_rp, ieee_quiet_nan)
        call gek%setup(kappa_list, grad_list, func_list, hess, .false., error)
        if (error /= 1) then
            write(stderr, *) "test_setup_gek failed: Covariance matrix that is not "// &
                "positive definite not detected."
            test_setup_gek = .false.
        end if

    end function test_setup_gek

    logical(c_bool) function test_distances_gek() bind(C)
        !
        ! this function tests the subroutine which returns the squared scaled distances
        ! of the data points of the GEK from a point and their derivatives
        !
        use otr_s_gek, only: gek_type, distances_gek
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix

        type(gek_type) :: gek
        real(rp) :: kappa(n_dim), eig_kappa(n_dim), dist_sq(n_points), &
                    dist_sq_deriv(n_dim, n_points)
        integer(ip) :: k

        ! assume tests pass
        test_distances_gek = .true.

        ! random GEK and point
        gek%n_dim = n_dim
        gek%n_points = n_points
        gek%eigvecs = generate_random_orthogonal_matrix(n_dim)
        allocate(gek%length_scales(n_dim), gek%kappa_list(n_dim, n_points))
        call random_number(gek%length_scales)
        gek%length_scales = gek%length_scales + 0.5_rp
        call random_number(gek%kappa_list)
        call random_number(kappa)
        call distances_gek(gek, kappa, dist_sq, dist_sq_deriv)
        eig_kappa = matmul(transpose(gek%eigvecs), kappa)
        do k = 1, n_points
            if (.not. is_close(dist_sq(k), sum(((gek%kappa_list(:, k) - eig_kappa) / &
                                                gek%length_scales)**2))) then
                write(stderr, *) "test_distances_gek failed: Wrong distance."
                test_distances_gek = .false.
            end if
            if (.not. all(is_close(dist_sq_deriv(:, k), &
                                   2.0_rp * (gek%kappa_list(:, k) - eig_kappa) / &
                                   gek%length_scales**2))) then
                write(stderr, *) "test_distances_gek failed: Wrong derivative."
                test_distances_gek = .false.
            end if
        end do

    end function test_distances_gek

    logical(c_bool) function test_energy_gek() bind(C)
        !
        ! this function tests the function which predicts the energy of the GEK
        !
        use otr_s_gek, only: gek_type, func_regularization

        type(gek_type) :: gek
        real(rp) :: kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    func_list(n_points), point(n_dim)
        integer(ip) :: k

        ! assume tests pass
        test_energy_gek = .true.

        ! energy at random points and at the data points, which it reproduces up to the
        ! regularization times the solution vector
        call setup_random_gek(gek, kappa_list, grad_list, func_list)
        do k = 1, n_points
            call random_number(point)
            if (.not. is_close(gek%energy(point), ref_energy_gek(gek, point), &
                               roundoff_tolerance(gek))) then
                write(stderr, *) "test_energy_gek failed: Wrong energy."
                test_energy_gek = .false.
            end if
            if (.not. is_close(gek%energy(kappa_list(:, k)), &
                               func_list(k) - func_regularization * gek%weights(k), &
                               roundoff_tolerance(gek))) then
                write(stderr, *) "test_energy_gek failed: Energy of a data point "// &
                    "not reproduced."
                test_energy_gek = .false.
            end if
        end do

    end function test_energy_gek

    logical(c_bool) function test_gradient_gek() bind(C)
        !
        ! this function tests the function which predicts the gradient of the GEK
        !
        use otr_s_gek, only: gek_type, grad_regularization

        type(gek_type) :: gek
        real(rp) :: kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    func_list(n_points), point(n_dim)
        integer(ip) :: k

        ! assume tests pass
        test_gradient_gek = .true.

        ! gradient at random points and at the data points, which it reproduces up to
        ! the regularization times the solution vector
        call setup_random_gek(gek, kappa_list, grad_list, func_list)
        do k = 1, n_points
            call random_number(point)
            if (.not. all(is_close(gek%gradient(point), ref_gradient_gek(gek, point), &
                                   roundoff_tolerance(gek)))) then
                write(stderr, *) "test_gradient_gek failed: Wrong gradient."
                test_gradient_gek = .false.
            end if
            if (.not. all(is_close( &
                gek%gradient(kappa_list(:, k)), &
                grad_list(:, k) - grad_regularization * &
                matmul(gek%eigvecs, gek%weights(n_points + k::n_points)), &
                roundoff_tolerance(gek)))) then
                write(stderr, *) "test_gradient_gek failed: Gradient of a data "// &
                    "point not reproduced."
                test_gradient_gek = .false.
            end if
        end do

    end function test_gradient_gek

    logical(c_bool) function test_hessian_gek() bind(C)
        !
        ! this function tests the function which predicts the Hessian of the GEK
        !
        use otr_s_gek, only: gek_type

        type(gek_type) :: gek
        real(rp) :: kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    func_list(n_points), point(n_dim)
        integer(ip) :: k

        ! assume tests pass
        test_hessian_gek = .true.

        ! Hessian at random points and at a data point
        call setup_random_gek(gek, kappa_list, grad_list, func_list)
        do k = 1, n_points
            call random_number(point)
            if (.not. all(is_close(gek%hessian(point), ref_hessian_gek(gek, point), &
                                   roundoff_tolerance(gek)))) then
                write(stderr, *) "test_hessian_gek failed: Wrong Hessian."
                test_hessian_gek = .false.
            end if
        end do
        if (.not. all(is_close(gek%hessian(kappa_list(:, 1)), ref_hessian_gek( &
            gek, kappa_list(:, 1)), roundoff_tolerance(gek)))) then
            write(stderr, *) "test_hessian_gek failed: Wrong Hessian at a data point."
            test_hessian_gek = .false.
        end if

    end function test_hessian_gek

    logical(c_bool) function test_dispersion_gek() bind(C)
        !
        ! this function tests the function which returns the dispersion of the GEK
        !
        use otr_s_gek, only: gek_type, confidence_factor

        type(gek_type) :: gek
        real(rp) :: kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    func_list(n_points), point(n_dim), dispersion, sigma
        integer(ip) :: k

        ! assume tests pass
        test_dispersion_gek = .true.

        ! dispersion at random points, vanishing at the data points and approaching the
        ! confidence interval of the process far from them
        call setup_random_gek(gek, kappa_list, grad_list, func_list)
        sigma = confidence_factor * sqrt(gek%process_variance)
        do k = 1, n_points
            call random_number(point)
            if (.not. is_close(gek%dispersion(point), ref_dispersion_gek(gek, point), &
                               1e-8_rp * sigma)) then
                write(stderr, *) "test_dispersion_gek failed: Wrong dispersion."
                test_dispersion_gek = .false.
            end if
            if (.not. is_at_most(gek%dispersion(kappa_list(:, k)), 1e-3_rp * sigma)) &
                then
                write(stderr, *) "test_dispersion_gek failed: Dispersion does not "// &
                    "vanish at a data point."
                test_dispersion_gek = .false.
            end if
        end do
        point = 1e3_rp
        if (.not. is_close(gek%dispersion(point), sigma, 1e-8_rp * sigma)) then
            write(stderr, *) "test_dispersion_gek failed: Dispersion far from the "// &
                "data points does not approach the confidence interval of the process."
            test_dispersion_gek = .false.
        end if

        ! excess dispersion over a dispersion at the current point below and above that
        ! of the point
        call random_number(point)
        dispersion = ref_dispersion_gek(gek, point)
        gek%excess_dispersion = .true.
        gek%current_dispersion = 0.25_rp * dispersion
        if (.not. is_close(gek%dispersion(point), 0.75_rp * dispersion, &
                           1e-8_rp * sigma)) then
            write(stderr, *) "test_dispersion_gek failed: Wrong excess dispersion."
            test_dispersion_gek = .false.
        end if
        gek%current_dispersion = 2.0_rp * dispersion
        if (.not. is_at_most(abs(gek%dispersion(point)), 0.0_rp)) then
            write(stderr, *) "test_dispersion_gek failed: Negative excess "// &
                "dispersion not floored."
            test_dispersion_gek = .false.
        end if

    end function test_dispersion_gek

    logical(c_bool) function test_matern_kernel() bind(C)
        !
        ! this function tests the subroutine which evaluates the Matern 5/2 kernel and
        ! its first three derivatives with respect to the squared scaled distance
        !
        use otr_s_gek, only: matern_kernel

        integer(ip), parameter :: n_d = 4
        real(rp), parameter :: h = 1e-5_rp
        real(rp) :: d(n_d), r(n_d), k0(n_d), k1(n_d), k2(n_d), k3(n_d), k0p(n_d), &
                    k1p(n_d), k2p(n_d), k3p(n_d), k0m(n_d), k1m(n_d), k2m(n_d), &
                    k3m(n_d), z0, z1, z2, z3

        ! assume tests pass
        test_matern_kernel = .true.

        ! kernel at random distances
        call random_number(d)
        d = 0.1_rp + 4.0_rp * d
        call matern_kernel(d, k0, k1, k2, k3)
        call matern_kernel(d + h, k0p, k1p, k2p, k3p)
        call matern_kernel(d - h, k0m, k1m, k2m, k3m)
        r = sqrt(d)
        if (.not. all(is_close(k0, (1.0_rp + sqrt(5.0_rp) * r + 5.0_rp / 3.0_rp * d) * &
                               exp(-sqrt(5.0_rp) * r), 1e-14_rp))) then
            write(stderr, *) "test_matern_kernel failed: Wrong kernel."
            test_matern_kernel = .false.
        end if

        ! derivatives against central finite differences of the next lower one
        if (.not. all(is_close(k1, (k0p - k0m) / (2.0_rp * h), 1e-8_rp))) then
            write(stderr, *) "test_matern_kernel failed: Wrong first derivative."
            test_matern_kernel = .false.
        end if
        if (.not. all(is_close(k2, (k1p - k1m) / (2.0_rp * h), 1e-8_rp))) then
            write(stderr, *) "test_matern_kernel failed: Wrong second derivative."
            test_matern_kernel = .false.
        end if
        if (.not. all(is_close(k3, (k2p - k2m) / (2.0_rp * h), 1e-6_rp))) then
            write(stderr, *) "test_matern_kernel failed: Wrong third derivative."
            test_matern_kernel = .false.
        end if

        ! at a vanishing distance, where the third derivative is set to zero
        call matern_kernel(0.0_rp, z0, z1, z2, z3)
        if (.not. (is_close(z0, 1.0_rp) .and. is_close(z1, -5.0_rp / 6.0_rp) .and. &
                   is_close(z2, 25.0_rp / 12.0_rp) .and. is_at_most(abs(z3), 0.0_rp))) &
            then
            write(stderr, *) "test_matern_kernel failed: Wrong kernel at a "// &
                "vanishing distance."
            test_matern_kernel = .false.
        end if

    end function test_matern_kernel

    logical(c_bool) function test_minimize_surrogate_otr() bind(C)
        !
        ! this function tests the subroutine which minimizes the surrogate with the
        ! OpenTrustRegion solver
        !
        use otr_s_gek, only: gek_type, minimize_surrogate_otr
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        type(gek_type), target :: gek
        real(rp) :: current_grad(n_dim), minimum(n_dim), kappa(n_dim), &
                    free_kappa(n_dim), perturbation(n_dim), beta, draw, expected_draw
        integer, allocatable :: rng_state(:)
        integer :: rng_size
        logical :: restricted, converged

        ! assume tests pass
        test_minimize_surrogate_otr = .true.

        ! without a restriction, the solve ends at a stationary point of the surrogate
        call setup_quadratic_gek(gek, current_grad, minimum)
        call minimize_surrogate_otr(gek, current_grad, 0.0_rp * current_grad, 1e10_rp, &
                                    free_kappa, restricted, converged)
        if (.not. converged .or. &
            .not. is_below(norm2(ref_gradient_gek(gek, free_kappa)), 1e-6_rp)) then
            write(stderr, *) "test_minimize_surrogate_otr failed: Unrestricted "// &
                "solve does not end at a stationary point."
            test_minimize_surrogate_otr = .false.
        end if

        ! a variance threshold below the dispersion there restricts the point, which
        ! still lowers the surrogate energy
        beta = 0.5_rp * ref_dispersion_gek(gek, free_kappa)
        call minimize_surrogate_otr(gek, current_grad, 0.0_rp * current_grad, beta, &
                                    kappa, restricted, converged)
        if (.not. restricted .or. converged .or. &
            .not. is_below(ref_dispersion_gek(gek, kappa), beta) .or. .not. is_below( &
                ref_energy_gek(gek, kappa), ref_energy_gek(gek, 0.0_rp * kappa))) then
            write(stderr, *) "test_minimize_surrogate_otr failed: Restricted solve "// &
                "does not end within the variance threshold."
            test_minimize_surrogate_otr = .false.
        end if

        ! a perturbation of the gradient moves the stationary point to that of the
        ! perturbed surrogate, and the solve, which reseeds the random number
        ! generator, leaves it as it found it
        call random_number(perturbation)
        perturbation = 0.1_rp * (perturbation - 0.5_rp)
        call random_seed(size=rng_size)
        allocate(rng_state(rng_size))
        call random_seed(get=rng_state)
        call random_number(expected_draw)
        call random_seed(put=rng_state)
        call minimize_surrogate_otr(gek, current_grad, perturbation, 1e10_rp, kappa, &
                                    restricted, converged)
        call random_number(draw)
        if (.not. converged .or. .not. is_below( &
            norm2(ref_gradient_gek(gek, kappa) + perturbation), 1e-6_rp)) then
            write(stderr, *) "test_minimize_surrogate_otr failed: Solve with a "// &
                "perturbed gradient does not end at a stationary point of the "// &
                "perturbed surrogate."
            test_minimize_surrogate_otr = .false.
        end if
        if (.not. is_close(draw, expected_draw, 0.0_rp)) then
            write(stderr, *) "test_minimize_surrogate_otr failed: Random number "// &
                "generator not left as it was found."
            test_minimize_surrogate_otr = .false.
        end if

        ! a perturbation that is not a number lets the solve fail at once, which ends
        ! it unconverged and unrestricted at the current point
        perturbation = ieee_value(1.0_rp, ieee_quiet_nan)
        call minimize_surrogate_otr(gek, current_grad, perturbation, 1e10_rp, kappa, &
                                    restricted, converged)
        if (converged .or. restricted .or. .not. all(is_at_most(abs(kappa), 0.0_rp))) &
            then
            write(stderr, *) "test_minimize_surrogate_otr failed: Failed solve not "// &
                "ended unconverged at the current point."
            test_minimize_surrogate_otr = .false.
        end if

    end function test_minimize_surrogate_otr

    logical(c_bool) function test_minimize_surrogate_rs_rfo() bind(C)
        !
        ! this function tests the subroutine which minimizes the surrogate with the
        ! micro iterations of restricted-step RFO
        !
        use otr_s_gek, only: gek_type, minimize_surrogate_rs_rfo, &
                             step_bound_increment, step_bound_tol
        use opentrustregion_unit_tests, only: identity_matrix

        type(gek_type), target :: gek
        real(rp) :: current_grad(n_dim), minimum(n_dim), kappa(n_dim), &
                    free_kappa(n_dim), perturbation(n_dim), beta, &
                    kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    func_list(n_points)
        integer(ip) :: error, k
        logical :: restricted, converged

        ! assume tests pass
        test_minimize_surrogate_rs_rfo = .true.

        ! without a restriction, the micro iterations end at a stationary point of the
        ! surrogate
        call setup_quadratic_gek(gek, current_grad, minimum)
        call minimize_surrogate_rs_rfo(gek, current_grad, 0.0_rp * current_grad, &
                                       1e10_rp, 1.0_rp, free_kappa, restricted, &
                                       converged, error)
        if (error /= 0 .or. .not. converged .or. &
            .not. is_below(norm2(ref_gradient_gek(gek, free_kappa)), 1e-6_rp)) then
            write(stderr, *) "test_minimize_surrogate_rs_rfo failed: Unrestricted "// &
                "micro iterations do not end at a stationary point."
            test_minimize_surrogate_rs_rfo = .false.
        end if

        ! a variance threshold below the dispersion there restricts the point, which
        ! still lowers the surrogate energy
        beta = 0.5_rp * ref_dispersion_gek(gek, free_kappa)
        call minimize_surrogate_rs_rfo(gek, current_grad, 0.0_rp * current_grad, beta, &
                                       1.0_rp, kappa, restricted, converged, error)
        if (error /= 0 .or. .not. restricted .or. converged .or. &
            .not. is_below(ref_dispersion_gek(gek, kappa), beta) .or. .not. is_below( &
                ref_energy_gek(gek, kappa), ref_energy_gek(gek, 0.0_rp * kappa))) then
            write(stderr, *) "test_minimize_surrogate_rs_rfo failed: Restricted "// &
                "micro iterations do not end within the variance threshold."
            test_minimize_surrogate_rs_rfo = .false.
        end if

        ! step bounds too small to reach the stationary point within the micro
        ! iterations end them without convergence
        call minimize_surrogate_rs_rfo(gek, current_grad, 0.0_rp * current_grad, &
                                       1e10_rp, 1e-3_rp, kappa, restricted, converged, &
                                       error)
        if (error /= 0 .or. converged .or. restricted) then
            write(stderr, *) "test_minimize_surrogate_rs_rfo failed: Micro "// &
                "iterations without convergence reported as converged or restricted."
            test_minimize_surrogate_rs_rfo = .false.
        end if

        ! a perturbation of the gradient moves the stationary point to that of the
        ! perturbed surrogate, and a variance threshold far below the dispersion there
        ! restricts already the first step, which ends the micro iterations
        call random_number(perturbation)
        perturbation = 0.1_rp * (perturbation - 0.5_rp)
        call minimize_surrogate_rs_rfo(gek, current_grad, perturbation, 1e10_rp, &
                                       1.0_rp, kappa, restricted, converged, error)
        if (error /= 0 .or. .not. converged .or. .not. is_below( &
            norm2(ref_gradient_gek(gek, kappa) + perturbation), 1e-6_rp)) then
            write(stderr, *) "test_minimize_surrogate_rs_rfo failed: Micro "// &
                "iterations with a perturbed gradient do not end at a stationary "// &
                "point of the perturbed surrogate."
            test_minimize_surrogate_rs_rfo = .false.
        end if
        beta = 0.01_rp * beta
        call minimize_surrogate_rs_rfo(gek, current_grad, 0.0_rp * current_grad, beta, &
                                       1.0_rp, free_kappa, restricted, converged, error)
        call minimize_surrogate_rs_rfo(gek, current_grad, perturbation, beta, 1.0_rp, &
                                       kappa, restricted, converged, error)
        if (error /= 0 .or. .not. restricted .or. converged .or. &
            all(is_close(kappa, free_kappa, 0.0_rp))) then
            write(stderr, *) "test_minimize_surrogate_rs_rfo failed: Perturbation "// &
                "does not move the first step, which the variance threshold restricts."
            test_minimize_surrogate_rs_rfo = .false.
        end if
        if (.not. level_shifted_newton_step(ref_hessian_gek(gek, 0.0_rp * kappa), &
                                            current_grad + perturbation, kappa)) then
            write(stderr, *) "test_minimize_surrogate_rs_rfo failed: Micro "// &
                "iterations not ended after their first step, which the variance "// &
                "threshold restricts."
            test_minimize_surrogate_rs_rfo = .false.
        end if

        ! negative curvature of the surrogate of a concave quadratic is flipped, which
        ! ends the micro iterations after their first step, a step of the flipped
        ! Hessian within the first step bound
        call random_number(minimum)
        minimum = minimum - 0.5_rp
        call random_number(kappa_list)
        kappa_list = 0.4_rp * kappa_list - 0.2_rp
        kappa_list(:, n_points) = 0.0_rp
        do k = 1, n_points
            grad_list(:, k) = minimum - kappa_list(:, k)
            func_list(k) = -0.5_rp * sum((kappa_list(:, k) - minimum)**2)
        end do
        call ref_setup_gek(gek, kappa_list, grad_list, func_list, &
                           -identity_matrix(n_dim))
        current_grad = grad_list(:, n_points)
        call minimize_surrogate_rs_rfo(gek, current_grad, 0.0_rp * current_grad, &
                                       1e10_rp, 1.0_rp, kappa, restricted, converged, &
                                       error)
        if (error /= 0 .or. converged .or. restricted .or. &
            .not. is_at_most(norm2(kappa), step_bound_increment + step_bound_tol) .or. &
            .not. level_shifted_newton_step(ref_hessian_gek(gek, 0.0_rp * kappa), &
                                            current_grad, kappa)) then
            write(stderr, *) "test_minimize_surrogate_rs_rfo failed: Negative "// &
                "curvature not flipped or micro iterations not ended after its step."
            test_minimize_surrogate_rs_rfo = .false.
        end if

    contains

        function level_shifted_newton_step(hess, grad, step) result(is_step)
            !
            ! this function returns whether a step is a level-shifted Newton step for a
            ! gradient with a level shift below the lowest eigenvalue of a Hessian
            ! whose negative eigenvalues are flipped, as a single RS-RFO step is
            !
            real(rp), intent(in) :: hess(:, :), grad(:), step(:)
            logical :: is_step

            real(rp) :: eigvecs(n_dim, n_dim), eigvals(n_dim), work(10 * n_dim), &
                        flipped(n_dim, n_dim), mu
            integer(ip) :: info, j
            external :: dsyev

            eigvecs = 0.5_rp * (hess + transpose(hess))
            call dsyev("V", "U", n_dim, eigvecs, n_dim, eigvals, work, &
                       size(work, kind=ip), info)
            flipped = 0.0_rp
            do j = 1, n_dim
                flipped = flipped + abs(eigvals(j)) * spread( &
                    eigvecs(:, j), 2, n_dim) * spread(eigvecs(:, j), 1, n_dim)
            end do
            mu = dot_product(step, matmul(flipped, step) + grad) / &
                 dot_product(step, step)
            is_step = all(is_close(matmul(flipped, step) + grad, mu * step, &
                                   1e-8_rp)) .and. is_below(mu, minval(abs(eigvals)))

        end function level_shifted_newton_step

    end function test_minimize_surrogate_rs_rfo

    logical(c_bool) function test_surrogate_of() bind(C)
        !
        ! this function tests the function which returns the surrogate the context of
        ! the inner solve points at
        !
        use otr_s_gek, only: surrogate_type, surrogate_of
        use test_reference, only: host_context_type

        type(surrogate_type), target :: surrogate
        type(host_context_type), target :: other
        class(*), pointer :: context

        ! assume tests pass
        test_surrogate_of = .true.

        ! a surrogate and another context
        context => surrogate
        if (.not. associated(surrogate_of(context), surrogate)) then
            write(stderr, *) "test_surrogate_of failed: Surrogate not returned."
            test_surrogate_of = .false.
        end if
        context => other
        if (associated(surrogate_of(context))) then
            write(stderr, *) "test_surrogate_of failed: Surrogate returned for "// &
                "another context."
            test_surrogate_of = .false.
        end if

    end function test_surrogate_of

    logical(c_bool) function test_update_orbs_surrogate() bind(C)
        !
        ! this function tests the subroutine which moves the inner solve on the
        ! surrogate
        !
        use otr_s_gek, only: gek_type, surrogate_type, update_orbs_surrogate, &
                             hess_x_surrogate

        type(gek_type), target :: gek
        type(surrogate_type), target :: surrogate
        class(*), pointer :: context
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp) :: kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    func_list(n_points), start(n_dim), kappa(n_dim), func, &
                    grad(n_dim), h_diag(n_dim), hess(n_dim, n_dim), expected(n_dim)
        integer(ip) :: error, i, i_case

        ! assume tests pass
        test_update_orbs_surrogate = .true.

        ! surrogate at the current point with a real gradient that differs from the
        ! surrogate's and a random perturbation of the gradient
        call setup_random_gek(gek, kappa_list, grad_list, func_list)
        surrogate%gek => gek
        allocate(surrogate%kappa(n_dim), surrogate%step(n_dim), &
                 surrogate%hess(n_dim, n_dim), surrogate%perturbation(n_dim))
        allocate(surrogate%current_grad(n_dim))
        call random_number(surrogate%current_grad)
        call random_number(surrogate%perturbation)
        call random_number(surrogate%current_func)
        context => surrogate

        ! a vanishing step at the current point returns the real gradient there, a
        ! random step from a random position and a vanishing one at a random position
        ! the surrogate's gradient at the position it moves to, both perturbed, and the
        ! energy with the linear term of the perturbation there
        do i_case = 1, 3
            start = 0.0_rp
            kappa = 0.0_rp
            if (i_case > 1) call random_number(start)
            if (i_case == 2) call random_number(kappa)
            surrogate%kappa = start
            call update_orbs_surrogate(kappa, func, grad, h_diag, hess_x_funptr, &
                                       error, context)
            if (error /= 0) then
                write(stderr, *) "test_update_orbs_surrogate failed: Produced error."
                test_update_orbs_surrogate = .false.
            end if
            kappa = start + kappa
            if (.not. all(is_close(surrogate%kappa, kappa)) .or. &
                .not. all(is_close(surrogate%step, kappa - start))) then
                write(stderr, *) "test_update_orbs_surrogate failed: Position or "// &
                    "step not recorded."
                test_update_orbs_surrogate = .false.
            end if
            expected = surrogate%current_grad + surrogate%perturbation
            if (i_case > 1) &
                expected = ref_gradient_gek(gek, kappa) + surrogate%perturbation
            if (.not. all(is_close(grad, expected, roundoff_tolerance(gek))) .or. &
                .not. all(is_close(surrogate%grad, expected, &
                                   roundoff_tolerance(gek)))) then
                write(stderr, *) "test_update_orbs_surrogate failed: Wrong gradient."
                test_update_orbs_surrogate = .false.
            end if
            hess = ref_hessian_gek(gek, kappa)
            if (.not. all(is_close(surrogate%hess, hess, roundoff_tolerance(gek)))) then
                write(stderr, *) "test_update_orbs_surrogate failed: Wrong Hessian."
                test_update_orbs_surrogate = .false.
            end if
            do i = 1, n_dim
                if (.not. is_close(h_diag(i), hess(i, i), roundoff_tolerance(gek))) then
                    write(stderr, *) "test_update_orbs_surrogate failed: Wrong "// &
                        "Hessian diagonal."
                    test_update_orbs_surrogate = .false.
                end if
            end do
            if (.not. is_close(func, ref_energy_gek( &
                gek, kappa) - surrogate%current_func + dot_product( &
                    surrogate%perturbation, kappa), roundoff_tolerance(gek))) then
                write(stderr, *) "test_update_orbs_surrogate failed: Wrong energy "// &
                    "with the perturbation relative to the current point."
                test_update_orbs_surrogate = .false.
            end if
            if (.not. associated(hess_x_funptr, hess_x_surrogate)) then
                write(stderr, *) "test_update_orbs_surrogate failed: Hessian "// &
                    "linear transformation of the surrogate not returned."
                test_update_orbs_surrogate = .false.
            end if
        end do

    end function test_update_orbs_surrogate

    logical(c_bool) function test_hess_x_surrogate() bind(C)
        !
        ! this function tests the subroutine which applies the surrogate's Hessian
        !
        use otr_s_gek, only: surrogate_type, hess_x_surrogate

        type(surrogate_type), target :: surrogate
        class(*), pointer :: context
        real(rp) :: x(n_dim), hess_x(n_dim)
        integer(ip) :: error

        ! assume tests pass
        test_hess_x_surrogate = .true.

        ! random Hessian and vector
        allocate(surrogate%hess(n_dim, n_dim))
        call random_number(surrogate%hess)
        call random_number(x)
        context => surrogate
        call hess_x_surrogate(x, hess_x, error, context)
        if (error /= 0 .or. .not. all(is_close(hess_x, matmul(surrogate%hess, x)))) then
            write(stderr, *) "test_hess_x_surrogate failed: Wrong Hessian linear "// &
                "transformation."
            test_hess_x_surrogate = .false.
        end if

    end function test_hess_x_surrogate

    logical(c_bool) function test_obj_func_surrogate() bind(C)
        !
        ! this function tests the function which returns the surrogate energy at a
        ! trial point of the inner solve
        !
        use otr_s_gek, only: gek_type, surrogate_type, obj_func_surrogate
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        type(gek_type), target :: gek
        type(surrogate_type), target :: surrogate
        class(*), pointer :: context
        real(rp) :: kappa_list(n_dim, n_points), grad_list(n_dim, n_points), &
                    func_list(n_points), kappa(n_dim), dispersion, func
        integer(ip) :: error

        ! assume tests pass
        test_obj_func_surrogate = .true.

        ! surrogate at a random position with a random perturbation of the gradient
        call setup_random_gek(gek, kappa_list, grad_list, func_list)
        surrogate%gek => gek
        allocate(surrogate%kappa(n_dim), surrogate%perturbation(n_dim))
        call random_number(surrogate%kappa)
        call random_number(surrogate%perturbation)
        call random_number(surrogate%current_func)
        context => surrogate
        call random_number(kappa)
        dispersion = ref_dispersion_gek(gek, surrogate%kappa + kappa)

        ! a trial point within the variance threshold returns its energy with the
        ! linear term of the perturbation relative to the current point
        surrogate%variance_threshold = 2.0_rp * dispersion
        func = obj_func_surrogate(kappa, error, context)
        if (error /= 0 .or. surrogate%rejected .or. .not. is_close( &
            func, &
            ref_energy_gek(gek, surrogate%kappa + kappa) - surrogate%current_func + &
            dot_product(surrogate%perturbation, surrogate%kappa + kappa), &
            roundoff_tolerance(gek))) then
            write(stderr, *) "test_obj_func_surrogate failed: Wrong energy with "// &
                "the perturbation within the variance threshold."
            test_obj_func_surrogate = .false.
        end if

        ! a trial point beyond it is rejected and recorded
        surrogate%variance_threshold = 0.5_rp * dispersion
        func = obj_func_surrogate(kappa, error, context)
        if (error /= 0 .or. .not. surrogate%rejected .or. ieee_is_finite(func) .or. &
            .not. func > 0.0_rp) then
            write(stderr, *) "test_obj_func_surrogate failed: Trial point beyond "// &
                "the variance threshold not rejected."
            test_obj_func_surrogate = .false.
        end if

    end function test_obj_func_surrogate

    logical(c_bool) function test_conv_check_surrogate() bind(C)
        !
        ! this function tests the function which stops the inner solve
        !
        use otr_s_gek, only: surrogate_type, conv_check_surrogate, inner_conv_tol

        type(surrogate_type), target :: surrogate
        class(*), pointer :: context
        integer(ip) :: error, i_case
        logical :: expected
        character(len=40) :: case_name

        ! assume tests pass
        test_conv_check_surrogate = .true.

        ! a rejection, which stops the inner solve whatever the gradient, and converged
        ! and unconverged gradients, whose convergence it returns
        allocate(surrogate%kappa(n_dim), surrogate%grad(n_dim), surrogate%step(n_dim))
        context => surrogate
        do i_case = 1, 3
            surrogate%rejected = i_case == 1
            surrogate%grad = merge(0.5_rp, 2.0_rp, i_case == 2) * inner_conv_tol
            surrogate%step = 0.0_rp
            expected = i_case <= 2
            select case (i_case)
            case (1)
                case_name = "rejection"
            case (2)
                case_name = "converged gradient"
            case (3)
                case_name = "unconverged gradient"
            end select
            if ((conv_check_surrogate(error, context) .neqv. expected) .or. &
                error /= 0) then
                write(stderr, *) "test_conv_check_surrogate failed: Wrong "// &
                    "convergence for "//trim(case_name)//"."
                test_conv_check_surrogate = .false.
            end if
        end do

    end function test_conv_check_surrogate

    logical(c_bool) function test_inner_converged() bind(C)
        !
        ! this function tests the function which returns whether an inner solve has
        ! converged
        !
        use otr_s_gek, only: inner_converged, inner_conv_tol

        real(rp) :: grad(n_dim), step(n_dim)
        integer(ip) :: i_case
        character(len=40) :: case_name

        ! assume tests pass
        test_inner_converged = .true.

        ! converged gradient and step, and a gradient, a step and a single element of
        ! the step above their thresholds, four and six times the gradient's for the
        ! step
        do i_case = 1, 4
            grad = 0.5_rp * inner_conv_tol
            step = 2.0_rp * inner_conv_tol
            select case (i_case)
            case (1)
                case_name = "converged gradient and step"
            case (2)
                grad = 2.0_rp * inner_conv_tol
                case_name = "gradient above threshold"
            case (3)
                step = 5.0_rp * inner_conv_tol
                case_name = "step above threshold"
            case (4)
                step = 0.0_rp
                step(1) = 6.5_rp * inner_conv_tol
                case_name = "step element above threshold"
            end select
            if (inner_converged(grad, step) .neqv. i_case == 1) then
                write(stderr, *) "test_inner_converged failed: Wrong convergence "// &
                    "for "//trim(case_name)//"."
                test_inner_converged = .false.
            end if
        end do

    end function test_inner_converged

    logical(c_bool) function test_rs_rfo() bind(C)
        !
        ! this function tests the subroutine which returns the restricted-step RFO step
        !
        use otr_s_gek, only: rs_rfo, step_bound_tol
        use otr_common_unit_tests, only: generate_random_orthogonal_matrix
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

        real(rp), parameter :: indefinite_eigenvalues(n_full) = &
            [-1.0_rp, -0.8_rp, -0.6_rp, -0.4_rp, -0.2_rp, 0.2_rp, 10.0_rp, 12.0_rp]
        real(rp) :: hess(n_full, n_full), grad(n_full), step(n_full), &
                    full_step(n_full), eigvecs(n_full, n_full), mu, lowest, &
                    aug(n_full + 1, n_full + 1), aug_eigvals(n_full + 1), &
                    work(10 * (n_full + 1)), shortest(n_full), trial(n_full)
        integer(ip) :: error, i, info
        external :: dsyev

        ! assume tests pass
        test_rs_rfo = .true.

        ! for a positive definite Hessian, the unrestricted step is that of the lowest
        ! eigenvector of the augmented Hessian, a level-shifted Newton step with a
        ! level shift below the lowest eigenvalue
        hess = generate_random_spd_matrix(n_full)
        call random_number(grad)
        grad = grad - 0.5_rp
        call rs_rfo(hess, grad, 1e3_rp, full_step, error)
        if (error /= 0) then
            write(stderr, *) "test_rs_rfo failed: Unrestricted step produced error."
            test_rs_rfo = .false.
        end if
        lowest = lowest_eigenvalue(hess)
        mu = dot_product(grad, full_step)
        if (.not. all(is_close(matmul(hess, full_step) + grad, mu * full_step, &
                               1e-10_rp)) .or. .not. is_below(mu, lowest)) then
            write(stderr, *) "test_rs_rfo failed: Unrestricted step is not the RFO "// &
                "step."
            test_rs_rfo = .false.
        end if

        ! a restricted step has the step bound as length and is a level-shifted Newton
        ! step with a level shift below the lowest eigenvalue
        call rs_rfo(hess, grad, 0.3_rp * norm2(full_step), step, error)
        if (error /= 0) then
            write(stderr, *) "test_rs_rfo failed: Restricted step produced error."
            test_rs_rfo = .false.
        end if
        if (.not. is_close(norm2(step), 0.3_rp * norm2(full_step), step_bound_tol)) then
            write(stderr, *) "test_rs_rfo failed: Restricted step does not have "// &
                "the length of the step bound."
            test_rs_rfo = .false.
        end if
        mu = dot_product(step, matmul(hess, step) + grad) / dot_product(step, step)
        if (.not. all(is_close(matmul(hess, step) + grad, mu * step, 1e-10_rp)) .or. &
            .not. is_below(mu, lowest)) then
            write(stderr, *) "test_rs_rfo failed: Restricted step is not a "// &
                "level-shifted Newton step."
            test_rs_rfo = .false.
        end if

        ! for an indefinite Hessian, the step is the shortest of the steps of the
        ! lowest seven eigenvectors of the augmented Hessian, whose eigenvalues lie
        ! between those of the Hessian, which leave the seventh the widest interval
        eigvecs = generate_random_orthogonal_matrix(n_full)
        hess = 0.0_rp
        do i = 1, n_full
            hess(:, i) = eigvecs(:, i) * indefinite_eigenvalues(i)
        end do
        hess = matmul(hess, transpose(eigvecs))

        ! shortest of the steps of the lowest seven eigenvectors of the augmented
        ! Hessian
        aug = 0.0_rp
        aug(:n_full, :n_full) = hess
        aug(:n_full, n_full + 1) = grad
        aug(n_full + 1, :n_full) = grad
        call dsyev("V", "U", n_full + 1, aug, n_full + 1, aug_eigvals, work, &
                   size(work, kind=ip), info)
        shortest = huge(1.0_rp)
        do i = 1, 7
            trial = aug(:n_full, i) / aug(n_full + 1, i)
            if (norm2(trial) < norm2(shortest)) shortest = trial
        end do
        call rs_rfo(hess, grad, 1e3_rp, step, error)
        if (error /= 0 .or. .not. all(is_close(step, shortest, 1e-10_rp))) then
            write(stderr, *) "test_rs_rfo failed: Step for an indefinite Hessian "// &
                "is not the shortest step of the lowest eigenvectors."
            test_rs_rfo = .false.
        end if

        ! a gradient that is not a number leaves no step
        grad(1) = ieee_value(1.0_rp, ieee_quiet_nan)
        call rs_rfo(hess, grad, 1e3_rp, step, error)
        if (error /= 1) then
            write(stderr, *) "test_rs_rfo failed: Missing step not detected."
            test_rs_rfo = .false.
        end if

    contains

        function lowest_eigenvalue(matrix) result(eigenvalue)
            !
            ! this function returns the lowest eigenvalue of a symmetric matrix
            !
            real(rp), intent(in) :: matrix(:, :)
            real(rp) :: eigenvalue

            real(rp) :: a(size(matrix, 1), size(matrix, 1)), eigvals(size(matrix, 1)), &
                        work(10 * size(matrix, 1))
            integer(ip) :: info
            external :: dsyev

            a = matrix
            call dsyev("N", "U", size(a, 1, kind=ip), a, size(a, 1, kind=ip), eigvals, &
                       work, size(work, kind=ip), info)
            eigenvalue = eigvals(1)

        end function lowest_eigenvalue

    end function test_rs_rfo

end module otr_s_gek_unit_tests
