! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek

    use opentrustregion_core, only: rp, ip, kw_len, settings_type, hess_x_type
    use otr_common, only: obj_func_grad_type

    implicit none

    ! interface for the chart transformation callback function, while otr_common holds
    ! that of the objective function that also returns the gradient at a trial point
    abstract interface
        subroutine change_reference_type(new_ref, kappa_list, local_grad_list, &
                                         grad_list, error, context)
            import :: rp, ip

            real(rp), intent(in), target :: new_ref(:)
            real(rp), intent(inout), target :: kappa_list(:, :), &
                                               local_grad_list(:, :), grad_list(:, :)
            integer(ip), intent(out) :: error
            class(*), intent(in), pointer :: context
        end subroutine change_reference_type
    end interface

    ! derived type for settings, whose optional callback functions are the objective
    ! function that also returns the gradient at a trial point and the static Hessian
    type, extends(settings_type) :: s_gek_settings_type
        procedure(obj_func_grad_type), pointer, nopass :: obj_func_grad => null()
        procedure(hess_x_type), pointer, nopass :: init_hess => null()
        character(len=kw_len) :: inner_solver = "otr", prior_update = "sr1", &
                                 acceptance = "trial"
        logical :: excess_dispersion = .true.
    end type

    type(s_gek_settings_type), parameter :: default_s_gek_settings = &
        s_gek_settings_type()

    ! options of the inner solver, of the update of the prior and of the acceptance
    character(len=kw_len), parameter :: &
        inner_solver_options(2) = [character(len=kw_len) :: "otr", "rs-rfo"], &
        prior_update_options(3) = [character(len=kw_len) :: "sr1", "bfgs", "none"], &
        acceptance_options(2) = [character(len=kw_len) :: "trial", "none"]

    ! warning for settings that were not initialized
    character(len=*), parameter :: s_gek_settings_uninitialized_warning_msg = &
        "S-GEK settings were not initialized. All S-GEK settings are set to "// &
        "default values"

    ! error origin of S-GEK, with which it tags the errors it detects itself, and of
    ! the chart transformation, with which its invoke routine tags the errors it returns
    integer(ip), parameter :: error_s_gek = 5500, error_change_reference = 7800

    ! maximum numbers of points of the history and of the window the GEK is fitted to,
    ! and number of micro iterations of the minimization of the surrogate
    integer(ip), parameter :: max_history_points = 30, max_gek_points = 20, n_micro = 50

    ! variance threshold: its floor and its cap, which the bound factor scales
    real(rp), parameter :: variance_threshold_floor = 5e-3_rp, &
                           variance_threshold_cap = 0.05_rp

    ! micro iterations of restricted-step rational function optimization: increment of
    ! the step bound per micro iteration, which the bound factor scales, tolerance of
    ! the length of a restricted step, scale of the step bound below which the variance
    ! threshold no longer shortens it, number of the lowest eigenvectors of the
    ! augmented Hessian whose steps are compared and maximum numbers of the doublings
    ! and the bisections of the scaling that shortens a step to the step bound
    real(rp), parameter :: step_bound_increment = 0.1_rp, step_bound_tol = 1e-7_rp, &
                           min_bound_scale = 1e-5_rp
    integer(ip), parameter :: n_rfo_eigvecs = 7, max_rfo_doublings = 100, &
                              max_rfo_bisections = 200

    ! convergence of the inner solves: threshold of the RMS of the gradient of the
    ! surrogate and the factors of it below which the RMS and every element of the step
    ! have to lie
    real(rp), parameter :: inner_conv_tol = 1e-7_rp, inner_step_rms_factor = 4.0_rp, &
                           inner_step_max_factor = 6.0_rp

    ! second-order range: norm of the gradient of the previous point below which it
    ! starts, and the bound factor outside it
    real(rp), parameter :: second_order_grad_norm = 0.075_rp, far_bound_factor = 1e4_rp

    ! stretching of the prior along the proposal while consecutive steps point in the
    ! same direction: growth of the stretch factor, the golden ratio, and the cosines
    ! of 5 degrees, above which the angle between consecutive steps lets it grow, and
    ! of 20 degrees, below which it resets it
    real(rp), parameter :: stretch_growth = 0.5_rp * (1.0_rp + sqrt(5.0_rp)), &
                           stretch_grow_cos = cos(atan(1.0_rp) / 9.0_rp), &
                           stretch_reset_cos = cos(4.0_rp * atan(1.0_rp) / 9.0_rp)

    ! GEK: bias of the trend above the highest energy, floor of the magnitudes of the
    ! prior's eigenvalues from which the length scales follow, regularizations of the
    ! energies and gradients and factor of the half width of the 95% confidence
    ! interval of a normal distribution
    real(rp), parameter :: trend_bias = 10.0_rp, prior_eigval_floor = 0.025_rp, &
                           func_regularization = 1e-13_rp, &
                           grad_regularization = 1e-10_rp, &
                           confidence_factor = 1.9599639845400542355_rp

    ! floors of the Hessian diagonal, of negative and of small elements
    real(rp), parameter :: h_diag_negative_floor = 1.0_rp, h_diag_floor = 0.05_rp

    ! gradient-enhanced Kriging model of the objective function: eigenvectors of the
    ! prior Hessian, in whose basis the kernel works, length scales, positions of the
    ! data points in this basis, Cholesky factor of the covariance matrix, weights of
    ! the energies and gradients of the data points, which the covariance matrix maps
    ! to them, trend, process variance and the dispersion at the current point, from
    ! which the excess dispersion is measured
    type :: gek_type
        integer(ip) :: n_dim = 0, n_points = 0
        real(rp), allocatable :: eigvecs(:, :), length_scales(:), kappa_list(:, :), &
                                 cov_chol(:, :), weights(:)
        real(rp) :: trend = 0.0_rp, process_variance = 0.0_rp, &
                    current_dispersion = 0.0_rp
        logical :: excess_dispersion = .false.
    contains
        procedure :: setup => setup_gek
        procedure :: energy => energy_gek
        procedure :: gradient => gradient_gek
        procedure :: hessian => hessian_gek
        procedure :: dispersion => dispersion_gek
    end type

    ! surrogate the inner solve minimizes, which its callback functions reach through
    ! the context of its settings: the GEK, the position, the step and the gradient of
    ! the last orbital update, the projected gradient at the current point, the Hessian
    ! at the position, the energy of the GEK at the current point, from which the
    ! objective function is measured, the variance threshold and whether a trial point
    ! reached it
    type :: surrogate_type
        type(gek_type), pointer :: gek => null()
        real(rp), allocatable :: kappa(:), step(:), grad(:), current_grad(:), &
                                 perturbation(:), hess(:, :)
        real(rp) :: current_func = 0.0_rp, variance_threshold = 0.0_rp
        logical :: rejected = .false.
    end type

contains

    subroutine s_gek_solver(update_orbs, obj_func, change_reference, n_param, error, &
                            settings, state, s_gek_settings)
        !
        ! this subroutine is the subspace gradient-enhanced Kriging optimizer of
        ! OpenMolcas (Lindh and co-workers) with the interface of the OpenTrustRegion
        ! solver: the orbital update moves to every point S-GEK takes, into whose chart
        ! the host's chart transformation, which receives the context of the S-GEK
        ! settings like their callback functions, re-expresses the history; from the
        ! solver settings it honours the convergence threshold or check, the number of
        ! macro iterations, the stability check at convergence with its settings and
        ! Hessian linear transformation, the projection, the gradient noise, and the
        ! logging, and it ends at the maximum number of macro iterations with the
        ! solver's error code
        !
        use opentrustregion_core, only: &
            update_orbs_type, obj_func_type, solver_settings_type, solver_state_type, &
            default_solver_settings, init_rng, init_defaults, solver_sanity_check, &
            invoke_update_orbs, invoke_obj_func, invoke_project, &
            check_stationary_point, add_error_origin, error_conv_check, &
            error_solver_max_iter, func_resolution_safety_factor, verbosity_info, &
            verbosity_warning, verbosity_error, perturb_vector, noise_level, &
            function_unchanged_warning_msg, stdout
        use otr_common, only: invoke_obj_func_grad
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        procedure(update_orbs_type), intent(in), pointer :: update_orbs
        procedure(obj_func_type), intent(in), pointer :: obj_func
        procedure(change_reference_type), intent(in), pointer :: change_reference
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        class(solver_settings_type), intent(in), optional :: settings
        type(solver_state_type), intent(out), target, optional :: state
        type(s_gek_settings_type), intent(in), optional :: s_gek_settings

        type(solver_settings_type) :: local_settings
        type(s_gek_settings_type) :: gek_settings
        type(solver_state_type), target :: local_state
        type(solver_state_type), pointer :: current_state
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), allocatable :: &
            kappa_list(:, :), local_grad_list(:, :), grad_list(:, :), func_list(:), &
            grad(:), h_diag(:), floored_diag(:), proposal(:), kappa(:), trial_grad(:), &
            escape(:), previous_step(:), zero(:), perturbation(:)
        logical, allocatable :: moved_list(:)
        integer(ip) :: n_points, n_eval, n_bfgs_points, n_gek_points, first, &
                       previous_point
        real(rp) :: func, grad_rms, func_resolution, bound_factor, beta, beta_factor, &
                    pred_func_diff, actual_func_diff, current_func, previous_func, &
                    trial_func, stretch_factor, stretch_rate, last_grad_norm, &
                    cos_angle, kappa_norm
        logical :: second_order_range, stretch_active, have_last_step, &
                   last_pred_small, last_unresolved, max_precision_reached, converged, &
                   stable, accepted, have_pred, macro_converged, local_region, &
                   restricted, inner_converged, trial_added
        character(len=300) :: msg
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        ! state in the given record or a local one, both from the default values
        if (present(state)) then
            current_state => state
        else
            current_state => local_state
        end if

        ! copies of the settings, or of the default settings, with their dependent
        ! defaults resolved and checked
        if (present(settings)) then
            local_settings = settings
        else
            local_settings = default_solver_settings
        end if
        call init_rng(local_settings%seed)
        call init_defaults(local_settings)
        call solver_sanity_check(local_settings, n_param, error)
        call add_error_origin(error, error_s_gek, local_settings)
        if (error /= 0) return
        if (present(s_gek_settings)) gek_settings = s_gek_settings
        call s_gek_sanity_check(gek_settings, error)
        call add_error_origin(error, error_s_gek, gek_settings)
        if (error /= 0) return

        ! history of up to its maximum number of points, the current one last
        allocate(kappa_list(n_param, max_history_points), &
                 local_grad_list(n_param, max_history_points), &
                 grad_list(n_param, max_history_points), &
                 func_list(max_history_points), moved_list(max_history_points), &
                 grad(n_param), h_diag(n_param), trial_grad(n_param), escape(n_param), &
                 zero(n_param))
        n_points = 0
        n_bfgs_points = 0
        zero = 0.0_rp

        ! print header
        call local_settings%log(repeat("-", 77), verbosity_info)
        call local_settings%log(" Iteration |     Objective function     | "// &
                                "Gradient RMS |   Step size  | Step", verbosity_info)
        call local_settings%log(repeat("-", 77), verbosity_info)

        ! starting point
        n_eval = 0
        call move(zero)
        if (error /= 0) return
        call append(zero, grad, grad, func, .true., .false.)
        call log_row(func, 0.0_rp, .false., .true., .false., grad)

        ! initialize the stretching, the second-order range, the factor of the variance
        ! threshold, the maximum precision and the region of the gradient noise
        stretch_factor = 1.0_rp
        stretch_rate = stretch_growth
        stretch_active = .false.
        second_order_range = .false.
        last_grad_norm = huge(1.0_rp)
        beta_factor = 1.0_rp
        have_last_step = .false.
        last_pred_small = .false.
        last_unresolved = .false.
        max_precision_reached = .false.
        converged = .false.
        macro_converged = .false.
        local_region = .false.

        do
            ! the start of the second-order range restarts the points of the BFGS
            ! proposal from the current point
            if (.not. second_order_range .and. &
                last_grad_norm < second_order_grad_norm) then
                second_order_range = .true.
                n_bfgs_points = 1
            end if

            ! gradient norm, floored Hessian diagonal, energy resolution and bound
            ! factor of the current point, and the state
            last_grad_norm = dnrm2(n_param, grad, 1_ip)
            floored_diag = floored_h_diag(h_diag)
            grad_rms = last_grad_norm / sqrt(real(n_param, kind=rp))
            func_resolution = func_resolution_safety_factor * epsilon(1.0_rp) * &
                              abs(func)
            bound_factor = merge(1.0_rp, far_bound_factor, second_order_range)
            current_state%macro_iteration = n_eval - 1
            current_state%func = func
            current_state%grad_rms = grad_rms

            ! maximum precision as for the OpenTrustRegion solver: neither the energy
            ! change of the last step nor its prediction can be resolved, or the step
            ! was rejected while its prediction could not be resolved
            max_precision_reached = have_last_step .and. last_pred_small .and. &
                                    last_unresolved

            ! an exactly vanishing gradient, which a symmetric problem can reach,
            ! leaves no step, so the point is taken as one at maximum precision
            if (last_grad_norm <= 0.0_rp) max_precision_reached = .true.
            current_state%max_precision_reached = max_precision_reached

            ! check for convergence, which takes precedence over maximum precision, and
            ! stability, whose escape step S-GEK takes
            if (associated(local_settings%conv_check)) then
                converged = local_settings%conv_check(error, local_settings%context)
                call add_error_origin(error, error_conv_check, local_settings)
                if (error /= 0) return
            else
                converged = grad_rms < local_settings%conv_tol
            end if
            if (converged) then
                max_precision_reached = .false.
                current_state%max_precision_reached = .false.
            else if (max_precision_reached) then
                call local_settings%log(function_unchanged_warning_msg, &
                                        verbosity_warning)
            end if
            if (converged .or. max_precision_reached) then
                if (stretch_active .and. .not. max_precision_reached) then
                    ! convergence while the prior is stretched switches the stretching
                    ! off
                    stretch_factor = 1.0_rp
                    stretch_rate = 1.0_rp
                    stretch_active = .false.
                else
                    ! the stability check is always performed at a stationary starting
                    ! point or at an exactly vanishing gradient, as by the solver
                    if (local_settings%stability .or. n_eval == 1 .or. &
                        last_grad_norm <= 0.0_rp) then
                        call check_stationary_point(func, grad, grad_rms, h_diag, &
                                                    obj_func, hess_x_funptr, n_eval, &
                                                    local_settings, current_state, &
                                                    stable, escape, kappa_norm, error)
                        call add_error_origin(error, error_s_gek, gek_settings)
                        if (error /= 0) return
                    else
                        stable = .true.
                    end if
                    if (stable .or. current_state%max_precision_reached) then
                        macro_converged = .true.
                        exit
                    end if
                    if (n_eval >= local_settings%n_macro) exit
                    previous_func = func
                    call move(escape)
                    if (error /= 0) return
                    call append(zero, grad, grad, func, .true., .false.)
                    call log_row(func, kappa_norm, .false., .true., .false., grad)
                    current_state%step_size = kappa_norm
                    have_last_step = .false.
                    beta_factor = 1.0_rp
                    cycle
                end if
            end if

            ! stop before the step of the last macro iteration, since the point it
            ! reaches would not be checked
            if (n_eval >= local_settings%n_macro) exit

            ! BFGS proposal over the points added since the start of the second-order
            ! range, and window of the GEK
            first = max(1_ip, n_points - n_bfgs_points + 1)
            proposal = inverse_bfgs_step(kappa_list(:, first:n_points), &
                                         grad_list(:, first:n_points), floored_diag)
            if (associated(local_settings%project)) then
                call invoke_project(proposal, local_settings, error)
                call add_error_origin(error, error_s_gek, gek_settings)
                if (error /= 0) return
            end if
            n_gek_points = min(n_points, max_gek_points)
            stretch_active = stretch_factor > 1.0_rp

            ! perturbation of the gradient of the surrogate as the solver's gradient
            ! noise, of the size of the local region after a step within its bound,
            ! which the first step, without a surrogate, does not take
            perturbation = 0.0_rp * grad
            if (local_settings%grad_noise .and. n_gek_points > 1) then
                perturbation = grad
                call perturb_vector(perturbation, noise_level( &
                    local_settings, local_region), local_settings, error)
                call add_error_origin(error, error_s_gek, gek_settings)
                if (error /= 0) return
                perturbation = perturbation - grad
            end if

            ! step: for a single point, the proposal, which is moved to without a
            ! prediction, else the S-GEK step
            allocate(kappa(n_param))
            if (n_gek_points == 1) then
                kappa = proposal
                restricted = .false.
                inner_converged = .true.
                have_pred = .false.
            else
                ! variance threshold of the last energy change, from the previous point
                ! moved to to the current point, scaled by its factor
                beta = beta_factor * &
                       variance_threshold_of(func - previous_func, bound_factor)
                call s_gek_step( &
                    gek_settings, kappa_list(:, n_points - n_gek_points + 1:n_points), &
                    grad_list(:, n_points - n_gek_points + 1:n_points), &
                    func_list(n_points - n_gek_points + 1:n_points), floored_diag, &
                    proposal, perturbation, stretch_factor, bound_factor, beta, kappa, &
                    restricted, inner_converged, pred_func_diff, error)
                call add_error_origin(error, error_s_gek, gek_settings)
                if (error /= 0) return
                have_pred = .true.
            end if
            local_region = have_pred .and. inner_converged

            ! stretch factor from the angle between the step and the previous step, the
            ! position of the previous point moved to, before the current point, last
            ! in the history, reversed, if the history still holds it
            previous_point = findloc(moved_list(:n_points - 1), .true., dim=1, &
                                     back=.true.)
            if (stretch_rate > 1.0_rp .and. n_gek_points > 1 .and. previous_point > 0) &
                then
                previous_step = -kappa_list(:, previous_point)
                cos_angle = ddot(n_param, previous_step, 1_ip, kappa, 1_ip) / max( &
                    dnrm2(n_param, previous_step, 1_ip) * dnrm2(n_param, kappa, 1_ip), &
                    tiny(1.0_rp))
                if (cos_angle < stretch_reset_cos .or. .not. inner_converged) then
                    stretch_factor = 1.0_rp
                else if (cos_angle > stretch_grow_cos) then
                    stretch_factor = stretch_factor * stretch_rate
                end if
            end if

            ! evaluate the step
            current_func = func
            trial_added = .false.
            if (gek_settings%acceptance == "trial" .and. have_pred) then
                ! the trial point is evaluated without moving, with its gradient by the
                ! objective function of the S-GEK settings if they provide one, which
                ! keeps a rejected one as data before the current point unless its
                ! energy is infinite, and otherwise by its energy alone with the
                ! objective function of the solver, which leaves the history unchanged;
                ! the step is taken only if it lowers the energy
                if (associated(gek_settings%obj_func_grad)) then
                    call invoke_obj_func_grad(gek_settings%obj_func_grad, kappa, &
                                              trial_func, trial_grad, gek_settings, &
                                              error_s_gek, error)
                    current_state%n_obj_func = current_state%n_obj_func + 1
                    call add_error_origin(error, error_s_gek, gek_settings)
                    if (error /= 0) return
                    trial_added = ieee_is_finite(trial_func)
                else
                    trial_func = invoke_obj_func(obj_func, kappa, local_settings, &
                                                 current_state%n_obj_func, error)
                    if (error /= 0) return
                end if
                n_eval = n_eval + 1
                accepted = ieee_is_finite(trial_func) .and. trial_func < current_func
                actual_func_diff = huge(1.0_rp)
                if (ieee_is_finite(trial_func)) &
                    actual_func_diff = trial_func - current_func
                if (trial_added .and. .not. accepted) then
                    call add_trial(kappa, trial_grad, trial_func)
                    if (error /= 0) return
                end if
                if (accepted) then
                    ! move to the trial point, which becomes the current point, last in
                    ! the history
                    previous_func = current_func
                    call invoke_update_orbs(update_orbs, kappa, func, grad, h_diag, &
                                            hess_x_funptr, local_settings, &
                                            current_state%n_update_orbs, error)
                    call add_error_origin(error, error_s_gek, gek_settings)
                    if (error /= 0) return
                    call invoke_change_reference( &
                        change_reference, kappa, kappa_list(:, :n_points), &
                        local_grad_list(:, :n_points), grad_list(:, :n_points), &
                        gek_settings, error)
                    call add_error_origin(error, error_s_gek, gek_settings)
                    if (error /= 0) return
                    call append(zero, grad, grad, func, .true., .false.)
                end if

                ! the row of the trial point, with its gradient where it is known, by
                ! the objective function of the S-GEK settings or by the move to it
                if (associated(gek_settings%obj_func_grad)) then
                    call log_row(trial_func, dnrm2(n_param, kappa, 1_ip), restricted, &
                                 inner_converged, .not. accepted, trial_grad)
                else if (accepted) then
                    call log_row(trial_func, dnrm2(n_param, kappa, 1_ip), restricted, &
                                 inner_converged, .false., grad)
                else
                    call log_row(trial_func, dnrm2(n_param, kappa, 1_ip), restricted, &
                                 inner_converged, .true.)
                end if
            else
                previous_func = current_func
                call move(kappa)
                if (error /= 0) return
                call append(zero, grad, grad, func, .true., .false.)
                call log_row(func, dnrm2(n_param, kappa, 1_ip), restricted, &
                             inner_converged, .false., grad)
                actual_func_diff = func - current_func
                accepted = .true.
            end if

            ! factor of the variance threshold
            beta_factor = updated_threshold_factor(beta_factor, accepted, trial_added)

            ! whether the energy change of the step and its prediction can be resolved
            have_last_step = have_pred
            if (have_pred) then
                last_pred_small = abs(pred_func_diff) <= func_resolution
                last_unresolved = abs(actual_func_diff) <= func_resolution .or. &
                                  actual_func_diff > 0.0_rp
            end if
            current_state%step_size = dnrm2(n_param, kappa, 1_ip)
            deallocate(kappa)
        end do

        ! stop if no convergence
        if (.not. macro_converged) then
            call local_settings%log("Orbital optimization has not converged!", &
                                    verbosity_error, .true.)
            error = error_solver_max_iter
            return
        end if

        ! finish logging as the solver does, with the evaluations S-GEK takes
        call local_settings%log(repeat("-", 77), verbosity_info)
        write(msg, '(A, I0)') "Total number of orbital updates: ", &
            current_state%n_update_orbs
        call local_settings%log(msg, verbosity_info)
        write(msg, '(A, I0)') "Total number of objective function evaluations: ", &
            current_state%n_obj_func
        call local_settings%log(msg, verbosity_info)

        ! flush output
        flush(stdout)

    contains

        subroutine move(point)
            !
            ! this subroutine moves to a point with the orbital update and re-expresses
            ! the history in its chart
            !
            real(rp), intent(in) :: point(:)

            ! orbital update
            call invoke_update_orbs(update_orbs, point, func, grad, h_diag, &
                                    hess_x_funptr, local_settings, &
                                    current_state%n_update_orbs, error)
            call add_error_origin(error, error_s_gek, gek_settings)
            if (error /= 0) return
            n_eval = n_eval + 1

            ! chart transformation of the history
            if (n_points > 0) call invoke_change_reference( &
                change_reference, point, kappa_list(:, :n_points), &
                local_grad_list(:, :n_points), grad_list(:, :n_points), gek_settings, &
                error)
            call add_error_origin(error, error_s_gek, gek_settings)

        end subroutine move

        subroutine append(new_kappa, new_local_grad, new_grad, new_func, moved, &
                          before_current)
            !
            ! this subroutine adds a point to the history, as its last point or before
            ! the current point, dropping the oldest one when it is full, and counts it
            ! among the points of the BFGS proposal
            !
            real(rp), intent(in) :: new_kappa(:), new_local_grad(:), new_grad(:), &
                                    new_func
            logical, intent(in) :: moved, before_current

            integer(ip) :: pos

            ! drop the oldest point of a full history
            if (n_points == max_history_points) then
                kappa_list(:, :n_points - 1) = kappa_list(:, 2:n_points)
                local_grad_list(:, :n_points - 1) = local_grad_list(:, 2:n_points)
                grad_list(:, :n_points - 1) = grad_list(:, 2:n_points)
                func_list(:n_points - 1) = func_list(2:n_points)
                moved_list(:n_points - 1) = moved_list(2:n_points)
                n_points = n_points - 1
            end if

            ! position of the new point, before which the current point moves to the end
            n_points = n_points + 1
            pos = n_points
            if (before_current) then
                pos = n_points - 1
                kappa_list(:, n_points) = kappa_list(:, pos)
                local_grad_list(:, n_points) = local_grad_list(:, pos)
                grad_list(:, n_points) = grad_list(:, pos)
                func_list(n_points) = func_list(pos)
                moved_list(n_points) = moved_list(pos)
            end if

            ! store the new point
            kappa_list(:, pos) = new_kappa
            local_grad_list(:, pos) = new_local_grad
            grad_list(:, pos) = new_grad
            func_list(pos) = new_func
            moved_list(pos) = moved
            n_bfgs_points = min(n_bfgs_points + 1, max_history_points)

        end subroutine append

        subroutine add_trial(trial_kappa, trial_local_grad, trial_point_func)
            !
            ! this subroutine adds a trial point, whose position is in the chart of the
            ! current point and whose gradient is local, before the current point, with
            ! its gradient in that chart
            !
            real(rp), intent(in) :: trial_kappa(:), trial_local_grad(:), &
                                    trial_point_func

            real(rp) :: kappa_col(n_param, 1), local_grad_col(n_param, 1), &
                        grad_col(n_param, 1)

            ! gradient of the trial point in the chart of the current point
            kappa_col(:, 1) = trial_kappa
            local_grad_col(:, 1) = trial_local_grad
            grad_col(:, 1) = trial_local_grad
            call invoke_change_reference(change_reference, zero, kappa_col, &
                                         local_grad_col, grad_col, gek_settings, error)
            call add_error_origin(error, error_s_gek, gek_settings)
            if (error /= 0) return

            ! add it before the current point
            call append(kappa_col(:, 1), local_grad_col(:, 1), grad_col(:, 1), &
                        trial_point_func, .false., .true.)

        end subroutine add_trial

        subroutine log_row(row_func, step_norm, row_restricted, row_converged, &
                           row_rejected, row_grad)
            !
            ! this subroutine prints the row of an evaluated point in the format of the
            ! OpenTrustRegion solver, with the RMS of its gradient if it is known and
            ! with "x" if its step was rejected, "*" if the variance threshold
            ! restricted it, "#" if its inner solve did not converge otherwise and " "
            ! else
            !
            real(rp), intent(in) :: row_func, step_norm
            logical, intent(in) :: row_restricted, row_converged, row_rejected
            real(rp), intent(in), optional :: row_grad(:)

            character :: row_marker
            character(len=9) :: grad_rms_text

            ! marker of the step and RMS of the gradient, then write and print the row
            row_marker = merge("*", merge(" ", "#", row_converged), row_restricted)
            if (row_rejected) row_marker = "x"
            grad_rms_text = ""
            if (present(row_grad)) write(grad_rms_text, '(1PE9.2)') &
                dnrm2(size(row_grad, kind=ip), row_grad, 1_ip) / &
                sqrt(real(size(row_grad), kind=rp))
            write(msg, '(4X, I4, 3X, "|", 4X, 1PE21.14, 3X, "|", 2X, A, 3X, "|", '// &
                  '2X, 1PE9.2, 3X, "|", 2X, A)') n_eval - 1, row_func, grad_rms_text, &
                step_norm, row_marker
            call local_settings%log(msg, verbosity_info)

        end subroutine log_row

    end subroutine s_gek_solver

    subroutine s_gek_sanity_check(settings, error)
        !
        ! this subroutine performs a sanity check for S-GEK input parameters, after
        ! converting the options to lowercase
        !
        use opentrustregion_core, only: verbosity_error, string_to_lowercase, &
                                        option_list

        type(s_gek_settings_type), intent(inout) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! convert strings to lowercase
        settings%inner_solver = string_to_lowercase(settings%inner_solver)
        settings%prior_update = string_to_lowercase(settings%prior_update)
        settings%acceptance = string_to_lowercase(settings%acceptance)

        ! check for character options
        if (.not. any(settings%inner_solver == inner_solver_options)) then
            call settings%log("Inner solver option unknown. Possible values are "// &
                              option_list(inner_solver_options)//".", verbosity_error, &
                              .true.)
            error = 1
            return
        end if
        if (.not. any(settings%prior_update == prior_update_options)) then
            call settings%log("Prior update option unknown. Possible values are "// &
                              option_list(prior_update_options)//".", verbosity_error, &
                              .true.)
            error = 1
            return
        end if
        if (.not. any(settings%acceptance == acceptance_options)) then
            call settings%log("Acceptance option unknown. Possible values are "// &
                              option_list(acceptance_options)//".", verbosity_error, &
                              .true.)
            error = 1
            return
        end if

    end subroutine s_gek_sanity_check

    subroutine invoke_change_reference(change_reference, new_ref, kappa_list, &
                                       local_grad_list, grad_list, settings, error)
        !
        ! this subroutine calls the host's chart transformation with the context of the
        ! S-GEK settings and tags its error with the origin of the chart
        ! transformation, while it returns an error of its own with the origin of S-GEK
        ! for positions or gradients that are not finite
        !
        use opentrustregion_core, only: verbosity_error, add_error_origin
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        procedure(change_reference_type), intent(in), pointer :: change_reference
        real(rp), intent(in), target :: new_ref(:)
        real(rp), intent(inout), target :: kappa_list(:, :), local_grad_list(:, :), &
                                           grad_list(:, :)
        type(s_gek_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        ! transform the history
        call change_reference(new_ref, kappa_list, local_grad_list, grad_list, error, &
                              settings%context)
        call add_error_origin(error, error_change_reference, settings)
        if (error /= 0) return

        ! check the result
        if (.not. (all(ieee_is_finite(kappa_list)) .and. &
                   all(ieee_is_finite(local_grad_list)) .and. &
                   all(ieee_is_finite(grad_list)))) then
            call settings%log("Chart transformation returned values that are not "// &
                              "finite.", verbosity_error, .true.)
            error = 1
            call add_error_origin(error, error_s_gek, settings)
        end if

    end subroutine invoke_change_reference

    function floored_h_diag(h_diag) result(floored)
        !
        ! this function floors a Hessian diagonal: negative elements to the larger of
        ! their magnitude and the floor of negative elements, small elements to the
        ! floor of small elements
        !
        real(rp), intent(in) :: h_diag(:)
        real(rp) :: floored(size(h_diag))

        ! floor negative and small elements
        where (h_diag < 0.0_rp)
            floored = max(h_diag_negative_floor, abs(h_diag))
        elsewhere (h_diag < h_diag_floor)
            floored = h_diag_floor
        elsewhere
            floored = h_diag
        end where

    end function floored_h_diag

    elemental function variance_threshold_of(func_diff, bound_factor) result(threshold)
        !
        ! this function returns the variance threshold of an energy change, its
        ! magnitude between the floor of the threshold and its cap, which the bound
        ! factor scales, equivalent to the variance threshold of OpenMolcas
        !
        real(rp), intent(in) :: func_diff, bound_factor
        real(rp) :: threshold

        threshold = min(variance_threshold_cap * bound_factor, &
                        max(variance_threshold_floor, abs(func_diff)))

    end function variance_threshold_of

    elemental function updated_threshold_factor(factor, accepted, trial_added) &
        result(updated)
        !
        ! this function returns the factor of the variance threshold after a step: one
        ! once the step is taken, half the factor after a rejected trial point that
        ! leaves the history unchanged, after which the same step would follow, and the
        ! factor unchanged after a rejected trial point that is kept as data
        !
        use opentrustregion_core, only: trust_radius_shrink_factor

        real(rp), intent(in) :: factor
        logical, intent(in) :: accepted, trial_added
        real(rp) :: updated

        if (accepted) then
            updated = 1.0_rp
        else if (trial_added) then
            updated = factor
        else
            updated = trust_radius_shrink_factor * factor
        end if

    end function updated_threshold_factor

    function inverse_bfgs_step(kappa_list, grad_list, h_diag) result(step)
        !
        ! this function returns the step -B^-1 g of the BFGS updates B of a positive
        ! Hessian diagonal, such as the floored one, over the pairs of consecutive
        ! points of a history, for the gradient g of the last point
        !
        use otr_common, only: bfgs_terms, woodbury_x

        real(rp), intent(in) :: kappa_list(:, :), grad_list(:, :), h_diag(:)
        real(rp) :: step(size(kappa_list, 1))

        real(rp), allocatable :: steps(:, :), grad_diffs(:, :), vectors(:, :), coeffs(:)
        logical, allocatable :: kept(:)
        integer(ip) :: n_points, n_terms

        ! terms of the BFGS updates of the pairs
        n_points = size(kappa_list, 2, kind=ip)
        steps = kappa_list(:, 2:) - kappa_list(:, :n_points - 1)
        grad_diffs = grad_list(:, 2:) - grad_list(:, :n_points - 1)
        allocate(vectors(size(kappa_list, 1), 2 * (n_points - 1)), &
                 coeffs(2 * (n_points - 1)), kept(n_points - 1))
        call bfgs_terms(steps, grad_diffs, spread(h_diag, 2, n_points - 1) * steps, &
                        vectors, coeffs, n_terms, kept)

        ! inverse of the updated diagonal applied to the gradient
        step = -grad_list(:, n_points) / h_diag
        call woodbury_x(vectors, coeffs, n_terms, &
                        vectors(:, :n_terms) / spread(h_diag, 2, n_terms), step)

    end function inverse_bfgs_step

    subroutine s_gek_step(settings, kappa_list, grad_list, func_list, h_diag, &
                          proposal, perturbation, stretch_factor, bound_factor, beta, &
                          kappa, restricted, converged, pred_func_diff, error)
        !
        ! this subroutine returns the S-GEK step of a window of points (columns of
        ! kappa_list and grad_list, energies func_list) in the chart of the current,
        ! last one: the minimum of the GEK of the data projected onto the subspace of
        ! the gradient and position differences, the current gradient, the proposal and
        ! the perturbation of the gradient, with the projected prior, the floored
        ! Hessian diagonal or the host's static Hessian, updated by SR1 or BFGS over
        ! the window if the settings ask for it and stretched along the proposal by the
        ! stretch factor, within the given variance threshold, for a gradient perturbed
        ! by the perturbation as by the solver's gradient noise, the bound factor
        ! scaling the step bound of restricted-step RFO; returns whether the variance
        ! threshold restricted the step and whether the inner solve converged within
        ! it, and the energy change the unperturbed GEK predicts
        !
        use opentrustregion_core, only: add_error_origin, verbosity_error, &
                                        orthogonalize_trial_vectors
        use otr_common, only: invoke_init_hess

        type(s_gek_settings_type), intent(in) :: settings
        real(rp), intent(in) :: kappa_list(:, :), grad_list(:, :), func_list(:), &
                                h_diag(:), proposal(:), perturbation(:), &
                                stretch_factor, bound_factor, beta
        real(rp), intent(out) :: kappa(:), pred_func_diff
        logical, intent(out) :: restricted, converged
        integer(ip), intent(out) :: error

        real(rp), allocatable :: basis(:, :), sub_kappa_list(:, :), &
                                 sub_grad_list(:, :), prior(:, :), direction(:), &
                                 stretch(:, :), sub_kappa(:), column(:), image(:), &
                                 sub_perturbation(:)
        type(gek_type), target :: gek
        integer(ip) :: n_points, n_sub, k, j
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! orthonormal basis of the subspace of the differences, the current gradient,
        ! the proposal and the perturbation, a vanishing one of which drops out, which
        ! are routinely linearly dependent
        n_points = size(func_list, kind=ip)
        allocate(basis(size(kappa_list, 1), 2 * n_points + 1))
        do k = 1, n_points - 1
            basis(:, 2 * k - 1) = grad_list(:, k + 1) - grad_list(:, k)
            basis(:, 2 * k) = kappa_list(:, k + 1) - kappa_list(:, k)
        end do
        basis(:, 2 * n_points - 1) = grad_list(:, n_points)
        basis(:, 2 * n_points) = proposal
        basis(:, 2 * n_points + 1) = perturbation
        call orthogonalize_trial_vectors(basis, settings, error, warn_lin_dep=.false.)
        if (error /= 0) then
            call add_error_origin(error, error_s_gek, settings)
            return
        end if
        n_sub = size(basis, 2, kind=ip)

        ! projected positions relative to the current point and projected gradients
        sub_kappa_list = matmul( &
            transpose(basis), kappa_list - spread(kappa_list(:, n_points), 2, n_points))
        sub_grad_list = matmul(transpose(basis), grad_list)

        ! projected prior: the floored Hessian diagonal or the host's static Hessian
        if (associated(settings%init_hess)) then
            allocate(prior(n_sub, n_sub), column(size(kappa_list, 1)), &
                     image(size(kappa_list, 1)))
            do j = 1, n_sub
                column = basis(:, j)
                call invoke_init_hess(settings%init_hess, column, image, settings, &
                                      error_s_gek, error)
                if (error /= 0) return
                prior(:, j) = matmul(transpose(basis), image)
            end do
            prior = 0.5_rp * (prior + transpose(prior))
        else
            prior = matmul(transpose(basis), spread(h_diag, 2, n_sub) * basis)
        end if

        ! prior updated by SR1 or BFGS over the window
        if (settings%prior_update /= "none") prior = &
            updated_prior(prior, sub_kappa_list, sub_grad_list, settings%prior_update)

        ! prior stretched along the proposal by the stretch factor, which scales it
        ! down along the proposal
        if (abs(stretch_factor - 1.0_rp) > 0.0_rp) then
            direction = matmul(transpose(basis), proposal) / &
                        dnrm2(size(proposal, kind=ip), proposal, 1_ip)
            stretch = spread(direction, 2, n_sub) * spread(direction, 1, n_sub) * &
                      (1.0_rp / stretch_factor - 1.0_rp)
            do j = 1, n_sub
                stretch(j, j) = stretch(j, j) + 1.0_rp
            end do
            prior = matmul(stretch, matmul(prior, stretch))
        end if

        ! GEK of the projected data
        call gek%setup(sub_kappa_list, sub_grad_list, func_list, prior, &
                       settings%excess_dispersion, error)
        if (error /= 0) then
            call settings%log("Covariance matrix of the GEK is not positive "// &
                              "definite.", verbosity_error, .true.)
            call add_error_origin(error, error_s_gek, settings)
            return
        end if

        ! minimum of the surrogate for a gradient perturbed by the projected
        ! perturbation within the variance threshold
        sub_perturbation = matmul(transpose(basis), perturbation)
        allocate(sub_kappa(n_sub))
        if (settings%inner_solver == "rs-rfo") then
            call minimize_surrogate_rs_rfo(gek, sub_grad_list(:, n_points), &
                                           sub_perturbation, beta, bound_factor, &
                                           sub_kappa, restricted, converged, error)
            if (error /= 0) then
                call settings%log("Restricted-step rational function optimization "// &
                                  "of the surrogate failed.", verbosity_error, .true.)
                call add_error_origin(error, error_s_gek, settings)
                return
            end if
        else
            call minimize_surrogate_otr(gek, sub_grad_list(:, n_points), &
                                        sub_perturbation, beta, sub_kappa, restricted, &
                                        converged)
        end if

        ! predicted energy change and step in the full space
        pred_func_diff = gek%energy(sub_kappa) - gek%energy(0.0_rp * sub_kappa)
        kappa = matmul(basis, sub_kappa)

    end subroutine s_gek_step

    function updated_prior(prior, kappa_list, grad_list, scheme) result(updated)
        !
        ! this function returns the symmetrized prior updated by the SR1 or BFGS
        ! updates, as the scheme names them, over the pairs of consecutive points of a
        ! window, from the oldest to the newest
        !
        use otr_common, only: sr1_terms, bfgs_terms

        real(rp), intent(in) :: prior(:, :), kappa_list(:, :), grad_list(:, :)
        character(len=*), intent(in) :: scheme
        real(rp) :: updated(size(prior, 1), size(prior, 2))

        real(rp), allocatable :: steps(:, :), grad_diffs(:, :), vectors(:, :), coeffs(:)
        logical, allocatable :: kept(:)
        integer(ip) :: n, n_points, n_terms, k

        ! symmetrized prior
        n = size(prior, 1, kind=ip)
        updated = 0.5_rp * (prior + transpose(prior))

        ! terms of the updates of the pairs
        n_points = size(kappa_list, 2, kind=ip)
        steps = kappa_list(:, 2:) - kappa_list(:, :n_points - 1)
        grad_diffs = grad_list(:, 2:) - grad_list(:, :n_points - 1)
        allocate(vectors(n, 2 * (n_points - 1)), coeffs(2 * (n_points - 1)), &
                 kept(n_points - 1))
        if (scheme == "sr1") then
            call sr1_terms(steps, grad_diffs, matmul(updated, steps), vectors, coeffs, &
                           n_terms, kept)
        else
            call bfgs_terms(steps, grad_diffs, matmul(updated, steps), vectors, &
                            coeffs, n_terms, kept)
        end if

        ! add the terms
        do k = 1, n_terms
            updated = updated + coeffs(k) * spread(vectors(:, k), 2, n) * &
                      spread(vectors(:, k), 1, n)
        end do

        ! symmetrize against rounding
        updated = 0.5_rp * (updated + transpose(updated))

    end function updated_prior

    subroutine setup_gek(self, kappa_list, grad_list, func_list, prior, &
                         excess_dispersion, error)
        !
        ! this subroutine sets up the GEK of the positions (columns of kappa_list),
        ! gradients and energies of data points with a prior Hessian: the length scales
        ! reproduce the magnitudes of the prior's eigenvalues, floored, as the
        ! curvature at a single data point, the trend lies the bias above the highest
        ! energy, and the energies and gradients are regularized; an error is returned
        ! when the covariance matrix is not positive definite
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        class(gek_type), intent(inout) :: self
        real(rp), intent(in) :: kappa_list(:, :), grad_list(:, :), func_list(:), &
                                prior(:, :)
        logical, intent(in) :: excess_dispersion
        integer(ip), intent(out) :: error

        integer(ip) :: n_dim, n_points, n_cov, i, j, k, mm, info, lwork
        real(rp), allocatable :: eigvals(:), work(:), scaled_diff(:, :, :), &
                                 dist_sq(:, :), k0(:, :), k1(:, :), k2(:, :), &
                                 k3(:, :), eig_grad_list(:, :), targets(:)
        real(rp) :: value
        real(rp), external :: ddot
        external :: dsyev, dpotrf, dpotrs

        ! initialize error flag
        error = 0

        ! dimension, number of data points and size of the covariance matrix
        n_dim = size(kappa_list, 1, kind=ip)
        n_points = size(kappa_list, 2, kind=ip)
        n_cov = n_points * (n_dim + 1)
        self%n_dim = n_dim
        self%n_points = n_points
        self%excess_dispersion = excess_dispersion

        ! eigenbasis of the symmetrized prior and length scales from its eigenvalues
        self%eigvecs = 0.5_rp * (prior + transpose(prior))
        allocate(eigvals(n_dim), work(1))
        call dsyev("V", "U", n_dim, self%eigvecs, n_dim, eigvals, work, -1_ip, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))
        call dsyev("V", "U", n_dim, self%eigvecs, n_dim, eigvals, work, lwork, info)
        if (info /= 0) then
            error = 1
            return
        end if
        self%length_scales = sqrt(5.0_rp / 3.0_rp * trend_bias / &
                                  max(abs(eigvals), prior_eigval_floor))

        ! positions and gradients of the data points in the eigenbasis
        self%kappa_list = matmul(transpose(self%eigvecs), kappa_list)
        eig_grad_list = matmul(transpose(self%eigvecs), grad_list)

        ! scaled differences and squared scaled distances between the data points and
        ! the kernel of these distances
        allocate(scaled_diff(n_dim, n_points, n_points), dist_sq(n_points, n_points))
        do j = 1, n_points
            do k = 1, n_points
                scaled_diff(:, k, j) = &
                    (self%kappa_list(:, k) - self%kappa_list(:, j)) / self%length_scales
                dist_sq(k, j) = ddot(n_dim, scaled_diff(:, k, j), 1_ip, &
                                     scaled_diff(:, k, j), 1_ip)
            end do
        end do
        allocate(k0(n_points, n_points), k1(n_points, n_points), &
                 k2(n_points, n_points), k3(n_points, n_points))
        call matern_kernel(dist_sq, k0, k1, k2, k3)

        ! covariance of the energies, of the energies with the gradients, whose rows
        ! and columns run over the data points within blocks of one component, and of
        ! the gradients
        if (allocated(self%cov_chol)) deallocate(self%cov_chol)
        allocate(self%cov_chol(n_cov, n_cov))
        self%cov_chol(:n_points, :n_points) = k0
        do j = 1, n_points
            do i = 1, n_dim
                do k = 1, n_points
                    value = -k1(k, j) * 2.0_rp * scaled_diff(i, k, j) / &
                            self%length_scales(i)
                    self%cov_chol(k, n_points + (i - 1) * n_points + j) = value
                    self%cov_chol(n_points + (i - 1) * n_points + j, k) = value
                end do
            end do
        end do
        do j = 1, n_points
            do mm = 1, n_dim
                do k = 1, n_points
                    do i = 1, n_dim
                        value = -k2(k, j) * (2.0_rp * scaled_diff( &
                            i, k, j) / self%length_scales(i)) * (2.0_rp * scaled_diff( &
                                mm, k, j) / self%length_scales(mm))
                        if (i == mm) &
                            value = value - k1(k, j) * 2.0_rp / self%length_scales(i)**2
                        self%cov_chol(n_points + (i - 1) * n_points + k, &
                                      n_points + (mm - 1) * n_points + j) = value
                    end do
                end do
            end do
        end do

        ! regularization of the energies and gradients
        do k = 1, n_points
            self%cov_chol(k, k) = self%cov_chol(k, k) + func_regularization
        end do
        do k = n_points + 1, n_cov
            self%cov_chol(k, k) = self%cov_chol(k, k) + grad_regularization
        end do

        ! Cholesky factorization, which fails for a matrix that is not positive
        ! definite, also one that is not finite, which not every LAPACK detects
        call dpotrf("U", n_cov, self%cov_chol, n_cov, info)
        if (info /= 0 .or. &
            .not. all(ieee_is_finite([(self%cov_chol(k, k), k=1, n_cov)]))) then
            error = 1
            return
        end if

        ! trend above the highest energy and weights of the energies relative to it and
        ! of the gradients, together with the process variance
        self%trend = maxval(func_list) + trend_bias
        allocate(targets(n_cov))
        targets(:n_points) = func_list - self%trend
        do i = 1, n_dim
            targets(n_points + (i - 1) * n_points + 1:n_points + i * n_points) = &
                eig_grad_list(i, :)
        end do
        self%weights = targets
        call dpotrs("U", n_cov, 1_ip, self%cov_chol, n_cov, self%weights, n_cov, info)
        if (info /= 0) then
            error = 1
            return
        end if
        self%process_variance = ddot(n_cov, targets, 1_ip, self%weights, 1_ip) / &
                                real(n_cov, kind=rp)

        ! dispersion at the newest data point, the current point, from which the excess
        ! dispersion is measured
        self%current_dispersion = 0.0_rp
        if (excess_dispersion) then
            self%excess_dispersion = .false.
            self%current_dispersion = self%dispersion(kappa_list(:, n_points))
            self%excess_dispersion = .true.
        end if

    end subroutine setup_gek

    subroutine distances_gek(self, kappa, dist_sq, dist_sq_deriv)
        !
        ! this subroutine returns the squared scaled distances of the data points of
        ! the GEK from a point kappa and their derivatives 2 (x_k - x) / l^2 with
        ! respect to the positions x_k of the data points in the eigenbasis of the
        ! prior, in which the point is x
        !
        class(gek_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:)
        real(rp), intent(out) :: dist_sq(:), dist_sq_deriv(:, :)

        real(rp) :: eig_kappa(self%n_dim), scaled_diff(self%n_dim)
        integer(ip) :: k
        real(rp), external :: ddot

        ! position in the eigenbasis of the prior
        eig_kappa = matmul(transpose(self%eigvecs), kappa)

        ! squared scaled distances and their derivatives
        do k = 1, self%n_points
            scaled_diff = (self%kappa_list(:, k) - eig_kappa) / self%length_scales
            dist_sq(k) = ddot(self%n_dim, scaled_diff, 1_ip, scaled_diff, 1_ip)
            dist_sq_deriv(:, k) = 2.0_rp * scaled_diff / self%length_scales
        end do

    end subroutine distances_gek

    function energy_gek(self, kappa) result(energy)
        !
        ! this function returns the energy the GEK predicts at a point kappa, the trend
        ! plus the covariances of the energy at the point with the energies and
        ! gradients of the data points times their weights
        !
        class(gek_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:)
        real(rp) :: energy

        real(rp) :: dist_sq(self%n_points), dist_sq_deriv(self%n_dim, self%n_points), &
                    k0(self%n_points), k1(self%n_points), k2(self%n_points), &
                    k3(self%n_points)
        integer(ip) :: i, n_points
        real(rp), external :: ddot

        ! kernel of the squared scaled distances of the data points
        n_points = self%n_points
        call distances_gek(self, kappa, dist_sq, dist_sq_deriv)
        call matern_kernel(dist_sq, k0, k1, k2, k3)

        ! trend and the weighted covariances with the energies and the gradients
        energy = self%trend + ddot(n_points, k0, 1_ip, self%weights, 1_ip)
        do i = 1, self%n_dim
            energy = energy + &
                     ddot(n_points, k1 * dist_sq_deriv(i, :), 1_ip, &
                          self%weights(n_points + (i - 1) * n_points + 1:), 1_ip)
        end do

    end function energy_gek

    function gradient_gek(self, kappa) result(gradient)
        !
        ! this function returns the gradient the GEK predicts at a point kappa, the
        ! derivatives of the covariances of the energy at the point times the weights
        !
        class(gek_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:)
        real(rp) :: gradient(self%n_dim)

        real(rp) :: dist_sq(self%n_points), dist_sq_deriv(self%n_dim, self%n_points), &
                    k0(self%n_points), k1(self%n_points), k2(self%n_points), &
                    k3(self%n_points), grad_weights(self%n_dim, self%n_points), &
                    deriv_weights(self%n_points), eig_gradient(self%n_dim)
        integer(ip) :: k, n_points
        real(rp), external :: ddot

        ! kernel of the squared scaled distances of the data points
        n_points = self%n_points
        call distances_gek(self, kappa, dist_sq, dist_sq_deriv)
        call matern_kernel(dist_sq, k0, k1, k2, k3)

        ! weights of the gradients per component and data point and their products with
        ! the derivatives of the distances
        grad_weights = transpose(reshape(self%weights(n_points + 1:), &
                                         [n_points, self%n_dim]))
        do k = 1, n_points
            deriv_weights(k) = ddot(self%n_dim, grad_weights(:, k), 1_ip, &
                                    dist_sq_deriv(:, k), 1_ip)
        end do

        ! derivatives of the weighted covariances in the eigenbasis of the prior
        eig_gradient = 0.0_rp
        do k = 1, n_points
            eig_gradient = eig_gradient - (self%weights(k) * k1(k) + k2(k) * &
                                           deriv_weights(k)) * dist_sq_deriv(:, k) - &
                           k1(k) * grad_weights(:, k) * 2.0_rp / self%length_scales**2
        end do

        ! gradient in the original basis
        gradient = matmul(self%eigvecs, eig_gradient)

    end function gradient_gek

    function hessian_gek(self, kappa) result(hessian)
        !
        ! this function returns the Hessian the GEK predicts at a point kappa, the
        ! second derivatives of the covariances of the energy at the point times the
        ! weights
        !
        class(gek_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:)
        real(rp) :: hessian(self%n_dim, self%n_dim)

        real(rp) :: dist_sq(self%n_points), dist_sq_deriv(self%n_dim, self%n_points), &
                    k0(self%n_points), k1(self%n_points), k2(self%n_points), &
                    k3(self%n_points), grad_weights(self%n_dim, self%n_points), &
                    deriv_weights(self%n_points), scaled_grad_weights(self%n_dim), &
                    eig_hessian(self%n_dim, self%n_dim), diagonal
        integer(ip) :: i, k, n_points
        real(rp), external :: ddot

        ! kernel of the squared scaled distances of the data points
        n_points = self%n_points
        call distances_gek(self, kappa, dist_sq, dist_sq_deriv)
        call matern_kernel(dist_sq, k0, k1, k2, k3)

        ! weights of the gradients per component and data point and their products with
        ! the derivatives of the distances
        grad_weights = transpose(reshape(self%weights(n_points + 1:), &
                                         [n_points, self%n_dim]))
        do k = 1, n_points
            deriv_weights(k) = ddot(self%n_dim, grad_weights(:, k), 1_ip, &
                                    dist_sq_deriv(:, k), 1_ip)
        end do

        ! second derivatives of the weighted covariances in the eigenbasis of the
        ! prior, of which the terms proportional to the identity add to the diagonal
        eig_hessian = 0.0_rp
        diagonal = 0.0_rp
        do k = 1, n_points
            scaled_grad_weights = grad_weights(:, k) * 2.0_rp / self%length_scales**2
            eig_hessian = eig_hessian + &
                          (self%weights(k) * k2(k) + k3(k) * deriv_weights(k)) * &
                          spread(dist_sq_deriv(:, k), 2, self%n_dim) * &
                          spread(dist_sq_deriv(:, k), 1, self%n_dim) + &
                          k2(k) * (spread(dist_sq_deriv(:, k), 2, self%n_dim) * &
                                   spread(scaled_grad_weights, 1, self%n_dim) + &
                                   spread(scaled_grad_weights, 2, self%n_dim) * &
                                   spread(dist_sq_deriv(:, k), 1, self%n_dim))
            diagonal = diagonal + self%weights(k) * k1(k) + k2(k) * deriv_weights(k)
        end do
        do i = 1, self%n_dim
            eig_hessian(i, i) = eig_hessian(i, i) + &
                                diagonal * 2.0_rp / self%length_scales(i)**2
        end do

        ! Hessian in the original basis
        hessian = matmul(self%eigvecs, matmul(eig_hessian, transpose(self%eigvecs)))

    end function hessian_gek

    function dispersion_gek(self, kappa) result(dispersion)
        !
        ! this function returns the dispersion of the GEK at a point kappa, the half
        ! width of the 95% confidence interval of the energy it predicts, or, for the
        ! excess dispersion, its increase over the dispersion at the current point,
        ! where the regularization keeps it from vanishing
        !
        class(gek_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:)
        real(rp) :: dispersion

        real(rp) :: dist_sq(self%n_points), dist_sq_deriv(self%n_dim, self%n_points), &
                    k0(self%n_points), k1(self%n_points), k2(self%n_points), &
                    k3(self%n_points), cov(self%n_points * (self%n_dim + 1)), &
                    solved(self%n_points * (self%n_dim + 1)), posterior_variance
        integer(ip) :: i, n_points, n_cov, info
        real(rp), external :: ddot
        external :: dpotrs

        ! kernel of the squared scaled distances of the data points
        n_points = self%n_points
        n_cov = n_points * (self%n_dim + 1)
        call distances_gek(self, kappa, dist_sq, dist_sq_deriv)
        call matern_kernel(dist_sq, k0, k1, k2, k3)

        ! covariances of the energy at the point with the energies and gradients of the
        ! data points
        cov(:n_points) = k0
        do i = 1, self%n_dim
            cov(n_points + (i - 1) * n_points + 1:n_points + i * n_points) = &
                k1 * dist_sq_deriv(i, :)
        end do

        ! posterior variance and the half width of the confidence interval, or its
        ! excess over that at the current point
        solved = cov
        call dpotrs("U", n_cov, 1_ip, self%cov_chol, n_cov, solved, n_cov, info)
        posterior_variance = 1.0_rp - ddot(n_cov, cov, 1_ip, solved, 1_ip)
        dispersion = confidence_factor * &
                     sqrt(max(posterior_variance, 0.0_rp) * self%process_variance)
        if (self%excess_dispersion) &
            dispersion = max(0.0_rp, dispersion - self%current_dispersion)

    end function dispersion_gek

    elemental subroutine matern_kernel(dist_sq, k0, k1, k2, k3)
        !
        ! this subroutine evaluates the Matern 5/2 kernel of a squared scaled distance
        ! d, (1 + sqrt(5 d) + 5 d / 3) exp(-sqrt(5 d)), and its first three derivatives
        ! with respect to d, the third set to zero at a vanishing distance, where it
        ! diverges as 1 / sqrt(d) while the terms it enters vanish as d^(3/2)
        !
        real(rp), intent(in) :: dist_sq
        real(rp), intent(out) :: k0, k1, k2, k3

        real(rp) :: dist, decay

        ! scaled distance and exponential decay
        dist = sqrt(dist_sq)
        decay = exp(-sqrt(5.0_rp) * dist)

        ! kernel and its derivatives
        k0 = (1.0_rp + sqrt(5.0_rp) * dist + 5.0_rp / 3.0_rp * dist_sq) * decay
        k1 = -5.0_rp / 6.0_rp * (1.0_rp + sqrt(5.0_rp) * dist) * decay
        k2 = 25.0_rp / 12.0_rp * decay
        if (dist > epsilon(1.0_rp)) then
            k3 = -25.0_rp * sqrt(5.0_rp) / 24.0_rp / dist * decay
        else
            k3 = 0.0_rp
        end if

    end subroutine matern_kernel

    subroutine minimize_surrogate_otr(gek, current_grad, perturbation, beta, kappa, &
                                      restricted, converged)
        !
        ! this subroutine minimizes the surrogate for a gradient perturbed by the
        ! perturbation, whose linear term the energy includes, with the OpenTrustRegion
        ! solver, whose objective function rejects points at or beyond the variance
        ! threshold beta and whose convergence check stops after such a rejection, and
        ! which leaves the random number generator, which the solver reseeds, as it
        ! found it; the solve ends at its last accepted point also when its macro
        ! iterations run out or it fails; it returns whether the variance threshold
        ! restricted the point and whether the solve converged within it, which a stop
        ! at maximum precision counts as, since the surrogate cannot resolve a lower
        ! energy there
        !
        use opentrustregion_core, only: solver, solver_settings_type, &
                                        update_orbs_type, obj_func_type

        type(gek_type), intent(in), target :: gek
        real(rp), intent(in) :: current_grad(:), perturbation(:), beta
        real(rp), intent(out) :: kappa(:)
        logical, intent(out) :: restricted, converged

        type(surrogate_type), target :: surrogate
        integer, allocatable :: rng_state(:)
        integer :: rng_size
        type(solver_settings_type) :: settings
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        integer(ip) :: n_dim, inner_error

        ! surrogate at the current point
        n_dim = gek%n_dim
        surrogate%gek => gek
        allocate(surrogate%kappa(n_dim), surrogate%step(n_dim), &
                 surrogate%hess(n_dim, n_dim))
        surrogate%kappa = 0.0_rp
        surrogate%step = 0.0_rp
        surrogate%current_grad = current_grad
        surrogate%perturbation = perturbation
        surrogate%current_func = gek%energy(surrogate%kappa)
        surrogate%variance_threshold = beta

        ! silent solve with the convergence check of the surrogate and as many macro
        ! iterations as the micro iterations of restricted-step RFO
        settings%conv_check => conv_check_surrogate
        settings%n_macro = n_micro
        settings%stability = .false.
        settings%verbose = 0
        settings%context => surrogate
        update_orbs_funptr => update_orbs_surrogate
        obj_func_funptr => obj_func_surrogate
        call random_seed(size=rng_size)
        allocate(rng_state(rng_size))
        call random_seed(get=rng_state)
        call solver(update_orbs_funptr, obj_func_funptr, n_dim, inner_error, settings)
        call random_seed(put=rng_state)

        ! restriction, convergence and last accepted point
        restricted = surrogate%rejected
        converged = inner_error == 0 .and. .not. restricted
        kappa = surrogate%kappa

    end subroutine minimize_surrogate_otr

    subroutine minimize_surrogate_rs_rfo(gek, current_grad, perturbation, beta, &
                                         bound_factor, kappa, restricted, converged, &
                                         error)
        !
        ! this subroutine minimizes the surrogate for a gradient perturbed by the
        ! perturbation, by micro iterations of restricted-step RFO, as OpenMolcas does:
        ! negative curvature of the surrogate is flipped, which ends the micro
        ! iterations after its step, and the step bound grows by its increment, scaled
        ! by the bound factor, per micro iteration and is halved together with the step
        ! until the dispersion of the new point is below the variance threshold beta,
        ! which ends the micro iterations after that step; it returns whether the
        ! variance threshold restricted a step and whether the micro iterations
        ! converged within it
        !
        type(gek_type), intent(in) :: gek
        real(rp), intent(in) :: current_grad(:), perturbation(:), beta, bound_factor
        real(rp), intent(out) :: kappa(:)
        logical, intent(out) :: restricted, converged
        integer(ip), intent(out) :: error

        real(rp) :: grad(size(current_grad)), step(size(current_grad)), &
                    new_kappa(size(current_grad)), &
                    hess(size(current_grad), size(current_grad)), &
                    eigvals(size(current_grad)), bound_scale, step_bound
        real(rp), allocatable :: eigvecs(:, :), work(:)
        integer(ip) :: micro, n, i, lwork, info
        logical :: terminate
        external :: dsyev

        ! initialize error flag
        error = 0

        ! micro iterations from the current point
        n = size(current_grad, kind=ip)
        kappa = 0.0_rp
        grad = current_grad + perturbation
        restricted = .false.
        converged = .false.
        terminate = .false.
        do micro = 1, n_micro
            bound_scale = 1.0_rp
            step_bound = step_bound_increment * bound_factor * real(micro, kind=rp)

            ! Hessian of the surrogate with its negative eigenvalues flipped
            hess = gek%hessian(kappa)
            eigvecs = 0.5_rp * (hess + transpose(hess))
            allocate(work(1))
            call dsyev("V", "U", n, eigvecs, n, eigvals, work, -1_ip, info)
            lwork = int(work(1), kind=ip)
            deallocate(work)
            allocate(work(lwork))
            call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
            deallocate(work)
            if (info /= 0) then
                error = 1
                return
            end if
            do i = 1, n
                if (eigvals(i) < 0.0_rp) then
                    terminate = .true.
                    hess = hess + 2.0_rp * abs(eigvals(i)) * &
                           spread(eigvecs(:, i), 2, n) * spread(eigvecs(:, i), 1, n)
                end if
            end do

            ! step restricted to points whose dispersion is below the variance threshold
            do
                call rs_rfo(hess, grad, step_bound, step, error)
                if (error /= 0) return
                new_kappa = kappa + step
                if (bound_scale < min_bound_scale .or. &
                    gek%dispersion(new_kappa) < beta) exit
                bound_scale = 0.5_rp * bound_scale
                step_bound = 0.5_rp * step_bound
                restricted = .true.
            end do

            ! convergence of the micro iterations within the variance threshold, which
            ! end after a restricted step, after flipping negative curvature or at
            ! their maximum number
            kappa = new_kappa
            grad = gek%gradient(kappa) + perturbation
            converged = .not. (restricted .or. terminate) .and. &
                        inner_converged(grad, step)
            if (converged .or. restricted .or. terminate) exit
        end do

    end subroutine minimize_surrogate_rs_rfo

    function surrogate_of(context) result(surrogate)
        !
        ! this function returns the surrogate the context of the inner solve points at
        !
        class(*), intent(in), pointer :: context
        type(surrogate_type), pointer :: surrogate

        ! select the surrogate
        surrogate => null()
        select type (context)
        type is (surrogate_type)
            surrogate => context
        end select

    end function surrogate_of

    subroutine update_orbs_surrogate(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                     context)
        !
        ! this subroutine moves the inner solve on the surrogate by kappa and returns
        ! the energy of the GEK relative to the current point, whose resolution then
        ! does not scale with the total energy, its gradient, or at the current point
        ! the projected gradient, which the GEK reproduces only up to its
        ! regularization, perturbed by the perturbation, whose linear term the energy
        ! includes, and its Hessian
        !
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(inout), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(surrogate_type), pointer :: surrogate
        integer(ip) :: i
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! move the position
        surrogate => surrogate_of(context)
        surrogate%kappa = surrogate%kappa + kappa
        surrogate%step = kappa

        ! gradient, Hessian and its diagonal, and energy relative to the current point
        if (sum(abs(surrogate%kappa)) > 0.0_rp) then
            surrogate%grad = surrogate%gek%gradient(surrogate%kappa) + &
                             surrogate%perturbation
        else
            surrogate%grad = surrogate%current_grad + surrogate%perturbation
        end if
        surrogate%hess = surrogate%gek%hessian(surrogate%kappa)
        grad = surrogate%grad
        do i = 1, size(h_diag)
            h_diag(i) = surrogate%hess(i, i)
        end do
        func = surrogate%gek%energy(surrogate%kappa) - surrogate%current_func + ddot( &
            size(kappa, kind=ip), surrogate%perturbation, 1_ip, surrogate%kappa, 1_ip)
        hess_x_funptr => hess_x_surrogate

    end subroutine update_orbs_surrogate

    subroutine hess_x_surrogate(x, hess_x, error, context)
        !
        ! this subroutine applies the Hessian of the surrogate at the position of the
        ! inner solve
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        type(surrogate_type), pointer :: surrogate

        ! initialize error flag
        error = 0

        ! apply the Hessian
        surrogate => surrogate_of(context)
        hess_x = matmul(surrogate%hess, x)

    end subroutine hess_x_surrogate

    function obj_func_surrogate(kappa, error, context) result(func)
        !
        ! this function returns the energy of the GEK relative to the current point,
        ! including the linear term of the perturbation, at a trial point of the inner
        ! solve, or an infinite value, which rejects it, if its dispersion reaches the
        ! variance threshold, which it records
        !
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf

        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        type(surrogate_type), pointer :: surrogate
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! reject a trial point at or beyond the variance threshold
        surrogate => surrogate_of(context)
        if (surrogate%gek%dispersion(surrogate%kappa + kappa) >= &
            surrogate%variance_threshold) then
            surrogate%rejected = .true.
            func = ieee_value(func, ieee_positive_inf)
        else
            func = surrogate%gek%energy(surrogate%kappa + kappa) - &
                   surrogate%current_func + &
                   ddot(size(kappa, kind=ip), surrogate%perturbation, 1_ip, &
                        surrogate%kappa + kappa, 1_ip)
        end if

    end function obj_func_surrogate

    function conv_check_surrogate(error, context) result(converged)
        !
        ! this function stops the inner solve at the first point accepted after a trial
        ! point reached the variance threshold, or once the gradient and the step of
        ! the surrogate have converged
        !
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        logical :: converged

        type(surrogate_type), pointer :: surrogate

        ! initialize error flag
        error = 0

        ! convergence after a rejection or of the gradient and the step
        surrogate => surrogate_of(context)
        converged = surrogate%rejected .or. &
                    inner_converged(surrogate%grad, surrogate%step)

    end function conv_check_surrogate

    function inner_converged(grad, step) result(converged)
        !
        ! this function returns whether an inner solve has converged: the RMS of the
        ! gradient of the surrogate lies below the convergence threshold of the inner
        ! solves, the RMS of the step below inner_step_rms_factor times it and every
        ! element of the step below inner_step_max_factor times it
        !
        real(rp), intent(in) :: grad(:), step(:)
        logical :: converged

        integer(ip) :: n
        real(rp), external :: dnrm2

        n = size(grad, kind=ip)
        converged = &
            dnrm2(n, grad, 1_ip) / sqrt(real(n, kind=rp)) < inner_conv_tol .and. &
            dnrm2(n, step, 1_ip) / sqrt(real(n, kind=rp)) < inner_step_rms_factor * &
            inner_conv_tol .and. &
            maxval(abs(step)) < inner_step_max_factor * inner_conv_tol

    end function inner_converged

    subroutine rs_rfo(hess, grad, step_bound, step, error)
        !
        ! this subroutine returns the restricted-step rational function optimization
        ! (RS-RFO) step of Besalu and Bofill: of the lowest n_rfo_eigvecs eigenvectors
        ! of the augmented Hessian, or all of fewer, the one with the shortest step,
        ! and, for a step longer than the step bound, the scaling alpha of the Hessian
        ! that shortens it to the step bound within its tolerance
        !
        real(rp), intent(in) :: hess(:, :), grad(:), step_bound
        real(rp), intent(out) :: step(:)
        integer(ip), intent(out) :: error

        real(rp) :: step_norm, lower, upper, middle, middle_norm
        integer(ip) :: iter

        ! initialize error flag
        error = 0

        ! unrestricted step
        call rfo_step(1.0_rp, step, step_norm)
        if (error /= 0) return
        if (step_norm <= step_bound) then
            step = -step
            return
        end if

        ! bracket the scaling that shortens the step to the step bound and bisect on it
        ! in logarithmic scale
        lower = 1.0_rp
        upper = 2.0_rp
        call rfo_step(upper, step, step_norm)
        do iter = 1, max_rfo_doublings
            if (step_norm <= step_bound .or. error /= 0) exit
            lower = upper
            upper = 2.0_rp * upper
            call rfo_step(upper, step, step_norm)
        end do
        if (error /= 0) return
        do iter = 1, max_rfo_bisections
            middle = sqrt(lower * upper)
            call rfo_step(middle, step, middle_norm)
            if (error /= 0) return
            if (abs(middle_norm - step_bound) <= step_bound_tol) exit
            if (middle_norm > step_bound) then
                lower = middle
            else
                upper = middle
            end if
        end do

        ! descent step
        step = -step

    contains

        subroutine rfo_step(alpha, rfo, rfo_norm)
            !
            ! this subroutine returns the RFO step of the augmented Hessian with the
            ! Hessian scaled by alpha and the negative gradient, and its norm
            !
            real(rp), intent(in) :: alpha
            real(rp), intent(out) :: rfo(:), rfo_norm

            integer(ip) :: n, i, lwork, info
            real(rp), allocatable :: aug_hess(:, :), eigvals(:), work(:)
            real(rp) :: last_component, trial(size(grad)), trial_norm
            logical :: found
            real(rp), external :: dnrm2
            external :: dsyev

            ! augmented Hessian and its eigenvectors
            n = size(grad, kind=ip)
            allocate(aug_hess(n + 1, n + 1), eigvals(n + 1), work(1))
            aug_hess = 0.0_rp
            aug_hess(:n, :n) = 0.5_rp * (hess + transpose(hess)) / alpha
            aug_hess(:n, n + 1) = -grad / sqrt(alpha)
            aug_hess(n + 1, :n) = -grad / sqrt(alpha)
            call dsyev("V", "U", n + 1, aug_hess, n + 1, eigvals, work, -1_ip, info)
            lwork = int(work(1), kind=ip)
            deallocate(work)
            allocate(work(lwork))
            call dsyev("V", "U", n + 1, aug_hess, n + 1, eigvals, work, lwork, info)
            if (info /= 0) then
                error = 1
                return
            end if

            ! shortest step of the lowest eigenvectors
            found = .false.
            rfo_norm = huge(1.0_rp)
            do i = 1, min(n_rfo_eigvecs, n + 1)
                last_component = aug_hess(n + 1, i)
                if (.not. abs(last_component) > 0.0_rp) cycle
                trial = aug_hess(:n, i) / sqrt(alpha) / last_component
                trial_norm = dnrm2(n, trial, 1_ip)
                if (trial_norm < rfo_norm) then
                    rfo = trial
                    rfo_norm = trial_norm
                    found = .true.
                end if
            end do
            if (.not. found) error = 1

        end subroutine rfo_step

    end subroutine rs_rfo

end module otr_s_gek
