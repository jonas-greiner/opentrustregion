! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion

    use, intrinsic :: iso_fortran_env, only: int32, int64, real64, &
                                             stdout => output_unit, stderr => error_unit

    implicit none

    integer, parameter :: rp = real64
#ifdef USE_ILP64
    integer, parameter :: ip = int64  ! 64-bit integers
#else
    integer, parameter :: ip = int32  ! 32-bit integers
#endif
    integer, parameter :: kw_len = 64

    ! mathematical constants
    real(rp), parameter :: pi = 4.0_rp * atan(1.0_rp)

    ! define trust region parameters
    real(rp), parameter :: default_spherical_trust_radius = 0.4_rp, &
                           default_ellipsoidal_trust_radius = 1.0_rp, &
                           trust_radius_shrink_ratio = 0.25_rp, &
                           trust_radius_expand_ratio = 0.75_rp, &
                           trust_radius_shrink_factor = 0.7_rp, &
                           trust_radius_expand_factor = 1.2_rp

    ! define error codes
    integer(ip), parameter :: error_solver = 100, error_stability_check = 200, &
                              error_obj_func = 1100, error_update_orbs = 1200, &
                              error_hess_x = 1300, error_precond = 1400, &
                              error_conv_check = 1500, error_project = 1600

    ! define useful parameters
    real(rp), parameter :: numerical_zero = 1e-14_rp, precond_floor = 1e-10_rp, &
                           hess_symm_thres = 1e-12_rp, residual_norm_floor = 1e-12_rp, &
                           level_shift_local_thres = 1e-12_rp

    ! define verbosity levels
    integer(ip), parameter :: verbosity_silent = 0, verbosity_error = 1, &
                              verbosity_warning = 2, verbosity_info = 3, &
                              verbosity_debug = 4

    ! define log messages
    character(len=*), parameter :: &
        random_trial_vector_warning_msg = &
            "Number of random trial vectors should be smaller than half the number "// &
            "of parameters.", &
        gram_schmidt_zero_vector_error_msg = &
            "Vector passed to Gram-Schmidt procedure is numerically zero.", &
        gram_schmidt_too_many_vectors_error_msg = &
            "Number of vectors in Gram-Schmidt procedure larger than dimension of "// &
            "vector space.", &
        project_warning_msg = &
            "Custom projection is provided. To optimize performance, OTR assumes "// &
            "that all other provided routines (update_orbs, hess_x, precond) are "// &
            "already projected onto the relevant orbital rotation subspace. If "// &
            "these routines are not self-projecting, redundant rotations may "// &
            "contaminate the trial space and cause convergence issues."

    ! interfaces for callback functions
    abstract interface
        subroutine hess_x_type(x, hess_x, error)
            import :: rp, ip

            real(rp), intent(in), target :: x(:)
            real(rp), intent(out), target :: hess_x(:)
            integer(ip), intent(out) :: error
        end subroutine hess_x_type
    end interface

    abstract interface
        subroutine update_orbs_type(kappa, func, grad, h_diag, hess_x_funptr, error)
            import :: rp, hess_x_type, ip

            real(rp), intent(in), target :: kappa(:)
            real(rp), intent(out) :: func
            real(rp), intent(out), target :: grad(:), h_diag(:)
            procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
            integer(ip), intent(out) :: error
        end subroutine update_orbs_type
    end interface

    abstract interface
        function obj_func_type(kappa, error) result(func)
            import :: rp, ip

            real(rp), intent(in), target :: kappa(:)
            integer(ip), intent(out) :: error
            real(rp) :: func
        end function obj_func_type
    end interface

    abstract interface
        subroutine precond_type(residual, mu, precond_residual, error)
            import :: rp, ip

            real(rp), intent(in), target :: residual(:)
            real(rp), intent(in) :: mu
            real(rp), intent(out), target :: precond_residual(:)
            integer(ip), intent(out) :: error
        end subroutine precond_type
    end interface

    abstract interface
        subroutine project_type(vector, error)
            import :: rp, ip

            real(rp), intent(inout), target :: vector(:)
            integer(ip), intent(out) :: error
        end subroutine project_type
    end interface

    abstract interface
        function conv_check_type(error) result(converged)
            import :: ip

            integer(ip), intent(out) :: error
            logical :: converged
        end function conv_check_type
    end interface

    abstract interface
        subroutine logger_type(message)
            character(*), intent(in) :: message
        end subroutine logger_type
    end interface

    ! derived type for solver settings
    type, abstract :: settings_type
        logical :: initialized = .false.
        integer(ip) :: verbose
        procedure(logger_type), pointer, nopass :: logger
    contains
        procedure :: log => print_message
        procedure(init_type), deferred :: init
    end type

    abstract interface
        subroutine init_type(self, error)
            import :: settings_type, ip

            class(settings_type), intent(out) :: self
            integer(ip), intent(out) :: error
        end subroutine init_type
    end interface

    type, abstract, extends(settings_type) :: optimizer_settings_type
        logical :: hess_symm
        real(rp) :: conv_tol
        integer(ip) :: n_random_trial_vectors, jacobi_davidson_start, seed
        procedure(precond_type), pointer, nopass :: precond
        procedure(project_type), pointer, nopass :: project
    end type

    type, extends(optimizer_settings_type) :: solver_settings_type
        logical :: stability, line_search
        real(rp) :: start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_macro, n_micro
        character(kw_len) :: subsystem_solver
        procedure(conv_check_type), pointer, nopass :: conv_check
    contains
        procedure :: init => init_solver_settings, print_results
    end type

    type, extends(optimizer_settings_type) :: stability_settings_type
        integer(ip) :: n_iter
        character(kw_len) :: diag_solver
    contains
        procedure :: init => init_stability_settings
    end type

    ! default settings
    type(solver_settings_type), parameter :: default_solver_settings = &
        solver_settings_type(precond = null(), project = null(), conv_check = null(), &
                             logger = null(), stability = .false., &
                             line_search = .false., hess_symm = .true., &
                             initialized = .true., conv_tol = 1e-5_rp, &
                             start_trust_radius = -1.0_rp, global_red_factor = 1e-3_rp, &
                             local_red_factor = 1e-4_rp, n_random_trial_vectors = 1, &
                             n_macro = 150, n_micro = 50, jacobi_davidson_start = 30, &
                             seed = 42, verbose = 0, subsystem_solver = "davidson")
    type(stability_settings_type), parameter :: default_stability_settings = &
        stability_settings_type(precond = null(), project = null(), logger = null(), &
                                hess_symm = .true., initialized = .true., &
                                conv_tol = 1e-8_rp, n_random_trial_vectors = 20, &
                                n_iter = 100, jacobi_davidson_start = 50, seed = 42, &
                                verbose = 0, diag_solver = "davidson")

    ! define global variables
    integer(ip) :: tot_orb_update = 0, tot_hess_x = 0

contains

    subroutine solver(update_orbs, obj_func, n_param, error, settings)
        !
        ! this subroutine is the main solver for orbital optimization
        !
        procedure(update_orbs_type), intent(in), pointer :: update_orbs
        procedure(obj_func_type), intent(in), pointer :: obj_func
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        type(solver_settings_type), intent(inout) :: settings

        type(stability_settings_type) :: stability_settings
        real(rp) :: trust_radius, func, grad_norm, grad_rms, mu, new_func, n_kappa, &
                    kappa_norm, lambda
        real(rp), allocatable :: kappa(:), grad(:), h_diag(:), precond_kappa(:)
        logical :: max_precision_reached, macro_converged, stable, &
                   jacobi_davidson_started, conv_check_passed
        integer(ip) :: imacro, imicro, imicro_jacobi_davidson, i
        character(300) :: msg
        integer(ip), parameter :: stability_n_points = 21
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = 0

        ! initialize settings
        if (.not. settings%initialized) then
            call settings%init(error)
            call add_error_origin(error, error_solver, settings)
            if (error /= 0) return
            call settings%log("Settings were not initialized. All settings are set "// &
                              "to default values", verbosity_warning)
        end if

        ! initialize maximum precision convergence
        max_precision_reached = .false.

        ! initialize macroiteration convergence
        macro_converged = .false.

        ! initialize stabilty boolean
        stable = .true.

        ! initialize random number generator
        call init_rng(settings%seed)

        ! initialize starting trust radius
        if (settings%start_trust_radius <= 0.0_rp) then
            if (settings%subsystem_solver == "davidson" .or. &
                settings%subsystem_solver == "jacobi-davidson") then
                trust_radius = default_spherical_trust_radius
            else
                trust_radius = default_ellipsoidal_trust_radius
            end if
        else
            trust_radius = settings%start_trust_radius
        end if

        ! print header
        call settings%log(repeat("-", 109), verbosity_info)
        call settings%log(" Iteration |     Objective function     | Gradient RMS |"// &
                          " Level shift |   Micro    | Trust radius | Step size ", &
                          verbosity_info)
        call settings%log("           |                            |              |"// &
                          "             | iterations |              |           ", &
                          verbosity_info)
        call settings%log(repeat("-", 109), verbosity_info)

        ! allocate arrays
        allocate(kappa(n_param), grad(n_param), h_diag(n_param), precond_kappa(n_param))

        ! initialize orbital rotation matrix
        kappa = 0.0_rp

        do imacro = 1, settings%n_macro
            if (.not. max_precision_reached) then
                ! calculate cost function, gradient and Hessian diagonal
                call update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error)
                call add_error_origin(error, error_update_orbs, settings)
                if (error /= 0) return

                ! perform sanity check
                if (imacro == 1) then
                    call solver_sanity_check(settings, n_param, grad, error)
                    call add_error_origin(error, error_solver, settings)
                    if (error /= 0) return
                end if

                ! calculate gradient norm
                grad_norm = dnrm2(n_param, grad, 1_ip)

                ! calculate RMS gradient
                grad_rms = grad_norm / sqrt(real(n_param, kind=rp))

                ! log results
                if (imacro == 1) then
                    call settings%print_results(imacro - 1, func, grad_rms)
                else
                    if (settings%subsystem_solver == "tcg") then
                        mu = 0.0_rp
                        jacobi_davidson_started = .false.
                    else if (settings%subsystem_solver == "gltr") then
                        mu = -lambda
                        jacobi_davidson_started = .false.
                    end if
                    if (.not. stable) then
                        call settings%print_results(imacro - 1, func, grad_rms, &
                                                    kappa_norm=kappa_norm)
                        stable = .true.
                    else if (jacobi_davidson_started) then
                        call settings%print_results(imacro - 1, func, grad_rms, &
                                                    level_shift=-mu, n_micro=imicro, &
                                                    imicro_jacobi_davidson=&
                                                    imicro_jacobi_davidson, &
                                                    trust_radius=trust_radius, &
                                                    kappa_norm=kappa_norm)
                    else
                        call settings%print_results(imacro - 1, func, grad_rms, &
                                                    level_shift=-mu, n_micro=imicro, &
                                                    trust_radius=trust_radius, &
                                                    kappa_norm=kappa_norm)
                    end if
                end if
            end if

            ! check for convergence and stability
            if (associated(settings%conv_check)) then
                conv_check_passed = settings%conv_check(error)
                call add_error_origin(error, error_conv_check, settings)
                if (error /= 0) return
            else
                conv_check_passed = .false.
            end if
            if (grad_rms < settings%conv_tol .or. max_precision_reached .or. &
                conv_check_passed) then
                ! always perform stability check if starting at stationary point
                if (settings%stability .or. imacro == 1) then
                    call stability_settings%init(error)
                    call add_error_origin(error, error_solver, settings)
                    if (error /= 0) return
                    stability_settings%precond => settings%precond
                    stability_settings%verbose = settings%verbose
                    stability_settings%logger => settings%logger
                    call stability_check(h_diag, hess_x_funptr, stable, error, &
                                         stability_settings, kappa=kappa)
                    call add_error_origin(error, error_stability_check, settings)
                    if (error /= 0) return
                    if (.not. stable) then
                        ! logarithmic line search
                        do i = 1, stability_n_points
                            n_kappa = 10.0_rp**(-(i - 1) / &
                                                real(stability_n_points - 1, rp) * &
                                                10.0_rp)
                            new_func = obj_func(n_kappa*kappa, error)
                            call add_error_origin(error, error_obj_func, settings)
                            if (error /= 0) return
                            if (new_func < func) then
                                kappa = n_kappa*kappa
                                exit
                            end if
                        end do
                        if (new_func >= func) then
                            call settings%log("Line search was unable to find "// &
                                              "lower objective function along "// &
                                              "unstable mode.", verbosity_error, .true.)
                            error = error_solver + 1
                            return
                        else if (imacro == 1) then
                            call settings%log("Started at saddle point. The "// &
                                              "algorithm will continue by moving "// &
                                              "along eigenvector direction "// &
                                              "corresponding to negative eigenvalue.", &
                                              verbosity_error, .true.)
                        else
                            call settings%log("Reached saddle point. This is "// &
                                              "likely due to symmetry and can be "// &
                                              "avoided by increasing the number of "// &
                                              "random trial vectors. The algorithm "// &
                                              "will continue by moving along "// &
                                              "eigenvector direction corresponding "// &
                                              "to negative eigenvalue.", &
                                              verbosity_error, .true.)
                        end if
                        max_precision_reached = .false.
                        cycle
                    else
                        macro_converged = .true.
                        exit
                    end if
                else
                    macro_converged = .true.
                    exit
                end if
            end if

            if (settings%subsystem_solver == "davidson" .or. &
                settings%subsystem_solver == "jacobi-davidson") then
                ! solve trust region subproblem with (Jacobi-)Davidson
                call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                            obj_func, hess_x_funptr, settings, &
                                            trust_radius, kappa, kappa_norm, mu, &
                                            imicro, imicro_jacobi_davidson, &
                                            jacobi_davidson_started, &
                                            max_precision_reached, error)
            else if (settings%subsystem_solver == "tcg") then
                ! solve trust region subproblem with truncated conjugate gradient
                call truncated_conjugate_gradient(func, grad, grad_norm, h_diag, &
                                                  n_param, obj_func, hess_x_funptr, &
                                                  settings, trust_radius, kappa, &
                                                  kappa_norm, imicro, &
                                                  max_precision_reached, error)
            else if (settings%subsystem_solver == "gltr") then
                ! solve trust region subproblem with generalized Lanczos
                call generalized_lanczos_trust_region(func, grad, grad_norm, h_diag, &
                                                      n_param, obj_func, &
                                                      hess_x_funptr, settings, &
                                                      trust_radius, kappa, kappa_norm, &
                                                      lambda, imicro, &
                                                      max_precision_reached, error)
            end if
            call add_error_origin(error, error_solver, settings)
            if (error /= 0) return

            ! perform line search
            if (max_precision_reached) then
                n_kappa = 0.0_rp
            else if (settings%line_search) then
                n_kappa = bracket(obj_func, kappa, 0.0_rp, 1.0_rp, settings, error)
                call add_error_origin(error, error_solver, settings)
                if (error /= 0) return
            else
                n_kappa = 1.0_rp
            end if

            ! set orbital rotation
            kappa = n_kappa * kappa
            kappa_norm = n_kappa * kappa_norm

            ! flush output
            flush(stdout)
            flush(stderr)

        end do

        ! deallocate arrays
        deallocate(kappa, grad, h_diag, precond_kappa)

        ! increment total number of orbital updates
        tot_orb_update = tot_orb_update + imacro

        ! stop if no convergence
        if (.not. macro_converged) then
            call settings%log("Orbital optimization has not converged!", &
                              verbosity_error, .true.)
            error = error_solver + 1
            return
        end if

        ! finish logging
        call settings%log(repeat("-", 109), verbosity_info)
        write (msg, '(A, I0)') "Total number of Hessian linear transformations: ", &
            tot_hess_x
        call settings%log(msg, verbosity_info)
        write (msg, '(A, I0)') "Total number of orbital updates: ", tot_orb_update
        call settings%log(msg, verbosity_info)

        ! reset global counter variables
        tot_orb_update = 0
        tot_hess_x = 0

        ! flush output
        flush (stdout)
        flush (stderr)

    end subroutine solver

    subroutine stability_check(h_diag, hess_x_funptr, stable, error, settings, kappa)
        !
        ! this subroutine performs a stability check
        !
        real(rp), intent(in) :: h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        logical, intent(out) :: stable
        integer(ip), intent(out) :: error
        type(stability_settings_type), intent(inout) :: settings
        real(rp), intent(out), optional :: kappa(:)

        integer(ip) :: n_param, n_trial, i, iter
        real(rp), allocatable :: solution(:), h_solution(:), residual(:), &
                                 basis_vec(:), h_basis_vec(:), red_space_basis(:, :), &
                                 h_basis(:, :), red_space_hess(:, :), &
                                 red_space_solution(:), row_vec(:), col_vec(:)
        real(rp) :: eigval, minres_tol, stability_rms
        character(300) :: msg
        real(rp), parameter :: stability_thresh = -1e-2_rp
        real(rp), external :: dnrm2, ddot
        external :: dgemm, dgemv

        ! initialize error flag
        error = 0

        ! initialize settings
        if (.not. settings%initialized) then
            call settings%init(error)
            call add_error_origin(error, error_stability_check, settings)
            if (error /= 0) return
            call settings%log("Settings were not initialized. All settings are set "// &
                              "to default values", verbosity_warning)
        end if

        ! initialize random number generator
        call init_rng(settings%seed)

        ! get number of parameters
        n_param = size(h_diag)

        ! perform sanity check
        call stability_sanity_check(settings, n_param, error)
        call add_error_origin(error, error_solver, settings)
        if (error /= 0) return

        ! generate trial vectors
        allocate(red_space_basis(n_param, 1 + settings%n_random_trial_vectors))
        red_space_basis(:, 1) = 0.0_rp
        red_space_basis(minloc(h_diag), 1) = 1.0_rp
        call generate_random_trial_vectors(red_space_basis, settings, error)
        call add_error_origin(error, error_stability_check, settings)
        if (error /= 0) return

        ! number of trial vectors
        n_trial = size(red_space_basis, 2)

        ! calculate linear transformations of basis vectors
        allocate(h_basis(n_param, n_trial))
        do i = 1, n_trial
            call hess_x_funptr(red_space_basis(:, i), h_basis(:, i), error)
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return
        end do

        ! increment number of Hessian linear transformations
        tot_hess_x = tot_hess_x + n_trial

        ! construct augmented Hessian in reduced space
        allocate(red_space_hess(n_trial, n_trial))
        call dgemm("T", "N", n_trial, n_trial, n_param, 1.0_rp, red_space_basis, &
                   n_param, h_basis, n_param, 0.0_rp, red_space_hess, n_trial)

        ! allocate arrays used throughout Davidson procedure
        allocate(red_space_solution(n_trial), solution(n_param), h_solution(n_param), &
                 residual(n_param), basis_vec(n_param), h_basis_vec(n_param))

        ! loop over iterations
        do iter = 1, settings%n_iter
            ! solve reduced space problem
            call mat_min_eig(red_space_hess, settings%hess_symm, eigval, &
                             red_space_solution, settings, error)
            call add_error_origin(error, error_stability_check, settings)
            if (error /= 0) return

            ! get full space solution
            call dgemv("N", n_param, n_trial, 1.0_rp, red_space_basis, n_param, &
                       red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)

            ! calculate Hessian linear transformation of solution
            call dgemv("N", n_param, n_trial, 1.0_rp, h_basis, n_param, &
                       red_space_solution, 1_ip, 0.0_rp, h_solution, 1_ip)

            ! calculate residual
            residual = h_solution - eigval*solution

            ! check convergence
            stability_rms = dnrm2(n_param, residual, 1_ip) / &
                sqrt(real(n_param, kind=rp))
            if (stability_rms < settings%conv_tol) exit

            if (settings%diag_solver == "davidson" .or. iter <= &
                settings%jacobi_davidson_start) then
                ! precondition residual
                call level_shifted_diag_precond(residual, 0.0_rp, h_diag, basis_vec, &
                                                settings, error)
                if (error /= 0) return

                ! orthonormalize to current orbital space to get new basis vector
                call gram_schmidt(basis_vec, red_space_basis, settings, error)
                call add_error_origin(error, error_stability_check, settings)
                if (error /= 0) return

                ! add linear transformation of new basis vector
                call hess_x_funptr(basis_vec, h_basis_vec, error)
                call add_error_origin(error, error_hess_x, settings)
                if (error /= 0) return

                ! increment Hessian linear transformations
                tot_hess_x = tot_hess_x + 1

            else
                ! solve Jacobi-Davidson correction equations
                minres_tol = 3.0_rp ** (-(iter - settings%jacobi_davidson_start - 1))
                call minres(-residual, hess_x_funptr, solution, eigval, minres_tol, &
                            basis_vec, h_basis_vec, settings, error)
                call add_error_origin(error, error_stability_check, settings)
                if (error /= 0) return

                ! orthonormalize to current orbital space to get new basis 
                ! vector
                call gram_schmidt(basis_vec, red_space_basis, settings, error, &
                                  lin_trans_vector=h_basis_vec, lin_trans_space=h_basis)
                call add_error_origin(error, error_stability_check, settings)
                if (error /= 0) return

                ! check if resulting linear transformation still respects Hessian 
                ! symmetry which can happen due to numerical noise accumulation
                if (abs(ddot(n_param, red_space_basis(:, n_trial), 1_ip, h_basis_vec, &
                             1_ip) &
                        - ddot(n_param, basis_vec, 1_ip, h_basis(:, n_trial), 1_ip)) > &
                    hess_symm_thres) then
                    call hess_x_funptr(basis_vec, h_basis_vec, error)
                    call add_error_origin(error, error_hess_x, settings)
                    if (error /= 0) return
                    tot_hess_x = tot_hess_x + 1
                end if
                
            end if

            ! increment trial vector count
            n_trial = n_trial + 1

            ! add new trial vector to orbital space
            call add_column(red_space_basis, basis_vec)

            ! add linear transformation of new basis vector
            call add_column(h_basis, h_basis_vec)

            ! construct new reduced space Hessian
            allocate(row_vec(n_trial))
            allocate(col_vec(n_trial))
            call dgemv("T", n_param, n_trial, 1.0_rp, red_space_basis, n_param, &
                       h_basis(:, n_trial), 1_ip, 0.0_rp, row_vec, 1_ip)
            if (settings%hess_symm) then
                col_vec = row_vec
            else
                call dgemv("T", n_param, n_trial, 1.0_rp, h_basis, n_param, &
                           red_space_basis(:, n_trial), 1_ip, 0.0_rp, col_vec, 1_ip)
            end if
            call extend_matrix(red_space_hess, row_vec, col_vec)
            deallocate(row_vec)
            deallocate(col_vec)

            ! reallocate reduced space solution
            deallocate(red_space_solution)
            allocate(red_space_solution(n_trial))

        end do

        ! check if stability check has converged
        if (stability_rms >= settings%conv_tol) &
            call settings%log("Stability check has not converged in the given "// &
                              "number of iterations.", verbosity_error, .true.)

        ! determine if saddle point
        stable = eigval > stability_thresh
        
        if (stable) then
            if (present(kappa)) kappa = 0.0_rp
        else
            if (present(kappa)) kappa = solution
            write (msg, '(A, F0.4)') "Solution not stable. Lowest eigenvalue: ", eigval
            call settings%log(msg, verbosity_error, .true.)
        end if

        ! deallocate quantities from Davidson iterations
        deallocate(solution, h_solution, residual, basis_vec, h_basis_vec, &
                   red_space_solution, red_space_hess, h_basis, red_space_basis)

        ! flush output
        flush (stdout)
        flush (stderr)

    end subroutine stability_check

    subroutine newton_step(aug_hess, grad_norm, red_space_basis, solution, &
                           red_space_solution, settings, error)
        !
        ! this subroutine performs a Newton step by solving the Newton equations in
        ! reduced space without a level shift
        !
        real(rp), intent(in) :: aug_hess(:, :), grad_norm, red_space_basis(:, :)
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(out) :: solution(:), red_space_solution(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, n_red, lwork, info
        integer(ip), allocatable :: ipiv(:)
        real(rp), allocatable :: red_hess(:, :), work(:)
        character(300) :: msg
        external :: dsysv, dgesv, dgemv

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(solution)

        ! reduced space size
        n_red = size(red_space_basis, 2)

        ! reduced space Hessian
        allocate(red_hess(n_red, n_red))
        red_hess = aug_hess(2:, 2:)

        ! set gradient
        red_space_solution = 0.0_rp
        red_space_solution(1) = -grad_norm

        ! allocate arrays
        allocate(ipiv(n_red))

        ! solve linear system
        if (settings%hess_symm) then
            ! query optimal workspace size
            lwork = -1
            allocate(work(1))
            call dsysv("U", n_red, 1_ip, red_hess, n_red, ipiv, red_space_solution, &
                       n_red, work, lwork, info)
            lwork = int(work(1), kind=ip)
            deallocate(work)
            allocate(work(lwork))

            ! solve linear system
            call dsysv("U", n_red, 1_ip, red_hess, n_red, ipiv, red_space_solution, &
                       n_red, work, lwork, info)

            ! deallocate work array
            deallocate(work)

            ! check for errors
            if (info /= 0) then
                write (msg, '(A, I0)') "Linear solver failed: Error in DSYSV, "// &
                                       "info = ", info
                call settings%log(msg, verbosity_error, .true.)
                error = 1
                return
            end if

        else
            ! general linear system solver since the Hessian can be non-symmetric for 
            ! Hessian approximations
            call dgesv(n_red, 1_ip, red_hess, n_red, ipiv, red_space_solution, n_red, &
                       info)

            ! check for errors
            if (info /= 0) then
                write (msg, '(A, I0)') "Linear solver failed: Error in DGESV, "// &
                                       "info = ", info
                call settings%log(msg, verbosity_error, .true.)
                error = 1
                return
            end if

        end if

        ! deallocate arrays
        deallocate(red_hess, ipiv)

        ! get solution in full space
        call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                   red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)

    end subroutine newton_step

    subroutine bisection(aug_hess, grad_norm, red_space_basis, trust_radius, solution, &
                         red_space_solution, mu, bracketed, settings, error)
        !
        ! this subroutine performs bisection to find the parameter alpha that matches
        ! the desired trust radius
        !
        real(rp), intent(inout) :: aug_hess(:, :)
        real(rp), intent(in) :: grad_norm, red_space_basis(:, :), trust_radius
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(out) :: solution(:), red_space_solution(:), mu
        logical, intent(out) :: bracketed
        integer(ip), intent(out) :: error

        real(rp) :: lower_alpha, middle_alpha, upper_alpha, lower_trust_dist, &
                    middle_trust_dist, upper_trust_dist
        integer(ip) :: n_param, n_red, iter
        real(rp), parameter :: lower_alpha_bound = 1e-4_rp, &
                               upper_alpha_bound = 1e6_rp, &
                               alpha_conv_factor = 1e-12_rp
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! initialize bracketing flag
        bracketed = .false.

        ! number of parameters
        n_param = size(solution)

        ! reduced space size
        n_red = size(red_space_basis, 2)

        ! lower and upper bracket for alpha
        lower_alpha = lower_alpha_bound
        upper_alpha = upper_alpha_bound

        ! solve reduced space problem with scaled gradient
        call get_ah_lowest_eigenvec(lower_alpha)
        if (error /= 0) return
        lower_trust_dist = dnrm2(n_param, solution, 1_ip) - trust_radius
        call get_ah_lowest_eigenvec(upper_alpha)
        if (error /= 0) return
        upper_trust_dist = dnrm2(n_param, solution, 1_ip) - trust_radius

        ! check if trust region is within bracketing range
        if ((lower_trust_dist*upper_trust_dist) > 0.0_rp) then
            solution = 0.0_rp
            red_space_solution = 0.0_rp
            mu = 0.0_rp
            return
        end if

        ! get middle alpha
        middle_alpha = sqrt(upper_alpha*lower_alpha)
        call get_ah_lowest_eigenvec(middle_alpha)
        if (error /= 0) return
        middle_trust_dist = dnrm2(n_param, solution, 1_ip) - trust_radius

        ! perform bisection to find root, converge to relative threshold to avoid 
        ! precision issues
        iter = 0
        do while (upper_alpha - lower_alpha > alpha_conv_factor * upper_alpha)
            ! targeted trust radius is in upper bracket
            if (lower_trust_dist*middle_trust_dist > 0.0_rp) then
                lower_alpha = middle_alpha
                lower_trust_dist = middle_trust_dist
                ! targeted trust radius is in lower bracket
            else
                upper_alpha = middle_alpha
                upper_trust_dist = middle_trust_dist
            end if
            ! get new middle alpha
            middle_alpha = sqrt(upper_alpha*lower_alpha)
            call get_ah_lowest_eigenvec(middle_alpha)
            if (error /= 0) return
            middle_trust_dist = dnrm2(n_param, solution, 1_ip) - trust_radius
            ! check if maximum number of iterations is reached
            iter = iter + 1
            if (iter > 100) then
                call settings%log("Maximum number of bisection iterations reached.", &
                                  verbosity_error, .true.)
                error = 1
                return
            end if
        end do

        bracketed = .true.

    contains

        subroutine get_ah_lowest_eigenvec(alpha)
            !
            ! this subroutine returns the lowest eigenvector for an augmented Hessian
            !
            real(rp), intent(in) :: alpha

            real(rp), allocatable :: eigvec(:)
            external :: dgemv

            ! finish construction of augmented Hessian
            aug_hess(1, 2) = alpha*grad_norm
            aug_hess(2, 1) = alpha*grad_norm

            ! allocate eigenvector
            allocate(eigvec(n_red + 1))

            ! perform eigendecomposition and get lowest eigenvalue and corresponding
            ! eigenvector
            call mat_min_eig(aug_hess, settings%hess_symm, mu, eigvec, settings, error)
            if (error /= 0) return

            ! check if eigenvector has level-shift component
            if (abs(eigvec(1)) <= numerical_zero) then
                call settings%log("Trial subspace too small. Increase "// &
                                  "n_random_trial_vectors.", verbosity_error, .true.)
                error = 1
                return
            end if

            ! scale eigenvector such that first element is equal to one and divide by
            ! alpha to get solution in reduced space
            red_space_solution = eigvec(2:)/eigvec(1)/alpha
            deallocate(eigvec)

            ! get solution in full space
            call dgemv("N", n_param, n_red, 1.0_rp, red_space_basis, n_param, &
                       red_space_solution, 1_ip, 0.0_rp, solution, 1_ip)

        end subroutine get_ah_lowest_eigenvec

    end subroutine bisection

    function bracket(obj_func, kappa, lower, upper, settings, error) result(n_kappa)
        !
        ! this function brackets a minimum (algorithm from numerical recipes)
        !
        procedure(obj_func_type), intent(in), pointer :: obj_func
        real(rp), intent(in) :: kappa(:), lower, upper
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp) :: n_kappa, f_upper, f_lower, n_a, n_b, n_c, n_u, f_a, f_b, f_c, f_u, &
                    n_u_lim, tmp1, tmp2, val, denom
        integer(ip) :: iter
        real(rp), parameter :: golden_ratio = (1.0_rp + sqrt(5.0_rp)) / 2.0_rp, &
                               grow_limit = 110.0_rp, denom_floor = 1e-21_rp

        ! initialize error flag
        error = 0

        ! initialize step length
        n_kappa = 0.0_rp

        ! evaluate function at upper and lower bounds
        f_lower = obj_func(lower*kappa, error)
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return
        f_upper = obj_func(upper*kappa, error)
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return

        ! ensure f_a > f_b
        if (f_upper > f_lower) then
            n_a = upper
            n_b = lower
            f_a = f_upper
            f_b = f_lower
        else
            n_a = lower
            n_b = upper
            f_a = f_lower
            f_b = f_upper
        end if

        ! default step
        n_c = n_b + golden_ratio*(n_b - n_a)
        f_c = obj_func(n_c*kappa, error)
        call add_error_origin(error, error_obj_func, settings)
        if (error /= 0) return

        ! continue looping until function at middle point is lower than at brackets
        iter = 0
        do while (f_c < f_b)
            ! compute value u by parabolic extrapolation
            tmp1 = (n_b - n_a)*(f_b - f_c)
            tmp2 = (n_b - n_c)*(f_b - f_a)
            val = tmp2 - tmp1
            denom = 2.0_rp*sign(max(abs(val), denom_floor), val)
            n_u = n_b - ((n_b - n_c)*tmp2 - (n_b - n_a)*tmp1)/denom

            ! maximum growth in parabolic fit
            n_u_lim = n_b + grow_limit*(n_c - n_b)

            ! check if u is between n_b and n_c
            if ((n_u - n_c)*(n_b - n_u) > 0.0) then
                ! evaluate function at n_u
                f_u = obj_func(n_u*kappa, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return

                ! check if minimum between n_b and n_c
                if (f_u < f_c) then
                    n_a = n_b
                    n_b = n_u
                    f_a = f_b
                    f_b = f_u
                    exit
                ! check if minium between n_a and n_u
                else if (f_u > f_b) then
                    n_c = n_u
                    f_c = f_u
                    exit
                end if

                ! parabolic fit did not help, default step
                n_u = n_c + golden_ratio*(n_c - n_b)
                f_u = obj_func(n_u*kappa, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return
            ! limit parabolic fit to its maximum allowed value
            else if ((n_u - n_u_lim)*(n_u_lim - n_c) >= 0.0_rp) then
                n_u = n_u_lim
                f_u = obj_func(n_u*kappa, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return
            ! parabolic fit is between n_c and its allowed limit
            else if ((n_u - n_u_lim)*(n_c - n_u) > 0.0_rp) then
                ! evaluate function at n_u
                f_u = obj_func(n_u*kappa, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return

                if (f_u < f_c) then
                    n_b = n_c
                    n_c = n_u
                    n_u = n_c + golden_ratio*(n_c - n_b)
                    f_b = f_c
                    f_c = f_u
                    f_u = obj_func(n_u*kappa, error)
                    call add_error_origin(error, error_obj_func, settings)
                    if (error /= 0) return
                end if
            ! reject parabolic fit and use default step
            else
                n_u = n_c + golden_ratio*(n_c - n_b)
                f_u = obj_func(n_u*kappa, error)
                call add_error_origin(error, error_obj_func, settings)
                if (error /= 0) return
            end if
            ! remove oldest point
            n_a = n_b
            n_b = n_c
            n_c = n_u
            f_a = f_b
            f_b = f_c
            f_c = f_u
            ! check if maximum number of iterations is reached
            iter = iter + 1
            if (iter > 100) then
                call settings%log("Maximum number of bracketing iterations reached.", &
                                  verbosity_error, .true.)
                error = 1
                return
            end if
        end do

        ! check if miniumum is bracketed
        if (.not. (((f_b < f_c .and. f_b <= f_a) .or. (f_b < f_a .and. f_b <= f_c)) &
                   .and. ((n_a < n_b .and. n_b < n_c) .or. &
                          (n_c < n_b .and. n_b < n_a)))) then
            call settings%log("Line search did not find minimum", verbosity_error, &
                              .true.)
            error = 1
            return
        end if

        ! set new multiplier
        n_kappa = n_b

    end function bracket

    subroutine extend_matrix(matrix, row, column)
        !
        ! this subroutine extends a matrix with a row and a column
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: row(:), column(:)

        real(rp), allocatable :: temp(:, :)

        ! copy data
        temp = matrix

        ! reallocate matrix
        deallocate(matrix)
        allocate(matrix(size(temp, 1) + 1, size(temp, 2) + 1))

        ! copy data back
        matrix(:size(temp, 1), :size(temp, 2)) = temp

        ! add new row and column
        matrix(:, size(temp, 2) + 1) = row
        matrix(size(temp, 1) + 1, :) = column

        ! deallocate temporary array
        deallocate(temp)

    end subroutine extend_matrix

    subroutine add_column(matrix, new_col)
        !
        ! this subroutine adds a new column to an array
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: new_col(:)

        real(rp), allocatable :: temp(:, :)

        ! copy data
        temp = matrix

        ! reallocate matrix
        deallocate(matrix)
        allocate(matrix(size(temp, 1), size(temp, 2) + 1))

        ! copy data back
        matrix(:, :size(temp, 2)) = temp

        ! add new column
        matrix(:, size(temp, 2) + 1) = new_col

        ! deallocate temporary array
        deallocate(temp)

    end subroutine add_column

    subroutine mat_min_eig(matrix, symm_matrix, lowest_eigval, lowest_eigvec, &
                           settings, error)
        !
        ! this subroutine returns the lowest eigenvalue and corresponding eigenvector 
        ! of a matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        logical, intent(in) :: symm_matrix
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: lowest_eigval, lowest_eigvec(:)
        integer(ip), intent(out) :: error

        if (symm_matrix) then
            call symm_mat_min_eig(matrix, lowest_eigval, lowest_eigvec, settings, error)
        else
            call general_mat_min_eig(matrix, lowest_eigval, lowest_eigvec, settings, &
                                     error)
        end if 

    end subroutine mat_min_eig

    subroutine symm_mat_min_eig(symm_matrix, lowest_eigval, lowest_eigvec, settings, &
                                error)
        !
        ! this subroutine returns the lowest eigenvalue and corresponding eigenvector 
        ! of a symmetric matrix
        !
        real(rp), intent(in) :: symm_matrix(:, :)
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: lowest_eigval, lowest_eigvec(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n
        real(rp), allocatable :: eigvals(:), eigvecs(:, :)

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(symm_matrix, 1)

        ! allocate arrays
        allocate(eigvals(n), eigvecs(n, n))

        ! diagonalize matrix
        call symm_diag(symm_matrix, eigvals, eigvecs, settings, error)
        if (error /= 0) return

        ! get lowest eigenvalue and corresponding eigenvector
        lowest_eigval = eigvals(1)
        lowest_eigvec = eigvecs(:, 1)

        ! deallocate eigenvalues and eigenvectors
        deallocate(eigvals, eigvecs)

    end subroutine symm_mat_min_eig

    subroutine symm_diag(symm_matrix, eigvals, eigvecs, settings, error)
        !
        ! this subroutine returns eigenvalues and eigenvectors of a symmetric matrix
        !
        real(rp), intent(in) :: symm_matrix(:, :)
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: eigvals(:), eigvecs(:, :)
        integer(ip), intent(out) :: error

        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:)
        character(300) :: msg
        external :: dsyev

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(symm_matrix, 1)

        ! copy matrix
        eigvecs = symm_matrix

        ! query optimal workspace size
        lwork = -1
        allocate(work(1))
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)

        ! deallocate work array
        deallocate(work)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine symm_diag

    subroutine general_mat_min_eig(matrix, lowest_eigval, lowest_eigvec, settings, &
                                   error)
        !
        ! this subroutine returns the lowest eigenvalue and corresponding eigenvector 
        ! of a square matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        class(settings_type), intent(in) :: settings
        real(rp), intent(out) :: lowest_eigval, lowest_eigvec(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n, lwork, info, min_idx
        real(rp), allocatable :: work(:), temp(:, :), eigvals(:), imag_eigvals(:), &
                                 left_eigvecs(:, :), right_eigvecs(:, :)
        character(300) :: msg
        external :: dgeev

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(matrix, 1)

        ! copy matrix to avoid modification of original matrix
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate(eigvals(n), imag_eigvals(n), left_eigvecs(n, n), right_eigvecs(n, n), &
                 work(1))
        call dgeev("N", "V", n, temp, n, eigvals, imag_eigvals, left_eigvecs, n, &
                   right_eigvecs, n, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! perform eigendecomposition
        call dgeev("N", "V", n, temp, n, eigvals, imag_eigvals, left_eigvecs, n, &
                   right_eigvecs, n, work, lwork, info)

        ! deallocate arrays
        deallocate(imag_eigvals, left_eigvecs, work)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DGEEV, "// &
                "info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! get lowest eigenvalue and corresponding eigenvector
        min_idx = minloc(eigvals, dim=1)
        lowest_eigval = eigvals(min_idx)
        lowest_eigvec = right_eigvecs(:, min_idx)

        ! deallocate eigenvalues and eigenvectors
        deallocate(eigvals, right_eigvecs)

    end subroutine general_mat_min_eig

    real(rp) function mat_min_eigval(matrix, symm_matrix, settings, error)
        !
        ! this function calculates the lowest eigenvalue of a matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        logical, intent(in) :: symm_matrix
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        if (symm_matrix) then
            mat_min_eigval = symm_mat_min_eigval(matrix, settings, error)
        else
            mat_min_eigval = general_mat_min_eigval(matrix, settings, error)
        end if

    end function mat_min_eigval

    real(rp) function symm_mat_min_eigval(matrix, settings, error)
        !
        ! this function calculates the lowest eigenvalue of a symmetric matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: eigvals(:), temp(:, :), work(:)
        integer(ip) :: n, lwork, info
        character(300) :: msg
        external :: dsyev

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(matrix, 1)

        ! copy matrix to avoid modification of original matrix
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate(eigvals(n), work(1))
        call dsyev("N", "U", n, temp, n, eigvals, work, lwork, info)
        lwork = int(work(1), kind=ip)
        deallocate(work)
        allocate(work(lwork))

        ! compute eigenvalues
        call dsyev("N", "U", n, temp, n, eigvals, work, lwork, info)

        ! deallocate temporary and work array
        deallocate(temp, work)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DSYEV, "// &
                "info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! get lowest eigenvalue
        symm_mat_min_eigval = eigvals(1)

        ! deallocate eigenvalues
        deallocate(eigvals)

    end function symm_mat_min_eigval

    real(rp) function general_mat_min_eigval(matrix, settings, error)
        !
        ! this function calculates the lowest eigenvalue of a square matrix
        !
        real(rp), intent(in) :: matrix(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: temp(:, :), eigvals(:), imag_eigvals(:), &
                                 left_eigvecs(:, :), right_eigvecs(:, :)
        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:)
        character(300) :: msg
        external :: dgeev

        ! initialize error flag
        error = 0

        ! size of matrix
        n = size(matrix, 1)

        ! copy matrix to avoid modification of original matrix
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate(eigvals(n), imag_eigvals(n), left_eigvecs(n, n), right_eigvecs(n, n), &
                 work(1))
        call dgeev("N", "N", n, temp, n, eigvals, imag_eigvals, left_eigvecs, n, &
                   right_eigvecs, n, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! compute eigenvalues
        call dgeev("N", "N", n, temp, n, eigvals, imag_eigvals, left_eigvecs, n, &
                   right_eigvecs, n, work, lwork, info)

        ! deallocate arrays
        deallocate(imag_eigvals, left_eigvecs, right_eigvecs, work)

        ! check for successful execution
        if (info /= 0) then
            write (msg, '(A, I0)') "Eigendecomposition failed: Error in DGEEV, "// &
                "info = ", info
            call settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! get lowest eigenvalue
        general_mat_min_eigval = minval(eigvals)

        ! deallocate eigenvalues and eigenvectors
        deallocate(eigvals)

    end function general_mat_min_eigval

    subroutine init_rng(seed)
        !
        ! this subroutine initializes the random number generator
        !
        integer(ip), intent(in) :: seed

        integer, allocatable :: seed_arr(:)
        integer :: n

        call random_seed(size=n)
        allocate(seed_arr(n))
        seed_arr = int(seed, kind=4)
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

    end subroutine init_rng

    function generate_trial_vectors(grad, grad_norm, h_diag, settings, error) &
        result(red_space_basis)
        !
        ! this function generates trial vectors
        !
        real(rp), intent(in) :: grad(:), grad_norm, h_diag(:)
        type(solver_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        real(rp), allocatable :: red_space_basis(:, :)

        integer(ip) :: min_idx, n_vectors
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! get minimum Hessian diagonal element
        min_idx = minloc(h_diag, dim=1)

        ! add direction if minimum Hessian diagonal element is negative
        if (h_diag(min_idx) < 0.0_rp .and. size(grad) > 2) then
            n_vectors = 2
            allocate(red_space_basis(size(grad), n_vectors + &
                     settings%n_random_trial_vectors))
            red_space_basis(:, 1) = grad/grad_norm
            red_space_basis(:, 2) = 0.0_rp
            red_space_basis(min_idx, 2) = 1.0_rp
            if (associated(settings%project)) then
                call settings%project(red_space_basis(:, 2), error)
                call add_error_origin(error, error_project, settings)
                if (error /= 0) return
            end if
            call gram_schmidt(red_space_basis(:, 2), &
                              reshape(red_space_basis(:, 1), &
                                      [size(red_space_basis, 1), 1]), &
                              settings, error)
            if (error /= 0) return
        else
            n_vectors = 1
            allocate(red_space_basis(size(grad), n_vectors + &
                     settings%n_random_trial_vectors))
            red_space_basis(:, 1) = grad/grad_norm
        end if

        call generate_random_trial_vectors(red_space_basis, settings, error)

    end function generate_trial_vectors

    subroutine generate_random_trial_vectors(red_space_basis, settings, error)
        !
        ! this subroutine generates random trial vectors
        !
        real(rp), intent(inout) :: red_space_basis(:, :)
        class(optimizer_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_param, n_trial, i
        real(rp), parameter :: rnd_vector_min_norm = 1e-3_rp
        real(rp), external :: dnrm2

        ! initialize error flag
        error = 0

        ! number of parameters
        n_param = size(red_space_basis, 1)

        ! number of trial vectors
        n_trial = size(red_space_basis, 2)

        do i = n_trial - settings%n_random_trial_vectors + 1, n_trial
            call random_number(red_space_basis(:, i))
            red_space_basis(:, i) = 2*red_space_basis(:, i) - 1
            do while (dnrm2(n_param, red_space_basis(:, i), 1_ip) < rnd_vector_min_norm)
                call random_number(red_space_basis(:, i))
                red_space_basis(:, i) = 2*red_space_basis(:, i) - 1
            end do
            if (associated(settings%project)) then
                call settings%project(red_space_basis(:, i), error)
                call add_error_origin(error, error_project, settings)
                if (error /= 0) return
            end if
            call gram_schmidt(red_space_basis(:, i), red_space_basis(:, :i - 1), &
                              settings, error)
            if (error /= 0) return
        end do

    end subroutine generate_random_trial_vectors

    subroutine gram_schmidt(vector, space, settings, error, lin_trans_vector, &
                            lin_trans_space)
        !
        ! this function orthonormalizes a vector with respect to a vector space
        ! this function can additionally also return a linear transformation of the 
        ! orthogonalized vector if the linear transformations of the vector and the 
        ! vector space are provided
        !
        real(rp), intent(inout) :: vector(:)
        real(rp), intent(in) :: space(:, :)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        real(rp), intent(inout), optional :: lin_trans_vector(:)
        real(rp), intent(in), optional :: lin_trans_space(:, :)

        real(rp), allocatable :: orth(:)
        real(rp) :: dot, norm
        integer(ip) :: n_param, n_vectors, iter, i
        real(rp), parameter :: zero_thres = 1e-16_rp, orth_thres = 1e-14_rp
        real(rp), external :: ddot, dnrm2
        external :: dgemv

        ! initialize error flag
        error = 0

        ! vector size
        n_param = size(vector)

        ! number of vectors
        n_vectors = size(space, 2)

        if (dnrm2(n_param, vector, 1_ip) < zero_thres) then
            call settings%log(gram_schmidt_zero_vector_error_msg, verbosity_error, &
                              .true.)
            error = 1
            return
        else if (n_vectors > n_param - 1) then
            call settings%log(gram_schmidt_too_many_vectors_error_msg, &
                              verbosity_error, .true.)
            error = 1
            return
        end if
        
        ! allocate array for orthogonalities
        allocate(orth(size(space, 2)))

        iter = 0
        if (.not. (present(lin_trans_vector) .and. present(lin_trans_space))) then
            do while (.true.)
                do i = 1, n_vectors
                    vector = orthogonal_projection(vector, space(:, i))
                end do
                vector = vector / dnrm2(n_param, vector, 1_ip)

                call dgemv("T", n_param, n_vectors, 1.0_rp, space, n_param, &
                           vector, 1_ip, 0.0_rp, orth, 1_ip)
                if (maxval(abs(orth)) < orth_thres) exit
                iter = iter + 1
                if (iter > 100) then
                    call settings%log("Maximum number of Gram-Schmidt iterations "// &
                                      "reached.", verbosity_error, .true.)
                    error = 1
                    return
                end if
            end do
        else
            do while (.true.)
                do i = 1, n_vectors
                    dot = ddot(n_param, vector, 1_ip, space(:, i), 1_ip)
                    vector = vector - dot * space(:, i)
                    lin_trans_vector = lin_trans_vector - dot * lin_trans_space(:, i)
                end do
                norm = dnrm2(n_param, vector, 1_ip)
                vector = vector / norm
                lin_trans_vector = lin_trans_vector / norm

                call dgemv("T", n_param, n_vectors, 1.0_rp, space, n_param, vector, &
                           1_ip, 0.0_rp, orth, 1_ip)
                if (maxval(abs(orth)) < orth_thres) exit
                iter = iter + 1
                if (iter > 100) then
                    call settings%log("Maximum number of Gram-Schmidt iterations "// &
                                      "reached.", verbosity_error, .true.)
                    error = 1
                    return
                end if
            end do
        end if

        ! allocate array for orthogonalities
        deallocate(orth)

    end subroutine gram_schmidt

    subroutine init_solver_settings(self, error)
        !
        ! this subroutine sets the optional settings to their default values
        !
        class(solver_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! ensure that class is actually solver_settings_type and not a subclass
        select type(settings => self)
        type is (solver_settings_type)
            settings = default_solver_settings
        class default
            call settings%log("Solver settings could not be initialized because "// &
                              "initialization routine received the wrong type. The "// &
                              "type solver_settings_type was likely subclassed "// &
                              "without providing an initialization routine.", &
                              verbosity_error, .true.)
            error = 1
        end select

    end subroutine init_solver_settings

    subroutine init_stability_settings(self, error)
        !
        ! this subroutine sets the optional settings to their default values
        !
        class(stability_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! ensure that class is actually stability_settings_type and not a subclass
        select type(settings => self)
        type is (stability_settings_type)
            settings = default_stability_settings
        class default
            call settings%log("Stability settings could not be initialized because "// &
                              "initialization routine received the wrong type. The "// &
                              "type stability_settings_type was likely subclassed "// &
                              "without providing an initialization routine.", &
                              verbosity_error, .true.)
            error = 1
        end select

    end subroutine init_stability_settings

    subroutine level_shifted_diag_precond(vector, mu, h_diag, precond_vector, &
                                          settings, error)
        !
        ! this function defines the default level-shifted diagonal preconditioner
        !
        real(rp), intent(in) :: vector(:), mu, h_diag(:)
        real(rp), intent(out) :: precond_vector(:)
        class(optimizer_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! check for user-defined preconditioner
        if (associated(settings%precond)) then
            call settings%precond(vector, mu, precond_vector, error)
            call add_error_origin(error, error_precond, settings)
            if (error /= 0) return
        ! construct level-shifted preconditioner
        else
            precond_vector = h_diag - mu
            where (abs(precond_vector) < precond_floor)
                precond_vector = precond_floor
            end where
            precond_vector = vector / precond_vector

            ! ensure basis vector stays in subspace
            if (associated(settings%project)) then
                call settings%project(precond_vector, error)
                call add_error_origin(error, error_project, settings)
                if (error /= 0) return
            end if
        end if
        
    end subroutine level_shifted_diag_precond

    subroutine abs_diag_precond(vector, h_diag, precond_vector, settings, error)
        !
        ! this function defines the default absolute diagonal preconditioner
        !
        real(rp), intent(in) :: vector(:), h_diag(:)
        real(rp), intent(out) :: precond_vector(:)
        class(optimizer_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! check for user-defined preconditioner
        if (associated(settings%precond)) then
            call settings%precond(vector, 0.0_rp, precond_vector, error)
            call add_error_origin(error, error_precond, settings)
            if (error /= 0) return
        ! construct positive-definite preconditioner
        else
            precond_vector = max(abs(h_diag), precond_floor)
            precond_vector = vector / precond_vector

            ! ensure basis vector stays in subspace
            if (associated(settings%project)) then
                call settings%project(precond_vector, error)
                call add_error_origin(error, error_project, settings)
                if (error /= 0) return
            end if
        end if
        
    end subroutine abs_diag_precond

    function orthogonal_projection(vector, direction) result(complement)
        !
        ! this function removes a certain direction from a vector
        !
        real(rp), intent(in) :: vector(:), direction(:)

        real(rp) :: complement(size(vector))

        real(rp), external :: ddot

        complement = vector - direction * &
                     ddot(int(size(vector), kind=ip), vector, 1_ip, direction, 1_ip)

    end function orthogonal_projection

    subroutine jacobi_davidson_correction(hess_x_funptr, vector, solution, eigval, &
                                          corr_vector, hess_vector, settings, error)
        !
        ! this subroutine performs the Jacobi-Davidson correction but also returns the 
        ! Hessian linear transformation since this can be reused
        !
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(in) :: vector(:), solution(:), eigval
        real(rp), intent(out) :: corr_vector(:), hess_vector(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0
        
        ! project solution out of vector
        corr_vector = orthogonal_projection(vector, solution)

        ! get Hessian linear transformation of projected vector
        call hess_x_funptr(corr_vector, hess_vector, error)
        call add_error_origin(error, error_hess_x, settings)
        if (error /= 0) return

        ! finish construction of correction
        corr_vector = orthogonal_projection(hess_vector - eigval * corr_vector, &
                                            solution)
    
    end subroutine jacobi_davidson_correction

    subroutine minres(rhs, hess_x_funptr, solution, eigval, r_tol, vec, hvec, &
                      settings, error, guess, max_iter)
        !
        ! this function uses the minimum residual method to iteratively solve the 
        ! linear system for the Jacobi-Davidson correction equation, modified from 
        ! SciPy implementation
        !
        real(rp), intent(in) :: rhs(:), r_tol, solution(:), eigval
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(out) :: vec(:), hvec(:)
        class(settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        real(rp), intent(in), optional :: guess(:)
        integer(ip), intent(in), optional :: max_iter

        integer(ip) :: n, max_iterations, iteration
        real(rp), parameter :: eps = epsilon(1.0_rp)
        real(rp) :: beta_start, beta, phi_bar, rhs1, old_beta, alfa, t_norm2, eps_ln, &
                    old_eps, cs, d_bar, sn, delta, g_bar, root, gamma, phi, g_max, &
                    g_min, tmp, rhs2, a_norm, vec_norm, qr_norm
        real(rp), allocatable :: matvec(:), r1(:), r2(:), y(:), w(:), hw(:), w1(:), &
                                 hw1(:), w2(:), hw2(:), v(:), hv(:)
        logical :: stop_iteration
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = 0

        ! initialze boolean to stop iterations
        stop_iteration = .false.

        ! size of problem
        n = size(rhs)

        ! allocate solution vector
        allocate(matvec(n), r1(n), r2(n), y(n), w(n), hw(n), w1(n), hw1(n), w2(n), &
                 hw2(n), v(n), hv(n))

        ! initial guess
        if (present(guess)) then
            vec = guess
            call jacobi_davidson_correction(hess_x_funptr, vec, solution, eigval, &
                                            matvec, hvec, settings, error)
            if (error /= 0) return
            tot_hess_x = tot_hess_x + 1
        else
            vec = 0.0_rp
            hvec = 0.0_rp
            matvec = 0.0_rp
        end if

        ! maximum number of iterations
        if (present(max_iter)) then
            max_iterations = max_iter
        else
            max_iterations = 5 * n
        end if

        r1 = rhs - matvec
        y = r1

        beta_start = dnrm2(n, r1, 1_ip)

        ! check if starting guess already describes solution
        if (beta_start < numerical_zero) return

        ! solution must be zero vector if rhs vanishes
        if (dnrm2(n, rhs, 1_ip) < numerical_zero) then
            vec = 0.0_rp
            hvec = 0.0_rp
            return
        end if

        ! initialize additional quantities
        beta = beta_start
        r2 = r1
        phi_bar = beta_start
        rhs1 = beta_start
        t_norm2 = 0.0_rp
        eps_ln = 0.0_rp
        cs = -1.0_rp
        d_bar = 0.0_rp
        sn = 0.0_rp
        w = 0.0_rp
        hw = 0.0_rp
        w2 = 0.0_rp
        hw2 = 0.0_rp
        g_max = 0.0_rp
        g_min = huge(1.0_rp)
        rhs2 = 0.0_rp

        iteration = 1
        do while (iteration <= max_iterations)
            ! scale trial vector
            v = y / beta

            ! apply Jacobi-Davidson projector to trial vector
            call jacobi_davidson_correction(hess_x_funptr, v, solution, eigval, y, hv, &
                                            settings, error)
            if (error /= 0) return
            tot_hess_x = tot_hess_x + 1

            ! get new trial vector
            if (iteration >= 2) y = y - (beta / old_beta) * r1
            alfa = ddot(n, v, 1_ip, y, 1_ip)
            y = y - (alfa / beta) * r2
            r1 = r2
            r2 = y
            old_beta = beta
            beta = dnrm2(n, r2, 1_ip)
            t_norm2 = t_norm2 + alfa**2 + old_beta**2 + beta**2

            if (iteration == 1 .and. beta / beta_start <= 10 * eps) &
                stop_iteration = .true.

            ! apply previous plane rotation
            ! [delta_k epsln_k+1] = [cs  sn][d_bar_k     0   ]
            ! [g_bar_k d_bar_k+1]   [sn -cs][alfa_k  beta_k+1]
            old_eps = eps_ln
            delta = cs * d_bar + sn * alfa
            g_bar = sn * d_bar - cs * alfa
            eps_ln = sn * beta
            d_bar = -cs * beta
            root = sqrt(g_bar ** 2 + d_bar ** 2)

            ! compute next plane rotation by calculating cs and sn
            gamma = max(sqrt(g_bar ** 2 + beta ** 2), eps)
            cs = g_bar / gamma
            sn = beta / gamma
            phi = cs * phi_bar
            phi_bar = sn * phi_bar

            ! update vec and hvec
            w1 = w2
            w2 = w
            w = (v - old_eps * w1 - delta * w2) / gamma
            hw1 = hw2
            hw2 = hw
            hw = (hv - old_eps * hw1 - delta * hw2) / gamma
            vec = vec + phi * w
            hvec = hvec + phi * hw

            ! round variables
            g_max = max(g_max, gamma)
            g_min = min(g_min, gamma)
            tmp = rhs1 / gamma
            rhs1 = rhs2 - delta * tmp
            rhs2 = -eps_ln * tmp

            ! estimate various norms and test for convergence
            a_norm = sqrt(t_norm2)
            vec_norm = dnrm2(n, vec, 1_ip)
            qr_norm = phi_bar

            ! check if rhs and initial vector are eigenvectors
            if (stop_iteration) then
                call settings%log("MINRES: beta2 = 0. If M = I, b and x are "// &
                                  "eigenvectors.", verbosity_debug)
                exit
            ! ||r||  / (||A|| ||x||)
            else if (vec_norm > 0.0_rp .and. a_norm > 0.0_rp .and. &
                qr_norm / (a_norm * vec_norm) <= r_tol) then
                    call settings%log("MINRES: A solution to Ax = b was found, "// &
                                      "given provided tolerance.", verbosity_debug)
                exit
            ! ||Ar|| / (||A|| ||r||)
            else if (a_norm < numerical_zero .and. root / a_norm <= r_tol) then
                call settings%log("MINRES: A least-squares solution was found, "// &
                                  "given provided tolerance.", verbosity_debug)
                exit
            ! check if reasonable accuracy is achieved with respect to machine precision
            else if (a_norm * vec_norm * eps >= beta_start) then
                call settings%log("MINRES: Reasonable accuracy achieved, given "// &
                                  "machine precision.", verbosity_debug)
                exit
            ! compare estimate of condition number of matrix to machine precision
            else if (g_max / g_min >= 0.1 / eps) then
                call settings%log("MINRES: x has converged to an eigenvector.", &
                                  verbosity_debug)
                exit
            ! check if maximum number of iterations has been reached
            else if (iteration == max_iterations) then
                call settings%log("MINRES: The iteration limit was reached.", &
                                  verbosity_error, .true.)
                error = 1
                return
            ! these tests ensure convergence is still achieved when r_tol 
            ! approaches machine precision
            else if (vec_norm > 0.0_rp .and. a_norm > 0.0_rp .and. &
                     1.0_rp + qr_norm / (a_norm * vec_norm) <= 1.0_rp) then
                call settings%log("MINRES: A solution to Ax = b was found, given "// &
                                  "provided tolerance.", verbosity_debug)
                exit
            else if (a_norm < numerical_zero .and. 1.0_rp + root / a_norm <= 1.0_rp) &
                then
                call settings%log("MINRES: A least-squares solution was found, "// &
                                  "given provided tolerance.", verbosity_debug)
                exit
            end if

            ! increment iteration counter
            iteration = iteration + 1

        end do

        deallocate(matvec, r1, r2, y, w, hw, w1, hw1, w2, hw2, v, hv)

    end subroutine minres

    subroutine print_results(self, iteration, func, grad_rms, level_shift, n_micro, &
                             imicro_jacobi_davidson, trust_radius, kappa_norm)
        !
        ! this function prints rows of the result table
        !
        class(solver_settings_type), intent(in) :: self
        integer(ip), intent(in) :: iteration
        real(rp), intent(in) :: func, grad_rms
        real(rp), intent(in), optional :: level_shift, trust_radius, kappa_norm
        integer(ip), intent(in), optional :: n_micro, imicro_jacobi_davidson

        character(11) :: iteration_str
        character(12) :: n_micro_str
        character(28) :: func_str
        character(14) :: grad_rms_str, trust_radius_str
        character(13) :: level_shift_str
        character(10) :: kappa_norm_str

        write(iteration_str, '(4X, I4, 3X)') iteration
        write(func_str, '(4X, 1PE21.14, 3X)') func
        write(grad_rms_str, '(2X, 1PE9.2, 3X)') grad_rms

        if (present(level_shift)) then
            write(level_shift_str, '(2X, 1PE9.2, 2X)') level_shift
        else
            write(level_shift_str, '(6X, "-", 6X)')
        end if

        if (present(n_micro)) then
            if (present(imicro_jacobi_davidson)) then
                write(n_micro_str, '(X, I3, X, "|", X, I3, 2X)') &
                    n_micro - imicro_jacobi_davidson + 1, imicro_jacobi_davidson - 1
            else
                write(n_micro_str, '(7X, I3, 2X)') n_micro
            end if
        else
            write(n_micro_str, '(8X, "-", 3X)')
        end if

        if (present(trust_radius)) then
            write(trust_radius_str, '(3X, 1PE9.2, 2X)') trust_radius
        else
            write(trust_radius_str, '(8X, "-", 5X)')
        end if

        if (present(kappa_norm)) then
            write(kappa_norm_str, '(X, 1PE9.2)') kappa_norm
        else
            write(kappa_norm_str, '(6X, "-", 3X)')
        end if

        call self%log(iteration_str // "|" // func_str // "|" // grad_rms_str &
                      // "|" // level_shift_str // "|" // n_micro_str // "|" // &
                      trust_radius_str // "|" // kappa_norm_str, verbosity_info)

    end subroutine print_results

    subroutine print_message(self, message, level, error)
        !
        ! this function performs logging
        !
        class(settings_type), intent(in) :: self
        character(*), intent(in) :: message
        integer(ip), intent(in) :: level
        logical, intent(in), optional :: error

        integer(ip), parameter :: max_length = 109
        integer(ip) :: i, out_unit
        character(:), dimension(:), allocatable :: substrings

        if (self%verbose >= level) then
            call split_string_by_space(message, max_length, substrings)
            if (associated(self%logger)) then
                do i = 1, size(substrings)
                    call self%logger(" " // substrings(i))
                end do
            else
                if (.not. present(error)) then
                    out_unit = stdout
                else if (.not. error) then
                    out_unit = stdout
                else
                    out_unit = stderr
                end if
                do i = 1, size(substrings)
                    write(out_unit, '(A)') " " // substrings(i)
                end do
            end if
            deallocate(substrings)
        end if

    end subroutine print_message

    subroutine split_string_by_space(input, max_length, substrings)
        !
        ! this function splits a string by spaces to produce substrings of a maximum 
        ! length
        !
        character(*), intent(in) :: input
        integer(ip), intent(in) :: max_length
        character(:), intent(out), allocatable :: substrings(:)
    
        integer(ip) :: i, len_input, start_pos, end_pos, space_pos
        
        character(max_length) :: temp_string
        character(max_length), allocatable :: temp_substrings(:)

        len_input = len_trim(input)
        start_pos = 1
    
        do while (start_pos <= len_input)
            end_pos = min(start_pos + max_length - 1, len_input)
            space_pos = 0
    
            ! find last space before max_length characters
            do i = end_pos, start_pos, -1
                if (input(i:i) == ' ') then
                    space_pos = i
                    exit
                end if
            end do
    
            ! force splitting if no space is found
            if (end_pos == len_input) then
                space_pos = len_input
            else if (space_pos == 0) then
                space_pos = end_pos
            end if
    
            ! add substring
            temp_string = " "
            temp_string(1:space_pos-start_pos+1) = input(start_pos:space_pos)
            if (allocated(substrings)) then
                allocate(temp_substrings(size(substrings) + 1))
                temp_substrings(1:size(substrings)) = substrings
                substrings = temp_substrings
                deallocate(temp_substrings)
            else
                allocate(character(max_length) :: substrings(1))
            end if
            substrings(size(substrings)) = temp_string
            start_pos = space_pos + 1
        end do

    end subroutine split_string_by_space

    subroutine level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                      obj_func, hess_x_funptr, settings, trust_radius, &
                                      solution, solution_norm, mu, imicro, &
                                      imicro_jacobi_davidson, jacobi_davidson_started, &
                                      max_precision_reached, error)
        !
        ! this subroutine performs level-shifted (Jacobi-)Davidson to solve the trust 
        ! region subproblem
        !
        real(rp), intent(in) :: func, grad(:), grad_norm, h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), pointer, intent(in) :: obj_func
        procedure(hess_x_type), pointer, intent(in) :: hess_x_funptr
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: trust_radius
        real(rp), intent(out) :: solution(:), solution_norm, mu
        integer(ip), intent(out) :: imicro, imicro_jacobi_davidson, error
        logical, intent(out) :: jacobi_davidson_started, max_precision_reached

        real(rp), allocatable :: red_space_basis(:, :), h_basis(:, :), aug_hess(:, :), &
                                 red_space_solution(:), basis_vec(:), h_basis_vec(:), &
                                 h_solution(:), residual(:), solution_normalized(:), &
                                 last_solution_normalized(:), row_vec(:), col_vec(:)
        integer(ip) :: n_trial, i, initial_imicro                          
        logical :: accept_step, micro_converged, newton, bracketed
        real(rp) :: aug_hess_min_eigval, residual_norm, red_factor, &
                    initial_residual_norm, new_func, ratio, minres_tol
        real(rp), parameter :: newton_eigval_thresh = -1e-5_rp, &
                               solution_overlap_thresh = 0.5_rp, &
                               residual_norm_max_red_factor = 0.8_rp
        real(rp), external :: dnrm2, ddot
        external :: dgemm, dgemv

        ! initialize error flag
        error = 0

        ! generate trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error /= 0) return

        ! number of trial vectors
        n_trial = size(red_space_basis, 2)

        ! increment number of Hessian linear transformations
        tot_hess_x = tot_hess_x + n_trial

        ! calculate linear transformations of basis vectors
        allocate(h_basis(n_param, n_trial))
        do i = 1, n_trial
            call hess_x_funptr(red_space_basis(:, i), h_basis(:, i), error)
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) return
        end do

        ! construct augmented Hessian in reduced space
        allocate(aug_hess(n_trial + 1, n_trial + 1))
        aug_hess = 0.0_rp
        call dgemm("T", "N", n_trial, n_trial, n_param, 1.0_rp, red_space_basis, &
                   n_param, h_basis, n_param, 0.0_rp, aug_hess(2, 2), n_trial + 1)

        ! allocate space for reduced space solution and Hessian linear transformation
        ! of basis vector
        allocate(red_space_solution(n_trial), h_solution(n_param), basis_vec(n_param), &
                 h_basis_vec(n_param), last_solution_normalized(n_param))

        ! decrease trust radius until micro iterations converge and step is accepted
        last_solution_normalized = 0.0_rp
        accept_step = .false.
        do while (.not. accept_step)
            micro_converged = .false.

            jacobi_davidson_started = .false.
            do imicro = 1, settings%n_micro
                ! do a Newton step if the model is positive definite and the step is 
                ! within the trust region
                newton = .false.
                aug_hess_min_eigval = mat_min_eigval(aug_hess(2:, 2:), &
                                                     settings%hess_symm, settings, &
                                                     error)
                if (error /= 0) return
                if (aug_hess_min_eigval > newton_eigval_thresh) then
                    call newton_step(aug_hess, grad_norm, red_space_basis, &
                                     solution, red_space_solution, settings, error)
                    if (error /= 0) return
                    mu = 0.0_rp
                    if (dnrm2(n_param, solution, 1_ip) < trust_radius) newton = .true.
                end if

                ! otherwise perform bisection to find the level shift
                if (.not. newton) then
                    call bisection(aug_hess, grad_norm, red_space_basis, trust_radius, &
                                   solution, red_space_solution, mu, bracketed, &
                                   settings, error)
                    if (error /= 0) return
                    if (.not. bracketed) exit
                end if

                ! calculate Hessian linear transformation of solution
                call dgemv("N", n_param, n_trial, 1.0_rp, h_basis, n_param, &
                           red_space_solution, 1_ip, 0.0_rp, h_solution, 1_ip)

                ! calculate residual
                residual = grad + h_solution - mu*solution

                ! calculate residual norm
                residual_norm = dnrm2(n_param, residual, 1_ip)

                ! determine reduction factor depending on whether local region is
                ! reached
                if (abs(mu) < level_shift_local_thres) then
                    red_factor = settings%local_red_factor
                else
                    red_factor = settings%global_red_factor
                end if

                ! get normalized solution vector
                solution_normalized = solution / dnrm2(n_param, solution, 1_ip)

                ! reset initial residual norm if solution changes
                if (ddot(n_param, last_solution_normalized, 1_ip, solution_normalized, &
                         1_ip)**2 < solution_overlap_thresh) then
                    initial_imicro = imicro
                    initial_residual_norm = residual_norm
                end if

                ! check if micro iterations have converged
                if (residual_norm < max(red_factor * grad_norm, residual_norm_floor)) &
                    then
                    micro_converged = .true.
                    exit
                ! check if Jacobi-Davidson is used and has not been started
                else if (settings%subsystem_solver == "jacobi-davidson" .and. .not. &
                         jacobi_davidson_started) then
                    ! check residual has not decreased sufficiently or if maximum of 
                    ! Davidson iterations has been reached
                    if ((imicro - initial_imicro >= 10 .and. residual_norm > &
                         residual_norm_max_red_factor * initial_residual_norm) .or. &
                        imicro > settings%jacobi_davidson_start) then

                        ! switch to Jacobi-Davidson
                        jacobi_davidson_started = .true.
                        imicro_jacobi_davidson = imicro
                    end if
                ! check if residual has not decreased sufficiently
                else if (imicro - initial_imicro >= 10 .and. residual_norm > &
                         residual_norm_max_red_factor * initial_residual_norm) then
                    exit
                end if

                ! save current solution
                last_solution_normalized = solution_normalized

                if (.not. jacobi_davidson_started) then
                    ! precondition residual
                    call level_shifted_diag_precond(residual, mu, h_diag, basis_vec, &
                                                    settings, error)
                    if (error /= 0) return

                    ! orthonormalize to current orbital space to get new basis vector
                    call gram_schmidt(basis_vec, red_space_basis, settings, error)
                    if (error /= 0) return

                    ! add linear transformation of new basis vector
                    call hess_x_funptr(basis_vec, h_basis_vec, error)
                    call add_error_origin(error, error_hess_x, settings)
                    if (error /= 0) return

                    ! increment Hessian linear transformations
                    tot_hess_x = tot_hess_x + 1

                else
                    ! solve Jacobi-Davidson correction equations
                    minres_tol = 3.0_rp ** (-(imicro - imicro_jacobi_davidson))
                    call minres(-residual, hess_x_funptr, solution_normalized, mu, &
                                minres_tol, basis_vec, h_basis_vec, settings, error)
                    if (error /= 0) return

                    ! orthonormalize to current orbital space to get new basis vector
                    call gram_schmidt(basis_vec, red_space_basis, settings, error, &
                                      lin_trans_vector=h_basis_vec, &
                                      lin_trans_space=h_basis)
                    if (error /= 0) return

                    ! check if resulting linear transformation still respects Hessian 
                    ! symmetry which can happen due to numerical noise accumulation
                    if (abs(ddot(n_param, red_space_basis(:, n_trial), 1_ip, &
                                 h_basis_vec, 1_ip) - &
                            ddot(n_param, basis_vec, 1_ip, h_basis(:, n_trial), 1_ip)) &
                        > hess_symm_thres) then
                        call hess_x_funptr(basis_vec, h_basis_vec, error)
                        call add_error_origin(error, error_hess_x, settings)
                        if (error /= 0) return
                    end if

                end if

                ! increment trial vector count
                n_trial = n_trial + 1

                ! add new trial vector to orbital space
                call add_column(red_space_basis, basis_vec)

                ! add linear transformation of new basis vector
                call add_column(h_basis, h_basis_vec)

                ! construct new augmented Hessian
                allocate(row_vec(n_trial + 1))
                row_vec(1) = 0.0_rp
                call dgemv("T", n_param, n_trial , 1.0_rp, red_space_basis, n_param, &
                           h_basis(:, n_trial ), 1_ip, 0.0_rp, row_vec(2:), 1_ip)
                col_vec = row_vec
                if (.not. settings%hess_symm) then
                    call dgemv("T", n_param, n_trial, 1.0_rp, h_basis, n_param, &
                               red_space_basis(:, n_trial ), 1_ip, 0.0_rp, &
                               col_vec(2:), 1_ip)
                end if
                call extend_matrix(aug_hess, row_vec, col_vec)
                deallocate(row_vec)
                deallocate(col_vec)

                ! reallocate reduced space solution
                deallocate(red_space_solution)
                allocate(red_space_solution(n_trial))
            end do

            ! evaluate function at predicted point
            new_func = obj_func(solution, error)
            call add_error_origin(error, error_obj_func, settings)
            if (error /= 0) return

            ! calculate ratio of evaluated function and predicted function
            ratio = (new_func - func) / ddot(n_param, solution, 1_ip, &
                                             grad + 0.5_rp * h_solution, 1_ip)

            ! decide whether to accept step and modify trust radius
            accept_step = accept_trust_region_step(solution, ratio, micro_converged, &
                                                   settings, trust_radius, &
                                                   max_precision_reached)
            if (max_precision_reached) exit
        end do

        ! deallocate quantities from microiterations
        deallocate(red_space_solution, aug_hess, red_space_basis, h_basis, h_solution, &
                   residual, basis_vec, h_basis_vec, solution_normalized, &
                   last_solution_normalized)

        ! get norm of orbital rotation
        solution_norm = dnrm2(n_param, solution, 1_ip)

    end subroutine level_shifted_davidson

    subroutine truncated_conjugate_gradient(func, grad, grad_norm, h_diag, n_param, &
                                            obj_func, hess_x_funptr, settings, &
                                            trust_radius, solution, solution_norm, &
                                            n_micro, max_precision_reached, error)
        !
        ! this subroutine performs truncated conjugate gradient to solve the trust 
        ! region subproblem, this implementation is a bit different from standard TCG 
        ! since it is not only checking whether the current direction has negative 
        ! curvature but whether the entire subspace does by performing an Cholesky 
        ! factorization of the tridiagonal Lanczos matrix on the fly, this 
        ! implementation is based on the implementation of the Steihaug-Toint method 
        ! in the GALAHAD library (https://github.com/ralna/GALAHAD)
        !
        real(rp), intent(in) :: func, grad(:), grad_norm, h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), intent(in), pointer :: obj_func
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: trust_radius
        real(rp), intent(out) :: solution(:), solution_norm
        integer(ip), intent(out) :: n_micro, error
        logical, intent(out) :: max_precision_reached
        
        real(rp), allocatable :: residual(:), vector(:), basis_vec(:)
        real(rp) :: ratio, new_func, pred_func, conv_tol, step_size, &
                    trial_solution_dot, basis_vec_dot, solution_dot, &
                    solution_basis_vec_dot, residual_dot, residual_dot_old, beta, &
                    lanczos_diag_elem, lanczos_off_diag_elem, curvature, &
                    lanczos_tridiag_chol_pivot
        integer(ip) :: imicro
        logical :: accept_step, micro_converged
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        ! initialize number of microiterations
        n_micro = 0

        ! compute the stopping tolerance
        conv_tol = max(settings%local_red_factor * grad_norm, residual_norm_floor)

        ! allocate space for vectors
        allocate(residual(n_param), vector(n_param), basis_vec(n_param))

        ! iterate until step is accepted
        accept_step = .false.
        do while (.not. accept_step)
            ! reset microiteration convergence threshold
            micro_converged = .false.

            ! reset problem
            residual = grad
            call perturb_vector(residual)
            solution = 0.0_rp

            ! initialize micro iteration convergence flag
            micro_converged = .false.

            ! initialize error flag
            error = 0

            ! initialize iteration counter
            imicro = 0

            ! assume solution reaches trust radius
            solution_norm = trust_radius

            ! initialize predicted function
            pred_func = func

            ! initialize Lanczos diagonal element
            lanczos_diag_elem = 0.0_rp

            ! initialize dot products
            solution_dot = 0.0_rp
            solution_basis_vec_dot = 0.0_rp

            ! start of microiterations
            do
                ! obtain the preconditioned residual
                call abs_diag_precond(residual, h_diag, vector, settings, error)
                call add_error_origin(error, error_precond, settings)
                if (error /= 0) exit

                ! obtain the preconditioned residual dot product
                residual_dot = get_preconditioned_residual_dot(residual, vector, &
                                                               settings, error)
                if (error /= 0) exit

                ! get coupling coefficient and Lanczos tridiagonal elements
                if (imicro > 0) then
                    beta = residual_dot / residual_dot_old
                    lanczos_diag_elem = beta / step_size
                    lanczos_off_diag_elem = sqrt(beta) / abs(step_size)
                end if

                ! test for an approximate solution
                if (sqrt(residual_dot) <= conv_tol) then
                    solution_norm = sqrt(solution_dot)
                    micro_converged = .true.
                    exit
                end if

                ! obtain the search direction and dot products
                if (imicro > 0) then
                    ! test to see if iteration limit has been exceeded
                    if (imicro >= settings%n_micro) then
                        solution_norm = sqrt(solution_dot)
                        exit
                    end if

                    basis_vec = -vector + beta * basis_vec
                    solution_basis_vec_dot = beta * &
                                            (solution_basis_vec_dot + step_size * &
                                             basis_vec_dot)
                    basis_vec_dot = residual_dot + basis_vec_dot * beta * beta
                else
                    basis_vec = -vector
                    basis_vec_dot = residual_dot
                end if
                residual_dot_old = residual_dot

                ! test for convergence
                if (dnrm2(n_param, basis_vec, 1_ip) <= 2.0_rp * epsilon(1.0_rp)) then
                    solution_norm = sqrt(solution_dot)
                    micro_converged = .true.
                    exit
                end if

                ! increment number of microiterations
                imicro = imicro + 1

                ! obtain the Hessian linear transformation of the new basis vector
                call hess_x_funptr(basis_vec, vector, error)
                call add_error_origin(error, error_hess_x, settings)
                if (error /= 0) exit
                tot_hess_x = tot_hess_x + 1

                ! obtain the curvature
                curvature = ddot(n_param, vector, 1_ip, basis_vec, 1_ip)

                ! obtain the stepsize and the new diagonal of the Lanczos tridiagonal
                if (abs(curvature) > 0.0_rp) then
                    step_size = residual_dot / curvature
                    lanczos_diag_elem = lanczos_diag_elem + 1.0_rp / step_size
                ! no curvature present so take an infinite step
                else
                    step_size = huge(1.0_rp) ** 0.25
                end if

                ! check that the Lanczos tridiagonal is still positive definite
                if (imicro > 1) then
                    lanczos_tridiag_chol_pivot = &
                        lanczos_diag_elem - &
                        (lanczos_off_diag_elem / lanczos_tridiag_chol_pivot) * &
                        lanczos_off_diag_elem
                else
                    lanczos_tridiag_chol_pivot = lanczos_diag_elem
                end if

                ! the matrix is indefinite
                if (lanczos_tridiag_chol_pivot <= 0.0_rp) then
                    ! find the appropriate point on the boundary
                    call find_point_on_boundary(basis_vec, trust_radius, &
                                                solution_basis_vec_dot, basis_vec_dot, &
                                                residual_dot, curvature, solution_dot, &
                                                solution, pred_func, step_size)
                    micro_converged = .true.
                    exit
                end if

                ! see if the new point is also interior
                trial_solution_dot = solution_dot + step_size * &
                                     (solution_basis_vec_dot + solution_basis_vec_dot &
                                      + step_size * basis_vec_dot)

                ! the new point is interior
                if (trial_solution_dot <= trust_radius ** 2) then
                    solution = solution + step_size * basis_vec
                    solution_dot = trial_solution_dot
                    pred_func = pred_func - 0.5_rp * step_size * step_size * curvature
                ! the new point is outside the trust region
                else
                    ! find the appropriate point on the boundary
                    call find_point_on_boundary(basis_vec, trust_radius, &
                                                solution_basis_vec_dot, basis_vec_dot, &
                                                residual_dot, curvature, solution_dot, &
                                                solution, pred_func, step_size)
                    micro_converged = .true.
                    exit
                end if

                ! update the residual
                residual = residual + step_size * vector
            end do
            if (error /= 0) exit

            ! add number of microiterations to total number of microiterations
            n_micro = n_micro + imicro

            ! evaluate function at predicted point
            new_func = obj_func(solution, error)
            call add_error_origin(error, error_obj_func, settings)
            if (error /= 0) exit

            ! calculate ratio of evaluated function and predicted function
            ratio = (new_func - func) / (pred_func - func)

            ! decide whether to accept step and modify trust radius
            accept_step = accept_trust_region_step(solution, ratio, micro_converged, &
                                                   settings, trust_radius, &
                                                   max_precision_reached)
            if (max_precision_reached) exit

        end do

        ! deallocate vectors
        deallocate(residual, vector, basis_vec)

    end subroutine truncated_conjugate_gradient

    subroutine generalized_lanczos_trust_region(func, grad, grad_norm, h_diag, &
                                                n_param, obj_func, hess_x_funptr, &
                                                settings, trust_radius, solution, &
                                                solution_norm, lambda, n_micro, &
                                                max_precision_reached, error)
        !
        ! this subroutine performs generalized lanczos trust region to solve the trust 
        ! region subproblem, this implementation is based on the implementation of GLTR 
        ! in the GALAHAD library (https://github.com/ralna/GALAHAD)
        !
        real(rp), intent(in) :: func, grad(:), grad_norm, h_diag(:)
        integer(ip), intent(in) :: n_param
        procedure(obj_func_type), intent(in), pointer :: obj_func
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: trust_radius
        real(rp), intent(out) :: solution(:), solution_norm, lambda
        integer(ip), intent(out) :: n_micro, error
        logical, intent(out) :: max_precision_reached
        
        real(rp), allocatable :: residual(:), eigenvec(:), lanczos_diag(:), &
                                 lanczos_off_diag(:), lanczos_diag_fact(:), &
                                 lanczos_off_diag_fact(:), red_space_rhs(:), &
                                 red_space_solution(:), red_space_eigenvec(:), &
                                 work(:), stepsize_list(:), residual_dot_list(:)
        real(rp) :: ratio, new_func, func_diff, lowest_eigval, hard_case_step_size, &
                    tau, pred_func, red_factor, conv_tol
        integer(ip) :: n_first_pass, n_second_pass, n_saved, n_red_space
        logical :: restart_lanczos, accept_step, micro_converged, hard_case, interior
        real(rp), external :: ddot, dnrm2

        ! initialize error flag
        error = 0

        ! initialize number of microiterations
        n_micro = 0
        
        ! set restart Lanczos boolean
        restart_lanczos = .false.

        ! allocate space for residual and for lowest eigenvector
        allocate(residual(n_param), eigenvec(n_param))

        ! allocate space for Lanczos tridiagonal
        allocate(lanczos_diag(settings%n_micro + 1), &
                 lanczos_off_diag(settings%n_micro), &
                 lanczos_diag_fact(settings%n_micro + 1), &
                 lanczos_off_diag_fact(settings%n_micro))

        ! allocate space for reduced space RHS, solution and eigenvector
        allocate(red_space_rhs(settings%n_micro + 1), &
                 red_space_solution(settings%n_micro + 1), &
                 red_space_eigenvec(settings%n_micro + 1))

        ! allocate work array for solving Lanczos subproblem
        allocate(work(settings%n_micro + 1))

        ! allocate space to store the stepsizes and residual norms which allows for 
        ! more efficient processing in the second pass
        allocate(stepsize_list(settings%n_micro), &
                 residual_dot_list(settings%n_micro))

        ! iterate until step is accepted
        accept_step = .false.
        do while (.not. accept_step)
            ! reset microiteration convergence threshold
            micro_converged = .false.

            ! reset problem
            residual = grad
            call perturb_vector(residual)
            solution = 0.0_rp
            eigenvec = 0.0_rp

            ! GLTR minimizer block
            gltr_minimizer: block
                ! check whether Lanczos is being run for the first time
                if (.not. restart_lanczos .or. n_red_space <= 0) then
                    ! perform first pass
                    call gltr_first_pass(func, grad_norm, h_diag, hess_x_funptr, &
                                         trust_radius, residual, solution, eigenvec, &
                                         lanczos_diag, lanczos_off_diag, &
                                         lanczos_diag_fact, lanczos_off_diag_fact, &
                                         red_space_rhs, red_space_solution, &
                                         red_space_eigenvec, work, stepsize_list, &
                                         residual_dot_list, pred_func, lambda, &
                                         solution_norm, lowest_eigval, tau, &
                                         micro_converged, interior, hard_case, &
                                         hard_case_step_size, n_first_pass, &
                                         n_red_space, n_saved, settings, error)
                    if (error /= 0) exit gltr_minimizer

                    ! check if number of micro iterations has exceeded limit or if 
                    ! interior solution has been found
                    if ((n_first_pass >= settings%n_micro .and. interior) .or. &
                         micro_converged) then
                        n_second_pass = 0
                        exit gltr_minimizer
                    end if

                ! repeated Lanczos solution with smaller trust-region radius
                else
                    ! no first pass necessary as we can reuse Lanczos factorization 
                    ! from first run for new trust radius
                    n_first_pass = 0

                    ! no vectors saved for this trust radius so full second pass 
                    ! necessary
                    n_saved = 0

                    ! find the solution to the Lanczos TR subproblem with this radius
                    call solve_tridiagonal_subproblem( &
                        n_red_space, lanczos_diag(:n_red_space), &
                        lanczos_off_diag(:n_red_space - 1), &
                        lanczos_diag_fact(:n_red_space), &
                        lanczos_off_diag_fact(:n_red_space - 1), &
                        red_space_rhs(:n_red_space), trust_radius, interior, .true., &
                        .false., lowest_eigval, lambda, func_diff, &
                        red_space_solution(:n_red_space), &
                        red_space_eigenvec(:n_red_space), work(:n_red_space), &
                        hard_case, hard_case_step_size)

                    ! record the optimal objective function value
                    pred_func = func + func_diff
                    
                    ! determine reduction factor depending on whether local region is 
                    ! reached
                    if (abs(lambda) < level_shift_local_thres) then
                        red_factor = settings%local_red_factor
                    else
                        red_factor = settings%global_red_factor
                    end if

                    ! compute the stopping tolerance
                    conv_tol = max(red_factor * grad_norm, residual_norm_floor)

                    ! check whether solution satisfies convergence criteria or whether 
                    ! trust radius needs to be decreased further
                    if (abs(lanczos_off_diag(n_red_space - 1) * &
                            red_space_solution(n_red_space)) > conv_tol) &
                        exit gltr_minimizer

                    ! intialize tau for second pass
                    tau = 1.0_rp

                    ! solution reaches trust radius
                    solution_norm = trust_radius

                end if

                ! second pass to obtain solution
                call gltr_second_pass(residual, tau, red_space_solution, &
                                      red_space_eigenvec, residual_dot_list, &
                                      stepsize_list, h_diag, hess_x_funptr, &
                                      n_red_space, n_saved, lambda, hard_case, &
                                      hard_case_step_size, solution, eigenvec, &
                                      n_second_pass, settings, error)
                if (error /= 0) exit gltr_minimizer

                ! check if number of micro iterations has exceeded limit
                if (n_first_pass >= settings%n_micro) exit gltr_minimizer

                ! successful return
                micro_converged = .true.

            end block gltr_minimizer
            if (error /= 0) exit

            ! add number of microiterations from first and second pass to total number 
            ! of microiterations
            n_micro = n_micro + n_first_pass

            ! evaluate function at predicted point
            new_func = obj_func(solution, error)
            call add_error_origin(error, error_obj_func, settings)
            if (error /= 0) exit

            ! calculate ratio of evaluated function and predicted function
            ratio = (new_func - func) / (pred_func - func)

            ! decide whether to accept step and modify trust radius
            accept_step = accept_trust_region_step(solution, ratio, micro_converged, &
                                                   settings, trust_radius, &
                                                   max_precision_reached)
            if (max_precision_reached) exit

            ! restart Lanczos with smaller trust region if step is not accepted
            restart_lanczos = .not. accept_step
            end do

        ! deallocate arrays
        deallocate(residual, eigenvec, lanczos_diag, lanczos_off_diag, &
                   lanczos_diag_fact, lanczos_off_diag_fact, red_space_rhs, &
                   red_space_solution, red_space_eigenvec, work, stepsize_list, &
                   residual_dot_list)

    end subroutine generalized_lanczos_trust_region

    logical function accept_trust_region_step(solution, ratio, micro_converged, &
                                              settings, trust_radius, &
                                              max_precision_reached)
        !
        ! this function checks whether the trust region step is accepted and modified 
        ! the trust region accordingly
        !
        real(rp), intent(in) :: solution(:), ratio
        logical, intent(in) :: micro_converged
        type(solver_settings_type), intent(in) :: settings
        real(rp), intent(inout) :: trust_radius
        logical, intent(out) :: max_precision_reached

        ! default to maximum precision not yet reached
        max_precision_reached = .false.

        ! decrease trust radius if micro iterations are unable to converge, if function 
        ! value has not decreased or if individual orbitals change too much
        if (.not. micro_converged .or. ratio < 0.0_rp .or. any(abs(solution) > pi/4)) &
            then
            trust_radius = trust_radius_shrink_factor * trust_radius
            accept_trust_region_step = .false.
            if (trust_radius < numerical_zero) then
                call settings%log("Trust radius too small. Convergence criterion "// &
                                  "is not fulfilled but calculation should be "// &
                                  "converged up to floating point precision.", &
                                  verbosity_error, .true.)
                max_precision_reached = .true.
                return
            end if
        ! check if step is too long
        else if (ratio < trust_radius_shrink_ratio) then
            trust_radius = trust_radius_shrink_factor * trust_radius
            accept_trust_region_step = .true.
        ! check if quadratic approximation is valid
        else if (ratio < trust_radius_expand_ratio) then
            accept_trust_region_step = .true.
        ! check if step is potentially too short
        else
            trust_radius = trust_radius_expand_factor * trust_radius
            accept_trust_region_step = .true.
        end if

    end function accept_trust_region_step

    subroutine solver_sanity_check(settings, n_param, grad, error)
        !
        ! this subroutine performs a sanity check for solver input parameters
        !
        type(solver_settings_type), intent(inout) :: settings
        integer(ip), intent(in) :: n_param
        real(rp), intent(in) :: grad(:)
        integer(ip), intent(out) :: error

        character(300) :: msg

        ! initialize error flag
        error = 0

        ! check that number of parameters is positive
        if (n_param < 1) then
            call settings%log("Number of parameters should be larger than 0.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! convert strings to lowercase
        settings%subsystem_solver = string_to_lowercase(settings%subsystem_solver)

        ! check that number of random trial vectors is below number of parameters
        if ((settings%subsystem_solver == "davidson" .or. &
             settings%subsystem_solver == "jacobi-davidson") .and. &
            settings%n_random_trial_vectors > n_param/2) then
            settings%n_random_trial_vectors = n_param/2
            write (msg, '(A, I0, A)') random_trial_vector_warning_msg//" Setting to ", &
                settings%n_random_trial_vectors, "."
            call settings%log(msg, verbosity_warning)
        end if

        ! sanity check for gradient size
        if (size(grad) /= n_param) then
            call settings%log("Size of gradient array returned by subroutine "// &
                              "update_orbs does not equal number of parameters.", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check for character options
        if (.not. (settings%subsystem_solver == "davidson" .or. &
                   settings%subsystem_solver == "jacobi-davidson" .or. &
                   settings%subsystem_solver == "tcg" .or. &
                   settings%subsystem_solver == "gltr")) then
            call settings%log("Subsystem solver option unknown. Possible values "// &
                              "are ""davidson"", ""jacobi-davidson"", ""tcg"" "// &
                              "(truncated conjugate gradient), and ""gltr"" "// &
                              "(generalized Lanczos trust region).", verbosity_error, &
                              .true.)
            error = 1
            return
        end if

        ! check whether projection functions is passed
        if (associated(settings%project)) call settings%log(project_warning_msg, &
                                                            verbosity_warning)

    end subroutine solver_sanity_check

    subroutine stability_sanity_check(settings, n_param, error)
        !
        ! this subroutine performs a sanity check for stability input parameters
        !
        type(stability_settings_type), intent(inout) :: settings
        integer(ip), intent(in) :: n_param
        integer(ip), intent(out) :: error
        character(300) :: msg

        ! initialize error flag
        error = 0

        ! convert strings to lowercase
        settings%diag_solver = string_to_lowercase(settings%diag_solver)

        ! check that number of random trial vectors is below number of parameters
        if (settings%n_random_trial_vectors > n_param/2) then
            settings%n_random_trial_vectors = n_param/2
            write (msg, '(A, I0, A)') random_trial_vector_warning_msg//" Setting to ", &
                settings%n_random_trial_vectors, "."
            call settings%log(msg, verbosity_warning)
        end if

        ! check for character options
        if (.not. (settings%diag_solver == "davidson" .or. &
                   settings%diag_solver == "jacobi-davidson")) then
            call settings%log("Diagonalization solver option unknown. Possible "// &
                              "values are ""davidson"" and ""jacobi-davidson""", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

        ! check whether projection functions is passed
        if (associated(settings%project)) call settings%log(project_warning_msg, &
                                                            verbosity_warning)

    end subroutine stability_sanity_check

    subroutine add_error_origin(error_code, error_origin, settings)
        !
        ! this function modifies the error code by adding the error's origin if it is 
        ! not already added
        !
        class(settings_type), intent(in) :: settings
        integer(ip), intent(inout) :: error_code
        integer(ip), intent(in) :: error_origin

        if (error_code > 0) then
            if (error_code < 100) error_code = error_origin + error_code
        else if (error_code < 0) then
            call settings%log("Negative error code encountered.", verbosity_error, &
                              .true.)
            error_code = error_origin + 1
        end if

    end subroutine add_error_origin

    function string_to_lowercase(str) result(lower_str)
        !
        ! this function converts a string to lower case
        !
        character(*), intent(in) :: str
        character(len(str)) :: lower_str
        integer(ip) :: i, code
        integer(ip), parameter :: ascii_upper_lower_diff = 32

        do i = 1, len(str)
            code = iachar(str(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower_str(i:i) = achar(code + ascii_upper_lower_diff)
            else
                lower_str(i:i) = str(i:i)
            end if
        end do
    
    end function string_to_lowercase

    subroutine gltr_first_pass(func, grad_norm, h_diag, hess_x_funptr, trust_radius, &
                               residual, solution, eigenvec, lanczos_diag, &
                               lanczos_off_diag, lanczos_diag_fact, &
                               lanczos_off_diag_fact, red_space_rhs, &
                               red_space_solution, red_space_eigenvec, work, &
                               stepsize_list, residual_dot_list, pred_func, lambda, &
                               solution_norm, lowest_eigval, tau, micro_converged, &
                               interior, hard_case, hard_case_step_size, imicro, &
                               n_red_space, n_saved, settings, error)
        !
        ! this subroutine performs the first Lanczos pass to compute the tridiagonal 
        ! matrix and the right hand side of the reduced problem
        !
        real(rp), intent(in) :: func, grad_norm, h_diag(:), trust_radius
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(rp), intent(inout) :: residual(:), solution(:), eigenvec(:)
        real(rp), intent(out) :: lanczos_diag(:), lanczos_off_diag(:), &
                                 lanczos_diag_fact(:), lanczos_off_diag_fact(:), &
                                 red_space_rhs(:), red_space_solution(:), &
                                 red_space_eigenvec(:), work(:), stepsize_list(:), &
                                 residual_dot_list(:), pred_func, lambda, &
                                 solution_norm, lowest_eigval, tau, hard_case_step_size
        logical, intent(out) :: micro_converged, interior, hard_case
        integer(ip), intent(out) :: imicro, n_red_space, n_saved
        type(solver_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        ! an estimate of the solution that gives at least fraction_opt times the 
        ! optimal objective value will be found, only saves computation time if second 
        ! pass is necessary
        real(rp), parameter :: fraction_opt = 1.0_rp

        real(rp), allocatable :: vector(:), basis_vec(:), min_func_list(:), &
                                 residual_start(:), residual_save(:), &
                                 basis_vec_save(:), precond_residuals_save(:, :)
        integer(ip) :: n_param, it, switch_iteration, extra_vectors, istat
        real(rp) :: step_size, f_tol, alpha, trial_solution_dot, u_norm, &
                    basis_vec_dot, solution_dot, solution_basis_vec_dot, residual_dot, &
                    residual_dot_old, beta, lanczos_diag_elem, lanczos_off_diag_elem, &
                    residual_norm, red_factor, conv_tol, last_red_space_solution, &
                    curvature, lanczos_tridiag_chol_pivot
        logical :: negative_curvature, try_warm, use_old
        real(rp), external :: dnrm2, ddot

        ! initialize error flag
        error = 0

        ! save starting residual for use in second pass
        residual_start = residual

        ! number of parameters
        n_param = size(residual)

        ! initialize micro iteration convergence flag
        micro_converged = .false.

        ! start interior to trust region
        interior = .true.

        ! no initial guess for Lagrange multiplier
        lambda = 0.0_rp

        ! initialize reduced space RHS
        red_space_rhs = 0.0_rp
        red_space_rhs(1) = 1.0_rp

        ! initialize coefficient for last basis vector
        last_red_space_solution = 0.0_rp

        ! initialize iteration counter
        imicro = 0

        ! initialize number of saved vectors
        n_saved = 0

        ! initialize hard case boolean
        hard_case = .false.

        ! initialize negative curvature boolean
        negative_curvature = .false.

        ! initialize warm start boolean which provides a potentially good guess for the 
        ! Lagrange multiplier
        try_warm = .false.

        ! initialize boolean which indicates that the lowest eigenvalue of the
        ! leading n-1 by n-1 block is given
        use_old = .false.

        ! assume solution reaches trust radius
        solution_norm = trust_radius

        ! initialize predicted function
        pred_func = func

        ! initialize Lanczos diagonal element
        lanczos_diag_elem = 0.0_rp

        ! initialize dot products
        solution_dot = 0.0_rp
        solution_basis_vec_dot = 0.0_rp

        ! allocate vector for the preconditioned residual and the Hessian vector 
        ! product and basis vector
        allocate(vector(n_param), basis_vec(n_param))

        ! allocate workspace for the sequence of smallest function values
        allocate(min_func_list(settings%n_micro + 1))

        ! check whether we can afford to store extra vectors to avoid recomputation in 
        ! the second pass
        extra_vectors = 0
        allocate(residual_save(n_param), basis_vec_save(n_param), stat=istat)
        if (istat == 0) then
            extra_vectors = settings%n_micro
            do
                allocate(precond_residuals_save(n_param, extra_vectors), &
                         stat=istat)
                if (istat == 0) exit
                extra_vectors = extra_vectors / 2
                if (extra_vectors == 0) then
                    deallocate(residual_save, basis_vec_save)
                    exit
                end if
            end do
        end if

        ! start of microiterations
        do
            ! obtain the preconditioned residual
            call abs_diag_precond(residual, h_diag, vector, settings, error)
            call add_error_origin(error, error_precond, settings)
            if (error /= 0) exit

            ! obtain the preconditioned residual dot product
            residual_dot = get_preconditioned_residual_dot(residual, vector, settings, &
                                                           error)
            if (error /= 0) exit
            residual_norm = sqrt(residual_dot)

            ! if the user has asked to save vectors, save preconditioned residual, 
            ! residual, and basis vector
            if (extra_vectors > 0) then
                if (imicro < extra_vectors) &
                    precond_residuals_save(:, imicro + 1) = vector
                if (imicro == extra_vectors) then
                    residual_save = residual
                    basis_vec_save = basis_vec
                end if
            end if

            ! get coupling coefficient and Lanczos tridiagonal elements
            if (imicro > 0) then
                beta = residual_dot / residual_dot_old
                lanczos_diag_elem = beta / step_size
                lanczos_off_diag_elem = sqrt(beta) / abs(step_size)
            ! set reduced space RHS for first iteration
            else
                red_space_rhs(1) = residual_norm
            end if

            ! determine reduction factor depending on whether local region is reached
            if (abs(lambda) < level_shift_local_thres) then
                red_factor = settings%local_red_factor
            else
                red_factor = settings%global_red_factor
            end if

            ! compute the stopping tolerance
            conv_tol = max(red_factor * grad_norm, residual_norm_floor)

            ! test for an interior approximate solution
            if (interior .and. residual_norm <= conv_tol) then
                solution_norm = sqrt(solution_dot)
                n_red_space = imicro
                if (n_red_space > 0) lambda = 0.0_rp
                micro_converged = .true.
                exit
            end if

            if (imicro > 0) then
                ! test to see if iteration limit has been exceeded
                if (imicro >= settings%n_micro .and. interior) then
                    solution_norm = sqrt(solution_dot)
                    n_red_space = imicro
                    exit
                end if

                ! obtain the search direction and dot products
                basis_vec = -vector + beta * basis_vec
                solution_basis_vec_dot = beta * &
                                         (solution_basis_vec_dot + step_size * &
                                          basis_vec_dot)
                basis_vec_dot = residual_dot + basis_vec_dot * beta * beta

                ! continue accumulating the Lanczos tridiagonal
                lanczos_diag(imicro + 1) = lanczos_diag_elem
                lanczos_off_diag(imicro) = lanczos_off_diag_elem

                ! check whether convergence on the trust region boundary has been 
                ! achieved or the iteration limit is reached
                if (imicro >= settings%n_micro .OR. (.NOT. interior .and. &
                    abs(lanczos_off_diag_elem * last_red_space_solution) <= conv_tol)) then
                    ! check whether any earlier point produces fraction_opt of the 
                    ! optimal solution and if yes, use this point to avoid iterations 
                    ! in second Lanczos pass
                    if (fraction_opt < 1.0_rp) then
                        f_tol = min_func_list(imicro) * fraction_opt
                        do n_red_space = 1, imicro
                            if (min_func_list(n_red_space) <= f_tol) exit
                        end do
                    else
                        n_red_space = imicro
                    end if

                    ! the required fraction of the optimal solution was achieved by an 
                    ! interior point, second pass is not needed
                    if (n_red_space <= switch_iteration) then
                        solution_norm = sqrt(solution_dot)
                        micro_converged = .true.
                        exit
                    end if

                    ! restore the solution to the Lanczos TR subproblem for this 
                    ! iteration
                    use_old = .false.
                    if (n_red_space <= 1 + switch_iteration) then
                        lambda = 0.0_rp
                        try_warm = .false.
                    end if

                    call solve_tridiagonal_subproblem( &
                        n_red_space, lanczos_diag(:n_red_space), &
                        lanczos_off_diag(:n_red_space - 1), &
                        lanczos_diag_fact(:n_red_space), &
                        lanczos_off_diag_fact(:n_red_space - 1), &
                        red_space_rhs(:n_red_space), trust_radius, interior, try_warm, &
                        use_old, lowest_eigval, lambda, min_func_list(n_red_space), &
                        red_space_solution(:n_red_space), &
                        red_space_eigenvec(:n_red_space), work(:n_red_space), &
                        hard_case, hard_case_step_size)

                    ! record the optimal objective function value and prepare to 
                    ! recover the approximate solution
                    pred_func = func + min_func_list(n_red_space)
                    tau = 1.0_rp

                    ! use saved vectors to start second pass
                    if (extra_vectors > 0) then
                        n_saved = min(n_red_space, extra_vectors)
                        do it = 1, n_saved
                            red_space_solution(it) = tau * &
                                                     (red_space_solution(it) / &
                                                      sqrt(residual_dot_list(it)))
                            if (hard_case) &
                                red_space_eigenvec(it) = tau * &
                                                         (red_space_eigenvec(it) / &
                                                          sqrt(residual_dot_list(it)))
                            tau = -sign(1.0_rp, stepsize_list(it)) * tau
                        end do

                        ! update the solution estimate using the saved vectors
                        solution = MATMUL(precond_residuals_save(:, :n_saved), &
                                          red_space_solution(1:n_saved))
                        if (hard_case) &
                            eigenvec = MATMUL(precond_residuals_save(:, :n_saved), &
                                              red_space_eigenvec(1:n_saved))

                        ! second pass not needed because number of saved vectors is 
                        ! sufficient to recover the solution
                        if (n_saved == n_red_space) then
                            ! if the hard case has occured, ensure that the recovered 
                            ! eigenvector has unit norm and compute the complete 
                            ! solution
                            if (hard_case) then
                                call hess_x_funptr(eigenvec, vector, error)
                                call add_error_origin(error, error_hess_x, settings)
                                if (error /= 0) exit
                                tot_hess_x = tot_hess_x + 1
                                u_norm = sqrt(-ddot(n_param, vector, 1_ip, eigenvec, &
                                                    1_ip) / lambda)
                                solution = solution + (hard_case_step_size / u_norm) * &
                                           eigenvec
                            end if
                            call abs_diag_precond(solution, h_diag, vector, settings, &
                                                  error)
                            call add_error_origin(error, error_precond, settings)
                            if (error /= 0) return
                            micro_converged = .true.
                            exit
                        end if
                        residual = residual_save
                        basis_vec = basis_vec_save
                        residual_dot_old = residual_dot_list(n_saved)
                    ! start second pass without any saved vectors, need to recompute 
                    ! the Lanczos factorization
                    else
                        residual = residual_start
                    end if
                    ! solution on boundary found
                    exit
                end if
            else
                ! obtain the search direction and dot products
                basis_vec = -vector
                basis_vec_dot = residual_dot
                lanczos_diag(1) = lanczos_diag_elem
            end if
            residual_dot_list(imicro + 1) = residual_dot
            residual_dot_old = residual_dot

            ! test for convergence
            if (interior .and. &
                dnrm2(n_param, basis_vec, 1_ip) <= 2.0_rp * epsilon(1.0_rp)) then
                solution_norm = sqrt(solution_dot)
                n_red_space = imicro
                micro_converged = .true.
                exit
            end if

            ! increment number of microiterations
            imicro = imicro + 1

            ! obtain the Hessian linear transformation of the new basis vector
            call hess_x_funptr(basis_vec, vector, error)
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) exit
            tot_hess_x = tot_hess_x + 1

            ! obtain the curvature
            curvature = ddot(n_param, vector, 1_ip, basis_vec, 1_ip)

            ! obtain the stepsize and the new diagonal of the Lanczos tridiagonal
            if (abs(curvature) > 0.0_rp) then
                step_size = residual_dot / curvature
                lanczos_diag_elem = lanczos_diag_elem + 1.0_rp / step_size
            ! no curvature present so take an infinite step
            else
                step_size = huge(1.0_rp) ** 0.25
            end if

            ! check that the Lanczos tridiagonal is still positive definite
            if (.NOT. negative_curvature) then
                if (imicro > 1) then
                    lanczos_tridiag_chol_pivot = &
                        lanczos_diag_elem - &
                        (lanczos_off_diag_elem / lanczos_tridiag_chol_pivot) * &
                        lanczos_off_diag_elem
                else
                    lanczos_tridiag_chol_pivot = lanczos_diag_elem
                end if
                negative_curvature = lanczos_tridiag_chol_pivot <= 0.0_rp
            end if

            ! the matrix is indefinite
            if (interior .and. negative_curvature) then
                ! find the appropriate point on the boundary
                call find_point_on_boundary(basis_vec, trust_radius, &
                                            solution_basis_vec_dot, basis_vec_dot, &
                                            residual_dot, curvature, solution_dot, &
                                            solution, pred_func, alpha)

                ! when the model has no curvature in the new basis vector direction 
                ! find the appropriate point and the gradient (residual) on the 
                ! boundary and stop
                if (abs(curvature) < numerical_zero) then
                    step_size = alpha
                    n_red_space = imicro
                    lambda = 0.0_rp
                    residual = residual + step_size * vector
                    micro_converged = .true.
                    exit
                ! if a more accurate solution is required, switch modes
                else
                    interior = .false.
                    switch_iteration = imicro - 1
                end if
            end if

            ! if the current estimate of the solution is interior, see if the new point
            ! is also interior
            if (interior) then
                trial_solution_dot = solution_dot + step_size * &
                                     (solution_basis_vec_dot + solution_basis_vec_dot &
                                      + step_size * basis_vec_dot)

                ! the new point is interior
                if (trial_solution_dot <= trust_radius ** 2) then
                    solution = solution + step_size * basis_vec
                    solution_dot = trial_solution_dot
                    pred_func = pred_func - 0.5_rp * step_size * step_size * curvature
                    min_func_list(imicro) = pred_func - func
                ! the new point is outside the trust region
                else
                    ! find the appropriate point on the boundary
                    call find_point_on_boundary(basis_vec, trust_radius, &
                                                solution_basis_vec_dot, basis_vec_dot, &
                                                residual_dot, curvature, solution_dot, &
                                                solution, pred_func, alpha)

                    ! switch modes
                    interior = .false.
                    switch_iteration = imicro - 1
                end if
            end if

            ! complete the new diagonal of the Lanczos tridiagonal matrix
            lanczos_diag(imicro) = lanczos_diag_elem
            stepsize_list(imicro) = step_size

            ! solve the subproblem if new point is not interior
            if (.NOT. interior) then
                call solve_tridiagonal_subproblem(imicro, lanczos_diag(:imicro), &
                                                  lanczos_off_diag(:imicro - 1), &
                                                  lanczos_diag_fact(:imicro), &
                                                  lanczos_off_diag_fact(:imicro - 1), &
                                                  red_space_rhs(:imicro), &
                                                  trust_radius, interior, try_warm, &
                                                  use_old, lowest_eigval, lambda, &
                                                  min_func_list(imicro), &
                                                  red_space_solution(:imicro), &
                                                  red_space_eigenvec(:imicro), &
                                                  work(:imicro), hard_case, &
                                                  hard_case_step_size)

                ! do a warm start next since we have a guess for the Lagrange multiplier
                try_warm = .true.

                use_old = lowest_eigval < 0.0_rp
                last_red_space_solution = red_space_solution(imicro)

            end if

            ! update the residual
            residual = residual + step_size * vector
        end do

        ! deallocate vectors
        deallocate(residual_start, vector, basis_vec, min_func_list)
        if (extra_vectors > 0) &
            deallocate(residual_save, basis_vec_save, &
                       precond_residuals_save)

    end subroutine gltr_first_pass
      
    subroutine gltr_second_pass(residual_start, tau_start, red_space_solution, &
                                red_space_eigenvec, residual_dot_list, stepsize_list, &
                                h_diag, hess_x_funptr, n_red_space, istart, lambda, &
                                hard_case, hard_case_step_size, solution, eigenvec, &
                                n_second_pass, settings, error)
        !
        ! this subroutine performs the second Lanczos pass to compute the approximate 
        ! solution and eigenvector if the hard case has occured, the starting point are 
        ! these vectors after the boundary point has been reached or after saved 
        ! vectors have been added if requested
        !
        real(rp), intent(in) :: residual_start(:), tau_start, red_space_solution(:), &
                                red_space_eigenvec(:), residual_dot_list(:), &
                                stepsize_list(:), h_diag(:)
        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        integer(ip), intent(in) :: n_red_space
        integer(ip), intent(in) :: istart
        real(rp), intent(in) :: lambda, hard_case_step_size
        logical, intent(in) :: hard_case
        real(rp), intent(inout) :: solution(:), eigenvec(:)
        integer(ip), intent(out) :: n_second_pass
        type(solver_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        
        integer(ip) :: n_param, imicro
        real(rp), allocatable :: vector(:), basis_vec(:), residual(:)
        real(rp) :: tau, step_size, beta, residual_dot, residual_dot_old, u_norm
        real(rp), external :: ddot

        ! initialize error flag
        error = 0

        ! initialize residual
        residual = residual_start

        ! number of parameters
        n_param = size(residual_start)

        ! initialize tau in case of saved vectors
        tau = tau_start

        ! initialize iteration counter in case of saved vectors
        imicro = istart

        ! allocate vector for the preconditioned residual and the Hessian vector product
        allocate(vector(n_param), basis_vec(n_param))

        ! start second pass loop
        do
            ! obtain the preconditioned residual
            call abs_diag_precond(residual, h_diag, vector, settings, error)
            call add_error_origin(error, error_precond, settings)
            if (error /= 0) exit

            ! obtain the scaled norm of the residual
            residual_dot = residual_dot_list(imicro + 1)

            ! update the solution estimate
            if (imicro /= 0) then
                solution = solution + tau * &
                           (red_space_solution(imicro + 1) / sqrt(residual_dot)) * &
                           vector
                if (hard_case) &
                    eigenvec = eigenvec + tau * &
                               (red_space_eigenvec(imicro + 1) / sqrt(residual_dot)) * &
                               vector
            else
                solution = tau * (red_space_solution(imicro + 1) / sqrt(residual_dot)) &
                           * vector
                if (hard_case) &
                    eigenvec = tau * &
                               (red_space_eigenvec(imicro + 1) / sqrt(residual_dot)) * &
                               vector
            end if

            ! if the approximate minimizer is complete, exit
            if (imicro + 1 == n_red_space) then
                n_second_pass = imicro - istart

                ! if the hard case has occured, ensure that the recovered eigenvector 
                ! has unit norm and compute the complete solution
                if (hard_case) then
                    call hess_x_funptr(eigenvec, vector, error)
                    call add_error_origin(error, error_hess_x, settings)
                    if (error /= 0) exit
                    tot_hess_x = tot_hess_x + 1
                    u_norm = sqrt(-ddot(n_param, vector, 1_ip, eigenvec, 1_ip) / lambda)
                    solution = solution + (hard_case_step_size / u_norm) * eigenvec
                end if

                exit
            end if

            ! get new basis vector
            if (imicro > 0) then
                beta = residual_dot / residual_dot_old
                basis_vec = - vector + beta * basis_vec
            else
                basis_vec = - vector
            end if
            residual_dot_old = residual_dot

            ! incremenet interation
            imicro = imicro + 1

            ! obtain the linear transformation of new basis vector
            call hess_x_funptr(basis_vec, vector, error)
            call add_error_origin(error, error_hess_x, settings)
            if (error /= 0) exit
            tot_hess_x = tot_hess_x + 1

            ! retreive the stepsize
            step_size = stepsize_list(imicro)

            ! update the residual
            residual = residual + step_size * vector
            tau = -sign(1.0_rp, step_size) * tau
        end do

        ! deallocate vectors
        deallocate(residual, vector, basis_vec)

    end subroutine gltr_second_pass

    subroutine perturb_vector(vector)
        !
        ! this subroutine perturbs the input vector by a random vector
        !
        real(rp), intent(inout) :: vector(:)

        real(rp), parameter :: random_noise_scale = 1e-4_rp
        integer(ip) :: n_param
        real(rp), allocatable :: random_vector(:)
        real(rp), external :: dnrm2

        n_param = size(vector)
        allocate(random_vector(n_param))
        random_vector = 0.0_rp
        do while (dnrm2(n_param, random_vector, 1_ip) < numerical_zero)
            call random_number(random_vector)
            random_vector = 2.0_rp * random_vector - 1.0_rp
        end do
        vector = vector + random_noise_scale * dnrm2(n_param, vector, 1_ip) * &
                 random_vector / dnrm2(n_param, random_vector, 1_ip)
        deallocate(random_vector)

    end subroutine perturb_vector

    function get_preconditioned_residual_dot(residual, precond_residual, settings, &
                                             error) result(residual_dot)
        !
        ! this function returns the dot product of the residual and the preconditioned 
        ! residual, throws an error if the preconditioner is not positive definite
        !
        real(rp), intent(in) :: residual(:), precond_residual(:)
        type(solver_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error
        real(rp) :: residual_dot
        real(rp), external :: ddot

        residual_dot = ddot(size(residual), residual, 1_ip, precond_residual, 1_ip)
        if (abs(residual_dot) < numerical_zero) residual_dot = 0.0_rp
        if (residual_dot < 0.0_rp) then
            if (maxval(abs(precond_residual)) < epsilon(1.0_rp) * &
                maxval(abs(residual))) then
                residual_dot = 0.0_rp
            else
                call settings%log("The passed preconditioner function is not "// &
                                  "positive definite.", verbosity_error, .true.)
                error = 1
            end if
        end if

    end function get_preconditioned_residual_dot

    subroutine find_point_on_boundary(basis_vec, trust_radius, solution_basis_vec_dot, &
                                      basis_vec_dot, residual_dot, curvature, &
                                      solution_dot, solution, pred_func, step_size)
        !
        ! this subroutine finds the appropriate point on the trust-region boundary in
        ! the basis vector direction, the predicted function value at this point is 
        ! updated to reflect the new point
        !
        real(rp), intent(in) :: basis_vec(:), trust_radius, solution_basis_vec_dot, &
                                basis_vec_dot, residual_dot, curvature
        real(rp), intent(inout) :: solution_dot, solution(:), pred_func
        real(rp), intent(out) :: step_size

        real(rp), parameter :: rel_tol = 1e-12_rp
        real(rp) :: other_root

        call quadratic_roots(solution_dot - trust_radius ** 2, &
                             2.0_rp * solution_basis_vec_dot, basis_vec_dot, &
                             rel_tol, other_root, step_size)
        solution_dot = solution_dot + step_size * &
                       (2.0_rp * solution_basis_vec_dot + step_size * basis_vec_dot)
        solution = solution + step_size * basis_vec
        pred_func = pred_func + step_size * &
                    (0.5_rp * step_size * curvature - residual_dot)
                                                
    end subroutine find_point_on_boundary

    subroutine solve_tridiagonal_subproblem(n_red_space, diagonal, off_diagonal, &
                                            diagonal_fact, off_diagonal_fact, linear, &
                                            trust_radius, interior, try_warm, use_old, &
                                            lowest_eigenval, lambda, func, solution, &
                                            eigenvec, work, hard_case, &
                                            hard_case_step_size)
        !
        ! this subroutine determines a vector x which approximately minimizes the 
        ! quadratic function
        !
        !   func(solution) = 1/2 <solution, tridiagonal solution> + <linear, solution>
        !
        ! subject to the Euclidean norm constraint ||solution|| <= trust_radius.

        ! - computes an approximate solution and a Lagrange multiplier lambda such that 
        !   either lambda is zero and ||solution|| <= (1+rtol)*trust_radius, or lambda 
        !   is positive and | ||solution|| - trust_radius | <= rtol * trust_radius
        ! - if solution_sol is the solution to the problem, the approximate solution 
        !   satisfies func(solution) <= func(solution_sol) * (1 - rtol) ** 2
        ! - diagonal and off_diagonal: tridiagonal matrix
        ! - diagonal_fact and off_diagonal_fact: LDL.T factorization of the tridiagonal 
        !   matrix shifted by lambda
        ! - try_warm is true: an initial estimate of lambda should be provided
        ! - use_old is true: the lowest eigenvalue of the leading n-1 by n-1 block 
        !   should be provided 
        ! - interior is true: an interior solution is possible (interior will be set to 
        !   true if an interior solution was found)
        !
        integer(ip), intent(in) :: n_red_space
        real(rp), intent(in) :: diagonal(n_red_space), off_diagonal(n_red_space - 1), &
                                linear(n_red_space)
        logical, intent(in) :: use_old, try_warm
        real(rp), intent(in) :: trust_radius
        logical, intent(inout) :: interior
        real(rp), intent(inout) :: lambda, lowest_eigenval
        real(rp), intent(out) :: func, hard_case_step_size, &
                                 off_diagonal_fact(n_red_space - 1), &
                                 diagonal_fact(n_red_space), solution(n_red_space), &
                                 eigenvec(n_red_space), work(n_red_space)
        logical, intent(out) :: hard_case

        real(rp), parameter :: rel_tol = 1e-12_rp, mach_eps = epsilon(1.0_rp)
        integer(ip), parameter :: iter_max = 100
        real(rp) :: solution_norm, func_linear_term, norm_eigenvec_solution_dot, dist, &
                    delta_lambda, pert_l
        integer(ip) :: iter, indefinite
        real(rp), external :: dnrm2, ddot
        external :: dpttrf

        ! initialize variables
        hard_case = .false.
        hard_case_step_size = 0.0_rp
        pert_l = mach_eps ** 0.75

        ! find a guess for lambda unless solution is interior
        find_lambda_guess: block
            ! try a warm start
            if (try_warm) then
                ! attempt the Cholesky factorization of lambda-shifted tridiagonal
                diagonal_fact = diagonal + lambda
                off_diagonal_fact = off_diagonal
                call dpttrf(n_red_space, diagonal_fact, off_diagonal_fact, indefinite)

                ! if shifted tridiagonal is positive definite, solve 
                ! (tridiagonal + lambda * I) solution = -linear
                if (indefinite == 0) then
                    work = -linear
                    call solve_tridiagonal(n_red_space, diagonal, off_diagonal, &
                                           lambda, diagonal_fact, off_diagonal_fact, &
                                           work, solution, eigenvec, func_linear_term)

                    ! if the solution lies outside the trust-region, it provides a good 
                    ! initial estimate of the solution to the TR problem
                    solution_norm = dnrm2(n_red_space, solution, 1_ip)
                    if (abs(solution_norm - trust_radius) <= rel_tol * trust_radius) &
                        then
                        func = -0.5_rp * &
                               (func_linear_term + lambda * solution_norm ** 2)
                        return
                    end if
                    if (solution_norm > trust_radius) exit find_lambda_guess
                end if
            end if

            ! if the warm start fails, check for an unconstrained solution
            if (interior) then
                ! attempt the Cholesky factorization of tridiagonal (no shifting 
                ! necessary since we are checking for an interior solution)
                diagonal_fact = diagonal
                off_diagonal_fact = off_diagonal
                call dpttrf(n_red_space, diagonal_fact, off_diagonal_fact, indefinite)

                ! if tridiagonal is positive definite, solve  T x = -linear
                if (indefinite == 0) then
                    solution = -linear
                    call inverse_iteration(n_red_space, diagonal_fact, &
                                           off_diagonal_fact, solution)

                    ! if the solution lies within the trust-region, it provides an 
                    ! interior solution to the TR problem
                    solution_norm = dnrm2(n_red_space, solution, 1_ip)
                    if (solution_norm <= trust_radius) then
                        lambda = 0.0_rp
                        func = 0.5_rp * ddot(n_red_space, linear, 1_ip, solution, 1_ip)
                        return
                    ! find optimal Lagrange multiplier with Newton's method
                    else
                        lambda = 0.0_rp
                    end if
                ! tridiagonal is indefinite
                else
                    interior = .false.
                end if
            ! no interior solution possible, tridiagonal must be indefinite
            else
                indefinite = 1
            end if

            ! the solution is not interior, compute the lowest eigenvalue
            if (indefinite > 0) then
                lowest_eigenval = &
                    get_tridiagonal_lowest_eigenvalue(n_red_space, diagonal, &
                                                      off_diagonal, rel_tol, use_old, &
                                                      lowest_eigenval)
                lowest_eigenval = min(lowest_eigenval, 0.0_rp)

                ! construct a Lagrange multiplier to ensure that shifted tridiagonal is 
                ! positive definite
                if (lowest_eigenval <= 0.0_rp) then
                    lambda = -lowest_eigenval * (1.0_rp + pert_l) + pert_l
                ! construct a Lagrange multiplier which can become negative for large 
                ! trust regions while keeping the shifted tridiagonal positive definite
                else
                    lambda = -lowest_eigenval * (1.0_rp - pert_l) + pert_l
                end if

                ! loop until Lagrange multiplier is found which makes lambda-shifted 
                ! tridiagonal positive definite
                do
                    ! attempt the Cholesky factorization of lambda-shifted tridiagonal
                    diagonal_fact = diagonal + lambda
                    off_diagonal_fact = off_diagonal
                    call dpttrf(n_red_space, diagonal_fact, off_diagonal_fact, &
                                indefinite)
                    if (indefinite == 0) exit

                    ! shifted tridiagonal is still numerically indefinite and must be 
                    ! perturbed a bit more
                    pert_l = 2.0_rp * pert_l
                    if (lowest_eigenval <= 0.0_rp) then
                        lambda = lambda * (1.0_rp + pert_l) + pert_l
                    else
                        lambda = lambda * (1.0_rp - pert_l) + pert_l
                    end if
                end do

                ! solve T x = -linear
                work = -linear
                call solve_tridiagonal(n_red_space, diagonal, off_diagonal, lambda, &
                                       diagonal_fact, off_diagonal_fact, work, &
                                       solution, eigenvec, func_linear_term)
                solution_norm = dnrm2(n_red_space, solution, 1_ip)

                ! if the step length stays below trust radius even as the exact 
                ! eigenvalue is approached by the Lagrange multiplier, the gradient 
                ! must be orthogonal to the corresponding eigenvector
                if (solution_norm < trust_radius) then
                    ! hard case is occuring
                    hard_case = .true.

                    ! lambda can be lowest eigenvalue
                    lambda = -lowest_eigenval

                    ! compute lowest eigenvector
                    call get_tridiagonal_lowest_eigenvector(n_red_space, &
                                                            lowest_eigenval, &
                                                            diagonal, off_diagonal, &
                                                            diagonal_fact, &
                                                            off_diagonal_fact, &
                                                            eigenvec)

                    ! get normalized overlap of eigenvector and solution (should be 
                    ! close to zero in hard case)
                    norm_eigenvec_solution_dot = &
                        ddot(n_red_space, eigenvec, 1_ip, solution, 1_ip) / trust_radius

                    ! get distance to trust-region boundary
                    dist = (trust_radius - solution_norm) * &
                           ((trust_radius + solution_norm) / trust_radius)
                    
                    ! compute the step size so that solution + hard_case_step_size * 
                    ! eigenvec lies on the trust-region boundary
                    hard_case_step_size = &
                        sign(dist / (abs(norm_eigenvec_solution_dot) + &
                                     sqrt(norm_eigenvec_solution_dot ** 2 + dist / &
                                          trust_radius)), norm_eigenvec_solution_dot)
                    func = -0.5_rp * (func_linear_term + lambda * trust_radius ** 2)
                    solution_norm = dnrm2(n_red_space, solution + hard_case_step_size &
                                          * eigenvec, 1_ip)
                    return
                end if
            else
                lowest_eigenval = 0.0_rp
            end if

        end block find_lambda_guess

        ! apply Newton's method starting from lambda guess
        do iter = 2, iter_max
            ! compute the Newton correction
            work = solution / solution_norm
            call forward_substitution(n_red_space, off_diagonal_fact, work)
            delta_lambda = ((solution_norm - trust_radius) / trust_radius) / &
                           ddot(n_red_space, work, 1_ip, work / diagonal_fact, 1_ip)

            ! check that the Newton correction is significant, otherwise return since 
            ! no further progress can be made
            if (abs(delta_lambda) < mach_eps * abs(lambda)) then
                func = -0.5_rp * (func_linear_term + lambda * solution_norm ** 2)
                return
            end if

            ! compute the new estimate of lambda
            lambda = lambda + delta_lambda

            ! find the Cholesky factorization of shifted tridiagonal
            diagonal_fact = diagonal + lambda
            off_diagonal_fact = off_diagonal
            call dpttrf(n_red_space, diagonal_fact, off_diagonal_fact, indefinite)

            ! solve the equation (tridiagonal + lambda * I) solution = -linear
            work = -linear
            call solve_tridiagonal(n_red_space, diagonal, off_diagonal, lambda, &
                                   diagonal_fact, off_diagonal_fact, work, solution, &
                                   eigenvec, func_linear_term)
            solution_norm = dnrm2(n_red_space, solution, 1_ip)

            ! test for convergence
            if (solution_norm - trust_radius <= rel_tol * trust_radius) then
                func = -0.5_rp * (func_linear_term + lambda * solution_norm ** 2)
                return
            end if
        end do

        ! could not converge in maximum number of iterations, return best available 
        ! approximation
        func = -0.5_rp * (func_linear_term + lambda * solution_norm ** 2)
        return

    end subroutine solve_tridiagonal_subproblem

    function get_tridiagonal_lowest_eigenvalue(n_elem, diagonal, off_diagonal, &
                                               rel_tol, use_old, old_lowest_eigenval) &
        result(lowest_eigenvalue)
        !
        ! this function computes the lowest eigenvalue of a symmetric tridiagonal 
        ! matrix
        !
        integer(ip), intent(in) :: n_elem
        real(rp), intent(in) :: rel_tol, old_lowest_eigenval
        logical, intent(in) :: use_old
        real(rp), intent(in) :: diagonal(n_elem), off_diagonal(n_elem - 1)

        real(rp), parameter :: perturb = 1e-6_rp, mach_eps = epsilon(mach_eps), &
                               tol = mach_eps ** 0.66
        integer(ip) :: i, n_neg_pivots
        real(rp) :: lower, upper, tol_interval, pivot, pivot_derivative, infinity, &
                    coeff_b, coeff_c, e_trial, root1, root2
        real(rp) :: lowest_eigenvalue

        ! special case: n_elem = 1
        if (n_elem == 1) then
            lowest_eigenvalue = diagonal(1)
            return
        end if

        ! initialize lower and upper bounds using Gersgorin bounds
        lower = min(diagonal(1) - abs(off_diagonal(1)), &
                    diagonal(n_elem) - abs(off_diagonal(n_elem - 1)), &
                    minval(diagonal(2:n_elem - 1) - abs(off_diagonal(:n_elem - 2)) - &
                           abs(off_diagonal(2:))))
        upper = max(diagonal(1) + abs(off_diagonal(1)), &
                    diagonal(n_elem) + abs(off_diagonal(n_elem - 1)), &
                    maxval(diagonal(2:n_elem - 1) + abs(off_diagonal(:n_elem - 2)) + &
                           abs(off_diagonal(2:))))

        ! initialize eigenvalue starting guess from guess or from bounds
        infinity = 1.0_rp + upper - lower
        if (use_old) then
            upper = min(upper, old_lowest_eigenval)
            lowest_eigenvalue = upper - perturb
        else
            lowest_eigenvalue = lower
        end if
        tol_interval = tol * (1.0_rp + 0.5_rp * abs(lower) + abs(upper))

        ! main iteration loop
        iter_loop: do
            ! compute the inertia of T - lowest_eigenvalue * I by implicitly factoring 
            ! the matrix
            n_neg_pivots = 0
            pivot = diagonal(1) - lowest_eigenvalue
            pivot_derivative = -1.0_rp

            ! if a zero pivot is encountered, reset the upper bound
            if (abs(pivot) < numerical_zero) then
                upper = lowest_eigenvalue
                lowest_eigenvalue = 0.5_rp * (lower + upper)
                cycle iter_loop
            ! if a negative pivot is encountered, exit
            else if (pivot < 0.0_rp) then
                n_neg_pivots = 1
            end if

            do i = 2, n_elem
                ! update the pivot
                pivot_derivative = -1.0_rp + pivot_derivative * &
                                   (off_diagonal(i - 1) / pivot) ** 2
                pivot = (diagonal(i) - (off_diagonal(i - 1) ** 2) / pivot) - &
                        lowest_eigenvalue

                ! check for zero pivot
                if (abs(pivot) < numerical_zero) then
                    ! return if a zero last pivot is encountered
                    if (n_neg_pivots == 0 .and. i == n_elem) then
                        return
                    ! reduce the upper bound if a zero pivot is encountered
                    else
                        upper = lowest_eigenvalue
                        lowest_eigenvalue = 0.5_rp * (lower + upper)
                        cycle iter_loop
                    end if
                ! increment number of negative pivots
                else if (pivot < 0.0_rp) then
                    n_neg_pivots = n_neg_pivots + 1
                    ! exit if more than one negative pivot is encountered
                    if (n_neg_pivots > 1) then
                        pivot = infinity
                        pivot_derivative = 1.0_rp
                        exit
                    end if
                end if
            end do

            ! increase the lower bound
            if (n_neg_pivots == 0) then
                lower = lowest_eigenvalue
            ! reduce the upper bound
            else
                upper = lowest_eigenvalue
            end if

            ! test for convergence
            if (abs(pivot) < tol .OR. upper - lower < tol_interval) then
                return
            end if

            ! compute the Newton step
            if (use_old) then
                coeff_b = 2.0_rp * lowest_eigenvalue + pivot + &
                          (lowest_eigenvalue - old_lowest_eigenval) * pivot_derivative
                coeff_c = -(lowest_eigenvalue - old_lowest_eigenval) * pivot + &
                          lowest_eigenvalue * coeff_b - lowest_eigenvalue ** 2
                call quadratic_roots(coeff_c, -coeff_b, 1.0_rp, rel_tol, root1, root2)
                e_trial = root1
            else
                e_trial = lowest_eigenvalue - pivot / pivot_derivative
            end if

            ! if the estimate lies in the interval (lower, e_2) and the Newton step
            ! continues to lie in [lower, upper], use the Newton step as the next 
            ! estimate
            if (n_neg_pivots <= 1 .and. (e_trial > lower .and. e_trial < upper)) then
                lowest_eigenvalue = e_trial
            ! otherwise bisect the bounds to get the new eigenvalue estimate
            else
                lowest_eigenvalue = 0.5_rp * (lower + upper)
            end if

        end do iter_loop

    end function get_tridiagonal_lowest_eigenvalue

    subroutine get_tridiagonal_lowest_eigenvector(n_elem, eigenval_est, diagonal, &
                                                  off_diagonal, diagonal_fact, &
                                                  off_diagonal_fact, eigenvec)
        !
        ! this subroutine computes an eigenvector corresponding to the lowest 
        ! eigenvalue of a symmetric tridiagonal matrix using inverse iteration
        !
        integer(ip), intent(in) :: n_elem
        real(rp), intent(in) :: eigenval_est
        real(rp), intent(in) :: diagonal(n_elem), off_diagonal(n_elem - 1)
        real(rp), intent(out):: diagonal_fact(n_elem), off_diagonal_fact(n_elem - 1), &
                                eigenvec(n_elem)

        real(rp), parameter :: perturb_factor = 1e-6_rp, conv_tol = 1e-8_rp
        integer(ip), parameter :: max_iter = 5
        integer(ip) :: indefinite, iter
        real(rp) :: wnorm, perturb
        real(rp), external :: dnrm2
        external :: dpttrf

        ! perturb eigenvalue estimate until shifted tridiagonal matrix becomes positive 
        ! definite
        perturb = perturb_factor * (1.0_rp - eigenval_est)
        do
            ! construct shifted tridiagonal matrix
            diagonal_fact = diagonal - (eigenval_est - perturb)
            if (n_elem > 1) off_diagonal_fact = off_diagonal

            ! attempt the Cholesky factorization of T - eigenval * I 
            call dpttrf(n_elem, diagonal_fact, off_diagonal_fact, indefinite)

            ! exit if shifted tridiagonal matrix is positive definite
            if (indefinite == 0) exit

            ! increase perturbation
            perturb = 10.0_rp * perturb
        end do

        ! initialize a random initial estimate of eigenvector
        call random_number(eigenvec)

        ! use inverse iteration to solve (T - eigenval_est I) eigenvec_new = eigenvec 
        ! which should converge very quickly since eigenvalue estimate is accurate
        do iter = 1, max_iter
            ! solve (T - eigenval_est I) eigenvec_new = eigenvec
            call inverse_iteration(n_elem, diagonal_fact, off_diagonal_fact, eigenvec)

            ! normalize eigenvector
            wnorm = 1.0_rp / dnrm2(n_elem, eigenvec, 1_ip)
            eigenvec = eigenvec * wnorm

            ! check for convergence
            if (abs(wnorm - perturb) <= conv_tol) exit
        end do

    end subroutine get_tridiagonal_lowest_eigenvector

    subroutine solve_tridiagonal(n_elem, diagonal, off_diagonal, lambda, &
                                 diagonal_fact, off_diagonal_fact, rhs, solution, &
                                 work, scaled_solution_norm_squared)
        !
        ! this subroutine solves the system (T + lambda I) x = b using iterative 
        ! refinement given the factors of the tridiagonal matrix T + lambda I
        !
        integer(ip), intent(in) :: n_elem
        real(rp), intent(in) :: lambda
        real(rp), intent(in) :: off_diagonal(n_elem - 1), &
                                off_diagonal_fact(n_elem - 1), diagonal(n_elem), &
                                diagonal_fact(n_elem), rhs(n_elem)
        real(rp), intent(out) :: solution(n_elem), work(n_elem), &
                                 scaled_solution_norm_squared

        integer(ip), parameter :: itmax = 1
        integer(ip) :: i, it
        real(rp), external :: ddot

        ! use inverse iteration to solve (T - lambda I) solution = rhs which should 
        ! converge very quickly since eigenvalue estimate is accurate
        work = rhs

        ! solve (T + lambda I) solution = rhs with inverse iteration
        call forward_substitution(n_elem, off_diagonal_fact, work)
        solution = work
        call diagonal_scaling(n_elem, diagonal_fact, solution)
        scaled_solution_norm_squared = ddot(n_elem, work, 1_ip, solution, 1_ip)
        call backward_substitution(n_elem, off_diagonal_fact, solution)

        ! start of iterative refinement, only a single iteration is used since this is 
        ! enough to achieve floating-point precision for well-conditioned problems
        do it = 1, itmax
            ! compute the residual r = b - (T + lambda I) x
            work = rhs - (diagonal + lambda) * solution
            do i = 1, n_elem - 1
                work(i) = work(i) - off_diagonal(i) * solution(i + 1)
                work(i + 1) = work(i + 1) - off_diagonal(i) * solution(i)
            end do

            ! solve (T + lambda I) dx = r with inverse iteration
            call inverse_iteration(n_elem, diagonal_fact, off_diagonal_fact, work)

            ! update solution
            solution = solution + work
        end do

    end subroutine solve_tridiagonal

    subroutine inverse_iteration(n, diagonal, off_diagonal, arr)
        !
        ! this subroutine performs in-place inverse iteration for a unit lower 
        ! bidiagonal matrix
        !
        integer(ip), intent(in) :: n
        real(rp), intent(in) :: diagonal(n), off_diagonal(n - 1)
        real(rp), intent(inout) :: arr(n)

        call forward_substitution(n, off_diagonal, arr)
        call diagonal_scaling(n, diagonal, arr)
        call backward_substitution(n, off_diagonal, arr)

    end subroutine inverse_iteration

    subroutine forward_substitution(n, off_diagonal, arr)
        !
        ! this subroutine performs in-place forward substitution L x = b for a unit 
        ! lower bidiagonal matrix 
        !
        integer(ip), intent(in) :: n
        real(rp), intent(in) :: off_diagonal(n - 1)
        real(rp), intent(inout) :: arr(n)

        integer(ip) :: i

        do i = 1, n - 1
            arr(i + 1) = arr(i + 1) - off_diagonal(i) * arr(i)
        end do

    end subroutine forward_substitution

    subroutine diagonal_scaling(n, diagonal, arr)
        !
        ! this subroutine performs in-place diagonal scaling D * x = b for a diagonal 
        ! matrix
        !
        integer(ip), intent(in) :: n
        real(rp), intent(in) :: diagonal(n)
        real(rp), intent(inout) :: arr(n)

        arr = arr / diagonal

    end subroutine diagonal_scaling

    subroutine backward_substitution(n, off_diagonal, arr)
        !
        ! this subroutine performs in-place backward substitution L.T x = b for a 
        ! unit lower bidiagonal matrix
        !
        integer(ip), intent(in) :: n
        real(rp), intent(in) :: off_diagonal(n - 1)
        real(rp), intent(inout) :: arr(n)

        integer(ip) :: i

        do i = n - 1, 1, -1
            arr(i) = arr(i) - off_diagonal(i) * arr(i + 1)
        end do

    end subroutine backward_substitution

    subroutine quadratic_roots(a0, a1, a2, tol, root1, root2)
        !
        ! this subroutine finds the number and values of real roots of an quadratic 
        ! equation (a2 * x**2 + a1 * x + a0 = 0) where a0, a1 and a2 are real
        !
        real(rp), intent(in) :: a2, a1, a0, tol
        real(rp), intent(out) :: root1, root2

        integer(ip) :: nroots
        real(rp) :: rhs, intermediate, quadratic, quadratic_derivative

        rhs = tol * a1 * a1
        ! function is quadratic
        if (abs(a0 * a2) > rhs) then
            root2 = a1 * a1 - 4.0_rp * a2 * a0
            ! numerical double root
            if (abs(root2) <= (epsilon(1.0_rp) * a1) ** 2) then
                nroots = 2
                root1 = -0.5_rp * a1 / a2
                root2 = root1
            ! complex not real roots
            else if (root2 < 0.0_rp) then
                nroots = 0
                root1 = 0.0_rp
                root2 = 0.0_rp
            ! distinct real roots
            else
                intermediate = -0.5_rp * (a1 + sign(sqrt(root2), a1))
                nroots = 2
                root1 = intermediate / a2
                root2 = a0 / intermediate
                if (root1 > root2) then
                    intermediate = root1
                    root1 = root2
                    root2 = intermediate
                end if
            end if
        ! function is lower-order polynomial
        else if (abs(a2) < numerical_zero) then
            if (abs(a1) < numerical_zero) then
                ! function is zero
                if (abs(a0) < numerical_zero) then
                    nroots = 1
                    root1 = 0.0_rp
                    root2 = 0.0_rp
                ! function is constant
                else
                    nroots = 0
                    root1 = 0.0_rp
                    root2 = 0.0_rp
                end if
            ! function is linear
            else
                nroots = 1
                root1 = -a0 / a1
                root2 = 0.0_rp
            end if
        ! function isvery ill-conditioned quadratic
        else
            nroots = 2
            if (-a1 / a2 > 0.0_rp) then
                root1 = 0.0_rp
                root2 = -a1 / a2
            else
                root1 = -a1 / a2
                root2 = 0.0_rp
            end if
        end if

        ! perform a Newton iteration to ensure that the roots are accurate
        if (nroots >= 1) then
            quadratic = (a2 * root1 + a1) * root1 + a0
            quadratic_derivative = 2.0_rp * a2 * root1 + a1
            if (abs(quadratic_derivative) > 0.0_rp) then
                root1 = root1 - quadratic / quadratic_derivative
                quadratic = (a2 * root1 + a1) * root1 + a0
            end if
            if (nroots == 2) then
                quadratic = (a2 * root2 + a1) * root2 + a0
                quadratic_derivative = 2.0_rp * a2 * root2 + a1
                if (abs(quadratic_derivative) > 0.0_rp) then
                    root2 = root2 - quadratic / quadratic_derivative
                    quadratic = (a2 * root2 + a1) * root2 + a0
                end if
            end if
        end if

    end subroutine quadratic_roots

end module opentrustregion
