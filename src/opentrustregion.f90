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

    ! mathematical constants
    real(rp), parameter :: pi = 4.d0*atan(1.d0)

    ! define default optional arguments
    logical, parameter :: solver_stability_default = .true., &
                          solver_line_search_default = .false.
    real(rp), parameter :: solver_conv_tol_default = 1d-5, &
                           solver_start_trust_radius_default = 0.4, &
                           solver_global_red_factor_default = 1d-3, &
                           solver_local_red_factor_default = 1d-4, &
                           stability_conv_tol_default = 1d-4
    integer(ip), parameter :: solver_n_random_trial_vectors_default = 1, &
                              solver_n_macro_default = 150, &
                              solver_n_micro_default = 100, &
                              solver_verbose_default = 0, solver_seed_default = 42, &
                              stability_n_random_trial_vectors_default = 20, &
                              stability_n_iter_default = 100, &
                              stability_verbose_default = 0

    ! derived type for solver settings
    type :: solver_settings_type
        logical :: stability, line_search
        real(rp) :: conv_tol, start_trust_radius, global_red_factor, local_red_factor
        integer(ip) :: n_random_trial_vectors, n_macro, n_micro, verbose, seed
    contains
        procedure :: init_solver_settings
    end type

    type :: stability_settings_type
        real(rp) :: conv_tol
        integer(ip) :: n_random_trial_vectors, n_iter, verbose
    contains
        procedure :: init_stability_settings
    end type

    ! interface for setting default values
    interface set_default
        module procedure set_default_real, set_default_integer, set_default_logical
    end interface

    ! define global variables
    integer(ip) :: tot_orb_update = 0, tot_hess_x = 0

    ! interfaces for callback functions
    abstract interface
        function hess_x_type(x) result(hess_x)
            import :: rp

            real(rp), intent(in) :: x(:)

            real(rp) :: hess_x(size(x))
        end function hess_x_type
    end interface

    abstract interface
        subroutine update_orbs_type(kappa, func, grad, h_diag, hess_x_funptr)
            import :: rp, hess_x_type

            real(rp), intent(in) :: kappa(:)

            real(rp), intent(out) :: func, grad(:), h_diag(:)
            procedure(hess_x_type), pointer, intent(out) :: hess_x_funptr
        end subroutine update_orbs_type
    end interface

    abstract interface
        function obj_func_type(kappa) result(func)
            import :: rp

            real(rp), intent(in) :: kappa(:)

            real(rp) :: func
        end function obj_func_type
    end interface

contains

    subroutine solver(update_orbs, obj_func, n_param, stability, line_search, &
                      conv_tol, n_random_trial_vectors, start_trust_radius, n_macro, &
                      n_micro, global_red_factor, local_red_factor, verbose, seed)
        !
        ! this subroutine is the main solver for orbital optimization
        !
        procedure(update_orbs_type), pointer, intent(in) :: update_orbs
        procedure(obj_func_type), intent(in), pointer :: obj_func
        integer(ip), intent(in) :: n_param
        logical, intent(in), optional :: stability, line_search
        real(rp), intent(in), optional :: conv_tol, start_trust_radius, &
                                          global_red_factor, local_red_factor
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_macro, n_micro, &
                                             verbose, seed

        type(solver_settings_type) :: settings
        real(rp) :: trust_radius, func, grad_norm, grad_rms, mu, residual_norm, &
                    red_factor, initial_residual_norm, new_func, ratio, n_kappa
        real(rp), dimension(n_param) :: kappa, grad, h_diag, solution, &
                                        last_solution_normalized, h_solution, &
                                        residual, solution_normalized, precond, &
                                        basis_vec
        real(rp), allocatable :: red_space_basis(:, :), h_basis(:, :), aug_hess(:, :), &
                                 red_space_solution(:), red_hess_vec(:)
        logical :: macro_converged = .false., stable = .true., accept_step, &
                   micro_converged, newton
        integer(ip) :: imacro, imicro, i, ntrial, initial_imicro
        integer(ip), parameter :: stability_n_points = 21
        procedure(hess_x_type), pointer :: hess_x_funptr
        character(:), allocatable :: error_flag
        real(rp), external :: dnrm2, ddot
        external :: dgemm, dgemv

        ! initialize settings
        call settings%init_solver_settings(stability, line_search, conv_tol, &
                                           n_random_trial_vectors, start_trust_radius, &
                                           n_macro, n_micro, global_red_factor, &
                                           local_red_factor, verbose, seed)

        ! initialize random number generator
        call init_rng(settings%seed)

        ! initialize orbital rotation matrix
        kappa = 0.d0

        ! check that number of random trial vectors is below number of parameters
        if (settings%n_random_trial_vectors > n_param/2) then
            settings%n_random_trial_vectors = n_param/2
            if (settings%verbose > 1) write (stderr, '(A, I0, A)') " Number of "// &
                "random trial vectors should be smaller than half the number of "// &
                "parameters. Setting to ", settings%n_random_trial_vectors, "."
        end if

        ! initialize starting trust radius
        trust_radius = settings%start_trust_radius

        ! print header
        if (settings%verbose > 2) then
            write (stdout, *) repeat("-", 109)
            write (stdout, *) " Iteration |     Objective function     |", &
                " Gradient RMS | Level shift |   Micro    | Trust radius | Step size "
            write (stdout, *) "           |                            |", &
                "              |             | iterations |              |           "
            write (stdout, *) repeat("-", 109)
        end if

        do imacro = 1, settings%n_macro
            ! calculate cost function, gradient and Hessian diagonal
            call update_orbs(kappa, func, grad, h_diag, hess_x_funptr)

            ! sanity check for array size
            if (imacro == 1 .and. size(grad) /= n_param) then
                error stop "Size of gradient array returned by subroutine "// &
                    "update_orbs does not equal number of parameters"
            end if

            ! calculate gradient norm
            grad_norm = dnrm2(n_param, grad, 1)

            ! calculate RMS gradient
            grad_rms = grad_norm/sqrt(real(n_param, kind=rp))

            ! log results
            if (settings%verbose > 2) then
                if (imacro == 1) then
                    write (stdout, '(6X, I3, 3X, "|", 4X, 1PE21.14, 3X, "|", 2X, '// &
                           '1PE9.2, 3X, "|", 6X, "-", 6X, "|", 8X, "-", 3X, "|", '// &
                           '8X, "-", 5X, "|", 6X, "-", 3X)') imacro - 1, func, grad_rms
                else if (.not. stable) then
                    write (stdout, '(6X, I3, 3X, "|", 4X, 1PE21.14, 3X, "|", 2X, '// &
                           '1PE9.2, 3X, "|", 6X, "-", 6X, "|", 8X, "-", 3X, "|", '// &
                           '8X, "-", 5X, "|", X, 1PE9.2)') imacro - 1, func, &
                        grad_rms, dnrm2(n_param, kappa, 1)
                    stable = .true.
                else
                    write (stdout, '(6X, I3, 3X, "|", 4X, 1PE21.14, 3X, "|", 2X, '// &
                           '1PE9.2, 3X, "|", 2X, 1PE9.2, 2X, "|", 6X, I3, 3X, "|", '// &
                           '3X, 1PE9.2, 2X, "|", X, 1PE9.2)') imacro - 1, func, &
                        grad_rms, -mu, imicro, trust_radius, dnrm2(n_param, kappa, 1)
                end if
            end if

            ! check for convergence and stability
            if (grad_rms < settings%conv_tol) then
                ! always perform stability check if starting at stationary point
                if (settings%stability .or. imacro == 1) then
                    call stability_check(grad, h_diag, hess_x_funptr, stable, kappa, &
                                         verbose=settings%verbose)
                    if (.not. stable) then
                        ! logarithmic line search
                        do i = 1, stability_n_points
                            n_kappa = 10.d0**(-(i - 1)/ &
                                              real(stability_n_points - 1, rp)*10.d0)
                            new_func = obj_func(n_kappa*kappa)
                            if (new_func < func) then
                                kappa = n_kappa*kappa
                                exit
                            end if
                        end do
                        if (new_func >= func) then
                            error stop "Line search was unable to find lower "// &
                                "objective function along unstable mode."
                        else if (imacro == 1) then
                            if (settings%verbose > 1) then
                                write (stderr, *) "Started at saddle point. The "// &
                                    "algorithm will continue by moving along "// &
                                    "eigenvector direction "
                                write (stderr, *) "corresponding to negative "// &
                                    "eigenvalue."
                            end if
                        else
                            if (settings%verbose > 1) then
                                write (stderr, *) "Reached saddle point. This is "// &
                                    "likely due to symmetry and can be avoided by "// &
                                    "increasing the number of random "
                                write (stderr, *) "trial vectors. The algorithm "// &
                                    "will continue by moving along eigenvector "
                                write (stderr, *) "direction corresponding to "// &
                                    "negative eigenvalue."
                            end if
                        end if
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

            ! generate trial vectors
            red_space_basis = generate_trial_vectors(settings%n_random_trial_vectors, &
                                                     grad, grad_norm, h_diag)

            ! number of trial vectors
            ntrial = size(red_space_basis, 2)

            ! calculate linear transformations of basis vectors
            allocate (h_basis(n_param, ntrial))
            do i = 1, ntrial
                h_basis(:, i) = hess_x_funptr(red_space_basis(:, i))
            end do

            ! construct augmented Hessian in reduced space
            allocate (aug_hess(ntrial + 1, ntrial + 1))
            aug_hess = 0.d0
            call dgemm("T", "N", ntrial, ntrial, n_param, 1.d0, red_space_basis, &
                       n_param, h_basis, n_param, 0.d0, aug_hess(2:, 2:), ntrial)

            ! allocate space for reduced space solution
            allocate (red_space_solution(size(red_space_basis, 2)))

            ! decrease trust radius until micro iterations converge and step is
            ! accepted
            last_solution_normalized = 0.d0
            accept_step = .false.
            do while (.not. accept_step)
                micro_converged = .false.
                do imicro = 1, settings%n_micro
                    ! do a Newton step if the model is positive definite and the step
                    ! is within the trust region
                    newton = .false.
                    if (min_eigval(aug_hess(2:, 2:)) > -1.d-5) then
                        call newton_step(aug_hess, grad_norm, red_space_basis, &
                                         solution, red_space_solution)
                        mu = 0.d0
                        if (dnrm2(n_param, solution, 1) < trust_radius) newton = .true.
                    end if

                    ! otherwise perform bisection to find the level shift
                    if (.not. newton) then
                        call bisection(aug_hess, grad_norm, red_space_basis, &
                                       trust_radius, solution, red_space_solution, mu, &
                                       error_flag)
                        if (error_flag == "Target trust region outside of "// &
                            "bracketing range.") exit
                    end if

                    ! calculate Hessian linear transformation of solution
                    call dgemv("N", n_param, size(h_basis, 2), 1.d0, h_basis, n_param, &
                               red_space_solution, 1, 0.d0, h_solution, 1)

                    ! calculate residual
                    residual = grad + h_solution - mu*solution

                    ! calculate residual norm
                    residual_norm = dnrm2(n_param, residual, 1)

                    ! determine reduction factor depending on whether local region is
                    ! reached
                    if (abs(mu) < 1.d-12) then
                        red_factor = settings%local_red_factor
                    else
                        red_factor = settings%global_red_factor
                    end if

                    ! get normalized solution vector
                    solution_normalized = solution/dnrm2(n_param, solution, 1)

                    ! reset initial residual norm if solution changes
                    if (ddot(n_param, last_solution_normalized, 1, &
                             solution_normalized, 1)**2 < 0.5d0) then
                        initial_imicro = imicro
                        initial_residual_norm = residual_norm
                    end if

                    ! check if micro iterations have converged
                    if (residual_norm < max(red_factor*grad_norm, 1d-12)) then
                        micro_converged = .true.
                        exit
                    else if (imicro - initial_imicro >= 10 .and. &
                            residual_norm > 0.8*initial_residual_norm) then
                        exit
                    end if

                    ! save current solution
                    last_solution_normalized = solution_normalized

                    ! precondition residual
                    precond = h_diag - mu
                    where (abs(precond) < 1e-10)
                        precond = 1e-10
                    end where
                    residual = residual/precond

                    ! orthonormalize to current orbital space to get new basis vector
                    basis_vec = gram_schmidt(residual, red_space_basis)

                    ! add new trial vector to orbital space
                    call add_column(red_space_basis, basis_vec)

                    ! add linear transformation of new basis vector
                    call add_column(h_basis, hess_x_funptr(basis_vec))

                    ! construct new augmented Hessian
                    allocate (red_hess_vec(size(red_space_basis, 2) + 1))
                    red_hess_vec(1) = 0.d0
                    call dgemv("T", n_param, size(red_space_basis, 2), 1.d0, &
                               red_space_basis, n_param, &
                               h_basis(:, size(red_space_basis, 2)), 1, 0.d0, &
                               red_hess_vec(2:), 1)
                    call extend_symm_matrix(aug_hess, red_hess_vec)
                    deallocate (red_hess_vec)

                    ! reallocate reduced space solution
                    deallocate (red_space_solution)
                    allocate (red_space_solution(size(red_space_basis, 2)))
                end do

                ! evaluate function at predicted point
                new_func = obj_func(solution)

                ! calculate ratio of evaluated function and predicted function
                ratio = (new_func - func)/ddot(n_param, solution, 1, &
                                               grad + 0.5*h_solution, 1)

                ! decrease trust radius if micro iterations are unable to converge, if
                ! function value has not decreased or if individual orbitals change too
                ! much
                if (.not. micro_converged .or. ratio < 0.d0 .or. &
                    any(abs(solution) > pi/4)) then
                    trust_radius = 0.7d0*trust_radius
                    accept_step = .false.
                    if (trust_radius < 1.d-10) error stop "Trust radius too small."
                    ! check if step is too long
                else if (ratio < 0.25d0) then
                    trust_radius = 0.7d0*trust_radius
                    accept_step = .true.
                    ! check if quadratic approximation is valid
                else if (ratio < 0.75d0) then
                    accept_step = .true.
                    ! check if step is potentially too short
                else
                    trust_radius = 1.2d0*trust_radius
                    accept_step = .true.
                end if

                ! flush output
                flush(stdout)
                flush(stderr)
            end do

            ! increment total number of Hessian linear transformations
            tot_hess_x = tot_hess_x + size(red_space_basis, 2)

            ! deallocate quantities from microiterations
            deallocate (red_space_solution)
            deallocate (aug_hess)
            deallocate (h_basis)
            deallocate (red_space_basis)

            ! perform line search
            if (settings%line_search) then
                n_kappa = bracket(obj_func, solution, 0.d0, 1.d0)
            else
                n_kappa = 1.d0
            end if

            ! set orbital rotation
            kappa = n_kappa*solution

        end do

        ! increment total number of orbital updates
        tot_orb_update = tot_orb_update + imacro

        ! stop if no convergence
        if (.not. macro_converged) error stop "Orbital optimization has not converged!"

        ! finish logging
        if (settings%verbose > 2) then
            write (stdout, *) repeat("-", 109)
            write (stdout, '(A, I0)') " Total number of Hessian linear "// &
                "transformations: ", tot_hess_x
            write (stdout, '(A, I0)') " Total number of orbital updates: ", &
                tot_orb_update
        end if

        ! flush output
        flush (stdout)
        flush (stderr)

    end subroutine solver

    subroutine stability_check(grad, h_diag, hess_x_funptr, stable, kappa, conv_tol, &
                               n_random_trial_vectors, n_iter, verbose)
        !
        ! this subroutine performs a stability check
        !
        real(rp), intent(in) :: grad(:), h_diag(:)
        procedure(hess_x_type), pointer, intent(in) :: hess_x_funptr
        logical, intent(out) :: stable
        real(rp), intent(out) :: kappa(:)
        real(rp), intent(in), optional :: conv_tol
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_iter, verbose

        type(stability_settings_type) :: settings
        integer(ip) :: n_param, ntrial, i, iter
        real(rp), dimension(size(grad)) :: full_grad, full_h_diag, solution, &
                                           h_solution, residual, precond, basis_vec
        real(rp) :: full_grad_norm, eigval
        real(rp), allocatable :: red_space_basis(:, :), h_basis(:, :), &
                                 red_space_hess(:, :), red_space_solution(:), &
                                 red_space_hess_vec(:)
        real(rp), external :: dnrm2
        external :: dgemm, dgemv

        ! initialize settings
        call settings%init_stability_settings(conv_tol, n_random_trial_vectors, &
                                              n_iter, verbose)

        ! get number of parameters
        n_param = size(grad)

        ! check that number of random trial vectors is below number of parameters
        if (settings%n_random_trial_vectors > n_param/2) then
            settings%n_random_trial_vectors = n_param/2
            if (settings%verbose > 1) write (stderr, '(A, I0, A)') " Number of "// &
                "random trial vectors should be smaller than half the number of "// &
                "parameters. Setting to ", settings%n_random_trial_vectors, "."
        end if

        ! get quantities for full Hessian
        full_grad = 2.d0*grad
        full_h_diag = 2.d0*h_diag
        full_grad_norm = dnrm2(n_param, full_grad, 1)

        ! generate trial vectors
        red_space_basis = generate_trial_vectors(settings%n_random_trial_vectors, &
                                                 full_grad, full_grad_norm, full_h_diag)

        ! number of trial vectors
        ntrial = size(red_space_basis, 2)

        ! calculate linear transformations of basis vectors
        allocate (h_basis(n_param, ntrial))
        do i = 1, ntrial
            h_basis(:, i) = full_hess_x(red_space_basis(:, i))
        end do

        ! construct augmented Hessian in reduced space
        allocate (red_space_hess(ntrial, ntrial))
        call dgemm("T", "N", ntrial, ntrial, n_param, 1.d0, red_space_basis, &
                   n_param, h_basis, n_param, 0.d0, red_space_hess, ntrial)

        ! allocate space for reduced space solution
        allocate (red_space_solution(size(red_space_basis, 2)))

        ! loop over iterations
        do iter = 1, settings%n_iter
            ! solve reduced space problem
            call symm_mat_min_eig(red_space_hess, eigval, red_space_solution)

            ! get full space solution
            call dgemv("N", size(red_space_basis, 1), size(red_space_basis, 2), 1.d0, &
                       red_space_basis, size(red_space_basis, 1), red_space_solution, &
                       1, 0.d0, solution, 1)

            ! calculate Hessian linear transformation of solution
            call dgemv("N", n_param, size(h_basis, 2), 1.d0, h_basis, n_param, &
                       red_space_solution, 1, 0.d0, h_solution, 1)

            ! calculate residual
            residual = h_solution - eigval*solution

            ! check convergence
            if (dnrm2(n_param, residual, 1) < settings%conv_tol) exit

            ! precondition residual
            precond = full_h_diag - eigval
            where (abs(precond) < 1e-10)
                precond = 1e-10
            end where
            residual = residual/precond

            ! orthonormalize to current orbital space to get new basis vector
            basis_vec = gram_schmidt(residual, red_space_basis)

            ! add new trial vector to orbital space
            call add_column(red_space_basis, basis_vec)

            ! add linear transformation of new basis vector
            call add_column(h_basis, full_hess_x(basis_vec))

            ! construct new reduced space Hessian
            allocate (red_space_hess_vec(size(red_space_basis, 2)))
            call dgemv("T", n_param, size(red_space_basis, 2), 1.d0, &
                       red_space_basis, n_param, &
                       h_basis(:, size(red_space_basis, 2)), 1, 0.d0, &
                       red_space_hess_vec, 1)
            call extend_symm_matrix(red_space_hess, red_space_hess_vec)
            deallocate (red_space_hess_vec)

            ! reallocate reduced space solution
            deallocate (red_space_solution)
            allocate (red_space_solution(size(red_space_basis, 2)))

        end do

        ! increment total number of Hessian linear transformations
        tot_hess_x = tot_hess_x + size(red_space_basis, 2)

        ! deallocate quantities from Davidson iterations
        deallocate (red_space_solution)
        deallocate (red_space_hess)
        deallocate (h_basis)
        deallocate (red_space_basis)

        ! determine if saddle point
        stable = eigval > -1.d-3
        if (stable) then
            kappa = 0.d0
        else
            kappa = solution
            if (settings%verbose > 1) then
                write (stderr, '(A, F0.2)') " Solution not stable. Lowest "// &
                    "eigenvalue: ", eigval
            end if
        end if

        ! flush output
        flush (stdout)
        flush (stderr)

    contains

        function full_hess_x(x)
            !
            ! this function calculates the full Hessian linear transformation
            !
            real(rp), intent(in) :: x(:)

            real(rp) :: full_hess_x(size(x))

            full_hess_x = 2.d0*hess_x_funptr(x)

        end function full_hess_x

    end subroutine stability_check

    subroutine newton_step(aug_hess, grad_norm, red_space_basis, solution, &
                           red_space_solution)
        !
        ! this subroutine performs a Newton step by solving the Newton equations in
        ! reduced space without a level shift
        !
        real(rp), intent(in) :: aug_hess(:, :), grad_norm, red_space_basis(:, :)
        real(rp), intent(out) :: solution(:), red_space_solution(:)

        integer(ip) :: nred, lwork, info
        real(rp) :: red_hess(size(red_space_basis, 2), size(red_space_basis, 2)), &
                    ipiv(size(red_space_basis, 2))
        real(rp), allocatable :: work(:)
        external :: dsysv, dgemv

        ! reduced space size
        nred = size(red_space_basis, 2)

        ! reduced space Hessian
        red_hess = aug_hess(2:, 2:)

        ! set gradient
        red_space_solution = 0.d0
        red_space_solution(1) = -grad_norm

        ! query optimal workspace size
        lwork = -1
        allocate (work(1))
        call dsysv("U", nred, 1, red_hess, nred, ipiv, red_space_solution, nred, work, &
                   lwork, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))

        ! solve linear system
        call dsysv("U", nred, 1, red_hess, nred, ipiv, red_space_solution, nred, work, &
                   lwork, info)

        ! deallocate work array
        deallocate (work)

        ! check for errors
        if (info /= 0) then
            write (stderr, '(A, I0)') " Error in DSYSV, info = ", info
            error stop "Linear solver failed"
        end if

        ! get solution in full space
        call dgemv("N", size(red_space_basis, 1), size(red_space_basis, 2), 1.d0, &
                   red_space_basis, size(red_space_basis, 1), red_space_solution, 1, &
                   0.d0, solution, 1)

    end subroutine newton_step

    subroutine bisection(aug_hess, grad_norm, red_space_basis, trust_radius, solution, &
                         red_space_solution, mu, error_flag)
        !
        ! this subroutine performs bisection to find the parameter alpha that matches
        ! the desired trust radius
        !
        real(rp), intent(inout) :: aug_hess(:, :)
        real(rp), intent(in) :: grad_norm, red_space_basis(:, :), trust_radius
        real(rp), intent(out) :: solution(:), red_space_solution(:), mu
        character(:), intent(out), allocatable, optional :: error_flag

        real(rp) :: lower_alpha, middle_alpha, upper_alpha, lower_trust_dist, &
                    middle_trust_dist, upper_trust_dist
        real(rp), external :: dnrm2

        ! lower and upper bracket for alpha
        lower_alpha = 1.d-4
        upper_alpha = 1.d6

        ! solve reduced space problem with scaled gradient
        call get_ah_lowest_eigenvec(lower_alpha)
        lower_trust_dist = dnrm2(size(solution), solution, 1) - trust_radius
        call get_ah_lowest_eigenvec(upper_alpha)
        upper_trust_dist = dnrm2(size(solution), solution, 1) - trust_radius

        ! check if trust region is within bracketing range
        if ((lower_trust_dist*upper_trust_dist) > 0.d0) then
            error_flag = raise_error("Target trust region outside of bracketing "// &
                                     "range.", present(error_flag))
            solution = 0.d0
            red_space_solution = 0.d0
            mu = 0.d0
            return
        end if

        ! get middle alpha
        middle_alpha = sqrt(upper_alpha*lower_alpha)
        call get_ah_lowest_eigenvec(middle_alpha)
        middle_trust_dist = dnrm2(size(solution), solution, 1) - trust_radius

        ! perform bisection to find root
        do while (upper_alpha - lower_alpha > 1.d-10)
            ! targeted trust radius is in upper bracket
            if (lower_trust_dist*middle_trust_dist > 0.d0) then
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
            middle_trust_dist = dnrm2(size(solution), solution, 1) - trust_radius
        end do

    contains

        subroutine get_ah_lowest_eigenvec(alpha)
            !
            ! this subroutine returns the lowest eigenvector for an augmented Hessian
            !
            real(rp), intent(in) :: alpha

            real(rp) :: eigvec(size(aug_hess, 1))
            external :: dgemv

            ! finish construction of augmented Hessian
            aug_hess(1, 2) = alpha*grad_norm
            aug_hess(2, 1) = alpha*grad_norm

            ! perform eigendecomposition and get lowest eigenvalue and corresponding
            ! eigenvector
            call symm_mat_min_eig(aug_hess, mu, eigvec)

            ! scale eigenvector such that first element is equal to one and divide by
            ! alpha to get solution in reduced space
            red_space_solution = eigvec(2:)/eigvec(1)/alpha

            ! get solution in full space
            call dgemv("N", size(red_space_basis, 1), size(red_space_basis, 2), 1.d0, &
                       red_space_basis, size(red_space_basis, 1), red_space_solution, &
                       1, 0.d0, solution, 1)

        end subroutine get_ah_lowest_eigenvec

    end subroutine bisection

    function bracket(obj_func, kappa, lower, upper) result(n_kappa)
        !
        ! this function brackets a minimum (algorithm from numerical recipes)
        !
        procedure(obj_func_type), intent(in), pointer :: obj_func
        real(rp), intent(in) :: kappa(:), lower, upper

        real(rp) :: n_kappa, f_upper, f_lower, n_a, n_b, n_c, n_u, f_a, f_b, f_c, f_u, &
                    n_u_lim, tmp1, tmp2, val, denom
        real(rp), parameter :: golden_ratio = (1.d0 + sqrt(5.d0))/2.d0, &
                               grow_limit = 110.d0

        ! evaluate function at upper and lower bounds
        f_lower = obj_func(lower*kappa)
        f_upper = obj_func(upper*kappa)

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
        f_c = obj_func(n_c*kappa)

        ! continue looping until function at middle point is lower than at brackets
        do while (f_c < f_b)
            ! compute value u by parabolic extrapolation
            tmp1 = (n_b - n_a)*(f_b - f_c)
            tmp2 = (n_b - n_c)*(f_b - f_a)
            val = tmp2 - tmp1
            denom = 2.d0*sign(max(abs(val), 1.d-21), val)
            n_u = n_b - ((n_b - n_c)*tmp2 - (n_b - n_a)*tmp1)/denom

            ! maximum growth in parabolic fit
            n_u_lim = n_b + grow_limit*(n_c - n_b)

            ! check if u is between n_b and n_c
            if ((n_u - n_c)*(n_b - n_u) > 0.0) then
                ! evaluate function at n_u
                f_u = obj_func(n_u*kappa)

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
                f_u = obj_func(n_u*kappa)
                ! limit parabolic fit to its maximum allowed value
            else if ((n_u - n_u_lim)*(n_u_lim - n_c) >= 0.d0) then
                n_u = n_u_lim
                f_u = obj_func(n_u*kappa)
                ! parabolic fit is between n_c and its allowed limit
            else if ((n_u - n_u_lim)*(n_c - n_u) > 0.d0) then
                ! evaluate function at n_u
                f_u = obj_func(n_u*kappa)

                if (f_u < f_c) then
                    n_b = n_c
                    n_c = n_u
                    n_u = n_c + golden_ratio*(n_c - n_b)
                    f_b = f_c
                    f_c = f_u
                    f_u = obj_func(n_u*kappa)
                end if
                ! reject parabolic fit and use default step
            else
                n_u = n_c + golden_ratio*(n_c - n_b)
                f_u = obj_func(n_u*kappa)
            end if
            ! remove oldest point
            n_a = n_b
            n_b = n_c
            n_c = n_u
            f_a = f_b
            f_b = f_c
            f_c = f_u

        end do

        ! check if miniumum is bracketed
        if (.not. (((f_b < f_c .and. f_b <= f_a) .or. (f_b < f_a .and. f_b <= f_c)) &
                   .and. ((n_a < n_b .and. n_b < n_c) .or. &
                          (n_c < n_b .and. n_b < n_a)))) &
            error stop "Line search did not find minimum"

        ! set new multiplier
        n_kappa = n_b

    end function bracket

    subroutine extend_symm_matrix(matrix, vector)
        !
        ! this subroutine extends a symmetric matrix
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: vector(:)

        real(rp), allocatable :: temp(:, :)

        ! allocate temporary array
        allocate (temp(size(matrix, 1), size(matrix, 2)))

        ! copy data
        temp = matrix

        ! reallocate matrix
        deallocate (matrix)
        allocate (matrix(size(temp, 1) + 1, size(temp, 2) + 1))

        ! copy data back
        matrix(:size(temp, 1), :size(temp, 2)) = temp

        ! add new row and column
        matrix(:, size(temp, 2) + 1) = vector
        matrix(size(temp, 1) + 1, :) = vector

        ! deallocate temporary array
        deallocate (temp)

    end subroutine extend_symm_matrix

    subroutine add_column(matrix, new_col)
        !
        ! this subroutine adds a new column to an array
        !
        real(rp), allocatable, intent(inout) :: matrix(:, :)
        real(rp), intent(in) :: new_col(:)

        real(rp), allocatable :: temp(:, :)

        ! allocate temporary array
        allocate (temp(size(matrix, 1), size(matrix, 2)))

        ! copy data
        temp = matrix

        ! reallocate matrix
        deallocate (matrix)
        allocate (matrix(size(temp, 1), size(temp, 2) + 1))

        ! copy data back
        matrix(:, :size(temp, 2)) = temp

        ! add new column
        matrix(:, size(temp, 2) + 1) = new_col

        ! deallocate temporary array
        deallocate (temp)

    end subroutine add_column

    subroutine symm_mat_min_eig(symm_matrix, lowest_eigval, lowest_eigvec)
        !
        ! this function returns the lowest eigenvalue and corresponding eigenvector of
        ! a symmetric matrix
        !
        real(rp), intent(in) :: symm_matrix(:, :)
        real(rp), intent(out) :: lowest_eigval, lowest_eigvec(:)

        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:)
        real(rp) :: eigvals(size(symm_matrix, 1)), &
                    eigvecs(size(symm_matrix, 1), size(symm_matrix, 2))
        external :: dsyev

        ! size of matrix
        n = size(symm_matrix, 1)

        ! copy matrix
        eigvecs = symm_matrix

        ! query optimal workspace size
        lwork = -1
        allocate (work(1))
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))

        ! perform eigendecomposition
        call dsyev("V", "U", n, eigvecs, n, eigvals, work, lwork, info)

        ! deallocate work array
        deallocate (work)

        ! check for successful execution
        if (info /= 0) then
            write (stderr, '(A, I0)') " Error in DSYEV, info = ", info
            error stop "Eigendecomposition failed"
        end if

        ! get lowest eigenvalue and corresponding eigenvector
        lowest_eigval = eigvals(1)
        lowest_eigvec = eigvecs(:, 1)

    end subroutine symm_mat_min_eig

    real(rp) function min_eigval(matrix)
        !
        ! this function calculates the lowest eigenvalue of a symmetric matrix
        !
        real(rp), intent(in) :: matrix(:, :)

        real(rp) :: eigvals(size(matrix, 1)), temp(size(matrix, 1), size(matrix, 2))
        integer(ip) :: n, lwork, info
        real(rp), allocatable :: work(:)
        external :: dsyev

        ! size of matrix
        n = size(matrix, 1)

        ! copy matrix to avoid modification of original matrix
        temp = matrix

        ! query optimal workspace size
        lwork = -1
        allocate (work(1))
        call dsyev("N", "U", n, temp, n, eigvals, work, lwork, info)
        lwork = int(work(1))
        deallocate (work)
        allocate (work(lwork))

        ! compute eigenvalues
        call dsyev("N", "U", n, temp, n, eigvals, work, lwork, info)

        ! deallocate work array
        deallocate (work)

        ! check for successful execution
        if (info /= 0) then
            write (stderr, '(A, I0)') " Error in DSYEV, info = ", info
            error stop "Eigendecomposition failed"
        end if

        ! get lowest eigenvalue
        min_eigval = eigvals(1)

    end function min_eigval

    subroutine init_rng(seed)
        !
        ! this subroutine initializes the random number generator
        !
        integer(ip), intent(in) :: seed

        integer, allocatable :: seed_arr(:)
        integer :: n

        call random_seed(size=n)
        allocate (seed_arr(n))
        seed_arr = int(seed, kind=4)
        call random_seed(put=seed_arr)
        deallocate (seed_arr)

    end subroutine init_rng

    function generate_trial_vectors(n_random_trial_vectors, grad, grad_norm, h_diag) &
        result(red_space_basis)
        !
        ! this function generates trial vectors
        !
        integer(ip), intent(in) :: n_random_trial_vectors
        real(rp), intent(in) :: grad(:), grad_norm, h_diag(:)

        real(rp), allocatable :: red_space_basis(:, :)

        integer(ip) :: min_idx, n_vectors, i
        real(rp) :: trial(size(grad))
        real(rp), external :: dnrm2

        min_idx = minloc(h_diag, dim=1)

        if (h_diag(min_idx) < 0.d0) then
            n_vectors = 2
            allocate (red_space_basis(size(grad), n_vectors + n_random_trial_vectors))
            red_space_basis(:, 1) = grad/grad_norm
            trial = 0.d0
            trial(min_idx) = 1.d0
            red_space_basis(:, 2) = gram_schmidt(trial, &
                                                 reshape(red_space_basis(:, 1), &
                                                         [size(red_space_basis, 1), 1]))
        else
            n_vectors = 1
            allocate (red_space_basis(size(grad), n_vectors + n_random_trial_vectors))
            red_space_basis(:, 1) = grad/grad_norm
        end if

        do i = 1, n_random_trial_vectors
            call random_number(trial)
            trial = 2*trial - 1
            do while (dnrm2(size(grad), trial, 1) < 1.d-3)
                call random_number(trial)
                trial = 2*trial - 1
            end do
            red_space_basis(:, n_vectors + i) = &
                gram_schmidt(trial, red_space_basis(:, :n_vectors + i - 1))
        end do

    end function generate_trial_vectors

    function gram_schmidt(vector, space, error_flag) result(orth_vector)
        !
        ! this function orthonormalizes a vector with respect to a vector space
        !
        real(rp), intent(in) :: vector(:), space(:, :)
        character(:), intent(out), allocatable, optional :: error_flag

        real(rp) :: orth_vector(size(vector)), orth(size(space, 2))

        integer(ip) :: i
        real(rp), external :: ddot, dnrm2

        if (dnrm2(size(vector), vector, 1) < 1.d-12) then
            error_flag = raise_error("Vector passed to Gram-Schmidt procedure is "// &
                                     "numerically zero.", present(error_flag))
            orth_vector = 0.d0
            return
        else if (size(space, 2) > size(space, 1) - 1) then
            error_flag = raise_error("Number of vectors in Gram-Schmidt procedure "// &
                                     "larger than dimension of vector space.", &
                                     present(error_flag))
            orth_vector = 0.d0
            return
        end if

        orth_vector = vector
        do while (.true.)
            do i = 1, size(space, 2)
                orth_vector = orth_vector - ddot(size(vector), orth_vector, 1, &
                                                 space(:, i), 1)*space(:, i)
            end do
            orth_vector = orth_vector/dnrm2(size(vector), orth_vector, 1)

            call dgemv("T", size(vector), size(space, 2), 1.d0, space, size(vector), &
                       orth_vector, 1, 0.d0, orth, 1)
            if (maxval(abs(orth)) < 1.d-14) exit
        end do

    end function gram_schmidt

    subroutine init_solver_settings(self, stability, line_search, conv_tol, &
                                    n_random_trial_vectors, start_trust_radius, &
                                    n_macro, n_micro, global_red_factor, &
                                    local_red_factor, verbose, seed)
        !
        ! this subroutine sets the optional settings to their default values
        !
        class(solver_settings_type), intent(inout) :: self
        logical, intent(in), optional :: stability, line_search
        real(rp), intent(in), optional :: conv_tol, start_trust_radius, &
                                          global_red_factor, local_red_factor
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_macro, n_micro, &
                                             verbose, seed

        self%stability = set_default(stability, solver_stability_default)
        self%line_search = set_default(line_search, solver_line_search_default)
        self%conv_tol = set_default(conv_tol, solver_conv_tol_default)
        self%n_random_trial_vectors = set_default(n_random_trial_vectors, &
                                                  solver_n_random_trial_vectors_default)
        self%start_trust_radius = set_default(start_trust_radius, &
                                              solver_start_trust_radius_default)
        self%n_macro = set_default(n_macro, solver_n_macro_default)
        self%n_micro = set_default(n_micro, solver_n_micro_default)
        self%global_red_factor = set_default(global_red_factor, &
                                             solver_global_red_factor_default)
        self%local_red_factor = set_default(local_red_factor, &
                                            solver_local_red_factor_default)
        self%verbose = set_default(verbose, solver_verbose_default)
        self%seed = set_default(seed, solver_seed_default)

    end subroutine init_solver_settings

    subroutine init_stability_settings(self, conv_tol, n_random_trial_vectors, n_iter, &
                                       verbose)
        !
        ! this subroutine sets the optional settings to their default values
        !
        class(stability_settings_type), intent(inout) :: self
        real(rp), intent(in), optional :: conv_tol
        integer(ip), intent(in), optional :: n_random_trial_vectors, n_iter, verbose

        self%conv_tol = set_default(conv_tol, stability_conv_tol_default)
        self%n_random_trial_vectors = set_default(n_random_trial_vectors, &
                                               stability_n_random_trial_vectors_default)
        self%n_iter = set_default(n_iter, stability_n_iter_default)
        self%verbose = set_default(verbose, stability_verbose_default)

    end subroutine init_stability_settings

    function set_default_real(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for reals
        !
        real(rp), intent(in)  :: default_value
        real(rp), intent(in), optional :: optional_value

        real(rp) :: variable

        if (present(optional_value)) then
            variable = optional_value
        else
            variable = default_value
        end if

    end function set_default_real

    function set_default_integer(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for integers
        !
        integer(ip), intent(in)  :: default_value
        integer(ip), intent(in), optional :: optional_value

        integer(ip) :: variable

        if (present(optional_value)) then
            variable = optional_value
        else
            variable = default_value
        end if

    end function set_default_integer

    function set_default_logical(optional_value, default_value) result(variable)
        !
        ! this function sets a default value for logicals
        !
        logical, intent(in)  :: default_value
        logical, intent(in), optional :: optional_value

        logical :: variable

        if (present(optional_value)) then
            variable = optional_value
        else
            variable = default_value
        end if

    end function set_default_logical

    function raise_error(error_msg, return_error)
        !
        ! this function raises an error
        !
        character(*), intent(in) :: error_msg
        logical, intent(in) :: return_error
        character(:), allocatable :: raise_error

        raise_error = error_msg
        if (.not. return_error) error stop error_msg

    end function raise_error

end module opentrustregion
