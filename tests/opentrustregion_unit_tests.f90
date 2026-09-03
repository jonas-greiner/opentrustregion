! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module opentrustregion_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! parameters for 6D Hartmann function
    real(rp), parameter :: alpha(4) = [1.0_rp, 1.2_rp, 3.0_rp, 3.2_rp]
    real(rp), parameter :: A(4, 6) = reshape([10.0_rp, 0.05_rp, 3.0_rp, 17.0_rp, &
                                              3.0_rp, 10.0_rp, 3.5_rp, 8.0_rp, &
                                              17.0_rp, 17.0_rp, 1.7_rp, 0.05_rp, &
                                              3.5_rp, 0.1_rp, 10.0_rp, 10.0_rp, &
                                              1.7_rp, 8.0_rp, 17.0_rp, 0.1_rp, &
                                              8.0_rp, 14.0_rp, 8.0_rp, 14.0_rp], [4, 6])
    real(rp), parameter :: P(4, 6) = reshape([0.1312_rp, 0.2329_rp, 0.2348_rp, &
                                              0.4047_rp, 0.1696_rp, 0.4135_rp, &
                                              0.1451_rp, 0.8828_rp, 0.5569_rp, &
                                              0.8307_rp, 0.3522_rp, 0.8732_rp, &
                                              0.0124_rp, 0.3736_rp, 0.2883_rp, &
                                              0.5743_rp, 0.8283_rp, 0.1004_rp, &
                                              0.3047_rp, 0.1091_rp, 0.5886_rp, &
                                              0.9991_rp, 0.6650_rp, 0.0381_rp], [4, 6])

    ! stationary points of 6D Hartmann function
    real(rp), parameter :: minimum1(6) = [0.20168951_rp, 0.15001069_rp, 0.47687398_rp, &
                                          0.27533243_rp, 0.31165162_rp, 0.65730053_rp]
    real(rp), parameter :: minimum2(6) = [0.40465313_rp, 0.88244493_rp, 0.84610160_rp, &
                                          0.57398969_rp, 0.13892673_rp, 0.03849589_rp]
    real(rp), parameter :: saddle_point(6) = [0.35278250_rp, 0.59374767_rp, &
                                              0.47631257_rp, 0.40058250_rp, &
                                              0.31111531_rp, 0.32397158_rp]

    ! define global current variable and Hessian for 6D Hartmann model and larger mock 
    ! Hessian so that these can be accessed by procedure pointers
    integer(ip), parameter :: n_param_mock_hess = 12
    real(rp) :: curr_vars(6), hess(6, 6), &
                mock_hess_mat(n_param_mock_hess, n_param_mock_hess)

    ! global log message
    character(:), allocatable :: log_message

contains

    ! 6D Hartmann function definition

    function hartmann6d_func(vars) result(f)
        !
        ! this function defines the Hartmann 6D function
        !
        real(rp), intent(in) :: vars(:)
        real(rp) :: f, exp_term(4)
        integer(ip) :: i

        do i = 1, 4
            exp_term(i) = exp(-sum(A(i, :)*(vars - P(i, :))**2))
        end do

        f = -sum(alpha*exp_term)

    end function hartmann6d_func

    subroutine hartmann6d_gradient(vars, grad)
        !
        ! this subroutine defines the Hartmann 6D function's gradient
        !
        real(rp), intent(in) :: vars(:)
        real(rp), intent(out) :: grad(:)
        real(rp) :: exp_term(4)
        integer(ip) :: i, j

        do i = 1, 4
            exp_term(i) = exp(-sum(A(i, :)*(vars - P(i, :))**2))
        end do

        do j = 1, size(vars)
            grad(j) = sum(2.0_rp*alpha*A(:, j)*(vars(j) - P(:, j))*exp_term)
        end do

    end subroutine hartmann6d_gradient

    subroutine hartmann6d_hessian(vars)
        !
        ! this subroutine defines the Hartmann 6D function's Hessian
        !
        real(rp), intent(in) :: vars(:)
        real(rp) :: exp_term(4)
        integer(ip) :: i, j

        do i = 1, 4
            exp_term(i) = exp(-sum(A(i, :)*(vars - P(i, :))**2))
        end do

        do i = 1, size(vars)
            hess(i, i) = 2.0_rp*sum(alpha*A(:, i)*exp_term* &
                                    (1.0_rp - 2.0_rp*A(:, i)*(vars(i) - P(:, i))**2))
            do j = 1, i - 1
                hess(i, j) = -4.0_rp*sum(alpha*A(:, i)*A(:, j)*(vars(i) - P(:, i))* &
                                         (vars(j) - P(:, j))*exp_term)
                hess(j, i) = hess(i, j)
            end do
        end do

    end subroutine hartmann6d_hessian

    function hartmann6d_hess_x(x)
        !
        ! this function defines the Hessian linear transformation operation for the
        ! Hartmann 6D function
        !
        real(rp), intent(in) :: x(:)
        real(rp) :: hartmann6d_hess_x(size(x))

        hartmann6d_hess_x = matmul(hess, x)

    end function hartmann6d_hess_x

    subroutine hess_x_fun(x, hess_x, error)
        !
        ! this function describes the Hessian linear transformation operation for the
        ! Hartmann 6D function
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        hess_x = hartmann6d_hess_x(x)

    end subroutine hess_x_fun

    subroutine mock_hess_x_fun(x, hess_x, error)
        !
        ! this subroutine describes the Hessian linear transformation for a larger 
        ! mock Hessian
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        hess_x = matmul(mock_hess_mat, x)

    end subroutine mock_hess_x_fun

    function obj_func(delta_vars, error) result(func)
        !
        ! this function describes the objective function evaluation for the Hartmann
        ! 6D function
        !
        real(rp), intent(in), target :: delta_vars(:)
        integer(ip), intent(out) :: error
        real(rp) :: func

        ! initialize error flag
        error = 0

        func = hartmann6d_func(curr_vars + delta_vars)

    end function obj_func

    subroutine update_orbs(delta_vars, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function describes the orbital update equivalent for the Hartmann 6D
        ! function
        !
        use opentrustregion, only: hess_x_type

        real(rp), intent(in), target :: delta_vars(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: i

        ! initialize error flag
        error = 0

        ! update variables
        curr_vars = curr_vars + delta_vars

        ! evaluate function, calculate gradient and Hessian diagonal and define
        ! Hessian linear transformation
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

    end subroutine update_orbs

    subroutine mock_precond(residual, mu, precond_residual, error)
        !
        ! this subroutine is a test subroutine for the preconditioner subroutine
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        precond_residual = mu * residual

        error = 0

    end subroutine mock_precond

    subroutine mock_precond_pd(residual, precond_residual, error)
        !
        ! this subroutine is a test subroutine for the positive-definite preconditioner 
        ! subroutine
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        precond_residual = 3.0_rp * residual

        error = 0

    end subroutine mock_precond_pd

    subroutine mock_project(vector, error)
        !
        ! this subroutine is a test projection subroutine that projects onto the
        ! subspace where the first two components are equal
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error

        vector(1) = 0.5_rp * (vector(1) + vector(2))
        vector(2) = vector(1)

        error = 0

    end subroutine mock_project

    subroutine mock_approx_hess_x(x, hess_x, error)
        !
        ! this subroutine is a test subroutine for the approximate Hessian linear 
        ! transformation subroutine
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: i

        hess_x = [(hess(i, i) * x(i), i = 1, size(x))]

        error = 0

    end subroutine mock_approx_hess_x

    subroutine mock_init_trial_space(trial_space, error)
        !
        ! this subroutine is a test subroutine for the trial space initialization 
        ! subroutine
        !
        real(rp), intent(out), target :: trial_space(:, :)
        integer(ip), intent(out) :: error

        integer(ip) :: i

        trial_space = 0.0_rp
        do i = 1, size(trial_space, 2)
            trial_space(i, i) = 1.0_rp
        end do

        error = 0

    end subroutine mock_init_trial_space

    subroutine logger(message)
        !
        ! this subroutine is a mock logging subroutine
        !
        character(*), intent(in) :: message

        log_message = log_message//trim(message)

    end subroutine logger

    subroutine setup_settings(settings)
        !
        ! this subroutine sets up a settings object for tests
        !
        use opentrustregion, only: settings_type

        class(settings_type), intent(inout) :: settings
        integer(ip) :: error

        call settings%init(error)
        settings%verbose = 3
        settings%logger => logger
        log_message = ""

    end subroutine setup_settings

    logical(c_bool) function test_solver() bind(C)
        !
        ! this function tests the solver subroutine
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, &
                                   solver_settings_type, solver, &
                                   default_settings => default_solver_settings

        integer(ip), parameter :: n_param = 6
        real(rp), parameter :: var_thres = 1e-6_rp
        integer(ip) :: error
        real(rp), allocatable :: final_grad(:)
        procedure(update_orbs_type), pointer :: update_orbs_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr
        type(solver_settings_type) :: settings

        ! assume tests pass
        test_solver = .true.

        ! start in quadratic region near minimum
        curr_vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! initialize settings
        call settings%init(error)

        ! allocate space for the final gradient
        allocate(final_grad(n_param))

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            default_settings%conv_tol) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > var_thres)) then
            write (stderr, *) "test_solver failed: Solver did not find correct minimum."
            test_solver = .false.
        end if

        ! start near saddle point
        curr_vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            default_settings%conv_tol) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > var_thres) .and. &
            any(abs(curr_vars - minimum2) > var_thres)) then
            write (stderr, *) "test_solver failed: Solver did not find correct minimum."
            test_solver = .false.
        end if

        ! start at saddle point
        curr_vars = saddle_point
        update_orbs_funptr => update_orbs
        obj_func_funptr => obj_func

        ! run solver, check if error has occured and check whether gradient is zero and 
        ! agrees with correct minimum
        call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
        if (error /= 0) then
            write (stderr, *) "test_solver failed: Produced error."
            test_solver = .false.
        end if
        call hartmann6d_gradient(curr_vars, final_grad)
        if (norm2(final_grad)/sqrt(real(n_param, kind=rp)) > &
            default_settings%conv_tol) then
            write (stderr, *) "test_solver failed: Solver did not find stationary "// &
                "point."
            test_solver = .false.
        end if
        if (any(abs(curr_vars - minimum1) > var_thres) .and. &
            any(abs(curr_vars - minimum2) > var_thres)) then
            write (stderr, *) "test_solver failed: Solver did not find minimum."
            test_solver = .false.
        end if

        ! deallocate space for the gradient
        deallocate(final_grad)

    end function test_solver

    logical(c_bool) function test_block_davidson() bind(C)
        !
        ! this function tests the Davidson procedure for both a single tracked
        ! eigenpair and a larger block
        !
        use opentrustregion, only: hess_x_type, stability_settings_type, block_davidson

        integer(ip), parameter :: n_param = n_param_mock_hess, n_block = 3, &
                                  lwork = 3 * n_param
        real(rp), parameter :: res_tol = 1e-8_rp, asymm_perturb = 0.1_rp
        real(rp) :: ref_mat(n_param, n_param), ref_eigvals(n_param), work(lwork), &
                    eigvals(n_block), eigvecs(n_param, n_block), &
                    ref_eigvals_im(n_param), dummy_vecs(1, 1), swap
        real(rp), allocatable :: trial_space(:, :)
        integer(ip) :: error, i, j, k, min_idx, info
        logical :: converged
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(stability_settings_type) :: settings
        external :: dsyev, dgeev

        ! assume tests pass
        test_block_davidson = .true.

        ! construct symmetric matrix with well separated eigenvalues
        do i = 1, n_param
            do j = 1, n_param
                mock_hess_mat(i, j) = 1.0_rp / real(i + j, kind=rp)
            end do
            mock_hess_mat(i, i) = mock_hess_mat(i, i) + real(i, kind=rp)
        end do

        ! obtain reference eigenvalues
        ref_mat = mock_hess_mat
        call dsyev("N", "U", n_param, ref_mat, n_param, ref_eigvals, work, lwork, info)
        if (info /= 0) then
            write (stderr, *) "test_block_davidson failed: Reference "// &
                "eigendecomposition failed."
            test_block_davidson = .false.
            return
        end if

        ! setup settings object
        call setup_settings(settings)

        ! set Hessian linear transformation function pointer
        hess_x_funptr => mock_hess_x_fun

        ! check block of eigenpairs
        allocate(trial_space(n_param, n_block))
        trial_space = 0.0_rp
        do i = 1, n_block
            trial_space(i, i) = 1.0_rp
        end do
        call block_davidson(hess_x_funptr, [(mock_hess_mat(i, i), i=1, n_param)], &
                            trial_space, .true., n_block, 200_ip, .false., eigvals, &
                            eigvecs, converged, settings, error, res_tol=res_tol)
        if (error /= 0) then
            write (stderr, *) "test_block_davidson failed: Produced error for "// &
                "block of eigenpairs."
            test_block_davidson = .false.
            return
        end if
        if (.not. converged) then
            write (stderr, *) "test_block_davidson failed: Block of eigenpairs did "// &
                "not converge."
            test_block_davidson = .false.
        end if
        do i = 1, n_block
            if (abs(eigvals(i) - ref_eigvals(i)) > tol) then
                write (stderr, *) "test_block_davidson failed: Block of eigenpairs "// &
                    "does not return the lowest eigenvalues."
                test_block_davidson = .false.
            end if
            if (norm2(matmul(mock_hess_mat, eigvecs(:, i)) - &
                      eigvals(i) * eigvecs(:, i)) > res_tol) then
                write (stderr, *) "test_block_davidson failed: Block of eigenpairs "// &
                    "does not return converged eigenvectors."
                test_block_davidson = .false.
            end if
            if (abs(norm2(eigvecs(:, i)) - 1.0_rp) > tol) then
                write (stderr, *) "test_block_davidson failed: Block of eigenpairs "// &
                    "does not return normalized eigenvectors."
                test_block_davidson = .false.
            end if
        end do
        deallocate(trial_space)

        ! check single eigenpair
        settings%conv_tol = res_tol
        allocate(trial_space(n_param, 1))
        trial_space = 0.0_rp
        trial_space(1, 1) = 1.0_rp
        call block_davidson(hess_x_funptr, [(mock_hess_mat(i, i), i=1, n_param)], &
                            trial_space, .true., 1_ip, 200_ip, .false., eigvals(:1), &
                            eigvecs(:, :1), converged, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_block_davidson failed: Produced error for "// &
                "single eigenpair."
            test_block_davidson = .false.
            return
        end if
        if (.not. converged) then
            write (stderr, *) "test_block_davidson failed: Single eigenpair did "// &
                "not converge."
            test_block_davidson = .false.
        end if
        if (abs(eigvals(1) - ref_eigvals(1)) > tol) then
            write (stderr, *) "test_block_davidson failed: Single eigenpair does "// &
                "not return the lowest eigenvalue."
            test_block_davidson = .false.
        end if
        deallocate(trial_space)

        ! construct non-symmetric transformation
        do i = 1, n_param
            do j = 1, n_param
                mock_hess_mat(i, j) = 1.0_rp/real(i + j, kind=rp)
            end do
            mock_hess_mat(i, i) = mock_hess_mat(i, i) + real(i, kind=rp)
        end do
        do i = 1, n_param
            do j = i + 1, n_param
                mock_hess_mat(i, j) = mock_hess_mat(i, j) + &
                                        asymm_perturb/real(i + j, kind=rp)
                mock_hess_mat(j, i) = mock_hess_mat(j, i) - &
                                        asymm_perturb/real(i + j, kind=rp)
            end do
        end do

        ! obtain reference eigenvalues which are then sorted by real part with a 
        ! selection sort since DGEEV returns them in no particular order
        ref_mat = mock_hess_mat
        call dgeev("N", "N", n_param, ref_mat, n_param, ref_eigvals, ref_eigvals_im, &
                   dummy_vecs, 1_ip, dummy_vecs, 1_ip, work, lwork, info)
        do i = 1, n_param - 1
            min_idx = i
            do k = i + 1, n_param
                if (ref_eigvals(k) < ref_eigvals(min_idx)) min_idx = k
            end do
            swap = ref_eigvals(i)
            ref_eigvals(i) = ref_eigvals(min_idx)
            ref_eigvals(min_idx) = swap
        end do

        ! check block of eigenpairs
        allocate(trial_space(n_param, n_block))
        trial_space = 0.0_rp
        do i = 1, n_block
            trial_space(i, i) = 1.0_rp
        end do
        call block_davidson(hess_x_funptr, [(mock_hess_mat(i, i), i=1, n_param)], &
                            trial_space, .false., n_block, 200_ip, .false., eigvals, &
                            eigvecs, converged, settings, error, res_tol=res_tol)
        if (error /= 0) then
            write (stderr, *) "test_block_davidson failed: Produced error for "// &
                "non-symmetric block of eigenpairs."
            test_block_davidson = .false.
            return
        end if
        if (.not. converged) then
            write (stderr, *) "test_block_davidson failed: Non-symmetric block of "// &
                "eigenpairs did not converge."
            test_block_davidson = .false.
        end if
        do i = 1, n_block
            if (abs(eigvals(i) - ref_eigvals(i)) > tol) then
                write (stderr, *) "test_block_davidson failed: Non-symmetric block "// &
                    "of eigenpairs does not return the eigenvalues with the lowest "// &
                    "real parts."
                test_block_davidson = .false.
            end if
            if (norm2(matmul(mock_hess_mat, eigvecs(:, i)) - &
                      eigvals(i) * eigvecs(:, i)) > res_tol) then
                write (stderr, *) "test_block_davidson failed: Non-symmetric block "// &
                    "of eigenpairs does not return converged eigenvectors."
                test_block_davidson = .false.
            end if
        end do
        deallocate(trial_space)

        ! test behavior when stopping at instability by building the well-separated 
        ! symmetric matrix and forcing its first diagonal entry to be clearly unstable, 
        ! letting a single iteration with a tight residual tolerance distinguish the 
        ! two settings unambiguously
        do i = 1, n_param
            do j = 1, n_param
                mock_hess_mat(i, j) = 1.0_rp/real(i + j, kind=rp)
            end do
            mock_hess_mat(i, i) = mock_hess_mat(i, i) + real(i, kind=rp)
        end do
        mock_hess_mat(1, 1) = -1.0_rp

        call setup_settings(settings)
        settings%conv_tol = 1e-8_rp
        allocate(trial_space(n_param, 1))
        trial_space = 0.0_rp
        trial_space(1, 1) = 1.0_rp

        ! when not stopping at instability, a single iteration must not report
        ! convergence, since the residual is far from the tight tolerance and the
        ! eigenvalue is not consulted on its own
        call block_davidson(hess_x_funptr, [(mock_hess_mat(i, i), i=1, n_param)], &
                            trial_space, .true., 1_ip, 1_ip, .false., eigvals(:1), &
                            eigvecs(:, :1), converged, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_block_davidson failed: Produced error when not "// &
                "stopping at instability."
            test_block_davidson = .false.
            return
        end if
        if (converged) then
            write (stderr, *) "test_block_davidson failed: Converged after a "// &
                "single iteration with an unconverged residual although "// &
                "not stopping at instability."
            test_block_davidson = .false.
        end if
        deallocate(trial_space)

        ! when stopping at instability, the same single iteration must instead report 
        ! convergence immediately, since the (exactly known) eigenvalue is already 
        ! below the instability threshold
        settings%stop_on_instability = .true.
        allocate(trial_space(n_param, 1))
        trial_space = 0.0_rp
        trial_space(1, 1) = 1.0_rp
        call block_davidson(hess_x_funptr, [(mock_hess_mat(i, i), i=1, n_param)], &
                            trial_space, .true., 1_ip, 1_ip, .false., eigvals(:1), &
                            eigvecs(:, :1), converged, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_block_davidson failed: Produced error when "// &
                "stopping at instability."
            test_block_davidson = .false.
            return
        end if
        if (.not. converged) then
            write (stderr, *) "test_block_davidson failed: Did not converge on an "// &
                "established instability when stopping at instability."
            test_block_davidson = .false.
        end if
        if (abs(eigvals(1) + 1.0_rp) > tol) then
            write (stderr, *) "test_block_davidson failed: Incorrect Ritz value "// &
                "when stopping at instability."
            test_block_davidson = .false.
        end if
        deallocate(trial_space)

    end function test_block_davidson

    logical(c_bool) function test_get_stability_trial_space() bind(C)
        !
        ! this function tests the subroutine which generates the initial trial space
        ! for the stability check
        !
        use opentrustregion, only: stability_settings_type, get_stability_trial_space

        integer(ip), parameter :: n_param = 4

        type(stability_settings_type) :: settings
        real(rp) :: h_diag(n_param), overlap, hess_diag_copy(n_param), hv(n_param), &
                    lambda, ref_eigvals(2)
        real(rp), allocatable :: red_space_basis(:, :)
        integer(ip) :: error, i, j, k, ref_idx(2)

        ! assume tests pass
        test_get_stability_trial_space = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_trial_vectors = 2
        settings%n_random_trial_vectors = 1
        settings%project => mock_project

        ! construct Hessian diagonal whose two lowest elements are the first two, so 
        ! the leading block is built from the first two directions, which the projector 
        ! above maps onto the same direction
        h_diag = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]

        ! generate trial space and check size and orthonormality
        call get_stability_trial_space(h_diag, red_space_basis, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_stability_trial_space failed: Produced error."
            test_get_stability_trial_space = .false.
            return
        end if
        if (size(red_space_basis, 2) /= settings%n_trial_vectors + &
            settings%n_random_trial_vectors) then
            write (stderr, *) "test_get_stability_trial_space failed: Incorrect "// &
                "number of trial vectors when the projector introduces a linear "// &
                "dependency."
            test_get_stability_trial_space = .false.
            return
        end if
        do i = 1, size(red_space_basis, 2)
            do j = 1, size(red_space_basis, 2)
                overlap = dot_product(red_space_basis(:, i), red_space_basis(:, j))
                if (i == j) overlap = overlap - 1.0_rp
                if (abs(overlap) > tol) then
                    write (stderr, *) "test_get_stability_trial_space failed: "// &
                        "Returned trial space is not orthonormal."
                    test_get_stability_trial_space = .false.
                end if
            end do
        end do

        ! use an approximate Hessian linear transformation to construct leading block
        call setup_settings(settings)
        settings%n_trial_vectors = 2
        settings%n_random_trial_vectors = 1
        call hartmann6d_hessian(minimum1)

        ! the Hessian diagonal preconditioner is given the opposite ranking of the
        ! approximate Hessian's true diagonal, so a silent fallback to the unit-vector 
        ! leading block would pick the wrong two directions entirely, letting the 
        ! eigenvector check below also confirm the correct branch ran
        h_diag = [(-hess(i, i), i = 1, n_param)]
        settings%approx_hess_x => mock_approx_hess_x

        ! generate trial space and check size, whether leading vectors are the lowest 
        ! eigenvectors of the approximate Hessian and orthonormality
        call get_stability_trial_space(h_diag, red_space_basis, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_stability_trial_space failed: Produced "// &
                "error with approximate Hessian linear transformation."
            test_get_stability_trial_space = .false.
            return
        end if
        if (size(red_space_basis, 2) /= settings%n_trial_vectors + &
            settings%n_random_trial_vectors) then
            write (stderr, *) "test_get_stability_trial_space failed: Incorrect "// &
                "number of trial vectors with approximate Hessian linear "// &
                "transformation."
            test_get_stability_trial_space = .false.
            return
        end if
        hess_diag_copy = [(hess(i, i), i = 1, n_param)]
        do k = 1, 2
            ref_idx(k) = minloc(hess_diag_copy, dim=1)
            hess_diag_copy(ref_idx(k)) = huge(1.0_rp)
        end do
        ref_eigvals = [(hess(ref_idx(k), ref_idx(k)), k = 1, 2)]
        do i = 1, settings%n_trial_vectors
            hv = [(hess(k, k) * red_space_basis(k, i), k = 1, n_param)]
            lambda = dot_product(red_space_basis(:, i), hv)
            if (norm2(hv - lambda * red_space_basis(:, i)) > tol) then
                write (stderr, *) "test_get_stability_trial_space failed: Leading "// &
                    "vector is not an eigenvector of approximate Hessian linear "// &
                    "transformation."
                test_get_stability_trial_space = .false.
            end if
            if (all(abs(lambda - ref_eigvals) > tol)) then
                write (stderr, *) "test_get_stability_trial_space failed: Leading "// &
                    "vector does not correspond to one of the two lowest "// &
                    "eigenvalues of approximate Hessian linear transformation."
                test_get_stability_trial_space = .false.
            end if
        end do
        do i = 1, size(red_space_basis, 2)
            do j = 1, size(red_space_basis, 2)
                overlap = dot_product(red_space_basis(:, i), red_space_basis(:, j))
                if (i == j) overlap = overlap - 1.0_rp
                if (abs(overlap) > tol) then
                    write (stderr, *) "test_get_stability_trial_space failed: "// &
                        "Returned trial space with approximate Hessian linear "// &
                        "transformation is not orthonormal."
                    test_get_stability_trial_space = .false.
                end if
            end do
        end do

        ! a supplied trial space initialization callback fully replaces the leading
        ! block and short-circuits the routine, so the random fill count has to be
        ! ignored entirely rather than added to the requested number of trial vectors
        call setup_settings(settings)
        settings%n_trial_vectors = 3
        settings%n_random_trial_vectors = 2
        settings%init_trial_space => mock_init_trial_space

        ! generate trial space and check size, whether the leading block is exactly
        ! the callback output and orthonormality
        call get_stability_trial_space(h_diag, red_space_basis, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_get_stability_trial_space failed: Produced "// &
                "error with trial space initialization."
            test_get_stability_trial_space = .false.
            return
        end if
        if (size(red_space_basis, 2) /= settings%n_trial_vectors) then
            write (stderr, *) "test_get_stability_trial_space failed: Trial space "// &
                "initialization callback result was padded with random vectors "// &
                "instead of being used as is."
            test_get_stability_trial_space = .false.
            return
        end if
        do i = 1, size(red_space_basis, 2)
            do k = 1, n_param
                if (abs(red_space_basis(k, i) - merge(1.0_rp, 0.0_rp, k == i)) > tol) &
                    then
                    write (stderr, *) "test_get_stability_trial_space failed: "// &
                        "Returned trial space does not match the trial space "// &
                        "initialization callback output."
                    test_get_stability_trial_space = .false.
                end if
            end do
        end do
        do i = 1, size(red_space_basis, 2)
            do j = 1, size(red_space_basis, 2)
                overlap = dot_product(red_space_basis(:, i), red_space_basis(:, j))
                if (i == j) overlap = overlap - 1.0_rp
                if (abs(overlap) > tol) then
                    write (stderr, *) "test_get_stability_trial_space failed: "// &
                        "Returned trial space with trial space initialization is "// &
                        "not orthonormal."
                    test_get_stability_trial_space = .false.
                end if
            end do
        end do

    end function test_get_stability_trial_space

    logical(c_bool) function test_stability_check() bind(C)
        !
        ! this function tests the stability check subroutine
        !
        use opentrustregion, only: hess_x_type, stability_settings_type, stability_check
        use test_reference, only: n_trial_vectors

        real(rp) :: vars(6), h_diag(6), direction(6), min_eigval, eigval_vec(6)
        procedure(hess_x_type), pointer :: hess_x_funptr
        logical :: stable
        integer(ip) :: error, i
        type(stability_settings_type) :: settings

        ! assume tests pass
        test_stability_check = .true.

        ! start at minimum and determine Hessian diagonal and define Hessian linear
        ! transformation
        vars = minimum1
        call hartmann6d_hessian(vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

        ! initialize settings
        call settings%init(error)

        ! run stability, check if error has occured check and determine whether minimum 
        ! is stable and returned direction and eigenvalue are correct
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             direction, min_eigval)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error."
            test_stability_check = .false.
        end if
        if (.not. stable) then
            write (stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of minimum."
            test_stability_check = .false.
        end if
        if (min_eigval <= 0.0_rp) then
            write (stderr, *) "test_stability_check failed: Stability check "// &
                "does not return correct eigenvalue for minimum."
            test_stability_check = .false.
        end if
        call hess_x_funptr(direction, eigval_vec, error)
        if (dot_product(direction, eigval_vec) - min_eigval > tol) then
            write (stderr, *) "test_stability_check failed: Stability check does "// &
                "not return correct eigenvector direction for minimum."
            test_stability_check = .false.
        end if

        ! initialize settings
        call settings%init(error)

        ! set initial trial space
        settings%n_trial_vectors = n_trial_vectors
        settings%init_trial_space => mock_init_trial_space

        ! run stability, check if error has occured check and determine whether minimum 
        ! is stable and the returned direction vanishes
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             direction, min_eigval)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error with "// &
                "custom trial space initialization."
            test_stability_check = .false.
        end if
        if (.not. stable) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "custom trial space initialization incorrectly classifies "// &
                "stability of minimum."
            test_stability_check = .false.
        end if
        if (min_eigval <= 0.0_rp) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "custom trial space initialization does not return correct "// &
                "eigenvalue for minimum."
            test_stability_check = .false.
        end if
        call hess_x_funptr(direction, eigval_vec, error)
        if (dot_product(direction, eigval_vec) - min_eigval > tol) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "custom trial space initialization does not return correct "// &
                "eigenvector direction for minimum."
            test_stability_check = .false.
        end if

        ! initialize settings
        call settings%init(error)

        ! set approximate Hessian linear transformation
        settings%approx_hess_x => mock_approx_hess_x

        ! run stability, check if error has occured check and determine whether minimum 
        ! is stable and the returned direction vanishes
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             direction, min_eigval)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error with "// &
                "approximate Hessian linear transformation."
            test_stability_check = .false.
        end if
        if (.not. stable) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "approximate Hessian linear transformation incorrectly classifies "// &
                "stability of minimum."
            test_stability_check = .false.
        end if
        if (min_eigval <= 0.0_rp) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "approximate Hessian linear transformation does not return correct "// &
                "eigenvalue for minimum."
            test_stability_check = .false.
        end if
        call hess_x_funptr(direction, eigval_vec, error)
        if (dot_product(direction, eigval_vec) - min_eigval > tol) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "approximate Hessian linear transformation does not return correct "// &
                "eigenvector direction for minimum."
            test_stability_check = .false.
        end if

        ! initialize settings
        call settings%init(error)

        ! start at saddle point and determine Hessian diagonal and define linear
        ! linear transformation
        vars = saddle_point
        call hartmann6d_hessian(vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        hess_x_funptr => hess_x_fun

        ! run stability check, check if error has occured and determine whether saddle 
        ! point is unstable and the returned direction is correct
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             direction, min_eigval)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error."
            test_stability_check = .false.
        end if
        if (stable) then
            write (stderr, *) "test_stability_check failed: Stability check "// &
                "incorrectly classifies stability of saddle point."
            test_stability_check = .false.
        end if
        if (min_eigval >= 0.0_rp) then
            write (stderr, *) "test_stability_check failed: Stability check "// &
                "does not return correct eigenvalue for saddle point."
            test_stability_check = .false.
        end if
        call hess_x_funptr(direction, eigval_vec, error)
        if (dot_product(direction, eigval_vec) - min_eigval > tol) then
            write (stderr, *) "test_stability_check failed: Stability check does "// &
                "not return correct eigenvector direction for saddle point."
            test_stability_check = .false.
        end if

        ! initialize settings
        call settings%init(error)

        ! set initial trial space
        settings%n_trial_vectors = n_trial_vectors
        settings%init_trial_space => mock_init_trial_space

        ! run stability check, check if error has occured and determine whether saddle 
        ! point is unstable and the returned direction is correct
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             direction, min_eigval)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error with "// &
                "custom trial space initialization."
            test_stability_check = .false.
        end if
        if (stable) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "custom trial space initialization incorrectly classifies "// &
                "stability of saddle point."
            test_stability_check = .false.
        end if
        if (min_eigval >= 0.0_rp) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "custom trial space initialization does not return correct "// &
                "eigenvalue for saddle point."
            test_stability_check = .false.
        end if
        call hess_x_funptr(direction, eigval_vec, error)
        if (dot_product(direction, eigval_vec) - min_eigval > tol) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "custom trial space initialization does not return correct "// &
                "eigenvector direction for saddle point."
            test_stability_check = .false.
        end if

        ! initialize settings
        call settings%init(error)

        ! set approximate Hessian linear transformation
        settings%approx_hess_x => mock_approx_hess_x

        ! run stability check, check if error has occured and determine whether saddle 
        ! point is unstable and the returned direction is correct
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             direction, min_eigval)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error with "// &
                "approximate Hessian linear transformation."
            test_stability_check = .false.
        end if
        if (stable) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "approximate Hessian linear transformation incorrectly classifies "// &
                "stability of saddle point."
            test_stability_check = .false.
        end if
        if (min_eigval >= 0.0_rp) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "approximate Hessian linear transformation does not return correct "// &
                "eigenvalue for saddle point."
            test_stability_check = .false.
        end if
        call hess_x_funptr(direction, eigval_vec, error)
        if (dot_product(direction, eigval_vec) - min_eigval > tol) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "approximate Hessian linear transformation does not return correct "// &
                "eigenvector direction for saddle point."
            test_stability_check = .false.
        end if

        ! initialize settings
        call settings%init(error)

        ! request a purely random trial space
        settings%n_trial_vectors = 0
        settings%n_random_trial_vectors = 3

        ! run stability check, check if error has occured and determine whether saddle 
        ! point is unstable and the returned direction is correct
        call stability_check(h_diag, hess_x_funptr, stable, error, settings, &
                             direction, min_eigval)
        if (error /= 0) then
            write (stderr, *) "test_stability_check failed: Produced error with "// &
                "purely random trial space."
            test_stability_check = .false.
        end if
        if (stable) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "purely random trial space incorrectly classifies stability of "// &
                "saddle point."
            test_stability_check = .false.
        end if
        if (min_eigval >= 0.0_rp) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "purely random trial space does not return correct eigenvalue for "// &
                "saddle point."
            test_stability_check = .false.
        end if
        call hess_x_funptr(direction, eigval_vec, error)
        if (dot_product(direction, eigval_vec) - min_eigval > tol) then
            write (stderr, *) "test_stability_check failed: Stability check with "// &
                "purely random trial space does not return correct eigenvector "// &
                "direction for saddle point."
            test_stability_check = .false.
        end if

    end function test_stability_check

    logical(c_bool) function test_newton_step() bind(C)
        !
        ! this function tests the Newton step subroutine
        !
        use opentrustregion, only: newton_step

        integer(ip), parameter :: lwork = 12
        real(rp) :: red_space_basis(6, 3), vars(6), grad(6), grad_norm, &
                    red_space_hess(3, 3), red_space_hess_left_eigvecs(3, 3), &
                    red_space_hess_right_eigvecs(3, 3), red_space_hess_eigvals(3), &
                    red_space_hess_imag_eigvals(3), solution(6), &
                    red_space_solution(3), work(lwork), temp(3, 3)
        integer(ip) :: i, j, info

        external :: dsyev, dgeev

        ! assume tests pass
        test_newton_step = .true.

        ! point in quadratic region near minimum
        vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]

        ! calculate gradient and Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)

        ! defined a reduced space basis
        red_space_basis(:, 1) = grad / grad_norm
        red_space_basis(:, 2) = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        red_space_basis(:, 3) = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]

        ! orthonormalize reduced space basis
        do i = 2, 3
            do j = 1, i - 1
                ! subtract projection onto previous vectors
                red_space_basis(:, i) = red_space_basis(:, i) - &
                    dot_product(red_space_basis(:, i), red_space_basis(:, j)) * &
                    red_space_basis(:, j)
            end do
            red_space_basis(:, i) = red_space_basis(:, i) / norm2(red_space_basis(:, i))
        end do

        ! define augmented Hessian
        do i = 1, 3
            do j = 1, 3
                red_space_hess(i, j) = &
                    dot_product(red_space_basis(:, i), &
                                matmul(hess, red_space_basis(:, j)))
            end do
        end do

        ! diagonalize Hessian
        red_space_hess_right_eigvecs = red_space_hess
        call dsyev("V", "U", 3_ip, red_space_hess_right_eigvecs, 3_ip, &
                   red_space_hess_eigvals, work, lwork, info)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform Newton step, determine whether resulting solution is correct in 
        ! reduced and full space
        call newton_step(grad_norm, red_space_basis, red_space_hess_eigvals, &
                         red_space_hess_right_eigvecs, red_space_hess_left_eigvecs, &
                         solution, red_space_solution)
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution)) &
                  + grad) > tol) then
            write (stderr, *) "test_newton_step failed: Newton step does not "// &
                "satisfy Newton equation for symmetric Hessian."
            test_newton_step = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_newton_step failed: Full space solution not "// &
                "correct for symmetric Hessian."
            test_newton_step = .false.
        end if

        ! make Hessian non-symmetric
        red_space_hess(1, 2) = red_space_hess(1, 2) + 0.5_rp

        ! diagonalize non-symmetric Hessian
        temp = red_space_hess
        call dgeev("V", "V", 3_ip, temp, 3_ip, red_space_hess_eigvals, &
                   red_space_hess_imag_eigvals, red_space_hess_left_eigvecs, 3_ip, &
                   red_space_hess_right_eigvecs, 3_ip, work, lwork, info)

        ! biorthonormalize eigenvectors
        do i = 1, 3
            red_space_hess_left_eigvecs(:, i) = &
                red_space_hess_left_eigvecs(:, i) / &
                dot_product(red_space_hess_left_eigvecs(:, i), &
                            red_space_hess_right_eigvecs(:, i))
        end do

        ! perform Newton step, determine whether resulting solution is correct in 
        ! reduced and full space
        call newton_step(grad_norm, red_space_basis, red_space_hess_eigvals,  &
                         red_space_hess_right_eigvecs, red_space_hess_left_eigvecs, &
                         solution, red_space_solution)
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution)) &
                  + grad) > tol) then
            write (stderr, *) "test_newton_step failed: Newton step does not "// &
                "satisfy Newton equation for asymmetric Hessian."
            test_newton_step = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_newton_step failed: Full space solution not "// &
                "correct for asymmetric Hessian."
            test_newton_step = .false.
        end if

    end function test_newton_step

    logical(c_bool) function test_bisection_ah() bind(C)
        !
        ! this function tests the augmented Hessian bisection subroutine
        !
        use opentrustregion, only: solver_settings_type, bisection_ah

        type(solver_settings_type) :: settings
        integer(ip), parameter :: lwork = 12
        real(rp) :: red_space_basis(6, 3), vars(6), grad(6), grad_norm, &
                    red_space_hess(3, 3), red_space_hess_eigvals(3), &
                    red_space_hess_left_eigvecs(3, 3), &
                    red_space_hess_right_eigvecs(3, 3), solution(6), &
                    red_space_solution(3), trust_radius, mu, work(lwork)
        integer(ip) :: i, j, error, info

        external :: dsyev

        ! assume tests pass
        test_bisection_ah = .true.

        ! setup settings object
        call setup_settings(settings)

        ! choose target trust radius
        trust_radius = 1.0_rp

        ! point with strong negative curvature
        vars = [0.29_rp, 0.47_rp, 0.66_rp, 0.41_rp, 0.23_rp, 0.26_rp]

        ! calculate gradient and Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)

        ! defined a reduced space basis
        red_space_basis(:, 1) = grad / grad_norm
        red_space_basis(:, 2) = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        red_space_basis(:, 3) = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]

        ! orthonormalize reduced space basis
        do i = 2, 3
            do j = 1, i - 1
                ! subtract projection onto previous vectors
                red_space_basis(:, i) = red_space_basis(:, i) - &
                    dot_product(red_space_basis(:, i), red_space_basis(:, j)) * &
                    red_space_basis(:, j)
            end do
            red_space_basis(:, i) = red_space_basis(:, i) / norm2(red_space_basis(:, i))
        end do

        ! define reduced space Hessian
        red_space_hess = matmul(transpose(red_space_basis), matmul(hess, &
                                                                   red_space_basis))

        ! diagonalize Hessian
        red_space_hess_right_eigvecs = red_space_hess
        call dsyev("V", "U", 3_ip, red_space_hess_right_eigvecs, 3_ip, &
                   red_space_hess_eigvals, work, 9_ip, info)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection, check whether error has occured and determine whether 
        ! resulting solution is correct in reduced and full space and respects target 
        ! trust radius
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection_ah failed: Produced error."
            test_bisection_ah = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write (stderr, *) "test_bisection_ah failed: Solution does not respect "// &
                "trust radius."
            test_bisection_ah = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write (stderr, *) "test_bisection_ah failed: Step does not satisfy "// &
                "level-shifted Newton equation."
            test_bisection_ah = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_bisection_ah failed: Full space solution not "// &
                              "correct for symmetric Hessian."
            test_bisection_ah = .false.
        end if
        
        ! construct custom augmented Hessian with negative eigenvalue and no coupling 
        ! to the gradient
        red_space_hess = 0.0_rp
        red_space_hess(1, 1) = 1.0_rp   ! Mode 1 (coupled to gradient)
        red_space_hess(2, 2) = 2.0_rp   ! Mode 2 (positive curvature)
        red_space_hess(1, 2) = 0.5_rp   ! coupling between mode 1 and 2
        red_space_hess(2, 1) = 0.5_rp
        red_space_hess(3, 3) = -2.0_rp  ! Mode 3 (lowest eigenvalue, isolated from 
                                        ! gradient)
        settings%hess_symm = .true.
        
        ! diagonalize artificial Hessian
        red_space_hess_right_eigvecs = red_space_hess
        call dsyev("V", "U", 3_ip, red_space_hess_right_eigvecs, 3_ip, &
                   red_space_hess_eigvals, work, 9_ip, info)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection, check whether error has occured and determine whether 
        ! resulting solution respects target trust radius, contains lowest eigenvector 
        ! component and has level shift equal to lowest eigenvalue
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection_ah failed: Produced error for hard case."
            test_bisection_ah = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write (stderr, *) "test_bisection_ah failed: Hard case solution does "// &
                "not respect trust radius."
            test_bisection_ah = .false.
        end if
        if (abs(mu - red_space_hess_eigvals(1)) > tol) then
            write (stderr, *) "test_bisection_ah failed: Level shift does not "// &
                "equal the most negative eigenvalue and hard case."
            test_bisection_ah = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write (stderr, *) "test_bisection_ah failed: Hard case solution does "// &
                "not contain the level-shifted Newton component."
            test_bisection_ah = .false.
        end if
        if (abs(dot_product(red_space_solution, red_space_hess_right_eigvecs(:, 1))) < &
            tol) then
            write (stderr, *) "test_bisection_ah failed: Hard case solution does "// &
                "not contain a component of the lowest eigenvector."
            test_bisection_ah = .false.
        end if

        ! point in quadratic region near minimum
        vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]

        ! calculate gradient and Hessian to define augmented Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)
        red_space_hess = matmul(transpose(red_space_basis), matmul(hess, &
                                                                   red_space_basis))

        ! diagonalize Hessian
        red_space_hess_right_eigvecs = red_space_hess
        call dsyev("V", "U", 3_ip, red_space_hess_right_eigvecs, 3_ip, &
                   red_space_hess_eigvals, work, 9_ip, info)
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection and determine whether routine correctly throws error since
        ! minimum is closer than target trust radius and no level shift is necessary
        call bisection_ah(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_left_eigvecs, &
                          red_space_hess_right_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error == 0) then
            write (stderr, *) "test_bisection_ah failed: Failed to produce error."
            test_bisection_ah = .false.
        end if

    end function test_bisection_ah

    logical(c_bool) function test_bisection_mu() bind(C)
        !
        ! this function tests the level-shift bisection subroutine
        !
        use opentrustregion, only: solver_settings_type, bisection_mu

        type(solver_settings_type) :: settings
        integer(ip), parameter :: lwork = 12
        real(rp) :: red_space_basis(6, 3), vars(6), grad(6), grad_norm, &
                    red_space_hess(3, 3), red_space_hess_left_eigvecs(3, 3), &
                    red_space_hess_right_eigvecs(3, 3), solution(6), &
                    red_space_solution(3), trust_radius, mu, work(lwork), temp(3, 3), &
                    red_space_hess_eigvals_re(3), red_space_hess_eigvals_im(3)
        complex(rp) :: red_space_hess_eigvals(3)
        integer(ip) :: i, j, error, info

        external :: dsyev

        ! assume tests pass
        test_bisection_mu = .true.

        ! setup settings object
        call setup_settings(settings)

        ! choose target trust radius
        trust_radius = 1.0_rp

        ! point with strong negative curvature
        vars = [0.29_rp, 0.47_rp, 0.66_rp, 0.41_rp, 0.23_rp, 0.26_rp]

        ! calculate gradient and Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)

        ! defined a reduced space basis
        red_space_basis(:, 1) = grad / grad_norm
        red_space_basis(:, 2) = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        red_space_basis(:, 3) = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]

        ! orthonormalize reduced space basis
        do i = 2, 3
            do j = 1, i - 1
                ! subtract projection onto previous vectors
                red_space_basis(:, i) = red_space_basis(:, i) - &
                    dot_product(red_space_basis(:, i), red_space_basis(:, j)) * &
                    red_space_basis(:, j)
            end do
            red_space_basis(:, i) = red_space_basis(:, i) / norm2(red_space_basis(:, i))
        end do

        ! define reduced space Hessian
        red_space_hess = matmul(transpose(red_space_basis), matmul(hess, &
                                                                   red_space_basis))

        ! diagonalize Hessian
        red_space_hess_right_eigvecs = red_space_hess
        call dsyev("V", "U", 3_ip, red_space_hess_right_eigvecs, 3_ip, &
                   red_space_hess_eigvals_re, work, 9_ip, info)
        red_space_hess_eigvals%re = red_space_hess_eigvals_re
        red_space_hess_eigvals%im = 0.0_rp
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection, check whether error has occured and determine whether 
        ! resulting solution is correct in reduced and full space and respects target 
        ! trust radius
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection_mu failed: Produced error for "// &
                "symmetric, indefinite Hessian."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write (stderr, *) "test_bisection_mu failed: Solution does not respect "// &
                "trust radius for symmetric, indefinite Hessian."
            test_bisection_mu = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) &
                         - mu * red_space_solution) + grad) > tol) then
            write (stderr, *) "test_bisection_mu failed: Step does not satisfy "// &
                "level-shifted Newton equation for symmetric, indefinite Hessian."
            test_bisection_mu = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_bisection_mu failed: Full space solution not "// &
                              "correct for symmetric, indefinite Hessian."
            test_bisection_mu = .false.
        end if

        ! make Hessian non-symmetric
        red_space_hess(1, 2) = red_space_hess(1, 2) + 0.5_rp

        ! diagonalize non-symmetric Hessian
        temp = red_space_hess
        call dgeev("V", "V", 3_ip, temp, 3_ip, red_space_hess_eigvals_re, &
                   red_space_hess_eigvals_im, red_space_hess_left_eigvecs, 3_ip, &
                   red_space_hess_right_eigvecs, 3_ip, work, lwork, info)
        red_space_hess_eigvals%re = red_space_hess_eigvals_re
        red_space_hess_eigvals%im = red_space_hess_eigvals_im

        ! perform bisection, check whether error has occured and determine whether 
        ! resulting solution is correct in reduced and full space and respects target 
        ! trust radius
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection_mu failed: Produced error for "// &
                "asymmetric, indefinite Hessian."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write (stderr, *) "test_bisection_mu failed: Solution does not respect "// &
                "trust radius for asymmetric, indefinite Hessian."
            test_bisection_mu = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write (stderr, *) "test_bisection_mu failed: Step does not satisfy "// &
                "level-shifted Newton equation for asymmetric, indefinite Hessian."
            test_bisection_mu = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_bisection_mu failed: Full space solution not "// &
                              "correct for asymmetric, indefinite Hessian."
            test_bisection_mu = .false.
        end if
        
        ! construct custom augmented Hessian with negative eigenvalue and no coupling 
        ! to the gradient
        red_space_hess = 0.0_rp
        red_space_hess(1, 1) = 1.0_rp   ! Mode 1 (coupled to gradient)
        red_space_hess(2, 2) = 2.0_rp   ! Mode 2 (positive curvature)
        red_space_hess(1, 2) = 0.5_rp   ! coupling between mode 1 and 2
        red_space_hess(2, 1) = 0.5_rp
        red_space_hess(3, 3) = -2.0_rp  ! Mode 3 (lowest eigenvalue, isolated from 
                                        ! gradient)
        
        ! diagonalize artificial Hessian
        red_space_hess_right_eigvecs = red_space_hess
        call dsyev("V", "U", 3_ip, red_space_hess_right_eigvecs, 3_ip, &
                   red_space_hess_eigvals_re, work, 9_ip, info)
        red_space_hess_eigvals%re = red_space_hess_eigvals_re
        red_space_hess_eigvals%im = 0.0_rp
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection, check whether error has occured and determine whether 
        ! resulting solution respects target trust radius, contains lowest eigenvector 
        ! component and has level shift equal to lowest eigenvalue
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection_mu failed: Produced error for "// &
                "symmetric Hessian and hard case."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write (stderr, *) "test_bisection_mu failed: Hard case solution does "// &
                "not respect trust radius for symmetric Hessian."
            test_bisection_mu = .false.
        end if
        if (abs(mu - red_space_hess_eigvals(1)) > tol) then
            write (stderr, *) "test_bisection_mu failed: Level shift does not equal "// &
                "the most negative eigenvalue for symmetric Hessian and hard case."
            test_bisection_mu = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write (stderr, *) "test_bisection_mu failed: Hard case solution does "// &
                "not contain the level-shifted Newton component for symmetric Hessian."
            test_bisection_mu = .false.
        end if
        if (abs(dot_product(red_space_solution, red_space_hess_right_eigvecs(:, 1))) < &
            tol) then
            write (stderr, *) "test_bisection_mu failed: Hard case solution does "// &
                "not contain a component of the lowest eigenvector for symmetric "// &
                "Hessian."
            test_bisection_mu = .false.
        end if

        ! make Hessian non-symmetric
        red_space_hess(1, 2) = red_space_hess(1, 2) + 0.5_rp

        ! diagonalize non-symmetric Hessian
        temp = red_space_hess
        call dgeev("V", "V", 3_ip, temp, 3_ip, red_space_hess_eigvals_re, &
                   red_space_hess_eigvals_im, red_space_hess_left_eigvecs, 3_ip, &
                   red_space_hess_right_eigvecs, 3_ip, work, lwork, info)
        red_space_hess_eigvals%re = red_space_hess_eigvals_re
        red_space_hess_eigvals%im = red_space_hess_eigvals_im

        ! perform bisection, check whether error has occured and determine whether 
        ! resulting solution respects target trust radius, contains lowest eigenvector 
        ! component and has level shift equal to lowest eigenvalue
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection_mu failed: Produced error for "// &
                "asymmetric Hessian and hard case."
            test_bisection_mu = .false.
        end if
        if (abs(norm2(solution) - trust_radius) > tol) then
            write (stderr, *) "test_bisection_mu failed: Hard case solution does "// &
                "not respect trust radius for asymmetric Hessian."
            test_bisection_mu = .false.
        end if
        if (abs(mu - minval(red_space_hess_eigvals%re)) > tol) then
            write (stderr, *) "test_bisection_mu failed: Level shift does not "// &
                "equal the most negative eigenvalue for asymmetric Hessian and "// &
                "hard case."
            test_bisection_mu = .false.
        end if
        if (norm2(matmul(red_space_basis, matmul(red_space_hess, red_space_solution) - &
                         mu * red_space_solution) + grad) > tol) then
            write (stderr, *) "test_bisection_mu failed: Hard case solution does "// &
                "not contain the level-shifted Newton component for asymmetric Hessian."
            test_bisection_mu = .false.
        end if
        if (abs(dot_product(red_space_solution, red_space_hess_right_eigvecs(:, 1))) < &
            tol) then
            write (stderr, *) "test_bisection_mu failed: Hard case solution does "// &
                "not contain a component of the lowest eigenvector for asymmetric "// &
                "Hessian."
            test_bisection_mu = .false.
        end if

        ! point in quadratic region near minimum
        vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]

        ! calculate gradient and Hessian to define reduced Hessian
        call hartmann6d_gradient(vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(vars)
        red_space_hess = matmul(transpose(red_space_basis), matmul(hess, &
                                                                   red_space_basis))

        ! defined a reduced space basis
        red_space_basis(:, 1) = grad / grad_norm
        red_space_basis(:, 2) = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        red_space_basis(:, 3) = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]

        ! orthonormalize reduced space basis
        do i = 2, 3
            do j = 1, i - 1
                ! subtract projection onto previous vectors
                red_space_basis(:, i) = red_space_basis(:, i) - &
                    dot_product(red_space_basis(:, i), red_space_basis(:, j)) * &
                    red_space_basis(:, j)
            end do
            red_space_basis(:, i) = red_space_basis(:, i) / norm2(red_space_basis(:, i))
        end do

        ! diagonalize Hessian
        red_space_hess_right_eigvecs = red_space_hess
        call dsyev("V", "U", 3_ip, red_space_hess_right_eigvecs, 3_ip, &
                   red_space_hess_eigvals_re, work, 9_ip, info)
        red_space_hess_eigvals%re = red_space_hess_eigvals_re
        red_space_hess_eigvals%im = 0.0_rp
        red_space_hess_left_eigvecs = red_space_hess_right_eigvecs

        ! perform bisection and determine whether routine correctly throws error since
        ! minimum is closer than target trust radius and no level shift is necessary
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection_mu failed: Produced error for "// &
                              "symmetric Hessian."
            test_bisection_mu = .false.
        end if
        if (norm2(solution) > trust_radius) then
            write (stderr, *) "test_bisection_mu failed: Solution does not respect "// &
                "trust radius for symmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if
        if (abs(mu) > tol) then
            write (stderr, *) "test_bisection_mu failed: Level shift does not "// &
                "vanish for symmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if
        if (norm2(matmul(red_space_basis, &
                         matmul(red_space_hess, red_space_solution)) + grad) > tol) then
            write (stderr, *) "test_bisection_mu failed: Step does not satisfy "// &
                "Newton equation for symmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_bisection_mu failed: Full space solution not "// &
                "correct for symmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if

        ! make Hessian non-symmetric
        red_space_hess(1, 2) = red_space_hess(1, 2) + 0.5_rp

        ! diagonalize non-symmetric Hessian
        temp = red_space_hess
        call dgeev("V", "V", 3_ip, temp, 3_ip, red_space_hess_eigvals_re, &
                   red_space_hess_eigvals_im, red_space_hess_left_eigvecs, 3_ip, &
                   red_space_hess_right_eigvecs, 3_ip, work, lwork, info)
        red_space_hess_eigvals%re = red_space_hess_eigvals_re
        red_space_hess_eigvals%im = red_space_hess_eigvals_im

        ! perform bisection and determine whether routine correctly throws error since
        ! minimum is closer than target trust radius and no level shift is necessary
        call bisection_mu(red_space_hess, grad_norm, red_space_basis, &
                          red_space_hess_eigvals, red_space_hess_right_eigvecs, &
                          red_space_hess_left_eigvecs, trust_radius, solution, &
                          red_space_solution, mu, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bisection_mu failed: Produced error for "// &
                              "asymmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if
        if (norm2(solution) > trust_radius) then
            write (stderr, *) "test_bisection_mu failed: Solution does not respect "// &
                "trust radius for asymmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if
        if (abs(mu) > tol) then
            write (stderr, *) "test_bisection_mu failed: Level shift does not "// &
                "vanish for asymmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if
        if (norm2(matmul(red_space_basis, &
                         matmul(red_space_hess, red_space_solution)) + grad) > tol) then
            write (stderr, *) "test_bisection_mu failed: Step does not satisfy "// &
                "Newton equation for asymmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if
        if (any(abs(solution - matmul(red_space_basis, red_space_solution)) > tol)) then
            write (stderr, *) "test_bisection_mu failed: Full space solution not "// &
                "correct for asymmetric, positive-definite Hessian."
            test_bisection_mu = .false.
        end if

    end function test_bisection_mu

    logical(c_bool) function test_bracket() bind(C)
        !
        ! this function tests the bisection subroutine
        !
        use opentrustregion, only: solver_settings_type, obj_func_type, bracket

        type(solver_settings_type) :: settings
        procedure(obj_func_type), pointer :: obj_func_funptr
        real(rp) :: vars(6), lower, upper, n
        integer(ip) :: error

        ! assume tests pass
        test_bracket = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define procedure pointer
        obj_func_funptr => obj_func

        ! define direction
        vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]

        ! define lower and upper bound
        lower = 0.0_rp
        upper = 0.5_rp

        ! perform bracket and determine if new point decreases objective function in
        ! comparison to lower and upper bound
        n = bracket(obj_func_funptr, vars, lower, upper, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bracket failed: Produced error."
            test_bracket = .false.
        end if
        if (hartmann6d_func(n*vars) >= hartmann6d_func(lower*vars) .and. &
            hartmann6d_func(n*vars) >= hartmann6d_func(upper*vars)) then
            write (stderr, *) "test_bracket failed: Line search does not produce "// &
                "lower function value than starting points."
            test_bracket = .false.
        end if

        ! try different order
        n = bracket(obj_func_funptr, vars, upper, lower, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_bracket failed: Produced error."
            test_bracket = .false.
        end if
        if (hartmann6d_func(n*vars) >= hartmann6d_func(lower*vars) .and. &
            hartmann6d_func(n*vars) >= hartmann6d_func(upper*vars)) then
            write (stderr, *) "test_bracket failed: Line search does not produce "// &
                "lower function value than starting points."
            test_bracket = .false.
        end if

    end function test_bracket

    logical(c_bool) function test_extend_matrix() bind(C)
        !
        ! this function tests the subroutine for extending a matrix
        !
        use opentrustregion, only: extend_matrix

        real(rp), allocatable :: matrix(:, :)
        real(rp) :: expected(3, 3), vector1(3), vector2(3)

        ! assume tests pass
        test_extend_matrix = .true.

        ! allocate and initialize symmetric matrix and vector to be added
        allocate(matrix(2, 2))
        matrix = reshape([1.0_rp, 2.0_rp, &
                          2.0_rp, 3.0_rp], [2, 2])
        vector1 = [4.0_rp, 5.0_rp, 6.0_rp]
        vector2 = [7.0_rp, 8.0_rp, 6.0_rp]

        ! initialize expected matrix
        expected = reshape([1.0_rp, 2.0_rp, 7.0_rp, &
                            2.0_rp, 3.0_rp, 8.0_rp, &
                            4.0_rp, 5.0_rp, 6.0_rp], [3, 3])

        ! call routine and determine if dimensions and values of resulting matrix match
        call extend_matrix(matrix, vector1, vector2)
        if (size(matrix, 1) /= 3 .or. size(matrix, 2) /= 3) then
            write (stderr, *) "test_extend_matrix failed: Incorrect matrix "// &
                "dimensions after extending."
            test_extend_matrix = .false.
        end if
        if (norm2(matrix - expected) > tol) then
            write (stderr, *) "test_extend_matrix failed: Incorrect matrix values "// &
                "after extending."
            test_extend_matrix = .false.
        end if

        ! deallocate matrix
        deallocate(matrix)

    end function test_extend_matrix

    logical(c_bool) function test_add_column() bind(C)
        !
        ! this function tests the subroutine for adding a column to a matrix
        !
        use opentrustregion, only: add_column

        real(rp), allocatable :: matrix(:, :)
        real(rp) :: expected(3, 3), new_col(3)

        ! assume tests pass
        test_add_column = .true.

        ! allocate and initialize matrix and vector to be added
        allocate(matrix(3, 2))
        matrix = reshape([1.0_rp, 2.0_rp, 3.0_rp, &
                          4.0_rp, 5.0_rp, 6.0_rp], [3, 2])
        new_col = [7.0_rp, 8.0_rp, 9.0_rp]

        ! initialize expected matrix
        expected = reshape([1.0_rp, 2.0_rp, 3.0_rp, &
                            4.0_rp, 5.0_rp, 6.0_rp, &
                            7.0_rp, 8.0_rp, 9.0_rp], [3, 3])

        ! call routine and determine if dimensions and values of resulting matrix match
        call add_column(matrix, new_col)
        if (size(matrix, 1) /= 3 .or. size(matrix, 2) /= 3) then
            write (stderr, *) "test_add_column failed: Incorrect matrix dimensions "// &
                "after adding column."
            test_add_column = .false.
        end if
        if (norm2(matrix - expected) > tol) then
            write (stderr, *) "test_add_column failed: Incorrect matrix values "// &
                "after adding column."
            test_add_column = .false.
        end if

        ! deallocate matrix
        deallocate(matrix)

    end function test_add_column

    logical(c_bool) function test_mat_min_eig() bind(C)
        !
        ! this function tests the subroutine for determining the minimum eigenvalue and
        ! corresponding eigenvector of a matrix, for both the symmetric and the
        ! non-symmetric case
        !
        use opentrustregion, only: solver_settings_type, mat_min_eig

        type(solver_settings_type) :: settings
        real(rp) :: symm_matrix(3, 3), general_matrix(3, 3)
        real(rp) :: eigval, eigvec(3)
        integer(ip) :: error

        ! assume tests pass
        test_mat_min_eig = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize symmetric matrix
        symm_matrix = reshape([3.0_rp, 1.0_rp, 1.0_rp, &
                               1.0_rp, 4.0_rp, 2.0_rp, &
                               1.0_rp, 2.0_rp, 5.0_rp], [3, 3])

        ! call routine for the symmetric case and determine if lowest eigenvalue and
        ! corresponding eigenvector are found
        call mat_min_eig(symm_matrix, .true., eigval, eigvec, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_mat_min_eig failed: Produced error for "// &
                "symmetric matrix."
            test_mat_min_eig = .false.
        end if
        if (abs(eigval - 2.30797852837_rp) > tol) then
            write (stderr, *) "test_mat_min_eig failed: Incorrect minimum "// &
                "eigenvalue for symmetric matrix."
            test_mat_min_eig = .false.
        end if
        if (norm2(matmul(symm_matrix, eigvec) - eigval * eigvec) > tol) then
            write (stderr, *) "test_mat_min_eig failed: Incorrect eigenvector "// &
                "corresponding to minimum eigenvalue for symmetric matrix."
            test_mat_min_eig = .false.
        end if

        ! initialize non-symmetric matrix
        general_matrix = reshape([3.0_rp, 2.0_rp, 4.0_rp, &
                                  1.0_rp, 4.0_rp, 6.0_rp, &
                                  3.0_rp, 5.0_rp, 5.0_rp], [3, 3])

        ! call routine for the non-symmetric case and determine if lowest eigenvalue
        ! and corresponding right eigenvector are found
        call mat_min_eig(general_matrix, .false., eigval, eigvec, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_mat_min_eig failed: Produced error for "// &
                "non-symmetric matrix."
            test_mat_min_eig = .false.
        end if
        if (abs(eigval + 1.43567185527_rp) > tol) then
            write (stderr, *) "test_mat_min_eig failed: Incorrect minimum "// &
                "eigenvalue for non-symmetric matrix."
            test_mat_min_eig = .false.
        end if
        if (norm2(matmul(general_matrix, eigvec) - eigval * eigvec) > tol) then
            write (stderr, *) "test_mat_min_eig failed: Incorrect eigenvector "// &
                "corresponding to minimum eigenvalue for non-symmetric matrix."
            test_mat_min_eig = .false.
        end if

    end function test_mat_min_eig

    logical(c_bool) function test_symm_mat_diag() bind(C)
        !
        ! this function tests the function for determining the eigenvalues and 
        ! eigenvectors for a symmetric matrix
        !
        use opentrustregion, only: solver_settings_type, symm_mat_diag

        type(solver_settings_type) :: settings
        real(rp) :: matrix(3, 3), eigvals(3), eigvecs(3, 3)
        integer(ip) :: error

        ! assume tests pass
        test_symm_mat_diag = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize symmetric matrix
        matrix = reshape([3.0_rp, 1.0_rp, 1.0_rp, &
                          1.0_rp, 4.0_rp, 2.0_rp, &
                          1.0_rp, 2.0_rp, 5.0_rp], [3, 3])

        ! call function and determine if lowest eigenvalue is found
        call symm_mat_diag(matrix, eigvals, eigvecs, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_symm_mat_diag failed: Produced error."
            test_symm_mat_diag = .false.
        end if
        if (norm2(matmul(matrix, eigvecs) - eigvecs * spread(eigvals, dim=1, &
            ncopies=size(eigvecs, 1))) > tol) then
            write (stderr, *) "test_symm_mat_diag failed: Incorrect eigenvectors "// &
                "and eigenvalues for matrix."
            test_symm_mat_diag = .false.
        end if

    end function test_symm_mat_diag

    logical(c_bool) function test_general_mat_diag() bind(C)
        !
        ! this function tests the function for determining the eigenvalues and left and 
        ! right eigenvectors for a square matrix
        !
        use opentrustregion, only: solver_settings_type, general_mat_diag

        type(solver_settings_type) :: settings
        real(rp) :: matrix(3, 3), right_eigvecs(3, 3), left_eigvecs(3, 3), diag(3, 3), &
                    eye(3, 3)
        complex(rp) :: eigvals(3)
        integer(ip) :: error, i

        ! assume tests pass
        test_general_mat_diag = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize symmetric matrix
        matrix = reshape([3.0_rp, 2.0_rp, 4.0_rp, &
                          1.0_rp, 4.0_rp, 6.0_rp, &
                          3.0_rp, 5.0_rp, 5.0_rp], [3, 3])

        ! call function and determine if lowest eigenvalue is found
        call general_mat_diag(matrix, eigvals, right_eigvecs, left_eigvecs, settings, &
                              error)
        diag = 0.0_rp
        eye  = 0.0_rp
        do i = 1, 3
            diag(i, i) = eigvals(i)%re
            eye(i, i)  = 1.0_rp
        end do
        if (error /= 0) then
            write (stderr, *) "test_general_mat_diag failed: Produced error."
            test_general_mat_diag = .false.
        end if
        if (norm2(matmul(transpose(left_eigvecs), right_eigvecs) - eye) > tol) then
            write(stderr, *) "test_general_mat_diag failed: Left and right "// &
                "eigenvectors are not biorthonormal."
            test_general_mat_diag = .false.
        end if
        if (norm2(matmul(matrix, right_eigvecs) - matmul(right_eigvecs, diag)) > tol) &
            then
            write(stderr, *) "test_general_mat_diag failed: Incorrect right "// &
                "eigenvectors or eigenvalues for matrix."
            test_general_mat_diag = .false.
        end if
        if (norm2(matmul(transpose(left_eigvecs), matrix) - &
                  matmul(diag, transpose(left_eigvecs))) > tol) then
            write(stderr, *) "test_general_mat_diag failed: Incorrect left "// &
                "eigenvectors or eigenvalues for matrix."
            test_general_mat_diag = .false.
        end if

    end function test_general_mat_diag

    logical(c_bool) function test_mat_lowest_eigpairs() bind(C)
        !
        ! this function tests the subroutine which returns the lowest eigenpairs of a
        ! matrix
        !
        use opentrustregion, only: solver_settings_type, mat_lowest_eigpairs

        integer(ip), parameter :: n = 3, n_target = 2

        real(rp) :: symm_matrix(n, n), general_matrix(n, n), eigvals(n_target), &
                    eigvecs(n, n_target), expected_symm(n_target), &
                    expected_general(n_target), residual(n)
        integer(ip) :: error, i
        type(solver_settings_type) :: settings

        ! assume tests pass
        test_mat_lowest_eigpairs = .true.

        ! setup settings object
        call setup_settings(settings)

        ! symmetric diagonal matrix with well separated eigenvalues
        symm_matrix = 0.0_rp
        symm_matrix(1, 1) = 5.0_rp
        symm_matrix(2, 2) = 1.0_rp
        symm_matrix(3, 3) = 3.0_rp
        expected_symm = [1.0_rp, 3.0_rp]

        ! call routine for the symmetric case and determine whether the lowest
        ! eigenvalues are returned in ascending order
        call mat_lowest_eigpairs(symm_matrix, .true., eigvals, eigvecs, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_mat_lowest_eigpairs failed: Produced error for "// &
                "symmetric matrix."
            test_mat_lowest_eigpairs = .false.
            return
        end if
        if (norm2(eigvals - expected_symm) > tol) then
            write (stderr, *) "test_mat_lowest_eigpairs failed: Incorrect lowest "// &
                "eigenvalues for symmetric matrix."
            test_mat_lowest_eigpairs = .false.
        end if
        do i = 1, n_target
            residual = matmul(symm_matrix, eigvecs(:, i)) - eigvals(i) * eigvecs(:, i)
            if (norm2(residual) > tol) then
                write (stderr, *) "test_mat_lowest_eigpairs failed: Returned "// &
                    "eigenvectors for symmetric matrix are not eigenvectors."
                test_mat_lowest_eigpairs = .false.
            end if
        end do

        ! upper triangular, and therefore non-symmetric, matrix whose eigenvalues are
        ! its diagonal entries
        general_matrix = 0.0_rp
        general_matrix(1, 1) = 1.0_rp
        general_matrix(2, 2) = 5.0_rp
        general_matrix(3, 3) = 3.0_rp
        general_matrix(1, 2) = 2.0_rp
        general_matrix(1, 3) = 3.0_rp
        general_matrix(2, 3) = 4.0_rp
        expected_general = [1.0_rp, 3.0_rp]

        ! call routine for the non-symmetric case and determine whether the
        ! eigenvalues with the lowest real parts are returned in ascending order
        call mat_lowest_eigpairs(general_matrix, .false., eigvals, eigvecs, settings, &
                                 error)
        if (error /= 0) then
            write (stderr, *) "test_mat_lowest_eigpairs failed: Produced error for "// &
                "non-symmetric matrix."
            test_mat_lowest_eigpairs = .false.
            return
        end if
        if (norm2(eigvals - expected_general) > tol) then
            write (stderr, *) "test_mat_lowest_eigpairs failed: Incorrect "// &
                "eigenvalues with lowest real parts for non-symmetric matrix."
            test_mat_lowest_eigpairs = .false.
        end if
        do i = 1, n_target
            residual = matmul(general_matrix, &
                              eigvecs(:, i)) - eigvals(i) * eigvecs(:, i)
            if (norm2(residual) > tol) then
                write (stderr, *) "test_mat_lowest_eigpairs failed: Returned right "// &
                    "eigenvectors for non-symmetric matrix are not eigenvectors."
                test_mat_lowest_eigpairs = .false.
            end if
        end do

    end function test_mat_lowest_eigpairs

    logical(c_bool) function test_init_rng() bind(C)
        !
        ! this function tests the initialization subroutine for the random number
        ! generator
        !
        use opentrustregion, only: init_rng

        integer(ip) :: seed1, seed2, i
        real(rp) :: rand_seq1(5), rand_seq2(5), rand_seq3(5)

        ! assume tests pass
        test_init_rng = .true.

        ! define seeds
        seed1 = 12345
        seed2 = 67890

        ! call rng with first seed
        call init_rng(seed1)
        do i = 1, 5
            call random_number(rand_seq1(i))
        end do

        ! call rng with first seed
        call init_rng(seed1)
        do i = 1, 5
            call random_number(rand_seq2(i))
        end do

        ! call rng with second seed
        call init_rng(seed2)
        do i = 1, 5
            call random_number(rand_seq3(i))
        end do

        ! check reproducibility
        if (any(abs(rand_seq1 - rand_seq2) > tol)) then
            write (stderr, *) "test_init_rng failed: RNG does not produce "// &
                "consistent sequences for the same seed."
            test_init_rng = .false.
        end if

        ! check variation
        if (all(abs(rand_seq1 - rand_seq3) < tol)) then
            write (stderr, *) "test_init_rng failed: RNG produces identical "// &
                "sequences for different seeds."
            test_init_rng = .false.
        end if

    end function test_init_rng

    logical(c_bool) function test_generate_trial_vectors() bind(C)
        !
        ! this function tests the function which generates trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: solver_settings_type, generate_trial_vectors

        type(solver_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        real(rp) :: grad(4), h_diag(4), grad_norm
        integer(ip) :: error, i, j

        ! assume tests pass
        test_generate_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_random_trial_vectors = 2

        ! define gradient
        grad = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        grad_norm = norm2(grad)

        ! define all positive Hessian diagonal elements
        h_diag = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]

        ! generate trial vectors and determine whether function returns the correct
        ! number of orthonormal trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error /= 0) then
            write (stderr, *) "test_generate_trial_vectors failed: Produced error."
            test_generate_trial_vectors = .false.
        end if
        if (.not. allocated(red_space_basis)) then
            write (stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not allocated."
            test_generate_trial_vectors = .false.
            return
        end if
        if (size(red_space_basis, 2) /= 1 + settings%n_random_trial_vectors) then
            write (stderr, *) "test_generate_trial_vectors failed: Incorrect "// &
                "number of vectors for Hessian with only positive diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        do i = 1, size(red_space_basis, 2)
            do j = i + 1, size(red_space_basis, 2)
                if (abs(dot_product(red_space_basis(:, i), red_space_basis(:, j))) > &
                    tol) then
                    write (stderr, *) "test_generate_trial_vectors failed: "// &
                        "Generated vectors are not orthonormal for Hessian with "// &
                        "only positive diagonal elements."
                    test_generate_trial_vectors = .false.
                end if
            end do
        end do

        ! deallocate reduced space basis
        deallocate(red_space_basis)

        ! define Hessian diagonal with negative elements
        h_diag = [-1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]

        ! generate trial vectors and determine whether function returns the correct
        ! number of orthonormal trial vectors
        red_space_basis = generate_trial_vectors(grad, grad_norm, h_diag, settings, &
                                                 error)
        if (error /= 0) then
            write (stderr, *) "test_generate_trial_vectors failed: Produced error."
            test_generate_trial_vectors = .false.
        end if
        if (.not. allocated(red_space_basis)) then
            write (stderr, *) "test_generate_trial_vectors failed: Reduced space "// &
                "basis not allocated."
            test_generate_trial_vectors = .false.
            return
        end if
        if (size(red_space_basis, 2) /= 2 + settings%n_random_trial_vectors) then
            write (stderr, *) "test_generate_trial_vectors failed: Incorrect "// &
                "number of vectors for Hessian with diagonal elements."
            test_generate_trial_vectors = .false.
        end if
        do i = 1, size(red_space_basis, 2)
            do j = i + 1, size(red_space_basis, 2)
                if (abs(dot_product(red_space_basis(:, i), red_space_basis(:, j))) > &
                    tol) then
                    write (stderr, *) "test_generate_trial_vectors failed: "// &
                        "Generated vectors are not orthonormal for Hessian with "// &
                        "diagonal elements."
                    test_generate_trial_vectors = .false.
                end if
            end do
        end do

        ! deallocate reduced space basis
        deallocate(red_space_basis)

    end function test_generate_trial_vectors

    logical(c_bool) function test_generate_random_trial_vectors() bind(C)
        !
        ! this function tests the function which generates random trial vectors for the
        ! Davidson procedure
        !
        use opentrustregion, only: solver_settings_type, generate_random_trial_vectors

        integer(ip), parameter :: n_fill = 2

        type(solver_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        integer(ip) :: error, i, j

        ! assume tests pass
        test_generate_random_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings)

        ! allocate reduced space basis and set first normalized basis vector, which
        ! is the leading column the routine has to leave alone and orthogonalize the
        ! randomly filled trailing columns against
        allocate(red_space_basis(4, n_fill + 1))
        red_space_basis(:, 1) = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        red_space_basis(:, 1) = red_space_basis(:, 1) / norm2(red_space_basis(:, 1))

        ! generate trial vectors and determine whether function returns orthonormal
        ! trial vectors
        call generate_random_trial_vectors(red_space_basis, n_fill, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_generate_random_trial_vectors failed: Produced "// &
                "error."
            test_generate_random_trial_vectors = .false.
        end if
        if (abs(dot_product(red_space_basis(:, 1), &
                            [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp] / &
                            norm2([1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp])) - 1.0_rp) > &
            tol) then
            write (stderr, *) "test_generate_random_trial_vectors failed: Leading "// &
                "column was modified."
            test_generate_random_trial_vectors = .false.
        end if
        do i = 2, size(red_space_basis, 2)
            if (abs(norm2(red_space_basis(:, i)) - 1.0_rp) > tol) then
                write (stderr, *) "test_generate_random_trial_vectors failed: "// &
                    "Generated vectors are not normalized."
                test_generate_random_trial_vectors = .false.
            end if
        end do
        do i = 1, size(red_space_basis, 2)
            do j = i + 1, size(red_space_basis, 2)
                if (abs(dot_product(red_space_basis(:, i), red_space_basis(:, j))) > &
                    tol) then
                    write (stderr, *) "test_generate_random_trial_vectors failed: "// &
                        "Generated vectors are not orthonormal."
                    test_generate_random_trial_vectors = .false.
                end if
            end do
        end do

        ! deallocate reduced space basis
        deallocate(red_space_basis)

    end function test_generate_random_trial_vectors

    logical(c_bool) function test_orthogonalize_trial_vectors() bind(C)
        !
        ! this function tests the subroutine which orthogonalizes trial vectors
        !
        use opentrustregion, only: stability_settings_type, orthogonalize_trial_vectors

        type(stability_settings_type) :: settings
        real(rp), allocatable :: red_space_basis(:, :)
        integer(ip) :: error

        ! assume tests pass
        test_orthogonalize_trial_vectors = .true.

        ! setup settings object
        call setup_settings(settings)

        ! allocate reduced space basis with 3 vectors, where the third vector is 
        ! linearly dependent
        allocate(red_space_basis(4, 3))
        red_space_basis(:, 1) = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]
        red_space_basis(:, 2) = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
        red_space_basis(:, 3) = [1.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]

        ! orthogonalize trial vectors
        call orthogonalize_trial_vectors(red_space_basis, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_orthogonalize_trial_vectors failed: Produced error."
            test_orthogonalize_trial_vectors = .false.
        end if
        if (size(red_space_basis, 2) /= 2) then
            write (stderr, *) "test_orthogonalize_trial_vectors failed: Linearly "// &
                "dependent vector was not removed."
            test_orthogonalize_trial_vectors = .false.
        end if
        if (any(abs(norm2(red_space_basis, dim=1) - 1.0_rp) > tol)) then
            write (stderr, *) "test_orthogonalize_trial_vectors failed: Vectors "// &
                "are not normalized."
            test_orthogonalize_trial_vectors = .false.
        end if
        if (abs(dot_product(red_space_basis(:, 1), red_space_basis(:, 2))) > tol) then
            write (stderr, *) "test_orthogonalize_trial_vectors failed: Vectors "// &
                "are not orthogonal."
            test_orthogonalize_trial_vectors = .false.
            end if
        deallocate(red_space_basis)

    end function test_orthogonalize_trial_vectors

    logical(c_bool) function test_gram_schmidt() bind(C)
        !
        ! this function tests the Gram-Schmidt subroutine which orthonormalizes a 
        ! vector to a given basis
        !
        use opentrustregion, only: solver_settings_type, gram_schmidt, &
                                   gram_schmidt_zero_vector_error_msg, &
                                   gram_schmidt_lin_dep_error_msg, &
                                   gram_schmidt_too_many_vectors_error_msg

        type(solver_settings_type) :: settings
        real(rp) :: vector(4), lin_trans_vector(4), vector_small(2), space(4, 2), &
                    symm_matrix(4, 4), lin_trans_space(4, 2), space_small(2, 2)
        integer(ip) :: error

        ! assume tests pass
        test_gram_schmidt = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define vector to be orthogonalized and space
        vector = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        space(:, 1) = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
        space(:, 2) = [0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp]

        ! perform Gram-Schmidt orthogonalization and determine whether added vector is
        ! orthonormalized
        call gram_schmidt(vector, space, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_gram_schmidt failed: Produced error."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not orthogonal."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.0_rp) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not normalized."
            test_gram_schmidt = .false.
        end if

        ! define vector to be orthogonalized and space
        vector = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        space(:, 1) = [0.0_rp, 1.0_rp, 0.0_rp, 0.0_rp]
        space(:, 2) = [0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp]

        ! define symmetric linear transformation and corresponding vector and space
        symm_matrix = reshape([ 1.0_rp, -5.0_rp,  8.0_rp,  0.0_rp, &
                               -5.0_rp,  2.0_rp, -6.0_rp,  9.0_rp, &
                                8.0_rp, -6.0_rp,  3.0_rp, -7.0_rp, &
                                0.0_rp,  9.0_rp, -7.0_rp,  4.0_rp], &
                              shape(symm_matrix), order=[2,1])
        lin_trans_vector = matmul(symm_matrix, vector)
        lin_trans_space = matmul(symm_matrix, space)

        ! perform Gram-Schmidt orthogonalization and determine whether added vector is
        ! orthonormalized and linear transformation is correct
        call gram_schmidt(vector, space, settings, error, lin_trans_vector, &
                          lin_trans_space)
        if (error /= 0) then
            write (stderr, *) "test_gram_schmidt failed: Produced error."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not orthogonal."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.0_rp) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector not normalized."
            test_gram_schmidt = .false.
        end if
        if (sum(abs(lin_trans_vector - matmul(symm_matrix, vector))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added linear "// &
                "transformation not correct."
            test_gram_schmidt = .false.
        end if

        ! define zero vector
        vector = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! throws error
        call gram_schmidt(vector, space, settings, error)
        if ((error == 0) .or. &
            (adjustl(log_message) /= gram_schmidt_zero_vector_error_msg)) then
            write (stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization for zero vector."
            test_gram_schmidt = .false.
        end if

        ! define linearly dependent vector
        vector = space(:, 1)

        ! reset log message
        log_message = ""

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! throws error
        call gram_schmidt(vector, space, settings, error)
        if ((error /= 2) .or. &
            (adjustl(log_message) /= gram_schmidt_lin_dep_error_msg)) then
            write (stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization for linearly dependent vector."
            test_gram_schmidt = .false.
        end if

        ! define vector and space such that projecting after orthogonalizing
        ! reintroduces overlap with the space, forcing repeated passes through the
        ! loop and exercising the projector inside it, not just once beforehand
        vector = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        settings%project => mock_project

        ! perform Gram-Schmidt orthogonalization and determine whether the projector's
        ! invariant, orthogonality and normalization are all satisfied simultaneously
        call gram_schmidt(vector, space, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_gram_schmidt failed: Produced error with "// &
                "custom projection function."
            test_gram_schmidt = .false.
        end if
        if (abs(dot_product(vector, space(:, 1))) > tol .or. &
            abs(dot_product(vector, space(:, 2))) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector with custom "// &
                "projection function not orthogonal."
            test_gram_schmidt = .false.
        end if
        if (abs(norm2(vector) - 1.0_rp) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector with custom "// &
                "projection function not normalized."
            test_gram_schmidt = .false.
        end if
        if (abs(vector(1) - vector(2)) > tol) then
            write (stderr, *) "test_gram_schmidt failed: Added vector with custom "// &
                "projection function does not satisfy projector."
            test_gram_schmidt = .false.
        end if

        ! define vector in space that is already complete
        vector_small = [1.0_rp, 2.0_rp]
        space_small(:, 1) = [1.0_rp, 0.0_rp]
        space_small(:, 2) = [0.0_rp, 1.0_rp]

        ! reset log message
        log_message = ""

        ! perform Gram-Schmidt orthogonalization and determine if function correctly
        ! throws error
        call gram_schmidt(vector_small, space_small, settings, error)
        if (error == 0 .or. &
            (adjustl(log_message) /= gram_schmidt_too_many_vectors_error_msg)) then
            write (stderr, *) "test_gram_schmidt failed: No error returned during "// &
                "orthogonalization when number of vectors is larger than dimension "// &
                "of vector space."
            test_gram_schmidt = .false.
        end if

    end function test_gram_schmidt

    logical(c_bool) function test_init_solver_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the solver settings
        !
        use opentrustregion, only: solver_settings_type, &
                                   default_settings => default_solver_settings
        use test_reference, only: test_associated_solver_funptr, operator(/=)

        type(solver_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_init_solver_settings = .true.

        ! initialize settings
        call settings%init(error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_init_solver_settings failed: Function raised "// &
                "error."
            test_init_solver_settings = .false.
        end if

        ! check that callback function pointers are not associated
        if (.not. test_associated_solver_funptr(settings, &
                                                "test_init_solver_settings")) &
            test_init_solver_settings = .false.

        ! check settings
        if (settings /= default_settings) then
            write (stderr, *) "test_init_solver_settings failed: Settings not "// &
                "initialized correctly."
            test_init_solver_settings = .false.
        end if

    end function test_init_solver_settings

    logical(c_bool) function test_init_stability_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the stability check
        ! settings
        !
        use opentrustregion, only: stability_settings_type, &
                                   default_settings => default_stability_settings
        use test_reference, only: test_associated_stability_funptr, operator(/=)

        type(stability_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_init_stability_settings = .true.

        ! initialize settings
        call settings%init(error)

        ! check for error
        if (error /= 0) then
            write (stderr, *) "test_init_stability_settings failed: Function "// &
                "raised error."
            test_init_stability_settings = .false.
        end if

        ! check that callback function pointers are not associated
        if (.not. test_associated_stability_funptr(settings, &
                                                   "test_init_stability_settings")) &
            test_init_stability_settings = .false.

        ! check settings
        if (settings /= default_settings) then
            write (stderr, *) "test_init_stability_settings failed: Settings not "// &
                "initialized correctly."
            test_init_stability_settings = .false.
        end if

    end function test_init_stability_settings

    logical(c_bool) function test_level_shifted_diag_precond() bind(C)
        !
        ! this function tests the subroutine that constructs the level-shifted diagonal 
        ! preconditioner
        !
        use opentrustregion, only: solver_settings_type, level_shifted_diag_precond

        real(rp) :: vector(3), mu, h_diag(3), precond_vector(3), merged_val
        type(solver_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_level_shifted_diag_precond = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize quantities
        vector = [1.0_rp, 1.0_rp, 1.0_rp]
        mu = -2.0_rp
        h_diag = [-2.0_rp, 1.0_rp, 2.0_rp]

        ! call subroutine and check if results match
        call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, settings, &
                                        error)
        if (error /= 0) then
            write (stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "error for default preconditioner."
            test_level_shifted_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [1e10_rp, 1.0_rp / 3, 0.25_rp]) > tol)) then
            write (stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "preconditioned vector not correct for default preconditioner."
            test_level_shifted_diag_precond = .false.
        end if

        ! test custom projector
        settings%project => mock_project
        merged_val = 0.5_rp * (1e10_rp + 1.0_rp / 3)

        ! call subroutine and check if results match
        call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, settings, &
                                        error)
        if (error /= 0) then
            write (stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "error for custom projection function."
            test_level_shifted_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [merged_val, merged_val, 0.25_rp]) > tol)) then
            write (stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom projection function."
            test_level_shifted_diag_precond = .false.
        end if

        ! test custom preconditioner
        settings%precond => mock_precond

        ! call subroutine and check if results match
        call level_shifted_diag_precond(vector, mu, h_diag, precond_vector, settings, &
                                        error)
        if (error /= 0) then
            write (stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "error for custom preconditioner."
            test_level_shifted_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [-2.0_rp, -2.0_rp, -2.0_rp]) > tol)) then
            write (stderr, *) "test_level_shifted_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom preconditioner."
            test_level_shifted_diag_precond = .false.
        end if

    end function test_level_shifted_diag_precond

    logical(c_bool) function test_rel_floor_diag_precond() bind(C)
        !
        ! this function tests the subroutine that constructs the relative-floor
        ! absoulte diagonal preconditioner used to define the ellipsoidal trust-region 
        ! metricfor TCG/GLTR
        !
        use opentrustregion, only: solver_settings_type, rel_floor_diag_precond

        real(rp) :: vector(3), h_diag(3), precond_vector(3)
        type(solver_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_rel_floor_diag_precond = .true.

        ! setup settings object
        call setup_settings(settings)

        ! initialize quantities; the first element is far below the relative floor
        ! (precond_rel_floor_factor * max|h_diag| = 1e-2 * 4 = 0.04) and should get
        ! floored, the second is negative and should be floored (or not) based on its
        ! absolute value rather than being shifted away, the third is unaffected
        vector = [1.0_rp, 1.0_rp, 1.0_rp]
        h_diag = [0.001_rp, -1.0_rp, 4.0_rp]

        ! call subroutine and check if results match
        call rel_floor_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_rel_floor_diag_precond failed: Returned error "// &
                "for default preconditioner."
            test_rel_floor_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [25.0_rp, 1.0_rp, 0.25_rp]) > tol)) then
            write (stderr, *) "test_rel_floor_diag_precond failed: Returned "// &
                "preconditioned vector not correct for default preconditioner."
            test_rel_floor_diag_precond = .false.
        end if

        ! check that the relative floor falls back to the absolute floor when h_diag
        ! vanishes everywhere, rather than leaving the floor itself at zero
        h_diag = [0.0_rp, 0.0_rp, 0.0_rp]
        call rel_floor_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_rel_floor_diag_precond failed: Returned error "// &
                "for vanishing Hessian diagonal."
            test_rel_floor_diag_precond = .false.
        end if
        if (any(abs(precond_vector - 1e10_rp) > tol)) then
            write (stderr, *) "test_rel_floor_diag_precond failed: Returned "// &
                "preconditioned vector not correct for vanishing Hessian diagonal."
            test_rel_floor_diag_precond = .false.
        end if
        h_diag = [0.001_rp, -1.0_rp, 4.0_rp]

        ! test custom projector
        settings%project => mock_project

        ! call subroutine and check if results match
        call rel_floor_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_rel_floor_diag_precond failed: Returned error "// &
                "for custom projection function."
            test_rel_floor_diag_precond = .false.
        end if
        if (any(abs(precond_vector - [13.0_rp, 13.0_rp, 0.25_rp]) > tol)) then
            write (stderr, *) "test_rel_floor_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom projection function."
            test_rel_floor_diag_precond = .false.
        end if

        ! test custom positive-definite preconditioner
        settings%precond_pd => mock_precond_pd

        ! call subroutine and check if results match
        call rel_floor_diag_precond(vector, h_diag, precond_vector, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_rel_floor_diag_precond failed: Returned error "// &
                "for custom positive-definite preconditioner."
            test_rel_floor_diag_precond = .false.
        end if
        if (any(abs(precond_vector - 3.0_rp) > tol)) then
            write (stderr, *) "test_rel_floor_diag_precond failed: Returned "// &
                "preconditioned vector not correct for custom positive-definite "// &
                "preconditioner."
            test_rel_floor_diag_precond = .false.
        end if

    end function test_rel_floor_diag_precond

    logical(c_bool) function test_get_precond_level_shift() bind(C)
        !
        ! this function tests the function that computes the level shift for the 
        ! default level-shifted diagonal preconditioner
        !
        use opentrustregion, only: get_precond_level_shift, precond_factor

        real(rp) :: h_diag(3)

        ! assume tests pass
        test_get_precond_level_shift = .true.

        ! initialize quantities
        h_diag = [-2.0_rp, 1.0_rp, 2.0_rp]

        ! call function and check if results match for negative minimum diagonal element
        if (abs(get_precond_level_shift(h_diag) - &
                (-2.0_rp - precond_factor * 5.0_rp / 3)) > tol) then
            write (stderr, *) "test_get_precond_level_shift failed: Returned level "// &
                "shift not correct for negative minimum diagonal element."
            test_get_precond_level_shift = .false.
        end if

        ! call function and check if results match for positive minimum diagonal element
        h_diag(1) = 1.0_rp
        if (abs(get_precond_level_shift(h_diag) - 0.0_rp) > tol) then
            write (stderr, *) "test_get_precond_level_shift failed: Returned level "// &
                "shift not correct for positive minimum diagonal element."
            test_get_precond_level_shift = .false.
        end if

    end function test_get_precond_level_shift

    logical(c_bool) function test_orthogonal_projection() bind(C)
        !
        ! this function tests the orthogonal projection function which removes a 
        ! certain direction from a vector
        !
        use opentrustregion, only: orthogonal_projection

        real(rp), dimension(4) :: vector, direction

        ! assume tests pass
        test_orthogonal_projection = .true.

        ! define vector and direction to be projected out, the latter needs to be 
        ! normalized
        vector = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp]
        direction = [0.0_rp, 1.0_rp, 2.0_rp, 0.0_rp] / sqrt(5.0_rp)

        ! perform orthogonal projection and determine whether vector contains direction
        vector = orthogonal_projection(vector, direction)
        if (abs(dot_product(vector, direction)) > tol) then
            write (stderr, *) "test_orthogonal_projection failed: Vector contains "// &
                "component from direction to be projected out."
            test_orthogonal_projection = .false.
        end if

    end function test_orthogonal_projection

    logical(c_bool) function test_jacobi_davidson_correction() bind(C)
        !
        ! this function tests the Jacobi-Davidson correction subroutine
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, &
                                   jacobi_davidson_correction

        type(solver_settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(6) :: vars, vector, solution, corr_vector, hess_vector
        integer(ip) :: error

        ! assume tests pass
        test_jacobi_davidson_correction = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define point near saddle point, define trial vector, and solution to be 
        ! projected out
        vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
        vector = [0.1_rp, 0.2_rp, 0.3_rp, 0.4_rp, 0.5_rp, 0.6_rp]
        solution = [1.0_rp, -2.0_rp, 2.0_rp, -1.0_rp, 1.0_rp, -2.0_rp]

        ! generate Hessian
        call hartmann6d_hessian(vars)

        ! define Hessian linear transformation
        hess_x_funptr => hess_x_fun

        ! calculate Jacobi-Davidson correction and compare values
        call jacobi_davidson_correction(hess_x_funptr, vector, solution, 0.5_rp, &
                                        corr_vector, hess_vector, settings, error)
        if (error /= 0) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned error."
            test_jacobi_davidson_correction = .false.
        end if
        if (sum(abs(corr_vector - [-96.940677944_rp, 203.929698480_rp, &
                                   -216.199768920_rp, 100.656941418_rp, &
                                   -90.624469448_rp, 212.045768918_rp])) > 1e-8_rp) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "correction vector wrong."
            test_jacobi_davidson_correction = .false.
        end if
        if (sum(abs(hess_vector - [14.407362159_rp, -18.566381727_rp, &
                                   6.546311286_rp, -10.441098685_rp, &
                                   20.923570656_rp, -10.250311288_rp])) > 1e-8_rp) then
            write (stderr, *) "test_jacobi_davidson_correction failed: Returned "// &
                "Hessian linear transformation wrong."
            test_jacobi_davidson_correction = .false.
        end if

    end function test_jacobi_davidson_correction

    logical(c_bool) function test_minres() bind(C)
        !
        ! this function tests the minimum residual method subroutine
        !
        use opentrustregion, only: solver_settings_type, hess_x_type, minres

        type(solver_settings_type) :: settings
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp), dimension(6) :: vars, rhs, solution, vector, hess_vector, corr_vector
        real(rp) :: mu
        real(rp), parameter :: rtol = 1e-14_rp
        integer(ip) :: error

        ! assume tests pass
        test_minres = .true.

        ! setup settings object
        call setup_settings(settings)

        ! define point near saddle point
        vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]

        ! define solution to be projected out
        solution = [1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp, 5.0_rp, 6.0_rp] / sqrt(91.0_rp)

        ! generate Hessian
        call hartmann6d_hessian(vars)

        ! define Hessian linear transformation
        hess_x_funptr => hess_x_fun

        ! define Rayleigh quotient
        mu = dot_product(solution, hartmann6d_hess_x(solution))

        ! define right hand side based on residual, this ensures rhs is orthogonal to 
        ! solution if mu describes the Rayleigh quotient and solution is normalized
        rhs = hartmann6d_hess_x(solution) - mu * solution

        ! run minimum residual method, check if Jacobi-Davidson correction equation is 
        ! solved and whether Hessian linear transformation is correct, if rhs and 
        ! solution are orthogonal (as in Jacobi-Davidson), the final vector will be 
        ! orthogonal to the solution vector and consequently the Hessian linear 
        ! transformation of the projected vector is equivalent to the Hessian linear 
        ! transformation of the vector itself
        call minres(-rhs, hess_x_funptr, solution, mu, rtol, vector, hess_vector, &
                    settings, error)
        if (error /= 0) then
            write (stderr, *) "test_minres failed: Returned error."
            test_minres = .false.
        end if
        corr_vector = vector - dot_product(vector, solution) * solution
        corr_vector = hartmann6d_hess_x(corr_vector) - mu * corr_vector
        corr_vector = corr_vector - dot_product(corr_vector, solution) * solution
        if (sum(abs(corr_vector + rhs)) > tol) then
            write (stderr, *) "test_minres failed: Returned solution does not "// & 
                "solve Jacobi-Davidson correction equation."
            test_minres = .false.
        end if
        if (sum(abs(hess_vector + dot_product(vector, solution) * &
                    hartmann6d_hess_x(solution) - hartmann6d_hess_x(vector))) > tol) &
            then
            write (stderr, *) "test_minres failed: Returned Hessian linear "// &
                "transformation wrong."
            test_minres = .false.
        end if

        ! run minimum residual method for vanishing right hand side
        rhs = 0.0_rp
        call minres(-rhs, hess_x_funptr, solution, mu, rtol, vector, hess_vector, &
                    settings, error)
        if (error /= 0) then
            write (stderr, *) "test_minres failed: Returned error."
            test_minres = .false.
        end if
        if (sum(abs(vector)) > tol) then
            write (stderr, *) "test_minres failed: Returned solution is not zero "// &
                "for a vanishing rhs."
            test_minres = .false.
        end if
        if (sum(abs(hess_vector)) > tol) then
            write (stderr, *) "test_minres failed: Returned Hessian linear "// &
                "transformation is not zero for a vanishing rhs."
            test_minres = .false.
        end if

    end function test_minres

    logical(c_bool) function test_print_results() bind(C)
        !
        ! this function tests the subroutine that prints the result table
        !
        use opentrustregion, only: solver_settings_type

        type(solver_settings_type) :: settings

        ! assume tests pass
        test_print_results = .true.

        ! setup settings object
        call setup_settings(settings)

        ! print row of results table without optional arguments and check if row is 
        ! correct
        call settings%print_results(1_ip, 2.0_rp, 3.0_rp)
        if (log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |      -      |        -   |        -     |      -   ") then
            write (stderr, *) "test_print_results failed: Printed row without "// &
                "optional arguments not correct."
            test_print_results = .false.
        end if

        ! reset log message
        log_message = ""

        ! print row of results table with optional arguments
        call settings%print_results(1_ip, 2.0_rp, 3.0_rp, 4.0_rp, 5_ip, 6_ip, 7.0_rp, &
                                    8.0_rp)
        if (log_message /= "        1   |     2.00000000000000E+00   "// &
            "|   3.00E+00   |   4.00E+00  |   0 |   5  |    7.00E+00  |  8.00E+00") then
            write (stderr, *) "test_print_results failed: Printed row with "// &
                "optional arguments not correct."
            test_print_results = .false.
        end if

    end function test_print_results

    logical(c_bool) function test_print_message() bind(C)
        !
        ! this function tests the logging subroutine
        !
        use opentrustregion, only: solver_settings_type, print_message, &
                                   verbosity_debug, verbosity_warning

        type(solver_settings_type) :: settings

        ! assume tests pass
        test_print_message = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if logging is correctly performed according to verbosity level when 
        ! logger is provided
        call print_message(settings, "This is a test message.", verbosity_warning)
        if (trim(log_message) /= " This is a test message.") then
            write (stderr, *) "test_print_message failed: Log message is not "// &
                "printed correctly even though it should be according to verbosity "// &
                "level."
            test_print_message = .false.
        end if
        call print_message(settings, "This is another test message.", verbosity_debug)
        if (log_message == " This is another test message.") then
            write (stderr, *) "test_print_message failed: Log message is printed "// &
                "even though it should not be according to verbosity level."
            test_print_message = .false.
        end if

    end function test_print_message

    logical(c_bool) function test_split_string_by_space() bind(C)
        !
        ! this function tests the subroutine which splits strings after a space if 
        ! they exceed a given maximum length
        !
        use opentrustregion, only: split_string_by_space

        character(23), parameter :: message = "This is a test message."
        character(:), allocatable :: substrings(:)

        ! assume tests pass
        test_split_string_by_space = .true.

        ! check if strings are split correctly on spaces
        call split_string_by_space(message, 8_ip, substrings)
        if (size(substrings) == 3) then
            if (trim(substrings(1)) /= "This is" .or. trim(substrings(2)) /= "a test" &
                .or. trim(substrings(3)) /= "message.") then
                write (stderr, *) "test_split_string_by_space failed: Split "// &
                    "strings incorrect."
                test_split_string_by_space = .false.
            end if
        else
            write (stderr, *) "test_split_string_by_space failed: Number of "// &
                "substrings incorrect."
            test_split_string_by_space = .false.
        end if

        ! check if strings are split correctly if splitting on spaces is not possible
        call split_string_by_space(message, 5_ip, substrings)
        if (size(substrings) == 5) then
            if (trim(substrings(1)) /= "This" .or. trim(substrings(2)) /= "is a" .or. &
                trim(substrings(3)) /= "test" .or. trim(substrings(4)) /= "messa" .or. &
                trim(substrings(5)) /= "ge.") then
                write (stderr, *) "test_split_string_by_space failed: Split "// &
                    "strings incorrect."
                test_split_string_by_space = .false.
            end if
        else
            write (stderr, *) "test_split_string_by_space failed: Number of "// &
                "substrings incorrect."
            test_split_string_by_space = .false.
        end if

    end function test_split_string_by_space

    logical(c_bool) function test_accept_trust_region_step() bind(C)
        !
        ! this function tests the subroutine which determines whether to accept a 
        ! trust-region step
        !
        use opentrustregion, only: solver_settings_type, accept_trust_region_step, &
                                   trust_radius_shrink_ratio, &
                                   trust_radius_expand_ratio, &
                                   trust_radius_shrink_factor, &
                                   trust_radius_expand_factor

        logical :: accept_step, max_precision_reached
        real(rp) :: solution(3), trust_radius, start_trust_radius
        type(solver_settings_type) :: settings

        ! assume tests pass
        test_accept_trust_region_step = .true.

        ! setup settings object
        call setup_settings(settings)

        solution = [0.3_rp, 0.3_rp, 0.3_rp]

        ! check if step is rejected and trust radius is correctly reduced if micro 
        ! iterations have not converged
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, norm2(solution), 1.0_rp, &
                                               1.0_rp, .false., settings, trust_radius, &
                                               max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// &
                "or trust radius not correctly reduced when micro iterations have "// &
                "not converged."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and trust radius is correctly reduced if ratio is 
        ! negative
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, norm2(solution), -1.0_rp, &
                                               1.0_rp, .true., settings, trust_radius, &
                                               max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// & 
                "or trust radius not correctly reduced when ratio is negative."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is rejected and trust radius is correctly reduced if 
        ! individual rotations are too large
        trust_radius = 1.0_rp
        solution(1) = 1.0_rp
        accept_step = accept_trust_region_step(solution, norm2(solution), 1.0_rp, &
                                               1.0_rp, .true., settings, trust_radius, &
                                               max_precision_reached)
        if (accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step accepted "// &
                "or trust radius not correctly reduced when individual rotations "// &
                "are too large."
            test_accept_trust_region_step = .false.
        end if
        solution(1) = 0.3_rp

        ! check if step is accepted and trust radius is correctly reduced if ratio is 
        ! too small
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, norm2(solution), &
                                               0.9_rp * trust_radius_shrink_ratio, &
                                               1.0_rp, .true., settings, trust_radius, &
                                               max_precision_reached)
        if (.not. accept_step .or. abs(trust_radius - trust_radius_shrink_factor) > &
            tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted or trust radius not correctly reduced when ratio is too "// &
                "small."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is accepted and trust radius is correctly reduced if ratio is 
        ! ok
        start_trust_radius = 1.0_rp
        trust_radius = start_trust_radius
        accept_step = accept_trust_region_step(solution, norm2(solution), &
                                               0.5_rp * (trust_radius_shrink_ratio + &
                                               trust_radius_expand_ratio), 1.0_rp, &
                                               .true., settings, trust_radius, &
                                               max_precision_reached)
        if (.not. accept_step .or. abs(trust_radius - start_trust_radius) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted or trust radius changed when ratio is acceptable."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is accepted and trust radius is correctly expanded if ratio is
        ! too large and the step reaches the trust region boundary
        start_trust_radius = 0.5_rp
        trust_radius = start_trust_radius
        accept_step = accept_trust_region_step(solution, norm2(solution), &
                                               1.1_rp * trust_radius_expand_ratio, &
                                               1.0_rp, .true., settings, trust_radius, &
                                               max_precision_reached)
        if (.not. accept_step .or. abs(trust_radius - trust_radius_expand_factor * &
            start_trust_radius) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted or trust radius not correctly expanded when ratio is too "// &
                "large."
            test_accept_trust_region_step = .false.
        end if

        ! check if step is accepted but trust radius is left unchanged if ratio is too
        ! large but the step does not reach the trust region boundary
        start_trust_radius = 1.0_rp
        trust_radius = start_trust_radius
        accept_step = accept_trust_region_step(solution, norm2(solution), &
                                               1.1_rp * trust_radius_expand_ratio, &
                                               1.0_rp, .true., settings, trust_radius, &
                                               max_precision_reached)
        if (.not. accept_step .or. abs(trust_radius - start_trust_radius) > tol) then
            write(stderr, *) "test_accept_trust_region_step failed: Step not "// &
                "accepted or trust radius incorrectly expanded when ratio is too "// &
                "large but step does not reach the trust region boundary."
            test_accept_trust_region_step = .false.
        end if

        ! check if maximum precision is correctly handled for numerically vanishing 
        ! function improvement
        trust_radius = 1.0_rp
        accept_step = accept_trust_region_step(solution, norm2(solution), 0.0_rp, &
                                               0.0_rp, .true., settings, trust_radius, &
                                               max_precision_reached)
        if (.not. max_precision_reached) then
            write(stderr, *) "test_accept_trust_region_step failed: Maximum "// &
                "precision not correctly handled for numerically vanishing "// &
                "function improvement."
            test_accept_trust_region_step = .false.
        end if

        ! check if maximum precision is correctly handled for numerically vanishing 
        ! trust radius
        trust_radius = 0.0_rp
        accept_step = accept_trust_region_step(solution, norm2(solution), 1.0_rp, &
                                               1.0_rp, .false., settings, trust_radius, &
                                               max_precision_reached)
        if (.not. max_precision_reached) then
            write(stderr, *) "test_accept_trust_region_step failed: Maximum "// &
                "precision not correctly handled for numerically vanishing trust "// &
                "radius."
            test_accept_trust_region_step = .false.
        end if

    end function test_accept_trust_region_step

    logical(c_bool) function test_init_defaults() bind(C)
        !
        ! this function tests the subroutine that resolves settings defaults that
        ! depend on other settings
        !
        use opentrustregion, only: solver_settings_type, init_defaults, &
                                   default_spherical_trust_radius, &
                                   default_ellipsoidal_trust_radius

        type(solver_settings_type) :: settings

        ! assume tests pass
        test_init_defaults = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check that an unset trust region shape resolves to spherical and the
        ! starting trust radius resolves to the spherical default for a
        ! Davidson-based subsystem solver
        settings%subsystem_solver = "davidson_ls"
        settings%trust_region_shape = "none"
        settings%start_trust_radius = 0.0_rp
        call init_defaults(settings)
        if (settings%trust_region_shape /= "spherical") then
            write (stderr, *) "test_init_defaults failed: Did not resolve unset "// &
                "trust region shape to spherical for Davidson-based subsystem solver."
            test_init_defaults = .false.
        end if
        if (abs(settings%start_trust_radius - default_spherical_trust_radius) > tol) &
            then
            write (stderr, *) "test_init_defaults failed: Did not resolve unset "// &
                "starting trust radius to spherical default for Davidson-based "// &
                "subsystem solver."
            test_init_defaults = .false.
        end if

        ! check that an unset trust region shape resolves to ellipsoidal and the
        ! starting trust radius resolves to the ellipsoidal default for a
        ! non-Davidson-based subsystem solver
        settings%subsystem_solver = "tcg"
        settings%trust_region_shape = "none"
        settings%start_trust_radius = 0.0_rp
        call init_defaults(settings)
        if (settings%trust_region_shape /= "ellipsoidal") then
            write (stderr, *) "test_init_defaults failed: Did not resolve unset "// &
                "trust region shape to ellipsoidal for non-Davidson-based "// &
                "subsystem solver."
            test_init_defaults = .false.
        end if
        if (abs(settings%start_trust_radius - default_ellipsoidal_trust_radius) > &
            tol) then
            write (stderr, *) "test_init_defaults failed: Did not resolve unset "// &
                "starting trust radius to ellipsoidal default for "// &
                "non-Davidson-based subsystem solver."
            test_init_defaults = .false.
        end if

        ! check that an explicitly set trust region shape and starting trust radius
        ! are left untouched
        settings%subsystem_solver = "tcg"
        settings%trust_region_shape = "spherical"
        settings%start_trust_radius = 0.7_rp
        call init_defaults(settings)
        if (settings%trust_region_shape /= "spherical") then
            write (stderr, *) "test_init_defaults failed: Explicitly set trust "// &
                "region shape was overwritten."
            test_init_defaults = .false.
        end if
        if (abs(settings%start_trust_radius - 0.7_rp) > tol) then
            write (stderr, *) "test_init_defaults failed: Explicitly set starting "// &
                "trust radius was overwritten."
            test_init_defaults = .false.
        end if

    end function test_init_defaults

    logical(c_bool) function test_solver_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the 
        ! solver
        !
        use opentrustregion, only: solver_settings_type, solver_sanity_check, &
                                   project_warning_msg, subsystem_solvers, &
                                   trust_region_shapes

        type(solver_settings_type) :: settings
        real(rp) :: grad(3)
        integer(ip) :: error, i

        ! assume tests pass
        test_solver_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! set trust region shape to a valid value so that the checks below do not
        ! trip over the unresolved default; the trust region shape checks themselves
        ! are further down
        settings%trust_region_shape = "spherical"

        ! check if error is incorrectly thrown for finite and non-negative number of 
        ! parameters
        settings%n_random_trial_vectors = 0
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "non-negative and non-vanishing number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if error is correctly thrown for vanishing number of parameters
        call solver_sanity_check(settings, 0_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for vanishing number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if error is correctly thrown for negative number of parameters
        call solver_sanity_check(settings, -1_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for negative number of parameters."
            test_solver_sanity_check = .false.
        end if

        ! check if number of random trial vectors is reduced correctly
        settings%n_random_trial_vectors = 3
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (settings%n_random_trial_vectors /= 1) then
            write(stderr, *) "test_solver_sanity_check failed: Number of random "// &
                "trial vectors not correctly set."
            test_solver_sanity_check = .false.
        end if

        ! check if gradient size is treated correctly
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error thrown for "// &
                "gradient size."
            test_solver_sanity_check = .false.
        end if
        call solver_sanity_check(settings, 4_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for correct incorrect gradient size."
            test_solver_sanity_check = .false.
        end if

        ! check if subsystem solver is correctly checked
        do i = 1, size(subsystem_solvers)
            settings%subsystem_solver = subsystem_solvers(i)
            call solver_sanity_check(settings, 3_ip, grad, error)
            if (error /= 0) then
                write(stderr, *) "test_solver_sanity_check failed: Error thrown "// &
                    "for " // trim(subsystem_solvers(i)) // " subsystem solver."
                test_solver_sanity_check = .false.
            end if
        end do
        settings%subsystem_solver = "unknown"
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for unknown subsystem solver."
            test_solver_sanity_check = .false.
        end if
        ! check if trust region shape is correctly checked; use a subsystem solver that 
        ! supports both shapes so this loop is not confounded by the shape/subsystem 
        ! solver compatibility check tested further below
        settings%subsystem_solver = "tcg"
        do i = 1, size(trust_region_shapes)
            settings%trust_region_shape = trust_region_shapes(i)
            call solver_sanity_check(settings, 3_ip, grad, error)
            if (error /= 0) then
                write(stderr, *) "test_solver_sanity_check failed: Error thrown "// &
                    "for " // trim(trust_region_shapes(i)) // " trust region shape."
                test_solver_sanity_check = .false.
            end if
        end do
        settings%trust_region_shape = "unknown"
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for unknown trust region shape."
            test_solver_sanity_check = .false.
        end if

        ! check that an ellipsoidal trust region is correctly rejected for a
        ! Davidson-based subsystem solver, which does not support it
        settings%subsystem_solver = "davidson_ls"
        settings%trust_region_shape = "ellipsoidal"
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (error == 0) then
            write(stderr, *) "test_solver_sanity_check failed: Error not thrown "// &
                "for ellipsoidal trust region with Davidson-based subsystem solver."
            test_solver_sanity_check = .false.
        end if
        settings%trust_region_shape = "spherical"

        ! check if Hessian symmetry is correctly handled
        settings%hess_symm = .false.
        do i = 1, size(subsystem_solvers)
            if (index(subsystem_solvers(i), "ls") > 0) cycle
            settings%subsystem_solver = subsystem_solvers(i)
            call solver_sanity_check(settings, 3_ip, grad, error)
            if (error == 0) then
                write(stderr, *) "test_solver_sanity_check failed: Error not "// &
                    "thrown for non-symmetric Hessian with " // &
                    trim(subsystem_solvers(i)) // " solver."
                test_solver_sanity_check = .false.
            end if
        end do

        ! reset log message
        log_message = ""

        ! check if preconditioner and projector are correctly checked
        settings%subsystem_solver = "davidson_ls"
        settings%project => mock_project
        call solver_sanity_check(settings, 3_ip, grad, error)
        if (adjustl(log_message) /= project_warning_msg) then
            write(stderr, *) "test_solver_sanity_check failed: Warning message "// &
                "not correctly printed when custom projecting function is set."
            test_solver_sanity_check = .false.
        end if

    end function test_solver_sanity_check

    logical(c_bool) function test_stability_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the 
        ! stability check
        !
        use opentrustregion, only: stability_settings_type, stability_sanity_check, &
                                   project_warning_msg

        type(stability_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_stability_sanity_check = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if number of random trial vectors is reduced correctly
        settings%n_random_trial_vectors = 3
        call stability_sanity_check(settings, 3_ip, error)
        if (settings%n_random_trial_vectors /= 1) then
            write(stderr, *) "test_stability_sanity_check failed: Number of random "// &
                "trial not correctly set."
            test_stability_sanity_check = .false.
        end if

        ! check if subsystem solver is correctly checked
        settings%diag_solver = "davidson"
        call stability_sanity_check(settings, 3_ip, error)
        if (error /= 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error thrown for "// &
                "davidson diagonalization solver."
            test_stability_sanity_check = .false.
        end if
        settings%diag_solver = "jacobi-davidson"
        call stability_sanity_check(settings, 3_ip, error)
        if (error /= 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error thrown for "// &
                "jacobi-davidson diagonalization solver."
            test_stability_sanity_check = .false.
        end if
        settings%diag_solver = "unknown"
        call stability_sanity_check(settings, 3_ip, error)
        if (error == 0) then
            write(stderr, *) "test_stability_sanity_check failed: Error not thrown "// &
                "for unknown diagonalization solver."
            test_stability_sanity_check = .false.
        end if

        ! reset log message
        log_message = ""

        ! check if preconditioner and projector are correctly checked
        settings%diag_solver = "davidson"
        settings%project => mock_project
        call stability_sanity_check(settings, 3_ip, error)
        if (adjustl(log_message) /= project_warning_msg) then
            write(stderr, *) "test_stability_sanity_check failed: Warning message "// &
                "not correctly printed when custom projecting function is set."
            test_stability_sanity_check = .false.
        end if

    end function test_stability_sanity_check

    logical(c_bool) function test_level_shifted_davidson() bind(C)
        !
        ! this function tests the level-shifted Davidson subroutine
        !
        use opentrustregion, only: solver_settings_type, subsystem_solvers

        type(solver_settings_type) :: settings
        integer(ip) :: i

        ! assume tests pass
        test_level_shifted_davidson = .true.

        ! setup settings object
        call setup_settings(settings)

        do i = 1, size(subsystem_solvers)
            if (index(subsystem_solvers(i), "davidson") == 0) cycle
            settings%subsystem_solver = subsystem_solvers(i)
            call run_level_shifted_davidson_test("near minimum", subsystem_solvers(i))
            call run_level_shifted_davidson_test("near saddle point", &
                                                 subsystem_solvers(i))
        end do

        contains

        subroutine run_level_shifted_davidson_test(region_str, solver_name)
            !
            ! this subroutine runs a level-shifted Davidson test for a given solver
            !
            use opentrustregion, only: obj_func_type, hess_x_type, &
                                       solver_settings_type, level_shifted_davidson, &
                                       trust_radius_shrink_ratio, &
                                       trust_radius_expand_ratio, &
                                       trust_radius_shrink_factor, &
                                       trust_radius_expand_factor

            character(*), intent(in) :: region_str, solver_name

            integer(ip), parameter :: n_param = 6
            real(rp) :: func, grad_norm, trust_radius, mu, ratio, solution_norm
            real(rp), dimension(n_param) :: grad, h_diag, solution
            integer(ip) :: i, imicro, imicro_jacobi_davidson, error
            procedure(obj_func_type), pointer :: obj_func_funptr
            procedure(hess_x_type), pointer :: hess_x_funptr
            logical :: jacobi_davidson_started, max_precision_reached

            ! initialize variables
            trust_radius = 0.4_rp
            obj_func_funptr => obj_func
            hess_x_funptr => hess_x_fun

            if (region_str == "near minimum") then
                curr_vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
            else if (region_str == "near saddle point") then
                curr_vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
            end if

            func = hartmann6d_func(curr_vars)
            call hartmann6d_gradient(curr_vars, grad)
            grad_norm = norm2(grad)
            call hartmann6d_hessian(curr_vars)
            h_diag = [(hess(i, i), i=1, size(h_diag))]

            ! run level-shifted Davidson
            call level_shifted_davidson(func, grad, grad_norm, h_diag, n_param, &
                                        obj_func_funptr, hess_x_funptr, settings, &
                                        trust_radius, solution, solution_norm, mu, &
                                        imicro, imicro_jacobi_davidson, &
                                        jacobi_davidson_started, &
                                        max_precision_reached, error)

            if (error /= 0) then
                write (stderr, *) "test_level_shifted_davidson failed: Produced "// &
                    "error " // trim(region_str) // " with " // trim(solver_name) // &
                    " solver."
                test_level_shifted_davidson = .false.
            end if
            if (region_str == "near minimum") then
                if (abs(mu) > tol) then
                    write (stderr, *) "test_level_shifted_davidson failed: Level "// &
                        "shift is not zero near minimum with " // trim(solver_name) // &
                        " solver."
                    test_level_shifted_davidson = .false.
                end if
                if (sum(abs(grad + hartmann6d_hess_x(solution))) > &
                    settings%local_red_factor * grad_norm) then
                    write (stderr, *) "test_level_shifted_davidson failed: "// &
                        "Solution does not describe Newton step near minimum with " // &
                        trim(solver_name) // " solver."
                    test_level_shifted_davidson = .false.
                end if
            else if (region_str == "near saddle point") then
                if (mu >= 0.0_rp) then
                    write (stderr, *) "test_level_shifted_davidson failed: Level "// &
                        "shift is not negative near saddle point with " // &
                        trim(solver_name) // " solver."
                    test_level_shifted_davidson = .false.
                end if
                if (sum(abs(grad + hartmann6d_hess_x(solution) - mu * solution)) > &
                    settings%global_red_factor * grad_norm) then
                    write (stderr, *) "test_level_shifted_davidson failed: "// &
                        "Solution does not describe level-shifted Newton step near "// &
                        "saddle point with "// trim(solver_name) // " solver."
                    test_level_shifted_davidson = .false.
                end if
            end if

            ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                    dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))

            if (region_str == "near minimum") then
                if ((ratio < trust_radius_shrink_ratio .and. &
                     solution_norm > trust_radius / trust_radius_shrink_factor) .or. &
                    (ratio > trust_radius_expand_ratio .and. &
                     solution_norm > trust_radius / trust_radius_expand_factor)) then
                    write (stderr, *) "test_level_shifted_davidson failed: "// &
                        "Solution does not stay within trust region near minimum " // &
                        "with " // trim(solver_name) // " solver."
                    test_level_shifted_davidson = .false.
                end if
            else if (region_str == "near saddle point") then
                if ((trust_radius_shrink_ratio > ratio .and. &
                     abs(solution_norm - trust_radius / trust_radius_shrink_factor) > &
                     tol) .or. &
                    (ratio > trust_radius_expand_ratio .and. &
                     abs(solution_norm - trust_radius / trust_radius_expand_factor) > &
                        tol)) then
                    write (stderr, *) "test_level_shifted_davidson failed: "// &
                        "Solution does not lie at trust region boundary near " // &
                        "saddle point with " // trim(solver_name) // " solver."
                    test_level_shifted_davidson = .false.
                end if
            end if

        end subroutine run_level_shifted_davidson_test

    end function test_level_shifted_davidson

    logical(c_bool) function test_truncated_conjugate_gradient() bind(C)
        !
        ! this function tests the truncated conjugate gradient subroutine
        !
        use opentrustregion, only: obj_func_type, hess_x_type, solver_settings_type, &
                                   truncated_conjugate_gradient, &
                                   trust_radius_shrink_ratio, &
                                   trust_radius_expand_ratio, &
                                   trust_radius_shrink_factor, &
                                   trust_radius_expand_factor

        integer(ip), parameter :: n_param = 6
        real(rp) :: func, grad_norm, trust_radius, ratio, solution_norm
        real(rp), dimension(n_param) :: grad, h_diag, solution
        integer(ip) :: i, imicro, error
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: settings
        logical :: max_precision_reached

        ! assume tests pass
        test_truncated_conjugate_gradient = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_micro = 50
        settings%n_random_trial_vectors = 1

        ! initialize variables
        trust_radius = 0.4_rp
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun

        ! start in quadratic region near minimum
        curr_vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]

        ! run truncated conjugate gradient, check whether the solution lies at the 
        ! trust region boundary and reduces the function value
        call truncated_conjugate_gradient(func, grad, grad_norm, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          trust_radius, solution, solution_norm, &
                                          imicro, max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_truncated_jacobi_davidson failed: Produced "// &
                "error near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if (ratio <= 0.0_rp) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near minimum."
            test_truncated_conjugate_gradient = .false.
        end if
        if ((ratio < trust_radius_shrink_ratio .and. solution_norm > trust_radius / &
             trust_radius_shrink_factor) .or. &
            (ratio > trust_radius_expand_ratio .and. solution_norm > trust_radius / &
             trust_radius_expand_factor)) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not stay within trust region near minimum."
            test_truncated_conjugate_gradient = .false.
        end if

        ! run truncated conjugate gradient with a spherical trust region and check that 
        ! the returned solution norm equals the true Euclidean norm of the solution
        settings%trust_region_shape = "spherical"
        call truncated_conjugate_gradient(func, grad, grad_norm, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          trust_radius, solution, solution_norm, &
                                          imicro, max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error near minimum for spherical trust region."
            test_truncated_conjugate_gradient = .false.
        end if
        if (abs(norm2(solution) - solution_norm) > tol) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Returned "// &
                "solution norm does not match Euclidean solution norm for "// &
                "spherical trust region."
            test_truncated_conjugate_gradient = .false.
        end if
        settings%trust_region_shape = "ellipsoidal"

        ! start near saddle point
        curr_vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        trust_radius = 0.4_rp

        ! run truncated conjugate gradient, check whether the solution lies at the 
        ! trust region boundary and reduces the function value
        call truncated_conjugate_gradient(func, grad, grad_norm, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          trust_radius, solution, solution_norm, &
                                          imicro, max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_truncated_jacobi_davidson failed: Produced "// &
                "error near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if (ratio <= 0.0_rp) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not reduce function value near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if
        if ((trust_radius_shrink_ratio > ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_shrink_factor) > tol) &
             .or. &
            (ratio > trust_radius_expand_ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_expand_factor) > tol)) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not lie at trust region boundary near saddle point."
            test_truncated_conjugate_gradient = .false.
        end if

        ! run truncated conjugate gradient with a spherical trust region near the
        ! saddle point and check that the returned solution norm equals the true
        ! Euclidean norm of the solution
        settings%trust_region_shape = "spherical"
        trust_radius = 0.4_rp
        call truncated_conjugate_gradient(func, grad, grad_norm, h_diag, n_param, &
                                          obj_func_funptr, hess_x_funptr, settings, &
                                          trust_radius, solution, solution_norm, &
                                          imicro, max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Produced "// &
                "error near saddle point for spherical trust region."
            test_truncated_conjugate_gradient = .false.
        end if
        if (abs(norm2(solution) - solution_norm) > tol) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Returned "// &
                "solution norm does not match Euclidean solution norm near saddle "// &
                "point for spherical trust region."
            test_truncated_conjugate_gradient = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if ((trust_radius_shrink_ratio > ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_shrink_factor) > tol) &
             .or. &
            (ratio > trust_radius_expand_ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_expand_factor) > tol)) then
            write (stderr, *) "test_truncated_conjugate_gradient failed: Solution "// &
                "does not lie at trust region boundary near saddle point for "// &
                "spherical trust region."
            test_truncated_conjugate_gradient = .false.
        end if
        settings%trust_region_shape = "ellipsoidal"

    end function test_truncated_conjugate_gradient

    logical(c_bool) function test_generalized_lanczos_trust_region() bind(C)
        !
        ! this function tests the generalized Lanczos trust region subroutine
        !
        use opentrustregion, only: obj_func_type, hess_x_type, solver_settings_type, &
                                   generalized_lanczos_trust_region, &
                                   trust_radius_shrink_ratio, &
                                   trust_radius_expand_ratio, &
                                   trust_radius_shrink_factor, &
                                   trust_radius_expand_factor, precond_rel_floor_factor

        integer(ip), parameter :: n_param = 6
        real(rp) :: func, grad_norm, trust_radius, lambda, ratio, solution_norm
        real(rp), dimension(n_param) :: grad, h_diag, solution, residual, precond
        integer(ip) :: i, imicro, error
        procedure(obj_func_type), pointer :: obj_func_funptr
        procedure(hess_x_type), pointer :: hess_x_funptr
        type(solver_settings_type) :: settings
        logical :: max_precision_reached

        ! assume tests pass
        test_generalized_lanczos_trust_region = .true.

        ! setup settings object
        call setup_settings(settings)
        settings%n_random_trial_vectors = 0

        ! initialize variables
        trust_radius = 0.4_rp
        obj_func_funptr => obj_func
        hess_x_funptr => hess_x_fun

        ! start in quadratic region near minimum
        curr_vars = [0.20_rp, 0.15_rp, 0.48_rp, 0.28_rp, 0.31_rp, 0.66_rp]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]

        ! run generalized Lanczos trust region without perturbation, check if error has 
        ! occured, whether the Lagrange multiplier shift vanishes and whether the 
        ! solution stays within trust region and describes the Newton step
        call generalized_lanczos_trust_region(func, grad, grad_norm, h_diag, n_param, &
                                              obj_func_funptr, hess_x_funptr, &
                                              settings, trust_radius, solution, &
                                              solution_norm, lambda, imicro, &
                                              max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Produced error near minimum."
            test_generalized_lanczos_trust_region = .false.
        end if
        if (abs(lambda) > tol) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Lagrange multiplier is not zero near minimum."
            test_generalized_lanczos_trust_region = .false.
        end if
        residual = grad + hartmann6d_hess_x(solution)
        if (sqrt(dot_product(residual, residual / h_diag)) > &
            settings%local_red_factor * grad_norm) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not describe Newton step near minimum."
            test_generalized_lanczos_trust_region = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if ((ratio < trust_radius_shrink_ratio .and. solution_norm > trust_radius &
             / trust_radius_shrink_factor) .or. &
            (ratio > trust_radius_expand_ratio .and. solution_norm > trust_radius &
            / trust_radius_expand_factor)) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not stay within trust region near minimum."
            test_generalized_lanczos_trust_region = .false.
        end if

        ! run generalized Lanczos trust region with perturbation, check if error has 
        ! occured, whether the Lagrange multiplier shift vanishes and whether the 
        ! solution stays within trust region and reduces the function value
        settings%n_random_trial_vectors = 1
        call generalized_lanczos_trust_region(func, grad, grad_norm, h_diag, n_param, &
                                              obj_func_funptr, hess_x_funptr, &
                                              settings, trust_radius, solution, &
                                              solution_norm, lambda, imicro, &
                                              max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Produced error near minimum for perturbed system."
            test_generalized_lanczos_trust_region = .false.
        end if
        if (abs(lambda) > tol) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Lagrange multiplier is not zero near minimum for perturbed system."
            test_generalized_lanczos_trust_region = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if (ratio <= 0.0_rp) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not reduce function value near minimum for "// &
                "perturbed system."
            test_generalized_lanczos_trust_region = .false.
        end if
        if ((ratio < trust_radius_shrink_ratio .and. solution_norm > trust_radius &
             / trust_radius_shrink_factor) .or. &
            (ratio > trust_radius_expand_ratio .and. solution_norm > trust_radius &
            / trust_radius_expand_factor)) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not stay within trust region near minimum for "// &
                "perturbed system."
            test_generalized_lanczos_trust_region = .false.
        end if

        ! run generalized Lanczos trust region with a spherical trust region and check
        ! that the returned solution norm equals the true Euclidean norm of the
        ! solution, unlike the default ellipsoidal trust region whose norm is measured
        ! in the preconditioner metric; GLTR falls back to no preconditioning to
        ! achieve this, so this also implicitly checks that this fallback still
        ! produces a valid, converged solution
        settings%trust_region_shape = "spherical"
        settings%n_random_trial_vectors = 0
        call generalized_lanczos_trust_region(func, grad, grad_norm, h_diag, n_param, &
                                              obj_func_funptr, hess_x_funptr, &
                                              settings, trust_radius, solution, &
                                              solution_norm, lambda, imicro, &
                                              max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Produced error near minimum for spherical trust region."
            test_generalized_lanczos_trust_region = .false.
        end if
        if (abs(norm2(solution) - solution_norm) > tol) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Returned solution norm does not match Euclidean solution norm "// &
                "for spherical trust region."
            test_generalized_lanczos_trust_region = .false.
        end if
        settings%trust_region_shape = "ellipsoidal"

        ! start near saddle point
        curr_vars = [0.35_rp, 0.59_rp, 0.48_rp, 0.40_rp, 0.31_rp, 0.32_rp]
        func = hartmann6d_func(curr_vars)
        call hartmann6d_gradient(curr_vars, grad)
        grad_norm = norm2(grad)
        call hartmann6d_hessian(curr_vars)
        h_diag = [(hess(i, i), i=1, size(h_diag))]
        trust_radius = 0.4_rp

        ! run generalized Lanczos trust region without perturbation, check if error has
        ! occured, whether the Lagrange multiplier is positive and whether the solution
        ! lies at the trust region boundary and describes a level-shifted Newton step
        settings%n_random_trial_vectors = 0
        call generalized_lanczos_trust_region(func, grad, grad_norm, h_diag, n_param, &
                                              obj_func_funptr, hess_x_funptr, &
                                              settings, trust_radius, solution, &
                                              solution_norm, lambda, imicro, &
                                              max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Produced error near saddle point."
            test_generalized_lanczos_trust_region = .false.
        end if
        if (lambda <= 0.0_rp) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Lagrange multiplier is not positive near saddle point."
            test_generalized_lanczos_trust_region = .false.
        end if
        precond = max(abs(h_diag), precond_rel_floor_factor * maxval(abs(h_diag)))
        residual = grad + hartmann6d_hess_x(solution) + lambda * solution * precond
        if (sqrt(dot_product(residual, residual / precond)) > &
            settings%global_red_factor * grad_norm) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not describe level-shifted Newton step near saddle "// &
                "point."
            test_generalized_lanczos_trust_region = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if ((trust_radius_shrink_ratio > ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_shrink_factor) > tol) &
             .or. &
            (ratio > trust_radius_expand_ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_expand_factor) > tol)) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not lie at trust region boundary near saddle point."
            test_generalized_lanczos_trust_region = .false.
        end if

        ! run generalized Lanczos trust region with a spherical trust region near the
        ! saddle point and check that the returned solution norm equals the true
        ! Euclidean norm of the solution and lies at the trust region boundary
        settings%trust_region_shape = "spherical"
        trust_radius = 0.4_rp
        call generalized_lanczos_trust_region(func, grad, grad_norm, h_diag, n_param, &
                                              obj_func_funptr, hess_x_funptr, &
                                              settings, trust_radius, solution, &
                                              solution_norm, lambda, imicro, &
                                              max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Produced error near saddle point for spherical trust region."
            test_generalized_lanczos_trust_region = .false.
        end if
        if (abs(norm2(solution) - solution_norm) > tol) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Returned solution norm does not match Euclidean solution norm "// &
                "near saddle point for spherical trust region."
            test_generalized_lanczos_trust_region = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if ((trust_radius_shrink_ratio > ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_shrink_factor) > tol) &
             .or. &
            (ratio > trust_radius_expand_ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_expand_factor) > tol)) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not lie at trust region boundary near saddle point "// &
                "for spherical trust region."
            test_generalized_lanczos_trust_region = .false.
        end if
        settings%trust_region_shape = "ellipsoidal"

        ! run generalized Lanczos trust region with perturbation, check if error has
        ! occured, whether the Lagrange multiplier is positive and whether the solution
        ! lies at the trust region boundary and reduces the function value
        settings%n_random_trial_vectors = 1
        call generalized_lanczos_trust_region(func, grad, grad_norm, h_diag, n_param, &
                                              obj_func_funptr, hess_x_funptr, &
                                              settings, trust_radius, solution, &
                                              solution_norm, lambda, imicro, &
                                              max_precision_reached, error)
        if (error /= 0) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Produced error near saddle point for perturbed system."
            test_generalized_lanczos_trust_region = .false.
        end if
        if (lambda <= 0.0_rp) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Lagrange multiplier is not positive near saddle point for "// &
                "perturbed system."
            test_generalized_lanczos_trust_region = .false.
        end if
        ratio = (hartmann6d_func(curr_vars + solution) - func) / &
                dot_product(solution, grad + 0.5_rp * hartmann6d_hess_x(solution))
        if (ratio <= 0.0_rp) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not reduce function value near saddle point for "// &
                "perturbed system."
            test_generalized_lanczos_trust_region = .false.
        end if
        if ((trust_radius_shrink_ratio > ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_shrink_factor) > tol) &
             .or. &
            (ratio > trust_radius_expand_ratio .and. &
             abs(solution_norm - trust_radius / trust_radius_expand_factor) > tol)) then
            write (stderr, *) "test_generalized_lanczos_trust_region failed: "// &
                "Solution does not lie at trust region boundary near saddle point "// &
                "for perturbed system."
            test_generalized_lanczos_trust_region = .false.
        end if

    end function test_generalized_lanczos_trust_region

    logical(c_bool) function test_add_error_origin() bind(C)
        !
        ! this function tests the subroutine that adds the error origin to an error 
        ! code
        !
        use opentrustregion, only: solver_settings_type, add_error_origin

        type(solver_settings_type) :: settings
        integer(ip) :: error

        ! assume tests pass
        test_add_error_origin = .true.

        ! setup settings object
        call setup_settings(settings)

        ! check if subroutine adds error origin correctly if no origin is present
        error = 1
        call add_error_origin(error, 100_ip, settings)
        if (error /= 101) then
            write (stderr, *) "test_add_error_origin failed: Error origin not "// &
                "correctly added."
            test_add_error_origin = .false.
        end if

        ! check if subroutine skips adding error origin if origin is already present
        call add_error_origin(error, 100_ip, settings)
        if (error /= 101) then
            write (stderr, *) "test_add_error_origin failed: Error code modified "// &
                "even though error origin is already present."
            test_add_error_origin = .false.
        end if

        ! check if subroutine does not modify error code when no error is encountered
        error = 0
        call add_error_origin(error, 100_ip, settings)
        if (error /= 0) then
            write (stderr, *) "test_add_error_origin failed: Error code modified "// &
                "even though error code of zero was passed."
            test_add_error_origin = .false.
        end if

        ! check if subroutine raises error for invalid error code
        error = -1
        call add_error_origin(error, 100_ip, settings)
        if (error /= 101) then
            write (stderr, *) "test_add_error_origin failed: Error code not "// &
                "correctly returned for invalid (negative) error code."
            test_add_error_origin = .false.
        end if

    end function test_add_error_origin

    logical(c_bool) function test_string_to_lowercase() bind(C)
        !
        ! this function tests the function that transfers strings to lowercase
        !
        use opentrustregion, only: string_to_lowercase

        character(*), parameter :: input  = 'OpenTrustRegion123!', &
                                   expect = 'opentrustregion123!'

        ! assume tests pass
        test_string_to_lowercase = .true.

        ! test transfer to lowercase
        if (string_to_lowercase(input) /= expect) then
            write (stderr, *) "test_string_to_lowercase failed: String not "// &
                "correctly transferred to lowercase."
            test_string_to_lowercase = .false.
        end if

    end function test_string_to_lowercase

    logical(c_bool) function test_string_in() bind(C)
        !
        ! this function tests the function that checks whether a string is contained in 
        ! another string
        !
        use opentrustregion, only: string_in

        character(*), parameter :: input  = 'OpenTrustRegion123!'

        ! assume tests pass
        test_string_in = .true.

        ! test if substring is in string
        if (.not. string_in("Trust", input)) then
            write (stderr, *) "test_string_in failed: Substring not correctly "// &
                "detected in string."
            test_string_in = .false.
        end if
        if (string_in("trust", input)) then
            write (stderr, *) "test_string_in failed: Substring detected in string "// &
                "even though it is not present."
            test_string_in = .false.
        end if

    end function test_string_in

end module opentrustregion_unit_tests
