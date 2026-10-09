! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip, update_orbs_c_type, hess_x_c_type, &
                           obj_func_c_type, precond_c_type, precond_pd_c_type, &
                           project_c_type, conv_check_c_type, logger_c_type
    use test_reference, only: tol, tol_c, n_param, n_param_c, check_host_context_c, &
                              arm_host_context_c, host_context_reached
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_loc, c_funptr, c_funloc, &
                                           c_char, c_associated, c_null_ptr, &
                                           c_null_char, c_null_funptr

    implicit none

    ! logical to test logging function
    logical :: test_logger

    ! a target of a type no callback bundle can hold, used to drive the guard the
    ! wrappers raise when they are handed a context they did not create
    real(rp), target :: foreign_context_target = 0.0_rp

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_c_type), pointer :: mock_update_orbs_ptr => mock_update_orbs
    procedure(update_orbs_c_type), pointer :: mock_update_orbs_no_hess_x_ptr => &
        mock_update_orbs_no_hess_x
    procedure(hess_x_c_type), pointer :: mock_hess_x_ptr => mock_hess_x
    procedure(obj_func_c_type), pointer :: mock_obj_func_ptr => mock_obj_func
    procedure(precond_c_type), pointer :: mock_precond_ptr => mock_precond
    procedure(precond_pd_c_type), pointer :: mock_precond_pd_ptr => mock_precond_pd
    procedure(project_c_type), pointer :: mock_project_ptr => mock_project
    procedure(conv_check_c_type), pointer :: mock_conv_check_ptr => mock_conv_check
    procedure(logger_c_type), pointer :: mock_logger_ptr => mock_logger

contains

    function mock_update_orbs(kappa, func, grad, h_diag, hess_x_c_funptr, context_c) &
        result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the orbital update C function
        !
        real(c_rp), intent(in), target :: kappa(*)
        real(c_rp), intent(out) :: func
        real(c_rp), intent(out), target :: grad(*), h_diag(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        func = sum(kappa(:n_param))

        grad(:n_param) = 2 * kappa(:n_param)

        h_diag(:n_param) = 3 * kappa(:n_param)

        hess_x_c_funptr = c_funloc(mock_hess_x)

        error = 0_c_ip

    end function mock_update_orbs

    function mock_update_orbs_no_hess_x(kappa, func, grad, h_diag, hess_x_c_funptr, &
                                        context_c) result(error) bind(C)
        !
        ! this function is a test function for an orbital update C function which
        ! reports success but does not provide a Hessian linear transformation
        !
        real(c_rp), intent(in), target :: kappa(*)
        real(c_rp), intent(out) :: func
        real(c_rp), intent(out), target :: grad(*), h_diag(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        error = mock_update_orbs(kappa, func, grad, h_diag, hess_x_c_funptr, context_c)
        hess_x_c_funptr = c_null_funptr

    end function mock_update_orbs_no_hess_x

    function mock_hess_x(x, hess_x, context_c) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the Hessian linear transformation C
        ! function
        !
        real(c_rp), intent(in), target :: x(*)
        real(c_rp), intent(out), target :: hess_x(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        hess_x(:n_param) = 4 * x(:n_param)

        error = 0

    end function mock_hess_x

    function mock_obj_func(kappa, func, context_c) result(error) bind(C)
        !
        ! this function is a test function for the C objective function
        !
        real(c_rp), intent(in), target :: kappa(*)
        real(c_rp), intent(out) :: func
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        func = sum(kappa(:n_param))

        error = 0

    end function mock_obj_func

    function mock_precond(residual, mu, precond_residual, context_c) result(error) &
        bind(C)
        !
        ! this function is a test function for the C preconditioner function
        !
        real(c_rp), intent(in), target :: residual(*)
        real(c_rp), intent(in) :: mu
        real(c_rp), intent(out), target :: precond_residual(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        precond_residual(:n_param) = mu * residual(:n_param)

        error = 0

    end function mock_precond

    function mock_precond_pd(residual, precond_residual, context_c) result(error) &
        bind(C)
        !
        ! this function is a test function for the C positive-definite preconditioner
        ! function
        !
        real(c_rp), intent(in), target :: residual(*)
        real(c_rp), intent(out), target :: precond_residual(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        precond_residual(:n_param) = 3 * residual(:n_param)

        error = 0

    end function mock_precond_pd

    function mock_project(vector, context_c) result(error) bind(C)
        !
        ! this function is a test function for the C projection function
        !
        real(c_rp), intent(inout), target :: vector(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        vector(:n_param) = 2 * vector(:n_param)

        error = 0

    end function mock_project

    function mock_conv_check(converged, context_c) result(error) bind(C)
        !
        ! this function is a test function for the convergence check function
        !
        logical(c_bool), intent(out) :: converged
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        converged = .true.

        error = 0

    end function mock_conv_check

    function mock_get_extra_trial_vectors(trial_vectors, n_extra_trial_vectors_c, &
                                          context_c) result(error) bind(C)
        !
        ! this function is a test function for the extra trial vector function
        !
        real(c_rp), intent(out), target :: trial_vectors(*)
        integer(c_ip), intent(in), value :: n_extra_trial_vectors_c
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: i, error

        ! check host context
        call check_host_context_c(context_c)

        do i = 1, n_extra_trial_vectors_c
            trial_vectors((i - 1) * n_param + 1:i * n_param) = real(i, kind=c_rp)
        end do

        error = 0

    end function mock_get_extra_trial_vectors

    function mock_conv_check_stability(residual, eigval, converged, context_c) &
        result(error) bind(C)
        !
        ! this function is a test function for the stability check convergence check
        ! function
        !
        real(c_rp), intent(in), target :: residual(*)
        real(c_rp), intent(in) :: eigval
        logical(c_bool), intent(out) :: converged
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        ! check host context
        call check_host_context_c(context_c)

        converged = abs(sum(residual(:n_param)) - eigval) < tol_c

        error = 0

    end function mock_conv_check_stability

    subroutine mock_logger(message_c, context_c) bind(C)
        !
        ! this function is a test function for the C logging function
        !
        character(kind=c_char), intent(in) :: message_c(*)
        type(c_ptr), intent(in), value :: context_c
        character(len=4) :: message

        ! check host context
        call check_host_context_c(context_c)

        message = transfer(message_c(1:4), message)
        if (message == "test") test_logger = .true.

    end subroutine mock_logger

    subroutine solver_set_mock(settings)
        !
        ! this function sets optional mock routines for the C solver settings object
        !
        use c_interface, only: solver_settings_type_c

        type(solver_settings_type_c) :: settings

        settings%precond = c_funloc(mock_precond)
        settings%precond_pd = c_funloc(mock_precond_pd)
        settings%project = c_funloc(mock_project)
        settings%get_extra_trial_vectors = c_funloc(mock_get_extra_trial_vectors)
        settings%conv_check = c_funloc(mock_conv_check)
        settings%stability_hess_x = c_funloc(mock_hess_x)
        settings%logger = c_funloc(mock_logger)
        call stability_set_mock(settings%stability_settings)

    end subroutine solver_set_mock

    subroutine stability_set_mock(settings)
        !
        ! this function sets optional mock routines for the C stability settings object
        !
        use c_interface, only: stability_settings_type_c

        type(stability_settings_type_c) :: settings

        settings%precond = c_funloc(mock_precond)
        settings%project = c_funloc(mock_project)
        settings%approx_hess_x = c_funloc(mock_hess_x)
        settings%get_extra_trial_vectors = c_funloc(mock_get_extra_trial_vectors)
        settings%conv_check = c_funloc(mock_conv_check_stability)
        settings%logger = c_funloc(mock_logger)

    end subroutine stability_set_mock

    subroutine solver_mock_callbacks(callbacks)
        !
        ! this subroutine bundles the mock routines that solver_set_mock sets, so that
        ! the Fortran wrappers of converted settings can call them through the context
        !
        use c_interface, only: c_callbacks_type

        type(c_callbacks_type), intent(inout) :: callbacks

        callbacks%precond => mock_precond
        callbacks%precond_pd => mock_precond_pd
        callbacks%project => mock_project
        callbacks%get_extra_trial_vectors => mock_get_extra_trial_vectors
        callbacks%conv_check => mock_conv_check
        callbacks%stability_hess_x => mock_hess_x
        callbacks%logger => mock_logger

    end subroutine solver_mock_callbacks

    subroutine stability_mock_callbacks(callbacks)
        !
        ! this subroutine bundles the mock routines that stability_set_mock sets, so
        ! that the Fortran wrappers of converted settings can call them through the
        ! context
        !
        use c_interface, only: c_callbacks_type

        type(c_callbacks_type), intent(inout) :: callbacks

        callbacks%precond => mock_precond
        callbacks%project => mock_project
        callbacks%approx_hess_x => mock_hess_x
        callbacks%get_extra_trial_vectors => mock_get_extra_trial_vectors
        callbacks%conv_check_stability => mock_conv_check_stability
        callbacks%logger => mock_logger

    end subroutine stability_mock_callbacks

    logical(c_bool) function test_solver_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the solver
        !
        use c_interface, only: solver_settings_type_c, solver, solver_c_wrapper
        use opentrustregion, only: standard_solver => solver
        use opentrustregion_mock, only: mock_solver, test_passed, mock_n_update_orbs, &
                                        mock_n_hess_x, mock_stability_n_hess_x
        use test_reference, only: assignment(=), ref_settings, stability_host_context

        type(c_funptr) :: update_orbs_c_funptr, obj_func_c_funptr
        type(solver_settings_type_c) :: settings
        integer(c_ip) :: error
        integer(ip) :: icase
        character(len=22), parameter :: case_names(2) = &
            [character(len=22) :: "without nested context", "with nested context"]

        ! assume tests pass
        test_solver_c_wrapper = .true.

        ! inject mock function
        solver => mock_solver

        ! get C function pointers to Fortran functions
        update_orbs_c_funptr = c_funloc(mock_update_orbs)
        obj_func_c_funptr = c_funloc(mock_obj_func)

        ! run once without and once with a context on the nested stability check
        ! settings, the callback functions of the internal stability check have to
        ! receive the solver's context in the former and the nested one in the latter
        do icase = 1, size(case_names)
            ! associate optional settings with values
            settings = ref_settings
            call solver_set_mock(settings)

            ! set host contexts
            call arm_host_context_c(settings%context)
            if (icase == 2) &
                settings%stability_settings%context = c_loc(stability_host_context)

            ! initialize logger logical
            test_logger = .false.

            ! call solver
            error = solver_c_wrapper(update_orbs_c_funptr, obj_func_c_funptr, &
                                     n_param_c, settings)

            ! check if logging subroutine was correctly called
            if (.not. test_logger) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Called logging "// &
                    "subroutine wrong "//trim(case_names(icase))//"."
            end if

            ! check if output variables are as expected
            if (error /= 0) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned error "// &
                    "boolean wrong "//trim(case_names(icase))//"."
            end if

            ! check if output fields are written back with the values set by the solver
            if (settings%max_precision_reached) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned maximum "// &
                    "precision reached flag wrong "//trim(case_names(icase))//"."
            end if
            if (settings%n_update_orbs /= int(mock_n_update_orbs, kind=c_ip)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned number of "// &
                    "orbital updates wrong "//trim(case_names(icase))//"."
            end if
            if (settings%n_hess_x /= int(mock_n_hess_x, kind=c_ip)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned number of "// &
                    "Hessian linear transformations wrong "//trim(case_names(icase))// &
                    "."
            end if
            if (settings%stability_settings%n_hess_x /= &
                int(mock_stability_n_hess_x, kind=c_ip)) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Returned number of "// &
                    "Hessian linear transformations of internal stability check "// &
                    "wrong "//trim(case_names(icase))//"."
            end if

            ! check that the callback functions of the internal stability check
            ! received the context of the nested settings when one was set
            if (icase == 2 .and. stability_host_context%n_calls == 0) then
                test_solver_c_wrapper = .false.
                write(stderr, *) "test_solver_c_wrapper failed: Callback functions "// &
                    "of internal stability check did not receive context of nested "// &
                    "settings."
            end if

            ! check that the host context reached the callback functions unchanged
            test_solver_c_wrapper = test_solver_c_wrapper .and. &
                                    host_context_reached("solver_c_wrapper")

            ! check if test has passed
            test_solver_c_wrapper = test_solver_c_wrapper .and. test_passed
        end do

        ! restore the procedure pointer so later tests do not inherit the mock
        solver => standard_solver

    end function test_solver_c_wrapper

    logical(c_bool) function test_stability_check_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the stability check
        !
        use c_interface, only: stability_settings_type_c, stability_check, &
                               stability_check_c_wrapper
        use opentrustregion, only: standard_stability_check => stability_check
        use opentrustregion_mock, only: mock_stability_check, test_passed
        use test_reference, only: assignment(=), ref_settings

        type(c_funptr) :: hess_x_c_funptr
        real(c_rp), allocatable :: h_diag(:)
        real(c_rp), allocatable, target :: kappa(:)
        real(c_rp), target :: min_eigval
        type(stability_settings_type_c) :: settings
        logical(c_bool) :: stable
        type(c_ptr) :: kappa_c_ptr, min_eigval_c_ptr
        integer(c_ip) :: error

        ! assume tests pass
        test_stability_check_c_wrapper = .true.

        ! inject mock function
        stability_check => mock_stability_check

        ! get C function pointers to Fortran functions
        hess_x_c_funptr = c_funloc(mock_hess_x)

        ! initialize Hessian diagonal and get C pointers
        allocate(h_diag(n_param))
        h_diag = 3.0_c_rp

        ! associate optional arguments with values
        settings = ref_settings
        call stability_set_mock(settings)

        ! set host context
        call arm_host_context_c(settings%context)

        ! unassociate returned direction pointer
        kappa_c_ptr = c_null_ptr
        min_eigval_c_ptr = c_null_ptr

        ! call stability check first without initialized returned direction
        error = stability_check_c_wrapper(h_diag, hess_x_c_funptr, n_param_c, stable, &
                                          settings, kappa_c_ptr, min_eigval_c_ptr)

        ! check that the host context reached the callback functions unchanged
        test_stability_check_c_wrapper = &
            test_stability_check_c_wrapper .and. &
            host_context_reached("stability_check_c_wrapper")

        ! check if test has passed
        test_stability_check_c_wrapper = test_stability_check_c_wrapper .and. &
                                         test_passed

        ! check if logging subroutine was correctly called
        if (.not. test_logger) then
            test_stability_check_c_wrapper = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Called "// &
                "logging subroutine wrong."
        end if

        ! check if output variables are as expected
        if (stable) then
            test_stability_check_c_wrapper = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "stability boolean wrong."
        end if

        if (error /= 0) then
            test_stability_check_c_wrapper = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "error code wrong."
        end if

        ! associate returned direction with value
        allocate(kappa(n_param))
        kappa_c_ptr = c_loc(kappa)
        min_eigval_c_ptr = c_loc(min_eigval)

        ! initialize logger logical
        test_logger = .false.

        ! call stability check with initilized returned direction
        error = stability_check_c_wrapper(h_diag, hess_x_c_funptr, n_param_c, stable, &
                                          settings, kappa_c_ptr, min_eigval_c_ptr)

        ! check if logging subroutine was correctly called
        if (.not. test_logger) then
            test_stability_check_c_wrapper = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Called "// &
                "logging subroutine wrong."
        end if

        ! check if output variables are as expected
        if (stable) then
            test_stability_check_c_wrapper = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "stability boolean wrong."
        end if

        if (any(abs(kappa - 1.0_c_rp) > tol_c)) then
            test_stability_check_c_wrapper = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "direction wrong."
        end if

        if (abs(min_eigval - (-1.0_c_rp)) > tol_c) then
            test_stability_check_c_wrapper = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "minimum eigenvalue wrong."
        end if

        if (error /= 0) then
            test_stability_check_c_wrapper = .false.
            write(stderr, *) "test_stability_check_c_wrapper failed: Returned "// &
                "error code wrong."
        end if
        deallocate(h_diag, kappa)

        ! check if test has passed
        test_stability_check_c_wrapper = test_passed .and. &
                                         test_stability_check_c_wrapper

        ! restore the procedure pointer so later tests do not inherit the mock
        stability_check => standard_stability_check

    end function test_stability_check_c_wrapper

    logical(c_bool) function test_store_optional_c_callbacks() bind(C)
        !
        ! this function tests the subroutine that stores the optional C callback
        ! functions and the host context of C settings in a callback bundle
        !
        use c_interface, only: c_callbacks_type, store_optional_c_callbacks
        use test_reference, only: host_context

        type(c_callbacks_type) :: callbacks

        ! assume test passes
        test_store_optional_c_callbacks = .true.

        ! store the callback functions and host context of initialized settings
        call store_optional_c_callbacks( &
            callbacks, .true._c_bool, c_funloc(mock_precond), &
            c_funloc(mock_precond_pd), c_funloc(mock_project), &
            c_funloc(mock_get_extra_trial_vectors), c_funloc(mock_conv_check), &
            c_funloc(mock_conv_check_stability), c_funloc(mock_hess_x), &
            c_funloc(mock_hess_x), c_funloc(mock_logger), c_loc(host_context))
        if (.not. all_stored()) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Callback "// &
                "functions of initialized settings not stored."
        end if
        if (.not. c_associated(callbacks%host_context, c_loc(host_context))) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Host context "// &
                "of initialized settings not stored."
        end if

        ! store initialized settings that provide nothing, the callback functions and
        ! host context already in the bundle have to be kept
        call store_optional_c_callbacks(callbacks, .true._c_bool, c_null_funptr, &
                                        c_null_funptr, c_null_funptr, c_null_funptr, &
                                        c_null_funptr, c_null_funptr, c_null_funptr, &
                                        c_null_funptr, c_null_funptr, c_null_ptr)
        if (.not. all_stored()) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Callback "// &
                "functions in the bundle replaced by ones that were not provided."
        end if
        if (.not. c_associated(callbacks%host_context, c_loc(host_context))) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Host context "// &
                "in the bundle replaced by one that was not provided."
        end if

        ! store the callback functions and host context of settings that were not
        ! initialized in an empty bundle, nothing may be stored
        callbacks = c_callbacks_type()
        call store_optional_c_callbacks( &
            callbacks, .false._c_bool, c_funloc(mock_precond), &
            c_funloc(mock_precond_pd), c_funloc(mock_project), &
            c_funloc(mock_get_extra_trial_vectors), c_funloc(mock_conv_check), &
            c_funloc(mock_conv_check_stability), c_funloc(mock_hess_x), &
            c_funloc(mock_hess_x), c_funloc(mock_logger), c_loc(host_context))
        if (associated(callbacks%precond) .or. associated(callbacks%precond_pd) .or. &
            associated(callbacks%project) .or. &
            associated(callbacks%get_extra_trial_vectors) .or. &
            associated(callbacks%conv_check) .or. &
            associated(callbacks%conv_check_stability) .or. &
            associated(callbacks%stability_hess_x) .or. &
            associated(callbacks%approx_hess_x) .or. associated(callbacks%logger)) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Callback "// &
                "functions of settings that were not initialized stored."
        end if
        if (c_associated(callbacks%host_context)) then
            test_store_optional_c_callbacks = .false.
            write(stderr, *) "test_store_optional_c_callbacks failed: Host context "// &
                "of settings that were not initialized stored."
        end if

    contains

        logical function all_stored()
            !
            ! this function checks that every optional callback function slot of the
            ! bundle holds its mock
            !
            all_stored = associated(callbacks%precond, mock_precond) .and. &
                         associated(callbacks%precond_pd, mock_precond_pd) .and. &
                         associated(callbacks%project, mock_project) .and. &
                         associated(callbacks%get_extra_trial_vectors, &
                                    mock_get_extra_trial_vectors) .and. &
                         associated(callbacks%conv_check, mock_conv_check) .and. &
                         associated(callbacks%conv_check_stability, &
                                    mock_conv_check_stability) .and. &
                         associated(callbacks%stability_hess_x, mock_hess_x) .and. &
                         associated(callbacks%approx_hess_x, mock_hess_x) .and. &
                         associated(callbacks%logger, mock_logger)

        end function all_stored

    end function test_store_optional_c_callbacks

    logical(c_bool) function test_update_orbs_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the orbital update
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use c_interface, only: c_callbacks_type, update_orbs_f_wrapper
        use test_reference, only: test_update_orbs_funptr

        procedure(update_orbs_type), pointer :: update_orbs_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        procedure(hess_x_type), pointer :: hess_x_funptr
        real(rp) :: kappa(n_param), func, grad(n_param), h_diag(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%update_orbs => mock_update_orbs
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        update_orbs_funptr => update_orbs_f_wrapper

        ! test orbital update wrapper
        test_update_orbs_f_wrapper = test_update_orbs_funptr( &
            update_orbs_funptr, "update_orbs_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_update_orbs_f_wrapper = test_update_orbs_f_wrapper .and. logical( &
            host_context_reached("update_orbs_f_wrapper"), kind=c_bool)

        ! an orbital update that succeeds without providing a Hessian linear
        ! transformation is an error
        callbacks%update_orbs => mock_update_orbs_no_hess_x
        kappa = 1.0_rp
        call update_orbs_f_wrapper(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                   context)
        if (error /= 1) then
            test_update_orbs_f_wrapper = .false.
            write(stderr, *) "test_update_orbs_f_wrapper failed: Did not report a "// &
                "missing Hessian linear transformation."
        end if

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        kappa = 1.0_rp
        call update_orbs_f_wrapper(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                   foreign_context)
        if (error /= 1) then
            test_update_orbs_f_wrapper = .false.
            write(stderr, *) "test_update_orbs_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_update_orbs_f_wrapper

    logical(c_bool) function test_hess_x_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the Hessian linear transformation
        !
        use opentrustregion, only: hess_x_type
        use c_interface, only: c_callbacks_type, hess_x_f_wrapper
        use test_reference, only: test_hess_x_funptr

        procedure(hess_x_type), pointer :: hess_x_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%hess_x => mock_hess_x
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        hess_x_funptr => hess_x_f_wrapper

        ! test Hessian linear transformation wrapper
        test_hess_x_f_wrapper = &
            test_hess_x_funptr(hess_x_funptr, "hess_x_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_hess_x_f_wrapper = test_hess_x_f_wrapper .and. logical( &
            host_context_reached("hess_x_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        x = 1.0_rp
        foreign_context => foreign_context_target
        call hess_x_f_wrapper(x, hess_x, error, foreign_context)
        if (error /= 1) then
            test_hess_x_f_wrapper = .false.
            write(stderr, *) "test_hess_x_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_hess_x_f_wrapper

    logical(c_bool) function test_obj_func_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the objective function
        !
        use opentrustregion, only: obj_func_type
        use c_interface, only: c_callbacks_type, obj_func_f_wrapper
        use test_reference, only: test_obj_func_funptr

        procedure(obj_func_type), pointer :: obj_func_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: kappa(n_param), func
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%obj_func => mock_obj_func
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        obj_func_funptr => obj_func_f_wrapper

        ! test objective function wrapper
        test_obj_func_f_wrapper = &
            test_obj_func_funptr(obj_func_funptr, "obj_func_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_obj_func_f_wrapper = test_obj_func_f_wrapper .and. logical( &
            host_context_reached("obj_func_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        kappa = 1.0_rp
        func = obj_func_f_wrapper(kappa, error, foreign_context)
        if (error /= 1) then
            test_obj_func_f_wrapper = .false.
            write(stderr, *) "test_obj_func_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_obj_func_f_wrapper

    logical(c_bool) function test_precond_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the preconditioner function
        !
        use opentrustregion, only: precond_type
        use c_interface, only: c_callbacks_type, precond_f_wrapper
        use test_reference, only: test_precond_funptr

        procedure(precond_type), pointer :: precond_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: residual(n_param), precond_residual(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%precond => mock_precond
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        precond_funptr => precond_f_wrapper

        ! test preconditioner wrapper
        test_precond_f_wrapper = &
            test_precond_funptr(precond_funptr, "precond_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_precond_f_wrapper = test_precond_f_wrapper .and. logical( &
            host_context_reached("precond_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        residual = 1.0_rp
        call precond_f_wrapper(residual, 1.0_rp, precond_residual, error, &
                               foreign_context)
        if (error /= 1) then
            test_precond_f_wrapper = .false.
            write(stderr, *) "test_precond_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_precond_f_wrapper

    logical(c_bool) function test_precond_pd_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the positive-definite
        ! preconditioner function
        !
        use opentrustregion, only: precond_pd_type
        use c_interface, only: c_callbacks_type, precond_pd_f_wrapper
        use test_reference, only: test_precond_pd_funptr

        procedure(precond_pd_type), pointer :: funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: residual(n_param), precond_residual(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%precond_pd => mock_precond_pd
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        funptr => precond_pd_f_wrapper

        ! test positive-definite preconditioner wrapper
        test_precond_pd_f_wrapper = &
            test_precond_pd_funptr(funptr, "precond_pd_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_precond_pd_f_wrapper = test_precond_pd_f_wrapper .and. logical( &
            host_context_reached("precond_pd_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        residual = 1.0_rp
        call precond_pd_f_wrapper(residual, precond_residual, error, foreign_context)
        if (error /= 1) then
            test_precond_pd_f_wrapper = .false.
            write(stderr, *) "test_precond_pd_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_precond_pd_f_wrapper

    logical(c_bool) function test_project_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the projection function
        !
        use opentrustregion, only: project_type
        use c_interface, only: c_callbacks_type, project_f_wrapper
        use test_reference, only: test_project_funptr

        procedure(project_type), pointer :: project_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: vector(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%project => mock_project
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        project_funptr => project_f_wrapper

        ! test projection wrapper
        test_project_f_wrapper = &
            test_project_funptr(project_funptr, "project_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_project_f_wrapper = test_project_f_wrapper .and. logical( &
            host_context_reached("project_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        vector = 1.0_rp
        call project_f_wrapper(vector, error, foreign_context)
        if (error /= 1) then
            test_project_f_wrapper = .false.
            write(stderr, *) "test_project_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_project_f_wrapper

    logical(c_bool) function test_conv_check_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the convergence check function
        !
        use opentrustregion, only: conv_check_type
        use c_interface, only: c_callbacks_type, conv_check_f_wrapper
        use test_reference, only: test_conv_check_funptr

        procedure(conv_check_type), pointer :: conv_check_funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        logical :: converged
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%conv_check => mock_conv_check
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        conv_check_funptr => conv_check_f_wrapper

        ! test convergence check wrapper
        test_conv_check_f_wrapper = test_conv_check_funptr( &
            conv_check_funptr, "conv_check_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_conv_check_f_wrapper = test_conv_check_f_wrapper .and. logical( &
            host_context_reached("conv_check_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        converged = conv_check_f_wrapper(error, foreign_context)
        if (error /= 1) then
            test_conv_check_f_wrapper = .false.
            write(stderr, *) "test_conv_check_f_wrapper failed: Did not report an "// &
                "invalid context."
        end if

    end function test_conv_check_f_wrapper

    logical(c_bool) function test_stability_hess_x_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the stability check Hessian
        ! linear transformation function
        !
        use opentrustregion, only: hess_x_type
        use c_interface, only: c_callbacks_type, stability_hess_x_f_wrapper
        use test_reference, only: test_hess_x_funptr

        procedure(hess_x_type), pointer :: funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%stability_hess_x => mock_hess_x
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        funptr => stability_hess_x_f_wrapper

        ! test stability check Hessian linear transformation wrapper
        test_stability_hess_x_f_wrapper = &
            test_hess_x_funptr(funptr, "stability_hess_x_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_stability_hess_x_f_wrapper = &
            test_stability_hess_x_f_wrapper .and. &
            logical(host_context_reached("stability_hess_x_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        x = 1.0_rp
        call stability_hess_x_f_wrapper(x, hess_x, error, foreign_context)
        if (error /= 1) then
            test_stability_hess_x_f_wrapper = .false.
            write(stderr, *) "test_stability_hess_x_f_wrapper failed: Did not "// &
                "report an invalid context."
        end if

    end function test_stability_hess_x_f_wrapper

    logical(c_bool) function test_approx_hess_x_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the approximate Hessian linear
        ! transformation function
        !
        use opentrustregion, only: hess_x_type
        use c_interface, only: c_callbacks_type, approx_hess_x_f_wrapper
        use test_reference, only: test_hess_x_funptr

        procedure(hess_x_type), pointer :: funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%approx_hess_x => mock_hess_x
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        funptr => approx_hess_x_f_wrapper

        ! test approximate Hessian linear transformation wrapper
        test_approx_hess_x_f_wrapper = &
            test_hess_x_funptr(funptr, "approx_hess_x_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_approx_hess_x_f_wrapper = test_approx_hess_x_f_wrapper .and. logical( &
            host_context_reached("approx_hess_x_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        x = 1.0_rp
        call approx_hess_x_f_wrapper(x, hess_x, error, foreign_context)
        if (error /= 1) then
            test_approx_hess_x_f_wrapper = .false.
            write(stderr, *) "test_approx_hess_x_f_wrapper failed: Did not report "// &
                "an invalid context."
        end if

    end function test_approx_hess_x_f_wrapper

    logical(c_bool) function test_get_extra_trial_vectors_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the extra trial vector
        ! function
        !
        use opentrustregion, only: get_extra_trial_vectors_type
        use c_interface, only: c_callbacks_type, get_extra_trial_vectors_f_wrapper
        use test_reference, only: test_get_extra_trial_vectors_funptr

        procedure(get_extra_trial_vectors_type), pointer :: funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: trial_vectors(n_param, 1)
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%get_extra_trial_vectors => mock_get_extra_trial_vectors
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        funptr => get_extra_trial_vectors_f_wrapper

        ! test extra trial vector wrapper
        test_get_extra_trial_vectors_f_wrapper = test_get_extra_trial_vectors_funptr( &
            funptr, "get_extra_trial_vectors_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_get_extra_trial_vectors_f_wrapper = &
            test_get_extra_trial_vectors_f_wrapper .and. logical( &
                host_context_reached("get_extra_trial_vectors_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        call get_extra_trial_vectors_f_wrapper(trial_vectors, error, foreign_context)
        if (error /= 1) then
            test_get_extra_trial_vectors_f_wrapper = .false.
            write(stderr, *) "test_get_extra_trial_vectors_f_wrapper failed: Did "// &
                "not report an invalid context."
        end if

    end function test_get_extra_trial_vectors_f_wrapper

    logical(c_bool) function test_conv_check_stability_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the stability check convergence
        ! check function
        !
        use opentrustregion, only: conv_check_stability_type
        use c_interface, only: c_callbacks_type, conv_check_stability_f_wrapper
        use test_reference, only: test_conv_check_stability_funptr

        procedure(conv_check_stability_type), pointer :: funptr
        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context
        real(rp) :: residual(n_param)
        logical :: converged
        integer(ip) :: error

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%conv_check_stability => mock_conv_check_stability
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! get pointer to subroutine
        funptr => conv_check_stability_f_wrapper

        ! test convergence check wrapper
        test_conv_check_stability_f_wrapper = test_conv_check_stability_funptr( &
            funptr, "conv_check_stability_f_wrapper", "", context)

        ! check that the wrapper handed the host context to the C function
        test_conv_check_stability_f_wrapper = &
            test_conv_check_stability_f_wrapper .and. &
            logical(host_context_reached("conv_check_stability_f_wrapper"), kind=c_bool)

        ! a context this module did not create is reported rather than dereferenced
        foreign_context => foreign_context_target
        residual = 1.0_rp
        converged = &
            conv_check_stability_f_wrapper(residual, 1.0_rp, error, foreign_context)
        if (error /= 1) then
            test_conv_check_stability_f_wrapper = .false.
            write(stderr, *) "test_conv_check_stability_f_wrapper failed: Did not "// &
                "report an invalid context."
        end if

    end function test_conv_check_stability_f_wrapper

    logical(c_bool) function test_logger_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the logging function
        !
        use c_interface, only: c_callbacks_type, logger_f_wrapper

        type(c_callbacks_type), target :: callbacks
        class(*), pointer :: context, foreign_context

        ! assume tests pass
        test_logger_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        callbacks%logger => mock_logger
        call arm_host_context_c(callbacks%host_context)
        context => callbacks

        ! call subroutine
        test_logger = .false.
        call logger_f_wrapper("test", context)

        ! check if logging test boolean is as expected
        if (.not. test_logger) then
            test_logger_f_wrapper = .false.
            write(stderr, *) "test_logger_f_wrapper failed: Returned logging "// &
                "subroutine wrong."
        end if

        ! check that the wrapper handed the host context to the C function
        test_logger_f_wrapper = test_logger_f_wrapper .and. logical( &
            host_context_reached("logger_f_wrapper"), kind=c_bool)

        ! a context this module did not create is dropped rather than dereferenced,
        ! the logger has no error channel so it must simply not be called
        foreign_context => foreign_context_target
        test_logger = .false.
        call logger_f_wrapper("test", foreign_context)
        if (test_logger) then
            test_logger_f_wrapper = .false.
            write(stderr, *) "test_logger_f_wrapper failed: Called logging "// &
                "subroutine with an invalid context."
        end if

    end function test_logger_f_wrapper

    logical(c_bool) function test_init_solver_settings_c() bind(C)
        !
        ! this function tests that the solver settings initialization routine correctly
        ! initializes all settings to their default values
        !
        use c_interface, only: solver_settings_type_c, init_solver_settings_c
        use opentrustregion, only: default_solver_settings
        use test_reference, only: test_associated_solver_c_funptr, operator(/=)

        type(solver_settings_type_c) :: settings

        ! assume test passes
        test_init_solver_settings_c = .true.

        ! initialize settings
        call init_solver_settings_c(settings)

        ! check that callback function pointers are not associated
        if (.not. test_associated_solver_c_funptr(settings, &
                                                  "test_init_solver_settings_c")) &
            test_init_solver_settings_c = .false.

        ! check settings
        if (settings /= default_solver_settings) then
            write(stderr, *) "test_init_solver_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_solver_settings_c = .false.
        end if

    end function test_init_solver_settings_c

    logical(c_bool) function test_init_stability_settings_c() bind(C)
        !
        ! this function tests that the stability check settings initialization routine
        ! correctly initializes all settings to their default values
        !
        use c_interface, only: stability_settings_type_c, init_stability_settings_c
        use opentrustregion, only: default_stability_settings
        use test_reference, only: test_associated_stability_c_funptr, operator(/=)

        type(stability_settings_type_c) :: settings

        ! assume test passes
        test_init_stability_settings_c = .true.

        ! initialize settings
        call init_stability_settings_c(settings)

        ! check that callback function pointers are not associated
        if (.not. test_associated_stability_c_funptr( &
            settings, "test_init_stability_settings_c")) &
            test_init_stability_settings_c = .false.

        ! check settings
        if (settings /= default_stability_settings) then
            write(stderr, *) "test_init_stability_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_stability_settings_c = .false.
        end if

    end function test_init_stability_settings_c

    logical(c_bool) function test_assign_solver_f_c() bind(C)
        !
        ! this function tests that the function that converts solver settings from C to
        ! Fortran correctly perform this conversion
        !
        use c_interface, only: solver_settings_type_c, c_callbacks_type, assignment(=)
        use opentrustregion, only: solver_settings_type, default_solver_settings
        use test_reference, only: assignment(=), ref_settings, &
                                  test_mock_solver_funptr, operator(/=)

        type(solver_settings_type_c) :: settings_c
        type(solver_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks, stability_callbacks

        ! assume test passes
        test_assign_solver_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_settings
        call solver_set_mock(settings_c)

        ! convert to Fortran settings and hand them the callbacks through context
        settings = settings_c
        call solver_mock_callbacks(callbacks)
        call stability_mock_callbacks(stability_callbacks)
        settings%context => callbacks
        settings%stability_settings%context => stability_callbacks

        ! check optional function pointers, which also calls the logging function
        test_logger = .false.
        if (.not. test_mock_solver_funptr(settings, "assign_solver_f_c", "")) &
            test_assign_solver_f_c = .false.
        if (.not. test_logger) then
            test_assign_solver_f_c = .false.
            write(stderr, *) "test_assign_solver_f_c failed: Called logging "// &
                "subroutine wrong."
        end if

        ! check against reference values
        if (settings /= ref_settings) then
            write(stderr, *) "test_assign_solver_f_c failed: Settings not "// &
                "converted correctly."
            test_assign_solver_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            write(stderr, *) "test_assign_solver_f_c failed: Settings not marked "// &
                "as initialized."
            test_assign_solver_f_c = .false.
        end if

        ! convert C settings that were not initialized, the custom values they still
        ! carry have to be replaced by the default settings
        settings_c%initialized = .false.
        settings = settings_c

        ! check that no callback functions were converted
        if (associated(settings%precond) .or. associated(settings%precond_pd) .or. &
            associated(settings%project) .or. &
            associated(settings%get_extra_trial_vectors) .or. &
            associated(settings%conv_check) .or. &
            associated(settings%stability_hess_x) .or. associated(settings%logger)) then
            write(stderr, *) "test_assign_solver_f_c failed: Function pointers "// &
                "converted for settings that were not initialized."
            test_assign_solver_f_c = .false.
        end if

        ! check against default values
        if (settings /= default_solver_settings) then
            write(stderr, *) "test_assign_solver_f_c failed: Settings that were "// &
                "not initialized not converted to default values."
            test_assign_solver_f_c = .false.
        end if

    end function test_assign_solver_f_c

    logical(c_bool) function test_assign_stability_f_c() bind(C)
        !
        ! this function tests that the function that converts stability check settings
        ! from C to Fortran correctly performs this conversion
        !
        use c_interface, only: stability_settings_type_c, c_callbacks_type, &
                               assignment(=)
        use opentrustregion, only: stability_settings_type, default_stability_settings
        use test_reference, only: assignment(=), ref_settings, &
                                  test_mock_stability_funptr, operator(/=)

        type(stability_settings_type_c) :: settings_c
        type(stability_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_stability_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_settings
        call stability_set_mock(settings_c)

        ! convert to Fortran settings and hand them the callbacks through context
        settings = settings_c
        call stability_mock_callbacks(callbacks)
        settings%context => callbacks

        ! check optional function pointers, which also calls the logging function
        test_logger = .false.
        if (.not. test_mock_stability_funptr(settings, "assign_stability_f_c", "")) &
            test_assign_stability_f_c = .false.
        if (.not. test_logger) then
            test_assign_stability_f_c = .false.
            write(stderr, *) "test_assign_stability_f_c failed: Called logging "// &
                "subroutine wrong."
        end if

        ! check against reference values
        if (settings /= ref_settings) then
            write(stderr, *) "test_assign_stability_f_c failed: Settings not "// &
                "converted correctly."
            test_assign_stability_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            test_assign_stability_f_c = .false.
            write(stderr, *) "test_assign_stability_f_c failed: Settings not "// &
                "marked as initialized."
        end if

        ! convert C settings that were not initialized, the custom values they still
        ! carry have to be replaced by the default settings
        settings_c%initialized = .false.
        settings = settings_c

        ! check that no callback functions were converted
        if (associated(settings%precond) .or. associated(settings%project) .or. &
            associated(settings%approx_hess_x) .or. &
            associated(settings%get_extra_trial_vectors) .or. &
            associated(settings%conv_check) .or. associated(settings%logger)) then
            test_assign_stability_f_c = .false.
            write(stderr, *) "test_assign_stability_f_c failed: Function pointers "// &
                "converted for settings that were not initialized."
        end if

        ! check against default values
        if (settings /= default_stability_settings) then
            test_assign_stability_f_c = .false.
            write(stderr, *) "test_assign_stability_f_c failed: Settings that were "// &
                "not initialized not converted to default values."
        end if

    end function test_assign_stability_f_c

    logical(c_bool) function test_assign_solver_c_f() bind(C)
        !
        ! this function tests that the function that converts solver settings from
        ! Fortran to C correctly performs this conversion
        !
        use opentrustregion, only: solver_settings_type
        use c_interface, only: solver_settings_type_c, assignment(=)
        use test_reference, only: ref_settings, assignment(=), &
                                  test_associated_solver_c_funptr, operator(/=)

        type(solver_settings_type) :: settings
        type(solver_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_solver_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_settings

        ! convert to C settings
        settings_c = settings

        ! check that callback function pointers are not associated
        if (.not. test_associated_solver_c_funptr(settings_c, "assign_solver_c_f")) &
            test_assign_solver_c_f = .false.

        ! check against reference values
        if (settings_c /= ref_settings) then
            write(stderr, *) "test_assign_solver_c_f failed: Settings not "// &
                "converted correctly."
            test_assign_solver_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_solver_c_f = .false.
            write(stderr, *) "test_assign_solver_c_f failed: Settings not marked "// &
                "as initialized."
        end if

    end function test_assign_solver_c_f

    logical(c_bool) function test_assign_stability_c_f() bind(C)
        !
        ! this function tests that the function that converts stability check settings
        ! from Fortran to C correctly performs this conversion
        !
        use opentrustregion, only: stability_settings_type
        use c_interface, only: stability_settings_type_c, assignment(=)
        use test_reference, only: ref_settings, assignment(=), &
                                  test_associated_stability_c_funptr, operator(/=)

        type(stability_settings_type) :: settings
        type(stability_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_stability_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_settings

        ! convert to C settings
        settings_c = settings

        ! check that callback function pointers are not associated
        if (.not. test_associated_stability_c_funptr( &
            settings_c, "assign_stability_c_f")) test_assign_stability_c_f = .false.

        ! check against reference values
        if (settings_c /= ref_settings) then
            write(stderr, *) "test_assign_stability_c_f failed: Settings not "// &
                "converted correctly."
            test_assign_stability_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_stability_c_f = .false.
            write(stderr, *) "test_assign_stability_c_f failed: Settings not "// &
                "marked as initialized."
        end if

    end function test_assign_stability_c_f

    logical(c_bool) function test_character_to_c() bind(C)
        !
        ! this function tests conversion of a Fortran character string to a C
        ! null-terminated character array
        !
        use c_interface, only: character_to_c

        character(len=*), parameter :: test_string = "test"
        character(kind=c_char), allocatable :: char_c(:)
        integer :: n, i

        ! assume test passes
        test_character_to_c = .true.

        ! perform conversion
        char_c = character_to_c(test_string)

        ! check length
        n = len_trim(test_string)
        if (size(char_c) /= n + 1) then
            write(stderr, *) "test_character_to_c failed: Character array has "// &
                "wrong size."
            test_character_to_c = .false.
        end if

        ! check characters
        do i = 1, n
            if (char_c(i) /= test_string(i:i)) then
                write(stderr, *) "test_character_to_c failed: Character array "// &
                    "mismatch at character ", i
                test_character_to_c = .false.
            end if
        end do

        ! check null terminator
        if (char_c(n + 1) /= c_null_char) then
            write(stderr, *) "test_character_to_c failed: Character array is "// &
                "missing null terminator."
            test_character_to_c = .false.
        end if

    end function test_character_to_c

    logical(c_bool) function test_character_from_c() bind(C)
        !
        ! this function tests conversion of a C null-terminated character array to a
        ! Fortran character string
        !
        use c_interface, only: character_from_c

        character(kind=c_char), parameter :: test_array(5) = &
            ["t", "e", "s", "t", c_null_char]
        character(len=:), allocatable :: char_f
        integer :: i

        ! assume test passes
        test_character_from_c = .true.

        ! perform conversion
        char_f = character_from_c(test_array)

        ! check equality
        if (len(char_f) /= size(test_array) - 1) then
            write(stderr, *) "test_character_from_c failed: Converted string has "// &
                "wrong length."
            test_character_from_c = .false.
            return
        end if

        ! check characters
        do i = 1, len(char_f)
            if (test_array(i) /= char_f(i:i)) then
                write(stderr, *) "test_character_from_c failed: String mismatch at "// &
                    "character ", i
                test_character_from_c = .false.
            end if
        end do

    end function test_character_from_c

end module c_interface_unit_tests
