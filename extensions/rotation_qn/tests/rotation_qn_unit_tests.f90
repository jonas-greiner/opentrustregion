! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_rotation_qn_unit_tests

    use opentrustregion_core, only: rp, ip, stderr, hess_x_type, precond_type, &
                                    project_type
    use test_reference, only: n_param, tol, is_close, is_at_most, host_context_type
    use otr_rotation_qn, only: transport_type, obj_func_grad_type, sr1_updating_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(transport_type), pointer :: mock_transport_ptr => mock_transport
    procedure(hess_x_type), pointer :: mock_init_hess_ptr => mock_init_hess
    procedure(precond_type), pointer :: mock_init_precond_ptr => mock_init_precond
    procedure(obj_func_grad_type), pointer :: mock_obj_func_grad_ptr => &
        mock_obj_func_grad
    procedure(project_type), pointer :: mock_project_ptr => mock_project

    ! diagonal of the initial Hessian the mock applies, indefinite so that the BFGS
    ! scheme meets a step along which the model so far has no positive curvature
    real(rp), parameter :: init_hess_diag(n_param) = [2.0_rp, 3.0_rp, -1.0_rp]

    ! host context of the host record of the solver level of the objects the tests set
    ! up, which the host's objective function and projection receive, distinct from
    ! that of the settings, which the transport function and the objective function,
    ! initial Hessian and preconditioner of the settings receive
    type(host_context_type), target :: level_host_context

    ! rotation QN object which records the point its orbital update adds and returns
    ! the error code it is set to, for the test of the orbital update
    type, extends(sr1_updating_type) :: recording_rotation_qn_type
        integer(ip) :: n_added = 0, error = 0
        real(rp), allocatable :: added_kappa(:), added_grad(:)
    contains
        procedure :: add => recording_add_rotation_qn
    end type

    ! rotation QN object which records the builds of its model and the number of pairs
    ! of the history when it was built and returns the error code it is set to, for the
    ! tests of the objective function and the addition
    type, extends(sr1_updating_type) :: recording_build_type
        integer(ip) :: n_built = 0, n_points_built = -1, error = 0
    contains
        procedure :: build_model => recording_build_model
    end type

contains

    subroutine mock_transport(geodesic, tangent_vector, error, context)
        !
        ! this subroutine is a mock of the host's transport function, which adds the
        ! geodesic to the tangent vector, so that the test sees which vectors were
        ! transported along which geodesic, records its call, checks the host context
        ! and returns the error code set in it or a value that is not a number
        !
        use otr_common_mock, only: record_host_call, inject_host_nan

        real(rp), intent(in), target :: geodesic(:)
        real(rp), intent(inout), target :: tangent_vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        error = record_host_call("transport", context)
        tangent_vector = tangent_vector + geodesic
        call inject_host_nan(context, 1_ip, tangent_vector(1))

    end subroutine mock_transport

    subroutine mock_init_hess(x, hess_x, error, context)
        !
        ! this subroutine is a mock of the initial Hessian of the settings, which
        ! applies the diagonal init_hess_diag, records its call, checks the host context
        ! and returns the error code set in it or a value that is not a number
        !
        use otr_common_mock, only: record_host_call, inject_host_nan

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        error = record_host_call("init_hess", context)
        hess_x = init_hess_diag * x
        call inject_host_nan(context, 1_ip, hess_x(1))

    end subroutine mock_init_hess

    subroutine mock_init_precond(residual, mu, precond_residual, error, context)
        !
        ! this subroutine is a mock of the preconditioner of the settings, the exact
        ! inverse of the level-shifted initial Hessian of the mock, which records its
        ! call, checks the host context and returns the error code set in it or a value
        ! that is not a number
        !
        use otr_common_mock, only: record_host_call, inject_host_nan

        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        error = record_host_call("precond", context)
        precond_residual = residual / (init_hess_diag - mu)
        call inject_host_nan(context, 1_ip, precond_residual(1))

    end subroutine mock_init_precond

    function mock_obj_func_grad(kappa, grad, error, context) result(func)
        !
        ! this function is a mock of the objective function of the settings, which
        ! returns the sum of the step as the energy and twice the step as the gradient,
        ! records its call, checks the host context and returns the error code set in
        ! it or a value that is not a number as the energy (position 1) or in the
        ! gradient (position 2)
        !
        use otr_common_mock, only: record_host_call, inject_host_nan

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out), target :: grad(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context
        real(rp) :: func

        error = record_host_call("obj_func_grad", context)
        func = sum(kappa)
        grad = 2.0_rp * kappa
        call inject_host_nan(context, 1_ip, func)
        call inject_host_nan(context, 2_ip, grad(1))

    end function mock_obj_func_grad

    subroutine mock_project(vector, error, context)
        !
        ! this subroutine is a mock of the host's projection, which removes the first
        ! component, records its call, checks the host context and returns the error
        ! code set in it
        !
        use otr_common_mock, only: record_host_call

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        error = record_host_call("project", context)
        vector(1) = 0.0_rp

    end subroutine mock_project

    subroutine recording_add_rotation_qn(self, kappa, grad, error)
        !
        ! this subroutine records the step and gradient of the point the recording
        ! rotation QN object is given and returns the error code it is set to
        !
        class(recording_rotation_qn_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), grad(:)
        integer(ip), intent(out) :: error

        self%n_added = self%n_added + 1
        self%added_kappa = kappa
        self%added_grad = grad
        error = self%error

    end subroutine recording_add_rotation_qn

    subroutine recording_build_model(self, error)
        !
        ! this subroutine records the build of the model of the recording rotation QN
        ! object and the number of pairs of its history and returns the error code it
        ! is set to
        !
        class(recording_build_type), intent(inout) :: self
        integer(ip), intent(out) :: error

        self%n_built = self%n_built + 1
        self%n_points_built = self%n_points
        error = self%error

    end subroutine recording_build_model

    subroutine generate_random_terms(object)
        !
        ! this subroutine gives an updating object a model of three random terms of
        ! either sign
        !
        use otr_rotation_qn, only: updating_type

        class(updating_type), intent(inout) :: object

        object%n_terms = 3
        call random_number(object%vectors(:, :3))
        call random_number(object%coeffs(:3))
        object%coeffs(:3) = object%coeffs(:3) - 0.5_rp

    end subroutine generate_random_terms

    subroutine setup_updating_object(object)
        !
        ! this subroutine sets up a new updating object as after three orbital updates,
        ! without the routines under test, with the default settings carrying the armed
        ! host context and the mock initial Hessian, the mock transport function, a
        ! host context of its own for the host record of its solver level, a random
        ! positive Hessian diagonal, two random points in its history, whose gradient
        ! differences are twice their steps, a random last gradient and a model without
        ! terms, whose arrays hold two terms per point
        !
        use otr_rotation_qn, only: updating_type, default_rotation_qn_settings, &
                                   error_rotation_qn
        use test_reference, only: arm_host_context

        class(updating_type), intent(out) :: object

        object%error_origin = error_rotation_qn
        object%n_param = n_param
        object%settings = default_rotation_qn_settings
        call arm_host_context(object%settings)
        object%settings%init_hess => mock_init_hess
        object%transport => mock_transport
        level_host_context = host_context_type()
        object%solver_level%host%context => level_host_context
        allocate(object%kappa_list(n_param, object%settings%max_points), &
                 object%grad_diff_list(n_param, object%settings%max_points), &
                 object%last_grad(n_param), object%h_diag(n_param), &
                 object%vectors(n_param, 2 * object%settings%max_points), &
                 object%coeffs(2 * object%settings%max_points))
        call random_number(object%h_diag)
        object%h_diag = object%h_diag + 1.0_rp
        object%kappa_list = 0.0_rp
        object%grad_diff_list = 0.0_rp
        object%n_points = 2
        call random_number(object%kappa_list(:, :2))
        object%grad_diff_list(:, :2) = 2.0_rp * object%kappa_list(:, :2)
        call random_number(object%last_grad)
        object%vectors = 0.0_rp
        object%coeffs = 0.0_rp
        object%n_terms = 0

    end subroutine setup_updating_object

    function model_matrix(object) result(matrix)
        !
        ! this function returns the matrix of the model of an updating object, the
        ! initial Hessian of the mock updated by the rank-1 terms of the object
        !
        use otr_rotation_qn, only: updating_type

        class(updating_type), intent(in) :: object
        real(rp) :: matrix(n_param, n_param)

        integer(ip) :: i, k

        matrix = 0.0_rp
        do i = 1, n_param
            matrix(i, i) = init_hess_diag(i)
        end do
        do k = 1, object%n_terms
            matrix = matrix + &
                     object%coeffs(k) * spread(object%vectors(:, k), 2, n_param) * &
                     spread(object%vectors(:, k), 1, n_param)
        end do

    end function model_matrix

    logical(c_bool) function test_rotation_qn_factory() bind(C)
        !
        ! this function tests the factory, which creates a rotation QN object of the
        ! updating scheme the settings select after their sanity check, initializes it
        ! with the settings and the host's transport function, wires the solver settings
        ! with it and the record of the state and lets the solver refresh the model
        ! after a rejected step if the settings have an objective function
        !
        use opentrustregion_core, only: update_orbs_type, obj_func_type, &
                                        solver_state_type
        use opentrustregion_unit_tests, only: setup_settings, test_context_type, logger
        use otr_common_test_reference, only: same_logger_and_context
        use otr_common, only: extension_solver_settings_type, extension_level_type, &
                              update_orbs_callback, obj_func_callback
        use otr_common_unit_tests, only: wired_level
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_rotation_qn, only: rotation_qn_factory, rotation_qn_settings_type, &
                                   updating_type, sr1_updating_type, &
                                   bfgs_updating_type, default_rotation_qn_settings, &
                                   error_rotation_qn
        use otr_rotation_qn_test_reference, only: ref_rotation_qn_settings, operator(/=)

        procedure(update_orbs_type), pointer :: update_orbs_funptr, &
                                                update_orbs_rotation_qn_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr, &
                                             obj_func_rotation_qn_funptr
        procedure(transport_type), pointer :: transport_funptr
        type(extension_solver_settings_type), allocatable, target :: solver_settings, &
                                                                     other_settings
        type(solver_state_type), target :: state, other_state
        type(rotation_qn_settings_type) :: settings, expected_settings
        type(test_context_type), target :: context
        type(extension_level_type), pointer :: level, other_level
        integer(ip) :: i_case, error
        logical :: sr1_expected, scheme_right, owned
        character(len=67), parameter :: case_names(3) = &
            [character(len=67) :: "for settings selecting SR1 in uppercase with an "// &
             "objective function", "for settings selecting BFGS in uppercase", &
             "without settings"]

        ! assume tests pass
        test_rotation_qn_factory = .true.

        ! set the host's callback functions
        update_orbs_funptr => mock_update_orbs
        obj_func_funptr => mock_obj_func
        transport_funptr => mock_transport

        ! create objects for settings selecting either updating scheme, the first with
        ! an objective function, and for no settings, each of which has to be the
        ! object of the selected scheme, initialized with the settings after their
        ! sanity check and the host's transport function, and to wire the solver
        ! settings, which own it, with the host's orbital update and objective function
        ! and the record of the state and let them refresh the model after a rejected
        ! step only with the objective function of the settings
        do i_case = 1, size(case_names)
            allocate(solver_settings)
            select case (i_case)
            case (1)
                settings = ref_rotation_qn_settings
                settings%qn_method = "SR1"
                settings%logger => logger
                settings%context => context
                settings%init_hess => mock_init_hess
                settings%init_precond => mock_init_precond
                settings%obj_func_grad => mock_obj_func_grad
                expected_settings = settings
                expected_settings%qn_method = "sr1"
            case (2)
                settings = ref_rotation_qn_settings
                settings%qn_method = "BFGS"
                expected_settings = ref_rotation_qn_settings
            case (3)
                expected_settings = default_rotation_qn_settings
            end select
            sr1_expected = expected_settings%qn_method == "sr1"
            if (i_case == 3) then
                call rotation_qn_factory( &
                    n_param, update_orbs_funptr, obj_func_funptr, transport_funptr, &
                    update_orbs_rotation_qn_funptr, obj_func_rotation_qn_funptr, &
                    solver_settings, state, error)
            else
                call rotation_qn_factory( &
                    n_param, update_orbs_funptr, obj_func_funptr, transport_funptr, &
                    update_orbs_rotation_qn_funptr, obj_func_rotation_qn_funptr, &
                    solver_settings, state, error, settings)
            end if
            if (error /= 0) then
                write(stderr, *) "test_rotation_qn_factory failed: Produced error "// &
                    trim(case_names(i_case))//"."
                test_rotation_qn_factory = .false.
                deallocate(solver_settings)
                cycle
            end if
            if (.not. associated(update_orbs_rotation_qn_funptr, &
                                 update_orbs_callback) .or. &
                .not. associated(obj_func_rotation_qn_funptr, obj_func_callback)) then
                write(stderr, *) "test_rotation_qn_factory failed: Role callback "// &
                    "functions not returned "//trim(case_names(i_case))//"."
                test_rotation_qn_factory = .false.
            end if
            level => wired_level(solver_settings)
            if (.not. associated(level)) then
                write(stderr, *) "test_rotation_qn_factory failed: Solver settings "// &
                    "not wired with an object "//trim(case_names(i_case))//"."
                test_rotation_qn_factory = .false.
                deallocate(solver_settings)
                cycle
            end if
            owned = .false.
            if (allocated(solver_settings%extensions)) then
                if (size(solver_settings%extensions) == 1) owned = &
                    associated(level%extension, solver_settings%extensions(1)%object)
            end if
            if (.not. owned) then
                write(stderr, *) "test_rotation_qn_factory failed: Object not "// &
                    "owned by the solver settings "//trim(case_names(i_case))//"."
                test_rotation_qn_factory = .false.
            end if
            select type (object => level%extension)
            type is (sr1_updating_type)
                scheme_right = sr1_expected
            type is (bfgs_updating_type)
                scheme_right = .not. sr1_expected
            class default
                scheme_right = .false.
            end select
            if (.not. scheme_right) then
                write(stderr, *) "test_rotation_qn_factory failed: Object of the "// &
                    "wrong updating scheme created "//trim(case_names(i_case))//"."
                test_rotation_qn_factory = .false.
            end if
            select type (object => level%extension)
            class is (updating_type)
                if (object%n_param /= n_param .or. &
                    object%settings /= expected_settings .or. &
                    .not. same_logger_and_context(object%settings, expected_settings)) &
                    then
                    write(stderr, *) "test_rotation_qn_factory failed: Object not "// &
                        "initialized with the settings "//trim(case_names(i_case))//"."
                    test_rotation_qn_factory = .false.
                end if
                if (.not. associated(object%transport, mock_transport)) then
                    write(stderr, *) "test_rotation_qn_factory failed: Host's "// &
                        "transport function not handed to the object "// &
                        trim(case_names(i_case))//"."
                    test_rotation_qn_factory = .false.
                end if
                if ((associated(object%settings%init_hess, mock_init_hess) .and. &
                     associated(object%settings%init_precond, mock_init_precond) .and. &
                     associated(object%settings%obj_func_grad, &
                                mock_obj_func_grad)) .neqv. i_case == 1) then
                    write(stderr, *) "test_rotation_qn_factory failed: Callback "// &
                        "functions of the settings not handed to the object "// &
                        trim(case_names(i_case))//"."
                    test_rotation_qn_factory = .false.
                end if
            end select
            if (.not. associated(level%host%update_orbs, mock_update_orbs) .or. &
                .not. associated(level%host%obj_func, mock_obj_func) .or. &
                .not. associated(level%extension%state, state)) then
                write(stderr, *) "test_rotation_qn_factory failed: Solver settings "// &
                    "not wired with the host's functions and the record of the "// &
                    "state "//trim(case_names(i_case))//"."
                test_rotation_qn_factory = .false.
            end if
            if (solver_settings%refresh_hess .neqv. i_case == 1) then
                write(stderr, *) "test_rotation_qn_factory failed: Refresh of the "// &
                    "model after a rejected step set wrong "// &
                    trim(case_names(i_case))//"."
                test_rotation_qn_factory = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create an object for parameters the sanity check rejects, which has to return
        ! its error without functions and leave the solver settings untouched
        allocate(solver_settings)
        call setup_settings(settings, context)
        update_orbs_rotation_qn_funptr => update_orbs_callback
        obj_func_rotation_qn_funptr => obj_func_callback
        call rotation_qn_factory(0_ip, update_orbs_funptr, obj_func_funptr, &
                                 transport_funptr, update_orbs_rotation_qn_funptr, &
                                 obj_func_rotation_qn_funptr, solver_settings, state, &
                                 error, settings)
        if (error /= error_rotation_qn + 1) then
            write(stderr, *) "test_rotation_qn_factory failed: Error not returned "// &
                "for parameters the sanity check rejects."
            test_rotation_qn_factory = .false.
        end if
        if (associated(update_orbs_rotation_qn_funptr) .or. &
            associated(obj_func_rotation_qn_funptr)) then
            write(stderr, *) "test_rotation_qn_factory failed: Functions returned "// &
                "for parameters the sanity check rejects."
            test_rotation_qn_factory = .false.
        end if
        if (associated(solver_settings%context) .or. &
            allocated(solver_settings%extensions)) then
            write(stderr, *) "test_rotation_qn_factory failed: Solver settings "// &
                "changed for parameters the sanity check rejects."
            test_rotation_qn_factory = .false.
        end if
        deallocate(solver_settings)

        ! create an object for solver settings a factory wired, given the host's own
        ! orbital update and objective function instead of those the factory returned,
        ! which the wiring refuses, so that the factory has to return its error without
        ! functions and leave the settings wired with the first object only
        allocate(solver_settings)
        call setup_settings(settings, context)
        call rotation_qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                                 transport_funptr, update_orbs_rotation_qn_funptr, &
                                 obj_func_rotation_qn_funptr, solver_settings, state, &
                                 error, settings)
        level => wired_level(solver_settings)
        call rotation_qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                                 transport_funptr, update_orbs_rotation_qn_funptr, &
                                 obj_func_rotation_qn_funptr, solver_settings, state, &
                                 error, settings)
        if (error /= error_rotation_qn + 1) then
            write(stderr, *) "test_rotation_qn_factory failed: Error not returned "// &
                "for the host's own functions given with solver settings a factory "// &
                "wired."
            test_rotation_qn_factory = .false.
        end if
        if (associated(update_orbs_rotation_qn_funptr) .or. &
            associated(obj_func_rotation_qn_funptr)) then
            write(stderr, *) "test_rotation_qn_factory failed: Functions returned "// &
                "for the host's own functions given with solver settings a factory "// &
                "wired."
            test_rotation_qn_factory = .false.
        end if
        owned = .false.
        if (associated(level) .and. allocated(solver_settings%extensions)) &
            owned = size(solver_settings%extensions) == 1 .and. &
                    associated(wired_level(solver_settings), level)
        if (.not. owned) then
            write(stderr, *) "test_rotation_qn_factory failed: Solver settings a "// &
                "factory wired changed for the host's own functions."
            test_rotation_qn_factory = .false.
        end if
        deallocate(solver_settings)

        ! create objects for two separate solver settings, which have to be distinct
        ! objects reading their own record of the state, so that calculations with
        ! separate settings can be used alternately
        allocate(solver_settings, other_settings)
        settings = ref_rotation_qn_settings
        call rotation_qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                                 transport_funptr, update_orbs_rotation_qn_funptr, &
                                 obj_func_rotation_qn_funptr, solver_settings, state, &
                                 error, settings)
        if (error /= 0) then
            write(stderr, *) "test_rotation_qn_factory failed: Produced error for "// &
                "the first of two separate solver settings."
            test_rotation_qn_factory = .false.
        end if
        call rotation_qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                                 transport_funptr, update_orbs_rotation_qn_funptr, &
                                 obj_func_rotation_qn_funptr, other_settings, &
                                 other_state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_rotation_qn_factory failed: Produced error for "// &
                "the second of two separate solver settings."
            test_rotation_qn_factory = .false.
        end if
        level => wired_level(solver_settings)
        other_level => wired_level(other_settings)
        if (.not. (associated(level) .and. associated(other_level))) then
            write(stderr, *) "test_rotation_qn_factory failed: Separate solver "// &
                "settings not wired."
            test_rotation_qn_factory = .false.
        else if (associated(level, other_level) .or. &
                 .not. associated(level%extension%state, state) .or. &
                 .not. associated(other_level%extension%state, other_state)) then
            write(stderr, *) "test_rotation_qn_factory failed: Separate solver "// &
                "settings not wired with distinct objects reading their own state."
            test_rotation_qn_factory = .false.
        end if
        deallocate(solver_settings, other_settings)

    end function test_rotation_qn_factory

    logical(c_bool) function test_rotation_qn_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! rotation QN input parameters
        !
        use opentrustregion_unit_tests, only: test_context_type, setup_error_logging
        use otr_rotation_qn, only: rotation_qn_settings_type, &
                                   rotation_qn_sanity_check, &
                                   default_rotation_qn_settings

        type(rotation_qn_settings_type) :: settings
        type(test_context_type), target :: log_context
        integer(ip) :: case_n_param, i_case, error
        character(len=35), parameter :: case_names(3) = &
            [character(len=35) :: "a vanishing number of parameters", &
             "a negative maximum number of points", "an unknown updating scheme"]

        ! assume tests pass
        test_rotation_qn_sanity_check = .true.

        ! check valid parameters, with no stored points, which keep no history, and the
        ! updating scheme in mixed case, which has to be converted to lowercase
        settings = default_rotation_qn_settings
        call setup_error_logging(settings, log_context)
        settings%max_points = 0
        settings%qn_method = "BfGs"
        call rotation_qn_sanity_check(settings, 1_ip, error)
        if (error /= 0 .or. log_context%log_message /= "") then
            write(stderr, *) "test_rotation_qn_sanity_check failed: Valid "// &
                "parameters rejected."
            test_rotation_qn_sanity_check = .false.
        end if
        if (settings%qn_method /= "bfgs") then
            write(stderr, *) "test_rotation_qn_sanity_check failed: Updating "// &
                "scheme not converted to lowercase."
            test_rotation_qn_sanity_check = .false.
        end if

        ! check that every invalid parameter is rejected with an error message
        do i_case = 1, size(case_names)
            settings = default_rotation_qn_settings
            call setup_error_logging(settings, log_context)
            case_n_param = 1
            select case (i_case)
            case (1)
                case_n_param = 0
            case (2)
                settings%max_points = -1
            case (3)
                settings%qn_method = "dfp"
            end select
            call rotation_qn_sanity_check(settings, case_n_param, error)
            if (error /= 1) then
                write(stderr, *) "test_rotation_qn_sanity_check failed: Error not "// &
                    "returned for "//trim(case_names(i_case))//"."
                test_rotation_qn_sanity_check = .false.
            end if
            if (log_context%log_message == "") then
                write(stderr, *) "test_rotation_qn_sanity_check failed: No error "// &
                    "message printed for "//trim(case_names(i_case))//"."
                test_rotation_qn_sanity_check = .false.
            end if
        end do

    end function test_rotation_qn_sanity_check

    logical(c_bool) function test_invoke_transport() bind(C)
        !
        ! this function tests the subroutine which calls the host's transport function
        !
        use otr_rotation_qn, only: rotation_qn_settings_type, invoke_transport
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              setup_error_logging
        use test_reference, only: arm_host_context, host_context_reached

        integer(ip), parameter :: error_origin = 4200
        type(rotation_qn_settings_type) :: settings
        type(test_context_type), target :: context
        procedure(transport_type), pointer :: transport_funptr
        real(rp) :: geodesic(n_param), vector(n_param), expected(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_invoke_transport = .true.

        ! call a transport function, which adds the geodesic to the tangent vector,
        ! with a random geodesic and tangent vector and check that it receives the host
        ! context and that the transported tangent vector is returned
        call setup_settings(settings, context)
        call arm_host_context(settings, context)
        transport_funptr => mock_transport
        call random_number(geodesic)
        call random_number(vector)
        expected = vector + geodesic
        call invoke_transport(transport_funptr, geodesic, vector, settings, &
                              error_origin, error)
        if (error /= 0) then
            write(stderr, *) "test_invoke_transport failed: Produced error."
            test_invoke_transport = .false.
        end if
        if (.not. all(is_close(vector, expected))) then
            write(stderr, *) "test_invoke_transport failed: Transported tangent "// &
                "vector not returned."
            test_invoke_transport = .false.
        end if
        if (.not. host_context_reached("invoke_transport", context)) &
            test_invoke_transport = .false.

        ! let the function fail and check that its error is passed on unchanged
        context%mock_error = 3
        call invoke_transport(transport_funptr, geodesic, vector, settings, &
                              error_origin, error)
        if (error /= 3) then
            write(stderr, *) "test_invoke_transport failed: Error of the function "// &
                "not passed on unchanged."
            test_invoke_transport = .false.
        end if
        context%mock_error = 0

        ! let the function return a value that is not a number, which has to be
        ! reported with the origin of the object and an error message
        call setup_error_logging(settings, context)
        context%nan_output = 1
        call invoke_transport(transport_funptr, geodesic, vector, settings, &
                              error_origin, error)
        context%nan_output = 0
        if (error /= error_origin + 1) then
            write(stderr, *) "test_invoke_transport failed: Tangent vector that is "// &
                "not finite not reported with the origin of the object."
            test_invoke_transport = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_invoke_transport failed: No error message "// &
                "printed for a tangent vector that is not finite."
            test_invoke_transport = .false.
        end if

    end function test_invoke_transport

    logical(c_bool) function test_invoke_init_hess() bind(C)
        !
        ! this function tests the subroutine which applies the initial Hessian of the
        ! settings
        !
        use otr_rotation_qn, only: rotation_qn_settings_type, invoke_init_hess
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              setup_error_logging
        use test_reference, only: arm_host_context, host_context_reached

        integer(ip), parameter :: error_origin = 4200
        type(rotation_qn_settings_type) :: settings
        type(test_context_type), target :: context
        procedure(hess_x_type), pointer :: init_hess_funptr
        real(rp) :: x(n_param), init_hess_x(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_invoke_init_hess = .true.

        ! apply an initial Hessian, which applies a diagonal, to a random vector and
        ! check that it receives the host context and that the result is returned
        call setup_settings(settings, context)
        call arm_host_context(settings, context)
        init_hess_funptr => mock_init_hess
        call random_number(x)
        call invoke_init_hess(init_hess_funptr, x, init_hess_x, settings, &
                              error_origin, error)
        if (error /= 0) then
            write(stderr, *) "test_invoke_init_hess failed: Produced error."
            test_invoke_init_hess = .false.
        end if
        if (.not. all(is_close(init_hess_x, init_hess_diag * x))) then
            write(stderr, *) "test_invoke_init_hess failed: Result not returned."
            test_invoke_init_hess = .false.
        end if
        if (.not. host_context_reached("invoke_init_hess", context)) &
            test_invoke_init_hess = .false.

        ! let the initial Hessian fail and check that its error is passed on unchanged
        context%mock_error = 3
        call invoke_init_hess(init_hess_funptr, x, init_hess_x, settings, &
                              error_origin, error)
        if (error /= 3) then
            write(stderr, *) "test_invoke_init_hess failed: Error of the initial "// &
                "Hessian not passed on unchanged."
            test_invoke_init_hess = .false.
        end if
        context%mock_error = 0

        ! let the initial Hessian return a value that is not a number, which has to be
        ! reported with the origin of the object and an error message
        call setup_error_logging(settings, context)
        context%nan_output = 1
        call invoke_init_hess(init_hess_funptr, x, init_hess_x, settings, &
                              error_origin, error)
        context%nan_output = 0
        if (error /= error_origin + 1) then
            write(stderr, *) "test_invoke_init_hess failed: Result that is not "// &
                "finite not reported with the origin of the object."
            test_invoke_init_hess = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_invoke_init_hess failed: No error message "// &
                "printed for a result that is not finite."
            test_invoke_init_hess = .false.
        end if

    end function test_invoke_init_hess

    logical(c_bool) function test_invoke_obj_func_grad() bind(C)
        !
        ! this function tests the subroutine which calls the objective function of the
        ! settings, which also returns the gradient at the trial point
        !
        use otr_rotation_qn, only: rotation_qn_settings_type, invoke_obj_func_grad
        use opentrustregion_unit_tests, only: test_context_type, setup_settings, &
                                              setup_error_logging
        use test_reference, only: arm_host_context, host_context_reached
        use, intrinsic :: ieee_arithmetic, only: ieee_class, ieee_positive_inf, &
                                                 operator(==)

        integer(ip), parameter :: error_origin = 4200
        type(rotation_qn_settings_type) :: settings
        type(test_context_type), target :: context
        procedure(obj_func_grad_type), pointer :: obj_func_grad_funptr
        real(rp) :: kappa(n_param), grad(n_param), func
        integer(ip) :: error

        ! assume tests pass
        test_invoke_obj_func_grad = .true.

        ! call an objective function, which returns the sum of the step as the energy
        ! and twice the step as the gradient, with a random step and check that it
        ! receives the host context and that both are returned
        call setup_settings(settings, context)
        call arm_host_context(settings, context)
        obj_func_grad_funptr => mock_obj_func_grad
        call random_number(kappa)
        call invoke_obj_func_grad(obj_func_grad_funptr, kappa, func, grad, settings, &
                                  error_origin, error)
        if (error /= 0) then
            write(stderr, *) "test_invoke_obj_func_grad failed: Produced error."
            test_invoke_obj_func_grad = .false.
        end if
        if (.not. is_close(func, sum(kappa)) .or. &
            .not. all(is_close(grad, 2.0_rp * kappa))) then
            write(stderr, *) "test_invoke_obj_func_grad failed: Energy or gradient "// &
                "not returned."
            test_invoke_obj_func_grad = .false.
        end if
        if (.not. host_context_reached("invoke_obj_func_grad", context)) &
            test_invoke_obj_func_grad = .false.

        ! let the function fail and check that its error is passed on unchanged
        context%mock_error = 3
        call invoke_obj_func_grad(obj_func_grad_funptr, kappa, func, grad, settings, &
                                  error_origin, error)
        if (error /= 3) then
            write(stderr, *) "test_invoke_obj_func_grad failed: Error of the "// &
                "function not passed on unchanged."
            test_invoke_obj_func_grad = .false.
        end if
        context%mock_error = 0

        ! let the function return an energy and a gradient that are not a number, the
        ! energy has to count as infinite and the gradient of an infinite energy is not
        ! checked
        context%nan_output = -1
        call invoke_obj_func_grad(obj_func_grad_funptr, kappa, func, grad, settings, &
                                  error_origin, error)
        context%nan_output = 0
        if (error /= 0 .or. .not. ieee_class(func) == ieee_positive_inf) then
            write(stderr, *) "test_invoke_obj_func_grad failed: Energy that is not "// &
                "a number not returned as infinite without an error."
            test_invoke_obj_func_grad = .false.
        end if

        ! let the function return a finite energy and a gradient that is not a number,
        ! which has to be reported with the origin of the object and an error message
        call setup_error_logging(settings, context)
        context%nan_output = 2
        call invoke_obj_func_grad(obj_func_grad_funptr, kappa, func, grad, settings, &
                                  error_origin, error)
        context%nan_output = 0
        if (error /= error_origin + 1) then
            write(stderr, *) "test_invoke_obj_func_grad failed: Gradient that is "// &
                "not finite not reported with the origin of the object."
            test_invoke_obj_func_grad = .false.
        end if
        if (len_trim(context%log_message) == 0) then
            write(stderr, *) "test_invoke_obj_func_grad failed: No error message "// &
                "printed for a gradient that is not finite."
            test_invoke_obj_func_grad = .false.
        end if

    end function test_invoke_obj_func_grad

    logical(c_bool) function test_update_orbs_rotation_qn() bind(C)
        !
        ! this function tests the rotation QN orbital update, which forwards to the
        ! host's orbital update with the host context of the solver level, returns its
        ! results and error, keeps the Hessian diagonal, returns the quasi-Newton
        ! Hessian linear transformation in place of the host's and adds the point for a
        ! nonzero step and for the starting point, whose gradient the first step pair
        ! needs, passing on the error of the addition
        !
        use opentrustregion_core, only: hess_x_type
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: hess_x_callback
        use otr_common_mock, only: mock_update_orbs
        use otr_rotation_qn, only: default_rotation_qn_settings, error_rotation_qn
        use test_reference, only: arm_host_context, host_context_reached

        type(recording_rotation_qn_type), target :: object
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: kappa(n_param), func, grad(n_param), h_diag(n_param)
        procedure(hess_x_type), pointer :: hess_x_funptr
        integer(ip) :: error, i_case, i
        logical :: kept
        character(len=31), parameter :: case_names(5) = &
            [character(len=31) :: "the starting point", "a zero step after the start", &
             "a nonzero step", "a failing host's orbital update", "a failing addition"]
        logical, parameter :: case_adds(5) = [.true., .false., .true., .false., .true.]
        integer(ip), parameter :: case_errors(5) = [0, 0, 0, 3, 4]

        ! assume tests pass
        test_update_orbs_rotation_qn = .true.

        ! set up an object whose solver level forwards to the host's mock orbital update
        ! with a host context distinct from that of its settings
        object%error_origin = error_rotation_qn
        object%n_param = n_param
        object%settings = default_rotation_qn_settings
        object%solver_level%host%update_orbs => mock_update_orbs
        object%solver_level%host%context => contexts(1)
        context => object%solver_level

        ! call the orbital update with a zero step before and after the starting point
        ! was added, with a nonzero step, whose entries sum to a negative number, so
        ! that only their magnitudes show it is nonzero, and with a failing host's
        ! orbital update or addition, and determine if it forwards with the host
        ! context of the solver level, returns the host's results and error, keeps the
        ! Hessian diagonal, returns the role callback function of the Hessian linear
        ! transformation and adds the point exactly where it has to
        do i_case = 1, size(case_names)
            call arm_host_context(object%settings, contexts(2))
            if (allocated(object%h_diag)) deallocate(object%h_diag)
            contexts(1)%n_calls = 0
            contexts(1)%mock_error = 0
            if (i_case == 4) contexts(1)%mock_error = case_errors(i_case)
            object%error = 0
            if (i_case == 5) object%error = case_errors(i_case)
            object%n_added = 0
            if (allocated(object%last_grad)) deallocate(object%last_grad)
            if (i_case > 1) object%last_grad = [(0.1_rp * i, i=1, n_param)]
            kappa = 0.0_rp
            if (i_case > 2) then
                call random_number(kappa)
                kappa = kappa - 1.0_rp
            end if
            call object%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                    context)
            if (error /= case_errors(i_case)) then
                write(stderr, *) "test_update_orbs_rotation_qn failed: Returned "// &
                    "error code wrong for "//trim(case_names(i_case))//"."
                test_update_orbs_rotation_qn = .false.
            end if
            if (.not. is_close(func, sum(kappa)) .or. &
                .not. all(is_close(grad, 2.0_rp * kappa)) .or. &
                .not. all(is_close(h_diag, 3.0_rp * kappa))) then
                write(stderr, *) "test_update_orbs_rotation_qn failed: Results of "// &
                    "the host's orbital update not returned for "// &
                    trim(case_names(i_case))//"."
                test_update_orbs_rotation_qn = .false.
            end if
            if (contexts(1)%n_calls /= 1 .or. contexts(2)%n_calls /= 0) then
                write(stderr, *) "test_update_orbs_rotation_qn failed: Host's "// &
                    "orbital update not called with the host context of the solver "// &
                    "level for "//trim(case_names(i_case))//"."
                test_update_orbs_rotation_qn = .false.
            end if
            if (.not. host_context_reached("update_orbs_rotation_qn", contexts(1))) &
                test_update_orbs_rotation_qn = .false.
            if (i_case /= 4) then
                if (.not. associated(hess_x_funptr, hess_x_callback)) then
                    write(stderr, *) "test_update_orbs_rotation_qn failed: Role "// &
                        "callback function of the Hessian linear transformation "// &
                        "not returned for "//trim(case_names(i_case))//"."
                    test_update_orbs_rotation_qn = .false.
                end if
                kept = allocated(object%h_diag)
                if (kept) kept = all(is_close(object%h_diag, h_diag, 0.0_rp))
                if (.not. kept) then
                    write(stderr, *) "test_update_orbs_rotation_qn failed: Hessian "// &
                        "diagonal not kept for "//trim(case_names(i_case))//"."
                    test_update_orbs_rotation_qn = .false.
                end if
            end if
            if (object%n_added /= merge(1, 0, case_adds(i_case))) then
                write(stderr, *) "test_update_orbs_rotation_qn failed: Point added "// &
                    "wrongly for "//trim(case_names(i_case))//"."
                test_update_orbs_rotation_qn = .false.
            else if (case_adds(i_case)) then
                if (.not. all(is_close(object%added_kappa, kappa)) .or. &
                    .not. all(is_close(object%added_grad, grad))) then
                    write(stderr, *) "test_update_orbs_rotation_qn failed: Step or "// &
                        "gradient of the added point wrong for "// &
                        trim(case_names(i_case))//"."
                    test_update_orbs_rotation_qn = .false.
                end if
            end if
        end do

    end function test_update_orbs_rotation_qn

    logical(c_bool) function test_obj_func_rotation_qn() bind(C)
        !
        ! this function tests the rotation QN objective function, which evaluates the
        ! host's objective function with the host context of the level without an
        ! objective function of the settings, and otherwise evaluates that with the host
        ! context of the settings, adds the pair of the trial step and the gradient
        ! difference at the trial point, transported back along the reversed step, to
        ! the history where a pair can be formed, remembers its step and rebuilds the
        ! model, passing on the errors of the functions it calls and of the build
        !
        use otr_common_mock, only: mock_obj_func
        use otr_common_unit_tests, only: level_context
        use otr_rotation_qn, only: error_rotation_qn
        use test_reference, only: host_context, host_context_reached
        use, intrinsic :: ieee_arithmetic, only: ieee_class, ieee_positive_inf, &
                                                 operator(==)

        type(recording_build_type), target :: object
        real(rp) :: kappa(n_param), func
        real(rp), allocatable :: kappa_list(:, :), grad_diff_list(:, :), last_grad(:)
        integer(ip) :: error, i_case, i_call, i_kind
        logical :: pair_expected
        character(len=8) :: call_number
        character(len=44), parameter :: case_names(5) = &
            [character(len=44) :: "a trial step", "an energy that is not a number", &
             "a vanishing step", "the starting point", &
             "a history without room for pairs"]

        ! assume tests pass
        test_obj_func_rotation_qn = .true.

        ! evaluate an object without an objective function of the settings, which has
        ! to forward to the host's objective function with the host context of the
        ! level and leave the history unchanged
        call setup_updating_object(object)
        object%solver_level%host%obj_func => mock_obj_func
        call random_number(kappa)
        func = object%obj_func(kappa, error, level_context(object, 1_ip))
        if (error /= 0 .or. .not. is_close(func, sum(kappa))) then
            write(stderr, *) "test_obj_func_rotation_qn failed: Host's objective "// &
                "function not evaluated without an objective function of the settings."
            test_obj_func_rotation_qn = .false.
        end if
        if (level_host_context%n_calls /= 1 .or. host_context%n_calls /= 0 .or. &
            object%n_points /= 2 .or. object%n_built /= 0) then
            write(stderr, *) "test_obj_func_rotation_qn failed: Host's objective "// &
                "function not called with the host context of the level or history "// &
                "changed without an objective function of the settings."
            test_obj_func_rotation_qn = .false.
        end if
        if (.not. host_context_reached("obj_func_rotation_qn", level_host_context)) &
            test_obj_func_rotation_qn = .false.

        ! evaluate an object with an objective function of the settings at a trial
        ! step, which has to add the pair, whose gradient the mock transport moves back
        ! by adding the reversed step, at an energy that is not a number, a vanishing
        ! step, the starting point and for a history without room for pairs, which have
        ! to leave the history unchanged
        do i_case = 1, size(case_names)
            call setup_updating_object(object)
            object%settings%obj_func_grad => mock_obj_func_grad
            call random_number(kappa)
            select case (i_case)
            case (2)
                host_context%nan_output = 1
            case (3)
                kappa = 0.0_rp
            case (4)
                object%n_points = 0
                deallocate(object%last_grad)
            case (5)
                object%settings%max_points = 0
                object%n_points = 0
                deallocate(object%kappa_list, object%grad_diff_list)
                allocate(object%kappa_list(n_param, 0), &
                         object%grad_diff_list(n_param, 0))
            end select
            kappa_list = object%kappa_list(:, :object%n_points)
            grad_diff_list = object%grad_diff_list(:, :object%n_points)
            if (allocated(object%last_grad)) last_grad = object%last_grad
            func = object%obj_func(kappa, error, level_context(object, 1_ip))
            pair_expected = i_case == 1
            if (error /= 0) then
                write(stderr, *) "test_obj_func_rotation_qn failed: Produced error "// &
                    "for "//trim(case_names(i_case))//"."
                test_obj_func_rotation_qn = .false.
            end if
            if (i_case == 2) then
                if (.not. ieee_class(func) == ieee_positive_inf) then
                    write(stderr, *) "test_obj_func_rotation_qn failed: Energy not "// &
                        "returned as infinite for "//trim(case_names(i_case))//"."
                    test_obj_func_rotation_qn = .false.
                end if
            else if (.not. is_close(func, sum(kappa))) then
                write(stderr, *) "test_obj_func_rotation_qn failed: Energy not "// &
                    "returned for "//trim(case_names(i_case))//"."
                test_obj_func_rotation_qn = .false.
            end if
            if (level_host_context%n_calls /= 0 .or. &
                host_context%n_calls /= merge(2, 1, pair_expected)) then
                write(stderr, *) "test_obj_func_rotation_qn failed: Objective and "// &
                    "transport functions not called exactly where needed with the "// &
                    "host context of the settings for "//trim(case_names(i_case))//"."
                test_obj_func_rotation_qn = .false.
            end if
            if (pair_expected) then
                if (object%n_points /= size(kappa_list, 2) + 1) then
                    write(stderr, *) "test_obj_func_rotation_qn failed: Pair not "// &
                        "added for "//trim(case_names(i_case))//"."
                    test_obj_func_rotation_qn = .false.
                else if (.not. all(is_close( &
                    object%kappa_list(:, :object%n_points), &
                    reshape([kappa_list, kappa], [n_param, object%n_points]))) .or. &
                         .not. all(is_close( &
                             object%grad_diff_list(:, :object%n_points), reshape( &
                                 [grad_diff_list, 2.0_rp * kappa - kappa - last_grad], &
                                 [n_param, object%n_points])))) then
                    write(stderr, *) "test_obj_func_rotation_qn failed: History "// &
                        "wrong for "//trim(case_names(i_case))//"."
                    test_obj_func_rotation_qn = .false.
                end if
                if (object%n_built /= 1 .or. object%n_points_built /= object%n_points) &
                    then
                    write(stderr, *) "test_obj_func_rotation_qn failed: Model not "// &
                        "built once for the new history for "// &
                        trim(case_names(i_case))//"."
                    test_obj_func_rotation_qn = .false.
                end if
            else if (object%n_points /= size(kappa_list, 2) .or. object%n_built /= 0) &
                then
                write(stderr, *) "test_obj_func_rotation_qn failed: History "// &
                    "changed for "//trim(case_names(i_case))//"."
                test_obj_func_rotation_qn = .false.
            end if
        end do

        ! let the objective function and the transport function fail in turn, by
        ! returning an error, which the routine has to pass on, and a value that is not
        ! a number, in the gradient for the objective function, which it has to report
        ! with the origin of the object, and let the build of the model fail, whose
        ! error it has to pass on
        do i_call = 1, 2
            write(call_number, '(I0)') i_call
            do i_kind = 1, 2
                call setup_updating_object(object)
                object%settings%obj_func_grad => mock_obj_func_grad
                host_context%failing_call = i_call
                if (i_kind == 1) host_context%mock_error = 7
                if (i_kind == 2) host_context%nan_output = merge(2, 1, i_call == 1)
                func = object%obj_func(kappa, error, level_context(object, 1_ip))
                if (i_kind == 1 .and. error /= 7) then
                    write(stderr, *) "test_obj_func_rotation_qn failed: Error of "// &
                        "call "//trim(call_number)//" not passed on."
                    test_obj_func_rotation_qn = .false.
                else if (i_kind == 2 .and. error /= error_rotation_qn + 1) then
                    write(stderr, *) "test_obj_func_rotation_qn failed: Value that "// &
                        "is not a number at call "//trim(call_number)// &
                        " not reported with the origin of the object."
                    test_obj_func_rotation_qn = .false.
                end if
            end do
        end do
        call setup_updating_object(object)
        object%settings%obj_func_grad => mock_obj_func_grad
        object%error = 5
        func = object%obj_func(kappa, error, level_context(object, 1_ip))
        if (error /= 5) then
            write(stderr, *) "test_obj_func_rotation_qn failed: Error of the build "// &
                "of the model not passed on."
            test_obj_func_rotation_qn = .false.
        end if

    end function test_obj_func_rotation_qn

    logical(c_bool) function test_add_rotation_qn() bind(C)
        !
        ! this function tests the subroutine which adds a new point to a rotation QN
        ! object, which transports the history along the step and, unless the settings
        ! have an objective function, which adds the pairs itself, the last gradient,
        ! appends the pair of the step and the gradient difference, takes the gradient
        ! as the last one and builds the model, and passes on the errors of the
        ! transport function and of the build
        !
        use otr_rotation_qn, only: error_rotation_qn
        use test_reference, only: host_context, host_context_reached

        type(recording_build_type) :: object
        real(rp) :: kappa(n_param), grad(n_param)
        real(rp), allocatable :: kappa_list(:, :), grad_diff_list(:, :), last_grad(:), &
                                 expected_kappa_list(:, :), &
                                 expected_grad_diff_list(:, :)
        integer(ip) :: error, n_points, n_calls, i_case, i_call, i_kind
        logical :: transported, appended
        character(len=40), parameter :: case_names(4) = &
            [character(len=40) :: "a history with room", "the starting point", &
             "a history without room for pairs", "settings with an objective function"]

        ! assume tests pass
        test_add_rotation_qn = .true.

        ! add a point to a history with room, at the starting point, which has no last
        ! gradient, to a history without room for pairs, where neither of the last two
        ! transports or appends anything, and with settings with an objective function,
        ! which only transports the history, and check that every vector of the history
        ! and the last gradient of an appended pair are transported along the step,
        ! which the mock transport adds to them, while the step itself is not, that the
        ! gradient becomes the last one and that the model is built once for the new
        ! history
        call random_number(kappa)
        call random_number(grad)
        do i_case = 1, size(case_names)
            call setup_updating_object(object)
            select case (i_case)
            case (2)
                object%n_points = 0
                deallocate(object%last_grad)
            case (3)
                object%settings%max_points = 0
                object%n_points = 0
                deallocate(object%kappa_list, object%grad_diff_list)
                allocate(object%kappa_list(n_param, 0), &
                         object%grad_diff_list(n_param, 0))
            case (4)
                object%settings%obj_func_grad => mock_obj_func_grad
            end select
            n_points = object%n_points
            kappa_list = object%kappa_list(:, :n_points)
            grad_diff_list = object%grad_diff_list(:, :n_points)
            if (i_case /= 2) last_grad = object%last_grad
            call object%add(kappa, grad, error)

            ! history the addition has to leave
            transported = i_case == 1 .or. i_case == 4
            appended = i_case == 1
            expected_kappa_list = kappa_list
            expected_grad_diff_list = grad_diff_list
            if (transported) then
                expected_kappa_list = kappa_list + spread(kappa, 2, n_points)
                expected_grad_diff_list = grad_diff_list + spread(kappa, 2, n_points)
            end if
            if (appended) then
                expected_kappa_list = reshape([expected_kappa_list, kappa], &
                                              [n_param, n_points + 1])
                expected_grad_diff_list = reshape( &
                    [expected_grad_diff_list, grad - (last_grad + kappa)], &
                    [n_param, n_points + 1])
            end if

            if (error /= 0) then
                write(stderr, *) "test_add_rotation_qn failed: Produced error for "// &
                    trim(case_names(i_case))//"."
                test_add_rotation_qn = .false.
            end if
            if (object%n_points /= size(expected_kappa_list, 2)) then
                write(stderr, *) "test_add_rotation_qn failed: Number of pairs "// &
                    "wrong for "//trim(case_names(i_case))//"."
                test_add_rotation_qn = .false.
            else if (.not. all(is_close(object%kappa_list(:, :object%n_points), &
                                        expected_kappa_list)) .or. &
                     .not. all(is_close(object%grad_diff_list(:, :object%n_points), &
                                        expected_grad_diff_list))) then
                write(stderr, *) "test_add_rotation_qn failed: History wrong for "// &
                    trim(case_names(i_case))//"."
                test_add_rotation_qn = .false.
            end if
            if (.not. all(is_close(object%last_grad, grad, 0.0_rp))) then
                write(stderr, *) "test_add_rotation_qn failed: Gradient not the "// &
                    "last one for "//trim(case_names(i_case))//"."
                test_add_rotation_qn = .false.
            end if
            if (object%n_built /= 1 .or. object%n_points_built /= object%n_points) then
                write(stderr, *) "test_add_rotation_qn failed: Model not built "// &
                    "once for the new history for "//trim(case_names(i_case))//"."
                test_add_rotation_qn = .false.
            end if
            if (host_context%n_calls /= &
                merge(2 * n_points, 0, transported) + merge(1, 0, appended)) then
                write(stderr, *) "test_add_rotation_qn failed: Transport function "// &
                    "not called with the host context for every vector of the "// &
                    "history and the last gradient of an appended pair for "// &
                    trim(case_names(i_case))//"."
                test_add_rotation_qn = .false.
            end if
            if (transported) then
                if (.not. host_context_reached("add_rotation_qn")) &
                    test_add_rotation_qn = .false.
            end if
        end do

        ! let every call of the transport function fail in turn, by returning an error,
        ! which the routine has to pass on, and a value that is not a number, which it
        ! has to report with the origin of the object
        call setup_updating_object(object)
        n_calls = 2 * object%n_points + 1
        do i_call = 1, n_calls
            do i_kind = 1, 2
                call setup_updating_object(object)
                host_context%failing_call = i_call
                if (i_kind == 1) host_context%mock_error = 7
                if (i_kind == 2) host_context%nan_output = 1
                call object%add(kappa, grad, error)
                if (i_kind == 1 .and. error /= 7) then
                    write(stderr, '(A, I0, A)') " test_add_rotation_qn failed: "// &
                        "Error of call ", i_call, &
                        " of the transport function not passed on."
                    test_add_rotation_qn = .false.
                else if (i_kind == 2 .and. error /= error_rotation_qn + 1) then
                    write(stderr, '(A, I0, A)') " test_add_rotation_qn failed: "// &
                        "Value that is not a number at call ", i_call, " of the "// &
                        "transport function not reported with the origin of the object."
                    test_add_rotation_qn = .false.
                end if
            end do
        end do

        ! let the build of the model fail, whose error has to be passed on
        call setup_updating_object(object)
        object%error = 5
        call object%add(kappa, grad, error)
        if (error /= 5) then
            write(stderr, *) "test_add_rotation_qn failed: Error of the build of "// &
                "the model not passed on."
            test_add_rotation_qn = .false.
        end if

    end function test_add_rotation_qn

    logical(c_bool) function test_add_pair_rotation_qn() bind(C)
        !
        ! this function tests the subroutine which appends a pair to the history of a
        ! rotation QN object, dropping the oldest pair of a full history
        !
        type(sr1_updating_type) :: object
        real(rp) :: kappa(n_param), grad_diff(n_param)
        real(rp), allocatable :: kappa_list(:, :), grad_diff_list(:, :)
        integer(ip) :: n_points, n_kept, i_case
        character(len=19), parameter :: case_names(2) = &
            [character(len=19) :: "a history with room", "a full history"]

        ! assume tests pass
        test_add_pair_rotation_qn = .true.

        ! append a random pair to a history with room and to a full history, whose
        ! oldest pair has to be dropped, while the other pairs keep their order
        call random_number(kappa)
        call random_number(grad_diff)
        do i_case = 1, size(case_names)
            call setup_updating_object(object)
            if (i_case == 2) then
                object%n_points = object%settings%max_points
                call random_number(object%kappa_list)
                call random_number(object%grad_diff_list)
            end if
            n_points = object%n_points
            kappa_list = object%kappa_list(:, :n_points)
            grad_diff_list = object%grad_diff_list(:, :n_points)
            call object%add_pair(kappa, grad_diff)
            n_kept = min(n_points, object%settings%max_points - 1)
            if (object%n_points /= n_kept + 1) then
                write(stderr, *) "test_add_pair_rotation_qn failed: Number of "// &
                    "pairs wrong for "//trim(case_names(i_case))//"."
                test_add_pair_rotation_qn = .false.
            else if ( &
                .not. all(is_close(object%kappa_list(:, :object%n_points), reshape( &
                    [kappa_list(:, n_points - n_kept + 1:), kappa], &
                    [n_param, n_kept + 1]), 0.0_rp)) .or. .not. all(is_close( &
                        object%grad_diff_list(:, :object%n_points), reshape( &
                            [grad_diff_list(:, n_points - n_kept + 1:), grad_diff], &
                            [n_param, n_kept + 1]), 0.0_rp))) then
                write(stderr, *) "test_add_pair_rotation_qn failed: History wrong "// &
                    "for "//trim(case_names(i_case))//"."
                test_add_pair_rotation_qn = .false.
            end if
        end do

    end function test_add_pair_rotation_qn

    logical(c_bool) function test_init_rotation_qn() bind(C)
        !
        ! this function tests the subroutine which initializes a rotation QN object
        !
        use otr_rotation_qn, only: rotation_qn_settings_type, sr1_updating_type, &
                                   error_rotation_qn
        use otr_rotation_qn_test_reference, only: ref_rotation_qn_settings, operator(/=)
        use otr_common_test_reference, only: same_logger_and_context
        use opentrustregion_unit_tests, only: test_context_type, logger

        type(sr1_updating_type) :: object
        type(rotation_qn_settings_type) :: settings
        type(test_context_type), target :: context

        ! assume tests pass
        test_init_rotation_qn = .true.

        ! initialize a new object, which has to perform the orbital update, the
        ! objective function, the Hessian linear transformation and the preconditioner
        ! itself, take over the settings and the number of parameters and allocate the
        ! arrays of the history and the model
        settings = ref_rotation_qn_settings
        settings%logger => logger
        settings%context => context
        call object%init(n_param, settings)
        if (.not. all([object%provides%update_orbs, object%provides%obj_func, &
                       object%provides%hess_x, object%provides%precond]) .or. &
            any([object%provides%precond_pd, object%provides%project, &
                 object%provides%get_extra_trial_vectors])) then
            write(stderr, *) "test_init_rotation_qn failed: Roles the object "// &
                "performs set wrong."
            test_init_rotation_qn = .false.
        end if
        if (object%error_origin /= error_rotation_qn) then
            write(stderr, *) "test_init_rotation_qn failed: Error origin of the "// &
                "object not set."
            test_init_rotation_qn = .false.
        end if
        if (object%n_param /= n_param .or. object%settings /= settings .or. &
            .not. same_logger_and_context(object%settings, settings)) then
            write(stderr, *) "test_init_rotation_qn failed: Settings or number of "// &
                "parameters not taken over."
            test_init_rotation_qn = .false.
        end if
        if (any(shape(object%kappa_list) /= [n_param, settings%max_points]) .or. &
            any(shape(object%grad_diff_list) /= [n_param, settings%max_points])) then
            write(stderr, *) "test_init_rotation_qn failed: Arrays of the history "// &
                "allocated wrong."
            test_init_rotation_qn = .false.
        end if
        if (any(shape(object%vectors) /= [n_param, 2 * settings%max_points]) .or. &
            any(shape(object%coeffs) /= [2 * settings%max_points])) then
            write(stderr, *) "test_init_rotation_qn failed: Arrays of the model, "// &
                "which hold two terms per pair, allocated wrong."
            test_init_rotation_qn = .false.
        end if

    end function test_init_rotation_qn

    logical(c_bool) function test_model_x_rotation_qn() bind(C)
        !
        ! this function tests the subroutine which applies the initial Hessian of the
        ! settings or, without one, the Hessian diagonal, updated by the first terms of
        ! the model of a rotation QN object, to a vector
        !
        use otr_rotation_qn, only: sr1_updating_type, error_rotation_qn
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type) :: object
        real(rp) :: x(n_param), model_x(n_param), expected(n_param)
        integer(ip) :: error, n_terms, k, i_kind
        character(len=8) :: terms

        ! assume tests pass
        test_model_x_rotation_qn = .true.

        ! apply the model of an object with three random terms of either sign, updated
        ! by none, some or all of them, which has to apply the initial Hessian once with
        ! the host context of the solver level and add exactly the first terms
        call setup_updating_object(object)
        call generate_random_terms(object)
        call random_number(x)
        do n_terms = 0, object%n_terms
            write(terms, '(I0)') n_terms
            host_context%n_calls = 0
            call object%model_x(x, model_x, n_terms, error)
            expected = init_hess_diag * x
            do k = 1, n_terms
                expected = expected + object%coeffs(k) * &
                           dot_product(object%vectors(:, k), x) * object%vectors(:, k)
            end do
            if (error /= 0) then
                write(stderr, *) "test_model_x_rotation_qn failed: Produced error "// &
                    "for "//trim(terms)//" terms."
                test_model_x_rotation_qn = .false.
            end if
            if (.not. all(is_close(model_x, expected))) then
                write(stderr, *) "test_model_x_rotation_qn failed: Model applied "// &
                    "wrong for "//trim(terms)//" terms."
                test_model_x_rotation_qn = .false.
            end if
            if (host_context%n_calls /= 1) then
                write(stderr, *) "test_model_x_rotation_qn failed: Initial Hessian "// &
                    "function not called once with the host context for "// &
                    trim(terms)//" terms."
                test_model_x_rotation_qn = .false.
            end if
        end do
        if (.not. host_context_reached("model_x_rotation_qn", host_context)) &
            test_model_x_rotation_qn = .false.

        ! apply the model without an initial Hessian of the settings, which has to take
        ! the Hessian diagonal as the initial Hessian
        object%settings%init_hess => null()
        call object%model_x(x, model_x, object%n_terms, error)
        expected = object%h_diag * x
        do k = 1, object%n_terms
            expected = expected + object%coeffs(k) * &
                       dot_product(object%vectors(:, k), x) * object%vectors(:, k)
        end do
        if (error /= 0 .or. .not. all(is_close(model_x, expected))) then
            write(stderr, *) "test_model_x_rotation_qn failed: Hessian diagonal "// &
                "not taken as the initial Hessian without one of the settings."
            test_model_x_rotation_qn = .false.
        end if

        ! let the initial Hessian fail, by returning an error, which the
        ! routine has to pass on, and a value that is not a number, which it has to
        ! report with the origin of the object
        do i_kind = 1, 2
            call setup_updating_object(object)
            if (i_kind == 1) host_context%mock_error = 7
            if (i_kind == 2) host_context%nan_output = 1
            call object%model_x(x, model_x, 0_ip, error)
            if (i_kind == 1 .and. error /= 7) then
                write(stderr, *) "test_model_x_rotation_qn failed: Error of the "// &
                    "initial Hessian not passed on."
                test_model_x_rotation_qn = .false.
            else if (i_kind == 2 .and. error /= error_rotation_qn + 1) then
                write(stderr, *) "test_model_x_rotation_qn failed: Value that is "// &
                    "not a number of the initial Hessian not reported with the "// &
                    "origin of the object."
                test_model_x_rotation_qn = .false.
            end if
        end do

    end function test_model_x_rotation_qn

    logical(c_bool) function test_hess_x_rotation_qn() bind(C)
        !
        ! this function tests the Hessian linear transformation of a rotation QN
        ! object, which applies the initial Hessian updated by every term of its model
        !
        use otr_common_unit_tests, only: level_context
        use otr_rotation_qn, only: sr1_updating_type, error_rotation_qn
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type), target :: object
        real(rp) :: x(n_param), hess_x(n_param), expected(n_param)
        integer(ip) :: error, k, i_kind

        ! assume tests pass
        test_hess_x_rotation_qn = .true.

        ! apply the transformation of an object with three random terms of either sign,
        ! which has to apply the initial Hessian once with the host context of the
        ! solver level and add every term
        call setup_updating_object(object)
        call generate_random_terms(object)
        call random_number(x)
        call object%hess_x(x, hess_x, error, level_context(object, 1_ip))
        expected = init_hess_diag * x
        do k = 1, object%n_terms
            expected = expected + object%coeffs(k) * &
                       dot_product(object%vectors(:, k), x) * object%vectors(:, k)
        end do
        if (error /= 0) then
            write(stderr, *) "test_hess_x_rotation_qn failed: Produced error."
            test_hess_x_rotation_qn = .false.
        end if
        if (.not. all(is_close(hess_x, expected))) then
            write(stderr, *) "test_hess_x_rotation_qn failed: Model not applied "// &
                "with every term."
            test_hess_x_rotation_qn = .false.
        end if
        if (host_context%n_calls /= 1) then
            write(stderr, *) "test_hess_x_rotation_qn failed: Initial Hessian "// &
                "function not called once with the host context."
            test_hess_x_rotation_qn = .false.
        end if
        if (.not. host_context_reached("hess_x_rotation_qn", host_context)) &
            test_hess_x_rotation_qn = .false.

        ! let the initial Hessian fail, by returning an error, which the
        ! routine has to pass on, and a value that is not a number, which it has to
        ! report with the origin of the object
        do i_kind = 1, 2
            call setup_updating_object(object)
            if (i_kind == 1) host_context%mock_error = 7
            if (i_kind == 2) host_context%nan_output = 1
            call object%hess_x(x, hess_x, error, level_context(object, 1_ip))
            if (i_kind == 1 .and. error /= 7) then
                write(stderr, *) "test_hess_x_rotation_qn failed: Error of the "// &
                    "initial Hessian not passed on."
                test_hess_x_rotation_qn = .false.
            else if (i_kind == 2 .and. error /= error_rotation_qn + 1) then
                write(stderr, *) "test_hess_x_rotation_qn failed: Value that is "// &
                    "not a number of the initial Hessian not reported with the "// &
                    "origin of the object."
                test_hess_x_rotation_qn = .false.
            end if
        end do

    end function test_hess_x_rotation_qn

    logical(c_bool) function test_init_precond_rotation_qn() bind(C)
        !
        ! this function tests the subroutine which applies the inverse of the
        ! level-shifted initial Hessian of a rotation QN object, the preconditioner of
        ! the settings or, without one, the level-shifted diagonal preconditioner of the
        ! Hessian diagonal followed by the host's projection of the level
        !
        use opentrustregion_core, only: precond_floor
        use opentrustregion_unit_tests, only: test_context_type, setup_error_logging
        use otr_common_unit_tests, only: level_context
        use otr_rotation_qn, only: sr1_updating_type, error_rotation_qn
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type), target :: object
        type(test_context_type), target :: log_context
        real(rp) :: x(n_param), precond_x(n_param), expected(n_param), divisors(n_param)
        real(rp), parameter :: mu = 0.3_rp
        integer(ip) :: error, i_kind

        ! assume tests pass
        test_init_precond_rotation_qn = .true.

        ! apply the preconditioner of the settings with their host context
        call setup_updating_object(object)
        object%settings%init_precond => mock_init_precond
        call random_number(x)
        call object%init_precond(x, mu, precond_x, error, level_context(object, 1_ip))
        if (error /= 0 .or. .not. all(is_close(precond_x, x / (init_hess_diag - mu)))) &
            then
            write(stderr, *) "test_init_precond_rotation_qn failed: Preconditioner "// &
                "of the settings not applied."
            test_init_precond_rotation_qn = .false.
        end if
        if (host_context%n_calls /= 1 .or. level_host_context%n_calls /= 0) then
            write(stderr, *) "test_init_precond_rotation_qn failed: Preconditioner "// &
                "of the settings not called once with their host context."
            test_init_precond_rotation_qn = .false.
        end if
        if (.not. host_context_reached("init_precond_rotation_qn")) &
            test_init_precond_rotation_qn = .false.

        ! apply the level-shifted diagonal preconditioner without a preconditioner of
        ! the settings, with a Hessian diagonal element at the level shift, whose
        ! divisor has to be floored, followed by the host's projection with the host
        ! context of the level
        call setup_updating_object(object)
        object%solver_level%host%project => mock_project
        object%h_diag(2) = mu
        call object%init_precond(x, mu, precond_x, error, level_context(object, 1_ip))
        divisors = object%h_diag - mu
        divisors(2) = precond_floor
        expected = x / divisors
        expected(1) = 0.0_rp
        if (error /= 0 .or. .not. all(is_close(precond_x, expected))) then
            write(stderr, *) "test_init_precond_rotation_qn failed: Diagonal "// &
                "preconditioner and host's projection not applied without a "// &
                "preconditioner of the settings."
            test_init_precond_rotation_qn = .false.
        end if
        if (level_host_context%n_calls /= 1) then
            write(stderr, *) "test_init_precond_rotation_qn failed: Host's "// &
                "projection not called once with the host context of the level."
            test_init_precond_rotation_qn = .false.
        end if

        ! let the preconditioner of the settings and the host's projection fail, whose
        ! errors have to be passed on
        do i_kind = 1, 2
            call setup_updating_object(object)
            if (i_kind == 1) then
                object%settings%init_precond => mock_init_precond
                host_context%mock_error = 7
            else
                object%solver_level%host%project => mock_project
                level_host_context%mock_error = 7
            end if
            call object%init_precond(x, mu, precond_x, error, &
                                     level_context(object, 1_ip))
            if (error /= 7) then
                write(stderr, *) "test_init_precond_rotation_qn failed: Error of "// &
                    "the "//trim(merge("preconditioner of the settings", &
                                       "host's projection             ", &
                                       i_kind == 1))//" not passed on."
                test_init_precond_rotation_qn = .false.
            end if
        end do

        ! let the preconditioner of the settings return a value that is not a number,
        ! which has to be reported with the origin of the object and an error message
        call setup_updating_object(object)
        call setup_error_logging(object%settings, log_context)
        object%settings%init_precond => mock_init_precond
        log_context%nan_output = 1
        call object%init_precond(x, mu, precond_x, error, level_context(object, 1_ip))
        if (error /= error_rotation_qn + 1) then
            write(stderr, *) "test_init_precond_rotation_qn failed: Result that is "// &
                "not finite not reported with the origin of the object."
            test_init_precond_rotation_qn = .false.
        end if
        if (len_trim(log_context%log_message) == 0) then
            write(stderr, *) "test_init_precond_rotation_qn failed: No error "// &
                "message printed for a result that is not finite."
            test_init_precond_rotation_qn = .false.
        end if

    end function test_init_precond_rotation_qn

    logical(c_bool) function test_precond_rotation_qn() bind(C)
        !
        ! this function tests the preconditioner of a rotation QN object, which applies
        ! the inverse of its level-shifted model by the Sherman-Morrison-Woodbury
        ! formula with the inverse of the level-shifted initial Hessian, or only that
        ! inverse where the level shift makes the model singular
        !
        use otr_common_unit_tests, only: level_context
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type), target :: object
        real(rp) :: residual(n_param), precond_residual(n_param), &
                    shifted_model(n_param, n_param)
        real(rp), parameter :: mu = 0.3_rp
        integer(ip) :: error, i, n_calls, i_call

        ! assume tests pass
        test_precond_rotation_qn = .true.

        ! apply the preconditioner of an object with three random terms and the exact
        ! inverse of the level-shifted initial Hessian as the preconditioner of the
        ! settings, whose result the level-shifted model has to map back to the
        ! residual, which has to call the preconditioner of the settings for the
        ! residual and every term with their host context
        call setup_updating_object(object)
        call generate_random_terms(object)
        object%settings%init_precond => mock_init_precond
        call random_number(residual)
        call object%precond(residual, mu, precond_residual, error, &
                            level_context(object, 1_ip))
        n_calls = host_context%n_calls
        shifted_model = model_matrix(object)
        do i = 1, n_param
            shifted_model(i, i) = shifted_model(i, i) - mu
        end do
        if (error /= 0) then
            write(stderr, *) "test_precond_rotation_qn failed: Produced error."
            test_precond_rotation_qn = .false.
        end if
        if (.not. all(is_close(matmul(shifted_model, precond_residual), residual))) then
            write(stderr, *) "test_precond_rotation_qn failed: Inverse of the "// &
                "level-shifted model not applied."
            test_precond_rotation_qn = .false.
        end if
        if (n_calls /= object%n_terms + 1) then
            write(stderr, *) "test_precond_rotation_qn failed: Preconditioner of "// &
                "the settings not called with their host context for the residual "// &
                "and every term."
            test_precond_rotation_qn = .false.
        end if
        if (.not. host_context_reached("precond_rotation_qn", host_context)) &
            test_precond_rotation_qn = .false.

        ! apply the preconditioner of an object without terms and of an object whose
        ! single term makes the level-shifted model singular, which have to apply only
        ! the inverse of the level-shifted initial Hessian
        call setup_updating_object(object)
        object%settings%init_precond => mock_init_precond
        call object%precond(residual, mu, precond_residual, error, &
                            level_context(object, 1_ip))
        if (error /= 0 .or. .not. all(is_close(precond_residual, &
                                               residual / (init_hess_diag - mu)))) then
            write(stderr, *) "test_precond_rotation_qn failed: Inverse of the "// &
                "level-shifted initial Hessian not applied without terms."
            test_precond_rotation_qn = .false.
        end if
        object%n_terms = 1
        object%vectors(:, 1) = 0.0_rp
        object%vectors(1, 1) = 1.0_rp
        object%coeffs(1) = mu - init_hess_diag(1)
        call object%precond(residual, mu, precond_residual, error, &
                            level_context(object, 1_ip))
        if (error /= 0 .or. .not. all(is_close(precond_residual, &
                                               residual / (init_hess_diag - mu)))) then
            write(stderr, *) "test_precond_rotation_qn failed: Inverse of the "// &
                "level-shifted initial Hessian not applied for a singular "// &
                "level-shifted model."
            test_precond_rotation_qn = .false.
        end if

        ! let every call of the preconditioner of the settings fail in turn, whose
        ! error has to be passed on
        do i_call = 1, n_calls
            call setup_updating_object(object)
            call generate_random_terms(object)
            object%settings%init_precond => mock_init_precond
            host_context%failing_call = i_call
            host_context%mock_error = 7
            call object%precond(residual, mu, precond_residual, error, &
                                level_context(object, 1_ip))
            if (error /= 7) then
                write(stderr, '(A, I0, A)') " test_precond_rotation_qn failed: "// &
                    "Error of call ", i_call, &
                    " of the host's preconditioner not passed on."
                test_precond_rotation_qn = .false.
            end if
        end do

    end function test_precond_rotation_qn

    logical(c_bool) function test_keep_pair() bind(C)
        !
        ! this function tests the subroutine which keeps a pair of the history of a
        ! rotation QN object as the next kept one
        !
        use otr_rotation_qn, only: sr1_updating_type, keep_pair

        type(sr1_updating_type) :: object
        real(rp), allocatable :: kappa_list(:, :), grad_diff_list(:, :)
        integer(ip) :: n_kept

        ! assume tests pass
        test_keep_pair = .true.

        ! keep the third of three random pairs as the second, which has to be copied
        ! there and counted, while the other pairs stay unchanged
        call setup_updating_object(object)
        object%n_points = 3
        call random_number(object%kappa_list(:, :3))
        call random_number(object%grad_diff_list(:, :3))
        kappa_list = object%kappa_list(:, :3)
        grad_diff_list = object%grad_diff_list(:, :3)
        n_kept = 1
        call keep_pair(object, 3_ip, n_kept)
        if (n_kept /= 2) then
            write(stderr, *) "test_keep_pair failed: Kept pair not counted."
            test_keep_pair = .false.
        end if
        if (.not. all(is_close(object%kappa_list(:, 2), kappa_list(:, 3), &
                               0.0_rp)) .or. &
            .not. all(is_close(object%grad_diff_list(:, 2), grad_diff_list(:, 3), &
                               0.0_rp))) then
            write(stderr, *) "test_keep_pair failed: Pair not kept as the next one."
            test_keep_pair = .false.
        end if
        if (.not. all(is_close(object%kappa_list(:, [1, 3]), kappa_list(:, [1, 3]), &
                               0.0_rp)) .or. &
            .not. all(is_close(object%grad_diff_list(:, [1, 3]), &
                               grad_diff_list(:, [1, 3]), 0.0_rp))) then
            write(stderr, *) "test_keep_pair failed: Other pairs changed."
            test_keep_pair = .false.
        end if

    end function test_keep_pair

    logical(c_bool) function test_build_model_sr1() bind(C)
        !
        ! this function tests the subroutine which builds the model of the SR1
        ! updating scheme from the history, by comparing with the recursive SR1 update
        ! of a dense matrix
        !
        use otr_rotation_qn, only: sr1_updating_type, error_rotation_qn, sr1_skip_thresh
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type) :: object
        real(rp) :: kappa_list(n_param, 4), grad_diff_list(n_param, 4), &
                    orthogonal(n_param), model(n_param, n_param), &
                    expected(n_param, n_param)
        integer(ip) :: error, n_calls, i_call, i_kind
        integer(ip), parameter :: kept(3) = [1, 2, 4]
        character(len=8) :: call_number

        ! assume tests pass
        test_build_model_sr1 = .true.

        ! build the model of a history of four random pairs, whose gradient
        ! differences, five times the steps plus noise, which keeps the denominators
        ! well away from zero, are not those of a quadratic function, so that only the
        ! recursive update reproduces the newest one, while the third has a residual
        ! orthogonal to its step, so that it has to be skipped and dropped from the
        ! history, whose other pairs keep their order
        call random_number(kappa_list)
        call random_number(grad_diff_list)
        grad_diff_list = 5.0_rp * kappa_list + grad_diff_list - 0.5_rp
        call random_number(orthogonal)
        orthogonal = orthogonal - dot_product(orthogonal, kappa_list(:, 3)) / &
                     dot_product(kappa_list(:, 3), kappa_list(:, 3)) * kappa_list(:, 3)
        model = ref_sr1_model(kappa_list(:, :2), grad_diff_list(:, :2))
        grad_diff_list(:, 3) = matmul(model, kappa_list(:, 3)) + orthogonal
        call setup_updating_object(object)
        object%n_points = 4
        object%kappa_list(:, :4) = kappa_list
        object%grad_diff_list(:, :4) = grad_diff_list
        call object%build_model(error)
        n_calls = host_context%n_calls
        if (error /= 0) then
            write(stderr, *) "test_build_model_sr1 failed: Produced error."
            test_build_model_sr1 = .false.
        end if
        if (object%n_points /= 3 .or. object%n_terms /= 3) then
            write(stderr, *) "test_build_model_sr1 failed: Pair with a residual "// &
                "orthogonal to its step not dropped."
            test_build_model_sr1 = .false.
        else if (.not. all(is_close(object%kappa_list(:, :3), kappa_list(:, kept), &
                                    0.0_rp)) .or. &
                 .not. all(is_close(object%grad_diff_list(:, :3), &
                                    grad_diff_list(:, kept), 0.0_rp))) then
            write(stderr, *) "test_build_model_sr1 failed: Kept pairs wrong."
            test_build_model_sr1 = .false.
        end if
        model = model_matrix(object)
        expected = ref_sr1_model(kappa_list, grad_diff_list)
        if (.not. all(is_close(model, expected))) then
            write(stderr, *) "test_build_model_sr1 failed: Model differs from the "// &
                "recursive SR1 update."
            test_build_model_sr1 = .false.
        end if
        if (.not. all(is_close(matmul(model, kappa_list(:, 4)), &
                               grad_diff_list(:, 4)))) then
            write(stderr, *) "test_build_model_sr1 failed: Model does not satisfy "// &
                "the secant condition of the newest pair."
            test_build_model_sr1 = .false.
        end if
        if (n_calls /= 4) then
            write(stderr, *) "test_build_model_sr1 failed: Initial Hessian "// &
                "function not called with the host context of the solver level for "// &
                "every pair."
            test_build_model_sr1 = .false.
        end if
        if (.not. host_context_reached("build_model_sr1", host_context)) &
            test_build_model_sr1 = .false.

        ! let every call of the initial Hessian fail in turn, by returning an
        ! error, which the routine has to pass on, and a value that is not a number,
        ! which it has to report with the origin of the object
        do i_call = 1, n_calls
            write(call_number, '(I0)') i_call
            do i_kind = 1, 2
                call setup_updating_object(object)
                object%n_points = 4
                object%kappa_list(:, :4) = kappa_list
                object%grad_diff_list(:, :4) = grad_diff_list
                host_context%failing_call = i_call
                if (i_kind == 1) host_context%mock_error = 7
                if (i_kind == 2) host_context%nan_output = 1
                call object%build_model(error)
                if (i_kind == 1 .and. error /= 7) then
                    write(stderr, *) "test_build_model_sr1 failed: Error of call "// &
                        trim(call_number)//" of the initial Hessian not passed on."
                    test_build_model_sr1 = .false.
                else if (i_kind == 2 .and. error /= error_rotation_qn + 1) then
                    write(stderr, *) "test_build_model_sr1 failed: Value that is "// &
                        "not a number at call "//trim(call_number)//" of the "// &
                        "initial Hessian not reported with the origin of the object."
                    test_build_model_sr1 = .false.
                end if
            end do
        end do

    contains

        function ref_sr1_model(steps, grad_diffs) result(matrix)
            !
            ! this function applies the SR1 update of every pair, skipping those whose
            ! denominator is too small, to the dense initial Hessian of the mock
            !
            real(rp), intent(in) :: steps(:, :), grad_diffs(:, :)
            real(rp) :: matrix(n_param, n_param)

            real(rp) :: residual(n_param), denominator
            integer(ip) :: i_param, j

            matrix = 0.0_rp
            do i_param = 1, n_param
                matrix(i_param, i_param) = init_hess_diag(i_param)
            end do
            do j = 1, size(steps, 2)
                residual = grad_diffs(:, j) - matmul(matrix, steps(:, j))
                denominator = dot_product(steps(:, j), residual)
                if (abs(denominator) <= &
                    sr1_skip_thresh * norm2(steps(:, j)) * norm2(residual)) cycle
                matrix = matrix + spread(residual, 2, n_param) * &
                         spread(residual, 1, n_param) / denominator
            end do

        end function ref_sr1_model

    end function test_build_model_sr1

    logical(c_bool) function test_build_model_bfgs() bind(C)
        !
        ! this function tests the subroutine which builds the model of the BFGS
        ! updating scheme from the history, by comparing with the recursive BFGS update
        ! of a dense matrix
        !
        use otr_rotation_qn, only: bfgs_updating_type, error_rotation_qn, &
                                   bfgs_skip_thresh
        use test_reference, only: host_context, host_context_reached

        type(bfgs_updating_type) :: object
        real(rp) :: kappa_list(n_param, 4), grad_diff_list(n_param, 4), &
                    model(n_param, n_param), expected(n_param, n_param)
        integer(ip) :: error, n_calls, i_call, i_kind
        integer(ip), parameter :: kept(2) = [1, 4]
        character(len=8) :: call_number

        ! assume tests pass
        test_build_model_bfgs = .true.

        ! build the model of a history of four pairs: the first with positive curvature
        ! and a step orthogonal to the last axis, along which the initial Hessian of the
        ! mock is negative, the second with negative curvature, the third with positive
        ! curvature along the last axis, along which the model so far is still
        ! negative, and the fourth with positive curvature and a step mostly within the
        ! first two axes, so that the second and the third have to be skipped and
        ! dropped from the history, whose other pairs keep their order
        call random_number(kappa_list)
        kappa_list = kappa_list + 0.5_rp
        kappa_list(n_param, 1) = 0.0_rp
        grad_diff_list(:, 1) = [4.0_rp, 5.0_rp, 6.0_rp] * kappa_list(:, 1)
        grad_diff_list(:, 2) = -kappa_list(:, 2)
        kappa_list(:, 3) = [0.0_rp, 0.0_rp, 1.0_rp]
        grad_diff_list(:, 3) = 2.0_rp * kappa_list(:, 3)
        kappa_list(n_param, 4) = 0.1_rp * kappa_list(n_param, 4)
        call random_number(grad_diff_list(:, 4))
        grad_diff_list(:, 4) = [4.0_rp, 5.0_rp, 6.0_rp] * kappa_list(:, 4) + &
                               0.1_rp * grad_diff_list(:, 4)
        call setup_updating_object(object)
        object%n_points = 4
        object%kappa_list(:, :4) = kappa_list
        object%grad_diff_list(:, :4) = grad_diff_list
        call object%build_model(error)
        n_calls = host_context%n_calls
        if (error /= 0) then
            write(stderr, *) "test_build_model_bfgs failed: Produced error."
            test_build_model_bfgs = .false.
        end if
        if (object%n_points /= 2 .or. object%n_terms /= 4) then
            write(stderr, *) "test_build_model_bfgs failed: Pairs without positive "// &
                "curvature in the history or the model not dropped."
            test_build_model_bfgs = .false.
        else if (.not. all(is_close(object%kappa_list(:, :2), kappa_list(:, kept), &
                                    0.0_rp)) .or. &
                 .not. all(is_close(object%grad_diff_list(:, :2), &
                                    grad_diff_list(:, kept), 0.0_rp))) then
            write(stderr, *) "test_build_model_bfgs failed: Kept pairs wrong."
            test_build_model_bfgs = .false.
        end if
        model = model_matrix(object)
        expected = ref_bfgs_model(kappa_list, grad_diff_list)
        if (.not. all(is_close(model, expected))) then
            write(stderr, *) "test_build_model_bfgs failed: Model differs from the "// &
                "recursive BFGS update."
            test_build_model_bfgs = .false.
        end if
        if (.not. all(is_close(matmul(model, kappa_list(:, 4)), &
                               grad_diff_list(:, 4)))) then
            write(stderr, *) "test_build_model_bfgs failed: Model does not satisfy "// &
                "the secant condition of the newest pair."
            test_build_model_bfgs = .false.
        end if
        if (n_calls /= 3) then
            write(stderr, *) "test_build_model_bfgs failed: Initial Hessian "// &
                "function not called with the host context for every pair with "// &
                "positive curvature."
            test_build_model_bfgs = .false.
        end if
        if (.not. host_context_reached("build_model_bfgs", host_context)) &
            test_build_model_bfgs = .false.

        ! let every call of the initial Hessian fail in turn, by returning an
        ! error, which the routine has to pass on, and a value that is not a number,
        ! which it has to report with the origin of the object
        do i_call = 1, n_calls
            write(call_number, '(I0)') i_call
            do i_kind = 1, 2
                call setup_updating_object(object)
                object%n_points = 4
                object%kappa_list(:, :4) = kappa_list
                object%grad_diff_list(:, :4) = grad_diff_list
                host_context%failing_call = i_call
                if (i_kind == 1) host_context%mock_error = 7
                if (i_kind == 2) host_context%nan_output = 1
                call object%build_model(error)
                if (i_kind == 1 .and. error /= 7) then
                    write(stderr, *) "test_build_model_bfgs failed: Error of call "// &
                        trim(call_number)//" of the initial Hessian not passed on."
                    test_build_model_bfgs = .false.
                else if (i_kind == 2 .and. error /= error_rotation_qn + 1) then
                    write(stderr, *) "test_build_model_bfgs failed: Value that is "// &
                        "not a number at call "//trim(call_number)//" of the "// &
                        "initial Hessian not reported with the origin of the object."
                    test_build_model_bfgs = .false.
                end if
            end do
        end do

    contains

        function ref_bfgs_model(steps, grad_diffs) result(matrix)
            !
            ! this function applies the BFGS update of every pair, skipping those
            ! without sufficiently positive curvature in the pair or in the model so
            ! far, to the dense initial Hessian of the mock
            !
            real(rp), intent(in) :: steps(:, :), grad_diffs(:, :)
            real(rp) :: matrix(n_param, n_param)

            real(rp) :: model_step(n_param), curvature, model_curvature
            integer(ip) :: i_param, j

            matrix = 0.0_rp
            do i_param = 1, n_param
                matrix(i_param, i_param) = init_hess_diag(i_param)
            end do
            do j = 1, size(steps, 2)
                curvature = dot_product(steps(:, j), grad_diffs(:, j))
                if (curvature <= bfgs_skip_thresh * norm2(steps(:, j)) * &
                    norm2(grad_diffs(:, j))) cycle
                model_step = matmul(matrix, steps(:, j))
                model_curvature = dot_product(steps(:, j), model_step)
                if (model_curvature <= &
                    bfgs_skip_thresh * norm2(steps(:, j)) * norm2(model_step)) cycle
                matrix = matrix - spread(model_step, 2, n_param) * &
                         spread(model_step, 1, n_param) / model_curvature + &
                         spread(grad_diffs(:, j), 2, n_param) * &
                         spread(grad_diffs(:, j), 1, n_param) / curvature
            end do

        end function ref_bfgs_model

    end function test_build_model_bfgs

end module otr_rotation_qn_unit_tests
