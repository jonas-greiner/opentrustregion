! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: n_param, tol
    use otr_qn, only: transport_type, init_hess_type, sr1_updating_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(transport_type), pointer :: mock_transport_ptr => mock_transport
    procedure(init_hess_type), pointer :: mock_init_hess_ptr => mock_init_hess

    ! quasi-Newton object which records the point its orbital update adds and returns
    ! the error code it is set to, for the test of the orbital update
    type, extends(sr1_updating_type) :: recording_qn_type
        integer(ip) :: n_added = 0, error = 0
        real(rp), allocatable :: added_kappa(:), added_grad(:)
    contains
        procedure :: add => recording_add_qn
    end type

contains

    subroutine mock_transport(geodesic, tangent_vector, error, context)
        !
        ! this subroutine is a mock of the host's transport function, which leaves the
        ! tangent vector unchanged, records its call, checks the host context and
        ! returns the error code set in it
        !
        use otr_common_mock, only: record_host_call

        real(rp), intent(in), target :: geodesic(:)
        real(rp), intent(inout), target :: tangent_vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        error = record_host_call("transport", context)

    end subroutine mock_transport

    subroutine mock_init_hess(vector, error, context)
        !
        ! this subroutine is a mock of the host's initial Hessian function, which
        ! applies the identity, records its call, checks the host context and returns
        ! the error code set in it
        !
        use otr_common_mock, only: record_host_call

        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        error = record_host_call("init_hess", context)

    end subroutine mock_init_hess

    subroutine recording_add_qn(self, kappa, grad, error)
        !
        ! this subroutine records the step and gradient of the point the recording
        ! quasi-Newton object is given and returns the error code it is set to
        !
        class(recording_qn_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), grad(:)
        integer(ip), intent(out) :: error

        self%n_added = self%n_added + 1
        self%added_kappa = kappa
        self%added_grad = grad
        error = self%error

    end subroutine recording_add_qn

    subroutine setup_updating_object(object)
        !
        ! this subroutine sets up a new updating object as after three orbital updates,
        ! without the routines under test, with the default settings carrying the armed
        ! host context, the mock transport and initial Hessian functions, two random
        ! points in its history, whose gradient differences are twice their steps, so
        ! that they have positive curvature along them, and a random last gradient, so
        ! that every loop over the history calls the host's functions more than once
        !
        use otr_qn, only: updating_type, default_qn_settings
        use test_reference, only: arm_host_context

        class(updating_type), intent(out) :: object

        object%n_param = n_param
        object%settings = default_qn_settings
        call arm_host_context(object%settings)
        object%transport => mock_transport
        object%init_hess => mock_init_hess
        allocate(object%kappa_list(n_param, object%settings%max_points), &
                 object%grad_diff_list(n_param, object%settings%max_points), &
                 object%last_grad(n_param))
        object%kappa_list = 0.0_rp
        object%grad_diff_list = 0.0_rp
        object%n_points = 2
        call random_number(object%kappa_list(:, :2))
        object%grad_diff_list(:, :2) = 2.0_rp * object%kappa_list(:, :2)
        call random_number(object%last_grad)

    end subroutine setup_updating_object

    logical(c_bool) function test_qn_factory() bind(C)
        !
        ! this function tests the factory, which creates a quasi-Newton object of the
        ! updating scheme the settings select after their sanity check, initializes it
        ! with the settings and the host's transport and initial Hessian functions and
        ! wires the solver settings with it and the record of the state
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, solver_state_type, &
                                   default_solver_settings
        use opentrustregion_unit_tests, only: setup_settings, test_context_type
        use otr_common, only: extension_solver_settings_type, extension_level_type, &
                              update_orbs_callback, obj_func_callback
        use otr_common_unit_tests, only: wired_level
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_qn, only: qn_factory, qn_settings_type, updating_type, &
                          sr1_updating_type, bfgs_updating_type, default_qn_settings, &
                          qn_settings_uninitialized_warning_msg
        use otr_qn_test_reference, only: ref_qn_settings, operator(/=)

        procedure(update_orbs_type), pointer :: update_orbs_funptr, &
                                                update_orbs_qn_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr, obj_func_qn_funptr
        procedure(transport_type), pointer :: transport_funptr
        procedure(init_hess_type), pointer :: init_hess_funptr
        type(extension_solver_settings_type), allocatable, target :: solver_settings, &
                                                                     other_settings
        type(solver_state_type), target :: state, other_state
        type(qn_settings_type) :: settings, expected_settings
        type(test_context_type), target :: context
        type(extension_level_type), pointer :: level, other_level
        integer(ip) :: i_case, error
        logical :: sr1_expected, scheme_right, owned
        character(len=41), parameter :: case_names(4) = &
            [character(len=41) :: "for settings selecting SR1 in uppercase", &
             "for settings selecting BFGS in uppercase", "without settings", &
             "for settings that were not initialized"]

        ! assume tests pass
        test_qn_factory = .true.

        ! set the host's callback functions
        update_orbs_funptr => mock_update_orbs
        obj_func_funptr => mock_obj_func
        transport_funptr => mock_transport
        init_hess_funptr => mock_init_hess

        ! create objects for settings selecting either updating scheme, for no settings
        ! and for settings that were not initialized, which are initialized after a
        ! warning through their logging function, each of which has to be the object of
        ! the selected scheme, initialized with the settings after their sanity check
        ! and the host's transport and initial Hessian functions, and to wire the
        ! solver settings, which own it, with the host's orbital update and objective
        ! function and the record of the state
        do i_case = 1, size(case_names)
            allocate(solver_settings)
            solver_settings%solver_settings_type = default_solver_settings
            select case (i_case)
            case (1)
                settings = ref_qn_settings
                settings%hess_update_scheme = "SR1"
                expected_settings = settings
                expected_settings%hess_update_scheme = "sr1"
            case (2)
                settings = ref_qn_settings
                settings%hess_update_scheme = "BFGS"
                expected_settings = ref_qn_settings
            case (3)
                expected_settings = default_qn_settings
            case (4)
                call setup_settings(settings, context)
                settings%initialized = .false.
                expected_settings = default_qn_settings
            end select
            sr1_expected = expected_settings%hess_update_scheme == "sr1"
            if (i_case == 3) then
                call qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                                transport_funptr, init_hess_funptr, &
                                update_orbs_qn_funptr, obj_func_qn_funptr, &
                                solver_settings, state, error)
            else
                call qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                                transport_funptr, init_hess_funptr, &
                                update_orbs_qn_funptr, obj_func_qn_funptr, &
                                solver_settings, state, error, settings)
            end if
            if (error /= 0) then
                write(stderr, *) "test_qn_factory failed: Produced error "// &
                    trim(case_names(i_case))//"."
                test_qn_factory = .false.
                deallocate(solver_settings)
                cycle
            end if
            if (.not. associated(update_orbs_qn_funptr, update_orbs_callback) .or. &
                .not. associated(obj_func_qn_funptr, obj_func_callback)) then
                write(stderr, *) "test_qn_factory failed: Role callback functions "// &
                    "not returned "//trim(case_names(i_case))//"."
                test_qn_factory = .false.
            end if
            level => wired_level(solver_settings)
            if (.not. associated(level)) then
                write(stderr, *) "test_qn_factory failed: Solver settings not "// &
                    "wired with an object "//trim(case_names(i_case))//"."
                test_qn_factory = .false.
                deallocate(solver_settings)
                cycle
            end if
            owned = .false.
            if (allocated(solver_settings%extensions)) then
                if (size(solver_settings%extensions) == 1) owned = &
                    associated(level%extension, solver_settings%extensions(1)%object)
            end if
            if (.not. owned) then
                write(stderr, *) "test_qn_factory failed: Object not owned by the "// &
                    "solver settings "//trim(case_names(i_case))//"."
                test_qn_factory = .false.
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
                write(stderr, *) "test_qn_factory failed: Object of the wrong "// &
                    "updating scheme created "//trim(case_names(i_case))//"."
                test_qn_factory = .false.
            end if
            select type (object => level%extension)
            class is (updating_type)
                if (object%n_param /= n_param .or. &
                    object%settings /= expected_settings) then
                    write(stderr, *) "test_qn_factory failed: Object not "// &
                        "initialized with the settings "//trim(case_names(i_case))//"."
                    test_qn_factory = .false.
                end if
                if (.not. associated(object%transport, mock_transport) .or. &
                    .not. associated(object%init_hess, mock_init_hess)) then
                    write(stderr, *) "test_qn_factory failed: Host's transport or "// &
                        "initial Hessian function not handed to the object "// &
                        trim(case_names(i_case))//"."
                    test_qn_factory = .false.
                end if
            end select
            if (.not. associated(level%host%update_orbs, mock_update_orbs) .or. &
                .not. associated(level%host%obj_func, mock_obj_func) .or. &
                .not. associated(level%extension%state, state)) then
                write(stderr, *) "test_qn_factory failed: Solver settings not "// &
                    "wired with the host's functions and the record of the state "// &
                    trim(case_names(i_case))//"."
                test_qn_factory = .false.
            end if
            if (i_case == 4 .and. index( &
                context%log_message, " "//qn_settings_uninitialized_warning_msg) == 0) &
                then
                write(stderr, *) "test_qn_factory failed: Warning not printed for "// &
                    "settings that were not initialized."
                test_qn_factory = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create an object for parameters the sanity check rejects, which has to return
        ! its error without functions and leave the solver settings untouched
        allocate(solver_settings)
        solver_settings%solver_settings_type = default_solver_settings
        call setup_settings(settings, context)
        update_orbs_qn_funptr => update_orbs_callback
        obj_func_qn_funptr => obj_func_callback
        call qn_factory(0_ip, update_orbs_funptr, obj_func_funptr, transport_funptr, &
                        init_hess_funptr, update_orbs_qn_funptr, obj_func_qn_funptr, &
                        solver_settings, state, error, settings)
        if (error /= 1) then
            write(stderr, *) "test_qn_factory failed: Error not returned for "// &
                "parameters the sanity check rejects."
            test_qn_factory = .false.
        end if
        if (associated(update_orbs_qn_funptr) .or. associated(obj_func_qn_funptr)) then
            write(stderr, *) "test_qn_factory failed: Functions returned for "// &
                "parameters the sanity check rejects."
            test_qn_factory = .false.
        end if
        if (associated(solver_settings%context) .or. &
            allocated(solver_settings%extensions)) then
            write(stderr, *) "test_qn_factory failed: Solver settings changed for "// &
                "parameters the sanity check rejects."
            test_qn_factory = .false.
        end if
        deallocate(solver_settings)

        ! create an object for solver settings a factory wired, given the host's own
        ! orbital update and objective function instead of those the factory returned,
        ! which the wiring refuses, so that the factory has to return its error without
        ! functions and leave the settings wired with the first object only
        allocate(solver_settings)
        solver_settings%solver_settings_type = default_solver_settings
        call setup_settings(settings, context)
        call qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                        transport_funptr, init_hess_funptr, update_orbs_qn_funptr, &
                        obj_func_qn_funptr, solver_settings, state, error, settings)
        level => wired_level(solver_settings)
        call qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                        transport_funptr, init_hess_funptr, update_orbs_qn_funptr, &
                        obj_func_qn_funptr, solver_settings, state, error, settings)
        if (error /= 1) then
            write(stderr, *) "test_qn_factory failed: Error not returned for the "// &
                "host's own functions given with solver settings a factory wired."
            test_qn_factory = .false.
        end if
        if (associated(update_orbs_qn_funptr) .or. associated(obj_func_qn_funptr)) then
            write(stderr, *) "test_qn_factory failed: Functions returned for the "// &
                "host's own functions given with solver settings a factory wired."
            test_qn_factory = .false.
        end if
        owned = .false.
        if (associated(level) .and. allocated(solver_settings%extensions)) &
            owned = size(solver_settings%extensions) == 1 .and. &
                    associated(wired_level(solver_settings), level)
        if (.not. owned) then
            write(stderr, *) "test_qn_factory failed: Solver settings a factory "// &
                "wired changed for the host's own functions."
            test_qn_factory = .false.
        end if
        deallocate(solver_settings)

        ! create objects for two separate solver settings, which have to be distinct
        ! objects reading their own record of the state, so that calculations with
        ! separate settings can be used alternately
        allocate(solver_settings, other_settings)
        solver_settings%solver_settings_type = default_solver_settings
        other_settings%solver_settings_type = default_solver_settings
        settings = ref_qn_settings
        call qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                        transport_funptr, init_hess_funptr, update_orbs_qn_funptr, &
                        obj_func_qn_funptr, solver_settings, state, error, settings)
        if (error /= 0) then
            write(stderr, *) "test_qn_factory failed: Produced error for the first "// &
                "of two separate solver settings."
            test_qn_factory = .false.
        end if
        call qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                        transport_funptr, init_hess_funptr, update_orbs_qn_funptr, &
                        obj_func_qn_funptr, other_settings, other_state, error, &
                        settings)
        if (error /= 0) then
            write(stderr, *) "test_qn_factory failed: Produced error for the "// &
                "second of two separate solver settings."
            test_qn_factory = .false.
        end if
        level => wired_level(solver_settings)
        other_level => wired_level(other_settings)
        if (.not. (associated(level) .and. associated(other_level))) then
            write(stderr, *) "test_qn_factory failed: Separate solver settings not "// &
                "wired."
            test_qn_factory = .false.
        else if (associated(level, other_level) .or. &
                 .not. associated(level%extension%state, state) .or. &
                 .not. associated(other_level%extension%state, other_state)) then
            write(stderr, *) "test_qn_factory failed: Separate solver settings not "// &
                "wired with distinct objects reading their own state."
            test_qn_factory = .false.
        end if
        deallocate(solver_settings, other_settings)

    end function test_qn_factory

    logical(c_bool) function test_qn_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! quasi-Newton input parameters
        !
        use opentrustregion_unit_tests, only: test_context_type, setup_error_logging
        use otr_qn, only: qn_settings_type, qn_sanity_check, default_qn_settings

        type(qn_settings_type) :: settings
        type(test_context_type), target :: log_context
        integer(ip) :: case_n_param, i_case, error
        character(len=35), parameter :: case_names(3) = &
            [character(len=35) :: "a vanishing number of parameters", &
             "a negative maximum number of points", "an unknown updating scheme"]

        ! assume tests pass
        test_qn_sanity_check = .true.

        ! check valid parameters, with no stored points, which keep no history, and the
        ! updating scheme in mixed case, which has to be converted to lowercase
        settings = default_qn_settings
        call setup_error_logging(settings, log_context)
        settings%max_points = 0
        settings%hess_update_scheme = "BfGs"
        call qn_sanity_check(settings, 1_ip, error)
        if (error /= 0 .or. log_context%log_message /= "") then
            write(stderr, *) "test_qn_sanity_check failed: Valid parameters rejected."
            test_qn_sanity_check = .false.
        end if
        if (settings%hess_update_scheme /= "bfgs") then
            write(stderr, *) "test_qn_sanity_check failed: Updating scheme not "// &
                "converted to lowercase."
            test_qn_sanity_check = .false.
        end if

        ! check that every invalid parameter is rejected with an error message
        do i_case = 1, size(case_names)
            settings = default_qn_settings
            call setup_error_logging(settings, log_context)
            case_n_param = 1
            select case (i_case)
            case (1)
                case_n_param = 0
            case (2)
                settings%max_points = -1
            case (3)
                settings%hess_update_scheme = "dfp"
            end select
            call qn_sanity_check(settings, case_n_param, error)
            if (error /= 1) then
                write(stderr, *) "test_qn_sanity_check failed: Error not returned "// &
                    "for "//trim(case_names(i_case))//"."
                test_qn_sanity_check = .false.
            end if
            if (log_context%log_message == "") then
                write(stderr, *) "test_qn_sanity_check failed: No error message "// &
                    "printed for "//trim(case_names(i_case))//"."
                test_qn_sanity_check = .false.
            end if
        end do

    end function test_qn_sanity_check

    logical(c_bool) function test_update_orbs_qn() bind(C)
        !
        ! this function tests the quasi-Newton orbital update, which forwards to the
        ! host's orbital update with the host context of the solver level, returns its
        ! results and error, returns the quasi-Newton Hessian linear transformation in
        ! place of the host's and adds the point for a nonzero step and for the
        ! starting point, whose gradient the first step pair needs, passing on the error
        ! of the addition
        !
        use opentrustregion, only: hess_x_type
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: hess_x_callback
        use otr_common_mock, only: mock_update_orbs
        use otr_qn, only: default_qn_settings
        use test_reference, only: tol, arm_host_context, host_context_reached

        type(recording_qn_type), target :: object
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: kappa(n_param), func, grad(n_param), h_diag(n_param)
        procedure(hess_x_type), pointer :: hess_x_funptr
        integer(ip) :: error, i_case, i
        character(len=31), parameter :: case_names(5) = &
            [character(len=31) :: "the starting point", "a zero step after the start", &
             "a nonzero step", "a failing host's orbital update", "a failing addition"]
        logical, parameter :: case_adds(5) = [.true., .false., .true., .false., .true.]
        integer(ip), parameter :: case_errors(5) = [0, 0, 0, 3, 4]

        ! assume tests pass
        test_update_orbs_qn = .true.

        ! set up an object whose solver level forwards to the host's mock orbital update
        ! with a host context distinct from that of its settings
        object%n_param = n_param
        object%settings = default_qn_settings
        object%solver_level%host%update_orbs => mock_update_orbs
        object%solver_level%host%context => contexts(1)
        context => object%solver_level

        ! call the orbital update with a zero step before and after the starting point
        ! was added, with a nonzero step, whose entries sum to a negative number, so
        ! that only their magnitudes show it is nonzero, and with a failing host's
        ! orbital update or addition, and determine if it forwards with the host
        ! context of the solver level, returns the host's results and error, returns
        ! the role callback function of the Hessian linear transformation and adds the
        ! point exactly where it has to
        do i_case = 1, size(case_names)
            call arm_host_context(object%settings, contexts(2))
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
                write(stderr, *) "test_update_orbs_qn failed: Returned error code "// &
                    "wrong for "//trim(case_names(i_case))//"."
                test_update_orbs_qn = .false.
            end if
            if (abs(func - sum(kappa)) > tol .or. &
                any(abs(grad - 2.0_rp * kappa) > tol) .or. &
                any(abs(h_diag - 3.0_rp * kappa) > tol)) then
                write(stderr, *) "test_update_orbs_qn failed: Results of the "// &
                    "host's orbital update not returned for "// &
                    trim(case_names(i_case))//"."
                test_update_orbs_qn = .false.
            end if
            if (contexts(1)%n_calls /= 1 .or. contexts(2)%n_calls /= 0) then
                write(stderr, *) "test_update_orbs_qn failed: Host's orbital "// &
                    "update not called with the host context of the solver level "// &
                    "for "//trim(case_names(i_case))//"."
                test_update_orbs_qn = .false.
            end if
            if (.not. host_context_reached("update_orbs_qn", contexts(1))) &
                test_update_orbs_qn = .false.
            if (i_case /= 4 .and. .not. associated(hess_x_funptr, hess_x_callback)) then
                write(stderr, *) "test_update_orbs_qn failed: Role callback "// &
                    "function of the Hessian linear transformation not returned "// &
                    "for "//trim(case_names(i_case))//"."
                test_update_orbs_qn = .false.
            end if
            if (object%n_added /= merge(1, 0, case_adds(i_case))) then
                write(stderr, *) "test_update_orbs_qn failed: Point added wrongly "// &
                    "for "//trim(case_names(i_case))//"."
                test_update_orbs_qn = .false.
            else if (case_adds(i_case)) then
                if (any(abs(object%added_kappa - kappa) > tol) .or. &
                    any(abs(object%added_grad - grad) > tol)) then
                    write(stderr, *) "test_update_orbs_qn failed: Step or gradient "// &
                        "of the added point wrong for "//trim(case_names(i_case))//"."
                    test_update_orbs_qn = .false.
                end if
            end if
        end do

    end function test_update_orbs_qn

    logical(c_bool) function test_init_qn_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the quasi-Newton
        ! settings
        !
        use otr_qn, only: qn_settings_type, default_qn_settings
        use otr_qn_test_reference, only: operator(/=)

        type, extends(qn_settings_type) :: extended_settings_type
        end type

        type(qn_settings_type) :: settings
        type(extended_settings_type) :: extended_settings
        integer(ip) :: error

        ! assume tests pass
        test_init_qn_settings = .true.

        ! initialize settings, which have to be set to the default values
        call settings%init(error)
        if (error /= 0) then
            write(stderr, *) "test_init_qn_settings failed: Produced error."
            test_init_qn_settings = .false.
        end if
        if (settings /= default_qn_settings) then
            write(stderr, *) "test_init_qn_settings failed: Settings not "// &
                "initialized correctly."
            test_init_qn_settings = .false.
        end if

        ! initialize settings of a type extending these settings, which cannot be set to
        ! the default values, and check that an error is returned
        call extended_settings%init(error)
        if (error /= 1) then
            write(stderr, *) "test_init_qn_settings failed: No error returned for "// &
                "extended settings type."
            test_init_qn_settings = .false.
        end if

    end function test_init_qn_settings

    logical(c_bool) function test_init_qn() bind(C)
        !
        ! this function tests the subroutine which initializes a quasi-Newton object
        !
        use otr_qn, only: qn_settings_type, sr1_updating_type
        use otr_qn_test_reference, only: ref_qn_settings, operator(/=)

        type(sr1_updating_type) :: object
        type(qn_settings_type) :: settings

        ! assume tests pass
        test_init_qn = .true.

        ! initialize a new object, which has to perform the orbital update and the
        ! Hessian linear transformation itself, take over the settings and the number
        ! of parameters and allocate an empty history
        settings = ref_qn_settings
        call object%init(n_param, settings)
        if (.not. all([object%provides%update_orbs, object%provides%hess_x]) .or. &
            any([object%provides%obj_func, object%provides%precond, &
                 object%provides%precond_pd, object%provides%project, &
                 object%provides%get_extra_trial_vectors])) then
            write(stderr, *) "test_init_qn failed: Roles the object performs set wrong."
            test_init_qn = .false.
        end if
        if (object%n_param /= n_param .or. object%settings /= settings) then
            write(stderr, *) "test_init_qn failed: Settings or number of "// &
                "parameters not taken over."
            test_init_qn = .false.
        end if
        if (object%n_points /= 0 .or. allocated(object%last_grad) .or. &
            any(shape(object%kappa_list) /= [n_param, settings%max_points]) .or. &
            any(shape(object%grad_diff_list) /= [n_param, settings%max_points])) then
            write(stderr, *) "test_init_qn failed: History not empty or arrays "// &
                "allocated wrong."
            test_init_qn = .false.
        end if

    end function test_init_qn

    logical(c_bool) function test_add_qn() bind(C)
        !
        ! this function tests the subroutine which adds a new point to the history of a
        ! quasi-Newton object, so far only that it hands every call of the host's
        ! transport function the host context of its settings and passes on the error
        ! of every call and that a pair the skipping criterion rejects leaves the
        ! history unchanged while the gradient becomes the last one, while its numerics
        ! are left to later tests
        !
        use otr_qn, only: bfgs_updating_type
        use test_reference, only: host_context, host_context_reached

        type(bfgs_updating_type) :: object
        real(rp) :: kappa(n_param), grad(n_param), orthogonal(n_param)
        integer(ip) :: error, n_calls, i_call
        character(len=8) :: call_number

        ! assume tests pass
        test_add_qn = .true.

        ! add a point to an object with two points in its history and a last gradient,
        ! for which the routine transports both history arrays of each point and the
        ! last gradient, of the scheme whose skipping criterion calls no host's
        ! function
        call setup_updating_object(object)
        call random_number(kappa)
        call random_number(grad)
        call object%add(kappa, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_add_qn failed: Produced error."
            test_add_qn = .false.
        end if
        n_calls = host_context%n_calls
        if (n_calls /= 5) then
            write(stderr, *) "test_add_qn failed: Transport function not called "// &
                "with the host context at every call site."
            test_add_qn = .false.
        end if
        if (.not. host_context_reached("add_qn")) test_add_qn = .false.

        ! let every call of the transport function fail in turn, whose error the
        ! routine has to pass on
        do i_call = 1, n_calls
            call setup_updating_object(object)
            host_context%failing_call = i_call
            host_context%mock_error = 7
            call object%add(kappa, grad, error)
            if (error /= 7) then
                write(call_number, '(I0)') i_call
                write(stderr, *) "test_add_qn failed: Error of call "// &
                    trim(call_number)//" of the transport function not passed on."
                test_add_qn = .false.
            end if
        end do

        ! add a point whose gradient differs from the last one, which the mock transport
        ! leaves unchanged, by a random vector orthogonal to the step, which the
        ! skipping criterion of the scheme rejects, so that the history has to keep its
        ! points and the gradient has to become the last one, from which the next pair
        ! starts
        call setup_updating_object(object)
        call random_number(orthogonal)
        orthogonal = orthogonal - &
                     dot_product(orthogonal, kappa) / dot_product(kappa, kappa) * kappa
        grad = object%last_grad + orthogonal
        call object%add(kappa, grad, error)
        if (error /= 0 .or. object%n_points /= 2) then
            write(stderr, *) "test_add_qn failed: Error produced or history "// &
                "changed for a pair the skipping criterion rejects."
            test_add_qn = .false.
        end if
        if (norm2(object%last_grad - grad) > tol) then
            write(stderr, *) "test_add_qn failed: Gradient not the last one after "// &
                "a pair the skipping criterion rejects."
            test_add_qn = .false.
        end if

    end function test_add_qn

    logical(c_bool) function test_sr1_skip_add() bind(C)
        !
        ! this function tests the function which decides whether the SR1 updating
        ! scheme skips a new point, so far only that it hands the host's initial
        ! Hessian function the host context of its settings and passes on the error of
        ! every call, while its numerics are left to later tests
        !
        use otr_qn, only: sr1_updating_type
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type) :: object
        real(rp) :: kappa(n_param), grad_diff(n_param)
        logical :: skip
        integer(ip) :: error, n_calls, i_call
        character(len=8) :: call_number

        ! assume tests pass
        test_sr1_skip_add = .true.

        ! decide on a new point
        call setup_updating_object(object)
        call random_number(kappa)
        grad_diff = 3.0_rp * kappa
        skip = object%skip_add(kappa, grad_diff, error)
        if (error /= 0) then
            write(stderr, *) "test_sr1_skip_add failed: Produced error."
            test_sr1_skip_add = .false.
        end if
        n_calls = host_context%n_calls
        if (n_calls /= 1) then
            write(stderr, *) "test_sr1_skip_add failed: Initial Hessian function "// &
                "not called with the host context."
            test_sr1_skip_add = .false.
        end if
        if (.not. host_context_reached("sr1_skip_add")) test_sr1_skip_add = .false.

        ! let every call of the initial Hessian function fail in turn, whose error the
        ! routine has to pass on
        do i_call = 1, n_calls
            call setup_updating_object(object)
            host_context%failing_call = i_call
            host_context%mock_error = 7
            skip = object%skip_add(kappa, grad_diff, error)
            if (error /= 7) then
                write(call_number, '(I0)') i_call
                write(stderr, *) "test_sr1_skip_add failed: Error of call "// &
                    trim(call_number)//" of the initial Hessian function not passed on."
                test_sr1_skip_add = .false.
            end if
        end do

    end function test_sr1_skip_add

    logical(c_bool) function test_hess_x_sr1() bind(C)
        !
        ! this function tests the Hessian linear transformation of the SR1 updating
        ! scheme, so far only that it hands every call of the host's initial Hessian
        ! function the host context of its settings and passes on the error of every
        ! call, while its numerics are left to later tests
        !
        use otr_common_unit_tests, only: level_context
        use otr_qn, only: sr1_updating_type
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type), target :: object
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error, n_calls, i_call
        character(len=8) :: call_number

        ! assume tests pass
        test_hess_x_sr1 = .true.

        ! apply the transformation of an object with two points in its history, for
        ! which the initial Hessian is applied to the trial vector and to the step of
        ! each point
        call setup_updating_object(object)
        call random_number(x)
        call object%hess_x(x, hess_x, error, level_context(object, 1_ip))
        if (error /= 0) then
            write(stderr, *) "test_hess_x_sr1 failed: Produced error."
            test_hess_x_sr1 = .false.
        end if
        n_calls = host_context%n_calls
        if (n_calls /= 3) then
            write(stderr, *) "test_hess_x_sr1 failed: Initial Hessian function not "// &
                "called with the host context at every call site."
            test_hess_x_sr1 = .false.
        end if
        if (.not. host_context_reached("hess_x_sr1")) test_hess_x_sr1 = .false.

        ! let every call of the initial Hessian function fail in turn, whose error the
        ! routine has to pass on
        do i_call = 1, n_calls
            call setup_updating_object(object)
            host_context%failing_call = i_call
            host_context%mock_error = 7
            call object%hess_x(x, hess_x, error, level_context(object, 1_ip))
            if (error /= 7) then
                write(call_number, '(I0)') i_call
                write(stderr, *) "test_hess_x_sr1 failed: Error of call "// &
                    trim(call_number)//" of the initial Hessian function not passed on."
                test_hess_x_sr1 = .false.
            end if
        end do

    end function test_hess_x_sr1

    logical(c_bool) function test_hess_x_bfgs() bind(C)
        !
        ! this function tests the Hessian linear transformation of the BFGS updating
        ! scheme, so far only that it hands every call of the host's initial Hessian
        ! function the host context of its settings and passes on the error of every
        ! call, while its numerics are left to later tests
        !
        use otr_common_unit_tests, only: level_context
        use otr_qn, only: bfgs_updating_type
        use test_reference, only: host_context, host_context_reached

        type(bfgs_updating_type), target :: object
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error, n_calls, i_call
        character(len=8) :: call_number

        ! assume tests pass
        test_hess_x_bfgs = .true.

        ! apply the transformation of an object with two points in its history, for
        ! which the initial Hessian is applied to the trial vector, to the step of each
        ! point when the compact scaling matrix is built and to each again when the
        ! transformation is assembled
        call setup_updating_object(object)
        call random_number(x)
        call object%hess_x(x, hess_x, error, level_context(object, 1_ip))
        if (error /= 0) then
            write(stderr, *) "test_hess_x_bfgs failed: Produced error."
            test_hess_x_bfgs = .false.
        end if
        n_calls = host_context%n_calls
        if (n_calls /= 5) then
            write(stderr, *) "test_hess_x_bfgs failed: Initial Hessian function "// &
                "not called with the host context at every call site."
            test_hess_x_bfgs = .false.
        end if
        if (.not. host_context_reached("hess_x_bfgs")) test_hess_x_bfgs = .false.

        ! let every call of the initial Hessian function fail in turn, whose error the
        ! routine has to pass on
        do i_call = 1, n_calls
            call setup_updating_object(object)
            host_context%failing_call = i_call
            host_context%mock_error = 7
            call object%hess_x(x, hess_x, error, level_context(object, 1_ip))
            if (error /= 7) then
                write(call_number, '(I0)') i_call
                write(stderr, *) "test_hess_x_bfgs failed: Error of call "// &
                    trim(call_number)//" of the initial Hessian function not passed on."
                test_hess_x_bfgs = .false.
            end if
        end do

    end function test_hess_x_bfgs

end module otr_qn_unit_tests
