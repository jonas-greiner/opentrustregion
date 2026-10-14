! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: n_param
    use otr_qn, only: transport_type, init_hess_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(transport_type), pointer :: mock_transport_ptr => mock_transport
    procedure(init_hess_type), pointer :: mock_init_hess_ptr => mock_init_hess

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

    subroutine setup_updating_object(object)
        !
        ! this subroutine sets up an updating object as after two orbital updates,
        ! without the routines under test, with the default settings carrying the armed
        ! host context, the mock transport and initial Hessian functions, one point in
        ! its history, whose gradient difference has positive curvature along its step,
        ! and the last gradient
        !
        use otr_qn, only: updating_type, default_qn_settings
        use test_reference, only: arm_host_context

        class(updating_type), intent(inout) :: object

        integer(ip) :: i

        object%n_param = n_param
        object%settings = default_qn_settings
        call arm_host_context(object%settings)
        object%transport => mock_transport
        object%init_hess => mock_init_hess
        allocate(object%kappa_list(n_param, object%settings%max_points), &
                 object%grad_diff_list(n_param, object%settings%max_points))
        object%kappa_list = 0.0_rp
        object%grad_diff_list = 0.0_rp
        object%n_points = 1
        object%kappa_list(:, 1) = [(0.1_rp * i, i=1, n_param)]
        object%grad_diff_list(:, 1) = 2.0_rp * object%kappa_list(:, 1)
        object%last_grad = [(0.2_rp * i - 0.3_rp, i=1, n_param)]

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
        use otr_qn_test_reference, only: ref_qn_settings, assignment(=), operator(/=)

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
                deallocate(solver_settings)
                cycle
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
        ! its error and leave the solver settings untouched
        allocate(solver_settings)
        solver_settings%solver_settings_type = default_solver_settings
        call setup_settings(settings, context)
        call qn_factory(0_ip, update_orbs_funptr, obj_func_funptr, transport_funptr, &
                        init_hess_funptr, update_orbs_qn_funptr, obj_func_qn_funptr, &
                        solver_settings, state, error, settings)
        if (error == 0 .or. associated(solver_settings%context) .or. &
            allocated(solver_settings%extensions)) then
            write(stderr, *) "test_qn_factory failed: Error not returned or solver "// &
                "settings changed for parameters the sanity check rejects."
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
        call qn_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                        transport_funptr, init_hess_funptr, update_orbs_qn_funptr, &
                        obj_func_qn_funptr, other_settings, other_state, error, &
                        settings)
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
            if (error == 0) then
                write(stderr, *) "test_qn_sanity_check failed: No error returned "// &
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
        use otr_qn_test_reference, only: ref_qn_settings, assignment(=), operator(/=)

        type(sr1_updating_type) :: object
        type(qn_settings_type) :: settings

        ! assume tests pass
        test_init_qn = .true.

        ! leave a history, a last gradient and arrays of another dimension from an
        ! earlier setup behind and initialize the object, which has to discard these,
        ! perform the orbital update and the Hessian linear transformation itself and
        ! take over the settings and the number of parameters
        object%n_points = 2
        allocate(object%kappa_list(n_param + 1, 2), &
                 object%grad_diff_list(n_param + 1, 2), object%last_grad(n_param + 1))
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
            write(stderr, *) "test_init_qn failed: History of an earlier setup "// &
                "kept or arrays allocated wrong."
            test_init_qn = .false.
        end if

    end function test_init_qn

    logical(c_bool) function test_add_qn() bind(C)
        !
        ! this function tests the subroutine which adds a new point to the history of a
        ! quasi-Newton object, so far only that it hands every call of the host's
        ! transport function the host context of its settings, while its numerics are
        ! left to later tests
        !
        use otr_qn, only: bfgs_updating_type
        use test_reference, only: host_context, host_context_reached

        type(bfgs_updating_type) :: object
        real(rp) :: kappa(n_param), grad(n_param)
        integer(ip) :: i, error

        ! assume tests pass
        test_add_qn = .true.

        ! add a point to an object with a point in its history and a last gradient, for
        ! which the routine transports both history arrays and the last gradient, of
        ! the scheme whose skipping criterion calls no host's function
        call setup_updating_object(object)
        kappa = [(0.05_rp * i, i=1, n_param)]
        grad = [(0.4_rp * i, i=1, n_param)]
        call object%add(kappa, grad, error)
        if (error /= 0) then
            write(stderr, *) "test_add_qn failed: Produced error."
            test_add_qn = .false.
        end if
        if (host_context%n_calls /= 3) then
            write(stderr, *) "test_add_qn failed: Transport function not called "// &
                "with the host context at every call site."
            test_add_qn = .false.
        end if
        if (.not. host_context_reached("add_qn")) test_add_qn = .false.

    end function test_add_qn

    logical(c_bool) function test_sr1_skip_add() bind(C)
        !
        ! this function tests the function which decides whether the SR1 updating
        ! scheme skips a new point, so far only that it hands the host's initial
        ! Hessian function the host context of its settings, while its numerics are
        ! left to later tests
        !
        use otr_qn, only: sr1_updating_type
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type) :: object
        real(rp) :: kappa(n_param), grad_diff(n_param)
        logical :: skip
        integer(ip) :: i, error

        ! assume tests pass
        test_sr1_skip_add = .true.

        ! decide on a new point
        call setup_updating_object(object)
        kappa = [(0.05_rp * i, i=1, n_param)]
        grad_diff = 3.0_rp * kappa
        skip = object%skip_add(kappa, grad_diff, error)
        if (error /= 0) then
            write(stderr, *) "test_sr1_skip_add failed: Produced error."
            test_sr1_skip_add = .false.
        end if
        if (host_context%n_calls /= 1) then
            write(stderr, *) "test_sr1_skip_add failed: Initial Hessian function "// &
                "not called with the host context."
            test_sr1_skip_add = .false.
        end if
        if (.not. host_context_reached("sr1_skip_add")) test_sr1_skip_add = .false.

    end function test_sr1_skip_add

    logical(c_bool) function test_hess_x_sr1() bind(C)
        !
        ! this function tests the Hessian linear transformation of the SR1 updating
        ! scheme, so far only that it hands every call of the host's initial Hessian
        ! function the host context of its settings, while its numerics are left to
        ! later tests
        !
        use otr_common_unit_tests, only: level_context
        use otr_qn, only: sr1_updating_type
        use test_reference, only: host_context, host_context_reached

        type(sr1_updating_type), target :: object
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_hess_x_sr1 = .true.

        ! apply the transformation of an object with a point in its history, for which
        ! the initial Hessian is applied to the trial vector and to the step of the
        ! point
        call setup_updating_object(object)
        x = 1.0_rp
        call object%hess_x(x, hess_x, error, level_context(object, 1_ip))
        if (error /= 0) then
            write(stderr, *) "test_hess_x_sr1 failed: Produced error."
            test_hess_x_sr1 = .false.
        end if
        if (host_context%n_calls /= 2) then
            write(stderr, *) "test_hess_x_sr1 failed: Initial Hessian function not "// &
                "called with the host context at every call site."
            test_hess_x_sr1 = .false.
        end if
        if (.not. host_context_reached("hess_x_sr1")) test_hess_x_sr1 = .false.

    end function test_hess_x_sr1

    logical(c_bool) function test_hess_x_bfgs() bind(C)
        !
        ! this function tests the Hessian linear transformation of the BFGS updating
        ! scheme, so far only that it hands every call of the host's initial Hessian
        ! function the host context of its settings, while its numerics are left to
        ! later tests
        !
        use otr_common_unit_tests, only: level_context
        use otr_qn, only: bfgs_updating_type
        use test_reference, only: host_context, host_context_reached

        type(bfgs_updating_type), target :: object
        real(rp) :: x(n_param), hess_x(n_param)
        integer(ip) :: error

        ! assume tests pass
        test_hess_x_bfgs = .true.

        ! apply the transformation of an object with a point in its history, for which
        ! the initial Hessian is applied to the trial vector, to the step of the point
        ! when the compact scaling matrix is built and to it again when the
        ! transformation is assembled
        call setup_updating_object(object)
        x = 1.0_rp
        call object%hess_x(x, hess_x, error, level_context(object, 1_ip))
        if (error /= 0) then
            write(stderr, *) "test_hess_x_bfgs failed: Produced error."
            test_hess_x_bfgs = .false.
        end if
        if (host_context%n_calls /= 3) then
            write(stderr, *) "test_hess_x_bfgs failed: Initial Hessian function "// &
                "not called with the host context at every call site."
            test_hess_x_bfgs = .false.
        end if
        if (.not. host_context_reached("hess_x_bfgs")) test_hess_x_bfgs = .false.

    end function test_hess_x_bfgs

end module otr_qn_unit_tests
