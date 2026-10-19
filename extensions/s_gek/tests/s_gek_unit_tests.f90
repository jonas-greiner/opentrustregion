! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol, n_param
    use otr_s_gek, only: change_reference_type, s_gek_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(change_reference_type), pointer :: mock_change_reference_ptr => &
        mock_change_reference

    ! S-GEK object which records the point its orbital update adds and returns the
    ! error code it is set to, for the test of the orbital update
    type, extends(s_gek_type) :: recording_s_gek_type
        integer(ip) :: n_added = 0, error = 0
        real(rp) :: added_func = 0.0_rp
        real(rp), allocatable :: added_kappa(:), added_grad(:), added_h_diag(:)
    contains
        procedure :: add => recording_add_s_gek
    end type

contains

    subroutine mock_change_reference(new_ref, n_points, kappa_list, local_grad_list, &
                                     grad_list, error, context)
        !
        ! this subroutine is a mock of the host's change of reference function, which
        ! leaves the parameters unchanged and takes the gradients over from the local
        ! gradients, as a change to a reference at which both coincide does, records
        ! its call, checks the host context and returns the error code set in it
        !
        use otr_common_mock, only: record_host_call

        real(rp), intent(in), target :: new_ref(:)
        integer(ip), intent(in) :: n_points
        real(rp), intent(inout), target :: kappa_list(:, :), local_grad_list(:, :), &
                                           grad_list(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        grad_list(:, :n_points) = local_grad_list(:, :n_points)
        error = record_host_call("change_reference", context)

    end subroutine mock_change_reference

    subroutine recording_add_s_gek(self, kappa, obj_func, grad, h_diag, error)
        !
        ! this subroutine records the step, objective function, gradient and Hessian
        ! diagonal of the point the recording S-GEK object is given and returns the
        ! error code it is set to
        !
        class(recording_s_gek_type), intent(inout) :: self
        real(rp), intent(in) :: kappa(:), obj_func, grad(:), h_diag(:)
        integer(ip), intent(out) :: error

        self%n_added = self%n_added + 1
        self%added_kappa = kappa
        self%added_func = obj_func
        self%added_grad = grad
        self%added_h_diag = h_diag
        error = self%error

    end subroutine recording_add_s_gek

    logical(c_bool) function test_s_gek_factory() bind(C)
        !
        ! this function tests the factory, which creates an S-GEK object after the
        ! sanity check of the settings, initializes it with the settings and the host's
        ! change of reference function and wires the solver settings with it and the
        ! record of the state
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, solver_state_type
        use opentrustregion_unit_tests, only: setup_settings, test_context_type, logger
        use otr_common_test_reference, only: same_logger_and_context
        use otr_common, only: extension_solver_settings_type, extension_level_type, &
                              update_orbs_callback, obj_func_callback
        use otr_common_unit_tests, only: wired_level
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_s_gek, only: s_gek_factory, s_gek_settings_type, s_gek_type, &
                             default_s_gek_settings, error_s_gek
        use otr_s_gek_test_reference, only: ref_s_gek_settings, operator(/=)

        procedure(update_orbs_type), pointer :: update_orbs_funptr, &
                                                update_orbs_s_gek_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr, obj_func_s_gek_funptr
        procedure(change_reference_type), pointer :: change_reference_funptr
        type(extension_solver_settings_type), allocatable, target :: solver_settings, &
                                                                     other_settings
        type(solver_state_type), target :: state, other_state
        type(s_gek_settings_type) :: settings, expected_settings
        type(test_context_type), target :: context
        type(extension_level_type), pointer :: level, other_level
        integer(ip) :: i_case, error
        logical :: owned
        character(len=16), parameter :: case_names(2) = &
            [character(len=16) :: "for settings", "without settings"]

        ! assume tests pass
        test_s_gek_factory = .true.

        ! set the host's callback functions
        update_orbs_funptr => mock_update_orbs
        obj_func_funptr => mock_obj_func
        change_reference_funptr => mock_change_reference

        ! create objects for settings and for no settings, each of which has to be
        ! initialized with the settings and the
        ! host's change of reference function and to wire the solver settings, which
        ! own it, with the host's orbital update and objective function and the record
        ! of the state
        do i_case = 1, size(case_names)
            allocate(solver_settings)
            select case (i_case)
            case (1)
                settings = ref_s_gek_settings
                settings%logger => logger
                settings%context => context
                expected_settings = settings
            case (2)
                expected_settings = default_s_gek_settings
            end select
            if (i_case == 2) then
                call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                                   change_reference_funptr, update_orbs_s_gek_funptr, &
                                   obj_func_s_gek_funptr, solver_settings, state, error)
            else
                call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                                   change_reference_funptr, update_orbs_s_gek_funptr, &
                                   obj_func_s_gek_funptr, solver_settings, state, &
                                   error, settings)
            end if
            if (error /= 0) then
                write(stderr, *) "test_s_gek_factory failed: Produced error "// &
                    trim(case_names(i_case))//"."
                test_s_gek_factory = .false.
                deallocate(solver_settings)
                cycle
            end if
            if (.not. associated(update_orbs_s_gek_funptr, update_orbs_callback) .or. &
                .not. associated(obj_func_s_gek_funptr, obj_func_callback)) then
                write(stderr, *) "test_s_gek_factory failed: Role callback "// &
                    "functions not returned "//trim(case_names(i_case))//"."
                test_s_gek_factory = .false.
            end if
            level => wired_level(solver_settings)
            if (.not. associated(level)) then
                write(stderr, *) "test_s_gek_factory failed: Solver settings not "// &
                    "wired with an object "//trim(case_names(i_case))//"."
                test_s_gek_factory = .false.
                deallocate(solver_settings)
                cycle
            end if
            owned = .false.
            if (allocated(solver_settings%extensions)) then
                if (size(solver_settings%extensions) == 1) owned = &
                    associated(level%extension, solver_settings%extensions(1)%object)
            end if
            if (.not. owned) then
                write(stderr, *) "test_s_gek_factory failed: Object not owned by "// &
                    "the solver settings "//trim(case_names(i_case))//"."
                test_s_gek_factory = .false.
            end if
            select type (object => level%extension)
            type is (s_gek_type)
                if (object%n_param /= n_param .or. &
                    object%settings /= expected_settings .or. &
                    .not. same_logger_and_context(object%settings, expected_settings)) &
                    then
                    write(stderr, *) "test_s_gek_factory failed: Object not "// &
                        "initialized with the settings "//trim(case_names(i_case))//"."
                    test_s_gek_factory = .false.
                end if
                if (.not. associated(object%change_reference, mock_change_reference)) &
                    then
                    write(stderr, *) "test_s_gek_factory failed: Host's change of "// &
                        "reference function not handed to the object "// &
                        trim(case_names(i_case))//"."
                    test_s_gek_factory = .false.
                end if
            class default
                write(stderr, *) "test_s_gek_factory failed: Object of the wrong "// &
                    "type created "//trim(case_names(i_case))//"."
                test_s_gek_factory = .false.
            end select
            if (.not. associated(level%host%update_orbs, mock_update_orbs) .or. &
                .not. associated(level%host%obj_func, mock_obj_func) .or. &
                .not. associated(level%extension%state, state)) then
                write(stderr, *) "test_s_gek_factory failed: Solver settings not "// &
                    "wired with the host's functions and the record of the state "// &
                    trim(case_names(i_case))//"."
                test_s_gek_factory = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create an object for parameters the sanity check rejects, which has to return
        ! its error without functions and leave the solver settings untouched
        allocate(solver_settings)
        call setup_settings(settings, context)
        update_orbs_s_gek_funptr => update_orbs_callback
        obj_func_s_gek_funptr => obj_func_callback
        call s_gek_factory(0_ip, update_orbs_funptr, obj_func_funptr, &
                           change_reference_funptr, update_orbs_s_gek_funptr, &
                           obj_func_s_gek_funptr, solver_settings, state, error, &
                           settings)
        if (error /= error_s_gek + 1) then
            write(stderr, *) "test_s_gek_factory failed: Error not returned for "// &
                "parameters the sanity check rejects."
            test_s_gek_factory = .false.
        end if
        if (associated(update_orbs_s_gek_funptr) .or. &
            associated(obj_func_s_gek_funptr)) then
            write(stderr, *) "test_s_gek_factory failed: Functions returned for "// &
                "parameters the sanity check rejects."
            test_s_gek_factory = .false.
        end if
        if (associated(solver_settings%context) .or. &
            allocated(solver_settings%extensions)) then
            write(stderr, *) "test_s_gek_factory failed: Solver settings changed "// &
                "for parameters the sanity check rejects."
            test_s_gek_factory = .false.
        end if
        deallocate(solver_settings)

        ! create an object for solver settings a factory wired, given the host's own
        ! orbital update and objective function instead of those the factory returned,
        ! which the wiring refuses, so that the factory has to return its error without
        ! functions and leave the settings wired with the first object only
        allocate(solver_settings)
        call setup_settings(settings, context)
        call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                           change_reference_funptr, update_orbs_s_gek_funptr, &
                           obj_func_s_gek_funptr, solver_settings, state, error, &
                           settings)
        level => wired_level(solver_settings)
        call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                           change_reference_funptr, update_orbs_s_gek_funptr, &
                           obj_func_s_gek_funptr, solver_settings, state, error, &
                           settings)
        if (error /= error_s_gek + 1) then
            write(stderr, *) "test_s_gek_factory failed: Error not returned for "// &
                "the host's own functions given with solver settings a factory wired."
            test_s_gek_factory = .false.
        end if
        if (associated(update_orbs_s_gek_funptr) .or. &
            associated(obj_func_s_gek_funptr)) then
            write(stderr, *) "test_s_gek_factory failed: Functions returned for "// &
                "the host's own functions given with solver settings a factory wired."
            test_s_gek_factory = .false.
        end if
        owned = .false.
        if (associated(level) .and. allocated(solver_settings%extensions)) &
            owned = size(solver_settings%extensions) == 1 .and. &
                    associated(wired_level(solver_settings), level)
        if (.not. owned) then
            write(stderr, *) "test_s_gek_factory failed: Solver settings a factory "// &
                "wired changed for the host's own functions."
            test_s_gek_factory = .false.
        end if
        deallocate(solver_settings)

        ! create objects for two separate solver settings, which have to be distinct
        ! objects reading their own record of the state, so that calculations with
        ! separate settings can be used alternately
        allocate(solver_settings, other_settings)
        settings = ref_s_gek_settings
        call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                           change_reference_funptr, update_orbs_s_gek_funptr, &
                           obj_func_s_gek_funptr, solver_settings, state, error, &
                           settings)
        if (error /= 0) then
            write(stderr, *) "test_s_gek_factory failed: Produced error for the "// &
                "first of two separate solver settings."
            test_s_gek_factory = .false.
        end if
        call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                           change_reference_funptr, update_orbs_s_gek_funptr, &
                           obj_func_s_gek_funptr, other_settings, other_state, error, &
                           settings)
        if (error /= 0) then
            write(stderr, *) "test_s_gek_factory failed: Produced error for the "// &
                "second of two separate solver settings."
            test_s_gek_factory = .false.
        end if
        level => wired_level(solver_settings)
        other_level => wired_level(other_settings)
        if (.not. (associated(level) .and. associated(other_level))) then
            write(stderr, *) "test_s_gek_factory failed: Separate solver settings "// &
                "not wired."
            test_s_gek_factory = .false.
        else if (associated(level, other_level) .or. &
                 .not. associated(level%extension%state, state) .or. &
                 .not. associated(other_level%extension%state, other_state)) then
            write(stderr, *) "test_s_gek_factory failed: Separate solver settings "// &
                "not wired with distinct objects reading their own state."
            test_s_gek_factory = .false.
        end if
        deallocate(solver_settings, other_settings)

    end function test_s_gek_factory

    logical(c_bool) function test_s_gek_sanity_check() bind(C)
        !
        ! this function tests the subroutine which performs a sanity check for the
        ! S-GEK input parameters
        !
        use opentrustregion_unit_tests, only: test_context_type, setup_error_logging
        use otr_s_gek, only: s_gek_settings_type, s_gek_sanity_check, &
                             default_s_gek_settings

        type(s_gek_settings_type) :: settings
        type(test_context_type), target :: log_context
        integer(ip) :: case_n_param, i_case, error
        character(len=36), parameter :: case_names(2) = &
            [character(len=36) :: "a vanishing number of parameters", &
             "a vanishing maximum number of points"]

        ! assume tests pass
        test_s_gek_sanity_check = .true.

        ! check valid parameters, with a single stored point, the current one
        settings = default_s_gek_settings
        call setup_error_logging(settings, log_context)
        settings%max_points = 1
        call s_gek_sanity_check(settings, 1_ip, error)
        if (error /= 0 .or. log_context%log_message /= "") then
            write(stderr, *) "test_s_gek_sanity_check failed: Valid parameters "// &
                "rejected."
            test_s_gek_sanity_check = .false.
        end if

        ! check that every invalid parameter is rejected with an error message
        do i_case = 1, size(case_names)
            settings = default_s_gek_settings
            call setup_error_logging(settings, log_context)
            case_n_param = 1
            select case (i_case)
            case (1)
                case_n_param = 0
            case (2)
                settings%max_points = 0
            end select
            call s_gek_sanity_check(settings, case_n_param, error)
            if (error /= 1) then
                write(stderr, *) "test_s_gek_sanity_check failed: Error not "// &
                    "returned for "//trim(case_names(i_case))//"."
                test_s_gek_sanity_check = .false.
            end if
            if (log_context%log_message == "") then
                write(stderr, *) "test_s_gek_sanity_check failed: No error message "// &
                    "printed for "//trim(case_names(i_case))//"."
                test_s_gek_sanity_check = .false.
            end if
        end do

    end function test_s_gek_sanity_check

    logical(c_bool) function test_update_orbs_s_gek() bind(C)
        !
        ! this function tests the S-GEK orbital update, which forwards to the host's
        ! orbital update with the host context of the solver level, returns its results
        ! and error, returns the S-GEK Hessian linear transformation in place of the
        ! host's and adds the point for a nonzero step and for the starting point, so
        ! that the S-GEK Hessian linear transformation never runs on a model without
        ! points, passing on the error of the addition
        !
        use opentrustregion, only: hess_x_type
        use opentrustregion_unit_tests, only: test_context_type
        use otr_common, only: hess_x_callback
        use otr_common_mock, only: mock_update_orbs
        use otr_s_gek, only: default_s_gek_settings, error_s_gek
        use test_reference, only: arm_host_context, host_context_reached

        type(recording_s_gek_type), target :: object
        type(test_context_type), target :: contexts(2)
        class(*), pointer :: context
        real(rp) :: kappa(n_param), func, grad(n_param), h_diag(n_param)
        procedure(hess_x_type), pointer :: hess_x_funptr
        integer(ip) :: error, i_case
        character(len=31), parameter :: case_names(5) = &
            [character(len=31) :: "the starting point", "a zero step after the start", &
             "a nonzero step", "a failing host's orbital update", "a failing addition"]
        logical, parameter :: case_adds(5) = [.true., .false., .true., .false., .true.]
        integer(ip), parameter :: case_errors(5) = [0, 0, 0, 3, 4]

        ! assume tests pass
        test_update_orbs_s_gek = .true.

        ! set up an object whose solver level forwards to the host's mock orbital update
        ! with a host context distinct from that of its settings
        object%n_param = n_param
        object%settings = default_s_gek_settings
        object%error_origin = error_s_gek
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
            object%n_points = merge(0, 1, i_case == 1)
            kappa = 0.0_rp
            if (i_case > 2) then
                call random_number(kappa)
                kappa = kappa - 1.0_rp
            end if
            call object%update_orbs(kappa, func, grad, h_diag, hess_x_funptr, error, &
                                    context)
            if (error /= case_errors(i_case)) then
                write(stderr, *) "test_update_orbs_s_gek failed: Returned error "// &
                    "code wrong for "//trim(case_names(i_case))//"."
                test_update_orbs_s_gek = .false.
            end if
            if (abs(func - sum(kappa)) > tol .or. &
                any(abs(grad - 2.0_rp * kappa) > tol) .or. &
                any(abs(h_diag - 3.0_rp * kappa) > tol)) then
                write(stderr, *) "test_update_orbs_s_gek failed: Results of the "// &
                    "host's orbital update not returned for "// &
                    trim(case_names(i_case))//"."
                test_update_orbs_s_gek = .false.
            end if
            if (contexts(1)%n_calls /= 1 .or. contexts(2)%n_calls /= 0) then
                write(stderr, *) "test_update_orbs_s_gek failed: Host's orbital "// &
                    "update not called with the host context of the solver level "// &
                    "for "//trim(case_names(i_case))//"."
                test_update_orbs_s_gek = .false.
            end if
            if (.not. host_context_reached("update_orbs_s_gek", contexts(1))) &
                test_update_orbs_s_gek = .false.
            if (i_case /= 4 .and. .not. associated(hess_x_funptr, hess_x_callback)) then
                write(stderr, *) "test_update_orbs_s_gek failed: Role callback "// &
                    "function of the Hessian linear transformation not returned "// &
                    "for "//trim(case_names(i_case))//"."
                test_update_orbs_s_gek = .false.
            end if
            if (object%n_added /= merge(1, 0, case_adds(i_case))) then
                write(stderr, *) "test_update_orbs_s_gek failed: Point added "// &
                    "wrongly for "//trim(case_names(i_case))//"."
                test_update_orbs_s_gek = .false.
            else if (case_adds(i_case)) then
                if (any(abs(object%added_kappa - kappa) > tol) .or. &
                    abs(object%added_func - func) > tol .or. &
                    any(abs(object%added_grad - grad) > tol) .or. &
                    any(abs(object%added_h_diag - h_diag) > tol)) then
                    write(stderr, *) "test_update_orbs_s_gek failed: Values of the "// &
                        "added point wrong for "//trim(case_names(i_case))//"."
                    test_update_orbs_s_gek = .false.
                end if
            end if
        end do

    end function test_update_orbs_s_gek

    logical(c_bool) function test_init_s_gek() bind(C)
        !
        ! this function tests the subroutine which initializes an S-GEK object
        !
        use otr_s_gek, only: s_gek_settings_type, s_gek_type, error_s_gek
        use otr_s_gek_test_reference, only: ref_s_gek_settings, operator(/=)
        use otr_common_test_reference, only: same_logger_and_context
        use opentrustregion_unit_tests, only: test_context_type, logger

        type(s_gek_type) :: object
        type(s_gek_settings_type) :: settings
        type(test_context_type), target :: context

        ! assume tests pass
        test_init_s_gek = .true.

        ! give the counters and the loosening factor, which a new object leaves
        ! undefined, other values and initialize the object, which has to set them,
        ! perform the orbital update and the Hessian linear transformation itself, take
        ! over the settings and the number of parameters and allocate empty lists
        object%n_points = 2
        object%n_covariance = 5
        object%loosen_factor = 2.0_rp
        settings = ref_s_gek_settings
        settings%logger => logger
        settings%context => context
        call object%init(n_param, settings)
        if (.not. all([object%provides%update_orbs, object%provides%hess_x]) .or. &
            any([object%provides%obj_func, object%provides%precond, &
                 object%provides%precond_pd, object%provides%project, &
                 object%provides%get_extra_trial_vectors])) then
            write(stderr, *) "test_init_s_gek failed: Roles the object performs "// &
                "set wrong."
            test_init_s_gek = .false.
        end if
        if (object%error_origin /= error_s_gek) then
            write(stderr, *) "test_init_s_gek failed: Error origin of the object "// &
                "not set."
            test_init_s_gek = .false.
        end if
        if (object%n_param /= n_param .or. object%settings /= settings .or. &
            .not. same_logger_and_context(object%settings, settings)) then
            write(stderr, *) "test_init_s_gek failed: Settings or number of "// &
                "parameters not taken over."
            test_init_s_gek = .false.
        end if
        if (object%n_points /= 0 .or. object%n_covariance /= 0 .or. &
            abs(object%loosen_factor - 1.0_rp) > tol) then
            write(stderr, *) "test_init_s_gek failed: Points, covariance size or "// &
                "loosening factor not initialized."
            test_init_s_gek = .false.
        end if
        if (any(shape(object%h_diag) /= [n_param]) .or. &
            any(shape(object%kappa_list) /= [n_param, 0_ip]) .or. &
            any(shape(object%disp_list) /= [n_param, 0_ip]) .or. &
            any(shape(object%obj_func_list) /= [0]) .or. &
            any(shape(object%grad_list) /= [n_param, 0_ip]) .or. &
            any(shape(object%local_grad_list) /= [n_param, 0_ip]) .or. &
            any(shape(object%y_list) /= [n_param, 0_ip])) then
            write(stderr, *) "test_init_s_gek failed: Arrays allocated wrong."
            test_init_s_gek = .false.
        end if

    end function test_init_s_gek

    logical(c_bool) function test_add_s_gek() bind(C)
        !
        ! this function tests the subroutine which adds a new point to an S-GEK object,
        ! so far only that it hands the host's change of reference function the host
        ! context of its settings and passes on its error, while its numerics are left
        ! to later tests
        !
        use otr_s_gek, only: s_gek_type, default_s_gek_settings, error_s_gek
        use test_reference, only: arm_host_context, host_context, host_context_reached

        type(s_gek_type) :: objects(2)
        real(rp) :: kappa(n_param), grad(n_param), h_diag(n_param)
        integer(ip) :: error, i_case

        ! assume tests pass
        test_add_s_gek = .true.

        ! add a random first point with a positive Hessian diagonal to objects as after
        ! their initialization, without the routine that initializes them, for which
        ! the routine moves the stored points to the new reference, with a change of
        ! reference function that succeeds and with one that fails, whose error the
        ! routine has to pass on
        call random_number(kappa)
        call random_number(grad)
        call random_number(h_diag)
        h_diag = h_diag + 1.0_rp
        do i_case = 1, size(objects)
            associate (object => objects(i_case))
                object%n_param = n_param
                object%settings = default_s_gek_settings
                object%error_origin = error_s_gek
                call arm_host_context(object%settings)
                object%change_reference => mock_change_reference
                object%n_points = 0
                object%n_covariance = 0
                object%loosen_factor = 1.0_rp
                allocate(object%h_diag(n_param), object%kappa_list(n_param, 0), &
                         object%disp_list(n_param, 0), object%obj_func_list(0), &
                         object%grad_list(n_param, 0), &
                         object%local_grad_list(n_param, 0), object%y_list(n_param, 0))
                if (i_case == 2) host_context%mock_error = 7
                call object%add(kappa, 1.0_rp, grad, h_diag, error)
            end associate
            if (i_case == 1) then
                if (error /= 0) then
                    write(stderr, *) "test_add_s_gek failed: Produced error."
                    test_add_s_gek = .false.
                end if
                if (host_context%n_calls /= 1) then
                    write(stderr, *) "test_add_s_gek failed: Change of reference "// &
                        "function not called with the host context."
                    test_add_s_gek = .false.
                end if
                if (.not. host_context_reached("add_s_gek")) test_add_s_gek = .false.
            else if (error /= 7) then
                write(stderr, *) "test_add_s_gek failed: Error of the change of "// &
                    "reference function not passed on."
                test_add_s_gek = .false.
            end if
        end do

    end function test_add_s_gek

end module otr_s_gek_unit_tests
