! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol, n_param
    use otr_s_gek, only: change_reference_type
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(change_reference_type), pointer :: mock_change_reference_ptr => &
        mock_change_reference

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

    logical(c_bool) function test_s_gek_factory() bind(C)
        !
        ! this function tests the factory, which creates an S-GEK object after the
        ! sanity check of the settings, initializes it with the settings and the host's
        ! change of reference function and wires the solver settings with it and the
        ! record of the state
        !
        use opentrustregion, only: update_orbs_type, obj_func_type, solver_state_type, &
                                   default_solver_settings
        use opentrustregion_unit_tests, only: setup_settings, test_context_type
        use otr_common, only: extension_solver_settings_type, extension_level_type, &
                              extension_level, update_orbs_callback, obj_func_callback
        use otr_common_mock, only: mock_update_orbs, mock_obj_func
        use otr_s_gek, only: s_gek_factory, s_gek_settings_type, s_gek_type, &
                             default_s_gek_settings, &
                             s_gek_settings_uninitialized_warning_msg
        use otr_s_gek_test_reference, only: ref_s_gek_settings, assignment(=), &
                                            operator(/=)

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
        character(len=38), parameter :: case_names(3) = &
            [character(len=38) :: "for settings", "without settings", &
             "for settings that were not initialized"]

        ! assume tests pass
        test_s_gek_factory = .true.

        ! set the host's callback functions
        update_orbs_funptr => mock_update_orbs
        obj_func_funptr => mock_obj_func
        change_reference_funptr => mock_change_reference

        ! create objects for settings, for no settings and for settings that were not
        ! initialized, which are initialized after a warning through their logging
        ! function, each of which has to be initialized with the settings and the
        ! host's change of reference function and to wire the solver settings, which
        ! own it, with the host's orbital update and objective function and the record
        ! of the state
        do i_case = 1, size(case_names)
            allocate(solver_settings)
            solver_settings%solver_settings_type = default_solver_settings
            select case (i_case)
            case (1)
                settings = ref_s_gek_settings
                expected_settings = settings
            case (2)
                expected_settings = default_s_gek_settings
            case (3)
                call setup_settings(settings, context)
                settings%initialized = .false.
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
            level => extension_level(solver_settings%context)
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
                deallocate(solver_settings)
                cycle
            end if
            select type (object => level%extension)
            type is (s_gek_type)
                if (object%n_param /= n_param .or. &
                    object%settings /= expected_settings) then
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
            if (i_case == 3 .and. &
                index(context%log_message, &
                      " "//s_gek_settings_uninitialized_warning_msg) == 0) then
                write(stderr, *) "test_s_gek_factory failed: Warning not printed "// &
                    "for settings that were not initialized."
                test_s_gek_factory = .false.
            end if
            deallocate(solver_settings)
        end do

        ! create an object for parameters the sanity check rejects, which has to return
        ! its error and leave the solver settings untouched
        allocate(solver_settings)
        solver_settings%solver_settings_type = default_solver_settings
        call setup_settings(settings, context)
        call s_gek_factory(0_ip, update_orbs_funptr, obj_func_funptr, &
                           change_reference_funptr, update_orbs_s_gek_funptr, &
                           obj_func_s_gek_funptr, solver_settings, state, error, &
                           settings)
        if (error == 0 .or. associated(solver_settings%context) .or. &
            allocated(solver_settings%extensions)) then
            write(stderr, *) "test_s_gek_factory failed: Error not returned or "// &
                "solver settings changed for parameters the sanity check rejects."
            test_s_gek_factory = .false.
        end if
        deallocate(solver_settings)

        ! create objects for two separate solver settings, which have to be distinct
        ! objects reading their own record of the state, so that calculations with
        ! separate settings can be used alternately
        allocate(solver_settings, other_settings)
        solver_settings%solver_settings_type = default_solver_settings
        other_settings%solver_settings_type = default_solver_settings
        settings = ref_s_gek_settings
        call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                           change_reference_funptr, update_orbs_s_gek_funptr, &
                           obj_func_s_gek_funptr, solver_settings, state, error, &
                           settings)
        call s_gek_factory(n_param, update_orbs_funptr, obj_func_funptr, &
                           change_reference_funptr, update_orbs_s_gek_funptr, &
                           obj_func_s_gek_funptr, other_settings, other_state, error, &
                           settings)
        level => extension_level(solver_settings%context)
        other_level => extension_level(other_settings%context)
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
            if (error == 0) then
                write(stderr, *) "test_s_gek_sanity_check failed: No error "// &
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

    logical(c_bool) function test_init_s_gek_settings() bind(C)
        !
        ! this function tests the subroutine which initializes the S-GEK settings
        !
        use otr_s_gek, only: s_gek_settings_type, default_s_gek_settings
        use otr_s_gek_test_reference, only: operator(/=)

        type, extends(s_gek_settings_type) :: extended_settings_type
        end type

        type(s_gek_settings_type) :: settings
        type(extended_settings_type) :: extended_settings
        integer(ip) :: error

        ! assume tests pass
        test_init_s_gek_settings = .true.

        ! initialize settings, which have to be set to the default values
        call settings%init(error)
        if (error /= 0) then
            write(stderr, *) "test_init_s_gek_settings failed: Produced error."
            test_init_s_gek_settings = .false.
        end if
        if (settings /= default_s_gek_settings) then
            write(stderr, *) "test_init_s_gek_settings failed: Settings not "// &
                "initialized correctly."
            test_init_s_gek_settings = .false.
        end if

        ! initialize settings of a type extending these settings, which cannot be set to
        ! the default values, and check that an error is returned
        call extended_settings%init(error)
        if (error /= 1) then
            write(stderr, *) "test_init_s_gek_settings failed: No error returned "// &
                "for extended settings type."
            test_init_s_gek_settings = .false.
        end if

    end function test_init_s_gek_settings

    logical(c_bool) function test_init_s_gek() bind(C)
        !
        ! this function tests the subroutine which initializes an S-GEK object
        !
        use otr_s_gek, only: s_gek_settings_type, s_gek_type
        use otr_s_gek_test_reference, only: ref_s_gek_settings, assignment(=), &
                                            operator(/=)

        type(s_gek_type) :: object
        type(s_gek_settings_type) :: settings

        ! assume tests pass
        test_init_s_gek = .true.

        ! leave points, a covariance size, a loosened factor and arrays of another
        ! dimension from an earlier setup behind and initialize the object, which has to
        ! discard these, perform the orbital update and the Hessian linear
        ! transformation itself and take over the settings and the number of parameters
        object%n_points = 2
        object%n_covariance = 5
        object%loosen_factor = 2.0_rp
        allocate(object%h_diag(n_param + 1), object%kappa_list(n_param + 1, 2), &
                 object%disp_list(n_param + 1, 1), object%obj_func_list(2), &
                 object%grad_list(n_param + 1, 2), &
                 object%local_grad_list(n_param + 1, 2), object%y_list(n_param + 1, 2))
        settings = ref_s_gek_settings
        call object%init(n_param, settings)
        if (.not. all([object%provides%update_orbs, object%provides%hess_x]) .or. &
            any([object%provides%obj_func, object%provides%precond, &
                 object%provides%precond_pd, object%provides%project, &
                 object%provides%get_extra_trial_vectors])) then
            write(stderr, *) "test_init_s_gek failed: Roles the object performs "// &
                "set wrong."
            test_init_s_gek = .false.
        end if
        if (object%n_param /= n_param .or. object%settings /= settings) then
            write(stderr, *) "test_init_s_gek failed: Settings or number of "// &
                "parameters not taken over."
            test_init_s_gek = .false.
        end if
        if (object%n_points /= 0 .or. object%n_covariance /= 0 .or. &
            abs(object%loosen_factor - 1.0_rp) > tol) then
            write(stderr, *) "test_init_s_gek failed: Points, covariance size or "// &
                "loosening factor of an earlier setup kept."
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
        ! context of its settings, while its numerics are left to later tests
        !
        use otr_s_gek, only: s_gek_type, default_s_gek_settings
        use test_reference, only: arm_host_context, host_context, host_context_reached

        type(s_gek_type) :: object
        real(rp) :: kappa(n_param), grad(n_param), h_diag(n_param)
        integer(ip) :: i, error

        ! assume tests pass
        test_add_s_gek = .true.

        ! add a first point to an object as after its initialization, without the
        ! routine that initializes it, for which the routine moves the stored points to
        ! the new reference
        object%n_param = n_param
        object%settings = default_s_gek_settings
        call arm_host_context(object%settings)
        object%change_reference => mock_change_reference
        object%n_points = 0
        object%n_covariance = 0
        object%loosen_factor = 1.0_rp
        allocate(object%h_diag(n_param), object%kappa_list(n_param, 0), &
                 object%disp_list(n_param, 0), object%obj_func_list(0), &
                 object%grad_list(n_param, 0), object%local_grad_list(n_param, 0), &
                 object%y_list(n_param, 0))
        kappa = [(0.1_rp * i, i=1, n_param)]
        grad = [(0.3_rp * i - 0.5_rp, i=1, n_param)]
        h_diag = [(1.0_rp + i, i=1, n_param)]
        call object%add(kappa, 1.0_rp, grad, h_diag, error)
        if (error /= 0) then
            write(stderr, *) "test_add_s_gek failed: Produced error."
            test_add_s_gek = .false.
        end if
        if (host_context%n_calls /= 1) then
            write(stderr, *) "test_add_s_gek failed: Change of reference function "// &
                "not called with the host context."
            test_add_s_gek = .false.
        end if
        if (.not. host_context_reached("add_s_gek")) test_add_s_gek = .false.

    end function test_add_s_gek

end module otr_s_gek_unit_tests
