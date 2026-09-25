! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol, tol_c
    use otr_oao_c_interface, only: evaluate_dm_c_type, get_response_c_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_funloc, c_associated

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(evaluate_dm_c_type), pointer :: &
        mock_evaluate_dm_cs_ptr => mock_evaluate_dm_cs, &
        mock_evaluate_dm_os_ptr => mock_evaluate_dm_os
    procedure(get_response_c_type), pointer :: &
        mock_get_response_cs_ptr => mock_get_response_cs, &
        mock_get_response_os_ptr => mock_get_response_os

contains

    function mock_evaluate_dm_cs(dm_ao, energy, fock, get_response_c_funptr) &
        result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function for 2D density matrices
        !
        use otr_oao_test_reference, only: n_ao, evaluate_dm_factors
        use otr_oao_unit_tests, only: record_mock_call

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*)
        type(c_funptr), intent(out), optional :: get_response_c_funptr
        integer(c_ip) :: error

        integer(ip) :: flat_len = n_ao**2

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(get_response_c_funptr)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) fock(:flat_len) = evaluate_dm_factors(1) * dm_ao(:flat_len)
        if (present(get_response_c_funptr)) &
            get_response_c_funptr = c_funloc(mock_get_response_cs)

        error = 0_c_ip

    end function mock_evaluate_dm_cs

    function mock_evaluate_dm_os(dm_ao, energy, fock, get_response_c_funptr) &
        result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function for 3D density matrices
        !
        use otr_oao_test_reference, only: n_ao, n_particle, evaluate_dm_factors
        use otr_oao_unit_tests, only: record_mock_call

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*)
        type(c_funptr), intent(out), optional :: get_response_c_funptr
        integer(c_ip) :: error

        integer(ip) :: flat_len = n_ao**2 * n_particle

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(get_response_c_funptr)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) fock(:flat_len) = evaluate_dm_factors(1) * dm_ao(:flat_len)
        if (present(get_response_c_funptr)) &
            get_response_c_funptr = c_funloc(mock_get_response_os)

        error = 0_c_ip

    end function mock_evaluate_dm_os

    function mock_get_response_cs(dm_ao, response) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the response C function for 2D
        ! density matrices
        !
        use otr_oao_test_reference, only: n_ao, evaluate_dm_factors

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out), target :: response(*)
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2

        response(:flat_len) = evaluate_dm_factors(2) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_get_response_cs

    function mock_get_response_os(dm_ao, response) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the response C function for 3D
        ! density matrices
        !
        use otr_oao_test_reference, only: n_ao, n_particle, evaluate_dm_factors

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out), target :: response(*)
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2 * n_particle

        response(:flat_len) = evaluate_dm_factors(2) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_get_response_os

    logical(c_bool) function test_oao_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO factory
        !
        use otr_oao_c_interface, only: oao_settings_type_c, oao_factory_cs, &
                                       oao_factory_os, oao_factory_c_wrapper
        use otr_oao_mock, only: mock_oao_factory_cs, mock_oao_factory_os, test_passed
        use otr_oao_test_reference, only: assignment(=), ref_oao_settings, n_ao, &
                                          n_particle, n_ao_c
        use c_interface_unit_tests, only: mock_logger, test_logger
        use test_reference, only: test_obj_func_c_funptr, test_update_orbs_c_funptr, &
                                  test_precond_c_funptr, test_precond_pd_c_funptr, &
                                  test_project_c_funptr, &
                                  test_get_extra_trial_vectors_c_funptr, ref_settings, &
                                  assignment(=), operator(/=)
        use c_interface, only: solver_settings_type_c

        real(c_rp), allocatable :: ao_overlap_c(:, :), dm_ao_2d_c(:, :), &
                                   dm_ao_3d_c(:, :, :)
        type(c_funptr) :: evaluate_dm_c_funptr, obj_func_oao_c_funptr, &
                          update_orbs_oao_c_funptr
        type(oao_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        integer(c_ip) :: n_particle_c, error_c

        ! assume tests pass
        test_oao_factory_c_wrapper = .true.

        ! number of particles
        n_particle_c = 1_c_ip

        ! inject mock functions
        oao_factory_cs => mock_oao_factory_cs
        oao_factory_os => mock_oao_factory_os

        ! allocate and initialize arrays
        allocate(dm_ao_2d_c(n_ao, n_ao), ao_overlap_c(n_ao, n_ao))
        dm_ao_2d_c = 1.0_c_rp
        ao_overlap_c = 2.0_c_rp

        ! get C function pointers to Fortran functions
        evaluate_dm_c_funptr = c_funloc(mock_evaluate_dm_cs)

        ! associate optional settings with values
        settings_c = ref_oao_settings
        settings_c%logger = c_funloc(mock_logger)

        ! initialize logger logical
        test_logger = .true.

        ! solver settings which are not yet initialized
        solver_settings_c%initialized = .false._c_bool

        ! call OAO orbital updating factory C wrapper for closed-shell case
        error_c = oao_factory_c_wrapper( &
            dm_ao_2d_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
            obj_func_oao_c_funptr, update_orbs_oao_c_funptr, solver_settings_c, &
            settings_c)

        ! check if logging subroutine was correctly called
        if (.not. test_logger) then
            test_oao_factory_c_wrapper = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Called logging "// &
                "subroutine wrong."
        end if

        ! check if output variables are as expected
        if (error_c /= 0) then
            test_oao_factory_c_wrapper = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Returned error "// &
                "code wrong."
        end if

        ! test returned objective function
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. &
            test_obj_func_c_funptr(obj_func_oao_c_funptr, "oao_factory_c_wrapper", &
                                   " by returned objective function")

        ! test returned orbital updating function
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. test_update_orbs_c_funptr( &
                update_orbs_oao_c_funptr, "oao_factory_c_wrapper", &
                " by returned orbital updating function")

        ! check if density matrix was updated
        if (any(abs(dm_ao_2d_c - 2.0_c_rp) > tol)) then
            test_passed = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Density matrix "// &
                "not updated correctly by returned orbital updating function."
        end if
        deallocate(dm_ao_2d_c)

        ! determine if the solver settings were initialized
        if (.not. solver_settings_c%initialized) then
            test_oao_factory_c_wrapper = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Solver settings "// &
                "not initialized."
        end if

        ! test functions wired into solver settings
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. &
            test_precond_c_funptr(solver_settings_c%precond, "oao_factory_c_wrapper", &
                                  " by wired level-shifted preconditioner function")
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. test_precond_pd_c_funptr( &
                solver_settings_c%precond_pd, "oao_factory_c_wrapper", &
                " by wired positive-definite preconditioner function")
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. &
            test_project_c_funptr(solver_settings_c%project, "oao_factory_c_wrapper", &
                                  " by wired projection function")
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. test_get_extra_trial_vectors_c_funptr( &
                solver_settings_c%get_extra_trial_vectors, "oao_factory_c_wrapper", &
                " by wired extra trial vector function")

        ! test functions wired into stability check settings
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. test_precond_c_funptr( &
                solver_settings_c%stability_settings%precond, "oao_factory_c_wrapper", &
                " by wired stability check level-shifted preconditioner function")
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. test_project_c_funptr( &
                solver_settings_c%stability_settings%project, "oao_factory_c_wrapper", &
                " by wired stability check projection function")
        test_oao_factory_c_wrapper = &
            test_oao_factory_c_wrapper .and. test_get_extra_trial_vectors_c_funptr( &
                solver_settings_c%stability_settings%get_extra_trial_vectors, &
                "oao_factory_c_wrapper", &
                " by wired stability check extra trial vector function")

        ! check if test has passed
        test_oao_factory_c_wrapper = test_oao_factory_c_wrapper .and. test_passed

        ! number of particles
        n_particle_c = 2_c_ip

        ! allocate and initialize 2D density matrix
        allocate(dm_ao_3d_c(n_ao, n_ao, n_particle))
        dm_ao_3d_c = 1.0_c_rp

        ! get C function pointers to Fortran functions
        evaluate_dm_c_funptr = c_funloc(mock_evaluate_dm_os)

        ! set the now initialized solver settings to reference values, which have to
        ! be kept
        solver_settings_c = ref_settings

        ! call OAO orbital updating factory C wrapper for open-shell case
        error_c = oao_factory_c_wrapper( &
            dm_ao_3d_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
            obj_func_oao_c_funptr, update_orbs_oao_c_funptr, solver_settings_c, &
            settings_c)

        ! deallocate arrays
        deallocate(dm_ao_3d_c, ao_overlap_c)

        ! determine if the initialized solver settings were kept
        if (solver_settings_c /= ref_settings) then
            test_oao_factory_c_wrapper = .false.
            write(stderr, *) "test_oao_factory_c_wrapper failed: Initialized "// &
                "solver settings not kept."
        end if

        ! check if tests for dm_ao_3d_c and evaluate_dm_c_funptr have passed
        test_oao_factory_c_wrapper = test_oao_factory_c_wrapper .and. test_passed

    end function test_oao_factory_c_wrapper

    logical(c_bool) function test_evaluate_dm_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating
        ! function
        !
        use otr_oao, only: evaluate_dm_cs_type, evaluate_dm_os_type
        use otr_oao_c_interface, only: evaluate_dm_before_wrapping, &
                                       evaluate_dm_cs_f_wrapper, &
                                       evaluate_dm_os_f_wrapper
        use otr_oao_test_reference, only: test_evaluate_dm_cs_funptr, &
                                          test_evaluate_dm_os_funptr
        use otr_oao_unit_tests, only: mock_requests

        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        integer(ip) :: request

        ! inject mock subroutine
        evaluate_dm_before_wrapping => mock_evaluate_dm_cs

        ! get pointer to subroutine
        evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) ::]
        test_evaluate_dm_f_wrapper = test_evaluate_dm_cs_funptr( &
            evaluate_dm_cs_funptr, "evaluate_dm_f_wrapper", "")
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request = 0, 3)])) then
            write (stderr, *) "test_evaluate_dm_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for closed-shell case "// &
                "wrong."
            test_evaluate_dm_f_wrapper = .false.
        end if

        ! inject mock subroutine
        evaluate_dm_before_wrapping => mock_evaluate_dm_os

        ! get pointer to subroutine
        evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) ::]
        test_evaluate_dm_f_wrapper = &
            test_evaluate_dm_f_wrapper .and. test_evaluate_dm_os_funptr( &
                evaluate_dm_os_funptr, "evaluate_dm_f_wrapper", "")
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request = 0, 3)])) then
            write (stderr, *) "test_evaluate_dm_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for open-shell case wrong."
            test_evaluate_dm_f_wrapper = .false.
        end if

    end function test_evaluate_dm_f_wrapper

    logical(c_bool) function test_get_response_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the response function
        !
        use otr_oao, only: get_response_cs_type, get_response_os_type
        use otr_oao_c_interface, only: get_response_before_wrapping, &
                                       get_response_cs_f_wrapper, &
                                       get_response_os_f_wrapper
        use otr_oao_test_reference, only: test_get_response_cs_funptr, &
                                          test_get_response_os_funptr

        procedure(get_response_cs_type), pointer :: get_response_cs_funptr
        procedure(get_response_os_type), pointer :: get_response_os_funptr

        ! inject mock subroutine
        get_response_before_wrapping => mock_get_response_cs

        ! get pointer to subroutine
        get_response_cs_funptr => get_response_cs_f_wrapper

        ! test response function wrapper
        test_get_response_f_wrapper = test_get_response_cs_funptr( &
            get_response_cs_funptr, "get_response_f_wrapper", "")

        ! inject mock subroutine
        get_response_before_wrapping => mock_get_response_os

        ! get pointer to subroutine
        get_response_os_funptr => get_response_os_f_wrapper

        ! test response function wrapper
        test_get_response_f_wrapper = &
            test_get_response_f_wrapper .and. test_get_response_os_funptr( &
                get_response_os_funptr, "get_response_f_wrapper", "")

    end function test_get_response_f_wrapper

    logical(c_bool) function test_obj_func_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO objective function
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_oao_c_interface, only: obj_func_oao_before_wrapping, &
                                       obj_func_oao_c_wrapper
        use otr_oao_mock, only: mock_obj_func_oao
        use test_reference, only: test_obj_func_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        obj_func_oao_before_wrapping => mock_obj_func_oao

        ! test objective function
        test_obj_func_oao_c_wrapper = test_obj_func_c_funptr( &
            c_funloc(obj_func_oao_c_wrapper), "obj_func_oao_c_wrapper", "")

    end function test_obj_func_oao_c_wrapper

    logical(c_bool) function test_update_orbs_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO orbital update
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_oao_c_interface, only: update_orbs_oao_before_wrapping, &
                                       update_orbs_oao_c_wrapper
        use otr_common_mock, only: mock_update_orbs
        use test_reference, only: test_update_orbs_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        update_orbs_oao_before_wrapping => mock_update_orbs

        ! test orbital updating function
        test_update_orbs_oao_c_wrapper = test_update_orbs_c_funptr( &
            c_funloc(update_orbs_oao_c_wrapper), "update_orbs_oao_c_wrapper", "")

    end function test_update_orbs_oao_c_wrapper

    logical(c_bool) function test_hess_x_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO Hessian linear transformation
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_oao_c_interface, only: hess_x_oao_before_wrapping, hess_x_oao_c_wrapper
        use otr_common_mock, only: mock_hess_x
        use test_reference, only: test_hess_x_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        hess_x_oao_before_wrapping => mock_hess_x

        ! test orbital updating function
        test_hess_x_oao_c_wrapper = test_hess_x_c_funptr( &
            c_funloc(hess_x_oao_c_wrapper), "hess_x_oao_c_wrapper", "")

    end function test_hess_x_oao_c_wrapper

    logical(c_bool) function test_precond_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO level-shifted preconditioner
        ! subroutine
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_oao_c_interface, only: precond_oao_before_wrapping, &
                                       precond_oao_c_wrapper
        use otr_oao_mock, only: mock_precond_oao
        use test_reference, only: test_precond_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        precond_oao_before_wrapping => mock_precond_oao

        ! test preconditioner function
        test_precond_oao_c_wrapper = test_precond_c_funptr( &
            c_funloc(precond_oao_c_wrapper), "precond_oao_c_wrapper", "")

    end function test_precond_oao_c_wrapper

    logical(c_bool) function test_precond_pd_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO positive-definite
        ! preconditioner subroutine
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_oao_c_interface, only: precond_pd_oao_before_wrapping, &
                                       precond_pd_oao_c_wrapper
        use otr_oao_mock, only: mock_precond_pd_oao
        use test_reference, only: test_precond_pd_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        precond_pd_oao_before_wrapping => mock_precond_pd_oao

        ! test preconditioner function
        test_precond_pd_oao_c_wrapper = test_precond_pd_c_funptr( &
            c_funloc(precond_pd_oao_c_wrapper), "precond_pd_oao_c_wrapper", "")

    end function test_precond_pd_oao_c_wrapper

    logical(c_bool) function test_project_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO projection subroutine
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_oao_c_interface, only: project_oao_before_wrapping, &
                                       project_oao_c_wrapper
        use otr_oao_mock, only: mock_project_oao
        use test_reference, only: test_project_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        project_oao_before_wrapping => mock_project_oao

        ! test orbital updating function
        test_project_oao_c_wrapper = test_project_c_funptr( &
            c_funloc(project_oao_c_wrapper), "project_oao_c_wrapper", "")

    end function test_project_oao_c_wrapper

    logical(c_bool) function test_get_extra_trial_vectors_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO extra trial vector subroutine
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_oao_c_interface, only: get_extra_trial_vectors_oao_before_wrapping, &
                                       get_extra_trial_vectors_oao_c_wrapper
        use otr_oao_mock, only: mock_get_extra_trial_vectors_oao
        use test_reference, only: test_get_extra_trial_vectors_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        get_extra_trial_vectors_oao_before_wrapping => mock_get_extra_trial_vectors_oao

        ! test extra trial vector function
        test_get_extra_trial_vectors_oao_c_wrapper = &
            test_get_extra_trial_vectors_c_funptr( &
                c_funloc(get_extra_trial_vectors_oao_c_wrapper), &
                "get_extra_trial_vectors_oao_c_wrapper", "")

    end function test_get_extra_trial_vectors_oao_c_wrapper

    logical(c_bool) function test_init_oao_settings_c() bind(C)
        !
        ! this function tests that the OAO settings initialization routine correctly
        ! initializes all settings to their default values
        !
        use otr_oao_c_interface, only: oao_settings_type_c, init_oao_settings_c
        use otr_oao, only: default_oao_settings
        use otr_oao_test_reference, only: operator(/=)

        type(oao_settings_type_c) :: settings

        ! assume test passes
        test_init_oao_settings_c = .true.

        ! initialize settings
        call init_oao_settings_c(settings)

        ! check function pointers
        if (c_associated(settings%logger)) then
            write(stderr, *) "test_init_oao_settings_c failed: Function pointers "// &
                "should not be initialized."
            test_init_oao_settings_c = .false.
        end if

        ! check settings
        if (settings /= default_oao_settings) then
            write(stderr, *) "test_init_oao_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_oao_settings_c = .false. 
        end if

    end function test_init_oao_settings_c

    logical(c_bool) function test_oao_deconstructor_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the OAO deconstructor
        !
        use otr_oao_c_interface, only: oao_deconstructor, oao_deconstructor_c_wrapper
        use otr_oao_mock, only: mock_oao_deconstructor, test_passed

        ! assume tests pass
        test_oao_deconstructor_c_wrapper = .true.

        ! inject mock functions
        oao_deconstructor => mock_oao_deconstructor

        ! initialize test logical
        test_passed = .false.

        ! call OAO orbital updating deconstructor C wrapper
        call oao_deconstructor_c_wrapper()

        ! check if test has passed
        test_oao_deconstructor_c_wrapper = test_passed

        ! check if test has passed
        if (.not. test_passed) then
            test_oao_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_oao_deconstructor_c_wrapper failed: "// &
                "Deconstructor called wrong."
        end if

    end function test_oao_deconstructor_c_wrapper

    logical(c_bool) function test_assign_oao_f_c() bind(C)
        !
        ! this function tests that the function that converts OAO settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao, only: oao_settings_type
        use otr_oao_test_reference, only: assignment(=), ref_oao_settings, operator(/=)
        use c_interface_unit_tests, only: mock_logger, test_logger

        type(oao_settings_type_c) :: settings_c
        type(oao_settings_type) :: settings

        ! assume test passes
        test_assign_oao_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_oao_settings
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings
        settings = settings_c

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_oao_f_c = .false.
            write(stderr, *) "test_assign_oao_f_c failed: Logging function not "// &
                "associated with value."
        else
            test_logger = .true.
            call settings%logger("test")
            if (.not. test_logger) then
                test_assign_oao_f_c = .false.
                write(stderr, *) "test_assign_oao_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
        end if

        ! check against reference values
        if (settings /= ref_oao_settings) then
            write(stderr, *) "test_assign_oao_f_c failed: Settings not converted "// &
                "correctly."
            test_assign_oao_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            write(stderr, *) "test_assign_oao_f_c failed: Settings not marked as "// &
                "initialized."
            test_assign_oao_f_c = .false.
        end if

    end function test_assign_oao_f_c

    logical(c_bool) function test_assign_oao_c_f() bind(C)
        !
        ! this function tests that the function that converts OAO settings from Fortran
        ! to C correctly performs this conversion
        !
        use otr_oao, only: oao_settings_type
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao_test_reference, only: ref_oao_settings, assignment(=), operator(/=)

        type(oao_settings_type)   :: settings
        type(oao_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_oao_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_oao_settings

        ! convert to C settings
        settings_c = settings

        ! check that callback function pointers are not associated
        if (c_associated(settings_c%logger)) then
            test_assign_oao_c_f = .false.
            write(stderr, *) "test_assign_oao_c_f failed: Logger function associated."
        end if

        ! check against reference values
        if (settings /= ref_oao_settings) then
            write(stderr, *) "test_assign_oao_c_f failed: Settings not converted "// &
                "correctly."
            test_assign_oao_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_oao_c_f = .false.
            write(stderr, *) "test_assign_oao_c_f failed: Settings not marked as "// &
                "initialized."
        end if

    end function test_assign_oao_c_f

end module otr_oao_c_interface_unit_tests
