! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol, tol_c
    use otr_arh_c_interface, only: evaluate_dm_os_c_type, evaluate_dm_cs_c_type
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_funptr, c_funloc

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(evaluate_dm_cs_c_type), pointer :: mock_arh_evaluate_dm_cs_ptr => &
        mock_arh_evaluate_dm_cs
    procedure(evaluate_dm_os_c_type), pointer :: mock_arh_evaluate_dm_os_ptr => &
        mock_arh_evaluate_dm_os

contains

    function mock_arh_evaluate_dm_cs(dm_ao, energy, fock, v_nonlinear) result(error) &
        bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function with a separate non-linear potential contribution for the
        ! closed-shell case
        !
        use otr_oao_test_reference, only: n_ao
        use otr_arh_test_reference, only: evaluate_dm_cs_factors
        use otr_oao_unit_tests, only: record_mock_call

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*), v_nonlinear(*)
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(v_nonlinear)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) &
            fock(:flat_len) = evaluate_dm_cs_factors(1) * dm_ao(:flat_len)
        if (present(v_nonlinear)) &
            v_nonlinear(:flat_len) = evaluate_dm_cs_factors(2) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_arh_evaluate_dm_cs

    function mock_arh_evaluate_dm_os(dm_ao, energy, fock, v_same_spin, &
                                     v_opposite_spin, v_nonlinear) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function with separate same- and opposite-spin potential contributions and a
        ! non-linear potential contribution for the open-shell case
        !
        use otr_oao_test_reference, only: n_ao, n_particle
        use otr_arh_test_reference, only: evaluate_dm_os_factors
        use otr_oao_unit_tests, only: record_mock_call

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*), v_same_spin(*), &
                                             v_opposite_spin(*), v_nonlinear(*)
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2 * n_particle

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + &
                              merge(2_ip, 0_ip, present(v_same_spin)) + &
                              merge(4_ip, 0_ip, present(v_opposite_spin)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) &
            fock(:flat_len) = evaluate_dm_os_factors(1) * dm_ao(:flat_len)
        if (present(v_same_spin)) &
            v_same_spin(:flat_len) = evaluate_dm_os_factors(2) * dm_ao(:flat_len)
        if (present(v_opposite_spin)) &
            v_opposite_spin(:flat_len) = evaluate_dm_os_factors(3) * dm_ao(:flat_len)
        if (present(v_nonlinear)) &
            v_nonlinear(:flat_len) = evaluate_dm_os_factors(4) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_arh_evaluate_dm_os

    logical(c_bool) function test_arh_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH factory
        !
        use otr_arh_c_interface, only: arh_settings_type_c, arh_factory_cs, &
                                       arh_factory_os, arh_factory_c_wrapper
        use otr_arh_mock, only: mock_arh_factory_cs, mock_arh_factory_os, test_passed
        use otr_arh_test_reference, only: assignment(=), ref_arh_settings
        use otr_oao_test_reference, only: n_ao, n_particle, n_ao_c
        use c_interface_unit_tests, only: mock_logger, test_logger
        use test_reference, only: test_obj_func_c_funptr, test_update_orbs_c_funptr, &
                                  test_precond_c_funptr, test_precond_pd_c_funptr, &
                                  test_project_c_funptr

        real(c_rp), allocatable :: ao_overlap_c(:, :), dm_ao_2d_c(:, :), &
                                   dm_ao_3d_c(:, :, :)
        type(c_funptr) :: evaluate_dm_c_funptr, obj_func_arh_c_funptr, &
                          update_orbs_arh_c_funptr, precond_arh_c_funptr, &
                          precond_pd_arh_c_funptr, project_arh_c_funptr
        type(arh_settings_type_c) :: settings_c
        integer(c_ip) :: n_particle_c, error_c

        ! assume tests pass
        test_arh_factory_c_wrapper = .true.

        ! number of particles
        n_particle_c = 1_c_ip

        ! inject mock functions
        arh_factory_cs => mock_arh_factory_cs
        arh_factory_os => mock_arh_factory_os

        ! allocate and initialize arrays
        allocate(dm_ao_2d_c(n_ao, n_ao), ao_overlap_c(n_ao, n_ao))
        dm_ao_2d_c = 1.0_c_rp
        ao_overlap_c = 2.0_c_rp

        ! get C function pointers to Fortran functions
        evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_cs)

        ! associate optional settings with values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger)

        ! initialize logger logical
        test_logger = .true.

        ! call ARH orbital updating factory C wrapper for closed-shell case
        error_c = arh_factory_c_wrapper( &
            dm_ao_2d_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
            obj_func_arh_c_funptr, update_orbs_arh_c_funptr, precond_arh_c_funptr, &
            precond_pd_arh_c_funptr, project_arh_c_funptr, settings_c)

        ! check if logging subroutine was correctly called
        if (.not. test_logger) then
            test_arh_factory_c_wrapper = .false.
            write(stderr, *) "test_arh_factory_c_wrapper failed: Called logging "// &
                "subroutine wrong."
        end if

        ! check if output variables are as expected
        if (error_c /= 0) then
            test_arh_factory_c_wrapper = .false.
            write(stderr, *) "test_arh_factory_c_wrapper failed: Returned error "// &
                "code wrong."
        end if

        ! test returned objective function
        test_arh_factory_c_wrapper = &
            test_arh_factory_c_wrapper .and. &
            test_obj_func_c_funptr(obj_func_arh_c_funptr, "arh_factory_c_wrapper", &
                                   " by returned objective function")

        ! test returned orbital updating function
        test_arh_factory_c_wrapper = &
            test_arh_factory_c_wrapper .and. test_update_orbs_c_funptr( &
                update_orbs_arh_c_funptr, "arh_factory_c_wrapper", &
                " by returned orbital updating function")

        ! check if density matrix was updated
        if (any(abs(dm_ao_2d_c - 2.0_c_rp) > tol)) then
            test_arh_factory_c_wrapper = .false.
            write(stderr, *) "test_arh_factory_c_wrapper failed: Density matrix "// &
                "not updated correctly by returned orbital updating function."
        end if
        deallocate(dm_ao_2d_c)

        ! test returned level-shifted preconditioner function
        test_arh_factory_c_wrapper = &
            test_arh_factory_c_wrapper .and. &
            test_precond_c_funptr(precond_arh_c_funptr, "arh_factory_c_wrapper", &
                                  " by returned level-shifted preconditioner function")

        ! test returned positive-definite preconditioner function
        test_arh_factory_c_wrapper = &
            test_arh_factory_c_wrapper .and. test_precond_pd_c_funptr( &
                precond_pd_arh_c_funptr, "arh_factory_c_wrapper", &
                " by returned positive-definite preconditioner function")

        ! test returned projection function
        test_arh_factory_c_wrapper = &
            test_arh_factory_c_wrapper .and. &
            test_project_c_funptr(project_arh_c_funptr, "arh_factory_c_wrapper", &
                                  " by returned projection function")

        ! check if test has passed
        test_arh_factory_c_wrapper = test_arh_factory_c_wrapper .and. test_passed

        ! number of particles
        n_particle_c = 2_c_ip

        ! allocate and initialize 2D density matrix
        allocate(dm_ao_3d_c(n_ao, n_ao, n_particle))
        dm_ao_3d_c = 1.0_c_rp

        ! get C function pointers to Fortran functions
        evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_os)

        ! call ARH orbital updating factory C wrapper for open-shell case
        error_c = arh_factory_c_wrapper( &
            dm_ao_3d_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
            obj_func_arh_c_funptr, update_orbs_arh_c_funptr, precond_arh_c_funptr, &
            precond_pd_arh_c_funptr, project_arh_c_funptr, settings_c)

        ! deallocate arrays
        deallocate(dm_ao_3d_c, ao_overlap_c)

        ! check if tests have passed
        test_arh_factory_c_wrapper = test_arh_factory_c_wrapper .and. test_passed

    end function test_arh_factory_c_wrapper

    logical(c_bool) function test_evaluate_dm_cs_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating C
        ! function with a separate non-linear potential contribution for the
        ! closed-shell case
        !
        use otr_arh, only: evaluate_dm_cs_type
        use otr_arh_c_interface, only: evaluate_dm_cs_before_wrapping, &
                                       evaluate_dm_cs_f_wrapper
        use otr_arh_test_reference, only: test_evaluate_dm_cs_funptr
        use otr_oao_unit_tests, only: mock_requests

        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        integer(ip) :: request

        ! inject mock subroutine
        evaluate_dm_cs_before_wrapping => mock_arh_evaluate_dm_cs

        ! get pointer to subroutine
        evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) ::]
        test_evaluate_dm_cs_f_wrapper = test_evaluate_dm_cs_funptr( &
            evaluate_dm_cs_funptr, "evaluate_dm_cs_f_wrapper", "")
        if (size(mock_requests) /= 4 .or. &
            .not. all([(any(mock_requests == request), request = 0, 3)])) then
            write (stderr, *) "test_evaluate_dm_cs_f_wrapper failed: Outputs "// &
                "passed on to density matrix evaluating C function for "// &
                "closed-shell case wrong."
            test_evaluate_dm_cs_f_wrapper = .false.
        end if

    end function test_evaluate_dm_cs_f_wrapper

    logical(c_bool) function test_evaluate_dm_os_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating C
        ! function with same- and opposite-spin potential contributions for the
        ! open-shell case
        !
        use otr_arh, only: evaluate_dm_os_type
        use otr_arh_c_interface, only: evaluate_dm_os_before_wrapping, &
                                       evaluate_dm_os_f_wrapper
        use otr_arh_test_reference, only: test_evaluate_dm_os_funptr
        use otr_oao_unit_tests, only: mock_requests

        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        integer(ip) :: request

        ! inject mock subroutine
        evaluate_dm_os_before_wrapping => mock_arh_evaluate_dm_os

        ! get pointer to subroutine
        evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) ::]
        test_evaluate_dm_os_f_wrapper = test_evaluate_dm_os_funptr( &
            evaluate_dm_os_funptr, "evaluate_dm_os_f_wrapper", "")
        if (size(mock_requests) /= 16 .or. &
            .not. all([(any(mock_requests == request), request = 0, 15)])) then
            write (stderr, *) "test_evaluate_dm_os_f_wrapper failed: Outputs "// &
                "passed on to density matrix evaluating C function for open-shell "// &
                "case wrong."
            test_evaluate_dm_os_f_wrapper = .false.
        end if

    end function test_evaluate_dm_os_f_wrapper

    logical(c_bool) function test_update_orbs_arh_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH orbital update
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_arh_c_interface, only: update_orbs_arh_before_wrapping, &
                                       update_orbs_arh_c_wrapper
        use otr_common_mock, only: mock_update_orbs
        use test_reference, only: test_update_orbs_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        update_orbs_arh_before_wrapping => mock_update_orbs

        ! test orbital updating function
        test_update_orbs_arh_c_wrapper = test_update_orbs_c_funptr( &
            c_funloc(update_orbs_arh_c_wrapper), "update_orbs_arh_c_wrapper", "")

    end function test_update_orbs_arh_c_wrapper

    logical(c_bool) function test_hess_x_arh_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH Hessian linear transformation
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_arh_c_interface, only: hess_x_arh_before_wrapping, hess_x_arh_c_wrapper
        use otr_common_mock, only: mock_hess_x
        use test_reference, only: test_hess_x_c_funptr, n_param

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        hess_x_arh_before_wrapping => mock_hess_x

        ! test orbital updating function
        test_hess_x_arh_c_wrapper = test_hess_x_c_funptr( &
            c_funloc(hess_x_arh_c_wrapper), "hess_x_arh_c_wrapper", "")

    end function test_hess_x_arh_c_wrapper

    logical(c_bool) function test_init_arh_settings_c() bind(C)
        !
        ! this function tests that the ARH settings initialization routine correctly
        ! initializes all settings to their default values
        !
        use otr_arh_c_interface, only: arh_settings_type_c, init_arh_settings_c
        use otr_arh, only: default_arh_settings
        use otr_arh_test_reference, only: operator(/=)

        type(arh_settings_type_c) :: settings

        ! assume test passes
        test_init_arh_settings_c = .true.

        ! initialize settings
        call init_arh_settings_c(settings)

        ! check settings
        if (settings /= default_arh_settings) then
            write(stderr, *) "test_init_arh_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_arh_settings_c = .false. 
        end if

    end function test_init_arh_settings_c

    logical(c_bool) function test_arh_deconstructor_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH deconstructor
        !
        use otr_arh_c_interface, only: arh_deconstructor, arh_deconstructor_c_wrapper
        use otr_arh_mock, only: mock_arh_deconstructor, test_passed

        ! assume tests pass
        test_arh_deconstructor_c_wrapper = .true.

        ! inject mock functions
        arh_deconstructor => mock_arh_deconstructor

        ! initialize test logical
        test_passed = .false.

        ! call ARH orbital updating deconstructor C wrapper
        call arh_deconstructor_c_wrapper()

        ! check if test has passed
        test_arh_deconstructor_c_wrapper = test_passed

        ! check if test has passed
        if (.not. test_passed) then
            test_arh_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_arh_deconstructor_c_wrapper failed: "// &
                "Deconstructor called wrong."
        end if

    end function test_arh_deconstructor_c_wrapper

    logical(c_bool) function test_assign_arh_f_c() bind(C)
        !
        ! this function tests that the function that converts ARH settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh, only: arh_settings_type
        use otr_arh_test_reference, only: assignment(=), operator(/=), ref_arh_settings
        use c_interface_unit_tests, only: mock_logger, test_logger

        type(arh_settings_type_c) :: settings_c
        type(arh_settings_type) :: settings

        ! assume test passes
        test_assign_arh_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger)

        ! convert to Fortran settings
        settings = settings_c

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_arh_f_c = .false.
            write(stderr, *) "test_assign_arh_f_c failed: Logging function not "// &
                "associated with value."
        else
            test_logger = .true.
            call settings%logger("test")
            if (.not. test_logger) then
                test_assign_arh_f_c = .false.
                write(stderr, *) "test_assign_arh_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
        end if

        ! check against reference values
        if (settings /= ref_arh_settings) then
            write(stderr, *) "test_assign_arh_f_c failed: Settings not converted "// &
                "correctly."
            test_assign_arh_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            write(stderr, *) "test_assign_arh_f_c failed: Settings not marked as "// &
                "initialized."
            test_assign_arh_f_c = .false.
        end if

    end function test_assign_arh_f_c

    logical(c_bool) function test_assign_arh_c_f() bind(C)
        !
        ! this function tests that the function that converts ARH settings from Fortran
        ! to C correctly performs this conversion
        !
        use otr_arh, only: arh_settings_type
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh_test_reference, only: assignment(=), operator(/=), ref_arh_settings

        type(arh_settings_type)   :: settings
        type(arh_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_arh_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_arh_settings

        ! convert to C settings
        settings_c = settings

        ! check logging function
        if (c_associated(settings_c%logger)) then
            test_assign_arh_c_f = .false.
            write(stderr, *) "test_assign_arh_c_f failed: Logger function associated."
        end if

        ! check against reference values
        if (settings /= ref_arh_settings) then
            write(stderr, *) "test_assign_arh_c_f failed: Settings not converted "// &
                "correctly."
            test_assign_arh_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_arh_c_f = .false.
            write(stderr, *) "test_assign_arh_c_f failed: Settings not marked as "// &
                "initialized."
        end if

    end function test_assign_arh_c_f

end module otr_arh_c_interface_unit_tests
