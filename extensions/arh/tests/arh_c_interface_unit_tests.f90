! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol, tol_c, n_param
    use otr_arh_test_reference, only: n_particle, n_ao, n_ao_c
    use otr_arh_c_interface, only: get_energy_c_type, get_fock_c_type, &
                                   get_fock_jk_c_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_funloc, c_associated, &
                                           c_f_procpointer

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(get_energy_c_type), pointer :: mock_get_energy_2d_ptr => &
                                             mock_get_energy_2d, &
                                             mock_get_energy_3d_ptr => &
                                             mock_get_energy_3d
    procedure(get_fock_c_type), pointer :: mock_get_fock_ptr => mock_get_fock
    procedure(get_fock_jk_c_type), pointer :: mock_get_fock_jk_ptr => mock_get_fock_jk

contains

    function mock_get_energy_2d(dm_ao, energy) result(error) bind(C)
        !
        ! this function is a test function for the energy C function for 2D density
        ! matrices
        !
        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        integer(c_ip) :: error

        energy = sum(dm_ao(:n_ao ** 2))

        error = 0

    end function mock_get_energy_2d

    function mock_get_energy_3d(dm_ao, energy) result(error) bind(C)
        !
        ! this function is a test function for the energy C function for 3D density
        ! matrices
        !
        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        integer(c_ip) :: error

        energy = sum(dm_ao(:n_ao ** 2 * n_particle))

        error = 0

    end function mock_get_energy_3d

    function mock_get_fock(dm_ao, energy, fock) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the Fock matrix C function
        !
        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), target :: fock(*)
        integer(c_ip) :: error

        energy = sum(dm_ao(:n_ao ** 2))

        fock(:n_ao ** 2) = 2 * dm_ao(:n_ao ** 2)

        error = 0_c_ip

    end function mock_get_fock

    function mock_get_fock_jk(dm_ao, energy, fock, coulomb, exchange) result(error) &
        bind(C)
        !
        ! this subroutine is a test subroutine for the Fock, Coulomb, and exchange 
        ! matrix C function
        !
        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), target :: fock(*), coulomb(*), exchange(*)
        integer(c_ip) :: error

        energy = sum(dm_ao(:n_ao ** 2 * n_particle))

        fock(:n_ao ** 2 * n_particle) = 2 * dm_ao(:n_ao ** 2 * n_particle)
        coulomb(:n_ao ** 2 * n_particle) = 3 * dm_ao(:n_ao ** 2 * n_particle)
        exchange(:n_ao ** 2 * n_particle) = 4 * dm_ao(:n_ao ** 2 * n_particle)

        error = 0_c_ip

    end function mock_get_fock_jk

    logical(c_bool) function test_arh_factory_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH factory
        !
        use otr_arh_c_interface, only: arh_settings_type_c, arh_factory_closed_shell, &
                                       arh_factory_open_shell, arh_factory_c_wrapper
        use otr_arh_mock, only: mock_arh_factory_closed_shell, &
                                mock_arh_factory_open_shell, test_passed
        use otr_arh_test_reference, only: assignment(=), ref_arh_settings
        use c_interface_unit_tests, only: mock_logger, test_logger
        use test_reference, only: test_obj_func_c_funptr, test_update_orbs_c_funptr, &
                                  test_project_c_funptr

        real(c_rp), allocatable :: ao_overlap_c(:, :), dm_ao_2d_c(:, :), &
                                   dm_ao_3d_c(:, :, :)
        type(c_funptr) :: get_energy_c_funptr, get_fock_c_funptr, &
                          obj_func_arh_c_funptr, update_orbs_arh_c_funptr, &
                          project_arh_c_funptr
        type(arh_settings_type_c) :: settings_c
        integer(c_ip) :: n_particle_c, error_c

        ! assume tests pass
        test_arh_factory_c_wrapper = .true.

        ! number of particles
        n_particle_c = 1_c_ip

        ! inject mock functions
        arh_factory_closed_shell => mock_arh_factory_closed_shell
        arh_factory_open_shell => mock_arh_factory_open_shell

        ! allocate and initialize arrays
        allocate(dm_ao_2d_c(n_ao, n_ao), ao_overlap_c(n_ao, n_ao))
        dm_ao_2d_c = 1.0_c_rp
        ao_overlap_c = 2.0_c_rp

        ! get C function pointers to Fortran functions
        get_energy_c_funptr = c_funloc(mock_get_energy_2d)
        get_fock_c_funptr = c_funloc(mock_get_fock)

        ! associate optional settings with values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger)

        ! initialize logger logical
        test_logger = .true.

        ! call ARH orbital updating factory C wrapper for closed-shell case
        error_c = arh_factory_c_wrapper(dm_ao_2d_c, ao_overlap_c, n_particle_c, &
                                        n_ao_c, get_energy_c_funptr, &
                                        get_fock_c_funptr, obj_func_arh_c_funptr, &
                                        update_orbs_arh_c_funptr, &
                                        project_arh_c_funptr, settings_c)

        ! deallocate 2D density matrix
        deallocate(dm_ao_2d_c)

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
        test_arh_factory_c_wrapper = test_arh_factory_c_wrapper .and. &
            test_obj_func_c_funptr(obj_func_arh_c_funptr, "arh_factory_c_wrapper", &
                                   " by returned objective function")

        ! test returned orbital updating function
        test_arh_factory_c_wrapper = test_arh_factory_c_wrapper .and. &
            test_update_orbs_c_funptr(update_orbs_arh_c_funptr, &
                                      "arh_factory_c_wrapper", &
                                      " by returned orbital updating function")

        ! test returned projection function
        test_arh_factory_c_wrapper = test_arh_factory_c_wrapper .and. &
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
        get_energy_c_funptr = c_funloc(mock_get_energy_3d)
        get_fock_c_funptr = c_funloc(mock_get_fock_jk)

        ! call ARH orbital updating factory C wrapper for open-shell case
        error_c = arh_factory_c_wrapper(dm_ao_3d_c, ao_overlap_c, n_particle_c, &
                                        n_ao_c, get_energy_c_funptr, &
                                        get_fock_c_funptr, obj_func_arh_c_funptr, &
                                        update_orbs_arh_c_funptr, &
                                        project_arh_c_funptr, settings_c)

        ! deallocate arrays
        deallocate(dm_ao_3d_c, ao_overlap_c)

        ! check if tests for dm_ao_3d_c, get_energy_c_funptr and get_fock_c_funptr have 
        ! passed
        test_arh_factory_c_wrapper = test_arh_factory_c_wrapper .and. test_passed

    end function test_arh_factory_c_wrapper

    logical(c_bool) function test_get_energy_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the energy function
        !
        use otr_arh, only: get_energy_closed_shell_type, get_energy_open_shell_type
        use otr_arh_c_interface, only: get_energy_before_wrapping, &
                                       get_energy_closed_shell_f_wrapper, &
                                       get_energy_open_shell_f_wrapper
        use otr_arh_test_reference, only: test_get_energy_closed_shell_funptr, &
                                          test_get_energy_open_shell_funptr

        procedure(get_energy_closed_shell_type), pointer :: &
            get_energy_closed_shell_funptr
        procedure(get_energy_open_shell_type), pointer :: get_energy_open_shell_funptr

        ! inject mock function
        get_energy_before_wrapping => mock_get_energy_2d

        ! get pointer to function
        get_energy_closed_shell_funptr => get_energy_closed_shell_f_wrapper

        ! test energy wrapper
        test_get_energy_f_wrapper = &
            test_get_energy_closed_shell_funptr(get_energy_closed_shell_funptr, &
                                                "get_energy_closed_shell_f_wrapper", &
                                                " for closed-shell case")

        ! inject mock function
        get_energy_before_wrapping => mock_get_energy_3d

        ! get pointer to function
        get_energy_open_shell_funptr => get_energy_open_shell_f_wrapper

        ! test energy wrapper
        test_get_energy_f_wrapper = test_get_energy_f_wrapper .and. &
            test_get_energy_open_shell_funptr(get_energy_open_shell_funptr, &
                                              "get_energy_open_shell_f_wrapper", &
                                              " for open-shell case")

    end function test_get_energy_f_wrapper

    logical(c_bool) function test_get_fock_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the Fock matrix function
        !
        use otr_arh, only: get_fock_type
        use otr_arh_c_interface, only: get_fock_before_wrapping, get_fock_f_wrapper
        use otr_arh_test_reference, only: test_get_fock_funptr

        procedure(get_fock_type), pointer :: get_fock_funptr

        ! inject mock subroutine
        get_fock_before_wrapping => mock_get_fock

        ! get pointer to subroutine
        get_fock_funptr => get_fock_f_wrapper

        ! test Fock matrix wrapper
        test_get_fock_f_wrapper = &
            test_get_fock_funptr(get_fock_funptr, "get_fock_f_wrapper", "")

    end function test_get_fock_f_wrapper

    logical(c_bool) function test_get_fock_jk_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the Fock, Coulomb and exchange 
        ! matrix function
        !
        use otr_arh, only: get_fock_jk_type
        use otr_arh_c_interface, only: get_fock_jk_before_wrapping, &
            get_fock_jk_f_wrapper
        use otr_arh_test_reference, only: test_get_fock_jk_funptr

        procedure(get_fock_jk_type), pointer :: get_fock_jk_funptr

        ! inject mock subroutine
        get_fock_jk_before_wrapping => mock_get_fock_jk

        ! get pointer to subroutine
        get_fock_jk_funptr => get_fock_jk_f_wrapper

        ! test Fock matrix wrapper
        test_get_fock_jk_f_wrapper = &
            test_get_fock_jk_funptr(get_fock_jk_funptr, "get_fock_jk_f_wrapper", "")
        
    end function test_get_fock_jk_f_wrapper

    logical(c_bool) function test_update_orbs_arh_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH orbital update
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_arh_c_interface, only: update_orbs_arh_before_wrapping, &
                                       update_orbs_arh_c_wrapper
        use otr_common_mock, only: mock_update_orbs
        use test_reference, only: test_update_orbs_c_funptr

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        update_orbs_arh_before_wrapping => mock_update_orbs

        ! test orbital updating function
        test_update_orbs_arh_c_wrapper = &
            test_update_orbs_c_funptr(c_funloc(update_orbs_arh_c_wrapper), &
                                      "update_orbs_arh_c_wrapper", "")

    end function test_update_orbs_arh_c_wrapper

    logical(c_bool) function test_hess_x_arh_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH Hessian linear transformation
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_arh_c_interface, only: hess_x_arh_before_wrapping, hess_x_arh_c_wrapper
        use otr_common_mock, only: mock_hess_x
        use test_reference, only: test_hess_x_c_funptr

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        hess_x_arh_before_wrapping => mock_hess_x

        ! test orbital updating function
        test_hess_x_arh_c_wrapper = &
            test_hess_x_c_funptr(c_funloc(hess_x_arh_c_wrapper), &
                                 "hess_x_arh_c_wrapper", "")

    end function test_hess_x_arh_c_wrapper

    logical(c_bool) function test_project_arh_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH projection subroutine
        !
        use otr_common_c_interface, only: n_param_global => n_param
        use otr_arh_c_interface, only: project_arh_before_wrapping, &
                                       project_arh_c_wrapper
        use otr_arh_mock, only: mock_project_arh
        use test_reference, only: test_project_c_funptr

        ! set global number of parameters for assumed size arrays
        n_param_global = n_param

        ! inject mock subroutine
        project_arh_before_wrapping => mock_project_arh

        ! test orbital updating function
        test_project_arh_c_wrapper = &
            test_project_c_funptr(c_funloc(project_arh_c_wrapper), &
                                  "project_arh_c_wrapper", "")

    end function test_project_arh_c_wrapper

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

        ! check function pointers
        if (c_associated(settings%logger)) then
            write(stderr, *) "test_init_arh_settings_c failed: Function pointers "// &
                "should not be initialized."
            test_init_arh_settings_c = .false.
        end if

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
        use otr_arh_c_interface, only: arh_deconstructor_closed_shell, &
                                       arh_deconstructor_open_shell, &
                                       arh_deconstructor_c_wrapper, &
                                       n_particle_global => n_particle, &
                                       n_ao_global => n_ao
        use otr_arh_mock, only: mock_arh_deconstructor_closed_shell, &
                                mock_arh_deconstructor_open_shell

        real(c_rp), allocatable :: dm_ao_2d_c(:, :), dm_ao_3d_c(:, :, :)
        integer(c_ip) :: error_c

        ! assume tests pass
        test_arh_deconstructor_c_wrapper = .true.

        ! set global number of particles and AOs for assumed size arrays
        n_particle_global = n_particle
        n_ao_global = n_ao

        ! inject mock functions
        arh_deconstructor_closed_shell => mock_arh_deconstructor_closed_shell
        arh_deconstructor_open_shell => mock_arh_deconstructor_open_shell

        ! allocate 2D density matrix
        allocate(dm_ao_2d_c(n_ao, n_ao))

        ! call ARH orbital updating deconstructor C wrapper
        error_c = arh_deconstructor_c_wrapper(dm_ao_2d_c)

        ! check if error is as expected
        if (error_c /= 0) then
            test_arh_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_arh_deconstructor_c_wrapper failed: Returned error."
        end if

        ! check if 2D density matrix is as expected
        if (any(abs(dm_ao_2d_c - 1.0_c_rp) > tol_c)) then
            test_arh_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_arh_deconstructor_c_wrapper failed: Returned AO "// &
                "density matrix for closed-shell case wrong."
        end if

        ! deallocate arrays
        deallocate(dm_ao_2d_c)

        ! allocate 3D density matrix
        allocate(dm_ao_3d_c(n_ao, n_ao, n_particle))

        ! call ARH orbital updating deconstructor C wrapper
        error_c = arh_deconstructor_c_wrapper(dm_ao_3d_c)

        ! check if 3D density matrix is as expected
        if (any(abs(dm_ao_3d_c - 1.0_c_rp) > tol_c)) then
            test_arh_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_arh_deconstructor_c_wrapper failed: Returned AO "// &
                "density matrix for open-shell case wrong."
        end if

        ! deallocate arrays
        deallocate(dm_ao_3d_c)

    end function test_arh_deconstructor_c_wrapper

    logical(c_bool) function test_assign_arh_f_c() bind(C)
        !
        ! this function tests that the function that converts ARH settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh, only: arh_settings_type
        use otr_arh_test_reference, only: assignment(=), ref_arh_settings, operator(/=)
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
        use otr_arh_test_reference, only: ref_arh_settings, assignment(=), operator(/=)

        type(arh_settings_type)   :: settings
        type(arh_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_arh_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_arh_settings

        ! convert to C settings
        settings_c = settings

        ! check that callback function pointers are not associated
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
