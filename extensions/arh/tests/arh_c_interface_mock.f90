! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip
    use otr_arh_test_reference, only: ref_arh_settings
    use otr_arh_c_interface, only: arh_factory_c_wrapper, init_arh_settings_c, &
                                   arh_deconstructor_c_wrapper
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char, c_f_pointer, c_loc

    implicit none

    logical(c_bool), bind(C) :: test_arh_factory_interface = .true._c_bool, &
                                test_arh_deconstructor_interface = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory_c_wrapper), pointer :: mock_arh_factory_c_wrapper_ptr => &
        mock_arh_factory_c_wrapper
    procedure(init_arh_settings_c), pointer :: mock_init_arh_settings_c_ptr => &
        mock_init_arh_settings_c
     procedure(arh_deconstructor_c_wrapper), pointer :: &
        mock_arh_deconstructor_c_wrapper_ptr => mock_arh_deconstructor_c_wrapper

contains

    function mock_arh_factory_c_wrapper(dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, &
                                        get_energy_c_funptr, update_dm_c_funptr, &
                                        obj_func_arh_c_funptr, &
                                        update_orbs_arh_c_funptr, &
                                        project_arh_c_funptr, settings_c) &
        result(error_c) bind(C, name="mock_arh_factory")
        !
        ! this subroutine is a mock routine for the ARH orbital updating factory C
        ! wrapper subroutine
        !
        use otr_arh_c_interface, only: arh_settings_type_c
        use c_interface, only: obj_func_c_type, update_orbs_c_type, project_c_type, &
                               logger_c_type
        use otr_oao_c_interface, only: dm_ao_3d_c
        use test_reference, only: tol_c
        use otr_oao_test_reference, only: test_get_energy_cs_c_funptr, &
                                          test_get_energy_os_c_funptr
        use otr_arh_test_reference, only: test_update_dm_os_c_funptr, &
                                          test_update_dm_cs_c_funptr, &
                                          operator(/=)
        use c_interface_unit_tests, only: mock_obj_func, mock_project
        use otr_oao_c_interface_mock, only: mock_update_orbs_oao

        real(c_rp), intent(in), target :: dm_ao_c(*), ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: get_energy_c_funptr, update_dm_c_funptr
        type(c_funptr), intent(out) :: obj_func_arh_c_funptr, &
                                       update_orbs_arh_c_funptr, project_arh_c_funptr
        type(arh_settings_type_c), intent(inout) :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message
        procedure(obj_func_c_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_c_type), pointer :: update_orbs_arh_funptr
        procedure(project_c_type), pointer :: project_arh_funptr

        ! set global pointer to density matrix so that it can be accessed in the mock 
        ! orbital updating function
        call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_3d_c, [n_ao_c, n_ao_c, n_particle_c])

        ! closed-shell case
        if (n_particle_c == 1) then
            ! check passed arrays
            if (any(abs(dm_ao_c(:n_ao_c ** 2) - 1.0_c_rp) > tol_c)) then
                write(stderr, *) "test_arh_factory_py_interface failed: Passed AO "// &
                    "density matrix for closed-shell case wrong."
                test_arh_factory_interface = .false.
            end if
            if (any(abs(ao_overlap_c(:n_ao_c ** 2) - 2.0_c_rp) > tol_c)) then
                write(stderr, *) "test_arh_factory_py_interface failed: Passed "// &
                    "AO overlap matrix wrong."
                test_arh_factory_interface = .false.
            end if

            ! test passed energy function
            test_arh_factory_interface = test_arh_factory_interface .and. &
                test_get_energy_cs_c_funptr(get_energy_c_funptr, &
                                            "arh_factory_py_interface", &
                                            " by given energy function")

            ! test passed density matrix updating function
            test_arh_factory_interface = test_arh_factory_interface .and. &
                test_update_dm_cs_c_funptr(update_dm_c_funptr, &
                                           "arh_factory_py_interface", " by given "// &
                                           "density matrix updating function with "// &
                                           "non-linear potential contribution")

            ! check if passed number of AOs is correct
            if (n_ao_c /= 3) then
                write (stderr, *) "test_arh_factory_py_interface failed: Passed "// &
                    "number of AOs wrong."
                test_arh_factory_interface = .false.
            end if

            ! get Fortran pointer to passed logging function and call it
            message = "test" // c_null_char
            call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
            call logger_funptr(message)

            ! check optional settings against reference values
            if (settings_c /= ref_arh_settings) then
                write(stderr, *) "test_arh_factory_py_interface failed: Passed "// &
                    "settings associated with wrong values."
                test_arh_factory_interface = .false.
            end if

        ! open-shell case
        else if (n_particle_c == 2) then
            ! check passed arrays
            if (any(abs(dm_ao_c(:n_ao_c ** 2 * n_particle_c) - 1.0_c_rp) > tol_c)) then
                write(stderr, *) "test_arh_factory_py_interface failed: Passed AO "// &
                    "density matrix for open-shell case wrong."
                test_arh_factory_interface = .false.
            end if

            ! test passed energy function
            test_arh_factory_interface = test_arh_factory_interface .and. &
                test_get_energy_os_c_funptr(get_energy_c_funptr, &
                                            "arh_factory_py_interface", &
                                            " by given energy function")

            ! test passed density matrix updating function
            test_arh_factory_interface = test_arh_factory_interface .and. &
                test_update_dm_os_c_funptr(update_dm_c_funptr, &
                                           "arh_factory_py_interface", " by given "// &
                                           "density matrix updating function with "// &
                                           "separate same- and opposite-spin "// &
                                           "potential contributions")

        ! number of particles is not correct
        else
            write (stderr, *) "test_arh_factory_py_interface failed: Passed number "// &
                "of particles wrong."
            test_arh_factory_interface = .false.

        end if

        ! set function pointers to mock to ARH mock functions
        obj_func_arh_c_funptr = c_funloc(mock_obj_func)
        update_orbs_arh_c_funptr = c_funloc(mock_update_orbs_oao)
        project_arh_c_funptr = c_funloc(mock_project)

        ! set return arguments
        error_c = 0

    end function mock_arh_factory_c_wrapper

    subroutine mock_init_arh_settings_c(settings) &
        bind(C, name="mock_init_arh_settings")
        !
        ! this subroutine is a mock routine for the C ARH setting initialization 
        ! subroutine
        !
        use otr_arh_c_interface, only: arh_settings_type_c
        use otr_arh_test_reference, only: assignment(=)

        type(arh_settings_type_c), intent(inout) :: settings

        ! set reference values
        settings = ref_arh_settings

    end subroutine mock_init_arh_settings_c

    subroutine mock_arh_deconstructor_c_wrapper() bind(C, name="mock_arh_deconstructor")
        !
        ! this subroutine is a mock routine for the C ARH deconstructor subroutine
        !
        test_arh_deconstructor_interface = .true._c_bool

    end subroutine mock_arh_deconstructor_c_wrapper

end module otr_arh_c_interface_mock
