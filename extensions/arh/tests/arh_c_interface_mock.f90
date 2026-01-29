! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip
    use otr_arh_c_interface, only: arh_factory_c_wrapper, init_arh_settings_c, &
                                   arh_deconstructor_c_wrapper
    use otr_arh_test_reference, only: ref_arh_settings, n_particle, n_ao
    use test_reference, only: n_param
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char

    implicit none

    logical(c_bool), bind(C) :: test_arh_factory_interface = .true._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory_c_wrapper), pointer :: mock_arh_factory_c_wrapper_ptr => &
        mock_arh_factory_c_wrapper
    procedure(init_arh_settings_c), pointer :: mock_init_arh_settings_c_ptr => &
        mock_init_arh_settings_c
    procedure(arh_deconstructor_c_wrapper), pointer :: &
        mock_arh_deconstructor_c_wrapper_2d_ptr => &
        mock_arh_deconstructor_c_wrapper_2d, &
        mock_arh_deconstructor_c_wrapper_3d_ptr => mock_arh_deconstructor_c_wrapper_3d

contains

    function mock_arh_factory_c_wrapper(dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, &
                                        get_energy_c_funptr, get_fock_c_funptr, &
                                        obj_func_arh_c_funptr, &
                                        update_orbs_arh_c_funptr, &
                                        precond_arh_c_funptr, settings_c) &
        result(error_c) bind(C, name="mock_arh_factory")
        !
        ! this subroutine is a mock routine for the ARH orbital updating factory C
        ! wrapper subroutine
        !
        use otr_arh_c_interface, only: arh_settings_type_c
        use c_interface, only: obj_func_c_type, update_orbs_c_type, precond_c_type, &
                               logger_c_type
        use test_reference, only: tol_c
        use otr_arh_test_reference, only: test_get_energy_closed_shell_c_funptr, &
                                          test_get_energy_open_shell_c_funptr, &
                                          test_get_fock_c_funptr, &
                                          test_get_fock_jk_c_funptr, operator(/=)
        use c_interface_unit_tests, only: mock_obj_func, mock_update_orbs, mock_precond

        real(c_rp), intent(in), target :: dm_ao_c(*), ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: get_energy_c_funptr, get_fock_c_funptr
        type(c_funptr), intent(out) :: obj_func_arh_c_funptr, &
                                       update_orbs_arh_c_funptr, precond_arh_c_funptr
        type(arh_settings_type_c), intent(in), value :: settings_c
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message
        procedure(obj_func_c_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_c_type), pointer :: update_orbs_arh_funptr
        procedure(precond_c_type), pointer :: precond_arh_funptr
        real(c_rp), allocatable :: test(:,:,:)

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
                test_get_energy_closed_shell_c_funptr(get_energy_c_funptr, &
                                                      "arh_factory_py_interface", &
                                                      " by given energy function "// &
                                                      "for closed-shell case")

            ! test passed Fock matrix function
            test_arh_factory_interface = test_arh_factory_interface .and. &
                test_get_fock_c_funptr(get_fock_c_funptr, "arh_factory_py_interface", &
                                       " by given Fock matrix function")

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

            ! set function pointers to mock to ARH mock functions
            obj_func_arh_c_funptr = c_funloc(mock_obj_func)
            update_orbs_arh_c_funptr = c_funloc(mock_update_orbs)
            precond_arh_c_funptr = c_funloc(mock_precond)

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
                test_get_energy_open_shell_c_funptr(get_energy_c_funptr, &
                                                    "arh_factory_py_interface", &
                                                    " by given energy function for "// &
                                                    "for open-shell case")
            ! test passed Fock matrix function
            test_arh_factory_interface = test_arh_factory_interface .and. &
                test_get_fock_jk_c_funptr(get_fock_c_funptr, &
                                          "arh_factory_py_interface", " by given "// &
                                          "Fock, Coulomb, and exchange matrix function")

        ! number of particles is not correct
        else
            write (stderr, *) "test_arh_factory_py_interface failed: Passed number "// &
                "of particles wrong."
            test_arh_factory_interface = .false.

        end if

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

    function mock_arh_deconstructor_c_wrapper_2d(dm_ao_c) result(error_c) &
        bind(C, name="mock_arh_deconstructor_2d")
        !
        ! this subroutine is a mock routine for the C ARH deconstructor subroutine for 
        ! 2D density matrices
        !
        real(c_rp), intent(out), target :: dm_ao_c(*)

        integer(c_ip) :: error_c

        ! set return arguments
        dm_ao_c(:n_ao ** 2) = 1.0_c_rp
        error_c = 0

    end function mock_arh_deconstructor_c_wrapper_2d

    function mock_arh_deconstructor_c_wrapper_3d(dm_ao_c) result(error_c) &
        bind(C, name="mock_arh_deconstructor_3d")
        !
        ! this subroutine is a mock routine for the C ARH deconstructor subroutine for 
        ! 3D density matrices
        !
        real(c_rp), intent(out), target :: dm_ao_c(*)

        integer(c_ip) :: error_c

        ! set return arguments
        dm_ao_c(:n_ao ** 2 * n_particle) = 1.0_c_rp
        error_c = 0

    end function mock_arh_deconstructor_c_wrapper_3d

end module otr_arh_c_interface_mock
