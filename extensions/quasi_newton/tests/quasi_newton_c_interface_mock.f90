! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_c_interface_mock

    use opentrustregion, only: stderr
    use c_interface, only: c_rp, c_ip
    use otr_qn_c_interface, only: update_orbs_qn_factory_c_wrapper, &
                                  init_qn_settings_c, &
                                  update_orbs_qn_deconstructor_c_wrapper
    use otr_qn_test_reference, only: ref_qn_settings
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_f_procpointer, &
                                           c_funloc, c_null_char

    implicit none

    logical(c_bool), bind(C) :: test_update_orbs_qn_factory_interface = .true._c_bool, &
        test_update_orbs_qn_deconstructor_interface = .false._c_bool

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_qn_factory_c_wrapper), pointer :: &
        mock_update_orbs_qn_factory_c_wrapper_ptr => &
        mock_update_orbs_qn_factory_c_wrapper
    procedure(init_qn_settings_c), pointer :: mock_init_qn_settings_c_ptr => &
        mock_init_qn_settings_c
    procedure(update_orbs_qn_deconstructor_c_wrapper), pointer :: &
        mock_update_orbs_qn_deconstructor_c_wrapper_ptr => &
        mock_update_orbs_qn_deconstructor_c_wrapper

contains

    function mock_update_orbs_qn_factory_c_wrapper(update_orbs_orig_c_funptr, &
                                                   n_param_c, settings_c, &
                                                   update_orbs_qn_c_funptr) &
        result(error_c) bind(C, name="mock_update_orbs_qn_factory")
        !
        ! this subroutine is a mock routine for the quasi-Newton orbital updating 
        ! factory C wrapper subroutine
        !
        use otr_qn_c_interface, only: qn_settings_type_c
        use c_interface, only: update_orbs_c_type, hess_x_c_type, logger_c_type
        use test_reference, only: test_update_orbs_c_funptr
        use otr_qn_test_reference, only: operator(/=)
        use c_interface_unit_tests, only: mock_update_orbs

        type(c_funptr), intent(in), value :: update_orbs_orig_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(qn_settings_type_c), intent(in), value :: settings_c
        type(c_funptr), intent(out) :: update_orbs_qn_c_funptr
        integer(c_ip) :: error_c

        procedure(logger_c_type), pointer :: logger_funptr
        character(:), allocatable, target :: message

        ! test passed orbital update function
        test_update_orbs_qn_factory_interface = &
            test_update_orbs_qn_factory_interface .and. &
            test_update_orbs_c_funptr(update_orbs_orig_c_funptr, &
                                      "update_orbs_qn_factory_py_interface", &
                                      " by given orbital updating function")

        ! check if passed number of parameters is correct
        if (n_param_c /= 3) then
            write (stderr, *) "test_update_orbs_qn_factory_interface failed: "// &
                "Passed number of parameters wrong."
            test_update_orbs_qn_factory_interface = .false.
        end if

        ! get Fortran pointer to passed logging function and call it
        message = "test" // c_null_char
        call c_f_procpointer(cptr=settings_c%logger, fptr=logger_funptr)
        call logger_funptr(message)

        ! check optional settings against reference values
        if (settings_c /= ref_qn_settings) then
            write(stderr, *) "test_update_orbs_qn_factory_py_interface failed: "// &
                "Passed settings associated with wrong values."
            test_update_orbs_qn_factory_interface = .false.
        end if

        ! set function pointer to mock quasi-Newton orbital updating function
        update_orbs_qn_c_funptr = c_funloc(mock_update_orbs)

        ! set return arguments
        error_c = 0

    end function mock_update_orbs_qn_factory_c_wrapper

    subroutine mock_init_qn_settings_c(settings) &
        bind(C, name="mock_init_qn_settings")
        !
        ! this subroutine is a mock routine for the C quasi-Newton setting 
        ! initialization subroutine
        !
        use otr_qn_c_interface, only: qn_settings_type_c
        use otr_qn_test_reference, only: assignment(=)

        type(qn_settings_type_c), intent(inout) :: settings

        ! set reference values
        settings = ref_qn_settings

    end subroutine mock_init_qn_settings_c

    subroutine mock_update_orbs_qn_deconstructor_c_wrapper() &
        bind(C, name="mock_update_orbs_qn_deconstructor")
        !
        ! this subroutine is a mock routine for the C quasi-Newton deconstructor 
        ! subroutine
        !
        test_update_orbs_qn_deconstructor_interface = .true._c_bool

    end subroutine mock_update_orbs_qn_deconstructor_c_wrapper

end module otr_qn_c_interface_mock
