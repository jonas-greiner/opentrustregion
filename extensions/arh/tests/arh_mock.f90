! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_mock

    use opentrustregion, only: rp, ip, stderr, obj_func_type, precond_type
    use otr_arh, only: arh_factory, arh_deconstructor
    use test_reference, only: tol
    use otr_arh_test_reference, only: ref_arh_settings, operator(/=)

    implicit none

    logical :: test_passed

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_factory), pointer :: mock_arh_factory_ptr => mock_arh_factory
    procedure(arh_deconstructor), pointer :: mock_arh_deconstructor_ptr => &
        mock_arh_deconstructor
    procedure(obj_func_type), pointer :: mock_obj_func_arh_ptr => mock_obj_func_arh
    procedure(precond_type), pointer ::  mock_precond_arh_ptr => mock_precond_arh

contains

    subroutine mock_arh_factory(dm_ao, ao_overlap, n_ao, get_energy_funptr, &
                                get_fock_funptr, obj_func_arh_funptr, &
                                update_orbs_arh_funptr, precond_arh_funptr, error, &
                                settings)
        !
        ! this function is a test function for the function which returns a modified
        ! orbital updating function
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use otr_arh, only: get_energy_type, get_fock_type, arh_settings_type
        use otr_arh_test_reference, only: test_get_energy_funptr, test_get_fock_funptr
        use otr_common_mock, only: mock_update_orbs

        real(rp), intent(in) :: dm_ao(:, :), ao_overlap(:, :)
        integer(ip), intent(in) :: n_ao
        procedure(get_energy_type), intent(in), pointer :: get_energy_funptr
        procedure(get_fock_type), intent(in), pointer :: get_fock_funptr
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        procedure(precond_type), intent(out), pointer :: precond_arh_funptr
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(in) :: settings

        ! initialize logical
        test_passed = .true.

        ! check passed arrays
        if (any(abs(dm_ao - 1.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed AO "// &
                "density matrix wrong."
        end if
        if (any(abs(ao_overlap - 2.0_rp) > tol)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed AO "// &
                "overlap matrix wrong."
        end if

        ! check number of AOs
        if (n_ao /= 3) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed number of "// &
                "AOs wrong."
        end if

        ! test passed energy function
        test_passed = test_passed .and. &
            test_get_energy_funptr(get_energy_funptr, "arh_factory_c_wrapper", &
                                   " by passed energy function")

        ! test passed Fock matrix function
        test_passed = test_passed .and. &
            test_get_fock_funptr(get_fock_funptr, "arh_factory_c_wrapper", &
                                 " by passed Fock matrix function")

        ! check if optional logging function is correctly passed
        if (.not. associated(settings%logger)) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: Passed logging "// &
                "function not associated with value."
        else
            call settings%logger("test")
        end if

        ! check if optional settings are correctly passed
        if (settings /= ref_arh_settings) then
            test_passed = .false.
            write (stderr, *) "test_arh_factory_c_wrapper failed: "// &
                "Passed optional settings associated with wrong values."
        end if

        ! set output quantities
        error = 0
        obj_func_arh_funptr => mock_obj_func_arh
        update_orbs_arh_funptr => mock_update_orbs
        precond_arh_funptr => mock_precond_arh

    end subroutine mock_arh_factory

    subroutine mock_arh_deconstructor(dm_ao, error)
        !
        ! this subroutine is a test function for the ARH deconstructor
        !
        real(rp), intent(out) :: dm_ao(:, :)
        integer(ip), intent(out) :: error

        dm_ao = 1.0_rp

        error = 0

    end subroutine mock_arh_deconstructor

    function mock_obj_func_arh(kappa, error) result(func)
        !
        ! this function is a test function for the ARH C objective function
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: func

        func = sum(kappa)

        error = 0

    end function mock_obj_func_arh

    subroutine mock_precond_arh(residual, mu, precond_residual, error)
        !
        ! this function is a test function for the ARH C preconditioner function
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        precond_residual = mu * residual

        error = 0

    end subroutine mock_precond_arh

end module otr_arh_mock
