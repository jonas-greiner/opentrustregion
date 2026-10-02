! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface_unit_tests

    use opentrustregion, only: ip
    use c_interface, only: c_rp, c_ip
    use otr_common_c_interface, only: evaluate_dm_c_type, get_response_c_type
    use, intrinsic :: iso_c_binding, only: c_funptr, c_funloc

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
        ! this function is a test function for the density matrix evaluating C function
        ! for 2D density matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao
        use otr_common_unit_tests, only: record_mock_call

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
        ! this function is a test function for the density matrix evaluating C function
        ! for 3D density matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao, n_particle
        use otr_common_unit_tests, only: record_mock_call

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
        ! this function is a test function for the response C function for 2D density
        ! matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out), target :: response(*)
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2

        response(:flat_len) = evaluate_dm_factors(2) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_get_response_cs

    function mock_get_response_os(dm_ao, response) result(error) bind(C)
        !
        ! this function is a test function for the response C function for 3D density
        ! matrices
        !
        use otr_common_test_reference, only: evaluate_dm_factors, n_ao, n_particle

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out), target :: response(*)
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2 * n_particle

        response(:flat_len) = evaluate_dm_factors(2) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_get_response_os

end module otr_common_c_interface_unit_tests
