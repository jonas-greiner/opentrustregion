! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_test_reference

    use opentrustregion, only: ip
    use c_interface, only: c_ip, c_rp
    use otr_common_test_reference, only: n_particle, n_occ, ref_orbital_settings_type, &
                                         ref_orbital_settings

    implicit none

    ! number of MOs for orbitals parameterized in the MO basis, fewer than the shared
    ! number of AOs, and the numbers of parameters for the shared occupations
    integer(ip), parameter :: n_mo = 4_ip
    integer(ip), parameter :: n_param_cs = n_occ(1) * (n_mo - n_occ(1)), &
                              n_param_os = sum(n_occ * (n_mo - n_occ))

    ! dimensions of the MO coefficients passed through the Python interface
    integer(c_ip), protected, bind(C, name="test_n_mo") :: n_mo_c = int(n_mo, kind=c_ip)
    integer(c_ip), protected, bind(C, name="test_n_occ") :: n_occ_c(n_particle) = &
        int(n_occ, kind=c_ip)

    ! occupation cases for routines acting on every particle channel: closed-shell,
    ! open-shell, and open-shell with an empty occupied or virtual block
    integer(ip), parameter :: n_cases = 4_ip
    integer(ip), parameter :: case_n_particle(n_cases) = &
        [1_ip, n_particle, n_particle, n_particle]
    integer(ip), parameter :: case_n_occ(n_particle, n_cases) = &
        reshape([n_occ(1), 0_ip, n_occ(1), n_occ(2), n_occ(1), 0_ip, n_mo, n_occ(2)], &
                [n_particle, n_cases])
    character(33), parameter :: case_names(n_cases) = &
        [character(33) :: "closed-shell", "open-shell", &
         "open-shell empty occupied channel", "open-shell empty virtual channel"]

    ! derived types for MO settings
    type, extends(ref_orbital_settings_type) :: ref_mo_settings_type
    end type

    type, bind(C) :: ref_mo_settings_type_c
        integer(c_ip) :: verbose
    end type

    ! general reference parameters
    type(ref_mo_settings_type), parameter :: ref_mo_settings = &
        ref_mo_settings_type(ref_orbital_settings_type = ref_orbital_settings)

    interface assignment(=)
        module procedure assign_ref_to_ref_c
        module procedure assign_ref_to_mo
        module procedure assign_ref_to_mo_c
    end interface

    interface operator(==)
        module procedure equal_mo_to_ref
        module procedure equal_mo_c_to_ref
        module procedure equal_mo
        module procedure equal_mo_c
    end interface

    interface operator(/=)
        module procedure not_equal_mo_to_ref
        module procedure not_equal_mo_c_to_ref
        module procedure not_equal_mo
        module procedure not_equal_mo_c
    end interface

contains

    function mo_coeff_pattern(n_rows, n_cols, n_channels, offset) result(pattern)
        !
        ! this function returns MO coefficients with the given numbers of AOs (rows),
        ! MOs (columns) and particle channels whose values encode their particle
        ! channel, AO and MO index (all counted from one) together with an offset
        !
        integer(c_ip), intent(in) :: n_rows, n_cols, n_channels
        real(c_rp), intent(in) :: offset
        real(c_rp) :: pattern(n_rows, n_cols, n_channels)

        integer(c_ip) :: i, j, k

        do k = 1, n_channels
            do j = 1, n_cols
                do i = 1, n_rows
                    pattern(i, j, k) = offset + real(100 * k + 10 * i + j, kind=c_rp)
                end do
            end do
        end do

    end function mo_coeff_pattern

    subroutine get_reference_mo_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the MO reference values for tests
        !
        type(ref_mo_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_mo_settings

    end subroutine get_reference_mo_values

    subroutine get_default_mo_values(default_values_out) bind(C)
        !
        ! this subroutine exports the default values for tests
        !
        use otr_mo, only: default_mo_settings
        use otr_mo_c_interface, only: mo_settings_type_c, assignment(=)

        type(mo_settings_type_c), intent(out) :: default_values_out

        default_values_out = default_mo_settings

    end subroutine get_default_mo_values

    subroutine assign_ref_to_mo(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set MO settings to
        ! reference values
        !
        use otr_mo, only: mo_settings_type

        type(mo_settings_type), intent(out) :: lhs
        type(ref_mo_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%logger => null()

        ! set reference values
        lhs%verbose = rhs%verbose

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_mo

    subroutine assign_ref_to_mo_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C MO settings to
        ! reference values
        !
        use otr_mo_c_interface, only: mo_settings_type_c, assignment(=)
        use otr_mo, only: mo_settings_type

        type(mo_settings_type_c), intent(out) :: lhs_c
        type(ref_mo_settings_type), intent(in) :: rhs

        type(mo_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_mo_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        type(ref_mo_settings_type_c), intent(out) :: lhs
        type(ref_mo_settings_type), intent(in) :: rhs

        lhs%verbose = int(rhs%verbose, kind=c_ip)

    end subroutine assign_ref_to_ref_c

    logical function equal_mo_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare MO settings to
        ! reference values
        !
        use otr_mo, only: mo_settings_type

        type(mo_settings_type), intent(in) :: lhs
        type(ref_mo_settings_type), intent(in) :: rhs

        equal_mo_to_ref = lhs%verbose == rhs%verbose

    end function equal_mo_to_ref

    logical function not_equal_mo_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare MO
        ! settings to reference values
        !
        use otr_mo, only: mo_settings_type

        type(mo_settings_type), intent(in) :: lhs
        type(ref_mo_settings_type), intent(in) :: rhs

        not_equal_mo_to_ref = .not. (lhs == rhs)

    end function not_equal_mo_to_ref

    logical function equal_mo_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare MO settings to
        ! reference values
        !
        use otr_mo_c_interface, only: mo_settings_type_c, assignment(=)
        use otr_mo, only: mo_settings_type

        type(mo_settings_type_c), intent(in) :: lhs_c
        type(ref_mo_settings_type), intent(in) :: rhs

        type(mo_settings_type) :: lhs

        lhs = lhs_c
        equal_mo_c_to_ref = lhs == rhs

    end function equal_mo_c_to_ref

    logical function not_equal_mo_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare MO
        ! settings to reference values
        !
        use otr_mo_c_interface, only: mo_settings_type_c

        type(mo_settings_type_c), intent(in) :: lhs
        type(ref_mo_settings_type), intent(in) :: rhs

        not_equal_mo_c_to_ref = .not. (lhs == rhs)

    end function not_equal_mo_c_to_ref

    logical function equal_mo(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare MO settings to
        ! different MO settings
        !
        use otr_mo, only: mo_settings_type

        type(mo_settings_type), intent(in) :: lhs, rhs

        equal_mo = lhs%verbose == rhs%verbose

    end function equal_mo

    logical function not_equal_mo(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare MO
        ! settings to different MO settings
        !
        use otr_mo, only: mo_settings_type

        type(mo_settings_type), intent(in) :: lhs, rhs

        not_equal_mo = .not. (lhs == rhs)

    end function not_equal_mo

    logical function equal_mo_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare MO settings to
        ! different MO settings
        !
        use otr_mo_c_interface, only: mo_settings_type_c, assignment(=)
        use otr_mo, only: mo_settings_type

        type(mo_settings_type_c), intent(in) :: lhs_c
        type(mo_settings_type), intent(in) :: rhs

        type(mo_settings_type) :: lhs

        lhs = lhs_c
        equal_mo_c = lhs == rhs

    end function equal_mo_c

    logical function not_equal_mo_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare MO
        ! settings to different MO settings
        !
        use otr_mo_c_interface, only: mo_settings_type_c
        use otr_mo, only: mo_settings_type

        type(mo_settings_type_c), intent(in) :: lhs_c
        type(mo_settings_type), intent(in) :: rhs

        not_equal_mo_c = .not. (lhs_c == rhs)

    end function not_equal_mo_c

end module otr_mo_test_reference
