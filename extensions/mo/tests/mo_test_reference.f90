! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_test_reference

    use opentrustregion, only: ip, kw_len
    use c_interface, only: c_ip, c_rp
    use otr_mo, only: mo_settings_type, mo_settings_uninitialized_warning_msg
    use otr_common_test_reference, only: n_particle, n_occ, ref_orbital_settings
    use, intrinsic :: iso_c_binding, only: c_size_t, c_sizeof, c_bool, c_char, &
                                           c_intptr_t, c_null_char, c_null_funptr, &
                                           c_null_ptr, c_associated

    implicit none

    ! number of MOs for orbitals parameterized in the MO basis, fewer than the shared
    ! number of AOs, and the number of parameters for the shared occupations
    integer(ip), parameter :: n_mo = 4_ip
    integer(ip), parameter :: n_param = sum(n_occ * (n_mo - n_occ))

    ! dimensions of the MO coefficients passed through the Python interface
    integer(c_ip), protected, bind(C, name="test_n_mo") :: n_mo_c = int(n_mo, kind=c_ip)
    integer(c_ip), protected, bind(C, name="test_n_occ") :: n_occ_c(n_particle) = &
        int(n_occ, kind=c_ip)

    ! occupation cases for routines acting on every particle channel: closed-shell,
    ! open-shell, open-shell with an empty occupied or virtual block, and closed- and
    ! open-shell with MOs of different irreps, chosen so that wrongly sliced irreps,
    ! packing in the wrong order or swapped channels give different parameters, whose
    ! number differs between the open-shell channels
    integer(ip), parameter :: n_cases = 6_ip
    integer(ip), parameter :: case_n_particle(n_cases) = &
        [1_ip, n_particle, n_particle, n_particle, 1_ip, n_particle]
    integer(ip), parameter :: case_n_occ(n_particle, n_cases) = &
        reshape([n_occ(1), 0_ip, n_occ(1), n_occ(2), n_occ(1), 0_ip, n_mo, n_occ(2), &
                 n_occ(1), 0_ip, n_occ(1), n_occ(2)], [n_particle, n_cases])
    integer(ip), parameter :: case_irreps(n_mo, n_particle, n_cases) = reshape( &
        [0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, &
         0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, &
         0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 0_ip, 1_ip, 1_ip, 1_ip, 0_ip, 0_ip, 0_ip, &
         0_ip, 1_ip, 0_ip, 0_ip, 1_ip, 0_ip, 0_ip, 1_ip, 1_ip], &
        [n_mo, n_particle, n_cases])
    character(len=33), parameter :: case_names(n_cases) = &
        [character(len=33) :: "closed-shell", "open-shell", &
         "open-shell empty occupied channel", "open-shell empty virtual channel", &
         "closed-shell with symmetry", "open-shell with symmetry"]

    ! irreps of the MOs passed through the Python interface, those of the last
    ! occupation case, open-shell with symmetry, whose first channel serves the
    ! closed-shell case
    integer(c_ip), protected, bind(C, name="test_orbsym") :: &
        orbsym_c(n_mo, n_particle) = int(case_irreps(:, :, n_cases), kind=c_ip)

    ! reference MO settings, those shared by the orbital bases
    type(mo_settings_type), parameter :: ref_mo_settings = &
        mo_settings_type(orbital_settings_type=ref_orbital_settings)

    ! warning the C factory prints for MO settings that were not initialized, as a
    ! null-terminated string for the C integration tests
    integer(ip), parameter :: settings_warning_len = &
        len(mo_settings_uninitialized_warning_msg) + 1
    character(kind=c_char), protected, &
        bind(C, name="test_mo_settings_uninitialized_warning") :: &
        settings_uninitialized_warning_c(settings_warning_len) = &
        transfer(mo_settings_uninitialized_warning_msg//c_null_char, "a", &
                 settings_warning_len)

    interface operator(==)
        module procedure equal_mo
        module procedure equal_mo_c
    end interface

    interface operator(/=)
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

    subroutine get_reference_mo_values(values_out) bind(C)
        !
        ! this subroutine exports the reference MO settings as C settings whose fields
        ! are set one by one without the conversion routines under test, the logging
        ! function and the host context are set to distinct addresses
        !
        use otr_mo_c_interface, only: mo_settings_type_c

        type(mo_settings_type_c), intent(out) :: values_out

        ! set logging function and host context to distinct addresses
        values_out%logger = transfer(1_c_intptr_t, c_null_funptr)
        values_out%context = transfer(2_c_intptr_t, c_null_ptr)

        ! set logical and integer
        values_out%initialized = ref_mo_settings%initialized
        values_out%verbose = int(ref_mo_settings%verbose, kind=c_ip)

    end subroutine get_reference_mo_values

    subroutine reference_mo_field(name_c, field_value, keyword_c) bind(C)
        !
        ! this subroutine returns the value of a field of the C MO settings exported by
        ! get_reference_mo_values by its name, a numeric, logical (as 1 or 0) or
        ! pointer (as its address) field in field_value and a keyword field as a
        ! null-terminated character array in keyword_c, so that the Python tests need no
        ! values of their own
        !
        use otr_mo_c_interface, only: mo_settings_type_c
        use test_reference, only: ref_character_from_c

        character(kind=c_char), intent(in) :: name_c(*)
        real(c_rp), intent(out) :: field_value
        character(kind=c_char), intent(out) :: keyword_c(kw_len + 1)

        type(mo_settings_type_c) :: values

        ! export reference values and select field by name
        field_value = 0.0_c_rp
        keyword_c = c_null_char
        call get_reference_mo_values(values)
        select case (ref_character_from_c(name_c))
        case ("logger")
            field_value = real(transfer(values%logger, 0_c_intptr_t), kind=c_rp)
        case ("context")
            field_value = real(transfer(values%context, 0_c_intptr_t), kind=c_rp)
        case ("initialized")
            field_value = merge(1.0_c_rp, 0.0_c_rp, values%initialized)
        case ("verbose")
            field_value = values%verbose
        end select

    end subroutine reference_mo_field

    logical(c_bool) function is_default_mo_settings(settings_c) bind(C)
        !
        ! this function checks whether C MO settings hold the default values without
        ! logging function and host context, so that a C test can compare against them
        ! without values of its own
        !
        use otr_mo, only: default_mo_settings
        use otr_mo_c_interface, only: mo_settings_type_c

        type(mo_settings_type_c), intent(in) :: settings_c

        is_default_mo_settings = settings_c == default_mo_settings .and. &
                                 .not. c_associated(settings_c%logger) .and. &
                                 .not. c_associated(settings_c%context)

    end function is_default_mo_settings

    integer(c_size_t) function mo_settings_size() bind(C)
        !
        ! this function returns the size of the C MO settings, so that the layout
        ! tests detect a field missing from the C or Python settings
        !
        use otr_mo_c_interface, only: mo_settings_type_c

        type(mo_settings_type_c) :: settings_c

        mo_settings_size = c_sizeof(settings_c)

    end function mo_settings_size

    logical function equal_mo(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare MO settings to
        ! different MO settings
        !
        use otr_common_test_reference, only: operator(==)

        type(mo_settings_type), intent(in) :: lhs, rhs

        equal_mo = lhs%orbital_settings_type == rhs%orbital_settings_type

    end function equal_mo

    logical function not_equal_mo(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare MO
        ! settings to different MO settings
        !
        type(mo_settings_type), intent(in) :: lhs, rhs

        not_equal_mo = .not. (lhs == rhs)

    end function not_equal_mo

    logical function equal_mo_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare C MO settings to
        ! Fortran MO settings field by field, independently of the conversion routines
        ! under test
        !
        use otr_mo_c_interface, only: mo_settings_type_c

        type(mo_settings_type_c), intent(in) :: lhs_c
        type(mo_settings_type), intent(in) :: rhs

        equal_mo_c = &
            (lhs_c%initialized .eqv. logical(rhs%initialized, kind=c_bool)) .and. &
            lhs_c%verbose == rhs%verbose

    end function equal_mo_c

    logical function not_equal_mo_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare C MO
        ! settings to Fortran MO settings
        !
        use otr_mo_c_interface, only: mo_settings_type_c

        type(mo_settings_type_c), intent(in) :: lhs_c
        type(mo_settings_type), intent(in) :: rhs

        not_equal_mo_c = .not. (lhs_c == rhs)

    end function not_equal_mo_c

end module otr_mo_test_reference
