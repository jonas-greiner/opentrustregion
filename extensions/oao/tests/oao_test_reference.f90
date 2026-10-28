! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_test_reference

    use opentrustregion_core, only: ip, kw_len
    use c_interface, only: c_ip, c_rp
    use otr_oao, only: oao_settings_type, oao_settings_uninitialized_warning_msg
    use otr_common_test_reference, only: n_particle, n_ao, ref_orbital_settings
    use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_intptr_t, c_null_char, &
                                           c_null_funptr, c_null_ptr, c_associated

    implicit none

    ! number of parameters of the OAO parameterization
    integer(ip), parameter :: n_param = n_particle * n_ao * (n_ao - 1) / 2

    ! reference OAO settings, those shared by the orbital bases
    type(oao_settings_type), parameter :: ref_oao_settings = &
        oao_settings_type(orbital_settings_type=ref_orbital_settings)

    ! warning the C factory prints for OAO settings that were not initialized, as a
    ! null-terminated string for the C integration tests
    integer(ip), parameter :: settings_warning_len = &
        len(oao_settings_uninitialized_warning_msg) + 1
    character(kind=c_char), protected, &
        bind(C, name="test_oao_settings_uninitialized_warning") :: &
        settings_uninitialized_warning_c(settings_warning_len) = &
        transfer(oao_settings_uninitialized_warning_msg//c_null_char, "a", &
                 settings_warning_len)

    interface operator(==)
        module procedure equal_oao
        module procedure equal_oao_c
    end interface

    interface operator(/=)
        module procedure not_equal_oao
        module procedure not_equal_oao_c
    end interface

contains

    function density_pattern(n_rows, n_channels, offset) result(pattern)
        !
        ! this function returns density matrices of the given numbers of AOs and
        ! particle channels whose values encode their particle channel and the pair of
        ! their AO indices regardless of their order (all counted from one) together
        ! with an offset, so that they are symmetric, as the Python interface expects
        ! of density matrices it passes without transposition, while their values
        ! reveal matrices of the wrong shape or channel
        !
        integer(c_ip), intent(in) :: n_rows, n_channels
        real(c_rp), intent(in) :: offset
        real(c_rp) :: pattern(n_rows, n_rows, n_channels)

        integer(c_ip) :: i, j, k

        do k = 1, n_channels
            do j = 1, n_rows
                do i = 1, n_rows
                    pattern(i, j, k) = &
                        offset + real(100 * k + 10 * max(i, j) + min(i, j), kind=c_rp)
                end do
            end do
        end do

    end function density_pattern

    subroutine get_reference_oao_values(values_out) bind(C)
        !
        ! this subroutine exports the reference OAO settings as C settings whose fields
        ! are set one by one without the conversion routines under test, the logging
        ! function and the host context are set to distinct addresses
        !
        use otr_oao_c_interface, only: oao_settings_type_c

        type(oao_settings_type_c), intent(out) :: values_out

        ! set logging function and host context to distinct addresses
        values_out%logger = transfer(1_c_intptr_t, c_null_funptr)
        values_out%context = transfer(2_c_intptr_t, c_null_ptr)

        ! set logical and integer
        values_out%initialized = .true.
        values_out%verbose = int(ref_oao_settings%verbose, kind=c_ip)

    end subroutine get_reference_oao_values

    subroutine reference_oao_field(name_c, field_value, keyword_c) bind(C)
        !
        ! this subroutine returns the value of a field of the C OAO settings exported by
        ! get_reference_oao_values by its name, a numeric, logical (as 1 or 0) or
        ! pointer (as its address) field in field_value and a keyword field as a
        ! null-terminated character array in keyword_c, so that the Python tests need no
        ! values of their own
        !
        use otr_oao_c_interface, only: oao_settings_type_c
        use test_reference, only: ref_character_from_c

        character(kind=c_char), intent(in) :: name_c(*)
        real(c_rp), intent(out) :: field_value
        character(kind=c_char), intent(out) :: keyword_c(kw_len + 1)

        type(oao_settings_type_c) :: values

        ! export reference values and select field by name
        field_value = 0.0_c_rp
        keyword_c = c_null_char
        call get_reference_oao_values(values)
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

    end subroutine reference_oao_field

    logical(c_bool) function is_default_oao_settings(settings_c) bind(C)
        !
        ! this function checks whether C OAO settings hold the default values without
        ! logging function and host context, so that a C test can compare against them
        ! without values of its own
        !
        use otr_oao, only: default_oao_settings
        use otr_oao_c_interface, only: oao_settings_type_c

        type(oao_settings_type_c), intent(in) :: settings_c

        is_default_oao_settings = settings_c == default_oao_settings .and. &
                                  .not. c_associated(settings_c%logger) .and. &
                                  .not. c_associated(settings_c%context)

    end function is_default_oao_settings

    logical function equal_oao(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to
        ! different OAO settings
        !
        use otr_common_test_reference, only: operator(==)

        type(oao_settings_type), intent(in) :: lhs, rhs

        equal_oao = lhs%orbital_settings_type == rhs%orbital_settings_type

    end function equal_oao

    logical function not_equal_oao(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare OAO
        ! settings to different OAO settings
        !
        type(oao_settings_type), intent(in) :: lhs, rhs

        not_equal_oao = .not. (lhs == rhs)

    end function not_equal_oao

    logical function equal_oao_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare C OAO settings to
        ! Fortran OAO settings field by field, independently of the conversion routines
        ! under test
        !
        use otr_oao_c_interface, only: oao_settings_type_c

        type(oao_settings_type_c), intent(in) :: lhs_c
        type(oao_settings_type), intent(in) :: rhs

        equal_oao_c = lhs_c%initialized .and. lhs_c%verbose == rhs%verbose

    end function equal_oao_c

    logical function not_equal_oao_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare C OAO
        ! settings to Fortran OAO settings
        !
        use otr_oao_c_interface, only: oao_settings_type_c

        type(oao_settings_type_c), intent(in) :: lhs_c
        type(oao_settings_type), intent(in) :: rhs

        not_equal_oao_c = .not. (lhs_c == rhs)

    end function not_equal_oao_c

end module otr_oao_test_reference
