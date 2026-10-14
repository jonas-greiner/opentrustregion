! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_test_reference

    use opentrustregion, only: ip, kw_len
    use c_interface, only: c_ip, c_rp
    use otr_oao, only: oao_settings_type
    use otr_common_test_reference, only: n_particle, n_ao, ref_orbital_settings
    use, intrinsic :: iso_c_binding, only: c_size_t, c_sizeof, c_bool, c_char, &
                                           c_intptr_t, c_null_char, c_null_funptr, &
                                           c_null_ptr, c_associated

    implicit none

    ! number of parameters of the OAO parameterization
    integer(ip), parameter :: n_param = n_particle * n_ao * (n_ao - 1) / 2

    ! reference OAO settings, those shared by the orbital bases
    type(oao_settings_type), parameter :: ref_oao_settings = &
        oao_settings_type(orbital_settings_type=ref_orbital_settings)

    interface operator(==)
        module procedure equal_oao
        module procedure equal_oao_c
    end interface

    interface operator(/=)
        module procedure not_equal_oao
        module procedure not_equal_oao_c
    end interface

contains

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
        values_out%initialized = ref_oao_settings%initialized
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

    integer(c_size_t) function oao_settings_size() bind(C)
        !
        ! this function returns the size of the C OAO settings, so that the layout
        ! tests detect a field missing from the C or Python settings
        !
        use otr_oao_c_interface, only: oao_settings_type_c

        type(oao_settings_type_c) :: settings_c

        oao_settings_size = c_sizeof(settings_c)

    end function oao_settings_size

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

        equal_oao_c = &
            (lhs_c%initialized .eqv. logical(rhs%initialized, kind=c_bool)) .and. &
            lhs_c%verbose == rhs%verbose

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
