! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_test_reference

    use opentrustregion, only: ip
    use c_interface, only: c_ip
    use otr_common_test_reference, only: n_particle, n_ao, ref_orbital_settings_type, &
                                         ref_orbital_settings

    implicit none

    ! number of parameters of the OAO parameterization
    integer(ip), parameter :: n_param = n_particle * n_ao * (n_ao - 1) / 2

    ! derived types for OAO settings
    type, extends(ref_orbital_settings_type) :: ref_oao_settings_type
    end type

    type, bind(C) :: ref_oao_settings_type_c
        integer(c_ip) :: verbose
    end type

    ! general reference parameters
    type(ref_oao_settings_type), parameter :: ref_oao_settings = &
        ref_oao_settings_type(ref_orbital_settings_type=ref_orbital_settings)

    interface assignment(=)
        module procedure assign_ref_to_ref_c
        module procedure assign_ref_to_oao
        module procedure assign_ref_to_oao_c
    end interface

    interface operator(==)
        module procedure equal_oao_to_ref
        module procedure equal_oao_c_to_ref
        module procedure equal_oao
        module procedure equal_oao_c
    end interface

    interface operator(/=)
        module procedure not_equal_oao_to_ref
        module procedure not_equal_oao_c_to_ref
        module procedure not_equal_oao
        module procedure not_equal_oao_c
    end interface

contains

    subroutine get_reference_oao_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the OAO reference values for tests
        !
        type(ref_oao_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_oao_settings

    end subroutine get_reference_oao_values

    subroutine get_default_oao_values(default_values_out) bind(C)
        !
        ! this subroutine exports the default values for tests
        !
        use otr_oao, only: default_oao_settings
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)

        type(oao_settings_type_c), intent(out) :: default_values_out

        default_values_out = default_oao_settings

    end subroutine get_default_oao_values

    subroutine assign_ref_to_oao(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set OAO settings to
        ! reference values
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type), intent(out) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%logger => null()

        ! set reference values
        lhs%verbose = rhs%verbose

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_oao

    subroutine assign_ref_to_oao_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C OAO settings to
        ! reference values
        !
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(out) :: lhs_c
        type(ref_oao_settings_type), intent(in) :: rhs

        type(oao_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_oao_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        type(ref_oao_settings_type_c), intent(out) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs

        lhs%verbose = int(rhs%verbose, kind=c_ip)

    end subroutine assign_ref_to_ref_c

    logical function equal_oao_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to
        ! reference values
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type), intent(in) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs

        equal_oao_to_ref = lhs%verbose == rhs%verbose

    end function equal_oao_to_ref

    logical function not_equal_oao_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare OAO
        ! settings to reference values
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type), intent(in) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs

        not_equal_oao_to_ref = .not. (lhs == rhs)

    end function not_equal_oao_to_ref

    logical function equal_oao_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to
        ! reference values
        !
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(in) :: lhs_c
        type(ref_oao_settings_type), intent(in) :: rhs

        type(oao_settings_type) :: lhs

        lhs = lhs_c
        equal_oao_c_to_ref = lhs == rhs

    end function equal_oao_c_to_ref

    logical function not_equal_oao_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare OAO
        ! settings to reference values
        !
        use otr_oao_c_interface, only: oao_settings_type_c

        type(oao_settings_type_c), intent(in) :: lhs
        type(ref_oao_settings_type), intent(in) :: rhs

        not_equal_oao_c_to_ref = .not. (lhs == rhs)

    end function not_equal_oao_c_to_ref

    logical function equal_oao(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to
        ! different OAO settings
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type), intent(in) :: lhs, rhs

        equal_oao = lhs%verbose == rhs%verbose

    end function equal_oao

    logical function not_equal_oao(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare OAO
        ! settings to different OAO settings
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type), intent(in) :: lhs, rhs

        not_equal_oao = .not. (lhs == rhs)

    end function not_equal_oao

    logical function equal_oao_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare OAO settings to
        ! different OAO settings
        !
        use otr_oao_c_interface, only: oao_settings_type_c, assignment(=)
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(in) :: lhs_c
        type(oao_settings_type), intent(in) :: rhs

        type(oao_settings_type) :: lhs

        lhs = lhs_c
        equal_oao_c = lhs == rhs

    end function equal_oao_c

    logical function not_equal_oao_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare OAO
        ! settings to different OAO settings
        !
        use otr_oao_c_interface, only: oao_settings_type_c
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(in) :: lhs_c
        type(oao_settings_type), intent(in) :: rhs

        not_equal_oao_c = .not. (lhs_c == rhs)

    end function not_equal_oao_c

end module otr_oao_test_reference
