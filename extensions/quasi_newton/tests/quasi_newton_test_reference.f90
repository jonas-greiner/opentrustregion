! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_test_reference

    use opentrustregion, only: ip, kw_len
    use c_interface, only: c_ip
    use, intrinsic :: iso_c_binding, only: c_char

    implicit none

    ! derived types for quasi-Newton settings
    type ref_qn_settings_type
        integer(ip) :: verbose
        character(kw_len, c_char) :: hess_update_scheme
    end type

    type, bind(C) :: ref_qn_settings_type_c
        integer(c_ip) :: verbose
        character(c_char) :: hess_update_scheme(kw_len + 1)
    end type

    ! general reference parameters
    type(ref_qn_settings_type) :: ref_qn_settings = &
        ref_qn_settings_type(verbose = 3, hess_update_scheme = "bfgs")

    interface assignment(=)
        module procedure assign_ref_to_qn
        module procedure assign_ref_to_qn_c
        module procedure assign_ref_to_ref_c
    end interface

    interface operator(==)
        module procedure equal_qn_to_ref
        module procedure equal_qn_c_to_ref
        module procedure equal_qn
        module procedure equal_qn_c
    end interface

    interface operator(/=)
        module procedure not_equal_qn_to_ref
        module procedure not_equal_qn_c_to_ref
        module procedure not_equal_qn
        module procedure not_equal_qn_c
    end interface

contains

    subroutine get_reference_qn_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the quasi-Newton reference values for tests
        !
        type(ref_qn_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_qn_settings

    end subroutine get_reference_qn_values

    subroutine assign_ref_to_qn(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set quasi-Newton 
        ! settings to reference values
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(out) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%logger => null()

        ! set reference values
        lhs%verbose = rhs%verbose
        lhs%hess_update_scheme = rhs%hess_update_scheme

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_qn

    subroutine assign_ref_to_qn_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C quasi-Newton 
        ! settings to reference values
        !
        use otr_qn_c_interface, only: qn_settings_type_c, assignment(=)
        use otr_qn, only: qn_settings_type

        type(qn_settings_type_c), intent(out) :: lhs_c
        type(ref_qn_settings_type), intent(in) :: rhs

        type(qn_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_qn_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        use c_interface, only: character_to_c

        type(ref_qn_settings_type_c), intent(out) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs

        lhs%verbose = int(rhs%verbose, kind=c_ip)
        lhs%hess_update_scheme = character_to_c(rhs%hess_update_scheme)

    end subroutine assign_ref_to_ref_c

    logical function equal_qn_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare quasi-Newton 
        ! settings to reference values
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(in) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs

        equal_qn_to_ref = lhs%verbose == rhs%verbose .and. lhs%hess_update_scheme == &
            rhs%hess_update_scheme

    end function equal_qn_to_ref

    logical function not_equal_qn_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare 
        ! quasi-Newton settings to reference values
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(in) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs

        not_equal_qn_to_ref = .not. (lhs == rhs)

    end function not_equal_qn_to_ref

    logical function equal_qn_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare quasi-Newton 
        ! settings to reference values
        !
        use otr_qn_c_interface, only: qn_settings_type_c, assignment(=)
        use otr_qn, only: qn_settings_type

        type(qn_settings_type_c), intent(in) :: lhs_c
        type(ref_qn_settings_type), intent(in) :: rhs

        type(qn_settings_type) :: lhs

        lhs = lhs_c
        equal_qn_c_to_ref = lhs == rhs

    end function equal_qn_c_to_ref

    logical function not_equal_qn_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare 
        ! quasi-Newton settings to reference values
        !
        use otr_qn_c_interface, only: qn_settings_type_c

        type(qn_settings_type_c), intent(in) :: lhs
        type(ref_qn_settings_type), intent(in) :: rhs
        
        not_equal_qn_c_to_ref = .not. (lhs == rhs)

    end function not_equal_qn_c_to_ref

    logical function equal_qn(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare quasi-Newton 
        ! settings to different quasi-Newton settings
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(in) :: lhs, rhs
        
        equal_qn = lhs%verbose == rhs%verbose .and. lhs%hess_update_scheme == &
            rhs%hess_update_scheme

    end function equal_qn

    logical function not_equal_qn(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare 
        ! quasi-Newton settings to different quasi-Newton settings
        !
        use otr_qn, only: qn_settings_type

        type(qn_settings_type), intent(in) :: lhs, rhs
        
        not_equal_qn = .not. (lhs == rhs)

    end function not_equal_qn

    logical function equal_qn_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare quasi-Newton 
        ! settings to different quasi-Newton settings
        !
        use otr_qn_c_interface, only: qn_settings_type_c, assignment(=)
        use otr_qn, only: qn_settings_type

        type(qn_settings_type_c), intent(in) :: lhs_c
        type(qn_settings_type), intent(in) :: rhs
        
        type(qn_settings_type) :: lhs

        lhs = lhs_c
        equal_qn_c = lhs == rhs

    end function equal_qn_c

    logical function not_equal_qn_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare 
        ! quasi-Newton settings to different quasi-Newton settings
        !
        use otr_qn_c_interface, only: qn_settings_type_c
        use otr_qn, only: qn_settings_type

        type(qn_settings_type_c), intent(in) :: lhs_c
        type(qn_settings_type), intent(in) :: rhs
        
        not_equal_qn_c = .not. (lhs_c == rhs)

    end function not_equal_qn_c

end module otr_qn_test_reference
