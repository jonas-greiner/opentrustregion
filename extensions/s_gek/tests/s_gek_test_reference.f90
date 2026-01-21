! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_test_reference

    use opentrustregion, only: ip
    use c_interface, only: c_ip
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

    ! derived types for S-GEK settings
    type ref_s_gek_settings_type
        logical :: use_subspace
        integer(ip) :: verbose, max_points
    end type

    type, bind(C) :: ref_s_gek_settings_type_c
        logical(c_bool) :: use_subspace
        integer(c_ip) :: verbose, max_points
    end type

    ! general reference parameters
    type(ref_s_gek_settings_type) :: ref_s_gek_settings = &
        ref_s_gek_settings_type(use_subspace = .false., verbose = 3, max_points = 10)

    interface assignment(=)
        module procedure assign_ref_to_s_gek
        module procedure assign_ref_to_s_gek_c
        module procedure assign_ref_to_ref_c
    end interface

    interface operator(==)
        module procedure equal_s_gek_to_ref
        module procedure equal_s_gek_c_to_ref
        module procedure equal_s_gek
        module procedure equal_s_gek_c
    end interface

    interface operator(/=)
        module procedure not_equal_s_gek_to_ref
        module procedure not_equal_s_gek_c_to_ref
        module procedure not_equal_s_gek
        module procedure not_equal_s_gek_c
    end interface

contains

    subroutine get_reference_s_gek_values(ref_settings_out) bind(C)
        !
        ! this subroutine exports the S-GEK reference values for tests
        !
        type(ref_s_gek_settings_type_c), intent(out) :: ref_settings_out

        ref_settings_out = ref_s_gek_settings

    end subroutine get_reference_s_gek_values

    subroutine assign_ref_to_s_gek(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to set S-GEK settings to 
        ! reference values
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(out) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs

        ! unassociate function pointers
        lhs%logger => null()

        ! set reference values
        lhs%use_subspace = rhs%use_subspace
        lhs%verbose = rhs%verbose
        lhs%max_points = rhs%max_points

        ! set initialization logical
        lhs%initialized = .true.

    end subroutine assign_ref_to_s_gek

    subroutine assign_ref_to_s_gek_c(lhs_c, rhs)
        !
        ! this subroutine overloads the assignment operator to set C S-GEK settings to 
        ! reference values
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(out) :: lhs_c
        type(ref_s_gek_settings_type), intent(in) :: rhs

        type(s_gek_settings_type) :: lhs

        lhs = rhs
        lhs_c = lhs

    end subroutine assign_ref_to_s_gek_c

    subroutine assign_ref_to_ref_c(lhs, rhs)
        !
        ! this subroutine overloads the assignment operator to convert reference values
        ! to their C counterpart
        !
        use c_interface, only: character_to_c

        type(ref_s_gek_settings_type_c), intent(out) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs

        lhs%use_subspace = logical(rhs%use_subspace, kind=c_bool)
        lhs%verbose = int(rhs%verbose, kind=c_ip)
        lhs%max_points = int(rhs%max_points, kind=c_ip)

    end subroutine assign_ref_to_ref_c

    logical function equal_s_gek_to_ref(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to 
        ! reference values
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(in) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs

        equal_s_gek_to_ref = (lhs%use_subspace .eqv. rhs%use_subspace) .and. &
                             lhs%verbose == rhs%verbose .and. lhs%max_points == &
                             rhs%max_points

    end function equal_s_gek_to_ref

    logical function not_equal_s_gek_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK 
        ! settings to reference values
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(in) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs

        not_equal_s_gek_to_ref = .not. (lhs == rhs)

    end function not_equal_s_gek_to_ref

    logical function equal_s_gek_c_to_ref(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to
        ! reference values
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(ref_s_gek_settings_type), intent(in) :: rhs

        type(s_gek_settings_type) :: lhs

        lhs = lhs_c
        equal_s_gek_c_to_ref = lhs == rhs

    end function equal_s_gek_c_to_ref

    logical function not_equal_s_gek_c_to_ref(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK
        ! settings to reference values
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c

        type(s_gek_settings_type_c), intent(in) :: lhs
        type(ref_s_gek_settings_type), intent(in) :: rhs
        
        not_equal_s_gek_c_to_ref = .not. (lhs == rhs)

    end function not_equal_s_gek_c_to_ref

    logical function equal_s_gek(lhs, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to 
        ! different S-GEK settings
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(in) :: lhs, rhs
        
        equal_s_gek = (lhs%use_subspace .eqv. rhs%use_subspace) .and. lhs%verbose == &
                      rhs%verbose .and. lhs%max_points == rhs%max_points

    end function equal_s_gek

    logical function not_equal_s_gek(lhs, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK
        ! settings to different S-GEK settings
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type), intent(in) :: lhs, rhs
        
        not_equal_s_gek = .not. (lhs == rhs)

    end function not_equal_s_gek

    logical function equal_s_gek_c(lhs_c, rhs)
        !
        ! this function overloads the comparison operator to compare S-GEK settings to
        ! settings to different S-GEK settings
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c, assignment(=)
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(s_gek_settings_type), intent(in) :: rhs
        
        type(s_gek_settings_type) :: lhs

        lhs = lhs_c
        equal_s_gek_c = lhs == rhs

    end function equal_s_gek_c

    logical function not_equal_s_gek_c(lhs_c, rhs)
        !
        ! this function overloads the negated comparison operator to compare S-GEK
        ! settings to different S-GEK settings
        !
        use otr_s_gek_c_interface, only: s_gek_settings_type_c
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(in) :: lhs_c
        type(s_gek_settings_type), intent(in) :: rhs
        
        not_equal_s_gek_c = .not. (lhs_c == rhs)

    end function not_equal_s_gek_c

end module otr_s_gek_test_reference
