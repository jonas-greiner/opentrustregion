! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_c_interface

    use opentrustregion, only: ip, rp, kw_len, update_orbs_type, hess_x_type
    use c_interface, only: c_ip, c_rp, update_orbs_c_type, hess_x_c_type
    use otr_common_c_interface, only: n_param
    use otr_qn, only: transport_type, init_hess_type, standard_update_orbs_qn_factory &
                      => update_orbs_qn_factory, standard_update_orbs_qn_deconstructor &
                      => update_orbs_qn_deconstructor
    use, intrinsic :: iso_c_binding, only: c_bool, c_char, c_funptr, c_funloc, &
                                           c_f_procpointer, c_associated, c_null_funptr

    implicit none

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function transport_c_type(geodesic_c, tangent_vector_c) &
            result(error_c) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in), target :: geodesic_c(*)
            real(c_rp), intent(inout), target :: tangent_vector_c(*)
            integer(c_ip) :: error_c
        end function transport_c_type
    end interface

    abstract interface
        function init_hess_c_type(vector_c) result(error_c) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(inout), target :: vector_c(*)
            integer(c_ip) :: error_c
        end function init_hess_c_type
    end interface

    ! define procedure pointer which will point to the Fortran procedures
    procedure(update_orbs_c_type), pointer :: update_orbs_orig_qn_before_wrapping => &
        null()
    procedure(transport_c_type), pointer :: transport_before_wrapping => null()
    procedure(init_hess_c_type), pointer :: init_hess_before_wrapping => null()
    procedure(update_orbs_type), pointer :: update_orbs_qn_before_wrapping => null()
    procedure(hess_x_type), pointer :: hess_x_qn_before_wrapping => null()

    ! derived type for quasi-Newton settings
    type, bind(C) :: qn_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose, max_points
        character(c_char) :: hess_update_scheme(kw_len + 1)
    end type

    procedure(standard_update_orbs_qn_factory), pointer :: update_orbs_qn_factory => &
        standard_update_orbs_qn_factory
    procedure(standard_update_orbs_qn_deconstructor), pointer :: &
        update_orbs_qn_deconstructor => standard_update_orbs_qn_deconstructor

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_type), pointer :: update_orbs_orig_qn_f_wrapper_ptr => &
        update_orbs_orig_qn_f_wrapper
    procedure(transport_type), pointer :: transport_f_wrapper_ptr => transport_f_wrapper
    procedure(init_hess_type), pointer :: init_hess_f_wrapper_ptr => init_hess_f_wrapper
    procedure(update_orbs_c_type), pointer :: update_orbs_qn_c_wrapper_ptr => &
        update_orbs_qn_c_wrapper
    procedure(hess_x_c_type), pointer :: hess_x_qn_c_wrapper_ptr => &
        hess_x_qn_c_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_qn_f_c
        module procedure assign_qn_c_f
    end interface

contains

    function update_orbs_qn_factory_c_wrapper(update_orbs_orig_c_funptr, &
                                              transport_c_funptr, init_hess_c_funptr, &
                                              n_param_c, settings_c, &
                                              update_orbs_qn_c_funptr) result(error_c) &
        bind(C, name="update_orbs_qn_factory")
        !
        ! this subroutine wraps the factory function for the subroutine to convert C 
        ! variables to Fortran variables
        !
        use otr_qn, only: qn_settings_type

        type(c_funptr), intent(in), value :: update_orbs_orig_c_funptr, &
                                             transport_c_funptr, init_hess_c_funptr
        integer(c_ip), intent(in), value :: n_param_c
        type(qn_settings_type_c), intent(inout) :: settings_c
        type(c_funptr), intent(out) :: update_orbs_qn_c_funptr
        integer(c_ip) :: error_c

        procedure(update_orbs_type), pointer :: update_orbs_funptr, &
                                                update_orbs_qn_funptr
        procedure(transport_type), pointer :: transport_funptr
        procedure(init_hess_type), pointer :: init_hess_funptr
        type(qn_settings_type) :: settings
        integer(ip) :: error

        ! associate the input C pointer to update_orbs subroutine to a Fortran
        ! procedure pointer
        call c_f_procpointer(cptr=update_orbs_orig_c_funptr, &
                             fptr=update_orbs_orig_qn_before_wrapping)
        call c_f_procpointer(cptr=transport_c_funptr, fptr=transport_before_wrapping)
        call c_f_procpointer(cptr=init_hess_c_funptr, fptr=init_hess_before_wrapping)

        ! associate procedure pointer to wrapper function
        update_orbs_funptr => update_orbs_orig_qn_f_wrapper
        transport_funptr => transport_f_wrapper
        init_hess_funptr => init_hess_f_wrapper

        ! convert number of parameters to Fortran kind and store globally to access 
        ! assumed size arrays passed from C to Fortran
        n_param = int(n_param_c, kind=ip)

        ! convert settings
        settings = settings_c

        ! associate the global procedure pointer to the update_orbs function
        call update_orbs_qn_factory(update_orbs_funptr, transport_funptr, &
                                    init_hess_funptr, n_param, error, settings, &
                                    update_orbs_qn_funptr)
        update_orbs_qn_before_wrapping => update_orbs_qn_funptr

        ! get a C function pointer to the update_orbs wrapper function
        update_orbs_qn_c_funptr = c_funloc(update_orbs_qn_c_wrapper)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function update_orbs_qn_factory_c_wrapper

    subroutine update_orbs_orig_qn_f_wrapper(kappa, func, grad, h_diag, hess_x, error)
        !
        ! this subroutine wraps the orbital update subroutine to convert Fortran 
        ! variables to C variables
        !
        use opentrustregion, only: hess_x_type
        use c_interface, only: update_orbs_f_wrapper_impl

        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x
        integer(ip), intent(out) :: error

        call update_orbs_f_wrapper_impl(update_orbs_orig_qn_before_wrapping, kappa, &
                                        func, grad, h_diag, hess_x, error)

    end subroutine update_orbs_orig_qn_f_wrapper

    subroutine transport_f_wrapper(geodesic, tangent_vector, error)
        !
        ! this subroutine wraps the transport subroutine to convert Fortran variables 
        ! to C variables
        !
        real(rp), intent(in), target :: geodesic(:)
        real(rp), intent(inout), target :: tangent_vector(:)
        integer(ip), intent(out) :: error

        real(c_rp), pointer :: geodesic_c(:), tangent_vector_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            geodesic_c => geodesic
            tangent_vector_c => tangent_vector
        else
            allocate(geodesic_c(size(geodesic, 1)))
            allocate(tangent_vector_c(size(tangent_vector, 1)))
            geodesic_c = real(geodesic, kind=c_rp)
            tangent_vector_c = real(tangent_vector, kind=c_rp)
        end if

        ! call transport C function
        error_c = transport_before_wrapping(geodesic_c, tangent_vector_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            tangent_vector = real(tangent_vector_c, kind=rp)
            deallocate(geodesic_c)
            deallocate(tangent_vector_c)
        end if

    end subroutine transport_f_wrapper

    subroutine init_hess_f_wrapper(vector, error)
        !
        ! this subroutine wraps the initial Hessian subroutine to convert Fortran 
        ! variables to C variables
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error

        real(c_rp), pointer :: vector_c(:)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            vector_c => vector
        else
            allocate(vector_c(size(vector, 1)))
            vector_c = real(vector, kind=c_rp)
        end if

        ! call initial Hessian C function
        error_c = init_hess_before_wrapping(vector_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            vector = real(vector_c, kind=rp)
            deallocate(vector_c)
        end if

    end subroutine init_hess_f_wrapper

    function update_orbs_qn_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                         hess_x_c_funptr) result(error_c) bind(C)
        !
        ! this function wraps the orbital update subroutine to convert Fortran 
        ! variables to C variables
        !
        use otr_common_c_interface, only: update_orbs_c_wrapper_impl
        
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(out) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        error_c = update_orbs_c_wrapper_impl(update_orbs_qn_before_wrapping, &
                                             hess_x_qn_before_wrapping, &
                                             hess_x_qn_c_wrapper, kappa_c, func_c, &
                                             grad_c, h_diag_c, hess_x_c_funptr)

    end function update_orbs_qn_c_wrapper

    function hess_x_qn_c_wrapper(x_c, hess_x_c) result(error_c) bind(C)
        !
        ! this function wraps the Hessian linear transformation to convert Fortran 
        ! variables to C variables
        !
        use otr_common_c_interface, only: hess_x_c_wrapper_impl

        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        integer(c_ip) :: error_c

        error_c = hess_x_c_wrapper_impl(hess_x_qn_before_wrapping, x_c, hess_x_c)

    end function hess_x_qn_c_wrapper

    subroutine init_qn_settings_c(settings_c) bind(C, name="init_qn_settings")
        !
        ! this subroutine initializes the quasi-Newton solver settings
        !
        use otr_qn, only: default_qn_settings

        type(qn_settings_type_c), intent(inout) :: settings_c
    
        settings_c = default_qn_settings

    end subroutine init_qn_settings_c

    subroutine update_orbs_qn_deconstructor_c_wrapper() &
        bind(C, name="update_orbs_qn_deconstructor")
        !
        ! this subroutine deallocates the quasi-Newton objects
        !
        call update_orbs_qn_deconstructor()

    end subroutine update_orbs_qn_deconstructor_c_wrapper

    subroutine assign_qn_f_c(settings, settings_c)
        !
        ! this subroutine converts quasi-Newton settings from C to Fortran
        !
        use otr_qn, only: qn_settings_type
        use c_interface, only: logger_before_wrapping, logger_f_wrapper, &
                               character_from_c

        type(qn_settings_type), intent(out) :: settings
        type(qn_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%logger)) then
                call c_f_procpointer(cptr=settings_c%logger, &
                                     fptr=logger_before_wrapping)
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert integers
            settings%verbose = int(settings_c%verbose, kind=ip)
            settings%max_points = int(settings_c%max_points, kind=ip)

            ! convert characters
            settings%hess_update_scheme = &
                character_from_c(settings_c%hess_update_scheme)

            ! set settings to initialized
            settings%initialized = .true.
        end if

    end subroutine assign_qn_f_c

    subroutine assign_qn_c_f(settings_c, settings)
        !
        ! this subroutine converts quasi-Newton settings from Fortran to C
        !
        use otr_qn, only: qn_settings_type
        use c_interface, only: character_to_c

        type(qn_settings_type_c), intent(out) :: settings_c
        type(qn_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%logger = c_null_funptr

            ! convert integers
            settings_c%verbose = int(settings%verbose, kind=c_ip)
            settings_c%max_points = int(settings%max_points, kind=c_ip)

            ! convert characters
            settings_c%hess_update_scheme = character_to_c(settings%hess_update_scheme)

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_qn_c_f

end module otr_qn_c_interface
