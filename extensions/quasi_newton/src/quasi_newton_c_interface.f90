! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_qn_c_interface

    use opentrustregion, only: ip, rp, kw_len, update_orbs_type, obj_func_type
    use c_interface, only: c_ip, c_rp, c_callbacks_type
    use otr_common_c_interface, only: extension_holder_type
    use otr_qn, only: &
        transport_type, init_hess_type, standard_qn_factory => qn_factory, &
        standard_update_orbs_qn_deconstructor => update_orbs_qn_deconstructor
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_char, c_funptr, c_funloc, &
                                           c_f_procpointer, c_associated, &
                                           c_null_funptr, c_null_ptr

    implicit none

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function transport_c_type(geodesic_c, tangent_vector_c, context_c) &
            result(error_c) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: geodesic_c(*)
            real(c_rp), intent(inout), target :: tangent_vector_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error_c
        end function transport_c_type
    end interface

    abstract interface
        function init_hess_c_type(vector_c, context_c) result(error_c) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(inout), target :: vector_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error_c
        end function init_hess_c_type
    end interface

    ! bundle of the C callback functions of quasi-Newton, which extends that of the
    ! core by the transport and initial Hessian functions
    type, extends(c_callbacks_type) :: qn_c_callbacks_type
        procedure(transport_c_type), pointer, nopass :: transport => null()
        procedure(init_hess_c_type), pointer, nopass :: init_hess => null()
    end type

    ! callback bundle of the transport and initial Hessian functions, and the logging
    ! function and host context of the quasi-Newton settings, held here since the
    ! settings outlive the factory call that sets it
    type(qn_c_callbacks_type), target :: qn_callbacks

    ! holder of the quasi-Newton object for C hosts, whose levels are the contexts of
    ! the C callback functions the factory returns and wires into the C solver settings
    type(extension_holder_type), target :: qn_holder

    ! derived type for quasi-Newton settings
    type, bind(C) :: qn_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose, max_points
        character(kind=c_char) :: hess_update_scheme(kw_len + 1)
        type(c_ptr) :: context
    end type

    procedure(standard_qn_factory), pointer :: qn_factory => standard_qn_factory
    procedure(standard_update_orbs_qn_deconstructor), pointer :: &
        update_orbs_qn_deconstructor => standard_update_orbs_qn_deconstructor

    ! create function pointers to ensure that routines comply with interface
    procedure(transport_type), pointer :: transport_f_wrapper_ptr => transport_f_wrapper
    procedure(init_hess_type), pointer :: init_hess_f_wrapper_ptr => init_hess_f_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_qn_f_c
        module procedure assign_qn_c_f
    end interface

contains

    function qn_factory_c_wrapper( &
        n_param_c, update_orbs_c_funptr, obj_func_c_funptr, transport_c_funptr, &
        init_hess_c_funptr, update_orbs_qn_c_funptr, obj_func_qn_c_funptr, &
        solver_settings_c, settings_c) result(error_c) bind(C, name="qn_factory")
        !
        ! this function wraps the quasi-Newton factory for C hosts, whose orbital
        ! update, objective function and C solver settings it hands to the holder of
        ! the quasi-Newton object, and writes the solver settings the factory wired
        ! back into the C solver settings
        !
        use opentrustregion, only: solver_settings_type
        use c_interface, only: solver_settings_type_c, update_orbs_f_wrapper, &
                               obj_func_f_wrapper
        use otr_qn, only: qn_settings_type, update_object
        use otr_common_c_interface, only: setup_extension_holder, &
                                          wire_solver_settings_c, &
                                          update_orbs_c_wrapper, obj_func_c_wrapper

        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             transport_c_funptr, init_hess_c_funptr
        type(c_funptr), intent(out) :: update_orbs_qn_c_funptr, obj_func_qn_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(qn_settings_type_c), intent(in) :: settings_c
        integer(c_ip) :: error_c

        procedure(update_orbs_type), pointer :: update_orbs_funptr, &
                                                update_orbs_qn_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr, obj_func_qn_funptr
        procedure(transport_type), pointer :: transport_funptr
        procedure(init_hess_type), pointer :: init_hess_funptr
        type(solver_settings_type) :: solver_settings
        type(qn_settings_type) :: settings
        type(qn_c_callbacks_type), target :: callbacks
        type(extension_holder_type) :: previous_holder
        type(qn_c_callbacks_type) :: previous_callbacks
        integer(ip) :: error

        ! no functions are returned unless the factory succeeds
        update_orbs_qn_c_funptr = c_null_funptr
        obj_func_qn_c_funptr = c_null_funptr

        ! convert settings and hand the transport and initial Hessian functions and the
        ! logging function with the host context of initialized settings to their
        ! wrappers in a bundle, which the settings point to
        settings = settings_c
        call c_f_procpointer(cptr=transport_c_funptr, fptr=callbacks%transport)
        call c_f_procpointer(cptr=init_hess_c_funptr, fptr=callbacks%init_hess)
        if (settings_c%initialized) then
            if (c_associated(settings_c%logger)) &
                call c_f_procpointer(cptr=settings_c%logger, fptr=callbacks%logger)
            callbacks%host_context = settings_c%context
        end if
        settings%context => callbacks

        ! keep the holder and the bundle held by the module, which are restored if the
        ! factory fails after the holder was set up, so that the functions an earlier
        ! call returned keep working
        previous_holder = qn_holder
        previous_callbacks = qn_callbacks

        ! set up the holder with the host's orbital update and objective function and
        ! convert the solver settings, unless the holder was already wired into these
        call setup_extension_holder(qn_holder, solver_settings_c, solver_settings, &
                                    settings, error, update_orbs_c_funptr, &
                                    obj_func_c_funptr)
        if (error /= 0) then
            error_c = int(error, kind=c_ip)
            return
        end if

        ! the bundle held by the module, in which the wrappers of the callback functions
        ! find them, takes over the bundle only once the holder accepted the solver
        ! settings, so that refused settings leave the callback functions of an object
        ! in use unchanged
        qn_callbacks = callbacks
        settings%context => qn_callbacks

        ! the factory reaches the host's callback functions through the core's
        ! wrappers, which find them in the bundles of the holder
        update_orbs_funptr => update_orbs_f_wrapper
        obj_func_funptr => obj_func_f_wrapper
        transport_funptr => transport_f_wrapper
        init_hess_funptr => init_hess_f_wrapper

        ! call factory function
        call qn_factory(int(n_param_c, kind=ip), update_orbs_funptr, obj_func_funptr, &
                        transport_funptr, init_hess_funptr, update_orbs_qn_funptr, &
                        obj_func_qn_funptr, solver_settings, error, settings)

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update and the objective function of the object, or restore the holder and
        ! the bundle held by the module if the factory failed
        if (error == 0) then
            call wire_solver_settings_c(qn_holder, update_object, solver_settings, &
                                        solver_settings_c)
            update_orbs_qn_c_funptr = c_funloc(update_orbs_c_wrapper)
            obj_func_qn_c_funptr = c_funloc(obj_func_c_wrapper)
        else
            qn_holder = previous_holder
            qn_callbacks = previous_callbacks
        end if

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function qn_factory_c_wrapper

    subroutine transport_f_wrapper(geodesic, tangent_vector, error, context)
        !
        ! this subroutine wraps the transport subroutine to convert Fortran variables
        ! to C variables
        !
        real(rp), intent(in), target :: geodesic(:)
        real(rp), intent(inout), target :: tangent_vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

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
        error_c = qn_callbacks%transport(geodesic_c, tangent_vector_c, &
                                         qn_callbacks%host_context)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            tangent_vector = real(tangent_vector_c, kind=rp)
            deallocate(geodesic_c)
            deallocate(tangent_vector_c)
        end if

    end subroutine transport_f_wrapper

    subroutine init_hess_f_wrapper(vector, error, context)
        !
        ! this subroutine wraps the initial Hessian subroutine to convert Fortran
        ! variables to C variables
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

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
        error_c = qn_callbacks%init_hess(vector_c, qn_callbacks%host_context)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            vector = real(vector_c, kind=rp)
            deallocate(vector_c)
        end if

    end subroutine init_hess_f_wrapper

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
        use c_interface, only: character_from_c, logger_f_wrapper

        type(qn_settings_type), intent(out) :: settings
        type(qn_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%logger)) then
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

            ! the host context cannot be converted
            settings_c%context = c_null_ptr

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_qn_c_f

end module otr_qn_c_interface
