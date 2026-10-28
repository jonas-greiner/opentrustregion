! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_rotation_qn_c_interface

    use opentrustregion_core, only: ip, rp, kw_len, update_orbs_type, obj_func_type
    use c_interface, only: c_ip, c_rp, c_callbacks_type
    use otr_rotation_qn, only: transport_type, init_hess_type, &
                               standard_rotation_qn_factory => rotation_qn_factory
    use, intrinsic :: iso_c_binding, only: c_size_t, c_sizeof, c_bool, c_ptr, c_char, &
                                           c_funptr, c_f_procpointer, c_associated, &
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

    ! bundle of the C callback functions of rotation QN, which extends that of the core
    ! by the transport and initial Hessian functions, held by the holder of the object,
    ! since the settings outlive the factory call, and reached through the context of
    ! the settings
    type, extends(c_callbacks_type) :: rotation_qn_c_callbacks_type
        procedure(transport_c_type), pointer, nopass :: transport => null()
        procedure(init_hess_c_type), pointer, nopass :: init_hess => null()
    end type

    ! derived type for rotation QN settings
    type, bind(C) :: rotation_qn_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose, max_points
        character(kind=c_char) :: qn_method(kw_len + 1)
        type(c_ptr) :: context
    end type

    procedure(standard_rotation_qn_factory), pointer :: rotation_qn_factory => &
        standard_rotation_qn_factory

    ! create function pointers to ensure that routines comply with interface
    procedure(transport_type), pointer :: transport_f_wrapper_ptr => transport_f_wrapper
    procedure(init_hess_type), pointer :: init_hess_f_wrapper_ptr => init_hess_f_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_rotation_qn_f_c
        module procedure assign_rotation_qn_c_f
    end interface

contains

    function rotation_qn_factory_c_wrapper( &
        n_param_c, update_orbs_c_funptr, obj_func_c_funptr, transport_c_funptr, &
        init_hess_c_funptr, update_orbs_rotation_qn_c_funptr, &
        obj_func_rotation_qn_c_funptr, extension_c, solver_settings_c, solver_state_c, &
        settings_c) result(error_c) bind(C, name="rotation_qn_factory")
        !
        ! this function wraps the rotation QN factory for C hosts, which creates a
        ! holder of the new rotation QN object, returned as the extension the host
        ! releases, and hands it the host's orbital update, objective function, C
        ! solver settings and record of the solver state, and writes the solver
        ! settings the factory wired back into the C solver settings
        !
        use opentrustregion_core, only: warn_uninitialized_settings, add_error_origin, &
                                        verbosity_error
        use c_interface, only: solver_settings_type_c, update_orbs_f_wrapper, &
                               obj_func_f_wrapper
        use otr_rotation_qn, only: &
            rotation_qn_settings_type, default_rotation_qn_settings, &
            rotation_qn_settings_uninitialized_warning_msg, error_rotation_qn
        use otr_common_c_interface, only: extension_holder_type, bundle_settings_c, &
                                          setup_extension_holder, &
                                          finish_extension_holder

        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             transport_c_funptr, init_hess_c_funptr
        type(c_funptr), intent(out) :: update_orbs_rotation_qn_c_funptr, &
                                       obj_func_rotation_qn_c_funptr
        type(c_ptr), intent(out) :: extension_c
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_ptr), intent(in), value :: solver_state_c
        type(rotation_qn_settings_type_c), intent(in), optional :: settings_c
        integer(c_ip) :: error_c

        procedure(update_orbs_type), pointer :: update_orbs_funptr, &
                                                update_orbs_rotation_qn_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr, &
                                             obj_func_rotation_qn_funptr
        procedure(transport_type), pointer :: transport_funptr
        procedure(init_hess_type), pointer :: init_hess_funptr
        type(rotation_qn_settings_type) :: settings
        type(rotation_qn_settings_type_c) :: local_settings_c
        type(extension_holder_type), pointer :: holder
        integer(ip) :: error

        ! no functions or extension are returned unless the factory succeeds
        update_orbs_rotation_qn_c_funptr = c_null_funptr
        obj_func_rotation_qn_c_funptr = c_null_funptr
        extension_c = c_null_ptr

        ! work on a copy of the given settings or, without any, of the default settings
        if (present(settings_c)) then
            local_settings_c = settings_c
        else
            local_settings_c = default_rotation_qn_settings
        end if

        ! convert the settings, which replaces settings that were not initialized by the
        ! default settings, and print the warning for these on stderr, since they
        ! provide no logging function
        settings = local_settings_c
        if (.not. local_settings_c%initialized) call warn_uninitialized_settings( &
            settings, rotation_qn_settings_uninitialized_warning_msg)

        ! create the holder, whose bundle hands the transport and initial Hessian
        ! functions and the logging function with the host context of initialized
        ! settings to their wrappers, which the settings point to
        allocate(holder)
        allocate(rotation_qn_c_callbacks_type :: holder%extension_callbacks)
        select type (callbacks => holder%extension_callbacks)
        type is (rotation_qn_c_callbacks_type)
            call c_f_procpointer(cptr=transport_c_funptr, fptr=callbacks%transport)
            call c_f_procpointer(cptr=init_hess_c_funptr, fptr=callbacks%init_hess)
            call bundle_settings_c(callbacks, local_settings_c%initialized, &
                                   local_settings_c%logger, local_settings_c%context)
        end select
        settings%context => holder%extension_callbacks

        ! set up the holder with the host's orbital update, objective function and
        ! record of the state and convert the solver settings
        call setup_extension_holder(holder, solver_settings_c, solver_state_c, &
                                    settings, error, update_orbs_c_funptr, &
                                    obj_func_c_funptr)

        ! the object can not run without the host's orbital update, objective
        ! function, transport and initial Hessian functions
        if (error == 0 .and. .not. (c_associated(update_orbs_c_funptr) .and. &
                                    c_associated(obj_func_c_funptr) .and. &
                                    c_associated(transport_c_funptr) .and. &
                                    c_associated(init_hess_c_funptr))) then
            call settings%log("Orbital update, objective function, transport and "// &
                              "initial Hessian functions have to be given.", &
                              verbosity_error, .true.)
            error = 1
        end if
        call add_error_origin(error, error_rotation_qn, settings)

        ! the factory wires the solver settings of the holder, which then own the
        ! object, reaches the host's callback functions through the core's wrappers,
        ! which find them in the bundles of the holder, and reads the state of the
        ! holder
        if (error == 0) then
            update_orbs_funptr => update_orbs_f_wrapper
            obj_func_funptr => obj_func_f_wrapper
            transport_funptr => transport_f_wrapper
            init_hess_funptr => init_hess_f_wrapper
            call rotation_qn_factory( &
                int(n_param_c, kind=ip), update_orbs_funptr, obj_func_funptr, &
                transport_funptr, init_hess_funptr, update_orbs_rotation_qn_funptr, &
                obj_func_rotation_qn_funptr, holder%solver_settings, holder%state, &
                error, settings)
        end if

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update and the objective function of the object and the holder, or free the
        ! holder if the factory failed
        call finish_extension_holder(holder, error, solver_settings_c, &
                                     update_orbs_rotation_qn_c_funptr, &
                                     obj_func_rotation_qn_c_funptr, extension_c)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function rotation_qn_factory_c_wrapper

    subroutine transport_f_wrapper(geodesic, tangent_vector, error, context)
        !
        ! this subroutine exposes the transport C function to Fortran, finding the C
        ! function and the host context in the bundle it receives as context
        !
        real(rp), intent(in), target :: geodesic(:)
        real(rp), intent(inout), target :: tangent_vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(c_ip) :: error_c

        ! the bundle holds the C function
        error = 1
        if (.not. associated(context)) return
        select type (callbacks => context)
        type is (rotation_qn_c_callbacks_type)
            ! call transport C function
            error_c = &
                callbacks%transport(geodesic, tangent_vector, callbacks%host_context)

            ! convert error code to Fortran kind
            error = int(error_c, kind=ip)
        end select

    end subroutine transport_f_wrapper

    subroutine init_hess_f_wrapper(vector, error, context)
        !
        ! this subroutine exposes the initial Hessian C function to Fortran, finding
        ! the C function and the host context in the bundle it receives as context
        !
        real(rp), intent(inout), target :: vector(:)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        integer(c_ip) :: error_c

        ! the bundle holds the C function
        error = 1
        if (.not. associated(context)) return
        select type (callbacks => context)
        type is (rotation_qn_c_callbacks_type)
            ! call initial Hessian C function
            error_c = callbacks%init_hess(vector, callbacks%host_context)

            ! convert error code to Fortran kind
            error = int(error_c, kind=ip)
        end select

    end subroutine init_hess_f_wrapper

    subroutine init_rotation_qn_settings_c(settings_c) &
        bind(C, name="init_rotation_qn_settings")
        !
        ! this subroutine initializes the rotation QN settings
        !
        use otr_rotation_qn, only: default_rotation_qn_settings

        type(rotation_qn_settings_type_c), intent(out) :: settings_c

        settings_c = default_rotation_qn_settings

    end subroutine init_rotation_qn_settings_c

    integer(c_size_t) function rotation_qn_settings_size() bind(C)
        !
        ! this function returns the size of the C rotation QN settings, which the
        ! Python interface compares with its own to detect a library of another version
        ! and the C layout tests with the header's
        !
        type(rotation_qn_settings_type_c) :: settings_c

        rotation_qn_settings_size = c_sizeof(settings_c)

    end function rotation_qn_settings_size

    subroutine assign_rotation_qn_f_c(settings, settings_c)
        !
        ! this subroutine converts rotation QN settings from C to Fortran
        !
        use otr_rotation_qn, only: rotation_qn_settings_type, &
                                   default_rotation_qn_settings
        use c_interface, only: character_from_c, logger_f_wrapper

        type(rotation_qn_settings_type), intent(out) :: settings
        type(rotation_qn_settings_type_c), intent(in) :: settings_c

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
            settings%qn_method = character_from_c(settings_c%qn_method)
        else
            ! C settings that were not initialized hold no values and are converted
            ! as the default settings
            settings = default_rotation_qn_settings
        end if

    end subroutine assign_rotation_qn_f_c

    subroutine assign_rotation_qn_c_f(settings_c, settings)
        !
        ! this subroutine converts rotation QN settings from Fortran to C
        !
        use otr_rotation_qn, only: rotation_qn_settings_type
        use c_interface, only: character_to_c

        type(rotation_qn_settings_type_c), intent(out) :: settings_c
        type(rotation_qn_settings_type), intent(in) :: settings

        ! callback functions cannot be converted
        settings_c%logger = c_null_funptr

        ! convert integers
        settings_c%verbose = int(settings%verbose, kind=c_ip)
        settings_c%max_points = int(settings%max_points, kind=c_ip)

        ! convert characters
        settings_c%qn_method = character_to_c(settings%qn_method)

        ! the host context cannot be converted
        settings_c%context = c_null_ptr

        ! Fortran settings are always initialized
        settings_c%initialized = .true._c_bool

    end subroutine assign_rotation_qn_c_f

end module otr_rotation_qn_c_interface
