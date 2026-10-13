! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface

    use opentrustregion, only: ip, rp, update_orbs_type, obj_func_type
    use c_interface, only: c_ip, c_rp, c_callbacks_type
    use otr_common_c_interface, only: extension_holder_type
    use otr_s_gek, only: &
        change_reference_type, standard_s_gek_factory => s_gek_factory, &
        standard_update_orbs_s_gek_deconstructor => update_orbs_s_gek_deconstructor
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_funloc, &
                                           c_f_procpointer, c_associated, &
                                           c_null_funptr, c_null_ptr

    implicit none

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function change_reference_c_type(new_ref_c, n_points_c, kappa_list_c, &
                                         local_grad_list_c, grad_list_c, context_c) &
            result(error_c) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: new_ref_c(*)
            integer(c_ip), intent(in), value :: n_points_c
            real(c_rp), intent(inout), target :: kappa_list_c(*), &
                                                 local_grad_list_c(*), grad_list_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error_c
        end function change_reference_c_type
    end interface

    ! bundle of the C callback functions of S-GEK, which extends that of the core by
    ! the change of reference function
    type, extends(c_callbacks_type) :: s_gek_c_callbacks_type
        procedure(change_reference_c_type), pointer, nopass :: change_reference => &
            null()
    end type

    ! callback bundle of the change of reference function, and the logging function and
    ! host context of the S-GEK settings, held here since the settings outlive the
    ! factory call that sets it
    type(s_gek_c_callbacks_type), target :: s_gek_callbacks

    ! holder of the S-GEK object for C hosts, whose levels are the contexts of the C
    ! callback functions the factory returns and wires into the C solver settings
    type(extension_holder_type), target :: s_gek_holder

    ! derived type for S-GEK settings
    type, bind(C) :: s_gek_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized, use_subspace
        integer(c_ip) :: verbose, max_points
        type(c_ptr) :: context
    end type

    procedure(standard_s_gek_factory), pointer :: s_gek_factory => &
        standard_s_gek_factory
    procedure(standard_update_orbs_s_gek_deconstructor), pointer :: &
        update_orbs_s_gek_deconstructor => standard_update_orbs_s_gek_deconstructor

    ! create function pointers to ensure that routines comply with interface
    procedure(change_reference_type), pointer :: change_reference_f_wrapper_ptr => &
        change_reference_f_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_s_gek_f_c
        module procedure assign_s_gek_c_f
    end interface

contains

    function s_gek_factory_c_wrapper( &
        n_param_c, update_orbs_c_funptr, obj_func_c_funptr, change_reference_c_funptr, &
        update_orbs_s_gek_c_funptr, obj_func_s_gek_c_funptr, solver_settings_c, &
        settings_c) result(error_c) bind(C, name="s_gek_factory")
        !
        ! this function wraps the S-GEK factory for C hosts, whose orbital update,
        ! objective function and C solver settings it hands to the holder of the S-GEK
        ! object, and writes the solver settings the factory wired back into the C
        ! solver settings
        !
        use opentrustregion, only: solver_settings_type
        use c_interface, only: solver_settings_type_c, update_orbs_f_wrapper, &
                               obj_func_f_wrapper
        use otr_s_gek, only: s_gek_settings_type, s_gek_object
        use otr_common_c_interface, only: setup_extension_holder, &
                                          wire_solver_settings_c, &
                                          update_orbs_c_wrapper, obj_func_c_wrapper

        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             change_reference_c_funptr
        type(c_funptr), intent(out) :: update_orbs_s_gek_c_funptr, &
                                       obj_func_s_gek_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(s_gek_settings_type_c), intent(in) :: settings_c
        integer(c_ip) :: error_c

        procedure(update_orbs_type), pointer :: update_orbs_funptr, &
                                                update_orbs_s_gek_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr, obj_func_s_gek_funptr
        procedure(change_reference_type), pointer :: change_reference_funptr
        type(solver_settings_type) :: solver_settings
        type(s_gek_settings_type) :: settings
        type(s_gek_c_callbacks_type), target :: callbacks
        type(extension_holder_type) :: previous_holder
        type(s_gek_c_callbacks_type) :: previous_callbacks
        integer(ip) :: error

        ! no functions are returned unless the factory succeeds
        update_orbs_s_gek_c_funptr = c_null_funptr
        obj_func_s_gek_c_funptr = c_null_funptr

        ! convert settings and hand the change of reference function and the logging
        ! function with the host context of initialized settings to their wrappers in a
        ! bundle, which the settings point to
        settings = settings_c
        call c_f_procpointer(cptr=change_reference_c_funptr, &
                             fptr=callbacks%change_reference)
        if (settings_c%initialized) then
            if (c_associated(settings_c%logger)) &
                call c_f_procpointer(cptr=settings_c%logger, fptr=callbacks%logger)
            callbacks%host_context = settings_c%context
        end if
        settings%context => callbacks

        ! keep the holder and the bundle held by the module, which are restored if the
        ! factory fails after the holder was set up, so that the functions an earlier
        ! call returned keep working
        previous_holder = s_gek_holder
        previous_callbacks = s_gek_callbacks

        ! set up the holder with the host's orbital update and objective function and
        ! convert the solver settings, unless the holder was already wired into these
        call setup_extension_holder(s_gek_holder, solver_settings_c, solver_settings, &
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
        s_gek_callbacks = callbacks
        settings%context => s_gek_callbacks

        ! the factory reaches the host's callback functions through the core's
        ! wrappers, which find them in the bundles of the holder
        update_orbs_funptr => update_orbs_f_wrapper
        obj_func_funptr => obj_func_f_wrapper
        change_reference_funptr => change_reference_f_wrapper

        ! call factory function
        call s_gek_factory(int(n_param_c, kind=ip), update_orbs_funptr, &
                           obj_func_funptr, change_reference_funptr, &
                           update_orbs_s_gek_funptr, obj_func_s_gek_funptr, &
                           solver_settings, error, settings)

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update and the objective function of the object, or restore the holder and
        ! the bundle held by the module if the factory failed
        if (error == 0) then
            call wire_solver_settings_c(s_gek_holder, s_gek_object, solver_settings, &
                                        solver_settings_c)
            update_orbs_s_gek_c_funptr = c_funloc(update_orbs_c_wrapper)
            obj_func_s_gek_c_funptr = c_funloc(obj_func_c_wrapper)
        else
            s_gek_holder = previous_holder
            s_gek_callbacks = previous_callbacks
        end if

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function s_gek_factory_c_wrapper

    subroutine change_reference_f_wrapper(new_ref, n_points, kappa_list, &
                                          local_grad_list, grad_list, error, context)
        !
        ! this subroutine wraps the change reference subroutine to convert Fortran
        ! variables to C variables
        !
        real(rp), intent(in), target :: new_ref(:)
        integer(ip), intent(in) :: n_points
        real(rp), intent(inout), target :: kappa_list(:, :), local_grad_list(:, :), &
                                           grad_list(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp), pointer :: new_ref_c(:), kappa_list_c(:, :), &
                               local_grad_list_c(:, :), grad_list_c(:, :)
        integer(c_ip) :: n_points_c, error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            new_ref_c => new_ref
            n_points_c = n_points
            kappa_list_c => kappa_list
            local_grad_list_c => local_grad_list
            grad_list_c => grad_list
        else
            allocate(new_ref_c(size(new_ref, 1)))
            allocate(kappa_list_c(size(kappa_list, 1), size(kappa_list, 2)))
            allocate( &
                local_grad_list_c(size(local_grad_list, 1), size(local_grad_list, 2)))
            allocate(grad_list_c(size(grad_list, 1), size(grad_list, 2)))
            new_ref_c = real(new_ref, kind=c_rp)
            n_points_c = int(n_points, kind=c_ip)
            kappa_list_c = real(kappa_list, kind=c_rp)
            local_grad_list_c = real(local_grad_list, kind=c_rp)
            grad_list_c = real(grad_list, kind=c_rp)
        end if

        ! call change reference C function
        error_c = s_gek_callbacks%change_reference( &
            new_ref_c, n_points_c, kappa_list_c, local_grad_list_c, grad_list_c, &
            s_gek_callbacks%host_context)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            kappa_list = real(kappa_list_c, kind=rp)
            local_grad_list = real(local_grad_list_c, kind=rp)
            grad_list = real(grad_list_c, kind=rp)
            deallocate(new_ref_c)
            deallocate(kappa_list_c)
            deallocate(local_grad_list_c)
            deallocate(grad_list_c)
        end if

    end subroutine change_reference_f_wrapper

    subroutine init_s_gek_settings_c(settings_c) bind(C, name="init_s_gek_settings")
        !
        ! this subroutine initializes the S-GEK solver settings
        !
        use otr_s_gek, only: default_s_gek_settings

        type(s_gek_settings_type_c), intent(inout) :: settings_c

        settings_c = default_s_gek_settings

    end subroutine init_s_gek_settings_c

    subroutine update_orbs_s_gek_deconstructor_c_wrapper() &
        bind(C, name="update_orbs_s_gek_deconstructor")
        !
        ! this subroutine deallocates the S-GEK objects
        !
        call update_orbs_s_gek_deconstructor()

    end subroutine update_orbs_s_gek_deconstructor_c_wrapper

    subroutine assign_s_gek_f_c(settings, settings_c)
        !
        ! this subroutine converts S-GEK settings from C to Fortran
        !
        use otr_s_gek, only: s_gek_settings_type
        use c_interface, only: logger_f_wrapper

        type(s_gek_settings_type), intent(out) :: settings
        type(s_gek_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%logger)) then
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert logicals
            settings%use_subspace = logical(settings_c%use_subspace)

            ! convert integers
            settings%verbose = int(settings_c%verbose, kind=ip)
            settings%max_points = int(settings_c%max_points, kind=ip)

            ! set settings to initialized
            settings%initialized = .true.
        end if

    end subroutine assign_s_gek_f_c

    subroutine assign_s_gek_c_f(settings_c, settings)
        !
        ! this subroutine converts S-GEK settings from Fortran to C
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(out) :: settings_c
        type(s_gek_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%logger = c_null_funptr

            ! convert logicals
            settings_c%use_subspace = logical(settings%use_subspace, kind=c_bool)

            ! convert integers
            settings_c%verbose = int(settings%verbose, kind=c_ip)
            settings_c%max_points = int(settings%max_points, kind=c_ip)

            ! the host context cannot be converted
            settings_c%context = c_null_ptr

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_s_gek_c_f

end module otr_s_gek_c_interface
