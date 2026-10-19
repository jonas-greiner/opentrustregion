! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_s_gek_c_interface

    use opentrustregion, only: ip, rp, update_orbs_type, obj_func_type
    use c_interface, only: c_ip, c_rp, c_callbacks_type
    use otr_s_gek, only: change_reference_type, standard_s_gek_factory => s_gek_factory
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_f_procpointer, &
                                           c_associated, c_null_funptr, c_null_ptr

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
    ! the change of reference function, held by the holder of the object, since the
    ! settings outlive the factory call, and reached through the context of the
    ! settings
    type, extends(c_callbacks_type) :: s_gek_c_callbacks_type
        procedure(change_reference_c_type), pointer, nopass :: change_reference => &
            null()
    end type

    ! derived type for S-GEK settings
    type, bind(C) :: s_gek_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized, use_subspace
        integer(c_ip) :: verbose, max_points
        type(c_ptr) :: context
    end type

    procedure(standard_s_gek_factory), pointer :: s_gek_factory => &
        standard_s_gek_factory

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
        update_orbs_s_gek_c_funptr, obj_func_s_gek_c_funptr, extension_c, &
        solver_settings_c, solver_state_c, settings_c) result(error_c) &
        bind(C, name="s_gek_factory")
        !
        ! this function wraps the S-GEK factory for C hosts, which creates a holder of
        ! the new S-GEK object, returned as the extension the host releases, and hands
        ! it the host's orbital update, objective function, C solver settings and record
        ! of the solver state, and writes the solver settings the factory wired back
        ! into the C solver settings
        !
        use opentrustregion, only: warn_uninitialized_settings, add_error_origin
        use otr_common, only: error_s_gek
        use c_interface, only: solver_settings_type_c, update_orbs_f_wrapper, &
                               obj_func_f_wrapper
        use otr_s_gek, only: s_gek_settings_type, default_s_gek_settings, &
                             s_gek_settings_uninitialized_warning_msg
        use otr_common_c_interface, only: extension_holder_type, bundle_settings_c, &
                                          setup_extension_holder, &
                                          finish_extension_holder

        integer(c_ip), intent(in), value :: n_param_c
        type(c_funptr), intent(in), value :: update_orbs_c_funptr, obj_func_c_funptr, &
                                             change_reference_c_funptr
        type(c_funptr), intent(out) :: update_orbs_s_gek_c_funptr, &
                                       obj_func_s_gek_c_funptr
        type(c_ptr), intent(out) :: extension_c
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_ptr), intent(in), value :: solver_state_c
        type(s_gek_settings_type_c), intent(in), optional :: settings_c
        integer(c_ip) :: error_c

        procedure(update_orbs_type), pointer :: update_orbs_funptr, &
                                                update_orbs_s_gek_funptr
        procedure(obj_func_type), pointer :: obj_func_funptr, obj_func_s_gek_funptr
        procedure(change_reference_type), pointer :: change_reference_funptr
        type(s_gek_settings_type) :: settings
        type(s_gek_settings_type_c) :: local_settings_c
        type(extension_holder_type), pointer :: holder
        integer(ip) :: error

        ! no functions or extension are returned unless the factory succeeds
        update_orbs_s_gek_c_funptr = c_null_funptr
        obj_func_s_gek_c_funptr = c_null_funptr
        extension_c = c_null_ptr

        ! work on a copy of the given settings or, without any, of the default settings
        if (present(settings_c)) then
            local_settings_c = settings_c
        else
            local_settings_c = default_s_gek_settings
        end if

        ! convert the settings, which replaces settings that were not initialized by the
        ! default settings, and print the warning for these on stderr, since they
        ! provide no logging function
        settings = local_settings_c
        if (.not. local_settings_c%initialized) call warn_uninitialized_settings( &
            settings, s_gek_settings_uninitialized_warning_msg)

        ! create the holder, whose bundle hands the change of reference function and the
        ! logging function with the host context of initialized settings to their
        ! wrappers, which the settings point to
        allocate(holder)
        allocate(s_gek_c_callbacks_type :: holder%extension_callbacks)
        select type (callbacks => holder%extension_callbacks)
        type is (s_gek_c_callbacks_type)
            call c_f_procpointer(cptr=change_reference_c_funptr, &
                                 fptr=callbacks%change_reference)
            call bundle_settings_c(callbacks, local_settings_c%initialized, &
                                   local_settings_c%logger, local_settings_c%context)
        end select
        settings%context => holder%extension_callbacks

        ! set up the holder with the host's orbital update, objective function and
        ! record of the state and convert the solver settings
        call setup_extension_holder(holder, solver_settings_c, solver_state_c, &
                                    settings, error, update_orbs_c_funptr, &
                                    obj_func_c_funptr)
        call add_error_origin(error, error_s_gek, settings)

        ! the factory wires the solver settings of the holder, which then own the
        ! object, reaches the host's callback functions through the core's wrappers,
        ! which find them in the bundles of the holder, and reads the state of the
        ! holder
        if (error == 0) then
            update_orbs_funptr => update_orbs_f_wrapper
            obj_func_funptr => obj_func_f_wrapper
            change_reference_funptr => change_reference_f_wrapper
            call s_gek_factory(int(n_param_c, kind=ip), update_orbs_funptr, &
                               obj_func_funptr, change_reference_funptr, &
                               update_orbs_s_gek_funptr, obj_func_s_gek_funptr, &
                               holder%solver_settings, holder%state, error, settings)
        end if

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update and the objective function of the object and the holder, or free the
        ! holder if the factory failed
        call finish_extension_holder(holder, error, solver_settings_c, &
                                     update_orbs_s_gek_c_funptr, &
                                     obj_func_s_gek_c_funptr, extension_c)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function s_gek_factory_c_wrapper

    subroutine change_reference_f_wrapper(new_ref, n_points, kappa_list, &
                                          local_grad_list, grad_list, error, context)
        !
        ! this subroutine wraps the change reference subroutine to convert Fortran
        ! variables to C variables, finding the C function and the host context in the
        ! bundle it receives as context
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

        ! the bundle holds the C function
        error = 1
        if (.not. associated(context)) return
        select type (callbacks => context)
        type is (s_gek_c_callbacks_type)
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
                allocate(local_grad_list_c(size(local_grad_list, 1), &
                                           size(local_grad_list, 2)))
                allocate(grad_list_c(size(grad_list, 1), size(grad_list, 2)))
                new_ref_c = real(new_ref, kind=c_rp)
                n_points_c = int(n_points, kind=c_ip)
                kappa_list_c = real(kappa_list, kind=c_rp)
                local_grad_list_c = real(local_grad_list, kind=c_rp)
                grad_list_c = real(grad_list, kind=c_rp)
            end if

            ! call change reference C function
            error_c = callbacks%change_reference(new_ref_c, n_points_c, kappa_list_c, &
                                                 local_grad_list_c, grad_list_c, &
                                                 callbacks%host_context)

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
        end select

    end subroutine change_reference_f_wrapper

    subroutine init_s_gek_settings_c(settings_c) bind(C, name="init_s_gek_settings")
        !
        ! this subroutine initializes the S-GEK settings
        !
        use otr_s_gek, only: default_s_gek_settings

        type(s_gek_settings_type_c), intent(inout) :: settings_c

        settings_c = default_s_gek_settings

    end subroutine init_s_gek_settings_c

    subroutine assign_s_gek_f_c(settings, settings_c)
        !
        ! this subroutine converts S-GEK settings from C to Fortran
        !
        use otr_s_gek, only: s_gek_settings_type, default_s_gek_settings
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
        else
            ! C settings that were not initialized hold no values and are converted
            ! as the default settings
            settings = default_s_gek_settings
        end if

    end subroutine assign_s_gek_f_c

    subroutine assign_s_gek_c_f(settings_c, settings)
        !
        ! this subroutine converts S-GEK settings from Fortran to C
        !
        use otr_s_gek, only: s_gek_settings_type

        type(s_gek_settings_type_c), intent(out) :: settings_c
        type(s_gek_settings_type), intent(in) :: settings

        ! callback functions cannot be converted
        settings_c%logger = c_null_funptr

        ! convert logicals
        settings_c%use_subspace = logical(settings%use_subspace, kind=c_bool)

        ! convert integers
        settings_c%verbose = int(settings%verbose, kind=c_ip)
        settings_c%max_points = int(settings%max_points, kind=c_ip)

        ! the host context cannot be converted
        settings_c%context = c_null_ptr

        ! Fortran settings are always initialized
        settings_c%initialized = .true._c_bool

    end subroutine assign_s_gek_c_f

end module otr_s_gek_c_interface
