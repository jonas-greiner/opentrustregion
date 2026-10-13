! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_c_interface

    use opentrustregion, only: ip, rp, obj_func_type, update_orbs_type
    use c_interface, only: c_ip, c_rp, update_orbs_c_type
    use otr_oao, only: standard_oao_factory_cs => oao_factory_cs, &
                       standard_oao_factory_os => oao_factory_os, &
                       standard_oao_deconstructor => oao_deconstructor
    use otr_common, only: evaluate_dm_cs_type, evaluate_dm_os_type, &
                          get_response_cs_type, get_response_os_type
    use otr_common_c_interface, only: orbital_c_callbacks_type, extension_holder_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_loc, &
                                           c_f_pointer, c_funloc, c_f_procpointer, &
                                           c_associated, c_null_funptr, c_null_ptr

    implicit none

    ! callback bundle of the density matrix evaluating function, the response function
    ! it returned last, and the logging function and host context of the OAO settings,
    ! held here since the settings outlive the factory call that sets it; the wrappers
    ! of the density matrix evaluating and response functions find it here rather than
    ! through the context they receive, since ARH in the OAO basis replaces the
    ! settings of the shared OAO object, context included
    type(orbital_c_callbacks_type), target :: oao_callbacks

    ! holder of the OAO object for C hosts, whose levels are the contexts of the C
    ! callback functions the factory returns and wires into the C solver settings
    type(extension_holder_type), target :: oao_holder

    ! density matrix of the C host, which the OAO object updates a copy of if the
    ! factory could not hand it the host's array in place, which the C wrapper of the
    ! orbital update then copies back
    real(c_rp), pointer :: dm_ao_3d_c(:, :, :) => null()

    ! derived type for OAO settings
    type, bind(C) :: oao_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose
        type(c_ptr) :: context
    end type

    procedure(standard_oao_factory_cs), pointer :: oao_factory_cs => &
        standard_oao_factory_cs
    procedure(standard_oao_factory_os), pointer :: oao_factory_os => &
        standard_oao_factory_os
    procedure(standard_oao_deconstructor), pointer :: oao_deconstructor => &
        standard_oao_deconstructor

    ! create function pointers to ensure that routines comply with interface
    procedure(update_orbs_c_type), pointer :: update_orbs_oao_c_wrapper_ptr => &
        update_orbs_oao_c_wrapper
    procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_oao_cs_f_wrapper_ptr => &
        evaluate_dm_oao_cs_f_wrapper
    procedure(evaluate_dm_os_type), pointer :: evaluate_dm_oao_os_f_wrapper_ptr => &
        evaluate_dm_oao_os_f_wrapper
    procedure(get_response_cs_type), pointer :: get_response_oao_cs_f_wrapper_ptr => &
        get_response_oao_cs_f_wrapper
    procedure(get_response_os_type), pointer :: get_response_oao_os_f_wrapper_ptr => &
        get_response_oao_os_f_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_oao_f_c
        module procedure assign_oao_c_f
    end interface

contains

    function oao_factory_c_wrapper( &
        dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
        update_orbs_oao_c_funptr, obj_func_oao_c_funptr, solver_settings_c, &
        settings_c) result(error_c) bind(C, name="oao_factory")
        !
        ! this function wraps the OAO factory for C hosts, whose C solver settings it
        ! hands to the holder of the OAO object, and writes the solver settings the
        ! factory wired back into the C solver settings
        !
        use opentrustregion, only: solver_settings_type
        use c_interface, only: solver_settings_type_c
        use otr_oao, only: oao_settings_type, oao_object
        use otr_common_c_interface, only: setup_extension_holder, &
                                          wire_solver_settings_c, obj_func_c_wrapper

        real(c_rp), intent(inout), target :: dm_ao_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_oao_c_funptr, obj_func_oao_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(oao_settings_type_c), intent(in) :: settings_c
        integer(c_ip) :: error_c

        real(rp), pointer, contiguous :: dm_ao_2d(:, :)
        real(rp), pointer, contiguous :: dm_ao_3d(:, :, :)
        real(rp), pointer :: ao_overlap(:, :)
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        type(solver_settings_type) :: solver_settings
        type(oao_settings_type) :: settings
        type(orbital_c_callbacks_type), target :: callbacks
        type(extension_holder_type) :: previous_holder
        type(orbital_c_callbacks_type) :: previous_callbacks
        integer(ip) :: n_particle, n_ao, error

        ! no functions are returned unless the factory succeeds
        update_orbs_oao_c_funptr = c_null_funptr
        obj_func_oao_c_funptr = c_null_funptr

        ! convert dimensions to Fortran kind
        n_particle = int(n_particle_c, kind=ip)
        n_ao = int(n_ao_c, kind=ip)

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            if (n_particle == 1) then
                call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_2d, [n_ao, n_ao])
            else
                call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_3d, [n_ao, n_ao, n_particle])
            end if
            call c_f_pointer(c_loc(ao_overlap_c(1)), ao_overlap, [n_ao, n_ao])
        else
            if (n_particle == 1) then
                allocate(dm_ao_2d(n_ao, n_ao))
                dm_ao_2d = reshape(real(dm_ao_c(:n_ao_c**2), kind=rp), [n_ao, n_ao])
            else
                allocate(dm_ao_3d(n_ao, n_ao, n_particle))
                dm_ao_3d = reshape(real(dm_ao_c(:(n_ao_c**2 * n_particle_c)), &
                                        kind=rp), [n_ao, n_ao, n_particle])
            end if
            allocate(ao_overlap(n_ao, n_ao))
            ao_overlap = reshape(real(ao_overlap_c(:n_ao**2), kind=rp), [n_ao, n_ao])
        end if

        ! associate procedure pointer to wrapper function
        if (n_particle == 1) then
            evaluate_dm_cs_funptr => evaluate_dm_oao_cs_f_wrapper
        else
            evaluate_dm_os_funptr => evaluate_dm_oao_os_f_wrapper
        end if

        ! convert settings and hand the density matrix evaluating function and the
        ! logging function with the host context of initialized settings to their
        ! wrappers in a bundle, which the settings point to
        settings = settings_c
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=callbacks%evaluate_dm)
        if (settings_c%initialized) then
            if (c_associated(settings_c%logger)) &
                call c_f_procpointer(cptr=settings_c%logger, fptr=callbacks%logger)
            callbacks%host_context = settings_c%context
        end if
        settings%context => callbacks

        ! keep the holder and the bundle held by the module, which are restored if the
        ! factory fails after the holder was set up, so that the functions an earlier
        ! call returned keep working
        previous_holder = oao_holder
        previous_callbacks = oao_callbacks

        ! set up the holder and convert the solver settings, unless the holder was
        ! already wired into these
        call setup_extension_holder(oao_holder, solver_settings_c, solver_settings, &
                                    settings, error)
        if (error /= 0) then
            error_c = int(error, kind=c_ip)
            return
        end if

        ! the bundle held by the module, in which the wrappers of the callback functions
        ! find them, takes over the bundle only once the holder accepted the solver
        ! settings, so that refused settings leave the callback functions of an object
        ! in use unchanged
        oao_callbacks = callbacks
        settings%context => oao_callbacks

        ! call factory function
        if (n_particle == 1) then
            call oao_factory_cs(dm_ao_2d, ao_overlap, n_particle, n_ao, &
                                evaluate_dm_cs_funptr, update_orbs_oao_funptr, &
                                obj_func_oao_funptr, solver_settings, error, settings)
        else
            call oao_factory_os(dm_ao_3d, ao_overlap, n_particle, n_ao, &
                                evaluate_dm_os_funptr, update_orbs_oao_funptr, &
                                obj_func_oao_funptr, solver_settings, error, settings)
        end if

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update, which copies the orbitals back to the host's array if the object
        ! updates a copy of them, and of the objective function of the object, or
        ! restore the holder and the bundle held by the module if the factory failed
        if (error == 0) then
            call wire_solver_settings_c(oao_holder, oao_object, solver_settings, &
                                        solver_settings_c)
            update_orbs_oao_c_funptr = c_funloc(update_orbs_oao_c_wrapper)
            obj_func_oao_c_funptr = c_funloc(obj_func_c_wrapper)
            if (rp /= c_rp) call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_3d_c, &
                                             [n_ao, n_ao, n_particle])
        else
            oao_holder = previous_holder
            oao_callbacks = previous_callbacks
        end if

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function oao_factory_c_wrapper

    recursive function update_orbs_oao_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                                 hess_x_c_funptr, context_c) &
        result(error_c) bind(C)
        !
        ! this function exposes the orbital update of the OAO object, and of ARH in the
        ! OAO basis, to C hosts through the shared C wrapper of the orbital update,
        ! after which it copies the density matrix the object updated back to the array
        ! of the C host if the factory could not hand it the host's array in place
        !
        use otr_common_c_interface, only: update_orbs_c_wrapper
        use otr_oao, only: oao_object

        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(inout) :: hess_x_c_funptr
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error_c

        ! update the orbitals
        error_c = update_orbs_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                        hess_x_c_funptr, context_c)

        ! copy the updated density matrix back to the host's array if it is a copy
        if (rp /= c_rp) dm_ao_3d_c = real(oao_object%dm_ao, kind=c_rp)

    end function update_orbs_oao_c_wrapper

    subroutine evaluate_dm_oao_cs_f_wrapper(dm, energy, fock, get_response_cs_funptr, &
                                            error, context)
        !
        ! this subroutine wraps the density matrix evaluating subroutine to convert
        ! Fortran variables to C variables for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target, contiguous :: fock(:, :)
        procedure(get_response_cs_type), intent(out), optional, pointer :: &
            get_response_cs_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp), pointer :: dm_3d(:, :, :), fock_3d(:, :, :)
        procedure(get_response_os_type), pointer :: get_response_funptr

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        nullify(fock_3d, get_response_funptr)
        if (present(fock)) fock_3d(1:size(fock, 1), 1:size(fock, 2), 1:1) => fock
        if (present(get_response_cs_funptr)) then
            call evaluate_dm_oao_os_f_wrapper(dm_3d, energy, fock_3d, &
                                              get_response_funptr, error, context)
            get_response_cs_funptr => null()
            if (error == 0) get_response_cs_funptr => get_response_oao_cs_f_wrapper
        else
            call evaluate_dm_oao_os_f_wrapper(dm_3d, energy, fock_3d, error=error, &
                                              context=context)
        end if

    end subroutine evaluate_dm_oao_cs_f_wrapper

    subroutine evaluate_dm_oao_os_f_wrapper(dm, energy, fock, get_response_funptr, &
                                            error, context)
        !
        ! this subroutine wraps the density matrix evaluating subroutine to convert
        ! Fortran variables to C variables for the open-shell case
        !
        use otr_common_c_interface, only: evaluate_dm_f_wrapper_impl

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target :: fock(:, :, :)
        procedure(get_response_os_type), intent(out), optional, pointer :: &
            get_response_funptr
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        ! call density matrix evaluating C function, which keeps a returned response C
        ! function in the bundle, and associate a Fortran procedure pointer with the
        ! wrapper of the response function
        call evaluate_dm_f_wrapper_impl(oao_callbacks, dm, energy, &
                                        present(get_response_funptr), error, fock)
        if (present(get_response_funptr)) then
            get_response_funptr => null()
            if (error == 0) get_response_funptr => get_response_oao_os_f_wrapper
        end if

    end subroutine evaluate_dm_oao_os_f_wrapper

    subroutine get_response_oao_cs_f_wrapper(dm, response, error, context)
        !
        ! this subroutine wraps the response subroutine to convert Fortran variables to
        ! C variables
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out), target, contiguous :: response(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp), pointer :: dm_3d(:, :, :), response_3d(:, :, :)

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        response_3d(1:size(response, 1), 1:size(response, 2), 1:1) => response
        call get_response_oao_os_f_wrapper(dm_3d, response_3d, error, context)

    end subroutine get_response_oao_cs_f_wrapper

    subroutine get_response_oao_os_f_wrapper(dm, response, error, context)
        !
        ! this subroutine wraps the response subroutine to convert Fortran variables to
        ! C variables
        !
        use otr_common_c_interface, only: get_response_f_wrapper_impl

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out), target :: response(:, :, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        call get_response_f_wrapper_impl(oao_callbacks, dm, response, error)

    end subroutine get_response_oao_os_f_wrapper

    subroutine init_oao_settings_c(settings_c) bind(C, name="init_oao_settings")
        !
        ! this subroutine initializes the OAO settings
        !
        use otr_oao, only: default_oao_settings

        type(oao_settings_type_c), intent(inout) :: settings_c

        settings_c = default_oao_settings

    end subroutine init_oao_settings_c

    subroutine oao_deconstructor_c_wrapper() bind(C, name="oao_deconstructor")
        !
        ! this subroutine deallocates the OAO objects
        !
        call oao_deconstructor()

    end subroutine oao_deconstructor_c_wrapper

    subroutine assign_oao_f_c(settings, settings_c)
        !
        ! this subroutine converts OAO settings from C to Fortran
        !
        use otr_oao, only: oao_settings_type
        use c_interface, only: logger_f_wrapper

        type(oao_settings_type), intent(out) :: settings
        type(oao_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%logger)) then
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert integers
            settings%verbose = int(settings_c%verbose, kind=ip)

            ! set settings to initialized
            settings%initialized = .true.
        end if

    end subroutine assign_oao_f_c

    subroutine assign_oao_c_f(settings_c, settings)
        !
        ! this subroutine converts OAO settings from Fortran to C
        !
        use otr_oao, only: oao_settings_type

        type(oao_settings_type_c), intent(out) :: settings_c
        type(oao_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%logger = c_null_funptr

            ! convert integers
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! the host context cannot be converted
            settings_c%context = c_null_ptr

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_oao_c_f

end module otr_oao_c_interface
