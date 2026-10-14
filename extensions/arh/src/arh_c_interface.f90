! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface

    use opentrustregion, only: ip, rp, kw_len, obj_func_type, update_orbs_type
    use c_interface, only: c_ip, c_rp, c_callbacks_type
    use otr_common_c_interface, only: extension_holder_type
    use otr_arh, only: standard_arh_factory_mo_cs => arh_factory_mo_cs, &
                       standard_arh_factory_mo_os => arh_factory_mo_os, &
                       standard_arh_factory_oao_cs => arh_factory_oao_cs, &
                       standard_arh_factory_oao_os => arh_factory_oao_os, &
                       standard_arh_deconstructor => arh_deconstructor, &
                       evaluate_dm_os_type, evaluate_dm_cs_type
    use, intrinsic :: iso_c_binding, only: &
        c_bool, c_ptr, c_funptr, c_loc, c_f_pointer, c_funloc, c_f_procpointer, &
        c_associated, c_char, c_null_funptr, c_null_ptr

    implicit none

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function evaluate_dm_c_type(dm_ao_c, energy_c, fock_c, v_coulomb_c, &
                                    v_exchange_c, v_nonlinear_c, context_c) &
            result(error_c) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            real(c_rp), intent(out), optional :: fock_c(*), v_coulomb_c(*), &
                                                 v_exchange_c(*), v_nonlinear_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error_c
        end function evaluate_dm_c_type
    end interface

    ! bundle of the C callback functions of ARH, which extends that of the core by the
    ! density matrix evaluating function
    type, extends(c_callbacks_type) :: arh_c_callbacks_type
        procedure(evaluate_dm_c_type), pointer, nopass :: evaluate_dm => null()
    end type

    ! callback bundle of the density matrix evaluating function, and the logging
    ! function and host context of the ARH settings, held here since the settings
    ! outlive the factory call that sets it
    type(arh_c_callbacks_type), target :: arh_callbacks

    ! holder of the ARH object for C hosts, whose levels are the contexts of the C
    ! callback functions the factories return and wire into the C solver settings
    type(extension_holder_type), target :: arh_holder

    ! derived type for ARH settings
    type, bind(C) :: arh_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose
        character(kind=c_char) :: arh_type(kw_len + 1)
        type(c_ptr) :: context
    end type

    procedure(standard_arh_factory_mo_cs), pointer :: arh_factory_mo_cs => &
        standard_arh_factory_mo_cs
    procedure(standard_arh_factory_mo_os), pointer :: arh_factory_mo_os => &
        standard_arh_factory_mo_os
    procedure(standard_arh_factory_oao_cs), pointer :: arh_factory_oao_cs => &
        standard_arh_factory_oao_cs
    procedure(standard_arh_factory_oao_os), pointer :: arh_factory_oao_os => &
        standard_arh_factory_oao_os
    procedure(standard_arh_deconstructor), pointer :: arh_deconstructor => &
        standard_arh_deconstructor

    ! create function pointers to ensure that routines comply with interface
    procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_f_wrapper_ptr => &
        evaluate_dm_os_f_wrapper
    procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_f_wrapper_ptr => &
        evaluate_dm_cs_f_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_arh_f_c
        module procedure assign_arh_c_f
    end interface

contains

    function arh_factory_mo_c_wrapper(mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, &
                                      n_ao_c, n_mo_c, evaluate_dm_c_funptr, &
                                      update_orbs_arh_c_funptr, obj_func_arh_c_funptr, &
                                      solver_settings_c, settings_c, orbsym_c) &
        result(error_c) bind(C, name="arh_factory_mo")
        !
        ! this function wraps the ARH factory for orbitals parameterized in the MO
        ! basis for C hosts, whose C solver settings it hands to the holder of the ARH
        ! object, and writes the solver settings the factory wired back into the C
        ! solver settings
        !
        use c_interface, only: solver_settings_type_c
        use otr_arh, only: arh_settings_type
        use otr_common_c_interface, only: bundle_settings_c, setup_extension_holder, &
                                          wire_solver_settings_c, obj_func_c_wrapper
        use otr_mo_c_interface, only: mo_coeff_3d_c, update_orbs_mo_c_wrapper

        real(c_rp), intent(inout), target :: mo_coeff_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in) :: n_occ_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c, n_mo_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_arh_c_funptr, obj_func_arh_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(arh_settings_type_c), intent(in) :: settings_c
        integer(c_ip), intent(in), optional :: orbsym_c(*)
        integer(c_ip) :: error_c

        real(rp), pointer, contiguous :: mo_coeff_2d(:, :)
        real(rp), pointer, contiguous :: mo_coeff_3d(:, :, :)
        real(rp), pointer :: ao_overlap(:, :)
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        type(arh_settings_type) :: settings
        type(arh_c_callbacks_type), target :: callbacks
        type(extension_holder_type) :: previous_holder
        type(arh_c_callbacks_type) :: previous_callbacks
        integer(ip) :: n_particle, n_ao, n_mo, error
        integer(ip), allocatable :: n_occ(:), orbsym(:, :), orbsym_1d(:)

        ! no functions are returned unless the factory succeeds
        update_orbs_arh_c_funptr = c_null_funptr
        obj_func_arh_c_funptr = c_null_funptr

        ! convert dimensions to Fortran kind
        n_particle = int(n_particle_c, kind=ip)
        n_ao = int(n_ao_c, kind=ip)
        n_mo = int(n_mo_c, kind=ip)
        n_occ = int(n_occ_c(:n_particle), kind=ip)

        ! convert the irreps of the MOs of every particle channel if they are given and
        ! the dimensions are valid, since the conversion would otherwise read out of
        ! bounds; the sanity check rejects invalid dimensions
        if (present(orbsym_c) .and. n_particle >= 1 .and. n_particle <= 2 .and. &
            n_mo >= 1) then
            orbsym = reshape(int(orbsym_c(:n_mo * n_particle), kind=ip), &
                             [n_mo, n_particle])
            if (n_particle == 1) orbsym_1d = orbsym(:, 1)
        end if

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            if (n_particle == 1) then
                call c_f_pointer(c_loc(mo_coeff_c(1)), mo_coeff_2d, [n_ao, n_mo])
            else
                call c_f_pointer(c_loc(mo_coeff_c(1)), mo_coeff_3d, &
                                 [n_ao, n_mo, n_particle])
            end if
            call c_f_pointer(c_loc(ao_overlap_c(1)), ao_overlap, [n_ao, n_ao])
        else
            if (n_particle == 1) then
                allocate(mo_coeff_2d(n_ao, n_mo))
                mo_coeff_2d = reshape(real(mo_coeff_c(:n_ao * n_mo), kind=rp), &
                                      [n_ao, n_mo])
            else
                allocate(mo_coeff_3d(n_ao, n_mo, n_particle))
                mo_coeff_3d = reshape(real(mo_coeff_c(:n_ao * n_mo * n_particle), &
                                           kind=rp), [n_ao, n_mo, n_particle])
            end if
            allocate(ao_overlap(n_ao, n_ao))
            ao_overlap = reshape(real(ao_overlap_c(:n_ao**2), kind=rp), [n_ao, n_ao])
        end if

        ! associate procedure pointer to wrapper function
        if (n_particle == 1) then
            evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper
        else
            evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper
        end if

        ! convert settings and hand the density matrix evaluating function and the
        ! logging function with the host context of initialized settings to their
        ! wrappers in a bundle, which the settings point to
        settings = settings_c
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=callbacks%evaluate_dm)
        call bundle_settings_c(callbacks, settings_c%initialized, settings_c%logger, &
                               settings_c%context)
        settings%context => callbacks

        ! keep the holder and the bundle held by the module, which are restored if the
        ! factory fails after the holder was set up, so that the functions an earlier
        ! call returned keep working
        previous_holder = arh_holder
        previous_callbacks = arh_callbacks

        ! set up the holder and convert the solver settings, unless the holder was
        ! already wired into these
        call setup_extension_holder(arh_holder, solver_settings_c, settings, error)
        if (error /= 0) then
            error_c = int(error, kind=c_ip)
            return
        end if

        ! the bundle held by the module, in which the wrappers of the callback functions
        ! find them, takes over the bundle only once the holder accepted the solver
        ! settings, so that refused settings leave the callback functions of an object
        ! in use unchanged
        arh_callbacks = callbacks
        settings%context => arh_callbacks

        ! call factory function
        if (n_particle == 1) then
            call arh_factory_mo_cs(mo_coeff_2d, ao_overlap, n_occ(1), n_particle, &
                                   n_ao, n_mo, evaluate_dm_cs_funptr, &
                                   update_orbs_arh_funptr, obj_func_arh_funptr, &
                                   arh_holder%solver_settings%solver_settings_type, &
                                   error, settings, orbsym_1d)
        else
            call arh_factory_mo_os(mo_coeff_3d, ao_overlap, n_occ, n_particle, n_ao, &
                                   n_mo, evaluate_dm_os_funptr, &
                                   update_orbs_arh_funptr, obj_func_arh_funptr, &
                                   arh_holder%solver_settings%solver_settings_type, &
                                   error, settings, orbsym)
        end if

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update of the orbital basis, which copies the orbitals back to the host's
        ! array if the object updates a copy of them, and of the objective function of
        ! the object, or restore the holder and the bundle held by the module if the
        ! factory failed
        if (error == 0) then
            call wire_solver_settings_c(arh_holder, solver_settings_c)
            update_orbs_arh_c_funptr = c_funloc(update_orbs_mo_c_wrapper)
            obj_func_arh_c_funptr = c_funloc(obj_func_c_wrapper)
            if (rp /= c_rp) call c_f_pointer(c_loc(mo_coeff_c(1)), mo_coeff_3d_c, &
                                             [n_ao, n_mo, n_particle])
        else
            arh_holder = previous_holder
            arh_callbacks = previous_callbacks
        end if

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function arh_factory_mo_c_wrapper

    function arh_factory_oao_c_wrapper( &
        dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
        update_orbs_arh_c_funptr, obj_func_arh_c_funptr, solver_settings_c, &
        settings_c) result(error_c) bind(C, name="arh_factory_oao")
        !
        ! this function wraps the ARH factory for orbitals parameterized in the OAO
        ! basis for C hosts, whose C solver settings it hands to the holder of the ARH
        ! object, and writes the solver settings the factory wired back into the C
        ! solver settings
        !
        use c_interface, only: solver_settings_type_c
        use otr_arh, only: arh_settings_type
        use otr_common_c_interface, only: bundle_settings_c, setup_extension_holder, &
                                          wire_solver_settings_c, obj_func_c_wrapper
        use otr_oao_c_interface, only: dm_ao_3d_c, update_orbs_oao_c_wrapper

        real(c_rp), intent(inout), target :: dm_ao_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_arh_c_funptr, obj_func_arh_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(arh_settings_type_c), intent(in) :: settings_c
        integer(c_ip) :: error_c

        real(rp), pointer, contiguous :: dm_ao_2d(:, :)
        real(rp), pointer, contiguous :: dm_ao_3d(:, :, :)
        real(rp), pointer :: ao_overlap(:, :)
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        type(arh_settings_type) :: settings
        type(arh_c_callbacks_type), target :: callbacks
        type(extension_holder_type) :: previous_holder
        type(arh_c_callbacks_type) :: previous_callbacks
        integer(ip) :: n_particle, n_ao, error

        ! no functions are returned unless the factory succeeds
        update_orbs_arh_c_funptr = c_null_funptr
        obj_func_arh_c_funptr = c_null_funptr

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
            evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper
        else
            evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper
        end if

        ! convert settings and hand the density matrix evaluating function and the
        ! logging function with the host context of initialized settings to their
        ! wrappers in a bundle, which the settings point to
        settings = settings_c
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=callbacks%evaluate_dm)
        call bundle_settings_c(callbacks, settings_c%initialized, settings_c%logger, &
                               settings_c%context)
        settings%context => callbacks

        ! keep the holder and the bundle held by the module, which are restored if the
        ! factory fails after the holder was set up, so that the functions an earlier
        ! call returned keep working
        previous_holder = arh_holder
        previous_callbacks = arh_callbacks

        ! set up the holder and convert the solver settings, unless the holder was
        ! already wired into these
        call setup_extension_holder(arh_holder, solver_settings_c, settings, error)
        if (error /= 0) then
            error_c = int(error, kind=c_ip)
            return
        end if

        ! the bundle held by the module, in which the wrappers of the callback functions
        ! find them, takes over the bundle only once the holder accepted the solver
        ! settings, so that refused settings leave the callback functions of an object
        ! in use unchanged
        arh_callbacks = callbacks
        settings%context => arh_callbacks

        ! call factory function
        if (n_particle == 1) then
            call arh_factory_oao_cs( &
                dm_ao_2d, ao_overlap, n_particle, n_ao, evaluate_dm_cs_funptr, &
                update_orbs_arh_funptr, obj_func_arh_funptr, &
                arh_holder%solver_settings%solver_settings_type, error, settings)
        else
            call arh_factory_oao_os( &
                dm_ao_3d, ao_overlap, n_particle, n_ao, evaluate_dm_os_funptr, &
                update_orbs_arh_funptr, obj_func_arh_funptr, &
                arh_holder%solver_settings%solver_settings_type, error, settings)
        end if

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update of the orbital basis, which copies the orbitals back to the host's
        ! array if the object updates a copy of them, and of the objective function of
        ! the object, or restore the holder and the bundle held by the module if the
        ! factory failed
        if (error == 0) then
            call wire_solver_settings_c(arh_holder, solver_settings_c)
            update_orbs_arh_c_funptr = c_funloc(update_orbs_oao_c_wrapper)
            obj_func_arh_c_funptr = c_funloc(obj_func_c_wrapper)
            if (rp /= c_rp) call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_3d_c, &
                                             [n_ao, n_ao, n_particle])
        else
            arh_holder = previous_holder
            arh_callbacks = previous_callbacks
        end if

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function arh_factory_oao_c_wrapper

    subroutine evaluate_dm_cs_f_wrapper(dm, energy, fock, v_coulomb, v_exchange, &
                                        v_nonlinear, error, context)
        !
        ! this subroutine wraps the density matrix evaluating subroutine to convert
        ! Fortran variables to C variables for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target, contiguous :: &
            fock(:, :), v_coulomb(:, :), v_exchange(:, :), v_nonlinear(:, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(rp), pointer :: dm_3d(:, :, :), fock_3d(:, :, :), v_coulomb_3d(:, :, :), &
                             v_exchange_3d(:, :, :), v_nonlinear_3d(:, :, :)

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        nullify(fock_3d, v_coulomb_3d, v_exchange_3d, v_nonlinear_3d)
        if (present(fock)) fock_3d(1:size(fock, 1), 1:size(fock, 2), 1:1) => fock
        if (present(v_coulomb)) &
            v_coulomb_3d(1:size(v_coulomb, 1), 1:size(v_coulomb, 2), 1:1) => v_coulomb
        if (present(v_exchange)) &
            v_exchange_3d(1:size(v_exchange, 1), 1:size(v_exchange, 2), 1:1) => &
            v_exchange
        if (present(v_nonlinear)) &
            v_nonlinear_3d(1:size(v_nonlinear, 1), 1:size(v_nonlinear, 2), 1:1) => &
            v_nonlinear
        call evaluate_dm_os_f_wrapper(dm_3d, energy, fock_3d, v_coulomb_3d, &
                                      v_exchange_3d, v_nonlinear_3d, error, context)

    end subroutine evaluate_dm_cs_f_wrapper

    subroutine evaluate_dm_os_f_wrapper(dm, energy, fock, v_coulomb, v_exchange, &
                                        v_nonlinear, error, context)
        !
        ! this subroutine wraps the density matrix evaluating subroutine to convert
        ! Fortran variables to C variables for the open-shell case
        !
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), optional, target :: &
            fock(:, :, :), v_coulomb(:, :, :), v_exchange(:, :, :), v_nonlinear(:, :, :)
        integer(ip), intent(out) :: error
        class(*), intent(in), pointer :: context

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :, :), fock_c(:, :, :), v_coulomb_c(:, :, :), &
                               v_exchange_c(:, :, :), v_nonlinear_c(:, :, :)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        nullify(fock_c, v_coulomb_c, v_exchange_c, v_nonlinear_c)
        if (rp == c_rp) then
            dm_c => dm
            if (present(fock)) fock_c => fock
            if (present(v_coulomb)) v_coulomb_c => v_coulomb
            if (present(v_exchange)) v_exchange_c => v_exchange
            if (present(v_nonlinear)) v_nonlinear_c => v_nonlinear
        else
            allocate(dm_c, source=real(dm, kind=c_rp))
            if (present(fock)) allocate(fock_c, mold=real(fock, kind=c_rp))
            if (present(v_coulomb)) &
                allocate(v_coulomb_c, mold=real(v_coulomb, kind=c_rp))
            if (present(v_exchange)) &
                allocate(v_exchange_c, mold=real(v_exchange, kind=c_rp))
            if (present(v_nonlinear)) &
                allocate(v_nonlinear_c, mold=real(v_nonlinear, kind=c_rp))
        end if

        ! call density matrix evaluating C function
        error_c = arh_callbacks%evaluate_dm(dm_c, energy_c, fock_c, v_coulomb_c, &
                                            v_exchange_c, v_nonlinear_c, &
                                            arh_callbacks%host_context)

        ! convert arguments to Fortran kind
        energy = real(energy_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            if (present(fock)) then
                fock = real(fock_c, kind=rp)
                deallocate(fock_c)
            end if
            if (present(v_coulomb)) then
                v_coulomb = real(v_coulomb_c, kind=rp)
                deallocate(v_coulomb_c)
            end if
            if (present(v_exchange)) then
                v_exchange = real(v_exchange_c, kind=rp)
                deallocate(v_exchange_c)
            end if
            if (present(v_nonlinear)) then
                v_nonlinear = real(v_nonlinear_c, kind=rp)
                deallocate(v_nonlinear_c)
            end if
            deallocate(dm_c)
        end if

    end subroutine evaluate_dm_os_f_wrapper

    subroutine init_arh_settings_c(settings_c) bind(C, name="init_arh_settings")
        !
        ! this subroutine initializes the ARH solver settings
        !
        use otr_arh, only: default_arh_settings

        type(arh_settings_type_c), intent(inout) :: settings_c

        settings_c = default_arh_settings

    end subroutine init_arh_settings_c

    subroutine arh_deconstructor_c_wrapper() bind(C, name="arh_deconstructor")
        !
        ! this subroutine deallocates the ARH objects
        !
        call arh_deconstructor()

    end subroutine arh_deconstructor_c_wrapper

    subroutine assign_arh_f_c(settings, settings_c)
        !
        ! this subroutine converts ARH settings from C to Fortran
        !
        use otr_arh, only: arh_settings_type
        use c_interface, only: character_from_c, logger_f_wrapper

        type(arh_settings_type), intent(out) :: settings
        type(arh_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%logger)) then
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert integers
            settings%verbose = int(settings_c%verbose, kind=ip)

            ! convert characters
            settings%arh_type = character_from_c(settings_c%arh_type)

            ! set settings to initialized
            settings%initialized = .true.
        end if

    end subroutine assign_arh_f_c

    subroutine assign_arh_c_f(settings_c, settings)
        !
        ! this subroutine converts ARH settings from Fortran to C
        !
        use otr_arh, only: arh_settings_type
        use c_interface, only: character_to_c

        type(arh_settings_type_c), intent(out) :: settings_c
        type(arh_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%logger = c_null_funptr

            ! convert integers
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! convert characters
            settings_c%arh_type = character_to_c(settings%arh_type)

            ! the host context cannot be converted
            settings_c%context = c_null_ptr

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_arh_c_f

end module otr_arh_c_interface
