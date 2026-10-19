! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_c_interface

    use opentrustregion, only: ip, rp, obj_func_type, update_orbs_type
    use c_interface, only: c_ip, c_rp
    use otr_mo, only: standard_mo_factory_cs => mo_factory_cs, &
                      standard_mo_factory_os => mo_factory_os
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_loc, &
                                           c_f_pointer, c_f_procpointer, c_associated, &
                                           c_null_funptr, c_null_ptr

    implicit none

    ! derived type for MO settings
    type, bind(C) :: mo_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose
        type(c_ptr) :: context
    end type

    procedure(standard_mo_factory_cs), pointer :: mo_factory_cs => &
        standard_mo_factory_cs
    procedure(standard_mo_factory_os), pointer :: mo_factory_os => &
        standard_mo_factory_os

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_mo_f_c
        module procedure assign_mo_c_f
    end interface

contains

    function mo_factory_c_wrapper( &
        mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
        evaluate_dm_c_funptr, update_orbs_mo_c_funptr, obj_func_mo_c_funptr, &
        extension_c, solver_settings_c, solver_state_c, settings_c, orbsym_c) &
        result(error_c) bind(C, name="mo_factory")
        !
        ! this function wraps the MO factory for C hosts, which creates a holder of the
        ! new MO object, returned as the extension the host releases, and hands it the
        ! host's C solver settings and record of the solver state, and writes the
        ! solver settings the factory wired back into the C solver settings
        !
        use opentrustregion, only: warn_uninitialized_settings, add_error_origin
        use c_interface, only: solver_settings_type_c
        use otr_mo, only: mo_settings_type, default_mo_settings, &
                          mo_settings_uninitialized_warning_msg
        use otr_common, only: evaluate_dm_cs_type, evaluate_dm_os_type, error_mo
        use otr_common_c_interface, only: &
            orbital_c_callbacks_type, extension_holder_type, bundle_settings_c, &
            setup_extension_holder, finish_extension_holder, evaluate_dm_cs_f_wrapper, &
            evaluate_dm_os_f_wrapper

        real(c_rp), intent(inout), target :: mo_coeff_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in) :: n_occ_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c, n_mo_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_mo_c_funptr, obj_func_mo_c_funptr
        type(c_ptr), intent(out) :: extension_c
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_ptr), intent(in), value :: solver_state_c
        type(mo_settings_type_c), intent(in), optional :: settings_c
        integer(c_ip), intent(in), optional :: orbsym_c(*)
        integer(c_ip) :: error_c

        real(rp), pointer, contiguous :: mo_coeff_2d(:, :), mo_coeff_3d(:, :, :)
        real(rp), pointer :: ao_overlap(:, :)
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
        procedure(obj_func_type), pointer :: obj_func_mo_funptr
        type(mo_settings_type) :: settings
        type(mo_settings_type_c) :: local_settings_c
        type(extension_holder_type), pointer :: holder
        integer(ip) :: n_particle, n_ao, n_mo, error
        integer(ip), allocatable :: n_occ(:), orbsym(:, :), orbsym_1d(:)

        ! no functions or extension are returned unless the factory succeeds
        update_orbs_mo_c_funptr = c_null_funptr
        obj_func_mo_c_funptr = c_null_funptr
        extension_c = c_null_ptr

        ! convert dimensions to Fortran kind and the occupations of every particle
        ! channel if the number of particle channels is valid, since the conversion
        ! would otherwise read out of bounds; the sanity check rejects an invalid number
        n_particle = int(n_particle_c, kind=ip)
        n_ao = int(n_ao_c, kind=ip)
        n_mo = int(n_mo_c, kind=ip)
        if (n_particle >= 1 .and. n_particle <= 2) then
            n_occ = int(n_occ_c(:n_particle), kind=ip)
        else
            allocate(n_occ(0))
        end if

        ! work on a copy of the given settings or, without any, of the default settings
        if (present(settings_c)) then
            local_settings_c = settings_c
        else
            local_settings_c = default_mo_settings
        end if

        ! convert the settings, which replaces settings that were not initialized by the
        ! default settings, and print the warning for these on stderr, since they
        ! provide no logging function
        settings = local_settings_c
        if (.not. local_settings_c%initialized) call warn_uninitialized_settings( &
            settings, mo_settings_uninitialized_warning_msg)

        ! convert the irreps of the MOs of every particle channel if they are given and
        ! the dimensions are valid, since the conversion would otherwise read out of
        ! bounds; the sanity check rejects invalid dimensions
        if (present(orbsym_c) .and. n_particle >= 1 .and. n_particle <= 2 .and. &
            n_mo >= 1) then
            orbsym = reshape(int(orbsym_c(:n_mo * n_particle), kind=ip), &
                             [n_mo, n_particle])
            if (n_particle == 1) orbsym_1d = orbsym(:, 1)
        end if

        ! the object rotates the host's MO coefficients in place or, if the real kind
        ! differs from C's, a copy of them, which the orbital update copies back
        if (rp == c_rp) then
            call c_f_pointer(c_loc(mo_coeff_c(1)), mo_coeff_3d, &
                             [n_ao, n_mo, n_particle])
            call c_f_pointer(c_loc(ao_overlap_c(1)), ao_overlap, [n_ao, n_ao])
        else
            allocate(mo_coeff_3d(n_ao, n_mo, n_particle))
            mo_coeff_3d = reshape(real(mo_coeff_c(:n_ao * n_mo * n_particle), &
                                       kind=rp), [n_ao, n_mo, n_particle])
            allocate(ao_overlap(n_ao, n_ao))
            ao_overlap = reshape(real(ao_overlap_c(:n_ao**2), kind=rp), [n_ao, n_ao])
        end if

        ! create the holder, whose bundle hands the density matrix evaluating function,
        ! if one is given, and the logging function with the host context of
        ! initialized settings to their wrappers, which the settings point to
        allocate(holder)
        allocate(orbital_c_callbacks_type :: holder%extension_callbacks)
        select type (callbacks => holder%extension_callbacks)
        type is (orbital_c_callbacks_type)
            if (c_associated(evaluate_dm_c_funptr)) call c_f_procpointer( &
                cptr=evaluate_dm_c_funptr, fptr=callbacks%evaluate_dm)
            call bundle_settings_c(callbacks, local_settings_c%initialized, &
                                   local_settings_c%logger, local_settings_c%context)
        end select
        settings%context => holder%extension_callbacks

        ! set up the holder with the host's record of the state and convert the solver
        ! settings
        call setup_extension_holder(holder, solver_settings_c, solver_state_c, &
                                    settings, error)
        call add_error_origin(error, error_mo, settings)

        ! the factory wires the solver settings of the holder, which then own the
        ! object, reaches the density matrix evaluating function, if one is given,
        ! through the shared wrappers, which find it in the bundle of the holder, and
        ! reads the state of the holder
        if (error == 0) then
            evaluate_dm_cs_funptr => null()
            evaluate_dm_os_funptr => null()
            if (n_particle == 1) then
                if (c_associated(evaluate_dm_c_funptr)) &
                    evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper
                mo_coeff_2d => mo_coeff_3d(:, :, 1)
                call mo_factory_cs( &
                    mo_coeff_2d, ao_overlap, n_occ(1), n_particle, n_ao, n_mo, &
                    evaluate_dm_cs_funptr, update_orbs_mo_funptr, obj_func_mo_funptr, &
                    holder%solver_settings, holder%state, error, settings, orbsym_1d)
            else
                if (c_associated(evaluate_dm_c_funptr)) &
                    evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper
                call mo_factory_os(mo_coeff_3d, ao_overlap, n_occ, n_particle, n_ao, &
                                   n_mo, evaluate_dm_os_funptr, update_orbs_mo_funptr, &
                                   obj_func_mo_funptr, holder%solver_settings, &
                                   holder%state, error, settings, orbsym)
            end if
        end if

        ! the orbital update copies the MO coefficients back to the host's array if the
        ! object rotates a copy of them
        if (error == 0 .and. rp /= c_rp) then
            holder%solver_level%orbitals => mo_coeff_3d
            call c_f_pointer(c_loc(mo_coeff_c(1)), holder%solver_level%orbitals_c, &
                             [n_ao, n_mo, n_particle])
        end if

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update and the objective function of the object and the holder, or free the
        ! holder if the factory failed
        call finish_extension_holder(holder, error, solver_settings_c, &
                                     update_orbs_mo_c_funptr, obj_func_mo_c_funptr, &
                                     extension_c)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function mo_factory_c_wrapper

    subroutine init_mo_settings_c(settings_c) bind(C, name="init_mo_settings")
        !
        ! this subroutine initializes the MO settings
        !
        use otr_mo, only: default_mo_settings

        type(mo_settings_type_c), intent(inout) :: settings_c

        settings_c = default_mo_settings

    end subroutine init_mo_settings_c

    subroutine assign_mo_f_c(settings, settings_c)
        !
        ! this subroutine converts MO settings from C to Fortran
        !
        use otr_mo, only: mo_settings_type, default_mo_settings
        use c_interface, only: logger_f_wrapper

        type(mo_settings_type), intent(out) :: settings
        type(mo_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%logger)) then
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert integers
            settings%verbose = int(settings_c%verbose, kind=ip)
        else
            ! C settings that were not initialized hold no values and are converted
            ! as the default settings
            settings = default_mo_settings
        end if

    end subroutine assign_mo_f_c

    subroutine assign_mo_c_f(settings_c, settings)
        !
        ! this subroutine converts MO settings from Fortran to C
        !
        use otr_mo, only: mo_settings_type

        type(mo_settings_type_c), intent(out) :: settings_c
        type(mo_settings_type), intent(in) :: settings

        ! callback functions cannot be converted
        settings_c%logger = c_null_funptr

        ! convert integers
        settings_c%verbose = int(settings%verbose, kind=c_ip)

        ! the host context cannot be converted
        settings_c%context = c_null_ptr

        ! Fortran settings are always initialized
        settings_c%initialized = .true._c_bool

    end subroutine assign_mo_c_f

end module otr_mo_c_interface
