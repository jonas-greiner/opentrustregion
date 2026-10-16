! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface

    use opentrustregion, only: ip, rp, kw_len, obj_func_type, update_orbs_type
    use c_interface, only: c_ip, c_rp, c_callbacks_type
    use otr_arh, only: standard_arh_factory_common => arh_factory_common, &
                       arh_evaluate_dm_os_type, arh_evaluate_dm_cs_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_ptr, c_funptr, c_f_procpointer, &
                                           c_associated, c_char, c_null_funptr, &
                                           c_null_ptr

    implicit none

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function arh_evaluate_dm_c_type(dm_ao_c, energy_c, fock_c, v_coulomb_c, &
                                        v_exchange_c, v_nonlinear_c, context_c) &
            result(error_c) bind(C)
            import :: c_rp, c_ip, c_ptr

            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            real(c_rp), intent(out), optional :: fock_c(*), v_coulomb_c(*), &
                                                 v_exchange_c(*), v_nonlinear_c(*)
            type(c_ptr), intent(in), value :: context_c
            integer(c_ip) :: error_c
        end function arh_evaluate_dm_c_type
    end interface

    ! bundle of the C callback functions of ARH, which extends that of the core by the
    ! density matrix evaluating function, held by the holder of the object, since the
    ! settings outlive the factory call, and reached through the context of the
    ! settings
    type, extends(c_callbacks_type) :: arh_c_callbacks_type
        procedure(arh_evaluate_dm_c_type), pointer, nopass :: evaluate_dm => null()
    end type

    ! derived type for ARH settings
    type, bind(C) :: arh_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose
        character(kind=c_char) :: arh_type(kw_len + 1)
        type(c_ptr) :: context
    end type

    procedure(standard_arh_factory_common), pointer :: arh_factory_common => &
        standard_arh_factory_common

    ! create function pointers to ensure that routines comply with interface
    procedure(arh_evaluate_dm_cs_type), pointer :: arh_evaluate_dm_cs_f_wrapper_ptr => &
        arh_evaluate_dm_cs_f_wrapper
    procedure(arh_evaluate_dm_os_type), pointer :: arh_evaluate_dm_os_f_wrapper_ptr => &
        arh_evaluate_dm_os_f_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_arh_f_c
        module procedure assign_arh_c_f
    end interface

contains

    function arh_factory_c_wrapper(evaluate_dm_c_funptr, update_orbs_arh_c_funptr, &
                                   obj_func_arh_c_funptr, extension_c, &
                                   solver_settings_c, solver_state_c, settings_c) &
        result(error_c) bind(C, name="arh_factory")
        !
        ! this function wraps the ARH factory for C hosts, which creates a holder of the
        ! new ARH object on the orbital basis object of the holder whose level is the
        ! context of the C solver settings, which the MO or OAO factory wired last, and
        ! hands it the host's C solver settings and record of the solver state, and
        ! writes the solver settings the factory wired back into the C solver settings,
        ! returning the holder as the extension the host releases before that of the
        ! orbital basis object
        !
        use opentrustregion, only: warn_uninitialized_settings, add_error_origin
        use c_interface, only: solver_settings_type_c
        use otr_arh, only: arh_settings_type, default_arh_settings, &
                           arh_settings_uninitialized_warning_msg
        use otr_common, only: extension_level_type, orbital_basis_type, error_arh
        use otr_common_c_interface, only: &
            extension_holder_type, holder_level_type, holder_level, bundle_settings_c, &
            setup_extension_holder, finish_extension_holder

        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(c_funptr), intent(out) :: update_orbs_arh_c_funptr, obj_func_arh_c_funptr
        type(c_ptr), intent(out) :: extension_c
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(c_ptr), intent(in), value :: solver_state_c
        type(arh_settings_type_c), intent(in), optional :: settings_c
        integer(c_ip) :: error_c

        procedure(arh_evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(arh_evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        type(arh_settings_type) :: settings
        type(arh_settings_type_c) :: local_settings_c
        type(extension_holder_type), pointer :: holder
        type(holder_level_type), pointer :: orbitals_held
        type(extension_level_type), pointer :: orbitals_level
        integer(ip) :: n_particle, error

        ! no functions or extension are returned unless the factory succeeds
        update_orbs_arh_c_funptr = c_null_funptr
        obj_func_arh_c_funptr = c_null_funptr
        extension_c = c_null_ptr

        ! work on a copy of the given settings or, without any, of the default settings
        if (present(settings_c)) then
            local_settings_c = settings_c
        else
            local_settings_c = default_arh_settings
        end if

        ! convert the settings and initialize them if they were not, after printing a
        ! warning, as the factory would, which would otherwise discard the context
        ! through which the wrappers of the host's callback functions find them
        settings = local_settings_c
        if (.not. settings%initialized) then
            call warn_uninitialized_settings(settings, &
                                             arh_settings_uninitialized_warning_msg)
            call settings%init(error)
        end if

        ! the context of the C solver settings is trusted to be the level of the holder
        ! of an extension object, whose object the Fortran factory checks to be an
        ! orbital basis object, whose number of particles selects the wrapper of the
        ! density matrix evaluating function
        orbitals_held => holder_level(solver_settings_c%context)
        orbitals_level => null()
        if (associated(orbitals_held)) orbitals_level => orbitals_held%level
        n_particle = 0
        if (associated(orbitals_level)) then
            select type (orbitals => orbitals_level%extension)
            class is (orbital_basis_type)
                n_particle = orbitals%n_particle
            end select
        end if

        ! create the holder, whose bundle hands the density matrix evaluating function
        ! and the logging function with the host context of initialized settings to
        ! their wrappers, which the settings point to
        allocate(holder)
        allocate(arh_c_callbacks_type :: holder%extension_callbacks)
        select type (callbacks => holder%extension_callbacks)
        type is (arh_c_callbacks_type)
            call c_f_procpointer(cptr=evaluate_dm_c_funptr, fptr=callbacks%evaluate_dm)
            call bundle_settings_c(callbacks, local_settings_c%initialized, &
                                   local_settings_c%logger, local_settings_c%context)
        end select
        settings%context => holder%extension_callbacks

        ! set up the holder with the host's record of the state and convert the solver
        ! settings
        call setup_extension_holder(holder, solver_settings_c, solver_state_c, &
                                    settings, error)
        call add_error_origin(error, error_arh, settings)

        ! the factory wires the solver settings of the holder, which then own the
        ! object, on the orbital basis object of the level, reaches the density matrix
        ! evaluating function through the wrapper of the shell of the orbitals, which
        ! finds it in the bundle of the holder, and reads the state of the holder
        if (error == 0) then
            evaluate_dm_cs_funptr => null()
            evaluate_dm_os_funptr => null()
            if (n_particle == 1) then
                evaluate_dm_cs_funptr => arh_evaluate_dm_cs_f_wrapper
            else
                evaluate_dm_os_funptr => arh_evaluate_dm_os_f_wrapper
            end if
            call arh_factory_common(orbitals_level, evaluate_dm_cs_funptr, &
                                    evaluate_dm_os_funptr, update_orbs_arh_funptr, &
                                    obj_func_arh_funptr, holder%solver_settings, &
                                    holder%state, error, settings)
        end if

        ! the orbital update copies the orbitals back to the host's array if the object
        ! moves a copy of them, as that of the orbital basis object does
        if (error == 0 .and. associated(orbitals_held)) then
            holder%solver_level%orbitals => orbitals_held%orbitals
            holder%solver_level%orbitals_c => orbitals_held%orbitals_c
        end if

        ! write the wired solver settings back and return the C wrappers of the orbital
        ! update and the objective function of the object and the holder, or free the
        ! holder if the factory failed
        call finish_extension_holder(holder, error, solver_settings_c, &
                                     update_orbs_arh_c_funptr, obj_func_arh_c_funptr, &
                                     extension_c)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function arh_factory_c_wrapper

    subroutine arh_evaluate_dm_cs_f_wrapper(dm, energy, fock, v_coulomb, v_exchange, &
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
        call arh_evaluate_dm_os_f_wrapper(dm_3d, energy, fock_3d, v_coulomb_3d, &
                                          v_exchange_3d, v_nonlinear_3d, error, context)

    end subroutine arh_evaluate_dm_cs_f_wrapper

    subroutine arh_evaluate_dm_os_f_wrapper(dm, energy, fock, v_coulomb, v_exchange, &
                                            v_nonlinear, error, context)
        !
        ! this subroutine wraps the density matrix evaluating subroutine to convert
        ! Fortran variables to C variables for the open-shell case, finding the C
        ! function and the host context in the bundle it receives as context
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

        ! the bundle holds the C function
        energy = 0.0_rp
        error = 1
        if (.not. associated(context)) return
        select type (callbacks => context)
        type is (arh_c_callbacks_type)
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
            error_c = callbacks%evaluate_dm(dm_c, energy_c, fock_c, v_coulomb_c, &
                                            v_exchange_c, v_nonlinear_c, &
                                            callbacks%host_context)

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
        end select

    end subroutine arh_evaluate_dm_os_f_wrapper

    subroutine init_arh_settings_c(settings_c) bind(C, name="init_arh_settings")
        !
        ! this subroutine initializes the ARH solver settings
        !
        use otr_arh, only: default_arh_settings

        type(arh_settings_type_c), intent(inout) :: settings_c

        settings_c = default_arh_settings

    end subroutine init_arh_settings_c

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
