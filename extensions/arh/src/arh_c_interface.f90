! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface

    use opentrustregion, only: ip, rp, kw_len, update_orbs_type, hess_x_type, &
                               project_type
    use c_interface, only: c_ip, c_rp, update_orbs_c_type, hess_x_c_type
    use otr_oao_c_interface, only: n_particle, n_ao
    use otr_arh, only: standard_arh_factory_closed_shell => arh_factory_closed_shell, &
                       standard_arh_factory_open_shell => arh_factory_open_shell, &
                       standard_arh_deconstructor => arh_deconstructor, &
                       update_dm_spin_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_loc, c_f_pointer, &
                                           c_funloc, c_f_procpointer, c_associated, &
                                           c_char, c_null_funptr

    implicit none

    ! define procedure pointer which will point to the Fortran procedures
    procedure(update_dm_spin_c_type), pointer :: update_dm_spin_before_wrapping => &
                                                 null()
    procedure(update_orbs_type), pointer :: update_orbs_arh_before_wrapping => null()
    procedure(hess_x_type), pointer :: hess_x_arh_before_wrapping => null()

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function update_dm_spin_c_type(dm_ao_c, energy_c, fock_c, v_same_spin_c, &
                                       v_opposite_spin_c, get_response_c_funptr) &
            result(error_c) bind(C)
            import :: c_rp, c_ip, c_funptr
    
            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            real(c_rp), intent(out), target :: fock_c(*), v_same_spin_c(*), &
                                               v_opposite_spin_c(*)
            type(c_funptr), intent(out) :: get_response_c_funptr
            integer(c_ip) :: error_c
        end function update_dm_spin_c_type
    end interface

    ! derived type for ARH settings
    type, bind(C) :: arh_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized, restricted
        integer(c_ip) :: verbose
        character(c_char) :: arh_type(kw_len + 1)
    end type

    procedure(standard_arh_factory_closed_shell), pointer :: arh_factory_closed_shell &
        => standard_arh_factory_closed_shell
    procedure(standard_arh_factory_open_shell), pointer :: arh_factory_open_shell &
        => standard_arh_factory_open_shell
    procedure(standard_arh_deconstructor), pointer :: arh_deconstructor => &
        standard_arh_deconstructor

    ! create function pointers to ensure that routines comply with interface
    procedure(update_dm_spin_type), pointer :: update_dm_spin_f_wrapper_ptr => &
        update_dm_spin_f_wrapper
    procedure(update_orbs_c_type), pointer :: update_orbs_arh_c_wrapper_ptr => &
        update_orbs_arh_c_wrapper
    procedure(hess_x_c_type), pointer :: hess_x_arh_c_wrapper_ptr => &
        hess_x_arh_c_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_arh_f_c
        module procedure assign_arh_c_f
    end interface

contains

    function arh_factory_c_wrapper(dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, &
                                   get_energy_c_funptr, update_dm_c_funptr, &
                                   obj_func_arh_c_funptr, update_orbs_arh_c_funptr, &
                                   project_arh_c_funptr, settings_c) result(error_c) &
        bind(C, name="arh_factory")
        !
        ! this subroutine wraps the factory function for the subroutine to convert C 
        ! variables to Fortran variables
        !
        use otr_arh, only: arh_settings_type
        use otr_oao, only: get_energy_2d_type, get_energy_3d_type, update_dm_2d_type, &
                           obj_func_type
        use otr_oao_c_interface, only: dm_ao_3d_c, get_energy_before_wrapping, &
                                       update_dm_before_wrapping, &
                                       get_energy_2d_f_wrapper, &
                                       update_dm_2d_f_wrapper, &
                                       get_energy_3d_f_wrapper, &
                                       obj_func_oao_before_wrapping, &
                                       project_oao_before_wrapping, &
                                       obj_func_oao_c_wrapper, project_oao_c_wrapper
        use otr_common_c_interface, only: n_param

        real(c_rp), intent(in), target :: dm_ao_c(*), ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: get_energy_c_funptr, update_dm_c_funptr
        type(arh_settings_type_c), intent(inout) :: settings_c
        type(c_funptr), intent(out) :: obj_func_arh_c_funptr, &
                                       update_orbs_arh_c_funptr, project_arh_c_funptr
        integer(c_ip) :: error_c

        real(rp), pointer, contiguous :: dm_ao_2d(:, :)
        real(rp), pointer :: dm_ao_3d(:, :, :), ao_overlap(:, :)
        procedure(get_energy_2d_type), pointer :: get_energy_2d_funptr
        procedure(get_energy_3d_type), pointer :: get_energy_3d_funptr
        procedure(update_dm_2d_type), pointer :: update_dm_funptr
        procedure(update_dm_spin_type), pointer :: update_dm_spin_funptr
        procedure(obj_func_type), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), pointer :: update_orbs_arh_funptr
        procedure(project_type), pointer :: project_arh_funptr
        type(arh_settings_type) :: settings
        integer(ip) :: error

        ! convert number of AOs to Fortran kind, calculate number of parameters and 
        ! store globally to access assumed size arrays passed from C to Fortran
        n_particle = int(n_particle_c, kind=ip)
        n_ao = int(n_ao_c, kind=ip)
        n_param = n_ao * (n_ao - 1) / 2
        if (.not. settings_c%restricted) n_param = n_particle * n_param

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            if (n_particle == 1) then
                call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_2d, [n_ao, n_ao])
            else
                call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_3d, [n_ao, n_ao, n_particle])
            end if
            call c_f_pointer(c_loc(ao_overlap_c(1)), ao_overlap, [n_ao, n_ao])
        else
            call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_3d_c, [n_ao, n_ao, n_particle])
            if (n_particle == 1) then
                allocate(dm_ao_2d(n_ao, n_ao))
                dm_ao_2d = reshape(real(dm_ao_c(:n_ao_c ** 2), kind=rp), [n_ao, n_ao])
            else
                allocate(dm_ao_3d(n_ao, n_ao, n_particle))
                dm_ao_3d = reshape(real(dm_ao_c(:(n_ao_c ** 2 * n_particle_c)), &
                                        kind=rp), [n_ao, n_ao, n_particle])
            end if
            allocate(ao_overlap(n_ao, n_ao))
            ao_overlap = reshape(real(ao_overlap_c(:n_ao ** 2), kind=rp), [n_ao, n_ao])
        end if

        ! associate the input C pointers to Fortran procedure pointers
        call c_f_procpointer(cptr=get_energy_c_funptr, fptr=get_energy_before_wrapping)
        if (n_particle == 1) then
            call c_f_procpointer(cptr=update_dm_c_funptr, &
                                 fptr=update_dm_before_wrapping)
        else
            call c_f_procpointer(cptr=update_dm_c_funptr, &
                                 fptr=update_dm_spin_before_wrapping)
        end if

        ! associate procedure pointer to wrapper function
        if (n_particle == 1) then
            get_energy_2d_funptr => get_energy_2d_f_wrapper
            update_dm_funptr => update_dm_2d_f_wrapper
        else
            get_energy_3d_funptr => get_energy_3d_f_wrapper
            update_dm_spin_funptr => update_dm_spin_f_wrapper
        end if

        ! convert settings
        settings = settings_c

        ! call factory function
        if (n_particle == 1) then
            call arh_factory_closed_shell(dm_ao_2d, ao_overlap, n_particle, n_ao, &
                                          get_energy_2d_funptr, update_dm_funptr, &
                                          obj_func_arh_funptr, update_orbs_arh_funptr, &
                                          project_arh_funptr, error, settings)
        else
            call arh_factory_open_shell(dm_ao_3d, ao_overlap, n_particle, n_ao, &
                                        get_energy_3d_funptr, update_dm_spin_funptr, &
                                        obj_func_arh_funptr, update_orbs_arh_funptr, &
                                        project_arh_funptr, error, settings)
        end if

        ! associate the global procedure pointers to the Fortran function pointers
        obj_func_oao_before_wrapping => obj_func_arh_funptr
        update_orbs_arh_before_wrapping => update_orbs_arh_funptr
        project_oao_before_wrapping => project_arh_funptr

        ! get a C function pointer to the C wrapper functions
        obj_func_arh_c_funptr = c_funloc(obj_func_oao_c_wrapper)
        update_orbs_arh_c_funptr = c_funloc(update_orbs_arh_c_wrapper)
        project_arh_c_funptr = c_funloc(project_oao_c_wrapper)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function arh_factory_c_wrapper

    subroutine update_dm_spin_f_wrapper(dm, energy, fock, v_same_spin, &
                                        v_opposite_spin, get_response_funptr, error)
        !
        ! this subroutine wraps the density matrix updating subroutine to convert 
        ! Fortran variables to C variables
        !
        use otr_oao, only: get_response_3d_type
        use otr_oao_c_interface, only: get_response_before_wrapping, &
                                       get_response_3d_f_wrapper

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target :: fock(:, :, :), v_same_spin(:, :, :), &
                                         v_opposite_spin(:, :, :)
        procedure(get_response_3d_type), intent(out), pointer :: get_response_funptr
        integer(ip), intent(out) :: error

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :, :), fock_c(:, :, :), v_same_spin_c(:, :, :), &
                               v_opposite_spin_c(:, :, :)
        type(c_funptr) :: get_response_c_funptr
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            dm_c => dm
            fock_c => fock
            v_same_spin_c => v_same_spin
            v_opposite_spin_c => v_opposite_spin
        else
            allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)), &
                     fock_c(size(dm, 1), size(dm, 2), size(dm, 3)), &
                     v_same_spin_c(size(dm, 1), size(dm, 2), size(dm, 3)), &
                     v_opposite_spin_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            dm_c = real(dm, kind=c_rp)
        end if

        ! call density matrix updating C function
        error_c = update_dm_spin_before_wrapping(dm_c, energy_c, fock_c, v_same_spin_c, &
                                                 v_opposite_spin_c, &
                                                 get_response_c_funptr)

        ! convert arguments to Fortran kind
        energy = real(energy_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            fock = real(fock_c, kind=rp)
            v_same_spin = real(v_same_spin_c, kind=rp)
            v_opposite_spin = real(v_opposite_spin_c, kind=rp)
            deallocate(dm_c, fock_c, v_same_spin_c, v_opposite_spin_c)
        end if

        ! associate the input C pointer to get_response function to a Fortran procedure
        ! pointer
        call c_f_procpointer(cptr=get_response_c_funptr, &
                             fptr=get_response_before_wrapping)

        ! associate procedure pointer to wrapper function
        get_response_funptr => get_response_3d_f_wrapper

    end subroutine update_dm_spin_f_wrapper

    function update_orbs_arh_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                       hess_x_c_funptr) result(error_c) bind(C)
        !
        ! this function wraps the orbital update subroutine to convert Fortran 
        ! variables to C variables
        !
        use otr_common_c_interface, only: update_orbs_c_wrapper_impl
        use otr_oao_c_interface, only: dm_ao_3d_c
        use otr_oao, only: oao_object

        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(out) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        error_c = update_orbs_c_wrapper_impl(update_orbs_arh_before_wrapping, &
                                             hess_x_arh_before_wrapping, &
                                             hess_x_arh_c_wrapper, kappa_c, func_c, &
                                             grad_c, h_diag_c, hess_x_c_funptr)

        if (rp /= c_rp) dm_ao_3d_c = real(oao_object%dm_ao, kind=c_rp)

    end function update_orbs_arh_c_wrapper

    function hess_x_arh_c_wrapper(x_c, hess_x_c) result(error_c) bind(C)
        !
        ! this function wraps the Hessian linear transformation to convert Fortran 
        ! variables to C variables
        !
        use otr_common_c_interface, only: hess_x_c_wrapper_impl

        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        integer(c_ip) :: error_c

        error_c = hess_x_c_wrapper_impl(hess_x_arh_before_wrapping, x_c, hess_x_c)

    end function hess_x_arh_c_wrapper

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
        use c_interface, only: logger_before_wrapping, logger_f_wrapper, &
                               character_from_c

        type(arh_settings_type), intent(out) :: settings
        type(arh_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%logger)) then
                call c_f_procpointer(cptr=settings_c%logger, &
                                     fptr=logger_before_wrapping)
                settings%logger => logger_f_wrapper
            else
                settings%logger => null()
            end if

            ! convert logicals
            settings%restricted = logical(settings_c%restricted)

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

            ! convert logicals
            settings_c%restricted = logical(settings%restricted, kind=c_bool)

            ! convert integers
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! convert characters
            settings_c%arh_type = character_to_c(settings%arh_type)

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_arh_c_f

end module otr_arh_c_interface
