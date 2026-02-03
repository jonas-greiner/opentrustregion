! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface

    use opentrustregion, only: ip, rp, obj_func_type, update_orbs_type, hess_x_type, &
                               project_type
    use c_interface, only: c_ip, c_rp, obj_func_c_type, update_orbs_c_type, &
                           hess_x_c_type, project_c_type
    use otr_common_c_interface, only: n_param
    use otr_arh, only: standard_arh_factory_closed_shell => arh_factory_closed_shell, &
                       standard_arh_factory_open_shell => arh_factory_open_shell, &
                       standard_arh_deconstructor_closed_shell => &
                       arh_deconstructor_closed_shell, &
                       standard_arh_deconstructor_open_shell => &
                       arh_deconstructor_open_shell, get_energy_closed_shell_type, &
                       get_energy_open_shell_type, get_fock_type, &
                       get_fock_jk_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_loc, c_f_pointer, &
                                           c_funloc, c_f_procpointer, c_associated, &
                                           c_null_funptr

    implicit none

    ! define procedure pointer which will point to the Fortran procedures
    procedure(get_energy_c_type), pointer :: get_energy_before_wrapping => null()
    procedure(get_fock_c_type), pointer :: get_fock_before_wrapping => null()
    procedure(get_fock_jk_c_type), pointer :: get_fock_jk_before_wrapping => null()
    procedure(obj_func_type), pointer :: obj_func_arh_before_wrapping => null()
    procedure(update_orbs_type), pointer :: update_orbs_arh_before_wrapping => null()
    procedure(hess_x_type), pointer :: hess_x_arh_before_wrapping => null()
    procedure(project_type), pointer :: project_arh_before_wrapping => null()

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function get_energy_c_type(dm_ao_c, energy_c) result(error_c) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            integer(c_ip) :: error_c
        end function get_energy_c_type
    end interface

    abstract interface
        function get_fock_c_type(dm_ao_c, energy_c, fock_c) result(error_c) bind(C)
            import :: c_rp, c_ip
    
            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            real(c_rp), intent(out), target :: fock_c(*)
            integer(c_ip) :: error_c
        end function get_fock_c_type
    end interface

    abstract interface
        function get_fock_jk_c_type(dm_ao_c, energy_c, fock_c, coulomb_c, exchange_c) &
            result(error_c) bind(C)
            import :: c_rp, c_ip
    
            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            real(c_rp), intent(out), target :: fock_c(*), coulomb_c(*), exchange_c(*)
            integer(c_ip) :: error_c
        end function get_fock_jk_c_type
    end interface

    ! derived type for ARH settings
    type, bind(C) :: arh_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized, restricted
        integer(c_ip) :: verbose
    end type

    ! global variables
    integer(ip) :: n_particle, n_ao

    procedure(standard_arh_factory_closed_shell), pointer :: arh_factory_closed_shell &
        => standard_arh_factory_closed_shell
    procedure(standard_arh_factory_open_shell), pointer :: arh_factory_open_shell &
        => standard_arh_factory_open_shell
    procedure(standard_arh_deconstructor_closed_shell), pointer :: &
        arh_deconstructor_closed_shell => standard_arh_deconstructor_closed_shell
    procedure(standard_arh_deconstructor_open_shell), pointer :: &
        arh_deconstructor_open_shell => standard_arh_deconstructor_open_shell

    ! create function pointers to ensure that routines comply with interface
    procedure(get_energy_closed_shell_type), pointer :: &
        get_energy_closed_shell_f_wrapper_ptr => get_energy_closed_shell_f_wrapper
    procedure(get_energy_open_shell_type), pointer :: &
        get_energy_open_shell_f_wrapper_ptr => get_energy_open_shell_f_wrapper
    procedure(get_fock_type), pointer :: get_fock_f_wrapper_ptr => &
        get_fock_f_wrapper
    procedure(get_fock_jk_type), pointer :: get_fock_jk_f_wrapper_ptr => &
        get_fock_jk_f_wrapper
    procedure(obj_func_c_type), pointer :: obj_func_arh_c_wrapper_ptr => &
        obj_func_arh_c_wrapper
    procedure(update_orbs_c_type), pointer :: update_orbs_arh_c_wrapper_ptr => &
        update_orbs_arh_c_wrapper
    procedure(hess_x_c_type), pointer :: hess_x_arh_c_wrapper_ptr => &
        hess_x_arh_c_wrapper
    procedure(project_c_type), pointer :: project_arh_c_wrapper_ptr => &
        project_arh_c_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_arh_f_c
        module procedure assign_arh_c_f
    end interface

contains

    function arh_factory_c_wrapper(dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, &
                                   get_energy_c_funptr, get_fock_c_funptr, &
                                   obj_func_arh_c_funptr, update_orbs_arh_c_funptr, &
                                   project_arh_c_funptr, settings_c) result(error_c) &
        bind(C, name="arh_factory")
        !
        ! this subroutine wraps the factory function for the subroutine to convert C 
        ! variables to Fortran variables
        !
        use otr_arh, only: arh_settings_type

        real(c_rp), intent(in), target :: dm_ao_c(*), ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: get_energy_c_funptr, get_fock_c_funptr
        type(arh_settings_type_c), intent(in), value :: settings_c
        type(c_funptr), intent(out) :: obj_func_arh_c_funptr, &
                                       update_orbs_arh_c_funptr, project_arh_c_funptr
        integer(c_ip) :: error_c

        real(rp), pointer :: dm_ao_2d(:, :), dm_ao_3d(:, :, :), ao_overlap(:, :)
        procedure(get_energy_closed_shell_type), pointer :: &
            get_energy_closed_shell_funptr
        procedure(get_energy_open_shell_type), pointer :: get_energy_open_shell_funptr
        procedure(get_fock_type), pointer :: get_fock_funptr
        procedure(get_fock_jk_type), pointer :: get_fock_jk_funptr
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
            call c_f_procpointer(cptr=get_fock_c_funptr, fptr=get_fock_before_wrapping)
        else
            call c_f_procpointer(cptr=get_fock_c_funptr, &
                                 fptr=get_fock_jk_before_wrapping)
        end if

        ! associate procedure pointer to wrapper function
        if (n_particle == 1) then
            get_energy_closed_shell_funptr => get_energy_closed_shell_f_wrapper
            get_fock_funptr => get_fock_f_wrapper
        else
            get_energy_open_shell_funptr => get_energy_open_shell_f_wrapper
            get_fock_jk_funptr => get_fock_jk_f_wrapper
        end if

        ! convert settings
        settings = settings_c

        ! call factory function
        if (n_particle == 1) then
            call arh_factory_closed_shell(dm_ao_2d, ao_overlap, n_particle, n_ao, &
                                          get_energy_closed_shell_funptr, &
                                          get_fock_funptr, obj_func_arh_funptr, &
                                          update_orbs_arh_funptr, project_arh_funptr, &
                                          error, settings)
        else
            call arh_factory_open_shell(dm_ao_3d, ao_overlap, n_particle, n_ao, &
                                        get_energy_open_shell_funptr, &
                                        get_fock_jk_funptr, obj_func_arh_funptr, &
                                        update_orbs_arh_funptr, project_arh_funptr, &
                                        error, settings)
        end if

        ! associate the global procedure pointers to the Fortran function pointers
        obj_func_arh_before_wrapping => obj_func_arh_funptr
        update_orbs_arh_before_wrapping => update_orbs_arh_funptr
        project_arh_before_wrapping => project_arh_funptr

        ! get a C function pointer to the C wrapper functions
        obj_func_arh_c_funptr = c_funloc(obj_func_arh_c_wrapper)
        update_orbs_arh_c_funptr = c_funloc(update_orbs_arh_c_wrapper)
        project_arh_c_funptr = c_funloc(project_arh_c_wrapper)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function arh_factory_c_wrapper

    function get_energy_closed_shell_f_wrapper(dm, error) result(energy)
        !
        ! this subroutine wraps the energy subroutine to convert Fortran variables to C 
        ! variables for the closed-shell case
        !
        real(rp), intent(in), target :: dm(:, :)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            dm_c => dm
        else
            allocate(dm_c(size(dm, 1), size(dm, 2)))
            dm_c = real(dm, kind=c_rp)
        end if

        ! call energy C function
        error_c = get_energy_before_wrapping(dm_c, energy_c)

        ! convert arguments to Fortran kind
        energy = real(energy_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            deallocate(dm_c)
        end if

    end function get_energy_closed_shell_f_wrapper

    function get_energy_open_shell_f_wrapper(dm, error) result(energy)
        !
        ! this subroutine wraps the energy subroutine to convert Fortran variables to C 
        ! variables for the open-shell case
        !
        real(rp), intent(in), target :: dm(:, :, :)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :, :)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            dm_c => dm
        else
            allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            dm_c = real(dm, kind=c_rp)
        end if

        ! call energy C function
        error_c = get_energy_before_wrapping(dm_c, energy_c)

        ! convert arguments to Fortran kind
        energy = real(energy_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            deallocate(dm_c)
        end if

    end function get_energy_open_shell_f_wrapper

    subroutine get_fock_f_wrapper(dm, energy, fock, error)
        !
        ! this subroutine wraps the energy and Fock matrix subroutine to convert 
        ! Fortran variables to C variables
        !
        real(rp), intent(in), target :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target :: fock(:, :)
        integer(ip), intent(out) :: error

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :), fock_c(:, :)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            dm_c => dm
            fock_c => fock
        else
            allocate(dm_c(size(dm, 1), size(dm, 2)))
            allocate(fock_c(size(dm, 1), size(dm, 2)))
            dm_c = real(dm, kind=c_rp)
        end if

        ! call Fock matrix C function
        error_c = get_fock_before_wrapping(dm_c, energy_c, fock_c)

        ! convert arguments to Fortran kind
        energy = real(energy_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            fock = real(fock_c, kind=rp)
            deallocate(dm_c)
            deallocate(fock_c)
        end if

    end subroutine get_fock_f_wrapper

    subroutine get_fock_jk_f_wrapper(dm, energy, fock, coulomb, exchange, error)
        !
        ! this subroutine wraps the energy, Fock, Coulomb and exchange matrix 
        ! subroutine to convert Fortran variables to C variables
        !
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target :: fock(:, :, :), coulomb(:, :, :), &
                                         exchange(:, :, :)
        integer(ip), intent(out) :: error

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :, :), fock_c(:, :, :), coulomb_c(:, :, :), &
                               exchange_c(:, :, :)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            dm_c => dm
            fock_c => fock
            coulomb_c => coulomb
            exchange_c => exchange
        else
            allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            allocate(fock_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            allocate(coulomb_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            allocate(exchange_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            dm_c = real(dm, kind=c_rp)
        end if

        ! call Fock matrix C function
        error_c = get_fock_jk_before_wrapping(dm_c, energy_c, fock_c, coulomb_c, &
                                              exchange_c)

        ! convert arguments to Fortran kind
        energy = real(energy_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            fock = real(fock_c, kind=rp)
            coulomb = real(coulomb_c, kind=rp)
            exchange = real(exchange_c, kind=rp)
            deallocate(dm_c)
            deallocate(fock_c)
            deallocate(coulomb_c)
            deallocate(exchange_c)
        end if

    end subroutine get_fock_jk_f_wrapper

    function obj_func_arh_c_wrapper(kappa_c, func_c) result(error_c) bind(C)
        !
        ! this function wraps the objective function subroutine to convert Fortran 
        ! variables to C variables
        !
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        integer(c_ip) :: error_c

        real(rp) :: func
        real(rp), pointer :: kappa(:)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            kappa => kappa_c(:n_param)
        else
            allocate(kappa(n_param))
            kappa = real(kappa_c(:n_param), kind=rp)
        end if

        ! call obj_func Fortran function
        func = obj_func_arh_before_wrapping(kappa, error)

        ! convert arguments to Fortran kind
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            deallocate(kappa)
        end if

    end function obj_func_arh_c_wrapper

    function update_orbs_arh_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
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

        error_c = update_orbs_c_wrapper_impl(update_orbs_arh_before_wrapping, &
                                             hess_x_arh_before_wrapping, &
                                             hess_x_arh_c_wrapper, kappa_c, func_c, &
                                             grad_c, h_diag_c, hess_x_c_funptr)

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

    function project_arh_c_wrapper(vector_c) result(error_c) bind(C)
        !
        ! this function wraps the projection subroutine to convert Fortran variables to 
        ! C variables
        !
        real(c_rp), intent(inout), target :: vector_c(*)
        integer(c_ip) :: error_c

        real(rp), pointer :: vector(:)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            vector => vector_c(:n_param)
        else
            allocate(vector(n_param))
            vector = real(vector_c(:n_param), kind=rp)
        end if

        ! call project Fortran subroutine
        call project_arh_before_wrapping(vector, error)

        ! convert arguments to Fortran kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            vector_c(:n_param) = real(vector, kind=c_rp)
            deallocate(vector)
        end if

    end function project_arh_c_wrapper

    subroutine init_arh_settings_c(settings_c) bind(C, name="init_arh_settings")
        !
        ! this subroutine initializes the ARH solver settings
        !
        use otr_arh, only: default_arh_settings

        type(arh_settings_type_c), intent(inout) :: settings_c
    
        settings_c = default_arh_settings

    end subroutine init_arh_settings_c

    function arh_deconstructor_c_wrapper(dm_ao_c) result(error_c) &
        bind(C, name="arh_deconstructor")
        !
        ! this subroutine deallocates the ARH objects
        !
        real(c_rp), intent(out), target :: dm_ao_c(*)
        integer(c_ip) :: error_c

        real(rp), pointer :: dm_ao_2d(:, :), dm_ao_3d(:, :, :)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            if (n_particle == 1) then
                call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_2d, [n_ao, n_ao])
            else
                call c_f_pointer(c_loc(dm_ao_c(1)), dm_ao_3d, [n_ao, n_ao, n_particle])
            end if
        else
            if (n_particle == 1) then
                allocate(dm_ao_2d(n_ao, n_ao))
            else
                allocate(dm_ao_3d(n_ao, n_ao, n_particle))
            end if
        end if

        ! call deconstructor Fortran function
        if (n_particle == 1) then
            call arh_deconstructor_closed_shell(dm_ao_2d, error)
        else
            call arh_deconstructor_open_shell(dm_ao_3d, error)
        end if

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            if (n_particle == 1) then
                dm_ao_c(:n_ao ** 2) = reshape(real(dm_ao_2d, kind=c_rp), [n_ao ** 2])
                deallocate(dm_ao_2d)
            else
                dm_ao_c(:(n_ao ** 2 * n_particle)) = &
                    reshape(real(dm_ao_3d, kind=c_rp), [n_ao ** 2 * n_particle])
                deallocate(dm_ao_3d)
            end if
        end if

    end function arh_deconstructor_c_wrapper

    subroutine assign_arh_f_c(settings, settings_c)
        !
        ! this subroutine converts ARH settings from C to Fortran
        !
        use otr_arh, only: arh_settings_type
        use c_interface, only: logger_before_wrapping, logger_f_wrapper

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

            ! set settings to initialized
            settings%initialized = .true.
        end if

    end subroutine assign_arh_f_c

    subroutine assign_arh_c_f(settings_c, settings)
        !
        ! this subroutine converts ARH settings from Fortran to C
        !
        use otr_arh, only: arh_settings_type

        type(arh_settings_type_c), intent(out) :: settings_c
        type(arh_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%logger = c_null_funptr

            ! convert logicals
            settings_c%restricted = logical(settings%restricted, kind=c_bool)

            ! convert integers
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_arh_c_f

end module otr_arh_c_interface
