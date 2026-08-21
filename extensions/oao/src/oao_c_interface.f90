! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_oao_c_interface

    use opentrustregion, only: ip, rp, obj_func_type, update_orbs_type, hess_x_type, &
                               precond_type, precond_pd_type, project_type
    use c_interface, only: c_ip, c_rp, obj_func_c_type, update_orbs_c_type, &
                           hess_x_c_type, precond_c_type, precond_pd_c_type, &
                           project_c_type
    use otr_oao, only: standard_oao_factory_cs => oao_factory_cs, &
                       standard_oao_factory_os => oao_factory_os, &
                       standard_oao_deconstructor => oao_deconstructor, &
                       get_energy_cs_type, get_energy_os_type, update_dm_cs_type, &
                       update_dm_os_type, get_response_cs_type, get_response_os_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_loc, c_f_pointer, &
                                           c_funloc, c_f_procpointer, c_associated, &
                                           c_null_funptr

    implicit none

    ! define procedure pointer which will point to the Fortran procedures
    procedure(get_energy_c_type), pointer :: get_energy_before_wrapping => null()
    procedure(update_dm_c_type), pointer :: update_dm_before_wrapping => null()
    procedure(get_response_c_type), pointer :: get_response_before_wrapping => null()
    procedure(obj_func_type), pointer :: obj_func_oao_before_wrapping => null()
    procedure(update_orbs_type), pointer :: update_orbs_oao_before_wrapping => null()
    procedure(hess_x_type), pointer :: hess_x_oao_before_wrapping => null()
    procedure(precond_type), pointer :: precond_oao_before_wrapping => null()
    procedure(precond_pd_type), pointer :: precond_pd_oao_before_wrapping => null()
    procedure(project_type), pointer :: project_oao_before_wrapping => null()

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
        function update_dm_c_type(dm_ao_c, energy_c, fock_c, get_response_c_funptr) &
            result(error_c) bind(C)
            import :: c_rp, c_ip, c_funptr
    
            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            real(c_rp), intent(out), target :: fock_c(*)
            type(c_funptr), intent(out) :: get_response_c_funptr
            integer(c_ip) :: error_c
        end function update_dm_c_type
    end interface

    abstract interface
        function get_response_c_type(dm_ao_c, response_c) result(error_c) bind(C)
            import :: c_rp, c_ip
    
            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out), target :: response_c(*)
            integer(c_ip) :: error_c
        end function get_response_c_type
    end interface

    ! derived type for OAO settings
    type, bind(C) :: oao_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose
    end type

    ! global variables
    integer(ip) :: n_particle, n_ao
    real(c_rp), pointer :: dm_ao_3d_c(:, :, :)

    procedure(standard_oao_factory_cs), pointer :: oao_factory_cs => &
        standard_oao_factory_cs
    procedure(standard_oao_factory_os), pointer :: oao_factory_os => &
        standard_oao_factory_os
    procedure(standard_oao_deconstructor), pointer :: oao_deconstructor => &
        standard_oao_deconstructor

    ! create function pointers to ensure that routines comply with interface
    procedure(get_energy_cs_type), pointer :: get_energy_cs_f_wrapper_ptr => &
        get_energy_cs_f_wrapper
    procedure(get_energy_os_type), pointer :: get_energy_os_f_wrapper_ptr => &
        get_energy_os_f_wrapper
    procedure(update_dm_cs_type), pointer :: update_dm_cs_f_wrapper_ptr => &
        update_dm_cs_f_wrapper
    procedure(update_dm_os_type), pointer :: update_dm_os_f_wrapper_ptr => &
        update_dm_os_f_wrapper
    procedure(get_response_cs_type), pointer :: get_response_cs_f_wrapper_ptr => &
        get_response_cs_f_wrapper
    procedure(get_response_os_type), pointer :: get_response_os_f_wrapper_ptr => &
        get_response_os_f_wrapper
    procedure(obj_func_c_type), pointer :: obj_func_oao_c_wrapper_ptr => &
        obj_func_oao_c_wrapper
    procedure(update_orbs_c_type), pointer :: update_orbs_oao_c_wrapper_ptr => &
        update_orbs_oao_c_wrapper
    procedure(hess_x_c_type), pointer :: hess_x_oao_c_wrapper_ptr => &
        hess_x_oao_c_wrapper
    procedure(precond_c_type), pointer :: precond_oao_c_wrapper_ptr => &
        precond_oao_c_wrapper
    procedure(precond_pd_c_type), pointer :: precond_pd_oao_c_wrapper_ptr => &
        precond_pd_oao_c_wrapper
    procedure(project_c_type), pointer :: project_oao_c_wrapper_ptr => &
        project_oao_c_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_oao_f_c
        module procedure assign_oao_c_f
    end interface

contains

    function oao_factory_c_wrapper(dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, &
                                   get_energy_c_funptr, update_dm_c_funptr, &
                                   obj_func_oao_c_funptr, update_orbs_oao_c_funptr, &
                                   precond_oao_c_funptr, precond_pd_oao_c_funptr, &
                                   project_oao_c_funptr, settings_c) result(error_c) &
        bind(C, name="oao_factory")
        !
        ! this subroutine wraps the factory function for the subroutine to convert C
        ! variables to Fortran variables
        !
        use otr_oao, only: oao_settings_type
        use otr_common_c_interface, only: n_param

        real(c_rp), intent(in), target :: dm_ao_c(*), ao_overlap_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c
        type(c_funptr), intent(in), value :: get_energy_c_funptr, update_dm_c_funptr
        type(oao_settings_type_c), intent(inout) :: settings_c
        type(c_funptr), intent(out) :: obj_func_oao_c_funptr, &
                                       update_orbs_oao_c_funptr, precond_oao_c_funptr, &
                                       project_oao_c_funptr, precond_pd_oao_c_funptr
        integer(c_ip) :: error_c

        real(rp), pointer, contiguous :: dm_ao_2d(:, :)
        real(rp), pointer, contiguous :: dm_ao_3d(:, :, :)
        real(rp), pointer :: ao_overlap(:, :)
        procedure(get_energy_cs_type), pointer :: get_energy_cs_funptr
        procedure(get_energy_os_type), pointer :: get_energy_os_funptr
        procedure(update_dm_cs_type), pointer :: update_dm_cs_funptr
        procedure(update_dm_os_type), pointer :: update_dm_os_funptr
        procedure(obj_func_type), pointer :: obj_func_oao_funptr
        procedure(update_orbs_type), pointer :: update_orbs_oao_funptr
        procedure(precond_type), pointer :: precond_oao_funptr
        procedure(precond_pd_type), pointer :: precond_pd_oao_funptr
        procedure(project_type), pointer :: project_oao_funptr
        type(oao_settings_type) :: settings
        integer(ip) :: error

        ! convert number of AOs to Fortran kind, calculate number of parameters and 
        ! store globally to access assumed size arrays passed from C to Fortran
        n_particle = int(n_particle_c, kind=ip)
        n_ao = int(n_ao_c, kind=ip)
        n_param = n_particle * n_ao * (n_ao - 1) / 2

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
        call c_f_procpointer(cptr=update_dm_c_funptr, fptr=update_dm_before_wrapping)

        ! associate procedure pointer to wrapper function
        if (n_particle == 1) then
            get_energy_cs_funptr => get_energy_cs_f_wrapper
            update_dm_cs_funptr => update_dm_cs_f_wrapper
        else
            get_energy_os_funptr => get_energy_os_f_wrapper
            update_dm_os_funptr => update_dm_os_f_wrapper
        end if

        ! convert settings
        settings = settings_c

        ! call factory function
        if (n_particle == 1) then
            call oao_factory_cs(dm_ao_2d, ao_overlap, n_particle, n_ao, &
                                get_energy_cs_funptr, update_dm_cs_funptr, &
                                obj_func_oao_funptr, update_orbs_oao_funptr, &
                                precond_oao_funptr, precond_pd_oao_funptr, &
                                project_oao_funptr, error, settings)
        else
            call oao_factory_os(dm_ao_3d, ao_overlap, n_particle, n_ao, &
                                get_energy_os_funptr, update_dm_os_funptr, &
                                obj_func_oao_funptr, update_orbs_oao_funptr, &
                                precond_oao_funptr, precond_pd_oao_funptr, &
                                project_oao_funptr, error, settings)
        end if

        ! associate the global procedure pointers to the Fortran function pointers
        obj_func_oao_before_wrapping => obj_func_oao_funptr
        update_orbs_oao_before_wrapping => update_orbs_oao_funptr
        precond_oao_before_wrapping => precond_oao_funptr
        precond_pd_oao_before_wrapping => precond_pd_oao_funptr
        project_oao_before_wrapping => project_oao_funptr

        ! get a C function pointer to the C wrapper functions
        obj_func_oao_c_funptr = c_funloc(obj_func_oao_c_wrapper)
        update_orbs_oao_c_funptr = c_funloc(update_orbs_oao_c_wrapper)
        precond_oao_c_funptr = c_funloc(precond_oao_c_wrapper)
        precond_pd_oao_c_funptr = c_funloc(precond_pd_oao_c_wrapper)
        project_oao_c_funptr = c_funloc(project_oao_c_wrapper)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function oao_factory_c_wrapper

    function get_energy_cs_f_wrapper(dm, error) result(energy)
        !
        ! this subroutine wraps the energy subroutine to convert Fortran variables to C 
        ! variables for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        real(rp), pointer :: dm_3d(:, :, :)

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        energy = get_energy_os_f_wrapper(dm_3d, error)

    end function get_energy_cs_f_wrapper

    function get_energy_os_f_wrapper(dm, error) result(energy)
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

    end function get_energy_os_f_wrapper

    subroutine update_dm_cs_f_wrapper(dm, energy, fock, get_response_cs_funptr, error)
        !
        ! this subroutine wraps the density matrix updating subroutine to convert 
        ! Fortran variables to C variables for the closed-shell case
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target, contiguous :: fock(:, :)
        procedure(get_response_cs_type), intent(out), pointer :: get_response_cs_funptr

        integer(ip), intent(out) :: error

        real(rp), pointer :: dm_3d(:, :, :), fock_3d(:, :, :)
        procedure(get_response_os_type), pointer :: get_response_funptr

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        fock_3d(1:size(fock, 1), 1:size(fock, 2), 1:1) => fock
        get_response_funptr => null()
        call update_dm_os_f_wrapper(dm_3d, energy, fock_3d, get_response_funptr, error)

        ! associate procedure pointer to wrapper function
        get_response_cs_funptr => get_response_cs_f_wrapper

    end subroutine update_dm_cs_f_wrapper

    subroutine update_dm_os_f_wrapper(dm, energy, fock, get_response_funptr, error)
        !
        ! this subroutine wraps the density matrix updating subroutine to convert 
        ! Fortran variables to C variables for the open-shell case
        !
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        real(rp), intent(out), target :: fock(:, :, :)
        procedure(get_response_os_type), intent(out), pointer :: get_response_funptr
        integer(ip), intent(out) :: error

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :, :), fock_c(:, :, :)
        type(c_funptr) :: get_response_c_funptr
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            dm_c => dm
            fock_c => fock
        else
            allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)), &
            fock_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            dm_c = real(dm, kind=c_rp)
        end if

        ! call Fock matrix C function
        error_c = update_dm_before_wrapping(dm_c, energy_c, fock_c, &
                                            get_response_c_funptr)

        ! convert arguments to Fortran kind
        energy = real(energy_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            fock = real(fock_c, kind=rp)
            deallocate(dm_c, fock_c)
        end if

        ! associate the input C pointer to get_response function to a Fortran procedure
        ! pointer
        call c_f_procpointer(cptr=get_response_c_funptr, &
                             fptr=get_response_before_wrapping)

        ! associate procedure pointer to wrapper function
        get_response_funptr => get_response_os_f_wrapper

    end subroutine update_dm_os_f_wrapper

    subroutine get_response_cs_f_wrapper(dm, response, error)
        !
        ! this subroutine wraps the response subroutine to convert Fortran variables 
        ! to C variables
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out), target, contiguous :: response(:, :)
        integer(ip), intent(out) :: error

        real(rp), pointer :: dm_3d(:, :, :), response_3d(:, :, :)

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        response_3d(1:size(response, 1), 1:size(response, 2), 1:1) => response
        call get_response_os_f_wrapper(dm_3d, response_3d, error)

    end subroutine get_response_cs_f_wrapper

    subroutine get_response_os_f_wrapper(dm, response, error)
        !
        ! this subroutine wraps the response subroutine to convert Fortran variables 
        ! to C variables
        !
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out), target :: response(:, :, :)
        integer(ip), intent(out) :: error

        real(c_rp), pointer :: dm_c(:, :, :), response_c(:, :, :)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        if (rp == c_rp) then
            dm_c => dm
            response_c => response
        else
            allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)), &
            response_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            dm_c = real(dm, kind=c_rp)
        end if

        ! call response C function
        error_c = get_response_before_wrapping(dm_c, response_c)

        ! convert arguments to Fortran kind
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            response = real(response_c, kind=rp)
            deallocate(dm_c, response_c)
        end if

    end subroutine get_response_os_f_wrapper

    function obj_func_oao_c_wrapper(kappa_c, func_c) result(error_c) bind(C)
        !
        ! this function wraps the objective function subroutine to convert Fortran 
        ! variables to C variables
        !
        use otr_common_c_interface, only: n_param

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
        func = obj_func_oao_before_wrapping(kappa, error)

        ! convert arguments to Fortran kind
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            deallocate(kappa)
        end if

    end function obj_func_oao_c_wrapper

    function update_orbs_oao_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                       hess_x_c_funptr) result(error_c) bind(C)
        !
        ! this function wraps the orbital update subroutine to convert Fortran 
        ! variables to C variables
        !
        use otr_common_c_interface, only: update_orbs_c_wrapper_impl
        use otr_oao, only: oao_object

        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(out) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        error_c = update_orbs_c_wrapper_impl(update_orbs_oao_before_wrapping, &
                                             hess_x_oao_before_wrapping, &
                                             hess_x_oao_c_wrapper, kappa_c, func_c, &
                                             grad_c, h_diag_c, hess_x_c_funptr)

        if (rp /= c_rp) dm_ao_3d_c = real(oao_object%dm_ao, kind=c_rp)

    end function update_orbs_oao_c_wrapper

    function hess_x_oao_c_wrapper(x_c, hess_x_c) result(error_c) bind(C)
        !
        ! this function wraps the Hessian linear transformation to convert Fortran 
        ! variables to C variables
        !
        use otr_common_c_interface, only: hess_x_c_wrapper_impl

        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        integer(c_ip) :: error_c

        error_c = hess_x_c_wrapper_impl(hess_x_oao_before_wrapping, x_c, hess_x_c)

    end function hess_x_oao_c_wrapper

    function precond_oao_c_wrapper(residual_c, mu_c, precond_residual_c) &
        result(error_c) bind(C)
        !
        ! this function wraps the level-shifted preconditioner subroutine to convert
        ! Fortran variables to C variables
        !
        use otr_common_c_interface, only: n_param

        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(in) :: mu_c
        real(c_rp), intent(out), target :: precond_residual_c(*)
        integer(c_ip) :: error_c

        real(rp) :: mu
        real(rp), pointer :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        mu = real(mu_c, kind=rp)
        if (rp == c_rp) then
            residual => residual_c(:n_param)
            precond_residual => precond_residual_c(:n_param)
        else
            allocate(residual(n_param), precond_residual(n_param))
            residual = real(residual_c(:n_param), kind=rp)
        end if

        ! call preconditioner Fortran subroutine
        call precond_oao_before_wrapping(residual, mu, precond_residual, error)

        ! convert arguments to Fortran kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            precond_residual_c(:n_param) = real(precond_residual, kind=c_rp)
            deallocate(residual, precond_residual)
        end if

    end function precond_oao_c_wrapper

    function precond_pd_oao_c_wrapper(residual_c, precond_residual_c) result(error_c) &
        bind(C)
        !
        ! this function wraps the positive-definite preconditioner subroutine to
        ! convert Fortran variables to C variables
        !
        use otr_common_c_interface, only: n_param

        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(out), target :: precond_residual_c(*)
        integer(c_ip) :: error_c

        real(rp), pointer :: residual(:), precond_residual(:)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            residual => residual_c(:n_param)
            precond_residual => precond_residual_c(:n_param)
        else
            allocate(residual(n_param), precond_residual(n_param))
            residual = real(residual_c(:n_param), kind=rp)
        end if

        ! call preconditioner Fortran subroutine
        call precond_pd_oao_before_wrapping(residual, precond_residual, error)

        ! convert arguments to Fortran kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            precond_residual_c(:n_param) = real(precond_residual, kind=c_rp)
            deallocate(residual, precond_residual)
        end if

    end function precond_pd_oao_c_wrapper

    function project_oao_c_wrapper(vector_c) result(error_c) bind(C)
        !
        ! this function wraps the projection subroutine to convert Fortran variables to 
        ! C variables
        !
        use otr_common_c_interface, only: n_param

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
        call project_oao_before_wrapping(vector, error)

        ! convert arguments to Fortran kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            vector_c(:n_param) = real(vector, kind=c_rp)
            deallocate(vector)
        end if

    end function project_oao_c_wrapper

    subroutine init_oao_settings_c(settings_c) bind(C, name="init_oao_settings")
        !
        ! this subroutine initializes the OAO solver settings
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
        use c_interface, only: logger_before_wrapping, logger_f_wrapper

        type(oao_settings_type), intent(out) :: settings
        type(oao_settings_type_c), intent(in) :: settings_c

        if (settings_c%initialized) then
            ! convert callback functions
            if (c_associated(settings_c%logger)) then
                call c_f_procpointer(cptr=settings_c%logger, &
                                     fptr=logger_before_wrapping)
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

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_oao_c_f

end module otr_oao_c_interface
