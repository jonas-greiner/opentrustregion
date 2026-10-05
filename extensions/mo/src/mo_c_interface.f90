! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_c_interface

    use opentrustregion, only: ip, rp, obj_func_type, update_orbs_type, hess_x_type, &
                               precond_type, precond_pd_type, &
                               get_extra_trial_vectors_type
    use c_interface, only: c_ip, c_rp, obj_func_c_type, update_orbs_c_type, &
                           hess_x_c_type, precond_c_type, precond_pd_c_type, &
                           get_extra_trial_vectors_c_type
    use otr_mo, only: standard_mo_factory_cs => mo_factory_cs, &
                      standard_mo_factory_os => mo_factory_os, &
                      standard_mo_deconstructor => mo_deconstructor
    use otr_common, only: evaluate_dm_cs_type, evaluate_dm_os_type, &
                          get_response_cs_type, get_response_os_type
    use otr_common_c_interface, only: evaluate_dm_c_type, get_response_c_type
    use, intrinsic :: iso_c_binding, only: c_bool, c_funptr, c_loc, c_f_pointer, &
                                           c_funloc, c_f_procpointer, c_associated, &
                                           c_null_funptr

    implicit none

    ! define procedure pointer which will point to the Fortran procedures
    procedure(evaluate_dm_c_type), pointer :: evaluate_dm_mo_before_wrapping => null()
    procedure(get_response_c_type), pointer :: get_response_mo_before_wrapping => null()
    procedure(obj_func_type), pointer :: obj_func_mo_before_wrapping => null()
    procedure(update_orbs_type), pointer :: update_orbs_mo_before_wrapping => null()
    procedure(hess_x_type), pointer :: hess_x_mo_before_wrapping => null()
    procedure(precond_type), pointer :: precond_mo_before_wrapping => null()
    procedure(precond_pd_type), pointer :: precond_pd_mo_before_wrapping => null()
    procedure(get_extra_trial_vectors_type), pointer :: &
        get_extra_trial_vectors_mo_before_wrapping => null()

    ! derived type for MO settings
    type, bind(C) :: mo_settings_type_c
        type(c_funptr) :: logger
        logical(c_bool) :: initialized
        integer(c_ip) :: verbose
    end type

    ! MO coefficients passed from C, which are updated from their Fortran copy after
    ! every orbital update if the real kinds differ
    real(c_rp), pointer :: mo_coeff_3d_c(:, :, :) => null()

    procedure(standard_mo_factory_cs), pointer :: mo_factory_cs => &
        standard_mo_factory_cs
    procedure(standard_mo_factory_os), pointer :: mo_factory_os => &
        standard_mo_factory_os
    procedure(standard_mo_deconstructor), pointer :: mo_deconstructor => &
        standard_mo_deconstructor

    ! create function pointers to ensure that routines comply with interface
    procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_mo_cs_f_wrapper_ptr => &
        evaluate_dm_mo_cs_f_wrapper
    procedure(evaluate_dm_os_type), pointer :: evaluate_dm_mo_os_f_wrapper_ptr => &
        evaluate_dm_mo_os_f_wrapper
    procedure(get_response_cs_type), pointer :: get_response_mo_cs_f_wrapper_ptr => &
        get_response_mo_cs_f_wrapper
    procedure(get_response_os_type), pointer :: get_response_mo_os_f_wrapper_ptr => &
        get_response_mo_os_f_wrapper
    procedure(obj_func_c_type), pointer :: obj_func_mo_c_wrapper_ptr => &
        obj_func_mo_c_wrapper
    procedure(update_orbs_c_type), pointer :: update_orbs_mo_c_wrapper_ptr => &
        update_orbs_mo_c_wrapper
    procedure(hess_x_c_type), pointer :: hess_x_mo_c_wrapper_ptr => hess_x_mo_c_wrapper
    procedure(precond_c_type), pointer :: precond_mo_c_wrapper_ptr => &
        precond_mo_c_wrapper
    procedure(precond_pd_c_type), pointer :: precond_pd_mo_c_wrapper_ptr => &
        precond_pd_mo_c_wrapper
    procedure(get_extra_trial_vectors_c_type), pointer :: &
        get_extra_trial_vectors_mo_c_wrapper_ptr => get_extra_trial_vectors_mo_c_wrapper

    ! interfaces for converting C settings to Fortran settings
    interface assignment(=)
        module procedure assign_mo_f_c
        module procedure assign_mo_c_f
    end interface

contains

    function mo_factory_c_wrapper(mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, &
                                  n_ao_c, n_mo_c, evaluate_dm_c_funptr, &
                                  obj_func_mo_c_funptr, update_orbs_mo_c_funptr, &
                                  solver_settings_c, settings_c, orbsym_c) &
        result(error_c) bind(C, name="mo_factory")
        !
        ! this subroutine wraps the factory function for the subroutine to convert C
        ! variables to Fortran variables
        !
        use opentrustregion, only: solver_settings_type
        use c_interface, only: solver_settings_type_c
        use otr_mo, only: mo_settings_type, count_mo_params
        use otr_common_c_interface, only: n_param

        real(c_rp), intent(inout), target :: mo_coeff_c(*)
        real(c_rp), intent(in), target :: ao_overlap_c(*)
        integer(c_ip), intent(in) :: n_occ_c(*)
        integer(c_ip), intent(in), value :: n_particle_c, n_ao_c, n_mo_c
        type(c_funptr), intent(in), value :: evaluate_dm_c_funptr
        type(solver_settings_type_c), intent(inout) :: solver_settings_c
        type(mo_settings_type_c), intent(inout) :: settings_c
        type(c_funptr), intent(out) :: obj_func_mo_c_funptr, update_orbs_mo_c_funptr
        integer(c_ip), intent(in), optional :: orbsym_c(*)
        integer(c_ip) :: error_c

        real(rp), pointer, contiguous :: mo_coeff_2d(:, :)
        real(rp), pointer, contiguous :: mo_coeff_3d(:, :, :)
        real(rp), pointer :: ao_overlap(:, :)
        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        procedure(obj_func_type), pointer :: obj_func_mo_funptr
        procedure(update_orbs_type), pointer :: update_orbs_mo_funptr
        type(solver_settings_type) :: solver_settings
        type(mo_settings_type) :: settings
        integer(ip) :: n_particle, n_ao, n_mo, error
        integer(ip), allocatable :: n_occ(:), orbsym(:, :), orbsym_1d(:)

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
            call c_f_pointer(c_loc(mo_coeff_c(1)), mo_coeff_3d_c, &
                             [n_ao, n_mo, n_particle])
            if (n_particle == 1) then
                allocate(mo_coeff_2d(n_ao, n_mo))
                mo_coeff_2d = real(mo_coeff_3d_c(:, :, 1), kind=rp)
            else
                allocate(mo_coeff_3d(n_ao, n_mo, n_particle))
                mo_coeff_3d = real(mo_coeff_3d_c, kind=rp)
            end if
            allocate(ao_overlap(n_ao, n_ao))
            ao_overlap = reshape(real(ao_overlap_c(:n_ao**2), kind=rp), [n_ao, n_ao])
        end if

        ! associate the input C pointers to Fortran procedure pointers
        call c_f_procpointer(cptr=evaluate_dm_c_funptr, &
                             fptr=evaluate_dm_mo_before_wrapping)

        ! associate procedure pointer to wrapper function
        if (n_particle == 1) then
            evaluate_dm_cs_funptr => evaluate_dm_mo_cs_f_wrapper
        else
            evaluate_dm_os_funptr => evaluate_dm_mo_os_f_wrapper
        end if

        ! convert settings
        settings = settings_c

        ! call factory function
        if (n_particle == 1) then
            call mo_factory_cs(mo_coeff_2d, ao_overlap, n_occ(1), n_particle, n_ao, &
                               n_mo, evaluate_dm_cs_funptr, obj_func_mo_funptr, &
                               update_orbs_mo_funptr, solver_settings, error, &
                               settings, orbsym_1d)
        else
            call mo_factory_os(mo_coeff_3d, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                               evaluate_dm_os_funptr, obj_func_mo_funptr, &
                               update_orbs_mo_funptr, solver_settings, error, &
                               settings, orbsym)
        end if

        ! calculate the number of parameters of a successful setup and store it
        ! globally to access assumed size arrays passed from C to Fortran
        if (error == 0) n_param = count_mo_params(n_occ, n_mo, orbsym)

        ! associate the global procedure pointers to the Fortran function pointers
        obj_func_mo_before_wrapping => obj_func_mo_funptr
        update_orbs_mo_before_wrapping => update_orbs_mo_funptr

        ! get a C function pointer to the C wrapper functions
        obj_func_mo_c_funptr = c_funloc(obj_func_mo_c_wrapper)
        update_orbs_mo_c_funptr = c_funloc(update_orbs_mo_c_wrapper)

        ! copy the solver settings
        if (error == 0) &
            call mo_set_solver_settings_c(solver_settings, solver_settings_c)

        ! convert return arguments to C kind
        error_c = int(error, kind=c_ip)

    end function mo_factory_c_wrapper

    subroutine mo_set_solver_settings_c(solver_settings, solver_settings_c)
        !
        ! this subroutine wires the C wrappers of the MO-basis preconditioners and
        ! extra trial vectors into C solver settings, leaving the projection to the
        ! caller
        !
        use opentrustregion, only: solver_settings_type, default_solver_settings
        use c_interface, only: solver_settings_type_c, assignment(=)

        type(solver_settings_type), intent(in) :: solver_settings
        type(solver_settings_type_c), intent(inout) :: solver_settings_c

        ! initialize settings
        if (.not. solver_settings_c%initialized) &
            solver_settings_c = default_solver_settings

        ! associate the global procedure pointers to the Fortran function pointers and
        ! get C function pointers to the C wrapper functions
        if (associated(solver_settings%precond)) then
            precond_mo_before_wrapping => solver_settings%precond
            solver_settings_c%precond = c_funloc(precond_mo_c_wrapper)
        end if
        if (associated(solver_settings%precond_pd)) then
            precond_pd_mo_before_wrapping => solver_settings%precond_pd
            solver_settings_c%precond_pd = c_funloc(precond_pd_mo_c_wrapper)
        end if
        if (associated(solver_settings%get_extra_trial_vectors)) then
            get_extra_trial_vectors_mo_before_wrapping => &
                solver_settings%get_extra_trial_vectors
            solver_settings_c%get_extra_trial_vectors = &
                c_funloc(get_extra_trial_vectors_mo_c_wrapper)
        end if
        if (associated(solver_settings%stability_settings%precond, &
                       solver_settings%precond)) &
            solver_settings_c%stability_settings%precond = &
            c_funloc(precond_mo_c_wrapper)
        if (associated(solver_settings%stability_settings%get_extra_trial_vectors, &
                       solver_settings%get_extra_trial_vectors)) &
            solver_settings_c%stability_settings%get_extra_trial_vectors = &
            c_funloc(get_extra_trial_vectors_mo_c_wrapper)

    end subroutine mo_set_solver_settings_c

    subroutine evaluate_dm_mo_cs_f_wrapper(dm, energy, fock, get_response_cs_funptr, &
                                           error)
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

        real(rp), pointer :: dm_3d(:, :, :), fock_3d(:, :, :)
        procedure(get_response_os_type), pointer :: get_response_funptr

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        nullify(fock_3d, get_response_funptr)
        if (present(fock)) fock_3d(1:size(fock, 1), 1:size(fock, 2), 1:1) => fock
        if (present(get_response_cs_funptr)) then
            call evaluate_dm_mo_os_f_wrapper(dm_3d, energy, fock_3d, &
                                             get_response_funptr, error)
            get_response_cs_funptr => null()
            if (error == 0) get_response_cs_funptr => get_response_mo_cs_f_wrapper
        else
            call evaluate_dm_mo_os_f_wrapper(dm_3d, energy, fock_3d, error=error)
        end if

    end subroutine evaluate_dm_mo_cs_f_wrapper

    subroutine evaluate_dm_mo_os_f_wrapper(dm, energy, fock, get_response_funptr, error)
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

        type(c_funptr) :: get_response_c_funptr

        ! call density matrix evaluating C function, and associate a returned C pointer
        ! to the response function with a Fortran procedure pointer
        if (present(get_response_funptr)) then
            call evaluate_dm_f_wrapper_impl(evaluate_dm_mo_before_wrapping, dm, &
                                            energy, error, fock, get_response_c_funptr)
            get_response_funptr => null()
            if (error == 0) then
                call c_f_procpointer(cptr=get_response_c_funptr, &
                                     fptr=get_response_mo_before_wrapping)
                get_response_funptr => get_response_mo_os_f_wrapper
            end if
        else
            call evaluate_dm_f_wrapper_impl(evaluate_dm_mo_before_wrapping, dm, &
                                            energy, error, fock)
        end if

    end subroutine evaluate_dm_mo_os_f_wrapper

    subroutine get_response_mo_cs_f_wrapper(dm, response, error)
        !
        ! this subroutine wraps the response subroutine to convert Fortran variables to
        ! C variables
        !
        real(rp), intent(in), target, contiguous :: dm(:, :)
        real(rp), intent(out), target, contiguous :: response(:, :)
        integer(ip), intent(out) :: error

        real(rp), pointer :: dm_3d(:, :, :), response_3d(:, :, :)

        dm_3d(1:size(dm, 1), 1:size(dm, 2), 1:1) => dm
        response_3d(1:size(response, 1), 1:size(response, 2), 1:1) => response
        call get_response_mo_os_f_wrapper(dm_3d, response_3d, error)

    end subroutine get_response_mo_cs_f_wrapper

    subroutine get_response_mo_os_f_wrapper(dm, response, error)
        !
        ! this subroutine wraps the response subroutine to convert Fortran variables to
        ! C variables
        !
        use otr_common_c_interface, only: get_response_f_wrapper_impl

        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out), target :: response(:, :, :)
        integer(ip), intent(out) :: error

        call get_response_f_wrapper_impl(get_response_mo_before_wrapping, dm, &
                                         response, error)

    end subroutine get_response_mo_os_f_wrapper

    function obj_func_mo_c_wrapper(kappa_c, func_c) result(error_c) bind(C)
        !
        ! this function wraps the objective function subroutine to convert Fortran
        ! variables to C variables
        !
        use otr_common_c_interface, only: obj_func_c_wrapper_impl

        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        integer(c_ip) :: error_c

        error_c = obj_func_c_wrapper_impl(obj_func_mo_before_wrapping, kappa_c, func_c)

    end function obj_func_mo_c_wrapper

    function update_orbs_mo_c_wrapper(kappa_c, func_c, grad_c, h_diag_c, &
                                      hess_x_c_funptr) result(error_c) bind(C)
        !
        ! this function wraps the orbital update subroutine to convert Fortran
        ! variables to C variables
        !
        use otr_common_c_interface, only: update_orbs_c_wrapper_impl
        use otr_mo, only: mo_object

        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(out) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        error_c = update_orbs_c_wrapper_impl( &
            update_orbs_mo_before_wrapping, hess_x_mo_before_wrapping, &
            hess_x_mo_c_wrapper, kappa_c, func_c, grad_c, h_diag_c, hess_x_c_funptr)

        ! update the rotated orbitals passed from C if they could not be rotated in
        ! place
        if (rp /= c_rp) mo_coeff_3d_c = real(mo_object%mo_coeff, kind=c_rp)

    end function update_orbs_mo_c_wrapper

    function hess_x_mo_c_wrapper(x_c, hess_x_c) result(error_c) bind(C)
        !
        ! this function wraps the Hessian linear transformation to convert Fortran
        ! variables to C variables
        !
        use otr_common_c_interface, only: hess_x_c_wrapper_impl

        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        integer(c_ip) :: error_c

        error_c = hess_x_c_wrapper_impl(hess_x_mo_before_wrapping, x_c, hess_x_c)

    end function hess_x_mo_c_wrapper

    function precond_mo_c_wrapper(residual_c, mu_c, precond_residual_c) &
        result(error_c) bind(C)
        !
        ! this function wraps the level-shifted preconditioner subroutine to convert
        ! Fortran variables to C variables
        !
        use otr_common_c_interface, only: precond_c_wrapper_impl

        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(in) :: mu_c
        real(c_rp), intent(out), target :: precond_residual_c(*)
        integer(c_ip) :: error_c

        error_c = precond_c_wrapper_impl(precond_mo_before_wrapping, residual_c, mu_c, &
                                         precond_residual_c)

    end function precond_mo_c_wrapper

    function precond_pd_mo_c_wrapper(residual_c, precond_residual_c) result(error_c) &
        bind(C)
        !
        ! this function wraps the positive-definite preconditioner subroutine to
        ! convert Fortran variables to C variables
        !
        use otr_common_c_interface, only: precond_pd_c_wrapper_impl

        real(c_rp), intent(in), target :: residual_c(*)
        real(c_rp), intent(out), target :: precond_residual_c(*)
        integer(c_ip) :: error_c

        error_c = precond_pd_c_wrapper_impl(precond_pd_mo_before_wrapping, residual_c, &
                                            precond_residual_c)

    end function precond_pd_mo_c_wrapper

    function get_extra_trial_vectors_mo_c_wrapper( &
        trial_vectors_c, n_extra_trial_vectors_c) result(error_c) bind(C)
        !
        ! this function wraps the extra trial vector subroutine to convert Fortran
        ! variables to C variables
        !
        use otr_common_c_interface, only: get_extra_trial_vectors_c_wrapper_impl

        real(c_rp), intent(out), target :: trial_vectors_c(*)
        integer(c_ip), intent(in), value :: n_extra_trial_vectors_c
        integer(c_ip) :: error_c

        error_c = get_extra_trial_vectors_c_wrapper_impl( &
            get_extra_trial_vectors_mo_before_wrapping, trial_vectors_c, &
            n_extra_trial_vectors_c)

    end function get_extra_trial_vectors_mo_c_wrapper

    subroutine init_mo_settings_c(settings_c) bind(C, name="init_mo_settings")
        !
        ! this subroutine initializes the MO settings
        !
        use otr_mo, only: default_mo_settings

        type(mo_settings_type_c), intent(inout) :: settings_c

        settings_c = default_mo_settings

    end subroutine init_mo_settings_c

    subroutine mo_deconstructor_c_wrapper() bind(C, name="mo_deconstructor")
        !
        ! this subroutine deallocates the MO objects
        !
        call mo_deconstructor()

    end subroutine mo_deconstructor_c_wrapper

    subroutine assign_mo_f_c(settings, settings_c)
        !
        ! this subroutine converts MO settings from C to Fortran
        !
        use otr_mo, only: mo_settings_type
        use c_interface, only: logger_before_wrapping, logger_f_wrapper

        type(mo_settings_type), intent(out) :: settings
        type(mo_settings_type_c), intent(in) :: settings_c

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

    end subroutine assign_mo_f_c

    subroutine assign_mo_c_f(settings_c, settings)
        !
        ! this subroutine converts MO settings from Fortran to C
        !
        use otr_mo, only: mo_settings_type

        type(mo_settings_type_c), intent(out) :: settings_c
        type(mo_settings_type), intent(in) :: settings

        if (settings%initialized) then
            ! callback functions cannot be converted
            settings_c%logger = c_null_funptr

            ! convert integers
            settings_c%verbose = int(settings%verbose, kind=c_ip)

            ! set settings to initialized
            settings_c%initialized = .true._c_bool
        end if

    end subroutine assign_mo_c_f

end module otr_mo_c_interface
