! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_c_interface

    use opentrustregion, only: ip, rp
    use c_interface, only: c_ip, c_rp
    use, intrinsic :: iso_c_binding, only: c_funptr, c_funloc, c_loc, c_f_pointer

    implicit none

    ! C-interoperable interfaces for the callback functions
    abstract interface
        function evaluate_dm_c_type(dm_ao_c, energy_c, fock_c, get_response_c_funptr) &
            result(error_c) bind(C)
            import :: c_rp, c_ip, c_funptr

            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out) :: energy_c
            real(c_rp), intent(out), optional :: fock_c(*)
            type(c_funptr), intent(out), optional :: get_response_c_funptr
            integer(c_ip) :: error_c
        end function evaluate_dm_c_type
    end interface

    abstract interface
        function get_response_c_type(dm_ao_c, response_c) result(error_c) bind(C)
            import :: c_rp, c_ip

            real(c_rp), intent(in), target :: dm_ao_c(*)
            real(c_rp), intent(out), target :: response_c(*)
            integer(c_ip) :: error_c
        end function get_response_c_type
    end interface

    ! global variables
    integer(ip) :: n_param

contains

    function update_orbs_c_wrapper_impl( &
        update_orbs_before_wrapping, hess_x_before_wrapping_funptr, &
        hess_x_c_wrapper_funptr, kappa_c, func_c, grad_c, h_diag_c, hess_x_c_funptr) &
        result(error_c)
        !
        ! this function wraps the orbital update subroutine to convert Fortran 
        ! variables to C variables
        !
        use opentrustregion, only: update_orbs_type, hess_x_type
        use c_interface, only: hess_x_c_type

        procedure(update_orbs_type), intent(in), pointer :: update_orbs_before_wrapping
        procedure(hess_x_type), intent(out), pointer :: hess_x_before_wrapping_funptr
        procedure(hess_x_c_type) :: hess_x_c_wrapper_funptr
        real(c_rp), intent(in), target :: kappa_c(*)
        real(c_rp), intent(out) :: func_c
        real(c_rp), intent(out), target :: grad_c(*), h_diag_c(*)
        type(c_funptr), intent(out) :: hess_x_c_funptr
        integer(c_ip) :: error_c

        real(rp) :: func
        real(rp), pointer :: kappa(:), grad(:), h_diag(:)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            kappa => kappa_c(:n_param)
            grad => grad_c(:n_param)
            h_diag => h_diag_c(:n_param)
        else
            allocate(kappa(n_param))
            allocate(grad(n_param))
            allocate(h_diag(n_param))
            kappa = real(kappa_c(:n_param), kind=rp)
        end if

        ! call update_orbs Fortran function
        call update_orbs_before_wrapping(kappa, func, grad, h_diag, &
                                         hess_x_before_wrapping_funptr, error)

        ! convert arguments to Fortran kind
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            grad_c(:n_param) = real(grad, kind=c_rp)
            h_diag_c(:n_param) = real(h_diag, kind=c_rp)
            deallocate(kappa)
            deallocate(grad)
            deallocate(h_diag)
        end if

        ! get a C function pointer to the hess_x wrapper function
        hess_x_c_funptr = c_funloc(hess_x_c_wrapper_funptr)

    end function update_orbs_c_wrapper_impl

    function hess_x_c_wrapper_impl(hess_x_funptr, x_c, hess_x_c) result(error_c)
        !
        ! this function wraps the Hessian linear transformation to convert Fortran 
        ! variables to C variables
        !
        use opentrustregion, only: hess_x_type

        procedure(hess_x_type), intent(in), pointer :: hess_x_funptr
        real(c_rp), intent(in), target :: x_c(*)
        real(c_rp), intent(out), target :: hess_x_c(*)
        integer(c_ip) :: error_c

        real(rp), pointer :: x(:), hess_x(:)
        integer(ip) :: error

        ! convert arguments to Fortran kind
        if (rp == c_rp) then
            x => x_c(:n_param)
            hess_x => hess_x_c(:n_param)
        else
            allocate(x(n_param))
            allocate(hess_x(n_param))
            x = real(x_c(:n_param), kind=rp)
        end if

        ! call Fortran function
        call hess_x_funptr(x, hess_x, error)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            hess_x_c(:n_param) = real(hess_x, kind=c_rp)
            deallocate(x)
            deallocate(hess_x)
        end if

    end function hess_x_c_wrapper_impl

    function obj_func_c_wrapper_impl(obj_func_before_wrapping, kappa_c, func_c) &
        result(error_c)
        !
        ! this function wraps the objective function subroutine to convert Fortran
        ! variables to C variables
        !
        use opentrustregion, only: obj_func_type

        procedure(obj_func_type), intent(in), pointer :: obj_func_before_wrapping
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
        func = obj_func_before_wrapping(kappa, error)

        ! convert arguments to Fortran kind
        func_c = real(func, kind=c_rp)
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            deallocate(kappa)
        end if

    end function obj_func_c_wrapper_impl

    function precond_c_wrapper_impl(precond_before_wrapping, residual_c, mu_c, &
                                    precond_residual_c) result(error_c)
        !
        ! this function wraps the level-shifted preconditioner subroutine to convert
        ! Fortran variables to C variables
        !
        use opentrustregion, only: precond_type

        procedure(precond_type), intent(in), pointer :: precond_before_wrapping
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
        call precond_before_wrapping(residual, mu, precond_residual, error)

        ! convert arguments to Fortran kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            precond_residual_c(:n_param) = real(precond_residual, kind=c_rp)
            deallocate(residual, precond_residual)
        end if

    end function precond_c_wrapper_impl

    function precond_pd_c_wrapper_impl(precond_pd_before_wrapping, residual_c, &
                                       precond_residual_c) result(error_c)
        !
        ! this function wraps the positive-definite preconditioner subroutine to
        ! convert Fortran variables to C variables
        !
        use opentrustregion, only: precond_pd_type

        procedure(precond_pd_type), intent(in), pointer :: precond_pd_before_wrapping
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
        call precond_pd_before_wrapping(residual, precond_residual, error)

        ! convert arguments to Fortran kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            precond_residual_c(:n_param) = real(precond_residual, kind=c_rp)
            deallocate(residual, precond_residual)
        end if

    end function precond_pd_c_wrapper_impl

    function get_extra_trial_vectors_c_wrapper_impl( &
        get_extra_trial_vectors_before_wrapping, trial_vectors_c, &
        n_extra_trial_vectors_c) result(error_c)
        !
        ! this function wraps the extra trial vector subroutine to convert Fortran
        ! variables to C variables
        !
        use opentrustregion, only: get_extra_trial_vectors_type

        procedure(get_extra_trial_vectors_type), intent(in), pointer :: &
            get_extra_trial_vectors_before_wrapping
        real(c_rp), intent(out), target :: trial_vectors_c(*)
        integer(c_ip), intent(in) :: n_extra_trial_vectors_c
        integer(c_ip) :: error_c

        real(rp), pointer :: trial_vectors(:, :)
        integer(ip) :: n_extra_trial_vectors, error

        ! convert arguments to Fortran kind
        n_extra_trial_vectors = int(n_extra_trial_vectors_c, kind=ip)
        if (rp == c_rp) then
            call c_f_pointer(c_loc(trial_vectors_c(1)), trial_vectors, &
                             [n_param, n_extra_trial_vectors])
        else
            allocate(trial_vectors(n_param, n_extra_trial_vectors))
        end if

        ! call extra trial vector Fortran subroutine
        call get_extra_trial_vectors_before_wrapping(trial_vectors, error)

        ! convert arguments to C kind
        error_c = int(error, kind=c_ip)
        if (rp /= c_rp) then
            trial_vectors_c(:n_param * n_extra_trial_vectors) = real( &
                reshape(trial_vectors, [n_param * n_extra_trial_vectors]), kind=c_rp)
            deallocate(trial_vectors)
        end if

    end function get_extra_trial_vectors_c_wrapper_impl

    subroutine evaluate_dm_f_wrapper_impl(evaluate_dm_before_wrapping, dm, energy, &
                                          error, fock, get_response_c_funptr)
        !
        ! this subroutine wraps the density matrix evaluating C function to convert
        ! Fortran variables to C variables, returning the C function pointer to the
        ! response function if one is requested
        !
        procedure(evaluate_dm_c_type), intent(in), pointer :: &
            evaluate_dm_before_wrapping
        real(rp), intent(in), target :: dm(:, :, :)
        real(rp), intent(out) :: energy
        integer(ip), intent(out) :: error
        real(rp), intent(out), optional, target :: fock(:, :, :)
        type(c_funptr), intent(out), optional :: get_response_c_funptr

        real(c_rp) :: energy_c
        real(c_rp), pointer :: dm_c(:, :, :), fock_c(:, :, :)
        integer(c_ip) :: error_c

        ! convert arguments to C kind
        nullify(fock_c)
        if (rp == c_rp) then
            dm_c => dm
            if (present(fock)) fock_c => fock
        else
            allocate(dm_c(size(dm, 1), size(dm, 2), size(dm, 3)))
            dm_c = real(dm, kind=c_rp)
            if (present(fock)) allocate(fock_c(size(dm, 1), size(dm, 2), size(dm, 3)))
        end if

        ! call density matrix evaluating C function
        if (present(get_response_c_funptr)) then
            error_c = evaluate_dm_before_wrapping(dm_c, energy_c, fock_c, &
                                                  get_response_c_funptr)
        else
            error_c = evaluate_dm_before_wrapping(dm_c, energy_c, fock_c)
        end if

        ! convert arguments to Fortran kind
        energy = real(energy_c, kind=rp)
        error = int(error_c, kind=ip)
        if (rp /= c_rp) then
            if (present(fock)) then
                fock = real(fock_c, kind=rp)
                deallocate(fock_c)
            end if
            deallocate(dm_c)
        end if

    end subroutine evaluate_dm_f_wrapper_impl

    subroutine get_response_f_wrapper_impl(get_response_before_wrapping, dm, response, &
                                           error)
        !
        ! this subroutine wraps the response C function to convert Fortran variables to
        ! C variables
        !
        procedure(get_response_c_type), intent(in), pointer :: &
            get_response_before_wrapping
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

    end subroutine get_response_f_wrapper_impl

end module otr_common_c_interface
