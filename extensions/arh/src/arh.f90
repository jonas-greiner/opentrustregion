! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh

    use opentrustregion, only: rp, ip, kw_len, obj_func_type, update_orbs_type, &
                               hess_x_type, precond_type, precond_pd_type, &
                               get_extra_trial_vectors_type, solver_settings_type
    use otr_common, only: orbital_settings_type, default_orbital_settings, &
                          orbital_basis_type, channel_rows
    use otr_oao, only: oao_type
    use otr_mo, only: mo_channel_type, mo_type

    implicit none

    ! define useful parameters
    real(rp), parameter :: history_step_ratio = 1e2_rp, ms_sr1_skip_thresh = 0.05_rp, &
                           eig_val_noise_factor = 10.0_rp

    type, extends(orbital_settings_type) :: arh_settings_type
        character(len=kw_len) :: arh_type
    contains
        procedure :: init => init_arh_settings
    end type

    type(arh_settings_type), parameter :: default_arh_settings = arh_settings_type( &
        orbital_settings_type=default_orbital_settings, arh_type="ms_sr1")

    ! define setting options
    character(len=kw_len), parameter :: arh_types(5) = &
        [character(len=kw_len) :: "arh", "symm_arh", "ms_psb", "ms_sp", "ms_sr1"]

    ! micro iteration limit of the subsystem solver, which the linear transformations
    ! of the approximate Hessian make cheap compared to the function evaluations
    integer(ip), parameter :: arh_n_micro = 300

    abstract interface
        subroutine evaluate_dm_cs_type(dm, energy, fock, v_coulomb, v_exchange, &
                                       v_nonlinear, error)
            import :: rp, ip

            real(rp), intent(in), target, contiguous :: dm(:, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), optional, target, contiguous :: &
                fock(:, :), v_coulomb(:, :), v_exchange(:, :), v_nonlinear(:, :)
            integer(ip), intent(out) :: error
        end subroutine evaluate_dm_cs_type

        subroutine evaluate_dm_os_type(dm, energy, fock, v_coulomb, v_exchange, &
                                       v_nonlinear, error)
            import :: rp, ip

            real(rp), intent(in), target :: dm(:, :, :)
            real(rp), intent(out) :: energy
            real(rp), intent(out), optional, target :: &
                fock(:, :, :), v_coulomb(:, :, :), v_exchange(:, :, :), &
                v_nonlinear(:, :, :)
            integer(ip), intent(out) :: error
        end subroutine evaluate_dm_os_type
    end interface

    ! ARH object, holding the approximate Hessian model built from the history of the
    ! orbitals it points to; the operations on the orbital basis which only ARH needs
    ! are type-bound procedures of the types extending it for the MO and the OAO basis,
    ! which point to the quantities of the object holding the orbitals that these
    ! operations need
    type, abstract :: arh_type
        type(arh_settings_type) :: settings
        logical :: evaluation_stale = .true., model_stale = .false.
        class(orbital_basis_type), pointer :: orbitals => null()
        real(rp), allocatable :: &
            v_coulomb(:, :, :), v_exchange(:, :, :), v_nonlinear(:, :, :), &
            a_sym(:, :), a_sym_nonlinear(:, :), a_inv(:, :), a_inv_comb(:, :), &
            dm_list(:, :, :, :), v_coulomb_list(:, :, :, :), &
            v_exchange_list(:, :, :, :), v_nonlinear_list(:, :, :, :), &
            linear_potential_dirs(:, :), nonlinear_potential_dirs(:, :), &
            dm_dirs(:, :), dm_dirs_nonlinear(:, :), expansion_dirs(:, :), &
            projection_dirs(:, :), coupling_matrix(:, :)
        procedure(evaluate_dm_os_type), pointer, nopass :: evaluate_dm_os => null()
        procedure(evaluate_dm_cs_type), pointer, nopass :: evaluate_dm_cs => null()
    contains
        procedure(rotate_trial_arh_type), deferred :: rotate_trial
        procedure(history_dm_arh_type), deferred :: history_dm
        procedure(to_history_basis_arh_type), deferred :: to_history_basis
        procedure(history_columns_arh_type), deferred :: history_columns
        procedure(channel_rows_arh_type), deferred :: history_channel_rows
        procedure(channel_rows_arh_type), deferred :: packed_channel_rows
        procedure(hess_x_static_arh_type), deferred :: hess_x_static
    end type

    abstract interface
        subroutine rotate_trial_arh_type(self, kappa, rot_dm_ao, rot_dm_hist, error)
            import :: arh_type, rp, ip

            class(arh_type), intent(in) :: self
            real(rp), intent(in) :: kappa(:)
            real(rp), intent(out) :: rot_dm_ao(:, :, :), rot_dm_hist(:, :, :)
            integer(ip), intent(out) :: error
        end subroutine rotate_trial_arh_type

        function history_dm_arh_type(self) result(dm)
            import :: arh_type, rp

            class(arh_type), intent(in) :: self
            real(rp), allocatable :: dm(:, :, :)
        end function history_dm_arh_type

        function to_history_basis_arh_type(self, matrix_ao) result(matrix)
            import :: arh_type, rp

            class(arh_type), intent(in) :: self
            real(rp), intent(in) :: matrix_ao(:, :, :)
            real(rp), allocatable :: matrix(:, :, :)
        end function to_history_basis_arh_type

        subroutine history_columns_arh_type(self, diff, cols, packed)
            import :: arh_type, rp

            class(arh_type), intent(in) :: self
            real(rp), intent(in) :: diff(:, :, :, :)
            real(rp), intent(out), allocatable :: cols(:, :), packed(:, :)
        end subroutine history_columns_arh_type

        function channel_rows_arh_type(self) result(rows)
            import :: arh_type, ip

            class(arh_type), intent(in) :: self
            integer(ip), allocatable :: rows(:, :)
        end function channel_rows_arh_type

        function hess_x_static_arh_type(self, x) result(hess_x)
            import :: arh_type, rp

            class(arh_type), intent(in) :: self
            real(rp), intent(in) :: x(:)
            real(rp), allocatable :: hess_x(:)
        end function hess_x_static_arh_type
    end interface

    ! ARH object for the MO basis, whose quantities the MO object holds
    type, extends(arh_type) :: arh_mo_type
        real(rp), pointer, contiguous :: mo_coeff(:, :, :) => null()
        real(rp), pointer :: ao_overlap(:, :) => null()
        type(mo_channel_type), pointer :: mo_channels(:) => null()
    contains
        procedure :: rotate_trial => rotate_trial_arh_mo
        procedure :: history_dm => history_dm_arh_mo
        procedure :: to_history_basis => to_history_basis_arh_mo
        procedure :: history_columns => history_columns_arh_mo
        procedure :: history_channel_rows => history_channel_rows_arh_mo
        procedure :: packed_channel_rows => packed_channel_rows_arh_mo
        procedure :: hess_x_static => hess_x_static_arh_mo
    end type

    ! ARH object for the OAO basis, whose quantities the OAO object holds
    type, extends(arh_type) :: arh_oao_type
        real(rp), pointer :: s_inv_sqrt(:, :) => null(), dm_oao(:, :, :) => null(), &
                             fock_oo(:, :, :) => null(), fock_vv(:, :, :) => null()
    contains
        procedure :: rotate_trial => rotate_trial_arh_oao
        procedure :: history_dm => history_dm_arh_oao
        procedure :: to_history_basis => to_history_basis_arh_oao
        procedure :: history_columns => history_columns_arh_oao
        procedure :: history_channel_rows => history_channel_rows_arh_oao
        procedure :: packed_channel_rows => packed_channel_rows_arh_oao
        procedure :: hess_x_static => hess_x_static_arh_oao
    end type

    ! define the constructors for these types
    interface arh_mo_type
        module procedure construct_arh_mo
    end interface

    interface arh_oao_type
        module procedure construct_arh_oao
    end interface

    ! global variables
    class(arh_type), allocatable, target :: arh_object

    ! create function pointers to ensure that routines comply with interface
    procedure(obj_func_type), pointer :: obj_func_arh_callback_ptr => &
        obj_func_arh_callback
    procedure(update_orbs_type), pointer :: update_orbs_arh_callback_ptr => &
        update_orbs_arh_callback
    procedure(hess_x_type), pointer :: hess_x_arh_callback_ptr => hess_x_arh_callback
    procedure(precond_type), pointer :: precond_arh_callback_ptr => precond_arh_callback
    procedure(precond_pd_type), pointer :: precond_pd_arh_callback_ptr => &
        precond_pd_arh_callback
    procedure(get_extra_trial_vectors_type), pointer :: &
        get_extra_trial_vectors_arh_callback_ptr => get_extra_trial_vectors_arh_callback

    ! define module procedures for different spin cases in every orbital basis, which
    ! share their call structure with the C and Python interfaces
    interface arh_factory_mo
        module procedure arh_factory_mo_cs, arh_factory_mo_os
    end interface

    interface arh_factory_oao
        module procedure arh_factory_oao_cs, arh_factory_oao_os
    end interface

    ! define module procedures for different spin cases and orbital bases
    interface arh_factory
        module procedure arh_factory_mo_cs, arh_factory_mo_os, arh_factory_oao_cs, &
            arh_factory_oao_os
    end interface

contains

    subroutine arh_factory_mo_cs(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                                 evaluate_dm_cs, obj_func_arh_funptr, &
                                 update_orbs_arh_funptr, solver_settings, error, &
                                 settings, orbsym)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! closed-shell case, with the orbitals parameterized in the MO basis, and wires
        ! the ARH preconditioners and extra trial vectors into the solver settings,
        ! which it also asks to rebuild the Hessian information of the subsystem solver
        ! after rejected steps and tells whether the approximate Hessian is symmetric;
        ! the MO coefficients are rotated in place, so they have to outlive the
        ! calculation; if the irreps of the MOs are given, only the occupied-virtual
        ! rotations between orbitals of the same irrep are parameters
        !
        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ, n_particle, n_ao, n_mo
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_cs
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings
        integer(ip), intent(in), optional :: orbsym(:)

        real(rp), pointer, contiguous :: mo_coeff_3d(:, :, :)
        integer(ip), allocatable :: orbsym_2d(:, :)

        ! initialize error flag
        error = 0

        ! call common setup, passing the irreps of the MOs only if they are given
        mo_coeff_3d(1:size(mo_coeff, 1), 1:size(mo_coeff, 2), 1:1) => mo_coeff
        if (present(orbsym)) orbsym_2d = reshape(orbsym, [size(orbsym, kind=ip), 1_ip])
        call arh_factory_mo_common(mo_coeff_3d, ao_overlap, [n_occ], n_particle, n_ao, &
                                   n_mo, error, settings, orbsym_2d)
        if (error /= 0) return
        nullify(mo_coeff_3d)

        ! set pointers to functions
        arh_object%evaluate_dm_cs => evaluate_dm_cs

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh_callback
        update_orbs_arh_funptr => update_orbs_arh_callback

        ! wire the remaining ARH routines into the solver settings
        call arh_set_solver_settings(solver_settings, arh_object, error)

    end subroutine arh_factory_mo_cs

    subroutine arh_factory_mo_os(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                                 evaluate_dm_os, obj_func_arh_funptr, &
                                 update_orbs_arh_funptr, solver_settings, error, &
                                 settings, orbsym)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! open-shell case, with the orbitals parameterized in the MO basis, and wires
        ! the ARH preconditioners and extra trial vectors into the solver settings,
        ! which it also asks to rebuild the Hessian information of the subsystem solver
        ! after rejected steps and tells whether the approximate Hessian is symmetric;
        ! the MO coefficients are rotated in place, so they have to outlive the
        ! calculation; if the irreps of the MOs of every particle channel are given,
        ! only the occupied-virtual rotations between orbitals of the same irrep are
        ! parameters
        !
        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings
        integer(ip), intent(in), optional :: orbsym(:, :)

        ! initialize error flag
        error = 0

        ! call common setup
        call arh_factory_mo_common(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, &
                                   n_mo, error, settings, orbsym)
        if (error /= 0) return

        ! set pointers to functions
        arh_object%evaluate_dm_os => evaluate_dm_os

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh_callback
        update_orbs_arh_funptr => update_orbs_arh_callback

        ! wire the remaining ARH routines into the solver settings
        call arh_set_solver_settings(solver_settings, arh_object, error)

    end subroutine arh_factory_mo_os

    subroutine arh_factory_mo_common(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, &
                                     n_mo, error, settings, orbsym)
        !
        ! this subroutine performs common ARH initialization operations for orbitals
        ! parameterized in the MO basis
        !
        use otr_mo, only: mo_factory_common, mo_object

        real(rp), intent(inout), target, contiguous :: mo_coeff(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_occ(:), n_particle, n_ao, n_mo
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings
        integer(ip), intent(in), optional :: orbsym(:, :)

        ! perform ARH sanity check
        call arh_sanity_check(settings, error)
        if (error /= 0) return

        ! call common MO setup, which performs the sanity check of the MO basis
        call mo_factory_common(mo_coeff, ao_overlap, n_occ, n_particle, n_ao, n_mo, &
                               error, settings, orbsym)
        if (error /= 0) return

        ! discard any state (in particular history and derived quantities) from a
        ! previous calculation
        if (allocated(arh_object)) deallocate(arh_object)

        ! parameterize the orbitals in the MO basis
        allocate(arh_object, source=arh_mo_type(mo_object))

        ! set (potentially new) settings
        arh_object%settings = settings

    end subroutine arh_factory_mo_common

    subroutine arh_factory_oao_cs(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_cs, &
                                  obj_func_arh_funptr, update_orbs_arh_funptr, &
                                  solver_settings, error, settings)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! closed-shell case, with the orbitals parameterized in the OAO basis, and
        ! wires the ARH preconditioners, projection and extra trial vectors into the
        ! solver settings, which it also asks to rebuild the Hessian information of the
        ! subsystem solver after rejected steps and tells whether the approximate
        ! Hessian is symmetric; the density matrix is updated in place, so it has to
        ! outlive the calculation
        !
        use otr_oao, only: project_oao_callback

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_cs_type), intent(in), pointer :: evaluate_dm_cs
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        real(rp), pointer, contiguous :: dm_ao_3d(:, :, :)

        ! initialize error flag
        error = 0

        ! call common setup
        dm_ao_3d(1:n_ao, 1:n_ao, 1:1) => dm_ao
        call arh_factory_oao_common(dm_ao_3d, ao_overlap, n_particle, n_ao, error, &
                                    settings)
        if (error /= 0) return
        nullify(dm_ao_3d)

        ! set pointers to functions
        arh_object%evaluate_dm_cs => evaluate_dm_cs

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh_callback
        update_orbs_arh_funptr => update_orbs_arh_callback

        ! wire the remaining ARH routines and the projection of the OAO basis into the
        ! solver settings
        call arh_set_solver_settings(solver_settings, arh_object, error, &
                                     project_oao_callback)

    end subroutine arh_factory_oao_cs

    subroutine arh_factory_oao_os(dm_ao, ao_overlap, n_particle, n_ao, evaluate_dm_os, &
                                  obj_func_arh_funptr, update_orbs_arh_funptr, &
                                  solver_settings, error, settings)
        !
        ! this function returns a modified ARH orbital updating function for the
        ! open-shell case, with the orbitals parameterized in the OAO basis, and wires
        ! the ARH preconditioners, projection and extra trial vectors into the solver
        ! settings, which it also asks to rebuild the Hessian information of the
        ! subsystem solver after rejected steps and tells whether the approximate
        ! Hessian is symmetric; the density matrix is updated in place, so it has to
        ! outlive the calculation
        !
        use otr_oao, only: project_oao_callback

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        procedure(evaluate_dm_os_type), intent(in), pointer :: evaluate_dm_os
        procedure(obj_func_type), intent(out), pointer :: obj_func_arh_funptr
        procedure(update_orbs_type), intent(out), pointer :: update_orbs_arh_funptr
        type(solver_settings_type), intent(inout) :: solver_settings
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! initialize error flag
        error = 0

        ! call common setup
        call arh_factory_oao_common(dm_ao, ao_overlap, n_particle, n_ao, error, &
                                    settings)
        if (error /= 0) return

        ! set pointers to functions
        arh_object%evaluate_dm_os => evaluate_dm_os

        ! get pointers to modified function
        obj_func_arh_funptr => obj_func_arh_callback
        update_orbs_arh_funptr => update_orbs_arh_callback

        ! wire the remaining ARH routines and the projection of the OAO basis into the
        ! solver settings
        call arh_set_solver_settings(solver_settings, arh_object, error, &
                                     project_oao_callback)

    end subroutine arh_factory_oao_os

    subroutine arh_factory_oao_common(dm_ao, ao_overlap, n_particle, n_ao, error, &
                                      settings)
        !
        ! this subroutine performs common ARH initialization operations for orbitals
        ! parameterized in the OAO basis
        !
        use otr_oao, only: oao_factory_common, oao_object

        real(rp), intent(inout), target, contiguous :: dm_ao(:, :, :)
        real(rp), intent(in) :: ao_overlap(:, :)
        integer(ip), intent(in) :: n_particle, n_ao
        integer(ip), intent(out) :: error
        type(arh_settings_type), intent(inout) :: settings

        ! perform ARH sanity check
        call arh_sanity_check(settings, error)
        if (error /= 0) return

        ! call common OAO setup, which performs the sanity check of the OAO basis
        call oao_factory_common(dm_ao, ao_overlap, n_particle, n_ao, error, settings)
        if (error /= 0) return

        ! discard any state (in particular history and derived quantities) from a
        ! previous calculation: the ARH factories are only ever called to start a new
        ! calculation, never to resume one, so state from an unrelated trajectory must
        ! never be reused, regardless of whether the new dimensions happen to match the
        ! old ones
        if (allocated(arh_object)) deallocate(arh_object)

        ! parameterize the orbitals in the OAO basis
        allocate(arh_object, source=arh_oao_type(oao_object))

        ! set (potentially new) settings
        arh_object%settings = settings

    end subroutine arh_factory_oao_common

    subroutine arh_sanity_check(settings, error)
        !
        ! this subroutine performs a sanity check for ARH input parameters
        !
        use opentrustregion, only: verbosity_error, string_to_lowercase

        type(arh_settings_type), intent(inout) :: settings
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! convert strings to lowercase
        settings%arh_type = string_to_lowercase(settings%arh_type)

        ! check for character options
        if (.not. any(settings%arh_type == arh_types)) then
            call settings%log("ARH type option unknown. Possible values are "// &
                              """arh"" (standard ARH), ""symm_arh"" (symmetrized "// &
                              "ARH), ""ms_sp"" (subspace-projected multisecant "// &
                              "method), ""ms_psb"" (multisecant PSB), and "// &
                              """ms_sr1"" (multisecant SR1 version).", &
                              verbosity_error, .true.)
            error = 1
            return
        end if

    end subroutine arh_sanity_check

    subroutine arh_set_solver_settings(solver_settings, arh, error, project)
        !
        ! this subroutine wires the ARH preconditioners and extra trial vectors and, if
        ! given, the projection of the orbital basis into the solver settings and those
        ! of its stability check, which only the OAO basis passes, since the parameters
        ! of the MO basis are non-redundant, so that the MO basis leaves the projection
        ! to the caller; it also asks the solver to rebuild the Hessian information of
        ! the subsystem solver after rejected steps, since the objective function adds
        ! the rejected point to the history, which changes the approximate Hessian,
        ! raises the micro iteration limit, since the linear transformations of the
        ! approximate Hessian are cheap, and tells it whether the approximate Hessian
        ! of the given ARH type is symmetric
        !
        use opentrustregion, only: project_type

        type(solver_settings_type), intent(inout) :: solver_settings
        class(arh_type), intent(in) :: arh
        integer(ip), intent(out) :: error
        procedure(project_type), optional :: project

        ! initialize error flag
        error = 0

        ! initialize settings
        if (.not. solver_settings%initialized) then
            call solver_settings%init(error)
            if (error /= 0) return
        end if

        ! set callback functions of the solver and its stability check
        solver_settings%precond => precond_arh_callback
        solver_settings%precond_pd => precond_pd_arh_callback
        solver_settings%get_extra_trial_vectors => get_extra_trial_vectors_arh_callback
        solver_settings%stability_settings%precond => precond_arh_callback
        solver_settings%stability_settings%get_extra_trial_vectors => &
            get_extra_trial_vectors_arh_callback
        if (present(project)) then
            solver_settings%project => project
            solver_settings%stability_settings%project => project
        end if

        ! rebuild the Hessian information after rejected steps
        solver_settings%refresh_hess = .true.

        ! converge the subsystem solver with the cheap approximate Hessian
        solver_settings%n_micro = arh_n_micro

        ! only the standard ARH type breaks the symmetry of the approximate Hessian
        solver_settings%hess_symm = arh%settings%arh_type /= "arh"

    end subroutine arh_set_solver_settings

    function obj_func_arh_callback(kappa, error) result(energy)
        !
        ! this function defines the energy evaluation, which also adds the evaluated
        ! point with its Coulomb, exact-exchange and non-linear potentials to the
        ! history, unless it is already there
        !
        real(rp), intent(in), target :: kappa(:)
        integer(ip), intent(out) :: error
        real(rp) :: energy

        integer(ip) :: n_ao, n_particle
        real(rp), allocatable :: rot_dm_ao(:, :, :), rot_dm_hist(:, :, :), &
                                 v_coulomb_ao(:, :, :), v_exchange_ao(:, :, :), &
                                 v_nonlinear_ao(:, :, :)

        ! initialize energy in case of error
        energy = 0.0_rp

        ! number of AOs and number of particles
        n_ao = arh_object%orbitals%n_ao
        n_particle = arh_object%orbitals%n_particle

        ! get rotated density matrix in the AO and history basis
        allocate(rot_dm_ao(n_ao, n_ao, n_particle), &
                 rot_dm_hist(n_ao, n_ao, n_particle), &
                 v_coulomb_ao(n_ao, n_ao, n_particle), &
                 v_exchange_ao(n_ao, n_ao, n_particle), &
                 v_nonlinear_ao(n_ao, n_ao, n_particle))
        call arh_object%rotate_trial(kappa, rot_dm_ao, rot_dm_hist, error)
        if (error /= 0) return

        ! calculate mean-field energy
        if (associated(arh_object%evaluate_dm_cs)) then
            call arh_object%evaluate_dm_cs( &
                rot_dm_ao(:, :, 1), energy, v_coulomb=v_coulomb_ao(:, :, 1), &
                v_exchange=v_exchange_ao(:, :, 1), &
                v_nonlinear=v_nonlinear_ao(:, :, 1), error=error)
        else
            call arh_object%evaluate_dm_os(rot_dm_ao, energy, v_coulomb=v_coulomb_ao, &
                                           v_exchange=v_exchange_ao, &
                                           v_nonlinear=v_nonlinear_ao, error=error)
        end if
        if (error /= 0) return

        ! update list of density and potential matrices
        if (allocated(arh_object%dm_list)) then
            if (.not. density_in_history(rot_dm_hist)) then
                call prepend(arh_object%dm_list, rot_dm_hist)
                call prepend(arh_object%v_coulomb_list, &
                             arh_object%to_history_basis(v_coulomb_ao))
                call prepend(arh_object%v_exchange_list, &
                             arh_object%to_history_basis(v_exchange_ao))
                call prepend(arh_object%v_nonlinear_list, &
                             arh_object%to_history_basis(v_nonlinear_ao))
                arh_object%model_stale = .true.
            end if
        end if

    end function obj_func_arh_callback

    subroutine update_orbs_arh_callback(kappa, func, grad, h_diag, hess_x_funptr, error)
        !
        ! this function defines the energy, gradient, and Hessian diagonal evaluation
        ! and the Hessian linear transformation on the basis of augmented Roothaan-Hall
        !
        real(rp), intent(in), target :: kappa(:)
        real(rp), intent(out) :: func
        real(rp), intent(out), target :: grad(:), h_diag(:)
        procedure(hess_x_type), intent(out), pointer :: hess_x_funptr
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle
        real(rp), allocatable :: dm_current(:, :, :), fock_ao(:, :, :), &
                                 v_coulomb_ao(:, :, :), v_exchange_ao(:, :, :), &
                                 v_nonlinear_ao(:, :, :)

        ! initialize error flag
        error = 0

        ! evaluate at the rotated orbitals, or at the current ones if these have not
        ! been evaluated yet
        if ((sum(abs(kappa)) > 0.0_rp) .or. arh_object%evaluation_stale) then
            ! number of AOs
            n_ao = arh_object%orbitals%n_ao

            ! number of particles
            n_particle = arh_object%orbitals%n_particle

            ! add the current point to the list of density and potential matrices if it
            ! has been evaluated
            if (.not. allocated(arh_object%dm_list)) then
                allocate(arh_object%dm_list(n_ao, n_ao, n_particle, 0), &
                         arh_object%v_coulomb_list(n_ao, n_ao, n_particle, 0), &
                         arh_object%v_exchange_list(n_ao, n_ao, n_particle, 0), &
                         arh_object%v_nonlinear_list(n_ao, n_ao, n_particle, 0))
            else if (.not. arh_object%evaluation_stale) then
                dm_current = arh_object%history_dm()
                if (.not. density_in_history(dm_current)) then
                    call prepend(arh_object%dm_list, dm_current)
                    call prepend(arh_object%v_coulomb_list, arh_object%v_coulomb)
                    call prepend(arh_object%v_exchange_list, arh_object%v_exchange)
                    call prepend(arh_object%v_nonlinear_list, arh_object%v_nonlinear)
                end if
            end if

            ! rotate orbitals, which have not been evaluated until this succeeds
            arh_object%evaluation_stale = .true.
            call arh_object%orbitals%rotate_orbitals(kappa, arh_object%settings, error)
            if (error /= 0) return

            ! get energy, Fock matrix, Coulomb and exact-exchange potentials, and
            ! non-linear potential
            allocate(fock_ao(n_ao, n_ao, n_particle), &
                     v_coulomb_ao(n_ao, n_ao, n_particle), &
                     v_exchange_ao(n_ao, n_ao, n_particle), &
                     v_nonlinear_ao(n_ao, n_ao, n_particle))
            if (associated(arh_object%evaluate_dm_cs)) then
                call arh_object%evaluate_dm_cs( &
                    arh_object%orbitals%dm_ao(:, :, 1), arh_object%orbitals%energy, &
                    fock_ao(:, :, 1), v_coulomb_ao(:, :, 1), v_exchange_ao(:, :, 1), &
                    v_nonlinear_ao(:, :, 1), error)
            else
                call arh_object%evaluate_dm_os( &
                    arh_object%orbitals%dm_ao, arh_object%orbitals%energy, fock_ao, &
                    v_coulomb_ao, v_exchange_ao, v_nonlinear_ao, error)
            end if
            if (error /= 0) then
                deallocate(fock_ao, v_coulomb_ao, v_exchange_ao, v_nonlinear_ao)
                return
            end if

            ! calculate gradient, Hessian diagonal and static part of the Hessian from
            ! the Fock matrix in the history basis, the basis the orbital object
            ! represents the density matrix in
            call arh_object%orbitals%calculate_grad_h_diag( &
                arh_object%to_history_basis(fock_ao))
            deallocate(fock_ao)

            ! transform Coulomb, exact-exchange and non-linear potentials to the
            ! history basis
            arh_object%v_coulomb = arh_object%to_history_basis(v_coulomb_ao)
            arh_object%v_exchange = arh_object%to_history_basis(v_exchange_ao)
            arh_object%v_nonlinear = arh_object%to_history_basis(v_nonlinear_ao)
            deallocate(v_coulomb_ao, v_exchange_ao, v_nonlinear_ao)

            ! assemble the approximate Hessian model from the history for the closed- or
            ! open-shell case
            if (associated(arh_object%evaluate_dm_cs)) then
                call build_hess_model_cs(error)
            else
                call build_hess_model_os(error)
            end if
            if (error /= 0) return

            ! the rotated orbitals have been evaluated
            arh_object%evaluation_stale = .false.
        end if

        ! set outputs
        func = arh_object%orbitals%energy
        grad = arh_object%orbitals%grad
        h_diag = arh_object%orbitals%h_diag
        hess_x_funptr => hess_x_arh_callback

    end subroutine update_orbs_arh_callback

    subroutine build_hess_model_cs(error)
        !
        ! this subroutine assembles the closed-shell approximate Hessian model from the
        ! history relative to the current point
        !
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list, n_acc, n_acc_nonlinear
        real(rp) :: min_residual
        real(rp), allocatable :: &
            dm_current(:, :, :), dm_diff(:, :, :, :), v_coulomb_diff(:, :, :, :), &
            v_exchange_diff(:, :, :, :), v_nonlinear_diff(:, :, :, :), dm_cols(:, :), &
            v_coulomb_cols(:, :), v_exchange_cols(:, :), v_linear_cols(:, :), &
            v_nonlinear_cols(:, :), dm_packed(:, :), v_coulomb_packed(:, :), &
            v_exchange_packed(:, :), v_linear_packed(:, :), v_nonlinear_packed(:, :), &
            chol(:, :), chol_nonlinear(:, :)
        integer(ip), allocatable :: map(:), map_nonlinear(:)
        logical, allocatable :: keep_nonlinear(:)

        ! initialize error flag
        error = 0

        ! number of AOs and number of particles
        n_ao = arh_object%orbitals%n_ao
        n_particle = arh_object%orbitals%n_particle

        ! prepare the density matrix and potential differences: the linear potential is
        ! split into its Coulomb and exact-exchange parts, whose difference it is, and
        ! the non-linear (XC) potential is kept apart
        n_list = size(arh_object%dm_list, 4)
        dm_current = arh_object%history_dm()
        allocate(dm_diff(n_ao, n_ao, n_particle, n_list), &
                 v_coulomb_diff(n_ao, n_ao, n_particle, n_list), &
                 v_exchange_diff(n_ao, n_ao, n_particle, n_list), &
                 v_nonlinear_diff(n_ao, n_ao, n_particle, n_list))
        do i = 1, n_list
            dm_diff(:, :, :, i) = arh_object%dm_list(:, :, :, i) - dm_current
            v_coulomb_diff(:, :, :, i) = arh_object%v_coulomb_list(:, :, :, i) - &
                                         arh_object%v_coulomb
            v_exchange_diff(:, :, :, i) = arh_object%v_exchange_list(:, :, :, i) - &
                                          arh_object%v_exchange
            v_nonlinear_diff(:, :, :, i) = arh_object%v_nonlinear_list(:, :, :, i) - &
                                           arh_object%v_nonlinear
        end do

        ! express the differences as history and packed columns, and the linear
        ! potential differences as the difference of their Coulomb and exact-exchange
        ! parts, which the transformation to columns preserves
        call arh_object%history_columns(dm_diff, dm_cols, dm_packed)
        call arh_object%history_columns(v_coulomb_diff, v_coulomb_cols, &
                                        v_coulomb_packed)
        call arh_object%history_columns(v_exchange_diff, v_exchange_cols, &
                                        v_exchange_packed)
        call arh_object%history_columns(v_nonlinear_diff, v_nonlinear_cols, &
                                        v_nonlinear_packed)
        deallocate(dm_diff, v_coulomb_diff, v_exchange_diff, v_nonlinear_diff)
        v_linear_cols = v_coulomb_cols - v_exchange_cols
        v_linear_packed = v_coulomb_packed - v_exchange_packed

        ! factorize the density-matrix-difference history for the linear part which
        ! resolves linear dependencies in the history; the coupling formulas below
        ! never need an explicit metric inverse since everything is expressed in the
        ! resulting orthonormalized basis via a single triangular solve
        call factorize_history(dm_cols, chol, map, n_acc)

        ! factorize the same history for the non-linear part, which can only be done
        ! once its response is known since the error in that response sets the shortest
        ! residual a direction has to contribute
        keep_nonlinear = history_step_mask(dm_cols)
        min_residual = resolvable_residual(dm_cols, v_nonlinear_cols, keep_nonlinear)
        call factorize_history(dm_cols, chol_nonlinear, map_nonlinear, &
                               n_acc_nonlinear, keep_nonlinear, min_residual)

        ! MS-SR1
        if (arh_object%settings%arh_type == "ms_sr1") then
            ! linear part: the difference of a Coulomb and an exact-exchange
            ! multisecant term, each exact since Coulomb and exact exchange are linear
            ! in the density matrix, and each bounded by its kernel, which is positive
            ! semidefinite, unlike the linear kernel itself; this also caches the
            ! packed potential difference directions the low-rank Hessian factors are
            ! assembled from, rebased into the orthonormalized S-basis
            call get_ms_a_inv_jk_cs( &
                dm_cols, v_coulomb_cols, v_exchange_cols, v_coulomb_packed, &
                v_exchange_packed, map, chol, arh_object%a_inv, &
                arh_object%linear_potential_dirs, arh_object%settings, error)
            if (error /= 0) return
            ! non-linear part: kept on its own system so that it contracts against its
            ! own inverse and the exact linear secant relationship is not averaged with
            ! the approximate one
            call get_ms_a_inv(dm_cols, v_nonlinear_cols, map_nonlinear, &
                              chol_nonlinear, arh_object%a_inv_comb, &
                              arh_object%settings, error)
            if (error /= 0) return

            ! cache the packed non-linear potential difference directions, rebased into
            ! the orthonormalized S-basis
            arh_object%nonlinear_potential_dirs = &
                rebase_dirs(v_nonlinear_packed, map_nonlinear, chol_nonlinear)
        ! ARH and related methods
        else
            ! cache the packed history-projection directions the low-rank Hessian
            ! factors are assembled from, rebased into the orthonormalized S-basis
            ! belonging to the linear and non-linear systems
            arh_object%dm_dirs = rebase_dirs(dm_packed, map, chol)
            arh_object%dm_dirs_nonlinear = &
                rebase_dirs(dm_packed, map_nonlinear, chol_nonlinear)
            if (arh_object%settings%arh_type /= "ms_sp") then
                arh_object%linear_potential_dirs = &
                    rebase_dirs(v_linear_packed, map, chol)
                arh_object%nonlinear_potential_dirs = &
                    rebase_dirs(v_nonlinear_packed, map_nonlinear, chol_nonlinear)
            end if

            ! construct A = S^T Y for the linear and non-linear system,
            ! congruence-transformed into its own orthonormalized S-basis
            if (arh_object%settings%arh_type == "ms_sp" .or. &
                arh_object%settings%arh_type == "ms_psb") then
                arh_object%a_sym = &
                    build_a_transformed(dm_cols, v_linear_cols, map, chol)
                arh_object%a_sym_nonlinear = build_a_transformed( &
                    dm_cols, v_nonlinear_cols, map_nonlinear, chol_nonlinear)
            end if
        end if

        ! assemble the low-rank (response) part of the approximate Hessian
        call get_low_rank_hess_factors()

        ! the model now reflects the whole history
        arh_object%model_stale = .false.

    end subroutine build_hess_model_cs

    subroutine build_hess_model_os(error)
        !
        ! this subroutine assembles the open-shell approximate Hessian model from the
        ! history relative to the current point
        !
        integer(ip), intent(out) :: error

        integer(ip) :: n_ao, n_particle, i, n_list, n_acc1, n_acc2, n_acc_nl, &
                       n_acc1_nl, n_acc2_nl
        real(rp) :: min_residual
        integer(ip), allocatable :: map1(:), map2(:), map_comb(:), map_nl(:), &
                                    map1_nl(:), map2_nl(:), map_comb_nl(:), &
                                    history_rows(:, :), packed_rows(:, :)
        logical, allocatable :: keep_nonlinear(:)
        real(rp), allocatable :: &
            dm_current(:, :, :), dm_diff(:, :, :, :), v_coulomb_diff(:, :, :, :), &
            v_exchange_diff(:, :, :, :), v_nonlinear_diff(:, :, :, :), &
            v_coulomb_alpha_diff(:, :, :, :), v_coulomb_beta_diff(:, :, :, :), &
            v_same_spin_diff(:, :, :, :), v_opposite_spin_diff(:, :, :, :), &
            dm_cols(:, :), v_coulomb_alpha_cols(:, :), v_coulomb_beta_cols(:, :), &
            v_exchange_cols(:, :), v_same_spin_cols(:, :), v_opposite_spin_cols(:, :), &
            v_nonlinear_cols(:, :), dm_cols1(:, :), dm_cols2(:, :), &
            v_nonlinear_cols1(:, :), v_nonlinear_cols2(:, :), dm_packed(:, :), &
            v_coulomb_alpha_packed(:, :), v_coulomb_beta_packed(:, :), &
            v_exchange_packed(:, :), v_same_spin_packed(:, :), &
            v_opposite_spin_packed(:, :), v_nonlinear_packed(:, :), chol1(:, :), &
            chol2(:, :), chol1_nl(:, :), chol2_nl(:, :), chol_comb(:, :), &
            chol_comb_nl(:, :), chol_nl(:, :)

        ! initialize error flag
        error = 0

        ! number of AOs and number of particles
        n_ao = arh_object%orbitals%n_ao
        n_particle = arh_object%orbitals%n_particle

        ! prepare the density matrix and potential differences; the Coulomb and
        ! exact-exchange potential differences are kept separate from the non-linear one
        ! since the former are exact at any distance from the current density and
        ! therefore take the full history, while the non-linear one describes a
        ! drifting Hessian and is fitted only to the history entries the screening
        ! below admits
        n_list = size(arh_object%dm_list, 4)
        dm_current = arh_object%history_dm()
        allocate(dm_diff(n_ao, n_ao, n_particle, n_list), &
                 v_coulomb_diff(n_ao, n_ao, n_particle, n_list), &
                 v_exchange_diff(n_ao, n_ao, n_particle, n_list), &
                 v_nonlinear_diff(n_ao, n_ao, n_particle, n_list))
        do i = 1, n_list
            dm_diff(:, :, :, i) = arh_object%dm_list(:, :, :, i) - dm_current
            v_coulomb_diff(:, :, :, i) = arh_object%v_coulomb_list(:, :, :, i) - &
                                         arh_object%v_coulomb
            v_exchange_diff(:, :, :, i) = arh_object%v_exchange_list(:, :, :, i) - &
                                          arh_object%v_exchange
            v_nonlinear_diff(:, :, :, i) = arh_object%v_nonlinear_list(:, :, :, i) - &
                                           arh_object%v_nonlinear
        end do

        ! express the density matrix and non-linear potential differences as history
        ! and packed columns
        call arh_object%history_columns(dm_diff, dm_cols, dm_packed)
        call arh_object%history_columns(v_nonlinear_diff, v_nonlinear_cols, &
                                        v_nonlinear_packed)
        deallocate(dm_diff, v_nonlinear_diff)

        ! express the linear potential differences as history and packed columns:
        ! multisecant SR1 takes the Coulomb potential of the alpha and of the beta
        ! density differences, which acts on both channels, and the exact-exchange
        ! potential of every channel's density difference, which acts on that channel
        ! only, while the ARH family takes the same-spin potential of every channel,
        ! its own Coulomb minus exact-exchange potential, and the opposite-spin
        ! potential, the Coulomb potential of the other channel; the history basis is
        ! shared by both channels, so that a potential can be moved between them
        if (arh_object%settings%arh_type == "ms_sr1") then
            v_coulomb_alpha_diff = spread(v_coulomb_diff(:, :, 1, :), 3, n_particle)
            v_coulomb_beta_diff = spread(v_coulomb_diff(:, :, 2, :), 3, n_particle)
            call arh_object%history_columns( &
                v_coulomb_alpha_diff, v_coulomb_alpha_cols, v_coulomb_alpha_packed)
            call arh_object%history_columns(v_coulomb_beta_diff, v_coulomb_beta_cols, &
                                            v_coulomb_beta_packed)
            call arh_object%history_columns(v_exchange_diff, v_exchange_cols, &
                                            v_exchange_packed)
            deallocate(v_coulomb_alpha_diff, v_coulomb_beta_diff)
        else
            v_same_spin_diff = v_coulomb_diff - v_exchange_diff
            v_opposite_spin_diff = v_coulomb_diff(:, :, [2, 1], :)
            call arh_object%history_columns(v_same_spin_diff, v_same_spin_cols, &
                                            v_same_spin_packed)
            call arh_object%history_columns( &
                v_opposite_spin_diff, v_opposite_spin_cols, v_opposite_spin_packed)
            deallocate(v_same_spin_diff, v_opposite_spin_diff)
        end if
        deallocate(v_coulomb_diff, v_exchange_diff)
        history_rows = arh_object%history_channel_rows()
        packed_rows = arh_object%packed_channel_rows()

        ! copy the rows of the individual spin channels for the density matrix
        ! differences, which are factorized per channel
        dm_cols1 = dm_cols(history_rows(1, 1):history_rows(2, 1), :)
        dm_cols2 = dm_cols(history_rows(1, 2):history_rows(2, 2), :)

        ! factorize the density-matrix-difference history once per channel which
        ! resolves linear dependencies in the history; the coupling formulas below
        ! never need an explicit metric inverse since everything is expressed in the
        ! resulting orthonormalized basis via a single triangular solve; the two
        ! channels never mix in this factorization (the metric itself is block-diagonal
        ! across channels), so they are factorized independently and then combined into
        ! one block-diagonal (chol, map) pair wherever the coupled types or the linear
        ! multisecant system need to act on both channels together
        call factorize_history(dm_cols1, chol1, map1, n_acc1)
        call factorize_history(dm_cols2, chol2, map2, n_acc2)
        call combine_channels(chol1, map1, chol2, map2, n_list, chol_comb, map_comb)

        ! MS-SR1
        if (arh_object%settings%arh_type == "ms_sr1") then
            ! linear part: spin-separated, which is exact since Coulomb and exact
            ! exchange are linear in the density matrix, and the difference of a
            ! Coulomb and an exact-exchange multisecant term, each bounded by its
            ! kernel, which is positive semidefinite, unlike the linear kernel itself;
            ! this also caches the packed potential difference directions the low-rank
            ! Hessian factors are assembled from, rebased into the orthonormalized
            ! S-basis
            call get_ms_a_inv_jk_os( &
                dm_cols, v_coulomb_alpha_cols, v_coulomb_beta_cols, v_exchange_cols, &
                v_coulomb_alpha_packed, v_coulomb_beta_packed, v_exchange_packed, &
                history_rows, packed_rows, map_comb, chol_comb, arh_object%a_inv, &
                arh_object%linear_potential_dirs, arh_object%settings, error)
            if (error /= 0) return
            ! non-linear part: get spin-combined multisecant SR1 matrix; the non-linear
            ! response mixes both channels at once, so this needs its own, separate
            ! combined-flat factorization of the same history
            keep_nonlinear = history_step_mask(dm_cols)
            min_residual = &
                resolvable_residual(dm_cols, v_nonlinear_cols, keep_nonlinear)
            call factorize_history(dm_cols, chol_nl, map_nl, n_acc_nl, keep_nonlinear, &
                                   min_residual)
            call get_ms_a_inv(dm_cols, v_nonlinear_cols, map_nl, chol_nl, &
                              arh_object%a_inv_comb, arh_object%settings, error)
            if (error /= 0) return

            ! cache the packed non-linear potential difference directions, which need
            ! no channel-splitting, rebased into the orthonormalized S-basis
            arh_object%nonlinear_potential_dirs = &
                rebase_dirs(v_nonlinear_packed, map_nl, chol_nl)
        ! ARH and related methods
        else
            ! screened per-channel factorization for the non-linear system; both
            ! channels are screened on the step length of the whole density matrix
            ! rather than on their own channel's, since the non-linear potential of
            ! either channel is a functional of both spin densities and its staleness
            ! is therefore set by the total step
            v_nonlinear_cols1 = &
                v_nonlinear_cols(history_rows(1, 1):history_rows(2, 1), :)
            v_nonlinear_cols2 = &
                v_nonlinear_cols(history_rows(1, 2):history_rows(2, 2), :)
            keep_nonlinear = history_step_mask(dm_cols)
            min_residual = &
                resolvable_residual(dm_cols1, v_nonlinear_cols1, keep_nonlinear)
            call factorize_history(dm_cols1, chol1_nl, map1_nl, n_acc1_nl, &
                                   keep_nonlinear, min_residual)
            min_residual = &
                resolvable_residual(dm_cols2, v_nonlinear_cols2, keep_nonlinear)
            call factorize_history(dm_cols2, chol2_nl, map2_nl, n_acc2_nl, &
                                   keep_nonlinear, min_residual)
            call combine_channels(chol1_nl, map1_nl, chol2_nl, map2_nl, n_list, &
                                  chol_comb_nl, map_comb_nl)

            ! cache the packed history-projection directions the low-rank Hessian
            ! factors are assembled from, per channel for the density matrix and the
            ! non-linear potential, which has no opposite-spin counterpart, and
            ! channel-combined for the linear potential, each rebased into the S-basis
            ! of its own system
            call cache_channel_split_dirs(dm_packed, packed_rows, map_comb, chol_comb, &
                                          arh_object%dm_dirs)
            call cache_channel_split_dirs(dm_packed, packed_rows, map_comb_nl, &
                                          chol_comb_nl, arh_object%dm_dirs_nonlinear)
            if (arh_object%settings%arh_type /= "ms_sp") then
                call cache_combined_channel_dirs( &
                    v_same_spin_packed, v_opposite_spin_packed, packed_rows, map_comb, &
                    chol_comb, arh_object%linear_potential_dirs)
                call cache_channel_split_dirs(v_nonlinear_packed, packed_rows, &
                                              map_comb_nl, chol_comb_nl, &
                                              arh_object%nonlinear_potential_dirs)
            end if

            ! construct A = S^T Y for the linear and non-linear system,
            ! congruence-transformed into its own combined orthonormalized S-basis
            if (arh_object%settings%arh_type == "ms_sp" .or. &
                arh_object%settings%arh_type == "ms_psb") then
                arh_object%a_sym = build_a_block_linear_os( &
                    dm_cols, v_same_spin_cols, v_opposite_spin_cols, history_rows, &
                    map_comb, chol_comb)
                arh_object%a_sym_nonlinear = build_a_block_nonlinear_os( &
                    dm_cols, v_nonlinear_cols, history_rows, map_comb_nl, chol_comb_nl)
            end if
        end if

        ! assemble the low-rank (response) part of the approximate Hessian
        call get_low_rank_hess_factors()

        ! the model now reflects the whole history
        arh_object%model_stale = .false.

    end subroutine build_hess_model_os

    subroutine rebuild_stale_hess_model(error)
        !
        ! this subroutine rebuilds the approximate Hessian model if the objective
        ! function has added points to the history since it was last assembled, so that
        ! the Hessian linear transformation and the preconditioner always use the whole
        ! history
        !
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        ! nothing to rebuild if the model reflects the whole history
        if (.not. arh_object%model_stale) return

        ! rebuild model for the closed- or open-shell case
        if (associated(arh_object%evaluate_dm_cs)) then
            call build_hess_model_cs(error)
        else
            call build_hess_model_os(error)
        end if

    end subroutine rebuild_stale_hess_model

    subroutine hess_x_arh_callback(x, hess_x, error)
        !
        ! this function defines the Hessian linear transformation on the basis of
        ! augmented Roothaan-Hall and related methods
        !
        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: hess_x(:)
        integer(ip), intent(out) :: error

        integer(ip) :: n_dirs
        real(rp), allocatable :: projected_x(:), coupled_x(:)

        external :: dgemv

        ! initialize error flag
        error = 0

        ! rebuild the approximate Hessian model if it is stale
        call rebuild_stale_hess_model(error)
        if (error /= 0) return

        ! get static part
        hess_x = arh_object%hess_x_static(x)

        ! add the response part, which the shared low-rank factors express directly in
        ! the packed parameter space as expansion * coupling * projection^T, already
        ! carrying both the projection onto the non-redundant subspace (trivial in the
        ! MO basis) and the scaling applied to the static part
        if (allocated(arh_object%coupling_matrix)) then
            n_dirs = size(arh_object%expansion_dirs, 2)
            allocate(projected_x(n_dirs), coupled_x(n_dirs))
            call dgemv("T", size(x, kind=ip), n_dirs, 1.0_rp, &
                       arh_object%projection_dirs, size(x, kind=ip), x, 1_ip, 0.0_rp, &
                       projected_x, 1_ip)
            call dgemv("N", n_dirs, n_dirs, 1.0_rp, arh_object%coupling_matrix, &
                       n_dirs, projected_x, 1_ip, 0.0_rp, coupled_x, 1_ip)
            call dgemv("N", size(hess_x, kind=ip), n_dirs, 1.0_rp, &
                       arh_object%expansion_dirs, size(hess_x, kind=ip), coupled_x, &
                       1_ip, 1.0_rp, hess_x, 1_ip)
            deallocate(projected_x, coupled_x)
        end if

    end subroutine hess_x_arh_callback

    subroutine inv_hess_x_arh(x, inv_hess_x, error, level_shift)
        !
        ! this subroutine applies the exact inverse of the approximate Hessian to a
        ! vector, optionally level-shifted as (G - level_shift*I)^-1; G is the static
        ! part D plus the low-rank correction which is inverted using the
        ! Sherman-Morrison-Woodbury identity written as
        !
        !     (D + E C P^T)^-1 = D^-1 - D^-1 E (I + C P^T D^-1 E)^-1 C P^T D^-1
        !
        ! with E the expansion directions, P the projection directions and C the
        ! coupling matrix
        !
        use opentrustregion, only: verbosity_error
        use otr_common, only: level_shifted_divisors

        real(rp), intent(in), target :: x(:)
        real(rp), intent(out), target :: inv_hess_x(:)
        integer(ip), intent(out) :: error
        real(rp), intent(in), optional :: level_shift

        integer(ip) :: n_param, n_dirs, i, info
        real(rp) :: mu
        character(len=300) :: msg
        real(rp), allocatable :: rotated_x(:), divisors(:), scaled_x(:), &
                                 rotated_expansion(:, :), rotated_projection(:, :), &
                                 weighted_projection(:, :), dirs_overlap(:, :), &
                                 bracket_matrix(:, :), projected_x(:), bracket_rhs(:), &
                                 bracket_solution(:), correction(:)
        integer(ip), allocatable :: ipiv(:)
        external :: dgemv, dgemm, dgesv

        ! initialize error flag
        error = 0

        ! rebuild the approximate Hessian model if it is stale
        call rebuild_stale_hess_model(error)
        if (error /= 0) return

        ! an absent level shift gives the plain inverse of the approximate Hessian
        mu = 0.0_rp
        if (present(level_shift)) mu = level_shift

        ! refresh the eigendecomposition if the static Hessian part has changed
        call arh_object%orbitals%refresh_hess_eigen(arh_object%settings, error)
        if (error /= 0) return

        ! rotate x into the static-Hessian eigenbasis and apply the level-shifted
        ! diagonal D^-1
        rotated_x = arh_object%orbitals%rotate_to_hess_eigenbasis(x)
        divisors = level_shifted_divisors(arh_object%orbitals%get_hess_eigval_pairs(), &
                                          mu)
        scaled_x = rotated_x / divisors

        ! fall back to the static part alone if there is no low-rank correction
        if (.not. allocated(arh_object%coupling_matrix)) then
            inv_hess_x = arh_object%orbitals%rotate_from_hess_eigenbasis(scaled_x)
            return
        end if

        ! get parameters
        n_param = size(x)
        n_dirs = size(arh_object%expansion_dirs, 2)

        ! rotate both sets of directions into the same eigenbasis
        allocate(rotated_expansion(n_param, n_dirs), &
                 rotated_projection(n_param, n_dirs))
        do i = 1, n_dirs
            rotated_expansion(:, i) = arh_object%orbitals%rotate_to_hess_eigenbasis( &
                arh_object%expansion_dirs(:, i))
            rotated_projection(:, i) = arh_object%orbitals%rotate_to_hess_eigenbasis( &
                arh_object%projection_dirs(:, i))
        end do

        ! projected_x = P^T D^-1 x
        allocate(projected_x(n_dirs))
        call dgemv("T", n_param, n_dirs, 1.0_rp, rotated_projection, n_param, &
                   scaled_x, 1_ip, 0.0_rp, projected_x, 1_ip)

        ! dirs_overlap = P^T D^-1 E
        allocate(weighted_projection(n_param, n_dirs))
        do i = 1, n_param
            weighted_projection(i, :) = rotated_projection(i, :) / divisors(i)
        end do
        allocate(dirs_overlap(n_dirs, n_dirs))
        call dgemm("T", "N", n_dirs, n_dirs, n_param, 1.0_rp, weighted_projection, &
                   n_param, rotated_expansion, n_param, 0.0_rp, dirs_overlap, n_dirs)
        deallocate(weighted_projection)

        ! solve (I + C * P^T D^-1 E) bracket_solution = C * projected_x
        allocate(bracket_matrix(n_dirs, n_dirs))
        call dgemm("N", "N", n_dirs, n_dirs, n_dirs, 1.0_rp, &
                   arh_object%coupling_matrix, n_dirs, dirs_overlap, n_dirs, 0.0_rp, &
                   bracket_matrix, n_dirs)
        do i = 1, n_dirs
            bracket_matrix(i, i) = bracket_matrix(i, i) + 1.0_rp
        end do
        allocate(bracket_rhs(n_dirs))
        call dgemv("N", n_dirs, n_dirs, 1.0_rp, arh_object%coupling_matrix, n_dirs, &
                   projected_x, 1_ip, 0.0_rp, bracket_rhs, 1_ip)
        allocate(bracket_solution(n_dirs), ipiv(n_dirs))
        bracket_solution = bracket_rhs
        call dgesv(n_dirs, 1_ip, bracket_matrix, n_dirs, ipiv, bracket_solution, &
                   n_dirs, info)
        if (info /= 0) then
            write(msg, '(A, I0)') "Level-shifted approximate Hessian is singular: "// &
                "Error in DGESV, info = ", info
            call arh_object%settings%log(msg, verbosity_error, .true.)
            error = 1
            return
        end if

        ! correction = D^-1 E bracket_solution, result = D^-1 x - correction
        allocate(correction(n_param))
        call dgemv("N", n_param, n_dirs, 1.0_rp, rotated_expansion, n_param, &
                   bracket_solution, 1_ip, 0.0_rp, correction, 1_ip)
        correction = correction / divisors
        inv_hess_x = &
            arh_object%orbitals%rotate_from_hess_eigenbasis(scaled_x - correction)

    end subroutine inv_hess_x_arh

    subroutine precond_arh_callback(residual, mu, precond_residual, error)
        !
        ! this subroutine defines the preconditioner of the ARH approximate Hessian
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(in) :: mu
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        call inv_hess_x_arh(residual, precond_residual, error, mu)

    end subroutine precond_arh_callback

    subroutine precond_pd_arh_callback(residual, precond_residual, error)
        !
        ! this subroutine defines the positive-definite preconditioner based on the
        ! exact eigendecomposition of the static part of the Hessian
        !
        real(rp), intent(in), target :: residual(:)
        real(rp), intent(out), target :: precond_residual(:)
        integer(ip), intent(out) :: error

        call arh_object%orbitals%precond_pd(residual, precond_residual, &
                                            arh_object%settings, error)

    end subroutine precond_pd_arh_callback

    subroutine get_extra_trial_vectors_arh_callback(trial_vectors, error)
        !
        ! this subroutine returns the extra trial vectors of the orbital basis for the
        ! solver's initial trial space
        !
        real(rp), intent(out), target :: trial_vectors(:, :)
        integer(ip), intent(out) :: error

        call arh_object%orbitals%get_extra_trial_vectors(trial_vectors, &
                                                         arh_object%settings, error)

    end subroutine get_extra_trial_vectors_arh_callback

    subroutine init_arh_settings(self, error)
        !
        ! this subroutine initializes the ARH settings
        !
        use opentrustregion, only: verbosity_error

        class(arh_settings_type), intent(out) :: self
        integer(ip), intent(out) :: error

        ! initialize error flag
        error = 0

        select type (settings => self)
        type is (arh_settings_type)
            settings = default_arh_settings
        class default
            call settings%log("Augmented Roothaan-Hall settings could not be "// &
                              "initialized because initialization routine received "// &
                              "the wrong type. The type arh_settings_type was "// &
                              "likely subclassed without providing an "// &
                              "initialization routine.", verbosity_error, .true.)
            error = 1
        end select

    end subroutine init_arh_settings

    subroutine arh_deconstructor()
        !
        ! this subroutine deallocates the ARH objects
        !
        use otr_oao, only: oao_deconstructor
        use otr_mo, only: mo_deconstructor

        if (allocated(arh_object)) deallocate(arh_object)
        call mo_deconstructor()
        call oao_deconstructor()

    end subroutine arh_deconstructor

    function construct_arh_mo(mo) result(arh)
        !
        ! this function returns the ARH object for the MO basis pointing to the given
        ! MO object and to those of its quantities which are allocated
        !
        type(mo_type), intent(in), target :: mo
        type(arh_mo_type) :: arh

        ! point to the orbitals
        arh%orbitals => mo

        ! associate matrices and particle channels
        if (associated(mo%mo_coeff)) arh%mo_coeff => mo%mo_coeff
        if (allocated(mo%ao_overlap)) arh%ao_overlap => mo%ao_overlap
        if (allocated(mo%mo_channels)) arh%mo_channels => mo%mo_channels

    end function construct_arh_mo

    function construct_arh_oao(oao) result(arh)
        !
        ! this function returns the ARH object for the OAO basis pointing to the given
        ! OAO object and to those of its quantities which are allocated
        !
        type(oao_type), intent(in), target :: oao
        type(arh_oao_type) :: arh

        ! point to the orbitals
        arh%orbitals => oao

        ! associate matrices
        if (allocated(oao%s_inv_sqrt)) arh%s_inv_sqrt => oao%s_inv_sqrt
        if (allocated(oao%dm_oao)) arh%dm_oao => oao%dm_oao
        if (allocated(oao%fock_oo)) arh%fock_oo => oao%fock_oo
        if (allocated(oao%fock_vv)) arh%fock_vv => oao%fock_vv

    end function construct_arh_oao

    subroutine rotate_trial_arh_mo(self, kappa, rot_dm_ao, rot_dm_hist, error)
        !
        ! this subroutine returns the density matrix rotated by kappa in the AO basis
        ! and in the history basis (S D S) without moving the current orbitals
        !
        use otr_mo, only: rotate_mo_coeff, mo_transform

        class(arh_mo_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:)
        real(rp), intent(out) :: rot_dm_ao(:, :, :), rot_dm_hist(:, :, :)
        integer(ip), intent(out) :: error

        real(rp), allocatable :: rot_mo_coeff(:, :, :)
        integer(ip) :: i

        ! rotate MO coefficients and construct the density matrix
        allocate(rot_mo_coeff, mold=self%mo_coeff)
        call rotate_mo_coeff(kappa, self%mo_coeff, self%ao_overlap, self%mo_channels, &
                             rot_mo_coeff, rot_dm_ao, self%settings, error)
        if (error /= 0) return

        ! multiply the density matrix by the AO overlap matrix from both sides
        do i = 1, self%orbitals%n_particle
            rot_dm_hist(:, :, i) = mo_transform(self%ao_overlap, rot_dm_ao(:, :, i))
        end do

    end subroutine rotate_trial_arh_mo

    subroutine rotate_trial_arh_oao(self, kappa, rot_dm_ao, rot_dm_hist, error)
        !
        ! this subroutine returns the density matrix rotated by kappa in the AO basis
        ! and in the history basis (the OAO basis) without moving the current orbitals
        !
        use otr_oao, only: rotate_dm_ao

        class(arh_oao_type), intent(in) :: self
        real(rp), intent(in) :: kappa(:)
        real(rp), intent(out) :: rot_dm_ao(:, :, :), rot_dm_hist(:, :, :)
        integer(ip), intent(out) :: error

        ! rotate density matrix in the AO and OAO basis
        call rotate_dm_ao(kappa, self%dm_oao, self%s_inv_sqrt, rot_dm_ao, &
                          self%settings, error, rot_dm_hist)

    end subroutine rotate_trial_arh_oao

    function history_dm_arh_mo(self) result(dm)
        !
        ! this function returns the current density matrix in the history basis (the
        ! density matrix in the AO basis multiplied by the AO overlap matrix from both
        ! sides, S D S) which unlike D transforms to the MO basis like the potentials,
        ! as C^T M C
        !
        use otr_mo, only: mo_transform

        class(arh_mo_type), intent(in) :: self
        real(rp), allocatable :: dm(:, :, :)

        integer(ip) :: i

        ! multiply the density matrix by the AO overlap matrix from both sides
        allocate(dm, mold=self%orbitals%dm_ao)
        do i = 1, self%orbitals%n_particle
            dm(:, :, i) = mo_transform(self%ao_overlap, self%orbitals%dm_ao(:, :, i))
        end do

    end function history_dm_arh_mo

    function history_dm_arh_oao(self) result(dm)
        !
        ! this function returns the current density matrix in the history basis (the
        ! OAO basis)
        !
        class(arh_oao_type), intent(in) :: self
        real(rp), allocatable :: dm(:, :, :)

        dm = self%dm_oao

    end function history_dm_arh_oao

    function to_history_basis_arh_mo(self, matrix_ao) result(matrix)
        !
        ! this function transforms a potential given in the AO basis to the history
        ! basis, which for the MO basis, since it changes with every orbital update, is
        ! the AO basis itself
        !
        class(arh_mo_type), intent(in) :: self
        real(rp), intent(in) :: matrix_ao(:, :, :)
        real(rp), allocatable :: matrix(:, :, :)

        matrix = matrix_ao

    end function to_history_basis_arh_mo

    function to_history_basis_arh_oao(self, matrix_ao) result(matrix)
        !
        ! this function transforms a potential given in the AO basis to the history
        ! basis (the OAO basis)
        !
        use otr_oao, only: symmetric_transformation

        class(arh_oao_type), intent(in) :: self
        real(rp), intent(in) :: matrix_ao(:, :, :)
        real(rp), allocatable :: matrix(:, :, :)

        matrix = symmetric_transformation(self%s_inv_sqrt, matrix_ao)

    end function to_history_basis_arh_oao

    subroutine history_columns_arh_mo(self, diff, cols, packed)
        !
        ! this subroutine expresses a history of difference matrices in the AO basis,
        ! with density matrices stored as S D S, as history columns (the full
        ! difference matrices in the current MO basis whose inner products define the
        ! approximate Hessian) and as packed columns (the parameters of their
        ! occupied-virtual blocks)
        !
        use otr_mo, only: mo_transform, pack_ov

        class(arh_mo_type), intent(in) :: self
        real(rp), intent(in) :: diff(:, :, :, :)
        real(rp), intent(out), allocatable :: cols(:, :), packed(:, :)

        integer(ip) :: n_list, n_mo, n_occ, k, i
        integer(ip), allocatable :: history_rows(:, :), packed_rows(:, :)
        real(rp), allocatable :: diff_mo(:, :)

        ! number of history entries and MOs
        n_list = size(diff, 4, kind=ip)
        n_mo = size(self%mo_coeff, 2, kind=ip)

        ! rows of every particle channel in the history and packed columns
        history_rows = self%history_channel_rows()
        packed_rows = self%packed_channel_rows()

        ! transform every difference matrix to the current MO basis, keep it in full as
        ! history column and retain the parameters of its occupied-virtual block as
        ! packed column
        allocate(cols(history_rows(2, self%orbitals%n_particle), n_list), &
                 packed(self%orbitals%n_param, n_list))
        do i = 1, self%orbitals%n_particle
            n_occ = self%mo_channels(i)%n_occ
            do k = 1, n_list
                diff_mo = mo_transform(self%mo_coeff(:, :, i), diff(:, :, i, k))
                cols(history_rows(1, i):history_rows(2, i), k) = reshape(diff_mo, &
                                                                         [n_mo * n_mo])
                packed(packed_rows(1, i):packed_rows(2, i), k) = &
                    pack_ov(diff_mo(:n_occ, n_occ + 1:), self%mo_channels(i))
            end do
        end do

    end subroutine history_columns_arh_mo

    subroutine history_columns_arh_oao(self, diff, cols, packed)
        !
        ! this subroutine expresses a history of difference matrices in the OAO basis
        ! as history columns (the flattened difference matrices whose inner products
        ! define the approximate Hessian) and as packed columns (their projections onto
        ! the non-redundant parameter space)
        !
        use otr_oao, only: project_asymm, pack_asymm

        class(arh_oao_type), intent(in) :: self
        real(rp), intent(in) :: diff(:, :, :, :)
        real(rp), intent(out), allocatable :: cols(:, :), packed(:, :)

        integer(ip) :: n_list, n_param, k
        integer(ip), allocatable :: packed_rows(:, :)

        ! number of history entries and parameters
        n_list = size(diff, 4, kind=ip)
        packed_rows = self%packed_channel_rows()
        n_param = packed_rows(2, size(packed_rows, 2))

        ! flatten every difference matrix into a history column
        cols = reshape(diff, [size(diff, 1, kind=ip) * size(diff, 2, kind=ip) * &
                              size(diff, 3, kind=ip), n_list])

        ! project every difference matrix onto the occupied-virtual and
        ! virtual-occupied subspace and pack it
        allocate(packed(n_param, n_list))
        do k = 1, n_list
            packed(:, k) = pack_asymm(project_asymm(diff(:, :, :, k), self%dm_oao), &
                                      n_param)
        end do

    end subroutine history_columns_arh_oao

    function history_channel_rows_arh_mo(self) result(rows)
        !
        ! this function returns the rows every particle channel occupies in a history
        ! column, which holds the full difference matrix of every channel in the MO
        ! basis
        !
        class(arh_mo_type), intent(in) :: self
        integer(ip), allocatable :: rows(:, :)

        rows = channel_rows(spread(size(self%mo_coeff, 2, kind=ip)**2, 1, &
                                   self%orbitals%n_particle))

    end function history_channel_rows_arh_mo

    function history_channel_rows_arh_oao(self) result(rows)
        !
        ! this function returns the rows every particle channel occupies in a history
        ! column, which holds the full difference matrix of every channel in the OAO
        ! basis
        !
        class(arh_oao_type), intent(in) :: self
        integer(ip), allocatable :: rows(:, :)

        integer(ip) :: n_ao, n_particle

        ! number of AOs and number of particles
        n_ao = self%orbitals%n_ao
        n_particle = self%orbitals%n_particle

        ! every particle channel holds the same number of rows
        rows = channel_rows(spread(n_ao**2, 1, n_particle))

    end function history_channel_rows_arh_oao

    function packed_channel_rows_arh_mo(self) result(rows)
        !
        ! this function returns the rows every particle channel occupies in a packed
        ! column, which holds the occupied-virtual rotations of every channel in the MO
        ! basis between orbitals of the same irrep
        !
        use otr_mo, only: mo_param_rows

        class(arh_mo_type), intent(in) :: self
        integer(ip), allocatable :: rows(:, :)

        rows = mo_param_rows(self%mo_channels)

    end function packed_channel_rows_arh_mo

    function packed_channel_rows_arh_oao(self) result(rows)
        !
        ! this function returns the rows every particle channel occupies in a packed
        ! column, which holds the antisymmetric parameters of every channel in the OAO
        ! basis
        !
        class(arh_oao_type), intent(in) :: self
        integer(ip), allocatable :: rows(:, :)

        integer(ip) :: n_ao, n_particle

        ! number of AOs and number of particles
        n_ao = self%orbitals%n_ao
        n_particle = self%orbitals%n_particle

        ! every particle channel holds the same number of rows
        rows = channel_rows(spread(n_ao * (n_ao - 1) / 2, 1, n_particle))

    end function packed_channel_rows_arh_oao

    function hess_x_static_arh_mo(self, x) result(hess_x)
        !
        ! this function applies the static part of the Hessian in the MO basis (which
        ! is built from the occupied-occupied and virtual-virtual blocks of the Fock
        ! matrix) to a trial vector
        !
        use otr_mo, only: hess_x_static_mo

        class(arh_mo_type), intent(in) :: self
        real(rp), intent(in) :: x(:)
        real(rp), allocatable :: hess_x(:)

        hess_x = hess_x_static_mo(x, self%mo_channels)

    end function hess_x_static_arh_mo

    function hess_x_static_arh_oao(self, x) result(hess_x)
        !
        ! this function applies the static part of the Hessian in the OAO basis (which
        ! is built from the occupied-occupied and virtual-virtual blocks of the Fock
        ! matrix) to a trial vector
        !
        use otr_oao, only: unpack_asymm, hess_x_static_oao, project_asymm, pack_asymm

        class(arh_oao_type), intent(in) :: self
        real(rp), intent(in) :: x(:)
        real(rp), allocatable :: hess_x(:)

        integer(ip) :: n_ao, n_particle
        real(rp), allocatable :: x_full(:, :, :), hess_x_full(:, :, :)

        ! number of AOs and number of particles
        n_ao = self%orbitals%n_ao
        n_particle = self%orbitals%n_particle

        ! unpack trial vector
        x_full = unpack_asymm(x, n_particle, n_ao)

        ! get static part
        hess_x_full = hess_x_static_oao(x_full, self%fock_oo, self%fock_vv)

        ! project, scale and pack the static part; the static part is already confined
        ! to the occupied-virtual and virtual-occupied subspace in exact arithmetic,
        ! but projecting prevents numerical leakage into the redundant subspace
        hess_x_full = project_asymm(hess_x_full, self%dm_oao)
        if (n_particle == 1) then
            hess_x = 4.0_rp * pack_asymm(hess_x_full, size(x, kind=ip))
        else
            hess_x = 2.0_rp * pack_asymm(hess_x_full, size(x, kind=ip))
        end if
        deallocate(x_full, hess_x_full)

    end function hess_x_static_arh_oao

    subroutine cache_channel_split_dirs(packed, rows, map, chol, dirs)
        !
        ! this subroutine caches the open-shell history projections one spin channel at
        ! a time: column i holds channel 1 of packed column i, column n_list+i holds
        ! channel 2; the two are never summed, since the per-channel metric contraction
        ! does not mix channels; the concatenated columns are then rebased into the
        ! orthonormalized basis defined by map and chol
        !
        real(rp), intent(in) :: packed(:, :), chol(:, :)
        integer(ip), intent(in) :: rows(:, :), map(:)
        real(rp), intent(out), allocatable :: dirs(:, :)

        integer(ip) :: n_list
        real(rp), allocatable :: raw(:, :)

        ! place every channel in its own set of columns, zeroing the other channel
        n_list = size(packed, 2, kind=ip)
        allocate(raw(size(packed, 1), 2 * n_list))
        raw = 0.0_rp
        raw(rows(1, 1):rows(2, 1), :n_list) = packed(rows(1, 1):rows(2, 1), :)
        raw(rows(1, 2):rows(2, 2), n_list + 1:) = packed(rows(1, 2):rows(2, 2), :)

        ! rebase into the orthonormalized basis
        dirs = rebase_dirs(raw, map, chol)

    end subroutine cache_channel_split_dirs

    subroutine cache_combined_channel_dirs(same, opp, rows, map, chol, dirs)
        !
        ! this subroutine caches the open-shell history projections two channels at a
        ! time: column i combines channel 1 of the same-spin potential with channel 2
        ! of the opposite-spin potential, column n_list+i the mirror image; the
        ! combined columns are then rebased into the orthonormalized basis defined by
        ! map and chol
        !
        real(rp), intent(in) :: same(:, :), opp(:, :), chol(:, :)
        integer(ip), intent(in) :: rows(:, :), map(:)
        real(rp), intent(out), allocatable :: dirs(:, :)

        integer(ip) :: n_list
        real(rp), allocatable :: raw(:, :)

        ! combine the channels of the same-spin and opposite-spin potentials
        n_list = size(same, 2, kind=ip)
        allocate(raw(size(same, 1), 2 * n_list))
        raw = 0.0_rp
        raw(rows(1, 1):rows(2, 1), :n_list) = same(rows(1, 1):rows(2, 1), :)
        raw(rows(1, 2):rows(2, 2), :n_list) = opp(rows(1, 2):rows(2, 2), :)
        raw(rows(1, 1):rows(2, 1), n_list + 1:) = opp(rows(1, 1):rows(2, 1), :)
        raw(rows(1, 2):rows(2, 2), n_list + 1:) = same(rows(1, 2):rows(2, 2), :)

        ! rebase into the orthonormalized basis
        dirs = rebase_dirs(raw, map, chol)

    end subroutine cache_combined_channel_dirs

    subroutine get_low_rank_hess_factors()
        !
        ! this subroutine assembles the low-rank part of the approximate Hessian in the
        ! packed parameter space, as
        !
        !     G_low_rank = expansion_dirs * coupling_matrix * transpose(projection_dirs)
        !
        ! by constructing expansion_dirs, coupling_matrix, and projection_dirs
        !
        ! the response is defined through Frobenius inner products
        ! <history_k, delta_dm(x)> of a history matrix with the density response to the
        ! parameter vector x, yet the whole correction can be expressed on packed
        ! parameter vectors alone, since in both orbital bases this contraction is
        ! twice the plain dot product <dirs(:, k), x> of the packed history column with
        ! x:
        ! - in the OAO basis, the density response is project_symm of the unpacked x
        !   and the packed column packs project_asymm of the history matrix, and the
        !   two projections (on antisymmetric and on symmetric input, respectively) are
        !   adjoint with respect to the Frobenius inner product
        ! - in the MO basis, the density response has x as both of its off-diagonal
        !   blocks, which vanish between orbitals of different irreps, and the packed
        !   column holds the occupied-virtual block of the symmetric history matrix
        !   between orbitals of the same irrep, so both off-diagonal blocks contribute
        !   the same dot product
        ! packing and projecting are linear, so the output expansion of the response
        ! collapses the same way
        !
        ! every closed-shell coupling matrix therefore carries an overall factor of 8:
        ! the factor of 2 above, times a factor of 4 (or 2) for the closed- (or
        ! open-)shell scaling of the Hessian linear transformation
        !
        integer(ip) :: n_linear, n_nonlinear, n_total, i
        real(rp) :: shell_scale

        ! discard the factors assembled for the previous history
        if (allocated(arh_object%expansion_dirs)) deallocate(arh_object%expansion_dirs)
        if (allocated(arh_object%projection_dirs)) &
            deallocate(arh_object%projection_dirs)
        if (allocated(arh_object%coupling_matrix)) &
            deallocate(arh_object%coupling_matrix)

        ! set scaling factor for closed- and open-shell systems
        shell_scale = merge(1.0_rp, 0.5_rp, arh_object%orbitals%n_particle == 1)

        select case (arh_object%settings%arh_type)
        ! multisecant SR1: the linear and non-linear parts of the potential enter as
        ! independent blocks, each contracted against its own inverse rather than
        ! against a single inverse of the summed potential
        case ("ms_sr1")
            if (.not. allocated(arh_object%linear_potential_dirs) .or. &
                .not. allocated(arh_object%nonlinear_potential_dirs)) return
            n_linear = size(arh_object%linear_potential_dirs, 2)
            n_nonlinear = size(arh_object%nonlinear_potential_dirs, 2)
            n_total = n_linear + n_nonlinear
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%orbitals%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, n_linear + 1:) = &
                arh_object%nonlinear_potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_linear, :n_linear) = 8.0_rp * shell_scale * &
                                                               arh_object%a_inv
            arh_object%coupling_matrix(n_linear + 1:, n_linear + 1:) = &
                8.0_rp * shell_scale * arh_object%a_inv_comb

        ! subspace-projected multisecant: the response both expands in and contracts
        ! against the density difference history alone
        case ("ms_sp")
            if (.not. allocated(arh_object%dm_dirs) .or. &
                .not. allocated(arh_object%dm_dirs_nonlinear)) return
            n_linear = size(arh_object%dm_dirs, 2)
            n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
            n_total = n_linear + n_nonlinear
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%orbitals%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_linear + 1:) = arh_object%dm_dirs_nonlinear
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_linear, :n_linear) = 8.0_rp * shell_scale * &
                                                               arh_object%a_sym
            arh_object%coupling_matrix(n_linear + 1:, n_linear + 1:) = &
                8.0_rp * shell_scale * arh_object%a_sym_nonlinear

        ! symmetrized ARH: the density and potential difference histories couple in
        ! both directions, so the coupling matrix is purely off-diagonal; also
        ! introduces an additional factor of 1/2
        case ("symm_arh")
            if (.not. allocated(arh_object%dm_dirs) .or. &
                .not. allocated(arh_object%linear_potential_dirs)) return
            n_linear = size(arh_object%dm_dirs, 2)
            n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
            n_total = 2 * (n_linear + n_nonlinear)
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%orbitals%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_linear + 1:2 * n_linear) = &
                arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, 2 * n_linear + 1:2 * n_linear + &
                                      n_nonlinear) = arh_object%dm_dirs_nonlinear
            arh_object%expansion_dirs(:, 2 * n_linear + n_nonlinear + 1:) = &
                arh_object%nonlinear_potential_dirs
            ! the metric scaling collapses in the orthonormalized S-basis
            arh_object%coupling_matrix = 0.0_rp
            do i = 1, n_linear
                arh_object%coupling_matrix(i, n_linear + i) = 4.0_rp * shell_scale
                arh_object%coupling_matrix(n_linear + i, i) = 4.0_rp * shell_scale
            end do
            do i = 1, n_nonlinear
                arh_object%coupling_matrix(2 * n_linear + i, 2 * n_linear + &
                                           n_nonlinear + i) = 4.0_rp * shell_scale
                arh_object%coupling_matrix(2 * n_linear + n_nonlinear + i, &
                                           2 * n_linear + i) = 4.0_rp * shell_scale
            end do

        ! multisecant PSB: the symmetrized ARH coupling with an additional
        ! density-density block subtracting the doubly counted curvature
        case ("ms_psb")
            if (.not. allocated(arh_object%dm_dirs) .or. &
                .not. allocated(arh_object%linear_potential_dirs)) return
            n_linear = size(arh_object%dm_dirs, 2)
            n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
            n_total = 2 * (n_linear + n_nonlinear)
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%orbitals%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%dm_dirs
            arh_object%expansion_dirs(:, n_linear + 1:2 * n_linear) = &
                arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, 2 * n_linear + 1:2 * n_linear + &
                                      n_nonlinear) = arh_object%dm_dirs_nonlinear
            arh_object%expansion_dirs(:, 2 * n_linear + n_nonlinear + 1:) = &
                arh_object%nonlinear_potential_dirs
            arh_object%coupling_matrix = 0.0_rp
            arh_object%coupling_matrix(:n_linear, :n_linear) = -8.0_rp * shell_scale * &
                                                               arh_object%a_sym
            arh_object%coupling_matrix(2 * n_linear + 1:2 * n_linear + n_nonlinear, &
                                       2 * n_linear + 1:2 * n_linear + n_nonlinear) = &
                -8.0_rp * shell_scale * arh_object%a_sym_nonlinear
            ! the metric scaling in the off-diagonal blocks collapses in the
            ! orthonormalized S-basis
            do i = 1, n_linear
                arh_object%coupling_matrix(i, n_linear + i) = 8.0_rp * shell_scale
                arh_object%coupling_matrix(n_linear + i, i) = 8.0_rp * shell_scale
            end do
            do i = 1, n_nonlinear
                arh_object%coupling_matrix(2 * n_linear + i, 2 * n_linear + &
                                           n_nonlinear + i) = 8.0_rp * shell_scale
                arh_object%coupling_matrix(2 * n_linear + n_nonlinear + i, &
                                           2 * n_linear + i) = 8.0_rp * shell_scale
            end do

        ! standard ARH: the response expands in the potential difference history while
        ! contracting against the density difference history, so unlike every other
        ! type the two sets of directions differ
        case ("arh")
            if (.not. allocated(arh_object%dm_dirs) .or. &
                .not. allocated(arh_object%linear_potential_dirs)) return
            n_linear = size(arh_object%dm_dirs, 2)
            n_nonlinear = size(arh_object%dm_dirs_nonlinear, 2)
            n_total = n_linear + n_nonlinear
            if (n_total == 0) return

            allocate(arh_object%expansion_dirs(arh_object%orbitals%n_param, n_total), &
                     arh_object%projection_dirs(arh_object%orbitals%n_param, n_total), &
                     arh_object%coupling_matrix(n_total, n_total))
            arh_object%expansion_dirs(:, :n_linear) = arh_object%linear_potential_dirs
            arh_object%expansion_dirs(:, n_linear + 1:) = &
                arh_object%nonlinear_potential_dirs
            arh_object%projection_dirs(:, :n_linear) = arh_object%dm_dirs
            arh_object%projection_dirs(:, n_linear + 1:) = arh_object%dm_dirs_nonlinear
            arh_object%coupling_matrix = 0.0_rp
            ! the metric scaling collapses in the orthonormalized S-basis
            do i = 1, n_total
                arh_object%coupling_matrix(i, i) = 8.0_rp * shell_scale
            end do

        case default
            return
        end select

        ! every type except standard ARH expands in and contracts against the same set
        ! of directions
        if (allocated(arh_object%expansion_dirs) .and. &
            .not. allocated(arh_object%projection_dirs)) &
            arh_object%projection_dirs = arh_object%expansion_dirs

    end subroutine get_low_rank_hess_factors

    subroutine build_a_part(dm_cols, v_cols, a)
        !
        ! this subroutine constructs A = S^T Y for a single part of a coupled
        ! potential-difference response from its history columns and symmetrizes it
        !
        real(rp), intent(in) :: dm_cols(:, :), v_cols(:, :)
        real(rp), intent(out) :: a(:, :)

        integer(ip) :: n_dm, flat_len
        external :: dgemm

        ! nothing to build for an empty history
        n_dm = size(dm_cols, 2, kind=ip)
        if (n_dm == 0) return

        ! A = S^T Y on the full history
        flat_len = size(dm_cols, 1, kind=ip)
        call dgemm("T", "N", n_dm, n_dm, flat_len, 1.0_rp, dm_cols, flat_len, v_cols, &
                   flat_len, 0.0_rp, a, n_dm)

        ! symmetrize A
        a = 0.5_rp * (a + transpose(a))

    end subroutine build_a_part

    function build_a_transformed(dm_cols, v_cols, map, chol) result(a_t)
        !
        ! this function builds A = S^T Y for a single (linear or non-linear) part of
        ! the coupled potential-difference response and congruence-transforms it into
        ! the orthonormalized S-basis belonging to that part
        !
        real(rp), intent(in) :: dm_cols(:, :), v_cols(:, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: a_t(:, :)

        integer(ip) :: n_diff
        real(rp), allocatable :: a(:, :)

        n_diff = size(dm_cols, 2, kind=ip)
        allocate(a(n_diff, n_diff))
        call build_a_part(dm_cols, v_cols, a)
        a_t = congruence_transform(a, map, chol)
        deallocate(a)

    end function build_a_transformed

    function build_a_block_linear_os(dm_cols, v_same_cols, v_opp_cols, rows, map, &
                                     chol) result(a_block)
        !
        ! this function builds the dense, cross-channel-symmetrized open-shell A = S^T
        ! Y matrix of the linear response entering the MS-SP and MS-PSB contributions;
        ! the same-spin response fills the per-channel diagonal blocks and the
        ! opposite-spin response the off-diagonal ones, and the result is
        ! congruence-transformed to the orthonormalized, rank-independent S-basis
        !
        real(rp), intent(in) :: dm_cols(:, :), v_same_cols(:, :), v_opp_cols(:, :), &
                                chol(:, :)
        integer(ip), intent(in) :: rows(:, :), map(:)
        real(rp), allocatable :: a_block(:, :)

        integer(ip) :: n_diff, j, n_rows
        real(rp), allocatable :: a_same(:, :, :), a_opp(:, :, :), a_full(:, :), &
                                 dm_rows(:, :), v_same_rows(:, :), v_opp_rows(:, :)
        external :: dgemm

        n_diff = size(dm_cols, 2, kind=ip)
        allocate(a_same(n_diff, n_diff, 2), a_opp(n_diff, n_diff, 2))

        do j = 1, 2
            ! copy the rows of the spin channel
            n_rows = rows(2, j) - rows(1, j) + 1
            dm_rows = dm_cols(rows(1, j):rows(2, j), :)
            v_same_rows = v_same_cols(rows(1, j):rows(2, j), :)
            v_opp_rows = v_opp_cols(rows(1, j):rows(2, j), :)

            ! same-spin and opposite-spin contributions, nothing for an empty history
            call build_a_part(dm_rows, v_same_rows, a_same(:, :, j))
            if (n_diff > 0) &
                call dgemm("T", "N", n_diff, n_diff, n_rows, 1.0_rp, dm_rows, n_rows, &
                           v_opp_rows, n_rows, 0.0_rp, a_opp(:, :, j), n_diff)
        end do

        ! cross-symmetrize off-diagonal blocks exactly, since the opposite-spin
        ! potential difference is purely linear
        call cross_symmetrize(a_opp(:, :, 1), a_opp(:, :, 2))

        allocate(a_full(2 * n_diff, 2 * n_diff))
        a_full(1:n_diff, 1:n_diff) = a_same(:, :, 1)
        a_full(1:n_diff, n_diff + 1:2 * n_diff) = a_opp(:, :, 1)
        a_full(n_diff + 1:2 * n_diff, 1:n_diff) = a_opp(:, :, 2)
        a_full(n_diff + 1:2 * n_diff, n_diff + 1:2 * n_diff) = a_same(:, :, 2)
        deallocate(a_same, a_opp)
        a_block = congruence_transform(a_full, map, chol)
        deallocate(a_full)

    end function build_a_block_linear_os

    function build_a_block_nonlinear_os(dm_cols, v_nonlinear_cols, rows, map, chol) &
        result(a_block)
        !
        ! this function builds the open-shell A = S^T Y matrix of the non-linear
        ! response entering the MS-SP and MS-PSB contributions; the non-linear response
        ! has no opposite-spin counterpart to cross-symmetrize against, so it fills the
        ! per-channel diagonal blocks only, and the result is congruence-transformed to
        ! the orthonormalized, rank-independent S-basis
        !
        real(rp), intent(in) :: dm_cols(:, :), v_nonlinear_cols(:, :), chol(:, :)
        integer(ip), intent(in) :: rows(:, :), map(:)
        real(rp), allocatable :: a_block(:, :)

        integer(ip) :: n_diff, j
        real(rp), allocatable :: a_same(:, :, :), a_full(:, :), dm_rows(:, :), &
                                 v_nonlinear_rows(:, :)

        n_diff = size(dm_cols, 2, kind=ip)
        allocate(a_same(n_diff, n_diff, 2))
        do j = 1, 2
            ! copy the rows of the spin channel
            dm_rows = dm_cols(rows(1, j):rows(2, j), :)
            v_nonlinear_rows = v_nonlinear_cols(rows(1, j):rows(2, j), :)
            call build_a_part(dm_rows, v_nonlinear_rows, a_same(:, :, j))
        end do

        allocate(a_full(2 * n_diff, 2 * n_diff))
        a_full = 0.0_rp
        a_full(1:n_diff, 1:n_diff) = a_same(:, :, 1)
        a_full(n_diff + 1:2 * n_diff, n_diff + 1:2 * n_diff) = a_same(:, :, 2)
        deallocate(a_same)
        a_block = congruence_transform(a_full, map, chol)
        deallocate(a_full)

    end function build_a_block_nonlinear_os

    subroutine cross_symmetrize(a12, a21)
        !
        ! this subroutine cross-symmetrizes two related off-diagonal blocks a12(i, k)
        ! and a21(k, i) of a larger matrix by averaging each with the transpose of the
        ! other, which is appropriate here since a12 and transpose(a21) are equal in
        ! exact arithmetic and any observed mismatch is therefore numerical noise
        !
        real(rp), intent(inout) :: a12(:, :), a21(:, :)

        real(rp), allocatable :: averaged(:, :)

        averaged = 0.5_rp * (a12 + transpose(a21))
        a12 = averaged
        a21 = transpose(averaged)
        deallocate(averaged)

    end subroutine cross_symmetrize

    function truncated_eigval_inv(eig_vals, eig_val_thresh) result(eig_vals_inv)
        !
        ! this function returns a hard-truncated pseudoinverse: the exact inverse is
        ! kept for eigenvalues above the given threshold, while eigenvalues at or below
        ! it are discarded entirely rather than smoothly damped
        !
        real(rp), intent(in) :: eig_vals(:), eig_val_thresh
        real(rp) :: eig_vals_inv(size(eig_vals))

        where (abs(eig_vals) > eig_val_thresh)
            eig_vals_inv = 1.0_rp / eig_vals
        elsewhere
            eig_vals_inv = 0.0_rp
        end where

    end function truncated_eigval_inv

    function history_step_mask(dm_cols) result(keep)
        !
        ! this function returns, for every history column, whether the non-linear
        ! multisecant system may use it, depending on step length
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: dm_cols(:, :)
        logical, allocatable :: keep(:)

        integer(ip) :: n_list, flat_len, i
        real(rp) :: shortest, cutoff
        real(rp), allocatable :: steps(:)
        real(rp), external :: dnrm2

        n_list = size(dm_cols, 2, kind=ip)
        allocate(keep(n_list), steps(n_list))
        keep = .true.
        if (n_list < 2) return

        ! get step lengths
        flat_len = size(dm_cols, 1, kind=ip)
        do i = 1, n_list
            steps(i) = dnrm2(flat_len, dm_cols(:, i), 1_ip)
        end do

        ! get shortest numerically relevant step
        shortest = huge(1.0_rp)
        do i = 1, n_list
            if (steps(i) > numerical_zero * maxval(steps)) &
                shortest = min(shortest, steps(i))
        end do
        if (shortest >= huge(1.0_rp)) return

        ! define cut off and decide on steps to keep
        cutoff = history_step_ratio * shortest
        keep = steps <= cutoff
        deallocate(steps)

    end function history_step_mask

    function median(x) result(med)
        !
        ! this function returns the median of an array
        !
        real(rp), intent(in) :: x(:)
        real(rp) :: med

        integer(ip) :: i, j, n
        real(rp) :: val
        real(rp), allocatable :: sorted(:)

        n = size(x)
        if (n == 0) then
            med = 0.0_rp
            return
        end if

        ! sort a copy by insertion, which is short and fast enough for the few hundred
        ! history directions this is called with at most
        sorted = x
        do i = 2, n
            val = sorted(i)
            j = i - 1
            do while (j >= 1)
                if (sorted(j) <= val) exit
                sorted(j + 1) = sorted(j)
                j = j - 1
            end do
            sorted(j + 1) = val
        end do

        ! average the two central entries if there is no single one
        if (mod(n, 2_ip) == 1) then
            med = sorted((n + 1) / 2)
        else
            med = 0.5_rp * (sorted(n / 2) + sorted(n / 2 + 1))
        end if
        deallocate(sorted)

    end function median

    function resolvable_residual(steps, responses, keep) result(min_residual)
        !
        ! this function returns the shortest residual norm a history direction has to
        ! contribute for its response to describe curvature rather than the error in
        ! it: the response is an exactly symmetric function of the step only for a
        ! vanishing step, so the antisymmetric part of the transformed response matrix
        ! Q^T Z collects both the numerical noise of the response and the finite-step
        ! error of a curvature that changes along the step, neither of which a
        ! symmetric fit can use; every column of Z is divided by the residual norm of
        ! its own direction, so multiplying the typical size of column k of that part
        ! by that residual norm undoes the division, and the median over the directions
        ! makes a single error scale of it; a direction shorter than that scale divided
        ! by the typical curvature contributes more amplified error than curvature
        !
        real(rp), intent(in) :: steps(:, :), responses(:, :)
        logical, intent(in) :: keep(:)
        real(rp) :: min_residual

        integer(ip) :: n_accepted, i, k, m
        integer(ip), allocatable :: map(:)
        real(rp), allocatable :: chol(:, :), q(:, :), z(:, :), a(:, :), &
                                 pair_asymmetries(:), asymmetries(:), residuals(:), &
                                 curvatures(:)
        real(rp) :: response_error, curvature
        external :: dgemm

        ! a direction needs at least two partners for a typical asymmetry of its own
        ! column to mean anything
        min_residual = 0.0_rp
        call factorize_history(steps, chol, map, n_accepted, keep)
        if (n_accepted < 3) return

        ! transformed response matrix Q^T Z in the orthonormalized step basis
        q = rebase_dirs(steps, map, chol)
        z = rebase_dirs(responses, map, chol)
        allocate(a(n_accepted, n_accepted))
        call dgemm("T", "N", n_accepted, n_accepted, size(q, 1, kind=ip), 1.0_rp, q, &
                   size(q, 1, kind=ip), z, size(z, 1, kind=ip), 0.0_rp, a, n_accepted)

        ! per-direction asymmetry, residual norm and curvature
        allocate(pair_asymmetries(n_accepted - 1), asymmetries(n_accepted), &
                 residuals(n_accepted), curvatures(n_accepted))
        do k = 1, n_accepted
            m = 0
            do i = 1, n_accepted
                if (i == k) cycle
                m = m + 1
                pair_asymmetries(m) = abs(a(i, k) - a(k, i))
            end do
            asymmetries(k) = median(pair_asymmetries)
            residuals(k) = abs(chol(k, k))
            curvatures(k) = abs(a(k, k))
        end do

        ! undo the amplification of every column to get the error of the response
        ! itself, and take the curvature the matrix describes at its typical size
        response_error = median(asymmetries * residuals)
        curvature = median(curvatures)

        ! residual norm at which the amplified error reaches the curvature scale
        if (curvature > 0.0_rp) min_residual = response_error / curvature
        deallocate(chol, map, q, z, a, pair_asymmetries, asymmetries, residuals, &
                   curvatures)

    end function resolvable_residual

    subroutine factorize_history(vecs, chol, map, n_accepted, keep, min_residual)
        !
        ! this subroutine performs a pivoted, rank-revealing Cholesky factorization of
        ! the Gram matrix of a set of history vectors, taken on the Gram normalized to
        ! a unit diagonal so that the pivot order and the rank tolerance test genuine
        ! linear dependence rather than step length: the secant conditions are
        ! homogeneous, so rescaling a column with its response leaves the approximate
        ! Hessian unchanged, while a length-driven pivot would discard the newest,
        ! shortest steps as dependent purely for being small; the normalization is
        ! undone in the returned factor, so together with the map back to original
        ! indices any other history can be re-expressed in the same orthonormalized
        ! basis by a single triangular solve, without forming a metric inverse
        !
        ! keep optionally restricts which columns a caller trusts; omitting it uses the
        ! whole history, which is right for a response that is exact at any distance
        !
        ! min_residual optionally demands a shortest residual norm a column has to
        ! contribute to be accepted, which is right for a response whose own error
        ! limits how short a direction can be and still be resolved; the Gram is then
        ! left unnormalized so that its pivots are comparable to that absolute length,
        ! which also makes the pivoting rank the columns by the length they contribute
        ! rather than by independence, since that length is what the tolerance bounds,
        ! and leaves the returned factor already carrying the norms; a non-positive
        ! value demands nothing
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: vecs(:, :)
        real(rp), intent(out), allocatable :: chol(:, :)
        integer(ip), intent(out), allocatable :: map(:)
        integer(ip), intent(out) :: n_accepted
        logical, intent(in), optional :: keep(:)
        real(rp), intent(in), optional :: min_residual

        integer(ip) :: vec_len, n_vecs, i, j, info
        logical :: normalize
        real(rp) :: tol
        real(rp), allocatable :: metric(:, :), work(:), col_norm(:)
        logical, allocatable :: usable(:)
        integer(ip), allocatable :: piv(:)
        external :: dsyrk, dpstrf

        vec_len = size(vecs, 1)
        n_vecs = size(vecs, 2)

        n_accepted = 0
        if (n_vecs == 0) then
            allocate(chol(0, 0), map(0))
            return
        end if

        ! decide whether the pivoting ranks the columns by independence or by the
        ! length they contribute
        normalize = .true.
        if (present(min_residual)) normalize = min_residual <= 0.0_rp

        ! generate the Gram matrix
        allocate(metric(n_vecs, n_vecs))
        metric = 0.0_rp
        call dsyrk("U", "T", n_vecs, vec_len, 1.0_rp, vecs, vec_len, 0.0_rp, metric, &
                   n_vecs)

        ! determine the norm of every step
        allocate(col_norm(n_vecs), usable(n_vecs))
        do j = 1, n_vecs
            col_norm(j) = sqrt(metric(j, j))
        end do

        ! check whether any steps are numerically vanishing
        usable = col_norm > numerical_zero * maxval(col_norm)
        if (present(keep)) usable = usable .and. keep

        ! normalize the usable steps to a unit diagonal,when requested, so that the
        ! pivoting ranks them by independence rather than by magnitude, and zero out
        ! the rest so that they are rejected rather than admitted with a zero pivot
        do j = 1, n_vecs
            do i = 1, j - 1
                if (.not. (usable(i) .and. usable(j))) then
                    metric(i, j) = 0.0_rp
                else if (normalize) then
                    metric(i, j) = metric(i, j) / (col_norm(i) * col_norm(j))
                end if
            end do
            if (.not. usable(j)) then
                metric(j, j) = 0.0_rp
            else if (normalize) then
                metric(j, j) = 1.0_rp
            end if
        end do

        ! tolerance for linear dependencies, relative to the unit diagonal, or the
        ! squared shortest resolvable residual on the unnormalized Gram
        if (normalize) then
            tol = n_vecs * numerical_zero
        else
            tol = max(min_residual**2, n_vecs * numerical_zero * maxval(col_norm)**2)
        end if

        ! perform pivoted rank-revealing Cholesky
        allocate(piv(n_vecs), work(2 * n_vecs))
        call dpstrf("U", n_vecs, metric, n_vecs, piv, n_accepted, tol, work, info)

        ! shrink to the accepted block only and undo the normalization by scaling,
        ! where it was applied
        allocate(chol(n_accepted, n_accepted), map(n_accepted))
        chol = metric(1:n_accepted, 1:n_accepted)
        map = piv(1:n_accepted)
        if (normalize) then
            do j = 1, n_accepted
                chol(:, j) = chol(:, j) * col_norm(map(j))
            end do
        end if
        deallocate(metric, piv, work, col_norm, usable)

    end subroutine factorize_history

    function rebase_dirs(dirs, map, chol) result(rebased)
        !
        ! this function re-expresses a set of packed history-direction columns in the
        ! orthonormalized basis by selecting and reordering the accepted columns via
        ! map, then applying a single triangular solve by the Cholesky factor,
        ! equivalent to right-multiplying by the inverse of the (implicit)
        ! upper-triangular factor without ever forming that inverse explicitly
        !
        real(rp), intent(in) :: dirs(:, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: rebased(:, :)

        integer(ip) :: n_param, n_accepted, j
        external :: dtrsm

        n_param = size(dirs, 1)
        n_accepted = size(map)
        allocate(rebased(n_param, n_accepted))
        do j = 1, n_accepted
            rebased(:, j) = dirs(:, map(j))
        end do

        if (n_accepted > 0) call dtrsm("R", "U", "N", "N", n_param, n_accepted, &
                                       1.0_rp, chol, n_accepted, rebased, n_param)

    end function rebase_dirs

    function congruence_transform(a, map, chol) result(a_tilde)
        !
        ! this function applies the congruence transformation A -> R^-T (P A P^T) R^-1,
        ! where P selects and reorders the rows/columns of A according to map and R is
        ! the Cholesky factor
        !
        real(rp), intent(in) :: a(:, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: a_tilde(:, :)

        integer(ip) :: n_accepted, i, j
        external :: dtrsm

        n_accepted = size(map)
        allocate(a_tilde(n_accepted, n_accepted))
        do j = 1, n_accepted
            do i = 1, n_accepted
                a_tilde(i, j) = a(map(i), map(j))
            end do
        end do

        if (n_accepted > 0) then
            call dtrsm("R", "U", "N", "N", n_accepted, n_accepted, 1.0_rp, chol, &
                       n_accepted, a_tilde, n_accepted)
            call dtrsm("L", "U", "T", "N", n_accepted, n_accepted, 1.0_rp, chol, &
                       n_accepted, a_tilde, n_accepted)
        end if

    end function congruence_transform

    subroutine combine_channels(chol1, map1, chol2, map2, n_offset, chol_comb, map_comb)
        !
        ! this subroutine assembles a block-diagonal Cholesky factor and concatenated,
        ! offset index map from two independent per-channel history factorizations, so
        ! that a single rebase or congruence-transform call can act on data that
        ! combines both channels
        !
        real(rp), intent(in) :: chol1(:, :), chol2(:, :)
        integer(ip), intent(in) :: map1(:), map2(:), n_offset
        real(rp), intent(out), allocatable :: chol_comb(:, :)
        integer(ip), intent(out), allocatable :: map_comb(:)

        integer(ip) :: n1, n2, n_total

        n1 = size(map1)
        n2 = size(map2)
        n_total = n1 + n2
        allocate(chol_comb(n_total, n_total), map_comb(n_total))
        chol_comb = 0.0_rp
        chol_comb(1:n1, 1:n1) = chol1
        chol_comb(n1 + 1:n_total, n1 + 1:n_total) = chol2
        map_comb(1:n1) = map1
        map_comb(n1 + 1:n_total) = n_offset + map2

    end subroutine combine_channels

    subroutine get_ms_a_inv(dm_cols, v_cols, map, chol, a_inv, settings, error)
        !
        ! this subroutine computes the pseudoinverse multisecant SR1 matrix for a
        ! single part of a coupled potential-difference response, from the symmetrized
        ! A of that part congruence-transformed into its own orthonormalized S-basis
        !
        use opentrustregion, only: symm_mat_diag

        real(rp), intent(in) :: dm_cols(:, :), v_cols(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), intent(in) :: chol(:, :)
        real(rp), intent(out), allocatable :: a_inv(:, :)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, n_accepted
        real(rp), allocatable :: a(:, :), a_tilde(:, :), eig_vecs(:, :), eig_vals(:), &
                                 eig_vals_inv(:)
        real(rp), allocatable :: y_gram(:, :)
        real(rp) :: eig_val_thresh

        ! initialize error flag
        error = 0

        ! handle empty history
        n_dm = size(dm_cols, 2, kind=ip)
        n_accepted = size(map)
        if (n_dm == 0 .or. n_accepted == 0) then
            allocate(a_inv(n_accepted, n_accepted))
            a_inv = 0.0_rp
            return
        end if

        ! build and symmetrize A = S^T Y
        allocate(a(n_dm, n_dm))
        call build_a_part(dm_cols, v_cols, a)

        ! congruence-transform to the orthonormalized, rank-independent S-basis
        a_tilde = congruence_transform(a, map, chol)
        deallocate(a)

        ! perform spectral decomposition
        allocate(eig_vecs(n_accepted, n_accepted), eig_vals(n_accepted))
        call symm_mat_diag(a_tilde, eig_vals, eig_vecs, settings, error)
        if (error /= 0) return

        ! construct inverse, discarding only eigenvalues at the level of numerical
        ! noise relative to the largest one
        allocate(eig_vals_inv(n_accepted))
        eig_val_thresh = eig_val_noise_factor * maxval(abs(eig_vals)) * epsilon(1.0_rp)
        eig_vals_inv = truncated_eigval_inv(eig_vals, eig_val_thresh)

        ! discard directions failing the MS-SR1 skipping criterion
        y_gram = response_gram(v_cols, map, chol)
        call apply_ms_sr1_skip(eig_vals, eig_vecs, y_gram, eig_vals_inv)
        deallocate(y_gram)

        ! reassemble the pseudoinverse
        a_inv = spectral_to_dense(eig_vecs, eig_vals_inv)
        deallocate(a_tilde, eig_vecs, eig_vals, eig_vals_inv)

    end subroutine get_ms_a_inv

    subroutine get_ms_a_inv_jk_cs(dm_cols, v_coulomb_cols, v_exchange_cols, &
                                  v_coulomb_packed, v_exchange_packed, map, chol, &
                                  a_inv, dirs, settings, error)
        !
        ! this subroutine computes the closed-shell multisecant SR1 term of the linear
        ! potential as the difference Y_J A_J^+ Y_J^T - Y_K A_K^+ Y_K^T of a Coulomb
        ! and an exact-exchange term with A = S^T Y, instead of a single term of the
        ! indefinite linear kernel; returns the coupling and the matching packed
        ! potential difference directions rebased into the orthonormalized S-basis
        !
        real(rp), intent(in) :: dm_cols(:, :), v_coulomb_cols(:, :), &
                                v_exchange_cols(:, :), v_coulomb_packed(:, :), &
                                v_exchange_packed(:, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), intent(out), allocatable :: a_inv(:, :), dirs(:, :)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm
        real(rp), allocatable :: a_coulomb(:, :), a_exchange(:, :)

        ! A = S^T Y of both kernels on the full history
        n_dm = size(dm_cols, 2, kind=ip)
        allocate(a_coulomb(n_dm, n_dm), a_exchange(n_dm, n_dm))
        call build_a_part(dm_cols, v_coulomb_cols, a_coulomb)
        call build_a_part(dm_cols, v_exchange_cols, a_exchange)

        ! pseudoinverses and directions of both terms
        call get_ms_jk_inv( &
            a_coulomb, a_exchange, rebase_dirs(v_coulomb_packed, map, chol), &
            rebase_dirs(v_exchange_packed, map, chol), all(v_exchange_cols == 0.0_rp), &
            map, chol, a_inv, dirs, settings, error)

    end subroutine get_ms_a_inv_jk_cs

    subroutine get_ms_a_inv_jk_os(dm_cols, v_coulomb_alpha_cols, v_coulomb_beta_cols, &
                                  v_exchange_cols, v_coulomb_alpha_packed, &
                                  v_coulomb_beta_packed, v_exchange_packed, rows, &
                                  packed_rows, map, chol, a_inv, dirs, settings, error)
        !
        ! this subroutine computes the open-shell multisecant SR1 term of the linear
        ! potential as the difference of a Coulomb and an exact-exchange term on the
        ! spin-separated density difference directions, those of every channel with the
        ! other channel set to zero, whose per-channel factorizations map and chol
        ! combine; the Coulomb potential of the alpha and beta density differences acts
        ! on both channels, which couples the directions of the two channels, the
        ! exact-exchange potential acts on the channel of its density difference only;
        ! returns the coupling and the matching packed potential difference directions
        ! rebased into the orthonormalized S-basis
        !
        real(rp), intent(in) :: &
            dm_cols(:, :), v_coulomb_alpha_cols(:, :), v_coulomb_beta_cols(:, :), &
            v_exchange_cols(:, :), v_coulomb_alpha_packed(:, :), &
            v_coulomb_beta_packed(:, :), v_exchange_packed(:, :), chol(:, :)
        integer(ip), intent(in) :: rows(:, :), packed_rows(:, :), map(:)
        real(rp), intent(out), allocatable :: a_inv(:, :), dirs(:, :)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_dm, n1, n2, i, k
        real(rp), allocatable :: a_coulomb(:, :), a_exchange(:, :), raw(:, :), &
                                 dirs_exchange(:, :)
        real(rp), external :: ddot

        ! number of history entries and of rows of the two spin channels in the history
        ! columns
        n_dm = size(dm_cols, 2, kind=ip)
        n1 = rows(2, 1) - rows(1, 1) + 1
        n2 = rows(2, 2) - rows(1, 2) + 1

        ! A = S^T Y of both kernels on the spin-separated directions
        allocate(a_coulomb(2 * n_dm, 2 * n_dm), a_exchange(2 * n_dm, 2 * n_dm))
        a_exchange = 0.0_rp
        do k = 1, n_dm
            do i = 1, n_dm
                a_coulomb(i, k) = ddot(n1, dm_cols(rows(1, 1):rows(2, 1), i), 1_ip, &
                                       v_coulomb_alpha_cols(rows(1, 1):rows(2, 1), k), &
                                       1_ip)
                a_coulomb(i, n_dm + k) = ddot( &
                    n1, dm_cols(rows(1, 1):rows(2, 1), i), 1_ip, &
                    v_coulomb_beta_cols(rows(1, 1):rows(2, 1), k), 1_ip)
                a_coulomb(n_dm + i, k) = ddot( &
                    n2, dm_cols(rows(1, 2):rows(2, 2), i), 1_ip, &
                    v_coulomb_alpha_cols(rows(1, 2):rows(2, 2), k), 1_ip)
                a_coulomb(n_dm + i, n_dm + k) = ddot( &
                    n2, dm_cols(rows(1, 2):rows(2, 2), i), 1_ip, &
                    v_coulomb_beta_cols(rows(1, 2):rows(2, 2), k), 1_ip)
                a_exchange(i, k) = ddot(n1, dm_cols(rows(1, 1):rows(2, 1), i), 1_ip, &
                                        v_exchange_cols(rows(1, 1):rows(2, 1), k), 1_ip)
                a_exchange(n_dm + i, n_dm + k) = ddot( &
                    n2, dm_cols(rows(1, 2):rows(2, 2), i), 1_ip, &
                    v_exchange_cols(rows(1, 2):rows(2, 2), k), 1_ip)
            end do
        end do

        ! the Coulomb responses of the alpha and the beta directions fill both
        ! channels, the exact-exchange responses stay in their own channel
        allocate(raw(size(v_coulomb_alpha_packed, 1), 2 * n_dm))
        raw(:, :n_dm) = v_coulomb_alpha_packed
        raw(:, n_dm + 1:) = v_coulomb_beta_packed
        call cache_channel_split_dirs(v_exchange_packed, packed_rows, map, chol, &
                                      dirs_exchange)

        ! pseudoinverses and directions of both terms
        call get_ms_jk_inv(a_coulomb, a_exchange, rebase_dirs(raw, map, chol), &
                           dirs_exchange, all(v_exchange_cols == 0.0_rp), map, chol, &
                           a_inv, dirs, settings, error)

    end subroutine get_ms_a_inv_jk_os

    subroutine get_ms_jk_inv(a_coulomb, a_exchange, dirs_coulomb, dirs_exchange, &
                             no_exchange, map, chol, a_inv, dirs, settings, error)
        !
        ! this subroutine returns the coupling and directions of the multisecant SR1
        ! term of the linear potential as the difference of a Coulomb and an
        ! exact-exchange term: both kernels are positive semidefinite, so that each
        ! term lies between zero and its kernel and a small eigenvalue of its A only
        ! comes with a small response, which a single term of the indefinite linear
        ! kernel cannot guarantee; both A are congruence-transformed into the
        ! orthonormalized S-basis and inverted on their eigenvalues above a noise floor
        ! relative to the larger of both spectra, so that an exact-exchange A of
        ! rounding noise is never inverted; the coupling is diag(A_J^+, -A_K^+) for the
        ! directions [Y_J, Y_K], or A_J^+ for Y_J alone without exact exchange
        !
        use opentrustregion, only: symm_mat_diag

        real(rp), intent(in) :: a_coulomb(:, :), a_exchange(:, :), dirs_coulomb(:, :), &
                                dirs_exchange(:, :), chol(:, :)
        logical, intent(in) :: no_exchange
        integer(ip), intent(in) :: map(:)
        real(rp), intent(out), allocatable :: a_inv(:, :), dirs(:, :)
        type(arh_settings_type), intent(in) :: settings
        integer(ip), intent(out) :: error

        integer(ip) :: n_accepted, n_param, i
        real(rp) :: eig_val_thresh
        real(rp), allocatable :: a_tilde_coulomb(:, :), a_tilde_exchange(:, :), &
                                 eig_vals_coulomb(:), eig_vals_exchange(:), &
                                 eig_vecs_coulomb(:, :), eig_vecs_exchange(:, :), &
                                 eig_vals_inv_coulomb(:), eig_vals_inv_exchange(:)

        ! initialize error flag
        error = 0

        ! handle empty history
        n_accepted = size(map)
        n_param = size(dirs_coulomb, 1, kind=ip)
        if (n_accepted == 0) then
            allocate(a_inv(0, 0), dirs(n_param, 0))
            return
        end if

        ! congruence-transform both A to the orthonormalized, rank-independent S-basis,
        ! enforcing the symmetry they have in exact arithmetic
        a_tilde_coulomb = congruence_transform(a_coulomb, map, chol)
        a_tilde_coulomb = 0.5_rp * (a_tilde_coulomb + transpose(a_tilde_coulomb))
        a_tilde_exchange = congruence_transform(a_exchange, map, chol)
        a_tilde_exchange = 0.5_rp * (a_tilde_exchange + transpose(a_tilde_exchange))

        ! perform spectral decompositions
        allocate( &
            eig_vals_coulomb(n_accepted), eig_vecs_coulomb(n_accepted, n_accepted), &
            eig_vals_exchange(n_accepted), eig_vecs_exchange(n_accepted, n_accepted))
        call symm_mat_diag(a_tilde_coulomb, eig_vals_coulomb, eig_vecs_coulomb, &
                           settings, error)
        if (error /= 0) return
        call symm_mat_diag(a_tilde_exchange, eig_vals_exchange, eig_vecs_exchange, &
                           settings, error)
        if (error /= 0) return

        ! invert the eigenvalues above the noise floor, which positive semidefinite
        ! matrices only fall below by rounding
        eig_val_thresh = eig_val_noise_factor * epsilon(1.0_rp) * max( &
            maxval(abs(eig_vals_coulomb)), maxval(abs(eig_vals_exchange)))
        allocate(eig_vals_inv_coulomb(n_accepted), eig_vals_inv_exchange(n_accepted))
        eig_vals_inv_coulomb = 0.0_rp
        eig_vals_inv_exchange = 0.0_rp
        do i = 1, n_accepted
            if (eig_vals_coulomb(i) > eig_val_thresh) &
                eig_vals_inv_coulomb(i) = 1.0_rp / eig_vals_coulomb(i)
            if (eig_vals_exchange(i) > eig_val_thresh) &
                eig_vals_inv_exchange(i) = 1.0_rp / eig_vals_exchange(i)
        end do

        ! assemble the coupling and the directions, without the exact-exchange term if
        ! there is no exact exchange
        if (no_exchange) then
            a_inv = spectral_to_dense(eig_vecs_coulomb, eig_vals_inv_coulomb)
            dirs = dirs_coulomb
        else
            allocate(a_inv(2 * n_accepted, 2 * n_accepted), &
                     dirs(n_param, 2 * n_accepted))
            a_inv = 0.0_rp
            a_inv(:n_accepted, :n_accepted) = &
                spectral_to_dense(eig_vecs_coulomb, eig_vals_inv_coulomb)
            a_inv(n_accepted + 1:, n_accepted + 1:) = &
                -spectral_to_dense(eig_vecs_exchange, eig_vals_inv_exchange)
            dirs(:, :n_accepted) = dirs_coulomb
            dirs(:, n_accepted + 1:) = dirs_exchange
        end if

    end subroutine get_ms_jk_inv

    function response_gram(v_cols, map, chol) result(y_gram)
        !
        ! this function returns the Gram matrix of the response history rebased into
        ! the orthonormalized S-basis, so that the squared norm of the response along
        ! an eigenvector v of the congruence-transformed A is v^T y_gram v
        !
        real(rp), intent(in) :: v_cols(:, :), chol(:, :)
        integer(ip), intent(in) :: map(:)
        real(rp), allocatable :: y_gram(:, :)

        integer(ip) :: n_dm, flat_len
        real(rp), allocatable :: gram(:, :)
        external :: dgemm

        n_dm = size(v_cols, 2, kind=ip)
        flat_len = size(v_cols, 1, kind=ip)
        allocate(gram(n_dm, n_dm))
        call dgemm("T", "N", n_dm, n_dm, flat_len, 1.0_rp, v_cols, flat_len, v_cols, &
                   flat_len, 0.0_rp, gram, n_dm)
        y_gram = congruence_transform(gram, map, chol)
        deallocate(gram)

    end function response_gram

    subroutine apply_ms_sr1_skip(eig_vals, eig_vecs, y_gram, eig_vals_inv)
        !
        ! this subroutine applies the multisecant analogue of the SR1 skipping
        ! criterion: an eigenvalue of the congruence-transformed A plays the role of
        ! the classical denominator s^T (y - H s), and a direction is discarded when it
        ! is small against the norm of the response it divides; the direction has unit
        ! norm in the orthonormalized S-basis, so the ratio is a dimensionless angle
        !
        real(rp), intent(in) :: eig_vals(:), eig_vecs(:, :), y_gram(:, :)
        real(rp), intent(inout) :: eig_vals_inv(:)

        integer(ip) :: i, n
        real(rp) :: y_norm
        real(rp), allocatable :: temp(:)
        real(rp), external :: ddot
        external :: dgemv

        n = size(eig_vals)
        allocate(temp(n))
        do i = 1, n
            call dgemv("N", n, n, 1.0_rp, y_gram, n, eig_vecs(:, i), 1_ip, 0.0_rp, &
                       temp, 1_ip)
            y_norm = sqrt(max(ddot(n, eig_vecs(:, i), 1_ip, temp, 1_ip), 0.0_rp))
            if (abs(eig_vals(i)) < ms_sr1_skip_thresh * y_norm) eig_vals_inv(i) = 0.0_rp
        end do
        deallocate(temp)

    end subroutine apply_ms_sr1_skip

    function spectral_to_dense(eigvecs, inv_eigvals) result(mat)
        !
        ! this function reconstructs the dense symmetric matrix V * diag(lambda) * V^T
        ! from its eigenvectors and eigenvalues
        !
        real(rp), intent(in) :: eigvecs(:, :), inv_eigvals(:)
        real(rp), allocatable :: mat(:, :)

        integer(ip) :: n, i
        real(rp), allocatable :: scaled(:, :)
        external :: dgemm

        n = size(eigvecs, 1)
        allocate(scaled(n, n), mat(n, n))
        scaled = eigvecs
        do i = 1, n
            scaled(:, i) = scaled(:, i) * inv_eigvals(i)
        end do
        call dgemm("N", "T", n, n, n, 1.0_rp, scaled, n, eigvecs, n, 0.0_rp, mat, n)
        deallocate(scaled)

    end function spectral_to_dense

    function density_in_history(dm) result(in_history)
        !
        ! this function reports whether a density is already held in the history
        !
        use opentrustregion, only: numerical_zero

        real(rp), intent(in) :: dm(:, :, :)
        logical :: in_history

        integer(ip) :: i
        real(rp) :: dm_scale

        in_history = .false.
        if (.not. allocated(arh_object%dm_list)) return

        ! judge the difference against the size of the density, so that the test is a
        ! relative one
        dm_scale = max(maxval(abs(dm)), numerical_zero)
        do i = 1, size(arh_object%dm_list, 4)
            if (maxval(abs(arh_object%dm_list(:, :, :, i) - dm)) <= &
                numerical_zero * dm_scale) then
                in_history = .true.
                return
            end if
        end do

    end function density_in_history

    subroutine prepend(list, new_array)
        !
        ! this subroutine prepends an array to a list of arrays of equal dimension
        !
        real(rp), intent(inout), allocatable :: list(:, :, :, :)
        real(rp), intent(in) :: new_array(:, :, :)

        integer(ip) :: n1, n2, n3
        real(rp), allocatable :: temp(:, :, :, :)

        n1 = size(new_array, 1)
        n2 = size(new_array, 2)
        n3 = size(new_array, 3)

        allocate(temp(n1, n2, n3, size(list, 4) + 1))
        temp(:, :, :, 1) = new_array
        temp(:, :, :, 2:size(list, 4) + 1) = list
        deallocate(list)
        list = temp

    end subroutine prepend

end module otr_arh
