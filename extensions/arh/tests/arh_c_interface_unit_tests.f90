! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_c_interface_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use c_interface, only: c_rp, c_ip
    use test_reference, only: tol
    use otr_arh_c_interface, only: evaluate_dm_c_type
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_ptr, c_funptr, &
                                           c_funloc, c_f_procpointer, c_loc

    implicit none

    ! create function pointers to ensure that routines comply with interface
    procedure(evaluate_dm_c_type), pointer :: mock_arh_evaluate_dm_cs_ptr => &
        mock_arh_evaluate_dm_cs
    procedure(evaluate_dm_c_type), pointer :: mock_arh_evaluate_dm_os_ptr => &
        mock_arh_evaluate_dm_os

contains

    function mock_arh_evaluate_dm_cs(dm_ao, energy, fock, v_coulomb, v_exchange, &
                                     v_nonlinear, context_c) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function with separate Coulomb, exact-exchange and non-linear potential
        ! contributions for the closed-shell case
        !
        use otr_common_test_reference, only: n_ao
        use otr_arh_test_reference, only: arh_evaluate_dm_factors
        use otr_common_unit_tests, only: record_mock_call
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*), v_coulomb(*), v_exchange(*), &
                                             v_nonlinear(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2

        ! check host context
        call check_host_context_c(context_c)

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + merge( &
            2_ip, 0_ip, present(v_coulomb)) + merge(4_ip, 0_ip, present(v_exchange)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) &
            fock(:flat_len) = arh_evaluate_dm_factors(1) * dm_ao(:flat_len)
        if (present(v_coulomb)) &
            v_coulomb(:flat_len) = arh_evaluate_dm_factors(2) * dm_ao(:flat_len)
        if (present(v_exchange)) &
            v_exchange(:flat_len) = arh_evaluate_dm_factors(3) * dm_ao(:flat_len)
        if (present(v_nonlinear)) &
            v_nonlinear(:flat_len) = arh_evaluate_dm_factors(4) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_arh_evaluate_dm_cs

    function mock_arh_evaluate_dm_os(dm_ao, energy, fock, v_coulomb, v_exchange, &
                                     v_nonlinear, context_c) result(error) bind(C)
        !
        ! this subroutine is a test subroutine for the density matrix evaluating C
        ! function with separate Coulomb, exact-exchange and non-linear potential
        ! contributions for the open-shell case
        !
        use otr_common_test_reference, only: n_ao, n_particle
        use otr_arh_test_reference, only: arh_evaluate_dm_factors
        use otr_common_unit_tests, only: record_mock_call
        use test_reference, only: check_host_context_c

        real(c_rp), intent(in), target :: dm_ao(*)
        real(c_rp), intent(out) :: energy
        real(c_rp), intent(out), optional :: fock(*), v_coulomb(*), v_exchange(*), &
                                             v_nonlinear(*)
        type(c_ptr), intent(in), value :: context_c
        integer(c_ip) :: error

        integer(c_ip) :: flat_len = n_ao**2 * n_particle

        ! check host context
        call check_host_context_c(context_c)

        call record_mock_call(merge(1_ip, 0_ip, present(fock)) + merge( &
            2_ip, 0_ip, present(v_coulomb)) + merge(4_ip, 0_ip, present(v_exchange)) + &
                              merge(8_ip, 0_ip, present(v_nonlinear)))
        energy = sum(dm_ao(:flat_len))
        if (present(fock)) &
            fock(:flat_len) = arh_evaluate_dm_factors(1) * dm_ao(:flat_len)
        if (present(v_coulomb)) &
            v_coulomb(:flat_len) = arh_evaluate_dm_factors(2) * dm_ao(:flat_len)
        if (present(v_exchange)) &
            v_exchange(:flat_len) = arh_evaluate_dm_factors(3) * dm_ao(:flat_len)
        if (present(v_nonlinear)) &
            v_nonlinear(:flat_len) = arh_evaluate_dm_factors(4) * dm_ao(:flat_len)

        error = 0_c_ip

    end function mock_arh_evaluate_dm_os

    logical(c_bool) function test_arh_factory_mo_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH factory for orbitals
        ! parameterized in the MO basis for the closed- and the open-shell case,
        ! without and with irreps of the MOs
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, assignment(=)
        use otr_arh, only: arh_object, arh_n_micro
        use otr_arh_c_interface, only: arh_callbacks, arh_settings_type_c, &
                                       arh_factory_mo_cs, arh_factory_mo_os, &
                                       arh_factory_mo_c_wrapper, arh_holder
        use otr_arh_mock, only: mock_arh_factory_mo_cs, mock_arh_factory_mo_os, &
                                mock_arh_factory_error, test_passed
        use otr_mo_mock, only: orbsym_passed
        use otr_mo_c_interface, only: update_orbs_mo_c_wrapper
        use otr_arh_test_reference, only: assignment(=), ref_arh_settings
        use otr_mo_test_reference, only: mo_coeff_pattern, n_mo, case_irreps, case_names
        use otr_common_test_reference, only: n_ao, n_particle, n_occ, n_ao_c, &
                                             shell_names
        use otr_common_c_interface_unit_tests, only: check_c_factory_wiring
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: host_context_type, arm_host_context_c, host_context, &
                                  host_context_reached

        real(c_rp), allocatable :: ao_overlap_c(:, :), mo_coeff_c(:, :, :)
        type(c_funptr) :: evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr
        type(arh_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        integer(c_ip) :: n_particle_c, error_c, n_mo_c, n_occ_c(n_particle)
        integer(c_ip), allocatable :: orbsym_c(:)
        integer(ip) :: irreps(n_mo, n_particle), i_call, i_sym
        logical :: with_irreps
        character(len=:), allocatable :: case_name
        type(host_context_type), target :: other_host_context

        ! assume tests pass
        test_arh_factory_mo_c_wrapper = .true.

        ! inject mock functions
        arh_factory_mo_cs => mock_arh_factory_mo_cs
        arh_factory_mo_os => mock_arh_factory_mo_os

        ! dimensions differ between AOs and MOs and occupations between the spin
        ! channels, so that the shapes are checked
        n_mo_c = int(n_mo, kind=c_ip)
        n_occ_c = int(n_occ, kind=c_ip)

        ! initialize AO overlap matrix
        allocate(ao_overlap_c(n_ao, n_ao))
        ao_overlap_c = 2.0_c_rp

        ! associate optional settings with values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! both spin cases pass through the same C wrapper, first without and then with
        ! the irreps of the MOs of the symmetric occupation case of the shell
        do i_call = 1, 2 * n_particle
            n_particle_c = int(1 + mod(i_call - 1, n_particle), kind=c_ip)
            with_irreps = i_call > n_particle
            case_name = trim(shell_names(n_particle_c))
            if (with_irreps) case_name = case_name//" with symmetry"

            ! initialize MO coefficients with values encoding their indices
            mo_coeff_c = mo_coeff_pattern(n_ao_c, n_mo_c, n_particle_c, 0.0_c_rp)

            ! get C function pointers to Fortran functions
            if (n_particle_c == 1) then
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_cs)
            else
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_os)
            end if

            ! arm the host context of the ARH settings, which the logging and density
            ! matrix evaluating functions receive
            call arm_host_context_c(settings_c%context)
            solver_settings_c = default_solver_settings

            ! set the irreps of the MOs passed to the wrapper, which an unallocated
            ! array passes as absent
            irreps = 0
            if (allocated(orbsym_c)) deallocate(orbsym_c)
            if (with_irreps) then
                i_sym = findloc(case_names, case_name, dim=1, kind=ip)
                irreps(:, :n_particle_c) = case_irreps(:, :n_particle_c, i_sym)
                orbsym_c = int( &
                    reshape(irreps(:, :n_particle_c), [n_mo * n_particle_c]), kind=c_ip)
            end if

            ! call ARH factory C wrapper
            test_passed = .false.
            mock_arh_factory_error = 0
            error_c = arh_factory_mo_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
                evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr, &
                solver_settings_c, settings_c, orbsym_c)

            ! check if the mock factory received the correct input and called the
            ! logging function
            if (.not. test_passed) then
                test_arh_factory_mo_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Mock "// &
                    "factory received wrong input for the "//case_name//" case."
            end if
            if (.not. host_context%logger_called) then
                test_arh_factory_mo_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Called "// &
                    "logging subroutine wrong for the "//case_name//" case."
            end if
            if (orbsym_passed .neqv. with_irreps) then
                test_arh_factory_mo_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Irreps of "// &
                    "the MOs passed on wrongly for the "//case_name//" case."
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("arh_factory_mo_c_wrapper")) &
                test_arh_factory_mo_c_wrapper = .false.

            ! check if output variables are as expected
            if (error_c /= 0) then
                test_arh_factory_mo_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Returned "// &
                    "error code wrong for the "//case_name//" case."
            end if

            ! check that the factory received the MO coefficients of the caller in
            ! place, unless the kinds differ and the object updates a copy
            if (rp == c_rp .and. any(abs(mo_coeff_c - 2.0_c_rp) > tol)) then
                test_arh_factory_mo_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: MO "// &
                    "coefficients of the caller not passed in place for the "// &
                    case_name//" case."
            end if

            ! check that the factory returned the C wrappers of the object, wired the
            ! solver settings with its holder and adjusted them to the approximate
            ! Hessian
            if (.not. check_c_factory_wiring("arh_factory_mo_c_wrapper", arh_holder, &
                                             arh_object, update_orbs_c_funptr, &
                                             obj_func_c_funptr, solver_settings_c, &
                                             c_funloc(update_orbs_mo_c_wrapper))) &
                test_arh_factory_mo_c_wrapper = .false.
            if (.not. solver_settings_c%refresh_hess .or. &
                solver_settings_c%n_micro /= arh_n_micro .or. &
                solver_settings_c%hess_symm) then
                test_arh_factory_mo_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Solver "// &
                    "settings not adjusted to the approximate Hessian for the "// &
                    case_name//" case."
            end if

            ! call the factory with the solver settings it wired and another host
            ! context, which has to produce an error without calling the factory,
            ! returning functions or replacing the bundle holding the callback functions
            ! of the object
            settings_c%context = c_loc(other_host_context)
            error_c = arh_factory_mo_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
                evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr, &
                solver_settings_c, settings_c, orbsym_c)
            if (error_c == 0 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. .not. c_associated( &
                    arh_callbacks%host_context, c_loc(host_context))) then
                test_arh_factory_mo_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Solver "// &
                    "settings the holder was wired into not refused before the "// &
                    "object was changed for the "//case_name//" case."
            end if

            ! call the factory for the initial MO coefficients with a mock factory
            ! which produces an error, which has to be returned without functions or
            ! writing back the solver settings or changing the holder or the bundle of
            ! the object
            mo_coeff_c = mo_coeff_pattern(n_ao_c, n_mo_c, n_particle_c, 0.0_c_rp)
            solver_settings_c = default_solver_settings
            update_orbs_c_funptr = c_funloc(mock_logger_c)
            obj_func_c_funptr = c_funloc(mock_logger_c)
            mock_arh_factory_error = 1
            error_c = arh_factory_mo_c_wrapper( &
                mo_coeff_c, ao_overlap_c, n_occ_c, n_particle_c, n_ao_c, n_mo_c, &
                evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr, &
                solver_settings_c, settings_c, orbsym_c)
            mock_arh_factory_error = 0
            if (error_c /= 1 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. &
                c_associated(solver_settings_c%context)) then
                test_arh_factory_mo_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Error of "// &
                    "the factory not passed on for the "//case_name//" case."
            end if
            if (.not. associated(arh_holder%solver_level%level, &
                                 arh_object%solver_level) .or. &
                .not. c_associated(arh_callbacks%host_context, c_loc(host_context))) &
                then
                write(stderr, *) "test_arh_factory_mo_c_wrapper failed: Holder or "// &
                    "callback functions of the object not restored after an error "// &
                    "of the factory for the "//case_name//" case."
                test_arh_factory_mo_c_wrapper = .false.
            end if
        end do
        deallocate(ao_overlap_c)
        if (allocated(arh_object)) deallocate(arh_object)

    end function test_arh_factory_mo_c_wrapper

    logical(c_bool) function test_arh_factory_oao_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH factory for orbitals
        ! parameterized in the OAO basis for the closed- and the open-shell case
        !
        use opentrustregion, only: default_solver_settings
        use c_interface, only: solver_settings_type_c, assignment(=)
        use otr_arh, only: arh_object, arh_n_micro
        use otr_arh_c_interface, only: arh_callbacks, arh_settings_type_c, &
                                       arh_factory_oao_cs, arh_factory_oao_os, &
                                       arh_factory_oao_c_wrapper, arh_holder
        use otr_arh_mock, only: mock_arh_factory_oao_cs, mock_arh_factory_oao_os, &
                                mock_arh_factory_error, test_passed
        use otr_arh_test_reference, only: assignment(=), ref_arh_settings
        use otr_oao_c_interface, only: update_orbs_oao_c_wrapper
        use otr_common_test_reference, only: n_ao, n_particle, n_ao_c, shell_names
        use otr_common_c_interface_unit_tests, only: check_c_factory_wiring
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: host_context_type, arm_host_context_c, host_context, &
                                  host_context_reached

        real(c_rp), allocatable :: ao_overlap_c(:, :), dm_ao_c(:, :, :)
        type(c_funptr) :: evaluate_dm_c_funptr, update_orbs_c_funptr, obj_func_c_funptr
        type(arh_settings_type_c) :: settings_c
        type(solver_settings_type_c) :: solver_settings_c
        integer(c_ip) :: n_particle_c, error_c
        character(len=:), allocatable :: case_name
        type(host_context_type), target :: other_host_context

        ! assume tests pass
        test_arh_factory_oao_c_wrapper = .true.

        ! inject mock functions
        arh_factory_oao_cs => mock_arh_factory_oao_cs
        arh_factory_oao_os => mock_arh_factory_oao_os

        ! initialize AO overlap matrix
        allocate(ao_overlap_c(n_ao, n_ao))
        ao_overlap_c = 2.0_c_rp

        ! associate optional settings with values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! both spin cases pass through the same C wrapper
        do n_particle_c = 1, int(n_particle, kind=c_ip)
            case_name = trim(shell_names(n_particle_c))

            ! initialize density matrix and get C function pointer to the density
            ! matrix evaluating function of the spin case
            allocate(dm_ao_c(n_ao, n_ao, n_particle_c))
            dm_ao_c = 1.0_c_rp
            if (n_particle_c == 1) then
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_cs)
            else
                evaluate_dm_c_funptr = c_funloc(mock_arh_evaluate_dm_os)
            end if

            ! arm the host context of the ARH settings, which the logging and density
            ! matrix evaluating functions receive
            call arm_host_context_c(settings_c%context)
            solver_settings_c = default_solver_settings

            ! call ARH factory C wrapper
            test_passed = .false.
            mock_arh_factory_error = 0
            error_c = arh_factory_oao_c_wrapper( &
                dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
                update_orbs_c_funptr, obj_func_c_funptr, solver_settings_c, settings_c)

            ! check if the mock factory received the correct input and called the
            ! logging function
            if (.not. test_passed) then
                test_arh_factory_oao_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Mock "// &
                    "factory received wrong input for the "//case_name//" case."
            end if
            if (.not. host_context%logger_called) then
                test_arh_factory_oao_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Called "// &
                    "logging subroutine wrong for the "//case_name//" case."
            end if

            ! check that the host context reached the callback functions unchanged
            if (.not. host_context_reached("arh_factory_oao_c_wrapper")) &
                test_arh_factory_oao_c_wrapper = .false.

            ! check if output variables are as expected
            if (error_c /= 0) then
                test_arh_factory_oao_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Returned "// &
                    "error code wrong for the "//case_name//" case."
            end if

            ! check that the factory received the density matrix of the caller in place,
            ! unless the kinds differ and the object updates a copy
            if (rp == c_rp .and. any(abs(dm_ao_c - 2.0_c_rp) > tol)) then
                test_arh_factory_oao_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Density "// &
                    "matrix of the caller not passed in place for the "//case_name// &
                    " case."
            end if

            ! check that the factory returned the C wrappers of the object, wired the
            ! solver settings with its holder and adjusted them to the approximate
            ! Hessian
            if (.not. check_c_factory_wiring("arh_factory_oao_c_wrapper", arh_holder, &
                                             arh_object, update_orbs_c_funptr, &
                                             obj_func_c_funptr, solver_settings_c, &
                                             c_funloc(update_orbs_oao_c_wrapper))) &
                test_arh_factory_oao_c_wrapper = .false.
            if (.not. solver_settings_c%refresh_hess .or. &
                solver_settings_c%n_micro /= arh_n_micro .or. &
                solver_settings_c%hess_symm) then
                test_arh_factory_oao_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Solver "// &
                    "settings not adjusted to the approximate Hessian for the "// &
                    case_name//" case."
            end if

            ! call the factory with the solver settings it wired and another host
            ! context, which has to produce an error without calling the factory,
            ! returning functions or replacing the bundle holding the callback functions
            ! of the object
            settings_c%context = c_loc(other_host_context)
            error_c = arh_factory_oao_c_wrapper( &
                dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
                update_orbs_c_funptr, obj_func_c_funptr, solver_settings_c, settings_c)
            if (error_c == 0 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. .not. c_associated( &
                    arh_callbacks%host_context, c_loc(host_context))) then
                test_arh_factory_oao_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Solver "// &
                    "settings the holder was wired into not refused before the "// &
                    "object was changed for the "//case_name//" case."
            end if

            ! call the factory for the initial density matrix with a mock factory which
            ! produces an error, which has to be returned without functions or writing
            ! back the solver settings or changing the holder or the bundle of the
            ! object
            dm_ao_c = 1.0_c_rp
            solver_settings_c = default_solver_settings
            update_orbs_c_funptr = c_funloc(mock_logger_c)
            obj_func_c_funptr = c_funloc(mock_logger_c)
            mock_arh_factory_error = 1
            error_c = arh_factory_oao_c_wrapper( &
                dm_ao_c, ao_overlap_c, n_particle_c, n_ao_c, evaluate_dm_c_funptr, &
                update_orbs_c_funptr, obj_func_c_funptr, solver_settings_c, settings_c)
            mock_arh_factory_error = 0
            if (error_c /= 1 .or. c_associated(update_orbs_c_funptr) .or. &
                c_associated(obj_func_c_funptr) .or. &
                c_associated(solver_settings_c%context)) then
                test_arh_factory_oao_c_wrapper = .false.
                write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Error of "// &
                    "the factory not passed on for the "//case_name//" case."
            end if
            if (.not. associated(arh_holder%solver_level%level, &
                                 arh_object%solver_level) .or. &
                .not. c_associated(arh_callbacks%host_context, c_loc(host_context))) &
                then
                write(stderr, *) "test_arh_factory_oao_c_wrapper failed: Holder or "// &
                    "callback functions of the object not restored after an error "// &
                    "of the factory for the "//case_name//" case."
                test_arh_factory_oao_c_wrapper = .false.
            end if
            deallocate(dm_ao_c)
        end do
        deallocate(ao_overlap_c)
        if (allocated(arh_object)) deallocate(arh_object)

    end function test_arh_factory_oao_c_wrapper

    logical(c_bool) function test_evaluate_dm_cs_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating C
        ! function with Coulomb, exact-exchange and non-linear potential contributions
        ! for the closed-shell case
        !
        use otr_arh, only: evaluate_dm_cs_type
        use otr_arh_c_interface, only: arh_callbacks, evaluate_dm_cs_f_wrapper
        use otr_arh_test_reference, only: test_evaluate_dm_cs_funptr
        use otr_common_unit_tests, only: mock_requests
        use test_reference, only: arm_host_context_c, host_context_reached

        procedure(evaluate_dm_cs_type), pointer :: evaluate_dm_cs_funptr
        class(*), pointer :: context
        integer(ip) :: request

        ! assume tests pass
        test_evaluate_dm_cs_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        arh_callbacks%evaluate_dm => mock_arh_evaluate_dm_cs
        call arm_host_context_c(arh_callbacks%host_context)
        context => arh_callbacks

        ! get pointer to subroutine
        evaluate_dm_cs_funptr => evaluate_dm_cs_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        if (.not. test_evaluate_dm_cs_funptr(evaluate_dm_cs_funptr, &
                                             "evaluate_dm_cs_f_wrapper", "", context)) &
            test_evaluate_dm_cs_f_wrapper = .false.
        if (size(mock_requests) /= 16 .or. &
            .not. all([(any(mock_requests == request), request=0, 15)])) then
            write(stderr, *) "test_evaluate_dm_cs_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for closed-shell case "// &
                "wrong."
            test_evaluate_dm_cs_f_wrapper = .false.
        end if

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("evaluate_dm_cs_f_wrapper")) &
            test_evaluate_dm_cs_f_wrapper = .false.

    end function test_evaluate_dm_cs_f_wrapper

    logical(c_bool) function test_evaluate_dm_os_f_wrapper() bind(C)
        !
        ! this function tests the Fortran wrapper for the density matrix evaluating C
        ! function with Coulomb, exact-exchange and non-linear potential contributions
        ! for the open-shell case
        !
        use otr_arh, only: evaluate_dm_os_type
        use otr_arh_c_interface, only: arh_callbacks, evaluate_dm_os_f_wrapper
        use otr_arh_test_reference, only: test_evaluate_dm_os_funptr
        use otr_common_unit_tests, only: mock_requests
        use test_reference, only: arm_host_context_c, host_context_reached

        procedure(evaluate_dm_os_type), pointer :: evaluate_dm_os_funptr
        class(*), pointer :: context
        integer(ip) :: request

        ! assume tests pass
        test_evaluate_dm_os_f_wrapper = .true.

        ! inject mock C function and hand the bundle to the wrapper as its context
        arh_callbacks%evaluate_dm => mock_arh_evaluate_dm_os
        call arm_host_context_c(arh_callbacks%host_context)
        context => arh_callbacks

        ! get pointer to subroutine
        evaluate_dm_os_funptr => evaluate_dm_os_f_wrapper

        ! test density matrix evaluating wrapper, which is called once for every
        ! combination of requested outputs and has to pass each on to the C function
        ! unchanged
        mock_requests = [integer(ip) :: ]
        if (.not. test_evaluate_dm_os_funptr(evaluate_dm_os_funptr, &
                                             "evaluate_dm_os_f_wrapper", "", context)) &
            test_evaluate_dm_os_f_wrapper = .false.
        if (size(mock_requests) /= 16 .or. &
            .not. all([(any(mock_requests == request), request=0, 15)])) then
            write(stderr, *) "test_evaluate_dm_os_f_wrapper failed: Outputs passed "// &
                "on to density matrix evaluating C function for open-shell case wrong."
            test_evaluate_dm_os_f_wrapper = .false.
        end if

        ! check that the wrapper handed the host context to the C function
        if (.not. host_context_reached("evaluate_dm_os_f_wrapper")) &
            test_evaluate_dm_os_f_wrapper = .false.

    end function test_evaluate_dm_os_f_wrapper

    logical(c_bool) function test_init_arh_settings_c() bind(C)
        !
        ! this function tests that the ARH settings initialization routine correctly
        ! initializes all settings to their default values
        !
        use otr_arh_c_interface, only: arh_settings_type_c, init_arh_settings_c
        use otr_arh, only: default_arh_settings
        use otr_arh_test_reference, only: operator(/=)

        type(arh_settings_type_c) :: settings

        ! assume test passes
        test_init_arh_settings_c = .true.

        ! initialize settings
        call init_arh_settings_c(settings)

        ! check settings
        if (settings /= default_arh_settings) then
            write(stderr, *) "test_init_arh_settings_c failed: Settings not "// &
                "initialized correctly."
            test_init_arh_settings_c = .false.
        end if

    end function test_init_arh_settings_c

    logical(c_bool) function test_arh_deconstructor_c_wrapper() bind(C)
        !
        ! this function tests the C wrapper for the ARH deconstructor
        !
        use otr_arh_c_interface, only: arh_deconstructor, arh_deconstructor_c_wrapper
        use otr_arh_mock, only: mock_arh_deconstructor, test_passed

        ! assume tests pass
        test_arh_deconstructor_c_wrapper = .true.

        ! inject mock functions
        arh_deconstructor => mock_arh_deconstructor

        ! initialize test logical
        test_passed = .false.

        ! call ARH orbital updating deconstructor C wrapper
        call arh_deconstructor_c_wrapper()

        ! check if test has passed
        test_arh_deconstructor_c_wrapper = test_passed

        ! check if test has passed
        if (.not. test_passed) then
            test_arh_deconstructor_c_wrapper = .false.
            write(stderr, *) "test_arh_deconstructor_c_wrapper failed: "// &
                "Deconstructor called wrong."
        end if

    end function test_arh_deconstructor_c_wrapper

    logical(c_bool) function test_assign_arh_f_c() bind(C)
        !
        ! this function tests that the function that converts ARH settings from C to
        ! Fortran correctly perform this conversion
        !
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh, only: arh_settings_type
        use otr_arh_test_reference, only: assignment(=), operator(/=), ref_arh_settings
        use c_interface, only: c_callbacks_type
        use otr_common_mock, only: mock_logger_c
        use test_reference, only: arm_host_context_c, host_context, host_context_reached

        type(arh_settings_type_c) :: settings_c
        type(arh_settings_type) :: settings
        type(c_callbacks_type), target :: callbacks

        ! assume test passes
        test_assign_arh_f_c = .true.

        ! initialize the C settings with custom values
        settings_c = ref_arh_settings
        settings_c%logger = c_funloc(mock_logger_c)

        ! convert to Fortran settings and hand them the C logging function in a bundle
        ! through context, as the factory does
        settings = settings_c
        callbacks%logger => mock_logger_c
        settings%context => callbacks

        ! check logging function
        if (.not. associated(settings%logger)) then
            test_assign_arh_f_c = .false.
            write(stderr, *) "test_assign_arh_f_c failed: Logging function not "// &
                "associated with value."
        else
            call arm_host_context_c(callbacks%host_context)
            call settings%logger("test", settings%context)
            if (.not. host_context%logger_called) then
                test_assign_arh_f_c = .false.
                write(stderr, *) "test_assign_arh_f_c failed: Called logging "// &
                    "subroutine wrong."
            end if
            if (.not. host_context_reached("assign_arh_f_c")) &
                test_assign_arh_f_c = .false.
        end if

        ! check against reference values
        if (settings /= ref_arh_settings) then
            write(stderr, *) "test_assign_arh_f_c failed: Settings not converted "// &
                "correctly."
            test_assign_arh_f_c = .false.
        end if

        ! check initialization flag
        if (.not. settings%initialized) then
            write(stderr, *) "test_assign_arh_f_c failed: Settings not marked as "// &
                "initialized."
            test_assign_arh_f_c = .false.
        end if

        ! convert C settings that were not initialized and check that they are left
        ! uninitialized without a logging function, so that the factory initializes
        ! them after its warning
        settings_c%initialized = .false.
        settings = settings_c
        if (settings%initialized .or. associated(settings%logger)) then
            write(stderr, *) "test_assign_arh_f_c failed: Settings that were not "// &
                "initialized converted."
            test_assign_arh_f_c = .false.
        end if

    end function test_assign_arh_f_c

    logical(c_bool) function test_assign_arh_c_f() bind(C)
        !
        ! this function tests that the function that converts ARH settings from Fortran
        ! to C correctly performs this conversion
        !
        use otr_arh, only: arh_settings_type
        use otr_arh_c_interface, only: arh_settings_type_c, assignment(=)
        use otr_arh_test_reference, only: assignment(=), operator(/=), ref_arh_settings

        type(arh_settings_type) :: settings
        type(arh_settings_type_c) :: settings_c

        ! assume test passes
        test_assign_arh_c_f = .true.

        ! initialize Fortran settings with reference values
        settings = ref_arh_settings

        ! convert to C settings
        settings_c = settings

        ! check that no callback function pointers and host contexts are associated
        if (c_associated(settings_c%logger)) then
            test_assign_arh_c_f = .false.
            write(stderr, *) "test_assign_arh_c_f failed: Logger function associated."
        end if
        if (c_associated(settings_c%context)) then
            write(stderr, *) "test_assign_arh_c_f failed: Host context associated."
            test_assign_arh_c_f = .false.
        end if

        ! check against reference values
        if (settings_c /= ref_arh_settings) then
            write(stderr, *) "test_assign_arh_c_f failed: Settings not converted "// &
                "correctly."
            test_assign_arh_c_f = .false.
        end if

        ! check initialization flag
        if (.not. settings_c%initialized) then
            test_assign_arh_c_f = .false.
            write(stderr, *) "test_assign_arh_c_f failed: Settings not marked as "// &
                "initialized."
        end if

    end function test_assign_arh_c_f

end module otr_arh_c_interface_unit_tests
