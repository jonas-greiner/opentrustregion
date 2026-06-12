! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_arh_unit_tests

    use opentrustregion, only: rp, ip, stderr
    use test_reference, only: tol
    use, intrinsic :: iso_c_binding, only: c_bool

    implicit none

contains

    function generate_random_density_matrix(n_ao, n_electrons) result(dm)
        !
        ! this function generates a random valid density matrix by first generating a 
        ! random set of orthonormal basis vectors and then summing the outer products 
        ! of these
        !
        use opentrustregion, only: numerical_zero

        integer(ip), intent(in) :: n_ao, n_electrons
        real(rp) :: dm(n_ao, n_ao)

        real(rp) :: vec(n_ao), val, U(n_ao, n_electrons)
        integer(ip) :: i, j

        ! generate random orthonormal basis
        do j = 1, n_electrons
            do
                call random_number(vec)
                do i = 1, j - 1
                    vec = vec - sum(vec * U(:, i)) * U(:, i)
                end do
                val = sqrt(sum(vec**2))
                if (val >= numerical_zero) then
                    U(:, j) = vec / val
                    exit
                end if
            end do
        end do

        ! construct valid density matrix by summing outer products of basis vectors
        dm = matmul(U, transpose(U))

    end function generate_random_density_matrix

    logical(c_bool) function test_get_response_contribution_closed_shell() bind(C)
        !
        ! this function tests the function that computes the response contribution to 
        ! the ARH Hessian for the closed-shell case
        !
        use otr_arh, only: get_response_contribution_closed_shell, arh_settings_type
        use opentrustregion, only: numerical_zero

        integer(ip), parameter :: n_ao = 6, n_electrons = 3, n_particle = 1, &
                                  n_diff = 2, lwork = 4 * n_diff
        real(rp) :: vec(n_ao), val, U(n_ao, n_electrons), &
                    dm_oao(n_ao, n_ao, n_particle), &
                    dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    x(n_ao, n_ao, n_particle), y(n_ao, n_ao, n_particle), &
                    hess(n_ao, n_ao, n_ao, n_ao), &
                    fock_diff(n_ao, n_ao, n_particle, n_diff), &
                    metric_eigvecs(n_diff, n_diff, 1), metric_eigvals(n_diff, 1), &
                    g_x(n_ao, n_ao, n_particle), g_y(n_ao, n_ao, n_particle), &
                    coeff(n_diff), dm_target(n_ao, n_ao, n_particle), &
                    fock_target(n_ao, n_ao, n_particle), &
                    expected_response(n_ao, n_ao, n_particle), proj_v(n_ao, n_ao), &
                    response(n_ao, n_ao, n_particle), work(lwork)
        integer(ip) :: i, j, k, l, m, info
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_response_contribution_closed_shell = .true.

        ! generate random density matrix
        dm_oao(:, :, 1) = generate_random_density_matrix(n_ao, n_electrons)

        ! generate valid density matrix differences
        do i = 1, n_diff
            call random_number(x(:, :, 1))
            x(:, :, 1) = x(:, :, 1) - transpose(x(:, :, 1))
            dm_diff(:, :, 1, i) = matmul(dm_oao(:, :, 1), x(:, :, 1))
            dm_diff(:, :, 1, i) = dm_diff(:, :, 1, i) + transpose(dm_diff(:, :, 1, i))
        end do

        ! create a mock Hessian operator (assuming array is only read for 
        ! first index >= second index)
        do l = 1, n_ao
            do k = 1, l
                do j = 1, l
                    do i = 1, merge(k, j, j == l)
                        call random_number(val)
                        hess(i, j, k, l) = val
                        hess(i, j, l, k) = val
                        hess(k, l, i, j) = val
                        hess(k, l, j, i) = val
                    end do
                end do
            end do
        end do

        ! generate corresponding Fock matrix differences
        do k = 1, n_diff
            do j = 1, n_ao
                do i = 1, j
                    fock_diff(i, j, 1, k) = sum(hess(i, j, :, :) * dm_diff(:, :, 1, k))
                    fock_diff(j, i, 1, k) = fock_diff(i, j, 1, k)
                end do
            end do
        end do

        ! compute and diagonalize metric
        metric_eigvecs = 0.0_rp
        do j = 1, n_diff
            do i = 1, n_diff
                metric_eigvecs(i, j, 1) = sum(dm_diff(:, :, 1, i) * dm_diff(:, :, 1, j))
            end do
        end do
        call dsyev("V", "U", 2_ip, metric_eigvecs(:, :, 1), 2_ip, &
                   metric_eigvals(:, 1), work, lwork, info)

        ! generate two random antisymmetric trial matrices
        call random_number(x)
        call random_number(y)
        x(:, :, 1) = x(:, :, 1) - transpose(x(:, :, 1))
        y(:, :, 1) = y(:, :, 1) - transpose(y(:, :, 1))

        ! compute the Hessian actions G(x) and G(y) and test symmetry of Hessian 
        ! operator for symmetrized ARH method
        settings%arh_type = "symmetric"
        g_x = get_response_contribution_closed_shell(dm_oao, x, dm_diff, fock_diff, &
                                                     metric_eigvals, metric_eigvecs, &
                                                     n_ao, n_particle, settings)
        g_y = get_response_contribution_closed_shell(dm_oao, y, dm_diff, fock_diff, &
                                                     metric_eigvals, metric_eigvecs, &
                                                     n_ao, n_particle, settings)
        if ((sum(y * g_x) - sum(x * g_y)) > tol) then
            write (stderr, *) "test_get_response_contribution_closed_shell failed: "// &
                "Hessian operator is not symmetric for symmetrized ARH method."
            test_get_response_contribution_closed_shell = .false.
        end if

        ! compute the Hessian actions G(x) and G(y) and test symmetry of Hessian 
        ! operator for multisecant PSB method
        settings%arh_type = "multisecant_psb"
        g_x = get_response_contribution_closed_shell(dm_oao, x, dm_diff, fock_diff, &
                                                     metric_eigvals, metric_eigvecs, &
                                                     n_ao, n_particle, settings)
        g_y = get_response_contribution_closed_shell(dm_oao, y, dm_diff, fock_diff, &
                                                     metric_eigvals, metric_eigvecs, &
                                                     n_ao, n_particle, settings)
        if ((sum(y * g_x) - sum(x * g_y)) > tol) then
            write (stderr, *) "test_get_response_contribution_closed_shell failed: "// &
                "Hessian operator is not symmetric for multisecant_psb ARH method."
            test_get_response_contribution_closed_shell = .false.
        end if

        ! create a linear combination of previous directions
        dm_target = 0.0_rp
        fock_target = 0.0_rp
        call random_number(coeff)
        do i = 1, n_diff
            dm_target(:, :, 1) = dm_target(:, :, 1) + coeff(i) * dm_diff(:, :, 1, i)
            fock_target(:, :, 1) = fock_target(:, :, 1) + coeff(i) * &
                                   fock_diff(:, :, 1, i)
        end do

        ! construct corresponding trial vector
        x(:, :, 1) = matmul(dm_oao(:, :, 1), dm_target(:, :, 1))
        x(:, :, 1) = x(:, :, 1) - transpose(x(:, :, 1))

        ! project out redundant components of expected response contribution
        proj_v = 0.0_rp
        do i = 1, n_ao
            proj_v(i, i) = 1.0_rp
        end do
        proj_v = proj_v - dm_oao(:, :, 1)
        expected_response(:, :, 1) = matmul(dm_oao(:, :, 1), &
                                            matmul(fock_target(:, :, 1), proj_v))
        expected_response(:, :, 1) = expected_response(:, :, 1) - &
                                     transpose(expected_response(:, :, 1))

        ! check if multisecant conditions are fulfilled for standard ARH method
        settings%arh_type = "standard"
        response = &
            get_response_contribution_closed_shell(dm_oao, x, dm_diff, fock_diff, &
                                                   metric_eigvals, metric_eigvecs, &
                                                   n_ao, n_particle, settings)
        if (norm2(response - expected_response) > tol) then
            write (stderr, *) "test_get_response_contribution_closed_shell failed: "// &
                "Hessian operator does not fulfill multisecant conditions for "// &
                "standard ARH method."
            test_get_response_contribution_closed_shell = .false.
        end if

        ! check if multisecant conditions are fulfilled for multisecant PSB method
        settings%arh_type = "multisecant_psb"
        response = &
            get_response_contribution_closed_shell(dm_oao, x, dm_diff, fock_diff, &
                                                   metric_eigvals, metric_eigvecs, &
                                                   n_ao, n_particle, settings)
        if (norm2(response - expected_response) > tol) then
            write (stderr, *) "test_get_response_contribution_closed_shell failed: "// &
                "Hessian operator does not fulfill multisecant conditions for "// &
                "multisecant_psb ARH method."
            test_get_response_contribution_closed_shell = .false.
        end if

    end function test_get_response_contribution_closed_shell


    logical(c_bool) function test_get_response_contribution_open_shell() bind(C)
        !
        ! this function tests the function that computes the response contribution to 
        ! the ARH Hessian for the open-shell case
        !
        use otr_arh, only: get_response_contribution_open_shell, arh_settings_type
        use opentrustregion, only: numerical_zero

        integer(ip), parameter :: n_ao = 6, n_electrons = 3, n_particle = 2, &
                                  n_diff = 2, lwork = 4 * n_diff
        real(rp) :: vec(n_ao), val, U(n_ao, n_electrons), &
                    dm_oao(n_ao, n_ao, n_particle), &
                    dm_diff(n_ao, n_ao, n_particle, n_diff), &
                    x(n_ao, n_ao, n_particle), y(n_ao, n_ao, n_particle), &
                    same_hess(n_ao, n_ao, n_ao, n_ao), &
                    opposite_hess(n_ao, n_ao, n_ao, n_ao), &
                    same_v_diff(n_ao, n_ao, n_particle, n_diff), &
                    opposite_v_diff(n_ao, n_ao, n_particle, n_diff), &
                    metric_eigvecs(n_diff, n_diff, n_particle), &
                    metric_eigvals(n_diff, n_particle), &
                    g_x(n_ao, n_ao, n_particle), g_y(n_ao, n_ao, n_particle), &
                    coeff(n_diff), dm_target(n_ao, n_ao, n_particle), &
                    v_target(n_ao, n_ao, n_particle), &
                    expected_response(n_ao, n_ao, n_particle), proj_v(n_ao, n_ao), &
                    response(n_ao, n_ao, n_particle), work(lwork)
        integer(ip) :: i, j, k, l, m, n, info
        type(arh_settings_type) :: settings

        ! assume tests pass
        test_get_response_contribution_open_shell = .true.

        do n = 1, n_particle
            ! generate random density matrix
            dm_oao(:, :, n) = generate_random_density_matrix(n_ao, n_electrons)

            ! generate valid density matrix differences
            do i = 1, n_diff
                call random_number(x(:, :, n))
                x(:, :, n) = x(:, :, n) - transpose(x(:, :, n))
                dm_diff(:, :, n, i) = matmul(dm_oao(:, :, n), x(:, :, n))
                dm_diff(:, :, n, i) = dm_diff(:, :, n, i) + &
                                      transpose(dm_diff(:, :, n, i))
            end do
        end do

        ! create mock Hessian operators (assuming arrays are only read for 
        ! first index >= second index)
        do l = 1, n_ao
            do k = 1, l
                do j = 1, l
                    do i = 1, merge(k, j, j == l)
                        call random_number(val)
                        same_hess(i, j, k, l) = val
                        same_hess(i, j, l, k) = val
                        same_hess(k, l, i, j) = val
                        same_hess(k, l, j, i) = val
                        call random_number(val)
                        opposite_hess(i, j, k, l) = val
                        opposite_hess(i, j, l, k) = val
                        opposite_hess(k, l, i, j) = val
                        opposite_hess(k, l, j, i) = val
                    end do
                end do
            end do
        end do

        do n = 1, n_particle
            ! generate corresponding spin-resolved potential matrix differences
            do k = 1, n_diff
                do j = 1, n_ao
                    do i = 1, j
                        same_v_diff(i, j, n, k) = sum(same_hess(i, j, :, :) * &
                                                      dm_diff(:, :, n, k))
                        same_v_diff(j, i, n, k) = same_v_diff(i, j, n, k)
                        opposite_v_diff(i, j, n, k) = sum(opposite_hess(i, j, :, :) * &
                                                          dm_diff(:, :, n, k))
                        opposite_v_diff(j, i, n, k) = opposite_v_diff(i, j, n, k)
                    end do
                end do
            end do

            ! compute and diagonalize metric
            metric_eigvecs(:, :, n) = 0.0_rp
            do j = 1, n_diff
                do i = 1, n_diff
                    metric_eigvecs(i, j, n) = sum(dm_diff(:, :, n, i) * &
                                                  dm_diff(:, :, n, j))
                end do
            end do
            call dsyev("V", "U", n_diff, metric_eigvecs(:, :, n), n_diff, &
                       metric_eigvals(:, n), work, lwork, info)
        end do

        ! generate two random antisymmetric trial matrices
        call random_number(x)
        call random_number(y)
        do n = 1, n_particle
            x(:, :, n) = x(:, :, n) - transpose(x(:, :, n))
            y(:, :, n) = y(:, :, n) - transpose(y(:, :, n))
        end do

        ! compute the Hessian actions G(x) and G(y) and test symmetry of Hessian 
        ! operator for symmetrized ARH method
        settings%arh_type = "symmetric"
        g_x = get_response_contribution_open_shell(dm_oao, x, dm_diff, same_v_diff, &
                                                   opposite_v_diff, metric_eigvals, &
                                                   metric_eigvecs, n_ao, n_particle, &
                                                   settings)
        g_y = get_response_contribution_open_shell(dm_oao, y, dm_diff, same_v_diff, &
                                                   opposite_v_diff, metric_eigvals, &
                                                   metric_eigvecs, n_ao, n_particle, &
                                                   settings)
        if ((sum(y * g_x) - sum(x * g_y)) > tol) then
            write (stderr, *) "test_get_response_contribution_open_shell failed: "// &
                "Hessian operator is not symmetric for symmetrized ARH method."
            test_get_response_contribution_open_shell = .false.
        end if

        ! compute the Hessian actions G(x) and G(y) and test symmetry of Hessian 
        ! operator for multisecant PSB method
        settings%arh_type = "multisecant_psb"
        g_x = get_response_contribution_open_shell(dm_oao, x, dm_diff, same_v_diff, &
                                                 opposite_v_diff, metric_eigvals, &
                                                 metric_eigvecs, n_ao, n_particle, &
                                                 settings)
        g_y = get_response_contribution_open_shell(dm_oao, y, dm_diff, same_v_diff, &
                                                 opposite_v_diff, metric_eigvals, &
                                                 metric_eigvecs, n_ao, n_particle, &
                                                 settings)
        if ((sum(y * g_x) - sum(x * g_y)) > tol) then
            write (stderr, *) "test_get_response_contribution_open_shell failed: "// &
                "Hessian operator is not symmetric for multisecant_PSB ARH method."
            test_get_response_contribution_open_shell = .false.
        end if

        dm_target = 0.0_rp
        v_target = 0.0_rp
        call random_number(coeff)
        do n = 1, n_particle
            ! create a linear combination of previous directions
            do j = 1, n_diff
                dm_target(:, :, n) = dm_target(:, :, n) + coeff(j) * dm_diff(:, :, n, j)
                v_target(:, :, n) = v_target(:, :, n) + coeff(j) * &
                                    same_v_diff(:, :, n, j)
                do i = 1, n_particle
                    if (i == n) cycle
                    v_target(:, :, n) = v_target(:, :, n) + coeff(j) * &
                                        opposite_v_diff(:, :, i, j)
                end do
            end do

            ! construct corresponding trial vector
            x(:, :, n) = matmul(dm_oao(:, :, n), dm_target(:, :, n))
            x(:, :, n) = x(:, :, n) - transpose(x(:, :, n))

            ! project out redundant components of expected response contribution
            proj_v = 0.0_rp
            do i = 1, n_ao
                proj_v(i, i) = 1.0_rp
            end do
            proj_v = proj_v - dm_oao(:, :, n)
            expected_response(:, :, n) = matmul(dm_oao(:, :, n), &
                                                matmul(v_target(:, :, n), proj_v))
            expected_response(:, :, n) = expected_response(:, :, n) - &
                                         transpose(expected_response(:, :, n))
        end do

        ! check if multisecant conditions are fulfilled for standard ARH method
        settings%arh_type = "standard"
        response = &
            get_response_contribution_open_shell(dm_oao, x, dm_diff, same_v_diff, &
                                                 opposite_v_diff, metric_eigvals, &
                                                 metric_eigvecs, n_ao, n_particle, &
                                                 settings)
        if (norm2(response - expected_response) > tol) then
            write (stderr, *) "test_get_response_contribution_open_shell failed: "// &
                "Hessian operator does not fulfill multisecant conditions for "// &
                "standard ARH method."
            test_get_response_contribution_open_shell = .false.
        end if

        ! check if multisecant conditions are fulfilled for multisecant PSB method
        settings%arh_type = "multisecant_psb"
        response = &
            get_response_contribution_open_shell(dm_oao, x, dm_diff, same_v_diff, &
                                                 opposite_v_diff, metric_eigvals, &
                                                 metric_eigvecs, n_ao, n_particle, &
                                                 settings)
        if (norm2(response - expected_response) > tol) then
            write (stderr, *) "test_get_response_contribution_open_shell failed: "// &
                "Hessian operator does not fulfill multisecant conditions for "// &
                "multisecant_psb ARH method."
            test_get_response_contribution_open_shell = .false.
        end if

    end function test_get_response_contribution_open_shell

end module otr_arh_unit_tests
