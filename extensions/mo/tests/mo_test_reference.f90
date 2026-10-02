! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_mo_test_reference

    use opentrustregion, only: ip
    use otr_common_test_reference, only: n_particle, n_occ

    implicit none

    ! number of MOs for orbitals parameterized in the MO basis, fewer than the shared
    ! number of AOs, and the numbers of parameters for the shared occupations
    integer(ip), parameter :: n_mo = 4_ip
    integer(ip), parameter :: n_param_cs = n_occ(1) * (n_mo - n_occ(1)), &
                              n_param_os = sum(n_occ * (n_mo - n_occ))

    ! occupation cases for routines acting on every particle channel: closed-shell,
    ! open-shell, and open-shell with an empty occupied or virtual block
    integer(ip), parameter :: n_cases = 4_ip
    integer(ip), parameter :: case_n_particle(n_cases) = &
        [1_ip, n_particle, n_particle, n_particle]
    integer(ip), parameter :: case_n_occ(n_particle, n_cases) = &
        reshape([n_occ(1), 0_ip, n_occ(1), n_occ(2), n_occ(1), 0_ip, n_mo, n_occ(2)], &
                [n_particle, n_cases])
    character(len=33), parameter :: case_names(n_cases) = &
        [character(len=33) :: "closed-shell", "open-shell", &
         "open-shell empty occupied channel", "open-shell empty virtual channel"]

end module otr_mo_test_reference
