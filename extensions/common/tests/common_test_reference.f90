! Copyright (C) 2025- Jonas Greiner
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

module otr_common_test_reference

    use opentrustregion, only: ip
    use c_interface, only: c_ip

    implicit none

    ! numbers of particle channels, number of AOs and number of the occupied orbitals
    ! of every channel, which differ between the channels
    integer(ip), parameter :: n_particle = 2_ip, n_ao = 5_ip, &
                              n_occ(n_particle) = [2_ip, 1_ip]
    integer(c_ip), protected, bind(C, name="test_n_ao") :: n_ao_c = int(n_ao, kind=c_ip)

contains

    function request_label(outputs, request) result(label)
        !
        ! this function describes which optional outputs of a density matrix evaluating
        ! function a test requests besides the energy, given as the bits of an integer
        ! in the order of the outputs, for the failure messages of tests which go
        ! through every combination of them
        !
        character(*), intent(in) :: outputs(:)
        integer(ip), intent(in) :: request
        character(:), allocatable :: label

        integer(ip) :: i, n_requested, n_listed

        ! list the energy and every requested output
        n_requested = count([(btest(request, i - 1), i = 1, size(outputs))])
        n_listed = 0
        label = " when requesting the energy"
        do i = 1, size(outputs)
            if (.not. btest(request, i - 1)) cycle
            n_listed = n_listed + 1
            if (n_listed == n_requested) then
                label = label//" and the "//trim(outputs(i))
            else
                label = label//", the "//trim(outputs(i))
            end if
        end do

    end function request_label

    function capitalized(text) result(capital)
        !
        ! this function returns a text with its first letter capitalized, so that the
        ! name of an output can start a failure message
        !
        character(*), intent(in) :: text
        character(len(text)) :: capital

        capital = text
        if (lge(text(1:1), "a") .and. lle(text(1:1), "z")) &
            capital(1:1) = achar(iachar(text(1:1)) - iachar("a") + iachar("A"))

    end function capitalized

end module otr_common_test_reference
