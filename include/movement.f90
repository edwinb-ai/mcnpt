module movement
    use types
    use parameters
    use energies

    implicit none

    public mcmove, adjust, mcvolume
contains
    ! This subroutine displace the system to a new configuration
    subroutine mcmove(x, y, z, ener, nattemp, nacc, del)
        real(dp), intent(in) :: del
        real(dp), intent(inout) :: ener
        integer, intent(inout) :: nattemp, nacc
        real(dp), intent(inout) :: x(:), y(:), z(:)

        ! Local variables
        integer :: no
        real(dp) :: xo, yo, zo, enero, enern, dener
        real(dp) :: rng

        nattemp = nattemp + 1

        call random_number(rng)
        no = 1 + floor(np * rng)
        call denergy(x, y, z, no, enero)

        xo = x(no)
        yo = y(no)
        zo = z(no)

        call random_number(rng)
        x(no) = x(no) + (rng - 0.5_dp) * del
        call random_number(rng)
        y(no) = y(no) + (rng - 0.5_dp) * del
        call random_number(rng)
        z(no) = z(no) + (rng - 0.5_dp) * del

        ! periodic boundary conditions
        x(no) = x(no) - (boxl * nint(x(no) / boxl))
        y(no) = y(no) - (boxl * nint(y(no) / boxl))
        z(no) = z(no) - (boxl * nint(z(no) / boxl))

        call denergy(x, y, z, no, enern)

        dener = enern - enero
        call random_number(rng)
        if (rng < exp(-dener / ktemp)) then
            ener = ener + dener
            nacc = nacc + 1
        else
            x(no) = xo
            y(no) = yo
            z(no) = zo
        end if
    end subroutine mcmove

    subroutine mcvolume(x, y, z, rhoave, ener, vattemp, vacc)
        real(dp), intent(inout) :: rhoave
        real(dp), intent(in) :: ener
        integer, intent(inout) :: vattemp, vacc
        real(dp), intent(inout) :: x(:), y(:), z(:)

        ! Local variables
        integer :: i
        real(dp) :: enern, dener
        real(dp) :: rng, volold, volnew, lnvolnew
        real(dp) :: voladjust, boxlnew, rhold, denpt

        ! Count as a movement always
        vattemp = vattemp + 1

        ! Estimate the new volume
        volold = boxl**3.0_dp
        call random_number(rng)
        lnvolnew = log(volold) + (dispvol * (rng - 0.5_dp))
        volnew = exp(lnvolnew)
        boxlnew = volnew**(1.0_dp / 3.0_dp)

        ! Adjust the particles to the new box
        voladjust = boxlnew / boxl
        boxl = boxlnew
        rc = boxl / 2.0_dp
        do i = 1, np
            x(i) = x(i) * voladjust
            y(i) = y(i) * voladjust
            z(i) = z(i) * voladjust
        end do

        ! Compute the new density
        rhold = rho
        rho = real(np, dp) / volnew

        ! Compute the energy after adjusting the box
        call energy(x, y, z, enern)
        dener = enern - ener
        ! Compute the full exponential term for the NPT ensemble
        denpt = pressure * (volnew - volold) + dener
        denpt = denpt - real(np + 1, dp) * log(volnew / volold) * ktemp

        ! Apply Metropolis criteria
        call random_number(rng)
        if (rng < exp(-denpt / ktemp)) then
            rhoave = rhoave + rho
            vacc = vacc + 1
        else
            do i = 1, np
                x(i) = x(i) / voladjust
                y(i) = y(i) / voladjust
                z(i) = z(i) / voladjust
            end do

            boxl = volold**(1.0_dp/3.0_dp)
            rho = rhold
            rc = boxl / 2.0_dp
        end if
    end subroutine mcvolume

    ! This subroutine adjusts the displacement of particles
    subroutine adjust(nattemp, nacc, del, tol)
            integer, intent(in) :: nattemp, nacc
            real(dp), intent(in) :: tol
            real(dp), intent(inout) :: del
            ! Local variables
            real(dp) :: ratio

            if (mod(nattemp, nacc) == 0) then
                ratio = real(nacc, dp) / real(nattemp, dp)
                if (ratio > tol) then
                    del = del * 1.05_dp
                else
                    del = del * 0.95_dp
                end if
            end if
    end subroutine adjust
end module movement