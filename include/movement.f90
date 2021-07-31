module movement
    use types, only: dp
    use parameters
    use energies, only: denergy
    use observables, only: rdf, sq

    implicit none

    public mcmove, adjust, average
contains
    ! This subroutine displace the system to a new configuration
    subroutine mcmove(x, y, z, ener, nattemp, nacc, del)
    real(dp), intent(in) :: del
    real(dp), intent(inout) :: ener
    integer, intent(inout) :: nattemp, nacc
    real(dp), intent(inout) :: x(:), y(:), z(:)
    ! it is because are evaluate

    ! Local variables
    integer :: no
    real(dp) :: xo, yo, zo, enero, enern, dener
    real(dp) :: rng

    nattemp = nattemp + 1

    call random_number(rng)
    no = int(rng*np)+1
    call denergy(x, y, z, no, enero) !intruduzco el valor de energia de conf inicial

    xo = x(no) !coordenadas de de la particula aleatoriamente escogida
    yo = y(no)
    zo = z(no)

    call random_number(rng)
    x(no) = x(no)+(rng-0.5_dp)*del
    call random_number(rng)
    y(no) = y(no)+(rng-0.5_dp)*del
    call random_number(rng)
    z(no) = z(no)+(rng-0.5_dp)*del

    ! periodic boundary conditions
    x(no) = x(no)-boxl*dnint(x(no)/boxl)
    y(no) = y(no)-boxl*dnint(y(no)/boxl)
    z(no) = z(no)-boxl*dnint(z(no)/boxl)

    call denergy(x, y, z, no, enern)

    dener = enern-enero
    call random_number(rng)
    if (rng < exp(-dener / ktemp)) then
        ener = ener+dener
        nacc = nacc+1
    else
        x(no) = xo
        y(no) = yo
        z(no) = zo
    end if
    end subroutine mcmove

    subroutine average(x, y, z, g, s, ener, nattemp, nacc, ng, naveg, del, dr, pbc)
    real(dp), intent(inout) :: x(:), y(:), z(:)
    real(dp), intent(inout) :: g(:), s(:)

    real(dp), intent(in) :: del, dr
    integer, intent(inout) :: nattemp, nacc, naveg, ng
    real(dp), intent(inout) :: ener
    integer, intent(in) :: pbc

    ! Local variables
    integer :: no
    real(dp) :: rng
    real(dp) :: xo, yo, zo, enero, enern, dener
    
    nattemp = nattemp + 1

    call random_number(rng)
    no = int(rng * np) + 1

    xo = x(no)
    yo = y(no)
    zo = z(no)

    call denergy(x, y, z, no, enero)

    call random_number(rng)
    x(no) = x(no)+(rng-0.5_dp)*del
    call random_number(rng)
    y(no) = y(no)+(rng-0.5_dp)*del
    call random_number(rng)
    z(no) = z(no)+(rng-0.5_dp)*del

    !periodic boundary conditions
    x(no) = x(no)-boxl*dnint(x(no)/boxl)
    y(no) = y(no)-boxl*dnint(y(no)/boxl)
    z(no) = z(no)-boxl*dnint(z(no)/boxl)

    call denergy(x, y, z, no, enern)

    dener = enern-enero

    call random_number(rng)
    if (rng < exp(-dener / ktemp)) then
        ener = ener + dener
        nacc = nacc + 1
        ng = ng + 1
        ! Calcular observables
        if (mod(ng, np) == 0) then
            naveg = naveg + 1
            call rdf(x, y, z, dr, g, pbc)
            call sq(x, y, z, s, pbc)
        end if
    else
        x(no) = xo
        y(no) = yo
        z(no) = zo
    end if
    end subroutine average

    subroutine mcvolume(x, y, z, ener, nattemp, nacc, del)
    real(dp), intent(in) :: del
    real(dp), intent(inout) :: ener
    integer, intent(inout) :: nattemp, nacc
    real(dp), intent(inout) :: x(:), y(:), z(:)

    ! Local variables
    integer :: no, i
    real(dp) :: xo, yo, zo, enero, enern, dener
    real(dp) :: rng, volold, volnew, lnvolold, lnvolnew
    real(dp) :: adjust, dispvol, boxlnew, rhold, denpt

    ! Count as a movement always
    nattemp = nattemp + 1

    ! Estimate the new volume
    dispvol = 0.001
    volold = boxl**3
    lnvolold = log(volold)
    call random_number(rng)
    lnvolnew = lnvolold + dispvol * (rng - 0.5_dp)
    boxlnew = exp(lnvolnew)**(1.0_dp / 3.0_dp)

    ! Adjust the particles to the new box
    adjust = boxlnew / boxl
    boxl = boxlnew
    do i = 1, np
        x(i) = x(i) * adjust
        y(i) = y(i) * adjust
        z(i) = z(i) * adjust
    end do

    ! Compute the new density
    rhold = rho
    rho = np / volnew

    call random_number(rng)
    no = int(rng*np) + 1
    call denergy(x, y, z, no, enero)
    xo = x(no)
    yo = y(no)
    zo = z(no)

    ! periodic boundary conditions
    x(no) = x(no)-boxl*dnint(x(no)/boxl)
    y(no) = y(no)-boxl*dnint(y(no)/boxl)
    z(no) = z(no)-boxl*dnint(z(no)/boxl)

    ! Compute the full energy
    call denergy(x, y, z, no, enern)
    dener = enern - enero
    denpt = pressure * (volnew - volold) + dener
    denpt = denpt + (np + 1) * (lnvolnew - lnvolold)

    ! Apply Metropolis criteria
    call random_number(rng)
    if (rng <= exp(-denpt / ktemp)) then
        ener = ener + dener
        nacc = nacc + 1
    else
        x(no) = xo
        y(no) = yo
        z(no) = zo

        do i = 1, np
            x(i) = x(i) / adjust
            y(i) = y(i) / adjust
            z(i) = z(i) / adjust
        end do

        boxl = volold**(1/3)
        rho = rhold
    end if
    end subroutine mcvolume

    ! This subroutine adjusts the displacement of particles
    subroutine adjust(nattemp, nacc, del)
        integer, intent(in) :: nattemp, nacc
        real(dp), intent(inout) :: del
        ! Local variables
        real(dp) :: ratio

        if (mod(nattemp, nacc) == 0) then
            ratio = real(nacc, dp)/real(nattemp, dp)
            if (ratio > 0.5_dp) then
                del = del*1.05_dp
            else
                del = del*0.95_dp
            end if
        end if
    end subroutine adjust
end module movement