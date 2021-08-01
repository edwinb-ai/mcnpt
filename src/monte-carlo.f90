program main
    use types, only: dp
    use parameters
    use utils
    use energies
    use movement

    implicit none

    ! Local variables, note that somes variables are initialized
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: del = 0.1_dp, ener, dv, dt2, dphi, sumq, qs
    real(dp) :: rng, d, dr, dq, rhoave, volratio
    real(dp) :: rhoaverage, rhosq, rhoprom, rhomean, rhodev
    integer :: nattemp = 0
    integer :: nacc = 1, nacco, nav, i, j, ncq = 0
    integer :: ng = 0, naveg = 0
    integer :: limT, u, nptvol, nptvolfreq, vacc, vattemp
    integer :: vacco, v, avevolfreq
    real(dp), allocatable :: rhoacc(:)
    ! Condiciones periódicas a la frontera
    integer :: pbc = 1

    ! Inicializar el RNG
    call random_seed()
    ! Read an input file that contains all the necessary information
    call parse_input('input.in', limT)
    ! Update the simulation parameters with this information
    rho = 6.0_dp * real(phi) / pi
    boxl = (np / rho)**(1.0_dp/3.0_dp)
    rc = boxl * 0.5_dp
    d = (1.0_dp/rho)**(1.0_dp/3.0_dp)
    nptvolfreq = np * 2
    avevolfreq = 25000
    rhoave = 0.0_dp
    nptvol = 1
    vacc = 1
    vattemp = 0
    j = 1
    rhoprom = 0.0_dp
    rhosq = 0.0_dp

    print*, 'rc = ', rc
    print*, 'dr = ', dr
    print*, 'Mean interparticle distance: ', d
    print*, 'Pressure = ', pressure
    print*, 'Reference density = ', rho

    ! Allocate memory for arrays
    allocate(x(np), y(np), z(np), rhoacc((limT / 2) / avevolfreq))
    rhoacc = 0.0_dp

    ! initial configuration and initial energy
    call iniconfig(x, y, z, d)
    call energy(x, y, z, ener)

    print*, 'E/N for the initial configuration:', ener/np

    ! MC cycle to thermalize the system
    open(newunit=u, file = 'energy.dat', status = 'unknown')
    open(newunit=v, file = 'density.dat', status = 'unknown')
    do i = 1, limT
        call mcmove(x, y, z, ener, nattemp, nacc, del)
        call adjust(nattemp, nacc, del, 0.5_dp)

        ! Adjust the box if asked for
        if ((nptvol == 1) .and. (mod(i, nptvolfreq)) == 0) then
            call mcvolume(x, y, z, rhoave, ener, vattemp, vacc)
            call adjust(vattemp, vacc, dispvol, 0.2_dp)
        end if
        
        if (mod(i, 100) == 0) then
            write(u, '(2f15.7)') i*1._dp, ener/np
        end if
        
        if (mod(i, 500000) == 0) then
            print*, 'MC Step, Particle disp, Energy / N'
            print*, i, del, ener/np
            print*, 'MC Step, Density average, box size, Vol ratio, Vol disp'
            volratio = real(vacc, dp) / real(vattemp, dp)
            print*, i, rhoave / vacc, boxl, volratio, dispvol
        end if

        ! Start accumulating results
        if (i > limT / 2) then
            if (mod(i, avevolfreq) == 0) then
                print*, 'Accumulating results...'
                print*, 'MC Step, Density average, box size, Vol ratio, Vol disp'
                volratio = real(vacc, dp) / real(vattemp, dp)
                rhoaverage = rhoave / vacc
                rhoacc(j) = rhoaverage
                rhoprom = rhoprom + rhoaverage
                rhosq = rhosq + rhoaverage**2
                print*, i, rhoaverage, boxl, volratio, dispvol
                write(v, *) i, rhoaverage
                j = j + 1
            end if
        end if
    end do

    ! print*, 'The system has thermalized'
    close(u)
    close(v)
    call block_average(rhoacc)
    ! Do some averaging
    rhoprom = rhoprom / real(j, dp)
    rhosq = rhosq / real(j, dp)
    rhodev = sqrt(rhosq - rhoprom**2)
    print*, 'Average, std deviation'
    print*, rhoprom, rhodev
    ! write the final configuration and the energy
    open(newunit=u, file = 'finalconf.dat', status = 'unknown')
    do i = 1, np
        write(u, '(3f15.7)') x(i), y(i), z(i)
    end do
    close(u)

    deallocate(x, y, z, rhoacc)

end program main