program main
    use types, only: dp
    use parameters
    use utils
    use energies
    use movement
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none

    ! Local variables, note that somes variables are initialized
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: del = 0.1_dp, ener
    real(dp) :: d, rhoave, volratio
    real(dp) :: rng
    real(dp) :: rhoaverage, rhosq, rhoprom, rhodev
    integer :: nattemp = 0
    integer :: rngint
    integer :: nacc = 1, i, j
    integer :: thermsteps, eqsteps, u, nptvolfreq, vacc, vattemp
    integer :: v, avevolfreq
    real(dp), allocatable :: rhoacc(:)

    ! Initialize the RNG
    call random_seed()

    ! Read an input file that contains all the necessary information
    call parse_input('input.in', eqsteps)
    
    ! Update the simulation parameters with this information
    boxl = (np / rho)**(1.0_dp/3.0_dp)
    rc = boxl * 0.5_dp
    d = (1.0_dp/rho)**(1.0_dp/3.0_dp)
    nptvolfreq = np * 2
    avevolfreq = 1000
    rhoave = 0.0_dp
    vacc = 1
    vattemp = 0
    j = 1
    rhoprom = 0.0_dp
    rhosq = 0.0_dp
    thermsteps = 2e8

    print*, 'rc = ', rc
    print*, 'Mean interparticle distance: ', d
    print*, 'Pressure = ', pressure
    print*, 'Reference density = ', rho

    ! Allocate memory for arrays
    allocate(x(np), y(np), z(np), rhoacc(eqsteps / avevolfreq))
    rhoacc = 0.0_dp

    if (from_file) then
        write(unit=output_unit, fmt='(a)') 'Reading from positions file...'
        open(newunit=u, file = 'configuration.dat', status = 'unknown')
            do i = 1, np
                read(u, *) x(i), y(i), z(i)
            end do
        close(u)
    end if

    ! initial configuration and initial energy
    call iniconfig(x, y, z, d)
    call energy(x, y, z, ener)

    print*, 'E/N for the initial configuration:', ener/np

    ! MC cycle to thermalize the system
    open(newunit=u, file = 'energy.dat', status = 'unknown')
    open(newunit=v, file = 'density.dat', status = 'unknown')
    
    do i = 1, thermsteps
        call random_number(rng)
        rngint = floor((np + 1) * rng)

        if (rngint < np) then
            call mcmove(x, y, z, ener, nattemp, nacc, del)
            call adjust(nattemp, nacc, del, 0.5_dp)
        else
            call mcvolume(x, y, z, rhoave, ener, vattemp, vacc)
            call adjust(vattemp, vacc, dispvol, 0.4_dp)
        end if
        
        if (mod(i, 100) == 0) then
            write(u, '(2f15.7)') real(i, dp), ener / real(np, dp)
        end if
        
        if (mod(i, 1000000) == 0) then
            write(unit=output_unit, fmt='(a)') 'MC Step, Particle disp, Energy / N'
            print*, i, del, ener/real(np, dp), real(nacc, dp) / real(nattemp, dp)
            write(unit=output_unit, fmt='(a)') 'MC Step, Density average, box size, Vol ratio, Vol disp'
            volratio = real(vacc, dp) / real(vattemp, dp)
            print*, i, rhoave / vacc, boxl, volratio, dispvol
        end if
    end do

    ! Start accumulating results
    write(unit=output_unit, fmt='(a)') 'Averaging starts...'
    do i = 1, eqsteps
        call random_number(rng)
        rngint = floor((np + 1) * rng)

        if (rngint < np) then
            call mcmove(x, y, z, ener, nattemp, nacc, del)
            call adjust(nattemp, nacc, del, 0.5_dp)
        else
            call mcvolume(x, y, z, rhoave, ener, vattemp, vacc)
            call adjust(vattemp, vacc, dispvol, 0.4_dp)
        end if
        
        if (mod(i, avevolfreq) == 0) then
            write(unit=u, fmt='(2f15.12)') real(i, dp), ener / real(np, dp)
        end if

        if (mod(i, avevolfreq) == 0) then
            volratio = real(vacc, dp) / real(vattemp, dp)
            rhoaverage = rhoave / vacc
            rhoacc(j) = rhoaverage
            rhoprom = rhoprom + rhoaverage
            rhosq = rhosq + rhoaverage**2
            write(unit=v, fmt='(f12.10)') rhoaverage
            j = j + 1
        end if
    end do

    close(u)
    close(v)
    call block_average(rhoacc)
    ! Do some averaging
    rhoprom = rhoprom / real(j, dp)
    rhosq = rhosq / real(j, dp)
    rhodev = sqrt(rhosq - rhoprom**2)
    write(unit=output_unit, fmt='(a)') 'Average, std deviation'
    write(unit=output_unit, fmt='(2f15.7)') rhoprom, rhodev
    ! write the final configuration and the energy
    open(newunit=u, file = 'configuration.dat', status = 'unknown')
    do i = 1, np
        write(unit=u, fmt='(3f15.7)') x(i), y(i), z(i)
    end do
    close(u)

    deallocate(x, y, z, rhoacc)

end program main