program main
    use types, only: dp
    use parameters
    use utils
    use energies
    use movement
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none

    ! Local scalar variables
    real(dp) :: d, volratio, del, rng, ener
    real(dp) :: rhoaverage, rhoave, current_volume, rhoprom
    integer :: rngint, i, j, k, nattemp, nacc
    integer :: thermsteps, eqsteps, u, vacc, vattemp
    integer :: v, avevolfreq, accsize
    ! Local arrays
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp), allocatable :: rhoacc(:), volacc(:), volsqacc(:)

    ! Initialize the RNG
    call random_seed()

    ! Read an input file that contains all the necessary information
    call parse_input('input.in', eqsteps, thermsteps, avevolfreq)
    
    ! Update the simulation parameters with this information
    boxl = (np / rho)**(1.0_dp/3.0_dp)
    rc = boxl / 2.0_dp
    d = (1.0_dp / rho)**(1.0_dp/3.0_dp)

    ! Initialization of variables
    del = 2.0_dp
    nattemp = 0
    nacc = 1
    rhoaverage = 0.0_dp
    rhoave = 0.0_dp
    rhoprom = 0.0_dp
    vacc = 1
    vattemp = 0
    j = 0
    current_volume = 0.0_dp

    print*, 'rc = ', rc
    print*, 'Mean interparticle distance: ', d
    print*, 'Pressure = ', pressure
    print*, 'Reference density = ', rho

    ! Allocate memory for arrays
    allocate(x(np), y(np), z(np))
    accsize = eqsteps / avevolfreq
    allocate(rhoacc(accsize), volacc(accsize), volsqacc(accsize))
    rhoacc = 0.0_dp
    volacc = 0.0_dp
    volsqacc = 0.0_dp

    ! Either read a configuration file or generate a new one
    if (from_file) then
        write(unit=output_unit, fmt='(a)') 'Reading from positions file...'
        open(newunit=u, file = 'configuration.dat', status = 'unknown')
            do i = 1, np
                read(u, *) x(i), y(i), z(i)
            end do
        close(u)
    else
        ! initial configuration as a simple lattice
        call iniconfig(x, y, z, d)
    end if

    ! Initial configuration energy, regardless of how it was created
    call energy(x, y, z, ener)

    print*, 'E/N for the initial configuration:', ener/np

    ! MC cycle to thermalize the system
    do i = 1, thermsteps
        do k = 1, np
            call random_number(rng)
            rngint = 1 + floor((np + 1) * rng)

            if (rngint <= np) then
                call mcmove(x, y, z, ener, nattemp, nacc, del)
                call adjust(nattemp, nacc, del, 0.35_dp)
            else
                call mcvolume(x, y, z, rhoave, ener, vattemp, vacc)
                call adjust(vattemp, vacc, dispvol, 0.25_dp)
            end if
        end do
            
        if (mod(i, 10000) == 0) then
            write(unit=output_unit, fmt='(a)') 'MC Step, Particle disp, Energy / N, Disp ratio'
            print*, i, del, ener/real(np, dp), real(nacc, dp) / real(nattemp, dp)
            write(unit=output_unit, fmt='(a)') 'MC Step, Density average, box size, Vol ratio, Vol disp'
            volratio = real(vacc, dp) / real(vattemp, dp)
            print*, i, rhoave / vacc, boxl, volratio, dispvol
        end if
    end do

    ! Reset accumulation variables
    nattemp = 0
    nacc = 1
    vattemp = 0
    vacc = 1
    ! Start accumulating results
    write(unit=output_unit, fmt='(a)') 'Averaging starts...'
    ! Open the necessary files for saving thermodynamical quantities
    open(newunit=u, file='energy.dat', status='unknown')
    open(newunit=v, file='density.dat', status='unknown')
    ! Production cyle
    do i = 1, eqsteps
        do k = 1, np
            call random_number(rng)
            rngint = 1 + floor((np + 1) * rng)

            if (rngint <= np) then
                call mcmove(x, y, z, ener, nattemp, nacc, del)
                ! call adjust(nattemp, nacc, del, 0.3_dp)
            else
                call mcvolume(x, y, z, rhoave, ener, vattemp, vacc)
                ! call adjust(vattemp, vacc, dispvol, 0.2_dp)
            end if
        end do
        
        if (mod(i, avevolfreq) == 0) then
            ! Save the value for the energy
            write(unit=u, fmt='(i12.2, f15.10)') i, ener/real(np, dp)
            
            ! Update the accumulation index
            j = j + 1

            ! Accumulate the results for the density
            current_volume = real(np, dp) / rho
            rhoaverage = rhoaverage + rho
            rhoprom = rhoaverage / real(j, dp)
            rhoacc(j) = rho
            volacc(j) = current_volume
            volsqacc(j) = current_volume**2.0_dp

            ! Save all results to file
            write(unit=v, fmt='(2f17.10)') rhoprom, current_volume
        end if
    end do

    ! Close off files that were opened for saving information
    close(u)
    close(v)

    ! Do some averaging for the density
    call calc_variable(rhoacc, 'Density', 'average_density.dat')
    ! Do some averaging for the volume
    call calc_variable(volacc, 'Volume', 'average_volume.dat')
    ! Do some averaging for the squared volume
    call calc_variable(volsqacc, 'Squared Volume', 'average_sqvolume.dat')

    ! write the final configuration to file
    open(newunit=u, file = 'configuration.dat', status = 'unknown')
    do i = 1, np
        write(unit=u, fmt='(3f14.8)') x(i), y(i), z(i)
    end do
    close(u)
end program main