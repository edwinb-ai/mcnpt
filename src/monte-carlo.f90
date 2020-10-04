program main
    use parameters
    use utils
    use energy
    use movement

    implicit none

    ! Local variables, note that somes variables was initialized
    real(dp), :: x(np), y(np), z(np) ! three vector of mp dimension
    real(dp), :: r(mr), g(mr), q(mr), s(mr) ! four vector of mr dimension
    real(dp), parameter :: dr = rc/mr, dq = pi/rc
    real(dp) :: del = 0.1_dp, ener , dv, dt2, dphi, sumq, qs ! variables used
    integer :: nattemp = 0
    integer :: nacc = 1, nacco, nav, i, j, ncq = 0
    integer :: ng = 0, naveg = 0
    integer, parameter :: limT = 100000, limG = 100000

    print*, 'rc = ',rc, 'dr = ', dr , 'dq = ', dq, 'boxl =', boxl
    print*
    print*, 'Mean interparticle distance: ', rho**(-1./3.)

    ! Gives values to q vector
    do i=1, mr
        q(i) = (i-1)*dq
    enddo

    ! Gives values to qx, qy, qz matrices and saving them
    open(121, file = 'qvectors_N.dat', status = 'unknown')
    do i = 1, mr
        do j = 1, nvq
            dt2 = pi*ranf(iseed)
            dphi = 2*pi*ranf(iseed)

            qx(i,j) = q(i)*dcos(dphi)*dsin(dt)
            qy(i,j) = q(i)*dsin(dphi)*dsin(dt)
            qz(i,j) = q(i)*dcos(dt)

            ncq = ncq+1
            qs = dsqrt(qx(i,j)*qx(i,j)+qy(i,j)*qy(i,j)+qz(i,j)*qz(i,j))

            write(121, '(3f15.7)') dfloat(ncq), qs
        end do
    end do
    close(121)

    !initial configuration
    call iniconfig(x, y, z) ! this subroutine create the data

    open(10, file = 'iniconf_N.dat', status = 'unknown')
    do i = 1, np
        write(10, '(3f15.7)') x(i), y(i), z(i)
    end do
    close(10)

    call energy(x, y, z, ener)

    print*, 'Energy per particle of the initial configuration:', ener/np

! MC cycle to thermalize the system
    open(30, file = 'energy_N.dat', status = 'unknown')
    do i = 1, limT
        call mcmove(x, y, z, ener, nattemp, nacc, del, iseed)
        call adjust(nattemp, nacc, del)
        ! mod(n,m) gives the remainder when n is divided by m
        if (mod(i, 100) == 0) write(30, '(3f15.7)') i*1._dp, ener/np
        if (mod(i, 10000) == 0) then
            print*, i, del, ener/np
            !pause
        end if
    end do

    print*, 'The system has thermalized'
 ! write the final configuration and the energy
    open(20, file = 'finalconf_N.dat', status = 'unknown')
    do i = 1, np
        write(20, '(3f15.7)') x(i), y(i), z(i) !guardo mi foto final
    end do

    close(20)
    close(30)

    !stop 'final  de primeras subrutinas'

    !MC cycle to calculate the g(r)
    nacco = nacc
    g(:) = 0.0_dp

    print*, iseed
    do i = 1, limG
        call average(x, y, z, g, s, ener, nattemp, nacc, ng, naveg, del, dr, iseed)
        call adjust(nattemp, nacc, del)
        if (mod(i, 10000) == 0) print*, i, 'calculating g(r) and S(q)'
    end do

    nav = nacc-nacco

    print*,'Average number for energy: ', nav ! cual es la diferencia nav y naveg
    print*,'Average number for g(r): ', naveg

!    stop 'final de 2da subrutinas'

! This is the radial distribution function
    open(50, file = 'gr_o_N.dat', status = 'unknown')
    do i = 2, mr
        r(i) = (i-1)*dr
        dv = (4._dp*pi*r(i)**2._dp*dr)*rho
        g(i) = g(i) / (np*naveg*dv)
        write(50, '(2f15.8)') r(i), g(i)
    end do
    close(50)

! This is the structure factor from the definition
    open(51, file = 'sq_o_N.dat', status = 'unknown')
    do i = 3, mr ! why i start in 3
        s(i) = s(i)/naveg
        if (q(i) < 30._dp) then
            write(51, '(3f15.7)') q(i), s(i)
        end if
    end do
    close(51)

end program main