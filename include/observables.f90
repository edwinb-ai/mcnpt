module observables
    use parameters
    implicit none
    save
    public rdf, sq
contains
    ! This subroutine calculates the pair potential between particles i & j
    subroutine rdf(x, y, z, dr, g)
    real(dp), intent(in) :: x(:), y(:), z(:)
    real(dp), intent(inout) :: g(:)

    real(dp), intent(in) :: dr

    ! Local variables
    integer :: i, j, nbin
    real(dp) :: rij, xij, yij, zij

    do i = 1, np - 1
        do j = i + 1, np
            xij = x(j)-x(i)
            yij = y(j)-y(i)
            zij = z(j)-z(i)
            
            xij = xij-boxl*dnint(xij/boxl)
            yij = yij-boxl*dnint(yij/boxl)
            zij = zij-boxl*dnint(zij/boxl)
            rij = norm2([xij, yij, zij])

            nbin = int(rij/dr) + 1
            if (nbin <= mr) then
                g(nbin) = g(nbin)+2._dp
            end if
        end do
    end do
    end subroutine rdf ! out g

    subroutine sq(x, y, z, s)
    real(dp), intent(in) :: x(:), y(:), z(;)
    real(dp), intent(inout) :: s(:)

    ! Local variables
    real(dp) :: auxc(:), auxs(:)
    integer :: i, kq, k, j
    real(dp) :: xaux, yaux, zaux, rij, arg, sum, parti, auxsq

    do i = 2, mr ! why i star in 2
        ! creating two zero vectors with leng nvq
        do kq = 1, nvq
            auxc(kq) = 0._dp
            auxs(kq) = 0._dp
        end do

        parti = 0._dp
        do k = 1, np
            xaux = x(k)-boxl*nint(x(k)/boxl)
            yaux = y(k)-boxl*nint(y(k)/boxl)
            zaux = z(k)-boxl*nint(z(k)/boxl)
            rij = norm2([xaux, yaux, zaux])

            if (rij < rc) then
                parti = parti+1._dp

                do j = 1, nvq
                arg = qx(i, j)*xaux+qy(i, j)*yaux+qz(i, j)*zaux
                auxc(j) = auxc(j)+dcos(arg)
                auxs(j) = auxs(j)+dsin(arg)
                end do
            end if
        end do

        sum = 0._dp
        do kq = 1, nvq
            sum = sum+auxc(kq)*auxc(kq)+auxs(kq)*auxs(kq)
        end do

        auxsq = sum/(nvq*parti)
        s(i) = s(i)+auxsq
    end do
    end subroutine sq ! out s
end module observables