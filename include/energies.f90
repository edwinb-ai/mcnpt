module energies
    use types, only: dp
    use parameters
    implicit none
    save
    public energy, denergy, potential
contains
    ! This configuration calculates the energy of a given configuration
    subroutine energy(x, y, z, ener)
        real(dp), intent(in) :: x(:), y(:), z(:)
        real(dp), intent(out) :: ener

        ! Local variables
        integer :: i, j
        real(dp) :: rij, xij, yij, zij, dB, uij
        
        ener = 0._dp

        do i = 1, np - 1
            do j = i + 1, np

                xij = x(j)-x(i)
                yij = y(j)-y(i)
                zij = z(j)-z(i)

                ! Minimum image convention
                xij = xij-boxl*nint(xij/boxl)
                yij = yij-boxl*nint(yij/boxl)
                zij = zij-boxl*nint(zij/boxl)

                rij = norm2([xij, yij, zij])

                if (rij < rc) then
                    if (rij < bpot) then
                        call potential(rij, uij)
                    else
                        uij = 0._dp
                    end if
                    ener = ener + uij
                end if
            end do
        end do
    end subroutine energy

    ! This subroutine calculates the difference in energy when a particle is displaced
    subroutine denergy(x, y, z, no, dener)
        real(dp), intent(in) :: x(:), y(:), z(:)
        real(dp), intent(out) :: dener
        integer, intent(in) :: no
        ! Local variables
        integer :: i
        real(dp) :: rij, xij, yij, zij, dB, uij

        dener = 0._dp ! initializing
        do i = 1, np
            if ( i == no ) cycle

            xij = x(no)-x(i)
            yij = y(no)-y(i)
            zij = z(no)-z(i)
            
            ! Minimum image convention
            xij = xij-boxl*dnint(xij/boxl)
            yij = yij-boxl*dnint(yij/boxl)
            zij = zij-boxl*dnint(zij/boxl)
            
            rij = norm2([xij, yij, zij])

            if (rij < rc) then
                if (rij < bpot) then
                    call potential(rij, uij)
                else
                    uij = 0._dp
                end if
                dener = dener + uij
            end if
        end do
    end subroutine denergy

    ! This configuration calculates the pair potential between particles i & j
    subroutine potential(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij
        ! Local
        real(dp) :: dA

        uij = (a2/dT)*((1._dp/rij)**dlr-(1._dp/rij)**dla)
       uij = uij + 1._dp/dT

    end subroutine potential
end module energies