module energies
    use ieee_arithmetic, only: ieee_positive_inf, ieee_value
    use types, only: dp
    use parameters
    implicit none
    save
    public energy, denergy, pseudohs, hardsphere
contains
    ! This configuration calculates the energy of a given configuration
    subroutine energy(x, y, z, ener)
        real(dp), intent(in) :: x(:), y(:), z(:)
        real(dp), intent(out) :: ener

        ! Local variables
        integer :: i, j
        real(dp) :: rij, xij, yij, zij, uij
        
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
                    ! call pseudohs(rij, uij)
                    call hardsphere(rij, uij)
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
        real(dp) :: rij, xij, yij, zij, uij

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
                call hardsphere(rij, uij)
                ! call pseudohs(rij, uij)
                dener = dener + uij
            end if
        end do
    end subroutine denergy

    ! This configuration calculates the pair potential between particles i & j
    subroutine pseudohs(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        if (rij < bpot) then
            uij = (a2/dT)*((1._dp/rij)**dlr-(1._dp/rij)**dla)
            uij = uij + 1._dp/dT
        else
            uij = 0.0_dp
        end if

    end subroutine pseudohs

    subroutine hardsphere(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij
        real(dp) :: rinf

        if (rij < 1.0_dp) then
            uij = ieee_value(rinf, ieee_positive_inf)
        else
            uij = 0.0_dp
        end if

    end subroutine hardsphere
end module energies