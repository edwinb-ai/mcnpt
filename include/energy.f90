module energy
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
        
        ! Inicializar variables
        ener = 0._dp

        do i = 1, np - 1
            do j = i + 1, np
                ! Distancias
                xij = x(j)-x(i)
                yij = y(j)-y(i)
                zij = z(j)-z(i)
                ! Condiciones periódicas a la frontera  
                xij = xij-boxl*dnint(xij/boxl)
                yij = yij-boxl*dnint(yij/boxl)
                zij = zij-boxl*dnint(zij/boxl)
                ! Distancia Euclidiana
                rij = norm2([xij, yij, zij])
                ! Aplicar el valor del potencial
                if (rij <. rc) then
                    dB = dl / (dl-1._dp)
                    if (rij <. dB) then
                        call potential(rij, uij)uij
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
            ! Si las partículas son iguales, que se salte la iteración
            if ( i == no ) cycle
            ! Distancias
            xij = x(no)-x(i)
            yij = y(no)-y(i)
            zij = z(no)-z(i)
            ! Condiciones periódicas a la frontera
            xij = xij-boxl*dnint(xij/boxl)
            yij = yij-boxl*dnint(yij/boxl)
            zij = zij-boxl*dnint(zij/boxl)
            ! Distancia Euclidiana
            rij = norm2([xij, yij, zij])
            ! Aplicar el valor del potencial
            if (rij < rc) then
                dB = dl/(dl-1._dp)
                if (rij < dB) then
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

        dA = dl*(dl/(dl-1._dp))**(dl-1._dp)
        uij = (dA/dT)*((1._dp/rij)**dl-(1._dp/rij)**(dl-1._dp))
        uij = uij+1._dp/dT

    end subroutine potential
end module energy