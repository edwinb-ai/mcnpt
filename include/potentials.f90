module potentials

use ieee_arithmetic, only: ieee_positive_inf, ieee_value
use types

! Intrinsic pseudohs potential arguments that need to be computed
! only once
real(dp), parameter :: dlr = 50.0_dp, dT = 1.4737_dp
real(dp), parameter :: dla = 49.0_dp
real(dp), parameter :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))
real(dp), parameter :: bpot = (dlr/dla)**(1.0_dp/(dlr-dla))

implicit none

public pseudohs, hardsphere, lennardjones, hertzian, &
softsphere, yukawa_attr, gaussian

contains
    subroutine pseudohs(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        if (rij < bpot) then
            uij = (a2/dT)*((1.0_dp/rij)**dlr-(1.0_dp/rij)**dla)
            uij = uij + 1.0_dp/dT
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

    subroutine lennardjones(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        real(dp) :: temp

        temp = 1.0_dp / (rij**6)

        uij = 4.0_dp * (temp**2 - temp)
        
    end subroutine lennardjones

    subroutine softsphere(rij, uij, n)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij
        integer, intent(in) :: n

        if (rij <= 1.0) then
            uij = (1.0_dp / rij)**n
        else
            uij = 0.0_dp
        end if

    end subroutine softsphere

    subroutine yukawa_attr(rij, uij, k)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij, k

        real(dp) :: rinf

        if (rij < 1.0) then
            uij = ieee_value(rinf, ieee_positive_inf)
        else
            uij = -(1.0_dp / rij) * exp(-k * (rij - 1.0_dp))
        end if

    end subroutine yukawa_attr

    subroutine gaussian(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        uij = exp(-0.5_dp * rij**2)
    end subroutine gaussian

    subroutine hertzian(rij, uij)
        real(dp), intent(inout) :: uij
        real(dp), intent(in) :: rij

        if (rij < 1.0_dp) then
            uij = (1.0_dp - rij)**(5.0_dp / 2.0_dp)
        else
            uij = 0.0_dp
        end if

    end subroutine hertzian

end module potentials