module parameters
    use types, only: dp
    implicit none
    save
! CONSTANT PARAMETERS
    

    ! constant values
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: diam = 1.0_dp

    ! User dependant parameters
    real(dp) :: phi, rho, boxl, rc, ktemp, pressure, dispvol
    integer :: np
end module parameters