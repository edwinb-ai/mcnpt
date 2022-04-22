module parameters
    use types, only: dp
    implicit none
    save
! CONSTANT PARAMETERS
    

    ! constant values
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: diam = 1.0_dp

    ! User dependant parameters
    real(dp) :: rho, boxl, rc, ktemp, pressure, del, dispvol
    integer :: np
    logical :: from_file
end module parameters