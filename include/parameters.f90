module parameters
    use types, only: dp
    implicit none
    save
! CONSTANT PARAMETERS
    ! pot arguments
    real(dp), parameter :: dlr = 50.0_dp, dT = 1.4737_dp
    real(dp), parameter :: dla = 49.0_dp
    real(dp), parameter :: a2 = (dlr/(dlr-dla))*(dlr/dla)**(dla/(dlr-dla))
    real(dp), parameter :: bpot = (dlr/dla)**(1.0_dp/(dlr-dla))

    ! parameters arguments
    real(dp), parameter :: phi = 0.45_dp
    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp), parameter :: rho = 6.0_dp*phi/pi
    real(dp), parameter :: diam = 1.0_dp

    ! box arguments
    integer, parameter :: np = 6**3! number of particles

    real(dp), parameter :: boxl = (np/rho)**(1.0_dp/3.0_dp)
    real(dp), parameter :: rc = boxl/2.0_dp

    ! mp and mr arguments
    integer, parameter :: mp = 1024, mr = 2**9
    integer, parameter :: nvq = 20
    real(dp), allocatable :: qx(:, :), qy(:, :), qz(:, :)
end module parameters