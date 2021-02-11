module utils
    use types, only: dp
    use parameters, only: np, rc
    implicit none
    save
    public iniconfig
contains
    subroutine iniconfig(x, y, z, d)
    ! defining three vector of mp dimension, it indicate that only are out variables
        real(dp), intent(inout) :: x(:), y(:), z(:)
        real(dp), intent(in) :: d

        ! Local variables
        integer :: i ! it is neccesary for restart i
        x(1) = -rc + d/2._dp
        y(1) = -rc + d/2._dp
        z(1) = -rc + d/2._dp

        do i = 1, np-1
            x(i+1) = x(i)+d
            y(i+1) = y(i)
            z(i+1) = z(i)
            if (x(i+1) > rc) then
                x(i+1) = -rc+d/2._dp
                y(i+1) = y(i+1)+d
                z(i+1) = z(i)
                if (y(i+1) > rc) then
                    x(i+1) = -rc+d/2._dp
                    y(i+1) = -rc+d/2._dp
                    z(i+1) = z(i+1)+d
                end if
            end if
        end do
    end subroutine iniconfig ! note that the out is x, y, z vector
end module utils