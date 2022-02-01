module utils
    use types
    use parameters
    use, intrinsic :: iso_fortran_env, only: output_unit

    implicit none

    public iniconfig, parse_input, calc_variable

contains
    subroutine iniconfig(x, y, z, d)
    ! defining three vector of mp dimension, it indicate that only are out variables
        real(dp), intent(inout) :: x(:), y(:), z(:)
        real(dp), intent(in) :: d

        ! Local variables
        integer :: i ! it is neccesary for restart i
        real(dp) :: half_box

        half_box = boxl / 2.0_dp

        x(1) = -half_box + (d / 2.0_dp)
        y(1) = -half_box + (d / 2.0_dp)
        z(1) = -half_box + (d / 2.0_dp)

        do i = 1, np-1
            x(i+1) = x(i)+d
            y(i+1) = y(i)
            z(i+1) = z(i)

            if (x(i+1) > half_box) then
                x(i+1) = -half_box + (d / 2.0_dp)
                y(i+1) = y(i+1)+d
                z(i+1) = z(i)

                if (y(i+1) > half_box) then
                    x(i+1) = -half_box + (d / 2.0_dp)
                    y(i+1) = -half_box + (d / 2.0_dp)
                    z(i+1) = z(i+1)+d
                end if
            end if
        end do
    end subroutine iniconfig ! note that the out is x, y, z vector

    subroutine parse_input(filein, limg, limterm, avefreq)
        character(len = *), intent(in) :: filein
        integer, intent(inout) :: limg, limterm, avefreq

        ! Local variables
        integer :: u

        open(newunit=u, file=filein, status='old')
        ! Read in the local variables
        read(u, *) rho
        read(u, *) ktemp
        read(u, *) pressure
        read(u, *) dispvol
        read(u, *) np
        read(u, *) limg
        read(u, *) limterm
        read(u, *) avefreq
        read(u, *) from_file
        close(u)
    end subroutine parse_input

    subroutine block_average(x, filein)
        real(dp), intent(in) :: x(:)
        character(len=*), intent(in) :: filein

        integer :: u
        integer :: nstep, tblock, nblock, blk, stp1, stp2, trun
        real(dp) :: a_blk, a_run, a_var, a_var_1, a_err, si, a_avg
        real(dp), allocatable :: xcentre(:)

        nstep = size(x)
        allocate(xcentre(nstep))
        a_avg   = sum(x) / real(nstep, dp)      ! Sample average
        xcentre = x - a_avg                    ! Centre the data
        a_var_1 = sum(xcentre**2) / real(nstep-1, dp) ! Bias-corrected sample variance
        print*, a_avg, sqrt(a_var_1)

        ! Open file to save information
        open(newunit=u, file=filein, status='unknown')

        do nblock = 20, 4, -1 ! Loop over number, and hence length, of blocks

            tblock = nstep / nblock              ! Block length in steps (rounded down)
            trun   = nblock*tblock               ! Run length in steps, accounting for rounding
            a_run  = sum(x(1:trun)) / real(trun, dp) ! Average of data
            a_var  = 0.0_dp                       ! Zero mean-square block average accumulator

            do blk = 1, nblock ! Loop over blocks
                stp1  = (blk-1)*tblock+1                            ! Start of block
                stp2  = blk*tblock                                  ! End of block
                a_blk = sum(x(stp1:stp2) - a_run)/real(tblock, dp) ! Block average of deviation
                a_var = a_var + a_blk**2                            ! Mean-square block average
            end do ! End loop over blocks

            a_var = a_var / real(nblock-1, dp)    ! Bias-corrected variance of block averages
            a_err = sqrt(a_var/real(nblock, dp))  ! Estimate of error from block-average variance
            si = tblock * a_var / a_var_1  ! Statistical inefficiency

            write(unit=u, fmt='(2f15.6)' ) a_run, a_err
        end do ! End loop over number, and hence length, of blocks

    ! WRITE ( unit=output_unit, fmt='(a)' ) 'Plateau at large tblock (small nblock)'
    ! WRITE ( unit=output_unit, fmt='(a)' ) 'should agree quite well with exact error estimate'
    ! WRITE ( unit=output_unit, fmt='(a)' ) 'Can plot SI or error**2 against 1/tblock'
    close(u)
    deallocate(xcentre)

    end subroutine block_average

    subroutine calc_variable(avevar, avename, fname)
        character(len=*), intent(in) :: avename, fname
        real(dp), intent(in) :: avevar(:)

        write(unit=output_unit, fmt='(a)') avename
        write(unit=output_unit, fmt='(a)') 'Average, std deviation'
        call block_average(avevar, fname)
    end subroutine calc_variable
end module utils