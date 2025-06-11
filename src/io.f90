module iomod
    use gridmod
    use windmod
    implicit none

contains

    subroutine read_wind(filename)
        character(len=*), intent(in) :: filename
        integer :: i, j, n, ios
        open(unit=10, file=filename, status='old', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening wind file: ', filename
            stop
        end if

        do n = 1, nt
            do j = 1, ny
                do i = 1, nx
                read(10,*) u10(n, i, j)
                end do
            end do
        end do

        close(10)
    end subroutine read_wind

    subroutine read_boundary_N(N, filename)
        real(8), intent(out) :: N(:,:,:,:)
        character(len=*), intent(in) :: filename
        integer :: i, j, k, l, ios
        open(unit=11, file=filename, status='old', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening boundary file: ', filename
            stop
        end if

        do l = 1, ntheta
            do k = 1, nsigma
                do j = 1, ny
                    do i = 1, nx
                        read(11,*) N(i,j,k,l)
                    end do
                end do
            end do
        end do
        close(11)
    end subroutine read_boundary_N

end module iomod
