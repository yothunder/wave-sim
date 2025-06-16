module iomod
    use gridmod
    use windmod
    implicit none

contains

    subroutine read_wind(filename, wind_arr)
        character(len=*), intent(in) :: filename
        real(8), intent(out) :: wind_arr(:,:,:)
        integer :: i, j, it, ios
        open(unit=10, file=filename, status='old', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening wind file: ', filename
            stop
        end if

        do it = 1, nt
            do j = 1, ny
                do i = 1, nx
                read(10,*) wind_arr(it, i, j)
                end do
            end do
        end do

        close(10)
    end subroutine read_wind

    ! subroutine read_boundary_file(N_bdy, filename)
    !     use gridmod
    !     implicit none
    !     real(8), intent(out) :: N_bdy(:,:,:,:,:)
    !     character(len=*), intent(in) :: filename
    !     integer :: it, i, j, k, l, ios
    !     open(unit=12, file=filename, status='old', iostat=ios)
    !     if (ios /= 0) then
    !         print *, 'Error opening boundary file: ', filename
    !         stop
    !     end if

    !     do it = 1, nt
    !         do l = 1, ntheta
    !             do k = 1, nsigma
    !                 do j = 1, ny
    !                     do i = 1, nx
    !                         read(12,*) N_bdy(it,i,j,k,l)
    !                     end do
    !                 end do
    !             end do
    !         end do
    !     end do

    !     close(12)
    ! end subroutine read_boundary_file

    ! subroutine read_boundary_file(N_bdy, filename)
    !     use gridmod
    !     implicit none
    !     real(8), intent(out) :: N_bdy(:,:,:,:,:)
    !     character(len=*), intent(in) :: filename
    !     integer :: ios

    !     open(unit=12, file=filename, status='old', iostat=ios)
    !     if (ios /= 0) then
    !         print *, 'Error opening boundary file: ', filename
    !         stop
    !     end if

    !     read(12, *) N_bdy

    !     close(12)
    ! end subroutine read_boundary_file

    subroutine read_boundary_file(N_bdy_t, filename, it)
        use gridmod
        implicit none
        real(8), intent(out) :: N_bdy_t(nx, ny, nsigma, ntheta)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: it
        integer :: i, j, k, l, n, ios
        real(8) :: tmp
        integer :: skip

        skip = (it-1)*nx*ny*nsigma*ntheta

        open(unit=12, file=filename, status='old', iostat=ios)
        if (ios /= 0) then
            print *, 'Error opening boundary file: ', filename
            stop
        end if

        ! Skip to the desired time step
        do n = 1, skip
            read(12,*) tmp
        end do

        ! Read the current time step
        do l = 1, ntheta
            do k = 1, nsigma
                do j = 1, ny
                    do i = 1, nx
                        read(12,*) N_bdy_t(i,j,k,l)
                    end do
                end do
            end do
        end do

        close(12)
    end subroutine read_boundary_file

    ! subroutine apply_boundary(N, N_bdy, it)
    !     use gridmod
    !     integer, intent(in) :: it
    !     real(8), intent(inout) :: N(:,:,:,:)
    !     real(8), intent(in) :: N_bdy(:,:,:,:,:)

    !     integer :: i, j, k, l

    !     do l = 1, ntheta
    !         do k = 1, nsigma
    !             do j = 1, ny
    !                 ! i=1 and i=nx boundaries
    !                 N(1, j, k, l) = N_bdy(it, 1, j, k, l)
    !                 N(nx, j, k, l) = N_bdy(it, nx, j, k, l)
    !             end do
    !             do i = 1, nx
    !                 ! j=1 and j=ny boundaries
    !                 N(i, 1, k, l) = N_bdy(it, i, 1, k, l)
    !                 N(i, ny, k, l) = N_bdy(it, i, ny, k, l)
    !             end do
    !         end do
    !     end do
    ! end subroutine apply_boundary

    subroutine apply_boundary(N, N_bdy, it)
        use gridmod
        real(8), intent(inout) :: N(:,:,:,:)
        real(8), intent(in) :: N_bdy(nx, ny, nsigma, ntheta)
        integer, intent(in) :: it
        integer :: i, j, k, l

        do l = 1, ntheta
            do k = 1, nsigma
                do j = 1, ny
                    N(1, j, k, l) = N_bdy(1, j, k, l)
                    N(nx, j, k, l) = N_bdy(nx, j, k, l)
                end do
                do i = 1, nx
                    N(i, 1, k, l) = N_bdy(i, 1, k, l)
                    N(i, ny, k, l) = N_bdy(i, ny, k, l)
                end do
            end do
        end do
    end subroutine apply_boundary

    subroutine write_diagnostics(Hs, Tm, n)
        use gridmod
        character(len=500) :: fname
        integer :: i, j
        real(8), intent(in) :: Hs(nx, ny), Tm(nx, ny)
        integer, intent(in) :: n

        write(fname, '("/home/yothunder/fort/wvmod/output/diag/diagnostics_", I4.4, ".txt")') n
        open(unit=20, file=fname, status='replace')
        do j = 1, ny
            do i = 1, nx
                write(20, '(F15.4, 1X, F15.4)') Hs(i,j), Tm(i,j)
            end do
        end do
        close(20)
    end subroutine write_diagnostics

    subroutine write_spectrum(N, it, i, j)
        use gridmod, only: nsigma, ntheta, dtheta, sigma
        implicit none
        integer, intent(in) :: it, i, j
        real(8), intent(in) :: N(:,:,:,:)
        
        real(8) :: E_sigma(nsigma)
        character(len=100) :: filename
        integer :: k, l, unit

        E_sigma = 0.0d0
        do k = 1, nsigma
            do l = 1, ntheta
                E_sigma(k) = E_sigma(k) + N(i, j, k, l)
            end do
            E_sigma(k) = E_sigma(k) * dtheta  ! approximate integral over theta
        end do

        write(filename, '("/home/yothunder/fort/wvmod/output/spc/spectrum_", I4.4, ".txt")') it
        open(newunit=unit, file=filename, status="replace")
        do k = 1, nsigma
            write(unit, '(1PE12.4, 1x, 1PE12.4)') sigma(k), E_sigma(k)
        end do
        close(unit)
    end subroutine write_spectrum


end module iomod
