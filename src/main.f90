program wavemod
    use gridmod
    use windmod
    use wavepropmod
    use physicsmod
    use iomod
    use paramsmod
    implicit none
    
    real(8), allocatable :: N(:,:,:,:), N_new(:,:,:,:)
    real(8), allocatable :: Hs(:,:), Tm(:,:)
    real(8) :: t
    real(8), parameter :: N_seed = 1.0d-6
    integer :: it, jsea, ispec, i, j, k, l

    call grid_init()

    allocate(N(nx, ny, nsigma, ntheta))
    allocate(N_new(nx, ny, nsigma, ntheta))
    allocate(N_bdy(nx, ny, nsigma, ntheta))
    allocate(u10(nt,nx,ny), v10(nt,nx,ny))
    allocate(Hs(nx, ny))
    allocate(Tm(nx, ny))

    print *, 'Starting wave model simulation...'
    print *, 'Grid size: ', nx, 'x', ny, 'x', nsigma, 'x', ntheta
    print *, 'Total time steps: ', nt
    print *, 'Time step (s): ', dt
    print *, 'Total simulation time (hr): ', tmax / 3600.0d0

    print *, 'Opening wind forcing...'
    call read_wind('/home/yothunder/fort/wvmod/input/u10.txt', u10)
    call read_wind('/home/yothunder/fort/wvmod/input/v10.txt', v10)

    N(:,:,:,:) = N_seed ! Initialize with a small background spectrum
    do it = 1, nt
        t = it * dt
        !$omp parallel do collapse(2) private(i, j, k, l, jsea, ispec)
        do j = 1, ny
            do i = 1, nx
                jsea = (j - 1) * nx + i
                do l = 1, ntheta
                    do k = 1, nsigma
                        ispec = (l-1) * nsigma + k
                        call propagate_wave(N, N_new, i, j, k, l, it)      ! Compute spatial propagation (via upwind scheme)
                        call calc_source_terms(N, N_new, it, i, j, k, l)   ! Apply source terms (wind input and dissipation)
                        N_new(i,j,k,l) = max(N_new(i,j,k,l), 0.0d0)
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        N = N_new
        call compute_params(N_new, Hs, Tm)
        call write_diagnostics(Hs, Tm, it)
        call write_spectrum(N, it, nx/2, ny/2)
        
        if (mod(it, 10) == 0) then
            write(*, '(A,I5,A,F10.2,A)') 'Step: ', it, ', Time (hr): ', t/3600, new_line('a')
        end if
    end do
    
    write(*, *) 'Simulation completed successfully.'

end program wavemod