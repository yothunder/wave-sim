program wavemod
    use gridmod
    use windmod
    use wavepropmod
    use physicsmod
    use iomod
    implicit none
    
    real(8), allocatable :: N(:,:,:,:), N_new(:,:,:,:)
    real(8) :: t
    integer :: it, jsea, ispec, i, j, k, l

    ! Allocate core fields
    allocate(N(nx, ny, nsigma, ntheta))
    allocate(N_new(nx, ny, nsigma, ntheta))
    allocate(u10(nt, nx, ny))

    ! Initialize grid and spectra
    call grid_init()

    ! Read wind field
    call read_wind('/home/yothunder/fort/wvmod/input/u10.txt')

    ! Read initial boundary spectrum
    call read_boundary_N(N, '/home/yothunder/fort/wvmod/input/bound.txt')

    N = 0.0d0
    N_new = 0.0d0

    do it = 1, nt
        t = it * dt

        do j = 1, ny
            do i = 1, nx
                jsea = (j - 1) * nx + i
                do l = 1, ntheta
                    do k = 1, nsigma
                        ispec = (l-1) * nsigma + k
                        ! Compute spatial propagation (e.g., via upwind scheme)
                        call propagate_wave(N, N_new, i, j, k, l)
                        ! Apply source terms (e.g., wind input and dissipation)
                        call calc_source_terms(N, N_new, it, i, j, k, l)
                    end do
                end do
            end do
        end do
        
        N = N_new

        if (mod(it, 24) == 0) then
            write(*, '(A,I5,A,F10.2,A)') 'Step: ', it, ', Time (hr): ', t/3600, new_line('a')
        end if
    end do
    
    write(*, *) 'Simulation completed successfully.'

end program wavemod