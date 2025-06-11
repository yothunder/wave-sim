module gridmod
    implicit none

    integer :: nx, ny, nsigma, ntheta, nt
    real(8) :: dx, dy, dt, tmax
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real(8), allocatable :: x(:), y(:), sigma(:), theta(:)

contains

    subroutine grid_init()
        integer :: i
        nx = 10
        ny = 10
        nsigma = 30
        ntheta = 36
        dx = 10.0d0
        dy = 10.0d0
        dt = 300
        tmax = 3600.0d0 * 6.0d0
        nt = int(tmax / dt)

        allocate(x(nx), y(ny), sigma(nsigma), theta(ntheta))
        x = [(i * dx, i = 0, nx - 1)]
        y = [(i * dy, i = 0, ny - 1)]

        call grid_sigma_init()
        theta = [( (i-1) * 2.0d0 * pi / ntheta, i = 1, ntheta)]
    end subroutine grid_init

    subroutine grid_sigma_init()
        integer :: i
        real(8) :: sigma_min, sigma_max, dlsigma

        sigma_min = 0.040d0
        sigma_max = 1.0d0
        dlsigma   = log(sigma_max / sigma_min) / (nsigma - 1)

        do i = 1, nsigma
            sigma(i) = sigma_min * exp(dlsigma * (i - 1))
        end do
    end subroutine grid_sigma_init

end module gridmod