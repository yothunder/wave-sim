module gridmod
    implicit none

    integer :: nx, ny, nsigma, ntheta, nt
    real(8) :: dx, dy, dt, tmax
    real(8), parameter :: pi = 4.0d0 * atan(1.0d0)
    real(8), allocatable :: x(:), y(:), sigma(:), theta(:)
    real(8) :: dtheta
    real(8), allocatable :: dsigma(:)

contains

    subroutine grid_init()
        integer :: i
        nx = 30
        ny = 30
        nsigma = 30
        ntheta = 36
        dx = 1.0d0
        dy = 1.0d0
        dt = 10.0d0
        tmax = 3600.0d0 * 1.0d0
        nt = int(tmax / dt)

        allocate(x(nx), y(ny), sigma(nsigma), theta(ntheta))
        x = [(i * dx, i = 0, nx - 1)]
        y = [(i * dy, i = 0, ny - 1)]

        call grid_sigma_init()
        dtheta = 2.0d0 * pi / ntheta
        theta = [( (i-1) * dtheta, i = 1, ntheta)]
    end subroutine grid_init

    subroutine grid_sigma_init()
        integer :: i
        real(8) :: sigma_min, sigma_max, dlsigma

        sigma_min = 0.040d0
        sigma_max = 1.0d0
        dlsigma   = log(sigma_max / sigma_min) / (nsigma - 1)
        allocate(dsigma(nsigma))

        do i = 1, nsigma
            sigma(i) = sigma_min * exp(dlsigma * (i - 1))
        end do

        do i = 1, nsigma - 1
            dsigma(i) = sigma(i+1) - sigma(i)
        end do
        dsigma(nsigma) = dsigma(nsigma - 1)

    end subroutine grid_sigma_init

end module gridmod