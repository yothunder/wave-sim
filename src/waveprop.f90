module wavepropmod
    use gridmod
    use windmod
    implicit none

contains

    subroutine propagate_wave(N, N_new, i, j, k, l, it)
        use gridmod
        implicit none

        integer, intent(in) :: i, j, k, l
        real(8), intent(inout) :: N(:,:,:,:), N_new(:,:,:,:)
        real(8) :: cg, cg_x, cg_y, cg_sigma, cg_theta, c
        real(8) :: adv_x, adv_y, adv_sigma, adv_theta
        real(8) :: N0, Nup
        integer :: l_up, it
        real(8), parameter :: g = 9.81d0

        N0 = N(i, j, k, l)

        ! Group velocity and spatial components
        cg = 0.5d0 * g / sigma(k)
        cg_x = cg * cos(theta(l))
        cg_y = cg * sin(theta(l))

        ! Spatial X advection (upwind)
        if (cg_x >= 0.0d0) then
            if (i > 1) then
                Nup = N(i-1,j,k,l)
            else
                Nup = N0
            end if
        else
            if (i < nx) then
                Nup = N(i+1,j,k,l)
            else
                Nup = N0
            end if
        end if
        adv_x = -cg_x * dt / dx * (N0 - Nup)

        ! Spatial Y advection (upwind)
        if (cg_y >= 0.0d0) then
            if (j > 1) then
                Nup = N(i,j-1,k,l)
            else
                Nup = N0
            end if
        else
            if (j < ny) then
                Nup = N(i,j+1,k,l)
            else
                Nup = N0
            end if
        end if
        adv_y = -cg_y * dt / dy * (N0 - Nup)

        ! Spectral frequency advection (σ)
        c = g / sigma(k)
        cg_sigma = 0.05d0 * max(u10(it,i,j) / c - 0.5d0, 0.0d0)
        if (cg_sigma >= 0.0d0) then
            if (k > 1) then
                Nup = N(i,j,k-1,l)
            else
                Nup = N0
            end if
        else
            if (k < nsigma) then
                Nup = N(i,j,k+1,l)
            else
                Nup = N0
            end if
        end if
        adv_sigma = -cg_sigma * dt / dsigma(k) * (N0 - Nup)

        ! Spectral directional advection (θ) with cyclic boundary
        cg_theta = 0.05d0 * sigma(k)
        if (cg_theta >= 0.0d0) then
            l_up = mod(l - 2 + ntheta, ntheta) + 1
        else
            l_up = mod(l, ntheta) + 1
        end if
        Nup = N(i,j,k,l_up)
        adv_theta = -cg_theta * dt / dtheta * (N0 - Nup)

        ! --- Final update ---
        N_new(i,j,k,l) = N0 + adv_x + adv_y + adv_sigma + adv_theta

    end subroutine propagate_wave

end module wavepropmod