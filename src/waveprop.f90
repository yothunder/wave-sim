module wavepropmod
    use gridmod
    implicit none

contains

    subroutine propagate_wave(N, N_new, i, j, k, l)
        implicit none
        integer, intent(in) :: i, j, k, l
        real(8), intent(inout) :: N(:,:,:,:), N_new(:,:,:,:)
        real(8) :: cg, cg_x, cg_y, cg_sigma, cg_theta
        real(8) :: adv_x, adv_y, adv_sigma, adv_theta
        real(8) :: N0, Nup
        real(8), parameter :: g = 9.81d0

        N0 = N(i, j, k, l)

        ! --- Dynamic group velocity ---
        cg = 0.5d0 * g / sigma(k)
        cg_x = cg * cos(theta(l))
        cg_y = cg * sin(theta(l))
        
        ! Spatial propagation using an upwind scheme
        if (i > 1) then
            adv_x = N(i-1,j,k,l)
        else
            adv_x = N0
        end if
        adv_x = -cg_x * dt / dx * (N0 - Nup)

        if (j > 1) then
            adv_y = N(i,j-1,k,l)
        else
            adv_y = N0
        end if
        adv_y = -cg_y * dt / dy * (N0 - Nup)

        ! --- Spectral propagation (frequency σ) ---
        cg_sigma = 0.0d0   ! Assume no frequency shifting for now
        adv_sigma = 0.0d0

        ! --- Spectral propagation (direction θ) ---
        cg_theta = 0.0d0   ! Assume no refraction yet
        adv_theta = 0.0d0

        ! --- Final update ---
        N_new(i,j,k,l) = N0 + adv_x + adv_y + adv_sigma + adv_theta

    end subroutine propagate_wave

end module wavepropmod