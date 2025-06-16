module physicsmod
    use gridmod
    use windmod
    implicit none
    real(8), parameter :: A_in = 0.250d0, B_ds = 0.0002d0, m_ds = 3.30d0
    real(8), parameter :: sigma_peak = 0.20d0
    real(8), parameter :: g = 9.81d0

contains
    subroutine calc_source_terms(N, N_new, it, i, j, k, l)
        implicit none
        integer, intent(in) :: it, i, j, k, l
        real(8), intent(inout) :: N(:,:,:,:), N_new(:,:,:,:)
        real(8) :: N0, S_tot, S_in, S_ds, beta_in, beta_ds, c
        real(8) :: u_local, v_local, theta_wind, cosfac

        N0 = N(i, j, k, l)

        u_local = u10(it, i, j)
        v_local = v10(it, i, j)
        theta_wind = atan2(v_local, u_local)
        cosfac = cos(theta(l) - theta_wind)
        c = g / max(sigma(k), 1.0d-6)

        beta_in = A_in * max((sqrt(u_local**2 + v_local**2) / c - 0.50d0), 0.0d0)**2
        S_in = beta_in * N0 * cosfac / max(sigma(k), 1.0d-6)

        beta_ds = B_ds * (sigma(k) / sigma_peak)**m_ds
        S_ds = beta_ds * N0 / max(sigma(k), 1.0d-6)

        S_tot = (S_in - S_ds)
        if (abs(S_tot * dt) > 0.5 * N0) then
            S_tot = sign(1.0d0, S_tot) * 0.5 * N0 / dt
        end if

        N_new(i, j, k, l) = N(i,j,k,l) + S_tot * dt
    end subroutine calc_source_terms
end module physicsmod