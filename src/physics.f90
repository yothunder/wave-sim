module physicsmod
    use gridmod
    use windmod
    implicit none
    real(8), parameter :: A_in = 0.001d0, B_ds = 0.0005d0, m_ds = 4.0d0
    real(8), parameter :: sigma_peak = 1.0d0 ! modify as needed

contains
    subroutine calc_source_terms(N, N_new, it, i, j, k, l)
        implicit none
        integer, intent(in) :: it, i, j, k, l
        real(8), intent(inout) :: N(:,:,:,:), N_new(:,:,:,:)
        real(8) :: N0, S_tot, S_in, S_ds, beta_in, beta_ds, c, u_local

        N0 = N(i, j, k, l)
        ! S_in = 0.001d0
        ! S_ds = 0.002d0

        u_local = u10(it, i, j)

        beta_in = A_in * max((u_local / c - 1.0d0), 0.0d0)**2
        S_in = beta_in * N0 / sigma(k)

        beta_ds = B_ds * (sigma(k) / sigma_peak)**m_ds
        S_ds = beta_ds * N0 / sigma(k)

        ! S_tot = (S_in - S_ds) * N0 / sigma(k)
        S_tot = (S_in - S_ds)

        N_new(i, j, k, l) = N_new(i, j, k, l) + S_tot * dt
    end subroutine calc_source_terms
end module physicsmod