module paramsmod
    use gridmod
    implicit none

contains

    subroutine compute_params(N, Hs, Tm)
        real(8), intent(in) :: N(:,:,:,:)
        real(8), intent(out) :: Hs(nx, ny), Tm(nx, ny)
        real(8) :: E, m0, m1
        integer :: i, j, k, l
        real(8), parameter :: pi = 4.0d0 * atan(1.0d0)

        do j = 1, ny
            do i = 1, nx
                m0 = 0.0d0
                m1 = 0.0d0

                do l = 1, ntheta
                    do k = 1, nsigma
                        E = N(i,j,k,l) * sigma(k)
                        m0 = m0 + E * dtheta * dsigma(k)
                        m1 = m1 + E * dtheta * dsigma(k) / max(sigma(k), 1.0d-6)
                    end do
                end do

                Hs(i,j) = 4.0d0 * sqrt(m0) 
                Tm(i,j) = 2.0d0 * pi * m0 / max(m1, 1.0d-8)

            end do
        end do
    end subroutine compute_params
end module paramsmod