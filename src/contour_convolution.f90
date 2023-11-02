!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine contour_convolution(chainlen, nx, ns, coeff, q1, q2, phi)
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: nx, ns
integer             :: tt, kk

real(8), intent(in)                       :: chainlen
real(8), intent(in), dimension(0:ns)      :: coeff
real(8), intent(in), dimension(0:nx,0:ns) :: q1, q2
real(8), intent(out), dimension(0:nx)     :: phi
real(8)                                   :: summ
!----------------------------------------------------------------------------------------------------------!
do kk = 0, nx
    summ = 0.d0
    do tt = 0, ns
        summ = summ + coeff(tt) * q1(kk,tt) * q2(kk,ns-tt)
    end do
    phi(kk) = summ / chainlen
enddo

return
!----------------------------------------------------------------------------------------------------------!
end subroutine contour_convolution
