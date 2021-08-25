!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

function get_nchains(coeff_x, nx, layer_area, phi, rho_seg_bulk, chainlen)
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: nx
integer             :: kk

real(8), intent(in)                  :: rho_seg_bulk, chainlen
real(8), intent(in), dimension(0:nx) :: coeff_x, layer_area, phi
real(8)                              :: summ, get_nchains
!----------------------------------------------------------------------------------------------------------!
summ = 0.d0
do kk = 0, nx
    summ = summ  + coeff_x(kk)* phi(kk) * layer_area(kk)
enddo

get_nchains = summ * 1.0d-30 * rho_seg_bulk / chainlen

return
!----------------------------------------------------------------------------------------------------------!
end function get_nchains
