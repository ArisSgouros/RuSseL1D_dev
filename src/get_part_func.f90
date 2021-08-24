!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

function get_part_func(nx, ns, layer_area, volume, coeff_x, q_final)
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: nx, ns
integer             :: kk

real(8), intent(in), dimension(0:nx)      :: coeff_x, layer_area
real(8), intent(in), dimension(0:nx,0:ns) :: q_final
real(8)                                   :: volume, summer, get_part_func
!----------------------------------------------------------------------------------------------------------!
summer = 0.d0
do kk = 0, nx
    summer = summer + coeff_x(kk) * q_final(kk,ns) * layer_area(kk)
enddo
get_part_func = summer/volume

return
!----------------------------------------------------------------------------------------------------------!
end function get_part_func
