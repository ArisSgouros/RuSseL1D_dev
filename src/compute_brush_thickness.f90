subroutine compute_brush_thickness(nx, layer_area, rx, coeff_x, phi, side, chain_type)
!----------------------------------------------------------------------------------------------------------!
use constants
use flags
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: nx, side
integer             :: kk

real(8), intent(in), dimension(0:nx) :: coeff_x, phi, layer_area, rx

real(8), dimension(0:nx) :: phi_side

real(8) :: br_thick, M_99, h_99, sum1, sum2, summer

character(6), intent(in) :: chain_type

character(40)            :: filename = ''
!----------------------------------------------------------------------------------------------------------!

if (side.eq.F_lo) then
    phi_side = phi
elseif (side.eq.F_hi) then
    do kk=0,nx
        phi_side(kk) = phi(nx-kk)
    enddo
endif

sum1 = 0.d00
sum2 = 0.d00

do kk = 0, nx
    sum1 = sum1 + coeff_x(kk) * layer_area(kk) * phi_side(kk) * rx(kk)**2
    sum2 = sum2 + coeff_x(kk) * layer_area(kk) * phi_side(kk)
end do

br_thick = sqrt(sum1/sum2)

summer = 0.d00

do kk = 0, nx
    summer = summer  + coeff_x(kk) * phi_side(kk) * layer_area(kk)
end do

M_99 = 0.99e00 * summer

summer = 0.d00

do kk = 1, nx
    summer = summer + coeff_x(kk) * phi_side(kk) * layer_area(kk)
    if (summer.gt.M_99) then
        h_99 = rx(kk-1)
        exit
    endif
enddo

write(filename,'(''o.brush_thickness_'',A6,''.out.txt'')')chain_type

open (file=filename,unit=32)
write(32,'(A30,E16.9)') 'brush thickness', br_thick
write(32,'(A30,E16.9)') 'h99%', h_99
close(32)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_brush_thickness
