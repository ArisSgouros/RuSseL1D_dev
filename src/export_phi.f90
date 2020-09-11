subroutine export_phi(rx, coeff_nx, layer_area, phi_matrix, phi_gr_lo, phi_gr_hi, phi_total)
!----------------------------------------------------------------------------------------------------------!
use parser_vars
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: ii, kk

real(8), intent(in), dimension(0:nx)              :: rx, coeff_nx, layer_area
real(8), intent(in), dimension(0:nx)              :: phi_matrix, phi_gr_lo, phi_gr_hi, phi_total

real(8), dimension(0:nx) :: cumul_matrix, cumul_gr_lo, cumul_gr_hi
!----------------------------------------------------------------------------------------------------------!

!export density profiles and fields
cumul_matrix = 0.d0
cumul_gr_lo = 0.d0
cumul_gr_hi = 0.d0
do ii=0,nx
    do kk=1,ii
        cumul_matrix(ii) = cumul_matrix(ii) + coeff_nx(kk) * phi_matrix(kk) * layer_area(kk)
        cumul_gr_lo(ii)  = cumul_gr_lo(ii)  + coeff_nx(kk) * phi_gr_lo(kk)  * layer_area(kk)
        cumul_gr_hi(ii)  = cumul_gr_hi(ii)  + coeff_nx(kk) * phi_gr_hi(kk)  * layer_area(kk)
    enddo
enddo
do ii=0,nx
    if (cumul_matrix(ii) > 0.d0) cumul_matrix(ii) = cumul_matrix(ii) / cumul_matrix(nx)
    if (cumul_gr_lo(ii) > 0.d0)  cumul_gr_lo(ii)  = cumul_gr_lo(ii)  / cumul_gr_lo(nx)
    if (cumul_gr_hi(ii) > 0.d0)  cumul_gr_hi(ii)  = cumul_gr_hi(ii)  / cumul_gr_hi(nx)
enddo

open(unit=120, file='o.phi.out.txt')
write(120,'(8(A14,2X))') 'r', 'phi_matrix', 'phi_gr_lo', 'phi_gr_hi', 'phi_tot', &
&                         'cumul_matrix', 'cumul_gr_lo', 'cumul_gr_hi'
do ii = 0, nx
    write (120,'(8(E14.7,2X))') rx(ii), phi_matrix(ii), phi_gr_lo(ii), phi_gr_hi(ii), phi_total(ii), &
&                                cumul_matrix(ii), cumul_gr_lo(ii), cumul_gr_hi(ii)
enddo
close(120)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_phi
