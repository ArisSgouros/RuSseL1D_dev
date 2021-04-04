subroutine export_phi(rx, phi_matrix, phi_gr_lo, phi_gr_hi, phi_total)
!----------------------------------------------------------------------------------------------------------!
use parser_vars, only: nx
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: ii

real(8), intent(in), dimension(0:nx) :: rx
real(8), intent(in), dimension(0:nx) :: phi_matrix, phi_gr_lo, phi_gr_hi, phi_total
!----------------------------------------------------------------------------------------------------------!
!export density profiles and fields

open(unit=120, file="o.phi")
write(120,'(5(A14,2X))')'r', "phi_matrix", "phi_gr_lo", "phi_gr_hi", "phi_tot"
do ii = 0, nx
    write (120,'(5(E14.7E3,2X))') rx(ii), phi_matrix(ii), phi_gr_lo(ii), phi_gr_hi(ii), phi_total(ii)
enddo
close(120)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_phi
