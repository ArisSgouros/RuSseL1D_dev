!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_field(rx, wa, wa_new)
!----------------------------------------------------------------------------------------------------------!
use parser_vars, only: nx
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: ii

real(8), intent(in), dimension(0:nx) :: rx
real(8), intent(in), dimension(0:nx) :: wa, wa_new
!----------------------------------------------------------------------------------------------------------!
open(unit=120, file="o.field")
write(120,'(3(A14,2X))') 'r', "wa_old", "wa_new"
do ii = 0, nx
    write (120,'(3(E14.7,2X))') rx(ii), wa(ii), wa_new(ii)
enddo
close(120)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_field
