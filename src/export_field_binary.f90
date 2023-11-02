!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_field_binary(wa, nx)
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: nx

real(8), intent(in), dimension(0:nx) :: wa
!----------------------------------------------------------------------------------------------------------!
character(len=50) :: field_filename_out = ''

write(field_filename_out,'("o.field.bin")')
open(unit=36, file = field_filename_out, form = "unformatted")
write(36) wa
close(36)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_field_binary
