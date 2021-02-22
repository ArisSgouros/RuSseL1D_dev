subroutine export_field_binary(wa, nx)
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: nx

real(8), intent(in), dimension(0:nx) :: wa
!----------------------------------------------------------------------------------------------------------!
character(len=50) :: field_filename_out = ''

write(field_filename_out,'(''field.out.bin'')')
open(unit=36, file = field_filename_out, form = "unformatted")
write(36) wa
close(36)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_field_binary