!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_field_binary(wa, nx, field_filename_out)
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  integer, intent(in) :: nx

  character(len=*), intent(in) :: field_filename_out
  real(8), intent(in), dimension(0:nx) :: wa
!----------------------------------------------------------------------------------------------------------!

  open (unit=36, file=field_filename_out, form="unformatted")
  write (36) wa
  close (36)

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_field_binary
