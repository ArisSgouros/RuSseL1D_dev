!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_field(rx, wa_kd1, wa_new_kd1, wa_kd2, wa_new_kd2, iter)
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: nx, export_multi
  use write_helper, only: adjl
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  integer :: ii

  real(8), intent(in), dimension(0:nx) :: rx
  real(8), intent(in), dimension(0:nx) :: wa_kd1, wa_new_kd1, wa_kd2, wa_new_kd2
  character(200) :: filename
  character(200) :: string_iter
  integer, intent(in) :: iter
!----------------------------------------------------------------------------------------------------------!

  string_iter = ""
  if (export_multi) then
    write (string_iter, '(I10)') iter
  end if
  write (filename, '("o.field",A10)') adjl(string_iter,10)

  open (unit=120, file=filename)
  write (120, '(5(A14,2X))') 'r', "wa_kd1", "wa_new_kd1", "wa_kd2", "wa_new_kd2"
  do ii = 0, nx
    write (120, '(5(E14.7,2X))') rx(ii), wa_kd1(ii), wa_new_kd1(ii), wa_kd2(ii), wa_new_kd2(ii)
  end do
  close (120)

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_field
