!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_field(rx, wa_kd1, wa_new_kd1, wa_kd2, wa_new_kd2)
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: nx
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  integer :: ii

  real(8), intent(in), dimension(0:nx) :: rx
  real(8), intent(in), dimension(0:nx) :: wa_kd1, wa_new_kd1, wa_kd2, wa_new_kd2
!----------------------------------------------------------------------------------------------------------!
  open (unit=120, file="o.field")
  write (120, '(5(A14,2X))') 'r', "wa_kd1", "wa_new_kd1", "wa_kd2", "wa_new_kd2"
  do ii = 0, nx
    write (120, '(5(E14.7,2X))') rx(ii), wa_kd1(ii), wa_new_kd1(ii), wa_kd2(ii), wa_new_kd2(ii)
  end do
  close (120)

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_field
