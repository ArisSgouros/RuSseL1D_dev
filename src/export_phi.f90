!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_phi(rx, phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_kd1, phi_kd2, phi_tot, iter)
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: nx, export_multi
  use write_helper, only: adjl
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  integer :: ii

  real(8), intent(in), dimension(0:nx) :: rx
  real(8), intent(in), dimension(0:nx) :: phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_kd1, phi_kd2, phi_tot
  character(200) :: filename
  character(200) :: string_iter
  integer, intent(in) :: iter
!----------------------------------------------------------------------------------------------------------!
!export density profiles and fields

  string_iter = ""
  if (export_multi) then
    write (string_iter, '(I10)') iter
  end if
  write (filename, '("o.phi",A10)') adjl(string_iter,10)

  open (unit=120, file=filename)
  write (120, '(8(A14,2X))') 'r', "phi_mxa", "phi_mxb", "phi_glo", "phi_ghi", "phi_kd1", "phi_kd2", "phi_tot"
  do ii = 0, nx
    write (120, '(8(E14.7E3,2X))') rx(ii), phi_mxa(ii), phi_mxb(ii), phi_glo(ii), phi_ghi(ii), phi_kd1(ii), phi_kd2(ii), phi_tot(ii)
  end do
  close (120)

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_phi
