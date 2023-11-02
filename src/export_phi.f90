!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_phi(rx, phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_tot)
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: nx
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  integer :: ii

  real(8), intent(in), dimension(0:nx) :: rx
  real(8), intent(in), dimension(0:nx) :: phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_tot
!----------------------------------------------------------------------------------------------------------!
!export density profiles and fields

  open (unit=120, file="o.phi")
  write (120, '(6(A14,2X))') 'r', "phi_mxa", "phi_mxb", "phi_glo", "phi_ghi", "phi_tot"
  do ii = 0, nx
    write (120, '(6(E14.7E3,2X))') rx(ii), phi_mxa(ii), phi_mxb(ii), phi_glo(ii), phi_ghi(ii), phi_tot(ii)
  end do
  close (120)

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_phi
