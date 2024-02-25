!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module arrays
!----------------------------------------------------------------------------------------------------------!
  integer                              :: n_dir_nodes
  integer, allocatable, dimension(:)   :: dir_nodes_id

  real(8)                              :: volume, surface_area, wa_bulk
  real(8), allocatable, dimension(:)   :: wall_custom_vars
  real(8), allocatable, dimension(:)   :: rr, irr, layer_area
  real(8), allocatable, dimension(:)   :: wa, wa_ifc, wa_ifc_new, wa_ifc_backup, Ufield, df_drho
  real(8), allocatable, dimension(:)   :: wa_mxa, wa_mxb, wa_glo, wa_ghi
  real(8), allocatable, dimension(:)   :: wa_ifc_mxa, wa_ifc_mxb, wa_ifc_glo, wa_ifc_ghi
  real(8), allocatable, dimension(:)   :: wa_ifc_new_mxa, wa_ifc_new_mxb, wa_ifc_new_glo, wa_ifc_new_ghi
  real(8), allocatable, dimension(:)   :: dx, rx, coeff_nx
  real(8), allocatable, dimension(:)   :: ds_mxa, ds_mxb, ds_glo, ds_ghi

  real(8), allocatable, dimension(:)   :: coeff_ns_mxa, coeff_ns_mxb, coeff_ns_glo, coeff_ns_ghi
  real(8), allocatable, dimension(:)   :: rs_mxa, rs_mxb, rs_glo, rs_ghi
  real(8), allocatable, dimension(:)   :: dir_nodes_rdiag
  real(8), allocatable, dimension(:)   :: dphi_dr, d2phi_dr2
  real(8), allocatable, dimension(:)   :: phi_tot, phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_kd1, phi_kd2
  real(8), allocatable, dimension(:, :) :: qmxa, qmxb, qglo, qghi, qglo_aux, qghi_aux
  real(8), allocatable, dimension(:, :) :: qfinal_mxa, qfinal_mxb, qfinal_glo, qfinal_ghi, qfinal_glo_aux, qfinal_ghi_aux
!----------------------------------------------------------------------------------------------------------!
end module arrays
