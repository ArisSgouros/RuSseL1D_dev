!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_arrays()
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: ns_mxa, ns_mxb, ns_glo, ns_ghi, nx
  use arrays
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!

  allocate (ds_mxa(0:ns_mxa), ds_glo(0:ns_glo), ds_ghi(0:ns_ghi))
  allocate (ds_mxb(0:ns_mxb))
  allocate (rs_mxa(0:ns_mxa), rs_glo(0:ns_glo), rs_ghi(0:ns_ghi))
  allocate (rs_mxb(0:ns_mxb))
  allocate (Ufield(0:nx), df_drho(0:nx))
  allocate (wa_ifc_mxa(0:nx), wa_ifc_mxb(0:nx), wa_ifc_glo(0:nx), wa_ifc_ghi(0:nx))
  allocate (wa_kd1(0:nx), wa_kd2(0:nx))
  allocate (wa_ifc_kd1(0:nx), wa_ifc_kd2(0:nx))
  allocate (wa_ifc_new_kd1(0:nx), wa_ifc_new_kd2(0:nx))
  allocate (wa_ifc_backup_kd1(0:nx), wa_ifc_backup_kd2(0:nx))
  allocate (rr(0:nx), irr(0:nx), layer_area(0:nx))
  allocate (phi_tot(0:nx), phi_mxa(0:nx), phi_mxb(0:nx), phi_glo(0:nx), phi_ghi(0:nx))
  allocate (phi_kd1(0:nx), phi_kd2(0:nx), phi_new_kd1(0:nx), phi_new_kd2(0:nx))
  allocate (qmxa(0:nx, 2))
  allocate (qmxb(0:nx, 2))
  allocate (qfinal_mxa(0:nx, 0:ns_mxa))
  allocate (qfinal_mxb(0:nx, 0:ns_mxb))
  allocate (qglo(0:nx, 2), qglo_aux(0:nx, 2), qghi(0:nx, 2), qghi_aux(0:nx, 2))
  allocate (qfinal_glo(0:nx, 0:ns_glo), qfinal_glo_aux(0:nx, 0:ns_glo))
  allocate (qfinal_ghi(0:nx, 0:ns_ghi), qfinal_ghi_aux(0:nx, 0:ns_ghi))
  allocate (coeff_nx(0:nx))
  allocate (coeff_ns_mxa(0:ns_mxa), coeff_ns_glo(0:ns_glo), &
  &                                                                          coeff_ns_ghi(0:ns_ghi))
  allocate (coeff_ns_mxb(0:ns_mxb))
  allocate (dir_nodes_id(0:nx), dir_nodes_rdiag(0:nx))
  allocate (dphi_dr(0:nx), d2phi_dr2(0:nx))
  allocate (dx(0:nx), rx(0:nx))

  dir_nodes_id = 0

  wa_ifc_mxa = 0.d0
  wa_ifc_mxb = 0.d0
  wa_ifc_glo = 0.d0
  wa_ifc_ghi = 0.d0
  wa_kd1 = 0.d0
  wa_kd2 = 0.d0
  wa_ifc_kd1 = 0.d0
  wa_ifc_kd2 = 0.d0
  wa_ifc_new_kd1 = 0.d0
  wa_ifc_new_kd2 = 0.d0
  wa_ifc_backup_kd1 = 0.d0
  wa_ifc_backup_kd2 = 0.d0
  dir_nodes_rdiag = 0.d0
  rr = 0.d0
  layer_area = 0.d0
  Ufield = 0.d0
  phi_tot = 0.d0

  phi_mxa = 0.d0
  phi_mxb = 0.d0
  phi_glo = 0.d0
  phi_ghi = 0.d0
  phi_kd1 = 0.d0
  phi_kd2 = 0.d0
  phi_new_kd1 = 0.d0
  phi_new_kd2 = 0.d0
  qmxa = 0.d0
  qmxb = 0.d0
  qglo = 0.d0
  qglo_aux = 0.d0
  qghi = 0.d0
  qghi_aux = 0.d0
  qfinal_mxa = 0.d0
  qfinal_mxb = 0.d0
  qfinal_glo = 0.d0
  qfinal_glo_aux = 0.d0
  qfinal_ghi = 0.d0
  qfinal_ghi_aux = 0.d0
  dx = 0.d0
  ds_mxa = 0.d0
  ds_mxb = 0.d0
  ds_glo = 0.d0
  ds_ghi = 0.d0
  rx = 0.d0
  rs_mxa = 0.d0
  rs_mxb = 0.d0
  rs_glo = 0.d0
  rs_ghi = 0.d0
  coeff_nx = 0.d0
  coeff_ns_mxa = 0.d0
  coeff_ns_mxb = 0.d0
  coeff_ns_glo = 0.d0
  coeff_ns_ghi = 0.d0
  df_drho = 0.d0
  dphi_dr = 0.d0
  d2phi_dr2 = 0.d0
!----------------------------------------------------------------------------------------------------------!
end subroutine init_arrays
