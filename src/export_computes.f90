!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_computes(qinit_lo, qinit_hi, iter)
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: geometry, &
                  & ns_mxa, ns_mxb, ns_glo, ns_ghi, &
                  & gnode_lo, gnode_hi, out_phi,   &
                  & Rg2_per_mon_mxa, Rg2_per_mon_mxb, Rg2_per_mon_glo, Rg2_per_mon_ghi, &
                  & out_field, out_q, out_phi_seg, out_ads_free, out_chainshape, out_brush_thickness,        &
                  & out_equimolar,                                                                           &
                  & exist_mxa, exist_mxb, exist_glo, exist_ghi, &
                  & rho_seg_bulk,&
                  & chainlen_mxa, chainlen_mxb, chainlen_glo, chainlen_ghi,    &
                  & bc_lo_mxa, bc_lo_mxb, bc_lo_grafted, &
                  & bc_hi_mxa, bc_hi_mxb, bc_hi_grafted, &
                  & nx, export_phi_seg_id,         &
                  & r_ads_lo, r_ads_hi, linear_solver, edwards_solver
  use flags, only: F_lo, F_hi
  use arrays, only: rx, coeff_nx, layer_area, &
                  & phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_kd1, phi_kd2, phi_tot, &
                  & qfinal_mxa, qfinal_mxb, qfinal_glo, qfinal_ghi, qfinal_glo_aux, qfinal_ghi_aux, &
                  & wa_ifc_kd1, wa_ifc_new_kd1, wa_ifc_kd2, wa_ifc_new_kd2, &
                  & rx, dx, ds_mxa, ds_glo, rr,    &
                  & ds_mxa, ds_mxb, ds_glo, ds_ghi, &
                  & coeff_ns_mxa, coeff_ns_mxb, coeff_ns_glo, coeff_ns_ghi, &
                  & rs_mxa, rs_mxb, rs_glo, rs_ghi
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  real(8), intent(in) :: qinit_lo, qinit_hi
  integer, intent(in) :: iter
!----------------------------------------------------------------------------------------------------------!
  if (out_phi) call export_phi(rx, phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_kd1, phi_kd2, phi_tot, iter)
  if (out_field) call export_field(rx, wa_ifc_kd1, wa_ifc_new_kd1, wa_ifc_kd2, wa_ifc_new_kd2, iter)

  if (exist_mxa) then
    if (out_q) call export_q(qfinal_mxa, ns_mxa, nx, rs_mxa, rx, "matrix")
    if (out_phi_seg) call compute_phi_seg(export_phi_seg_id, chainlen_mxa, coeff_ns_mxa, ns_mxa, nx,  &
&                                            rx, qfinal_mxa, qfinal_mxa, "matrix")
    if (out_ads_free) call compute_phi_ads_states(coeff_ns_mxa, rr, rx, dx, ds_mxa, wa_ifc_kd1, phi_mxa,   &
&                                                   qfinal_mxa, bc_lo_mxa, bc_hi_mxa, geometry, r_ads_lo,&
&                                                    r_ads_hi, chainlen_mxa, Rg2_per_mon_mxa, ns_mxa, "matrix")
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon_mxa, geometry, gnode_lo, edwards_solver,&
&                                               linear_solver, bc_lo_mxa, bc_hi_mxa, qinit_lo, coeff_nx,    &
&                                               rr, layer_area, rx, nx, dx, chainlen_mxa, ns_mxa, ds_mxa,&
&                                               wa_ifc_kd1, phi_mxa, qfinal_mxa, "matrix")
  end if

  if (exist_mxb) then
    if (out_q) call export_q(qfinal_mxb, ns_mxb, nx, rs_mxb, rx, "mxb___")
    if (out_phi_seg) call compute_phi_seg(export_phi_seg_id, chainlen_mxb, coeff_ns_mxb, ns_mxb, nx,  &
&                                            rx, qfinal_mxb, qfinal_mxb, "mxb___")
    if (out_ads_free) call compute_phi_ads_states(coeff_ns_mxb, rr, rx, dx, ds_mxb, wa_ifc_kd1, phi_mxb,   &
&                                                   qfinal_mxb, bc_lo_mxb, bc_hi_mxb, geometry, r_ads_lo,&
&                                                    r_ads_hi, chainlen_mxb, Rg2_per_mon_mxb, ns_mxb, "mxb___")
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon_mxb, geometry, gnode_lo, edwards_solver,&
&                                               linear_solver, bc_lo_mxb, bc_hi_mxb, qinit_lo, coeff_nx,    &
&                                               rr, layer_area, rx, nx, dx, chainlen_mxb, ns_mxb, ds_mxb,&
&                                               wa_ifc_kd1, phi_mxb, qfinal_mxb, "mxb___")
  end if

  if (exist_glo) then
    if (out_q) call export_q(qfinal_glo, ns_glo, nx, rs_glo, rx, "gra_lo")
    if (out_phi_seg) call compute_phi_seg(export_phi_seg_id, chainlen_glo, coeff_ns_glo, &
                                          ns_glo, nx, rx, qfinal_glo, qfinal_glo_aux, "gra_lo")
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_glo, F_lo, "gra_lo")
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon_glo, geometry, gnode_lo, edwards_solver,&
&                                             linear_solver, bc_lo_grafted, bc_hi_grafted, qinit_lo, coeff_nx,&
&                                             rr, layer_area, rx, nx, dx, chainlen_glo, ns_glo, &
&                                             ds_glo, wa_ifc_kd1, phi_glo, qfinal_glo, "gra_lo")
  end if

  if (exist_ghi) then
    if (out_q) call export_q(qfinal_ghi, ns_ghi, nx, rs_ghi, rx, "gra_hi")
    if (out_phi_seg) call compute_phi_seg(export_phi_seg_id, chainlen_ghi, coeff_ns_ghi,   &
&                                                 ns_ghi, nx, rx, qfinal_ghi, qfinal_ghi_aux, "gra_hi")
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_ghi, F_hi, "gra_hi")
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon_ghi, geometry, gnode_hi, edwards_solver, &
&                                             linear_solver, bc_lo_grafted, bc_hi_grafted, qinit_hi, coeff_nx,&
&                                             rr, layer_area, rx, nx, dx, chainlen_ghi, ns_ghi, &
&                                             ds_ghi, wa_ifc_kd1, phi_ghi, qfinal_ghi, "gra_hi")
  end if

  if (out_equimolar) call compute_equimolar(nx, rx, coeff_nx, layer_area, phi_tot, "total_")

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_computes
