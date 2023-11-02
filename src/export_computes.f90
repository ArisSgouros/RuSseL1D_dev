!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_computes(qinit_lo, qinit_hi)
!----------------------------------------------------------------------------------------------------------!
use parser_vars, only: geometry, ns_glo, ns_ghi, gnode_lo, gnode_hi, out_phi,   &
                & Rg2_per_mon_mxa, Rg2_per_mon_mxb, Rg2_per_mon_glo, Rg2_per_mon_ghi, &
                & out_field, out_q, out_phi_seg, out_ads_free, out_chainshape, out_brush_thickness,        &
                & out_equimolar,                                                                           &
                & exist_glo, exist_ghi, exist_mxa, ns_mxa, rho_seg_bulk,&
                & edwards_solver, rho_seg_bulk, chainlen_mxa, exist_ghi, exist_glo,       &
                & chainlen_mxb, &
                & exist_mxa, ns_glo, ns_ghi, chainlen_glo, chainlen_ghi,    &
                & exist_mxb, &
                & ns_mxb, &
                & bc_lo_mxa, bc_hi_mxa, bc_lo_grafted, bc_hi_grafted, nx, export_phi_seg_id,         &
                & bc_lo_mxb, bc_hi_mxb, &
                & r_ads_lo, r_ads_hi, linear_solver
use flags,  only: F_lo, F_hi
use arrays, only: rx, coeff_nx, layer_area, phi_mxa, phi_glo, phi_ghi, phi_tot, qfinal_mxa,    &
                & phi_mxb, qfinal_mxb,    &
                & qfinal_glo, qfinal_ghi, wa_ifc, wa_ifc_new, rx, dx, ds_mxa, ds_glo, rr,    &
                & ds_mxb, &
                & ds_ghi, coeff_ns_mxa, rs_glo, rs_ghi, rs_mxa,                &
                & rs_mxb, &
                & coeff_ns_mxb, &
                & qfinal_glo_aux, qfinal_ghi_aux,                                                       &
                & coeff_ns_glo, coeff_ns_ghi
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
real(8), intent(in) :: qinit_lo, qinit_hi
!----------------------------------------------------------------------------------------------------------!
if (out_phi)   call export_phi(rx, phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_tot)
if (out_field) call export_field(rx, wa_ifc, wa_ifc_new)

if (exist_mxa) then
    if (out_q)          call export_q(qfinal_mxa, ns_mxa, nx, rs_mxa, rx, "matrix")
    if (out_phi_seg)    call compute_phi_seg(export_phi_seg_id, chainlen_mxa, coeff_ns_mxa, ns_mxa, nx,  &
&                                            rx, qfinal_mxa, qfinal_mxa,"matrix")
    if (out_ads_free)   call compute_phi_ads_states(coeff_ns_mxa, rr, rx, dx, ds_mxa, wa_ifc, phi_mxa,   &
&                                                   qfinal_mxa, bc_lo_mxa, bc_hi_mxa, geometry, r_ads_lo,&
&                                                    r_ads_hi, chainlen_mxa, Rg2_per_mon_mxa, ns_mxa, "matrix")
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon_mxa, geometry, gnode_lo, edwards_solver,&
&                                               linear_solver, bc_lo_mxa, bc_hi_mxa, qinit_lo, coeff_nx,    &
&                                               rr, layer_area, rx, nx, dx, chainlen_mxa, ns_mxa, ds_mxa,&
&                                               wa_ifc, phi_mxa, qfinal_mxa, "matrix")
endif

if (exist_mxb) then
    if (out_q)          call export_q(qfinal_mxb, ns_mxb, nx, rs_mxb, rx, "mxb___")
    if (out_phi_seg)    call compute_phi_seg(export_phi_seg_id, chainlen_mxb, coeff_ns_mxb, ns_mxb, nx,  &
&                                            rx, qfinal_mxb, qfinal_mxb,"mxb___")
    if (out_ads_free)   call compute_phi_ads_states(coeff_ns_mxb, rr, rx, dx, ds_mxb, wa_ifc, phi_mxb,   &
&                                                   qfinal_mxb, bc_lo_mxb, bc_hi_mxb, geometry, r_ads_lo,&
&                                                    r_ads_hi, chainlen_mxb, Rg2_per_mon_mxb, ns_mxb, "mxb___")
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon_mxb, geometry, gnode_lo, edwards_solver,&
&                                               linear_solver, bc_lo_mxb, bc_hi_mxb, qinit_lo, coeff_nx,    &
&                                               rr, layer_area, rx, nx, dx, chainlen_mxb, ns_mxb, ds_mxb,&
&                                               wa_ifc, phi_mxb, qfinal_mxb, "mxb___")
endif

if (exist_glo) then
    if (out_q)               call export_q(qfinal_glo,  ns_glo, nx, rs_glo, rx, "gra_lo")
    if (out_phi_seg)         call compute_phi_seg(export_phi_seg_id, chainlen_glo, coeff_ns_glo,   &
                                                  ns_glo, nx, rx, qfinal_glo, qfinal_glo_aux,"gra_lo")
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_glo, F_lo, "gra_lo")
    if (out_chainshape)    call compute_chainshape(rho_seg_bulk, Rg2_per_mon_glo, geometry, gnode_lo, edwards_solver,&
&                                             linear_solver, bc_lo_grafted, bc_hi_grafted, qinit_lo, coeff_nx,&  
&                                             rr, layer_area, rx, nx, dx, chainlen_glo, ns_glo, &
&                                             ds_glo, wa_ifc, phi_glo, qfinal_glo, "gra_lo")
endif

if (exist_ghi) then
    if (out_q)               call export_q(qfinal_ghi,  ns_ghi, nx, rs_ghi, rx, "gra_hi")
    if (out_phi_seg)         call compute_phi_seg(export_phi_seg_id, chainlen_ghi, coeff_ns_ghi,   &
&                                                 ns_ghi, nx, rx, qfinal_ghi, qfinal_ghi_aux ,"gra_hi")
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_ghi, F_hi, "gra_hi")
    if (out_chainshape)    call compute_chainshape(rho_seg_bulk, Rg2_per_mon_ghi, geometry, gnode_hi, edwards_solver, &
&                                             linear_solver, bc_lo_grafted, bc_hi_grafted, qinit_hi, coeff_nx,&
&                                             rr, layer_area, rx, nx, dx, chainlen_ghi, ns_ghi, &
&                                             ds_ghi, wa_ifc, phi_ghi, qfinal_ghi, "gra_hi")
endif

if (out_equimolar) call compute_equimolar(nx, rx, coeff_nx, layer_area, phi_tot, "total_")

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_computes
