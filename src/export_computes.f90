!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_computes(qinit_lo, qinit_hi)
!----------------------------------------------------------------------------------------------------------!
use parser_vars, only: geometry, Rg2_per_mon, ns_grafted_lo, ns_grafted_hi, gnode_lo, gnode_hi, out_phi,   &
                & out_field, out_q, out_phi_seg, out_ads_free, out_chainshape, out_brush_thickness,        &
                & out_equimolar,                                                                           &
                & grafted_lo_exist, grafted_hi_exist, matrixA_exist, ns_matrixA, ns_matrixA_aux, rho_seg_bulk,&
                & edwards_solver, rho_seg_bulk, chainlen_matrixA, grafted_hi_exist, grafted_lo_exist,       &
                & matrixA_exist, ns_grafted_lo, ns_grafted_hi, chainlen_grafted_lo, chainlen_grafted_hi,    &
                & bc_lo_matrixA, bc_hi_matrixA, bc_lo_grafted, bc_hi_grafted, nx, export_phi_seg_id,         &
                & r_ads_lo, r_ads_hi, linear_solver
use flags,  only: F_lo, F_hi
use arrays, only: rx, coeff_nx, layer_area, phi_matrixA, phi_gr_lo, phi_gr_hi, phi_total, qmatrixA_final,    &
                & qgr_final_lo, qgr_final_hi, wa_ifc, wa_ifc_new, rx, dx, ds_matrixA, ds_grafted_lo, rr,    &
                & ds_grafted_hi, coeff_ns_matrixA, rs_grafted_lo, rs_grafted_hi, rs_matrixA, rs_matrixA_aux,  &
                & coeff_ns_grafted_lo, coeff_ns_grafted_hi
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
real(8), intent(in) :: qinit_lo, qinit_hi
!----------------------------------------------------------------------------------------------------------!
if (out_phi)   call export_phi(rx, phi_matrixA, phi_gr_lo, phi_gr_hi, phi_total)
if (out_field) call export_field(rx, wa_ifc, wa_ifc_new)
if (out_q)     call export_q(qmatrixA_final, ns_matrixA_aux, nx, rs_matrixA_aux, rx, "mataux")

if (matrixA_exist) then
    if (out_q)          call export_q(qmatrixA_final, ns_matrixA, nx, rs_matrixA, rx, "matrixA")
    if (out_phi_seg)    call compute_phi_seg(export_phi_seg_id, chainlen_matrixA, coeff_ns_matrixA, ns_matrixA, nx,  &
&                                            rx, qmatrixA_final, qmatrixA_final,"matrixA")
    if (out_ads_free)   call compute_phi_ads_states(coeff_ns_matrixA, rr, rx, dx, ds_matrixA, wa_ifc, phi_matrixA,   &
&                                                   qmatrixA_final, bc_lo_matrixA, bc_hi_matrixA, geometry, r_ads_lo,&
&                                                    r_ads_hi, chainlen_matrixA, ns_matrixA, "matrixA")
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_lo, edwards_solver,    &
&                                               linear_solver, bc_lo_matrixA, bc_hi_matrixA, qinit_lo, coeff_nx,    &
&                                               rr, layer_area, rx, nx, dx, chainlen_matrixA, ns_matrixA, ds_matrixA,&
&                                               wa_ifc, phi_matrixA, qmatrixA_final, "matrixA")
endif

if (grafted_lo_exist) then
    if (out_q)               call export_q(qgr_final_lo,  ns_grafted_lo, nx, rs_grafted_lo, rx, "gra_lo")
    if (out_phi_seg)         call compute_phi_seg(export_phi_seg_id, chainlen_grafted_lo, coeff_ns_grafted_lo,   &
                                                  ns_grafted_lo, nx, rx, qgr_final_lo, qmatrixA_final,"gra_lo")
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_gr_lo, F_lo, "gra_lo")
    if (out_chainshape)    call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_lo, edwards_solver,&
&                                             linear_solver, bc_lo_grafted, bc_hi_grafted, qinit_lo, coeff_nx,&  
&                                             rr, layer_area, rx, nx, dx, chainlen_grafted_lo, ns_grafted_lo, &
&                                             ds_grafted_lo, wa_ifc, phi_gr_lo, qgr_final_lo, "gra_lo")
endif

if (grafted_hi_exist) then
    if (out_q)               call export_q(qgr_final_hi,  ns_grafted_hi, nx, rs_grafted_hi, rx, "gra_hi")
    if (out_phi_seg)         call compute_phi_seg(export_phi_seg_id, chainlen_grafted_hi, coeff_ns_grafted_hi,   &
&                                                 ns_grafted_hi, nx, rx, qgr_final_hi, qmatrixA_final,"gra_hi")
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_gr_hi, F_hi, "gra_hi")
    if (out_chainshape)    call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_hi, edwards_solver, &
&                                             linear_solver, bc_lo_grafted, bc_hi_grafted, qinit_hi, coeff_nx,&
&                                             rr, layer_area, rx, nx, dx, chainlen_grafted_hi, ns_grafted_hi, &
&                                             ds_grafted_hi, wa_ifc, phi_gr_hi, qgr_final_hi, "gra_hi")
endif

if (out_equimolar) call compute_equimolar(nx, rx, coeff_nx, layer_area, phi_total, "total_")

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_computes
