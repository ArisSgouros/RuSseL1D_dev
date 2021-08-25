!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_computes(qinit_lo, qinit_hi)
!----------------------------------------------------------------------------------------------------------!
use parser_vars, only: geometry, Rg2_per_mon, ns_grafted_lo, ns_grafted_hi, gnode_lo, gnode_hi, out_phi,   &
                & out_field, out_q, out_phi_seg, out_ads_free, out_chainshape, out_brush_thickness,        &
                & grafted_lo_exist, grafted_hi_exist, matrix_exist, ns_matrix, ns_matrix_aux, rho_seg_bulk,&
                & edwards_solver, rho_seg_bulk, chainlen_matrix, grafted_hi_exist, grafted_lo_exist,       &
                & matrix_exist, ns_grafted_lo, ns_grafted_hi, chainlen_grafted_lo, chainlen_grafted_hi,    &
                & bc_lo_matrix, bc_hi_matrix, bc_lo_grafted, bc_hi_grafted, nx, export_phi_seg_id,         &
                & r_ads_lo, r_ads_hi, linear_solver
use flags,  only: F_lo, F_hi
use arrays, only: rx, coeff_nx, layer_area, phi_matrix, phi_gr_lo, phi_gr_hi, phi_total, qmatrix_final,    &
                & qgr_final_lo, qgr_final_hi, wa_ifc, wa_ifc_new, rx, dx, ds_matrix, ds_grafted_lo, rr,    &
                & ds_grafted_hi, coeff_ns_matrix, rs_grafted_lo, rs_grafted_hi, rs_matrix, rs_matrix_aux,  &
                & coeff_ns_grafted_lo, coeff_ns_grafted_hi
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
real(8), intent(in) :: qinit_lo, qinit_hi
!----------------------------------------------------------------------------------------------------------!
if (out_phi)   call export_phi(rx, phi_matrix, phi_gr_lo, phi_gr_hi, phi_total)
if (out_field) call export_field(rx, wa_ifc, wa_ifc_new)
if (out_q)     call export_q(qmatrix_final, ns_matrix_aux, nx, rs_matrix_aux, rx, "mataux")

if (matrix_exist) then
    if (out_q)          call export_q(qmatrix_final, ns_matrix, nx, rs_matrix, rx, "matrix")
    if (out_phi_seg)    call compute_phi_seg(export_phi_seg_id, chainlen_matrix, coeff_ns_matrix, ns_matrix, nx,  &
&                                            rx, qmatrix_final, qmatrix_final,"matrix")
    if (out_ads_free)   call compute_phi_ads_states(coeff_ns_matrix, rr, rx, dx, ds_matrix, wa_ifc, phi_matrix,   &
&                                                   qmatrix_final, bc_lo_matrix, bc_hi_matrix, geometry, r_ads_lo,&
&                                                    r_ads_hi, chainlen_matrix, ns_matrix, "matrix")
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_lo, edwards_solver,    &
&                                               linear_solver, bc_lo_matrix, bc_hi_matrix, qinit_lo, coeff_nx,    &
&                                               rr, layer_area, rx, nx, dx, chainlen_matrix, ns_matrix, ds_matrix,&
&                                               wa_ifc, phi_matrix, qmatrix_final, "matrix")
endif

if (grafted_lo_exist) then
    if (out_q)               call export_q(qgr_final_lo,  ns_grafted_lo, nx, rs_grafted_lo, rx, "gra_lo")
    if (out_phi_seg)         call compute_phi_seg(export_phi_seg_id, chainlen_grafted_lo, coeff_ns_grafted_lo,   &
                                                  ns_grafted_lo, nx, rx, qgr_final_lo, qmatrix_final,"gra_lo")
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_gr_lo, F_lo, "gra_lo")
    if (out_chainshape)    call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_lo, edwards_solver,&
&                                             linear_solver, bc_lo_grafted, bc_hi_grafted, qinit_lo, coeff_nx,&  
&                                             rr, layer_area, rx, nx, dx, chainlen_grafted_lo, ns_grafted_lo, &
&                                             ds_grafted_lo, wa_ifc, phi_gr_lo, qgr_final_lo, "gra_lo")
endif

if (grafted_hi_exist) then
    if (out_q)               call export_q(qgr_final_hi,  ns_grafted_hi, nx, rs_grafted_hi, rx, "gra_hi")
    if (out_phi_seg)         call compute_phi_seg(export_phi_seg_id, chainlen_grafted_hi, coeff_ns_grafted_hi,   &
&                                                 ns_grafted_hi, nx, rx, qgr_final_hi, qmatrix_final,"gra_hi")
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_gr_hi, F_hi, "gra_hi")
    if (out_chainshape)    call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_hi, edwards_solver, &
&                                             linear_solver, bc_lo_grafted, bc_hi_grafted, qinit_hi, coeff_nx,&
&                                             rr, layer_area, rx, nx, dx, chainlen_grafted_hi, ns_grafted_hi, &
&                                             ds_grafted_hi, wa_ifc, phi_gr_hi, qgr_final_hi, "gra_hi")
endif

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_computes
