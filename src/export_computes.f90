subroutine export_computes(qinit_lo, qinit_hi)
!----------------------------------------------------------------------------------------------------------!
use parser_vars
use constants
use flags
use arrays
use eos
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
real(8),intent(in) :: qinit_lo, qinit_hi
!----------------------------------------------------------------------------------------------------------!

if (out_phi)   call export_phi(rx, coeff_nx, layer_area, phi_matrix, phi_gr_lo, phi_gr_hi, phi_total)
if (out_field) call export_field(rx, wa_ifc, wa_ifc_new)
if (out_q)     call export_q(qmatrix_final, ns_matrix_aux, nx, rs_matrix_aux, rx, 'mataux')

if (matrix_exist) then
    if (out_q)          call export_q(qmatrix_final, ns_matrix, nx, rs_matrix, rx, 'matrix')
    if (out_end_middle) call compute_phi_end_middle(ns_matrix, nx, rx, qmatrix_final, qmatrix_final,'matrix')
    if (out_ads_free)   call compute_matrix_ads_free(coeff_ns_matrix, rr, rx, dx, ds_matrix, wa_ifc, phi_matrix, qmatrix_final)
    if (out_chainshape) call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_lo, edwards_solver, &
&             bc_lo_matrix, bc_hi_matrix, qinit_lo, coeff_nx, rr, layer_area, &
&             rx, nx, dx, chainlen_matrix, ns_matrix, ds_matrix, wa_ifc, phi_matrix, qmatrix_final, 'matrix')
endif

if (grafted_lo_exist) then
    if (out_q)               call export_q(qgr_final_lo,  ns_grafted_lo, nx, rs_grafted_lo, rx, 'gra_lo')
    if (out_end_middle)      call compute_phi_end_middle(ns_grafted_lo, nx, rx, qgr_final_lo, qmatrix_final,'gra_lo')
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_gr_lo, F_lo, 'gra_lo')
    if (out_chainshape)      call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_lo, edwards_solver, &
&             bc_lo_grafted, bc_hi_grafted, qinit_lo, coeff_nx, &
&             rr, layer_area, rx, nx, dx, chainlen_grafted_lo, ns_grafted_lo, ds_grafted_lo, wa_ifc, phi_gr_lo, qgr_final_lo, 'gra_lo')
endif

if (grafted_hi_exist) then
    if (out_q)               call export_q(qgr_final_hi,  ns_grafted_hi, nx, rs_grafted_hi, rx, 'gra_hi')
    if (out_end_middle)      call compute_phi_end_middle(ns_grafted_hi, nx, rx, qgr_final_hi, qmatrix_final,'gra_hi')
    if (out_brush_thickness) call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_gr_hi, F_hi, 'gra_hi')
    if (out_chainshape)      call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_hi, edwards_solver, &
&             bc_lo_grafted, bc_hi_grafted, qinit_hi, coeff_nx, &
&             rr, layer_area, rx, nx, dx, chainlen_grafted_hi, ns_grafted_hi, ds_grafted_hi, wa_ifc, phi_gr_hi, qgr_final_hi, 'gra_hi')
endif

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_computes
