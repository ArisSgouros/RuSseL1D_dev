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

call export_phi(rx, coeff_nx, layer_area, phi_matrix, phi_gr_lo, phi_gr_hi, phi_total)
call export_field(rx, wa, wa_new)
call export_q(qmatrix_final, ns_matrix_aux, nx, rs_matrix_aux, rx, 'mataux')

if (matrix_exist) then
    call export_q(qmatrix_final, ns_matrix,     nx, rs_matrix,     rx, 'matrix')
    call compute_phi_end_middle(ns_matrix, nx, rx, qmatrix_final, qmatrix_final,'matrix')
    call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_lo, edwards_solver, &
&             bc_lo_matrix, bc_hi_matrix, qinit_lo, coeff_nx, &
&             rr, layer_area, rx, nx, dx, chainlen_matrix, ns_matrix, ds_matrix, wa, phi_matrix, qmatrix_final, 'matrix')
    call compute_matrix_ads_free(coeff_ns_matrix, rr, rx, dx, ds_matrix, wa, phi_matrix, qmatrix_final)
endif

if (grafted_lo_exist) then
    call export_q(qgr_final_lo,  ns_grafted_lo, nx, rs_grafted_lo, rx, 'gra_lo')
    call compute_phi_end_middle(ns_grafted_lo, nx, rx, qgr_final_lo, qmatrix_final,'gra_lo')
    call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_lo, edwards_solver, &
&             bc_lo_grafted, bc_hi_grafted, qinit_lo, coeff_nx, &
&             rr, layer_area, rx, nx, dx, chainlen_grafted_lo, ns_grafted_lo, ds_grafted_lo, wa, phi_gr_lo, qgr_final_lo, 'gra_lo')
    call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_gr_lo, F_lo, 'gra_lo')
endif

if (grafted_hi_exist) then
    call export_q(qgr_final_hi,  ns_grafted_hi, nx, rs_grafted_hi, rx, 'gra_hi')
    call compute_phi_end_middle(ns_grafted_hi, nx, rx, qgr_final_hi, qmatrix_final,'gra_hi')
    call compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode_hi, edwards_solver, &
&             bc_lo_grafted, bc_hi_grafted, qinit_hi, coeff_nx, &
&             rr, layer_area, rx, nx, dx, chainlen_grafted_hi, ns_grafted_hi, ds_grafted_hi, wa, phi_gr_hi, qgr_final_hi, 'gra_hi')
    call compute_brush_thickness(nx, layer_area, rx, coeff_nx, phi_gr_hi, F_hi, 'gra_hi')
endif

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_computes
