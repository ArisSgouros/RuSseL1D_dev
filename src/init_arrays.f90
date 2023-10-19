!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_arrays()
!----------------------------------------------------------------------------------------------------------!
use parser_vars, only: ns_matrixA, ns_matrix_auxA, ns_grafted_lo, ns_grafted_hi, nx
use arrays
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!

allocate(ds_matrixA(0:ns_matrixA),ds_matrix_auxA(0:ns_matrix_auxA),ds_grafted_lo(0:ns_grafted_lo),ds_grafted_hi(0:ns_grafted_hi))
allocate(rs_matrixA(0:ns_matrixA),rs_matrix_auxA(0:ns_matrix_auxA),rs_grafted_lo(0:ns_grafted_lo),rs_grafted_hi(0:ns_grafted_hi))
allocate(wa(0:nx),wa_ifc(0:nx),wa_ifc_new(0:nx),wa_ifc_backup(0:nx), Ufield(0:nx),df_drho(0:nx))
allocate(rr(0:nx),irr(0:nx),layer_area(0:nx))
allocate(phi_total(0:nx),phi_matrixA(0:nx),phi_gr_lo(0:nx),phi_gr_hi(0:nx))
allocate(qmatrixA(0:nx,2),qgr_lo(0:nx,2),qgr_hi(0:nx,2))
allocate(qmatrix_finalA(0:nx,0:ns_matrix_auxA),qgr_final_lo(0:nx,0:ns_grafted_lo),qgr_final_hi(0:nx,0:ns_grafted_hi))
allocate(coeff_nx(0:nx))
allocate(coeff_ns_matrixA(0:ns_matrixA),coeff_ns_matrix_auxA(0:ns_matrix_auxA),coeff_ns_grafted_lo(0:ns_grafted_lo), &
&                                                                          coeff_ns_grafted_hi(0:ns_grafted_hi))
allocate(dir_nodes_id(0:nx), dir_nodes_rdiag(0:nx))
allocate(dphi_dr(0:nx), d2phi_dr2(0:nx))
allocate(dx(0:nx),rx(0:nx))

dir_nodes_id = 0

wa_bulk             = 0.d0
wa                  = 0.d0
wa_ifc              = 0.d0
wa_ifc_new          = 0.d0
wa_ifc_backup       = 0.d0
dir_nodes_rdiag     = 0.d0
rr                  = 0.d0
layer_area          = 0.d0
Ufield              = 0.d0
phi_total           = 0.d0

phi_matrixA          = 0.d0
phi_gr_lo           = 0.d0
phi_gr_hi           = 0.d0
qmatrixA             = 0.d0
qgr_lo              = 0.d0
qgr_hi              = 0.d0
qmatrix_finalA       = 0.d0
qgr_final_lo        = 0.d0
qgr_final_hi        = 0.d0
dx                  = 0.d0
ds_matrixA           = 0.d0
ds_matrix_auxA       = 0.d0
ds_grafted_lo       = 0.d0
ds_grafted_hi       = 0.d0
rx                  = 0.d0
rs_matrixA           = 0.d0
rs_matrix_auxA       = 0.d0
rs_grafted_lo       = 0.d0
rs_grafted_hi       = 0.d0
coeff_nx            = 0.d0
coeff_ns_matrixA     = 0.d0
coeff_ns_matrix_auxA = 0.d0
coeff_ns_grafted_lo = 0.d0
coeff_ns_grafted_hi = 0.d0
df_drho             = 0.d0
dphi_dr             = 0.d0
d2phi_dr2           = 0.d0
!----------------------------------------------------------------------------------------------------------!
end subroutine init_arrays
