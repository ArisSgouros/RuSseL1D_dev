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
real(8), allocatable, dimension(:)   :: dx, rx, coeff_nx
real(8), allocatable, dimension(:)   :: ds_matrixA,ds_matrixB, ds_matrixA_aux,ds_matrixB_aux, ds_grafted_lo, ds_grafted_hi

real(8), allocatable, dimension(:)   :: coeff_ns_matrixA,coeff_ns_matrixB, coeff_ns_matrixA_aux,coeff_ns_matrixB_aux, coeff_ns_grafted_lo, coeff_ns_grafted_hi
real(8), allocatable, dimension(:)   :: rs_matrixA,rs_matrixB, rs_matrixA_aux,rs_matrixB_aux, rs_grafted_lo, rs_grafted_hi
real(8), allocatable, dimension(:)   :: dir_nodes_rdiag
real(8), allocatable, dimension(:)   :: dphi_dr, d2phi_dr2
real(8), allocatable, dimension(:)   :: phi_total, phi_matrixA,phi_matrixB, phi_gr_lo, phi_gr_hi
real(8), allocatable, dimension(:,:) :: qmatrixA,qmatrixB, qgr_lo, qgr_hi
real(8), allocatable, dimension(:,:) :: qmatrixA_final,qmatrixB_final, qgr_final_lo, qgr_final_hi
!----------------------------------------------------------------------------------------------------------!
end module arrays
