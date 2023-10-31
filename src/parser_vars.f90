!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module parser_vars
!----------------------------------------------------------------------------------------------------------!
logical :: grafted_lo_exist, grafted_hi_exist, matrixA_exist,matrixB_exist
logical :: read_field
logical :: square_gradient
logical :: out_phi, out_q, out_field, out_phi_seg
logical :: out_chainshape, out_ads_free, out_brush_thickness
logical :: out_equimolar
logical :: wall_auto
logical :: wall_hamaker, wall_square_well, wall_ramp, wall_vacuum, wall_hybrid
logical :: wall_custom, wall_table, wall_hamaker_well

integer :: bc_hi_matrixA,bc_hi_matrixB, bc_lo_matrixA,bc_lo_matrixB, bc_hi_grafted, bc_lo_grafted
integer :: edwards_solver, linear_solver
integer :: contour_discret_scheme, spatial_discret_scheme
integer :: contour_integr_scheme, spatial_integr_scheme
integer :: nx, ns_matrixA,ns_matrixB, ns_matrixA_aux,ns_matrixB_aux, ns_grafted_lo, ns_grafted_hi
integer :: gnode_lo, gnode_hi, geometry, max_iter
integer :: thermo_every, field_every, compute_every, check_stability_every
integer :: export_phi_seg_id
integer :: wall_type, wall_side
integer :: n_wall_custom_vars

real(8) :: max_wa_error, r_ads_lo, r_ads_hi
real(8) :: chainlen_matrixA,chainlen_matrixB, chainlen_matrixA_aux,chainlen_matrixB_aux, chainlen_grafted_lo, chainlen_grafted_hi, chainlen_bulk
real(8) :: ds_ave_matrixA,ds_ave_matrixB, ds_ave_grafted_lo, ds_ave_grafted_hi
real(8) :: lx, dx_ave
real(8) :: graft_pos, gdens_lo, gdens_hi, delta
real(8) :: Temp, beta, Pressure, CN, bond_length, frac, Rg2_per_mon
real(8) :: rho_mol_bulk, rho_seg_bulk, rho_mass_bulk
real(8) :: mon_mass, sphere_radius
real(8) :: Apol, Asolid, sig_pol, sig_solid
real(8) :: A_sq_well, sigma_sq_well
real(8) :: A_ramp, sigma_ramp
real(8) :: hamaker_well_constant, hamaker_well_rc
real(8) :: wall_pos, E_wall_target
real(8) :: k_gr, k_gr_tilde

character(100) :: field_in_filename
!----------------------------------------------------------------------------------------------------------!
end module parser_vars
