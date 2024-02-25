!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module parser_vars
!----------------------------------------------------------------------------------------------------------!
  logical :: exist_glo, exist_ghi, exist_mxa, exist_mxb
  logical :: read_field
  logical :: square_gradient
  logical :: out_phi, out_q, out_field, out_phi_seg
  logical :: out_chainshape, out_ads_free, out_brush_thickness
  logical :: out_equimolar
  logical :: wall_auto
  logical :: wall_hamaker, wall_square_well, wall_ramp, wall_vacuum, wall_hybrid
  logical :: wall_custom, wall_table, wall_hamaker_well

  integer :: bc_lo_mxa, bc_lo_mxb, bc_lo_grafted
  integer :: bc_hi_mxa, bc_hi_mxb, bc_hi_grafted
  integer :: edwards_solver, linear_solver
  integer :: contour_discret_scheme, spatial_discret_scheme
  integer :: contour_integr_scheme, spatial_integr_scheme
  integer :: nx
  integer :: ns_mxa, ns_mxb, ns_glo, ns_ghi
  integer :: gnode_lo, gnode_hi, geometry, max_iter
  integer :: thermo_every, field_every, compute_every, check_stability_every
  integer :: export_phi_seg_id
  integer :: wall_type, wall_side
  integer :: n_wall_custom_vars
  integer :: mxa_kind, mxb_kind, glo_kind, ghi_kind

  real(8) :: max_wa_error, r_ads_lo, r_ads_hi
  real(8) :: chainlen_mxa, chainlen_mxb, chainlen_glo, chainlen_ghi, chainlen_bulk
  real(8) :: ds_ave_mxa, ds_ave_mxb, ds_ave_glo, ds_ave_ghi
  real(8) :: lx, dx_ave
  real(8) :: graft_pos, gdens_lo, gdens_hi, delta
  real(8) :: Temp, beta, Pressure, frac
  real(8) :: CN_mxa, CN_mxb, CN_glo, CN_ghi
  real(8) :: bond_length_mxa, bond_length_mxb, bond_length_glo, bond_length_ghi
  real(8) :: Rg2_per_mon_mxa, Rg2_per_mon_mxb, Rg2_per_mon_glo, Rg2_per_mon_ghi
  real(8) :: rho_mol_bulk, rho_seg_bulk, rho_mass_bulk
  real(8) :: mon_mass, sphere_radius
  real(8) :: Apol, Asolid, sig_pol, sig_solid
  real(8) :: A_sq_well, sigma_sq_well
  real(8) :: A_ramp, sigma_ramp
  real(8) :: hamaker_well_constant, hamaker_well_rc
  real(8) :: wall_pos, E_wall_target
  real(8) :: k_gr, k_gr_tilde
  real(8) :: chi12

  character(100) :: field_in_filename
!----------------------------------------------------------------------------------------------------------!
end module parser_vars
