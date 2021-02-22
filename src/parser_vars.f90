module parser_vars
!----------------------------------------------------------------------------------------------------------!
logical :: grafted_lo_exist, grafted_hi_exist, matrix_exist
logical :: read_field
logical :: wall_auto
logical :: square_gradient
logical :: out_phi, out_q, out_field, out_end_middle
logical :: out_chainshape, out_ads_free, out_brush_thickness
logical :: wall_hamaker, wall_square_well, wall_ramp, wall_vacuum, wall_hybrid

integer :: bc_hi_matrix, bc_lo_matrix, bc_hi_grafted, bc_lo_grafted
integer :: edwards_solver
integer :: contour_discret_scheme, spatial_discret_scheme
integer :: contour_integr_scheme, spatial_integr_scheme
integer :: nx, ns_matrix, ns_matrix_aux, ns_grafted_lo, ns_grafted_hi
integer :: gnode_lo, gnode_hi, geometry, max_iter
integer :: thermo_every, field_every, compute_every, check_stability_every

integer :: wall_type, wall_side

real(8) :: max_wa_error, r_critical
real(8) :: chainlen_matrix, chainlen_matrix_aux, chainlen_grafted_lo, chainlen_grafted_hi
real(8) :: ds_ave_matrix, ds_ave_grafted_lo, ds_ave_grafted_hi
real(8) :: lx, dx_ave
real(8) :: graft_pos, gdens_lo, gdens_hi, delta
real(8) :: Temp, beta, Pressure, CN, bond_length, frac, Rg2_per_mon
real(8) :: rho_mol_bulk, rho_seg_bulk, rho_mass_bulk
real(8) :: mon_mass, sphere_radius
real(8) :: Apol, Asolid, sig_pol, sig_solid
real(8) :: A_sq_well, sigma_sq_well
real(8) :: A_ramp, sigma_ramp
real(8) :: wall_pos, E_wall_target
real(8) :: k_gr, k_gr_tilde
!----------------------------------------------------------------------------------------------------------!
end module parser_vars