!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_scf_params()
!----------------------------------------------------------------------------------------------------------!
  use eos, only: rsl_N, T_tilde, P_tilde, rho_tilde_bulk, T_star, P_star, V_star, rho_star,             &
                        & eos_type, eos_rho_tilde_0
  use flags, only: F_sanchez_lacombe
  use parser_vars, only: chainlen_mxa, chainlen_mxb, chainlen_glo, chainlen_ghi, &
                        & ds_ave_mxa, ds_ave_mxb,  &
                        & ds_ave_glo, ds_ave_ghi, ns_mxa, ns_mxb, ns_glo, ns_ghi,         &
                        & exist_mxa, exist_mxb, exist_glo, exist_ghi, rho_seg_bulk,       &
                        & rho_mol_bulk, rho_mass_bulk, pressure, Temp, k_gr, k_gr_tilde, mon_mass,               &
                        & square_gradient, gdens_hi, gdens_lo, chainlen_bulk, &
                        & bond_length_mxa, bond_length_mxb, bond_length_glo, bond_length_ghi,       &
                        & CN_mxa, CN_mxb, CN_glo, CN_ghi,                                           &
                        & Rg2_per_mon_mxa, Rg2_per_mon_mxb, Rg2_per_mon_glo, Rg2_per_mon_ghi
  use constants, only: N_avog, boltz_const_Joule_molK, boltz_const_Joule_K, gr_cm3_to_kg_m3, iow, tol
  use write_helper, only: adjl
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  real(8) :: SL_kappa_T
!----------------------------------------------------------------------------------------------------------!
  write (iow, '(A85)') adjl("-----------------------------INITIALIZE THE SCF PARAMETERS---------------------------", 85)
  write (*, '(A85)') adjl("-----------------------------INITIALIZE THE SCF PARAMETERS---------------------------", 85)

  if (exist_mxa) then
    Rg2_per_mon_mxa = bond_length_mxa**2*CN_mxa/6.d00
    ns_mxa = 2*nint(0.5d0*chainlen_mxa/ds_ave_mxa)
    write (iow, '(3X,A45,F16.4," Angstrom")') adjl("mxa radious of gyration:", 45), sqrt(Rg2_per_mon_mxa*chainlen_mxa)
    write (*, '(3X,A45,F16.4," Angstrom")') adjl("mxa radious of gyration:", 45), sqrt(Rg2_per_mon_mxa*chainlen_mxa)
    write (iow, '(3X,A45,I16," nodes")') adjl("mxa nodes along chain contour:", 45), ns_mxa
    write (*, '(3X,A45,I16," nodes")') adjl("mxa nodes along chain contour:", 45), ns_mxa
  end if
  if (exist_mxb) then
    Rg2_per_mon_mxb = bond_length_mxb**2*CN_mxb/6.d00
    ns_mxb = 2*nint(0.5d0*chainlen_mxb/ds_ave_mxb)
    write (iow, '(3X,A45,F16.4," Angstrom")') adjl("mxb radious of gyration:", 45), sqrt(Rg2_per_mon_mxb*chainlen_mxb)
    write (*, '(3X,A45,F16.4," Angstrom")') adjl("mxb radious of gyration:", 45), sqrt(Rg2_per_mon_mxb*chainlen_mxb)
    write (iow, '(3X,A45,I16," nodes")') adjl("mxb nodes along chain contour:", 45), ns_mxb
    write (*, '(3X,A45,I16," nodes")') adjl("mxb nodes along chain contour:", 45), ns_mxb
  end if
  if (exist_glo) then
    Rg2_per_mon_glo = bond_length_glo**2*CN_glo/6.d00
    ns_glo = 2*nint(0.5d0*chainlen_glo/ds_ave_glo)
    write (iow, '(3X,A45,F16.4," Angstrom")') adjl("glo radious of gyration:", 45), sqrt(Rg2_per_mon_glo*chainlen_glo)
    write (*, '(3X,A45,F16.4," Angstrom")') adjl("glo radious of gyration:", 45), sqrt(Rg2_per_mon_glo*chainlen_glo)
    write (iow, '(3X,A45,I16," nodes")') adjl("glo nodes along chain contour:", 45), ns_glo
    write (*, '(3X,A45,I16," nodes")') adjl("glo nodes along chain contour:", 45), ns_glo
  end if
  if (exist_ghi) then
    Rg2_per_mon_ghi = bond_length_ghi**2*CN_ghi/6.d00
    ns_ghi = 2*nint(0.5d0*chainlen_ghi/ds_ave_ghi)
    write (iow, '(3X,A45,F16.4," Angstrom")') adjl("ghi radious of gyration:", 45), sqrt(Rg2_per_mon_ghi*chainlen_ghi)
    write (*, '(3X,A45,F16.4," Angstrom")') adjl("ghi radious of gyration:", 45), sqrt(Rg2_per_mon_ghi*chainlen_ghi)
    write (iow, '(3X,A45,I16," nodes")') adjl("ghi nodes along chain contour:", 45), ns_ghi
    write (*, '(3X,A45,I16," nodes")') adjl("ghi nodes along chain contour:", 45), ns_ghi
  end if

  if (eos_type .eq. F_sanchez_lacombe) then
    write (iow, '(3X,A45)') adjl("Computation of the mass density from SL EoS..", 45)
    write (*, '(3X,A45)') adjl("Computation of the mass density from SL EoS..", 45)
    V_star = boltz_const_Joule_K*T_star/P_star
    T_tilde = Temp/T_star
    P_tilde = Pressure/P_star
    rsl_N = (mon_mass*P_star)/(rho_star*1.d03*boltz_const_Joule_molK*T_star)
    rho_tilde_bulk = eos_rho_tilde_0(T_tilde, P_tilde, rsl_N*chainlen_bulk)
    rho_mass_bulk = rho_tilde_bulk*rho_star
    rho_mass_bulk = rho_mass_bulk
    write (iow, '(3X,A45,F16.4," g/cm3")') adjl("mass density was recomputed as:", 45), rho_mass_bulk/gr_cm3_to_kg_m3
    write (*, '(3X,A45,F16.4," g/cm3")') adjl("mass density was recomputed as:", 45), rho_mass_bulk/gr_cm3_to_kg_m3

    if (square_gradient .and. k_gr_tilde > tol) k_gr = 2.d0*P_star*rsl_N**2*V_star**(8.d0/3.d0)*k_gr_tilde

    write (iow, '(3X,A45,E16.8," J m^5")') adjl("Sanchez-Lacombe influence parameter:", 45), k_gr
    write (*, '(3X,A45,E16.8," J m^5")') adjl("Sanchez-Lacombe influence parameter:", 45), k_gr

    SL_kappa_T = 1.0d0/(P_star*T_tilde*rho_tilde_bulk*(1.d0/(1.d0/rho_tilde_bulk - 1.d0) + &
&                                          1.d0/(rsl_N*chainlen_bulk)) - 2.d0*rho_tilde_bulk**2*P_star)
    write (iow, '(3X,A45,E16.4," Pa-1")') adjl("Sanchez-Lacombe isothermal compressibility:", 45), SL_kappa_T
    write (*, '(3X,A45,E16.4," Pa-1")') adjl("Sanchez-Lacombe isothermal compressibility:", 45), SL_kappa_T
  end if

  rho_mol_bulk = rho_mass_bulk/mon_mass*gr_cm3_to_kg_m3
  rho_seg_bulk = rho_mol_bulk*N_avog

  write (iow, '(3X,A45,F16.4," mol/m3")') adjl("molar density in bulk", 45), rho_mol_bulk
  write (*, '(3X,A45,F16.4," mol/m3")') adjl("molar density in bulk", 45), rho_mol_bulk
!----------------------------------------------------------------------------------------------------------!
end subroutine init_scf_params
