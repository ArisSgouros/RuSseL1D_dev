!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_energies_incompressible(free_energy)
!----------------------------------------------------------------------------------------------------------!
  use arrays, only: ufield, layer_area, coeff_nx, &
                       & qfinal_mxa, qfinal_mxb, qfinal_glo, qfinal_ghi, qfinal_glo_aux, qfinal_ghi_aux, &
                       & wa_ifc_new_kd1, wa_kd1, wa_bulk_kd1, wa_ifc_new_kd2, wa_kd2, wa_bulk_kd2, &
                       & phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_tot, &
                       & surface_area, volume, rx, &
                       & phi_kd1, phi_kd2
  use constants, only: pi
  use flags, only: F_both
  use parser_vars, only: wall_hamaker, rho_seg_bulk, lx, &
                       & ns_mxa, ns_mxb, ns_glo, ns_ghi, &
                       & exist_mxa, exist_mxb, exist_glo, exist_ghi, &
                       & k_gr, lx, nx, &
                       & chainlen_mxa, chainlen_mxb, chainlen_glo, chainlen_ghi, chainlen_bulk, &
                       & gnode_lo, gnode_hi, &
                       & gdens_lo, gdens_hi, beta, sig_solid, wall_pos, wall_side, asolid, &
                       & chi12, fh_rho_bulk, fh_rho_bulk_kd2, fh_press_bulk, fh_fraction
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  integer :: kk

  real(8), intent(out)     :: free_energy
  real(8), dimension(0:nx) :: prof_solid
  real(8), dimension(0:nx) :: prof_Flory, prof_fh_pressure
  real(8)                  :: get_nchains, get_part_func
  real(8)                  :: nchglo, nchghi, nchmxa, nchmxb, nchmxa_b1, nchmxa_b2, nchmxb_b1, nchmxb_b2
  real(8)                  :: part_func_mxa, part_func_mxb
  real(8)                  :: E_nkTlnQm, E_solid, E_NLnQ_mxa, E_NLnQ_mxb,  &
                       &      E_solid_solid, rho_seg_bulk_mx12, &
                       &      E_Flory, E_Flory_b1, E_Flory_b2, E_fh_pressure, &
                       &      E_fact, E_fact_b1, E_fact_b2, &
                       &      phi_kd1_b1, phi_kd1_b2, phi_kd2_b1, phi_kd2_b2
!----------------------------------------------------------------------------------------------------------!
  E_NLnQ_mxa = 0.d0
  E_NLnQ_mxb = 0.d0
  E_nkTlnQm = 0.d0
  E_solid = 0.d0
  E_solid_solid = 0.d0
  E_Flory = 0.d0
  E_Flory_b1 = 0.d0
  E_Flory_b2 = 0.d0
  E_fh_pressure = 0.d0
  E_fact = 0.d0
  E_fact_b1 = 0.d0
  E_fact_b2 = 0.d0
  nchglo = 0.d0
  nchghi = 0.d0
  nchmxa = 0.d0
  nchmxb = 0.d0
  nchmxa_b1 = 0.d0
  nchmxa_b2 = 0.d0
  nchmxb_b1 = 0.d0
  nchmxb_b2 = 0.d0
  phi_kd1_b1 = fh_rho_bulk
  phi_kd2_b1 = 1.0 - fh_rho_bulk
  phi_kd2_b2 = fh_rho_bulk_kd2
  phi_kd1_b2 = 1.0 - fh_rho_bulk_kd2

  prof_solid = 0.d0
  prof_Flory = 0.d0
  prof_fh_pressure = 0.d0
  do kk = 0, nx
    ! total contributions
    prof_solid(kk) = phi_tot(kk)*Ufield(kk)/beta
    ! flory huggins
    !prof_Flory(kk) = -chi12*(phi_kd1(kk)*phi_kd2(kk) - fh_rho_bulk*(1.d0 - fh_rho_bulk))
    prof_Flory(kk) = -chi12*(phi_kd1(kk)*phi_kd2(kk))
    prof_fh_pressure(kk) = -(0.5*(wa_ifc_new_kd1(kk) + wa_ifc_new_kd2(kk)) - 0.5d0*chi12)
  end do

  do kk = 0, nx
    ! total contributions
    E_solid = E_solid + coeff_nx(kk)*prof_solid(kk)*layer_area(kk)
    ! Flory Huggins
    E_Flory = E_Flory + coeff_nx(kk)*prof_Flory(kk)*layer_area(kk)
    E_fh_pressure = E_fh_pressure + coeff_nx(kk)*prof_fh_pressure(kk)*layer_area(kk)
  end do

! total contributions
  E_solid = E_solid*1.d-30*rho_seg_bulk
  E_Flory = E_Flory*1.d-30*rho_seg_bulk/beta
  E_Flory_b1 = -rho_seg_bulk*(fh_fraction*volume*1.d-30)*chi12*phi_kd1_b1*phi_kd2_b1
  E_Flory_b1 = E_Flory_b1/beta
  E_Flory_b2 = -rho_seg_bulk*((1.0 - fh_fraction)*volume*1.d-30)*chi12*phi_kd1_b2*phi_kd2_b2
  E_Flory_b2 = E_Flory_b2/beta

  E_fh_pressure = E_fh_pressure*1.d-30*rho_seg_bulk/beta

  
  E_fact = 0.d0
  E_fact_b1 = 0.d0
  E_fact_b2 = 0.d0

  ! mxa
  nchmxa = get_nchains(coeff_nx, nx, layer_area, phi_mxa, rho_seg_bulk, chainlen_mxa)
  part_func_mxa = get_part_func(nx, ns_mxa, layer_area, volume, coeff_nx, qfinal_mxa)
  E_NLnQ_mxa = -nchmxa/beta*log(part_func_mxa)
  E_fact = E_fact + (nchmxa*log(nchmxa) - nchmxa)/beta

  rho_seg_bulk_mx12 = rho_seg_bulk*phi_kd1_b1*fh_fraction
  nchmxa_b1 = get_nchains(coeff_nx, nx, layer_area, phi_tot, rho_seg_bulk_mx12, chainlen_mxa)
  E_fact_b1 = E_fact_b1 - (nchmxa_b1*log(nchmxa_b1) - nchmxa_b1)/beta

  rho_seg_bulk_mx12 = rho_seg_bulk*phi_kd1_b2*(1.0 - fh_fraction)
  nchmxa_b2 = get_nchains(coeff_nx, nx, layer_area, phi_tot, rho_seg_bulk_mx12, chainlen_mxa)
  E_fact_b2 = E_fact_b2 - (nchmxa_b2*log(nchmxa_b2) - nchmxa_b2)/beta

  ! mxb
  nchmxb = get_nchains(coeff_nx, nx, layer_area, phi_mxb, rho_seg_bulk, chainlen_mxb)
  part_func_mxb = get_part_func(nx, ns_mxb, layer_area, volume, coeff_nx, qfinal_mxb)
  E_NLnQ_mxb = -nchmxb/beta*log(part_func_mxb)
  E_fact = E_fact + (nchmxb*log(nchmxb) - nchmxb)/beta

  rho_seg_bulk_mx12 = rho_seg_bulk*phi_kd2_b1*(fh_fraction)
  nchmxb_b1 = get_nchains(coeff_nx, nx, layer_area, phi_tot, rho_seg_bulk_mx12, chainlen_mxb)
  E_fact_b1 = E_fact_b1 - (nchmxb_b1*log(nchmxb_b1) - nchmxb_b1)/beta

  rho_seg_bulk_mx12 = rho_seg_bulk*phi_kd2_b2*(1.0 - fh_fraction)
  nchmxb_b2 = get_nchains(coeff_nx, nx, layer_area, phi_tot, rho_seg_bulk_mx12, chainlen_mxb)
  E_fact_b2 = E_fact_b2 - (nchmxb_b2*log(nchmxb_b2) - nchmxb_b2)/beta

  E_nkTlnQm = 0.d0
  if (exist_glo) then
    nchglo = get_nchains(coeff_nx, nx, layer_area, phi_glo, rho_seg_bulk, chainlen_glo)
    E_nkTlnQm = E_nkTlnQm - nchglo/beta*log(qfinal_glo_aux(gnode_lo, ns_glo))
  end if
  if (exist_ghi) then
    nchghi = get_nchains(coeff_nx, nx, layer_area, phi_ghi, rho_seg_bulk, chainlen_ghi)
    E_nkTlnQm = E_nkTlnQm - nchghi/beta*log(qfinal_ghi_aux(gnode_hi, ns_ghi))
  end if

  if (wall_side .eq. F_both .and. wall_hamaker) then
    E_solid_solid = (surface_area*1.e-20)*(Asolid/PI)*                              &
 &                  (+(sig_solid*1.e-10)**6/(360.d0*((lx + 2.d0*wall_pos)*1.0d-10)**8) &
 &                    - 1.d0/(12.d0*((lx + 2.d0*wall_pos)*1.0d-10)**2))
  end if

! cast the energy contributions in units mJ/m^2
! total contributions
  E_solid = E_solid/(surface_area*1.e-20)*1e+3
  E_NLnQ_mxa = E_NLnQ_mxa/(surface_area*1.e-20)*1e+3
  E_NLnQ_mxb = E_NLnQ_mxb/(surface_area*1.e-20)*1e+3
  E_nkTlnQm = E_nkTlnQm/(surface_area*1.e-20)*1e+3
  E_Flory = E_Flory/(surface_area*1.e-20)*1e+3
  E_Flory_b1 = E_Flory_b1/(surface_area*1.e-20)*1e+3
  E_Flory_b2 = E_Flory_b2/(surface_area*1.e-20)*1e+3
  E_fh_pressure = E_fh_pressure/(surface_area*1.e-20)*1e+3
  E_solid_solid = E_solid_solid/(surface_area*1.e-20)*1e+3
  E_fact = E_fact/(surface_area*1.e-20)*1e+3
  E_fact_b1 = E_fact_b1/(surface_area*1.e-20)*1e+3
  E_fact_b2 = E_fact_b2/(surface_area*1.e-20)*1e+3
!estimate the free energy in mJ/m^2
  free_energy = E_fact + E_fact_b1 + E_fact_b2 + E_Flory + E_Flory_b1 + E_Flory_b2 + E_fh_pressure + E_NLnQ_mxa + E_NLnQ_mxb + E_nkTlnQm + E_solid + E_solid_solid

  open (unit=777, file="o.energies")
  write (777, '(21(A16))') "free_energy", "E_fact", "E_fact_b1", "E_fact_b2", "E_Flory", "E_Flory_b1","E_Flory_b2","E_fh_pressure", "NLnQmxa", "NLnQmxb", "nkTlnQm_ns", &
  &                        "N_ch1", "N_ch2", "N_ch1_b1", "N_ch1_b2", "N_ch2_b1", "N_ch2_b2", &
  &                        "E_solid", "solid_solid", "Qmxa", "Qmxb"

  write (777, '(21(E16.7))') free_energy, E_fact, E_fact_b1, E_fact_b2,  E_Flory, E_Flory_b1, E_Flory_b2, E_fh_pressure, E_NLnQ_mxa, E_NLnQ_mxb, E_nkTlnQm, &
  &                        nchmxa, nchmxb, nchmxa_b1, nchmxa_b2, nchmxb_b1, nchmxb_b2, &
  &                        E_solid, E_solid_solid, part_func_mxa, part_func_mxb
  close (777)

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_energies_incompressible
