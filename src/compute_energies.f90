!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_energies(free_energy)
!----------------------------------------------------------------------------------------------------------!
use eos,         only: eos_ff, eos_df_drho
use arrays,      only: ufield, layer_area, coeff_nx, qmatrixA_final, wa_ifc_new, wa, wa_bulk, phi_matrixA,   & 
                     & phi_gr_hi, phi_gr_lo, phi_total, dphi_dr, d2phi_dr2, qgr_final_lo, qgr_final_hi,    &
                     & surface_area, volume, rx, qgr_final_lo_aux, qgr_final_hi_aux
use constants,   only: pi
use flags,       only: F_both
use parser_vars, only: wall_hamaker, rho_seg_bulk, lx, ns_grafted_lo, ns_grafted_hi, ns_matrixA, wall_side, &
                     & grafted_lo_exist, grafted_hi_exist, matrixA_exist, k_gr, lx, nx, chainlen_matrixA,    &
                     & chainlen_grafted_lo, chainlen_grafted_hi, chainlen_bulk, gnode_lo, gnode_hi,        &
                     & gdens_lo, gdens_hi, beta, Rg2_per_mon, sig_solid, wall_pos, asolid
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: kk

real(8), intent(out)     :: free_energy
real(8), dimension(0:nx) :: prof_eos_f, prof_eos_rdfdr, prof_sgt_f, prof_sgt_rdfdr, prof_field, prof_solid
real(8), dimension(0:nx) :: prof_rhoglo_wifc, prof_rhoghi_wifc, prof_rhom_wifc, prof_solid_glo, prof_solid_ghi, prof_solid_m
real(8), dimension(0:nx) :: phi_end, rho_end, A_stretch
real(8)                  :: get_nchains, get_part_func, nchgr_lo, nchgr_hi, part_func_matrixA
real(8)                  :: E_eos_f, E_eos_rdfdr, E_rhoVkTQ, E_nkTlnQm, E_sgt_f, E_sgt_rdfdr, E_field, E_solid
real(8)                  :: E_rhoglo_wifc, E_rhoghi_wifc, E_rhomx_wifc, E_solid_glo, E_solid_ghi, E_solid_m, E_solid_solid
real(8)                  :: E_stretch_glo, E_stretch_add_glo, E_stretch_ghi, E_stretch_add_ghi, dz
!----------------------------------------------------------------------------------------------------------!
E_eos_f           = 0.d0
E_eos_rdfdr       = 0.d0
E_rhoVkTQ         = 0.d0
E_nkTlnQm         = 0.d0
E_sgt_f           = 0.d0
E_sgt_rdfdr       = 0.d0
E_field           = 0.d0
E_solid           = 0.d0
E_rhoglo_wifc     = 0.d0
E_rhoghi_wifc     = 0.d0
E_rhomx_wifc      = 0.d0
E_solid_glo       = 0.d0
E_solid_ghi       = 0.d0
E_solid_m         = 0.d0
E_solid_solid     = 0.d0
E_stretch_glo     = 0.d0
E_stretch_ghi     = 0.d0
E_stretch_add_glo = 0.d0
E_stretch_add_ghi = 0.d0

prof_eos_f       = 0.d0
prof_eos_rdfdr   = 0.d0
prof_sgt_f       = 0.d0
prof_sgt_rdfdr   = 0.d0
prof_field       = 0.d0
prof_solid       = 0.d0
prof_rhoglo_wifc = 0.d0
prof_rhoghi_wifc = 0.d0
prof_rhom_wifc   = 0.d0
prof_solid_glo   = 0.d0
prof_solid_ghi   = 0.d0
prof_solid_m     = 0.d0


do kk = 0, nx
    ! total contributions
    prof_eos_f(kk)       =  eos_ff(phi_total(kk))
    prof_eos_rdfdr(kk)   = -phi_total(kk) * eos_df_drho(phi_total(kk))
    prof_field(kk)       = -phi_total(kk) * wa(kk)
    prof_solid(kk)       =  phi_total(kk) * Ufield(kk) / beta
    prof_sgt_f(kk)       =  0.5d0 * k_gr * (rho_seg_bulk * dphi_dr(kk))**2
    prof_sgt_rdfdr(kk)   =  k_gr * (rho_seg_bulk**2 * phi_total(kk) * d2phi_dr2(kk))
    ! partial contributions
    prof_rhoglo_wifc(kk) = -phi_gr_lo(kk)  * wa_ifc_new(kk)
    prof_rhoghi_wifc(kk) = -phi_gr_hi(kk)  * wa_ifc_new(kk)
    prof_rhom_wifc(kk)   = -phi_matrixA(kk) * wa_ifc_new(kk)
    prof_solid_glo(kk)   =  phi_gr_lo(kk)  * Ufield(kk) / beta
    prof_solid_ghi(kk)   =  phi_gr_hi(kk)  * Ufield(kk) / beta
    prof_solid_m(kk)     =  phi_matrixA(kk) * Ufield(kk) / beta
end do

if (matrixA_exist) then       
    do kk = 0, nx
        prof_eos_f(kk)     = prof_eos_f(kk)     - eos_ff(1.d0)
        prof_eos_rdfdr(kk) = prof_eos_rdfdr(kk) + 1.d0*eos_df_drho(1.d0)
        prof_field(kk)     = prof_field(kk)     + 1.d0*wa_bulk
    end do
endif

do kk = 0, nx
    ! total contributions
    E_eos_f     = E_eos_f     + coeff_nx(kk) * prof_eos_f(kk)     * layer_area(kk)
    E_eos_rdfdr = E_eos_rdfdr + coeff_nx(kk) * prof_eos_rdfdr(kk) * layer_area(kk)
    E_field     = E_field     + coeff_nx(kk) * prof_field(kk)     * layer_area(kk)
    E_sgt_f     = E_sgt_f     + coeff_nx(kk) * prof_sgt_f(kk)     * layer_area(kk)
    E_sgt_rdfdr = E_sgt_rdfdr + coeff_nx(kk) * prof_sgt_rdfdr(kk) * layer_area(kk)
    E_solid     = E_solid     + coeff_nx(kk) * prof_solid(kk)     * layer_area(kk)
    ! partial contributions
    E_rhoglo_wifc = E_rhoglo_wifc + coeff_nx(kk) * prof_rhoglo_wifc(kk) * layer_area(kk)
    E_rhoghi_wifc = E_rhoghi_wifc + coeff_nx(kk) * prof_rhoghi_wifc(kk) * layer_area(kk)
    E_rhomx_wifc  = E_rhomx_wifc  + coeff_nx(kk) * prof_rhom_wifc(kk)   * layer_area(kk)
    E_solid_glo   = E_solid_glo   + coeff_nx(kk) * prof_solid_glo(kk)   * layer_area(kk)
    E_solid_ghi   = E_solid_ghi   + coeff_nx(kk) * prof_solid_ghi(kk)   * layer_area(kk)
    E_solid_m     = E_solid_m     + coeff_nx(kk) * prof_solid_m(kk)     * layer_area(kk)
end do

! total contributions
E_eos_f       = E_eos_f     * 1.d-30
E_eos_rdfdr   = E_eos_rdfdr * 1.d-30 * rho_seg_bulk
E_field       = E_field     * 1.d-30 * rho_seg_bulk / beta
E_solid       = E_solid     * 1.d-30 * rho_seg_bulk
E_sgt_f       = E_sgt_f     * 1.d-30
E_sgt_rdfdr   = E_sgt_rdfdr * 1.d-30
! partial contributions
E_rhoglo_wifc = E_rhoglo_wifc * 1.d-30 * rho_seg_bulk / beta
E_rhoghi_wifc = E_rhoghi_wifc * 1.d-30 * rho_seg_bulk / beta
E_rhomx_wifc  = E_rhomx_wifc  * 1.d-30 * rho_seg_bulk / beta
E_solid_glo   = E_solid_glo   * 1.d-30 * rho_seg_bulk
E_solid_ghi   = E_solid_ghi   * 1.d-30 * rho_seg_bulk
E_solid_m     = E_solid_m     * 1.d-30 * rho_seg_bulk

if (matrixA_exist) then       
   part_func_matrixA = get_part_func(nx, ns_matrixA, layer_area, volume, coeff_nx, qmatrixA_final)
   E_rhoVkTQ = volume*1.0d-30*rho_seg_bulk/beta/chainlen_matrixA*(1.d0-part_func_matrixA)
endif

E_nkTlnQm = 0.d0
if (grafted_lo_exist) then
   nchgr_lo  = get_nchains(coeff_nx, nx, layer_area, phi_gr_lo, rho_seg_bulk, chainlen_grafted_lo)
   E_nkTlnQm = E_nkTlnQm -nchgr_lo / beta * log(qgr_final_lo_aux(gnode_lo,ns_grafted_lo))
endif
if (grafted_hi_exist) then
   nchgr_hi  = get_nchains(coeff_nx, nx, layer_area, phi_gr_hi, rho_seg_bulk, chainlen_grafted_hi)
   E_nkTlnQm = E_nkTlnQm -nchgr_hi / beta * log(qgr_final_hi_aux(gnode_hi,ns_grafted_hi))
endif

if (wall_side.eq.F_both.and.wall_hamaker) then
   E_solid_solid = (surface_area * 1.e-20) * (Asolid/PI) *                              &
&                  ( + (sig_solid*1.e-10)**6 / (360.d0*((lx+2.d0*wall_pos)*1.0d-10)**8) &
&                    -         1.d0          / (12.d0 *((lx+2.d0*wall_pos)*1.0d-10)**2) )
endif
!
! estimate the contribution due to chain stretching
!
if (grafted_lo_exist) then
   !calculate the profile of the chain ends
   do kk = 0, nx
       phi_end(kk) = 1.d0 / chainlen_grafted_lo * (qgr_final_lo(kk,ns_grafted_lo) * qgr_final_lo_aux(kk,0))
   enddo
   rho_end = phi_end * rho_seg_bulk*1.d-30
   ! calculate the stretching free energy for every chain end abstained by dz from the grafting point
   do kk = 0, nx
       dz = rx(kk)-rx(gnode_lo)
       A_stretch(kk) = 3.0/(2.d0 * beta * (Rg2_per_mon*chainlen_grafted_lo*6.d0)) * dz**2
   end do
   ! calculate the total stretching free energy from all chain ends
   do kk = 0, nx
       E_stretch_glo = E_stretch_glo + A_stretch(kk) * rho_end(kk) * coeff_nx(kk) * layer_area(kk)
   end do
   ! calculate the additive constant
   E_stretch_add_glo = -0.5d0 * surface_area * gdens_lo / beta
   ! Test: The integral of phi_end must equal the number of grafted chains
   !rho_end_int = 0.d0
   !do kk = 0, nx
   !    rho_end_int = rho_end_int + coeff_nx(kk)* rho_end(kk) * layer_area(kk)
   !end do
endif
if (grafted_hi_exist) then
   !calculate the profile of the chain ends
   do kk = 0, nx
       phi_end(kk) = 1.d0 / chainlen_grafted_hi * (qgr_final_hi(kk,ns_grafted_hi) * qgr_final_hi_aux(kk,0))
   enddo
   rho_end = phi_end * rho_seg_bulk*1.d-30
   ! calculate the stretching free energy for every chain end abstained by dz from the grafting point
   do kk = 0, nx
       dz = rx(kk)-rx(gnode_hi)
       A_stretch(kk) = 3.0/(2.d0 * beta * (Rg2_per_mon*chainlen_grafted_hi*6.d0)) * dz**2
   end do
   ! calculate the total stretching free energy from all chain ends
   do kk = 0, nx
       E_stretch_ghi = E_stretch_ghi + A_stretch(kk) * rho_end(kk) * coeff_nx(kk) * layer_area(kk)
   end do
   ! calculate the additive constant
   E_stretch_add_ghi = -0.5d0 * surface_area * gdens_hi / beta
   ! Test: The integral of phi_end must equal the number of grafted chains
   !rho_end_int = 0.d0
   !do kk = 0, nx
   !    rho_end_int = rho_end_int + coeff_nx(kk)* rho_end(kk) * layer_area(kk)
   !end do
endif

! cast the energy contributions in units mJ/m^2
! total contributions
E_eos_f     = E_eos_f     / (surface_area * 1.e-20) * 1e+3
E_eos_rdfdr = E_eos_rdfdr / (surface_area * 1.e-20) * 1e+3
E_field     = E_field     / (surface_area * 1.e-20) * 1e+3
E_solid     = E_solid     / (surface_area * 1.e-20) * 1e+3
E_rhoVkTQ   = E_rhoVkTQ   / (surface_area * 1.e-20) * 1e+3
E_nkTlnQm   = E_nkTlnQm   / (surface_area * 1.e-20) * 1e+3
E_sgt_f     = E_sgt_f     / (surface_area * 1.e-20) * 1e+3
E_sgt_rdfdr = E_sgt_rdfdr / (surface_area * 1.e-20) * 1e+3
! partial contributions
E_rhoglo_wifc     = E_rhoglo_wifc     / (surface_area * 1.e-20) * 1e+3
E_rhoghi_wifc     = E_rhoghi_wifc     / (surface_area * 1.e-20) * 1e+3
E_rhomx_wifc      = E_rhomx_wifc      / (surface_area * 1.e-20) * 1e+3
E_solid_glo       = E_solid_glo       / (surface_area * 1.e-20) * 1e+3
E_solid_ghi       = E_solid_ghi       / (surface_area * 1.e-20) * 1e+3
E_solid_m         = E_solid_m         / (surface_area * 1.e-20) * 1e+3
E_stretch_glo     = E_stretch_glo     / (surface_area * 1.e-20) * 1e+3
E_stretch_ghi     = E_stretch_ghi     / (surface_area * 1.e-20) * 1e+3
E_stretch_add_glo = E_stretch_add_glo / (surface_area * 1.e-20) * 1e+3
E_stretch_add_ghi = E_stretch_add_ghi / (surface_area * 1.e-20) * 1e+3
E_solid_solid     = E_solid_solid     / (surface_area * 1.e-20) * 1e+3
!estimate the free energy in mJ/m^2
free_energy = E_eos_f + E_sgt_f + E_eos_rdfdr + E_sgt_rdfdr + E_rhoVkTQ + E_nkTlnQm + E_solid_solid

open(unit=777, file = "o.energies")
write(777,'(20(A16))')  "eos_f", "eos_rdfdr", "sgt_f", "sgt_rdfdr", "rhoVkTQ", "nkTlnQm_ns",       &
&                       "field", "solid", "free_energy", "E_rhoglo_wifc", "E_rhoghi_wifc",         &
&                       "E_rhomx_wifc", "E_str_glo", "E_str_add_glo", "E_str_ghi", "E_str_add_ghi",&
&                       "solid_glo", "solid_ghi", "solid_m", "solid_solid"

write(777,'(20(E16.7))') E_eos_f, E_eos_rdfdr, E_sgt_f, E_sgt_rdfdr, E_rhoVkTQ, E_nkTlnQm, E_field,       &
&                        E_solid, free_energy, E_rhoglo_wifc, E_rhoghi_wifc, E_rhomx_wifc, E_stretch_glo, &
&                        E_stretch_add_glo, E_stretch_ghi, E_stretch_add_ghi, E_solid_glo, E_solid_ghi,   &
&                        E_solid_m, E_solid_solid
close(777)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_energies
