!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_energies(free_energy)
!----------------------------------------------------------------------------------------------------------!
use eos,         only: eos_ff, eos_df_drho
use arrays,      only: ufield, layer_area, coeff_nx, &
                     & qfinal_mxa, qfinal_mxb, qfinal_glo, qfinal_ghi, qfinal_glo_aux, qfinal_ghi_aux, &
                     & wa_ifc_new, wa, wa_bulk,  & 
                     & phi_mxa, phi_mxb, phi_glo, phi_ghi, phi_tot, &
                     & dphi_dr, d2phi_dr2, &
                     & surface_area, volume, rx
use constants,   only: pi
use flags,       only: F_both
use parser_vars, only: wall_hamaker, rho_seg_bulk, lx, &
                     & ns_mxa, ns_mxb, ns_glo, ns_ghi, &
                     & exist_mxa, exist_mxb, exist_glo, exist_ghi, &
                     & k_gr, lx, nx, &
                     & chainlen_mxa, chainlen_mxb, chainlen_glo, chainlen_ghi, chainlen_bulk, &
                     & gnode_lo, gnode_hi, &
                     & gdens_lo, gdens_hi, beta, sig_solid, wall_pos, wall_side, asolid, &
                     & Rg2_per_mon_glo, Rg2_per_mon_ghi
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: kk

real(8), intent(out)     :: free_energy
real(8), dimension(0:nx) :: prof_eos_f, prof_eos_rdfdr, prof_sgt_f, prof_sgt_rdfdr, prof_field, prof_solid
real(8), dimension(0:nx) :: prof_rho_wifc_mxa, prof_rho_wifc_mxb, prof_rho_wifc_glo, prof_rho_wifc_ghi, &
                     &      prof_solid_mxa, prof_solid_mxb, prof_solid_glo, prof_solid_ghi
real(8), dimension(0:nx) :: phi_end, rho_end, A_stretch
real(8)                  :: get_nchains, get_part_func, nchglo, nchghi
real(8)                  :: part_func_mxa, part_func_mxb
real(8)                  :: E_eos_f, E_eos_rdfdr, E_nkTlnQm, E_sgt_f, E_sgt_rdfdr, E_field, E_solid, &
                     &      E_rhoVkTQ_mxa, E_rhoVkTQ_mxb
real(8)                  :: E_rho_wifc_mxa, E_rho_wifc_mxb, E_rho_wifc_glo, E_rho_wifc_ghi, &
                     &      E_solid_mxa, E_solid_mxb, E_solid_glo, E_solid_ghi, E_solid_solid
real(8)                  :: E_stretch_glo, E_stretch_add_glo, E_stretch_ghi, E_stretch_add_ghi, dz
!----------------------------------------------------------------------------------------------------------!
E_eos_f           = 0.d0
E_eos_rdfdr       = 0.d0
E_rhoVkTQ_mxa         = 0.d0
E_rhoVkTQ_mxb         = 0.d0
E_nkTlnQm         = 0.d0
E_sgt_f           = 0.d0
E_sgt_rdfdr       = 0.d0
E_field           = 0.d0
E_solid           = 0.d0
E_rho_wifc_glo     = 0.d0
E_rho_wifc_ghi     = 0.d0
E_rho_wifc_mxa      = 0.d0
E_rho_wifc_mxb      = 0.d0
E_solid_glo       = 0.d0
E_solid_ghi       = 0.d0
E_solid_mxa         = 0.d0
E_solid_mxb         = 0.d0
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
prof_rho_wifc_glo = 0.d0
prof_rho_wifc_ghi = 0.d0
prof_rho_wifc_mxa   = 0.d0
prof_rho_wifc_mxb   = 0.d0
prof_solid_glo   = 0.d0
prof_solid_ghi   = 0.d0
prof_solid_mxa     = 0.d0
prof_solid_mxb     = 0.d0


do kk = 0, nx
    ! total contributions
    prof_eos_f(kk)       =  eos_ff(phi_tot(kk))
    prof_eos_rdfdr(kk)   = -phi_tot(kk) * eos_df_drho(phi_tot(kk))
    prof_field(kk)       = -phi_tot(kk) * wa(kk)
    prof_solid(kk)       =  phi_tot(kk) * Ufield(kk) / beta
    prof_sgt_f(kk)       =  0.5d0 * k_gr * (rho_seg_bulk * dphi_dr(kk))**2
    prof_sgt_rdfdr(kk)   =  k_gr * (rho_seg_bulk**2 * phi_tot(kk) * d2phi_dr2(kk))
    ! partial contributions
    prof_rho_wifc_glo(kk) = -phi_glo(kk)  * wa_ifc_new(kk)
    prof_rho_wifc_ghi(kk) = -phi_ghi(kk)  * wa_ifc_new(kk)
    prof_rho_wifc_mxa(kk)   = -phi_mxa(kk) * wa_ifc_new(kk)
    prof_rho_wifc_mxb(kk)   = -phi_mxb(kk) * wa_ifc_new(kk)
    prof_solid_glo(kk)   =  phi_glo(kk)  * Ufield(kk) / beta
    prof_solid_ghi(kk)   =  phi_ghi(kk)  * Ufield(kk) / beta
    prof_solid_mxa(kk)     =  phi_mxa(kk) * Ufield(kk) / beta
    prof_solid_mxb(kk)     =  phi_mxb(kk) * Ufield(kk) / beta
end do

if (exist_mxa.or.exist_mxb) then       
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
    E_rho_wifc_glo = E_rho_wifc_glo + coeff_nx(kk) * prof_rho_wifc_glo(kk) * layer_area(kk)
    E_rho_wifc_ghi = E_rho_wifc_ghi + coeff_nx(kk) * prof_rho_wifc_ghi(kk) * layer_area(kk)
    E_rho_wifc_mxa  = E_rho_wifc_mxa  + coeff_nx(kk) * prof_rho_wifc_mxa(kk)   * layer_area(kk)
    E_rho_wifc_mxb  = E_rho_wifc_mxb  + coeff_nx(kk) * prof_rho_wifc_mxb(kk)   * layer_area(kk)
    E_solid_glo   = E_solid_glo   + coeff_nx(kk) * prof_solid_glo(kk)   * layer_area(kk)
    E_solid_ghi   = E_solid_ghi   + coeff_nx(kk) * prof_solid_ghi(kk)   * layer_area(kk)
    E_solid_mxa     = E_solid_mxa     + coeff_nx(kk) * prof_solid_mxa(kk)     * layer_area(kk)
    E_solid_mxb     = E_solid_mxb     + coeff_nx(kk) * prof_solid_mxb(kk)     * layer_area(kk)
end do

! total contributions
E_eos_f       = E_eos_f     * 1.d-30
E_eos_rdfdr   = E_eos_rdfdr * 1.d-30 * rho_seg_bulk
E_field       = E_field     * 1.d-30 * rho_seg_bulk / beta
E_solid       = E_solid     * 1.d-30 * rho_seg_bulk
E_sgt_f       = E_sgt_f     * 1.d-30
E_sgt_rdfdr   = E_sgt_rdfdr * 1.d-30
! partial contributions
E_rho_wifc_glo = E_rho_wifc_glo * 1.d-30 * rho_seg_bulk / beta
E_rho_wifc_ghi = E_rho_wifc_ghi * 1.d-30 * rho_seg_bulk / beta
E_rho_wifc_mxa  = E_rho_wifc_mxa  * 1.d-30 * rho_seg_bulk / beta
E_rho_wifc_mxb  = E_rho_wifc_mxb  * 1.d-30 * rho_seg_bulk / beta
E_solid_glo   = E_solid_glo   * 1.d-30 * rho_seg_bulk
E_solid_ghi   = E_solid_ghi   * 1.d-30 * rho_seg_bulk
E_solid_mxa     = E_solid_mxa     * 1.d-30 * rho_seg_bulk
E_solid_mxb     = E_solid_mxb     * 1.d-30 * rho_seg_bulk

if (exist_mxa) then
   part_func_mxa = get_part_func(nx, ns_mxa, layer_area, volume, coeff_nx, qfinal_mxa)
   E_rhoVkTQ_mxa = volume*1.0d-30*rho_seg_bulk/beta/chainlen_mxa*(1.d0-part_func_mxa)
endif
if (exist_mxb) then
   part_func_mxb = get_part_func(nx, ns_mxb, layer_area, volume, coeff_nx, qfinal_mxb)
   E_rhoVkTQ_mxb = volume*1.0d-30*rho_seg_bulk/beta/chainlen_mxb*(1.d0-part_func_mxb)
endif

E_nkTlnQm = 0.d0
if (exist_glo) then
   nchglo  = get_nchains(coeff_nx, nx, layer_area, phi_glo, rho_seg_bulk, chainlen_glo)
   E_nkTlnQm = E_nkTlnQm -nchglo / beta * log(qfinal_glo_aux(gnode_lo,ns_glo))
endif
if (exist_ghi) then
   nchghi  = get_nchains(coeff_nx, nx, layer_area, phi_ghi, rho_seg_bulk, chainlen_ghi)
   E_nkTlnQm = E_nkTlnQm -nchghi / beta * log(qfinal_ghi_aux(gnode_hi,ns_ghi))
endif

if (wall_side.eq.F_both.and.wall_hamaker) then
   E_solid_solid = (surface_area * 1.e-20) * (Asolid/PI) *                              &
&                  ( + (sig_solid*1.e-10)**6 / (360.d0*((lx+2.d0*wall_pos)*1.0d-10)**8) &
&                    -         1.d0          / (12.d0 *((lx+2.d0*wall_pos)*1.0d-10)**2) )
endif
!
! estimate the contribution due to chain stretching
!
if (exist_glo) then
   !calculate the profile of the chain ends
   do kk = 0, nx
       phi_end(kk) = 1.d0 / chainlen_glo * (qfinal_glo(kk,ns_glo) * qfinal_glo_aux(kk,0))
   enddo
   rho_end = phi_end * rho_seg_bulk*1.d-30
   ! calculate the stretching free energy for every chain end abstained by dz from the grafting point
   do kk = 0, nx
       dz = rx(kk)-rx(gnode_lo)
       A_stretch(kk) = 3.0/(2.d0 * beta * (Rg2_per_mon_glo*chainlen_glo*6.d0)) * dz**2
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
if (exist_ghi) then
   !calculate the profile of the chain ends
   do kk = 0, nx
       phi_end(kk) = 1.d0 / chainlen_ghi * (qfinal_ghi(kk,ns_ghi) * qfinal_ghi_aux(kk,0))
   enddo
   rho_end = phi_end * rho_seg_bulk*1.d-30
   ! calculate the stretching free energy for every chain end abstained by dz from the grafting point
   do kk = 0, nx
       dz = rx(kk)-rx(gnode_hi)
       A_stretch(kk) = 3.0/(2.d0 * beta * (Rg2_per_mon_ghi*chainlen_ghi*6.d0)) * dz**2
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
E_rhoVkTQ_mxa   = E_rhoVkTQ_mxa   / (surface_area * 1.e-20) * 1e+3
E_rhoVkTQ_mxb   = E_rhoVkTQ_mxb   / (surface_area * 1.e-20) * 1e+3
E_nkTlnQm   = E_nkTlnQm   / (surface_area * 1.e-20) * 1e+3
E_sgt_f     = E_sgt_f     / (surface_area * 1.e-20) * 1e+3
E_sgt_rdfdr = E_sgt_rdfdr / (surface_area * 1.e-20) * 1e+3
! partial contributions
E_rho_wifc_glo     = E_rho_wifc_glo     / (surface_area * 1.e-20) * 1e+3
E_rho_wifc_ghi     = E_rho_wifc_ghi     / (surface_area * 1.e-20) * 1e+3
E_rho_wifc_mxa      = E_rho_wifc_mxa      / (surface_area * 1.e-20) * 1e+3
E_rho_wifc_mxb      = E_rho_wifc_mxb      / (surface_area * 1.e-20) * 1e+3
E_solid_glo       = E_solid_glo       / (surface_area * 1.e-20) * 1e+3
E_solid_ghi       = E_solid_ghi       / (surface_area * 1.e-20) * 1e+3
E_solid_mxa         = E_solid_mxa         / (surface_area * 1.e-20) * 1e+3
E_solid_mxb         = E_solid_mxb         / (surface_area * 1.e-20) * 1e+3
E_stretch_glo     = E_stretch_glo     / (surface_area * 1.e-20) * 1e+3
E_stretch_ghi     = E_stretch_ghi     / (surface_area * 1.e-20) * 1e+3
E_stretch_add_glo = E_stretch_add_glo / (surface_area * 1.e-20) * 1e+3
E_stretch_add_ghi = E_stretch_add_ghi / (surface_area * 1.e-20) * 1e+3
E_solid_solid     = E_solid_solid     / (surface_area * 1.e-20) * 1e+3
!estimate the free energy in mJ/m^2
free_energy = E_eos_f + E_sgt_f + E_eos_rdfdr + E_sgt_rdfdr + E_rhoVkTQ_mxa + E_rhoVkTQ_mxb + E_nkTlnQm + E_solid_solid

open(unit=777, file = "o.energies")
write(777,'(23(A16))')  "eos_f", "eos_rdfdr", "sgt_f", "sgt_rdfdr", "rhoVkTQmxa", "rhoVkTQmxb", "nkTlnQm_ns",       &
&                       "field", "solid", "free_energy", "E_rho_wifc_glo", "E_rho_wifc_ghi",         &
&                       "E_rho_wifc_mxa", "E_rho_wifc_mxb", "E_str_glo", "E_str_add_glo", "E_str_ghi", "E_str_add_ghi",&
&                       "solid_glo", "solid_ghi", "solid_mxa", "solid_mxb", "solid_solid"

write(777,'(23(E16.7))') E_eos_f, E_eos_rdfdr, E_sgt_f, E_sgt_rdfdr, E_rhoVkTQ_mxa, E_rhoVkTQ_mxb, E_nkTlnQm, E_field,       &
&                        E_solid, free_energy, E_rho_wifc_glo, E_rho_wifc_ghi, E_rho_wifc_mxa, E_rho_wifc_mxb, E_stretch_glo, &
&                        E_stretch_add_glo, E_stretch_ghi, E_stretch_add_ghi, E_solid_glo, E_solid_ghi,   &
&                        E_solid_mxa, E_solid_mxb, E_solid_solid
close(777)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_energies
