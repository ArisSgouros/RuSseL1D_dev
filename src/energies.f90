subroutine energies(free_energy)
!----------------------------------------------------------------------------------------------------------!
use eos
use arrays
use parser_vars
use constants
use flags
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: kk

real(8)                                          :: nchgr_lo, nchgr_hi
real(8), intent(out)                             :: free_energy
real(8)                                          :: get_nchains, get_part_func

real(8)                                          :: E_eos_f, E_eos_rdfdr, E_rhoVkTQ, E_nkTlnQm, E_grad_A, E_grad_B
real(8), dimension(0:nx)                         :: prof_eos_f, prof_eos_rdfdr, prof_term_3, prof_term_4, prof_grad_A, prof_grad_B
! TEMPORARY VARIABLES
integer                                          :: scheme = 1
real(8)                                          :: rho_tilde, part_func_matrix
real(8)                                          :: term1, term2, term3, term4, term4b
real(8), dimension(0:nx)                         :: prof_term_1, prof_term_2, prof_term_4b


real(8) :: sum1, sum2

!----------------------------------------------------------------------------------------------------------!
E_eos_f     = 0.d0
E_eos_rdfdr = 0.d0
E_rhoVkTQ   = 0.d0
E_nkTlnQm   = 0.d0
E_grad_A    = 0.d0
E_grad_B    = 0.d0

prof_eos_f = 0.d0
prof_eos_rdfdr = 0.d0
prof_grad_A = 0.d0
prof_grad_B = 0.d0

if (scheme.eq.1) then

  do kk = 0, nx
      prof_eos_f(kk)     =  eos_ff(phi_total(kk))! - eos_ff(1.d0)
      prof_eos_rdfdr(kk) = - phi_total(kk) * eos_df_drho(phi_total(kk))! +1.d0*eos_df_drho(1.d0)
      prof_grad_A(kk)    = -0.5d0 * k_gr * (rho_seg_bulk * dphi_dr(kk))**2
      prof_grad_B(kk)    =  0.5d0 * k_gr * (rho_seg_bulk**2 * phi_total(kk) * d2phi_dr2(kk))
  end do

  if (matrix_exist) then       
      do kk = 0, nx
          prof_eos_f(kk)     = prof_eos_f(kk) - eos_ff(1.d0)
          prof_eos_rdfdr(kk) = prof_eos_rdfdr(kk) +1.d0*eos_df_drho(1.d0)
      end do
  endif

  do kk = 0, nx
      E_eos_f     = E_eos_f + coeff_nx(kk)* prof_eos_f(kk)*layer_area(kk)
      E_eos_rdfdr = E_eos_rdfdr + coeff_nx(kk)* prof_eos_rdfdr(kk)*layer_area(kk)
      E_grad_A      = E_grad_A + coeff_nx(kk)* prof_grad_A(kk)*layer_area(kk)
      E_grad_B      = E_grad_B + coeff_nx(kk)* prof_grad_B(kk)*layer_area(kk)
  end do

  E_eos_f     = E_eos_f * 1.d-30
  E_eos_rdfdr = E_eos_rdfdr * 1.d-30 * rho_seg_bulk
  E_grad_A      = E_grad_A * 1.d-30
  E_grad_B      = E_grad_B * 1.d-30

  if (matrix_exist) then       
     part_func_matrix = get_part_func(nx, ns_matrix, layer_area, volume, coeff_nx, qmatrix_final)
     E_rhoVkTQ = volume*1.0d-30*rho_seg_bulk*boltz_const_Joule_K*Temp/chainlen_matrix*(1.d0- part_func_matrix)
  endif

  nchgr_lo   = get_nchains(coeff_nx, nx, layer_area, phi_gr_lo, rho_seg_bulk, chainlen_grafted_lo)
  nchgr_hi   = get_nchains(coeff_nx, nx, layer_area, phi_gr_hi, rho_seg_bulk, chainlen_grafted_hi)
  E_nkTlnQm = 0.d0
  if (grafted_lo_exist) E_nkTlnQm = E_nkTlnQm -nchgr_lo * boltz_const_Joule_K * Temp * log(qmatrix_final(gnode_lo,ns_grafted_lo))
  if (grafted_hi_exist) E_nkTlnQm = E_nkTlnQm -nchgr_hi * boltz_const_Joule_K * Temp * log(qmatrix_final(gnode_hi,ns_grafted_hi))

  !cast the terms 1-4 in units mJ/m^2
  E_eos_f     = E_eos_f / (surface_area * 1.e-20) * 1e+3
  E_eos_rdfdr = E_eos_rdfdr / (surface_area * 1.e-20) * 1e+3
  E_rhoVkTQ   = E_rhoVkTQ / (surface_area * 1.e-20) * 1e+3
  E_nkTlnQm   = E_nkTlnQm / (surface_area * 1.e-20) * 1e+3
  E_grad_A    = E_grad_A / (surface_area * 1.e-20) * 1e+3
  E_grad_B    = E_grad_B / (surface_area * 1.e-20) * 1e+3

  !estimate the free energy in mJ/m^2
  free_energy = E_eos_f + E_eos_rdfdr + E_rhoVkTQ + E_nkTlnQm + E_grad_A

  open(unit=777, file = 'o.energies.out.txt')
  write(777,'(7(A16))')"eos_f", "eos_rho_df_drho", "rhoVkTQ", "nkTlnQm_ns", "grad_A", "grad_B", "free_energy"
  write(777,'(7(E16.7))')E_eos_f, E_eos_rdfdr, E_rhoVkTQ, E_nkTlnQm, E_grad_A, E_grad_B, free_energy

elseif (scheme.eq.2) then

  free_energy = 0.d0

  do kk = 0, nx
      rho_tilde = rho_tilde_bulk * phi_total(kk)

      prof_term_2(kk)  = rho_tilde**2
      prof_term_3(kk)  = T_tilde*log(1-rho_tilde)
      prof_term_4(kk)  = T_tilde * rho_tilde
      prof_term_4b(kk) = T_tilde * rho_tilde * (-1.d0/(rsl_N*chainlen_matrix))
      ! Approach A: gradient
      prof_grad_A(kk)    = -0.5d0 * k_gr / P_star * (rho_seg_bulk * dphi_dr(kk))**2
      ! Approach B: Laplacian
      prof_grad_B(kk)    = 0.5d0 * k_gr / P_star * (rho_seg_bulk**2 * phi_total(kk) * d2phi_dr2(kk))

      sum1 = sum1 + dphi_dr(kk)**2
      sum2 = sum2 + phi_total(kk) * d2phi_dr2(kk)
  end do

  do kk = 0, nx
      term2   = term2 + coeff_nx(kk)* prof_term_2(kk)*layer_area(kk)
      term3   = term3 + coeff_nx(kk)* prof_term_3(kk)*layer_area(kk)
      term4   = term4 + coeff_nx(kk)* prof_term_4(kk)*layer_area(kk)
      term4b  = term4b + coeff_nx(kk)* prof_term_4b(kk)*layer_area(kk)
      E_grad_A= E_grad_A + coeff_nx(kk)* prof_grad_A(kk)*layer_area(kk)
      E_grad_B= E_grad_B + coeff_nx(kk)* prof_grad_B(kk)*layer_area(kk)
  end do

  term2    = term2*1.d-30*P_star
  term3    = term3*1.d-30*P_star
  term4    = term4*1.d-30*P_star
  term4b   = term4b*1.d-30*P_star
  E_grad_A = E_grad_A*1.d-30*P_star
  E_grad_B = E_grad_B*1.d-30*P_star

  term2    = term2 / (surface_area * 1.e-20) * 1e+3
  term3    = term3 / (surface_area * 1.e-20) * 1e+3
  term4    = term4 / (surface_area * 1.e-20) * 1e+3
  term4b   = term4b / (surface_area * 1.e-20) * 1e+3
  E_grad_A = E_grad_A / (surface_area * 1.e-20) * 1e+3
  E_grad_B = E_grad_B / (surface_area * 1.e-20) * 1e+3

  if (.not.matrix_exist) term4b = 0.d0

  nchgr_lo   = get_nchains(coeff_nx, nx, layer_area, phi_gr_lo, rho_seg_bulk, chainlen_grafted_lo)
  nchgr_hi   = get_nchains(coeff_nx, nx, layer_area, phi_gr_hi, rho_seg_bulk, chainlen_grafted_hi)
  E_nkTlnQm = 0.d0
  if (grafted_lo_exist) E_nkTlnQm = E_nkTlnQm -nchgr_lo * boltz_const_Joule_K * Temp * log(qmatrix_final(gnode_lo,ns_grafted_lo))
  if (grafted_hi_exist) E_nkTlnQm = E_nkTlnQm -nchgr_hi * boltz_const_Joule_K * Temp * log(qmatrix_final(gnode_hi,ns_grafted_hi))
  E_nkTlnQm   = E_nkTlnQm / (surface_area * 1.e-20) * 1e+3

  free_energy = term2 + term3 + term4 + term4b + E_grad_A + E_nkTlnQm

open(unit=777, file = 'o.energies.out.txt')
write(777,'(7(A16))') "term2", "term3", "term4", "term4b", "E_grad_A", "E_grad_B", "free_energy"
write(777,'(7(E16.7))') term2, term3, term4, term4b, E_grad_A, E_grad_B, free_energy

endif

! The section below is useful for debugging purposes
!write(777,*)
!write(777,'(3(A16))')'r','term1','term2'
!do kk = 0, nx
!    write(777,'(3(E16.7))')rx(kk), prof_eos_f(kk), prof_eos_rdfdr(kk)
!enddo
close(777)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine energies
