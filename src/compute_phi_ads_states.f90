!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_phi_ads_states(coeff, rr, rx, dx, ds, wa, phi_mx, qfinal_mx, &
&                                 bc_lo_type, bc_hi_type, geometry, r_ads_lo, r_ads_hi, &
&                                 chainlen, Rg2_per_mon, ns, chain_type)
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: nx, edwards_solver, linear_solver, lx, wall_pos
  use flags, only: F_bc_dirichlet_eq_0, F_bc_dirichlet_eq_1, F_sphere
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  real(8), intent(in)        :: r_ads_lo, r_ads_hi, chainlen, Rg2_per_mon
  integer, dimension(0:nx + 4) :: dir_nodes_id
  integer                    :: tt, kk, ii, counter, n_dir_nodes

  integer, intent(in)        :: bc_lo_type, bc_hi_type, geometry, ns

  real(8), intent(in), dimension(0:nx)      :: dx, rx
  real(8), intent(in), dimension(0:ns)      :: ds
  real(8), intent(in), dimension(0:ns)      :: coeff
  real(8), intent(in), dimension(0:nx)      :: wa, phi_mx, rr
  real(8), intent(in), dimension(0:nx, 0:ns) :: qfinal_mx
  real(8), dimension(0:nx + 4)                :: dir_nodes_rdiag
  real(8), dimension(0:nx, 0:ns)             :: qfinal_not_alo, qfinal_full_alo, qfinal_alo_part, qf_full_final, &
  &                                            qfinal_not_ahi, qfinal_full_ahi, qfinal_ahi_part
  real(8), dimension(0:nx)                  :: phi_not_alo, phi_alo, phi_alo_full, phi_alo_part, phi_loop_flo,   &
  &                                            phi_tail_flo, phi_loop_alo, phi_tail_alo, phi_not_ahi, phi_ahi,   &
  &                                            phi_ahi_full, phi_ahi_part, phi_loop_fhi, phi_tail_fhi,           &
  &                                            phi_loop_ahi, phi_tail_ahi, phi_free, phi_bridge
  real(8), dimension(0:nx, 1:2)              :: q_part
  real(8)                                   :: distance

  character(6), intent(in) :: chain_type
  character(40)            :: filename = ''
!----------------------------------------------------------------------------------------------------------!
  phi_free = 0.d0
  phi_bridge = 0.d0
  phi_not_alo = 0.d0
  phi_alo = 0.d0
  phi_alo_full = 0.d0
  phi_alo_part = 0.d0
  phi_tail_flo = 0.d0
  phi_loop_flo = 0.d0
  phi_tail_alo = 0.d0
  phi_loop_alo = 0.d0

  phi_not_ahi = 0.d0
  phi_ahi = 0.d0
  phi_ahi_full = 0.d0
  phi_ahi_part = 0.d0
  phi_tail_fhi = 0.d0
  phi_loop_fhi = 0.d0
  phi_tail_ahi = 0.d0
  phi_loop_ahi = 0.d0

  qf_full_final = 0.d0
  qfinal_not_alo = 0.d0
  qfinal_full_alo = 0.d0
  qfinal_alo_part = 0.d0
  qfinal_not_ahi = 0.d0
  qfinal_full_ahi = 0.d0
  qfinal_ahi_part = 0.d0

!
! Partition function of chains that are fully adsorbed to the lo boundary
!
!set the dirichlet boundary conditions
  q_part = 0.d0
  n_dir_nodes = 0
  dir_nodes_id = 0
  dir_nodes_rdiag = 0.d0

!dirichlet lower bound
  if (bc_lo_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_lo_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    n_dir_nodes = n_dir_nodes + 1
  end if
!dirichlet upper bound
  if (bc_hi_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_hi_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    dir_nodes_rdiag(n_dir_nodes) = dir_nodes_rdiag(n_dir_nodes)
    n_dir_nodes = n_dir_nodes + 1
  end if

!dirichlet BC on extra nodes for chainshape
  do ii = 0, nx
    distance = rx(ii)
    if (distance .ge. r_ads_lo) then
      counter = ii
      q_part(ii, 1) = 0.d00
      qfinal_full_alo(ii, 0) = 0.d00
      dir_nodes_id(n_dir_nodes) = ii
      dir_nodes_rdiag(n_dir_nodes) = 0.d0
      n_dir_nodes = n_dir_nodes + 1
    else
      q_part(ii, 1) = 1.d00
      qfinal_full_alo(ii, 0) = 1.d00
    end if
  end do

  if (geometry .eq. F_sphere) then
    do ii = 0, nx
      q_part(ii, 1) = q_part(ii, 1)*rr(ii)
      qfinal_full_alo(ii, 0) = qfinal_full_alo(ii, 0)*rr(ii)
    end do
    do ii = 0, n_dir_nodes - 1
      dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
    end do
  end if

  call solver_edwards(bc_lo_type, bc_hi_type, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
             & Rg2_per_mon, nx, ns, dx, ds, edwards_solver, linear_solver, wa, q_part, qfinal_full_alo)

  if (geometry .eq. F_sphere) then
    do tt = 0, ns
      do ii = 0, nx
        qfinal_full_alo(ii, tt) = qfinal_full_alo(ii, tt)/rr(ii)
      end do
    end do
  end if

!
! Partition function of chains that are not adsorbed to the lo boundary
!
!set the dirichlet boundary conditions
  q_part = 0.d0
  n_dir_nodes = 0
  dir_nodes_id = 0
  dir_nodes_rdiag = 0.d0

!dirichlet lower bound
  if (bc_lo_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_lo_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    n_dir_nodes = n_dir_nodes + 1
  end if
!dirichlet upper bound
  if (bc_hi_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_hi_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    dir_nodes_rdiag(n_dir_nodes) = dir_nodes_rdiag(n_dir_nodes)
    n_dir_nodes = n_dir_nodes + 1
  end if

!dirichlet BC on extra nodes for chainshape
  do ii = 0, nx
    distance = rx(ii)
    if (distance .lt. r_ads_lo) then
      counter = ii
      q_part(ii, 1) = 0.d00
      qfinal_not_alo(ii, 0) = 0.d00
      dir_nodes_id(n_dir_nodes) = ii
      dir_nodes_rdiag(n_dir_nodes) = 0.d0
      n_dir_nodes = n_dir_nodes + 1
    else
      q_part(ii, 1) = 1.d00
      qfinal_not_alo(ii, 0) = 1.d00
    end if
  end do

  if (geometry .eq. F_sphere) then
    do ii = 0, nx
      q_part(ii, 1) = q_part(ii, 1)*rr(ii)
      qfinal_not_alo(ii, 0) = qfinal_not_alo(ii, 0)*rr(ii)
    end do
    do ii = 0, n_dir_nodes - 1
      dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
    end do
  end if

  call solver_edwards(bc_lo_type, bc_hi_type, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
             & Rg2_per_mon, nx, ns, dx, ds, edwards_solver, linear_solver, wa, q_part, qfinal_not_alo)

  if (geometry .eq. F_sphere) then
    do tt = 0, ns
      do ii = 0, nx
        qfinal_not_alo(ii, tt) = qfinal_not_alo(ii, tt)/rr(ii)
      end do
    end do
  end if

!
! Partition function of free chains (not adsorbed either to the lo or hi boundary)
!
!set the dirichlet boundary conditions
  q_part = 0.d0
  n_dir_nodes = 0
  dir_nodes_id = 0
  dir_nodes_rdiag = 0.d0

!dirichlet lower bound
  if (bc_lo_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_lo_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    n_dir_nodes = n_dir_nodes + 1
  end if
!dirichlet upper bound
  if (bc_hi_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_hi_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    dir_nodes_rdiag(n_dir_nodes) = dir_nodes_rdiag(n_dir_nodes)
    n_dir_nodes = n_dir_nodes + 1
  end if

!dirichlet BC on extra nodes for chainshape
  do ii = 0, nx
    if (rx(ii) .lt. r_ads_lo .or. lx + 2.0*wall_pos - rx(ii) .lt. r_ads_hi) then
      counter = ii
      q_part(ii, 1) = 0.d00
      qf_full_final(ii, 0) = 0.d00
      dir_nodes_id(n_dir_nodes) = ii
      dir_nodes_rdiag(n_dir_nodes) = 0.d0
      n_dir_nodes = n_dir_nodes + 1
    else
      q_part(ii, 1) = 1.d00
      qf_full_final(ii, 0) = 1.d00
    end if
  end do

  if (geometry .eq. F_sphere) then
    do ii = 0, nx
      q_part(ii, 1) = q_part(ii, 1)*rr(ii)
      qf_full_final(ii, 0) = qf_full_final(ii, 0)*rr(ii)
    end do
    do ii = 0, n_dir_nodes - 1
      dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
    end do
  end if

  call solver_edwards(bc_lo_type, bc_hi_type, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
             & Rg2_per_mon, nx, ns, dx, ds, edwards_solver, linear_solver, wa, q_part, qf_full_final)

  if (geometry .eq. F_sphere) then
    do tt = 0, ns
      do ii = 0, nx
        qf_full_final(ii, tt) = qf_full_final(ii, tt)/rr(ii)
      end do
    end do
  end if

!
! Partition function of chains that are not adsorbed to the hi boundary
!
!set the dirichlet boundary conditions
  q_part = 0.d0
  n_dir_nodes = 0
  dir_nodes_id = 0
  dir_nodes_rdiag = 0.d0

!dirichlet lower bound
  if (bc_lo_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_lo_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    n_dir_nodes = n_dir_nodes + 1
  end if
!dirichlet upper bound
  if (bc_hi_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_hi_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    dir_nodes_rdiag(n_dir_nodes) = dir_nodes_rdiag(n_dir_nodes)
    n_dir_nodes = n_dir_nodes + 1
  end if

!dirichlet BC on extra nodes for chainshape
  do ii = 0, nx
    distance = lx + 2.0*wall_pos - rx(ii) !! hard sphere wall?!?
    if (distance .lt. r_ads_hi) then
      counter = ii
      q_part(ii, 1) = 0.d00
      qfinal_not_ahi(ii, 0) = 0.d00
      dir_nodes_id(n_dir_nodes) = ii
      dir_nodes_rdiag(n_dir_nodes) = 0.d0
      n_dir_nodes = n_dir_nodes + 1
    else
      q_part(ii, 1) = 1.d00
      qfinal_not_ahi(ii, 0) = 1.d00
    end if
  end do

  if (geometry .eq. F_sphere) then
    do ii = 0, nx
      q_part(ii, 1) = q_part(ii, 1)*rr(ii)
      qfinal_not_ahi(ii, 0) = qfinal_not_ahi(ii, 0)*rr(ii)
    end do
    do ii = 0, n_dir_nodes - 1
      dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
    end do
  end if

  call solver_edwards(bc_lo_type, bc_hi_type, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
             & Rg2_per_mon, nx, ns, dx, ds, edwards_solver, linear_solver, wa, q_part, qfinal_not_ahi)

  if (geometry .eq. F_sphere) then
    do tt = 0, ns
      do ii = 0, nx
        qfinal_not_ahi(ii, tt) = qfinal_not_ahi(ii, tt)/rr(ii)
      end do
    end do
  end if

!
! Partition function of chains that are fully adsorbed to the hi boundary
!
!set the dirichlet boundary conditions
  q_part = 0.d0
  n_dir_nodes = 0
  dir_nodes_id = 0
  dir_nodes_rdiag = 0.d0

!dirichlet lower bound
  if (bc_lo_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_lo_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    n_dir_nodes = n_dir_nodes + 1
  end if
!dirichlet upper bound
  if (bc_hi_type .eq. F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1
  else if (bc_hi_type .eq. F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    dir_nodes_rdiag(n_dir_nodes) = dir_nodes_rdiag(n_dir_nodes)
    n_dir_nodes = n_dir_nodes + 1
  end if

!dirichlet BC on extra nodes for chainshape
  do ii = 0, nx
    distance = lx + 2.0*wall_pos - rx(ii)
    if (distance .ge. r_ads_hi) then
      counter = ii
      q_part(ii, 1) = 0.d00
      qfinal_full_ahi(ii, 0) = 0.d00
      dir_nodes_id(n_dir_nodes) = ii
      dir_nodes_rdiag(n_dir_nodes) = 0.d0
      n_dir_nodes = n_dir_nodes + 1
    else
      q_part(ii, 1) = 1.d00
      qfinal_full_ahi(ii, 0) = 1.d00
    end if
  end do

  if (geometry .eq. F_sphere) then
    do ii = 0, nx
      q_part(ii, 1) = q_part(ii, 1)*rr(ii)
      qfinal_full_ahi(ii, 0) = qfinal_full_ahi(ii, 0)*rr(ii)
    end do
    do ii = 0, n_dir_nodes - 1
      dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
    end do
  end if

  call solver_edwards(bc_lo_type, bc_hi_type, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
             & Rg2_per_mon, nx, ns, dx, ds, edwards_solver, linear_solver, wa, q_part, qfinal_full_ahi)

  if (geometry .eq. F_sphere) then
    do tt = 0, ns
      do ii = 0, nx
        qfinal_full_ahi(ii, tt) = qfinal_full_ahi(ii, tt)/rr(ii)
      end do
    end do
  end if

  qfinal_alo_part = qfinal_mx - qfinal_not_alo - qfinal_full_alo
  qfinal_ahi_part = qfinal_mx - qfinal_not_ahi - qfinal_full_ahi

  call contour_convolution(chainlen, nx, ns, coeff, qf_full_final, qf_full_final, phi_free)
  call contour_convolution(chainlen, nx, ns, coeff, qfinal_full_alo, qfinal_full_alo, phi_alo_full)
  call contour_convolution(chainlen, nx, ns, coeff, qfinal_full_ahi, qfinal_full_ahi, phi_ahi_full)
  call contour_convolution(chainlen, nx, ns, coeff, qfinal_not_alo, qfinal_not_alo, phi_not_alo)
  call contour_convolution(chainlen, nx, ns, coeff, qfinal_not_ahi, qfinal_not_ahi, phi_not_ahi)

  call contour_convolution(chainlen, nx, ns, coeff, qfinal_not_alo, qfinal_alo_part, phi_tail_flo)
  call contour_convolution(chainlen, nx, ns, coeff, qfinal_full_alo, qfinal_alo_part, phi_tail_alo)
  call contour_convolution(chainlen, nx, ns, coeff, qfinal_not_ahi, qfinal_ahi_part, phi_tail_fhi)
  call contour_convolution(chainlen, nx, ns, coeff, qfinal_full_ahi, qfinal_ahi_part, phi_tail_ahi)

  call contour_convolution(chainlen, nx, ns, coeff, qfinal_alo_part, qfinal_alo_part, phi_loop_flo)
  call contour_convolution(chainlen, nx, ns, coeff, qfinal_ahi_part, qfinal_ahi_part, phi_loop_fhi)

  phi_alo = phi_mx - phi_not_alo
  phi_ahi = phi_mx - phi_not_ahi
  phi_tail_flo = 2.d0*phi_tail_flo
  phi_tail_fhi = 2.d0*phi_tail_fhi
  phi_tail_alo = 2.d0*phi_tail_alo
  phi_tail_ahi = 2.d0*phi_tail_ahi
  phi_alo_part = phi_alo - phi_alo_full
  phi_ahi_part = phi_ahi - phi_ahi_full
  phi_bridge = max(phi_alo + phi_ahi - phi_mx + phi_free, 0.d0)

  write (filename, '("o.phi_states_",A6)') chain_type
  open (unit=37, file=filename)
  write (37, '(A14,17(2X,A20))') 'r', "phi_m", "f", &
  &                                     "!a+", "a-", "afull-", &
  &                                     "apart-", "loop_f-", "tail_f-", &
  &                                     "tail_a-", "!a-", "a+", &
  &                                     "afull+", "apart+", "loop_f+", &
  &                                     "tail_f+", "tail_a+", "bridge"
  do kk = 0, nx
    write (37, '(F14.7,17(2X,F20.11))') rx(kk), phi_mx(kk), phi_free(kk),     &
&                                     phi_not_ahi(kk), phi_alo(kk), phi_alo_full(kk),     &
&                                     phi_alo_part(kk), phi_loop_flo(kk), phi_tail_flo(kk), &
&                                     phi_tail_alo(kk), phi_not_alo(kk), phi_ahi(kk),        &
&                                     phi_ahi_full(kk), phi_ahi_part(kk), phi_loop_fhi(kk),  &
&                                     phi_tail_fhi(kk), phi_tail_ahi(kk), phi_bridge(kk)
  end do
  close (37)

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_phi_ads_states
