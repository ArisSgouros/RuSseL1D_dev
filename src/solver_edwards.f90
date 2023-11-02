!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine solver_edwards(bc_lo, bc_hi, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, Rg2_per_mon, &
                 & nx, ns, dx, ds, edwards_solver, linear_solver, wa, q, q_final)
!----------------------------------------------------------------------------------------------------------!
  use flags, only: F_semi_implicit, F_implicit, F_bc_neuman, F_bc_periodic, F_tridag, F_gelim
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  integer, intent(in)                    :: bc_lo, bc_hi, n_dir_nodes, nx, ns, edwards_solver, linear_solver
  integer, intent(in), dimension(0:nx + 4) :: dir_nodes_id
  integer                                :: tt, ii, id

  real(8), intent(in)                          :: Rg2_per_mon
  real(8), intent(in), dimension(0:nx + 4)       :: dir_nodes_rdiag
  real(8), intent(in), dimension(0:nx)         :: dx
  real(8), intent(in), dimension(0:nx)         :: wa
  real(8), intent(in), dimension(0:ns)         :: ds
  real(8), intent(inout), dimension(0:nx, 1:2)  :: q
  real(8), intent(inout), dimension(0:nx, 0:ns) :: q_final
  real(8), dimension(0:nx)                     :: adiag, bdiag, cdiag, rdiag, u, diff_number
  real(8), dimension(0:nx, 0:nx)                :: stiff
  real(8)                                      :: temp1
!----------------------------------------------------------------------------------------------------------!
  do tt = 1, ns
    do ii = 1, nx
      diff_number(ii) = 0.5d0*Rg2_per_mon*ds(tt)/dx(ii)**2
    end do
    diff_number(0) = 0.5d0*Rg2_per_mon*ds(tt)/dx(1)**2  !because dx(0) is not defined

    if (edwards_solver .eq. F_semi_implicit) then
      !interior nodes
      do ii = 1, nx - 1
        adiag(ii) = -diff_number(ii)
        cdiag(ii) = -diff_number(ii)
        bdiag(ii) = 1.d0 + 2.d0*diff_number(ii) + 0.5d0*ds(tt)*wa(ii)

        temp1 = (q(ii + 1, 1) - 2.d0*q(ii, 1) + q(ii - 1, 1))*diff_number(ii)
        rdiag(ii) = temp1 + (1 - 0.5d0*ds(tt)*wa(ii))*q(ii, 1)
      end do
      !Neumann boundary conditions
      if (bc_lo .eq. F_bc_neuman) then
        cdiag(0) = -diff_number(0)
        bdiag(0) = 1.d0 + 1.d0*diff_number(0) + 0.5d0*ds(tt)*wa(0)
        temp1 = (-1.d0*q(0, 1) + q(1, 1))*diff_number(0)
        rdiag(0) = temp1 + q(0, 1) - 0.5*ds(tt)*wa(0)*q(0, 1)
      end if
      if (bc_hi .eq. F_bc_neuman) then
        adiag(nx) = -diff_number(nx)
        bdiag(nx) = 1.d0 + 1.d0*diff_number(nx) + 0.5d0*ds(tt)*wa(nx)
        temp1 = (-1.d0*q(nx, 1) + q(nx - 1, 1))*diff_number(nx)
        rdiag(nx) = temp1 + q(nx, 1) - 0.5*ds(tt)*wa(nx)*q(nx, 1)
      end if
      !periodic boundary conditions
      if (bc_lo .eq. F_bc_periodic) then
        adiag(0) = -diff_number(0)
        cdiag(0) = -diff_number(0)
        bdiag(0) = 1.d0 + 2.d0*diff_number(0) + 0.5d0*ds(tt)*wa(0)

        temp1 = (q(0 + 1, 1) - 2.d0*q(0, 1) + q(nx, 1))*diff_number(0)
        rdiag(0) = temp1 + (1 - 0.5d0*ds(tt)*wa(0))*q(0, 1)
      end if
      if (bc_hi .eq. F_bc_periodic) then
        adiag(nx) = -diff_number(nx)
        cdiag(nx) = -diff_number(nx)
        bdiag(nx) = 1.d0 + 2.d0*diff_number(nx) + 0.5d0*ds(tt)*wa(nx)

        temp1 = (q(0, 1) - 2.d0*q(nx, 1) + q(nx - 1, 1))*diff_number(nx)
        rdiag(nx) = temp1 + (1 - 0.5d0*ds(tt)*wa(nx))*q(nx, 1)
      end if
    elseif (edwards_solver .eq. F_implicit) then
      !interior nodes
      do ii = 1, nx - 1
        adiag(ii) = -2.d0*diff_number(ii)
        cdiag(ii) = -2.d0*diff_number(ii)
        bdiag(ii) = 1.d0 + 4.d0*diff_number(ii) + ds(tt)*wa(ii)
        rdiag(ii) = q(ii, 1)
      end do
      !Neumann boundary conditions
      if (bc_lo .eq. F_bc_neuman) then
        cdiag(0) = -4.d0*diff_number(0)
        bdiag(0) = 1.d0 + 4.d0*diff_number(0) + ds(tt)*wa(0)
        rdiag(0) = q(0, 1)
      end if
      if (bc_hi .eq. F_bc_neuman) then
        adiag(nx) = -4.d0*diff_number(nx)
        bdiag(nx) = 1.d0 + 4.d0*diff_number(nx) + ds(tt)*wa(nx)
        rdiag(nx) = q(nx, 1)
      end if
      !periodic boundary conditions
      if (bc_lo .eq. F_bc_periodic) then
        adiag(0) = -2.d0*diff_number(0)
        cdiag(0) = -2.d0*diff_number(0)
        bdiag(0) = 1.d0 + 4.d0*diff_number(0) + ds(tt)*wa(0)
        rdiag(0) = q(0, 1)
      end if
      if (bc_hi .eq. F_bc_periodic) then
        adiag(nx) = -2.d0*diff_number(nx)
        cdiag(nx) = -2.d0*diff_number(nx)
        bdiag(nx) = 1.d0 + 4.d0*diff_number(nx) + ds(tt)*wa(nx)
        rdiag(nx) = q(nx, 1)
      end if
    end if

    !Dirichlet boundary conditions
    do ii = 0, n_dir_nodes - 1
      id = dir_nodes_id(ii)
      cdiag(id) = 0.0d0
      adiag(id) = 0.0d0
      bdiag(id) = 1.0d0
      rdiag(id) = dir_nodes_rdiag(ii)
    end do

    if (linear_solver .eq. F_tridag) then
      call solver_tridag(nx + 1, adiag, bdiag, cdiag, rdiag, u)
    elseif (linear_solver .eq. F_gelim) then
      !assembly stiffness matrix for gelim
      stiff = 0.d0
      do ii = 1, nx - 1
        stiff(ii, ii - 1) = adiag(ii)
        stiff(ii, ii) = bdiag(ii)
        stiff(ii, ii + 1) = cdiag(ii)
      end do
      stiff(0, 0) = bdiag(0)
      stiff(0, 1) = cdiag(0)
      stiff(nx, nx) = bdiag(nx)
      stiff(nx, nx - 1) = adiag(nx)

      if ((bc_lo .eq. F_bc_periodic) .and. (bc_hi .eq. F_bc_periodic)) then
        stiff(0, nx) = adiag(0)
        stiff(nx, 0) = cdiag(nx)
      end if

      !do ii = 0,nx
      !   write(123,'(51(E16.4E3,2X))') (stiff(ii,jj),jj=0,nx)
      !   write(456,'(2(E16.4E3,2X))') u(ii), rdiag(ii)
      !enddo

      call solver_gelim(nx + 1, stiff, rdiag, u)
    end if

    do ii = 0, nx
      q(ii, 2) = u(ii)
    end do

    do ii = 0, nx
      q_final(ii, tt) = q(ii, 2)
      q(ii, 1) = q(ii, 2)
    end do
  end do

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine solver_edwards
