subroutine solver_edwards(bc_lo, bc_hi, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, Rg2_per_mon, &
                 & nx, ns, dx, ds, edwards_solver, wa, q, q_final)
!----------------------------------------------------------------------------------------------------------!
use flags, only: F_semi_implicit, F_implicit, F_bc_neuman
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in)                  :: bc_lo, bc_hi, n_dir_nodes, nx, ns, edwards_solver
integer, intent(in), dimension(0:nx) :: dir_nodes_id
integer                              :: tt, ii, id

real(8), intent(in)                          :: Rg2_per_mon
real(8), intent(in), dimension(0:nx)         :: dir_nodes_rdiag
real(8), intent(in), dimension(0:nx)         :: dx
real(8), intent(in), dimension(0:nx)         :: wa
real(8), intent(in), dimension(0:ns)         :: ds
real(8), intent(inout), dimension(0:nx,1:2)  :: q
real(8), intent(inout), dimension(0:nx,0:ns) :: q_final
real(8), dimension(0:nx)                     :: adiag, bdiag, cdiag, rdiag, u, diff_number
real(8)                                      :: temp1
!----------------------------------------------------------------------------------------------------------!
do tt = 1, ns
    do ii = 1, nx
        diff_number(ii)  = 0.5d0*Rg2_per_mon*ds(tt)/dx(ii)**2
    enddo

    if (edwards_solver.eq.F_semi_implicit) then
        do ii = 1, nx-1
            adiag(ii) = -diff_number(ii)
            cdiag(ii) = -diff_number(ii)
            bdiag(ii) = 1.d0 + 2.d0*diff_number(ii) + 0.5d0*ds(tt)*wa(ii)

            temp1     = (q(ii+1,1) - 2.d0*q(ii,1) + q(ii-1,1))*diff_number(ii)
            rdiag(ii) = temp1 + (1 - 0.5d0*ds(tt)*wa(ii)) * q(ii,1)
        enddo
        if (bc_lo.eq.F_bc_neuman) then
            adiag(0) = -diff_number(0)
            bdiag(0) = 1.d0 + 1.d0*diff_number(0) + 0.5d0*ds(tt)*wa(0)
            temp1     = (-1.d0*q(0,1) + q(1,1)) * diff_number(0)
            rdiag(0) = temp1 + q(0,1) - 0.5*ds(tt)*wa(0) * q(0,1)
        endif
        if (bc_hi.eq.F_bc_neuman) then
            adiag(nx) = -diff_number(nx)
            bdiag(nx) = 1.d0 + 1.d0*diff_number(nx) + 0.5d0*ds(tt)*wa(nx)
            temp1     = (-1.d0*q(nx,1) + q(nx-1,1)) * diff_number(nx)
            rdiag(nx) = temp1 + q(nx,1) - 0.5*ds(tt)*wa(nx) * q(nx,1)
        endif
    elseif (edwards_solver.eq.F_implicit) then
        do ii = 1, nx-1
            adiag(ii) = -2.d0*diff_number(ii)
            cdiag(ii) = -2.d0*diff_number(ii)
            bdiag(ii) = 1.d0 + 4.d0*diff_number(ii) + ds(tt)*wa(ii)
            rdiag(ii) = q(ii,1)
        enddo
        if (bc_lo.eq.F_bc_neuman) then
            adiag(0) = -4.d0*diff_number(0)
            bdiag(0) = 1.d0 + 4.d0*diff_number(0) + ds(tt)*wa(0)
            rdiag(0) = q(0,1)
        endif
        if (bc_hi.eq.F_bc_neuman) then
            adiag(nx) = -4.d0*diff_number(nx)
            bdiag(nx) = 1.d0 + 4.d0*diff_number(nx) + ds(tt)*wa(nx)
            rdiag(nx) = q(nx,1)
        endif
    endif

    do ii = 0, n_dir_nodes-1
        id = dir_nodes_id(ii)
        cdiag(id) = 0.0d0
        adiag(id) = 0.0d0
        bdiag(id) = 1.0d0
        rdiag(id) = dir_nodes_rdiag(ii)
    enddo

    call solver_tridag(nx+1, adiag, bdiag, cdiag, rdiag, u)

    do ii = 0,nx
        q(ii,2) = u(ii)
    enddo

    do ii = 0, nx
        q_final(ii,tt) = q(ii,2)
        q(ii,1)        = q(ii,2)
    enddo
enddo

return
!----------------------------------------------------------------------------------------------------------!
end subroutine solver_edwards
