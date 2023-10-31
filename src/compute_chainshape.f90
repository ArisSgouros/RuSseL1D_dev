!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_chainshape(rho_seg_bulk, Rg2_per_mon, geometry, gnode, edwards_solver,       &
                            & linear_solver, bc_lower_type, bc_upper_type, qinit, coeff_x, rr,  &
                            & layer_area, rx, nx, dx, chainlen, ns, ds, wa, phi, q_final, chain_type)
!----------------------------------------------------------------------------------------------------------!
use flags, only: F_bc_dirichlet_eq_0, F_bc_dirichlet_eq_1, F_sphere
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in)      :: bc_lower_type, bc_upper_type, geometry, edwards_solver, linear_solver, ns, nx
integer, intent(in)      :: gnode
integer, dimension(0:nx) :: dir_nodes_id
integer                  :: tt, kk, ii, mm, n_dir_nodes

character(6), intent(in) :: chain_type
character(30)            :: filename = ''

real(8), intent(in), dimension(0:nx)      :: dx, rx
real(8), intent(in), dimension(0:ns)      :: ds
real(8), intent(in)                       :: qinit
real(8), intent(in)                       :: chainlen, Rg2_per_mon, rho_seg_bulk
real(8), intent(in), dimension(0:nx)      :: coeff_x, rr, layer_area, wa, phi
real(8), intent(in), dimension(0:nx,0:ns) :: q_final
real(8), dimension(0:nx)                  :: dir_nodes_rdiag
real(8), dimension(0:nx,0:ns)             :: qshape_final
real(8), dimension(0:nx,1:2)              :: qshape
real(8), dimension(0:nx)                  :: n_shape, p_cross
real(8)                                   :: sum_final, sum_q, sum_phi
!----------------------------------------------------------------------------------------------------------!
p_cross = 0.d0
n_shape = 0.d0

do kk = 1, nx-1

    if (chain_type.eq."matrixA") then
        do ii = 0, nx
            qshape(ii,1)       = 1.d0
            qshape_final(ii,0) = 1.d0
        enddo
    elseif (chain_type.eq."gra_lo".or.chain_type.eq."gra_hi") then
        do ii=0,nx
            qshape(ii,1)       = 0.d00
            qshape_final(ii,0) = 0.d00
        enddo
        qshape(gnode,1)       = qinit
        qshape_final(gnode,0) = qinit
    endif

    !set the dirichlet boundary conditions
    n_dir_nodes = 0
    dir_nodes_id = 0
    dir_nodes_rdiag = 0.d0

    !dirichlet lower bound
    if (bc_lower_type.eq.F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1   
    else if (bc_lower_type.eq.F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = 0
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1   
    endif
    !dirichlet upper bound
    if (bc_upper_type.eq.F_bc_dirichlet_eq_0) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 0.d0
        n_dir_nodes = n_dir_nodes + 1   
    else if (bc_upper_type.eq.F_bc_dirichlet_eq_1) then
        dir_nodes_id(n_dir_nodes) = nx
        dir_nodes_rdiag(n_dir_nodes) = 1.0d0
        n_dir_nodes = n_dir_nodes + 1   
    endif

    !dirichlet BC on extra nodes for chainshape
    dir_nodes_id(n_dir_nodes) = kk
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1   

    if (geometry.eq.F_sphere) then
        do ii = 0, nx
            qshape(ii,1)       = qshape(ii,1) * rr(ii)
            qshape_final(ii,0) = qshape_final(ii,0) * rr(ii)
        enddo
        do ii = 0, n_dir_nodes-1
            dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
        enddo
    endif
    call solver_edwards(bc_lower_type, bc_upper_type, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, Rg2_per_mon, &
               & nx, ns, dx, ds, edwards_solver, linear_solver, wa, qshape, qshape_final)

    if (geometry.eq.F_sphere) then
        do tt = 0, ns
            do ii = 0, nx
                qshape_final(ii,tt) = qshape_final(ii,tt)/rr(ii)
           enddo
        enddo
    endif

    !integrate
    sum_final = 0.d0
    sum_Q= 0.d0
    do mm = 0, nx
        sum_final = sum_final + coeff_x(mm) * qshape_final(mm,ns) * layer_area(mm)
        sum_Q = sum_Q + coeff_x(mm) * q_final(mm,ns) * layer_area(mm)
    enddo

    p_cross(kk)  = 1.0 - sum_final / sum_Q
enddo!kk

sum_phi = 0.0d0
do kk = 0, nx
    sum_phi = sum_phi + phi(kk) * coeff_x(kk) * layer_area(kk)
enddo

do ii=0,nx
    n_shape(ii) = p_cross(ii) * rho_seg_bulk / chainlen * sum_phi * 1.e-30 /layer_area(ii)
enddo

write(filename,'("o.chain_shape_",A6)')chain_type
open(unit=38, file=filename)
write(38,'(A14,6(2X,A14))')'z', "z_Rg", "phi", "pcross", "phi/n_shape", "q1/p_cross", "n_shape"
do ii = 1, nx-1
    write(38,'(E14.7,6(2X,E14.7))') rx(ii), rx(ii)/dsqrt(Rg2_per_mon*chainlen), phi(ii), p_cross(ii), &
                                    phi(ii)/n_shape(ii), q_final(ii,ns)/p_cross(ii), n_shape(ii)
enddo
close(38)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_chainshape
