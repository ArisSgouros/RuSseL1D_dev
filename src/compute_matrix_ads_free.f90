subroutine compute_matrix_ads_free(coeff, rr, rx, dx, ds, wa, phi_matrix, qmatrix_final)
!----------------------------------------------------------------------------------------------------------!
use parser_vars, only: bc_lo_matrix, bc_hi_matrix, chainlen_matrix, ns_matrix, nx, edwards_solver,&
                     & geometry, Rg2_per_mon, r_critical
use flags, only: F_bc_dirichlet_eq_0, F_bc_dirichlet_eq_1, F_sphere
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, dimension(0:nx) :: dir_nodes_id
integer                  :: tt, kk, ii, counter, n_dir_nodes

real(8), intent(in), dimension(0:nx)             :: dx, rx
real(8), intent(in), dimension(0:ns_matrix)      :: ds
real(8), intent(in), dimension(0:ns_matrix)      :: coeff
real(8), intent(in), dimension(0:nx)             :: wa, phi_matrix, rr
real(8), intent(in), dimension(0:nx,0:ns_matrix) :: qmatrix_final
real(8), dimension(0:nx)                         :: dir_nodes_rdiag
real(8), dimension(0:nx,0:ns_matrix)             :: qfree_final
real(8), dimension(0:nx)                         :: phi_ads, phi_free, phi_loop, phi_tail
real(8), dimension(0:nx,0:ns_matrix)             :: qads_final
real(8), dimension(0:nx,1:2)                     :: qfree
real(8)                                          :: distance
!----------------------------------------------------------------------------------------------------------!
phi_free    = 0.d0
phi_ads     = 0.d0
phi_tail    = 0.d0
phi_loop    = 0.d0
qfree_final = 0.d0

do ii = 0, nx
    distance = rx(ii)
    if (distance.lt.r_critical) then
        counter = ii
        qfree(ii,1)       = 0.d00
        qfree_final(ii,0) = 0.d00
    else
        qfree(ii,1)       = 1.d00
        qfree_final(ii,0) = 1.d00
    endif
enddo

!set the dirichlet boundary conditions
n_dir_nodes = 0
dir_nodes_id = 0
dir_nodes_rdiag = 0.d0

!dirichlet lower bound
if (bc_lo_matrix.eq.F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1   
else if (bc_lo_matrix.eq.F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = 0
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    n_dir_nodes = n_dir_nodes + 1   
endif
!dirichlet upper bound
if (bc_hi_matrix.eq.F_bc_dirichlet_eq_0) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1   
else if (bc_hi_matrix.eq.F_bc_dirichlet_eq_1) then
    dir_nodes_id(n_dir_nodes) = nx
    dir_nodes_rdiag(n_dir_nodes) = 1.0d0
    dir_nodes_rdiag(n_dir_nodes) = dir_nodes_rdiag(n_dir_nodes)
    n_dir_nodes = n_dir_nodes + 1   
endif

!dirichlet BC on extra nodes for chainshape
do ii = 0, counter
    dir_nodes_id(n_dir_nodes) = ii
    dir_nodes_rdiag(n_dir_nodes) = 0.d0
    n_dir_nodes = n_dir_nodes + 1   
enddo

if (geometry.eq.F_sphere) then
    do ii = 0, nx
        qfree(ii,1)       = qfree(ii,1) * rr(ii)
        qfree_final(ii,0) = qfree_final(ii,0) * rr(ii)
    enddo
    do ii = 0, n_dir_nodes-1
        dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
    enddo
endif

call solver_edwards(bc_lo_matrix, bc_hi_matrix, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
           & Rg2_per_mon, nx, ns_matrix, dx, ds, edwards_solver, wa, qfree, qfree_final)

if (geometry.eq.F_sphere) then
    do tt = 0, ns_matrix
        do ii = 0, nx
            qfree_final(ii,tt) = qfree_final(ii,tt)/rr(ii)
        enddo
    enddo
endif

do tt = 0, ns_matrix
    do ii = 0, nx
        qads_final(ii,tt) = qmatrix_final(ii,tt) - qfree_final(ii,tt)
    enddo
enddo

call contour_convolution(chainlen_matrix, nx, ns_matrix, coeff, qfree_final, qfree_final, phi_free)
call contour_convolution(chainlen_matrix, nx, ns_matrix, coeff, qads_final, qads_final, phi_loop)
call contour_convolution(chainlen_matrix, nx, ns_matrix, coeff, qfree_final, qads_final, phi_tail)

do kk = 0, nx
    phi_ads(kk) = phi_matrix(kk) - phi_free(kk)
enddo

open (unit=37, file="o.ads_free_loop_tail")
write(37,'(A14,5(2X,A14))') 'r', "phi_ads", "phi_free", "phi_loop", "phi_tail", "phi_matrix"
do kk = 0, nx
    write(37,'(F14.7,5(2X,F14.7))') rx(kk), phi_ads(kk), phi_free(kk), phi_loop(kk), 2.d0 * phi_tail(kk), phi_matrix(kk)
enddo
close (37)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_matrix_ads_free
