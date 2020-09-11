program fd_1d
!----------------------------------------------------------------------------------------------------------!
use parser_vars
use constants
use flags
use arrays
use eos
use write_helper
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
logical :: convergence = .false.

integer :: iter, jj, ii, tt

real(8) :: get_nchains
real(8) :: wa_error = 1.d10, en_error = 1.d10, free_energy_prev = 1.d10
real(8) :: qinit_lo = 0.d0, qinit_hi = 0.d0
real(8) :: free_energy = 0.d0
real(8) :: nchgr_lo = 0.d0, nchgr_hi = 0.d0, nch_matrix = 0.d0
!----------------------------------------------------------------------------------------------------------!
open(unit=iow, file = 'o.scft.out.txt')

call parser
call init_scf_params
call init_arrays
call init_mesh
call init_geom
call init_field

write(iow,'(A85)')adjl('---------------------------------BEGIN THE SIMULATION--------------------------------',85)
write(*  ,'(A85)')adjl('---------------------------------BEGIN THE SIMULATION--------------------------------',85)

write(*,  '(2X,A15,4(2X,A15))') 'Iteration', 'energy (mN/m)','error (k_B T)', 'gdens_lo', 'gdens_hi'
write(iow,'(2X,A15,4(2X,A15))') 'Iteration', 'energy (mN/m)','error', 'gdens_lo', 'gdens_hi'

do iter = 0, max_iter

    !matrix chains
    do ii = 0, nx
        qmatrix(ii,1)       = 1.d0
        qmatrix_final(ii,0) = 1.d0
    enddo
  
    !set the dirichlet boundary conditions for matrix chains
    n_dir_nodes = 0
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
        n_dir_nodes = n_dir_nodes + 1   
    endif

    if (geometry.eq.F_sphere) then
        do ii = 0, nx
            qmatrix(ii,1)       = qmatrix(ii,1) * rr(ii)
            qmatrix_final(ii,0) = qmatrix_final(ii,0) * rr(ii)
        enddo
        do ii = 0, n_dir_nodes-1
            dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
        enddo
    endif
 
    call edwards(bc_lo_matrix, bc_hi_matrix, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&           Rg2_per_mon, nx, ns_matrix_aux, dx, ds_matrix_aux, edwards_solver, wa, qmatrix, qmatrix_final)

    if (geometry.eq.F_sphere) then
        do tt = 0, ns_matrix_aux
            do ii = 0, nx
                qmatrix_final(ii,tt) = qmatrix_final(ii,tt)*irr(ii)
            enddo
        enddo
    endif

    call convolution(chainlen_matrix, nx, ns_matrix, coeff_ns_matrix, qmatrix_final, qmatrix_final, phi_matrix)

    !edwards diffusion for grafted chains in the lower boundary
    if (grafted_lo_exist) then
        do ii = 0, nx
            qgr_lo(ii,1)       = 0.d0
            qgr_final_lo(ii,0) = 0.d0
        enddo

        delta = (1.0/( dx(gnode_lo) * layer_area(gnode_lo))) * 1e+30 !1/m3
        qinit_lo = chainlen_grafted_lo * delta * (gdens_lo * surface_area) &
&                     / (rho_seg_bulk * qmatrix_final(gnode_lo,ns_grafted_lo))

        qgr_lo(gnode_lo,1)       = qinit_lo
        qgr_final_lo(gnode_lo,0) = qinit_lo

        !set the dirichlet boundary conditions for grafted chains
        n_dir_nodes = 0
        !dirichlet lower bound
        if (bc_lo_grafted.eq.F_bc_dirichlet_eq_0) then
            dir_nodes_id(n_dir_nodes) = 0
            dir_nodes_rdiag(n_dir_nodes) = 0.d0
            n_dir_nodes = n_dir_nodes + 1   
        else if (bc_lo_grafted.eq.F_bc_dirichlet_eq_1) then
            dir_nodes_id(n_dir_nodes) = 0
            dir_nodes_rdiag(n_dir_nodes) = 1.0d0
            n_dir_nodes = n_dir_nodes + 1   
        endif
        !dirichlet upper bound
        if (bc_hi_grafted.eq.F_bc_dirichlet_eq_0) then
            dir_nodes_id(n_dir_nodes) = nx
            dir_nodes_rdiag(n_dir_nodes) = 0.d0
            n_dir_nodes = n_dir_nodes + 1   
        else if (bc_hi_grafted.eq.F_bc_dirichlet_eq_1) then
            dir_nodes_id(n_dir_nodes) = nx
            dir_nodes_rdiag(n_dir_nodes) = 1.0d0
            n_dir_nodes = n_dir_nodes + 1   
        endif

        if (geometry.eq.F_sphere) then
            qgr_lo(gnode_lo,1)       = qgr_lo(gnode_lo,1) * rr(gnode_lo)
            qgr_final_lo(gnode_lo,0) = qgr_final_lo(gnode_lo,0) * rr(gnode_lo)
            do ii = 0, n_dir_nodes-1
                dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
            enddo
        endif

        call edwards(bc_lo_grafted, bc_hi_grafted, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&                Rg2_per_mon, nx, ns_grafted_lo, dx, ds_grafted_lo, edwards_solver, wa, qgr_lo, qgr_final_lo)

        if (geometry.eq.F_sphere) then
            do tt = 0, ns_grafted_lo
                do ii = 0, nx
                    qgr_final_lo(ii,tt) = qgr_final_lo(ii,tt) * irr(ii)
                enddo
            enddo
        endif

        call convolution(chainlen_grafted_lo, nx, ns_grafted_lo, coeff_ns_grafted_lo, qmatrix_final, qgr_final_lo, phi_gr_lo)
    endif


    !edwards diffusion for grafted chains in the higher boundary
    if (grafted_hi_exist) then
        do ii = 0, nx
            qgr_hi(ii,1)       = 0.d0
            qgr_final_hi(ii,0) = 0.d0
        enddo

        delta = (1.0/( dx(gnode_hi) * layer_area(gnode_hi))) * 1e+30 !1/m3
        qinit_hi = chainlen_grafted_hi * delta * (gdens_hi * surface_area) &
&                     / (rho_seg_bulk * qmatrix_final(gnode_hi,ns_grafted_hi))

        qgr_hi(gnode_hi,1)       = qinit_hi
        qgr_final_hi(gnode_hi,0) = qinit_hi

        !set the dirichlet boundary conditions for grafted chains
        n_dir_nodes = 0
        !dirichlet lower bound
        if (bc_lo_grafted.eq.F_bc_dirichlet_eq_0) then
            dir_nodes_id(n_dir_nodes) = 0
            dir_nodes_rdiag(n_dir_nodes) = 0.d0
            n_dir_nodes = n_dir_nodes + 1   
        else if (bc_lo_grafted.eq.F_bc_dirichlet_eq_1) then
            dir_nodes_id(n_dir_nodes) = 0
            dir_nodes_rdiag(n_dir_nodes) = 1.0d0
            n_dir_nodes = n_dir_nodes + 1   
        endif
        !dirichlet upper bound
        if (bc_hi_grafted.eq.F_bc_dirichlet_eq_0) then
            dir_nodes_id(n_dir_nodes) = nx
            dir_nodes_rdiag(n_dir_nodes) = 0.d0
            n_dir_nodes = n_dir_nodes + 1   
        else if (bc_hi_grafted.eq.F_bc_dirichlet_eq_1) then
            dir_nodes_id(n_dir_nodes) = nx
            dir_nodes_rdiag(n_dir_nodes) = 1.0d0
            n_dir_nodes = n_dir_nodes + 1   
        endif

        if (geometry.eq.F_sphere) then
            qgr_hi(gnode_hi,1)       = qgr_hi(gnode_hi,1) * rr(gnode_hi)
            qgr_final_hi(gnode_hi,0) = qgr_final_hi(gnode_hi,0) * rr(gnode_hi)
            do ii = 0, n_dir_nodes-1
                dir_nodes_rdiag(ii) = dir_nodes_rdiag(ii)*rr(dir_nodes_id(ii))
            enddo
        endif

        call edwards(bc_hi_grafted, bc_hi_grafted, n_dir_nodes, dir_nodes_id, dir_nodes_rdiag, &
&                Rg2_per_mon, nx, ns_grafted_hi, dx, ds_grafted_hi, edwards_solver, wa, qgr_hi, qgr_final_hi)

        if (geometry.eq.F_sphere) then
            do tt = 0, ns_grafted_hi
                do ii = 0, nx
                    qgr_final_hi(ii,tt) = qgr_final_hi(ii,tt) * irr(ii)
                enddo
            enddo
        endif

        call convolution(chainlen_grafted_hi, nx, ns_grafted_hi, coeff_ns_grafted_hi, qmatrix_final, qgr_final_hi, phi_gr_hi)
    endif

    phi_total = 0.d0
    do jj = 0, nx
        if (matrix_exist)     phi_total(jj) = phi_total(jj) + phi_matrix(jj)
        if (grafted_lo_exist) phi_total(jj) = phi_total(jj) + phi_gr_lo(jj)
        if (grafted_hi_exist) phi_total(jj) = phi_total(jj) + phi_gr_hi(jj)
    enddo

    if (square_gradient) then
        ! Warning! this is not compatible with the nonuniforma spatial discretization scheme
        if (geometry.eq.F_sphere) then
            do jj = 1, nx-1
                d2phi_dr2(jj) = (phi_total(jj-1)*rr(jj-1) -2.d0*phi_total(jj)*rr(jj) + phi_total(jj+1)*rr(jj+1)) / (dx(jj)*dx(jj)*1.e-20*rr(jj))
                dphi_dr(jj)   = (phi_total(jj+1)*rr(jj+1) - phi_total(jj)*rr(jj)) / (dx(jj)*1.e-10*rr(jj))
            enddo
            d2phi_dr2(0)  = (phi_total(0)*(rr(0)-dx(1))    -2.d0*phi_total(0)*rr(0)  + phi_total(1)*rr(1))  / (dx(1)*dx(1)*1.e-20*rr(0))
            d2phi_dr2(nx) = (phi_total(nx-1)*rr(nx-1) -2.d0*phi_total(nx)*rr(nx) + phi_total(nx)*(rr(nx)+dx(1))) / (dx(nx)*dx(nx)*1.e-20*rr(nx))
        
            dphi_dr(0)  = (phi_total(1)*rr(1) - phi_total(0)*rr(0)) / (dx(1)*1.e-10*rr(0))
            dphi_dr(nx) = 0.d0
        else
            do jj = 1, nx-1
                d2phi_dr2(jj) = (phi_total(jj-1) -2.d0*phi_total(jj) + phi_total(jj+1)) / (dx(jj)*1.e-10)**2
                dphi_dr(jj)   = (phi_total(jj+1) - phi_total(jj)) / (dx(jj)*1.e-10)
            enddo
            d2phi_dr2(0)  = (phi_total(0)    -2.d0*phi_total(0)  + phi_total(1))  / (dx(1) *1.e-10)**2
            d2phi_dr2(nx) = (phi_total(nx-1) -2.d0*phi_total(nx) + phi_total(nx)) / (dx(nx)*1.e-10)**2
        
            dphi_dr(0)  = (phi_total(1) - phi_total(0)) / (dx(1)*1.e-10)
            dphi_dr(nx) = 0.d0
        endif
    endif

    !calculate new field and maximum absolute wa_error
    do jj = 0, nx
        wa_new(jj) = +(eos_df_drho(phi_total(jj)) - eos_df_drho(1.d0)) / (boltz_const_Joule_K * Temp) &
&                    -k_gr * (rho_seg_bulk * d2phi_dr2(jj)) / (boltz_const_Joule_K * Temp) &
&                  + Ufield(jj)
    enddo

    wa_error = 0.d0
    do jj = 0, nx
        wa_error = max(wa_error, dabs((wa_new(jj)-wa(jj))))
    end do

    !apply field mixing rule and update field
    do jj = 0, nx
        wa(jj) = (1.d0 - frac) * wa(jj) + frac * wa_new(jj)
    end do

    call energies(free_energy)
    en_error = abs(free_energy_prev - free_energy)
    free_energy_prev = free_energy

    convergence = (wa_error.lt.max_wa_error).or.(en_error.lt.max_en_error)

    if (mod(iter,field_every)  .eq.0.or.convergence) call export_field_binary(wa,nx)
    if (mod(iter,compute_every).eq.0.or.convergence) call export_computes(qinit_lo, qinit_hi)
    if (mod(iter,thermo_every) .eq.0.or.convergence) then
        nch_matrix = get_nchains(coeff_nx, nx, layer_area, phi_matrix, rho_seg_bulk, chainlen_matrix)
        nchgr_lo   = get_nchains(coeff_nx, nx, layer_area, phi_gr_lo,  rho_seg_bulk, chainlen_grafted_lo)
        nchgr_hi   = get_nchains(coeff_nx, nx, layer_area, phi_gr_hi,  rho_seg_bulk, chainlen_grafted_hi)
        write(*,  '(2X,I15,4(2X,E15.7))') iter, free_energy, wa_error, nchgr_lo / surface_area, nchgr_hi / surface_area
        write(iow,'(2X,I15,4(2X,E15.7))') iter, free_energy, wa_error, nchgr_lo / surface_area, nchgr_hi / surface_area
    endif
    if (convergence) exit
enddo

if (iter.eq.max_iter)         write(iow,'(/''Convergence of max iterations'',I16,'' iterations'')') iter
if (wa_error.lt.max_wa_error) write(iow,'(/''Convergence of max field error'',E16.7,'' k_B T'')') wa_error
if (en_error.lt.max_en_error) write(iow,'(/''Convergence of max energy error'',E16.7,'' mJ/m^2'')') en_error
if (iter.eq.max_iter)         write(*  ,'(/''Convergence of max iterations'',I16,'' iterations'')') iter
if (wa_error.lt.max_wa_error) write(*  ,'(/''Convergence of max field error'',E16.7,'' k_B T'')') wa_error
if (en_error.lt.max_en_error) write(*  ,'(/''Convergence of max energy error'',E16.7,'' mJ/m^2'')') en_error

write(iow,'(A85)')adjl('-------------------------------------FINAL OUTPUT------------------------------------',85)
write(*  ,'(A85)')adjl('-------------------------------------FINAL OUTPUT------------------------------------',85)
write(iow,'(3X,A17,E16.7,'' mN/m'')') adjl('Free energy:',17),free_energy
write(*  ,'(3X,A17,E16.7,'' mN/m'')') adjl('Free energy:',17),free_energy
if (matrix_exist) then
    write(iow,'(A85)')adjl('-------------------------------------MATRIX CHAINS-----------------------------------',85)
    write(*  ,'(A85)')adjl('-------------------------------------MATRIX CHAINS-----------------------------------',85)
    write(iow,'(3X,A17,E16.7,'' chains'')') adjl('matrix chains:',17),nch_matrix
    write(*  ,'(3X,A17,E16.7,'' chains'')') adjl('matrix chains:',17),nch_matrix
endif
if (grafted_lo_exist) then
    write(iow,'(A85)')adjl('---------------------------------GRAFTED LO CHAINS-----------------------------------',85)
    write(*  ,'(A85)')adjl('---------------------------------GRAFTED LO CHAINS-----------------------------------',85)
    write(iow,'(3X,A17,E16.7,'' chains'')')          adjl('grafted lo chains:',17),nchgr_lo
    write(*  ,'(3X,A17,E16.7,'' chains'')')          adjl('grafted lo chains:',17),nchgr_lo
    write(iow,'(3X,A17,E16.7)')                      adjl('q_f(r_gi,0):',17),qmatrix_final(gnode_lo,ns_grafted_lo)
    write(*  ,'(3X,A17,E16.7)')                      adjl('q_f(r_gi,0):',17),qmatrix_final(gnode_lo,ns_grafted_lo)
    write(iow,'(3X,A17,E16.7)')                      adjl('q_g(r_gi,N_g):',17),qgr_final_lo(gnode_lo,0)
    write(*  ,'(3X,A17,E16.7)')                      adjl('q_g(r_gi,N_g):',17),qgr_final_lo(gnode_lo,0)
    write(iow,'(3X,A17,E16.7,'' chains/angstrom'')') adjl('grafting density:',17),nchgr_lo / surface_area
    write(*  ,'(3X,A17,E16.7,'' chains/angstrom'')') adjl('grafting density:',17),nchgr_lo / surface_area
    write(iow,'(3X,A17,E16.7)'                 )     adjl('error:',17),(nchgr_lo / surface_area - gdens_lo) / gdens_lo * 100.0
    write(*  ,'(3X,A17,E16.7)'                 )     adjl('error:',17),(nchgr_lo / surface_area - gdens_lo) / gdens_lo * 100.0
endif
if (grafted_hi_exist) then
    write(iow,'(A85)')adjl('---------------------------------GRAFTED HI CHAINS-----------------------------------',85)
    write(*  ,'(A85)')adjl('---------------------------------GRAFTED HI CHAINS-----------------------------------',85)
    write(iow,'(3X,A17,E16.7,'' chains'')')          adjl('grafted hi chains:',17),nchgr_hi
    write(*  ,'(3X,A17,E16.7,'' chains'')')          adjl('grafted hi chains:',17),nchgr_hi
    write(iow,'(3X,A17,E16.7)')                      adjl('q_f(r_gi,0):',17),qmatrix_final(gnode_hi,ns_grafted_hi)
    write(*  ,'(3X,A17,E16.7)')                      adjl('q_f(r_gi,0):',17),qmatrix_final(gnode_hi,ns_grafted_hi)
    write(iow,'(3X,A17,E16.7)')                      adjl('q_g(r_gi,N_g):',17),qgr_final_hi(gnode_hi,0)
    write(*  ,'(3X,A17,E16.7)')                      adjl('q_g(r_gi,N_g):',17),qgr_final_hi(gnode_hi,0)
    write(iow,'(3X,A17,E16.7,'' chains/angstrom'')') adjl('grafting density:',17),nchgr_hi / surface_area
    write(*  ,'(3X,A17,E16.7,'' chains/angstrom'')') adjl('grafting density:',17),nchgr_hi / surface_area
    write(iow,'(3X,A17,E16.7)'                 )     adjl('error:',17),(nchgr_hi / surface_area - gdens_hi) / gdens_hi * 100.0
    write(*  ,'(3X,A17,E16.7)'                 )     adjl('error:',17),(nchgr_hi / surface_area - gdens_hi) / gdens_hi * 100.0
endif
!----------------------------------------------------------------------------------------------------------!
end program fd_1d
