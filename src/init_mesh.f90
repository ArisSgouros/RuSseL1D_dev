!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_mesh()
!----------------------------------------------------------------------------------------------------------!
use flags
use eos
use parser_vars, only: spatial_integr_scheme, contour_integr_scheme, spatial_discret_scheme, lx, nx,       &
                     & contour_discret_scheme, ns_matrixA, ns_matrix_auxA, ns_grafted_hi, ns_grafted_lo,     &
                     & chainlen_matrixA, chainlen_grafted_hi, chainlen_grafted_lo, matrix_existA,            &
                     & grafted_lo_exist, grafted_hi_exist, chainlen_matrix_auxA 
use constants
use arrays, only: coeff_nx, coeff_ns_matrixA, coeff_ns_matrix_auxA, coeff_ns_grafted_lo, coeff_ns_grafted_hi,&
                & rx, rs_matrixA, rs_matrix_auxA, rs_grafted_hi, rs_grafted_lo, dx, ds_matrixA, ds_matrix_auxA,&
                & ds_grafted_lo, ds_grafted_hi
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
logical :: symmetric
!----------------------------------------------------------------------------------------------------------!
!discretize the spatial dimension and retrieve dx, rx and coeff_nx
symmetric = .false.
call generate_mesh(spatial_integr_scheme, spatial_discret_scheme, symmetric, lx, nx, dx, rx, coeff_nx)

!discretize the chain contour retrieve ds, rs and coeff_ns
symmetric = .true.

call generate_mesh(contour_integr_scheme, contour_discret_scheme, symmetric, chainlen_matrix_auxA, ns_matrix_auxA, &
&                                                             ds_matrix_auxA, rs_matrix_auxA, coeff_ns_matrix_auxA)

if (matrix_existA)     call generate_mesh(contour_integr_scheme, contour_discret_scheme, symmetric, &
&                                        chainlen_matrixA, ns_matrixA, ds_matrixA, rs_matrixA, coeff_ns_matrixA)
if (grafted_lo_exist) call generate_mesh(contour_integr_scheme, contour_discret_scheme, symmetric, &
&                                        chainlen_grafted_lo, ns_grafted_lo, ds_grafted_lo,        &
&                                        rs_grafted_lo, coeff_ns_grafted_lo)
if (grafted_hi_exist) call generate_mesh(contour_integr_scheme, contour_discret_scheme, symmetric, &
&                                        chainlen_grafted_hi, ns_grafted_hi, ds_grafted_hi,        &
&                                        rs_grafted_hi, coeff_ns_grafted_hi)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_mesh
