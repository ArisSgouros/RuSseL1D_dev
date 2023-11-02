!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_mesh()
!----------------------------------------------------------------------------------------------------------!
use flags
use eos
use parser_vars, only: spatial_integr_scheme, contour_integr_scheme, spatial_discret_scheme, lx, nx,       &
                     & contour_discret_scheme, ns_mxa,ns_mxb, ns_ghi, ns_glo,     &
                     & chainlen_mxa,chainlen_mxb, chainlen_ghi, chainlen_glo, exist_mxa,exist_mxb,            &
                     & exist_glo, exist_ghi
use constants
use arrays, only: coeff_nx, coeff_ns_mxa,coeff_ns_mxb, coeff_ns_glo, coeff_ns_ghi,&
                & rx, rs_mxa,rs_mxb, rs_ghi, rs_glo, dx, ds_mxa,ds_mxb, &
                & ds_glo, ds_ghi
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

if (exist_mxa) call generate_mesh(contour_integr_scheme, contour_discret_scheme, symmetric, &
&                                        chainlen_mxa, ns_mxa, ds_mxa, rs_mxa, coeff_ns_mxa)
if (exist_mxb) call generate_mesh(contour_integr_scheme, contour_discret_scheme, symmetric, &
&                                        chainlen_mxb, ns_mxb, ds_mxb, rs_mxb, coeff_ns_mxb)
if (exist_glo) call generate_mesh(contour_integr_scheme, contour_discret_scheme, symmetric, &
&                                        chainlen_glo, ns_glo, ds_glo,        &
&                                        rs_glo, coeff_ns_glo)
if (exist_ghi) call generate_mesh(contour_integr_scheme, contour_discret_scheme, symmetric, &
&                                        chainlen_ghi, ns_ghi, ds_ghi,        &
&                                        rs_ghi, coeff_ns_ghi)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_mesh
