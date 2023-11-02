!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module flags
!----------------------------------------------------------------------------------------------------------!
! mesh discretization
integer, parameter :: F_uniform           = 0
integer, parameter :: F_nonuniform        = 1
! boundary conditions
integer, parameter :: F_bc_neuman         = -1
integer, parameter :: F_bc_dirichlet_eq_0 = 0
integer, parameter :: F_bc_dirichlet_eq_1 = 1
integer, parameter :: F_bc_periodic       = 2
! system geometry
integer, parameter :: F_film              = 0
integer, parameter :: F_cylinder          = -1
integer, parameter :: F_sphere            = 1
! EoS
integer, parameter :: F_helfand           = 0
integer, parameter :: F_sanchez_lacombe   = 1
! wall types
integer, parameter :: F_hybrid            = -1
integer, parameter :: F_vacuum            = 0
integer, parameter :: F_hamaker           = 1
integer, parameter :: F_square_well       = 2
integer, parameter :: F_ramp              = 3
integer, parameter :: F_hamaker_well      = 4
integer, parameter :: F_custom            = 9
integer, parameter :: F_table             = 10
! chain types
integer, parameter :: F_lo                = -1
integer, parameter :: F_both              = 0
integer, parameter :: F_hi                = 1
! edwards solution scheme
integer, parameter :: F_implicit          = 0
integer, parameter :: F_semi_implicit     = 1
! edwards solution scheme
integer, parameter :: F_rectangle_rule    = 0
integer, parameter :: F_simpson_rule      = 1
! linear system solver
integer, parameter :: F_tridag            = 0
integer, parameter :: F_gelim             = 1
!----------------------------------------------------------------------------------------------------------!
end module flags
