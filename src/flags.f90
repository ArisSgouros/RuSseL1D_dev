module flags
!----------------------------------------------------------------------------------------------------------!
! mesh discretization
integer, parameter :: F_uniform           = 0
integer, parameter :: F_nonuniform        = 1
! boundary conditions
integer, parameter :: F_bc_neuman         = -1
integer, parameter :: F_bc_dirichlet_eq_0 = 0
integer, parameter :: F_bc_dirichlet_eq_1 = 1
! system geometry
integer, parameter :: F_film              = 0
integer, parameter :: F_cylinder          = -1
integer, parameter :: F_sphere            = 1
! EoS
integer, parameter :: F_helfand           = 0
integer, parameter :: F_sanchez_labombe   = 1
! wall types
integer, parameter :: F_vacuum            = 0
integer, parameter :: F_hamaker           = 1
integer, parameter :: F_square_well       = 2
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
!----------------------------------------------------------------------------------------------------------!
end module flags
