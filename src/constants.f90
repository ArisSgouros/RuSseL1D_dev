!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module constants
!----------------------------------------------------------------------------------------------------------!
real(8), parameter :: N_avog                 = 6.02214085700E+23 !Avogadro's Constant
real(8), parameter :: pi                     = 3.14159265359
real(8), parameter :: boltz_const_Joule_molK = 8.3144598         !Boltzmann constant in J/mol/K.
real(8), parameter :: boltz_const_Joule_K    = 1.38064852E-23    !Boltzmann constant in J/K.
real(8), parameter :: atm_to_pa              = 101325            !pa
real(8), parameter :: gr_cm3_to_kg_m3        = 1000              !kg/m^3
real(8), parameter :: tol                    = 1.0E-13

integer, parameter :: iow = 10
!----------------------------------------------------------------------------------------------------------!
end module constants
