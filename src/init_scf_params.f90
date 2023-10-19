!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_scf_params()
!----------------------------------------------------------------------------------------------------------!
use eos,          only: rsl_N, T_tilde, P_tilde, rho_tilde_bulk, T_star, P_star, V_star, rho_star,             &
                      & eos_type, eos_rho_tilde_0
use flags,        only: F_sanchez_lacombe
use parser_vars,  only: bond_length, chainlen_matrixA, chainlen_grafted_lo, chainlen_grafted_hi, ds_ave_matrixA, &
                      & ds_ave_grafted_lo, ds_ave_grafted_hi, ns_matrixA, ns_grafted_lo, ns_grafted_hi,         &
                      & matrix_existA, grafted_lo_exist, grafted_hi_exist, Rg2_per_mon, CN, rho_seg_bulk,       &
                      & rho_mol_bulk, rho_mass_bulk, pressure, Temp, k_gr, k_gr_tilde, mon_mass,               &
                      & square_gradient, gdens_hi, gdens_lo, ns_matrix_auxA, chainlen_matrix_auxA, chainlen_bulk
use constants,    only: N_avog, boltz_const_Joule_molK, boltz_const_Joule_K, gr_cm3_to_kg_m3, iow, tol
use write_helper, only: adjl
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
real(8) :: ds_matrix_auxA, SL_kappa_T
!----------------------------------------------------------------------------------------------------------!
write(iow,'(A85)')adjl("-----------------------------INITIALIZE THE SCF PARAMETERS---------------------------",85)
write(*  ,'(A85)')adjl("-----------------------------INITIALIZE THE SCF PARAMETERS---------------------------",85)

Rg2_per_mon  = bond_length**2 * CN / 6.d00

if (matrix_existA) then
    ns_matrixA           = 2 * nint(0.5d0 * chainlen_matrixA / ds_ave_matrixA)
    ds_matrix_auxA       = ds_ave_matrixA
    write(iow,'(3X,A45,F16.4," Angstrom")')adjl("Matrix radious of gyration:",45), sqrt(Rg2_per_mon*chainlen_matrixA)
    write(*  ,'(3X,A45,F16.4," Angstrom")')adjl("Matrix radious of gyration:",45), sqrt(Rg2_per_mon*chainlen_matrixA)
    write(iow,'(3X,A45,I16," nodes")')     adjl("Matrix nodes along chain contour:",45), ns_matrixA
    write(*  ,'(3X,A45,I16," nodes")')     adjl("Matrix nodes along chain contour:",45), ns_matrixA
endif
if (grafted_lo_exist) then
    ns_grafted_lo       = 2 * nint(0.5d0 * chainlen_grafted_lo / ds_ave_grafted_lo)
    ds_matrix_auxA       = ds_ave_grafted_lo
    write(iow,'(3X,A45,F16.4," Angstrom")')adjl("grafted lo radious of gyration:",45), sqrt(Rg2_per_mon*chainlen_grafted_lo)
    write(*  ,'(3X,A45,F16.4," Angstrom")')adjl("grafted lo radious of gyration:",45), sqrt(Rg2_per_mon*chainlen_grafted_lo)
    write(iow,'(3X,A45,I16," nodes")')     adjl("grafted lo nodes along chain contour:",45), ns_grafted_lo
    write(*  ,'(3X,A45,I16," nodes")')     adjl("grafted lo nodes along chain contour:",45), ns_grafted_lo
endif
if (grafted_hi_exist) then
    ns_grafted_hi       = 2 * nint(0.5d0 * chainlen_grafted_hi / ds_ave_grafted_hi)
    ds_matrix_auxA       = ds_ave_grafted_hi
    write(iow,'(3X,A45,F16.4," Angstrom")')adjl("grafted hi radious of gyration:",45), sqrt(Rg2_per_mon*chainlen_grafted_hi)
    write(*  ,'(3X,A45,F16.4," Angstrom")')adjl("grafted hi radious of gyration:",45), sqrt(Rg2_per_mon*chainlen_grafted_hi)
    write(iow,'(3X,A45,I16," nodes")')     adjl("grafted hi nodes along chain contour:",45), ns_grafted_hi
    write(*  ,'(3X,A45,I16," nodes")')     adjl("grafted hi nodes along chain contour:",45), ns_grafted_hi
endif

ns_matrix_auxA       = max(ns_matrixA, ns_grafted_lo, ns_grafted_hi)
chainlen_matrix_auxA = max(chainlen_matrixA, chainlen_grafted_lo,chainlen_grafted_hi)
ds_matrix_auxA       = max(ds_ave_matrixA, ds_ave_grafted_lo, ds_ave_grafted_hi)

! miscellenious checks 
if ((matrix_existA    .and.abs(ds_matrix_auxA-ds_ave_matrixA)    >tol) .or. &
&   (grafted_lo_exist.and.abs(ds_matrix_auxA-ds_ave_grafted_lo)>tol) .or. &
&   (grafted_hi_exist.and.abs(ds_matrix_auxA-ds_ave_grafted_hi)>tol)) then
    write(iow,*) "Error: the nonconstant contour discret scheme does not work with different chain lengths."
    write(*  ,*) "Error: the nonconstant contour discret scheme does not work with different chain lengths."
    STOP
endif


if (eos_type.eq.F_sanchez_lacombe) then
    write(iow,'(3X,A45)')adjl("Computation of the mass density from SL EoS..",45)
    write(*  ,'(3X,A45)')adjl("Computation of the mass density from SL EoS..",45)
    V_star         = boltz_const_Joule_K * T_star / P_star
    T_tilde        = Temp  / T_star
    P_tilde        = Pressure / P_star
    rsl_N          = (mon_mass * P_star) / (rho_star * 1.d03 * boltz_const_Joule_molK * T_star)         
    rho_tilde_bulk = eos_rho_tilde_0(T_tilde, P_tilde, rsl_N*chainlen_bulk)
    rho_mass_bulk  = rho_tilde_bulk * rho_star
    rho_mass_bulk  = rho_mass_bulk
    write(iow,'(3X,A45,F16.4," g/cm3")')adjl("mass density was recomputed as:",45), rho_mass_bulk/gr_cm3_to_kg_m3
    write(*  ,'(3X,A45,F16.4," g/cm3")')adjl("mass density was recomputed as:",45), rho_mass_bulk/gr_cm3_to_kg_m3

    if (square_gradient.and.k_gr_tilde>tol) k_gr = 2.d0 * P_star * rsl_N**2 * V_star**(8.d0/3.d0) * k_gr_tilde

    write(iow,'(3X,A45,E16.8," J m^5")')adjl("Sanchez-Lacombe influence parameter:",45), k_gr
    write(*  ,'(3X,A45,E16.8," J m^5")')adjl("Sanchez-Lacombe influence parameter:",45), k_gr

    SL_kappa_T = 1.0d0 / (  P_star * T_tilde * rho_tilde_bulk * ( 1.d0 / (1.d0 / rho_tilde_bulk - 1.d0) + &
&                                          1.d0 / (rsl_N*chainlen_bulk)) - 2.d0*rho_tilde_bulk**2*P_star)
    write(iow,'(3X,A45,E16.4," Pa-1")')adjl("Sanchez-Lacombe isothermal compressibility:",45), SL_kappa_T
    write(*  ,'(3X,A45,E16.4," Pa-1")')adjl("Sanchez-Lacombe isothermal compressibility:",45), SL_kappa_T
endif

rho_mol_bulk = rho_mass_bulk/mon_mass*gr_cm3_to_kg_m3
rho_seg_bulk = rho_mol_bulk*N_avog

write(iow,'(3X,A45,F16.4," mol/m3")')adjl("molar density in bulk",45), rho_mol_bulk
write(*  ,'(3X,A45,F16.4," mol/m3")')adjl("molar density in bulk",45), rho_mol_bulk
!----------------------------------------------------------------------------------------------------------!
end subroutine init_scf_params
