!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_solid()
!----------------------------------------------------------------------------------------------------------!
use flags
use eos
use parser_vars
use constants
use arrays
use force_fields
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: kk, side, iot, io

logical :: file_exists

real(8), dimension(0:nx) :: u_table
real(8)                  :: r_pol, r_solid, h12, Urep, Uatt, iUrep, iUatt, xx, aux

character(20)  :: table_filename = "in.table"
character(200) :: line
!----------------------------------------------------------------------------------------------------------!
Urep = 0.d0
Uatt = 0.d0
h12  = 0.d0

iot = 35

if (wall_table) then
   inquire(file=table_filename, exist=file_exists)
   if (.not.file_exists) then
       write(iow,'("File ",A15," does not exist!")') table_filename
       write(*  ,'("File ",A15," does not exist!")') table_filename
       STOP
   endif
   u_table = 0.d0
   open(unit=iot, file = table_filename, form="FORMATTED", status="OLD", action="READ")
   do kk = 0, nx
      read(iot,"(A200)",iostat=io) line
      if (io/=0) exit
      READ(line, *) aux, u_table(kk)
      !write(*,*) kk, aux, u_table(kk)
   END DO
   CLOSE(iot)
endif

open(unit=211, file = "o.Usolid")
write(211,('(A16,5(2X,A19))')) "r12", "h12", "Ufield(kk)", "Uatt", "Urep", "Uatt+Urep"

do kk = 0, nx
    do side = -1, 1, 2
        if (side.eq.F_lo) xx = rx(kk)
        if (side.eq.F_hi) xx = rx(nx-kk)
        if (wall_side.ne.side.and.wall_side.ne.F_both) then
            cycle
        endif
        Urep = 0.d0
        Uatt = 0.d0

        if (wall_hamaker) then
            r_pol = (3./4./pi/(rho_mol_bulk*N_avog))**(1./3.)*1.e10
            r_solid = sphere_radius
            h12 = xx - r_pol

            if (geometry.eq.F_sphere) then
                call hamaker_sphere_sphere(h12*1.e-10, r_pol*1.e-10, r_solid*1.e-10, sig_pol*1.e-10, &
&                    sig_solid*1.e-10, Apol, Asolid, iUrep, iUatt)
            elseif (geometry.eq.F_film) then
                call hamaker_sphere_plate( h12*1.e-10, r_pol*1.e-10,                 sig_pol*1.e-10, &
&                    sig_solid*1.e-10, Apol, Asolid, iUrep, iUatt)
            endif
            Urep = Urep + iUrep
            Uatt = Uatt + iUatt

            if (h12 .le. 0.d0) then
                write(*,*)"ERROR: Hamaker polymer/solid distance < 0!"
                write(*,*)"       Increase the distance of the hard sphere wall."
                STOP
            endif
        endif
        if (wall_square_well) then
            if (xx < sigma_sq_well) Uatt = Uatt + A_sq_well
        endif
        if (wall_ramp) then
            if (xx < sigma_ramp) Uatt = Uatt + A_ramp * (sigma_ramp - xx) / sigma_ramp
        endif
        if (wall_custom) then
            Uatt = Uatt + exp(-xx/wall_custom_vars(2)) * ( wall_custom_vars(4) + wall_custom_vars(1) * &
&                         sin((2 * PI / wall_custom_vars(3)) * (xx+ wall_custom_vars(5)) ))
        endif
        if (wall_table) then
            Uatt = Uatt + u_table(kk)
        endif
        if (wall_vacuum) then
            continue
        endif
        Urep = Urep * beta
        Uatt = Uatt * beta
        Ufield(kk) = Ufield(kk) + Urep + Uatt
    enddo

    write(211,('(E16.9,5(2X,E19.9))')) rx(kk), h12, Ufield(kk), Uatt, Urep, Uatt+Urep
enddo
close(211)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_solid
