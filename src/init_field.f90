subroutine init_field()
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
integer :: kk, side
logical :: FILE_EXISTS
real(8) :: r_pol, r_solid, h12, Urep, Uatt, xx
character(20) :: field_filename = "field.in.bin"
!----------------------------------------------------------------------------------------------------------!
Urep = 0.d0
Uatt = 0.d0
h12 = 0.d0

open(unit=211, file = 'o.Usolid.out.txt')
write(211,('(A16,4(2X,A19))')) "h12", "Ufield(kk)", "Uatt", "Urep", "Uatt+Urep"

do kk = 0, nx
    do side = -1, 1, 2
        if (side.eq.F_lo) xx = rx(kk)
        if (side.eq.F_hi) xx = rx(nx-kk)
        if (wall_side.ne.side.and.wall_side.ne.F_both) then
            cycle
        endif
        Urep = 0.d0
        Uatt = 0.d0

        if (wall_type.eq.F_hamaker) then
            r_pol = (3./4./pi/(rho_mol_bulk*N_avog))**(1./3.)*1.e10
            r_solid = sphere_radius
            h12 = xx - r_pol

            if (geometry.eq.F_sphere) then
                call hamaker_sphere_sphere(h12*1.e-10, r_pol*1.e-10, r_solid*1.e-10, sig_pol*1.e-10, &
&                    sig_solid*1.e-10, Apol, Asolid, Urep, Uatt)
            elseif (geometry.eq.F_film) then
                call hamaker_sphere_plate( h12*1.e-10, r_pol*1.e-10,                 sig_pol*1.e-10, &
&                    sig_solid*1.e-10, Apol, Asolid, Urep, Uatt)
            endif

            if (h12 .le. 0.d0) then
                write(*,*)"ERROR: Hamaker polymer/solid distance < 0!"
                write(*,*)"       Increase the distance of the hard sphere wall."
                STOP
            endif
        elseif (wall_type.eq.F_square_well) then
            if (xx < sigma_sq_well) Uatt = -A_sq_well
        elseif (wall_type.eq.F_vacuum) then
            continue
        endif
        Urep = Urep / (boltz_const_Joule_K * Temp)
        Uatt = Uatt / (boltz_const_Joule_K * Temp)
        Ufield(kk) = Ufield(kk) + Urep + Uatt
    enddo

    write(211,('(E16.9,4(2X,E19.9))')) h12, Ufield(kk), Uatt, Urep, Uatt+Urep
enddo
close(211)

if (read_field) then

    INQUIRE(FILE=field_filename, EXIST=FILE_EXISTS)

    if (FILE_EXISTS) then
        open(unit=21, file = field_filename, form = 'unformatted')
        read (21) wa
        close(21)
    else
        write(*,'(''File '',A15,'' does not exist!'')')field_filename
        STOP
    endif
endif

return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_field
