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
real(8) :: r_pol, r_solid, h12, Urep, Uatt, iUrep, iUatt, xx
character(20) :: field_filename = "field.in.bin"
!----------------------------------------------------------------------------------------------------------!
Urep = 0.d0
Uatt = 0.d0
h12 = 0.d0

open(unit=211, file = 'o.Usolid.out.txt')
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

if (read_field) then

    INQUIRE(FILE=field_filename, EXIST=FILE_EXISTS)

    if (FILE_EXISTS) then
        open(unit=21, file = field_filename, form = 'unformatted')
        read (21) wa_ifc
        close(21)
    else
        write(*,'(''File '',A15,'' does not exist!'')')field_filename
        STOP
    endif
endif

! Store the initial field for backup
wa_ifc_backup = wa_ifc

return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_field
