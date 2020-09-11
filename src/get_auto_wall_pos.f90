function get_auto_wall_pos()
!----------------------------------------------------------------------------------------------------------!
use flags
use parser_vars
use constants
use force_fields
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
real(8) :: get_auto_wall_pos
real(8) :: a1, a2, h12, Urep, Uatt, E_wall, xlo, xhi, xx, error
!----------------------------------------------------------------------------------------------------------!
a1 = (3./4./pi/rho_seg_bulk)**(1./3.)
a2 = sphere_radius*1.e-10

xlo = 0.d0
xhi = min(lx*1.d-10, 3.d0 * a1)


error = 1.d10
do while (.true.)
    xx = 0.5d0 * (xlo + xhi)
    h12 = xx - a1

    if (h12 .le. 0.d0) then
       xlo = xx
       cycle
    endif

    if (geometry.eq.F_sphere) then
        call hamaker_sphere_sphere(h12, a1, a2, sig_pol*1.e-10, &
&            sig_solid*1.e-10, Apol, Asolid, Urep, Uatt)
    elseif (geometry.eq.F_film) then
        call hamaker_sphere_plate(h12, a1,     sig_pol*1.e-10, &
&            sig_solid*1.e-10, Apol, Asolid, Urep, Uatt)
    endif

    Urep = Urep / (boltz_const_Joule_K * Temp)
    Uatt = Uatt / (boltz_const_Joule_K * Temp)

    error = abs(E_wall - (Urep + Uatt))
    E_wall = Urep + Uatt

    if (error<1.d-8) then
        !write(*,'(A20, 5(E18.9))')"Convergenced", E_wall, Urep, Uatt, xlo*1.d10, xhi*1.d10
        exit
    endif

    if (E_wall .gt. E_wall_target) then
       xlo = xx
    else
       xhi = xx
    endif
enddo

get_auto_wall_pos = xx * 1.d10

return
!----------------------------------------------------------------------------------------------------------!
end function get_auto_wall_pos
