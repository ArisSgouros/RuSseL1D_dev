!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

function get_auto_wall_pos()
!----------------------------------------------------------------------------------------------------------!
use flags,        only: F_sphere, F_film
use parser_vars,  only: rho_seg_bulk, wall_hamaker, wall_ramp, wall_square_well, sphere_radius, sig_solid, &
                      & sig_pol, lx, sigma_ramp, sigma_sq_well, geometry, beta, E_wall_target, Apol,       &
                      & Asolid, A_sq_well, A_ramp
use constants,    only: pi
use force_fields, only: hamaker_sphere_sphere, hamaker_sphere_plate
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
real(8) :: get_auto_wall_pos
real(8) :: a1, a2, h12, Urep, Uatt, E_wall, xlo, xhi, xx, error, iUrep, iUatt
real(8) :: r_pol, r_solid
!----------------------------------------------------------------------------------------------------------!
a1 = (3./4./pi/rho_seg_bulk)**(1./3.)
a2 = sphere_radius*1.e-10

xlo = 0.d0
xhi = min(lx*1.d-10, 3.d0 * a1)

error = 1.d10
do while (.true.)
    xx = 0.5d0 * (xlo + xhi)

    Urep = 0.d0
    Uatt = 0.d0

    if (wall_hamaker) then
        r_pol = (3./4./pi/rho_seg_bulk)**(1./3.)
        r_solid = sphere_radius*1.e-10
        h12 = xx - r_pol

        if (h12 .le. 0.d0) then
            xlo = xx
            cycle
        endif
        if (geometry.eq.F_sphere) then
            call hamaker_sphere_sphere(h12, r_pol, r_solid, sig_pol*1.e-10, &
&                sig_solid*1.e-10, Apol, Asolid, iUrep, iUatt)
        elseif (geometry.eq.F_film) then
            call hamaker_sphere_plate( h12, r_pol,          sig_pol*1.e-10, &
&                sig_solid*1.e-10, Apol, Asolid, iUrep, iUatt)
        endif
        Urep = Urep + iUrep
        Uatt = Uatt + iUatt
    endif

    if (wall_square_well) then
        if (xx < sigma_sq_well*1.d-10) Uatt = Uatt + A_sq_well
    endif
    if (wall_ramp) then
        if (xx < sigma_ramp*1.d-10)    Uatt = Uatt + A_ramp * (sigma_ramp*1.d-10 - xx) / (sigma_ramp*1.d-10)
    endif

    Urep = Urep * beta
    Uatt = Uatt * beta

    error = abs(E_wall - (Urep + Uatt))
    E_wall = Urep + Uatt

    if (error<1.d-8) then
        !write(*,'(A20, 7(E18.9))')"Convergenced", E_wall_target, E_wall, error, Urep, Uatt, xlo*1.d10, xhi*1.d10
        exit
    endif

    if (E_wall .gt. E_wall_target) then
       xlo = xx
    else
       xhi = xx
    endif
enddo

if (E_wall_target - E_wall < 1.d-4) then
    write(*,*) "Automatic recalibration of wall position converged!"
else
    write(*,*) "Automatic recalibration of wall position could not converge!"
    STOP
endif

get_auto_wall_pos = xx * 1.d10

return
!----------------------------------------------------------------------------------------------------------!
end function get_auto_wall_pos
