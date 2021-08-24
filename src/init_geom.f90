!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_geom()
!----------------------------------------------------------------------------------------------------------!
use flags,        only: F_film, F_sphere
use parser_vars,  only: geometry, wall_pos, wall_auto, nx, sphere_radius
use constants,    only: pi, iow
use arrays,       only: rx, rr, irr, layer_area, coeff_nx, volume, surface_area
use write_helper, only: adjl
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: kk

real(8) :: get_auto_wall_pos
!----------------------------------------------------------------------------------------------------------!
if (wall_auto) then
    write(iow,'(3X,A45)')adjl("automatic recalibration of wall position..",45)
    write(*  ,'(3X,A45)')adjl("automatic recalibration of wall position..",45)
    wall_pos = get_auto_wall_pos()
    write(iow,'(3X,A45,F16.9,'' Angstrom'')')adjl("wall_pos was recalibrated to",45), wall_pos
    write(*  ,'(3X,A45,F16.9,'' Angstrom'')')adjl("wall_pos was recalibrated to",45), wall_pos
endif

if (geometry.eq.F_film) then
    do kk = 0, nx
        rx(kk)  = rx(kk) + wall_pos
        rr(kk)  = rx(kk)
        irr(kk) = 1.d0 / rr(kk)
        layer_area(kk) = 1.d0
    enddo
    surface_area = layer_area(0)
elseif (geometry.eq.F_sphere) then
    do kk = 0, nx
        rx(kk)  = rx(kk) + wall_pos
        rr(kk)  = rx(kk) + sphere_radius
        irr(kk) = 1.d0 / rr(kk)
        layer_area(kk) = 4.d0*pi*rr(kk)**2.
    enddo
    surface_area = 4.d0 * pi * sphere_radius**2
endif
volume = 0.d0
do kk = 0, nx
    volume = volume + coeff_nx(kk) * layer_area(kk)
enddo
return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_geom
