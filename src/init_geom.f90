subroutine init_geom()
!----------------------------------------------------------------------------------------------------------!
use flags
use eos
use parser_vars
use constants
use arrays
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: kk
!----------------------------------------------------------------------------------------------------------!
!initialize the distance of each node with the center of NP for each geometry
!geometry = 0 (film)
!geometry = 1 (sphere)
if (geometry.eq.F_film) then
    do kk = 0, nx
        rx(kk) = rx(kk) + wall_pos
        rr(kk) = rx(kk)
        irr(kk) = 1.d0 / rr(kk)
        layer_area(kk) = 1.d0
    enddo
    volume = lx * 1.d0**2
    surface_area = 1.d0
elseif (geometry.eq.F_sphere) then
    do kk = 0, nx
        rx(kk) = rx(kk) + wall_pos
        rr(kk) = rx(kk) + sphere_radius
        irr(kk) = 1.d0 / rr(kk)
        layer_area(kk) = 4.d0*pi*rr(kk)**2.
    enddo
    volume = 4./3.*pi*(rr(nx)**3.-rr(0)**3.)
    surface_area = 4.d0 * pi * sphere_radius**2
endif

return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_geom
