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
if (geometry.eq.F_film) then
    do kk = 0, nx
        rx(kk) = rx(kk) + wall_pos
        rr(kk) = rx(kk)
        irr(kk) = 1.d0 / rr(kk)
        layer_area(kk) = 1.d0
    enddo
    surface_area = layer_area(0)
elseif (geometry.eq.F_sphere) then
    do kk = 0, nx
        rx(kk) = rx(kk) + wall_pos
        rr(kk) = rx(kk) + sphere_radius
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
