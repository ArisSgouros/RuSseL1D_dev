!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module force_fields
use constants
!--------------------------------------------------------------------!
implicit none
!--------------------------------------------------------------------!
contains

subroutine hamaker_sphere_plate(h_12, r_pol, sigma1, sigma2, A1, A2, Urep, Uatt)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: h_12, r_pol, sigma1, sigma2, A1, A2
real(8), intent(out) :: Urep, Uatt
real(8)              :: sigma = 0.d0, hamaker_constant = 0.d0
!------------------------------------------------------------------------------------------------------!
!***SI units are used everywhere***
Urep = 0.d0
Uatt = 0.d0

hamaker_constant = sqrt(A1*A2)

sigma = (sigma1 + sigma2)/2.d0 

!Ruckenstein and Priere, 1976a; Feke et al., 1984
Uatt = - (hamaker_constant/6.d0) * (1.d0/(h_12/r_pol) + 1.d0/(2.d0+h_12/r_pol) + log((h_12/r_pol)/(2.d0+h_12/r_pol)))

Urep = (hamaker_constant/7.56d3)*(sigma/r_pol)**6.d0*((8.d0 + h_12/r_pol)/(2.d0 + h_12/r_pol)**7.d0 &
                                                           & + (6.d0 - h_12/r_pol)/(h_12/r_pol)**7.d0)
return
!------------------------------------------------------------------------------------------------------!
end subroutine hamaker_sphere_plate

subroutine hamaker_sphere_sphere(h_12, r1, r2, sigma1, sigma2, A1, A2, Urep, Uatt)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: h_12, r1, r2, sigma1, sigma2, A1, A2
real(8), intent(out) :: Urep, Uatt
real(8)              :: r12 = 0.d0, sigma = 0.d0, hamaker_constant = 0.d0
!------------------------------------------------------------------------------------------------------!
!***SI units are used everywhere***
Urep = 0.d0
Uatt = 0.d0

hamaker_constant = sqrt(A1*A2)

sigma = (sigma1 + sigma2)/2.d0

r12 = r1 + r2 + h_12  !=r_centers

Uatt = -hamaker_constant/(6.d00) * (                      &
       (2*r1*r2)/(r12**2.-(r1+r2)**2.) +                  &
       (2*r1*r2)/(r12**2.-(r1-r2)**2.) +                  &
        log((r12**2.-(r1+r2)**2.)/(r12**2.-(r1-r2)**2.))  )

Urep = hamaker_constant*sigma**6./(37800.*r12) * (                           &
       ( r12**2.-7.*r12*(r1+r2)+6*(r1**2.+7.*r1*r2+r2**2.))/(r12-r1-r2)**7.  &
      +( r12**2.+7.*r12*(r1+r2)+6*(r1**2.+7.*r1*r2+r2**2.))/(r12+r1+r2)**7.  &
      -( r12**2.+7.*r12*(r1-r2)+6*(r1**2.-7.*r1*r2+r2**2.))/(r12+r1-r2)**7.  &
      -( r12**2.-7.*r12*(r1-r2)+6*(r1**2.-7.*r1*r2+r2**2.))/(r12-r1+r2)**7.  )

return
!------------------------------------------------------------------------------------------------------!
end subroutine hamaker_sphere_sphere

subroutine hamaker_well_point_plate(h12, rc, cc, Uatt)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: h12, rc, cc
real(8), intent(out) :: Uatt
!------------------------------------------------------------------------------------------------------!
!***SI units are used everywhere***
if (h12 < rc) then
    Uatt = cc * ((h12**3 - rc**3) / 3.d0 + rc**2 * (rc - h12))
else
    Uatt = 0.d0
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine hamaker_well_point_plate

subroutine hamaker_well_point_sphere(h12, rc, cc, Rsol, Uatt)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: h12, rc, cc, Rsol
real(8), intent(out) :: Uatt
real(8)              :: r_center, Rrc
!------------------------------------------------------------------------------------------------------!
!***SI units are used everywhere***
r_center = h12 + Rsol
Rrc = MIN(r_center+Rsol, rc)

if (h12 < rc) then
    Uatt = cc / r_center * (                  &
               +0.5d0 * Rrc**2 * Rsol**2      &
               -0.5d0 * Rrc**2 * r_center**2  &
               +2.d0/3.d0 * Rrc**3 * r_center &
               -0.25d0 * Rrc**4               &
               -0.5d0 * h12**2 * Rsol**2      &
               +0.5d0 * h12**2 * r_center**2  &
               -2.d0/3.d0 * h12**3 * r_center &
               +0.25d0 * h12**4               &
             )
else
    Uatt = 0.d0
endif

return
!------------------------------------------------------------------------------------------------------!
end subroutine hamaker_well_point_sphere

subroutine square_well_potential(xx, sigma, Uwell, Uatt)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: xx, sigma, Uwell
real(8), intent(out) :: Uatt
!------------------------------------------------------------------------------------------------------!
if (xx < sigma) then
    Uatt = Uatt + Uwell
else
    Uatt = 0.d0
endif
end subroutine square_well_potential

subroutine ramp_potential(xx, sigma, Uwell, Uatt)
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
real(8), intent(in)  :: xx, sigma, Uwell
real(8), intent(out) :: Uatt
!------------------------------------------------------------------------------------------------------!
if (xx < sigma) then
    Uatt = Uwell * (sigma - xx) / sigma
else
    Uatt = 0.d0
endif
end subroutine ramp_potential

!--------------------------------------------------------------------!
end module force_fields
