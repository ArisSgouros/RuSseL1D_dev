!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine generate_mesh(integr_scheme, discret_scheme, symmetric, domain_size, nu, du, ru, coeff)
!------------------------------------------------------------------------------------------------------!
  use flags, only: F_uniform, F_nonuniform, F_simpson_rule, F_rectangle_rule
  use constants, only: pi
  use integr_coeffs, only: modsimpson_rule, rectangle_rule
!------------------------------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------------------------------!
  integer             :: ii
  integer, intent(in) :: integr_scheme, discret_scheme, nu
  logical, intent(in) :: symmetric

  real(8), intent(in)                   :: domain_size
  real(8), intent(out), dimension(0:nu) :: du, ru, coeff
!------------------------------------------------------------------------------------------------------!
  if (discret_scheme .eq. F_uniform) then
    du = domain_size/dble(nu)
    du(0) = 0.d0
    do ii = 0, nu
      ru(ii) = dble(ii)*du(ii)
    end do
  elseif (discret_scheme .eq. F_nonuniform) then
    du(0) = 0.d0
    do ii = 1, nu
      if (symmetric) then
        ru(ii) = domain_size*0.5d0*(1.d0 - DCOS(pi*(dble(ii))/dble(nu)))         !symmetric scheme
      else
        ru(ii) = domain_size*(1.d0 - DCOS(pi*(dble(ii))/(dble(nu)*2.d0))) !asymmetric scheme
      end if

      du(ii) = ru(ii) - ru(ii - 1)
    end do
  end if

  if (integr_scheme .eq. F_simpson_rule) call modsimpson_rule(du, nu, coeff)
  if (integr_scheme .eq. F_rectangle_rule) call rectangle_rule(du, nu, coeff)
!------------------------------------------------------------------------------------------------------!
end subroutine generate_mesh
