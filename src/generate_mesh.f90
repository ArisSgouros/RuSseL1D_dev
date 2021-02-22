subroutine generate_mesh(integr_scheme, discret_scheme, symmetric, domain_size, nu, du, ru, coeff)
!------------------------------------------------------------------------------------------------------!
use flags
use constants
use integr_coeffs
!------------------------------------------------------------------------------------------------------!
implicit none
!------------------------------------------------------------------------------------------------------!
integer             :: ii
integer, intent(in) :: integr_scheme, discret_scheme, nu
logical, intent(in) :: symmetric

real(8), intent(in)                   :: domain_size
real(8), intent(out), dimension(0:nu) :: du, ru, coeff
!------------------------------------------------------------------------------------------------------!
if (discret_scheme.eq.F_uniform) then
    du = domain_size / dble(nu)
    du(0) = 0.d0
    do ii = 0, nu
        ru(ii) = dble(ii) * du(ii)
    enddo
elseif (discret_scheme.eq.F_nonuniform) then
    du(0)=0.d0
    do ii = 1, nu
        if (symmetric) then
            ru(ii) = domain_size * 0.5d0 * (1.d0 - DCOS(pi * (dble(ii)) /  dble(nu)))         !symmetric scheme
        else
            ru(ii) = domain_size *         (1.d0 - DCOS(pi * (dble(ii)) / (dble(nu) * 2.d0))) !asymmetric scheme
        endif

        du(ii) = ru(ii) - ru(ii-1)
    enddo
endif

if (integr_scheme.eq.F_simpson_rule)   call modsimpson_rule(du, nu, coeff)
if (integr_scheme.eq.F_rectangle_rule) call rectangle_rule(du, nu, coeff)
!------------------------------------------------------------------------------------------------------!
end subroutine generate_mesh
