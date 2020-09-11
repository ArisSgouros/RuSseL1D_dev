module integr_coeffs
!----------------------------------------------------------------------------------------------------------!
contains
!----------------------------------------------------------------------------------------------------------!
    subroutine rectangle_rule(ds, ns, coeff)
    !------------------------------------------------------------------------------!
    implicit none
    !------------------------------------------------------------------------------!
    integer, intent(in) :: ns

    real(8), intent(in), dimension(0:ns)  :: ds
    real(8), intent(out), dimension(0:ns) :: coeff
    !------------------------------------------------------------------------------!
    coeff = ds
    coeff(0) = ds(1)
    return
    !------------------------------------------------------------------------------!
    end subroutine rectangle_rule

    !
    ! The Simpson's rule has been modified to work with nonuniform meshes
    !
    subroutine modsimpson_rule(ds, ns, coeff)
    !------------------------------------------------------------------------------!
    implicit none
    !------------------------------------------------------------------------------!
    integer, intent(in) :: ns
    integer             :: n

    real(8), intent(in), dimension(0:ns)  :: ds
    real(8), intent(out), dimension(0:ns) :: coeff
    real(8), dimension(0:ns)              :: x
    !------------------------------------------------------------------------------!
    x(0) = 0.d0
    do n = 1, ns
       x(n) = x(n-1) + ds(n)
    enddo

    coeff = 0.d0

    do n = 1, ns-1, 2
       coeff(n-1) = coeff(n-1) + ( x(n+1)-x(n-1) )*( 2*x(n-1)+x(n+1)-3*x(n) ) &
&                              / ( 6.d0 * ( x(n-1) - x(n) ))
    
       coeff(n)   = coeff(n)   + ( x(n-1) - x(n+1) )**3.d0 &
&                              / ( 6.d0 * ( x(n) - x(n-1) ) * (x(n) - x(n+1) ) )
    
       coeff(n+1) = coeff(n+1) + ( x(n+1)-x(n-1) )*( 2*x(n+1)+x(n-1)-3*x(n) ) &
&                              / ( 6.d0 * ( x(n+1) - x(n) ))
    enddo

    return
    !------------------------------------------------------------------------------!
    end subroutine modsimpson_rule
!------------------------------------------------------------------------------!
end module integr_coeffs
