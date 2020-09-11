subroutine tridag(n, a, b, c, r, u)
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: n
integer             :: k

real(8), intent(inout), dimension(n) :: a, b, c
real(8), intent(out), dimension(n)   :: r, u
!----------------------------------------------------------------------------------------------------------!
!Thomas Algorithm
!See Page 301 Numerical methods for engineers
!Thomas algorithm is a efficient method to solve tridiagonal system
!       /                                                  \ /      \      /       \
!      | b(1)    c(1)                                       || u(1)  |     | r(1)  |
!      | a(2)    b(2)    c(2)                               || u(2)  |     | r(2)  |
!      |         a(3)    b(3)    c(3)                       || u(3)  |     | r(3)  |
!      |                .      .      .                     ||   .   | =   | .     |
!      |                     a(n-2)  b(n-2)  c(n-2)         || u(n-2)|     | r(n-2)|
!      |                             a(n-1)  b(n-1)  c(n-1) || u(n-1)|     | r(n-1)|
!      |                                      a(n)     b(n) || u(n)  |     | r(n)  |
!      \                                                    /\       /      \      /
!As convectional LU decomposition the algorithm consists of three steps

!decomposition
do k = 2, n
    a(k) = a(k)/b(k-1)
    b(k) = b(k)-a(k)*c(k-1)
enddo

!forward substitution
do k = 2, n
    r(k) = r(k)-a(k)*r(k-1)
enddo

!back substitution
u(n) = r(n)/b(n)
do k = n-1, 1, -1
    u(k) = (r(k)-c(k)*u(k+1) )/b(k)
enddo

return
!----------------------------------------------------------------------------------------------------------!
end subroutine tridag
