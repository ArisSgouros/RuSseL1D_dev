!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine solver_gelim(nx, sk, r1, u)
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer, intent(in) :: nx
integer             :: i, j, k, l

real(8), dimension(nx, nx), intent(inout) :: sk
real(8), dimension(nx), intent(inout)     :: r1
real(8), dimension(nx), intent(out)       :: u
real(8), dimension(nx)                    :: s
real(8)                                   :: c, pivot, store
!----------------------------------------------------------------------------------------------------------!
! step 1: begin forward elimination
do k=1, nx-1

! step 2: "scaling"
! s(i) will have the largest element from row i 
  do i=k,nx                       ! loop over rows
    s(i) = 0.0
    do j=k,nx                    ! loop over elements of row i
      s(i) = max(s(i),abs(sk(i,j)))
    end do
  end do

! step 3: "pivoting 1" 
! find a row with the largest pivoting element
  pivot = abs(sk(k,k)/s(k))
  l = k
  do j=k+1,nx
    if(abs(sk(j,k)/s(j)) > pivot) then
      pivot = abs(sk(j,k)/s(j))
      l = j
    endif
  enddo

! Check if the system has a sigular matrix
  if(pivot == 0.0) then
    write(6,*) " The matrix is singular.. "
    return
  endif

! step 4: "pivoting 2" interchange rows k and l (if needed)
if (l /= k) then
  do j=k,nx
     store = sk(k,j)
     sk(k,j) = sk(l,j)
     sk(l,j) = store
  enddo
  store = r1(k)
  r1(k) = r1(l)
  r1(l) = store
endif

! step 5: the elimination (after scaling and pivoting)
   do i=k+1,nx
      c=sk(i,k)/sk(k,k)
      sk(i,k) = 0.0
      r1(i)=r1(i)- c*r1(k)
      do j=k+1,nx
         sk(i,j) = sk(i,j)-c*sk(k,j)
      enddo
   enddo
enddo

! step 6: back substiturion 
u(nx) = r1(nx)/sk(nx,nx)
do i=nx-1,1,-1
   c=0.0
   do j=i+1,nx
     c= c + sk(i,j)*u(j)
   enddo 
   u(i) = (r1(i)- c)/sk(i,i)
enddo

return
!----------------------------------------------------------------------------------------------------------!
end subroutine solver_gelim
