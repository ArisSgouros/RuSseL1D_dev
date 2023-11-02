!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine export_q(q_final, ns, nx, rs, rx, chain_type)
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: tt, ii

integer, intent(in) :: ns, nx

real(8), intent(in), dimension(0:ns)      :: rs
real(8), intent(in), dimension(0:nx)      :: rx
real(8), intent(in), dimension(0:nx,0:ns) :: q_final

character(6), intent(in) :: chain_type
character(30)            :: filename = ''
!----------------------------------------------------------------------------------------------------------!
write(filename,'("o.q_",A6)') chain_type
open (file=filename,unit=366)

write (366,'(A20)',advance='no') 'r'
do tt=0,ns
    write (366,'(E20.9)',advance='no') rs(tt)
enddo
write (366,*)
do ii=0,nx
    write (366,'(E20.9)',advance='no') rx(ii)
    do tt=0,ns
       write (366,'(E20.9)',advance='no') q_final(ii,tt)
    enddo
    write (366,*)
enddo

close(366)

return
!----------------------------------------------------------------------------------------------------------!
end subroutine export_q
