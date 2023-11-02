!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_phi_seg(iseg, chainlen, coeff_ns, ns, nx, rx, qfinal, qfinal_aux, chain_type)
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  integer, intent(in) :: iseg, ns, nx
  integer             :: kk, ii

  real(8), intent(in), dimension(0:ns)      :: coeff_ns
  real(8), intent(in), dimension(0:nx)      :: rx
  real(8), intent(in), dimension(0:nx, 0:ns) :: qfinal, qfinal_aux
  real(8), intent(in)                       :: chainlen
  real(8), dimension(0:nx)                  :: phi_seg

  character(6), intent(in) :: chain_type
  character(40)            :: filename = ''
!----------------------------------------------------------------------------------------------------------!
  phi_seg = 0.d0

  write (filename, '("o.phi_segs_",A6)') chain_type
  open (unit=122, file=filename)

  if (iseg .ne. -1) then
    !start profile post-processing
    do kk = 0, nx
      phi_seg(kk) = (qfinal(kk, ns - iseg)*qfinal_aux(kk, iseg))/chainlen
    end do

    write (122, '(2(2X,A16))') 'r', "phi_seg(r)"
    do kk = 0, nx
      write (122, '(3(2X,E16.9))') rx(kk), phi_seg(kk)
    end do
    close (unit=122)
  end if

  if (iseg .eq. -1) then
    write (122, '(2X,A16)', advance='no') 'r'
    do ii = 0, ns
      write (122, '(2X,I16)', advance='no') ii
    end do
    write (122, *)
    do ii = 0, ns
      write (122, '(2X,E16.9)', advance='no') coeff_ns(ii)
    end do
    write (122, *)
    do kk = 0, nx
      write (122, '(2X,E16.9)', advance='no') rx(kk)
      do ii = 0, ns
        write (122, '(2X,E16.9)', advance='no') (qfinal(kk, ii)*qfinal_aux(kk, ns - ii))/chainlen
      end do
      write (122, *)
    end do
  end if

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_phi_seg
