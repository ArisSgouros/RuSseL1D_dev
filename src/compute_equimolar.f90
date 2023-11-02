!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine compute_equimolar(nx, rx, coeff_x, layer_area, phi, chain_type)
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  real(8), parameter :: tol_phi_0 = 1e-8
  real(8), parameter :: tol_phi_bulk = 1e-3

  integer, intent(in) :: nx
  integer             :: kk, jj

  real(8), intent(in), dimension(0:nx) :: phi, rx, coeff_x, layer_area
  real(8)                              :: phi_a, phi_b, phi05
  real(8)                              :: rx_phi05_lo, rx_phi05_hi, rx_equi_lo, rx_equi_hi
  real(8)                              :: slope, intercept

  logical :: is_lo, is_hi, is_bulk

  real(8) :: integ_old, integ_new, phi_phase, phi_excess

  character(6), intent(in) :: chain_type
  character(40)            :: filename = ''
!----------------------------------------------------------------------------------------------------------!

  phi_b = 0.d0 ! gas phase

  phi_a = 0.d0 ! bulk phase
  do kk = 0, nx
    phi_a = max(phi_a, phi(kk))
  end do
! in case the layer does not feature a bulk region set it to the maximum phi
  phi_a = min(1.d0, phi_a)

! The phi05 is traced at the region where phi = 0.5 * phi_a
  phi05 = 0.5d0*phi_a

  is_lo = (dabs(phi(0)) < tol_phi_0)
  is_hi = (dabs(phi(nx)) < tol_phi_0)
  is_bulk = (dabs(phi_a - 1.d0) < tol_phi_bulk)

  rx_phi05_lo = 0.d0
  rx_equi_lo = 0.d0
  if (is_lo) then
    ! Calculate phi05_lo
    do kk = 1, nx
      if (phi(kk) > phi05 .and. phi05 > phi(kk - 1)) then
        slope = (phi(kk) - phi(kk - 1))/(rx(kk) - rx(kk - 1))
        intercept = -rx(kk)*slope + phi(kk)
        rx_phi05_lo = (phi05 - intercept)/slope
        exit
      end if
    end do

    ! Calculate the r_equi_lo if bulk region exists
    if (is_bulk) then
      integ_old = 0.d0
      do jj = 1, nx
        integ_new = 0.d0
        do kk = 1, nx
          if (rx(kk) < rx(jj)) then
            phi_phase = phi_b
          else
            phi_phase = phi_a
          end if
          phi_excess = phi(kk) - phi_phase

          ! address cumulative numerical errors for high curvatures (phi =0.999..)
          if (dabs(phi(kk) - phi_a) < tol_phi_bulk) exit

          ! surface excess up to phi = phi_a
          if (phi(kk) > phi_a) exit

          integ_new = integ_new + coeff_x(kk)*layer_area(kk)*phi_excess
        end do

        if (integ_new > 0.d0) then
          slope = (integ_new - integ_old)/(rx(jj) - rx(jj - 1))
          intercept = -rx(jj)*slope + integ_new
          rx_equi_lo = (0.d0 - intercept)/slope
          exit
        else
          integ_old = integ_new
        end if
      end do
    end if
  end if

! Calculate phi05_hi
  rx_phi05_hi = 0.d0
  rx_equi_hi = 0.d0
  if (is_hi) then
    do kk = nx - 1, 0, -1
      if (phi(kk) > phi05 .and. phi05 > phi(kk + 1)) then
        slope = (phi(kk) - phi(kk + 1))/(rx(kk) - rx(kk + 1))
        intercept = -rx(kk)*slope + phi(kk)
        rx_phi05_hi = (phi05 - intercept)/slope
        exit
      end if
    end do

    ! Calculate the r_equi_hi if bulk region exists
    if (is_bulk) then
      integ_old = 0.d0
      do jj = nx - 1, 1, -1
        integ_new = 0.d0
        do kk = nx - 1, 0, -1
          if (rx(kk) >= rx(jj)) then
            phi_phase = phi_b
          else
            phi_phase = phi_a
          end if
          phi_excess = phi(kk) - phi_phase

          ! address cumulative numerical errors for high curvatures (phi =0.999..)
          if (dabs(phi(kk) - phi_a) < tol_phi_bulk) exit

          ! surface excess up to phi = phi_a
          if (phi(kk) > phi_a) exit

          integ_new = integ_new + coeff_x(kk)*layer_area(kk)*phi_excess
        end do

        if (integ_new > 0.d0) then
          slope = (integ_new - integ_old)/(rx(jj - 1) - rx(jj))
          intercept = -rx(jj - 1)*slope + integ_new
          rx_equi_hi = (0.d0 - intercept)/slope
          exit
        else
          integ_old = integ_new
        end if
      end do
    end if

  end if

  write (filename, '("o.equimolar_",A6)') chain_type

  open (file=filename, unit=32)
  write (32, '(4(A16))') "side", "r_phi05", "r_phi_equi", "rel_diff"
  write (32, '(A16,3(E16.6))') "lo", rx_phi05_lo, rx_equi_lo, rx_equi_lo - rx_phi05_lo
  write (32, '(A16,3(E16.6))') "hi", rx_phi05_hi, rx_equi_hi, rx_phi05_hi - rx_equi_hi
  close (32)

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine compute_equimolar
