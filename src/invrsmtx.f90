    subroutine invrsmtx(U, n, Uinv)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
      ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
      integer, intent(in)::n
      integer :: i, j, k, l, m, irow
      real(8):: big, dum
      real(8), dimension(1:n, 1:n)::temp
      real(8), dimension(1:n, 1:n)::U
      real(8), dimension(1:n, 1:n)::Uinv
      do i = 1, n
      do j = 1, n
        temp(i, j) = U(i, j)

      end do
      end do
      !build the identity matrix
      do i = 1, n
        do j = 1, n
          Uinv(i, j) = 0.0
        end do
        Uinv(i, i) = 1.0
      end do

      do i = 1, n! this is the big loop over all the columns of a(n,n)
        ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
        ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
        big = temp(i, i)
        do j = i, n
        if (temp(j, i) .gt. big) then
          big = temp(j, i)

          irow = j
        end if
        end do
        ! interchange lines i with irow for both a() and b() matrices
        if (big .gt. temp(i, i)) then
          do k = 1, n
            dum = temp(i, k)
            temp(i, k) = temp(irow, k)
            temp(irow, k) = dum
            dum = Uinv(i, k) ! matrix b()
            Uinv(i, k) = Uinv(irow, k)
            Uinv(irow, k) = dum
          end do
        end if
        ! divide all entries in line i from a(i,j) by the value a(i,i);
        ! same operation for the identity matrix
        dum = temp(i, i)
        do j = 1, n
          temp(i, j) = temp(i, j)/dum
          Uinv(i, j) = Uinv(i, j)/dum
        end do
! make zero all entries in the column a(j,i); same operation for indent()
        do j = i + 1, n
          dum = temp(j, i)
          do k = 1, n
            temp(j, k) = temp(j, k) - dum*temp(i, k)
            Uinv(j, k) = Uinv(j, k) - dum*Uinv(i, k)
          end do
        end do
      end do

! substract appropiate multiple of row j from row j-1
      do i = 1, n - 1
      do j = i + 1, n
        dum = temp(i, j)
        do l = 1, n
          temp(i, l) = temp(i, l) - dum*temp(j, l)
          Uinv(i, l) = Uinv(i, l) - dum*Uinv(j, l)
        end do
      end do
      end do

    end subroutine invrsmtx

