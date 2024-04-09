pure function invrsmtx(U,U2) result(Uinv,Uinv2)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    complex(wp), intent(in) :: U(2,2),U2(2,2)   !! Matrix
    complex(wp),intent(out)             :: Uinv(2,2),Uinv2 !! Inverse matrix
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(U(1,1)*U(2,2) - U(1,2)*U(2,1))

    ! Calculate the inverse of the matrix
    Uinv(1,1) = +detinv * A(2,2)
    Uinv(2,1) = -detinv * A(2,1)
    Uinv(1,2) = -detinv * A(1,2)
    Uinv(2,2) = +detinv * A(1,1)
    ! Calculate the inverse determinant of the matrix2
    detinv2= 1/(U2(1,1)*U2(2,2) - U2(1,2)*U2(2,1))

    ! Calculate the inverse of the matrix
    Uinv2(1,1) = +detinv2 * U2(2,2)
    Uinv2(2,1) = -detinv2 * U2(2,1)
    Uinv2(1,2) = -detinv2 * U2(1,2)
    Uinv2(2,2) = +detinv2 * U2(1,1)
    
  end function
