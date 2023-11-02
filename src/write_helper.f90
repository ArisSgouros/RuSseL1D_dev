!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

module write_helper
!--------------------------------------------------------------------!
  Implicit none
!--------------------------------------------------------------------!
  integer, parameter :: PLOG = 0
  integer, parameter :: PSCREEN = 1
  integer, parameter :: PBOTH = 2

  character(len=100) :: MESSAGE
!--------------------------------------------------------------------!
contains

  subroutine output(TARGET, iow, MESSAGE)
    character(len=100) :: MESSAGE
    integer            :: TARGET
    integer            :: iow

    if (TARGET .eq. 0) then
      write (iow, *) MESSAGE
    elseif (TARGET .eq. 1) then
      write (6, *) MESSAGE
    elseif (TARGET .eq. 2) then
      write (iow, *) MESSAGE
      write (6, *) MESSAGE
    else
      write (6, *) "WRITE_HELPER: WRONG TARGET TO output FUNCTION"
      write (6, *) "    CHOOSE BETWEEN 0, 1 and 2"
    end if
  end subroutine output

  function adjl(string, length) result(r)
    character(len=*)      :: string
    integer               :: length
    character(len=length) :: r
    r = adjustl(string)
  end function adjl
!--------------------------------------------------------------------!
end module write_helper
