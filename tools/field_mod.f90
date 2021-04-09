program mod_field_bin
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer :: kk, nx_old, nx_mod, nx_max

logical :: FILE_EXISTS

character(200) :: field_old_filename = "field.out.bin"
character(200) :: field_mod_filename = "field.mod.bin"

real(8) :: lx_old, lx_mod
real(8) :: dx_old, dx_mod

real(8), allocatable, dimension(:) :: wa_old, wa_mod

character(200) :: argument
!----------------------------------------------------------------------------------------------------------!

call GET_COMMAND_ARGUMENT(1,argument)
read( argument, * ) lx_old
call GET_COMMAND_ARGUMENT(2,argument)
read( argument, * ) lx_mod
call GET_COMMAND_ARGUMENT(3,field_old_filename)
call GET_COMMAND_ARGUMENT(4,field_mod_filename)

write(*,*)'Old length of the domain:',lx_old
write(*,*)'New length of the domain:',lx_mod
write(*,*)'Old field filename:',field_old_filename
write(*,*)'New field filename:',field_mod_filename

!lx_old = 70.0
dx_old = 0.5
nx_old = 2 * int(0.5d0*lx_old/dx_old)

!lx_mod = 71.0
dx_mod = dx_old
nx_mod = 2 * int(0.5d0*lx_mod/dx_mod)

nx_max = max(nx_old,nx_mod)

allocate(wa_old(0:nx_old),wa_mod(0:nx_mod))

wa_old = 0.d0
wa_mod = 0.d0

INQUIRE(FILE=field_old_filename, EXIST=FILE_EXISTS)

if (FILE_EXISTS) then
    open(unit=21, file = field_old_filename, form = 'unformatted')
    read (21) wa_old
    close(21)
else
    write(*,'(''File '',A200,'' does not exist'')') field_old_filename
    STOP
endif

! increase/decrease domain
do kk = 0, nx_max / 2
   wa_mod(kk)        = wa_old(kk)
   wa_mod(nx_mod-kk) = wa_old(nx_old-kk)
end do

! set field to unity in a set range
!do kk = 0, nx_old
!   if ( 100 < kk .and. kk < 222) then
!      wa_mod(kk) = 0.0015
!   else
!      wa_mod(kk) = wa_old(kk)
!   endif
!end do

do kk = 0, nx_max
   write(*,*)kk,wa_old(kk),wa_mod(kk)
end do

open(unit=21, file = field_mod_filename, form = 'unformatted')
write (21) wa_mod
close(21)

!----------------------------------------------------------------------------------------------------------!
end program mod_field_bin
