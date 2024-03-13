!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_field()
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: read_field, field_in_filename
  use arrays, only: wa_ifc_kd1, wa_ifc_backup_kd1
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
  if (read_field) then
    open (unit=21, file=field_in_filename, form='unformatted')
    read (21) wa_ifc_kd1
    close (21)
  end if

!store the initial field for backup
  wa_ifc_backup_kd1 = wa_ifc_kd1

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_field
