!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_field()
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: read_field, field_in_filename, mxa_kind, mxb_kind, glo_kind, ghi_kind
  use arrays, only: wa_ifc_kd1, wa_ifc_backup_kd1, wa_ifc_kd2, wa_ifc_backup_kd2, &
                   & wa_ifc_mxa, wa_ifc_mxb, wa_ifc_glo, wa_ifc_ghi
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!
! TODO: add option for reading both kd1 and kd2 fields
  if (read_field) then
    open (unit=21, file=field_in_filename, form='unformatted')
    read (21) wa_ifc_kd1
    close (21)
  end if

!store the initial field for backup
  wa_ifc_backup_kd1 = wa_ifc_kd1

! TODO: revert tentative initial assignment from kd1 field
  wa_ifc_kd2 = wa_ifc_kd1
  wa_ifc_backup_kd2 = wa_ifc_kd1

! TODO: add option for randomized initial fields

  if (mxa_kind==1) wa_ifc_mxa = wa_ifc_kd1
  if (mxb_kind==1) wa_ifc_mxb = wa_ifc_kd1
  if (glo_kind==1) wa_ifc_glo = wa_ifc_kd1
  if (ghi_kind==1) wa_ifc_ghi = wa_ifc_kd1
  if (mxa_kind==2) wa_ifc_mxa = wa_ifc_kd2
  if (mxb_kind==2) wa_ifc_mxb = wa_ifc_kd2
  if (glo_kind==2) wa_ifc_glo = wa_ifc_kd2
  if (ghi_kind==2) wa_ifc_ghi = wa_ifc_kd2

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_field
