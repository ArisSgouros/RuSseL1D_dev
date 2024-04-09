!RuSseL1D - Copyright (C) 2021 C. J. Revelas, A. P. Sgouros, A. T. Lakkas
!
!See the LICENSE file in the root directory for license information.

subroutine init_field()
!----------------------------------------------------------------------------------------------------------!
  use parser_vars, only: read_field, random_field, random_field_magn, &
                        & mxa_kind, mxb_kind, glo_kind, ghi_kind, &
                        & exist_kd1, exist_kd2, &
                        & nx
  use arrays, only: wa_ifc_kd1, wa_ifc_backup_kd1, wa_ifc_kd2, wa_ifc_backup_kd2, &
                   & wa_ifc_mxa, wa_ifc_mxb, wa_ifc_glo, wa_ifc_ghi
!----------------------------------------------------------------------------------------------------------!
  implicit none
!----------------------------------------------------------------------------------------------------------!

  integer :: jj

  wa_ifc_kd1 = 0.d0
  wa_ifc_kd2 = 0.d0

  if (read_field) then
    if (exist_kd1) then
      open (unit=21, file="in.field.bin", form='unformatted')
      read (21) wa_ifc_kd1
      close (21)
    end if
    if (exist_kd2) then
      open (unit=21, file="in.field_kd2.bin", form='unformatted')
      read (21) wa_ifc_kd2
      close (21)
    end if
  end if

  if (random_field) then
    if (exist_kd1) then
      do jj = 0, nx
        wa_ifc_kd1(jj) = random_field_magn*(1.0d0 - 2.0d0*RAND())
      end do
    end if
    if (exist_kd2) then
      do jj = 0, nx
        wa_ifc_kd2(jj) = random_field_magn*(1.0d0 - 2.0d0*RAND())
      end do
    end if
  end if

!store the initial field for backup
  wa_ifc_backup_kd1 = wa_ifc_kd1
  wa_ifc_backup_kd2 = wa_ifc_kd2

  if (mxa_kind == 1) wa_ifc_mxa = wa_ifc_kd1
  if (mxb_kind == 1) wa_ifc_mxb = wa_ifc_kd1
  if (glo_kind == 1) wa_ifc_glo = wa_ifc_kd1
  if (ghi_kind == 1) wa_ifc_ghi = wa_ifc_kd1
  if (mxa_kind == 2) wa_ifc_mxa = wa_ifc_kd2
  if (mxb_kind == 2) wa_ifc_mxb = wa_ifc_kd2
  if (glo_kind == 2) wa_ifc_glo = wa_ifc_kd2
  if (ghi_kind == 2) wa_ifc_ghi = wa_ifc_kd2

  return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_field
