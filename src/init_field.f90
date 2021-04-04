subroutine init_field()
!----------------------------------------------------------------------------------------------------------!
use parser_vars, only: read_field, field_in_filename
use arrays,      only: wa_ifc, wa_ifc_backup
!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
if (read_field) then
    open(unit=21, file = field_in_filename, form = 'unformatted')
    read (21) wa_ifc
    close(21)
endif

!store the initial field for backup
wa_ifc_backup = wa_ifc

return
!----------------------------------------------------------------------------------------------------------!
end subroutine init_field
