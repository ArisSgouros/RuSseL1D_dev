subroutine parser()
!----------------------------------------------------------------------------------------------------------!
use flags
use parser_vars
use constants
use eos
use write_helper

!----------------------------------------------------------------------------------------------------------!
implicit none
!----------------------------------------------------------------------------------------------------------!
integer        :: ior = 638

character(200) :: line
character(15)  :: input_filename = 'input.in.txt'
character(12)  :: field_filename = 'field.in.bin'
character(100) :: ERROR_MESSAGE

logical :: FILE_EXISTS

integer :: Reason, wall_itype

logical :: log_system_geometry             = .false.
logical :: log_sphere_radius               = .false.
logical :: log_lx                          = .false.
logical :: log_dx                          = .false.
logical :: log_temperature                 = .false.
logical :: log_pressure                    = .false.

logical :: log_characteristic_ratio        = .false.
logical :: log_bond_length                 = .false.
logical :: log_monomer_mass                = .false.
logical :: log_mass_density                = .false.

logical :: log_number_of_iterations        = .false.
logical :: log_thermo_every                = .false.
logical :: log_compute_every               = .false.
logical :: log_field_every                 = .false.
logical :: log_check_stability_every       = .false.
logical :: log_max_wa_error                = .false.
logical :: log_fraction_of_new_field       = .false.
logical :: log_read_field                  = .false.

logical :: log_out_phi                     = .false.
logical :: log_out_field                   = .false.
logical :: log_out_q                       = .false.
logical :: log_out_chainshape              = .false.
logical :: log_out_ads_free                = .false.
logical :: log_out_end_middle              = .false.
logical :: log_out_brush_thickness         = .false.

logical :: log_edwards_solver              = .false.
logical :: log_spatial_discret_scheme      = .false.
logical :: log_contour_discret_scheme      = .false.
logical :: log_spatial_integr_scheme       = .false.
logical :: log_contour_integr_scheme       = .false.

logical :: log_wall_type                   = .false.
logical :: log_wall_coeffs_hamaker         = .false.
logical :: log_wall_coeffs_square_well     = .false.
logical :: log_wall_coeffs_ramp            = .false.
logical :: log_wall_pos                    = .false.
logical :: log_wall_side                   = .false.
logical :: log_wall_pos_auto               = .false.

logical :: log_matrix_exist                = .false.
logical :: log_chain_length_matrix         = .false.
logical :: log_ds_ave_matrix               = .false.
logical :: log_r_critical                  = .false.

logical :: log_grafted_lo_exist            = .false.
logical :: log_chain_length_grafted_lo     = .false.
logical :: log_ds_ave_grafted_lo           = .false.
logical :: log_gdens_lo                    = .false.

logical :: log_grafted_hi_exist            = .false.
logical :: log_chain_length_grafted_hi     = .false.
logical :: log_ds_ave_grafted_hi           = .false.
logical :: log_gdens_hi                    = .false.

logical :: log_position_of_grafted         = .false.

logical :: log_hi_BC_of_matrix             = .false.
logical :: log_lo_BC_of_matrix             = .false.
logical :: log_hi_BC_of_grafted            = .false.
logical :: log_lo_BC_of_grafted            = .false.

logical :: log_eos_type                    = .false.
logical :: log_eos_coeffs                  = .false.
logical :: log_influence_param             = .false.
logical :: log_real_influence_param        = .false.
!----------------------------------------------------------------------------------------------------------!
wall_hamaker     = .false.
wall_square_well = .false.
wall_ramp        = .false.
wall_vacuum      = .false.
wall_hybrid      = .false.
!----------------------------------------------------------------------------------------------------------!

INQUIRE(FILE=input_filename, EXIST=FILE_EXISTS)

if (FILE_EXISTS) then
    open(unit=ior, file = input_filename)
else
    write(ERROR_MESSAGE,'(''File '',A15,'' does not exist!'')')input_filename
    write(iow,*)ERROR_MESSAGE
    write(*  ,*)ERROR_MESSAGE
    STOP
endif

do
    read(ior,'(A100)',IOSTAT=Reason) line

    if (Reason > 0)  then
        write(iow,*)"Problem with the input file!"
        write(*  ,*)"Problem with the input file!"
        STOP
    elseif (Reason < 0) then
        write(*,*)"Input parameter file was read successfully!"
        exit
    else
        ! system setup
        if (index(line,'! domain geometry') > 0) then
            read(line,'(I9)') geometry
            log_system_geometry = .true.
        elseif (index(line,'! domain lx') > 0) then
            read(line,*) lx
            log_lx = .true.
        elseif (index(line,'! domain dx') > 0) then
            read(line,*) dx_ave
            log_dx = .true.
        elseif (index(line,'! domain sphere_radius') > 0) then
            read(line,'(E16.9)') sphere_radius
            log_sphere_radius = .true.
        elseif (index(line,'! system temperature') > 0) then
            read(line,'(E16.9)') Temp
            log_temperature = .true.
        elseif (index(line,'! system pressure') > 0) then
            read(line,'(E16.9)') Pressure
            log_pressure = .true.
    
        ! polymer parameters
        elseif (index(line,'! polymer C_inf') > 0) then
            read(line,'(E16.9)') CN
            log_characteristic_ratio = .true.
        elseif (index(line,'! polymer bond_length') > 0) then
            read(line,'(E16.9)') bond_length
            log_bond_length = .true.
        elseif (index(line,'! polymer monomer_mass') > 0) then
            read(line,'(E16.9)') mon_mass
            log_monomer_mass = .true.
        elseif (index(line,'! polymer mass_density') > 0) then
            read(line,'(E16.9)') rho_mass_bulk
            log_mass_density = .true.

        ! solution parameters
        elseif (index(line,'! field iterations') > 0) then
            read(line,'(I10)') max_iter
            log_number_of_iterations = .true.
        elseif (index(line,'! thermo every') > 0) then
            read(line,'(I10)') thermo_every
            log_thermo_every = .true.
        elseif (index(line,'! field every') > 0) then
            read(line,'(I10)') field_every
            log_field_every = .true.
        elseif (index(line,'! compute every') > 0) then
            read(line,'(I10)') compute_every
            log_compute_every = .true.
        elseif (index(line,'! check stability every') > 0) then
            read(line,'(I10)') check_stability_every
            log_check_stability_every = .true.
        elseif (index(line,'! field read') > 0) then
            read(line,'(L10)') read_field
            log_read_field = .true.
        elseif (index(line,'! field max_error') > 0) then
            read(line,'(E16.9)') max_wa_error
            log_max_wa_error = .true.
        elseif (index(line,'! field mixing_fraction') > 0) then
            read(line,'(E16.9)') frac
            log_fraction_of_new_field = .true.
        elseif (index(line,'! edwards solver') > 0) then
            read(line,*) edwards_solver
            log_edwards_solver = .true.
        elseif (index(line,"! discret contour") > 0) then
            read(line,*) contour_discret_scheme
            log_contour_discret_scheme = .true.
        elseif (index(line,"! discret spatial") > 0) then
            read(line,*) spatial_discret_scheme
            log_spatial_discret_scheme = .true.
        elseif (index(line,"! integr contour") > 0) then
            read(line,*) contour_integr_scheme
            log_contour_integr_scheme = .true.
        elseif (index(line,"! integr spatial") > 0) then
            read(line,*) spatial_integr_scheme
            log_spatial_integr_scheme = .true.

        ! computes
        elseif (index(line,'! export phi') > 0) then
            read(line,'(L10)') out_phi
            log_out_phi = .true.
        elseif (index(line,'! export field') > 0) then
            read(line,'(L10)') out_field
            log_out_field = .true.
        elseif (index(line,'! export q') > 0) then
            read(line,'(L10)') out_q
            log_out_q = .true.
        elseif (index(line,'! export shape') > 0) then
            read(line,'(L10)') out_chainshape
            log_out_chainshape = .true.
        elseif (index(line,'! export ads_free') > 0) then
            read(line,'(L10)') out_ads_free
            log_out_ads_free = .true.
        elseif (index(line,'! export brush') > 0) then
            read(line,'(L10)') out_brush_thickness
            log_out_brush_thickness = .true.
        elseif (index(line,'! export end_middle') > 0) then
            read(line,'(L10)') out_end_middle
            log_out_end_middle = .true.

        ! wall section
        elseif (index(line,'! wall type') > 0) then
            read(line,'(I10)') wall_type
            if (wall_type == F_vacuum)      wall_vacuum = .true.
            if (wall_type == F_hamaker)     wall_hamaker = .true.
            if (wall_type == F_square_well) wall_square_well = .true.
            if (wall_type == F_ramp)        wall_ramp = .true.
            if (wall_type == F_hybrid)      wall_hybrid = .true.
            log_wall_type = .true.
        elseif (index(line,'! wall coeffs') > 0) then
            if (wall_hybrid) then
                read(line,*) wall_itype
                if (wall_itype.eq.F_hamaker) then
                    read(line,*) wall_itype, sig_pol, sig_solid, Apol, Asolid
                    wall_hamaker = .true.
                    log_wall_coeffs_hamaker = .true.
                endif
                if (wall_itype.eq.F_square_well) then
                    read(line,*) wall_itype, sigma_sq_well, A_sq_well
                    wall_square_well = .true.
                    log_wall_coeffs_square_well = .true.
                endif
                if (wall_itype.eq.F_ramp) then
                    read(line,*) wall_itype, sigma_ramp, A_ramp
                    wall_ramp = .true.
                    log_wall_coeffs_ramp = .true.
                endif
            else
                if (wall_hamaker) then
                    read(line,*) sig_pol, sig_solid, Apol, Asolid
                    log_wall_coeffs_hamaker = .true.
                endif
                if (wall_square_well) then
                    read(line,*) sigma_sq_well, A_sq_well
                    log_wall_coeffs_square_well = .true.
                endif
                if (wall_ramp) then
                    read(line,*) sigma_ramp, A_ramp
                    log_wall_coeffs_ramp = .true.
                endif
            endif
        elseif (index(line,'! wall pos set') > 0) then
            read(line,'(E16.9)') wall_pos
            log_wall_pos = .true.
        elseif (index(line,'! wall pos auto') > 0) then
            read(line,'(E16.9)') E_wall_target
            log_wall_pos_auto = .true.
        elseif (index(line,'! wall side') > 0) then
            read(line,'(I10)') wall_side
            log_wall_side = .true.
        !matrix chains
        elseif (index(line,'! matrix set') > 0) then
            read(line,'(L10)') matrix_exist
            log_matrix_exist = .true.
        elseif (index(line,'! matrix chain_length') > 0) then
            read(line,*) chainlen_matrix
            log_chain_length_matrix = .true.
        elseif (index(line,'! matrix ds') > 0) then
            read(line,*) ds_ave_matrix
            log_ds_ave_matrix = .true.
        elseif (index(line,'matrix r_adsorbed') > 0) then
            read(line,'(E16.9)') r_critical
            log_r_critical = .true.
        ! grafted lo chains
        elseif (index(line,'! grafted lo set') > 0) then
            read(line,'(L10)') grafted_lo_exist
            log_grafted_lo_exist = .true.
        elseif (index(line,'! grafted lo chain_length') > 0) then
            read(line,*) chainlen_grafted_lo
            log_chain_length_grafted_lo = .true.
        elseif (index(line,'! grafted lo ds') > 0) then
            read(line,*) ds_ave_grafted_lo
            log_ds_ave_grafted_lo = .true.
        elseif (index(line,'! grafted lo grafting_density') > 0) then
            read(line,'(E16.9)') gdens_lo
            log_gdens_lo = .true.
        ! grafted hi chains
        elseif (index(line,'! grafted hi set') > 0) then
            read(line,'(L10)') grafted_hi_exist
            log_grafted_hi_exist = .true.
        elseif (index(line,'! grafted hi chain_length') > 0) then
            read(line,*) chainlen_grafted_hi
            log_chain_length_grafted_hi = .true.
        elseif (index(line,'! grafted hi ds') > 0) then
            read(line,*) ds_ave_grafted_hi
            log_ds_ave_grafted_hi = .true.
        elseif (index(line,'! grafted hi grafting_density') > 0) then
            read(line,'(E16.9)') gdens_hi
            log_gdens_hi = .true.
        ! grafted misc
        elseif (index(line,'! grafted distance_from_solid') > 0) then
            read(line,*) graft_pos
            log_position_of_grafted = .true.
        ! boundary condition
        elseif (index(line,'! boundary_condition lo matrix') > 0) then
            read(line,'(I10)') bc_lo_matrix
            log_lo_BC_of_matrix = .true.
        elseif (index(line,'! boundary_condition hi matrix') > 0) then
            read(line,'(I10)') bc_hi_matrix
            log_hi_BC_of_matrix = .true.
        elseif (index(line,'! boundary_condition lo grafted') > 0) then
            read(line,'(I10)') bc_lo_grafted
            log_lo_BC_of_grafted = .true.
        elseif (index(line,'! boundary_condition hi grafted') > 0) then
            read(line,'(I10)') bc_hi_grafted
            log_hi_BC_of_grafted = .true.
        ! EOS
        elseif (index(line,'! EOS type') > 0) then
            read(line,'(I6)') eos_type
            log_eos_type = .true.
        elseif (index(line,'! EOS coeffs') > 0) then
            if (eos_type.eq.0) then
                read(line,*) HF_kappa_T
            else if (eos_type.eq.1) then
                read(line,*) rho_star, T_star, P_star
            end if
            log_eos_coeffs = .true.
        elseif (index(line,'! EOS influence_parameter') > 0) then
            read(line,'(E16.9)') k_gr_tilde
            log_influence_param = .true.
        elseif (index(line,'! EOS real_influence_parameter') > 0) then
            read(line,'(E16.9)') k_gr
            log_real_influence_param = .true.
        endif
    endif
enddo

close(ior)

!*******************************************************************!
!              Check the inputs of the parameter file               !
!*******************************************************************!
write(iow,'(A85)')adjl('-----------------------------------SYSTEM PARAMETERS---------------------------------',85)
write(*  ,'(A85)')adjl('-----------------------------------SYSTEM PARAMETERS---------------------------------',85)
if (log_system_geometry) then
    if (geometry.eq.F_film) then
        write(iow,'(3X,A45,A16)')adjl('System geometry:',45),adjustl('film')
        write(*  ,'(3X,A45,A16)')adjl('System geometry:',45),adjustl('film')
    elseif (geometry.eq.F_sphere) then
        write(iow,'(3X,A45,A16)')adjl('System geometry:',45),adjustl('sphere')
        write(*  ,'(3X,A45,A16)')adjl('System geometry:',45),adjustl('sphere')
    else
        write(iow,'(3X,A80,I3,'' or '',I3)')'Error: Unknown system geometry. Choose either:',F_film,F_sphere
        STOP
    endif
else
    write(iow,'(3X,A45)') 'Error: System geometry was not set'
    write(*  ,'(3X,A45)') 'Error: System geometry was not set'
    STOP
endif

if (geometry.eq.F_sphere) then
if (log_sphere_radius) then
    write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('Sphere radius:',45),sphere_radius
    write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('Sphere radius:',45),sphere_radius
else
    write(iow,'(3X,A45)')'Error: Sphere radius not found'
    write(*  ,'(3X,A45)')'Error: Sphere radius not found'
    STOP
endif
endif

if (log_lx) then
    write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('Domain in Angstroms:',45),lx
    write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('Domain in Angstroms:',45),lx
else
    write(iow,'(3X,A45)')'Error: domain length was not set'
    write(*  ,'(3X,A45)')'Error: domain length was not set'
    STOP
endif

if (log_dx) then
    nx = 2 * int(0.5d0*lx/dx_ave)
    write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('Domain discretization:',45), dx_ave
    write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('Domain discretization:',45), dx_ave
    write(iow,'(3X,A45,I16,'' nodes'')')                adjl('Nodes across the domain:',45), nx
    write(*  ,'(3X,A45,I16,'' nodes'')')                adjl('Nodes across the domain:',45), nx
else
    write(iow,'(3X,A45,I16)')'Error: domain discretization was not set'
    write(*  ,'(3X,A45,I16)')'Error: domain discretization was not set'
    STOP
endif

if (log_temperature) then
    write(iow,'(3X,A45,F16.4,'' K'')')adjl('Temperature:',45),Temp
    write(*  ,'(3X,A45,F16.4,'' K'')')adjl('Temperature:',45),Temp
    beta = 1.d0 / (boltz_const_Joule_K * Temp)
else
    write(iow,'(3X,A45)')'Error: temperature not found..'
    write(*  ,'(3X,A45)')'Error: temperature not found..'
    STOP
endif

if (log_pressure) then
    write(iow,'(3X,A45,F16.4,'' atm'')')adjl('Pressure:',45), Pressure
    write(*  ,'(3X,A45,F16.4,'' atm'')')adjl('Pressure:',45), Pressure
    Pressure = Pressure * atm_to_pa
else
    Pressure = 0.d0
    write(iow,'(3X,A45,F16.4,'' atm'')')adjl('*Pressure not found. Auto:',45),Pressure
    write(iow,'(3X,A45,F16.4,'' atm'')')adjl('*Pressure not found. Auto:',45),Pressure
endif

write(iow,'(A85)')adjl('----------------------------------POLYMER PROPERTIES---------------------------------',85)
write(*  ,'(A85)')adjl('----------------------------------POLYMER PROPERTIES---------------------------------',85)

if (log_characteristic_ratio) then
    write(iow,'(3X,A45,F16.4)')adjl('Chain characteristic ratio C_infinity:',45) ,CN
    write(*  ,'(3X,A45,F16.4)')adjl('Chain characteristic ratio C_infinity:',45) ,CN
else
    write(iow,'(3X,A45)')adjl('Error: chain characteristic ratio not found..',45)
    write(*  ,'(3X,A45)')adjl('Error: chain characteristic ratio not found..',45)
    STOP
endif

if (log_bond_length) then
    write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('Bond length:',45), bond_length
    write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('Bond length:',45), bond_length
else
    write(iow,'(3X,A45)') 'Error: bond length not found..'
    write(*  ,'(3X,A45)') 'Error: bond length not found..'
    STOP
endif

if (log_monomer_mass) then
    write(iow,'(3X,A45,F16.4,'' g/mol'')')adjl('Monomer mass:',45), mon_mass
    write(*  ,'(3X,A45,F16.4,'' g/mol'')')adjl('Monomer mass:',45), mon_mass
else
    write(iow,'(3X,A45)') 'Error: monomer mass not found'
    write(*  ,'(3X,A45)') 'Error: monomer mass not found'
    STOP
endif

if (log_mass_density) then
    write(iow,'(3X,A45,F16.4,'' g/cm^3'')')adjl('Mass density:',45), rho_mass_bulk
    write(*  ,'(3X,A45,F16.4,'' g/cm^3'')')adjl('Mass density:',45), rho_mass_bulk
else
    write(iow,'(3X,A45)') 'Error: Mass density not found'
    write(*  ,'(3X,A45)') 'Error: Mass density not found'
    STOP
endif
rho_mass_bulk = rho_mass_bulk * gr_cm3_to_kg_m3 ! SI

write(iow,'(A85)')adjl('--------------------------------SIMULATION PARAMETERS--------------------------------',85)
write(*  ,'(A85)')adjl('--------------------------------SIMULATION PARAMETERS--------------------------------',85)

if (log_field_every) then
    write(iow,'(3X,A45,I16,A16)')adjl('Export field every:',45),field_every,adjustl('steps')
    write(*  ,'(3X,A45,I16,A16)')adjl('Export field every:',45),field_every,adjustl('steps')
    if (field_every.le.0) then
        write(iow,'(3X,A45)') '*set a positive field output value'
        write(*  ,'(3X,A45)') '*set a positive field output value'
        STOP
    endif  
else
    field_every = 1000
    continue
endif

if (log_thermo_every) then
    write(iow,'(3X,A45,I16,A16)')adjl('Output thermodynamics every:',45),thermo_every,adjustl('steps')
    write(*  ,'(3X,A45,I16,A16)')adjl('Output thermodynamics every:',45),thermo_every,adjustl('steps')
    if (thermo_every.le.0) then
        write(iow,'(3X,A16)') '*set a positive thermo output value'
        write(*  ,'(3X,A16)') '*set a positive thermo output value'
        STOP
    endif  
else
    thermo_every = 1000
    continue
endif

if (log_compute_every) then
    write(iow,'(3X,A45,I16,'' steps'')')adjl('Output computes every:',45),compute_every
    write(*  ,'(3X,A45,I16,'' steps'')')adjl('Output computes every:',45),compute_every
    if (compute_every.le.0) then
        write(iow,'(3X,A45)') '*Computes have been skipped'
        write(6  ,'(3X,A45)') '*Computes have been skipped'
    endif  
else
    compute_every = 50000
    continue
endif

if (log_out_phi) then
    write(iow,'(3X,A45,L16)')adjl('Output phi:',45),out_phi
    write(6  ,'(3X,A45,L16)')adjl('Output phi:',45),out_phi
else
    out_phi = .True.
    continue
endif

if (log_out_field) then
    write(iow,'(3X,A45,L16)')adjl('Output field:',45),out_field
    write(6  ,'(3X,A45,L16)')adjl('Output field:',45),out_field
else
    out_field = .True.
    continue
endif

if (log_out_q) then
    write(iow,'(3X,A45,L16)')adjl('Output q:',45),out_q
    write(6  ,'(3X,A45,L16)')adjl('Output q:',45),out_q
else
    out_q = .True.
    continue
endif

if (log_out_chainshape) then
    write(iow,'(3X,A45,L16)')adjl('Output chain shape:',45),out_chainshape
    write(6  ,'(3X,A45,L16)')adjl('Output chain shape:',45),out_chainshape
else
    out_chainshape = .True.
    continue
endif

if (log_out_ads_free) then
    write(iow,'(3X,A45,L16)')adjl('Output ads free:',45),out_ads_free
    write(6  ,'(3X,A45,L16)')adjl('Output ads free:',45),out_ads_free
else
    out_ads_free = .True.
    continue
endif

if (log_out_end_middle) then
    write(iow,'(3X,A45,L16)')adjl('Output end middle:',45),out_end_middle
    write(6  ,'(3X,A45,L16)')adjl('Output end middle:',45),out_end_middle
else
    out_end_middle = .True.
    continue
endif

if (log_out_brush_thickness) then
    write(iow,'(3X,A45,L16)')adjl('Output brush thickness:',45),out_brush_thickness
    write(6  ,'(3X,A45,L16)')adjl('Output brush thickness:',45),out_brush_thickness
else
    out_brush_thickness = .True.
    continue
endif

if (log_check_stability_every) then
    write(iow,'(3X,A45,I16,A16)')adjl('The stability will be checked every:',45),check_stability_every,adjustl('steps')
    write(6  ,'(3X,A45,I16,A16)')adjl('The stability will be checked every:',45),check_stability_every,adjustl('steps')
    if (check_stability_every.le.0) then
        write(iow,'(3X,A45)') '*stability check will not be performed'
        write(6  ,'(3X,A45)') '*stability check will not be performed'
    endif  
else
    check_stability_every = 0
    continue
endif

if (log_read_field) then
    if (read_field)then
        write(iow,'(3X,A45,A16)')adjl('Field will be read from file:',45),adjustl(field_filename)
        write(*  ,'(3X,A45,A16)')adjl('Field will be read from file:',45),adjustl(field_filename)
    else
        write(iow,'(3X,A45)')adjl('Field will be initialized to zero',45)
        write(*  ,'(3X,A45)')adjl('Field will be initialized to zero',45)
    end if
else
    write(iow,'(3X,A45,'' initialization to zero'')')adjl('Field read flag not found. Auto:',45)
    write(*  ,'(3X,A45,'' initialization to zero'')')adjl('Field read flag not found. Auto:',45)
    read_field=.false.
endif

if (log_max_wa_error) then
    write(iow,'(3X,A45,F16.4,'' J/k_BT'')')adjl('Field error tolerance:',45),max_wa_error
    write(*  ,'(3X,A45,F16.4,'' J/k_BT'')')adjl('Field error tolerance:',45),max_wa_error
else
    max_wa_error = 0.0D+00
    write(iow,'(3X,A45,E16.9,'' k_B T'')')adjl('*Maximum field error not found. Auto:',45),max_wa_error
    write(*  ,'(3X,A45,E16.9,'' k_B T'')')adjl('*Maximum field error not found. Auto:',45),max_wa_error
endif

if (log_fraction_of_new_field) then
    write(iow,'(3X,A45,F16.9)')adjl('Field mixing fraction:',45),frac
    write(*  ,'(3X,A45,F16.9)')adjl('Field mixing fraction:',45),frac
else
    write(iow,'(3X,A45)')'*Field mixing fraction not found..'
    write(*  ,'(3X,A45)')'*Field mixing fraction not found..'
    STOP
endif

if (log_edwards_solver) then
    if (edwards_solver.eq.F_implicit) then
        write(iow,'(3X,A45,A16)')adjl('Edwards solver',45),adjustl('implicit')
        write(*  ,'(3X,A45,A16)')adjl('Edwards solver',45),adjustl('implicit')
    elseif (edwards_solver.eq.F_semi_implicit) then
        write(iow,'(3X,A45,A16)')adjl('Edwards solver',45),adjustl('semi-implicit')
        write(*  ,'(3X,A45,A16)')adjl('Edwards solver',45),adjustl('semi-implicit')
    else
        write(iow,'(3X,A61)')adjl('*Edwards solver flag must be set to either 0 or 1..',45)
        write(*  ,'(3X,A61)')adjl('*Edwards solver flag must be set to either 0 or 1..',45)
        STOP
    endif
else
    edwards_solver = F_implicit
    write(iow,'(3X,A45)')adjl('*Edwards solver not set.. we will use the implicit..',45)
    write(*  ,'(3X,A45)')adjl('*Edwards solver not set.. we will use the implicit..',45)
endif

if (log_number_of_iterations) then
    write(iow,'(3X,A45,I16,'' iterations'')')adjl('Maximum iterations:',45),max_iter
    write(*  ,'(3X,A45,I16,'' iterations'')')adjl('Maximum iterations:',45),max_iter
else
    max_iter = 500000   
    write(iow,'(3X,A45)')adjl('*Maximum iterations not found. Auto:',45),max_iter
    write(*  ,'(3X,A45)')adjl('*Maximum iterations not found. Auto:',45),max_iter
endif

if (log_contour_discret_scheme) then
    if (contour_discret_scheme.eq.F_uniform) then
        write(iow,'(3X,A45,A16)')adjl("Spacing for chain contour discret:",45),adjustl('uniform')
        write(*  ,'(3X,A45,A16)')adjl("Spacing for chain contour discret:",45),adjustl('uniform')
    elseif ( contour_discret_scheme.eq.F_nonuniform) then
        write(iow,'(3X,A45,A16)')adjl("Spacing for chain contour discret:",45),adjustl('nonuniform')
        write(*  ,'(3X,A45,A16)')adjl("Spacing for chain contour discret:",45),adjustl('nonuniform')
    else
        write(iow,'(3X,A61,I3,'' or '',I3)')adjl('*Error: Spacing for chain contour discret must be set to',61),F_uniform,F_nonuniform
        write(*  ,'(3X,A61,I3,'' or '',I3)')adjl('*Error: Spacing for chain contour discret must be set to',61),F_uniform,F_nonuniform
    endif
else
    contour_discret_scheme = F_uniform
    write(iow,'(3X,A45,A16)')adjl("*Chain contour spacing discret not found. Auto:",45),adjustl('uniform')
    write(*  ,'(3X,A45,A16)')adjl("*Chain contour spacing discret not found. Auto:",45),adjustl('uniform')
endif


if (log_spatial_discret_scheme) then
    if (spatial_discret_scheme.eq.F_uniform) then
        write(iow,'(3X,A45,A16)')adjl("Spacing for spatial discret:",45),adjustl('uniform')
        write(*  ,'(3X,A45,A16)')adjl("Spacing for spatial discret:",45),adjustl('uniform')
    elseif ( spatial_discret_scheme.eq.F_nonuniform) then
        write(iow,'(3X,A45,A16)')adjl("Spacing for spatial discret:",45),adjustl('nonuniform')
        write(*  ,'(3X,A45,A16)')adjl("Spacing for spatial discret:",45),adjustl('nonuniform')
    else
        write(iow,'(3X,A61,I3,'' or '',I3)')adjl('*Error: Spacing for spatial discret must be set to',61),F_uniform,F_nonuniform
        write(*  ,'(3X,A61,I3,'' or '',I3)')adjl('*Error: Spacing for spatial discret must be set to',61),F_uniform,F_nonuniform
    endif
else
    spatial_discret_scheme = F_uniform
    write(iow,'(3X,A45,A16)')adjl("*Spatial spacing discret not found. Auto:",45),adjustl('uniform')
    write(*  ,'(3X,A45,A16)')adjl("*Spatial spacing discret not found. Auto:",45),adjustl('uniform')
endif

if (log_contour_integr_scheme) then
    if (contour_integr_scheme.eq.F_rectangle_rule) then
        write(iow,'(3X,A45,A16)')adjl("Chain contour integration rule:",45),adjustl('Recrangle rule')
        write(*  ,'(3X,A45,A16)')adjl("Chain contour integration rule:",45),adjustl('Recrangle rule')
    elseif ( contour_integr_scheme.eq.F_simpson_rule) then
        write(iow,'(3X,A45,A16)')adjl("Chain contour integration rule:",45),adjustl('Simpson rule')
        write(*  ,'(3X,A45,A16)')adjl("Chain contour integration rule:",45),adjustl('Simpson rule')
    else
        write(iow,'(3X,A61,I3,'' or '',I3)')adjl('*Error: Chain contour integr rule must be set to',61),F_rectangle_rule, F_simpson_rule
        write(*  ,'(3X,A61,I3,'' or '',I3)')adjl('*Error: Chain contour integr rule must be set to',61),F_rectangle_rule, F_simpson_rule
    endif
else
    contour_integr_scheme = F_simpson_rule
    write(iow,'(3X,A45,A16)')adjl("*Chain contour integr rule not found. Auto:",45),adjustl('Simpson rule')
    write(*  ,'(3X,A45,A16)')adjl("*Chain contour integr rule not found. Auto:",45),adjustl('Simpson rule')
endif

if (log_spatial_integr_scheme) then
    if (spatial_integr_scheme.eq.F_rectangle_rule) then
        write(iow,'(3X,A45,A16)')adjl("Spatial integration rule:",45),adjustl('Recrangle rule')
        write(*  ,'(3X,A45,A16)')adjl("Spatial integration rule:",45),adjustl('Recrangle rule')
    elseif ( spatial_integr_scheme.eq.F_simpson_rule) then
        write(iow,'(3X,A45,A16)')adjl("Spatial integration rule:",45),adjustl('Simpson rule')
        write(*  ,'(3X,A45,A16)')adjl("Spatial integration rule:",45),adjustl('Simpson rule')
    else
        write(iow,'(3X,A61,I3,'' or '',I3)')adjl('*Error: Spatial integr rule must be set to',61),F_rectangle_rule, F_simpson_rule
        write(*  ,'(3X,A61,I3,'' or '',I3)')adjl('*Error: Spatial integr rule must be set to',61),F_rectangle_rule, F_simpson_rule
    endif
else
    spatial_integr_scheme = F_simpson_rule
    write(iow,'(3X,A45,A16)')adjl("*Spatial integr rule not found. Auto:",45),adjustl('Simpson rule')
    write(*  ,'(3X,A45,A16)')adjl("*Spatial integr rule not found. Auto:",45),adjustl('Simpson rule')
endif

if (log_lo_BC_of_matrix) then
    if (bc_lo_matrix == F_bc_neuman)then
        write(iow,'(3X,A45,A16,'' (dq/dr=0)'')')adjl('matrix boundary condition for lo edge:',45),adjustl('Neumann')
        write(*  ,'(3X,A45,A16,'' (dq/dr=0)'')')adjl('matrix boundary condition for lo edge:',45),adjustl('Neumann')
    else if(bc_lo_matrix == F_bc_dirichlet_eq_0) then
        write(iow,'(3X,A45,A16,'' q=0'')')adjl('matrix boundary condition for lo edge:',45),adjustl('Dirichlet')
        write(*  ,'(3X,A45,A16,'' q=0'')')adjl('matrix boundary condition for lo edge:',45),adjustl('Dirichlet')
    else if(bc_lo_matrix == F_bc_dirichlet_eq_1) then
        write(iow,'(3X,A45,A16,'' q=1'')')adjl('matrix boundary condition for lo edge:',45),adjustl('Dirichlet')
        write(*  ,'(3X,A45,A16,'' q=1'')')adjl('matrix boundary condition for lo edge:',45),adjustl('Dirichlet')
    else
        write(iow,'(3X,A150)')adjl('Error: wrong matrix boundary condition for lo edge.. (choose between -1, 0 and 1',150)
        write(*  ,'(3X,A150)')adjl('Error: wrong matrix boundary condition for lo edge.. (choose between -1, 0 and 1',150)
        STOP
    end if 
else
    write(iow,'(3X,A150)')adjl('Error: matrix boundary condition for lo edge not found..',150)
    write(*  ,'(3X,A150)')adjl('Error: matrix boundary condition for lo edge not found..',150)
    STOP
endif

if (log_hi_BC_of_matrix) then
    if (bc_hi_matrix == F_bc_neuman)then
        write(iow,'(3X,A45,A16,'' (dq/dr=0)'')')adjl('matrix boundary condition for hi edge:',45),adjustl('Neumann')
        write(*  ,'(3X,A45,A16,'' (dq/dr=0)'')')adjl('matrix boundary condition for hi edge:',45),adjustl('Neumann')
    else if(bc_hi_matrix == F_bc_dirichlet_eq_0) then
        write(iow,'(3X,A45,A16,'' q=0'')')adjl('matrix boundary condition for hi edge:',45),adjustl('Dirichlet')
        write(*  ,'(3X,A45,A16,'' q=0'')')adjl('matrix boundary condition for hi edge:',45),adjustl('Dirichlet')
    else if(bc_hi_matrix == F_bc_dirichlet_eq_1) then
        write(iow,'(3X,A45,A16,'' q=1'')')adjl('matrix boundary condition for hi edge:',45),adjustl('Dirichlet')
        write(*  ,'(3X,A45,A16,'' q=1'')')adjl('matrix boundary condition for hi edge:',45),adjustl('Dirichlet')
    else
        write(iow,'(3X,A150)')adjl('Error: wrong matrix boundary condition for hi edge.. (choose between -1, 0 and 1',150)
        write(*  ,'(3X,A150)')adjl('Error: wrong matrix boundary condition for hi edge.. (choose between -1, 0 and 1',150)
        STOP
    end if 
else
    write(iow,'(3X,A150)')adjl('Error: matrix boundary condition for hi edge not found..',150)
    write(*  ,'(3X,A150)')adjl('Error: matrix boundary condition for hi edge not found..',150)
    STOP
endif

if (log_lo_BC_of_grafted) then
    if (bc_lo_grafted == F_bc_neuman)then
        write(iow,'(3X,A45,A16,'' (dq/dr=0)'')')adjl('grafted boundary condition for lo edge:',45),adjustl('Neumann')
        write(*  ,'(3X,A45,A16,'' (dq/dr=0)'')')adjl('grafted boundary condition for lo edge:',45),adjustl('Neumann')
    else if(bc_lo_grafted == F_bc_dirichlet_eq_0) then
        write(iow,'(3X,A45,A16,'' q=0'')')adjl('grafted boundary condition for lo edge:',45),adjustl('Dirichlet')
        write(*  ,'(3X,A45,A16,'' q=0'')')adjl('grafted boundary condition for lo edge:',45),adjustl('Dirichlet')
    else if(bc_lo_grafted == F_bc_dirichlet_eq_1) then
        write(iow,'(3X,A45,A16,'' q=1'')')adjl('grafted boundary condition for lo edge:',45),adjustl('Dirichlet')
        write(*  ,'(3X,A45,A16,'' q=1'')')adjl('grafted boundary condition for lo edge:',45),adjustl('Dirichlet')
    else
        write(iow,'(3X,A150)')adjl('Error: wrong grafted boundary condition for lo edge.. (choose between -1, 0 and 1',150)
        write(*  ,'(3X,A150)')adjl('Error: wrong grafted boundary condition for lo edge.. (choose between -1, 0 and 1',150)
        STOP
    end if 
else
    write(iow,'(3X,A150)')adjl('Error: grafted boundary condition for lo edge not found..',150)
    write(*  ,'(3X,A150)')adjl('Error: grafted boundary condition for lo edge not found..',150)
    STOP
endif

if (log_hi_BC_of_grafted) then
    if (bc_hi_grafted == F_bc_neuman)then
        write(iow,'(3X,A45,A16,'' (dq/dr=0)'')')adjl('grafted boundary condition for hi edge:',45),adjustl('Neumann')
        write(*  ,'(3X,A45,A16,'' (dq/dr=0)'')')adjl('grafted boundary condition for hi edge:',45),adjustl('Neumann')
    else if(bc_hi_grafted == F_bc_dirichlet_eq_0) then
        write(iow,'(3X,A45,A16,'' q=0'')')adjl('grafted boundary condition for hi edge:',45),adjustl('Dirichlet')
        write(*  ,'(3X,A45,A16,'' q=0'')')adjl('grafted boundary condition for hi edge:',45),adjustl('Dirichlet')
    else if(bc_hi_grafted == F_bc_dirichlet_eq_1) then
        write(iow,'(3X,A45,A16,'' q=1'')')adjl('grafted boundary condition for hi edge:',45),adjustl('Dirichlet')
        write(*  ,'(3X,A45,A16,'' q=1'')')adjl('grafted boundary condition for hi edge:',45),adjustl('Dirichlet')
    else
        write(iow,'(3X,A150)')adjl('Error: wrong grafted boundary condition for hi edge.. (choose between -1, 0 and 1',150)
        write(*  ,'(3X,A150)')adjl('Error: wrong grafted boundary condition for hi edge.. (choose between -1, 0 and 1',150)
        STOP
    end if 
else
    write(iow,'(3X,A150)')adjl('Error: grafted boundary condition for hi edge not found..',150)
    write(*  ,'(3X,A150)')adjl('Error: grafted boundary condition for hi edge not found..',150)
    STOP
endif

if (log_matrix_exist) then
    if (matrix_exist) then
        write(iow,'(A85)')adjl('-------------------------------------MATRIX CHAINS-----------------------------------',85)
        write(*  ,'(A85)')adjl('-------------------------------------MATRIX CHAINS-----------------------------------',85)
    else
        matrix_exist = .false.
    endif
else
    continue
endif

if (matrix_exist) then ! Start of matrix chains
if (log_ds_ave_matrix) then
    write(iow,'(3X,A45,E16.9,'' monomers'')')adjl('Chain contour discretization:',45), ds_ave_matrix
    write(*  ,'(3X,A45,E16.9,'' monomers'')')adjl('Chain contour discretization:',45), ds_ave_matrix
else
    write(iow,'(3X,A45)')'*Chain contour discretization was not set'
    write(iow,'(3X,A45)')'*Chain contour discretization was not set'
    STOP
endif

if (log_chain_length_matrix) then
    write(iow,'(3X,A45,F16.9,'' monomers'')')adjl('Chain length:',45), chainlen_matrix
    write(*  ,'(3X,A45,F16.9,'' monomers'')')adjl('Chain length:',45), chainlen_matrix
else
    write(iow,'(3X,A45)') 'Chain length of matrix chains was not set'
    write(*  ,'(3X,A45)') 'Chain length of matrix chains was not set'
    STOP
endif

if (log_r_critical) then
    if (r_critical > 0.d0)then
        write(iow,'(3X,A45,F16.9, '' Angstrom'')')adjl('Critical position of adsorbed segs',45), r_critical
        write(*  ,'(3X,A45,F16.9, '' Angstrom'')')adjl('Critical position of adsorbed segs',45), r_critical
    else
        write(iow,*)"r_critical < 0"
        write(*  ,*)"r_critical < 0"
        STOP
    end if
else
    r_critical = 10.0d0
    write(iow,'(3X,A45,E16.9, '' Angstrom'')')adjl('*Grafting position not found. It was set to',45),r_critical,adjustl('Angstrom')
    write(*  ,'(3X,A45,E16.9, '' Angstrom'')')adjl('*Grafting position not found. It was set to',45),r_critical,adjustl('Angstrom')
endif
else
    ds_ave_matrix = 0.d0
    chainlen_matrix = 0.d0
endif ! end of matrix chains


if (log_grafted_lo_exist) then
    if (grafted_lo_exist) then
        write(iow,'(A85)')adjl('----------------------------------GRAFTED LO CHAINS----------------------------------',85)
        write(*  ,'(A85)')adjl('----------------------------------GRAFTED LO CHAINS----------------------------------',85)
    else
        grafted_lo_exist = .false.
    endif
else
    grafted_lo_exist = .false.
endif

if (grafted_lo_exist) then ! start of grafted lo chains
if (log_ds_ave_grafted_lo) then
    write(iow,'(3X,A45,E16.9,'' monomers'')')adjl('Chain contour discretization:',45), ds_ave_grafted_lo
    write(*  ,'(3X,A45,E16.9,'' monomers'')')adjl('Chain contour discretization:',45), ds_ave_grafted_lo
else
    write(iow,'(3X,A45)')'*Chain contour discretization was not set'
    write(iow,'(3X,A45)')'*Chain contour discretization was not set'
    STOP
endif

if (log_chain_length_grafted_lo) then
    write(iow,'(3X,A45,F16.9,'' monomers'')')adjl('Chain length:',45), chainlen_grafted_lo
    write(*  ,'(3X,A45,F16.9,'' monomers'')')adjl('Chain length:',45), chainlen_grafted_lo
else
    write(iow,'(3X,A45)') 'Chain length of grafted_lo chains was not set'
    write(*  ,'(3X,A45)') 'Chain length of grafted_lo chains was not set'
    STOP
endif

if (log_gdens_lo) then
    write(iow,'(3X,A45,E16.9,'' chains/Angstrom^2'')')adjl('Grafting density:',45), gdens_lo
    write(*  ,'(3X,A45,E16.9,'' chains/Angstrom^2'')')adjl('Grafting density:',45), gdens_lo
    if (gdens_lo < 0.d0)then
        write(iow,*)"ERROR: gdens < 0"
        write(*  ,*)"ERROR: gdens < 0"
        STOP
    endif
else
    gdens_lo = 0.d0
    write(iow,'(3X,A45)') adjl('Grafting density not found.. It was set to zero.',45)
    write(*  ,'(3X,A45)') adjl('Grafting density not found.. It was set to zero.',45)
endif
else
    ds_ave_grafted_lo = 0.d0
    chainlen_grafted_lo = 0.d0
    gdens_lo = 0.d0
endif ! end of grafted lo chains

if (log_grafted_hi_exist) then
    if (grafted_hi_exist) then
        write(iow,'(A85)')adjl('----------------------------------GRAFTED HI CHAINS----------------------------------',85)
        write(*  ,'(A85)')adjl('----------------------------------GRAFTED HI CHAINS----------------------------------',85)
    else
        grafted_hi_exist = .false.
    endif
else
    grafted_hi_exist = .false.
endif

if (grafted_hi_exist) then ! start of grafted hi chains
if (log_ds_ave_grafted_hi) then
    write(iow,'(3X,A45,E16.9,'' monomers'')')adjl('Chain contour discretization:',45), ds_ave_grafted_hi
    write(*  ,'(3X,A45,E16.9,'' monomers'')')adjl('Chain contour discretization:',45), ds_ave_grafted_hi
else
    write(iow,'(3X,A45)')'*Chain contour discretization was not set'
    write(iow,'(3X,A45)')'*Chain contour discretization was not set'
    STOP
endif

if (log_chain_length_grafted_hi) then
    write(iow,'(3X,A45,F16.9,'' monomers'')')adjl('Chain length:',45), chainlen_grafted_hi
    write(*  ,'(3X,A45,F16.9,'' monomers'')')adjl('Chain length:',45), chainlen_grafted_hi
else
    write(iow,'(3X,A45)') 'Chain length of grafted_hi chains was not set'
    write(*  ,'(3X,A45)') 'Chain length of grafted_hi chains was not set'
    STOP
endif

if (log_gdens_hi) then
    write(iow,'(3X,A45,E16.9,'' chains/Angstrom^2'')')adjl('Grafting density:',45), gdens_hi
    write(*  ,'(3X,A45,E16.9,'' chains/Angstrom^2'')')adjl('Grafting density:',45), gdens_hi
    if (gdens_hi < 0.d0)then
        write(iow,*)"ERROR: gdens < 0"
        write(*  ,*)"ERROR: gdens < 0"
        STOP
    endif
else
    gdens_hi = 0.d0
    write(iow,'(3X,A45)') adjl('Grafting density not found.. It was set to zero.',45)
    write(*  ,'(3X,A45)') adjl('Grafting density not found.. It was set to zero.',45)
endif
else
    ds_ave_grafted_hi = 0.d0
    chainlen_grafted_hi = 0.d0
    gdens_hi = 0.d0
endif ! end of grafted hi chains

if (grafted_lo_exist.or.grafted_hi_exist) then
if (log_position_of_grafted) then
    if (abs(graft_pos)<tol)then
        gnode_lo = 1
        write(iow,'(3X,A45,I16,'' nodes'')') adjl('Solid-grafted point distance:',45), gnode_lo
        write(*  ,'(3X,A45,I16,'' nodes'')') adjl('Solid-grafted point distance:',45), gnode_lo
    elseif (graft_pos > 0.d0) then
        gnode_lo = int(anint((graft_pos - wall_pos) * dble(nx) / lx))
        write(*,*)graft_pos, wall_pos, nx, lx, gnode_lo
        if (gnode_lo .le. 1) then
            gnode_lo = 1
        endif
        if (graft_pos.lt.wall_pos) then
            write(iow,*)"*Error: graft_pos < wall_pos"
            write(*  ,*)"*Error: graft_pos < wall_pos"
            STOP
        endif
        write(iow,'(3X,A45,I16,'' nodes'')') adjl('Solid-grafted point distance:',45), gnode_lo
        write(*  ,'(3X,A45,I16,'' nodes'')') adjl('Solid-grafted point distance:',45), gnode_lo
    else
        write(iow,*)"graft_pos < 0"
        write(*  ,*)"graft_pos < 0"
        STOP
    end if
else
    gnode_lo = 1
    write(iow,'(3X,A45,I6,'' nodes'')')adjl('*Solid-grafted distance not found. Auto:',45), gnode_lo
    write(*  ,'(3X,A45,I6,'' nodes'')')adjl('*Solid-grafted distance not found. Auto:',45), gnode_lo
endif
gnode_hi = nx - gnode_lo
endif

write(iow,'(A85)')adjl('------------------------------SETUP THE EQUATION OF STATE----------------------------',85)
write(*  ,'(A85)')adjl('------------------------------SETUP THE EQUATION OF STATE----------------------------',85)
if (log_eos_type) then
    if (eos_type.eq.F_helfand) then
        write(iow,'(3X,A45)')adjl('The Helfand EoS was chosen with coeffs:',45)
        write(*  ,'(3X,A45)')adjl('The Helfand EoS was chosen with coeffs:',45)
    elseif (eos_type.eq.F_sanchez_labombe) then
        write(iow,'(3X,A45)')adjl('The Sanchez-Lacombe EoS was chosen with coeffs:',45)
        write(*  ,'(3X,A45)')adjl('The Sanchez-Lacombe EoS was chosen with coeffs:',45)
    else
        write(iow,'(A45,I16)') 'Error: EOS flag different than 0 (HF) or 1 (SL)',eos_type
        write(*  ,'(A45,I16)') 'Error: EOS flag different than 0 (HF) or 1 (SL)',eos_type
        STOP
    endif
else
    write(iow,'(3X,A45)')'Error: EOS flag not set'
    write(*  ,'(3X,A45)')'Error: EOS flag not set'
    STOP
endif

if (log_eos_coeffs) then
    if (eos_type.eq.F_helfand) then
        write(iow,'(3X,A45,1(E16.4),'' Pa^-1'')') adjl('*Isothermal compressibility:',45),HF_kappa_T
        write(*  ,'(3X,A45,1(E16.4),'' Pa^-1'')') adjl('*Isothermal compressibility:',45),HF_kappa_T
    elseif (eos_type.eq.F_sanchez_labombe) then
        write(iow,'(3X,A45,F16.4,'' kg/m^3'')')   adjl('*rho_star :',45),rho_star
        write(iow,'(3X,A45,F16.4,'' K'')')        adjl('*T_star   :',45),T_star
        write(iow,'(3X,A45,F16.4,'' Pa'')')       adjl('*P_star   :',45),P_star
        write(*  ,'(3X,A45,F16.4,'' kg/m^3'')')   adjl('*rho_star :',45),rho_star
        write(*  ,'(3X,A45,F16.4,'' K'')')        adjl('*T_star   :',45),T_star
        write(*  ,'(3X,A45,F16.4,'' Pa'')')       adjl('*P_star   :',45),P_star
    endif
else
    write(iow,'(3X,A45)') adjl('Error: EOS coeffs not set',45)
    write(*  ,'(3X,A45)') adjl('Error: EOS coeffs not set',45)
    STOP
endif

if (log_influence_param) then
    write(iow,'(3X,A45,F16.4)') adjl('Reduced influence parameter:',45),k_gr_tilde
    write(*  ,'(3X,A45,F16.4)') adjl('Reduced influence parameter:',45),k_gr_tilde
else
    k_gr_tilde = 0.d0
endif
if (log_real_influence_param) then
    write(iow,'(3X,A45,E16.9,'' J m^5/mol^2'')') adjl('Influence parameter:',45),k_gr
    write(*  ,'(3X,A45,E16.9,'' J m^5/mol^2'')') adjl('Influence parameter:',45),k_gr
else
    k_gr = 0.d0
endif

if (log_influence_param.and.log_real_influence_param) then
    write(iow,'(3X,A45)') adjl('Error: Set either the real OR the reduced infl. parameter',45)
    write(*  ,'(3X,A45)') adjl('Error: Set either the real OR the reduced infl. parameter',45)
    STOP
endif

if (log_influence_param.or.log_real_influence_param) then
    square_gradient = .true.
else
    square_gradient = .false.
endif


write(iow,'(A85)')adjl('------------------------------------SETUP THE WALLS----------------------------------',85)
write(*  ,'(A85)')adjl('------------------------------------SETUP THE WALLS----------------------------------',85)

if (log_wall_type) then
    if (wall_hybrid) then
        write(iow,'(3X,A45)')adjl('Hybrid wall potential was selected:',45)
        write(*  ,'(3X,A45)')adjl('Hybrid wall potential was selected:',45)
    endif
    if (wall_hamaker.and.log_wall_coeffs_hamaker) then
        write(iow,'(3X,A45)')adjl('Coefficients of the Hamaker potential:',45)
        write(*  ,'(3X,A45)')adjl('Coefficients of the Hamaker potential:',45)
        write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('*sigma polymer:',45),sig_pol
        write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('*sigma solid:',45),sig_solid
        write(iow,'(3X,A45,F16.4,'' 10^-20 J'')')adjl('*Hamaker constant of polymer:',45), Apol
        write(iow,'(3X,A45,F16.4,'' 10^-20 J'')')adjl('*Hamaker constant of solid:',45), Asolid
        write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('*sigma polymer:',45),sig_pol
        write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('*sigma solid:',45),sig_solid
        write(*  ,'(3X,A45,F16.4,'' 10^-20 J'')')adjl('*Hamaker constant of polymer:',45), Apol
        write(*  ,'(3X,A45,F16.4,'' 10^-20 J'')')adjl('*Hamaker constant of solid:',45), Asolid
        Asolid = Asolid * 1.e-20 ! SI
        Apol   = Apol   * 1.e-20 ! SI
    if (wall_hamaker.and..not.log_wall_coeffs_hamaker) then
            write(iow,'(3X,A150)')adjl('Error: The coefficients of the hamaker potential were not found!',150)
            write(*  ,'(3X,A150)')adjl('Error: The coefficients of the hamaker potential were not found!',150)
            STOP
        endif
    endif
    if (wall_square_well.and.log_wall_coeffs_square_well) then
        write(iow,'(3X,A45)')adjl('Coefficients of the square well potential:',45)
        write(*  ,'(3X,A45)')adjl('Coefficients of the square well potential:',45)
        write(iow,'(3X,A45,F16.4,'' 10^-20 J'')')adjl('*Energy barrier:',45), A_sq_well
        write(iow,'(3X,A45,F16.4,'' kBT'')')     adjl('*Energy barrier:',45), A_sq_well* 1.e-20 * beta
        write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('*sigma:',45)         , sigma_sq_well
        write(*  ,'(3X,A45,F16.4,'' 10^-20 J'')')adjl('*Energy barrier:',45), A_sq_well
        write(*  ,'(3X,A45,F16.4,'' kBT'')')     adjl('*Energy barrier:',45), A_sq_well* 1.e-20 * beta
        write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('*sigma:',45)         , sigma_sq_well
        A_sq_well = A_sq_well * 1.e-20 ! SI
    endif
    if (wall_square_well.and..not.log_wall_coeffs_square_well) then
        write(iow,'(3X,A150)')adjl('Error: The coefficients of the square well potential were not found!',150)
        write(*  ,'(3X,A150)')adjl('Error: The coefficients of the square well potential were not found!',150)
        STOP
    endif
    if (wall_ramp.and.log_wall_coeffs_ramp) then
        write(iow,'(3X,A45)')adjl('Coefficients of the ramp potential:',45)
        write(*  ,'(3X,A45)')adjl('Coefficients of the ramp potential:',45)
        write(iow,'(3X,A45,F16.4,'' 10^-20 J'')')adjl('*Energy barrier:',45), A_ramp
        write(iow,'(3X,A45,F16.4,'' kBT'')')     adjl('*Energy barrier:',45), A_ramp* 1.e-20 * beta
        write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('*sigma:',45),sigma_ramp
        write(*  ,'(3X,A45,F16.4,'' 10^-20 J'')')adjl('*Energy barrier:',45), A_ramp
        write(*  ,'(3X,A45,F16.4,'' kBT'')')     adjl('*Energy barrier:',45), A_ramp* 1.e-20 * beta
        write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('*sigma:',45),sigma_ramp
        A_ramp = A_ramp * 1.e-20 ! SI
    endif
    if (wall_ramp.and..not.log_wall_coeffs_ramp) then
        write(iow,'(3X,A150)')adjl('Error: The coefficients of the ramp potential were not found!',150)
        write(*  ,'(3X,A150)')adjl('Error: The coefficients of the ramp potential were not found!',150)
        STOP
    endif

    if (wall_vacuum) then
        write(iow,'(3X,A45,A16)')adjl('Wall type:',45),adjustl('vacuum')
        write(*  ,'(3X,A45,A16)')adjl('Wall type:',45),adjustl('vacuum')
    endif
else
    wall_type = F_vacuum
    wall_vacuum = .true.
    write(iow,'(3X,A45)')adjl('*wall type not found.. it was set to vacuum',45)
    write(*  ,'(3X,A45)')adjl('*wall type not found.. it was set to vacuum',45)
endif

if (log_wall_side) then
    if (wall_side.eq.F_lo) then
        write(iow,'(3X,A45,A16)')adjl('Side of the solid wall:',45),'lo'
        write(*  ,'(3X,A45,A16)')adjl('Side of the solid wall:',45),'lo'
    elseif (wall_side.eq.F_both) then
        write(iow,'(3X,A45,A16)')adjl('Side of the solid wall:',45),'lo and hi'
        write(*  ,'(3X,A45,A16)')adjl('Side of the solid wall:',45),'lo and hi'
    elseif (wall_side.eq.F_hi) then
        write(iow,'(3X,A45,A16)')adjl('Side of the solid wall:',45),'hi'
        write(*  ,'(3X,A45,A16)')adjl('Side of the solid wall:',45),'hi'
    else
        write(iow,'(3X,A61)')adjl('Error: Wrong value to wall side! (-1: lo, 0: both, 1:hi)',45)
        write(*  ,'(3X,A61)')adjl('Error: Wrong value to wall side! (-1: lo, 0: both, 1:hi)',45)
        STOP
    endif
else
    write(iow,'(3X,A45)')adjl('*Side of the solid wall not detected..',45)
    write(*  ,'(3X,A45)')adjl('*Side of the solid wall not detected..',45)
    if (wall_vacuum) then
        wall_side = F_both
        write(iow,'(3X,A45)')'    ..it will be set to both sides'
        write(*  ,'(3X,A45)')'    ..it will be set to both sides'
    else
        STOP
    endif
endif

if (log_wall_pos) then
    write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('Position of the hard-sphere wall:',45), wall_pos
    write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('Position of the hard-sphere wall:',45), wall_pos
else
    wall_pos = 0.d0
    write(iow,'(3X,A45,F16.4,'' Angstrom'')')adjl('*Hard sphere wall position not found. Auto:',45),wall_pos
    write(*  ,'(3X,A45,F16.4,'' Angstrom'')')adjl('*Hard sphere wall position not found. Auto:',45),wall_pos
endif

if (log_wall_pos_auto.and.(wall_hamaker.or.wall_ramp)) then
    wall_auto = .true.
    write(iow,'(3X,A45,F16.4,'' k_B T'')')adjl('Recalibration of hard-sphere wall. E-target:',45), E_wall_target
    write(*  ,'(3X,A45,F16.4,'' k_B T'')')adjl('Recalibration of hard-sphere wall. E-target:',45), E_wall_target
else
    wall_auto = .false.
endif

return
!----------------------------------------------------------------------------------------------------------!
end subroutine parser
