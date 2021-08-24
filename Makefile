.SUFFIXES:	(.SUFFIXES) .F .h .p .f90


#new change from cjrevelas

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  ##
# 1. PREPROCESSOR                                                             ##
# 1.1 System Specific C PreProcessor                                          ##
# Linux:                                                                      ## 
CPP = /usr/bin/cpp -P -traditional -t -W 
CPPFLAGS = -I./include/ 
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ##

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##
# 2. SHARED SOURCE                                                            ##
# Shared source files directories:
SHARDIR = ./include  
# gnu Make folders to search for specific files                               ##
vpath %.F $(SHARDIR)
vpath %.f90 $(SHARDIR)
vpath %.c $(SHARDIR)
vpath %.h $(SHARDIR)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ##



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##
# 3. COMPILERS AND FLAGS                                                      ##
# 3.1 Intel Compiler                                                				##
ifeq  ($(notdir $(shell which ifort 2>&1)),ifort)
# 3.2 Intel Fortran serial compilation                                        ##
FCNEW = ifort
FCOLD = ifort -fpp 
CCOMP = icc -std=gnu99 
CPPCOMP = icpc
FRFLAGS = -O2 -prec-sqrt -prec-div -align -static -ip -ipo -heap-arrays
FDFLAGS = -O0 -g -traceback -fpe0 -fp-stack-check -heap-arrays -ftrapuv\
          -check pointers -check bounds -warn all
CDFLAGS = -O0 -g -traceback 

#
else
# 3.3 gnuFortran serial compilation                                           ## 
FCNEW = gfortran
FCOLD = gfortran
CCOMP = gcc -std=c99
CPPCOMP = g++
#FRFLAGS = -O3  
FRFLAGS = -O3 -ffree-line-length-none -m64 -mtune=native -march=native\
#          -funroll-loops -finline-limit=600 -fwhole-program -fstack-arrays 
FDFLAGS = -m64 -g -O0 -pedantic-errors -frepack-arrays -fdump-core -fbounds-check\
           -fimplicit-none -fbacktrace -ffree-line-length-none -frange-check\
           -Wall -Waliasing -Wampersand\
           -Wsurprising -Wunderflow -W
#Wimplicit-interface\ -Wconversion\

CDFLAGS = -g -fbacktrace

endif
# ------------------------------------------------                            ##                      



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##
# 4. Choose whether we are building for DEBUG or RELEASE                      ##
# 4.1 Debug or release flags for Fortran files                                ##
FCFLAGS = $(FRFLAGS) 
FCFLAGS = $(FDFLAGS)
# 4.2 Debug or release flags for C files                                      ##
CCFLAGS = $(FRFLAGS)
#CCFLAGS = $(CDFLAGS)
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ##



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##
# 5. Linker and appropriate flags
LD = $(FCNEW)
LDFLAGS = $(FCFLAGS)
LIBS = #-lstdc++ # -lm 
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ##



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##
# 6. Rules for compiling file-types                                           ##
.f90.o:
	$(FCNEW) -c $(FCFLAGS) $(LIBS) $*.f90
#.cpp.o: 
#	$(CPPCOMP) -c $(CCFLAGS) $(LIBS) $*.cpp
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ##


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##
MODULES =   flags.o\
	    constants.o\
            parser_vars.o\
            write_helper.o\
            eos.o\
            force_fields.o\
            integr_coeffs.o\
	    arrays.o\

OBJECTS =   parser.o\
            init_scf_params.o\
            init_arrays.o\
            init_mesh.o\
            init_geom.o\
            init_field.o\
            init_solid.o\
            generate_mesh.o\
            solver_edwards.o\
            solver_tridag.o\
            contour_convolution.o\
	    compute_energies.o\
	    get_part_func.o\
            get_nchains.o\
            get_auto_wall_pos.o\
            compute_chainshape.o\
	    compute_phi_seg.o\
	    compute_phi_ads_states.o\
	    compute_brush_thickness.o\
	    export_q.o\
            export_phi.o\
            export_field.o\
            export_field_binary.o\
            export_computes.o\
           
MAIN_OBJECT = main.o
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ##

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##
# 8. Executable                                                               ##
CMD = RuSseL1D

$(CMD): $(MODULES) $(OBJECTS) $(MAIN_OBJECT) 
	$(LD) -o $(CMD) $(LDFLAGS) $(LIBS) $(MODULES) $(OBJECTS) $(MAIN_OBJECT) 

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ##


clean:
	rm -f *.o *.mod 

cleaner:
	rm -f *.o *.mod *.exe fort.* o.* *__genmod.f90

cleantest:
	rm -f *.o *.mod *.exe fort.* o.* *__genmod.f90 LOG.*

cleanout:
	rm -f fort.* o.* *__genmod.f90 LOG.*

test:
	./test_integrity/test_integrity.sh
