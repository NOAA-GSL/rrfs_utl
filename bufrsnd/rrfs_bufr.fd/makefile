################################################################################
# 
#     Makefile for WRFBUFR code
#
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#     The following macros will be of interest:
#
#         CMD      - name of the executable
#         FC       - name of Fortran compiler
#         CPP      - name of CPP
#         ARCH     - architecture
#         CPPFLAGS - CPP flags
#         OPTS     - compiler code optimizations
#         TRAPS    - runtime traps for floating point exceptions
#         DEBUG    - -g
#
#################################################################################
#
# Define the name of the executable
#
CMD = rrfs_bufr.x
#
# CPP, Compiler, and Linker Options
#

FC       = ${COMP_MP}

CPP      = /lib/cpp -P
ARCH     =
CPPFLAGS = 
OPTS     = -O -fp-model precise  -convert big_endian -assume noold_ldout_format
FREE     = -FR
TRAPS    = 
DEBUG    = -g -traceback  -fpe0 -ftrapuv  -check all 


#
# Assemble Options
#
FFLAGS   = $(OPTS) $(TRAPS)  -I$(NETCDF_INCLUDES)
FFLAGST  = $(OPTS) $(FREE) $(TRAPS)   -I$(NETCDF_INCLUDES)

FFLAGS_DEBUG   = $(OPTS) $(TRAPS) $(DEBUG) -I$(NETCDF_INCLUDES)
FFLAGST_DEBUG  = $(OPTS) $(FREE) $(TRAPS) $(DEBUG) -I$(NETCDF_INCLUDES)

LDFLAGS  =
LIBS     =  \
           $(W3EMC_LIB4) $(W3NCO_LIB4) $(BACIO_LIB4) \
           -L$(NETCDF_LIBRARIES) -lnetcdf -lnetcdff
#
#
# free-format object files
#
OBJST=  PARA_RANGE.o

OBJST_DEBUG= PARA_RANGE.o
#
#
# fixed-format object files
#
#
OBJS= VRBLS2D_mod.o VRBLS3D_mod.o MASKS_mod.o SOIL_mod.o \
WRFBUFR.o  \
PROF_FV3SAR_NET.o iw3jdn.o w3fs26.o w3movdat.o w3reddat.o 

OBJS_DEBUG= VRBLS2D_mod.o VRBLS3D_mod.o MASKS_mod.o SOIL_mod.o \
WRFBUFR.o  \
PROF_FV3SAR_NET.o iw3jdn.o w3fs26.o w3movdat.o w3reddat.o

#
# Includes
#
INCLUDES= parmsoil
#

DEPS= $(INCLUDES)

.SUFFIXES:	.F .f .o

.F.f:
	$(CPP) $(CPPFLAGS) $< > $*.f

$(CMD):	$(OBJST) $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJST) $(OBJS) $(LIBS) 

debug:	$(OBJST_DEBUG) $(OBJS_DEBUG)
	$(FC) $(LDFLAGS) -o rrfs_bufr.x $(OBJST_DEBUG) $(OBJS_DEBUG) $(LIBS) 

$(OBJS)::	$(DEPS)
	$(FC) -c $(FFLAGS) $(@:.o=.f)

$(OBJS_DEBUG)::	$(DEPS)
	$(FC) -c $(FFLAGS_DEBUG) $(@:.o=.f)

$(OBJST)::	$(DEPS)
	$(FC) -c $(FFLAGST) $(@:.o=.f90)

$(OBJST_DEBUG)::	$(DEPS)
	$(FC) -c $(FFLAGST_DEBUG) $(@:.o=.f90)

clean:	
	/bin/rm -f  *.lst *.o *.mod
