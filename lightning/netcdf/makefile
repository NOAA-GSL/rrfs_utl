# makefile for lightning_porcess
#! /bin/sh -v

# Define the fortran compiler and options

include ../../configure.process

INCLD=  -I${NETCDF}/include -I../../../include
LIBS =  -L${NETCDF}/lib -lnetcdf  -L${GSILIBS} \
        ../../../src/main/constants.o

OBJS = Check_Alaska.o Check_NLDN.o \
       netCDFsub_geo.o netCDFsub_lightning.o \
       cio.o \
       module_debug.o misc_definitions_module.o constants_module.o \
       module_map_utils.o gridmod_gsi_map_utils.o write_bufr_lght.o \
       process_Lightning.o

all: process_Lightning

process_Lightning: ${OBJS}
	${FC} -o process_Lightning.exe ${FLAGS} ${OBJS} ${LIBS} 

.SUFFIXES : .F90 .f90 .f .c .o

.F90.o :
	${FC} ${FLAGS} ${INCLD} -c $<

.f90.o :
	${FC} ${FLAGS} ${INCLD} -c $<

.c.o :
	${CC} ${CPPFLAGS} -c $<

clean:
	/bin/rm -f *.o *.mod process_Lightning.exe
