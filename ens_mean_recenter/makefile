# makefile for nssl_porcess
#! /bin/sh -v

# Define the fortran compiler and options

FLAGS= -O3 -r8 -warn all -implicitnone -g -traceback -fp-model strict -qopenmp

INCLD= -I${NETCDF}/include
LIBS = -L${NETCDF}/lib -lnetcdf -lnetcdff

OBJS = module_fv3lam_stats.o  module_mpi_arrange.o

all: ens_mean_recenter

ens_mean_recenter: ${OBJS} ens_mean_recenter.o
	mpif90 -f90=ifort -o ens_mean_recenter.exe ${FLAGS} ens_mean_recenter.o ${OBJS} ${LIBS} 

.SUFFIXES : .F90 .f90 .f .c .o

.F90.o :
	mpif90 -f90=ifort ${FLAGS} ${INCLD} -c $<

.f90.o :
	mpif90 -f90=ifort ${FLAGS} ${INCLD} -c $<

clean:
	/bin/rm -f *.o ens_mean_recenter.exe *.mod
