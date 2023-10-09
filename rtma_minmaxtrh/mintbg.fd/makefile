##############################################################
#     Makefile for urmamaxobs
#     Use:
#       make            - build the executable
#       make clean      - start with a clean slate
##############################################################
# Tunable parameters:
#   FC Name of the Fortran compiling system to use
#   LDFLAGS Options of the loader
#   FFLAGS Options of the compilier
#   DEBUG Options of the compilier included for debugging
#   LIBS List of libraries
#   CMD Name of the executable

FC = ifort
LDFLAGS= -convert big_endian
BINDIR= ../../exec
INC=param.incl
LIBS=
CMD=rtma_mintbg

DEBUG_OPT ?= 1
ifeq ($(DEBUG_OPT), 1)
  DEBUG= -g -traceback
else
  DEBUG=
endif

FFLAGS= -O3 -I $(INC) $(DEBUG)
SRCS= calcminbg.f90 writeminbg.f90 domain_dims.f maxmin.f maxmin_ak.f

# clear out all suffixes
 .SUFFIXES:
# # list only those we use
 .SUFFIXES: .o .f .f90
#

all: $(CMD)

$(CMD): ${SRCS}
	${FC} ${LDFLAGS} -o $(CMD) ${SRCS} ${FFLAGS}

clean:
	/bin/rm -f *.o $(CMD)

install:
	-mv $(CMD) ${BINDR}/
