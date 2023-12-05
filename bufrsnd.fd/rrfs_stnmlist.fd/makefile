SHELL=/bin/sh
#
# This makefile was produced by /usr/bin/fmgen at 11:24:59 AM on 09/25/96
# If it is invoked by the command line
#	make -f makefile
# it will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and produce an executable named a.out.
#
# If it is invoked by the command line
#	make -f makefile a.out.prof
# it will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and produce an executable which profiles
# named a.out.prof.
#
# To remove all the objects but leave the executables use the command line
#	make -f makefile clean
#
# The parameters SRCS and OBJS should not need to be changed.  If, however,
# you need to add a new module add the name of the source module to the
# SRCS parameter and add the name of the resulting object file to the OBJS
# parameter.  The new modules are not limited to fortran, but may be C, YACC,
# LEX, or CAL.  An explicit rule will need to be added for PASCAL modules.
#
SRCS=	stnmlist.f  

OBJS=	stnmlist.o 
OBJS_DEBUG=	stnmlist.o 

# Tunable parameters
#
# FC		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling
#
#FC =		${COMP}
FC =		${COMPILER}
FC_DEBUG =	${COMPILER}  -g -traceback  -ftrapuv  -check all  -fp-stack-check
LDFLAGS =	
LIBS = ${W3NCO_LIB4} ${W3EMC_LIB4} ${BUFR_LIB4}

CMD =	        rrfs_stnmlist.x
PROFLIB =


# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:		$(CMD)

$(CMD)::		$(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)
debug::		$(OBJS_DEBUG)
	$(FC_DEBUG) $(LDFLAGS) -o hireswfv3_stnmlist.x $(OBJS) $(LIBS)
$(OBJS_DEBUG):
	$(FC_DEBUG)  -c  $(@:.o=.f)


clean:
	-rm -f $(OBJS)
