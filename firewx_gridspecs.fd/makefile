LIBS =          ${W3NCO_LIB4} ${W3EMC_LIB4} 
FC=		ftn
LDFLAGS= -convert big_endian
ESSL=
SRCS=		firewx_gridspecs.f	
OBJS=		firewx_gridspecs.o 
CMD=		${OUTDIR}/firewx_gridspecs.x

all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS) $(ESSL) -o $(@) $(OBJS) $(LIBS)
clean:
	rm -rf $(OBJS) $(CMD)

