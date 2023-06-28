FC = ifort
FFLAGS = -fopenmp #-warn

all: remap_dwinds.so

remap_dwinds.so: remap_dwinds.f90
	$(FC) $(FFLAGS) -c $<
	f2py3.8 -DF2PY_REPORT_ATEXIT --build-dir . -c --f90flags='-fopenmp' -L/pe/intel/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin/libiomp5.so -liomp5 -m remap_dwinds remap_dwinds.f90

.PHONY: clean
clean:
	rm -f *.so *.mod *.o

