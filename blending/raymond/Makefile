FC = ifort

remap.so: raymond.f
	$(FC) $(FFLAGS) -c $<
	f2py3.8 -DF2PY_REPORT_ATEXIT --build-dir . -c -m raymond raymond.f

.PHONY: clean
clean:
	rm -f *.so *.o

