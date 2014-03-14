
prefix = /usr/local/bin
mpif90 = /opt/local/lib/openmpi/bin/mpif90
f90 = gfortran
flags = -ftree-vectorize -O3 -fopenmp -ffree-line-length-none
lapacklib = /usr/local/lib

fit : fit_charges.f90 dcd.f90
	$(f90) -c fit_charges.f90 dcd.f90  $(flags) -L$(lapacklib) -llapack -lblas
	$(f90)  fit_charges.o dcd.o -o fit_charges.omp.x  $(flags) -L$(lapacklib) -llapack -lblas
	cp fit_charges.omp.x /Users/mmccull/Dropbox/charge_fitting/test_systems/actin_monomer_atp_fil

clean:
	rm -f *.o *.mod *.x

