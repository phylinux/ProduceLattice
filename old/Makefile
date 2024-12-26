#FCC=ifort
#FC=-g -CU -CB

FCC=gfortran
FC=-O0 -fbounds-check -Wuninitialized

all:
	$(FCC) $(FC) produce_lattice.f90 -o produce_lattice

clean:
	rm -f *.mod produce_lattice
