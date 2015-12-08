
MFile=produce_lattice
FModule=lat_par
PFile=plotlat

$(MFile): $(FModule).o $(MFile).o
	ifort $(FModule).o $(MFile).o -o $(MFile)
	#ifort -g $(FModule).o $(MFile).o -o $(MFile)

$(MFile).pl: $(FModule).o $(PFile).o $(MFile).o
	   ifort $(FModule).o $(PFile).o $(MFile).o -o $(MFile) ~/bin/lib/libvogle.a -lX11

$(FModule).o: $(FModule).f90
	ifort -c $(FModule).f90
	#ifort -g -c $(FModule).f90

$(PFile).o: $(PFile).f90
	ifort -c $(PFile).f90
	#ifort -g -c $(PFile).f90

$(MFile).o: $(MFile).f90
	ifort -c $(MFile).f90
	#ifort -g -c $(MFile).f90

run:
	./$(MFile)

cpto:
	cp $(MFile) ~/Documents/Research/WORKS/Hard-Core-Bose-Hubbard-In-Pyrochlore/

.PHONY: clean

clean:
	-rm *.o *.mod *~ $(MFile)
	-rm *.dat *.list lattice
