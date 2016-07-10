FCC=ifort
PLIB=ar cr
LIBPATH=/home/cj/bin/lib/
MOD= lat_par.o
OBJS= chain.o \
	  square.o \
	  honeycomb.o \
	  triangular.o \
	  kagome.o \
	  cubic.o \
	  pyrochlore.o \
	  produce_lattice.o

DEG=
#DEG=-g

TAGE=produce_lattice.a
all: $(MOD) $(TAGE) INSTALL

$(MOD): lat_par.f90
	$(FCC) $(DEG) -c lat_par.f90

$(TAGE): $(OBJS)
	$(PLIB) $(TAGE) $(MOD) $(OBJS)

$(OBJS): %.o: %.f90
	$(FCC) $(DEG) -c $< -o $@

INSTALL:
	cp -f $(TAGE) $(LIBPATH)


#TAGE=produce_lattice
#DEG=-g
#all: $(MOD) $(TAGE)
#
#$(MOD): lat_par.f90
#	$(FCC) $(DEG) -c lat_par.f90
#$(OBJS): %.o: %.f90
#	$(FCC) $(DEG) -c $<
#$(TAGE): $(OBJS) $(MOD)
#	$(FCC) $(DEG) $(MOD) $(OBJS) -o $(TAGE)
#


clean:
	rm -f $(MOD) $(OBJS) *.mod

.PHONY: all clean
