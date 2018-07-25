#
FC=ifort
FF=-O3 -CB -traceback -module mod/
#
EZCDF=-module $(HOME)/TOOLS/sosie/mod $(HOME)/TOOLS/sosie/obj/io_ezcdf.o $(HOME)/TOOLS/sosie/obj/drwn.o
#
NCDF=/opt/netcdf-4.1.1
#
#
INC=-I$(NCDF)/include
LIB=-L$(NCDF)/lib -lnetcdf
#
#
All:  mnth_offset.x
#
mnth_offset.x: src/mnth_offset.f90
	$(FC) $(FF) $(INC) $(EZCDF) src/mnth_offset.f90 -o mnth_offset.x $(LIB)
#
#
clean:
	rm -f *.o *.nc *.dat
	rm -f *.x *.o*
#
