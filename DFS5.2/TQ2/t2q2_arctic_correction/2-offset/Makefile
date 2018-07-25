#
FC=ifort
FF=-O3 -CB -traceback -module mod/
#FF=-g -CB -traceback -module mod/
#
SOSIEROOT=/home/users/molines/sosie_sourceforge
STUFFROOT=/home/users/dussin/DFS5/DFS5.2/TQ2/t2q2_arctic_correction/stuff
#EZCDF=-module $(SOSIEROOT)/mod $(SOSIEROOT)/obj/io_ezcdf.o $(SOSIEROOT)/obj/drwn.o
EZCDF=-module $(SOSIEROOT)/mod $(STUFFROOT)/io_ezcdf.o $(SOSIEROOT)/obj/drwn.o
#
NCDF=/opt/netcdf-4.1.1-ifort
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
