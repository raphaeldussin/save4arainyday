#
FC=ifort
#FF=-O3 -CB -traceback
FF=-O3
#
EZCDF=-module $(HOME)/TOOLS/sosie/mod $(HOME)/TOOLS/sosie/obj/io_ezcdf.o
NCDF=/opt/netcdf-4.1.1
#
LIB=-L$(NCDF)/lib -lnetcdf
INC=-I$(NCDF)/include
#
#
#
All: mk_q2_dfs4.x mk_t2_dfs4.x
#All: mk_t2_dfs4.x
#
mk_q2_dfs4.x: src/mk_q2_dfs4.f90
	$(FC) $(FF) $(EZCDF) src/mk_q2_dfs4.f90 -o mk_q2_dfs4.x $(LIB)
#
mk_t2_dfs4.x: src/mk_t2_dfs4.f90
	$(FC) $(FF) $(EZCDF) src/mk_t2_dfs4.f90 -o mk_t2_dfs4.x $(LIB)
#	
#
clean:
	rm -f  *.nc out  *.x *.o* *.mod
#
