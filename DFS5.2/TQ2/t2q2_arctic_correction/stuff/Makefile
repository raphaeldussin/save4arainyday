# Fortran compiler :
FC = ifort
# #
# # Root directory for Netcdf :
NCDF = /opt/netcdf-4.1.1-ifort
# #
# # Fortran compilation flags :
# #FF= -O3 -fp-model precise -module mod/
FF= -O3 -fp-model precise -traceback -CB 
# #
#

NCDF_INC=$(NCDF)/include


all: io_ezcdf.o drwn.o

io_ezcdf.o: io_ezcdf.f90
	$(FC) $(FF) -I$(NCDF_INC) -c io_ezcdf.f90 -o io_ezcdf.o

drwn.o: drwn.F90
	$(FC) $(FF) -c drwn.F90 -o drwn.o

