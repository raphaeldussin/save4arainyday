### Makefile for GIG correction

# portable
#NCDF= -L/usr/lib/libnetcdf.a -L/usr/lib/libnetcdff.a -I/usr/include -lnetcdf -lnetcdff
#F90=gfortran
#FFLAGS= -O  $(NCDF) -ffree-line-length-none -fconvert=big-endian

NCROOT=/opt/netcdf-4.1.1-ifort
NCDF= -L$(NCROOT)/lib -I$(NCROOT)/include -lnetcdf
F90=ifort
FFLAGS= -O2 $(NCDF) 



EXEC = correct_t2

all: $(EXEC)

## rules
correct_t2: correct_t2.f90
	$(F90) correct_t2.f90 -o correct_t2 $(FFLAGS)

