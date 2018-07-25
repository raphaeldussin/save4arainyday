### Makefile for offline tools (corrections of G. Garric and A. Storto)

# gfortran options
#NCDF= -L/usr/lib/libnetcdf.a -L/usr/lib/libnetcdff.a -I/usr/include -lnetcdf -lnetcdff
#F90=gfortran
#FFLAGS= -O  $(NCDF) -ffree-line-length-none -fconvert=big-endian

# ifort options (meom3pc102)
NCDF= -L/opt/netcdf-4.1.1-ifort/lib -I/opt/netcdf-4.1.1-ifort/include -lnetcdf
F90=ifort
#FFLAGS= -O2 $(NCDF) 
#FFLAGS= -O2 $(NCDF) -openmp
FFLAGS= -O2 $(NCDF) -parallel

EXEC = correc_radsw_ggweights correc_radlw_ggweights correc_precip_ggweights \
       replace_meanfield correc_radsw_ggweights_nopolar correc_precip_astorto \
       correc_radsw_coef

all: $(EXEC)

## rules
correc_radsw_ggweights: correc_radsw_ggweights.f90
	$(F90) correc_radsw_ggweights.f90 -o correc_radsw_ggweights $(FFLAGS)

correc_radsw_ggweights_nopolar: correc_radsw_ggweights_nopolar.f90
	$(F90) correc_radsw_ggweights_nopolar.f90 -o correc_radsw_ggweights_nopolar $(FFLAGS)

correc_radlw_ggweights: correc_radlw_ggweights.f90
	$(F90) correc_radlw_ggweights.f90 -o correc_radlw_ggweights $(FFLAGS)

correc_precip_ggweights: correc_precip_ggweights.f90
	$(F90) correc_precip_ggweights.f90 -o correc_precip_ggweights $(FFLAGS)

correc_precip_astorto: correc_precip_astorto.f90
	$(F90) correc_precip_astorto.f90 -o correc_precip_astorto $(FFLAGS)

replace_meanfield: replace_meanfield.f90
	$(F90) replace_meanfield.f90 -o replace_meanfield $(FFLAGS)

correc_radsw_coef: correc_radsw_coef.f90
	$(F90) correc_radsw_coef.f90 -o correc_radsw_coef $(FFLAGS)

correc_radlw_coef: correc_radlw_coef.f90
	$(F90) correc_radlw_coef.f90 -o correc_radlw_coef $(FFLAGS)

