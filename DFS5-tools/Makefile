### Makefile for GIG correction

# portable
#NCDF= -L/usr/lib/libnetcdf.a -L/usr/lib/libnetcdff.a -I/usr/include -lnetcdf -lnetcdff
#F90=gfortran
#FFLAGS= -O  $(NCDF) -ffree-line-length-none -fconvert=big-endian

NCROOT=/opt/netcdf-4.1.1-ifort
NCDF= -L$(NCROOT)/lib -I$(NCROOT)/include -lnetcdf
F90=ifort
FFLAGS= -O2 $(NCDF) 



EXEC = correc_radsw_ggweights correc_radlw_ggweights correc_precip_ggweights correc_radsw_ggweights_nopolar correc_precip_astorto \
       detrend_precip detrend_precip_v2 rescale_precip blend_precip

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

detrend_precip: detrend_precip.f90
	$(F90) detrend_precip.f90 -o detrend_precip $(FFLAGS)

detrend_precip_v2: detrend_precip_v2.f90
	$(F90) detrend_precip_v2.f90 -o detrend_precip_v2 $(FFLAGS)

rescale_precip: rescale_precip.f90
	$(F90) rescale_precip.f90 -o rescale_precip $(FFLAGS)

blend_precip: blend_precip.f90
	$(F90) blend_precip.f90 -o blend_precip $(FFLAGS)

