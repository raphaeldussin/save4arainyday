#############################################################################
### Makefile for FARC compiled codes

include macro.meolkerg
#include macro.yourmachine


#############################################################################

EXEC = global_timeserie trend

all: $(EXEC)

clean:
	\rm global_timeserie trend trend.o

## rules
global_timeserie: global_timeserie.f90
	$(F90) global_timeserie.f90 -o global_timeserie $(FFLAGS)

trend: trend.o
	$(CC) trend.o $(LFLAGS) -o trend

trend.o: trend.c
	$(CC) $(CFLAGS) -c trend.c

