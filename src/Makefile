FC = gfortran
DBG := -O0 -g -fbacktrace -Wall -Wextra
FFLAGS := $(DBG)
OBJS = types.o share.o solver.o param.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: mora

mora: main.f90 $(OBJS)
	$(FC) $(FFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FFLAGS) $<

%.o: %.mod

clean:
	$(RM) mora *.o *.mod
