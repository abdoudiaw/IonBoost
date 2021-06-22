FC = gfortran
# flags
DBG := -O0 -g -fbacktrace -Wall -Wextra #-pedantic -fcheck=all -finit-real=snan -ffpe-trap=invalid -std=f2018
FFLAGS := $(DBG)
#
OBJS = types.o share.o solver.o param.o
#
.PHONY: all clean
.SUFFIXES: .f90 .o

all: Ionboost

Ionboost: IonBoost28.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

%.o: %.mod

clean:
	$(RM) IonBoost28 *.o *.mod
