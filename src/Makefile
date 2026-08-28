FC = gfortran
DBG := -O0 -g -fbacktrace -Wall -Wextra
FFLAGS := $(DBG)
OBJS = types.o config.o linear_solver.o electrons.o grid.o ions.o fields.o pusher.o diagnostics.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: ionboost

ionboost: main.f90 $(OBJS)
	$(FC) $(FFLAGS) $< $(OBJS) -o $@

types.o: types.f90
config.o: config.f90 types.o
linear_solver.o: linear_solver.f90 types.o
electrons.o: electrons.f90 types.o
grid.o: grid.f90 types.o config.o
ions.o: ions.f90 types.o grid.o
fields.o: fields.f90 types.o grid.o ions.o electrons.o linear_solver.o
pusher.o: pusher.f90 types.o grid.o ions.o fields.o linear_solver.o
diagnostics.o: diagnostics.f90 types.o grid.o ions.o fields.o

.f90.o:
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) ionboost *.o *.mod
