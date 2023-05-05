FC = gfortran
FCFLAGS = -O3

OBJS = #mod_share.o mod_math.o mod_io.o


.PHONY: all clean
.SUFFIXES: .f90 .o

all: sVLASOV

sVLASOV: IonBoost72.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFLAGS) $<

%.o: %.mod

clean:
	$(RM) IonBoost72 *.o *.mod
