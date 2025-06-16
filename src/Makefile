SRCDIR = .
FC = gfortran
FFLAGS = -O2 -fopenmp -fbounds-check
OBJS = grid.o waveprop.o windmod.o physics.o io.o params.o main.o 

all: wvmod

wvmod: $(OBJS) 
	$(FC) $(FFLAGS) -o wave_model $(OBJS)

grid.o: $(SRCDIR)/grid.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/grid.f90

waveprop.o: $(SRCDIR)/waveprop.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/waveprop.f90

windmod.o: $(SRCDIR)/windmod.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/windmod.f90

physics.o: $(SRCDIR)/physics.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/physics.f90

io.o: $(SRCDIR)/io.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/io.f90

params.o: $(SRCDIR)/params.f90
	$(FC) $(FFLAGS) -c $(SRCDIR)/params.f90

clean:
	rm -f $(OBJS) wave_model
	rm -rf $(SRCDIR)/../output/diag/*txt
	rm -rf $(SRCDIR)/../output/spc/*txt

main.o: $(SRCDIR)/main.f90 grid.o waveprop.o physics.o windmod.o io.o params.o
	$(FC) $(FFLAGS) -c $(SRCDIR)/main.f90