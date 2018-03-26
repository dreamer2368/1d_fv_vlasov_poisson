### Compilers & flags
F90=mpifort

TOPDIR = $(dir $(firstword $(MAKEFILE_LIST)))
SRCDIR = src
OBJDIR = obj
MODDIR = mod

FFTWLIBS=
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = 


EXE = debye
F90SRC = main.f90 \
		constants.f90 \
		modInputHelper.f90 \
		modMPI.f90 \
		MatrixVector.f90 \
		modPlasma.f90 \
		modCircuit.f90 \
		modCircuitBC.f90 \
		Limiter.f90 \
		modPlasmaBC.f90 \
		modSource.f90 \
		modQoI.f90 \
		modRecord.f90 \
		timeStep.f90 \
		init.f90 \
		testmodules.f90
F90OBJ = $(F90SRC:%.f90=$(OBJDIR)/%.o)

### Targets
all: dir $(EXE)
run: $(EXE) 
	./$(EXE)
dir: $(OBJDIR) $(MODDIR)

# Sub-directories
$(OBJDIR):
	mkdir $(TOPDIR)$(OBJDIR)
$(MODDIR):
	mkdir $(TOPDIR)$(MODDIR)

# Link object files to executables
$(EXE): $(F90OBJ)
	@echo $(dir $(firstword $(MAKEFILE_LIST)))
	$(F90) -o $(EXE) $(F90OBJ) $(LIBS)

# All .o files depend on the corresponding .f90 file
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(F90) -J$(MODDIR) -c -o $@ $<

# Dependencies
$(OBJDIR)/modMPI.o : $(OBJDIR)/constants.o
$(OBJDIR)/MatrixVector.o : $(OBJDIR)/constants.o
$(OBJDIR)/Limiter.o : $(OBJDIR)/constants.o
$(OBJDIR)/modInputHelper.o : $(OBJDIR)/modMPI.o
$(OBJDIR)/modPlasmaBC.o : $(OBJDIR)/Limiter.o
$(OBJDIR)/modCircuit.o : $(OBJDIR)/MatrixVector.o
$(OBJDIR)/modCircuitBC.o : $(OBJDIR)/modCircuit.o
$(OBJDIR)/modPlasma.o : $(OBJDIR)/MatrixVector.o \
						$(OBJDIR)/modPlasmaBC.o
$(OBJDIR)/modSource.o : $(OBJDIR)/modPlasma.o \
						$(OBJDIR)/modCircuit.o
$(OBJDIR)/modQoI.o : $(OBJDIR)/modPlasma.o \
					$(OBJDIR)/modCircuit.o
$(OBJDIR)/modRecord.o : $(OBJDIR)/modPlasma.o \
						$(OBJDIR)/modCircuit.o
$(OBJDIR)/timeStep.o : $(OBJDIR)/modQoI.o \
						$(OBJDIR)/modRecord.o \
						$(OBJDIR)/modCircuitBC.o \
						$(OBJDIR)/modSource.o
$(OBJDIR)/init.o : $(OBJDIR)/modPlasma.o \
					$(OBJDIR)/modCircuitBC.o \
					$(OBJDIR)/modRecord.o
$(OBJDIR)/testmodules.o : $(OBJDIR)/init.o \
							$(OBJDIR)/timeStep.o \
							$(OBJDIR)/modRecord.o
$(OBJDIR)/main.o: $(OBJDIR)/testmodules.o \
					$(OBJDIR)/modInputHelper.o

clean:
	rm -rf $(OBJDIR) $(MODDIR) $(EXE)

.PHONY: dir all run clean
