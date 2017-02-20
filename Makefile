### Compilers & flags
F90=gfortran

FFTWLIBS=
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = $(LAPACKLIB)


EXE = exec
F90SRC = main.f90 paramaters.f90 modPlasma.f90 init.f90 MatrixVector.f90 machinePrecision.f90 timeStep.f90 modRecord.f90 Limiter.f90
F90OBJ = main.o parameters.o modPlasma.o init.o MatrixVector.o machinePrecision.o timeStep.o modRecord.o Limiter.o

### Targets
all: $(EXE)
run: $(EXE) 
	./$(EXE)

# Link object files to executables
$(EXE): $(F90OBJ)
	$(F90) -o $(EXE) $(F90OBJ) $(LIBS)

# All .o files depend on the corresponding .f90 file
%.o: %.f90
	$(F90) -c $<

# Dependencies
parameters.o : machinePrecision.o modPlasma.o modRecord.o
main.o: init.o timeStep.o modRecord.o
modPlasma.o : machinePrecision.o
init.o : parameters.o MatrixVector.o modRecord.o
MatrixVector.o : machinePrecision.o
timeStep.o : parameters.o MatrixVector.o modPlasma.o modRecord.o Limiter.o
modRecord.o : machinePrecision.o modPlasma.o MatrixVector.o
Limiter.o : modPlasma.o

clean:
	rm *.o *.mod $(EXE)

.PHONY: all run clean


