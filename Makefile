### Compilers & flags
F90=gfortran

FFTWLIBS=
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = 


EXE = exec
F90SRC = main.f90 constants.f90 MatrixVector.f90 modPlasma.f90 Limiter.f90 modQoI.f90 modRecord.f90 timeStep.f90 init.f90 testmodules.f90
F90OBJ = main.o constants.o MatrixVector.o modPlasma.o Limiter.o modQoI.o modRecord.o timeStep.o init.o testmodules.o

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
MatrixVector.o : constants.o
modPlasma.o : constants.o
Limiter.o : constants.o
modQoI.o : modPlasma.o MatrixVector.o
modRecord.o : modPlasma.o MatrixVector.o
timeStep.o : modQoI.o modRecord.o Limiter.o
init.o : modPlasma.o modRecord.o
testmodules.o : init.o timeStep.o modRecord.o
main.o: testmodules.o


clean:
	rm *.o *.mod $(EXE)

.PHONY: all run clean


