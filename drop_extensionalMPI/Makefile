#include path for petsc
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

CLINKER = g++

extens: main.o ../sourceMPI/geometryMPI.o ../sourceMPI/combaseMPI.o ../sourceMPI/remeshMPI.o ../sourceMPI/matrixMPI.o chkopts
	-${CLINKER} -o $@ main.o ../sourceMPI/geometryMPI.o ../sourceMPI/combaseMPI.o ../sourceMPI/remeshMPI.o ../sourceMPI/matrixMPI.o ${PETSC_LIB} ../myLIB/alglib.a
	${RM}$
