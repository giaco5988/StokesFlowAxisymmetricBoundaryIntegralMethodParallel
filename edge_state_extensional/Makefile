# Makefile to compile Rising droplet, 26.01.2016
# The sharp (#) symbols indicate comment, which are not executed

# Please adapt the following lines 
# compiler of choice: g++, icc, etc.
CPP = g++

# compiler options for g++ with Xcode
OPTS = -lstdc++ -framework Accelerate -O3 -lpthread
#OPTS = -lstdc++ -framework Accelerate -lpthread

# compiler options for g++ with Lapack libraries
#OPTS = -fopenmp -lstdc++ -llapack -lblas -lpthread
#OPTS = -lstdc++ -llapack -lblas -lpthread

# compiler options for icc and mkl (Intel)
# OPTS = -lstdc++ -mkl -O3 -openmp 

# library path for instance for LAPACK, if needed
# LIBS = -L/opt/intel/composer_xe_2013.1.119/mkl/lib/

# include path for headers, if needed
# INCS = -I/Developer/SDKs/MacOSX10.7.sdk/usr/lib/gcc/i686-apple-darwin11/4.2.1/include

# That is all, only edit the lines below for new source files
# or to add user defined cases to the compiler (example: user1)

# source files except main classes
SRC = ../source/matrix.cpp ../source/remesh.cpp ../source/combase.cpp ../source/geometrydrop.cpp ../myLIB/alglib.a

# source files for compiler test
#TRC = combase_cote.cpp bndblock.cpp matrixblock.cpp bndelement.cpp matblc1.cpp -w

# Lines hereafter don’t need modification unless you want to add your own cases. 
 
#all: tutorial1 tutorial2 tutorial3 tutorial4 tutorial5

edge_state: main.cpp
	$(CPP) -o edge_state_extensional.o main.cpp $(LIBS) $(INCS) $(SRC) $(OPTS)
	@if [ $$? = 0 ]; then \
	  echo "  Compiling $@ done" ; \
	else \
	  echo " Compiling $@ failed" ; \
	fi
comp_test: main_channel.cpp
	$(CPP) -o ./test.o main.cpp $(TRC) $(INCS) $(LIBS) 
	@if [ $$? = 0 ]; then \
	  echo "  Compiling without OpenMP and Lapack successful "; \
	else \
	  echo " Compiling failed, check the error messages "; \
	fi
	rm -rf ./test.o
clean:
	rm -rf ../*.o
