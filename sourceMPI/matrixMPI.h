//
//  matrix.h
//
//  Created by Giacomo Gallino on 26.01.16.
//

#ifndef matrixMPI_h
#define matrixMPI_h

#include "geometryMPI.h"
#include "petscksp.h"

class MatrixMPI {

  //these are the main charecteristic of the object
  int type; //option for topology (drop, 2 drops, drop and walls etc etc)
  int BCkind; //option for kind of BC 
  int nG; //pointer to the elements of different parts
  int MeshKind; //decide which remesh to use
  //double *A; //matrix
  double *RHS; //rhs
  GeometryMPI *DropHere;
  double lambda; //global variable for the viscosity ratio
  double Ca; //global variable for the capillary number
  int loop; //number of loops
  double dt; //time step
  char *resdir; //directory where to save the files
  double CaExtensional; //capillary number for extensioanl flow
  int checkpoint; //decide when to save the data
  int VolCorr; //decide if to do apply volume correction
  int DoRemeshDistribution; //decide how often execute global remesh
  
  //PETSc create petsc complex objects like vectors matrices and oinear solver                                                                                        
  Vec sol, b, fake, rhs, transfer;   //vectors                                                                                                                          
  Mat A, U;   //matrices                                                                                                                                        
  KSP ksp;    //linear solver object                                                                                                                    
  PetscInt NODES;
  VecScatter scatter;
  PC pc; //preconditioner

public:

  //constructor
  MatrixMPI (GeometryMPI*, double, double, int);
  //get topology type
  int GetType();
  //build matrix
  void Build();
  //build matrix in parallel with petsc                                                                                                                                                                                                   
  void BuildAndSolvePETSc();
  //solve matrix
  void Solve(int);
  //set number of loops and time step
  void SetLoopDt(int, double);
  //dolve and advance in time with RK1
  void SolveRK1();
  //dolve and advance in time with RK2
  void SolveRK2();
  //set directory where to save the files
  void SetResdir(char*);
  //set number of threads in parallel computation
  void SetThreads(int);
  //set checkpoint
  void SetCheckpoint(int);
  //activate volume correction
  void SetVolumeCorr(int);
  //set do remesh distribution
  void SetDoRemeshDistribution(int);
  //set remesh kind
  void SetMeshKind(int);
    
};


#endif
