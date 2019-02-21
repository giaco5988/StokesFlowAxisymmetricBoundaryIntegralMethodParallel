/*
Rising droplet BEM Giacomo Gallino
Date: 26.01.2016
*/

static char help[] = "Basic vector routines.\n\n";

//libraries
#include "../sourceMPI/geometryMPI.h"
#include "../sourceMPI/remeshMPI.h"
#include "../sourceMPI/matrixMPI.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <petscksp.h>


// The program begins
int main (int argc, char ** argv) {

  time_t t0;

  time(&t0);  /* get current time; same as: now = time(NULL)  */

  //direcory where saving results
  char res [546] = "/Users/Giacomo/Documents/C++/drop_extensionalMPI/results/";

  //droplets parameters
  double x0 = 0;             //center of mass x coordinate
  double y0 = 0;             //center of mass y coordinate
  int n = 50;               //number of elements
  double D = 0.02;            //ellipticity
  double Ca = 0.05;             //capillary number
  double lambda = 1;         //viscosity ratio

  //simulations parameters
  double dt = 0.1;          //time step
  int loop = 1000;              //number of loop
  int RK = 2;                  //choose order of Runge Kutta
  //int PETSc = 1;               //choose to use PETSc
  int checkpoint = 1;         //when save data
  int VolCorr = 2;            //volume correction

  //remesh
  int UseRemesh = 1;         //if 0 use mesh stab, if 1 use remesh with distribution
  int DoRemesh = 10;         //how often to perform remesh

  //start petsc
  PetscInitialize(&argc,&argv,(char*)0,help);

  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  //PetscOptionsGetInt(NULL,"-view_randomvalues",&rank,NULL);
  PetscPrintf(PETSC_COMM_SELF,"%i \n",rank);

  //options, if 1 override data if directory already exists
  int override = 1;

  //output, starting message
  printf("Axisymmetric drop in extensional flow starts with %i elements, Ca=%f and lambda=%f\n",n,Ca,lambda);

  /////////////////////////////// SAVINGS DATA OPERATIONS //////////////////////////////////
  //create folder and save data to this folder
  char respath [246];
  sprintf(respath,"%sCa=%f_lambda=%f_delta=%f_elem=%i_dt=%f_loop=%i_RK=%i",res,Ca,lambda,D,n,dt,loop,RK);
  mkdir(respath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (mkdir(respath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)&(1-override)){
    printf("Saving directory already existing, kill simulation\n");
    exit(0);
  }

  //save options
  char saveTo[246];
  sprintf(saveTo,"%s/ParametersOptions.txt",respath);

  FILE* pfile;
  pfile = fopen(saveTo,"w");
  fprintf(pfile,"Ca  lambda  delta  elem  dt loop  RK checkpoint\n"); //header
  fprintf(pfile,"%f  %f  %f  %i  %f  %i %i %i \n",Ca, lambda, D, n, dt, loop, RK, checkpoint); //options data
  fclose(pfile);
  /////////////////////////////// END SAVINGS DATA OPERATIONS //////////////////////////////////

  ////////////////////////////// GEOMETRY OBJECT INITAILIZATIONS  /////////////////////////////
  //initialize drop object
  GeometryMPI drop(x0, y0, D, n);
  ////////////////////////////// END GEOMETRY OBJECT INITAILIZATIONS  /////////////////////////

  ////////////////////////////// REMESH OBJECT INITILAIZATION /////////////////////////////////

  //HERE REMESH

  ////////////////////////////// END REMESH OBEJECT INITIALIZATION /////////////////////////////

  //here starts the simulation
  /////////////////////////////  SOLUTION OBJECT INITIALIZATION AND RUN  /////////////////////////////
  //initialize solution object
  MatrixMPI solveDrop(&drop, lambda, Ca, 2);

  //set type of remesh
  solveDrop.SetMeshKind(UseRemesh);

  //set how often global remesh is done                                                                                                                                                                                                   
  solveDrop.SetDoRemeshDistribution(DoRemesh);

  //set checkpoint for data saving
  solveDrop.SetCheckpoint(checkpoint);

  //set volume correction
  solveDrop.SetVolumeCorr(VolCorr);

  //set time step and number of loops
  solveDrop.SetLoopDt(loop, dt);

  //choose output directory
  solveDrop.SetResdir(respath);

  //solveDrop.BuildAndSolvePETSc();
  //solveDrop.BuildAndSolvePETSc();
  if (RK==1){
    //solve the system and advance with RK1
    solveDrop.SolveRK1();
  } else if (RK==2){
    //solve the system and advance with RK2
    solveDrop.SolveRK2();
  }
  /////////////////////////////  END SOLUTION OBJECT OPERATIONS ///////////////////////////////

  //close PETSc
  PetscFinalize();

  //determine time of the simulation 
  time_t t1;

  time(&t1);  /* get current time; same as: now = time(NULL)  */

  double seconds = difftime(t1,t0);

  printf ("Time of the simulation is %.f seconds \n", seconds);

// End of program
    std::cout << "Finished and Exiting\n"; 
    return 0;
  
}
