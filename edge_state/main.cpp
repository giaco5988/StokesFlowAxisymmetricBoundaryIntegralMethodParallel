/*
Rising droplet BEM Giacomo Gallino
Date: 26.01.2016
*/

//libraries
#include "../source/geometrydrop.h"
#include "../source/remesh.h"
#include "../source/matrix.h"
#include <iostream>
//#include <string>
//#include <fstream>
#include <stdio.h>
#include <stdlib.h>
//#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>



// The program begins
int main (int argc, char * const argv[]) {

  time_t t0;

  time(&t0);  /* get current time; same as: now = time(NULL)  */

  //direcory where saving results
  char res [246] = "/Users/Giacomo/Documents/C++/edge_state/results/";

  //droplets parameters
  double x0 = 0;             //center of mass x coordinate
  double y0 = 0;             //center of mass y coordinate
  int n = 500;               //number of elements
  int IDshape = 1;            //ellipticity
  double Ca = 6;             //capillary number
  double lambda = 1;         //viscosity ratio

  //shape parameters
  int IDdelta = 1;           //counter for delta
  double StepDelta = 1e-6;   //step on search in delta
  double DeltaStart = 0.03;  //starting deformation for shape0
  int THISround = 1;              //round of smaller Delta
  int IteGetShape = 10;       //iteration of previous simulation form which I get the shape
  double Dphi = 0.5;

  //numericss parameters
  double dt = 0.001;          //time step
  int loop = 25000;              //number of loop
  int checkpoint = 100;        //data saving
  int DoRemesh = 100;		//how often do global remesh
  int RK = 2;                  //choose order of Runge Kutta
  int CPUs = 1;                //choose number of CPUs to use
  int VolCorr = 2;             //if volume correction is used
  int extrapolate = 1;         //decide if to interpolate or extrapolate for new shape
  double extraMulti = 3.0;      //safety risizing for eztrapolation

  //mesh control parameters                                                                                                                                                                                         
  int MeshType = 1;         // if 0 is mesh stabilization, if 1 is classic mesh control                                                                                                                                 
  double Adaptivity = 1.5;  //set mesh stab adaptivity                                                                                                                                                                    
  int Cluster = 3;          //FOR MESH STAB set clustering option: //1 is classic based on curvature, 2 clustering close to axis, 3 on the right part of the droplet                                                       
  int WhichDistr = 4;       //FOR GLOBAL REMESHING: 1 and 2 are on right part of droplet, 3 is uniform, 4 is very clustered to the right (exponential)                                                                           
  double AdaptDistribution = 0.1; //set distribution adaptivity 

  //options, if 1 override data if directory already exists
  int override = 1;

  if (argc>2){

    IDshape = atof(argv[1]);
    IDdelta = atof(argv[2]);
    THISround = atof(argv[3]);

  }

  //output, starting message
  printf("Edge state of rising droplet for IDshape=%i IDdelta=%i Round=%i starts with %i elements, Ca=%f and lambda=%f\n",IDshape,IDdelta,THISround,n,Ca,lambda);

  /////////////////////////////// SAVINGS DATA OPERATIONS //////////////////////////////////
  //create folder and save data to this folder
  char respath [246];
  
  //save path for loading
  char load1 [246]; char load2 [246];
  sprintf(load1,"%sCa=%f_lambda=%f_",res,Ca,lambda);
  sprintf(load2,"elem=%i_dt=%f_loop=%i_RK=%i_CPUs=%i",n,dt,loop,RK,CPUs);
  /////////////////////////////// END SAVINGS DATA OPERATIONS //////////////////////////////////

  ////////////////////////////// GEOMETRY OBJECT INITAILIZATIONS  /////////////////////////////
  //initialize drop object
  GeometryDrop drop(IDshape, IDdelta, DeltaStart, StepDelta, n);

  ////////////////////////////////////////////////////////////////////////////////////////////

  //remesh object                                                                                                                                                                                                    
  Remesh MeshDrop(&drop, MeshType);
  MeshDrop.SetAdaptDrop(Adaptivity);//adaptivity                                                                                                                                                                      
  MeshDrop.SetRK(RK); //I'm using RK2                                                                                                                                                                                 
  MeshDrop.SetOptionCluster(Cluster); //1 is classic based on curvature, 2 clustering close to axis, 3 on the right part of the droplet                                                                                         
  MeshDrop.SetTuneDistributionAdaptivity(AdaptDistribution);
  MeshDrop.SetOptionDistribution(WhichDistr);

  ///////////////////////////// SOLUTION OBJECT INITIALIZATION ////////////////////////////////
  //initialize solution object BCkind 0-normal, 1-gravity, 2-extensional flow
  Matrix solveDrop(&drop, &MeshDrop, lambda, Ca, 1);

  /////////////////////////////////////////////////////////////////////////////////////////////

  //create folder and save data to this folder                                                                                                                                                                                              
  sprintf(respath,"%sCa=%f_lambda=%f_IDshape=%i_IDdelta=%i_Round=%i_elem=%i_dt=%f_loop=%i_RK=%i_CPUs=%i",res,Ca,lambda,IDshape,IDdelta,THISround,n,dt,loop,RK,CPUs);
  if (mkdir(respath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)){

    if (1-override) {
      printf("Saving directory already existing, kill simulation\n");
      exit(0);
    }

  }

  /////////////////////////////////   UPLOAD PREVIOUS DATA   //////////////////////////////////
  //if first shape set round.txt=1
  if (IDshape==1){

    drop.SetRound(1,respath);

  } else if (IDshape>1) { //if it after first shape geometry is uplaoded from previous simulation

    //initialize temporary variables                                                                                                                                                          
    double *xStable, *yStable, *xUnstable, *yUnstable;
    xStable = new double [n+1];
    yStable = new double [n+1];
    xUnstable = new double [n+1];
    yUnstable = new double [n+1];

    printf("Shape is upload from previous simulation");

    //if round==1 I have to look at the previous shape
    if (THISround==1) {

      //check which is good round of previous shape
      int IDround = drop.SetPreviousRound(load1, load2, IDshape-1);
      //get first unstable shape
      int FirstUnstable = drop.GetFirstUnstable(IDshape-1, IDround, load1, load2);
      //check if step is not too big and eventually cahnge it
      IteGetShape = drop.CheckIteGetShape(IDshape-1, FirstUnstable, IDround, IteGetShape, load1, load2);
      //get shapes at a certain time, first stable and first unstable                                                                                                       
      drop.GetShapesAtTime(IDshape-1, FirstUnstable, IDround, IteGetShape, load1, load2, xStable, yStable, xUnstable, yUnstable);

    } else if (THISround>1) { //in this case I have to look at previous round

      //check if previous round was already fine enough                                                                                                                                                                   
      drop.CheckPreviousRound(IDshape, IDdelta, THISround, load1, load2);
      //get first unstable shape                                                                                                                                                                                                           
      int FirstUnstable = drop.GetFirstUnstable(IDshape, THISround-1, load1, load2);
      //get shapes at a certain time, first stable and first unstable                                                                                                                                         
      drop.GetShapesAtTime(IDshape, FirstUnstable, THISround-1, 0, load1, load2, xStable, yStable, xUnstable, yUnstable);

    }

  if (extrapolate==0) {
    //once the shapes is chosen, do weighted average of the last stable nad first unstable
    drop.AverageShapes(Dphi, IDdelta, xStable, yStable, xUnstable, yUnstable);
  } else if (extrapolate==1) {
    //extrapolate shape from stable one
    double deltaS = drop.ComputeDelta(1.0, 1, xStable, yStable, xUnstable, yUnstable);
    if (THISround>1){
      extraMulti = 1.0;
    }
    double deltaU = extraMulti*drop.ComputeDelta(1.0, 2, xStable, yStable, xUnstable, yUnstable);
    drop.ExtrapolateShapes(Dphi, IDdelta, xStable, yStable, deltaS, deltaU, xUnstable, yUnstable);
    drop.AverageShapes(Dphi ,IDdelta, xStable, yStable, xUnstable, yUnstable);
  }

  //compute the norm of the perturbation which is equal to the norm of R(theta)-1                                                                                                                                        
  double delta1  = drop.ComputeDelta(Dphi, 1, xStable, yStable, xUnstable, yUnstable);
  double delta2  = drop.ComputeDelta(Dphi, 2, xStable, yStable, xUnstable, yUnstable);

  //check if delta step is under the desired one                                                                                                                                                          
  drop.CheckDeltaStep(delta1, delta2, StepDelta, load1, load2, IDshape, IDdelta, THISround);

  printf("deltaStep=%1.16f \n",delta2-delta1);

  //free memory
  delete[] xStable;
  delete[] yStable;
  delete[] xUnstable;
  delete[] yUnstable;
  }

  //save options                                                                                                                                                                                                                          
  char saveTo[246];
  sprintf(saveTo,"%s/ParametersOptions.txt",respath);

  FILE* pfile;
  pfile = fopen(saveTo,"w");
  fprintf(pfile,"Ca  lambda  IDshape IDdelta DeltaStart StepDelta  elem  dt loop  RK CPUs checkpoint VolCorr StepEdge \n"); //header                                                                                                          
  fprintf(pfile,"%f  %f  %i %i  %1.16f %1.16f %i  %f  %i %i %i %i %i %i\n",Ca, lambda, IDshape, IDdelta, DeltaStart, StepDelta, n, dt, loop, RK, CPUs, checkpoint, VolCorr, IteGetShape); //options data                                              
  fclose(pfile);
  ////////////////////////////// END UPLOADING PREVIOUS DATA  /////////////////////////

  //here starts the simulation
  /////////////////////////////   RUN!!!!!!!!!!!!!!!!!!!  /////////////////////////////
  //set checkpoint for data saving
  solveDrop.SetCheckpoint(checkpoint);

  //set how often global remesh is done
  solveDrop.SetDoRemeshDistribution(DoRemesh);

  //set number of therads for parallel computation
  solveDrop.SetThreads(CPUs);

  //set time step and number of loops
  solveDrop.SetLoopDt(loop, dt);

  //choose output directory
  solveDrop.SetResdir(respath);

  //choose if activating volume correction
  solveDrop.SetVolumeCorr(VolCorr);

  if (RK==1){
  	//solve the system and advance with RK1
  	solveDrop.SolveRK1();
  	} else if (RK==2){
  	//solve the system and advance with RK2
  	solveDrop.SolveRK2();
  }
  /////////////////////////////  END RUN  ///////////////////////////////

  //determine time of the simulation 
  time_t t1;

  time(&t1);  /* get current time; same as: now = time(NULL)  */

  double seconds = difftime(t1,t0);

  printf ("Time of the simulation is %.f seconds \n", seconds);

// End of program
    std::cout << "Finished and Exiting\n"; 
    return 0;
  
}
