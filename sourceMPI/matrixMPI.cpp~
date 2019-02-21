//
//  matrix.cpp
//
//  Created by Giacomo Gallino on 26.01.16.
//

//static char help[] = "Basic vector routines.\n\n";

#include "matrixMPI.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include "combaseMPI.h"
#include <stdlib.h>
#include <fstream>
#include <string>
#include <assert.h>
#include "remeshMPI.h"
#include <petscksp.h>
//#include <petsctime.h>

using namespace std;

//constructor of the class
MatrixMPI::MatrixMPI (GeometryMPI *drop, double visc, double capillary, int BC){

  //get rank of this processor
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  //for private variables
  Ca = capillary; //capillary number
  CaExtensional = 0;  //capillary number for extensional flow
  lambda = visc; //viscosity ratio

  //type of topology
  //0 is only droplet
  type = 0;

  //display topology type
  if (type==0 && rank==0){
    printf("Topology is initialized as one drop \n");
  }
  
  //decide BC kind
  //0 is only surface tension
  //1 is also gravity
  //2 set Ca=1 by definition
  BCkind = BC;

  //dispaly BC kind
  if(BCkind==0 && rank==0){
    printf("BC are calssic surface tension \n");
  }
  else if(BCkind==1 && rank==0)
  {
    printf("BC are like Leal rising droplet \n");
  }
  else if(BCkind==2){
    if (rank==0) {
      printf("BC are like Stone extensional flow \n");
    }
    CaExtensional = Ca;
    Ca = 1;
  }
  
  nG = drop->getElem(); //drop number elements
  DropHere = drop; //drop object as local variable
  
  //create petsc complex objects like vectors matrices and oinear solver
  NODES = nG+1;

  //create vectors                                                                                                                  
  VecCreate(PETSC_COMM_WORLD,&sol);
  VecSetSizes(sol,PETSC_DECIDE,2*NODES);
  VecSetFromOptions(sol);
  VecDuplicate(sol,&b);
  VecDuplicate(sol,&rhs);
  
  //create scatter context
  //VecScatterCreate(sol,from,transfer,to,&scatter);
  VecScatterCreateToAll(sol,&scatter,&transfer);

  //create matrices
  //left hand side matrix
  //MatCreate(PETSC_COMM_WORLD,&A);
  MatCreateDense(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,2*NODES,2*NODES,NULL,&A);
  //MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,2*NODES,2*NODES);
  MatSetFromOptions(A);
  MatSetUp(A);
  
  //rhs matrix
  //MatCreate(PETSC_COMM_WORLD,&U);
  MatCreateDense(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,2*NODES,2*NODES,NULL,&U);
  //MatSetSizes(U,PETSC_DECIDE,PETSC_DECIDE,2*NODES,2*NODES);
  MatSetFromOptions(U);
  MatSetUp(U);

  //fake vector just to get ownership                                                                                                                                                                                                       
  VecCreate(PETSC_COMM_WORLD,&fake);
  VecSetSizes(fake,PETSC_DECIDE,NODES);
  VecSetFromOptions(fake);
  //VecGetOwnershipRange(fake,&Istart,&Iend);

  //linear solver
  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetOperators(ksp,A,A);
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCJACOBI);
  //KSPSetFromOptions(ksp);

  //allocate memory for rhs
  //RHS = new PetscScalar [2*NODES];
  PetscMalloc1(2*NODES,&RHS);

}

//set the directory where to save the files
void MatrixMPI::SetResdir(char* res){

  resdir = res;

  //get rank of this processor
  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  if (rank==0) {

    printf("Results are saved to %s \n",res);

  }

}

//get topology type
int MatrixMPI::GetType(){

  int n = type;
  return n;
  
}

//build matrix and solve it using parallel computing with PETSc library
void MatrixMPI::BuildAndSolvePETSc(){

  clock_t t0 = clock();

  //define variables
  PetscScalar fx, fy;
  PetscScalar K1, K2, K, Nx, Ny;
  PetscScalar fcapillary, fbuoyancy, f;

  //define PETSc variables
  PetscInt nodes, Istart, Iend, J, I;

  //set the values of the vector and matrices to 0
  VecSet(sol, 0);
  VecSet(b, 0);
  VecSet(rhs, 0);
  MatZeroEntries(A);
  MatZeroEntries(U);

  //determine the type of topology, here is only one droplet
  if (type==0){

    //NOT CLEAR WHAT IS DOING
    //PetscOptionsGetInt(NULL,"-NodesPetsc",&NodesPetsc,NULL);

    //petsc variable for vector and matrix size
    nodes = nG+1;

    //create vectors                                                                                                                             
    VecGetOwnershipRange(fake,&Istart,&Iend);

    //initialize variables
    PetscInt elem = nG; //number of elements
    //elements extremeties
    PetscScalar *x = DropHere->GetXCoord();
    PetscScalar *y = DropHere->GetYCoord();
    //nodes coordinates elements extremeties (here is the same)
    PetscScalar *x0 = DropHere->GetXCoord();
    PetscScalar *y0 = DropHere->GetYCoord();
    PetscScalar xcm, h0, h1;
    PetscScalar GXXa, GXYa, GYXa, GYYa, GXXb, GXYb, GYXb, GYYb;
    PetscScalar axHere, bxHere, cxHere, dxHere, ayHere, byHere, cyHere, dyHere;

    //memory allocation for spline coefficient
    PetscScalar *ax, *bx, *cx, *dx, *ay, *by, *cy, *dy;
    ax = new PetscScalar [elem];
    bx = new PetscScalar [elem];
    cx = new PetscScalar [elem];
    dx = new PetscScalar [elem];
    ay = new PetscScalar [elem];
    by = new PetscScalar [elem];
    cy = new PetscScalar [elem];
    dy = new PetscScalar [elem];
    
    //compute center of mass of the droplet
    xcm = CenterMassAxis(x,y,elem);

    //compute droplet curvature apprximating interface with a spline
    SplineSymmetric(x, y, ax, bx, cx, dx, ay, by, cy, dy, elem);

    //loop over elements
    for (J=Istart; J<Iend; J++){

      //from second to before last node
      if (J>0 && J<elem) {

      //splines coeff at this loop
      axHere = ax[J]; bxHere = bx[J]; cxHere = cx[J]; dxHere = dx[J];
      ayHere = ay[J]; byHere = by[J]; cyHere = cy[J]; dyHere = dy[J];

      //printf("%f %i \n",axHere,J);

      //in plane curvature on first node of the element
      K1 = CurvSplineR(bxHere, cxHere, byHere, cyHere);

      //compute normals to the node
      Nx = NormalVectorNxR(bxHere, byHere);
      Ny = NormalVectorNyR(bxHere, byHere);

      //for integration on splines
      h0 = sqrt(bxHere*bxHere+byHere*byHere);
      h1 =  sqrt((bxHere+2*cxHere+3*dxHere)*(bxHere+2*cxHere+3*dxHere)+(byHere+2*cyHere+3*dyHere)*(byHere+2*cyHere+3*dyHere));
      
      //azimuthal component of the curvature
      K2 = Ny/y[J];

      //total curvature
      K = K1+K2;

      //stresses on first element (Leal adimensionalization)
      //first node
      fcapillary = K/Ca;
      fbuoyancy = (x[J]-xcm)*3.0*(1+1.5*lambda)/(1+lambda)*(BCkind==1);
      f = fcapillary+fbuoyancy;
      fx = f*Nx;
      fy = f*Ny;

      //store stresses in vector
      PetscInt cols[2];
      cols[0] = 2*J; cols[1] = 2*J+1;
      PetscScalar valsVec[2];
      valsVec[0] = fx; valsVec[1] = fy;
      VecSetValues(b,2,cols,valsVec,ADD_VALUES);

      //build matrix for equal viscosity, terms on the diagonal
      PetscInt rows[2]; PetscScalar valsMat[4];
      //coordinates
      rows[0] = 2*J; rows[1] = 2*J+1;
      cols[0] = 2*J; cols[1] = 2*J+1;
      //values
      valsMat[0] = -4*M_PI*(1+lambda); valsMat[1] = 0;
      valsMat[2] = 0; valsMat[3] = -4*M_PI*(1+lambda);
      MatSetValues(A, 2, rows, 2, cols, valsMat, ADD_VALUES);

      //loop over nodes
      for (I=0; I<nodes; I++){

	//set to zero for new integration
	GXXa = 0.0;
	GXYa = 0.0;
	GYXa = 0.0;
	GYYa = 0.0;
        GXXb = 0.0;
	GXYb = 0.0;
	GYXb = 0.0;
	GYYb = 0.0;

	//compute and integrate green's functions single layer
	computeGTlinearSPlines1layer(x0[I],y0[I],h0,h1,axHere,bxHere,cxHere,dxHere,ayHere,byHere,cyHere,dyHere,GXXa,GXYa,GYXa,GYYa,GXXb,GXYb,GYXb,GYYb);
	
        //right hand side
	//PhiA contribution
	//coordinates
	rows[0] = 2*I; rows[1] = 2*I+1;
	cols[0] = 2*J; cols[1] = 2*J+1;
	//values
	valsMat[0] = GXXa; valsMat[1] = GXYa;
	valsMat[2] = GYXa; valsMat[3] = GYYa;
	MatSetValues(U, 2, rows, 2, cols, valsMat, ADD_VALUES);

	//PhiB contribution
	//coordinates
        rows[0] = 2*I; rows[1] = 2*I+1;
        cols[0] = 2*J+2; cols[1] = 2*J+3;
	//values
        valsMat[0] = GXXb; valsMat[1] = GXYb;
        valsMat[2] = GYXb; valsMat[3] = GYYb;
	MatSetValues(U, 2, rows, 2, cols, valsMat, ADD_VALUES);

      }

      //if it is first node
      } else if (J==0) {

	//splines coeff at this loop                                                                                                                                                                                               
	axHere = ax[J]; bxHere = bx[J]; cxHere = cx[J]; dxHere = dx[J];
	ayHere = ay[J]; byHere = by[J]; cyHere = cy[J]; dyHere = dy[J];

	//in plane curvature on first node of the element                                                                                                                                                                                       
	K1 = CurvSplineR(bxHere, cxHere, byHere, cyHere);

	//compute normals to the node                                                                                                                                                                                                           
	Nx = NormalVectorNxR(bxHere, byHere);
	Ny = NormalVectorNyR(bxHere, byHere);

	//for integration on splines                                                                                                                                                                                                            
	h0 = sqrt(bxHere*bxHere+byHere*byHere);
	h1 =  sqrt((bxHere+2*cxHere+3*dxHere)*(bxHere+2*cxHere+3*dxHere)+(byHere+2*cyHere+3*dyHere)*(byHere+2*cyHere+3*dyHere));

	//total curvature                                                                                                                                                                                                           
	K = 2*K1;

	//stresses on first element (Leal adimensionalization)                                                                                                                                                                           
	//first node                                                                                                                                                                                                                         
	fcapillary = K/Ca;
	fbuoyancy = (x[J]-xcm)*3.0*(1+1.5*lambda)/(1+lambda)*(BCkind==1);
	f = fcapillary+fbuoyancy;
	fx = f*Nx;
	fy = f*Ny;

	//store stresses in vector                                                                                                                                                                                                    
	PetscInt cols[2];
	cols[0] = 2*J; cols[1] = 2*J+1;
	PetscScalar valsVec[2];
	valsVec[0] = fx; valsVec[1] = fy;
	VecSetValues(b,2,cols,valsVec,ADD_VALUES);

	//build matrix for equal viscosity, terms on the diagonal                                                                                                                                                               
	PetscInt rows[2]; PetscScalar valsMat[4];
	//coordinates                                                                                                                                                                                                                 
	rows[0] = 2*J; rows[1] = 2*J+1;
	cols[0] = 2*J; cols[1] = 2*J+1;
	//values                                                                                                                                                                                                              
	valsMat[0] = -4*M_PI*(1+lambda); valsMat[1] = 0;
	valsMat[2] = 0; valsMat[3] = -4*M_PI*(1+lambda);
	MatSetValues(A, 2, rows, 2, cols, valsMat, ADD_VALUES);

	//loop over nodes                                                                                                                                                                                                 
	for (I=0; I<nodes; I++){

	  //set to zero for new integration                                                                                                                                                                               
	  GXXa = 0.0;
	  GXYa = 0.0;
	  GYXa = 0.0;
	  GYYa = 0.0;
	  GXXb = 0.0;
	  GXYb = 0.0;
	  GYXb = 0.0;
	  GYYb = 0.0;

	  //compute and integrate green's functions single layer                                                                                                                                                                  
	  computeGTlinearSPlines1layer(x0[I],y0[I],h0,h1,axHere,bxHere,cxHere,dxHere,ayHere,byHere,cyHere,dyHere,GXXa,GXYa,GYXa,GYYa,GXXb,GXYb,GYXb,GYYb);
	 
	  //right hand side                                                                                                                                                                                                           
	  //PhiA contribution                                                                                                                                                                                                        
	  //coordinates                                                                                                                                                                                                             
	  rows[0] = 2*I; rows[1] = 2*I+1;
	  cols[0] = 2*J; cols[1] = 2*J+1;
	  //values                                                                                                                                                                                                                 
	  valsMat[0] = GXXa; valsMat[1] = GXYa;
	  valsMat[2] = GYXa; valsMat[3] = GYYa;
	  MatSetValues(U, 2, rows, 2, cols, valsMat, ADD_VALUES);
	  
	  //PhiB contribution                                                                                                                                                                                                
	  //coordinates                                                                                                                                                                                                   
	  rows[0] = 2*I; rows[1] = 2*I+1;
	  cols[0] = 2*J+2; cols[1] = 2*J+3;
	  //values                                                                                                                                                                                                             
	  valsMat[0] = GXXb; valsMat[1] = GXYb;
	  valsMat[2] = GYXb; valsMat[3] = GYYb;
	  MatSetValues(U, 2, rows, 2, cols, valsMat, ADD_VALUES);

	}

      //if it is last node
      } else if (J==elem) {

	//splines coeff at this loop                                                                                                               
	axHere = ax[J-1]; bxHere = bx[J-1]; cxHere = cx[J-1]; dxHere = dx[J-1];
        ayHere = ay[J-1]; byHere = by[J-1]; cyHere = cy[J-1]; dyHere = dy[J-1];

	//LAST NODE                                                                                                                                 
	K1 = CurvSplineL(bxHere, cxHere, dxHere, byHere, cyHere, dyHere);                                                                                                                                                             
	//compute normals to the node                                                                                                                      
	Nx = NormalVectorNxL(bxHere, cxHere, dxHere, byHere, cyHere, dyHere);                                                                                         Ny = NormalVectorNyL(bxHere, cxHere, dxHere, byHere, cyHere, dyHere);                                                                                                                                                         
	//total curvature                                                                                                                          
	K = 2*K1;                                                                                                                                                                                                                                 
	//stresses on first element (Leal adimensionalization)                                                                                        
	fcapillary = K/Ca;                                                                                                                                            fbuoyancy = (x[elem+1]-xcm)*3.0*(1+1.5*lambda)/(1+lambda)*(BCkind==1);                                                                                        f = fcapillary+fbuoyancy;                                                                                                                                     fx = f*Nx;                                                                                                                                                    fy = f*Ny;
	
	//store stresses in vector                                                                                                                                                                                    
	PetscInt cols[2];
	cols[0] = 2*J; cols[1] = 2*J+1;
	PetscScalar valsVec[2];
	valsVec[0] = fx; valsVec[1] = fy;
	VecSetValues(b,2,cols,valsVec,ADD_VALUES);

	//diagonal terms in matrix A
	PetscInt rows[2]; PetscScalar valsMat[4];
        //coordinates                                                                                                                                                                                                 
	rows[0] = 2*J; rows[1] = 2*J+1;
        cols[0] = 2*J; cols[1] = 2*J+1;
        //values                                                                                                                                                                                         
	valsMat[0] = -4*M_PI*(1+lambda); valsMat[1] = 0;
        valsMat[2] = 0; valsMat[3] = -4*M_PI*(1+lambda);
        MatSetValues(A, 2, rows, 2, cols, valsMat, ADD_VALUES);

	}

    }

    //assembly
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    MatAssemblyBegin(U,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(U,MAT_FINAL_ASSEMBLY);

    //compute right hand side rhs = U*b
    MatMult(U, b, rhs);

    //check matrices and vectors
    //VecView(b,PETSC_VIEWER_STDOUT_WORLD);

    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    //MatView(A,PETSC_VIEWER_STDOUT_WORLD);
                                                                                                                            
    //MatView(U,PETSC_VIEWER_STDOUT_WORLD);

    //update curvature and normal vectors geometry object
    DropHere->UpdateNormal(ax, bx , cx, dx, ay, by, cy, dy, nG);
    DropHere->UpdateCurv(ax, bx , cx, dx, ay, by, cy, dy, nG);
    
    //evenually add extensional flow
    if (CaExtensional!=0) {

      PetscScalar G = CaExtensional/2;
      
      for (PetscInt i=Istart; i<Iend; i++){

	PetscInt cols[2];      PetscScalar valsVec[2];
	cols[0] = 2*i;         cols[1] = 2*i+1;
	valsVec[0] = -16*M_PI*G*x0[i];   valsVec[1] = 8*M_PI*G*y0[i];

	VecSetValues(rhs,2,cols,valsVec,ADD_VALUES);

      }

    }

    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    //determine time of the simulation 
    clock_t t1 = clock()-t0;
    double secondsMatrix = ((double)t1) / CLOCKS_PER_SEC * 1000;
    if (rank==0) {
       char saveTo[500];
       sprintf(saveTo,"%s/BuildMatrix.txt",resdir);
       FILE* pfile;
       //save geometrical data
       pfile = fopen(saveTo,"w");
       fprintf(pfile,"%f \n",secondsMatrix); //header
       fclose(pfile);
    }
    printf ("Only build matrix takes %.f milliseconds \n", secondsMatrix);

    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);
    //VecView(rhs,PETSC_VIEWER_STDOUT_WORLD);

    //solve the linear system Ax = rhs <--> Ax = Ub
    KSPSolve(ksp,rhs,sol);

    //check is saving time is correct
    //for (int i=0; i<10e7; i++){
      //double r = 1;
    //}

    VecAssemblyBegin(sol);
    VecAssemblyEnd(sol);
    //VecView(sol,PETSC_VIEWER_STDOUT_WORLD);

     //determine time of the simulation 
    clock_t t2 = clock() - t1 - t0;
    double secondsSolve = ((double)t2) / CLOCKS_PER_SEC * 1000;
     if (rank==0) {
       char saveTo[500];
       sprintf(saveTo,"%s/SolveMatrix.txt",resdir);
       FILE* pfile;
       //save geometrical data
       pfile = fopen(saveTo,"w");
       fprintf(pfile,"%f \n",secondsSolve); //header
       fclose(pfile);
    }
     printf ("Only solving linear system takes %.f milliseconds \n", secondsSolve);
     

    //push solution to global variable
    VecGetOwnershipRange(sol,&Istart,&Iend);

    VecScatterBegin(scatter,sol,transfer,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter,sol,transfer,INSERT_VALUES,SCATTER_FORWARD);
    //VecView(transfer,PETSC_VIEWER_STDOUT_WORLD);
    VecGetArray(transfer, &RHS);

    //free memory
    delete[] ax;
    delete[] bx;
    delete[] cx;
    delete[] dx;
    delete[] ay;
    delete[] by;
    delete[] cy;
    delete[] dy;

  }

}

//set volume correction
void MatrixMPI::SetVolumeCorr(int i){

  if (i==1){
    printf("Volume correction is activated \n");
  }

  VolCorr = i;

}

//compute solution and advance the interface with RK1
void MatrixMPI::SolveRK1(){

  //advance the solution with Runge Kutta order 1
  if (type==0){

    double x, y, u, v;
    int count = 0;
    //int CountRemesh = 0;
    PetscMPIInt rank;
    PetscInt Istart, Iend;

    //choose mesh stabilization for droplet and set mesh parameters
    RemeshMPI MeshDrop(DropHere);
    MeshDrop.SetAdaptDrop(1.5);
    MeshDrop.SetRK(1); //I'm using RK1
    
    for (int k=0; k<loop; k++){

      //get rank of this processor
      //PetscMPIInt rank;
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

      //volume correction
      DropHere->VolumeCorrection(VolCorr);

      if (rank==0) {

	//save drop data
	if (k==count||k==loop-1) {
	  printf("Save data \n");
	  DropHere->SaveDropData(resdir, count/checkpoint);
	  count += checkpoint;
	}
      }

      //build ans solve with PETSc
      BuildAndSolvePETSc();

      //force symmetry condition
      RHS[1] = 0;
      RHS[2*(nG+1)-1] = 0;

      //for (int i=0; i<2*nG+2; i++) {
      //printf("%f %i \n",RHS[i],i);
      //}

      if (rank==0){
	//print current status of the simulation
	printf("Loop %i of %i RK1 \n",k+1,loop);
      }

      //execute remesh
      //MeshDrop.ExecuteRemesh(DropHere, RHS);
      if (MeshKind==0){
	MeshDrop.PassiveMeshStabilization(DropHere,RHS);
      }

      //VecGetOwnershipRange(fake,&Istart,&Iend);

      if (MeshKind==0){
	//update interface position when using mesh stabilization
	for (int i=0; i<nG+1; i++){
	//for (PetscInt i=Istart; i<Iend; i++){

	  //interface coordinates
	  x = DropHere->GetOneXCoord(i);
	  y = DropHere->GetOneYCoord(i);

	  //get velocity from mesh stab
	  u = MeshDrop.GetOneVelU(i);
	  v = MeshDrop.GetOneVelV(i);
	  x += u*dt;
	  y += v*dt;

	  //update coordinates in droplet object
	  DropHere->SetOneCoord(x, y, i);

	  //printf("%f %f \n",b[2*i],b[2*i+1]);

	  //clean rhs
	  RHS[2*i] = 0;
	  RHS[2*i+1] = 0;
	  
	}
      } else {
	//update interface position
	for (int i=0; i<nG+1; i++){
	  //for (PetscInt i=Istart; i<Iend; i++){

	//interface coordinates
	x = DropHere->GetOneXCoord(i);
	y = DropHere->GetOneYCoord(i);
	  
	//adavance interface
	x += RHS[2*i]*dt;
	y += RHS[2*i+1]*dt;

	//update coordinates in droplet object
	DropHere->SetOneCoord(x, y, i);

	//printf("%f %f %i \n",RHS[2*i],RHS[2*i+1],i);

	//clean rhs
	RHS[2*i] = 0;
	RHS[2*i+1] = 0;
	}
	
      }

      //CheclBreakDomain();

    }

  }
  
}

void MatrixMPI::SetDoRemeshDistribution(int i){

  DoRemeshDistribution = i;

}

void MatrixMPI::SetCheckpoint(int i){

  checkpoint = i;

}

//compute solution and advance the interface with RK2
void MatrixMPI::SolveRK2(){

  //advance the solution with Runge Kutta order 2
  if (type==0){

    double x1, y1, x2, y2, *nxTemp, *nyTemp, *k1Temp, *k2Temp;
    int count = 0;
    int CountRemesh = 0;

    //initialize vectors
    k1Temp = new double [nG+1];
    k2Temp = new double [nG+1];
    nxTemp = new double [nG+1];
    nyTemp = new double [nG+1];

    //choose mesh stabilization for droplet and set mesh parameters
    RemeshMPI MeshDrop(DropHere);
    MeshDrop.SetAdaptDrop(1.5);//adaptivity
    MeshDrop.SetRK(2); //I'm using RK2
    MeshDrop.SetOptionCluster(1); //1 is classic based on curvature, 2 clustering close to axis, 3 on the right part of the droplet
    MeshDrop.SetTuneDistributionAdaptivity(10.0);
    MeshDrop.SetOptionDistribution(3); //1 is on the right, 2 is even more on the right, 3 is uniform
    
    for (int k=0; k<loop; k++){

      //get rank of this processor
      PetscMPIInt rank;
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

      //FIRST RK LOOP
      //remesh using distribution
      if (k==CountRemesh && MeshKind==1){

	int numLoop = 1;
	if (k==0){

	  numLoop = 0;

	}

	for (int k=0; k<numLoop; k++){
	  //remesh using distribution                                                                                                                                                                                                               
	  double *xStart, *yStart;
	  xStart = new double [nG+1];
	  yStart = new double [nG+1];

	  MeshDrop.RemeshDistribution(DropHere, xStart, yStart);

	  //push to coordiantes to drop object                                                                                                                                                                                 
	  for (int i=0; i<nG+1; i++){

	    DropHere->SetOneCoord(xStart[i], yStart[i], i);
	    
	  }

	  //free memory                                             
	  delete[] xStart;
	  delete[] yStart;

	}

	CountRemesh += DoRemeshDistribution;

      }

      //volume correction
      DropHere->VolumeCorrection(VolCorr);

      if (rank==0) {

      //save drop data
      if (k==count||k==loop-1) {
	printf("Save data \n");
        DropHere->SaveDropData(resdir, count/checkpoint);
	count += checkpoint;
      }

      }

      //build and solve using PETSc
      BuildAndSolvePETSc();

      //force symmetry condition
      RHS[1] = 0;
      RHS[2*(nG+1)-1] = 0;

      if (rank==0) {
	//print current status of the simulation
	printf("Loop %i of %i RK2 .",k+1,loop);
      }

      //old geometrical parameters
      nxTemp = DropHere->GetNormalX();
      nyTemp = DropHere->GetNormalY();
      k1Temp = DropHere->GetK1();
      k2Temp = DropHere->GetK2();

      for (int i=0; i<nG+1; i++){

	//store data before moving the interace
	//normal vector
	DropHere->SetOldNormal(nxTemp[i], nyTemp[i], i);
	DropHere->SetOldCurv(k1Temp[i], k2Temp[i], i);

	x1 = DropHere->GetOneXCoord(i);
	y1 = DropHere->GetOneYCoord(i);

	//printf("%f %f \n",x1,y1);

	DropHere->SetOneOldCoord(x1, y1, i);
	
	//adavance interface
	x1 += RHS[2*i]*dt/2;
	y1 += RHS[2*i+1]*dt/2;

	//printf("%f %f \n",b[2*i],b[2*i+1]);

	//update coordinates in droplet object
	DropHere->SetOneCoord(x1, y1, i);
	
	//clean rhs
	RHS[2*i] = 0;
	RHS[2*i+1] = 0;
       
      }

      //THE END OF FIRST RK LOOP

      //START SECOND RK LOOP
      //build and solve with petsc
      BuildAndSolvePETSc();

      //force symmetry condition
      RHS[1] = 0;
      RHS[2*(nG+1)-1] = 0;

      if (rank==0) {
	//print current status of the simulation
	printf(".\n");
      }

      //execute passive mesh stabilization
      //MeshDrop.ExecuteRemesh(DropHere, RHS);
      if (MeshKind==0) {
	MeshDrop.PassiveMeshStabilization(DropHere, RHS);
      }

      if (MeshKind==0){
	//update interface position when using mesh stabilization
	for (int i=0; i<nG+1; i++){

	  x2 = DropHere->GetOneOldXCoord(i);
	  y2 = DropHere->GetOneOldYCoord(i);

	  //get velocity from mesh stab
	  RHS[2*i] = MeshDrop.GetOneVelU(i);
	  RHS[2*i+1] = MeshDrop.GetOneVelV(i);
	  x2 += RHS[2*i]*dt;
	  y2 += RHS[2*i+1]*dt;

	  //update coordinates in droplet object
	  DropHere->SetOneCoord(x2, y2, i);

	  //clean rhs
	  RHS[2*i] = 0;
	  RHS[2*i+1] = 0;
	  
	}
      } else {
	//update interface position
	for (int i=0; i<nG+1; i++){

	  //upload old interafce position
	  x2 = DropHere->GetOneOldXCoord(i);
	  y2 = DropHere->GetOneOldYCoord(i);
	  
	  //advance interface
	  x2 += RHS[2*i]*dt;
	  y2 += RHS[2*i+1]*dt;

	  //clean rhs
	  RHS[2*i] = 0;
	  RHS[2*i+1] = 0;

	  //update coordinates in droplet object
	  DropHere->SetOneCoord(x2, y2, i);

	}
      }

      //break if the drop retuned sperical
      if (DropHere->CheckBreakStable()) {
	break;
      }

    }

    //free memory
    delete nxTemp;
    delete nyTemp;
    delete k1Temp;
    delete k2Temp;

  }
  
}

void MatrixMPI::SetLoopDt(int nLOOP, double deltaT){

  loop = nLOOP;
  dt = deltaT;

}

//set which remesh I'm using: 0 is mesh stabilzation, 1 is remesh with distribution
void MatrixMPI::SetMeshKind(int x){

  MeshKind = x;

}
