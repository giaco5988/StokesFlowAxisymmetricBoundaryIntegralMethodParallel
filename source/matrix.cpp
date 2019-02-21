//
//  matrix.cpp
//
//  Created by Giacomo Gallino on 26.01.16.
//

//static char help[] = "Basic vector routines.\n\n";

#include "matrix.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include "combase.h"
#include <stdlib.h>
#include <fstream>
#include <string>
#include <assert.h>
#include "remesh.h"
//#include <petscksp.h>

//this in the case openMP is not supported on the architecture
#if defined (_OPENMP)
    #include <omp.h>
#endif

using namespace std;

//constructor of the class
Matrix::Matrix (GeometryDrop *drop, Remesh *MeshDrop, double visc, double capillary, int BC){

  //for private variables
  Ca = capillary; //capillary number
  CaExtensional = 0;  //capillary number for extensional flow
  lambda = visc; //viscosity ratio

  //type of topology
  //0 is only droplet
  type = 0;

  //display topology type
  if (type==0){
  printf("Topology is initialized as one drop \n");
  }
  
  //decide BC kind
  //0 is only surface tension
  //1 is also gravity
  //2 set Ca=1 by definition
  BCkind = BC;

  //dispaly BC kind
  if(BCkind==0){
    printf("BC are calssic surface tension \n");
  }
  else if(BCkind==1)
  {
    printf("BC are like Leal rising droplet \n");
  }
  else if(BCkind==2){
    printf("BC are like Stone extensional flow \n");
    CaExtensional = Ca;
    Ca = 1;
  }
  
  nG = drop->getElem(); //drop number elements
  DropHere = drop; //drop object as global variable
  MeshHere = MeshDrop; // mesh object as global variable
  
  //allocate momory for matrix (linear elements)
  A = new double [4*(nG+1)*(nG+1)];

  //allocate memory for rhs
  RHS = new double [2*(nG+1)];

}

//set residuals for breaking
void Matrix::SetResiduals(double res){

  ResidualsBreak = res;

}

//set the number of threads
void Matrix::SetThreads(int i){

  // if parellel is on, set the number of threads
  #if defined (_OPENMP)
        omp_set_num_threads(i);
	printf("Number of threads is set to %i \n",i);
  #endif

}

//set the directory where to save the files
void Matrix::SetResdir(char* res){

  resdir = res;

  printf("Results are saved to %s \n",res);

}

//get topology type
int Matrix::GetType(){

  int n = type;
  return n;
  
}

//build matrix
void Matrix::Build(){

  //double Dx, Dy, xMid, yMid, dl;
  //double distXA, distYA, distA, distXB, distYB, distB;
  double f1x, f1y, f2x, f2y;
  double K1r, K1l, K2r, K2l, Kr, Kl, NxR, NyR, NxL, NyL;
  double f1capillary, f2capillary, f1buoyancy, f2buoyancy, f1, f2;

  //determine the type of topology, here is only one droplet
  if (type==0){

    //initialize variables
    int elem = nG; //number of elements
    int nodes = elem+1; //number of nodes
    //elements extremeties
    double *x = DropHere->GetXCoord();
    double *y = DropHere->GetYCoord();
    //nodes coordinates elements extremeties (here is the same)
    double *x0 = DropHere->GetXCoord();
    double *y0 = DropHere->GetYCoord();
    double xcm, h0, h1;
    double GXXa, GXYa, GYXa, GYYa, GXXb, GXYb, GYXb, GYYb;
    double axHere, bxHere, cxHere, dxHere, ayHere, byHere, cyHere, dyHere;

    //memory allocation for spline coefficient
    double *ax, *bx, *cx, *dx, *ay, *by, *cy, *dy;//, *x, *y, *x0, *y0;
    ax = new double [elem];
    bx = new double [elem];
    cx = new double [elem];
    dx = new double [elem];
    ay = new double [elem];
    by = new double [elem];
    cy = new double [elem];
    dy = new double [elem];
    
    //compute center of mass of the droplet
    xcm = CenterMassAxis(x,y,elem);

    //compute droplet curvature apprximating interface with a spline
    SplineSymmetric(x ,y, ax, bx, cx, dx, ay, by, cy, dy, nG);

    //loop over elements
    for (int j=0; j<elem; j++){

      //splines coeff at this loop
      axHere = ax[j]; bxHere = bx[j]; cxHere = cx[j]; dxHere = dx[j];
      ayHere = ay[j]; byHere = by[j]; cyHere = cy[j]; dyHere = dy[j];

      //in plane curvature on first node of the element
      K1r = CurvSplineR(bxHere, cxHere, byHere, cyHere);
      //in plane curvature on second node of the element
      K1l = CurvSplineL(bxHere, cxHere, dxHere, byHere, cyHere, dyHere);

      //compute normals to the node
      NxR = NormalVectorNxR(bxHere, byHere);
      NyR = NormalVectorNyR(bxHere, byHere);
      NxL = NormalVectorNxL(bxHere, cxHere, dxHere, byHere, cyHere, dyHere);
      NyL = NormalVectorNyL(bxHere, cxHere, dxHere, byHere, cyHere, dyHere);

      //for integration on splines
      h0 = sqrt(bxHere*bxHere+byHere*byHere);
      h1 =  sqrt((bxHere+2*cxHere+3*dxHere)*(bxHere+2*cxHere+3*dxHere)+(byHere+2*cyHere+3*dyHere)*(byHere+2*cyHere+3*dyHere));
      
      //azimuthal component of the curvature
      K2r = NyR/y[j]; K2l = NyL/y[j+1];
      //exception if it is first or last element
      if (j==0){
	K2r = K1r;
      }
      if (j==elem-1){
	K2l = K1l;
      }

      //total curvature
      Kr = K1r+K2r;
      Kl = K1l+K2l;

      //stresses on first element (Leal adimensionalization)
      //first node
      f1capillary = Kr/Ca;
      f1buoyancy = (x[j]-xcm)*3.0*(1+1.5*lambda)/(1+lambda)*(BCkind==1);
      f1 = f1capillary+f1buoyancy;
      f1x = f1*NxR;
      f1y = f1*NyR;
      //second node
      f2capillary = Kl/Ca;
      f2buoyancy = (x[j+1]-xcm)*3.0*(1+1.5*lambda)/(1+lambda)*(BCkind==1);
      f2 = f2capillary+f2buoyancy;
      f2x = f2*NxL;
      f2y = f2*NyL;

      //build matrix for equal viscosity, terms on the diagonal
      A[(2*(2*(elem+1)))*j+2*j] = -4*M_PI*(1+lambda);
      A[(2*(2*(elem+1)))*j+2*j+2*(elem+1)+1] = -4*M_PI*(1+lambda);

      //loop over nodes
      for (int i=0; i<nodes; i++){

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
	computeGTlinearSPlines1layer(x0[i],y0[i],h0,h1,axHere,bxHere,cxHere,dxHere,ayHere,byHere,cyHere,dyHere,GXXa,GXYa,GYXa,GYYa,GXXb,GXYb,GYXb,GYYb);
	
        //right hand side
	RHS[2*i] += GXXa*f1x + GXYa*f1y + GXXb*f2x + GXYb*f2y;
	RHS[2*i+1] += GYXa*f1x + GYYa*f1y + GYXb*f2x + GYYb*f2y;

      }

    }

    //fill last slot, not taken into account in the loop
    A[(2*(2*(nG+1)))*nG+2*nG] = -4*M_PI*(1+lambda);
    A[(2*(2*(nG+1)))*nG+2*nG+2*(nG+1)+1] = -4*M_PI*(1+lambda);

    //update curvature and normal vectors geometry object
    DropHere->UpdateNormal(ax, bx , cx, dx, ay, by, cy, dy, nG);
    DropHere->UpdateCurv(ax, bx , cx, dx, ay, by, cy, dy, nG);
    
    //evenually add extensional flow
    if (CaExtensional!=0) {

      double G = CaExtensional/2;
      double u, v;
      
      for (int i=0; i<nodes; i++){

	u = 2*G*x0[i];
        v = -G*y0[i];

	RHS[2*i] -= 8*M_PI*u;
	RHS[2*i+1] -= 8*M_PI*v;

      }

    }

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

//solve the matrix system
void Matrix::Solve(int ite){
  
  //tag breakOn
  DropHere->BreakOn(resdir, ite);

  //for (int i=0; i<nG+1;i++){

  //printf("%f %f %i %i\n",b[i],b[2*i+1],i,ite);

  //}

  //solve
  dgesv(A, RHS, 2*(nG+1));

  //tag breakOff
  DropHere->BreakOff(resdir, ite);
  
}

//set volume correction
void Matrix::SetVolumeCorr(int i){

  if (i==1){
    printf("Volume correction is activated \n");
  }

  VolCorr = i;

}

void Matrix::SaveResultsData (char *resu, int ite) {

  char saveTo[246];
  sprintf(saveTo,"%s/result%i.txt",resu,ite);

  FILE* pfile;
  //save results data                                                                                                                                                            
  pfile = fopen(saveTo,"w");
  fprintf(pfile,"u  v \n"); //header                                                                                                                                      
  for (int i=0; i<nG+1; i++){

    fprintf(pfile,"%1.16f  %1.16f\n",RHS[2*i],RHS[2*i+1]);

  }
  fclose(pfile);

}


//compute solution and advance the interface with RK1
void Matrix::SolveRK1(){

  //advance the solution with Runge Kutta order 1
  if (type==0){

    double x, y, u, v;
    int count = 0;
    int CountRemesh = 0;

    //choose mesh stabilization for droplet and set mesh parameters
    //Remesh MeshDrop(DropHere, 10);
    //MeshDrop.SetAdaptDrop(1.5);
    //MeshDrop.SetRK(1); //I'm using RK1
    
    for (int k=0; k<loop; k++){

      //volume correction
      DropHere->VolumeCorrection(VolCorr);

      //save drop data
      if (k==count||k==loop-1) {
        printf("Save data \n");
        DropHere->SaveDropData(resdir, count/checkpoint);
        count += checkpoint;
      }

      //build matrix
      Build();

      //solve linear system
      Solve(k/checkpoint);

      //force symmetry condition
      RHS[1] = 0;
      RHS[2*(nG+1)-1] = 0;

      //save drop geometry and results data                                                                                                                 
      if (k==count||k==loop-1) {

	SaveResultsData(resdir, count/checkpoint);
      }

      //print current status of the simulation
      printf("Loop %i of %i RK1 \n",k+1,loop);

      //execute remesh
      MeshHere->ExecuteRemesh(DropHere, RHS);

      if (MeshHere->GetMeshKind()==0){
	//update interface position when using mesh stabilization
	for (int i=0; i<nG+1; i++){

	  //interface coordinates
	  x = DropHere->GetOneXCoord(i);
	  y = DropHere->GetOneYCoord(i);

	  //get velocity from mesh stab
	  u = MeshHere->GetOneVelU(i);
	  v = MeshHere->GetOneVelV(i);
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

	//interface coordinates
	x = DropHere->GetOneXCoord(i);
	y = DropHere->GetOneYCoord(i);
	  
	//adavance interface
	x += RHS[2*i]*dt;
	y += RHS[2*i+1]*dt;

	//update coordinates in droplet object
	DropHere->SetOneCoord(x, y, i);

	//clean rhs
	RHS[2*i] = 0;
	RHS[2*i+1] = 0;
	}
      }

      //DropHere->CheckBreakDomain();

    }

  }
  
}

void Matrix::SetDoRemeshDistribution(int i){

  DoRemeshDistribution = i;

}

void Matrix::SetCheckpoint(int i){

  checkpoint = i;

}

//set how much remesh to do when I start the simulation
void Matrix::SetFirstRemeshLoop(int n){

  FirstRemeshLoop = n;

}

//compute solution and advance the interface with RK2
void Matrix::SolveRK2(){

  //advance the solution with Runge Kutta order 2
  if (type==0){

    double x1, y1, x2, y2, *nxTemp, *nyTemp, *k1Temp, *k2Temp, *Dellipse;
    int count = 0;
    int CountRemesh = 0;

    //initialize vectors
    k1Temp = new double [nG+1];
    k2Temp = new double [nG+1];
    nxTemp = new double [nG+1];
    nyTemp = new double [nG+1];
    Dellipse = new double [loop];

    //compute how often to check for residuals
    double Tscale = 1/Ca;
    if (BCkind == 2){
      Tscale = 1/CaExtensional;
    }
    int WhenResiduals = Tscale/dt;
    double Residuals = 100;
    //printf("%i \n",WhenResiduals);
    
    for (int k=0; k<loop; k++){

      //compute ellipsicity for residuals
      Dellipse[k] = DropHere->ComputeEllipse();

      //FIRST RK LOOP

      //remesh using distribution
      if (k==CountRemesh && MeshHere->GetMeshKind()==1){

	int numLoop = 1;
	if (k==0){

	  numLoop = FirstRemeshLoop;

	}

	for (int k=0; k<numLoop; k++){
	  //remesh using distribution                                                                                                                                                                                                               
	  double *xStart, *yStart;
	  xStart = new double [nG+1];
	  yStart = new double [nG+1];

	  MeshHere->RemeshDistribution(DropHere, xStart, yStart);

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

      //save drop data
      if (k==count||k==loop-1) {
	//printf("Ciao \n");
	printf("Save data \n");
	//printf("Ciao \n");
        DropHere->SaveDropData(resdir, count/checkpoint);
	count += checkpoint;
      }

      //build matrix
      Build();

      //solve linear system
      Solve(k/checkpoint);

      //force symmetry condition
      RHS[1] = 0;
      RHS[2*(nG+1)-1] = 0;

      //print current status of the simulation
      printf("Loop %i of %i RK2 .",k+1,loop);

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
      Build();

      //solve linear system
      Solve(k/checkpoint);

      //force symmetry condition
      RHS[1] = 0;
      RHS[2*(nG+1)-1] = 0;

      //save drop geometry and results data                                                                                                                       
      if (k==count||k==loop-1) {

        SaveResultsData(resdir, count/checkpoint);
      }

      //print current status of the simulation
      printf(". Residuals=%f \n",Residuals);

      //execute remesh
      MeshHere->ExecuteRemesh(DropHere, RHS);

      if (MeshHere->GetMeshKind()==0){
	//update interface position when using mesh stabilization
	for (int i=0; i<nG+1; i++){

	  x2 = DropHere->GetOneOldXCoord(i);
	  y2 = DropHere->GetOneOldYCoord(i);

	  //get velocity from mesh stab
	  RHS[2*i] = MeshHere->GetOneVelU(i);
	  RHS[2*i+1] = MeshHere->GetOneVelV(i);
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

      //break if a steady shape has been obtained
      if (k>WhenResiduals) {
	Residuals = sqrt(pow(Dellipse[k]-Dellipse[k-WhenResiduals],2));
	//printf("%f %f %f \n", Dellipse[k], Dellipse[k-WhenResiduals],Residuals);
	if (Residuals<ResidualsBreak & Dellipse[k]<0.3) {
	  printf("Steady shape has been reached \n");
	  break;
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
    delete Dellipse;

  }
  
}

void Matrix::SetLoopDt(int nLOOP, double deltaT){

  loop = nLOOP;
  dt = deltaT;

}
