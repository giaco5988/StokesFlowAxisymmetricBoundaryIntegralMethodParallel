//
//  geometrydrop.cpp
//
//  Created by Giacomo Gallino on 26.01.16.
//

#include "geometrydrop.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include "combase.h"
#include <stdlib.h>
#include <fstream>
#include <string>
#include <assert.h>
#include <time.h>
#include <algorithm>

using namespace std;

//constructor of the class
GeometryDrop::GeometryDrop (double x0, double y0, double D, int elem){

  //memory allocation
  xCoord = new double [elem+1];
  yCoord = new double [elem+1];
  k1 = new double [elem+1];
  k2 = new double [elem+1];
  nx = new double [elem+1];
  ny = new double [elem+1];
  OldxCoord = new double [elem+1];
  OldyCoord = new double [elem+1];
  nxOld = new double [elem+1];
  nyOld = new double [elem+1];
  k1Old = new double [elem+1];
  k2Old = new double [elem+1];

  //set starting volume of droplet
  V0drop = 4.0/3.0*M_PI;

  //in order to check running time
  time(&Tstart);

  //semi axis from D and volume of sphere with radius=1
  double b = pow((1.0-D)/(1.0+D),1.0/3);
  double a = 1.0/pow(b,2);
  
  //define number of elements in the private class
  DropElem = elem;

  //build drop, volume is the one of a sphere of radius 1
  double dtheta = M_PI/elem;
  double theta = 0.0;
  for(int i=0; i<DropElem+1; i++){
  xCoord[i] = x0 + a*cos(theta);
  yCoord[i] = y0 + b*sin(theta);

  theta = theta+dtheta;
  }

}

//overloading constructor of the class (upload data from prevoius simulation)
GeometryDrop::GeometryDrop (int IDshape, int IDdelta, double D, double StepDelta, int elem){

  //memory allocation
  xCoord = new double [elem+1];
  yCoord = new double [elem+1];
  k1 = new double [elem+1];
  k2 = new double [elem+1];
  nx = new double [elem+1];
  ny = new double [elem+1];
  OldxCoord = new double [elem+1];
  OldyCoord = new double [elem+1];
  nxOld = new double [elem+1];
  nyOld = new double [elem+1];
  k1Old = new double [elem+1];
  k2Old = new double [elem+1];

  //define number of elements in the private class
  DropElem = elem;

  //if (IDshape==1){

    //semi axis from D and volume of sphere with radius=1
    D += StepDelta*IDdelta;
    double b = pow((1.0-D)/(1.0+D),1.0/3);
    double a = 1.0/pow(b,2);

    printf("Edge state: first shape start as an ellipsoid of ellipsicity D=%1.16f \n",D);

    //build drop, volume is the one of a sphere of radius 1                                                                                                                                                                  
    double dtheta = M_PI/elem;
    double theta = 0.0;
    for(int i=0; i<DropElem+1; i++){
      xCoord[i] = a*cos(theta);
      yCoord[i] = b*sin(theta);

      //printf("%f %f %i \n",xCoord[i],yCoord[i],i);

      theta = theta+dtheta;
    }

    //}

  //set starting volume of droplet
  V0drop = 4.0/3.0*M_PI;

  //in order to check running time
  time(&Tstart);

}

//check if previous round was already fine enough
void GeometryDrop::CheckPreviousRound(int IDshape, int IDdelta, int IDround, char* resBefore, char* resAfter){

  //check if previous round was fine enough and eventually break                                                                                                                                                                          
  char LoadFrom[500];

  int BreakRound;
  sprintf(LoadFrom,"%sIDshape=%i_IDdelta=1_Round=%i_%s/round.txt",resBefore,IDshape,IDround-1,resAfter);
  fstream myfile(LoadFrom, std::ios_base::in);
  //printf("%s \n",LoadFrom);                                                                                                                                                                                                             
  myfile >> BreakRound;
  if (BreakRound==1) {
    //save round for other simulations                                                                                                                                                                                                    
    char saveTo[500];
    sprintf(saveTo,"%sIDshape=%i_IDdelta=%i_Round=%i_%s/round.txt",resBefore,IDshape,IDdelta,IDround,resAfter);                                                                                                                                                                                                                                                                 
    FILE* pfile;
    pfile = fopen(saveTo,"w");
    fprintf(pfile,"%i \n",1);
    fclose(pfile);
    printf("\n Previous round or even before was already fine enough, break! \n");
    exit(0);
  }
  myfile.close();

}

//check if delta step is under the desired one
void GeometryDrop::CheckDeltaStep(double delta1, double delta2, double StepDelta, char* resBefore, char*resAfter, int IDshape, int IDdelta, int IDround){

  if (fabs(delta2-delta1)<StepDelta) {

    //save round for other simulations                                                                                                                                                                                                      
    char saveTo[500];
    sprintf(saveTo,"%sIDshape=%i_IDdelta=%i_Round=%i_%s/round.txt",resBefore,IDshape,IDdelta,IDround,resAfter);

    //printf("%s",saveTo);                                                                                                            
                                                                                                                                                                                                       
    FILE* pfile;
    pfile = fopen(saveTo,"w");
    fprintf(pfile,"%i \n",1);
    fclose(pfile);
    printf("Desired delta step has been reached, next refinement will break \n");

  } else {

    //save round for other simulations                                                                                                                                                                                                     \
                                                                                                                                                                                                                                              
    char saveTo[500];
    sprintf(saveTo,"%sIDshape=%i_IDdelta=%i_Round=%i_%s/round.txt",resBefore,IDshape,IDdelta,IDround,resAfter);

    FILE* pfile;
    pfile = fopen(saveTo,"w");
    fprintf(pfile,"%i \n",0);
    fclose(pfile);
    printf("Desired delta step has not been reached yet \n");

  }

}

//compute the norm of the perturbation which is equal to the norm of R(theta)-1
double GeometryDrop::ComputeDelta(double Dphi, int IDdelta, double *xStable, double *yStable, double *xUnstable, double *yUnstable){

  //check delta step and eventually put a flag to break succesive simuylations                                                                                                                                                              
  //get coordinate of the closest to the stable one                                                                                                                                                                                         
  double *xHere, *yHere;
  xHere = new double [DropElem+1];
  yHere = new double [DropElem+1];
  double phi = Dphi*(IDdelta-1.0);
  //average the two shapes                                                                                                                                                                                                                    
  for (int k=0; k<DropElem+1; k++){

    xHere[k] = xStable[k]*(1.0-phi) + xUnstable[k]*phi;
    yHere[k] = yStable[k]*(1.0-phi) + yUnstable[k]*phi;                                                                                                                                                                                  
  }

  double xcm = CenterMassAxis(xHere,yHere,DropElem);
  double INT = 0;
  double dtheta, r1, r2, xTemp1, yTemp1, xTemp2, yTemp2, theta1, theta2;
  //integration with trapezi rule                                                                                                                                                                                                           
  for (int k=0; k<DropElem; k++) {

    xTemp1 = xHere[k]-xcm;
    yTemp1 = yHere[k];
    xTemp2 = xHere[k+1]-xcm;
    yTemp2 = yHere[k+1];
    
    //trick to avoid fake angle due to imprecise y coordinate
    theta1 = atan(yTemp1/xTemp1);
    theta2 = atan(yTemp2/xTemp2);
    
    if (xTemp1<0){

      theta1 += M_PI;

    }

    if (xTemp2<0){

      theta2 += M_PI;

    }

    dtheta = theta2-theta1;
    
    r1 = sqrt(pow(xTemp1,2)+pow(yTemp1,2))-1;
    r2 = sqrt(pow(xTemp2,2)+pow(yTemp2,2))-1;

    INT += (pow(r1,2)+pow(r2,2))*dtheta*0.5;
    //printf("%f %f %f %f %f \n",INT,xcm,dtheta,r1,r2);
    //printf("%f %f %f  %i \n",dtheta, theta2, xTemp2, k); 

  }
  
  //free memory
  delete[] xHere;
  delete[] yHere;

  //compute delta as ||R-1|| with ||f||=1                                                                                                                                                                                                   
  double delta = sqrt(INT);
  return delta;

}

//extrapolate shapes from stable one                                                                                                                                          
void GeometryDrop::ExtrapolateShapes(double Dphi, int IDdelta, double *xStable, double *yStable, double deltaStable, double deltaUnstable, double *xUnstable, double *yUnstable){

  printf("New shape is the results of an extrapolation \n");

  //double *xUnstable, *yUnstable;
  //xUnstable = new double [DropElem+1];
  //yUnstable = new double [DropElem+1];

  //compute "unstable" shape
  double xcm = CenterMassAxis(xStable,yStable,DropElem);
  double theta;
  for (int k=0; k<DropElem+1; k++) {

    xStable[k] -= xcm;
    theta = atan(yStable[k]/xStable[k]) + M_PI*(xStable[k]<0);
    xUnstable[k] = (xStable[k] - cos(theta))/deltaStable*deltaUnstable + cos(theta) + xcm;
    yUnstable[k] = (yStable[k] - sin(theta))/deltaStable*deltaUnstable + sin(theta);

    //printf("%f %f %f %i \n",theta,xUnstable[k],yUnstable[k],k);

  }

}

//weigth average of shapes
void GeometryDrop::AverageShapes(double Dphi, int IDdelta, double *xStable, double *yStable, double *xUnstable, double *yUnstable){

  printf("New shape is the result of an interpolation \n");

  double phi = Dphi*(IDdelta-1.0);
  //average the two shapes                                                                                                                                                                                                                  
  for (int k=0; k<DropElem+1; k++){

    xCoord[k] = xStable[k]*(1.0-phi) + xUnstable[k]*phi;
    yCoord[k] = yStable[k]*(1.0-phi) + yUnstable[k]*phi;                                                                                                                                                                                                          

  }

}

//checl if ite get shape is too large (simulation already break) and eventually reduce it
int GeometryDrop::CheckIteGetShape(int IDshape, int IDdelta, int IDround, int step, char* resBefore, char* resAfter){

  char LoadFrom[500];                                                                                                                                                                                                                       
  sprintf(LoadFrom,"%sIDshape=%i_IDdelta=%i_Round=%i_%s/break.txt",resBefore,IDshape,IDdelta,IDround,resAfter);
  fstream CheckBreak(LoadFrom, std::ios_base::in);
  
  int myBreak, IteBreak;
  
  if (CheckBreak.is_open()){

    CheckBreak >> myBreak;
    CheckBreak >> IteBreak;

  }
  CheckBreak.close();

  int NewStep = step;
  int count = 0;
  while (NewStep>IteBreak) {
  
    NewStep = 0.9*NewStep;

    //printf("\n %i \n",NewStep);

    //avoid inifinite loop
    if (count>100) {

      printf("No sufficiently small step has been found");
      break;

    }

    count += 1;

  }

  return NewStep;

}

//get shapes at a certain time, first stable and first unstable
void GeometryDrop::GetShapesAtTime(int IDshape, int IDdelta, int IDround, int step, char* resBefore, char* resAfter, double *xStable, double *yStable, double *xUnstable, double *yUnstable){

  char LoadFrom[500];

  printf(" IDshape=%i IDdelta=%i-IDdelta=%i Round=%i at step %i \n", IDshape, IDdelta, IDdelta-1, IDround, step);	

  //the IDdelta of the first unstable shape is given by count (I take an average of the last stable and first unstable)                                                                                                                   
  //load first unstable                                                                                                                                                                                                      
  sprintf(LoadFrom,"%sIDshape=%i_IDdelta=%i_Round=%i_%s/drop%i.txt",resBefore,IDshape,IDdelta,IDround,resAfter,step);
  fstream myUnstable(LoadFrom, std::ios_base::in);
  myUnstable.is_open();
  //load last stable   
  sprintf(LoadFrom,"%sIDshape=%i_IDdelta=%i_Round=%i_%s/drop%i.txt",resBefore,IDshape,IDdelta-1,IDround,resAfter,step);
  fstream myStable(LoadFrom, std::ios_base::in);
  myStable.is_open();

  char str1[500];
  char str2[500];

  if (myStable.is_open()&&myUnstable.is_open()){
    int count = 0;
    while (count<DropElem+1) {

      myStable.getline(str1,500);
      myUnstable.getline(str2,500);

      myStable >> xStable[count] >> yStable[count];
      myUnstable >> xUnstable[count] >> yUnstable[count];

      //printf("ciao");
      //printf("%f %f %i \n",xStable[count],yStable[count], count);                                                                                                                                                                       

      count += 1;
    }
  }
  myStable.close();
  myUnstable.close();

}

//get trough the folder and figure out which is the first that got unstable
int GeometryDrop::GetFirstUnstable(int IDshape, int CountRound, char*resBefore, char* resAfter){
                                                                                                                                                             
  int myBreak = 0;
  int count = 0;

  char LoadFrom[500];
  while (myBreak==0){

    //load from this path                                                                                                                                                                                                                 
    sprintf(LoadFrom,"%sIDshape=%i_IDdelta=%i_Round=%i_%s/break.txt",resBefore,IDshape,count+1,CountRound,resAfter);
    ifstream myfile(LoadFrom, std::ios_base::in);
    myfile >> myBreak;				

    //avoid infinite                                                                                                                                                                                                                      
    if (count>100){			

      //too many loops                                                                                                                                                                                                                    
      printf("\n No treshold has been found, kill simulation \n");
      exit(0);					

    }
    myfile.close();
    count += 1;

    //printf("%i \n",count);                                                                                                                                                                                                              
  }

  return count;

}

//set which is the first round of previous shape which is refined enough
int GeometryDrop::SetPreviousRound(char *resBefore, char* resAfter, int IDshape){

  int CountRound = 0;
  int LoadRound = 0;

  //search last round of previous shape                                                                                                                                                                                                   
  char LoadFrom[256];
  //load from this path                                                                                                                                                                                                                   
  while (LoadRound==0){

    sprintf(LoadFrom,"%sIDshape=%i_IDdelta=1_Round=%i_%s/round.txt",resBefore,IDshape,CountRound+1,resAfter);
    fstream myfile(LoadFrom, std::ios_base::in);
    //printf("%s \n",LoadFrom);                                                                                                                                                                                                           
    myfile >> LoadRound;

    CountRound += 1;

    //printf("\n %i \n",LoadRound);                                                                                                                                                                                                       
    myfile.close();

  }

  return CountRound;

}

//set round.txt
void GeometryDrop::SetRound(int i, char* resu){

  char saveTo[500];
  sprintf(saveTo,"%s/round.txt",resu);

  FILE* pfile;
  pfile = fopen(saveTo,"w");
  fprintf(pfile,"%i \n",i);			
  fclose(pfile);

}

//get initial volume
double GeometryDrop::GetV0drop(){

  double V = V0drop;
  return V;

}

//nbreak if I come back to the stable shape
int GeometryDrop::CheckBreakStable(){

  int BREAK = 0;
  double xcm = CenterMassAxis(xCoord,yCoord,DropElem);
  double check = fabs(xCoord[0]-1-xcm);
  if (check<1e-3){

    printf("Drop has returned to the spherical configuration \n");
    BREAK = 1;

  }

  return BREAK;
}

//compute doplet ellpsicity
double GeometryDrop::ComputeEllipse(){

  //define variables
  double a, b, D, xcm;
  xcm = CenterMassAxis(xCoord,yCoord,DropElem);

  //define major and minor axis
  a = xCoord[0]-xcm;
  b = yCoord[DropElem/2+1];
  
  //printf("%f %f \n", a, b);

  //ellipsicuty
  D = (a-b)/(a+b);
  return D;

}

void GeometryDrop::CheckBreakDomain(){

  int outDomain = 0;
  //check if the interface goes trough the axis
  for (int i=0; i<DropElem+1; i++){

    outDomain += (yCoord[i]<0);

  }

  if (outDomain>0){
    printf("Interface escaped the domain");
  }
  //assert when interface goes though the domain
  assert(outDomain==0);

}

//break on
void GeometryDrop::BreakOn(char *resu, int ite){

  if (ite>1){
  char saveTo[256];
  sprintf(saveTo,"%s/break.txt",resu);

  FILE* pfile;
  pfile = fopen(saveTo,"w");
  fprintf(pfile,"%i \n",1); //break On
  fprintf(pfile,"%i \n",ite); //break On
  fclose(pfile);
  }

}

//break off
void GeometryDrop::BreakOff(char *resu, int ite){

  char saveTo[256];
  sprintf(saveTo,"%s/break.txt",resu);

  FILE* pfile;
  pfile = fopen(saveTo,"w");
  fprintf(pfile,"%i \n",0); //break Off
  fprintf(pfile,"%i \n",ite); //break Off
  fclose(pfile);

}

//execute volume correction                                                                                                                                                   
void GeometryDrop::VolumeCorrection(int kkk){                                                                                                                                                                                                                 

  if(kkk==1){
    double A, V, xCoordHere, yCoordHere;

    for (int i=0; i<10; i++){
    //compute area and volume of the drop                                                                                                                                                                                   
    ComputeAreaVolume(A,V);

    double V0 = GetV0drop();
    //display the interface based on Pranay                                                                                                                                                                                     
    double displace = -1.0*(V-V0)/A;

    double *nx = GetNormalX();
    double *ny = GetNormalY();
    //displace interface points                                                                                                                                                                                                               
    for (int k=0; k<DropElem+1;k++){

      xCoordHere = GetOneXCoord(k);
      yCoordHere = GetOneYCoord(k);
      xCoordHere += displace*nx[k];
      yCoordHere += displace*ny[k];

      SetOneCoord(xCoordHere,yCoordHere,k);

    }
    }                                                                                                                                                                                   

  } else if(kkk==2){

    double A, V, xCoordHere, yCoordHere, xcm;

    for (int i=0; i<1; i++){
      //compute area and volume of the drop                                                                                                                                                                                                     
      ComputeAreaVolume(A,V);
      xcm = CenterMassAxis(xCoord,yCoord,DropElem);

      double V0 = GetV0drop();
      //display the interface based on Pranay                                                                                                                                                                                                   
      double alpha = pow(V0/V,1.0/3);
      //printf("%f\n",alpha);

      double *nx = GetNormalX();
      double *ny = GetNormalY();
      //displace interface points                                                                                                                                                                                             
      
      for (int k=0; k<DropElem+1;k++){

	xCoordHere = GetOneXCoord(k);
	yCoordHere = GetOneYCoord(k);
	xCoordHere = (xCoordHere-xcm)*alpha + xcm;
	yCoordHere = yCoordHere*alpha;

	SetOneCoord(xCoordHere,yCoordHere,k);

      }
    }

  }

}

void GeometryDrop::ComputeAreaVolume(double &A,double &V){

  double *ax, *bx, *cx, *dx, *ay, *by, *cy, *dy;
  ax = new double [DropElem];
  bx = new double [DropElem];
  cx = new double [DropElem];
  dx = new double [DropElem];
  ay = new double [DropElem];
  by = new double [DropElem];
  cy = new double [DropElem];
  dy = new double [DropElem];

  SplineSymmetric(xCoord ,yCoord, ax, bx, cx, dx, ay, by, cy, dy, DropElem);

  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};

  A = 0;
  V = 0;
  double GPhere, GWhere, beta, deta, dbeta;
  //loop on elements
  for (int i=0; i<DropElem; i++){

    //loop on gauss points
    for (int k=0; k<6; k++) {

      GPhere = (GP[k]+1)/2;
      GWhere = GW[k]/2;

      beta = ay[i]+by[i]*GPhere+cy[i]*pow(GPhere,2)+dy[i]*pow(GPhere,3);
      deta = bx[i]+2*cx[i]*GPhere+3*dx[i]*pow(GPhere,2);
      dbeta = by[i]+2*cy[i]*GPhere+3*dy[i]*pow(GPhere,2);

      //area
      A += (beta*sqrt(deta*deta+dbeta*dbeta)*GWhere)*2*M_PI;

      //volume
      V += -beta*beta*deta*M_PI*GWhere;

    }

  }

  delete []ax;
  delete []bx;
  delete []cx;
  delete []dx;
  delete []ay;
  delete []by;
  delete []cy;
  delete []dy;

}

//print coordinates to screen
void GeometryDrop::PrintCoord () {

  for (int i=0; i<DropElem+1; i++){

    printf("x=%f y=%f node %i\n",xCoord[i],yCoord[i],i+1);

  }

}

//save coordinates to folder
void GeometryDrop::SaveCoord (char *resu, int ite) {

  char saveTo[128];
  sprintf(saveTo,"%s/drop%i.txt",resu,ite);

  FILE* pfile;
  pfile = fopen(saveTo,"w");
  fprintf(pfile,"x  y\n"); //header
  for (int i=0; i<DropElem+1; i++){

    fprintf(pfile,"%f  %f\n",xCoord[i],yCoord[i]);
    
  }
  fclose(pfile);

}

//update curvature vector
void GeometryDrop::UpdateCurv(double* ax, double* bx, double* cx, double* dx, double* ay, double* by, double* cy, double* dy, int i) {

  for (int k=0; k<i; k++){

    k1[k] = CurvSplineR(bx[k], cx[k], by[k], cy[k]);
    k2[k] = ny[k]/yCoord[k];

  }

  //first and last point are on the axis
  k2[0] = k1[0];
  k1[i] = CurvSplineL(bx[i-1], cx[i-1], dx[i-1], by[i-1], cy[i-1], dy[i-1]);
  k2[i] = k1[i];

}

//update curvature vector
void GeometryDrop::UpdateNormal(double* ax, double* bx, double* cx, double* dx, double* ay, double* by, double* cy, double* dy, int i) {
  
  for (int k=0; k<i; k++){

    nx[k] = NormalVectorNxR(bx[k], by[k]);
    ny[k] = NormalVectorNyR(bx[k], by[k]);

    //printf("%f %f \n",nx[i],ny[i]);

  }

  //first and last point are on the axis
  nx[i] = NormalVectorNxL(bx[i-1], cx[i-1], dx[i-1], by[i-1], cy[i-1], dy[i-1]);
  ny[i] = NormalVectorNyL(bx[i-1], cx[i-1], dx[i-1], by[i-1], cy[i-1], dy[i-1]);

  //printf("%f %f \n",nx[i],ny[i]);

}

//save coordinates to folder
void GeometryDrop::SaveDropData (char *resu, int ite) {

  char saveTo[246];
  sprintf(saveTo,"%s/drop%i.txt",resu,ite);

  time_t Tnow;
  time(&Tnow);
  double seconds = difftime(Tnow,Tstart);

  FILE* pfile;
  //save geometry data
  pfile = fopen(saveTo,"w");
  fprintf(pfile,"x  y  k1  k2 \n"); //header
  for (int i=0; i<DropElem+1; i++){

    fprintf(pfile,"%1.16f  %1.16f  %1.9f  %1.9f \n",xCoord[i],yCoord[i],k1[i],k2[i]);
    
  }
  fclose(pfile);

  FILE *pfile2;
  //save running time
  sprintf(saveTo,"%s/time%i.txt",resu,ite);
  pfile2 = fopen(saveTo,"w");
  fprintf(pfile2,"time\n"); //header
  fprintf(pfile2,"%1.9f\n",seconds);
  fclose(pfile2);

}

int GeometryDrop::getElem () {

  int l = DropElem;
  return l;
} 

double* GeometryDrop::GetXCoord () {

  double *x = xCoord;
  return x;
}

double GeometryDrop::GetOneXCoord (int i) {

  double x = xCoord[i];
  return x;
}

double GeometryDrop::GetOneOldXCoord (int i) {

  double x = OldxCoord[i];
  return x;
}

double GeometryDrop::GetOneOldYCoord (int i) {

  double x = OldyCoord[i];
  return x;
}

double* GeometryDrop::GetOldXCoord () {

  double *x = OldxCoord;
  return x;
}

double* GeometryDrop::GetYCoord () {

  double *y = yCoord;
  return y;
}

double GeometryDrop::GetOneYCoord (int i) {

  double y = yCoord[i];
  return y;
}

double* GeometryDrop::GetOldYCoord () {

  double *y = OldyCoord;
  return y;
}

void GeometryDrop::GetCoord(double *x, double *y) {

  x = xCoord;
  y = yCoord;
}

void GeometryDrop::SetCoord(double *x, double *y) {

  xCoord = x;
  yCoord = y;
  
}

void GeometryDrop::SetOneCoord(double x, double y, int i) {

  xCoord[i] = x;
  yCoord[i] = y;
  
}

void GeometryDrop::SetOldCoord(double *x, double *y) {

  OldxCoord = x;
  OldyCoord = y;
  
}

void GeometryDrop::SetOneOldCoord(double x, double y, int i) {

  OldxCoord[i] = x;
  OldyCoord[i] = y;
  
}

double* GeometryDrop::GetK1() {

  double *x = k1;
  return x;
}

double* GeometryDrop::GetK2() {

  double *x = k2;
  return x;
}

double* GeometryDrop::GetOldK1() {

  double *x = k1Old;
  return x;
}

double* GeometryDrop::GetOldK2() {

  double *x = k2Old;
  return x;
}

void GeometryDrop::SetCurv(double x, double y, int point) {

  k1[point] = x;
  k2[point] = y;
  
}

void GeometryDrop::SetOldCurv(double x, double y, int point) {

  k1Old[point] = x;
  k2Old[point] = y;
  
}

double* GeometryDrop::GetNormalX() {

  double *x = nx;
  return x;
}

double* GeometryDrop::GetNormalY() {

  double *x = ny;
  return x;
}

double* GeometryDrop::GetOldNormalX() {

  double *x = nxOld;
  return x;
}

double* GeometryDrop::GetOldNormalY() {

  double *x = nyOld;
  return x;
}

void GeometryDrop::SetNormal(double x, double y, int point) {

  nx[point] = x;
  ny[point] = y;
  
}

void GeometryDrop::SetOldNormal(double x, double y, int point) {

  nxOld[point] = x;
  nyOld[point] = y;
  
}
