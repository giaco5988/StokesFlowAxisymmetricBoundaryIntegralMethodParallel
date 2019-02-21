/*
 *  remesh.cpp
 *
 *  Created by Mathias Nagel on 28.01.16.
 *
 */
#include "remesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "combase.h"
#include <numeric>
#include "../alglib/src/interpolation.h"

using namespace alglib;

//constructor of the class when doing mesh stabilization
Remesh::Remesh(GeometryDrop* geo, int kind){

  //drop here
  //GeoHere = drop;
  //u = b;
  MESHkind = kind;
  int n = geo->getElem();
  unew = new double [n+1];
  vnew = new double [n+1];

  if (MESHkind==0){
    printf("Passive mesh stabilization is activated \n");
  }
  
}

double* Remesh::GetVelU(){

double* x = unew;
return x;

}

double* Remesh::GetVelV(){

double* x = vnew;
return x;

}

void Remesh::SetVelU(double x, int point){

  unew[point] = x;

}

void Remesh::SetVelV(double x, int point){

  vnew[point] = x;

}

int Remesh::GetMeshKind(){

  int x = MESHkind;
  return x;

}

double Remesh::GetOneVelU(int takeVel){

double x = unew[takeVel];
return x;

}

double Remesh::GetOneVelV(int takeVel){

double x = vnew[takeVel];
return x;

}

void Remesh::SetAdaptDrop(double r){

  g = r;

}

void Remesh::SetRK(int i){

  meshRK = i;

}

void Remesh::ExecuteRemesh(GeometryDrop* GeoHere, double* u){

  //passive mesh stabilization
  if (MESHkind==0){

    //perform passive mesh stabilization
    PassiveMeshStabilization(GeoHere, u);
  
  }

}

//set option for clustering
void Remesh::SetOptionCluster(int i){

  OptionCluster = i;

}

//set option for choosing distribution                                                                                                                                                                                                      
void Remesh::SetOptionDistribution(int i){

  OptionDistribution = i;

}

//set tune remesh adaptivity
void Remesh::SetTuneDistributionAdaptivity(double adapt){

  TuneDistributionAdaptivity = adapt;

}

//remesh globally using a selected distribution
void Remesh::RemeshDistribution(GeometryDrop* GeoHere, double *xOut, double *yOut){

  if (OptionDistribution==1){
    printf("Remesh globally, more nodes on right part 1 \n");
  } else if (OptionDistribution==2){
    printf("Remesh globally, more nodes on right part 2 \n");
  } else if (OptionDistribution==3){
    printf("Remesh globally, uniform distribution \n");
  } else if (OptionDistribution==4){
    printf("Remesh globally, more nodes on right part 3 \n");
  }

  for (int ccc = 0; ccc<10; ccc++) {

  int elem = GeoHere->getElem();
  double *k1, *k2, *nx, *ny, *x, *y;

  //get data from drop object                                                                                                                                                                                                     
  k1 = GeoHere->GetK1();
  x = GeoHere->GetXCoord();
  y = GeoHere->GetYCoord();

  //allocate memory
  double *distribution, *ax, *bx, *cx, *dx, *ay, *by, *cy, *dy, *DS, *cumDS, *cumDistr;
  distribution = new double [elem];
  ax = new double [elem];
  bx = new double [elem];
  cx = new double [elem];
  dx = new double [elem];
  ay = new double [elem];
  by = new double [elem];
  cy = new double [elem];
  dy = new double [elem];
  DS = new double [elem];
  cumDS = new double [elem+1];
  cumDistr = new double [elem];

  //Compute distribution
  double ratio = 0.1;
  double dist = 1;
  //double opt = 2;
  ChooseDistributionNodes(x, y, k1, distribution, TuneDistributionAdaptivity, elem, elem-1, dist, OptionDistribution, ratio); // dist and ratio are USELESS at the moment, keep for future development
  
  //compute splines coefficient
  SplineSymmetric(x,y,ax,bx,cx,dx,ay,by,cy,dy,elem);

  //compute length of single elements
  double dt = 1.0/1000;
  double t1, t2, dsHere, hx1, hy1, hx2, hy2, h1, h2;
  //t2 = 0;
  for (int i=0; i<elem+1; i++){

    dsHere = 0;
    t2 = 0;
    for (int k=0; k<1000; k++) {
      
      t1 = t2;
      t2 += dt;
      //trapezi rule on the subinterval
      hx1 = bx[i]+2*cx[i]*t1+3*dx[i]*pow(t1,2);
      hx2 = bx[i]+2*cx[i]*t2+3*dx[i]*pow(t2,2);
      hy1 = by[i]+2*cy[i]*t1+3*dy[i]*pow(t1,2);
      hy2 = by[i]+2*cy[i]*t2+3*dy[i]*pow(t2,2);
      h1 = sqrt(pow(hx1,2)+pow(hy1,2));
      h2 = sqrt(pow(hx2,2)+pow(hy2,2));
      dsHere += (h1+h2)*dt/2;
      //printf("%f \n",dt);

    }
    DS[i] = dsHere;
    //printf("%f %i \n",dsHere,i);

  }
  
  cumDS[0] = 0.0;
  //compute cumulative sum of ds and distribution
  for (int i=0; i<elem; i++){
   
    cumDS[i+1] = std::accumulate(DS,DS+i+1,0.0);
    //printf("%f %i \n",cumDS[i],i);

    cumDistr[i] = std::accumulate(distribution,distribution+i+1,0.0);
    //printf("%f %f %i \n",distribution[i],cumDistr[i],i); 

  }

  

  //SPLINE FITTING                                                                                                                                                                                                                 
  real_1d_array xxx;
  xxx.setlength(elem+1);
  real_1d_array yyy;
  yyy.setlength(elem+1);
  ae_int_t info;
  double v;
  spline1dinterpolant spline1;
  spline1dinterpolant spline2;
  spline1dfitreport rep;
  double rho;

  //fill vectors
  for (int i=0; i<elem+1; i++){

    xxx[i] = cumDS[i];
    yyy[i] = x[i];

  }
  
  //compute spline
  rho = +3.0;
  //spline1dfitpenalized(xxx, yyy, 50, rho, info, spline1, rep);
  spline1dbuildcubic(xxx, yyy, spline1);
  //printf("%d\n", int(info)); // EXPECTED: 1                                                                                                                                                                                                 
  //fill vectors                                                                                                                                                                                                                              
  for (int i=0; i<elem+1; i++){

    xxx[i] = cumDS[i];
    yyy[i] = y[i];

  }

  //compute spline                                                                                                                                                                                                                            
  rho = +3.0;
  //spline1dfitpenalized(xxx, yyy, 50, rho, info, spline2, rep);
  spline1dbuildcubic(xxx, yyy, spline2);
  //printf("%d\n", int(info)); // EXPECTED: 1 
  
                                                                                                                                                                                                                        
  //evaluate spline at the desired point
  double spacing = 0;
  for (int i=0; i<elem; i++){

    xOut[i] = spline1dcalc(spline1,spacing);
    yOut[i] = spline1dcalc(spline2,spacing);
    spacing = cumDistr[i]/cumDistr[elem-1]*cumDS[elem];
    //printf("%f %f %f %i \n",spacing,xOut[i],yOut[i],i);

  }
  //last term
  xOut[elem] = spline1dcalc(spline1,spacing);
  yOut[elem] = spline1dcalc(spline2,spacing);

  //printf("%f \n",spacing)
  
  //free memory
  delete[] distribution;
  delete[] ax;
  delete[] bx;
  delete[] cx;
  delete[] dx;
  delete[] ay;
  delete[] by;
  delete[] cy;
  delete[] dy;
  delete[] DS;
  delete[] cumDS;
  delete[] cumDistr;

  }

}

//like in Zinchenko 1996 and 2013
void Remesh::PassiveMeshStabilization(GeometryDrop* GeoHere, double* u){

  printf("Passive mesh stabilization \n");
  
  int elem = GeoHere->getElem();
  double *k1, *k2, *nx, *ny, *x, *y;
  
  if (meshRK==1){
    //get data from drop object
    k1 = GeoHere->GetK1();
    k2 = GeoHere->GetK2();
    nx = GeoHere->GetNormalX();
    ny = GeoHere->GetNormalY();
    x = GeoHere->GetXCoord();
    y = GeoHere->GetYCoord();
  } else if (meshRK==2) {
    k1 = GeoHere->GetOldK1();
    k2 = GeoHere->GetOldK2();
    nx = GeoHere->GetOldNormalX();
    ny = GeoHere->GetOldNormalY();
    x = GeoHere->GetOldXCoord();
    y = GeoHere->GetOldYCoord();
  }

  //local variables
  double dx1, dy1, dl1, dx2, dy2, dl2, dlSqua1, dlSqua2, du, dv, duij, dvij;
  double uij1, uij2, vij1, vij2, H1, H2, fx, fy, fxij, fyij, a, b, xi, eta, tempFx, tempFy, Fold, F;
  double *G , *saveFX, *saveFY, *U, *u1, *saveDU, *saveDV, *distribution;
  double K = 0;
  double h21i, h21j, h22i, h22j, hij21, hij22;
  double fxBEF = 0;
  double fyBEF = 0;
  double duBEF = 0;
  double dvBEF = 0;
  double aCoeff = 0;
  double bCoeff = 0;
  double cCoeff = 0;
  double dCoeff = 0;
  double eCoeff = 0;

  //memory allocation
  G = new double [elem+1];
  saveFX = new double [elem+1];
  saveFY = new double [elem+1];
  saveDU = new double [elem+1];
  saveDV = new double [elem+1];
  U = new double [2*elem+2];
  u1 = new double [2*elem+2];
  distribution = new double [elem+1];

  //choose distribution
  if (OptionCluster==3){

    double ratio = 0.1;
    double dist = 1;
    //double OptionDistribution = 1;
    ChooseDistributionNodes(x, y, k1, distribution, TuneDistributionAdaptivity, elem, elem, dist, OptionDistribution, ratio); //dist and ratio are USELESS AT THE MOMENT, keep for future development

  }

  //apply ditribution
  for (int i=0; i<elem+1; i++){

    //distibution that we try to mantain                                                                                                                                                                                                         //classical based on curvature
    if (OptionCluster==1) {
      G[i] = k1[i]*k1[i] + k2[i]*k2[i] + 0.004;
    }else if (OptionCluster==2){//cluster close to axis (if ther is a tail they go there)
      G[i] = k1[i]*k1[i] + 1.0/pow(y[i],2);
      if (i==0){
	G[0] = k1[0]*k1[0] + 1.0/pow(y[1],2);
      }
      if (i==elem){
	G[elem] = k1[elem]*k1[elem] + 1.0/pow(y[elem-1],2);
      } 
    }else if (OptionCluster==3){//choose distribution with certain criteria

      G[i] = distribution[i];

    }
    
  }

  for (int i=1; i<elem; i++){

    //edge before
    dx1 = x[i]-x[i-1];
    dy1 = y[i]-y[i-1];
    dl1 = sqrt(dx1*dx1 + dy1*dy1);
    
    //edge after
    dx2 = x[i+1]-x[i];
    dy2 = y[i+1]-y[i];
    dl2 = sqrt(dx2*dx2 + dy2*dy2);

    K += pow(G[i],g)*(pow(dl1,2)+pow(dl2,2));

  }
  //first and last term are different because it is an open line
  dx1 = x[1]-x[0];
  dy1 = y[1]-y[0];
  dl1 = sqrt(dx1*dx1 + dy1*dy1);
  K += pow(G[0],g)*pow(dl1,2);
  
  dx2 = x[elem]-x[elem-1];
  dy2 = y[elem]-y[elem-1];
  dl2 = sqrt(dx2*dx2 + dy2*dy2);
  K += pow(G[elem],g)*pow(dl2,2);

  K = 0.5/elem*K;
  ////////////////////// STEEPEST DESCENT METHOD ////////////////////
  //loop on the nodes, the main one, I don't treat first and last one
  for (int i=1; i<elem; i++){

    //dx and dy for edge before
    dx1 = x[i]-x[i-1];
    dy1 = y[i]-y[i-1];
    dlSqua1 = dx1*dx1+dy1*dy1;

    //dx and dy for edge after
    dx2 = x[i+1]-x[i];
    dy2 = y[i+1]-y[i];
    dlSqua2 = dx2*dx2+dy2*dy2;
    
    //hij2 for edge before
    h21i = K*pow(G[i-1],-g);
    h21j = K*pow(G[i],-g);
    hij21 = 0.5*(h21i+h21j);

    //hij2 for edge after
    h22i = K*pow(G[i],-g);
    h22j = K*pow(G[i+1],-g);
    hij22 = 0.5*(h22i+h22j);

    //vij for edge before
    uij1 = u[2*i]-u[2*i-2];
    vij1 = u[2*i+1]-u[2*i-1];

    //vij for edge after
    uij2 = u[2*(i+1)]-u[2*i];
    vij2 = u[2*i+3]-u[2*i+1];

    //H for edge before
    H1 = 4*pow(1/hij21-hij21/pow(dlSqua1,2),2);

    //H for edge after
    H2 = 4*pow(1/hij22-hij22/pow(dlSqua2,2),2);
    
    //f
    fx = 2*H1*dx1*uij1*dx1 - 2*H2*dx2*uij2*dx2;
    fy = 2*H1*dy1*vij1*dy1 - 2*H2*dy2*vij2*dy2;

    //remove normal component
    tempFx = fx - (fx*nx[i]+fy*ny[i])*nx[i];
    tempFy = fy - (fx*nx[i]+fy*ny[i])*ny[i];
    fx = tempFx;
    fy = tempFy;
    
    fxij = fx-fxBEF;
    fyij = fy-fyBEF;

    //compute parabola coefficients
    aCoeff += H1*pow(dx1*fxij+dy1*fyij,2);
    bCoeff += 2*H1*((dx1*uij1+dy1*vij1) * (dx1*fxij+dy1*fyij));

    fxBEF = fx;
    fyBEF = fy;

    saveFX[i] = fx;
    saveFY[i] = fy;
    
  }
  //compute parabola coefficients, add last term
  fxij = 0.0-fxBEF;
  fyij = 0.0-fyBEF;
  aCoeff += H2*pow(dx2*fxij+dy2*fyij,2);
  bCoeff += 2*H2*((dx2*uij2+dy2*vij2) * (dx2*fxij+dy2*fyij));

  xi = -bCoeff/2/aCoeff;
  double uTemp, vTemp, fxTemp, fyTemp;

  //printf("%f %f\n",aCoeff,bCoeff);

  //new velocities as loacal and global variables
  for (int i=0; i<elem+1; i++){

    uTemp = u[2*i];
    vTemp = u[2*i+1];
    fxTemp = saveFX[i];
    fyTemp = saveFY[i];
    U[2*i] = uTemp + xi*fxTemp;
    U[2*i+1] = vTemp + xi*fyTemp;

  }

  //to avoid leckage NOT CLEAR WHY
  U[0] = u[0];
  U[1] = 0.0;
  U[2*elem] = u[2*elem];
  U[2*elem+1] = 0.0;

  //compute the function tha is being minimized
  Fold = computeF(x, y, G, u, K, elem);
  F = computeF(x, y, G, U, K, elem);
  
  //store old velocity for next iteration
    for (int i=0; i<2*elem+2;i++){
      u1[i] = u[i];
      u[i] = U[i];
    }

  ////////////////////// CONJUGATE GRADIENT METHOD ////////////////////
  int count = 0;
  
  while(fabs(F-Fold) > fabs(F)*1e-6) {

    Fold = F;

    //reset to zero
    aCoeff = 0;
    bCoeff = 0;
    cCoeff = 0;
    dCoeff = 0;
    eCoeff = 0;
    fxBEF = 0;
    fyBEF = 0;
    duBEF = 0;
    dvBEF = 0;

    //loop on the nodes, the main one, I don't treat first and last one
    for (int i=1; i<elem; i++){

      //dx and dy for edge before
      dx1 = x[i]-x[i-1];
      dy1 = y[i]-y[i-1];
      dlSqua1 = dx1*dx1+dy1*dy1;

      //dx and dy for edge after
      dx2 = x[i+1]-x[i];
      dy2 = y[i+1]-y[i];
      dlSqua2 = dx2*dx2+dy2*dy2;
    
      //hij2 for edge before
      h21i = K*pow(G[i-1],-g);
      h21j = K*pow(G[i],-g);
      hij21 = 0.5*(h21i+h21j);

      //hij2 for edge after
      h22i = K*pow(G[i],-g);
      h22j = K*pow(G[i+1],-g);
      hij22 = 0.5*(h22i+h22j);

      //vij for edge before
      uij1 = u[2*i]-u[2*i-2];
      vij1 = u[2*i+1]-u[2*i-1];

      //vij for edge after
      uij2 = u[2*(i+1)]-u[2*i];
      vij2 = u[2*i+3]-u[2*i+1];

      //H for edge before
      H1 = 4*pow(1/hij21-hij21/pow(dlSqua1,2),2);

      //H for edge after
      H2 = 4*pow(1/hij22-hij22/pow(dlSqua2,2),2);
    
      //f
      fx = 2*H1*dx1*uij1*dx1 - 2*H2*dx2*uij2*dx2;
      fy = 2*H1*dy1*vij1*dy1 - 2*H2*dy2*vij2*dy2;

      //compute du and dv
      du = u[2*i]-u1[2*i];
      dv = u[2*i+1]-u1[2*i+1];

      //compute duij and dvij
      duij = du-duBEF;
      dvij = dv-dvBEF;

      //remove normal component
      tempFx = fx - (fx*nx[i]+fy*ny[i])*nx[i];
      tempFy = fy - (fx*nx[i]+fy*ny[i])*ny[i];
      fx = tempFx;
      fy = tempFy;

      //compute fxij and fyij
      fxij = fx-fxBEF;
      fyij = fy-fyBEF;

      //compute paraboloid coefficients
      aCoeff += H1*pow(dx1*fxij+dy1*fyij,2);
      bCoeff += 2*H1*((dx1*fxij+dy1*fyij) * (dx1*duij+dy1*dvij));
      cCoeff += H1*pow(dx1*duij+dy1*dvij,2);
      dCoeff += 2*H1*((dx1*fxij+dy1*fyij) * (dx1*uij1+dy1*vij1));
      eCoeff += 2*H1*((dx1*uij1+dy1*vij1) * (dx1*duij+dy1*dvij));

      fxBEF = fx;
      fyBEF = fy;
      duBEF = du;
      dvBEF = dv;

      saveFX[i] = fx;
      saveFY[i] = fy;
      saveDU[i] = du;
      saveDV[i] = dv;

      //printf("%f %f \n",fxij,fyij);
    
    }
    //compute paraboloid coefficients, add last term
    fxij = 0.0-fxBEF;
    fyij = 0.0-fyBEF;
    duij = 0.0-duBEF;
    dvij = 0.0-dvBEF;
    aCoeff += H2*pow(dx2*fxij+dy2*fyij,2);
    bCoeff += 2*H2*((dx2*fxij+dy2*fyij) * (dx2*duij+dy2*dvij));
    cCoeff += H2*pow(dx2*duij+dy2*dvij,2);
    dCoeff += 2*H2*((dx2*fxij+dy2*fyij) * (dx2*uij2+dy2*vij2));
    eCoeff += 2*H2*((dx2*uij2+dy2*vij2) * (dx2*duij+dy2*dvij));
    
    //phi and eta
    xi = (-2*cCoeff*dCoeff+bCoeff*eCoeff)/(4*aCoeff*cCoeff-bCoeff*bCoeff);
    eta = (-2*aCoeff*eCoeff+bCoeff*dCoeff)/(4*aCoeff*cCoeff-bCoeff*bCoeff);

    //new velocities as local variables
    for (int i=0; i<elem+1; i++){

      uTemp = u[2*i];
      vTemp = u[2*i+1];
      fxTemp = saveFX[i];
      fyTemp = saveFY[i];
      U[2*i] = uTemp + xi*fxTemp + eta*saveDU[i];
      U[2*i+1] = vTemp + xi*fyTemp + eta*saveDV[i];

    }

    //to avoid leckage NOT CLEAR WHY
    U[0] = u[0];
    U[1] = 0.0;
    U[2*elem] = u[2*elem];
    U[2*elem+1] = 0.0;

    //compute function F, consition for convergence
    F = computeF(x, y, G, U, K, elem);

    //store old velocity for next iteration
    for (int i=0; i<2*elem+2;i++){
      u1[i] = u[i];
      u[i] = U[i];
    }

    //break if too many iteration
    //if (count>500) {

    //printf("Force break in mesh stabilization, too many iteration to find the minimum \n");
    //break;

    //}
    
    //printf("%i \n",count);
    count += 1;

  }
  //end of conjugate gradient method
  

  //send to global variables
  for (int i=0; i<elem+1; i++){

    uTemp = u[2*i];
    vTemp = u[2*i+1];
    fxTemp = saveFX[i];
    fyTemp = saveFY[i];
    unew[i] = uTemp + xi*fxTemp;
    vnew[i] = vTemp + xi*fyTemp;

  }

  //to avoid leckage NOT CLEAR WHY
  unew[0] = u[0];
  vnew[0] = 0.0;
  unew[elem] = u[2*elem];
  vnew[elem] = 0.0;
  
  //free memory
  delete[] u1;
  delete[] U;
  delete[] G;
  delete[] saveFX;
  delete[] saveFY;
  delete[] saveDU;
  delete[] saveDV;
  delete[] distribution;

}

double Remesh::computeF(double *x, double *y, double *G, double *u, double K, int elem){

  double h21, h22, hij2, uij, vij, dx, dy, dl2;
  double H = 0;
  double F = 0;

  for (int i=0; i<elem; i++){

    h21 = K*pow(G[i],-g);
    h22 = K*pow(G[i+1],-g);
    hij2 = 0.5*(h21+h22);
    uij = u[2*i+2]-u[2*i];
    vij = u[2*i+3]-u[2*i+1];
    dx = x[i+1]-x[i];
    dy = y[i+1]-y[i];
    dl2 = dx*dx+dy*dy;

    H = 4*pow(1/hij2-hij2/(pow(dl2,2)),2);

    F += H*pow((dx*uij+dy*vij),2);
    
  }

  return F;

}

