/*
 *  combase.cpp
 *  Ulambator v0.3
 *
 *  Created by Mathias Nagel on 29.01.10.
 *  Copyright 2010 EPFL. All rights reserved.
 *
 */

#include "combase.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../alglib/src/interpolation.h"

using namespace alglib;

void bessk01(double x, double &k0, double &k1)
{
    double t,t2,dtmp,dtmp1,i0,i1;
	
//    if (x == 0.0) {
//        k0 = 1e308;
//        k1 = 1e308;
//    }
    if (x < 3.75) {
        t = x/3.75;
        t2 = t*t;
        i0 = (((((0.0045813*t2+0.0360768)*t2+0.2659732)*t2+
				1.2067492)*t2+3.0899424)*t2+3.5156229)*t2+1.0;
        i1 = x*(((((0.00032411*t2+0.00301532)*t2+0.02658733*t2+
				   0.15084934)*t2+0.51498869)*t2+0.87890594)*t2+0.5);
    }
    else {
        t = 3.75/x;
        dtmp1 = exp(x)/sqrt(x);
        dtmp = (((((((0.00392377*t-0.01647633)*t+0.026355537)*t-0.02057706)*t+
				   0.00916281)*t-0.00157565)*t+0.00225319)*t+0.01328592)*t+0.39894228;
        i0 = dtmp*dtmp1;
        dtmp = (((((((-0.00420059*t+0.01787654)*t-0.02895312)*t+0.02282967)*t-
				   0.01031555)*t+0.00163801)*t-0.00362018)*t-0.03988024)*t+0.39894228;
        i1 = dtmp*dtmp1;
    }
    if (x < 2.0) {
        t = 0.5*x;
        t2 = t*t;    // already calculated above
        dtmp = (((((0.0000074*t2+0.0001075)*t2+0.00262698)*t2+0.0348859)*t2+
				 0.23069756)*t2+0.4227842)*t2-0.57721566;
        k0 = dtmp - i0*log(t);
        dtmp = (((((-0.00004686*t2-0.00110404)*t2-0.01919402)*t2-
				  0.18156897)*t2-0.67278578)*t2+0.15443144)*t2+1.0;
        k1 = dtmp/x + i1*log(t);
    }
    else {
        t = 2.0/x;
        dtmp1 = exp(-x)/sqrt(x);
        dtmp = (((((0.00053208*t-0.0025154)*t+0.00587872)*t-0.01062446)*t+
				 0.02189568)*t-0.07832358)*t+1.25331414;
        k0 = dtmp*dtmp1;
        dtmp = (((((-0.00068245*t+0.00325614)*t-0.00780353)*t+0.01504268)*t-
				 0.0365562)*t+0.23498619)*t+1.25331414;
        k1 = dtmp*dtmp1;
    }
}

void bessi12(double x,double &i1,double &i2)
{
    double t,t2,dtmp,dtmp1,i0;
	
    if (x == 0.0) {
        i0 = 1.0;
        i1 = 0.0;
		i2 = 0.0;
    }
    if (x < 3.75) {
        t = x/3.75;
        t2 = t*t;
        i0 = (((((0.0045813*t2+0.0360768)*t2+0.2659732)*t2+
				1.2067492)*t2+3.0899424)*t2+3.5156229)*t2+1.0;
        i1 = x*((((((0.00032411*t2+0.00301532)*t2+0.02658733)*t2+
				   0.15084934)*t2+0.51498869)*t2+0.87890594)*t2+0.5);
    }
    else {
        t = 3.75/x;
        dtmp1 = exp(x)/sqrt(x);
        dtmp = (((((((0.00392377*t-0.01647633)*t+0.026355537)*t-0.02057706)*t+
				   0.00916281)*t-0.00157565)*t+0.00225319)*t+0.01328592)*t+0.39894228;
        i0 = dtmp*dtmp1;
        dtmp = (((((((-0.00420059*t+0.01787654)*t-0.02895312)*t+0.02282967)*t-
				   0.01031555)*t+0.00163801)*t-0.00362018)*t-0.03988024)*t+0.39894228;
        i1 = dtmp*dtmp1;
    }

	i2 = i0-2.0*i1/x;
}

//integration on linear elements of green functions single layer axisymmetric
void computeGTlinear1layer(double x0, double y0, double xMid,double yMid,double dx,double dy,double &GXXa, double &GXYa,double &GYXa,double &GYYa, double &GXXb, double &GXYb, double &GYXb, double &GYYb){

  double x, y, dl;
  double SXX, SXY, SYX, SYY;
  double phiA, phiB;

  //check if I'm on the axis
  bool axis = (y0 < 10e-8);
          
  //decide if doing singular treatment
  double distXA = (xMid-dx/2)-x0;
  double distYA = (yMid-dy/2)-y0;
  double distA = sqrt(distXA*distXA+distYA*distYA);
  bool dolog1 = (distA>10e-8);
          
  //decide if doing singular treatment
  double distXB = (xMid+dx/2)-x0;
  double distYB = (yMid+dy/2)-y0;
  double distB = sqrt(distXB*distXB+distYB*distYB);
  bool dolog2 = (distB>10e-8);
          
  bool dolog = (dolog1+dolog2 > 1.9);

  dl = sqrt(dx*dx+dy*dy);
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};

  //gauss integration
  for (int l=0; l<6; l++){

    //hat function
    phiA = 1-(GP[l]+1)/2;
    phiB = (GP[l]+1)/2;

    //current gauss point
    x = GP[l]/2*dx+xMid;
    y = GP[l]/2*dy+yMid;

    //printf("%f\n",y0);

    //compute green functions
    computeS(x,y,x0,y0,dolog,axis,SXX,SXY,SYX,SYY);

    //printf("%f\n",y0);
    
    //printf("%f %f %f %f %f %i\n",x, y, x0, y0, SXX ,dolog1);

    //integration on phiA
    GXXa += phiA*GW[l]*dl*SXX/2;
    GXYa += phiA*GW[l]*dl*SXY/2;
    GYXa += phiA*GW[l]*dl*SYX/2;
    GYYa += phiA*GW[l]*dl*SYY/2;

    //integration on phiB
    GXXb += phiB*GW[l]*dl*SXX/2;
    GXYb += phiB*GW[l]*dl*SXY/2;
    GYXb += phiB*GW[l]*dl*SYX/2;
    GYYb += phiB*GW[l]*dl*SYY/2;

    //if (l==5){
    //printf("%f \n",GXXa);
    //}

  }

  //singular treatment: add part form analytical integration (as in research notes) IF I'M ON THE AXIS I DON'T NEED SINGULAR TREATMENT
  GXXa += (1.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
  GYYa += (1.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          
  GXXb += (0.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
  GYYb += (0.5*dl-dl*log(dl))*(1-dolog1)*(1-axis);
          
  GXXa += (0.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
  GYYa += (0.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
          
  GXXb += (1.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);
  GYYb += (1.5*dl-dl*log(dl))*(1-dolog2)*(1-axis);

  //printf("%f \n",GXYb);
  
}

//choose distribution floowing different criteria
void ChooseDistributionNodes(double *xxx, double *yyy, double *K, double *distrRES, double tune, int elem, int elemNEW, double dist, int opt, double ratio){

  //printf("ciao %i \n",opt);

  //allocate dynamic memory
  double* temp;
  temp = new double [elem+1];

  if (opt==1) {  //more mesh on the rigth part of the droplet

    //printf("Distribution number 1 on right part of the droplet \n");
    
    for (int i=0; i<elem+1; i++){
        
      //temp[i] = xxx[i]-xxx[elem] + (xxx[0]-xxx[elem])*tune; //avoid negative numbers
      temp[i] = pow(xxx[i]-xxx[elem] + (xxx[0]-xxx[elem]),tune);

    }
        
  } else if (opt==2){ //even more mesh on right part of the droplet

    double l, l0;
    l0 = 2*(xxx[0]-xxx[elem]);
    //printf("Distribution number 2 on right part of the droplet \n");

    for (int i=0; i<elem+1; i++){

      if (i==0){
	l = sqrt(pow(xxx[i]-xxx[elem],2)+pow(yyy[i]-yyy[i-1],2));
      } else {
	l += sqrt(pow(xxx[i]-xxx[i-1],2)+pow(yyy[i]-yyy[i-1],2));
      }
      temp[i] = 1/pow(l+l0,tune);

    }

  } else if (opt==3) {

    //printf("Use uniform ditribution \n");

    for (int i=0; i<elem+1; i++) {
      temp[i] = 1.0;
    }

  } else if (opt==4) { //a lot of clustering on the right part of the droplet

    double distRight;
    for (int i=0; i<elem+1; i++) {

      //distance from right tip
      distRight = xxx[0]-xxx[i];

      //distribution
      temp[i] = 1/(exp(distRight+tune)-1);
      //printf("%f %f %i \n",distRight,temp[i],i);

    }

  }

  //SPLINE FITTING
  real_1d_array x;
  x.setlength(elem+1);
  real_1d_array y;
  y.setlength(elem+1);
  ae_int_t info;
  double v;
  spline1dinterpolant s;
  spline1dfitreport rep;
  double rho;

  //fill vectors
  for (int i=0; i<elem+1; i++){

    x[i] = i;
    y[i] = temp[i];

    //printf("%f %f %i \n",xxx[0],y[i],i);

  }

  // In real life applications you may need some moderate degree of fitting,
  // so we try to fit once more with rho=3.0.
  rho = +3.0;
  spline1dfitpenalized(x, y, 50, rho, info, s, rep);
  //printf("%d\n", int(info)); // EXPECTED: 1

  double dq = elem/elemNEW;
  //printf("%f \n",dq);
  double interrogate;
  //evaluate spline at the desired point
  for (int i=0; i<elemNEW+1; i++){

    interrogate = i*dq;
    distrRES[i] = 1.0/spline1dcalc(s,interrogate);

  }
  
  delete[] temp;

}

//integration on linear elements of green functions single layer axisymmetric
void computeGTlinearSPlines1layer(double x0, double y0, double h0,double h1, double ax, double bx, double cx, double dx, double ay, double by, double cy, double dy, double &GXXa, double &GXYa,double &GYXa,double &GYYa, double &GXXb, double &GXYb, double &GYXb, double &GYYb){

  double SXX, SXY, SYX, SYY;
  double phiA, phiB, GPpar, h;
  double eta, beta, deta, dbeta;
  double nx, ny;

  //check if I'm on the axis
  bool axis = (y0 < 10e-8);
          
  //decide if doing singular treatment
  double distXA = ax-x0;
  double distYA = ay-y0;
  double distA = sqrt(distXA*distXA+distYA*distYA);
  bool dolog1 = (distA>10e-8);
          
  //decide if doing singular treatment
  double distXB = ax+bx+cx+dx-x0;
  double distYB = ay+by+cy+dy-y0;
  double distB = sqrt(distXB*distXB+distYB*distYB);
  bool dolog2 = (distB>10e-8);
  
  //Gauss points and Gauss weigths
  double GP[] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
  double GW[] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};

  //gauss integration
  for (int l=0; l<6; l++){

    //gauss point in parametric space
    GPpar = (GP[l]+1)/2;
                                
    //gauss points in physical space
    eta = ax+bx*GPpar+cx*GPpar*GPpar+dx*GPpar*GPpar*GPpar;
    beta = ay+by*GPpar+cy*GPpar*GPpar+dy*GPpar*GPpar*GPpar;
    deta = bx+2*cx*GPpar+3*dx*GPpar*GPpar;
    dbeta = by+2*cy*GPpar+3*dy*GPpar*GPpar;
                
    //normal to the spline
    nx = dbeta/sqrt(deta*deta+dbeta*dbeta);
    ny = -deta/sqrt(deta*deta+dbeta*dbeta);

    //metrics term
    h = sqrt(deta*deta+dbeta*dbeta);

    //hat function
    phiA = 1-(GP[l]+1)/2;
    phiB = (GP[l]+1)/2;

    //printf("%f\n",y0);

    //compute green functions
    computeSSplines(GPpar,h,h0,h1,eta,beta,x0,y0,axis,dolog1,dolog2,SXX,SXY,SYX,SYY);
    //printf("%f %f %f %f %f %i\n",x, y, x0, y0, SXX ,dolog1);

    //integration on phiA
    GXXa += phiA*GW[l]*SXX/2 * h;
    GXYa += phiA*GW[l]*SXY/2 * h;
    GYXa += phiA*GW[l]*SYX/2 * h;
    GYYa += phiA*GW[l]*SYY/2 * h;

    //integration on phiB
    GXXb += phiB*GW[l]*SXX/2 * h;
    GXYb += phiB*GW[l]*SXY/2 * h;
    GYXb += phiB*GW[l]*SYX/2 * h;
    GYYb += phiB*GW[l]*SYY/2 * h;

  }

  //singular treatment: add part form analytical integration (as in research notes) IF I'M ON THE AXIS I DON'T NEED SINGULAR TREATMENT
  GXXa += 1.5*h0*(1-dolog1)*(1-axis);
  GYYa += 1.5*h0*(1-dolog1)*(1-axis);
          
  GXXb += 0.5*h0*(1-dolog1)*(1-axis);
  GYYb += 0.5*h0*(1-dolog1)*(1-axis);
          
  GXXa += 0.5*h1*(1-dolog2)*(1-axis);
  GYYa += 0.5*h1*(1-dolog2)*(1-axis);
          
  GXXb += 1.5*h1*(1-dolog2)*(1-axis);
  GYYb += 1.5*h1*(1-dolog2)*(1-axis);

  //printf("%f \n",GXYb);
  
}

//compute greem's function single layer axisymmetric
void computeS(double x,double y,double x0,double y0,bool dolog,bool axis,double &SXX, double &SXY, double &SYX, double &SYY){

  //printf("%f %f %f %f\n",x, y, x0, y0);
  
  //compute green function physical variables
  double DX = x-x0;
  double DY = y-y0;
  double d2 = DX*DX+DY*DY;
  double d = sqrt(d2);
  double r0 = y0;
  double r = y;
  double k2 = (4*r0*r)/(DX*DX+pow((r+r0),2));
  double k = sqrt(k2);
  double r2 = r*r;
  double r3 = r*r*r;
  
  //dummy variables
  double k2p = 1-k2;
  double k4 = k2*k2;
  double k6 = k2*k2*k2;
  double k5 = k4*k;
  double yy5 = pow(r0*r,2.5);
  double DX2 = DX*DX;
  double DX3 = DX*DX*DX;
                
  //compute elliptic integral of first and second kind
  double tol = 0.0000000000001;
  double F = 0.5*M_PI;
  double E = 1.0;
  double G = 1.0;
  double B = k;
  double D = 10*tol;
  double C;

  //while (fabs(D) > tol) {
  while (sqrt(D*D) > tol) {
      C = sqrt(1.0-B*B);
      B = (1.0-C)/(1.0+C);
      D = F*B;
      F = F+D;
      G = 0.50*G*B;
      E = E+G;
  }

  E = F*(1.0-0.50*k2*E);
                
  //often used variables for single layer
  double I10 = 2.0*k/sqrt(r0*r)*F;
  double I11 = k/pow(r0*r,1.5)*((DX2+r0*r0+r2)*F-(DX2+r0*r0+r2+2*r*r0)*E);
  double I30 = 2*k/pow(r0*r,0.5)*E/d2;
  double I31 = k/pow(r0*r,1.5)*(-F+(r0*r0+r2+DX2)/d2*E);
  double I32 = k/pow(r0*r,2.5)*(-F*(DX2+r0*r0+r2) + (DX2*DX2+2.0*DX2*r2+2.0*DX2*r0*r0+r2*r2+r0*r0*r0*r0)*E/d2);
                
  if (axis) {
                    
      double d3  = d2*d;
      double d5  = d3*d2;
      I10 = 2.0*M_PI/d;
      I11 = 0.0;
      I30 = 2*M_PI/d3;
      I31 = 0.0;
      I32 = M_PI/d3;
      
  }
                
  //single layer
  SXX = r*(I10+DX*DX*I30) + 2*log(d)*(1-dolog)*(1-axis);
  SXY = r*DX*(r*I30-r0*I31);
  SYX = r*DX*(r*I31-r0*I30);
  SYY = r*(I11+(r*r+r0*r0)*I31-r*r0*(I30+I32)) + 2*log(d)*(1-dolog)*(1-axis);
  
}

//compute greem's function single layer axisymmetric
void computeSSplines(double GPpar, double h, double h0, double h1, double x, double y, double x0, double y0, bool axis, bool dolog1, bool dolog2, double &SXX, double &SXY, double &SYX, double &SYY){
  
  //compute green function physical variables
  double DX = x-x0;
  double DY = y-y0;
  double d2 = DX*DX+DY*DY;
  double d = sqrt(d2);
  double r0 = y0;
  double r = y;
  double k2 = (4*r0*r)/(DX*DX+pow((r+r0),2));
  double k = sqrt(k2);
  double r2 = r*r;
  double r3 = r*r*r;
  
  //dummy variables
  double k2p = 1-k2;
  double k4 = k2*k2;
  double k6 = k2*k2*k2;
  double k5 = k4*k;
  double yy5 = pow(r0*r,2.5);
  double DX2 = DX*DX;
  double DX3 = DX*DX*DX;
                
  //compute elliptic integral of first and second kind
  double tol = 0.0000000000001;
  double F = 0.5*M_PI;
  double E = 1.0;
  double G = 1.0;
  double B = k;
  double D = 10*tol;
  double C;

  //while (fabs(D) > tol) {
  while (sqrt(D*D) > tol) {
      C = sqrt(1.0-B*B);
      B = (1.0-C)/(1.0+C);
      D = F*B;
      F = F+D;
      G = 0.50*G*B;
      E = E+G;
  }

  E = F*(1.0-0.50*k2*E);
                
  //often used variables for single layer
  double I10 = 2.0*k/sqrt(r0*r)*F;
  double I11 = k/pow(r0*r,1.5)*((DX2+r0*r0+r2)*F-(DX2+r0*r0+r2+2*r*r0)*E);
  double I30 = 2*k/pow(r0*r,0.5)*E/d2;
  double I31 = k/pow(r0*r,1.5)*(-F+(r0*r0+r2+DX2)/d2*E);
  double I32 = k/pow(r0*r,2.5)*(-F*(DX2+r0*r0+r2) + (DX2*DX2+2.0*DX2*r2+2.0*DX2*r0*r0+r2*r2+r0*r0*r0*r0)*E/d2);
                
  if (axis) {
                    
      double d3  = d2*d;
      double d5  = d3*d2;
      I10 = 2.0*M_PI/d;
      I11 = 0.0;
      I30 = 2*M_PI/d3;
      I31 = 0.0;
      I32 = M_PI/d3;
      
  }
                
  //single layer
  SXX = r*(I10+DX*DX*I30);
  SXX += (2*log(d) - 2*log(d/GPpar) + log(GPpar)*(1-h0/h))*(1-dolog1)*(1-axis);     //remove singularity on the left
  SXX += (2*log(d) - 2*log(d/(1-GPpar)) + log(1-GPpar)*(1-h1/h))*(1-dolog2)*(1-axis);     //remove singularity on the right
                
  SXY = r*DX*(r*I30-r0*I31);
  SYX = r*DX*(r*I31-r0*I30);
                
  SYY = r*(I11+(r*r+r0*r0)*I31-r*r0*(I30+I32));
  SYY += (2*log(d) - 2*log(d/GPpar) + log(GPpar)*(1-h0/h))*(1-dolog1)*(1-axis);     //remove singularity on the left
  SYY += (2*log(d) - 2*log(d/(1-GPpar)) + log(1-GPpar)*(1-h1/h))*(1-dolog2)*(1-axis);     //remove singularity on the right
  
}

//compute center of mass of an axisymmetric droplet using trapezi rule for integration
double CenterMassAxis(double *x, double *y, int n){

  double integ;
  double dx;
  double first;
  double second;
  double num = 0;
  double den = 0;

  //perform integral with trapezi rule
  for(int i=0; i<n; i++){

    dx = x[i+1]-x[i];

    first = pow(y[i],2)*x[i];
    second = pow(y[i+1],2)*x[i+1];
    num += (first+second)*dx/2;

    first = pow(y[i],2);
    second = pow(y[i+1],2);
    den += (first+second)*dx/2;

  }

  integ = num/den;

  return integ;

}

//compute spline coefficient using splines perpendicular to the axis
void SplineSymmetric(double *x ,double *y, double *ax, double *bx, double *cx, double *dx, double *ay, double *by, double *cy, double *dy, int n){

  double *deltaX;
  double *deltaY;
  double *A;
  
  deltaX = (double*)calloc(n+1,sizeof(double));
  deltaY = (double*)calloc(n+1,sizeof(double));
  A = (double*)calloc((n+1)*(n+1),sizeof(double));
  
  //build deltaX and deltaY
  deltaX[0] = 0;
  deltaX[n] = 0;
  deltaY[0] = 3*(y[1]-y[0]);
  deltaY[n] = 3*(y[n]-y[n-1]);
  for (int k=0; k<n-1; k++){

    deltaX[k+1] = 3*(x[k+2]-x[k]);
    deltaY[k+1] = 3*(y[k+2]-y[k]);

    //fprintf(pFILE, "%f \n", deltaY[k+1]);

  }

  //build Ay
  for (int k=0; k<n; k++){

    A[k+k*(n+1)] += 2.0;
    A[k+1+k*(n+1)] = 1.0;
    A[k+n+1+k*(n+1)] = 1.0;
    A[k+2+n+k*(n+1)] = 2.0;

    //printf("%f %i \n", A[k+2+n+k*(n+1)], k);

  }

  //solve linear system
  dgesv(A, deltaY, n+1);

  //build Ax
  for (int k=0; k<n; k++){

    A[k+k*(n+1)] += 2.0;
    A[k+1+k*(n+1)] = 1.0;
    A[k+n+1+k*(n+1)] = 1.0;
    A[k+2+n+k*(n+1)] = 2.0;

    //printf("%f %i \n", A[k+2+n+k*(n+1)], k);

  }
  A[0] = 1; A[n+1] = 0; A[(n+1)*(n+1)-2-n] = 0; A[(n+1)*(n+1)-1] = 1;
  
  //solve linear system
  dgesv(A, deltaX, n+1);

  //FILE *pFILE;
  //pFILE = fopen ("coeff.txt","w");
  
  //write coefficents
  for (int k=0; k<n; k++){

    //x coeff
    ax[k] = x[k];
    bx[k] = deltaX[k];
    cx[k] = 3*(x[k+1]-x[k])-2*deltaX[k]-deltaX[k+1];
    dx[k] = 2*(x[k]-x[k+1])+deltaX[k]+deltaX[k+1];

    //y coeff
    ay[k] = y[k];
    by[k] = deltaY[k];
    cy[k] = 3*(y[k+1]-y[k])-2*deltaY[k]-deltaY[k+1];
    dy[k] = 2*(y[k]-y[k+1])+deltaY[k]+deltaY[k+1];

    //fprintf(pFILE, "%f %f %f %f %f %f %f %f\n", ax[k], bx[k], cx[k], dx[k], ay[k], by[k], cy[k], dy[k]);

  }

  //fclose(pFILE);

  //free memory
  free(A);
  free(deltaX);
  free(deltaY);
  
  }

//x component on first node of the element
double NormalVectorNxR(double& bx, double& by){

  double nxR = by/sqrt(bx*bx+by*by);

  return nxR;

}

//y component on first node of the element
double NormalVectorNyR(double& bx, double& by){

  double nyR = -bx/sqrt(bx*bx+by*by);

  return nyR;

}

//x component on second node of the element
double NormalVectorNxL(double &bx, double &cx, double &dx, double &by, double &cy, double &dy){

  double nxL = (by+2*cy+3*dy)/sqrt((bx+2*cx+3*dx)*(bx+2*cx+3*dx)+(by+2*cy+3*dy)*(by+2*cy+3*dy));

  return nxL;

}

//y component on second node of the element
double NormalVectorNyL(double &bx, double &cx, double &dx, double &by, double &cy, double &dy){

  double nyL = -(bx+2*cx+3*dx)/sqrt((bx+2*cx+3*dx)*(bx+2*cx+3*dx)+(by+2*cy+3*dy)*(by+2*cy+3*dy));

  return nyL;

}

//curvature on first node of the element
double CurvSplineR(double &bx, double &cx, double &by, double &cy){

  double K = (2*bx*cy-2*by*cx)/pow(bx*bx+by*by,1.5);

  return K;

}

//curvature on second node of the element
double CurvSplineL(double &bx, double &cx, double &dx, double &by, double &cy, double &dy){

  double num = (bx+2*cx+3*dx)*(2*cy+6*dy)-(by+2*cy+3*dy)*(2*cx+6*dx);
  double den = pow((bx+2*cx+3*dx)*(bx+2*cx+3*dx)+(by+2*cy+3*dy)*(by+2*cy+3*dy),1.5);
  
  double K = num/den;

  return K;

}

void dgemv(double *A, double *x, double *y, int m, int n, double alpha , double beta)
{
	int inc = 1;
	char trans = 'N';
	
	dgemv_(&trans, &m, &n, &alpha, A, &m, x, &inc, &beta, y, &inc);
}

void dgemm(double *A, double* B, double *C, int m, int n, int k, double alpha, double beta)
{
	char trans;
	
	trans = 'N';
	
	dgemm_(&trans, &trans, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
}

int dsysv(double *A, double *b, int n)
{
    int nrhs, lda, ldb, info;
	int *ipiv;
	int status = 0;
    char uplo = 'U';
	double *work;
    int lwork = int(n*n/10);  // TEST ceil()
    work = new double[lwork];
	nrhs = 1; /* The Fortran routine that I am calling will solve a
			   system of equations with multiple right hand sides.
			   nrhs gives the number of right hand sides, and then b
			   must be a matrix of dimension n by nrhs. */
	
	lda = n; // The leading dimension of A
	
	ipiv = new int[n]; // Allocate memory to store pivot information
	
	ldb = n;
	
	/* The Fortran routine replaces the input information in b with the
     results.  Since we are interested in the results, we put the
     initial conditions in the output array and use that in the
     call. */
	
    //	for (int i=0; i<n; i++)
	
    //	x[i] = b[i];
	
	// Now the function call
	info = 0;
    
    dsysv_(&uplo, &n, &nrhs, A, &lda, ipiv, b, &ldb, work, &lwork, &info);
	
	// Clean up the memory before returning to the calling program
    assert(info==0);
    assert(b[0]==b[0]);
	delete ipiv;
	return status;
}



int dgesv(double *A, double *b, int n)
{
	int nrhs, lda, ldb, info;
	int *ipiv;
	int status = 0;
	
	nrhs = 1; /* The Fortran routine that I am calling will solve a
			   system of equations with multiple right hand sides.
			   nrhs gives the number of right hand sides, and then b
			   must be a matrix of dimension n by nrhs. */
	
	lda = n; // The leading dimension of A
	
	ipiv = new int[n]; // Allocate memory to store pivot information
	
	ldb = n;
	
	/* The Fortran routine replaces the input information in b with the
	   results.  Since we are interested in the results, we put the
	   initial conditions in the output array and use that in the
	   call. */
	
	//for (int i=0; i<n; i++){ 
	//printf("%f \n", b[i]);
	//}
	
	// Now the function call
	info = 0;
	dgesv_(&n, &nrhs, A, &lda, ipiv, b, &ldb, &info);

	/*FILE* pfile;
	pfile = fopen("results.dat","w");
	for (int i=0; i<n; i++){
	  fprintf(pfile, "%f \n", b[i]);
	}
	fclose(pfile);*/
	
	// Clean up the memory before returning to the calling program
	assert(info==0);
	assert(b[0]==b[0]);
	delete ipiv;
	return status;
}



int dgsev2rhs(double *A, double *b1, double *b2, int n)
{
	int nrhs, lda, ldb, info;
    double *b;
	int *ipiv;
	int status = 0;
	
	nrhs = 2; /* The Fortran routine that I am calling will solve a
			   system of equations with multiple right hand sides.
			   nrhs gives the number of right hand sides, and then b
			   must be a matrix of dimension n by nrhs. */
	
	lda = n; // The leading dimension of A
	
	ipiv = new int[n]; // Allocate memory to store pivot information
	b = new double[2*n];
    
    for (int j=0; j<n; j++) {
        b[j] = b1[j];
        b[j+n] = b2[j];
    }
    
	ldb = n;
	
	/* The Fortran routine replaces the input information in b with the
	   results.  Since we are interested in the results, we put the
	   initial conditions in the output array and use that in the
	   call. */
	
    //	for (int i=0; i<n; i++) 
	
    //	x[i] = b[i];
	
	// Now the function call
	info = 0;
	dgesv_(&n, &nrhs, A, &lda, ipiv, b, &ldb, &info);
	
	// Clean up the memory before returning to the calling program
	assert(info==0);
	assert(b[0]==b[0]);
	for (int j=0; j<n; j++) {
        b1[j] = b[j];
        b2[j] = b[j+n];
    }
    
	delete ipiv;
	return status;
}

int dgels(double *A, double *b, int m, int n, double *x)
{
	int nrhs, lda, ldb, info;

	int lwork;
	int status =0;
	
	char trans = 'N';
	nrhs = 1; /* The Fortran routine that I am calling will solve a
			   system of equations with multiple right hand sides.
			   nrhs gives the number of right hand sides, and then b
			   must be a matrix of dimension n by nrhs. */
	
	lda = m; // The leading dimension of A

	lwork = 2*n*m;
	
	ldb = m;
	
	/* The Fortran routine replaces the input information in b with the
     results.  Since we are interested in the results, we put the
     initial conditions in the output array and use that in the
     call. */
	
	for (int i=0; i<m; i++) x[i] = b[i];
	
	// Now the function call
	
	dgels_(&trans, &m, &n, &nrhs, A, &lda, x, &ldb, b, &lwork, &info);
	
	// Clean up the memory before returning to the calling program
	assert(info==0);

	return status;
}

int dgeinv(double *A, int N){
    int info;
    int *ipiv;
    ipiv = new int[N]; // Allocate memory to store pivot
    double *work;
    int lwork = int(N*N*1.1);  // TEST ceil()
    work = new double[lwork];
    
    dgetrf_(&N, &N, A, &N, ipiv, &info);
    dgetri_(&N, A, &N, ipiv, work, &lwork, &info);
    
    delete work;
    return info;

}

int schur(double *A, double *B, double *C, double *D, double *f1, double *f2, int M, int N){
    double *work, *vork;
    int info;
    work = new double[M*N];
    vork = new double[M];
    
    // Matrix multiplications
    dgemm(A, B, work, M, N, M, 1.0, 0.0);
    dgemm(C, work, D, N, N, M, -1.0, 1.0);
    
    // Vector multiplications
    dgemv(A, f1, vork, M, M, 1.0, 0.0);
    dgemv(C, vork, f2, N, M, -1.0, 1.0);
    
    info = dgesv(D, f2, N);
    
    delete work;
    delete vork;
    return info;
}

int dgesch(double *A, double *f1, int M, int N){
    double *work, *vork;
    int info;
    double *B, *C, *D, *f2;
    work = new double[M*N];
    vork = new double[M];
    C = &A[M*M];
    B = &A[M*(N+M)];
    D = &A[M*(M+2*N)];
    f2 = &f1[M];
    
    // Matrix multiplications
    dgemm(A, B, work, M, N, M, 1.0, 0.0);
    dgemm(C, work, D, N, N, M, -1.0, 1.0);
    
    // Vector multiplications
    dgemv(A, f1, vork, M, M, 1.0, 0.0);
    dgemv(C, vork, f2, N, M, -1.0, 1.0);
    
    info = dgesv(D, f2, N);
    
    delete work;
    delete vork;
    return info;
}


void dg_ctof(double **in, double *lout, int rows, int cols)
{
	int i, j;
	
	for (i=0; i<rows; i++) {
		for (j=0; j<cols; j++) {
			lout[i+j*rows] = in[i][j];
		}
	}
}

