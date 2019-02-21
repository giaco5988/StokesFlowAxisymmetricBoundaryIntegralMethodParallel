/*
 *  remesh.h
 *
 *  Created by Giacomo Gallino on 28.01.16.
 *  
 */
#ifndef REMESH_H
#define REMESH_H
#include "geometrydrop.h"

class Remesh {

  //mesh management parameters
  double g; //adaptivity coefficient
  GeometryDrop* GeoHere;
  int MESHkind; //type of remeshing
  double *unew;
  double *vnew;
  int meshRK; //type of time marching
  int OptionCluster; //option for clustering
  int OptionDistribution; //option for choosing distribution
  double TuneDistributionAdaptivity; //set sharpeness of the distribution that is used for remesh of mesh stabilization
  
public:

  //constructor
  Remesh(GeometryDrop*,int);
  //do remesh
  void ExecuteRemesh(GeometryDrop*, double*);
  //passive mesh stabilization
  void PassiveMeshStabilization(GeometryDrop*, double*);
  //get output velocity
  double* GetVelU();
  double* GetVelV();
  //set output velocity
  void SetVelU(double,int);
  void SetVelV(double,int);
  //get which kind of remeshing
  int GetMeshKind();
  //get ONE output velocity
  double GetOneVelU(int);
  double GetOneVelV(int);
  //decide adaptivity for drop nodes
  void SetAdaptDrop(double);
  //set type of time marching
  void SetRK(int);
  //compute F for mesh stab
  double computeF(double *, double *, double *, double *, double, int);
  //set option for node clustering
  void SetOptionCluster(int);
  //set option for choosing distribution                                                                                                                                                                                                  
  void SetOptionDistribution(int);
  //remsh following distribution
  void RemeshDistribution(GeometryDrop*, double*, double*);
  //set tune distribution
  void SetTuneDistributionAdaptivity(double);
    
    
};
#endif
