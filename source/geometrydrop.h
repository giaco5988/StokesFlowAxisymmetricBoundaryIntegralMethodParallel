//
//  geometrydrop.h
//
//  Created by Giacomo Gallino on 26.01.16.
//

#ifndef geometrydrop_h
#define geometrydrop_h
#include <time.h>

//#include "bndelementaxis.h"

//#include <string>

class GeometryDrop {

  //these are the main charecteristic of the object
  double *xCoord, *yCoord; //drop coordinate
  double *k1, *k2; //drop curvature
  double *nx, *ny; //drop normal vector
  int DropElem; //number of elements
  double *OldxCoord, *OldyCoord, *nxOld, *nyOld, *k1Old, *k2Old; //useful for time marching
  time_t Tstart;  //time at which the simulation starts
  double V0drop;  //starting volume of the droplet

public:

  //constructor
  GeometryDrop (double, double, double, int);
  //constructor, upload data
  GeometryDrop (int, int, double, double, int);
  //print x and y coordiante to screen
  void PrintCoord ();
  //save x and y coordinates to folder
  void SaveCoord (char*, int);
  //save x and y coordinates k1 k2 nx ny to folder
  void SaveDropData (char*, int);
  //get number of elemnts
  int getElem ();
  //get x coordinates as pointer
  double* GetXCoord();
  //get y coordinates as pointer
  double* GetYCoord();
  //get ONE x coordinates as pointer
  double GetOneXCoord(int);
  //get ONE y coordinates as pointer
  double GetOneYCoord(int);
  //get ONE OLD x coordinates as pointer
  double GetOneOldXCoord(int);
  //get ONE OLD y coordinates as pointer
  double GetOneOldYCoord(int);
  //get x coordinates as pointer
  double* GetOldXCoord();
  //get y coordinates as pointer
  double* GetOldYCoord();
  //get x and y coordinates as pointer
  void GetCoord(double*,double*);
  //set interface coordinates
  void SetCoord(double*, double*);
  //set ONE interface coordinates
  void SetOneCoord(double, double, int);
  //set OLD interface coordinates
  void SetOldCoord(double*, double*);
  //set ONE OLD interface coordinates
  void SetOneOldCoord(double, double, int);
  //get curvature
  double* GetK1();
  double* GetK2();
  //get OLD curvature
  double* GetOldK1();
  double* GetOldK2();
  //set curvature
  void SetCurv(double, double, int);
  //update curvature vector
  void UpdateCurv(double*, double*, double*, double*, double*, double*, double*, double*, int);
  //set OLD curvature
  void SetOldCurv(double, double, int);
  //get normal vector
  double* GetNormalX();
  double* GetNormalY();
  //get OLD normal vector
  double* GetOldNormalX();
  double* GetOldNormalY();
  //set normal vector
  void SetNormal(double, double, int);
  //update normal vector
  void UpdateNormal(double*, double*, double*, double*, double*, double*, double*, double*, int);
  //set OLD normal vector
  void SetOldNormal(double, double, int);
  //compute area an volume of the droplet
  void ComputeAreaVolume(double&,double&);
  //get V0drop
  double GetV0drop();
  //volume correction
  void VolumeCorrection(int);
  //break if the interface escape the domain
  void CheckBreakDomain();
  //monitor if the simulation crash
  void BreakOn(char*, int);
  void BreakOff(char*, int);
  //set value of round.txt
  void SetRound(int, char*);
  //set which is the first round of previous shape which is refined enough                                                                                
  int SetPreviousRound(char*, char*, int);
  //get trough the folder and figure out which is the first that got unstable                                                                                                                                       
  int GetFirstUnstable(int, int, char*, char*);
  //get shapes at a certain time, first stable and first unstable                                                                                                                        
  void GetShapesAtTime(int, int, int, int, char*, char*, double *, double *, double *, double *);
  //weigth average of shapes                                                                                                                                                                                                    
  void AverageShapes(double, int, double*, double*, double*, double*);
  //compute the norm of the perturbation which is equal to the norm of R(theta)-1                                                                                                                                  
  double ComputeDelta(double , int , double *, double *, double *, double*);
  //check if previous round was already fine enough                                                    
  void CheckPreviousRound(int, int, int, char*, char*);
  //check if delta step is under the desired one                                                                                                                                       
  void CheckDeltaStep(double, double, double, char*, char*, int, int, int);
  //checl if ite get shape is too large (simulation already break) and eventually reduce it                                                                                                                          
  int CheckIteGetShape(int, int, int, int, char*, char*);
  //check if drop returned to sphere
  int CheckBreakStable();
  //compute ellipse
  double ComputeEllipse();
  //extrpolate shape
  void ExtrapolateShapes(double , int , double*, double*, double, double, double*, double*);
    
};


#endif
