/*
 *  Class for doing routine statistical tasks. 
 *
 *
 */

#ifndef    _STAT_H
#define    _STAT_H

#include "iostream"
#include "fstream"

using namespace std;

//read in a file with N rows and M items per row
double **RRfuleA(int N, int M, const char *fn); //with allocation
void RRfuleNA(int N, int M, double **mt, const char *fn); //without allocation


//select a tournament in the first tz spaces of dx; low fitness first
void tselect(double *fit,  //fitness
             int *dx,      //sorting index
             int tz,       //tournament size
             int pz);      //population size

//select a tournament in the first tz spaces of dx; highest fitness first
void Tselect(double *fit,  //fitness
             int *dx,      //sorting index
             int tz,       //tournament size
             int pz);      //population size

void shuffleDX(int *dx, int N); //shuffle a sorting index

void smallfirst(double *fit, int *dx, int N); //special size two tournaments

//fitness proportional selection  - 
//if ttl=0 it is computer otherwise user supplied
int
FPS(double *fit, double ttl, int n); //perform fitness proportional selection

double uniform(double z);  //uniform mutation of size z
double Gauss(double z);    //Gaussian with SD z

//compute the entropy of a sequence of integers
double Ientropy(int *dt, int n);

class dset {//data set manipulation class

 public:

  //constructors and destructors
  dset();    //create an empty data set
  ~dset();   //clear the data set

  //input methods
  void add(double v);
  void add(const double *v, int k);
  void add(int *v, int k);

  //use methods
  void compute();
  void clear();

  //reference methods
  double Rmu();
  double Rsg();
  double Rmax();
  double Rmin();
  double RCI95();

 private:

  int flag;
  double sum;
  double sums;
  int n;
  int set;
  double mu;
  double sg;
  double max;
  double min;
  double CI95;

};

#endif /* _STAT_H */
