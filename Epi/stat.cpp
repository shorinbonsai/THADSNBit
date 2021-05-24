/*
 *  Class for doing routine statistical tasks. 
 *
 *
 */

#include <cmath>
#include <cstdlib>

#include "stat.h"

//read in a file with N rows and M items per row
double **RRfuleA(int N, int M, const char *fn) {//with allocation

  double **mt;  //matrix
  int i;        //loop index

  mt = new double *[N];                   //allocate the matrix spine
  for (i = 0; i < N; i++)mt[i] = new double[M]; //allocate the rows
  RRfuleNA(N, M, mt, fn);  //call the reader
  return (mt);  //return the allocated matrix

}

void RRfuleNA(int N, int M, double **mt, const char *fn) {//without allocation

  fstream inp;    //input file
  int i, j, k;      //loop index variables
  char buf[1000]; //input buffer

  inp.open(fn, ios::in);  //open file
  for (i = 0; i < N; i++) {//loop over the rows of the file
    inp.getline(buf, 999); //get the line
    mt[i][0] = atof(buf);  //first entry of line
    k = 0;  //intialize the pointer
    for (j = 1; j < M; j++) {//loop over the rest of the row
      while (buf[k] != ' ')k++;
      while (buf[k] == ' ')k++; //find next entry
      mt[i][j] = atof(buf + k); //translate it
    }
  }
  inp.close();  //close file
}


//select a tournament in the first tz spaces of dx LOWEST FITNESS FIRST
void tselect(double *fit,  //fitness
             int *dx,      //sorting index
             int tz,       //tournament size
             int pz) {//population size

  int i, rp, sw;

  for (i = 0; i < tz; i++) {//select tournament members
    rp = lrand48() % pz; //select random member
    sw = dx[i];
    dx[i] = dx[rp];
    dx[rp] = sw; //swap it in
  }

  do {//sort the tournament into lowest fitness first order
    rp = 0;  //reset flag
    for (i = 0; i < tz - 1; i++)
      if (fit[dx[i]] > fit[dx[i + 1]]) {//if out of order
        sw = dx[i];
        dx[i] = dx[i + 1];
        dx[i + 1] = sw; //swap
        rp = 1;  //set flag
      }
  } while (rp == 1); //while flag set

}

//select a tournament in the first tz spaces of dx  HIGHEST FITNESS FIRST
void Tselect(double *fit,  //fitness
             int *dx,      //sorting index
             int tz,       //tournament size
             int pz) {//population size

  int i, rp, sw;

  for (i = 0; i < tz; i++) {//select tournament members
    rp = lrand48() % pz; //select random member
    sw = dx[i];
    dx[i] = dx[rp];
    dx[rp] = sw; //swap it in
  }

  do {//sort the tournament into highest fitness first order
    rp = 0;  //reset flag
    for (i = 0; i < tz - 1; i++)
      if (fit[dx[i]] < fit[dx[i + 1]]) {//if out of order
        sw = dx[i];
        dx[i] = dx[i + 1];
        dx[i + 1] = sw; //swap
        rp = 1;  //set flag
      }
  } while (rp == 1); //while flag set

}

//This routine assumes that the index is already filled with 1..N in some order
void shuffleDX(int *dx, int N) {//shuffle a sorting index

  int i, rp, sw;

  for (i = 0; i < N; i++) {//select tournament members
    rp = lrand48() % N; //select random member
    sw = dx[i];
    dx[i] = dx[rp];
    dx[rp] = sw; //swap it in
  }

}

void smallfirst(double *fit, int *dx, int N) {//special size two tournaments

  int i, sw;

  for (i = 0; i < N; i += 2) {//loop over pairs
    if (fit[dx[i]] > fit[dx[i + 1]]) {//if the pair is out of order
      sw = dx[i];
      dx[i] = dx[i + 1];
      dx[i + 1] = sw;  //put them in order
    }
  }
}


//fitness proportional selection  - 
//if ttl=0 it is computer otherwise user supplied
int
FPS(double *fit, double ttl, int n) {//perform fitness proportional selection

  int i, p;
  double dart;

  if (ttl == 0)for (i = 0; i < n; i++)ttl += fit[i];
  p = 0;
  dart = drand48() * ttl - fit[0];
  while (p < n) {
    if (dart < 0)return (p);
    dart -= fit[++p];
  }
  return (n - 1);
}

double uniform(double z) {//uniform mutation of size z

  return (2 * z * drand48() - z);

}

double Gauss(double z) {//Gaussian with SD z

  return (z * sqrt(-2 * log(drand48())) * cos(2 * M_PI * drand48()));

}

//compute the entropy of a sequence of integers
double Ientropy(int *dt, int n) {

  double val, ttl; //accumulators for n Log(n) and the total
  double ret;     //return value
  int i;          //loop index

  val = ttl = 0.0; //zero the accumulators
  for (i = 0; i < n; i++) {//loop over data
    if (dt[i] > 0) {
      ttl += dt[i];
      val += dt[i] * log(dt[i]);
    }
  }

  ret = log(ttl) - val / ttl;
  ret /= log(2.0);

  return (ret);

}

//constructors and destructors
dset::dset() {//create an empty data set

  clear();

}

dset::~dset() {//clear the data set

  //nothing to do

}

//input methods
void dset::add(double v) {//Add one number

  if (set == 0) {
    set = 1;
    max = min = v;
  }
  if (v > max)max = v;
  if (v < min)min = v;
  sum += v;
  sums += (v * v);
  n++;
  flag = 0;

}

void dset::add(const double v[], int k) {//Add k numbers

  int i;

  if (set == 0) {
    set = 1;
    max = min = v[0];
  }
  for (i = 0; i < k; i++) {
    sum += v[i];
    sums += (v[i] * v[i]);
    if (v[i] > max)max = v[i];
    if (v[i] < min)min = v[i];
  }
  n += k;
  flag = 0;

}

void dset::add(int *v, int k) {//Add k numbers

  int i;
  double cv;

  if (set == 0) {
    set = 1;
    max = min = ((double) v[0]);
  }
  for (i = 0; i < k; i++) {
    cv = ((double) v[i]);
    sum += cv;
    sums += (cv * cv);
    if (cv > max)max = cv;
    if (cv < min)min = cv;
  }
  n += k;
  flag = 0;

}

//use methods
void dset::compute() {//run the statistical computations

  double sc;

  if (n == 0)return;
  sc = ((double) n);
  mu = sum / sc;
  if (n > 1) {
    sg = sums / sc - mu * mu;
    if (sg > 0)sg = sqrt(sg); else sg = 0;
    CI95 = 1.96 * sg / sqrt(sc - 1.0);
  } else sg = CI95 = 0.0;
  flag = 1;

}

void dset::clear() {//clear out the local variables

  set = n = 0;
  sum = sums = 0.0;
  flag = 0;
  mu = sg = max = min = CI95 = 0.0;

}

//reference methods
double dset::Rmu() {//report the mean

  if (flag == 0)compute();
  return (mu);
}

double dset::Rsg() {//report the standard deviation

  if (flag == 0)compute();
  return (sg);

}

double dset::Rmax() {//report the maximum

  if (flag == 0)compute();
  return (max);

}

double dset::Rmin() {//report the minimum

  if (flag == 0)compute();
  return (min);

}

double dset::RCI95() {//report a 95% confidence interval

  if (flag == 0)compute();
  return (CI95);

}
 
