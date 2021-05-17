#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <ctime>

using namespace std;

#include "setu.h"
#include "./BitSprayer/bitspray.h"

#define Z 128

int main() {

  graph G(Z), H(Z);
  double g[Z], h[Z], q[Z];
  int i, j;
  fstream inp;
  int *D;
  int sz;

  inp.open("Graph0.grf", ios::in);
  G.read(inp);
  inp.close();

  sz = G.size();
  cout << sz << endl;
  D = new int[sz];
  G.dfrom(0, D);
  for (i = 0; i < sz; i++)cout << D[i] << " ";
  cout << endl;


  delete[] D;

  //G.write(cout);

}
 
