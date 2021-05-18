#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>

using namespace std;

#include "setu.h"
#include "stat.h"

/*************************algorithm controls******************************/
#define PL 16
#define NSE 50
#define alpha 0.5
#define mepl 3              //  Minimum epidemic length
#define rse 5               //  Re-try short epidemics
//#define ftl 50              //  Final test length
#define verbose true
#define runs 30
#define mevs 40000
#define RI ((long)mevs/100)
#define NmC (long)9
#define EDGB 2              //  Minimum degree for swap
#define popsize 1000
#define verts 160
#define GL 256
#define RNS 91207819
#define MAXL (long)pow(verts, 3)  //  Size of integer
#define MNM 3
#define tsize 7
#define omega 0.5
#define edgeAdd 1
#define triProb 0.5

/**************************Variable dictionary************************/
int pop[popsize][GL];           //  Population of command strings
double fit[popsize];            //  Fitness array
int dx[popsize];                //  Sorting index
double PD[PL + 1];               //  Profile dictionary
double CmD[NmC];                //  Command densities
int mode;                       //  Profile?
bool ringG;
graph iG;

/****************************Procedures********************************/
void initalg(const char *pLoc);         //initialize the algorithm
int validloci();                       //generate an acceptable large integer
void express(graph &G, const int *cmd);       //express a command string
void initpop();                        //initialize population
void matingevent();                    //run a mating event
void report(ostream &aus);             //make a statistical report
void reportbest(ostream &aus, ostream &difc);
void createReadMe(ostream &aus);
void cmdLineIntro(ostream &aus);
void cmdLineRun(int run, ostream &aus);

/****************************Main Routine*******************************/
int main(int argc, char *argv[]) {
  /**
   * Output Location, Profile Location, Command Densities
   */

  mode = 1;   //  0 - Epidemic Length, 1 - Profile Matching
  ringG = false;
  /*
   * Mode 0 -> Epidemic Length (w Densities)
   * Mode 1 -> Profile Matching (w Densities)
   * Mode 2 -> Profile Matching (w Bitsprayers)
   */

  fstream stat, best, dchar, readme, iGOut;   //statistics, best structures
  char fn[60];          //file name construction buffer
  char *outLoc = argv[1];
  char *pLoc = argv[2];

  initalg(pLoc);
  if (mode < 2) { // Densities
    int offset = 3;
    for (int cmd = 0; cmd < NmC; cmd++) {
      CmD[cmd] = strtod(argv[cmd + offset], nullptr);
    }
  }

  sprintf(fn, "%sbest.lint", outLoc);
  best.open(fn, ios::out);
  sprintf(fn, "%sdifc.dat", outLoc);
  dchar.open(fn, ios::out);
  sprintf(fn, "%sreadme.dat", outLoc);
  readme.open(fn, ios::out);
  createReadMe(readme);
  readme.close();
  if (!ringG) {
    sprintf(fn, "%sinitGraph.dat", outLoc);
    iGOut.open(fn, ios::out);
    iG.write(iGOut);
    iGOut.close();
  }
  if (verbose) {
    cmdLineIntro(cout);
  }

  for (int run = 0; run < runs; run++) {
    if (verbose) cmdLineRun(run, cout);
    sprintf(fn, "%srun%02d.dat", outLoc, run); // File name
    stat.open(fn, ios::out);
    initpop();
    report(stat); //report the statistics of the initial population
    for (int mev = 0; mev < mevs; mev++) {//do mating events
      matingevent();  //run a mating event
      if ((mev + 1) % RI == 0) {//Time for a report?
        if (verbose) {
          cout << left << setw(5) << run;
          cout << left << setw(4) << (mev + 1) / RI;
        }
        report(stat); //report statistics
      }
    }
    stat.close();
    reportbest(best, dchar);
    cout << "Done run " << run << " of " << runs - 1 << endl;
  }
  best.close();
  dchar.close();
  return (0);  //keep the system happy
}

void createReadMe(ostream &aus) {
  aus << "This file contains the info about the files in this folder." << endl;
  aus << "Graph Evolution Tool." << endl;
  aus << "Initial Graph: ";
  if (ringG) {
    aus << "Ring with m = 2" << endl;
  } else {
    aus << "Power Law Cluster Graph" << endl;
  }
  switch (mode) {
    case 0:
      aus << "Fitness Function: Epidemic Length" << endl;
      aus << "Representation: THADS-N with Operation Densities" << endl;
      aus << "Gene length: " << GL << endl;
      aus << "Densities: ";
      for (double i : CmD) {
        aus << i << " ";
      }
      aus << endl;
      break;
    case 1:
      aus << "Fitness Function: Profile Matching" << endl;
      aus << "Representation: THADS-N with Operation Densities" << endl;
      aus << "Gene length: " << GL << endl;
      aus << "Densities: ";
      for (double i : CmD) {
        aus << i << " ";
      }
      aus << endl;
      aus << "Profile: ";
      for (double i : PD) {
        aus << i << " ";
      }
      aus << endl;
      break;
//    case 2:
//      //TODO: complete
//      break;
    default:
      break;
  }
  aus << endl;
  aus << "The parameter settings are as follows: " << endl;
  aus << "Number of sample epidemics: " << NSE << endl;
  aus << "Alpha: " << alpha << endl;
  aus << "Minimum epidemic length: " << mepl << endl;
  aus << "Re-tries for short epidemics: " << rse << endl;
  aus << "Runs: " << runs << endl;
  aus << "Mating events: " << mevs << endl;
  aus << "Minimum degree for swap: " << EDGB << endl;
  aus << "Population size: " << popsize << endl;
  aus << "Number of vertices: " << verts << endl;
  aus << "Maximum number of mutations: " << MNM << endl;
  aus << "Tournament size: " << tsize << endl;
  aus << "Decay strength for diffusion characters: " << omega << endl;
  if (!ringG) {
    aus << "Edges added at each step: " << edgeAdd << endl;
    aus << "Triangle creation probability: " << triProb << endl;
  }
  aus << endl;
  aus << "The file descriptions are as follows: " << endl;
  aus << "best.lint -> the best fitness and it's associated gene for each run";
  aus << endl;
  aus << "difc.dat -> the diffusion characters of the best graph for each run";
  aus << endl;
  aus << "run##.dat -> population statistics for each run" << endl;
  if (!ringG) {
    aus << "initGraph.dat -> the initial power law clustering graph" << endl;
  }
}

void cmdLineIntro(ostream &aus) {
  aus << "Graph Evolution Tool." << endl;
  aus << "Initial Graph: ";
  if (ringG) {
    aus << "Ring with m = 2" << endl;
  } else {
    aus << "Power Law Cluster Graph" << endl;
  }
  switch (mode) {
    case 0:
      aus << "Fitness Function: Epidemic Length" << endl;
      aus << "Representation: THADS-N with Operation Densities" << endl;
      aus << "Densities: ";
      for (double i : CmD) {
        aus << i << " ";
      }
      aus << endl;
      break;
    case 1:
      aus << "Fitness Function: Profile Matching" << endl;
      aus << "Representation: THADS-N with Operation Densities" << endl;
      aus << "Densities: ";
      for (double i : CmD) {
        aus << i << " ";
      }
      aus << endl;
      aus << "Profile: ";
      for (double i : PD) {
        aus << i << " ";
      }
      aus << endl;
      break;
//    case 2:
//      //TODO: complete
//      break;
    default:
      break;
  }
  aus << "Check readme.dat for more information about parameters/output.";
  aus << endl;
}

void cmdLineRun(int run, ostream &aus) {
  aus << endl << "Beginning Run " << run << " of " << runs - 1 << endl;
  aus << left << setw(5) << "Run";
  aus << left << setw(4) << "RI";
  aus << left << setw(10) << "Mean";
  aus << left << setw(12) << "95% CI";
  aus << left << setw(10) << "SD";
  aus << left << setw(8) << "Best";
  aus << endl;
  aus << left << setw(5) << run;
  aus << left << setw(4) << "0";
}

void initalg(const char *pLoc) {//initialize the algorithm
  fstream inp;    //input file
  char buf[20];   //input buffer

  srand48(RNS);                 //read the random number seed
  if (!ringG) {
    iG.create(verts);
    iG.PCG(verts, edgeAdd, triProb);
  }
  if (mode == 1) {
    inp.open(pLoc, ios::in);      //open input file
    for (int i = 0; i < PL; i++) {
      PD[i] = 0;  //pre-fill missing values
    }
    PD[0] = 1;  //put in patient zero
    for (int i = 0; i < PL; i++) {      //loop over input values
      inp.getline(buf, 19);  //read in the number
      PD[i + 1] = strtod(buf, nullptr);    //translate the number
    }
    cout << "Profile: "; //Prints out the profile being tested
    for (int i = 0; i <= PL; i++) {
      cout << PD[i] << " ";
    }
    cout << endl;
    inp.close();
  }
}

//This routine generates valid loci for the expression routine
int validloci() {//generate an acceptable large integer
  // TODO: Insert SDA
  int cmd;       //command type generated
  double dart;   //Random command

  dart = drand48() - CmD[0];//throw the dart
  cmd = 0;                  //initialize the position on the dartboard
  while ((dart > 0) && (cmd < NmC - 1)) {
    dart -= CmD[++cmd]; //walk the board
  }
  cmd += (int) (NmC * (lrand48() % MAXL)); //add in the large integer part
  return (cmd);  //return the generated command
}

//This is expression of a large integer represenation
void express(graph &G, const int *cmd) {//express a command string
  int a, b, c;  //decoded values
  int cdv;    //command value
  int block;  //integer carving block

  if (ringG) {
    G.RNGnm(verts, 2);  //  Initial graph
  } else {
    G.copy(iG);
  }

  for (int i = 0; i < GL; i++) {//loop over the commands (genetic loci)
    block = cmd[i];  //get integer
    cdv = (int) (block % NmC); //slice of the command
    block /= NmC;    //clear command information from block
    switch (cdv) {  //What command is it?
      case 0://Toggle
        a = block % verts;
        b = (block / verts) % verts; //get vertex numbers
        G.toggle(a, b);    //toggle edge {a,b}
        break;
      case 1://Hop
        //get vertex numbers
        a = block % verts;
        b = (block / verts) % verts;
        c = (block / verts / verts) % verts;
        G.hop(a, b, c);
        break;
      case 2://Add
        a = block % verts;
        b = (block / verts) % verts; //get vertex numbers
        G.add(a, b);       //add edge {a,b}
        break;
      case 3://Delete
        a = block % verts;
        b = (block / verts) % verts; //get vertex numbers
        G.del(a, b);       //delete edge a,b
        break;
      case 4://Swap
        a = block % (verts * 10);
        b = (block / (verts * 10)) % (verts * 10); //get vertex numbers
        G.edgeswap(a, b, EDGB);
        break;
      case 5://Local Toggle
        //get vertex numbers
        a = block % verts;
        b = (block / verts) % verts;
        c = (block / verts / verts) % verts;
        G.loggle(a, b, c);
        break;
      case 6: // Local Add
        //get vertex numbers
        a = block % verts;
        b = (block / verts) % verts;
        c = (block / verts / verts) % verts;
        G.ladd(a, b, c);
        break;
      case 7: // Local Del
        //get vertex numbers
        a = block % verts;
        b = (block / verts) % verts;
        c = (block / verts / verts) % verts;
        G.ldel(a, b, c);
        break;
      case 8: // Null - Do nothing
        break;
      default:
        cout << "ERROR: default of switch";
        break;
    }
  }
}

double fitness(int *cmd) {//compute the epidemic length fitness
  graph G(verts);      //scratch graph
  int max, len, ttl;   //maximum, length, and total removed
  int cnt;             //counter for tries
  double prof[verts];  //profile variable
  double trials[NSE];  //stores squared error for each trial
  int en;              //epidemic number
  double delta;        //difference between profile and trial
  double accu = 0.0;         //accumulator
  express(G, cmd);               //create the graph

  if (mode == 0) {
    for (en = 0; en < NSE; en++) {
      cnt = 0;
      do {
        G.SIR(0, max, len, ttl, alpha);
        cnt++;
      } while (len < mepl && cnt < rse);
      trials[en] = len;
    }
    for (double trial : trials) {//loop over trials
      accu += trial;
    }
    accu = accu / NSE;
  } else if (mode == 1) {
    for (en = 0; en < NSE; en++) {//loop over epidemics
      cnt = 0;
      do {
        G.SIRProfile(0, max, len, ttl, alpha, prof);
        cnt++;
      } while (len < mepl && cnt < rse);
      trials[en] = 0;  //zero the current squared error
      if (len < PL + 1) {
        len = PL + 1;  //find length of epi/prof (longer)
      }
      for (int i = 0; i < len; i++) {//loop over time periods
        delta = prof[i] - PD[i];
        trials[en] += delta * delta;
      }
      trials[en] = sqrt(trials[en] / len); //convert to RMS error
    }

    for (double trial : trials) {
      accu += trial;
    }
    accu = accu / NSE;
  }
  return accu;  //return the fitness value
}

void initpop() {//initialize population
  for (int i = 0; i < popsize; i++) {    //loop over the population
    /***INITIALIZATION CODE***/
    for (int j = 0; j < GL; j++) {
      pop[i][j] = validloci(); //fill in the loci
    }
    fit[i] = fitness(pop[i]);  //compute its fitness
    dx[i] = i;                 //refresh the sorting index
  }
}

void matingevent() {//run a mating event
  int rp, sw;   //loop index, random position, swap variable
  int cp1, cp2;   //crossover points

  //perform tournament selection, highest fitness first
  if (mode == 0) {
    tselect(fit, dx, tsize, popsize);
  } else {
    Tselect(fit, dx, tsize, popsize);
  }

  //selection and crossover
  cp1 = (int) lrand48() % GL;
  cp2 = (int) lrand48() % GL;
  if (cp1 > cp2) {
    sw = cp1;
    cp1 = cp2;
    cp2 = sw;
  }
  for (int i = 0; i < cp1; i++) {
    pop[dx[0]][i] = pop[dx[tsize - 2]][i];
    pop[dx[1]][i] = pop[dx[tsize - 1]][i];
  }
  for (int i = cp1; i < cp2; i++) {
    pop[dx[0]][i] = pop[dx[tsize - 1]][i];
    pop[dx[1]][i] = pop[dx[tsize - 2]][i];
  }
  for (int i = cp2; i < GL; i++) {
    pop[dx[0]][i] = pop[dx[tsize - 2]][i];
    pop[dx[1]][i] = pop[dx[tsize - 1]][i];
  }

  //mutation
  rp = (int) lrand48() % MNM + 1;
  for (int i = 0; i < rp; i++) {
    pop[dx[0]][lrand48() % GL] = validloci();
  }
  rp = (int) lrand48() % MNM + 1;
  for (int i = 0; i < rp; i++) {
    pop[dx[1]][lrand48() % GL] = validloci();
  }

  //update fitness
  fit[dx[0]] = fitness(pop[dx[0]]);
  fit[dx[1]] = fitness(pop[dx[1]]);

  // Skeptical tournament selection
  if (mode == 0) {
    fit[dx[tsize - 1]] = fitness(pop[dx[tsize - 1]]);
    fit[dx[tsize - 2]] = fitness(pop[dx[tsize - 2]]);
  }
}

void report(ostream &aus) {//make a statistical report
  dset D;
  D.add(fit, popsize);  //load fitness

  //print report
  if (mode == 0) {
    aus << left << setw(10) << D.Rmu();
    aus << left << setw(12) << D.RCI95();
    aus << left << setw(10) << D.Rsg();
    aus << left << setw(8) << D.Rmax() << endl;
    if (verbose) {
      cout << left << setw(10) << D.Rmu();
      cout << left << setw(12) << D.RCI95();
      cout << left << setw(10) << D.Rsg();
      cout << left << setw(8) << D.Rmax() << endl;
    }
  } else if (mode == 1) {
    aus << left << setw(10) << D.Rmu();
    aus << left << setw(12) << D.RCI95();
    aus << left << setw(10) << D.Rsg();
    aus << left << setw(8) << D.Rmin() << endl;
    if (verbose) {
      cout << left << setw(10) << D.Rmu();
      cout << left << setw(12) << D.RCI95();
      cout << left << setw(10) << D.Rsg();
      cout << left << setw(8) << D.Rmin() << endl;
    }
  }
}

void reportbest(ostream &aus, ostream &difc) {//report the best graph
  int b;       //loop indices and best pointer
  graph G(verts);  //scratch graph
  double En;
  static double M[verts][verts];
  static double Ent[verts];
//  fstream initGOut;    //graph output file
//  char fn[60];     //file name construction buffer

  b = 0;
  if (mode == 0) {
    for (int i = 1; i < popsize; i++) {
      if (fit[i] > fit[b]) {
        b = i; //find best fitness
      }
    }
  } else {
    for (int i = 1; i < popsize; i++) {
      if (fit[i] < fit[b]) {
        b = i; //find best fitness
      }
    }
  }

  //output the fitness and the gene
  aus << fit[b] << " -fitness" << endl;
  // TODO: Write the SDA
  // Write the gene
  aus << "Gene" << endl;
  aus << pop[b][0];
  for (int i = 1; i < GL; i++) {
    aus << " " << pop[b][i];
  }
  aus << endl;

  G.empty(verts);
  express(G, pop[b]);
  aus << "Graph" << endl;
//  sprintf(fn, "%sGraph%d.dat", outLoc, run);
//  gout.open(fn, ios::out);
  G.write(aus);
  for (int i = 0; i < G.size(); i++) {
    G.DiffChar(i, omega, M[i]);
  }

  for (int i = 0; i < G.size(); i++) {//loop over vertices
    En = 0.0;  //prepare En for normalizing
    for (int j = 0; j < G.size(); j++) {
      En += M[i][j];  //total the amount of gas at vertex i
    }
    for (int j = 0; j < G.size(); j++) {
      M[i][j] /= En;  //normalize so sum(gas)=1
    }
    En = 0.0;  //zero the entropy accumulator
    for (int j = 0; j < G.size(); j++) {//build up the individual entropy terms
      if (M[i][j] > 0) {
        En += -M[i][j] * log(M[i][j]);  //this is entropy base E
      }
    }
    Ent[i] = En / log(2);  //convert entropy to Log base 2
  }

  //Now sort the entropy vector
  bool more = false;  //no swaps
  do {//swap out-of-order entries
    more = false;
    for (int i = 0; i < G.size() - 1; i++) {
      if (Ent[i] < Ent[i + 1]) {
        En = Ent[i];
        Ent[i] = Ent[i + 1];
        Ent[i + 1] = En;  //swap
        more = true;  //set the flag that a swap happened
      }
    }
  } while (more); //until in order

  difc << Ent[0];  //output first entropy value
  for (int i = 1; i < G.size(); i++) {
    difc << " " << Ent[i];  //output remaining values
  }
  difc << endl;  //end line
}

