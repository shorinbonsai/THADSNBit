/*
 * A bit sprayer is a Moore configuration self driving automata.
 *
 */

#ifndef    _BITSPRAY_H
#define    _BITSPRAY_H

using namespace std;

class bitspray {//self driving automata bit generator class

 public:

  bitspray();           //creates an unallocated bitspray
  bitspray(int S);      //create a bitspray with buffer S states
  bitspray(bitspray &other);  //copy constructor
  ~bitspray();                //destructor

  //utility routines
  void create(int S);         //allocate S states
  void randomize();           //fill in new random values
  void copy(bitspray &other); //copy routine
  void destroy();             //deallocate

  //use routines
  void reset(int *Q, int &H, int &T); //reset to initial state
//  int bit();                  //get the next bit
//  int val(int b);             //return an integer value based on b bits
//  double rel(int b);          //return a [0,1) value based on b bits

  //Genetics
  void tpc(bitspray &other);        //two point crossover
  void mutate(int nm);              //mutate a response or transition

  //Input/Output
  void print(ostream &aus);    //write in human readable form
  void write(ostream &aus);    //write to output stream
  void read(istream &inp);     //read from input stream
//  void getBuf(ostream &out);   //print out the current buffer from hd to tl
//  bool hasnext(int b) const;         //more room in buffer
  void next(int *Q, int &H, int &T, int qz); //add next value(s) to Q

 private:

  //Nb==0 is the unallocated tag
  int init;     //initial output
//  int Nb;       //size of the buffer for output that is not yet input
//  int *buf;     //buffer for output
  int St;       //number of states
  int cs;       //current state
//  int hd, tl;    //pointer to buffered output, head, tail
  int **trans;  //transitions  [n][2]
  int **resp;   //responses     [n][*]Nabcde...   N is length of response

};
#endif /*BITSPRAY.H*/
