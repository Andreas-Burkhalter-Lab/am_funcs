#include <iostream>
using namespace std;
#include "NeighborhoodStatistics.h"

// Compile with g++ -g -I.. nbrstats.cpp -o nbrstats

int main()
{
  NeighborhoodStatistics<float> nstats(0.001,Moments::Full);

  cout << nstats.isSignificant(5,1,3) << endl;
  cout << nstats.isSignificant(5,2,4) << endl;
  cout << nstats.isSignificant(5,2,5) << endl;
  cout << nstats.isSignificant(5,2,6) << endl;
  cout << nstats.isSignificant(5,2,7) << endl;
  cout << nstats.isSignificant(35,2,6) << endl;

  nstats.print();
}

  
  
