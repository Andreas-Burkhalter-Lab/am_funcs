#include <iostream>
using namespace std;
#include "OutlierStatistics.h"

// Compile with g++ -g -I.. outlierstats.cpp -o outlierstats

int main()
{
  OutlierStatistics<double> pstats(0.1);

  cout << pstats.isSignificant(5,1,3) << endl;
  cout << pstats.isSignificant(5,2,4) << endl;
  cout << pstats.isSignificant(5,2,5) << endl;
  cout << pstats.isSignificant(5,2,6) << endl;
  cout << pstats.isSignificant(5,2,7) << endl;
  cout << pstats.isSignificant(35,2,6) << endl;

  pstats.print();
}

  
  
