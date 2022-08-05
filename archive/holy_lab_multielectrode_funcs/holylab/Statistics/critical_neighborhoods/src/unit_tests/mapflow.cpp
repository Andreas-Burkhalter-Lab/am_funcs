#include <vector>
#include <algorithm>
#include <iostream>
#include "CriticalNeighborhood.h"

// Compile with g++ -g -I.. mapflow.cpp -o mapflow

using namespace std;

int main()
{
  int theMap[] = {1,2,1};
  int theN[] = {5,37,37};
  int len = 3;
  vector<int> mapsTo(theMap,theMap+len);
  vector<int> n(theN,theN+len);

  CriticalNeighborhood::flowMap(mapsTo,n);

  for (int i = 0; i < len; i++)
    cout << mapsTo[i] << ' ';
  cout << endl;
}
