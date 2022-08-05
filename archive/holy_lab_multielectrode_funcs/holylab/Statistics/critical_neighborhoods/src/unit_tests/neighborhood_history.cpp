#include <iostream>
#include "NeighborhoodHistory.h"

using namespace std;

int main()
{
  vector<int> sortOrder;
  int i;
  NeighborhoodHistory history;
  
  for (i = 0; i < 9; i++)
    sortOrder.push_back(i);
  history.add(sortOrder);
  cout << history.isCycling() << ' ' << history.isAtMax() << endl;
  sortOrder.clear();
  for (i = 0; i < 10; i++)
    sortOrder.push_back(i);
  history.add(sortOrder);
  cout << history.isCycling() << ' ' << history.isAtMax() << endl;
  sortOrder.clear();
  for (i = 0; i < 9; i++)
    sortOrder.push_back(i);
  sortOrder[1] = 2;
  sortOrder[2] = 1;
  history.add(sortOrder);
  cout << history.isCycling() << ' ' << history.isAtMax() << endl;
  sortOrder.clear();
  for (i = 0; i < 9; i++)
    sortOrder.push_back(i);
  history.add(sortOrder);
  cout << history.isCycling() << ' ' << history.isAtMax() << endl;
  sortOrder.clear();
  for (i = 0; i < 10; i++)
    sortOrder.push_back(i);
  history.add(sortOrder);
  cout << history.isCycling() << ' ' << history.isAtMax() << endl;
}
  
