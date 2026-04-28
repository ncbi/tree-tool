// speed_test.cpp

#include <iostream>
#include <chrono>
#include <cmath>
using namespace std;


int main() 
{
  const auto start = chrono::high_resolution_clock::now ();

  double d = 0.0;
  for (size_t i = 0; i < 100000000; i++) 
    d += sin (i);
  cout << d << endl;

  const auto end = chrono::high_resolution_clock::now ();
  chrono::duration<double> elapsed = end - start;
  cout << "Elapsed time: " << elapsed. count () << " seconds" << endl;
    
  return 0;
}
