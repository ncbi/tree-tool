// speed_test.cpp

#include <iostream>
#include <chrono>
#include <cmath>
using namespace std;

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    // --- The code you want to test goes here ---
    double d = 0.0;
    for(int i = 0; i < 100000000; ++i) 
      d += std::sin(i);
    cout << d << endl;

    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    return 0;
}
