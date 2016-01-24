#include <iostream>
#include <fstream>

#include "CIE_matching_functions.h"

using namespace std;

int NUMBER_OF_WAVELENGTHS = 40;

int main() {
  cout << "Hello, World!" << endl;

  // spectral refelectance distribution
  struct point *srd = new struct point[NUMBER_OF_WAVELENGTHS];
  string line;
  int i = 0;

  ifstream myfile("mysteryS16.spectrum");
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      cout << srd[i].x << srd[i].y << '\n';
      ++i;
    }
    myfile.close();
  }

  else cout << "Unable to open file";

  return 0;
}