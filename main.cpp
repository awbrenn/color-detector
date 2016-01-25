#include <iostream>
#include <fstream>

#include "CIE_matching_functions.h"

using namespace std;

int NUMBER_OF_WAVELENGTHS = 40;
struct point *srd;
struct point *spd;


float interpolate(float wavelength, point* array)
{
  int base_wavelength = (int) wavelength;
  int index = (base_wavelength - 380)/10;
  float result;

  if (base_wavelength == 770)
    result = array[index].y;
  else
    result = (array[index+1].y-array[index].y) * (wavelength-array[index].y) + array[index].y;

  return result;
}


float integrate(point* array1, point* array2)
{
  float delta_lambda = 0.1;
  float wavelength = 380.0;
  int iterations = (int)(((float)(NUMBER_OF_WAVELENGTHS*10))/delta_lambda);
  float result = 0.0;

  // Calculate N
  for (int i=0; i<iterations; ++i)
  {
    result += interpolate(wavelength, array1) * interpolate(wavelength, array2) * delta_lambda;
    wavelength += delta_lambda;
  }
}

void calculateXYZ(float* X, float* Y, float* Z)
{
  float N;

  N = integrate(spd, fy);
  *X = (1.0/N) * integrate(srd, fx);
  *Y = (1.0/N) * integrate(srd, fy);
  *Z = (1.0/N) * integrate(srd, fz);
}

int main() {
  float X,Y,Z;

  // spectral refelectance distribution
  srd = new struct point[NUMBER_OF_WAVELENGTHS];
  // spectral power distribution
  spd = new struct point[NUMBER_OF_WAVELENGTHS];

  for (int i=0; i<NUMBER_OF_WAVELENGTHS; ++i) {
    spd[i].x = 380.0 + i*10;
    spd[i].y = 1.0;
  }
  

  ifstream myfile("mysteryS16.spectrum");
  if (myfile.is_open()) {
    for (int i=0;i<NUMBER_OF_WAVELENGTHS; ++i) {
        myfile >> srd[i].x >> srd[i].y;
    }
    myfile.close();
  }
  else cout << "Unable to open file";

  cout << srd[10].x << " " << srd[10].y << endl;

  calculateXYZ(&X, &Y, &Z);

  cout << X << " " << Y << " " << Z << endl;

  return 0;
}