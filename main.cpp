#include <iostream>
#include <fstream>

#include "CIE_matching_functions.h"
#include "vecmat/Vector.h"
#include "vecmat/Matrix.h"

using namespace std;

int NUMBER_OF_WAVELENGTHS = 40;
struct point *srd;
struct point *spd;
double xw = 0.31271;
double yw = 0.32902;
double xr = 0.690;
double yr = 0.300;
double xg = 0.205;
double yg = 0.715;
double xb = 0.150;
double yb = 0.045;


double interpolate(double wavelength, point* array)
{
  int base_wavelength = (int) wavelength;
  int index = (base_wavelength - 380)/10;
  double result;

  if (base_wavelength == 770)
    result = array[index].y;
  else
    result = array[index+1].y-array[index].y * (wavelength-array[index].y) + array[index].y;

  return result;
}


double integrate(point* array1, point* array2)
{
  double delta_lambda = 0.01;
  double wavelength = 380.0;
  int iterations = (int)(((double)(NUMBER_OF_WAVELENGTHS*10))/delta_lambda);
  double result = 0.0;

  // Calculate N
  for (int i=0; i<iterations; ++i)
  {
    result += interpolate(wavelength, array1) * interpolate(wavelength, array2) * delta_lambda;
    wavelength += delta_lambda;
  }

  return result;
}

void calculateXYZ(double* X, double* Y, double* Z)
{
  double N;

  N = integrate(spd, fy);
  *X = (1.0/N) * integrate(srd, fx);
  *Y = (1.0/N) * integrate(srd, fy);
  *Z = (1.0/N) * integrate(srd, fz);
}

double calculatez(double x, double y)
{
  return (1.0-x-y)/y;
}

void recoverWhite(double* Xw, double* Yw, double* Zw)
{
  *Yw = 1.0;
  *Xw = xw/yw;
  *Zw = calculatez(xw, yw);
}


int main() {
  double X,Y,Z;
  double Xw, Yw, Zw;

  // spectral refelectance distribution
  srd = new struct point[NUMBER_OF_WAVELENGTHS];
  // spectral power distribution
  spd = new struct point[NUMBER_OF_WAVELENGTHS];

  for (int i=0; i<NUMBER_OF_WAVELENGTHS; ++i) {
    spd[i].x = 380.0 + i*10;
    spd[i].y = 1.0;
  }

  ifstream myfile("/home/awbrenn/Documents/workspace/6050/1/mysteryS15.spectrum");
  if (myfile.is_open()) {
    for (int i=0;i<NUMBER_OF_WAVELENGTHS; ++i) {
        myfile >> srd[i].x >> srd[i].y;
    }
    myfile.close();
  }
  else cout << "Unable to open file";

  // Step 0: compute four integrals for X,Y,Z
  calculateXYZ(&X, &Y, &Z);

  // Step 1: recover white Xw, Yw, Zw,
  recoverWhite(&Xw, &Yw, &Zw);

  // Step 2: solve matrix equation for Cr, Cg, Cb
  const double zr = calculatez(xr, yr);
  const double zg = calculatez(xg, yg);
  const double zb = calculatez(xb, yb);
  Vector3d XwYwZw(Xw, Yw, Zw);
  Matrix3x3 xyzrgb(xr, xg, xb, yr, yg, yb, zr, zg, zb);
  Vector3d CrCgCb = xyzrgb.inv() * XwYwZw;

  // Step 3: Solve matrix equation for R,G,B
  Vector3d XYZ(X,Y,Z);
  Matrix3x3 xyzrgbC(xr*CrCgCb[0], xg*CrCgCb[1], xb*CrCgCb[2],
                    yr*CrCgCb[0], yg*CrCgCb[1], yb*CrCgCb[2],
                    zr*CrCgCb[0], zg*CrCgCb[1], zb*CrCgCb[2]);
  Vector3d RGB = xyzrgbC.inv() * XYZ;

  cout << "xyzrgb\n" << xyzrgb << "\n" << endl;
  cout << "X, Y, Z\n" << XYZ << "\n" << endl;
  cout << "Cr, Cg, Cb\n" << CrCgCb << "\n" << endl;
  cout << "RGB\n" << RGB << "\n" << endl;

  return 0;
}