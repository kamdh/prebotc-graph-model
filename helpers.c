#include <stdio.h>
#include <math.h>
#include "helpers.h"

double infHN( double A, double B, double V ) {
  return 1.0/(1.0 + exp( (V - A)/B ));
}

double tau( double A, double B, double C, double V ) {
  return A / cosh( 0.5 * (V - B)/C );
}
    
double phi( double x, double kNa ) {
  return pow(x,3) / (pow(x, 3) + pow(kNa, 3));
}
