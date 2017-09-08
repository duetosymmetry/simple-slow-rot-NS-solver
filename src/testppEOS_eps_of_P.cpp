#include "ppEOS.hpp"
#include <cmath>
#include <iostream>

int main() {

  const unsigned int numSamples = 400;
  const double P_min = 1.e19;
  const double P_max = 1.e38;

  ppEOS SLy( pow(10., 34.384), // p1
	     3.005,           // Gamma1
	     2.988,           // Gamma2
	     2.851 );        // Gamma3

  // This will be log spaced
  const double log10_P_min = log10(P_min);
  const double log10_P_max = log10(P_max);
  const double delta = (log10_P_max-log10_P_min)/(numSamples-1);

  for (unsigned int i=0; i<numSamples; i++) {
    double P = pow( 10., log10_P_min + i*delta );
    std::cout << "{"
              << P << ", "
              << SLy.epsilonOfP(P)
              << "}," << std::endl;
  };

  return 0;

};
