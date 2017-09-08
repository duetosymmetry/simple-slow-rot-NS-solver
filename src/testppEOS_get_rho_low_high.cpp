#include "ppEOS.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>

int main( int argc, char *argv[] ) {

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " Gamma1." << std::endl;
    return 1;
  };

  double Gamma1 = atof(argv[1]);

  ppEOS SLy( pow(10., 34.384), // p1
	     Gamma1,          // Gamma1
	     2.988,           // Gamma2
	     2.851 );        // Gamma3

  std::cout << SLy.get_rho_low_high() << std::endl;

  return 0;

};
