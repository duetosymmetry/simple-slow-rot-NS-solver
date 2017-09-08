/*
 * ppEOS.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2013 Dec. 1
 *
 * This is the implementation for piecewise polyttropic equations of state
 * following Read, Lackey, Owen, and Friedman (2009) [arXiv:0812.2163]
 *
 * This parametrized EOS has parameters: Gamma1, Gamma2, Gamma3, and p1.
 * All numbers are in CGS (densities in g/cm^3, pressures in dyne/cm^2)
 *
 * The transition densities for the three power laws are
 *   rho12 = 5.01187e14; // 10^14.7 g/cm^3 = 1.85 rho_nuc
 *   rho23 = 1.e15;      // g/cm^3
 *
 * BELOW nuclear densities, we use a FIXED analytic fit that comes
 * from Table II: TODO ERROR IN PAPER
 *
 *   Ki          Γi      ρi
 * ----------- ------- -----------
 * 6.80110e-09 1.58425 2.44034e+07
 * 1.06186e-06 1.28733 3.78358e+11
 * 5.32697e+01 0.62223 2.62780e+12
 * 3.99874e-08 1.35692 –
 *
 */

#include "ppEOS.hpp"
#include <cmath>
#include <string>
#include <sstream>

#include "constants.hpp" /* for G and c */

// Eq. (1)
double ppEOS::pressure( double rho ) const
{
  int i=0;
  while ( (rho >= rho_trans[i]) && (i<6) ) i++;

  return K[i] * pow(rho, Gamma[i]);

};

// Eq. (6)
double ppEOS::epsilon( double rho ) const
{
  double P_g_cm_3 = pressure( rho ) / c2_cm_s;

  int i=0;
  while ( (rho >= rho_trans[i]) && (i<6) ) i++;
  
  return (1. + a[i])*rho + P_g_cm_3 / (Gamma[i] - 1.);

};

double ppEOS::epsilonOfP( double P ) const
{

  int i=0;
  while ( (P >= P_trans[i]) && (i<6) ) i++;
  
  return (1. + a[i]) * pow( P / K[i], 1/Gamma[i] )
    + P / c2_cm_s / (Gamma[i] - 1. );

};

double ppEOS::geomepsilonOfgeomP( double geomP ) const
{
  double f;

  /* convert p from cm^-2 into dyne/cm^2 */
  geomP /= G_cgs/(c_cm_s*c_cm_s*c_cm_s*c_cm_s);

  f = epsilonOfP( geomP );
  /* convert f from g/cm^3 to cm^-2 */
  f *= G_cgs/(c_cm_s*c_cm_s);

  return f;
};

ppEOS::ppEOS( double p1, double Gamma1, double Gamma2, double Gamma3 )
{

  setUpConsts();
  setModel(p1, Gamma1, Gamma2, Gamma3);

};

// extended version
ppEOS::ppEOS( double p1, double Gamma1, double Gamma2, double Gamma3, double rho1 )
{

  setUpConsts();
  setModel(p1, Gamma1, Gamma2, Gamma3, rho1);

};

/* Extended version of setModel -- changes rho_trans[4] (aka rho1) */
void ppEOS::setModel( double p1, double Gamma1, double Gamma2, double Gamma3, double rho1 )
{

  this->p1 = p1;
  this->rho_trans[4] = rho1;
  this->Gamma[4] = Gamma1;
  this->Gamma[5] = Gamma2;
  this->Gamma[6] = Gamma3;

  // Compute the K's via Eq. (11) [arXiv version]
  K[4] = p1 / pow( rho_trans[4], Gamma[4] );
  K[5] = K[4] * pow( rho_trans[4], Gamma[4]-Gamma[5] );
  K[6] = K[5] * pow( rho_trans[5], Gamma[5]-Gamma[6] );

  // The transition from low to high density EOS
  rho_trans[3] = pow( K[4]/K[3] , 1./(Gamma[3]-Gamma[4]) );

  // Compute the transition P's
  for (int i=0; i<6; i++)
    P_trans[i] = pressure(rho_trans[i]);

  // Compute the ai's for epsilon
  a[0]=0.;
  for (int i=1; i<=6; i++) {
    // Eq. (6) [arXiv]
    double eps_i_minus_1 = (1. + a[i-1])*rho_trans[i-1]
      + P_trans[i-1] / c2_cm_s / (Gamma[i-1] - 1.);
    // Eq. (7) [arXiv]
    a[i] = eps_i_minus_1 / rho_trans[i-1] - 1.
      - K[i] / (Gamma[i] - 1.) *
      pow(rho_trans[i-1], Gamma[i] - 1.) / c2_cm_s;
  };

};

void ppEOS::setModel( double p1, double Gamma1, double Gamma2, double Gamma3 )
{
  const double Reads_rho1 = 5.01187e14; // 10^14.7 g/cm^3 = 1.85 rho_nuc

  // Call the more general version of setModel, passing in Read's value for rho1
  setModel(p1, Gamma1, Gamma2, Gamma3, Reads_rho1);
};

void ppEOS::setUpConsts() {

  rho_trans[0] = 2.44034e+07;
  rho_trans[1] = 3.78358e+11;
  rho_trans[2] = 2.62780e+12;
  // rho_trans[3] is determined from p1 and Gamma1 (inputs)
  // rho_trans[4] used to be a constant, but is now set by setModel
  // rho_trans[4] = 5.01187e14; // 10^14.7 g/cm^3 = 1.85 rho_nuc
  rho_trans[5] = 1.e15;

  // There is an error in Read's paper that dropped a factor of c^2.
  K[0] = 6.80110e-09 * c2_cm_s;
  K[1] = 1.06186e-06 * c2_cm_s;
  K[2] = 5.32697e+01 * c2_cm_s;
  K[3] = 3.99874e-08 * c2_cm_s;

  Gamma[0] = 1.58425;
  Gamma[1] = 1.28733;
  Gamma[2] = 0.62223;
  Gamma[3] = 1.35692;

};

std::string ppEOS::summary() const
{
  std::ostringstream o;
  o << "EOS summary: "
    << "rho1=" << rho_trans[4]
    << ", p1=" << p1
    << ", Gamma1=" << Gamma[4]
    << ", Gamma2=" << Gamma[5]
    << ", Gamma3=" << Gamma[6];

  return o.str();
};
