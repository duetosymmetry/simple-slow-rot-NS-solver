/*
 * ppEOS.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2013 Dec. 1
 *
 * This is the interface for piecewise polyttropic equations of state
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

#pragma once

#include <string>

struct ppEOS {

  /* constructor -- just calls setModel. */
  ppEOS( double p1, double Gamma1, double Gamma2, double Gamma3 );

  /* Extended constructor -- allows rho_trans[4] (aka rho1) to be set.
     rho1 is in units of g/cm^3
     Calls extended version of setModel.
   */
  ppEOS( double p1, double Gamma1, double Gamma2, double Gamma3, double rho1 );

  /* to assign a new model, if you want that for some reason */
  void setModel( double p1, double Gamma1, double Gamma2, double Gamma3 );

  /* Extended version of setModel -- allows rho_trans[4] (aka rho1) to be set.
   */
  void setModel( double p1, double Gamma1, double Gamma2, double Gamma3, double rho1 );

  /* pressure (in units of dyne/cm^2) from 
   * a density (in units of g/cm^3)
   */
  double pressure(double rho) const;

  /* energy density (in units of g/cm^3) from
   * a density (in units of g/cm^3)
   */
  double epsilon(double rho) const;

  /* energy density (in units of g/cm^3) from
   * a pressure (in units of dyne/cm^2)
   */
  double epsilonOfP(double P) const;

  /* energy density (in units of cm^-2) from
   * a pressure (in units of cm^-2)
   */
  double geomepsilonOfgeomP(double geomP) const;

  /* get rho_low_high
   */
  double get_rho_low_high() const { return rho_trans[3]; };

  /* The minimum pressure that makes sense */
  double p_min() const { return 0.0; };

  std::string summary() const;

protected:

  /* set up constants (low-density numbers) */
  void setUpConsts();

  double K[7];
  double Gamma[7];
  double p1;
  /* transition densities and pressure: there are 6 transitions between the 7 regions */
  double rho_trans[6];
  double P_trans[6];
  /* the a_i's are used for calculating epsilon of rho or P */
  double a[7];

};
