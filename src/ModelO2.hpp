/*
 * ModelO2.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Feb
 *
 * This is the interface for a class which solves for the O(a^2) NS solution
 *
 */

#pragma once

#include "BackgroundModel.hpp"
#include "ModelO1.hpp"
#include <string>
#include <vector>
#include <memory>
#include <gsl/gsl_spline.h>

class ModelO2 {

 public:

  ModelO2( ModelO1 &model1 );

  ~ModelO2();

  void reset();
  void solve();

  /* returns if we currently have a solution for a model */
  bool isSolved() const { return solved; };

  /* physical solution */

  // Access to interpolants
  /* K2 (in cm^?) as a function of radius (in cm) */
  double K2Ofr ( double R_cm );
  /* h2 (in cm^?) as a function of radius (in cm) */
  double h2Ofr ( double R_cm );
  /* m2 (in cm^?) as a function of radius (in cm) */
  double m2Ofr ( double R_cm );
  /* xi2 (in cm^?) as a function of radius (in cm) */
  double xi2Ofr ( double R_cm );

  // Access to points in solution

  // The maximum index value which may be passed to f, q
  // i.e. the solutions are indexed with 0 <= i <= iMax().
  // This just passes through to the background solution
  long   iMax() const { return bg.iMax(); };

  double K2  ( unsigned long i ) const { return _K2 [i]; };
  double h2  ( unsigned long i ) const { return _h2 [i]; };

  // Derived quantities
  // Quadrupole
  double Q() const;
  // Reduced quadrupole
  double Qbar() const;

  // Calculation functions
  // (algebra) compute xi2 of r (in cm) and some value of h2
  double xi2Ofrh2( double R_cm, double h2 );
  /* (algebra) compute m2 of r (in cm) and some value of h2,
   * OPTIONALLY turning off the "inhomogeneous" term
   * (i.e. the part of m2 not proportional to K2 or h2).
   * By default inhomogeneous terms are INCLUDED.
   */
  double m2Ofrh2( double R_cm, double h2, bool inhomogeneous = true );

  // The bg and O(1) models we're based upon
  // Thanks to Curran Muhlberger for helping me realize that this can't be const.
  BackgroundModel &bg;
  ModelO1         &model1;

  std::string summary() const;

protected:

  void assertBGSolved(    const std::string & message ) const;
  void assertO1Solved(    const std::string & message ) const;
  void assertSolved(      const std::string & message ) const;
  void assertSplineBuilt( const std::string & message ) const;

  void safeDeallocSplines();
  void resetAccel();

  void buildSplines();

  /* Solves the homogeneous or inhomogeneous diffeq with
   * unphysical initial conditions, to be matched later.
   * Stores the solution into the vectors K2_sol, h2_sol.
   * The unphysical initial condition is h2_in, such that
   * h2 =  h2_in * r^2 + O(higher) as r->0.
   * K2 = -h2_in * r^2 + O(higher) as r->0.
   */
  void solveUnphysical( double * K2_sol, double * h2_sol,
			double h2_in, bool inhomogeneous );

  // Storage

  bool solved, splineBuilt;

  std::vector<double>_K2, _h2;

  double _Q;

  gsl_spline 
    *spline_K2,
    *spline_h2;

  std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>
    acc;

};
