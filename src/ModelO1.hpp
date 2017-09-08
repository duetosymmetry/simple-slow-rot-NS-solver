/*
 * ModelO1.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * This is the interface for a class which solves for the O(a^1) NS solution
 *
 */

#pragma once

#include "BackgroundModel.hpp"
#include <string>
#include <vector>
#include <memory>
#include <gsl/gsl_spline.h>

class ModelO1 {

 public:

  ModelO1( BackgroundModel &bg );

  ~ModelO1();

  void reset();
  void solve();

  /* returns if we currently have a solution for a model */
  bool isSolved() const { return solved; };

  /* physical solution */

  // Access to interpolants
  /* omega1 (in cm^-1) as a function of radius (in cm) */
  double omega1Ofr ( double R_cm );
  /* domega1dr (in cm^-2) as a function of radius (in cm) */
  double domega1drOfr ( double R_cm );

  // Access to points in solution

  // The maximum index value which may be passed to f, q
  // i.e. the solutions are indexed with 0 <= i <= iMax().
  // This just passes through to the background solution
  long   iMax() const { return bg.iMax(); };

  double omega1    ( unsigned long i ) const { return _omega1    [i]; };
  double domega1dr ( unsigned long i ) const { return _domega1dr [i]; };

  // Derived quantities
  // Moment of inertia
  double I() const;
  // Reduced moment of inertia
  double Ibar() const;

  // The bg model we're based upon
  // Thanks to Curran Muhlberger for helping me realize that this can't be const.
  BackgroundModel &bg;

  std::string summary() const;

protected:

  void assertBGSolved(    const std::string & message ) const;
  void assertSolved(      const std::string & message ) const;
  void assertSplineBuilt( const std::string & message ) const;

  void safeDeallocSplines();
  void resetAccel();

  void buildSplines();

  // Calculates initial conditions and stores them in
  // the [0] elements of _omega1, _domega1dr
  void initialConditions();

  // Storage

  bool solved, splineBuilt;

  std::vector<double> _omega1, _domega1dr;

  gsl_spline 
    *spline_omega1,
    *spline_domega1dr;

  std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>
    acc;

};
