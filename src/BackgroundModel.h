/*
 * BackgroundModel.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * This is the interface for a class which solves for a background
 * neutron star model.
 *
 */

#pragma once

#include "ppEOS.h"
#include <string>
#include <vector>
#include <memory>
#include <gsl/gsl_spline.h>

class BackgroundModel {

 public:

  /* takes pc in units of cm^-2 */
  BackgroundModel( const ppEOS &eos, double pc );

  ~BackgroundModel();

  void reset();
  void solve();
  /* takes pc in units of cm^-2, also resets. */
  void setpc( double new_pc ) { pc = new_pc; reset(); };

  /* returns if we currently have a solution for a model */
  bool isSolved() const { return solved; };

  /* physical solution */

  // Access to interpolants

  /* m (in cm^1) as a function of radius (in cm) */
  double mOfr    ( double R_cm );
  /* nu (in cm^0) as a function of radius (in cm) */
  double nuOfr   ( double R_cm );
  /* pressure (in cm^-2) as a function of radius (in cm) */
  double pOfr    ( double R_cm );
  /* rho (in cm^-2) as a function of radius (in cm) */
  double rhoOfr  ( double R_cm );
  /* dnudR (in cm^-1) as a function of radius (in cm) */
  double dnudROfr( double R_cm );

  // Access to points in solution

  // The maximum index value which may be passed to r, m, nu, p
  // i.e. the solutions are indexed with 0 <= i <= iMax().
  long   iMax() const { return i_max; };

  double r ( unsigned long i ) const { return _r [i]; };
  double m ( unsigned long i ) const { return _m [i]; };
  double nu( unsigned long i ) const { return _nu[i]; };
  double p ( unsigned long i ) const { return _p [i]; };

  // Derived quantities
  // The stellar radius (in cm)
  double R() const;
  // The total mass (in cm)
  double M() const;

  // This is public ... it's const so it'll never get modified anyway

  const ppEOS &eos;

  // Access to raw data. I will return this as const but
  // the caller may not respect the constness ... :(
  const double * r_raw() const { return _r.data(); };

  std::string summary() const;

protected:

  void assertSolved(      const std::string & message ) const;
  void assertSplineBuilt( const std::string & message ) const;

  void safeDeallocSplines();
  void resetAccel();

  void buildSplines();

  // Calculates initial conditions and stores them in
  // the [0] elements of _r, _m, _nu, _p
  void initialConditions();

  // Storage

  /* pc is in units of cm^-2 */
  double pc;

  bool solved, splineBuilt;

  long i_max;

  std::vector<double>
    _r,
    _m,
    _nu,
    _p;

  gsl_spline 
    *spline_m,
    *spline_nu,
    *spline_p;

  std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>
    acc;

};
