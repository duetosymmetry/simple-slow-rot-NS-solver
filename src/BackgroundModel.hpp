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

#include "ppEOS.hpp"
#include "Conformal.hpp"
#include <string>
#include <vector>
#include <memory>
#include <gsl/gsl_spline.h>

class BackgroundModel {

 public:

  /* takes pc in units of cm^-2 */
  BackgroundModel( const ppEOS &eos,
                   const Conformal &conf,
                   double pc, double phic);

  ~BackgroundModel();

  void reset();
  void solve();
  /* takes pc in units of cm^-2, also resets. */
  void setpc( double new_pc ) { pc = new_pc; reset(); };
  /* takes phic in units of cm^0, also resets. */
  void setphic( double new_phic ) { phic = new_phic; reset(); };

  /* returns if we currently have a solution for a model */
  bool isSolved() const { return solved; };

  /* physical solution */

  // Access to interpolants

  /* m (in cm^1) as a function of radius (in cm) */
  double muOfr   ( double R_cm );
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

  double r  ( unsigned long i ) const { return _r  [i]; };
  double mu ( unsigned long i ) const { return _mu [i]; };
  double nu ( unsigned long i ) const { return _nu [i]; };
  double phi( unsigned long i ) const { return _phi[i]; };
  double psi( unsigned long i ) const { return _psi[i]; };
  double p  ( unsigned long i ) const { return _p  [i]; };
  double mbar( unsigned long i ) const { return _mbar[i]; };

  // Derived quantities
  // The areal radius in Einstein frame (in cm)
  double R_areal_Einstein() const;
  // The areal radius in Jordan frame (in cm)
  double R_areal_Jordan() const;
  // The ADM mass (in cm)
  double M_ADM() const;
  // The baryonic mass (in cm)
  double M_b() const;
  // The asymptotic scalar field phi_0 [cm^0]
  double phi_0() const;
  // Scaled charge alpha [cm^0]
  double alpha() const;
  // Scalar charge omega [cm^1]
  double omega() const;

  // This is public ... it's const so it'll never get modified anyway

  const ppEOS &eos;
  const Conformal &conf;

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
  // the [0] elements of _r, _mu, _nu, _phi, _psi, _p
  void initialConditions();

  // Compute the aux surface/body quantities
  void computeSurfaceBodyQuantities(double r_s, double y_s[]);

  //// Storage

  // Parameters
  /* pc is in units of cm^-2 */
  double pc;
  double phic;

  // State
  bool solved, splineBuilt;

  long i_max;

  std::vector<double>
    _r,
    _mu,
    _nu,
    _phi,
    _psi,
    _p,
    _mbar;

  gsl_spline 
    *spline_mu,
    *spline_nu,
    *spline_phi,
    *spline_psi,
    *spline_p;

  std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)>
    acc;

  // Aux quantities for surface values, charge, phi_0, etc.
  double nu_prime_s;
  double _phi_0;
  double m_ADM;
  double alpha_body;
  double omega_body;

};
