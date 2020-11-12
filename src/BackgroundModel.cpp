/*
 * BackgroundModel.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * This is the implementation for a class which solves for a background
 * neutron star model.
 *
 */

#include "BackgroundModel.hpp"
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

// for integration magic numbers:
// BG_MAX_SIZE, EPSR, DR, EPSABS, EPSREL, P_RATIO
#include "integMagicNums.hpp"

//////////////////////////////////////////////////////////////////////
// Forward decls
//////////////////////////////////////////////////////////////////////

// For the gsl integration routine
int RHS_0_gsl(double r, const double y[], double f[], void *params);

//////////////////////////////////////////////////////////////////////
// Class member functions
//////////////////////////////////////////////////////////////////////

BackgroundModel::BackgroundModel( const ppEOS &eos, double pc )
  : eos(eos), pc(pc), solved(false), splineBuilt(false), i_max(0),
    _r(BG_MAX_SIZE), _m(BG_MAX_SIZE), _nu(BG_MAX_SIZE), _p(BG_MAX_SIZE),
    spline_m(nullptr), spline_nu(nullptr), spline_p(nullptr),
    acc(gsl_interp_accel_alloc(), &gsl_interp_accel_free)
{
};

BackgroundModel::~BackgroundModel()
{
  safeDeallocSplines();
};

void BackgroundModel::assertSolved( const std::string & message ) const
{
#ifndef NOASSERT
  if (!solved) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void BackgroundModel::assertSplineBuilt( const std::string & message ) const
{
#ifndef NOASSERT
  if (!splineBuilt) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void BackgroundModel::safeDeallocSplines()
{  

  if (spline_p)
    { gsl_spline_free(spline_p);  spline_p=nullptr; };
  if (spline_nu)
    { gsl_spline_free(spline_nu); spline_nu=nullptr; };
  if (spline_m)
    { gsl_spline_free(spline_m);  spline_m=nullptr; };

};

void BackgroundModel::resetAccel()
{  

  gsl_interp_accel_reset(acc.get());

};

void BackgroundModel::buildSplines()
{
  assertSolved("Tried to build splines when not solved.");

  // gsl spline interface is stupid ...
  // need to init with the size already known.

  safeDeallocSplines();

  spline_m   = gsl_spline_alloc(gsl_interp_cspline, i_max+1);
  spline_nu  = gsl_spline_alloc(gsl_interp_cspline, i_max+1);
  spline_p   = gsl_spline_alloc(gsl_interp_cspline, i_max+1);

  gsl_spline_init (spline_m  , _r.data(), _m.data(),  i_max+1);
  gsl_spline_init (spline_nu , _r.data(), _nu.data(), i_max+1);
  gsl_spline_init (spline_p  , _r.data(), _p.data(),  i_max+1);

  resetAccel();

  splineBuilt = true;
};

void BackgroundModel::reset()
{
  solved      = false;
  splineBuilt = false;

  i_max = 0;

  safeDeallocSplines();
  resetAccel();

  // TODO: do we want to clear anything?
};

//////////////////////////////////////////////////////////////////////
// Getting solution values
//////////////////////////////////////////////////////////////////////

/* m (in cm^1) as a function of radius (in cm) */
double BackgroundModel::mOfr( double R_cm )
{
  assertSplineBuilt("Tried to get m when splines not built.");

  return gsl_spline_eval (spline_m, R_cm, acc.get());
};

/* nu (in cm^0) as a function of radius (in cm) */
double BackgroundModel::nuOfr( double R_cm )
{
  assertSplineBuilt("Tried to get nu when splines not built.");

  return gsl_spline_eval (spline_nu, R_cm, acc.get());
};

/* pressure (in cm^-2) as a function of radius (in cm) */
double BackgroundModel::pOfr( double R_cm )
{
  assertSplineBuilt("Tried to get p when splines not built.");

  return gsl_spline_eval (spline_p, R_cm, acc.get());
};

/* rho (in cm^-2) as a function of radius (in cm) */
double BackgroundModel::rhoOfr( double R_cm )
{

  return eos.geomepsilonOfgeomP( pOfr(R_cm) );

};

/* dnudR (in cm^-1) as a function of radius (in cm) */
double BackgroundModel::dnudROfr( double R_cm )
{
  assertSplineBuilt("Tried to get dnudR when splines not built.");

  return gsl_spline_eval_deriv( spline_nu, R_cm, acc.get() );
};

// The stellar radius (in cm)
double BackgroundModel::R() const
{
  assertSolved("Tried to get R when not solved.");
  
  return _r[i_max];
};

// The total mass (in cm)
double BackgroundModel::M() const
{
  assertSolved("Tried to get M when not solved.");
  
  return _m[i_max];
};

std::string BackgroundModel::summary() const
{
  assertSolved("Tried to get summary when not solved.");

  const double m = M(), r = R();

  std::ostringstream o;

  o << eos.summary()
    << "; BG summary: M=" << m
    << ", R=" << r
    << ", C=" << m/r;

  return o.str();
};

//////////////////////////////////////////////////////////////////////
// Initial conditions
//////////////////////////////////////////////////////////////////////
// Calculates initial conditions and stores them in
// the [0] elements of _r, _m, _nu, _p
void BackgroundModel::initialConditions()
{

  const double rhoc = eos.geomepsilonOfgeomP( pc );
  const double rho2 = 0.;
  const double p2   = -(2.0/3.0*(rhoc+pc))*M_PI*(rhoc+3.0*pc);

  // We have to start slightly away from 0
  const double r  = EPSR;
  const double r2 = r*r;
  const double r3 = r2*r;
  const double r5 = r3*r2;

  _r [0]  = r;
  _p [0]  = pc+p2*r2;
  _m [0]  = (4.0/3.0)*M_PI*rhoc*r3 + (4.0/5.0)*M_PI*rho2*r5;
  _nu[0]  = (4.0/3.0)*M_PI*r2*(rhoc + 3.0*pc);

};

//////////////////////////////////////////////////////////////////////
// SOLVE
//////////////////////////////////////////////////////////////////////
void BackgroundModel::solve()
{

  // TODO: check if we've already solved?

  // Set up the GSL ODE solver
  gsl_odeiv2_system sys = {RHS_0_gsl, // function which computes d/dr of system
                           0,         // Jacobian -- we aren't specifying it
                           3,         // Number of equations in the system
                           this};     // params-- pointer to me.


  gsl_odeiv2_driver *ode_driver =
    gsl_odeiv2_driver_alloc_y_new (&sys,                  // system
                                   gsl_odeiv2_step_rkf45, // stepper: embedded RK45
                                   EPSR,                  // initial step size
                                   EPSABS,                // absolute error
                                   EPSREL);               // relative error

  /* State variables for integration */
  int i;
  double r;
  const double dr = DR;
  double y[3];

  /* Set up initial conditions in [0] elements of storage */
  initialConditions();

  /* Set state variables with initial conditions */
  r    = _r [0];

  y[0] = _p [0];
  y[1] = _m [0];
  y[2] = _nu[0];

  const double p_min = eos.p_min();

  /* Integration (zeroth-order in spin) */
  for(i=0; _p[i] > P_RATIO*pc && _p[i] > p_min && i<BG_MAX_SIZE-1; i++)
  {

      int status = gsl_odeiv2_driver_apply (ode_driver, &r, r+dr, y);

      if (status != GSL_SUCCESS)
      {
        printf ("error, return value=%d\n", status);
        break;
      }

    _p [i+1] = y[0];
    _m [i+1] = y[1];
    _nu[i+1] = y[2];
    _r [i+1] = r;

    i_max=i;

  }

  gsl_odeiv2_driver_free (ode_driver);

  /* Surface Values */
  double m_final = _m[i_max];
  double r_final = r;
  double nu_final = _nu[i_max];

  /* Normalization Values */
  double normalization_nu = log(1.0 - 2.0*m_final/r_final) - nu_final; 
  /* This is such that newnu = oldnu + normaliazation_nu = 1 - 2 M_f/R_f */
  for(i=0; i<=i_max; i++)
    _nu[i] += normalization_nu;

  solved = true;

  buildSplines();

};

//////////////////////////////////////////////////////////////////////
// helper functions
//////////////////////////////////////////////////////////////////////

// This is what gsl uses as source term for the first order ODE
int RHS_0_gsl(double r, const double y[], double f[], void *params)
{
  // y[0] -- p      (geometric units)
  // y[1] -- m      (geometric units)
  // y[2] -- nu     (geometric units)
  // f[0] -- dp/dr  (geometric units)
  // f[1] -- dm/dr  (geometric units)
  // f[2] -- dnu/dr (geometric units)

  const BackgroundModel *myModel = (BackgroundModel*) params;
  const double eps = myModel->eos.geomepsilonOfgeomP(y[0]);
  const double r2  = r*r;
  const double r3  = r2*r;

  /* dp/dr = */  f[0] = -((y[0]+eps)*(y[1]+4.0*M_PI*r3*y[0]))/(r*(r-2.0*y[1]));

  /* dm/dr = */  f[1] = 4.0*M_PI*r2*eps;
  /* dnu/dr = */ f[2] = 2.0*(y[1]+4.0*M_PI*r3*y[0])/(r*(r-2.0*y[1]));
  return GSL_SUCCESS;
}
