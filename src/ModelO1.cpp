/*
 * ModelO1.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * This is the implementation for a class which solves for the O(a^1) NS solution
 *
 */

#include "ModelO1.hpp"
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
int RHS_1_gsl(double r, const double y[], double f[], void *params);

//////////////////////////////////////////////////////////////////////
// Class member functions
//////////////////////////////////////////////////////////////////////

ModelO1::ModelO1( BackgroundModel &bg )
  : bg(bg), solved(false), splineBuilt(false),
    _omega1(BG_MAX_SIZE), _domega1dr(BG_MAX_SIZE),
    spline_omega1(nullptr), spline_domega1dr(nullptr),
    acc(gsl_interp_accel_alloc(), &gsl_interp_accel_free)
{
};

ModelO1::~ModelO1()
{
  safeDeallocSplines();
};

void ModelO1::assertBGSolved( const std::string & message ) const
{
#ifndef NOASSERT
  if (! (bg.isSolved() ) ) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void ModelO1::assertSolved( const std::string & message ) const
{
#ifndef NOASSERT
  if (!solved) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void ModelO1::assertSplineBuilt( const std::string & message ) const
{
#ifndef NOASSERT
  if (!splineBuilt) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void ModelO1::safeDeallocSplines()
{  

  if (spline_omega1)
    { gsl_spline_free(spline_omega1); spline_omega1=nullptr; };
  if (spline_domega1dr)
    { gsl_spline_free(spline_domega1dr); spline_domega1dr=nullptr; };

};

void ModelO1::resetAccel()
{

  gsl_interp_accel_reset(acc.get());

};

void ModelO1::buildSplines()
{

  assertBGSolved("Tried to build splines when O(0) not solved.");
  assertSolved("Tried to build splines when not solved.");

  // gsl spline interface is stupid ...
  // need to init with the size already known.

  safeDeallocSplines();

  spline_omega1    = gsl_spline_alloc(gsl_interp_cspline, iMax()+1);
  spline_domega1dr = gsl_spline_alloc(gsl_interp_cspline, iMax()+1);

  double * _r = (double *) bg.r_raw();

  gsl_spline_init (spline_omega1    , _r, _omega1.data(),    iMax()+1);
  gsl_spline_init (spline_domega1dr , _r, _domega1dr.data(), iMax()+1);

  resetAccel();

  splineBuilt = true;

};

void ModelO1::reset()
{

  solved      = false;
  splineBuilt = false;

  safeDeallocSplines();
  resetAccel();

};

//////////////////////////////////////////////////////////////////////
// Getting solution values
//////////////////////////////////////////////////////////////////////

/* omega1 (in cm^-1) as a function of radius (in cm) */
double ModelO1::omega1Ofr( double R_cm )
{
  assertSplineBuilt("Tried to get omega1 when splines not built.");
  assertBGSolved("Tried to get omega1 when O(0) not solved.");

  return gsl_spline_eval (spline_omega1, R_cm, acc.get());
};

/* domega1dr (in cm^-2) as a function of radius (in cm) */
double ModelO1::domega1drOfr( double R_cm )
{
  assertSplineBuilt("Tried to get domega1dr when splines not built.");
  assertBGSolved("Tried to get domega1dr when O(0) not solved.");

  return gsl_spline_eval (spline_domega1dr, R_cm, acc.get());
};

// Moment of inertia
double ModelO1::I() const
{
  assertBGSolved("Tried to get I when O(0) not solved.");
  assertSolved("Tried to get I when O(1) not solved.");

  // TODO: Document where this comes from.
  const double r = bg.R();
  return 1/6. * (r*r*r*r) * _domega1dr[iMax()];
};

// Reduced moment of inertia
double ModelO1::Ibar() const
{
  assertBGSolved("Tried to get Ibar when O(0) not solved.");
  assertSolved("Tried to get Ibar when O(1) not solved.");

  const double M = bg.M(), M3=M*M*M;
  // TODO: Document where this comes from.
  return I() / M3;
};

std::string ModelO1::summary() const
{
  assertBGSolved("Tried to get summary when O(0) not solved.");
  assertSolved("Tried to get summary when O(1) not solved.");

  const double i = I();

  std::ostringstream o;

  o << bg.summary()
    << "; O(1) summary: I=" << i;

  return o.str();
};

//////////////////////////////////////////////////////////////////////
// Initial conditions
//////////////////////////////////////////////////////////////////////
// Calculates initial conditions and stores them in
// the [0] elements of _omega1, _domega1dr
void ModelO1::initialConditions()
{

  assertBGSolved("Tried to initialize O(1) when O(0) not solved.");

  double r    = bg.r(0);
  double pc   = bg.p(0);
  double rhoc = bg.eos.geomepsilonOfgeomP( pc );

  /* This I.C. does not satisfy the B.C. at the surface. We need to rescale later. */
  double omega1_c0 = 1.0;

  _omega1[0]    = omega1_c0+(8.0/5.0)*M_PI*(rhoc+pc)*omega1_c0*pow(r,2.0);
  _domega1dr[0] = (16.0/5.0)*M_PI*(rhoc+pc)*omega1_c0*r;

};

//////////////////////////////////////////////////////////////////////
// SOLVE
//////////////////////////////////////////////////////////////////////
void ModelO1::solve()
{

  assertBGSolved("Tried to solve O(1) when O(0) not solved.");

  // Set up the GSL ODE solver
  gsl_odeiv2_system sys = {RHS_1_gsl, // function which computes d/dr of system
                           0,         // Jacobian -- we aren't specifying it
                           2,         // Number of equations in the system
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
  double y[2];

  /* Set up initial conditions in [0] elements of storage */
  initialConditions();

  /* Set state variables with initial conditions */
  r    = bg.r(0);

  y[0] = _omega1[0];
  y[1] = _domega1dr[0];

  long i_max = bg.iMax();

  /* Integration (first order in spin) */
  for( i = 1; i <= i_max; i++ )
  {

    int status = gsl_odeiv2_driver_apply (ode_driver, &r, bg.r(i), y);

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value=%d\n", status);
      break;
    };

    _omega1[i]    = y[0];
    _domega1dr[i] = y[1];

  }

  gsl_odeiv2_driver_free (ode_driver);

  /* Surface Values */
  double r_final         = bg.r(i_max);
  double omega1_final    = _omega1[i_max];
  double domega1dr_final = _domega1dr[i_max];

  double zeta = 1.0/(omega1_final + r_final/3.0*domega1dr_final);

  /* Normalization Values */
  for(i=0; i<=i_max; i++)
  {
    _omega1[i]    *= zeta;
    _domega1dr[i] *= zeta;
  };

  solved = true;

  buildSplines();

};

//////////////////////////////////////////////////////////////////////
// helper functions
//////////////////////////////////////////////////////////////////////

// This is what gsl uses as source term for the first order ODE
int RHS_1_gsl(double r, const double y[], double f[], void *params)
{

  // y[0] -- omega1          (geometric units)
  // y[1] -- domega1dr       (geometric units)
  // f[0] -- d(omega1)/dr    (geometric units)
  // f[1] -- d(domega1dr)/dr (geometric units)

  ModelO1 *myModel = (ModelO1*) params;

  const double m   = myModel->bg.mOfr(r);
  const double p   = myModel->bg.pOfr(r);
  const double rho = myModel->bg.eos.geomepsilonOfgeomP( p );

  // This is Eq. (19) of arXiv:1303.1528
  /* d(omega1)/dr = */	  f[0] = y[1];
  /* d(domega1dr)/dr = */ f[1] = ( p + rho )*4.0*M_PI/(1.0 - 2.0*m/r)*(y[1]*r+4.0*y[0]) - 4.0*y[1]/r;
  return GSL_SUCCESS;

}
