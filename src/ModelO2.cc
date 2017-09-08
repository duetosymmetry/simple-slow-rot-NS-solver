 /*
 * ModelO2.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Feb
 *
 * This is the implementation for a class which solves for the O(a^2) NS solution
 *
 */

#include "ModelO2.h"
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
#include "integMagicNums.h"

//////////////////////////////////////////////////////////////////////
// Forward decls
//////////////////////////////////////////////////////////////////////

// The gsl function will take two params: the model itself,
// and a boolean flag for whether or not to
// include the inhomogeneous terms in the source
struct rhsParams
{
  ModelO2 &model;
  bool inhomogeneous;
};

// Helper functions
// For the gsl integration routine
int RHS_2_gsl(double r, const double y[], double f[], void *params);

/* m2_impl implements the calculation of m2 when you already
 * have all the necessary ingredients.
 * This is purely a helper function since RHS_2_gsl
 * already gets all the inputs, so no need to do it again
 */
inline double m2_impl( double R_cm, double R4,
                       double h2,
                       double p, double rho,
                       double exp_minus_lambda,
                       double exp_minus_nu,
                       double omega1,
                       double domega1dr,
                       bool inhomogeneous
                       );

//////////////////////////////////////////////////////////////////////
// Class member functions
//////////////////////////////////////////////////////////////////////

ModelO2::ModelO2( ModelO1 &model1 )
  : bg(model1.bg), model1(model1), solved(false), splineBuilt(false),
    _K2(BG_MAX_SIZE), _h2(BG_MAX_SIZE), _Q(0),
    spline_K2(nullptr), spline_h2(nullptr),
    acc(gsl_interp_accel_alloc(), &gsl_interp_accel_free)
{
};

ModelO2::~ModelO2()
{
  safeDeallocSplines();
};

void ModelO2::assertBGSolved( const std::string & message ) const
{
#ifndef NOASSERT
  if (! (bg.isSolved() ) ) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void ModelO2::assertO1Solved( const std::string &  message ) const
{
#ifndef NOASSERT
  if (! (model1.isSolved() ) ) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void ModelO2::assertSolved( const std::string & message ) const
{
#ifndef NOASSERT
  if (!solved) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void ModelO2::assertSplineBuilt( const std::string & message ) const
{
#ifndef NOASSERT
  if (!splineBuilt) {
    std::cerr << message << std::endl;
    abort();
  };
#endif
};

void ModelO2::safeDeallocSplines()
{

  if (spline_K2)
    { gsl_spline_free(spline_K2); spline_K2=nullptr; };
  if (spline_h2)
    { gsl_spline_free(spline_h2); spline_h2=nullptr; };

};

void ModelO2::resetAccel()
{
  gsl_interp_accel_reset(acc.get());
};

void ModelO2::buildSplines()
{

  assertBGSolved("Tried to build splines when O(0) not solved.");
  assertO1Solved("Tried to build splines when O(1) not solved.");
  assertSolved("Tried to build splines when not solved.");

  // gsl spline interface is stupid ...
  // need to init with the size already known.

  safeDeallocSplines();

  spline_K2  = gsl_spline_alloc(gsl_interp_cspline, iMax()+1);
  spline_h2  = gsl_spline_alloc(gsl_interp_cspline, iMax()+1);

  double * _r = (double *) bg.r_raw();

  gsl_spline_init (spline_K2 , _r, _K2.data(),  iMax()+1);
  gsl_spline_init (spline_h2 , _r, _h2.data(),  iMax()+1);

  resetAccel();

  splineBuilt = true;

};

void ModelO2::reset()
{

  solved      = false;
  splineBuilt = false;

  safeDeallocSplines();
  resetAccel();

};

//////////////////////////////////////////////////////////////////////
// Getting solution values
//////////////////////////////////////////////////////////////////////

/* K2 (in cm^??) as a function of radius (in cm) */
double ModelO2::K2Ofr( double R_cm )
{
  assertSplineBuilt("Tried to get K2 when splines not built.");
  assertBGSolved("Tried to get K2 when O(0) not solved.");

  return gsl_spline_eval (spline_K2, R_cm, acc.get());
};

/* h2 (in cm^??) as a function of radius (in cm) */
double ModelO2::h2Ofr( double R_cm )
{
  assertSplineBuilt("Tried to get h2 when splines not built.");
  assertBGSolved("Tried to get h2 when O(0) not solved.");

  return gsl_spline_eval (spline_h2, R_cm, acc.get());
};

/* m2 (in cm^??) as a function of radius (in cm) */
double ModelO2::m2Ofr( double R_cm )
{
  assertSplineBuilt("Tried to get m2 when splines not built.");
  assertBGSolved("Tried to get m2 when O(0) not solved.");

  const double h2 = h2Ofr(R_cm);

  return m2Ofrh2( R_cm, h2, true );
};

/* xi2 (in cm^??) as a function of radius (in cm) */
double ModelO2::xi2Ofr( double R_cm )
{
  assertSplineBuilt("Tried to get xi2 when splines not built.");
  assertBGSolved("Tried to get xi2 when O(0) not solved.");

  const double h2 = h2Ofr(R_cm);

  return xi2Ofrh2( R_cm, h2 );
};

// (algebra) compute xi2 of r (in cm) and some value of h2
// This is Eq. (27) of Yagi and Yunes arXiv:1303.1528v2
double ModelO2::xi2Ofrh2( double R_cm, double h2 )
{
  assertSplineBuilt("Tried to get xi2(R,h2) when splines not built.");
  assertBGSolved(   "Tried to get xi2(R,h2) when O(0) not solved.");

  const double R2 = R_cm*R_cm;
  const double R3 = R2*R_cm;

  const double m                = bg.mOfr(R_cm);
  const double p                = bg.pOfr(R_cm);
  const double exp_minus_lambda = 1. - 2.* m/R_cm;
  const double exp_minus_nu     = exp( - bg.nuOfr(R_cm) );
  const double omega1           = model1.omega1Ofr(R_cm);

  return - R2*exp_minus_lambda*( 3.*h2 + exp_minus_nu*R2*omega1*omega1 )
           / ( 3. * ( m + 4.*M_PI*p*R3 ));
};

/* (algebra) compute m2 of r (in cm) and some value of h2,
 * OPTIONALLY turning off the "inhomogeneous" term
 * (i.e. the part of m2 not proportional to K2 or h2).
 * By default inhomogeneous terms are INCLUDED.
 */
double ModelO2::m2Ofrh2( double R_cm, double h2, bool inhomogeneous )
{
  assertBGSolved(   "Tried to get m2(R,h2) when O(0) not solved.");
  assertO1Solved(   "Tried to get m2(R,h2) when O(1) not solved.");

  const double R4 = R_cm*R_cm*R_cm*R_cm;

  const double m                = bg.mOfr(R_cm);
  const double p                = bg.pOfr(R_cm);
  const double rho              = bg.eos.geomepsilonOfgeomP(p);
  const double exp_minus_lambda = 1. - 2.* m/R_cm;
  const double exp_minus_nu     = exp( - bg.nuOfr(R_cm) );
  const double omega1           = model1.omega1Ofr(R_cm);
  const double domega1dr        = model1.domega1drOfr(R_cm);

  return m2_impl( R_cm, R4,
                  h2,
                  p, rho,
                  exp_minus_lambda,
                  exp_minus_nu,
                  omega1,
                  domega1dr,
                  inhomogeneous);

};

/* m2_impl implements the calculation of m2 when you already
 * have all the necessary ingredients.
 * This is purely a helper function since RHS_2_gsl
 * already gets all the inputs, so no need to do it again
 * This is Eq. (28) of Yagi and Yunes arXiv:1303.1528v2
 * NOTE: There is a typo in Eq. (28). The domega1dr is
 * supposed to be squared!!
 */
inline double m2_impl( double R_cm, double R4,
                       double h2,
                       double p, double rho,
                       double exp_minus_lambda,
                       double exp_minus_nu,
                       double omega1,
                       double domega1dr,
                       bool inhomogeneous
                       )
{

  const double m2_homog = - R_cm * exp_minus_lambda * h2;

  if (inhomogeneous)
  {
    return m2_homog
      + (1./6.) * R4 * exp_minus_lambda * exp_minus_nu
          * ( R_cm * exp_minus_lambda * domega1dr * domega1dr
              + 16. * M_PI * R_cm * omega1 * omega1 * (rho + p) );
  } else
    return m2_homog;

};

// Quadrupole moment
double ModelO2::Q() const
{
  assertSolved("Tried to get Q when not solved.");

  return _Q;
};

// Reduced quadrupole moment
double ModelO2::Qbar() const
{
  assertSolved(  "Tried to get Qbar when not solved.");
  assertBGSolved("Tried to get Qbar when O(0) not solved.");
  assertO1Solved("Tried to get Qbar when O(1) not solved.");

  const double I = model1.I();
  return - _Q * bg.M() / (I*I);
};

std::string ModelO2::summary() const
{
  assertBGSolved("Tried to get summary when O(0) not solved.");
  assertO1Solved("Tried to get summary when O(1) not solved.");
  assertSolved("Tried to get summary when O(2) not solved.");

  const double q = Q();

  std::ostringstream o;

  o << model1.summary()
    << "; O(2) summary: Q=" << q;

  return o.str();
};

//////////////////////////////////////////////////////////////////////
// SOLVE
//////////////////////////////////////////////////////////////////////

/* Solves the homogeneous or inhomogeneous diffeq with
 * unphysical initial conditions, to be matched later.
 * Stores the solution into the vectors K2_sol, h2_sol.
 * The unphysical initial condition is h2_in, such that
 * h2 =  h2_in * r^2 + O(higher) as r->0.
 * K2 = -h2_in * r^2 + O(higher) as r->0.
 */
void ModelO2::solveUnphysical( double * K2_sol, double * h2_sol,
                               double h2_in, bool inhomogeneous )
{
  assertBGSolved("Tried to solve O(2) when O(0) not solved.");
  assertO1Solved("Tried to solve O(2) when O(1) not solved.");

  // Set up the GSL ODE solver
  // Need a rhsParams struct
  rhsParams params = { *this,           // model
                       inhomogeneous }; // inhomogeneous parameter was passed in

  gsl_odeiv2_system sys = {RHS_2_gsl, // function which computes d/dr of system
                           0,         // Jacobian -- we aren't specifying it
                           2,         // Number of differential equations in the system
                           &params};  // params


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

  /* Set state variables with initial conditions */
  // This starts slightly away from the origin.
  r = bg.r(0);

  // Initial conditions: this is Eqs. (36) and (35)
  K2_sol[0] = y[0] = -h2_in * r*r;
  h2_sol[0] = y[1] =  h2_in * r*r;

  long i_max = bg.iMax();

  /* Integrate */
  for( i = 1; i <= i_max; i++ )
  {

    int status = gsl_odeiv2_driver_apply (ode_driver, &r, bg.r(i), y);

    if (status != GSL_SUCCESS)
    {
      printf ("error, return value=%d\n", status);
      break;
    };

    K2_sol[i] = y[0];
    h2_sol[i] = y[1];

  }

  gsl_odeiv2_driver_free (ode_driver);

};

void ModelO2::solve()
{

  assertBGSolved("Tried to solve O(2) when O(0) not solved.");
  assertO1Solved("Tried to solve O(2) when O(1) not solved.");

  ////////////////////////////////////////////////////////////
  // Approach:
  // Solve the system [Eqs. (29) and (30)] two times.
  // First the inhomogeneous solution (with bogus ICs),
  // then the homogeneous solution    (with bogus ICs).
  // Find the correct linear combination of the two solutions
  // by demanding continuity at the stellar surface (attaching
  // to the analytic exterior solution).

  std::vector<double>
    K2_homogeneous(bg.iMax() + 1),
    h2_homogeneous(bg.iMax() + 1),
    K2_particular( bg.iMax() + 1),
    h2_particular( bg.iMax() + 1);

  // WARNING! MAGIC NUMBER
  // I have checked that this choice is arbitrary and
  // doesn't affect anything.
  const double h2_in = 100.;

  solveUnphysical( K2_particular.data(), h2_particular.data(),
                   h2_in,
                   true ); // true means inhomogeneous

  solveUnphysical( K2_homogeneous.data(), h2_homogeneous.data(),
                   h2_in,
                   false ); // false means homogeneous

  ////////////////////////////////////////////////////////////
  // Match boundary conditions at surface
  // Below I took Kent's original code and replaced variable
  // names, and some pow() expressions, and
  // log(1-2M/R) with nu_surf.
  // I hope that I did this all correctly!

  const long i_max = bg.iMax();

  const double h2_final_h = h2_homogeneous[i_max];
  const double k2_final_h = K2_homogeneous[i_max];
  const double h2_final_p = h2_particular[i_max];
  const double k2_final_p = K2_particular[i_max];

  const double M = bg.M(), M2 = M*M, M3 = M2*M, M4 = M2*M2;
  const double R = bg.R(), R2 = R*R, R3 = R2*R, R4 = R2*R2;
  const double I = model1.I(), I2 = I*I;

  // In the exterior, nu = -lambda; thus at the surface, log(1-2M/R) = nu
  // See Eq. (13)
  // QUESTION: WHY is there a difference between nu(i_max) and log(1-2M/R) ?
  // const double nu_surf = bg.nu(i_max);
  const double nu_surf = log(1.-2.*M/R);

  // TODO WHERE DO THESE EXPRESSIONS COME FROM?
  const double h2_ext_p = (M+R)*I2/(R4*M);
  const double h2_ext_h = (1./2.)*(-3.*R2*pow(R-2.*M,2.)*nu_surf+18.*M2*R2-8.*M3*R-6.*M*R3-4.*M4)/(R*M2*(R-2.*M));
  const double k2_ext_p = (-2.*M-R)*I2/(R4*M);
  const double k2_ext_h = -(1./2.)*(4.*M3+6.*R*M2*nu_surf-6.*R*M2-6.*R2*M-3.*R3*nu_surf)/(R*M2);

  // TODO WHERE DO THESE EXPRESSIONS COME FROM?
  const double C2 = -(h2_final_h*k2_final_p+k2_final_h*h2_ext_p-h2_final_h*k2_ext_p-k2_final_h*h2_final_p)/(-h2_final_h*k2_ext_h+k2_final_h*h2_ext_h);
  const double C2_int = (h2_final_p*k2_ext_h-h2_ext_h*k2_final_p+h2_ext_h*k2_ext_p-h2_ext_p*k2_ext_h)/(-h2_final_h*k2_ext_h+k2_final_h*h2_ext_h);

  // Combine for the solution
  for (long i = 0; i <= i_max; i++)
  {
    _K2[i] = C2_int * K2_homogeneous[i] + K2_particular[i];
    _h2[i] = C2_int * h2_homogeneous[i] + h2_particular[i];
  };

  // Compute Q
  // TODO WHERE DOES THIS COME FROM?
  _Q = -I2/M - 8./5.*C2*M3;

  solved = true;

  buildSplines();

};

//////////////////////////////////////////////////////////////////////
// helper functions
//////////////////////////////////////////////////////////////////////

// This is what gsl uses as source term for the first order ODE
int RHS_2_gsl(double r, const double y[], double f[], void *vparams)
{

  // y[0] -- K2       (geometric units)
  // y[1] -- h2       (geometric units)
  // f[0] -- d(K2)/dr (geometric units)
  // f[1] -- d(h2)/dr (geometric units)

  // this will probably get optimized away,
  // and I don't care if it doesn't
  const double K2 = y[0], h2 = y[1];

  const rhsParams * const params = (rhsParams*) vparams;
  ModelO2 &my = params->model;

  const double r2         = r*r;
  const double r3         = r*r2;
  const double r4         = r2*r2;

  const double m          = my.bg.mOfr(r);
  const double p          = my.bg.pOfr(r);
  const double rho        = my.bg.eos.geomepsilonOfgeomP( p );
  const double exp_minus_nu = exp( - my.bg.nuOfr(r) );
  const double omega1     = my.model1.omega1Ofr(r);
  const double domega1dr  = my.model1.domega1drOfr(r);

  const double exp_minus_lambda = 1. - 2.* m/r;
  const double exp_lambda       = 1./exp_minus_lambda;

  // We need the "m2" function for this r and h2 (stored in y[1])
  // but remember that we need to evaluate it either with or
  // without the "homogeneous" terms included.
  // Here we use the m2_impl function because we already fetched
  // most of the things we need to calculate it.
  const double m2 = m2_impl( r, r4,
                             h2,
                             p, rho,
                             exp_minus_lambda,
                             exp_minus_nu,
                             omega1,
                             domega1dr,
                             params->inhomogeneous);

  ////////////////////////////////////////////////////////////
  // Eq. (29) and (30) of arXiv:1303.1528
  // are a system for K2' and h2', of the form
  //   K2' + h2' = rest29
  // C K2' + h2' = rest30
  // from which we get
  //  K2' = ( - rest29 + rest30) / (C - 1.)
  //  h2' = ( C*rest29 - rest30) / (C - 1.)
  // NOTE: For r/R << 1, C is very close to 1.
  //       This causes the system to be poorly conditioned.
  //       Therefore here we compute C-1, which can be found
  //       accurately.
  // const double C = (r - m + 4.*M_PI*p*r3)/r * exp_lambda;
  const double Cminus1 = (m + 4.*M_PI*p*r3)/(r - m - m);

  // The only inhomogeneous parts here are in m2,
  // so that is automatically taken care of
  const double rest29 = (r - 3.*m - 4*M_PI*p*r3)/r2 * exp_lambda * h2
                  + (r - m + 4.*M_PI*p*r3)/r3 * exp_lambda * exp_lambda * m2;

  // The inhomogeneous part of m2 is automatically taken care of
  double rest30 = (3. - 4.*M_PI*(p+rho)*r2)/r * exp_lambda * h2
    + 2./r * exp_lambda * K2
    + (1. + 8.*M_PI*p*r2)/r2 * exp_lambda * exp_lambda * m2;

  if (params->inhomogeneous) // include inhomogeneous terms
  {
    rest30 += r3/12. * exp_minus_nu * domega1dr * domega1dr
      - 4./3.*M_PI*(p+rho)*r3* omega1 * omega1 * exp_minus_nu * exp_lambda;
  };

  ////////////////////////////////////////////////////////////
  /* d(K2)/dr */ f[0] = (              - rest29 + rest30) / Cminus1;
  /* d(h2)/dr */ f[1] = ( (1. + Cminus1)*rest29 - rest30) / Cminus1;
  return GSL_SUCCESS;

}
