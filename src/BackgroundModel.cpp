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

BackgroundModel::BackgroundModel( const ppEOS &eos,
                                  const Conformal &conf,
                                  double pc, double phic)
  : eos(eos), conf(conf),
    pc(pc), phic(phic),
    solved(false), splineBuilt(false), i_max(0),
    _r(BG_MAX_SIZE), _mu(BG_MAX_SIZE), _nu(BG_MAX_SIZE),
    _phi(BG_MAX_SIZE), _psi(BG_MAX_SIZE), _p(BG_MAX_SIZE),
    spline_mu(nullptr), spline_nu(nullptr),
    spline_phi(nullptr), spline_psi(nullptr), spline_p(nullptr),
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
    { gsl_spline_free(spline_p);   spline_p=nullptr; };
  if (spline_psi)
    { gsl_spline_free(spline_psi); spline_psi=nullptr; };
  if (spline_phi)
    { gsl_spline_free(spline_phi); spline_phi=nullptr; };
  if (spline_nu)
    { gsl_spline_free(spline_nu);  spline_nu=nullptr; };
  if (spline_mu)
    { gsl_spline_free(spline_mu);  spline_mu=nullptr; };

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

  spline_mu  = gsl_spline_alloc(gsl_interp_cspline, i_max+1);
  spline_nu  = gsl_spline_alloc(gsl_interp_cspline, i_max+1);
  spline_phi = gsl_spline_alloc(gsl_interp_cspline, i_max+1);
  spline_psi = gsl_spline_alloc(gsl_interp_cspline, i_max+1);
  spline_p   = gsl_spline_alloc(gsl_interp_cspline, i_max+1);

  gsl_spline_init (spline_mu , _r.data(), _mu.data(),  i_max+1);
  gsl_spline_init (spline_nu , _r.data(), _nu.data(),  i_max+1);
  gsl_spline_init (spline_phi, _r.data(), _phi.data(), i_max+1);
  gsl_spline_init (spline_psi, _r.data(), _psi.data(), i_max+1);
  gsl_spline_init (spline_p  , _r.data(), _p.data(),   i_max+1);

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

/* mu (in cm^1) as a function of radius (in cm) */
double BackgroundModel::muOfr( double R_cm )
{
  assertSplineBuilt("Tried to get mu when splines not built.");

  return gsl_spline_eval (spline_mu, R_cm, acc.get());
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

  return eos.geomepsilonOfgeomP( p(R_cm) );

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
  
  return _mu[i_max];
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
// the [0] elements of _r, _mu, _nu, _phi, _psi, _p
void BackgroundModel::initialConditions()
{

  const double A_c = conf.A(phic);
  const double A_c4 = A_c*A_c*A_c*A_c;
  const double alpha_c = conf.alpha(phic);
  const double eps_c = eos.geomepsilonOfgeomP( pc );
  const double mu_3 = 8.*M_PI*A_c4*eps_c;
  const double nu_2 = 8.*M_PI*A_c4*pc + mu_3/3.;
  const double psi_1 = 4./3.*M_PI*A_c4*alpha_c*(eps_c-3.*pc);  
  const double p_2   = -(eps_c+pc)*(0.5*nu_2+alpha_c*psi_1);

  // We have to start slightly away from 0
  const double r  = EPSR;
  const double r2 = r*r;
  const double r3 = r2*r;

  _r [0]  = r;
  _mu[0]  = mu_3*r3/6.;
  _nu[0]  = 0.5*nu_2*r2;
  _phi[0] = phic + 0.5*psi_1*r2;
  _psi[0] = psi_1*r;
  _p [0]  = pc + 0.5*p_2*r2;

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
                           5,         // Number of equations in the system
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
  double y[5];

  /* Set up initial conditions in [0] elements of storage */
  initialConditions();

  /* Set state variables with initial conditions */
  r    = _r [0];

  y[0] = _mu[0];
  y[1] = _nu[0];
  y[2] = _phi[0];
  y[3] = _psi[0];
  y[4] = _p [0];

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

    _r [i+1] = r;

    _mu[i+1]  = y[0];
    _nu[i+1]  = y[1];
    _phi[i+1] = y[2];
    _psi[i+1] = y[3];
    _p [i+1]  = y[4];

    i_max=i;

  }

  gsl_odeiv2_driver_free (ode_driver);

  /* Surface Values */
  double r_final = r;
  double mu_final = _mu[i_max];
  
  // TODO
  // Compute surface values, save summary quantities like
  // \phi_0, \alpha_A, \omega_A, m_A
  
  solved = true;

  buildSplines();

};

//////////////////////////////////////////////////////////////////////
// helper functions
//////////////////////////////////////////////////////////////////////

// This is what gsl uses as source term for the first order ODE
int RHS_0_gsl(double r, const double y[], double f[], void *params)
{
  // Compiler should get rid of these variables
  const double mu  = y[0]; // cm^-1
  // nu does not appear in any RHS expression
  // const double nu  = y[1]; // cm^0
  const double phi = y[2]; // cm^0
  const double psi = y[3]; // cm^-1
  const double p   = y[4]; // ? geometric units

  const double psi2 = psi*psi;

  const double rMinusMu = r - mu;
  const double rMinus2Mu = r - 2.*mu;

  const BackgroundModel *myModel = (BackgroundModel*) params;

  const double A = myModel->conf.A(phi);
  const double A4 = A*A*A*A;
  const double alpha = myModel->conf.alpha(phi);

  const double eps = myModel->eos.geomepsilonOfgeomP(p);
  const double r2  = r*r;

  // f[0] -- dmu/dr  (cm^-2)
  // f[1] -- dnu/dr  (cm^-1)
  // f[2] -- dphi/dr (cm^-1)
  // f[3] -- dpsi/dr (cm^-2)
  // f[4] -- dp/dr   (?)

  /* dm/dr = */   f[0] = 4.*M_PI*r2*A4*eps + 0.5*r*rMinus2Mu*psi2;
  /* dnu/dr = */  f[1] = 8.*M_PI*r2*A4*p/rMinus2Mu + r*psi2 + 2.*mu/(r*rMinus2Mu);
  /* dphi/dr = */ f[2] = psi;
  /* dpsi/dr = */ f[3] = 4.*M_PI*r*A4/rMinus2Mu
                           * ( alpha*(eps-3.*p) + r*psi*(eps-p) )
                         - 2.*rMinusMu/(r*rMinus2Mu) * psi;
  /* dp/dr = */   f[4] = - (eps+p)*( 0.5*f[1] + alpha*psi );
  
  return GSL_SUCCESS;
}
