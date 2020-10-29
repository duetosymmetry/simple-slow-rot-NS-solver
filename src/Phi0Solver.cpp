/*
 * Phi0Solver.cpp
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2020 Oct
 *
 * Implementation for class Phi0Solver
 *
 */

#include "Phi0Solver.hpp"

#include <memory>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

double phi0_error_gsl(double phic, void *params);

const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;

Phi0Solver::Phi0Solver( double epsabs, double epsrel, int max_iter )
  : epsabs(epsabs), epsrel(epsrel), max_iter(max_iter),
    gsl_solver(gsl_root_fsolver_alloc(T), gsl_root_fsolver_free)
{
  F.function = &phi0_error_gsl;
};

struct Phi0Params {
  BackgroundModel *model0;
  double phi0_goal;
};

int Phi0Solver::solve( BackgroundModel &model0, double phi0_goal, double phic_low, double phic_high )
{
  Phi0Params params {&model0, phi0_goal};
  F.params = &params;

  gsl_root_fsolver_set(gsl_solver.get(), &F, phic_low, phic_high);

  int status;
  int iter = 0;

  do {
    iter++;
    status = gsl_root_fsolver_iterate (gsl_solver.get());
    // r = gsl_root_fsolver_root (gsl_solver);
    phic_low = gsl_root_fsolver_x_lower (gsl_solver.get());
    phic_high = gsl_root_fsolver_x_upper (gsl_solver.get());
    status = gsl_root_test_interval (phic_low, phic_high,
                                     epsabs, epsrel);
  } while (status == GSL_CONTINUE && iter < max_iter);

  return status;

};

double phi0_error_gsl(double phic, void *vparams)
{
  Phi0Params *params = (Phi0Params*) vparams;

  params->model0->setphic(phic);
  params->model0->solve();

  return params->model0->phi_0() - params->phi0_goal;
};
