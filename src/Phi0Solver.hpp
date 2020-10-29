/*
 * Phi0Solver.hpp
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2020 Oct
 *
 * Class for controlling solving a BG model with a desired phi0 instead of phic
 *
 */

#pragma once

#include <memory>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "BackgroundModel.hpp"

class Phi0Solver
{
public:

  Phi0Solver( double epsabs, double epsrel, int max_iter );
  int solve( BackgroundModel &model0, double phi0_goal,
             double phic_min, double phic_max );

  void setParams( double new_epsabs, double new_epsrel, int new_max_iter )
  {
    epsabs   = new_epsabs;
    epsrel   = new_epsrel;
    max_iter = new_max_iter;
  };

  double getEpsAbs() { return epsabs; };
  double getEpsRel() { return epsrel; };
  int getMaxIter() { return max_iter; };
  double getPhi0() { return phi0; };

private:

  double epsabs, epsrel;
  int max_iter;

  double phi0;

  std::unique_ptr<gsl_root_fsolver, decltype(&gsl_root_fsolver_free)> gsl_solver;
  gsl_function F;

};
