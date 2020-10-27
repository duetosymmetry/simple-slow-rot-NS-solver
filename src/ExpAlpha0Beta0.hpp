/*
 * ExpAlpha0Beta0.hpp
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2020 Oct
 *
 * Class deriving the interface of Conformal,
 * implementing A(phi) = exp(alpha_0*phi + 0.5*beta_0*phi^2)
 *
 */

#pragma once

#include "Conformal.hpp"

class ExpAlpha0Beta0 : public Conformal
{
public:

  const double alpha_0;
  const double beta_0;

  ExpAlpha0Beta0(double alpha_0, double beta_0)
    : alpha_0(alpha_0), beta_0(beta_0) {};

  double A(double phi) const;
  double alpha(double phi) const;
};
