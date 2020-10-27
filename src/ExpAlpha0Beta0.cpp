/*
 * ExpAlpha0Beta0.cpp
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2020 Oct
 *
 * Implementation of ExpAlpha0Beta0,
 * A(phi) = exp(alpha_0*phi + 0.5*beta_0*phi^2)
 *
 */

#include "ExpAlpha0Beta0.hpp"
#include <math.h>

double ExpAlpha0Beta0::A(double phi) const
{
  return exp( (alpha_0+0.5*beta_0*phi)*phi );
};

double ExpAlpha0Beta0::alpha(double phi) const
{
  return alpha_0 + beta_0*phi;
};
