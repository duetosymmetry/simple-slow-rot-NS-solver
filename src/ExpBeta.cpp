/*
 * ExpBeta.cpp
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2020 Oct
 *
 * Implementations of class Conformal, A(phi) = exp(0.5*beta*phi^2)
 *
 */

#include "ExpBeta.hpp"
#include <math.h>

double ExpBeta::A(double phi) const
{
  return exp(0.5*beta*phi*phi);
};

double ExpBeta::alpha(double phi) const
{
  return beta*phi;
};
