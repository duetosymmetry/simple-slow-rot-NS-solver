/*
 * ExpBeta.hpp
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2020 Oct
 *
 * Instance of class Conformal, implementing A(phi) = exp(0.5*beta*phi^2)
 *
 */

#pragma once

#include "Conformal.hpp"

class ExpBeta : public Conformal
{
public:

  const double beta;

  ExpBeta(double beta) : beta(beta) {};

  double A(double phi) const;
  double alpha(double phi) const;
};
