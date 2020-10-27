/*
 * Conformal.hpp
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2020 Oct
 *
 * Pure virtual base class that gives interface for conformal factor and log deriv
 *
 */

#pragma once

class Conformal
{
public:
  virtual double A(double phi) const = 0;
  virtual double alpha(double phi) const = 0;
};
