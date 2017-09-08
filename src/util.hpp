/*
 * util.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Feb
 *
 * Declarations of utility (helper) functions useful for running models
 *
 */

#pragma once

#include "BackgroundModel.hpp"
#include "ModelO1.hpp"
#include "ModelO2.hpp"

#include <gsl/gsl_spline.h>

#include <functional>
#include <vector>
#include <algorithm>
#include <queue>

// Convenience function: runs the BG, O(1), and O(2) models.
void runBG12( double pc, BackgroundModel &model0,
              ModelO1 &model1, ModelO2 &model2 );

// Determine the number of monotonically increasing entries in an
// array, starting at the beginning.
template< class T >
int numMonotonicallyIncreasing( const std::vector<T> & vals )
{
  auto firstGE =
    std::adjacent_find(vals.begin(), vals.end(), std::greater_equal<T>());
  return std::distance(vals.begin(), firstGE);
};

// This info struct is used for queueing the extra models
struct extra_info
{
  double pc;
  std::string prop_name;
  double prop_value;
};

// Add some extra runs into a queue
void enqueue_extra( std::queue< extra_info > & extra_q,
                    const std::vector<double> & extra_arg,
                    const std::vector<double> & extra_arg_lookup,
                    double arg_min, double arg_max,
                    gsl_spline * spline, gsl_interp_accel *acc,
                    const std::string &label );
