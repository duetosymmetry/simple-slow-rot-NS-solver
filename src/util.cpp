/*
 * util.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Feb
 *
 * Implementation for utility (helper) functions
 *
 */

#include "util.hpp"
#include "BackgroundModel.hpp"
#include "ModelO1.hpp"
#include "ModelO2.hpp"
#include <gsl/gsl_spline.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <queue>
#include <iostream>

// Runs the BG, O(1), and O(2) models.
void runBG12( double pc, BackgroundModel &model0,
              ModelO1 &model1, ModelO2 &model2 )
{
  //////////////////////////////
  // Background model
  // This also has the side effect of resetting model0
  model0.setpc( pc );
  std::cerr << "Solving at O(a^0) ... ";
  model0.solve();
  //////////////////////////////
  // O(1) model
  model1.reset();
  std::cerr << "solving at O(a^1) ... ";
  model1.solve();
  //////////////////////////////
  // O(2) model
  model2.reset();
  std::cerr << "solving at O(a^2) ... ";
  model2.solve();
};

////////////////////////////////////////////////////////////
// Enqueue array of extra args with a given spline
void enqueue_extra( std::queue< extra_info > & extra_q,
                    const std::vector<double> & extra_arg,
                    const std::vector<double> & extra_arg_lookup,
                    double arg_min, double arg_max,
                    gsl_spline * spline, gsl_interp_accel *acc,
                    const std::string &label )
{
  if (extra_arg.size() != extra_arg_lookup.size()) {
    std::cerr <<
      "WARNING: extra_arg and extra_arg_lookup are "
      "different lengths, ignoring everything!" << std::endl;
    return;
  };

  for (unsigned int i=0; i < extra_arg.size(); i++) {
    const double arg = extra_arg[i];
    const double arg_lookup = extra_arg_lookup[i];
    if ((arg_lookup < arg_min) || (arg_lookup > arg_max)) {
      std::cerr << "WARNING: extra " << label << " #" << i
                << " is out of range. Skipping." << std::endl;
      continue;
    };
    double pc;
    // interpolate, unless spline is null.
    pc = spline ? gsl_spline_eval(spline, arg_lookup, acc) : arg;
    extra_q.push((extra_info){ pc, label, arg });
  };
};
