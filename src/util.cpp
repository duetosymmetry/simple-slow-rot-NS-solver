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
#include <gsl/gsl_spline.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <queue>
#include <iostream>

// Runs the BG model
void runBG( double pc, double phic,
            BackgroundModel &model0 )
{
  //////////////////////////////
  // Background model
  // This has the side effect of resetting model0
  model0.setpc( pc );
  model0.setphic( phic );
  std::cerr << "Solving at O(a^0) ... ";
  model0.solve();
};

////////////////////////////////////////////////////////////
// Enqueue array of extra args with a given spline
void enqueue_extra( std::queue< extra_info > & extra_q,
                    const std::vector<double> & extra_arg,
                    const std::vector<double> & extra_arg_lookup,
                    double arg_min, double arg_max,
                    gsl_spline * pc_spline,
                    gsl_spline * phic_spline,
                    gsl_interp_accel *acc,
                    const std::string &label )
{
  if (extra_arg.size() != extra_arg_lookup.size()) {
    std::cerr <<
      "WARNING: extra_arg and extra_arg_lookup are "
      "different lengths, ignoring everything!" << std::endl;
    return;
  };

  for (unsigned int i=0; i < extra_arg.size(); i++) {
    // TODO CURRENTLY ONLY ONE ARG?
    const double arg = extra_arg[i];
    const double arg_lookup = extra_arg_lookup[i];
    if ((arg_lookup < arg_min) || (arg_lookup > arg_max)) {
      std::cerr << "WARNING: extra " << label << " #" << i
                << " is out of range. Skipping." << std::endl;
      continue;
    };
    double pc, phic;
    // interpolate, unless spline is null.
    pc = pc_spline ? gsl_spline_eval(pc_spline, arg_lookup, acc) : arg;
    // interpolate, unless spline is null.
    phic = phic_spline ? gsl_spline_eval(phic_spline, arg_lookup, acc) : arg;
    extra_q.push((extra_info){ pc, phic, label, arg });
  };
};
