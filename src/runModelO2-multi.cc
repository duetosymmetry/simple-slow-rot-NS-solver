/*
 * runModelO2-multi.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Feb
 *
 * This program runs multiple ModelO2 models
 *
 */

#include "constants.h" /* for G and c */
#include "ppEOSTable.h"
#include "BackgroundModel.h"
#include "ModelO1.h"
#include "ModelO2.h"
#include "writeModels.h"
#include "util.h"
#include <vector>
#include <gsl/gsl_spline.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>

#include "runModelO2-multi-cmdline.h"

int parseCmdlineConfigInto(struct gengetopt_args_info *args_info,
                           int argc, char *argv[]);

int main(int argc, char *argv[])
{

  // Useful logging
  std::cerr << "This is " __FILE__ ", compiled on " __DATE__ << std::endl;

  ////////////////////////////////////////////////////////////
  // Parse command line and config file
  ////////////////////////////////////////////////////////////

  static gengetopt_args_info args_info;

  if (0 != parseCmdlineConfigInto(&args_info, argc, argv) )
  {
    // Error messages should have already been emitted by parseCmdlineConfigInto
    cmdline_parser_free(&args_info);

    // Exit with error
    return 1;
  };

  ////////////////////////////////////////////////////////////
  // Put command line parameters into EoS, pc, etc.

  ppEOS EoS = findEOS("SLy");

  if (args_info.eos_name_given) {
    EoS = findEOS(args_info.eos_name_arg);
  } else {
    EoS = ppEOS( pow(10., args_info.log10p1_arg),
                 args_info.Gamma1_arg,
                 args_info.Gamma2_arg,
                 args_info.Gamma3_arg,
                 pow(10., args_info.log10rho1_arg) );
  };

  const int num = args_info.num_arg;
  if (num < 2)
  {
    std::cerr << "ERROR: num *must* be at least two!" << std::endl;
    return 1;
  };

  // These come in as dyne/cm^2, convert them to geometric units
  const double pc_low  = args_info.pc_low_arg  * G_cgs/(c_cm_s*c_cm_s*c_cm_s*c_cm_s);
  const double pc_high = args_info.pc_high_arg * G_cgs/(c_cm_s*c_cm_s*c_cm_s*c_cm_s);

  if (!(pc_low > 0) || !(pc_high > 0))
  {
    std::cerr << "ERROR: pressures must be positive!" << std::endl;
    return 1;
  };

  if (!(pc_low < pc_high))
  {
    std::cerr << "ERROR: pc_low must be smaller than pc_high!" << std::endl;
    return 1;
  };

  // This will be log spaced
  const double log10_pc_low  = log10(pc_low);
  const double log10_pc_high = log10(pc_high);
  const double delta = (log10_pc_high-log10_pc_low)/(num-1);

  // summary file
  std::ofstream o( args_info.out_arg );
  std::cerr << "Writing to file " << args_info.out_arg << std::endl;
  writeSummaryHeader2( o );

  ////////////////////////////////////////////////////////////
  // Storage for M/R/C interpolations
  // We invert R so it's monotonically increasing with pc
  std::vector<double>
    pcs(num),
    Ms(num),
    invRs(num),
    Cs(num),
    invQbars(num),
    invIbars(num);

  ////////////////////////////////////////////////////////////
  // Model variables
  double pc =  pc_low;
  BackgroundModel model0( EoS, pc );
  ModelO1 model1( model0 );
  ModelO2 model2( model1 );

  for ( unsigned int i=0; i < num; i++ ) {
    pc = pow( 10., log10_pc_low + i*delta );

    std::cerr << "Model #" << i+1 << "/" << num << ", pc=" << pc << " ... ";

    runBG12( pc, model0, model1, model2 );

    // Save this stuff for spline building later
    pcs  [i] = pc;
    Ms   [i] = model0.M();
    invRs[i] = 1. / model0.R();
    Cs   [i] = model0.M() /  model0.R();
    invQbars[i] = 1. / model2.Qbar();
    invIbars[i] = 1. / model1.Ibar();

    std::cerr << "writing summary line." << std::endl;
    writeSummaryLine( o, model0, model1, model2 );
  };

  o.close();

  ////////////////////////////////////////////////////////////
  // Running extra models

  // Build splines
  const int
    nIncM = numMonotonicallyIncreasing(Ms),
    nIncR = numMonotonicallyIncreasing(invRs),
    nIncC = numMonotonicallyIncreasing(Cs),
    nIncQ = numMonotonicallyIncreasing(invQbars),
    nIncI = numMonotonicallyIncreasing(invQbars);

  gsl_interp_accel * acc = gsl_interp_accel_alloc();
  gsl_spline
    * Mspline       = gsl_spline_alloc (gsl_interp_cspline, nIncM),
    * invRspline    = gsl_spline_alloc (gsl_interp_cspline, nIncR),
    * Cspline       = gsl_spline_alloc (gsl_interp_cspline, nIncC),
    * invQbarspline = gsl_spline_alloc (gsl_interp_cspline, nIncQ),
    * invIbarspline = gsl_spline_alloc (gsl_interp_cspline, nIncI);

  if (args_info.extra_M_given > 0) {
    std::cerr << "Building M spline ... ";
    gsl_spline_init (Mspline, Ms.data(), pcs.data(), nIncM);
  } else
    std::cerr << "No extra M models requested, skipping M spline ... "
              << std::endl;

  if (args_info.extra_R_given > 0) {
    std::cerr << "Building R spline ... ";
    gsl_spline_init (invRspline, invRs.data(), pcs.data(), nIncR);
  } else
    std::cerr << "No extra R models requested, skipping R spline ... "
              << std::endl;

  if (args_info.extra_C_given > 0) {
    std::cerr << "Building C spline ... ";
    gsl_spline_init (Cspline, Cs.data(), pcs.data(), nIncC);
  } else
    std::cerr << "No extra C models requested, skipping C spline ... "
              << std::endl;

  if (args_info.extra_Qbar_given > 0) {
    std::cerr << "Building Qbar spline ... ";
    gsl_spline_init (invQbarspline, invQbars.data(), pcs.data(), nIncQ);
  } else
    std::cerr << "No extra Qbar models requested, skipping Qbar spline ... "
              << std::endl;

  if (args_info.extra_Ibar_given > 0) {
    std::cerr << "Building Ibar spline ... ";
    gsl_spline_init (invIbarspline, invIbars.data(), pcs.data(), nIncI);
  } else
    std::cerr << "No extra Ibar models requested, skipping Qbar spline ... "
              << std::endl;

  // queue for the extra models
  std::queue< extra_info > extra_q;

  ////////////////////////////////////////////////////////////
  // Enqueue pc
  if (args_info.extra_pc_given) {
    std::vector<double>
      extra_pc_vec( args_info.extra_pc_arg,
                    args_info.extra_pc_arg + args_info.extra_pc_given);
    enqueue_extra( extra_q,
                   extra_pc_vec, extra_pc_vec,
                   0., INFINITY,
                   nullptr, nullptr, "pc");
  };

  ////////////////////////////////////////////////////////////
  // Enqueue M
  if (args_info.extra_M_given)
  {
  std::vector<double>
    extra_M_vec( args_info.extra_M_arg,
                 args_info.extra_M_arg + args_info.extra_M_given);
  // The argument comes in as solar masses.
  // Convert to cm.
  std::vector<double> extra_M_cm(args_info.extra_M_given);
  for (unsigned int i=0; i < args_info.extra_M_given; i++)
    extra_M_cm[i] = extra_M_vec[i] * GMsun_cm;

  enqueue_extra( extra_q,
                 extra_M_vec, extra_M_cm,
                 Ms[0], Ms[nIncM-1],
                 Mspline, acc, "M");
  };
  ////////////////////////////////////////////////////////////
  // Enqueue R
  if (args_info.extra_R_given)
  {
  std::vector<double>
    extra_R_vec( args_info.extra_R_arg,
                 args_info.extra_R_arg + args_info.extra_R_given);
  // The argument comes in as km.
  // Convert to cm, then take reciprocal to interpolate
  std::vector<double> extra_invR(args_info.extra_R_given);
  for (unsigned int i=0; i < args_info.extra_R_given; i++)
    extra_invR[i] = 1.e-5 / extra_R_vec[i];

  enqueue_extra( extra_q,
                 extra_R_vec,
                 extra_invR,
                 invRs[0], invRs[nIncR-1],
                 invRspline, acc, "R");
  };
  ////////////////////////////////////////////////////////////
  // Enqueue C
  if (args_info.extra_C_given) {
  std::vector<double>
    extra_C_vec( args_info.extra_C_arg,
                 args_info.extra_C_arg + args_info.extra_C_given);
  enqueue_extra( extra_q,
                 extra_C_vec, extra_C_vec,
                 Cs[0], Cs[nIncC-1],
                 Cspline, acc, "C");
  };

  ////////////////////////////////////////////////////////////
  // Enqueue Qbars
  if (args_info.extra_Qbar_given)
  {
  std::vector<double>
    extra_Qbar_vec( args_info.extra_Qbar_arg,
                    args_info.extra_Qbar_arg + args_info.extra_Qbar_given);
  // The argument comes in as Qbars.
  // Take reciprocal to interpolate
  std::vector<double> extra_invQ(args_info.extra_Qbar_given);
  for (unsigned int i=0; i < args_info.extra_Qbar_given; i++)
    extra_invQ[i] = 1. / extra_Qbar_vec[i];

  enqueue_extra( extra_q,
                 extra_Qbar_vec, extra_invQ,
                 invQbars[0], invQbars[nIncQ-1],
                 invQbarspline, acc, "Qbar");
  };
  ////////////////////////////////////////////////////////////
  // Enqueue Ibars
  if (args_info.extra_Ibar_given)
  {
  std::vector<double>
    extra_Ibar_vec( args_info.extra_Ibar_arg,
                    args_info.extra_Ibar_arg + args_info.extra_Ibar_given);
  // The argument comes in as Ibars.
  // Take reciprocal to interpolate
  std::vector<double> extra_invI(args_info.extra_Ibar_given);
  for (unsigned int i=0; i < args_info.extra_Ibar_given; i++)
    extra_invI[i] = 1. / extra_Ibar_vec[i];

  enqueue_extra( extra_q,
                 extra_Ibar_vec, extra_invI,
                 invIbars[0], invIbars[nIncI-1],
                 invIbarspline, acc, "Ibar");
  };
  ////////////////////////////////////////////////////////////
  // Run extra models
  unsigned int i=0;
  const unsigned int num_valid_extra = extra_q.size();
  while (!extra_q.empty()) {

    i++;

    extra_info info = extra_q.front();
    extra_q.pop();

    pc = info.pc;

    std::cerr << "Running extra model #" << i << "/" << num_valid_extra
        << ", " << info.prop_name << "=" << info.prop_value
        << std::endl;

    runBG12( pc, model0, model1, model2 );

    std::stringstream extra_filename;

    extra_filename << args_info.out_extra_basename_arg
       << "-" << info.prop_name
       << "-" << info.prop_value << ".dat";

    std::cerr << "Writing BG+O(1)+O(2) to file "
              << extra_filename.str() << std::endl;
    writeBG12( extra_filename.str().c_str(), model2 );

  };

  // Exit cleanly
  gsl_spline_free (invIbarspline);
  gsl_spline_free (invQbarspline);
  gsl_spline_free (Cspline);
  gsl_spline_free (invRspline);
  gsl_spline_free (Mspline);
  gsl_interp_accel_free (acc);

  cmdline_parser_free(&args_info);
  return 0;
}

//////////////////////////////////////////////////////////////////////
// Command line and config file parsing
//////////////////////////////////////////////////////////////////////
int parseCmdlineConfigInto(struct gengetopt_args_info *args_info,
          int argc, char *argv[])
{
  int result = 0;

  struct cmdline_parser_params *params;

  /* initialize the parameters structure */
  params = cmdline_parser_params_create();

  /*
     initialize args_info, but don't check for required options.
     override defaults.
  */
  params->initialize = 1;
  params->check_required = 0;
  params->override = 1;

  /* call the command line parser */
  if (cmdline_parser_ext (argc, argv, args_info, params) != 0) {
    result = 1;
    goto stop;
  };

  /*
     DO NOT override command line options.
     do not initialize args_info.
     check for required options.
     NOTICE: we must NOT skip the 0 assignment to initialize,
     since its default value is 1.
  */
  params->initialize = 0;
  params->override = 0;
  params->check_required = 1;

  /* call the config file parser */

  std::cerr << "Reading config file " << args_info->conf_file_arg << std::endl;

  if (cmdline_parser_config_file
      (args_info->conf_file_arg, args_info, params) != 0)
    {
      result = 1;
      goto stop;
    }

  /* make sure everything required was specified somewhere */
  if (cmdline_parser_required(args_info, argv[0]) != 0)
  {
    result = 1;
    goto stop;
  };

  std::cerr << "Arguments in config file and on command line:" << std::endl;
  cmdline_parser_dump(stderr, args_info);

 stop:
  free(params);
  return result;

};
