/*
 * runPcPhi0Grid.cpp
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2020 Oct
 *
 * This program runs a grid of stellar models at different values of phi0
 *
 */

#include "integMagicNums.hpp"
#include "constants.hpp" /* for G and c */
#include "ppEOSTable.hpp"
#include "BackgroundModel.hpp"
#include "ExpAlpha0Beta0.hpp"
#include "Phi0Solver.hpp"
#include "writeModels.hpp"
#include "util.hpp"
#include <vector>
#include <gsl/gsl_spline.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>

#include "runPcPhi0Grid-cmdline.h"

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
  const double alpha_0 = args_info.alpha_arg;
  const double beta_0 = args_info.beta_arg;
  ExpAlpha0Beta0 conf(alpha_0, beta_0);

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

  //////////////////////////////
  // pc grid
  const int n_p = args_info.n_p_arg;
  if (n_p < 1)
  {
    std::cerr << "ERROR: n-p *must* be at least one!" << std::endl;
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
  const double delta_pc = (n_p > 1) ? (log10_pc_high-log10_pc_low)/(n_p-1.) : 0.;

  //////////////////////////////
  // phic grid
  const int n_phi0 = args_info.n_phi0_arg;
  if (n_phi0 < 1)
  {
    std::cerr << "ERROR: n-phi0 *must* be at least one!" << std::endl;
    return 1;
  };

  const double phic_low  = args_info.phic_low_arg;
  const double phic_high = args_info.phic_high_arg;

  if (!(phic_low < phic_high))
  {
    std::cerr << "ERROR: phic_low must be smaller than phic_high!" << std::endl;
    return 1;
  };

  const double phi0_low  = args_info.phi0_low_arg;
  const double phi0_high = args_info.phi0_high_arg;

  if (!(phi0_low < phi0_high))
  {
    std::cerr << "ERROR: phi0_low must be smaller than phi0_high!" << std::endl;
    return 1;
  };

  // This will be linearly spaced
  const double delta_phi0 = (n_phi0 > 1) ? (phi0_high-phi0_low)/(n_phi0-1.) : 0.;

  //////////////////////////////
  // summary file
  std::ofstream o( args_info.out_arg );
  std::cerr << "Writing to file " << args_info.out_arg << std::endl;
  writeSummaryHeader( o );

  ////////////////////////////////////////////////////////////
  // Model variables
  double pc =  pc_low;
  double phi0 =  phi0_low;
  BackgroundModel model0( EoS, conf, pc, phi0 );
  Phi0Solver phi0solver( EPSABS, EPSREL, MAX_ITER );

  unsigned int i=0;
  const unsigned int tot = n_p*n_phi0;

  for ( unsigned int ip=0; ip < n_p; ip++ ) {
    pc = pow( 10., log10_pc_low + ip*delta_pc );
    model0.setpc(pc);

    for ( unsigned int iphi=0; iphi < n_phi0; iphi++, i++ ) {
      phi0 = phi0_low + iphi*delta_phi0;

      std::cerr << "Model #" << i+1 << "/" << tot
                << ", (" << ip+1 << "/" << n_p
                << ", " << iphi+1 << "/" << n_phi0
                << "), pc=" << pc
                << ", phi0=" << phi0 << " ... ";

      int ret = phi0solver.solve(model0, phi0, phic_low, phic_high);
      std::cerr << "Returned " << ret << " from phi0solver" << std::endl;

      std::cerr << "writing summary line." << std::endl;
      writeSummaryLine( o, model0 );
    };
  };

  o.close();

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
