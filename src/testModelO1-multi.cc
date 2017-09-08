/*
 * testModelO1-multi.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Feb
 *
 * This program run multiple ModelO1 models
 *
 */

#include "constants.h" /* for G and c */
#include "ppEOSTable.h"
#include "BackgroundModel.h"
#include "ModelO1.h"
#include "writeModels.h"
#include <cmath>
#include <iostream>
#include <fstream>

#include "testModelO1-multi-cmdline.h"

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
     args_info.Gamma3_arg );
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

  // This will be log spaced
  const double log10_pc_low  = log10(pc_low);
  const double log10_pc_high = log10(pc_high);
  const double delta = (log10_pc_high-log10_pc_low)/(num-1);

  // summary file
  std::ofstream o(args_info.out_arg);

  std::cerr << "Writing to file " << args_info.out_arg << std::endl;

  writeSummaryHeader1( o );

  ////////////////////////////////////////////////////////////
  // Main variables
  double pc =  pc_low;
  BackgroundModel model0( EoS, pc );
  ModelO1 model1( model0 );

  for ( unsigned int i=0; i < num; i++ )
  {
    pc = pow( 10., log10_pc_low + i*delta );

    //////////////////////////////
    // Background model
    // This also has the side effect of resetting model0
    model0.setpc( pc );

    std::cerr << "Model #" << i+1 << "/" << num << ", pc=" << pc << " ... ";

    model0.solve();

    //////////////////////////////
    // O(1) model
    model1.reset();
    std::cerr << "solving at O(a^1) ... ";

    model1.solve();

    std::cerr << "writing summary line." << std::endl;

    writeSummaryLine( o, model0, model1 );

  };


  o.close();

  // Exit cleanly
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
