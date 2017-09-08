/*
 * backgroundModelTest.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * This program just runs a BackgroundModel in order to test said class.
 *
 */

#include "constants.h" /* for G and c */
#include "ppEOSTable.h"
#include "BackgroundModel.h"
#include "writeModels.h"
#include <cmath>
#include <iostream>

#include "backgroundModelTest-cmdline.h"

#define CONFIG_FILENAME "backgroundModel.conf"

int main(int argc, char *argv[])
{

  ppEOS EoS = findEOS("SLy");
  double pc;
  
  ////////////////////////////////////////////////////////////
  // Parse command line and config file
  ////////////////////////////////////////////////////////////

  static gengetopt_args_info args_info;

  struct cmdline_parser_params *params;

  /* initialize the parameters structure */
  params = cmdline_parser_params_create();
  
  /* 
     initialize args_info, but don't check for required options
     NOTICE: the other fields are initialized to their default values
  */
  params->check_required = 0;

  /* call the config file parser */
  if (cmdline_parser_config_file
      (CONFIG_FILENAME, &args_info, params) != 0)
    exit(1);

  /* 
     override config file options,
     do not initialize args_info, don't check for required options.
  */
  params->initialize = 0;
  params->override = 1;
  params->check_required = 1;

  /* call the command line parser */
  if (cmdline_parser_ext (argc, argv, &args_info, params) != 0)
    exit(1);

  /* make sure everything required was specified somewhere */
  if (cmdline_parser_required(&args_info, argv[0]) != 0)
    exit(1);

  // Extra info
  std::cerr << "This is " __FILE__ ", compiled on " __DATE__ << std::endl;
  std::cerr << "Read arguments of " << CONFIG_FILENAME << std::endl;
  std::cerr << "Arguments in config file and on command line:" << std::endl;
  cmdline_parser_dump(stderr, &args_info);

  ////////////////////////////////////////////////////////////
  // Put command line parameters into EoS

  if (args_info.eos_name_given) {
    EoS = findEOS(args_info.eos_name_arg);
  } else {
    EoS = ppEOS( pow(10., args_info.log10p1_arg),
                 args_info.Gamma1_arg,
                 args_info.Gamma2_arg,
                 args_info.Gamma3_arg );
  };

  // And into pc
  // This starts out in dyne/cm^2
  pc = args_info.pc_arg;
  // Convert to geometric units
  pc *= G_cgs/(c_cm_s*c_cm_s*c_cm_s*c_cm_s);

  ////////////////////////////////////////////////////////////
  // Create the model and solve
  BackgroundModel model0( EoS, pc );
  
  model0.solve();
  
  writeBackgroundModel( args_info.out_arg, model0 );

  return 0;
}
