/*
 * testModelO1-single.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * This program runs a single ModelO1 in order to test said class.
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

#include "testModelO1-single-cmdline.h"

int parseCmdlineConfigInto(struct gengetopt_args_info *args_info,
			   int argc, char *argv[]);

int main(int argc, char *argv[])
{

  ppEOS EoS = findEOS("SLy");
  double pc;

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
  // Create the bg model and solve
  BackgroundModel model0( EoS, pc );

  std::cerr << "Solving background model" << std::endl;

  model0.solve();

  std::cerr << "Writing background model to file " << args_info.out_0_arg << std::endl;

  // Can now write bg model to file
  writeBackgroundModel( args_info.out_0_arg, model0 );

  ////////////////////////////////////////////////////////////
  // Creat the O(a^1) model and solve
  ModelO1 model1( model0 );

  std::cerr << "Solving model at O(a^1)" << std::endl;

  model1.solve();

  std::cerr << "Writing model at O(a^1) to file " << args_info.out_1_arg << std::endl;

  // Can now write bg model to file
  writeModelO1( args_info.out_1_arg, model1 );

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
