#include "ppEOS.h"
#include "ppEOSTable.h"
#include <cmath>
#include <iostream>
#include <fstream>

#include "testppEOS-cmdline.h"

int parseCmdlineInto(struct gengetopt_args_info *args_info,
			   int argc, char *argv[]);

int main(int argc, char *argv[]) {

  // Useful logging
  std::cerr << "This is " __FILE__ ", compiled on " __DATE__ << std::endl;

  ////////////////////////////////////////////////////////////
  // Parse command line and config file
  ////////////////////////////////////////////////////////////

  static gengetopt_args_info args_info;

  if (0 != parseCmdlineInto(&args_info, argc, argv) )
  {
    // Error messages should have already been emitted by parseCmdlineInto
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

  const int num_points = args_info.num_points_arg;

  if (num_points < 2)
  {
    std::cerr << "ERROR: num-points *must* be at least two!" << std::endl;
    return 1;
  };

  const double rho_low = args_info.rho_low_arg;
  const double rho_high = args_info.rho_high_arg;

  // This will be log spaced
  const double log10_rho_low = log10(rho_low);
  const double log10_rho_high = log10(rho_high);
  const double delta = (log10_rho_high-log10_rho_low)/(num_points-1);

  std::ofstream o(args_info.out_arg);

  std::cerr << "Writing to file " << args_info.out_arg << std::endl;

  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#    1  mass density   [g/cm^3]\n"
    "#    2  pressure       [dyne/cm^2]\n"
    "#    3  energy density [g/cm^3]\n";

  o << header;

  const char s = ' ';

  for (unsigned int i=0; i<num_points; i++) {
    double rho = pow( 10., log10_rho_low + i*delta );
    o << rho << s <<  EoS.pressure(rho)
      << s << EoS.epsilon(rho) << std::endl;
  };

  o.close();

  // Exit cleanly
  cmdline_parser_free(&args_info);
  return 0;
};

//////////////////////////////////////////////////////////////////////
// Command line and config file parsing
//////////////////////////////////////////////////////////////////////
int parseCmdlineInto(struct gengetopt_args_info *args_info,
			    int argc, char *argv[])
{
  int result = 0;

  /* call the command line parser */
  if (cmdline_parser (argc, argv, args_info) != 0) {
    result = 1;
    goto stop;
  };

  std::cerr << "Arguments on command line:" << std::endl;
  cmdline_parser_dump(stderr, args_info);

 stop:
  return result;

};
