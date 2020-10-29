/** @file runPhi0Family-cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.23
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt */

#ifndef RUNPHI0FAMILY_CMDLINE_H
#define RUNPHI0FAMILY_CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "runPhi0Family"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "runPhi0Family"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0.1"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * conf_file_arg;	/**< @brief Configuration file to read (default='phi0Family.conf').  */
  char * conf_file_orig;	/**< @brief Configuration file to read original value given at command line.  */
  const char *conf_file_help; /**< @brief Configuration file to read help description.  */
  double pc_low_arg;	/**< @brief Smallest central pressure [dyne/cm^2] (default='4.e+34').  */
  char * pc_low_orig;	/**< @brief Smallest central pressure [dyne/cm^2] original value given at command line.  */
  const char *pc_low_help; /**< @brief Smallest central pressure [dyne/cm^2] help description.  */
  double pc_high_arg;	/**< @brief Largest  central pressure [dyne/cm^2] (default='2.e+36').  */
  char * pc_high_orig;	/**< @brief Largest  central pressure [dyne/cm^2] original value given at command line.  */
  const char *pc_high_help; /**< @brief Largest  central pressure [dyne/cm^2] help description.  */
  int n_p_arg;	/**< @brief Number of pc values in grid (pressure is log-spaced) (default='16').  */
  char * n_p_orig;	/**< @brief Number of pc values in grid (pressure is log-spaced) original value given at command line.  */
  const char *n_p_help; /**< @brief Number of pc values in grid (pressure is log-spaced) help description.  */
  double phi0_arg;	/**< @brief Asymptotic scalar phi_0 (dimensionless) (default='0.0043').  */
  char * phi0_orig;	/**< @brief Asymptotic scalar phi_0 (dimensionless) original value given at command line.  */
  const char *phi0_help; /**< @brief Asymptotic scalar phi_0 (dimensionless) help description.  */
  double phic_low_arg;	/**< @brief Lower search bound for central value of phi (dimensionless) (default='-0.3').  */
  char * phic_low_orig;	/**< @brief Lower search bound for central value of phi (dimensionless) original value given at command line.  */
  const char *phic_low_help; /**< @brief Lower search bound for central value of phi (dimensionless) help description.  */
  double phic_high_arg;	/**< @brief Upper search bound for central value of phi (dimensionless) (default='0.3').  */
  char * phic_high_orig;	/**< @brief Upper search bound for central value of phi (dimensionless) original value given at command line.  */
  const char *phic_high_help; /**< @brief Upper search bound for central value of phi (dimensionless) help description.  */
  double alpha_arg;	/**< @brief alpha_0 (dimensionless) (default='0.').  */
  char * alpha_orig;	/**< @brief alpha_0 (dimensionless) original value given at command line.  */
  const char *alpha_help; /**< @brief alpha_0 (dimensionless) help description.  */
  double beta_arg;	/**< @brief alpha_0 (dimensionless) (default='-4.5').  */
  char * beta_orig;	/**< @brief alpha_0 (dimensionless) original value given at command line.  */
  const char *beta_help; /**< @brief alpha_0 (dimensionless) help description.  */
  char * eos_name_arg;	/**< @brief Name of the EoS model. (default='SLy').  */
  char * eos_name_orig;	/**< @brief Name of the EoS model. original value given at command line.  */
  const char *eos_name_help; /**< @brief Name of the EoS model. help description.  */
  double log10p1_arg;	/**< @brief The log (base 10) of (p1 in dyne/cm^2). (default='34.384').  */
  char * log10p1_orig;	/**< @brief The log (base 10) of (p1 in dyne/cm^2). original value given at command line.  */
  const char *log10p1_help; /**< @brief The log (base 10) of (p1 in dyne/cm^2). help description.  */
  double Gamma1_arg;	/**< @brief Power law index Gamma1 (default='3.005').  */
  char * Gamma1_orig;	/**< @brief Power law index Gamma1 original value given at command line.  */
  const char *Gamma1_help; /**< @brief Power law index Gamma1 help description.  */
  double Gamma2_arg;	/**< @brief Power law index Gamma2 (default='2.988').  */
  char * Gamma2_orig;	/**< @brief Power law index Gamma2 original value given at command line.  */
  const char *Gamma2_help; /**< @brief Power law index Gamma2 help description.  */
  double Gamma3_arg;	/**< @brief Power law index Gamma3 (default='2.851').  */
  char * Gamma3_orig;	/**< @brief Power law index Gamma3 original value given at command line.  */
  const char *Gamma3_help; /**< @brief Power law index Gamma3 help description.  */
  double log10rho1_arg;	/**< @brief The log (base 10) of (rho1 in g/cm^3). (default='14.7').  */
  char * log10rho1_orig;	/**< @brief The log (base 10) of (rho1 in g/cm^3). original value given at command line.  */
  const char *log10rho1_help; /**< @brief The log (base 10) of (rho1 in g/cm^3). help description.  */
  char * out_arg;	/**< @brief The file name to output the multiple model summary (default='phi0-family-summary.dat').  */
  char * out_orig;	/**< @brief The file name to output the multiple model summary original value given at command line.  */
  const char *out_help; /**< @brief The file name to output the multiple model summary help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int conf_file_given ;	/**< @brief Whether conf-file was given.  */
  unsigned int pc_low_given ;	/**< @brief Whether pc-low was given.  */
  unsigned int pc_high_given ;	/**< @brief Whether pc-high was given.  */
  unsigned int n_p_given ;	/**< @brief Whether n-p was given.  */
  unsigned int phi0_given ;	/**< @brief Whether phi0 was given.  */
  unsigned int phic_low_given ;	/**< @brief Whether phic-low was given.  */
  unsigned int phic_high_given ;	/**< @brief Whether phic-high was given.  */
  unsigned int alpha_given ;	/**< @brief Whether alpha was given.  */
  unsigned int beta_given ;	/**< @brief Whether beta was given.  */
  unsigned int eos_name_given ;	/**< @brief Whether eos-name was given.  */
  unsigned int log10p1_given ;	/**< @brief Whether log10p1 was given.  */
  unsigned int Gamma1_given ;	/**< @brief Whether Gamma1 was given.  */
  unsigned int Gamma2_given ;	/**< @brief Whether Gamma2 was given.  */
  unsigned int Gamma3_given ;	/**< @brief Whether Gamma3 was given.  */
  unsigned int log10rho1_given ;	/**< @brief Whether log10rho1 was given.  */
  unsigned int out_given ;	/**< @brief Whether out was given.  */

  int param_fit_mode_counter; /**< @brief Counter for mode param_fit */
  int param_manual_mode_counter; /**< @brief Counter for mode param_manual */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * The config file parser (deprecated version)
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_config_file() instead
 */
int cmdline_parser_configfile (const char *filename,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The config file parser
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_config_file (const char *filename,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);

extern const char *cmdline_parser_eos_name_values[];  /**< @brief Possible values for eos-name. */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* RUNPHI0FAMILY_CMDLINE_H */
