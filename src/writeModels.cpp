/*
 * writeModels.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * Definitions for functions which write models to files
 *
 */

#include "writeModels.hpp"
#include "BackgroundModel.hpp"
#include "constants.hpp"

#include <string>
#include <fstream>

void writeBackgroundModel( std::ostream &o, const BackgroundModel &model )
{

  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#   01  index of radial point in model\n"
    "#   02  radius         [cm^1]\n"
    "#   03  mu(r)          [cm^1]\n"
    "#   04  nu(r)          [cm^0]\n"
    "#   05  phi(r)         [cm^0]\n"
    "#   06  A(phi(r))      [cm^0]\n"
    "#   07  pressure       [cm^-2]\n"
    "#   08  energy density [cm^-2]\n"
    "#   09  rest mass dens [cm^-2]\n"
    "#   10  mbar(r)        [cm^1]\n";

  const char s = ' ';

  o << header;

  const long iMax = model.iMax();
  for ( int i=0; i <= iMax; i++)
  {
    const double r       = model.r(i);
    const double mu      = model.mu(i);
    const double nu      = model.nu(i);
    const double phi     = model.phi(i);
    const double A       = model.conf.A(phi);
    const double geomP   = model.p(i);
    const double geomeps = model.eos.geomepsilonOfgeomP( geomP );
    const double geomrho = model.eos.geomrhoOfgeomP( geomP );
    const double mbar    = model.mbar(i);
    o << i << s
      << r << s
      << mu << s
      << nu << s
      << phi << s
      << A << s
      << geomP << s
      << geomeps << s
      << geomrho << s
      << mbar << s
      << std::endl;
  };
};

void writeBackgroundModel( const std::string & filename, const BackgroundModel &model )
{
  std::ofstream of(filename);
  writeBackgroundModel( of, model );
  of.close();
};

void writeSummaryHeader( std::ostream &o )
{
  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#    1  central pressure  [cm^-2]\n"
    "#    2  central phi       [cm^0]\n"
    "#    3  M_ADM             [Msun]\n"
    "#    4  M_b               [Msun]\n"
    "#    5  R_E               [km^1]\n"
    "#    6  R_J               [km^1]\n"
    "#    7  phi_0             [cm^0]\n"
    "#    8  alpha             [cm^0]\n"
    "#    9  omega             [cm^1]\n";

  o << header;

};

void writeSummaryLine( std::ostream &o,
                       const BackgroundModel &model0)
{

  const char s = ' ';

  o << model0.p(0) << s
    << model0.phi(0) << s
    << model0.M_ADM()/GMsun_cm << s
    << model0.M_b()/GMsun_cm << s
    << model0.R_areal_Einstein()*1.e-5 << s
    << model0.R_areal_Jordan()*1.e-5 << s
    << model0.phi_0() << s
    << model0.alpha() << s
    << model0.omega() << s
    << std::endl;

};
