/*
 * writeModels.cc
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * Definitions for functions which write models to files
 *
 */

#include "writeModels.h"
#include "BackgroundModel.h"
#include "ModelO1.h"
#include "ModelO2.h"

#include <string>
#include <fstream>

void writeBackgroundModel( std::ostream &o, const BackgroundModel &model )
{

  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#    1  index of radial point in model\n"
    "#    2  radius         [cm]\n"
    "#    3  m(r)           [cm]\n"
    "#    4  nu(r)          [1]\n"
    "#    5  pressure       [cm^-2]\n"
    "#    6  energy density [cm^-2]\n";

  const char s = ' ';

  o << header;

  const long iMax = model.iMax();
  for ( int i=0; i <= iMax; i++)
  {
    const double r       = model.r(i);
    const double m       = model.m(i);
    const double nu      = model.nu(i);
    const double geomP   = model.p(i);
    const double geomeps = model.eos.geomepsilonOfgeomP( geomP );
    o << i << s
      << r << s
      << m << s
      << nu << s
      << geomP << s
      << geomeps << s
      << std::endl;
  };
};

void writeBackgroundModel( const std::string & filename, const BackgroundModel &model )
{
  std::ofstream of(filename);
  writeBackgroundModel( of, model );
  of.close();
};

void writeModelO1( std::ostream &o, const ModelO1 &model )
{

  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#    1  index of radial point in model\n"
    "#    2  radius         [cm]\n"
    "#    3  omega1         [cm^-1]\n";

  const char s = ' ';

  o << header;

  const long iMax = model.iMax();
  for ( int i=0; i <= iMax; i++)
  {
    double r       = model.bg.r(i);
    double omega1  = model.omega1(i);
    o << i << s << r << s << omega1
      << std::endl;
  };
};

void writeModelO1( const std::string & filename, const ModelO1 &model )
{
  std::ofstream of(filename);
  writeModelO1( of, model );
  of.close();
};

void writeModelO2( std::ostream &o, ModelO2 &model )
{

  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#    1  index of radial point in model\n"
    "#    2  radius [cm]\n"
    "#    3  K2     [cm^??]\n"
    "#    4  h2     [cm^??]\n"
    "#    5  xi2    [cm^??]\n";

  const char s = ' ';

  o << header;

  const long iMax = model.iMax();
  for ( int i=0; i <= iMax; i++)
  {
    const double r   = model.bg.r(i);
    const double K2  = model.K2(i);
    const double h2  = model.h2(i);
    const double xi2 = model.xi2Ofr(r);
    o << i << s
      << r << s
      << K2 << s
      << h2 << s
      << xi2
      << std::endl;
  };
};

void writeModelO2( const std::string & filename, ModelO2 &model )
{
  std::ofstream of(filename);
  writeModelO2( of, model );
  of.close();
};

void writeBG12( std::ostream &o, ModelO2 &model2 )
{

  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#    1  index of radial point in model\n"
    "#    2  radius         [cm]\n"
    "#    3  m(r)           [cm]\n"
    "#    4  nu(r)          [1]\n"
    "#    5  pressure       [cm^-2]\n"
    "#    6  energy density [cm^-2]\n"
    "#    7  omega1         [cm^-1]\n"
    "#    8  K2             [cm^??]\n"
    "#    9  h2             [cm^??]\n"
    "#    10 xi2            [cm^??]\n";

  const char s = ' ';

  o << "# " << model2.summary() << std::endl;
  o << header;

  const BackgroundModel &bg = model2.bg;
  const ModelO1 &model1     = model2.model1;

  const long iMax = bg.iMax();
  for ( int i=0; i <= iMax; i++)
  {
    const double r       = bg.r(i);
    const double m       = bg.m(i);
    const double nu      = bg.nu(i);
    const double geomP   = bg.p(i);
    const double geomeps = bg.eos.geomepsilonOfgeomP( geomP );
    const double omega1  = model1.omega1(i);
    const double K2      = model2.K2(i);
    const double h2      = model2.h2(i);
    const double xi2     = model2.xi2Ofr(r);
    o << i << s
      << r << s
      << m << s
      << nu << s
      << geomP << s
      << geomeps << s
      << omega1 << s
      << K2 << s
      << h2 << s
      << xi2
      << std::endl;
  };
};

void writeBG12( const std::string & filename, ModelO2 &model )
{
  std::ofstream of(filename);
  writeBG12( of, model );
  of.close();
};

void writeSummaryHeader1( std::ostream &o )
{
  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#    1  central pressure  [cm^-2]\n"
    "#    2  total mass        [cm]\n"
    "#    3  stellar radius    [cm]\n"
    "#    4  moment of inertia [cm^??]\n";

  o << header;

};

void writeSummaryHeader2( std::ostream &o )
{
  // Put a header at the beginning of the file with the format
  const std::string header =
    "# Format for this file:\n"
    "# col#  description\n"
    "#    1  central pressure  [cm^-2]\n"
    "#    2  total mass        [cm]\n"
    "#    3  stellar radius    [cm]\n"
    "#    4  moment of inertia [cm^??]\n"
    "#    5  quadrupole moment [cm^??]\n";

  o << header;

};

void writeSummaryLine( std::ostream &o,
                       const BackgroundModel &model0,
                       const ModelO1 &model1)
{

  const char s = ' ';

  o << model0.p(0) << s
    << model0.M() << s
    << model0.R() << s
    << model1.I() << std::endl;

};

void writeSummaryLine( std::ostream &o,
                       const BackgroundModel &model0,
                       const ModelO1 &model1,
                       const ModelO2 &model2)
{

  const char s = ' ';

  o << model0.p(0) << s
    << model0.M() << s
    << model0.R() << s
    << model1.I() << s
    << model2.Q() << std::endl;

};
