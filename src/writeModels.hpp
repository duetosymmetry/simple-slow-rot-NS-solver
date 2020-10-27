/*
 * writeModels.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * Declarations for functions which write models to files
 *
 */

#pragma once

#include "BackgroundModel.hpp"

#include <string>
#include <ostream>

void writeBackgroundModel( std::ostream &o,
                           const BackgroundModel &model );
void writeBackgroundModel( const std::string & filename,
                           const BackgroundModel &model );

void writeSummaryHeader( std::ostream &o );

void writeSummaryLine(    std::ostream &o,
                          const BackgroundModel &model0);
