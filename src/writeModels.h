/*
 * writeModels.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * Declarations for functions which write models to files
 *
 */

#pragma once

#include "BackgroundModel.h"
#include "ModelO1.h"
#include "ModelO2.h"

#include <string>
#include <ostream>

void writeBackgroundModel( std::ostream &o,
                           const BackgroundModel &model );
void writeBackgroundModel( const std::string & filename,
                           const BackgroundModel &model );

void writeModelO1( std::ostream &o,
                   const ModelO1 &model );
void writeModelO1( const std::string & filename,
                   const ModelO1 &model );

void writeModelO2( std::ostream &o,
                   ModelO2 &model );
void writeModelO2( const std::string & filename,
                   ModelO2 &model );

void writeBG12( std::ostream &o,
                ModelO2 &model );
void writeBG12( const std::string & filename,
                ModelO2 &model );

void writeSummaryHeader1( std::ostream &o );
void writeSummaryHeader2( std::ostream &o );
void writeSummaryLine(    std::ostream &o,
                          const BackgroundModel &model0,
                          const ModelO1 &model1);
void writeSummaryLine(    std::ostream &o,
                          const BackgroundModel &model0,
                          const ModelO1 &model1,
                          const ModelO2 &model2 );
