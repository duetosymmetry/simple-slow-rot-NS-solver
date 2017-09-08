/*
 * ppEOSTable.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2013 Dec. 17
 *
 * This is the interface for getting a parametrized EOS from
 * Table III of Read, Lackey, Owen, and Friedman (2009) [arXiv:0812.2163]
 *
 */

#pragma once

#include <string>
#include "ppEOS.h"

ppEOS findEOS(const std::string &name);
