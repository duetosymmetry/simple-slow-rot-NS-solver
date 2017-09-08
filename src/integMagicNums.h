/*
 * integMagicNums.h
 * Author: Leo C. Stein (leo.stein@gmail.com)
 * Date: 2014 Jan
 *
 * Magic numbers for integration
 *
 */

#pragma once

// Length of arrays
#define BG_MAX_SIZE 1000000
// The R epsilon (in cm) to start the integration away from r=0
#define EPSR 0.5
// The R step size (in cm)
#define DR 500.
// Absolute and relative error tolerances
#define EPSABS 3.e-5
#define EPSREL 0.0
// the pressure contrast to end integration
#define P_RATIO 1.e-10
