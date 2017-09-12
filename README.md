[//]: # (Distributed under the MIT License.)
[//]: # (See LICENSE for details.)

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](LICENSE)

# Simple slow-rotation neutron star structure solver
This started out as Kent's code for
[1303.1528](https://arxiv.org/abs/1303.1528), and was initially Nico's
code before that.  I rewrote it a lot.  You will thus notice wildly
varying programming styles throughout.  It can still be greatly
improved, but as we say, it works for me for now.

This code is capable of solving for a neutron star structure at orders
0, 1, and 2 in the slow rotation expansion.  Order 0 is of course the
"Tolman-Oppenheimer-Volkoff" (TOV) equation.  At orders 1 and 2 we can
extract the moment of inertia and quadrupole moment.

The piecewise polytropic EOS implements the prescription in Read,
Lackey, Owen, and Friedman
(2009) [0812.2163](https://arxiv.org/abs/0812.2163), along with their
fits for named EOSs.

## Dependencies
* make
* makedepend
* GSL
* gengetopt
* A C++ compiler. I am using clang and have not tested others.
* I'm probably forgetting something, add it here if you pick up on a
  missing dep.
