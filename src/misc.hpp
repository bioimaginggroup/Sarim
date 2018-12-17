////////////////////////////////////////////////////////////////////////////////
//                  _
//     ________ ___(_) ___
//    / __/ _  | __| |/   \
//    \__ \(_) | | | | Y Y |
//    /___/\_|_|_| |_|_|_|_|
//
//
// sarim is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// sarim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with sarim. If not, see <http://www.gnu.org/licenses/>.
//
// This file contains:
// -------------------
//
//   Definition for miscellaneous function need by sarim
//      Include:
//          - Incomplete Cholesky for Lanczos-Approximation
//          - Random number samples for gaussian, gamma and uniform distributions
//
// Written by:
//   Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MISC_H_
#define MISC_H_

#include <RcppEigen.h>
#include <random>

Eigen::SparseMatrix<double> ichol (const Eigen::SparseMatrix<double> & Q);

Eigen::VectorXd random_gamma (const int & n, const double & shape, const double & scale);

Eigen::VectorXd random_gauss (const int & n);

Eigen::VectorXd random_uniform (const int & n);

#endif // MISC_H_
