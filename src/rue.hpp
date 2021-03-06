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
//   Declaration for Rue (2001) algorithm to sample form gaussian distribution
//
// Written by:
//   Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////

#ifndef RUE_H_
#define RUE_H_

#include <RcppEigen.h>


// structure for multiple returns
struct RueSolv 
{
    Eigen::VectorXd ga;
    Eigen::VectorXd mu;
};


// Rue-algorithm
//
// first decompose Q by Cholesky, i.e. Q = LL'
// a sample x ~ N(0, Q^{-1}) can then be obtained by sampling z from N(0, I)
// and solving L'x = z
// A sample from N(mu, Q^{-1}) is solving Q mu = b; first by solving Lw = b 
// and subsequently solving L'mu = v
// Apply a linear constraint with  Q^{-1}*A'*(A*Q^{-1}*A')^{-1} * (Ax - e)
// where A a row vector of ones and e = 0
RueSolv algorithm (const Eigen::SparseMatrix<double> & Q, 
                   const Eigen::VectorXd & b,
                   const std::string & lin_con);


#endif // RUE_H_