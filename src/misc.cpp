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

#ifndef RUE_H_
#define RUE_H_

#include <RcppEigen.h>
#include <random>
// [[Rcpp::depends(dqrng)]]
#include <dqrng.h>
#include "misc.hpp"


// function to calculate the incomplete cholesky decomposition, return a sparse matrix
// need for lanczos algorithm
Eigen::SparseMatrix<double> ichol (const Eigen::SparseMatrix<double> & Q) {
    Eigen::IncompleteCholesky<double> icholOfQ(Q);
    return icholOfQ.matrixL();
};


// faster C++-sampling from a gamma distribution
Eigen::VectorXd random_gamma (const int & n, const double & shape, const double & scale) {
    

     Eigen::VectorXd e(n);
    
    std::random_device device_random;
    std::default_random_engine generator(device_random());
    std::gamma_distribution<> distribution(shape, scale);
    
    for (int i = 0; i < n; ++i)
    {
        e(i) = distribution(generator);
    }
      return e;
};



// faster C++-sampling from a gaussian distribution
Eigen::VectorXd random_gauss (const int & n)
{

  Eigen::VectorXd e(n);
  //dqrng::dqset_seed(time(0));
  dqrng::dqRNGkind("Xoroshiro128+");
  Rcpp::NumericVector e0 = dqrng::dqrnorm(n,0.0,1.0);
  for (int i = 0; i < n; ++i)
  {
    e(i) = e0(i);
  }
  return e;

};


// faster C++-sampling from a uniform distribution
Eigen::VectorXd random_uniform (const int & n) {
    Eigen::VectorXd e(n);
    
    std::random_device device_random;
    std::default_random_engine generator(device_random());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    
    for (int i = 0; i < n; ++i)
    {
        e(i) = distribution(generator);
    }
    return e;
};

#endif // RUE_H_