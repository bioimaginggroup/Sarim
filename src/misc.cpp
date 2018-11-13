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
#include <dqrng_distribution.h>
#include "misc.hpp"

class SplitMix : public dqrng::random_64bit_generator {
public:
  SplitMix (result_type seed) : state(seed) {};
  result_type operator() () {
    result_type z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  void seed(result_type seed) {state = seed;}
  
private:
  result_type state;
};

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
  //Test with dqrng
  
  /*
    std::random_device device_random;
    std::default_random_engine generator(device_random());
    std::normal_distribution<> distribution(0, 1);
    Eigen::VectorXd e(n);
     
    for (int i = 0; i < n; ++i)
    {
        e(i) = distribution(generator);
    }
    return e;
*/
  Rcpp::NumericVector f(n);
  auto rng = dqrng::generator<SplitMix>(42);
    dqrng::normal_distribution dist(0.0, 1.0);
    f = dqrng::generate<dqrng::normal_distribution, Rcpp::NumericVector>(n, rng, dist);
    Eigen::VectorXd e(n);
    
    for (int i = 0; i < n; ++i)
    {
      e(i) = f(i);
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