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
//   Definition for iterative sampling from a gaussian distribution by Lanczos
//
// Written by:
//   Christopher Küster
//
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "lanczos.hpp"
#include "misc.hpp"


// generate the requried tridiagonalmatrix 
// Input: two vectors alpha and beta, with alpa on main diagonal 
//        and beta on upper and lower one off diagonal
inline Eigen::MatrixXd tri_mat_creation (const int & n, 
                                         const Eigen::VectorXd & a, 
                                         const Eigen::VectorXd & b) 
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd BU = A;
    Eigen::MatrixXd BL = A;
    
    A.diagonal() = a.segment(0, n);
    BU.diagonal<1>() = b.segment(1, n);
    BL.diagonal<-1>() = b.segment(1, n);
    
    return A + BU + BL;
};



// lanczos algorithm
// Input: Q     precision matrix, is precision matrix from ~N(mu, Q)
//        m     number of iterations, can be set by user
//        F1,F2 incomplete cholesky factorization
//        thr   when the algorithm stop
Lanczos algorithm (const Eigen::SparseMatrix<double> & Q, 
                   const int & m, 
                   const Eigen::SparseMatrix<double> & F1, 
                   const Eigen::SparseMatrix<double> & F2, 
                   const double & thr) 
{
    Lanczos lanczos_solver;
//    int coco=0;
    int iter; 
    int n = Q.cols();

    Eigen::VectorXd z = random_gauss(n);

    Eigen::MatrixXd V (n, m + 1);

    V.col(0) = z / z.norm();

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(m + 1);
    Eigen::VectorXd alpha = Eigen::VectorXd::Zero(m);

    Eigen::VectorXd w, x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd x_old = x;
    Eigen::MatrixXd T;

    Eigen::VectorXd err = alpha;

    Eigen::VectorXd a = x;
    Eigen::VectorXd b = a;


    Eigen::SparseLU<Eigen::SparseMatrix<double> > LUF1(F1); 
    Eigen::SparseLU<Eigen::SparseMatrix<double> > LUF2(F2); 

    for (int j = 0; j < m; j++) {
      iter = j + 1;
      x_old = x;
//      double temp;
      a = LUF2.solve(V.col(j));
        b = Q * a;
        w = LUF1.solve(b);

        if(j > 0) {
            w -= beta(j) * V.col(j - 1);
        }
        
        alpha.row(j) = w.transpose() * V.col(j);
        
        w -= alpha(j) * V.col(j);
        
        beta(j + 1) = w.norm();
        
        V.col(j + 1) = w / w.norm();

        if (j > 1) {
            T = tri_mat_creation (j + 1, alpha, beta);
            Eigen::VectorXd e1 = Eigen::VectorXd::Zero(j + 1);
            e1(0) = 1;
            x = z.norm() * V.leftCols(j) * T.sqrt().lu().solve(e1);
            err(j) = (x - x_old).norm() / x.norm();
            if (err(j) < thr) {
                break;
            }
        }
    }

    lanczos_solver.x = LUF2.solve(x);
    lanczos_solver.error = err;
    lanczos_solver.Iteration = iter;
    return lanczos_solver;
}