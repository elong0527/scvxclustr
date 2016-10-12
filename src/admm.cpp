// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List admm_eigen(   Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w,
             double gamma1, double gamma2, Eigen::Map<Eigen::VectorXd> gamma2_weight,
             double nu, double tol_abs, double max_iter, double type, Function f){

            using namespace Eigen;

            int nrow = X.rows();
            int ncol = X.cols();
            int nk = nrow * ( nrow - 1 )/2;
            int n = nrow;
            MatrixXd Lambda(ncol, nk); Lambda.fill(0);
            MatrixXd V(ncol, nk); V.fill(0);

            MatrixXd I = MatrixXd::Identity(nrow, nrow);
            MatrixXd J(nrow, nrow); J.fill(1);

            MatrixXd N(nrow, nrow);
            MatrixXd N_inv(nrow, nrow);

            // defined for update A
            VectorXd e(nrow);
            VectorXd term1(nrow);
            VectorXd x0(nrow); x0.fill(0);
            MatrixXd Y(nrow, ncol); Y.fill(0);
            MatrixXd Z(nrow, ncol); Z.fill(0);

            MatrixXd A(nrow, ncol); A = X;
            MatrixXd A0(nrow, ncol); A0.fill(0);
            int l;


            // defined for update V
            MatrixXd Lambda0(ncol, nk); Lambda0.fill(0);
            MatrixXd V0(ncol, nk); V0.fill(0);
            double term2 = sqrt(1 + nrow*nu);
            VectorXd term3(ncol);
            double sigma;
            VectorXd pair(2);
            pair(0) = 0;

            // evaluate convergence
            VectorXd eva(3); eva.fill(1);

            int iter = 0;
            while( (iter < max_iter) & (eva.maxCoeff() > tol_abs)){
            iter++;
            A0 = A;
              // Update A

              V0 = V + Lambda / nu;
              N  = term2 * I - (term2 - 1)/n * J;
              N_inv = ( I + J * (term2 - 1) / n ) / term2;
              for(int i = 0; i < ncol; i++){

                term1.fill(0);
                l = 0;
                for(int l1 = 0; l1 < nrow; l1++ ){
                  for(int l2 = 0; l2 < nrow; l2++ ){
                    if(l1 < l2){
                      e.fill(0);
                      e(l1) = 1;
                      e(l2) = -1;
                      term1 += V0(i,l) * e;
                      l++;
                    }
                  }
                }
                Z.col(i) = X.col(i) + nu * term1;
              }

              Y = N_inv * Z;


              for(int ll =0; ll < ncol; ll++){
                NumericVector tmp = f(Y.col(ll), N, gamma2, gamma2_weight(ll));
                A0.col(ll) = as<Map<VectorXd> >(tmp);
              }

              eva(0) = (A - A0).norm();
              A = A0;

              // Update V

                l = 0;
                for(int l1 = 0; l1 < nrow; l1++ ){
                  for(int l2 = 0; l2 < nrow; l2++ ){
                    if(l1 < l2){

                      term3 = A.row(l1) - A.row(l2) ;
                      term3 -= Lambda.col(l) / nu;
                      sigma = gamma1 * w(l) / nu;
                      pair(1) = 1 - sigma/ term3.norm()  ;
                      V0.col(l) =  pair.maxCoeff() * term3;
                      Lambda0.col(l) =  - A.row(l1) + A.row(l2) ;
                      Lambda0.col(l) = Lambda.col(l) + nu * (V0.col(l) + Lambda0.col(l));
                      l++;

                    }
                  }
                }
                eva(1) = (V - V0).norm();
                eva(2) = (Lambda - Lambda0).norm();

                V = V0;
                Lambda = Lambda0;
            }


            // return wrap( A);

            return Rcpp::List::create(Rcpp::Named("A") = A,
                                      Rcpp::Named("V") = V,
                                      Rcpp::Named("Lambda") = Lambda,
                                      Rcpp::Named("eva") = eva,
                                      Rcpp::Named("iter") = iter
                                      );
}




