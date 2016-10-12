// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::interfaces(r,cpp)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List ama_eigen(   Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w, Eigen::Map<Eigen::MatrixXd> ix,
                   double gamma1, double gamma2, Eigen::Map<Eigen::VectorXd> gamma2_weight,
                   double nu, double tol_abs, double max_iter, double type){

    using namespace Eigen;

    int n = X.cols();
    int p = X.rows();
    int nk = ix.rows();
    MatrixXd Lambda(p, nk); Lambda.fill(0);
    MatrixXd V(p, nk);

    // defined for update A
    VectorXd e(n);
    VectorXd term1(n);
    VectorXd x0(n); x0.fill(0);
    VectorXd Z(n);
    MatrixXd A(p, n); A = X;
    MatrixXd A0(p, n); A0.fill(0);

    // defined for update V
    double norm_lambda;
    double norm_term3;
    VectorXd term3(p);
    MatrixXd Lambda0(p, nk); Lambda0.fill(0);

    // defined for while loop
    VectorXd eva(2); eva.fill(1); // evaluate convergence
    int l1;
    int l2;
    int iter = 0;

    MatrixXd X_T(n, p); X_T = X.transpose();

    while( (iter < max_iter) & (eva.maxCoeff() > tol_abs)){
        iter++;
      //update A
      for(int i = 0; i < p; i++){
        term1.fill(0);
        for(int l=0; l < nk; l++){
          l1 = ix(l,0);
          l2 = ix(l,1);

          term1(l1) += Lambda(i,l);
          term1(l2) -= Lambda(i,l);
          }

          Z = X_T.col(i) + term1;
          A0.row(i) = std::max(0.0, 1 - gamma2 * gamma2_weight(i) / Z.norm() ) * Z;

      }


      eva(0) = (A - A0).norm();
      //eva(0) = 1;
      A = A0;



      //update Lambda
      for(int l = 0; l < nk; l++){
        l1 = ix(l,0);
        l2 = ix(l,1);
        term3 = Lambda.col(l);
        term3 -= nu * (A.col(l1) - A.col(l2));
        norm_lambda = gamma1 * w(l);
        norm_term3 = term3.norm();


        if(type == 2){
                // L2 dual norm projection on to ball
          if(norm_term3 < norm_lambda){
            Lambda0.col(l) = term3;
          }else{
            Lambda0.col(l) = term3/norm_term3 * norm_lambda;
          }
        }
      }

      eva(1) = (Lambda0 - Lambda).norm();
      //eva(1) = 1;
      Lambda = Lambda0;

      }

    // update V after convergence
    for(int l = 0; l < nk; l++){
      l1 = ix(l,0);
      l2 = ix(l,1);
      V.col(l)  = A.col(l1) - A.col(l2);
    }


    return Rcpp::List::create(Rcpp::Named("A") = A,
                              Rcpp::Named("V") = V,
                              Rcpp::Named("Lambda") = Lambda,
                              Rcpp::Named("iter") = iter
    );
}


