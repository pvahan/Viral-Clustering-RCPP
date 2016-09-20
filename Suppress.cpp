#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec Suppress(arma::mat data, arma::vec y) {
  arma::vec uniqy = unique(y);
  int n = data.n_rows, k = uniqy.size();
  arma::mat centers(uniqy.size(),data.n_cols);
  
  for(int i = 0; i < k; i++) {
    centers.row(i) = mean( data.rows( find(y == uniqy (i) ) ), 0);
  }
  
  for(int i = 0; i < n; i++) {
    int minIndex = 0;
    double minVecDist = norm(centers.row(0) - data.row(i), 2);
    double curDist = minVecDist;
    for(int j = 0; j < k; j++) {
      curDist = norm(centers.row(j) - data.row(i), 2);
      if(curDist < minVecDist) {
        minVecDist = curDist;
        minIndex = j;
      }
    }
    y(i) = uniqy(minIndex);
  }
  
  return y;
}


