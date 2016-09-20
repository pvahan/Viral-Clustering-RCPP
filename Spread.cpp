#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec Spread(arma::vec y, arma::vec samp, arma::vec sampMat, arma::vec zz, int K, int KTrue, int n) {
	int j, jj, temp;
  for(int i = 0; i < n; i++) {
    j = samp(i) - 1;
    jj = sampMat(j)-1;
    zz(y(j)-1) = zz(y(j)-1) - 1;
    zz(y(jj)-1) = zz(y(jj)-1) + 1;

    temp = zz(y(j)-1);
    y(j) = y(jj);
    if(temp == 0) {
      K = K - 1;
      if(K <= KTrue) break;
    }

  }
  return y;
}


