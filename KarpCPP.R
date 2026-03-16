library(Rcpp)

Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List karp_cpp(NumericMatrix cij){

  int n = cij.nrow();

  NumericMatrix dkv(n, n+1);

  for(int i=0;i<n;i++){
    for(int j=0;j<n+1;j++){
      dkv(i,j) = R_PosInf;
    }
  }

  dkv(0,0) = 0.0;

  for(int k=1;k<=n;k++){

    for(int j=0;j<n;j++){

      double best = R_PosInf;

      for(int i=0;i<n;i++){

        double val = dkv(i,k-1) + cij(i,j);

        if(val < best)
          best = val;
      }

      dkv(j,k) = best;
    }
  }

  NumericMatrix dndk(n,n);

  for(int i=0;i<n;i++){
    for(int k=0;k<n;k++){

      dndk(i,k) = (dkv(i,n) - dkv(i,k)) / (n-k);

    }
  }

  NumericVector dmax(n);

  for(int i=0;i<n;i++){

    double best = dndk(i,0);

    for(int k=1;k<n;k++){

      if(dndk(i,k) > best)
        best = dndk(i,k);

    }

    dmax[i] = best;
  }

  double mu_star = dmax[0];

  for(int i=1;i<n;i++){
    if(dmax[i] < mu_star)
      mu_star = dmax[i];
  }

  NumericVector shortest(n);

  for(int i=0;i<n;i++){

    double best = R_PosInf;

    for(int k=0;k<=n;k++){

      double val = dkv(i,k) - mu_star*k;

      if(val < best)
        best = val;
    }

    shortest[i] = best;
  }

  return List::create(
    Named("mu_star") = mu_star,
    Named("shortest") = shortest
  );
}
')