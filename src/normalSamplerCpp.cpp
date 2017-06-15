#include "mleHeader.h"
using namespace Rcpp;


double sampleBoundedTruncNorm(double mu, double sd, double lower, double upper) {
  double u = runif(1)[0] ;
  double phiB = R::pnorm5(upper, mu, sd, 1, 0) ;
  double phiA = R::pnorm5(lower, mu, sd, 1, 0) ;
  double quantile = u * phiB + phiA * (1 - u) ;
  double sample = R::qnorm(quantile, mu, sd, 1, 0) ;
  return sample ;
}

// [[Rcpp::export]]
NumericVector mvtSampler(NumericVector y,
                         NumericVector mu,
                         IntegerVector selected,
                         NumericMatrix threshold,
                         NumericMatrix precision,
                         int nsamp, int burnin, int trim,
                         bool verbose) {
  int totalIter = burnin + (nsamp - 1) * trim ;
  int frac = std::round(totalIter / 5) ;
  int p = y.length() ;
  NumericMatrix samples(nsamp, p) ;
  NumericVector samp = clone(y) ;
  double condmean, condsd, pprob, nprob, u ;
  int row = 0;

  for(int i = 0 ; i < (burnin + trim * nsamp) ; i ++) {
    for(int j = 0 ; j < p ; j ++) {
      condmean = computeConditionalMean(mu, samp, precision, 1, j) ;
      condsd = 1 / std::sqrt(precision(j, j)) ;
      if(selected[j] == 1) {
        nprob = R::pnorm(threshold(j, 0), condmean, condsd, 1, 1) ;
        pprob = R::pnorm(threshold(j, 1), condmean, condsd, 0, 1) ;
        nprob = 1 / (1 + std::exp(pprob - nprob)) ;
        u = runif(1)[0] ;
        if(u < nprob) {
          samp[j] = sampleUnivTruncNorm(condmean, condsd, threshold(j, 0)) ;
        } else {
          samp[j] = sampleUnivTruncNorm(-condmean, condsd, -threshold(j, 1)) ;
          samp[j] = -samp[j] ;
        }
      } else {
        samp[j] = sampleBoundedTruncNorm(condmean, condsd, threshold(j, 0), threshold(j, 1)) ;
      }
    }

    if(verbose & ((((i + 1) % frac) == 0) | (i == totalIter - 1))) {
      int out = (i + 1.0) / (1.0 * totalIter) * 100 ;
      Rcpp::Rcout<<out<<"% " ;
    }

    if(i >= (burnin - 1) & (i - burnin + 1) % trim == 0) {
      for(int j = 0 ; j < samples.ncol() ; j++) {
        samples(row, j) = samp[j] ;
      }
      if(++row == nsamp) {
        break ;
      }
    }
  }

  if(verbose) Rcpp::Rcout<<"\n" ;

  return samples ;
}


