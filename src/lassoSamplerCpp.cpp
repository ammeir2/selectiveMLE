#include <Rcpp.h>
using namespace Rcpp;

# define M_PI           3.14159265358979323846  /* pi */
const double log2pi = log(2.0 * 3.1415926535897932384);

void printVec(NumericVector x) {
  for(int i = 0; i < x.length() ; i ++) {
    Rcpp::Rcout<<x[i]<<" ";
  }
  Rcpp::Rcout<<"\n" ;
}

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

void randomShuffle(IntegerVector a) {
  std::random_shuffle(a.begin(), a.end(), randWrapper);
}

double sampleExtreme(double mu, double sd, double threshold) {
  double sign = 1 ;
  double proposal ;
  double alpha ;
  double phi ;

  sign = -1 ;
  mu *= sign ;
  threshold = threshold * sign ;

  // rescaling
  threshold = (threshold - mu) / sd ;
  alpha = (threshold + sqrt(std::pow(threshold, 2) + 4)) / 2 ;

  bool reject = true ;
  int iter = 0;
  while(reject & (iter++ < 10000)) {
    proposal = threshold + R::rexp(1 / alpha) ;
    phi = exp(-std::pow(proposal - alpha, 2) / 2) ;
    if(runif(1)[0] < phi) {
      reject = false ;
    }
  }

  proposal = proposal * sd + mu ;
  return proposal * sign;
}

double sampleUnivTruncNorm(double mu, double sd, double threshold) {
  double u = runif(1)[0] ;
  double phiThreshold, sample ;

  // if(isnan(mu)) {
  //   Rcpp::Rcout<<"mu is nan \n" ;
  //   return 0 ;
  // }

  if((std::abs(mu - threshold) / sd) > 3) {
    return sampleExtreme(mu, sd, threshold) ;
  }

  phiThreshold = R::pnorm5(threshold, mu, sd, 1, 0) ;
  sample = R::qnorm5(u * phiThreshold, mu, sd, 1, 0) ;

  int tries = 0 ;
  while(isnan(sample) & tries ++ < 10) {
    sample =  sampleExtreme(mu, sd, threshold) ;
  }

  return sample ;
}

double computeConditionalMean(NumericVector mu,
                              NumericVector samp,
                              const NumericMatrix XmX,
                              double yvar,
                              int index) {
  double result = 0 ;

  for(int j = 0; j < mu.length() ; j ++) {
    if(j != index) {
      result += XmX(index, j) * (samp[j] - mu[j]) / yvar ;
    }
  }

  result = result / XmX(index, index) * yvar ;
  result = mu[index] + result ;
  return result ;
}

int innerZeroCondition(NumericVector samp,
                       const NumericVector l,
                       const NumericVector u) {
  int count = 1;
  double x ;

  for(int j = 0; j < samp.length() ; j++) {
    x = samp[j] ;
    if(x > u[j] | x < l[j]) {
      count-- ;
      break ;
    }
  }

  return count;
}

int innerCheckOneCondition(NumericVector samp,
                           NumericVector threshold) {
  double x ;
  for(int i = 0; i < samp.length() ; i++) {
    x = - std::abs(samp[i]) ;
    if(x > threshold[i]) {
      return 0 ;
    }
  }

  return 1;
}

void computeZeroThresholds(const NumericMatrix u0mat,
                           NumericVector signs,
                           NumericVector l,
                           NumericVector u) {
  for(int i = 0; i < u0mat.nrow() ; i ++) {
    u[i] = 1;
    l[i] = -1;
    for(int j = 0; j < signs.length() ; j ++) {
      u[i] -= u0mat(i, j) * signs[j] ;
      l[i] -= u0mat(i, j) * signs[j] ;
    }
  }
}

void computeOneThreshold(NumericVector signs,
                         double lambda,
                         const NumericMatrix XmXinv,
                         NumericVector u) {
  for(int i = 0 ; i < XmXinv.nrow() ; i++) {
    u[i] = 0 ;
    for(int j = 0 ; j < XmXinv.ncol() ; j++) {
      u[i] += XmXinv(i, j) * signs[j] ;
    }
    u[i] *= lambda ;
  }
}

double computeDiffThreshold(NumericVector signs,
                            double lambda,
                            const NumericMatrix XmXinv,
                            int coordinate) {
  double result = 0;
  for(int i = 0 ; i < XmXinv.ncol() ; i++) {
    if(i == coordinate) {
      result -= XmXinv(coordinate, i) * signs[i] ;
    } else {
      result += XmXinv(coordinate, i) * signs[i] ;
    }
  }

  return result * lambda ;
}

void copyVector(NumericVector to, NumericVector from) {
  for(int i = 0 ; i < to.length() ; i++) {
    to[i] = from[i] ;
  }
}

void aOneSampler(NumericVector &samp,
                 const NumericVector &signs,
                 const NumericVector &u,
                 const NumericVector &mean,
                 const NumericMatrix &XmX,
                 const NumericVector &condSigma,
                 double ysigsq,
                 int burnin, int maxiter) {
  double condMean ;
  double condSD ;
  int i, j ;

  int iterations = 0 ;
  while(iterations++ < maxiter) {
    for(int j = 0; j < samp.length() ; j++) {
      condMean = computeConditionalMean(mean, samp, XmX, ysigsq, j) ;
      condSD = std::sqrt(condSigma[j]) ;
      samp[j] = sampleUnivTruncNorm(-signs[j] * condMean, condSD, -signs[j] * u[j]) ;
      samp[j] *= -signs[j] ;
    }
  }
}

void zeroSampler(NumericVector &zsamp,
                 NumericVector &zeroSamp,
                 NumericVector zeroMean,
                 NumericMatrix sqrtMat,
                 NumericVector lzero,
                 NumericVector uzero,
                 int miniters, int maxiters) {
  double eps = 0.00000000001 ;

  double currentLower, currentUpper ;
  double threshold ;
  double adjustment ;
  double x, sigmaij ;
  int i, j ;

  for(int i = 0; i < zsamp.length() ; i++) zsamp[i] = 0 ;
  for(int i = 0; i < zeroMean.length() ; i++) zeroSamp[i] = zeroMean[i] ;

  for(int k = 0 ; k < maxiters ; k++) {
    for(int i = 0 ; i < zsamp.length() ; i ++) {
      // Finding constraints for next sample
      currentLower = -10000000000 ;
      currentUpper = 10000000000 ;
      for(int j = 0; j < zeroSamp.length() ; j ++) {
        sigmaij = sqrtMat(j, i) ;
        if(std::abs(sigmaij) < eps) continue ;
        adjustment = - zeroSamp[j] + sigmaij * zsamp[i] ;
        if(sigmaij > 0) {
          currentLower = std::max(currentLower, (lzero[j] + adjustment) / sigmaij) ;
          currentUpper = std::min(currentUpper, (uzero[j] + adjustment) / sigmaij);
        } else {
          currentLower = std::max(currentLower, (uzero[j] + adjustment) / sigmaij) ;
          currentUpper = std::min(currentUpper, (lzero[j] + adjustment) / sigmaij);
        }
      }

      // truncated sampling
      x = R::pnorm5(currentUpper, 0, 1, 1, 0) - R::pnorm5(currentLower, 0, 1, 1, 0);
      x *= runif(1)[0] ;
      x += R::pnorm5(currentLower, 0, 1, 1, 0) ;
      x = R::qnorm5(x, 0, 1, 1, 0) ;
      for(int j = 0; j < sqrtMat.nrow() ; j++) {
        sigmaij = sqrtMat(j, i) ;
        zeroSamp[j] += x * sigmaij - sigmaij * zsamp[i];
      }
      zsamp[i] = x ;
    }
    if(k > miniters) {
      if(innerZeroCondition(zeroSamp, lzero, uzero) == 1) {
        return ;
      }
    }
  }
}

double sign(double x) {
  double eps = 0.00000001 ;
  if(std::abs(x) < eps) {
    return 0.0 ;
  } else if(x < 0) {
    return -1.0 ;
  } else {
    return 1.0 ;
  }
}

void computeGradient(NumericVector gradient, NumericVector Xy,
                     NumericVector samp, NumericMatrix XmX,
                     double stepCoef, double stepRate,
                     double gradientBound,
                     int iter, int delay) {
  stepCoef = stepCoef / std::pow(std::max(1, iter + 1 - delay), stepRate) ;
  for(int i = 0 ; i < gradient.length() ; i++) {
    gradient[i] = Xy[i];
    for(int j = 0 ; j < XmX.ncol() ; j++) {
      gradient[i] -= XmX(i, j) * samp[j] ;
    }
    gradient[i] *= stepCoef ;
    if(std::abs(gradient[i]) > gradientBound) {
      gradient[i] = gradientBound * sign(gradient[i]) ;
    }
  }
}

void computeConditionalBeta(NumericVector &betaOut,
                            const NumericMatrix &estimateMat,
                            int iter) {
  int meanStart = std::max(0, iter - 300) ;
  double denominator = iter - meanStart + 1 ;
  for(int i = 0; i < estimateMat.ncol() ; i++) {
    betaOut[i] = 0;
    int denominator = 0;
    for(int j = meanStart ; j <= iter ; j++) {
      denominator++ ;
      betaOut[i] += estimateMat(j, i) ;
    }
    betaOut[i] /= denominator ;
  }
}

void boundBeta(NumericVector &estimate, NumericVector naive) {
  for(int i = 0 ; i < estimate.length() ; i++) {
    if(naive[i] < 0) {
      if(estimate[i] < naive[i]) {
        estimate[i] = naive[i] ;
      } else if(estimate[i] > 0) {
        estimate[i] = 0;
      }
    } else {
      if(estimate[i] > naive[i]) {
        estimate[i] = naive[i] ;
      } else if(estimate[i] < 0) {
        estimate[i] = 0 ;
      }
    }
  }
}

double mvtLogDens(NumericVector samp, NumericVector mean,
                  NumericMatrix XmX, double ysigsq) {
  NumericVector diff = NumericVector(samp.length()) ;
  for(int i = 0 ; i < diff.length() ; i++) {
    diff[i] = samp[i] - mean[i] ;
  }

  NumericVector inner = NumericVector(diff.length()) ;
  for(int i = 0 ; i < inner.length() ; i ++) {
    inner[i] = 0 ;
    for(int j = 0 ; j < inner.length() ; j ++) {
      inner[i] += XmX(i, j) * diff[j] ;
    }
  }

  double result = 0 ;
  for(int i = 0 ; i < inner.length() ; i ++) {
    result += inner[i] * diff[i] ;
  }

  return - 0.5 / ysigsq * result ;
}

void mhSampler(NumericVector samp, NumericVector oldSamp,
               NumericVector signs, NumericVector newsigns,
               NumericVector mean, NumericMatrix &sigma,
               NumericVector condSigma,
               double lambda,
               NumericVector u, NumericVector newuone,
               NumericMatrix &XmXinv,
               NumericMatrix &XmX, double ysigsq,
               NumericVector lZero, NumericVector uZero,
               NumericVector newlzero, NumericVector newuzero,
               NumericMatrix &u0mat, NumericVector zeroSamp,
               int maxiter, IntegerVector order, NumericVector initEst,
               bool methodExact) {
  int j ;
  double condSD, condMean ;
  double thres, psame, pdiff ;
  int zeroCondition, oneCondition ;
  double newsamp ;
  double forwardDens, reverseDens, sd;
  double oldMVTdens, newMVTdens, mhRatio ;

  for(int iter = 0 ; iter < maxiter ; iter++) {
    randomShuffle(order) ;
    for(int i = 0 ; i < samp.length() ; i ++) {
      // Deciding whether to change signs
      j = order[i] ;
      if(std::abs(mean[j]) > std::abs(initEst[j]) &
         signs[j] == sign(initEst[j])) continue ;
      newsigns[j] *= - 1 ;
      condMean = computeConditionalMean(mean, samp, XmX, ysigsq, j) ;
      condSD = std::sqrt(condSigma[j]) ;
      thres  = computeDiffThreshold(newsigns, lambda, XmXinv, j) ;
      psame = R::pnorm(-signs[j] * u[j], -signs[j] * condMean, condSD, 1, 1) ;
      pdiff = R::pnorm(signs[j] * thres, signs[j] * condMean, condSD, 1, 1) ;
      pdiff = 1.0 / (1.0 + std::exp(psame - pdiff)) ;
      if(runif(1)[0] > pdiff) {
        newsamp = sampleUnivTruncNorm(-signs[j] * condMean, condSD, - signs[j] * u[j]) ;
        samp[j] = -signs[j] * newsamp ;
        oldSamp[j] = samp[j] ;
        newsigns[j] *= -1 ;
        continue ;
      }

      // check zero condition
      computeZeroThresholds(u0mat, newsigns, newlzero, newuzero) ;
      if(methodExact) {
        zeroCondition = innerZeroCondition(zeroSamp, newlzero, newuzero) ;
      }

      if(methodExact & (zeroCondition == 0)) {
        newsamp = sampleUnivTruncNorm(-signs[j] * condMean, condSD, - signs[j] * u[j]) ;
        samp[j] = -signs[j] * newsamp ;
        oldSamp[j] = samp[j] ;
        newsigns[j] *= -1 ;
        continue ;
      }

      // Check if other coordinates need to be sampled
      computeOneThreshold(newsigns, lambda, XmXinv, newuone) ;
      oneCondition = innerCheckOneCondition(samp, newuone) ;
      if(oneCondition == 1) { // if no, just sample the one coordinate
        signs[j] *= - 1;
        copyVector(u, newuone) ;
        newsamp = sampleUnivTruncNorm(-signs[j] * condMean, condSD, -signs[j] * u[j]) ;
        newsamp *= -signs[j] ;
        samp[j] = newsamp ;
        oldSamp[j] = newsamp ;
        if(methodExact) {
          copyVector(lZero, newlzero) ;
          copyVector(uZero, newuzero) ;
        }
        continue ;
      }

      // If yes, then do MH step
      forwardDens = 0;
      reverseDens = 0;
      for(int k = 0; k < samp.length() ; k++) {
        if(k == j) {
          newsamp = sampleUnivTruncNorm(signs[j] * condMean, condSD, signs[j] * thres) ;
          newsamp *= signs[j] ;
          samp[j] = newsamp ;
          forwardDens += R::dnorm(newsamp, condMean, condSD, 1) ;
        } else {
          sd = std::sqrt(sigma(k, k)) ;
          newsamp = sampleUnivTruncNorm(-signs[k] * samp[k], sd, -signs[k] * newuone[j]) ;
          newsamp *= -signs[k] ;
          samp[k] = newsamp ;
          forwardDens += R::dnorm(newsamp, oldSamp[k], sd, 1) ;
          forwardDens -= R::pnorm(-signs[k] * newuone[k], -signs[k] * oldSamp[k], sd, 1, 1) ;
          reverseDens += R::dnorm(oldSamp[k], samp[k], sd, 1) ;
          reverseDens -= R::pnorm(-signs[k] * u[k], -signs[k] * samp[k], sd, 1, 1) ;
        }
      }
      condMean = computeConditionalMean(mean, samp, XmX, ysigsq, j) ;
      reverseDens += R::dnorm(oldSamp[j], condMean, condSD, 1) ;
      newMVTdens = mvtLogDens(samp, mean, XmX, ysigsq) ;
      oldMVTdens = mvtLogDens(oldSamp, mean, XmX, ysigsq) ;

     // Do we boldly go where no sampler has gone before?
      mhRatio = newMVTdens - oldMVTdens + reverseDens - forwardDens ;
      mhRatio = std::exp(mhRatio) ;
      //Rcpp::Rcout<<" MH "<<mhRatio ;
      if(runif(1)[0] < mhRatio) {
        signs[j] *= -1 ;
        copyVector(oldSamp, samp) ;
        copyVector(u, newuone) ;
        if(methodExact) {
          copyVector(lZero, newlzero) ;
          copyVector(uZero, newuzero) ;
        }
      } else {
        newsigns[j] *= -1 ;
        copyVector(samp, oldSamp) ;
      }
    }
  }
}


// [[Rcpp::export]]
void   lassoSampler(const NumericVector initEst,
                    const NumericVector initSamp,
                    NumericMatrix oneCov, NumericMatrix XmX,
                    NumericMatrix XmXinv,
                    NumericVector condSigma,
                    double lambda, double ysigsq,
                    NumericVector zeroMean, NumericMatrix sqrtZero,
                    NumericMatrix u0mat,
                    int n, int p,
                    int nsamp, int burnin,
                    NumericVector Xy,
                    NumericMatrix &estimateMat, NumericMatrix &sampMat,
                    int delay, double stepRate, double stepCoef,
                    double gradientBound, int assumeConvergence,
                    NumericVector naive, bool methodExact, bool verbose) {
  // Initializing Sampling order
  IntegerVector order = IntegerVector(initEst.length()) ;
  for(int i = 0; i < order.length() ; i++) order[i] = i ;

  // initalizing signs
  NumericVector samp = clone(initSamp) ;
  NumericVector oldsamp = clone(initSamp) ;
  NumericVector estimate = clone(initEst) ;
  NumericVector signs = NumericVector(samp.length()) ;
  for(int i = 0; i < signs.length() ; i++) {
    signs[i] = -1 ;
    if(samp[i] > 0) signs[i] = 1 ;
  }
  NumericVector newsigns = clone(signs) ;

  // Initializing zero sample
  int k = samp.length() ;
  NumericVector zsamp, zeroSamp, lzero, uzero, newlzero, newuzero ;
  if(methodExact) {
    int zeroRank = std::min(p - k, n - k) ;
    zsamp = NumericVector(zeroRank) ;
    zeroSamp = NumericVector(p - k) ;
    lzero = NumericVector(p - k) ;
    uzero = NumericVector(p - k) ;
    newlzero = NumericVector(p - k) ;
    newuzero = NumericVector(p - k) ;
    computeZeroThresholds(u0mat, signs, lzero, uzero) ;
  }

  // initializing one sampler
  double currentThreshold ;
  NumericVector uone = NumericVector(samp.length()) ;
  NumericVector newuone = NumericVector(samp.length()) ;
  computeOneThreshold(signs, lambda, XmXinv, uone) ;

  // initializing optimizer
  NumericVector gradient = NumericVector(samp.length()) ;

  for(int optimIter = 0 ; optimIter < sampMat.nrow() ; optimIter ++) {
    // SAMPLING STARTS
    for(int sampIter = 0 ; sampIter < nsamp ; sampIter++){
      // sampling A0y
      if(methodExact) {
        zeroSampler(zsamp, zeroSamp,
                    zeroMean, sqrtZero,
                    lzero, uzero, 5, 40) ;
      }

      // Sampling Signs
      mhSampler(samp, oldsamp, signs, newsigns,
                estimate, oneCov, condSigma,
                lambda, uone, newuone,
                XmXinv, XmX, ysigsq,
                lzero, uzero, newlzero, newuzero,
                u0mat, zeroSamp,
                1, order, initEst, methodExact) ;

      // Sampling regression Coefs
      aOneSampler(samp, signs, uone, estimate, XmX, condSigma, ysigsq, 10, 10) ;
      copyVector(oldsamp, samp) ;
    }

    // Rcpp::Rcout<<optimIter<<" ";
    // printVec(samp) ;

    // Checking that Sample is reasonable
    for(int l = 0 ; l < samp.length() ; l++) {
      if(std::abs((naive[l] - samp[l]) / std::sqrt(oneCov(l, l))) > 10) {
        for(int i = 0 ; i < samp.length() ; i ++) {
          samp[i] = naive[i] ;
          signs[i] = sign(naive[i]) ;
        }
        copyVector(oldsamp, samp) ;
        copyVector(newsigns, signs) ;
        computeOneThreshold(signs, lambda, XmXinv, uone) ;
        if(methodExact) {
          computeZeroThresholds(u0mat, signs, lzero, uzero) ;
        }
        break ;
      }
    }

    // SAMPLING ENDS

    // OPTIMIZATION STARTS
    // computing gradient and updating estimate
    if(optimIter < assumeConvergence) {
      // Beta Update
      computeGradient(gradient, Xy, samp, XmX,
                      stepCoef, stepRate, gradientBound,
                      optimIter, delay) ;
      for(int i = 0 ; i < gradient.length() ; i++) {
        estimate[i] += gradient[i] ;
      }
      boundBeta(estimate, naive) ;
    }

    // `reporting' sample and estimate
    for(int i = 0 ; i < sampMat.ncol() ; i++) {
      sampMat(optimIter, i) = samp[i] ;
      estimateMat(optimIter, i) = estimate[i] ;
    }

    if(optimIter == (assumeConvergence - 1)) {
      computeConditionalBeta(estimate, estimateMat, optimIter) ;
      burnin *= 2 ;
    }

    // OPTIMIZATION ENDS
    int frac = floor(sampMat.nrow() / 20) ;
    if((optimIter + 1) % frac == 0 & verbose) {
      double percent = round((optimIter + 1.0) / (sampMat.nrow() + 1.0) * 100) ;
      Rcpp::Rcout<<percent<<"\% " ;
      for(int i = 0 ; i < samp.length() ; i ++) {
        samp[i] = naive[i] ;
        signs[i] = sign(naive[i]) ;
      }
      copyVector(oldsamp, samp) ;
      copyVector(newsigns, signs) ;
      computeOneThreshold(signs, lambda, XmXinv, uone) ;
      if(methodExact) {
        computeZeroThresholds(u0mat, signs, lzero, uzero) ;
      }
    }
  }
  if(verbose) {
    Rcpp::Rcout<<"\n" ;
  }
}
