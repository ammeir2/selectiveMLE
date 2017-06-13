#include <Rcpp.h>
using namespace Rcpp;
double sampleUnivTruncNorm(double mu, double sd, double threshold);

double computeConditionalMean(NumericVector mu,
                              NumericVector samp,
                              const NumericMatrix XmX,
                              double yvar,
                              int index);
