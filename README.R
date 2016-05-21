################################################################################
# File: README.R
# Aim : A brief introduction about the usage for gCoda
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 2015-05-21
#-------------------------------------------------------------------------------
# Package required: 
#                   huge
# Files needed: 
#               gcoda.R for gcoda 
#-------------------------------------------------------------------------------
# Function parameter description:
# function: gcoda
#   Input:
#                 x ------ n x p data matrix (row/column is sample/variable)
#                          n samples & p compositional variables
#            counts ------ Is the compositional data matrix a count matrix? 
#                          Default: FALSE
#            pseudo ------ pseudo count if counts = TRUE
#                          Default: 0.5
#  lambda.min.ratio ------ lambda(max) / lambda(min)
#                          Default: 1e-4     
#           nlambda ------ number of tuning parameters
#                          Default: 15
#        ebic.gamma ------ gamma value of EBIC
#                          Default: 0.5
#   Output: 
#      A list structure contains:
#            lambda ------ lambda sequence for compuation of EBIC score
#           nloglik ------ negative log likelihood for lambda sequence
#                df ------ number of edges for lambda sequence
#              path ------ sparse pattern for lambda sequence
#              icov ------ inverse covariance matrix for lambda sequence
#        ebic.score ------ EBIC score for lambda sequence
#             refit ------ sparse pattern with best EBIC score 
#          opt.icov ------ inverse covariance matrix with best EBIC score
#        opt.lambda ------ lambda with best EBIC score
#-------------------------------------------------------------------------------
require(huge);
source("R/gCoda.R");
#-------------------------------------------------------------------------------
# 1 Basic example (no edges in the conditional dependence network)
# 1.1 Generate logistic normal variables
n <- 100;
p <- 20;
x <- matrix(rnorm(n * p), nrow = n); 
x.frac <- exp(x) / rowSums(exp((x)));
totCount <- round(runif(n = n,  min = 1000, max = 2000));
x.count <- x.frac * totCount;
# 1.2 Run gCoda 
# using fraction
res_gcoda_frac <- gcoda(x = x.frac, counts = F);
# using counts
res_gcoda_count <- gcoda(x = x.count, counts = T);
# 1.3 Get the estimation of the inverse covariance matrix
{
  cat("gCoda using fraction data:\n");
  print(round(res_gcoda_frac$opt.icov, 2));
  cat("gCoda using count data:\n");
  print(round(res_gcoda_count$opt.icov, 2));
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
