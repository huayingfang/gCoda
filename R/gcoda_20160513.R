################################################################################
# File: gcoda.R
# Aim : Conditional dependence network inference for compositional data
#-------------------------------------------------------------------------------
# Author: Fang Huaying (Peking University)
# Email : hyfang@pku.edu.cn
# Date  : 13MAY2016
#-------------------------------------------------------------------------------
require(huge);
#-------------------------------------------------------------------------------
gcoda <- function(x, counts = F, pseudo = 0.5, lambda.min.ratio = 1e-4, 
                  nlambda = 15, ebic.gamma = 0.5) {
    # Counts or fractions?
    if(counts) {
      x <- x + pseudo;
      x <- x / rowSums(x);
    }
    n <- nrow(x);
    p <- ncol(x);
    # Log transformation for compositional data
    S <- var(log(x) - rowMeans(log(x)));
    # Generate lambda via lambda.min.ratio and nlambda
    lambda.max <- max(max(S - diag(p)), -min(S - diag(p)));
    lambda.min <- lambda.min.ratio * lambda.max;
    lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda));
    # Store fit result for gcoda via a series of lambda
    fit <- list();
    fit$lambda <- lambda;
    fit$nloglik <- rep(0, nlambda);
    fit$df <- rep(0, nlambda);
    fit$path <- list();
    fit$icov <- list();
    # Compute solution paths for gcoda
    icov <- diag(p);
    for(i in 1:nlambda) {
        out.gcoda <- gcoda_sub(A = S, iSig = icov, lambda = lambda[i]);
        icov <- out.gcoda$iSig;
        fit$nloglik[i] <- out.gcoda$nloglik;
        fit$icov[[i]] <- icov;
        fit$path[[i]] <- 0 + (abs(icov) > 1e-6);
        diag(fit$path[[i]]) <- 0;
        fit$df[i] <- sum(fit$path[[i]]) / 2;
    }
    # Continue if edge density is too small
    imax <- nlambda + 15;
    emin <- p * (p - 1) / 2 * 0.618;
    while(fit$df[i] <= emin && i <= imax && lambda[i] > 1e-6) {
      cat(i, "Going on!\n");
      lambda <- c(lambda, lambda[i]/2);
      i <- i + 1;
      fit$lambda[i] <- lambda[i]; 
      out.gcoda <- gcoda_sub(A = S, iSig = icov, lambda = lambda[i]);
      icov <- out.gcoda$iSig;
      fit$nloglik[i] <- out.gcoda$nloglik;
      fit$icov[[i]] <- icov;
      fit$path[[i]] <- 0 + (abs(icov) > 1e-6);
      diag(fit$path[[i]]) <- 0;
      fit$df[i] <- sum(fit$path[[i]]) / 2;
    }
    # Compute EBIC score for lambda selection
    fit$ebic.score <- n * fit$nloglik + log(n) * fit$df + 
      4 * ebic.gamma * log(p) * fit$df;
    fit$opt.index <- which.min(fit$ebic.score);
    fit$refit <- fit$path[[fit$opt.index]];
    fit$opt.icov <- fit$icov[[fit$opt.index]];
    fit$opt.lambda <- fit$lambda[fit$opt.index];
    return(fit);
}
#-------------------------------------------------------------------------------
# Optimization for gcoda with given lambda
gcoda_sub <- function(A, iSig = NULL, lambda = 0.1, tol_err = 1e-4, 
                      k_max = 100) {
  p <- ncol(A);
  if(is.null(iSig)) {
    iSig <- diag(p);
  }

  err <- 1;
  k <- 0;
  fval_cur <- Inf;
  
  while(err > tol_err && k < k_max) {
    iSig_O <- rowSums(iSig);
    iS_iSig <- 1 / sum(iSig_O);
    iSig_O2 <- iSig_O * iS_iSig;
    A_iSig_O2 <- rowSums(A * rep(iSig_O2, each = p));    
    A2 <- A - A_iSig_O2 - rep(A_iSig_O2, each = p) + 
          sum(iSig_O2 * A_iSig_O2) + iS_iSig;
    iSig2 <- huge_glasso_mod(S = A2, lambda = lambda);

    fval_new <- obj_gcoda(iSig = iSig2, A = A, lambda = lambda);
    xerr <- max(abs(iSig2 - iSig) / (abs(iSig2) + 1));
    err <- min(xerr, abs(fval_cur - fval_new)/(abs(fval_new) + 1));

    k <- k + 1;
    iSig <- iSig2;
    fval_cur <- fval_new;
  }
  nloglik <- fval_cur - lambda * sum(abs(iSig));

  if(k >= k_max) {
    cat("WARNING of gcoda_sub:\n", "\tMaximum Iteration:", k_max, 
      "&& Relative error:", err, "!\n");
  }
  
  return(list(iSig = iSig, nloglik = nloglik));
}
#----------------------------------------
# Objective function value of gcoda (negative log likelihood + penalty)
obj_gcoda <- function(iSig, A, lambda) {
  p <- ncol(A);
  iSig_O <- rowSums(iSig);
  S_iSig <- sum(iSig_O);
  nloglik <- - log(det(iSig)) + sum(iSig * A) + log(S_iSig) - 
            sum(iSig_O * rowSums(A * rep(iSig_O, each = p))) / S_iSig;
  pen <- lambda * sum(abs(iSig));
  return(nloglik + pen);
}
#----------------------------------------
# Modified huge::huge.glasso for quick preparation
# Input S must be covariance matrix
require(huge);
huge_glasso_mod <- function(S, lambda) {
    icov <- diag(1/(diag(S) + lambda));

    z <- which(rowSums(abs(S) > lambda) > 1);
    q <- length(z);
    if (q > 0) {
        out.glasso = .C("hugeglasso", S = as.double(S[z, z]), 
                        W = as.double(S[z, z]), T = as.double(diag(q)), 
                        dd = as.integer(q), lambda = as.double(lambda), 
                        df = as.integer(0), PACKAGE = "huge");
        icov[z, z] = matrix(out.glasso$T, ncol = q);
    }

    return(icov);
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

