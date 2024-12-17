# Obtain statistics for linear regression

# we need to correct the F test to differentiate between the model with intercept + covariate and the model with intercept, covariate and rhythmicity.
do_regression_cov <- function(Y, X) {
  betas <- matrix(NA, nrow=ncol(X), ncol=nrow(Y))
    lower_bounds <- matrix(NA, nrow = ncol(X), ncol = nrow(Y))
    upper_bounds <- matrix(NA, nrow = ncol(X), ncol = nrow(Y))
  dfmod <- numeric(nrow(Y))
  dfres <- numeric(nrow(Y))
  sstot <- numeric(nrow(Y))
  ssres <- numeric(nrow(Y))
  rsqs  <- numeric(nrow(Y))
  fs    <- numeric(nrow(Y))

  nainds <- is.na(Y)
  groups <- rowSums(nainds)
  if(any(groups != 0)) {
    categ <- lapply(split(nainds[groups!=0,,drop=FALSE], row(nainds[groups!=0,,drop=FALSE])), which)
    groups[groups!=0] <- match(categ, unique(categ))
  }
  for(g in unique(groups)) {
    rowinds <- groups == g
    colinds <- !nainds[match(g, groups),]
      print(sum(rowinds))
      print(sum(colinds))
    y   <- Y[rowinds,colinds,drop=FALSE]
    res <- stats::.lm.fit(X[colinds,,drop=FALSE], t(y))

    betas[,rowinds] <- res$coefficients
    betas[,rowinds][abs(betas[,rowinds]) < .Machine$double.eps] <- 0

    dfres[rowinds] <- sum(colinds) - res$rank
    dfmod[rowinds] <- res$rank - 1

    sstot[rowinds] <- rowSums((y - rowMeans(y, na.rm=TRUE))^2, na.rm=TRUE)
    sstot[rowinds][sstot[rowinds] < .Machine$double.eps] <- 0
    ssres[rowinds] <- colSums(res$residuals^2, na.rm=TRUE)
    ssres[rowinds][ssres[rowinds] < .Machine$double.eps] <- 0
    isequal <- abs(ssres[rowinds]-sstot[rowinds]) < .Machine$double.eps^0.5
    ssres[rowinds][isequal] <- sstot[rowinds][isequal]
    rsqs[rowinds] <- 1 - (ssres[rowinds]/sstot[rowinds])

    ssmod <- sstot[rowinds] - ssres[rowinds]
    msres <- ssres[rowinds] / dfres[rowinds]
    msmod <- ssmod / dfmod[rowinds]
    fs[rowinds] <- msmod / msres
      
      # confidence intervals for coefficients
      #print(dim(X))
      #print(head(X))
      
      # the X^tX is gonna be the same for each gene, and has dimensions ncoefs X ncoefs
      XTX_inv <- solve(t(X[colinds,,drop=FALSE]) %*% X[colinds,,drop=FALSE])
      #print(dim(XTX_inv))
      #print(length(msres))
      
      # the error comes from here: we have to multiply XTX_inv by each element in mrses, which contains
      # one residual_variance  value for each estimation.
      # standard errors needs to be a ngenes times ncoefs list or matrix.
      
      # there are as many standard errors as beta coeffs. 
      # the vector of standard errors needs to be the same dimension as X, which is nrow = genes, ncols = betas
      standard_errors <- sapply(msres, function(x) t(sqrt(diag(x*XTX_inv))))
      t_critical <- qt(0.975, dfres)
      
      lower_bounds[, rowinds] <- betas[,rowinds] - t_critical * standard_errors
    upper_bounds[, rowinds] <- betas[,rowinds] + t_critical * standard_errors

    
      if (ncol(X) == 4){
          res_reduced = stats::.lm.fit(X[colinds,c(1, 4),drop=FALSE], t(y))
          ssres_reduced = colSums(res_reduced$residuals^2, na.rm=TRUE)
          
          ssmod =  ssres_reduced[rowinds] - ssres[rowinds]
          p = 2 # number of coefficients removed from the full model
          ssmod = ssres_reduced - ssres
          dfmod = 2
          msmod = ssmod/dfmod
          fs[rowinds] <- msmod / msres
      }
      
  }

  ps <- stats::pf(fs, dfmod, dfres, lower.tail=FALSE)

  list(betas=betas, lower_bounds = lower_bounds, upper_bounds = upper_bounds,
       stats=data.frame(dfmod=dfmod, dfres=dfres, sstot=sstot, ssres=ssres, rsq=rsqs, f=fs, p=ps))
}
