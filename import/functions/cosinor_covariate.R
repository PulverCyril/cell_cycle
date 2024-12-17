## originally from matrixTests/R/cosinor.R

row_cosinor_cov <- function(x, t, period=24, cov = NULL) {
  force(x)
  force(t)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)
  assert_numeric_vec_length(t, ncol(x))
  assert_numeric_vec_length(period, 1)
  assert_all_in_open_interval(period, 0, Inf)

  bad <- is.na(t)
  if(any(bad)) {
    warning(sum(bad), ' columns dropped due to missing time information')
    x <- x[,!bad, drop=FALSE]
    t <- t[!bad]
  }

  period <- rep_len(period, min(1, nrow(x)))

  hasinfx <- is.infinite(x)
  x[hasinfx] <- NA
  hasinfx <- rowSums(hasinfx) > 0


  nobs  <- rowSums(!is.na(x))

  b0 <- rep(1, ncol(x))
  b1 <- sinpi(2*t/period)
  b2 <- cospi(2*t/period)
  B  <- cbind(b0, b1, b2, cov)

    #print("B, mapping to X in do_regression_cov: ")
    #print(dim(B))
    #print(head(B))
  res <- do_regression_cov(x, B)

  mesor     <- res$betas[1,]
  phase     <- atan2(-res$betas[2,], res$betas[3,])
  acrophase <- (-1*phase / (2*pi) * period + period) %% period
  amplitude <- sqrt(res$betas[2,]^2 + res$betas[3,]^2)
    
    beta <- res$betas[2, ]
    gamma <- res$betas[3, ]
    delta = rep(NA, times = length(beta))
    if(!is.null(cov)) {
        delta <- res$betas[4, ]
        }

    beta_2.5 <- res$lower_bounds[2, ]
    beta_97.5 <- res$upper_bounds[2, ]
    
    gamma_2.5 <- res$lower_bounds[3, ]
    gamma_97.5 <- res$upper_bounds[3, ]


    phase_2.5     <- atan2(-res$lower_bounds[2,], res$lower_bounds[3,])
    acrophase_2.5 <- (phase_2.5 / (2*pi) * period + period) %% period

    phase_97.5     <- atan2(-res$upper_bounds[2,], res$upper_bounds[3,])
    acrophase_97.5 <- (phase_97.5 / (2*pi) * period + period) %% period
    
    phase_lower_upper <- atan2(-res$lower_bounds[2,], res$upper_bounds[3,])
    phase_upper_lower <- atan2(-res$upper_bounds[2,], res$lower_bounds[3,])

    amplitude_2.5 <- sqrt(res$lower_bounds[2,]^2 + res$lower_bounds[3,]^2)
    amplitude_97.5 <- sqrt(res$upper_bounds[2,]^2 + res$upper_bounds[3,]^2)



  w1 <- hasinfx
  showWarning(w1, 'had infinite observations that were removed')

  w2 <- nobs < 3
  showWarning(w2, 'had less than 3 complete observations: no p-values produced, amplitude and acrophase will be unreliable')

  w3 <- nobs == 3
  showWarning(w3, 'had exactly 3 complete observations: no p-values produced')

  w4 <- !w2 & !w3 & res$stats$dfmod < 2
  showWarning(w4, 'had less than 3 unique timepoints within the specified period: amplitude and acrophase will be unreliable')

  w5 <- !w2 & !w3 & res$stats$rsq == 1
  showWarning(w5, 'had essentially perfect fit')

  w6 <- !w2 & !w3 & res$stats$sstot == 0
  showWarning(w6, 'were essentially constant')

  res$stats[w2 | w3 | w6, c("f","p")] <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  res_df <- data.frame(obs=nobs, mesor=mesor, amplitude=amplitude, acrophase=acrophase,
             rsquared=res$stats$rsq, df.model=res$stats$dfmod,
             df.residual=res$stats$dfres, statistic=res$stats$f,
             pvalue=res$stats$p, period=period, row.names=rnames,
                       beta = beta, gamma = gamma, delta = delta,
                        beta_2.5 = beta_2.5,
                      beta_97.5 = beta_97.5,
                       gamma_2.5 = gamma_2.5,
                       gamma_97.5 = gamma_97.5,
                       amplitude_2.5 = amplitude_2.5,
                       amplitude_97.5 = amplitude_97.5,
                       acrophase_2.5 = acrophase_2.5,
                       acrophase_97.5 = acrophase_97.5,
                       phase = phase,
                       phase_2.5 = phase_2.5,
                       phase_97.5 = phase_97.5,
                       phase_lower_upper = phase_lower_upper,
                       phase_upper_lower = phase_upper_lower

             )
    return(res_df)
}
