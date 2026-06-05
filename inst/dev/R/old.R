## Old stuff


## Bayesian logistic regression
##------------------------------------------------------------------------------
logit.R <- function(y, X, n=rep(1, length(y)),
                    m0=rep(0, ncol(X)), P0=matrix(0, nrow=ncol(X), ncol=ncol(X)),
                    samp=1000, burn=500, verbose=500)
{
  ## X: n by p matrix
  ## y: n by 1 vector, avg response
  ## n: n by 1 vector, # of obs at distinct x

  ## Combine data.
  ## new.data = logit.combine(y, X, n);
  ## y = new.data$y;
  ## X = new.data$X;
  ## n = new.data$n;
  ## n.prior = 0.0;

  X = as.matrix(X);
  y = as.numeric(y)

  p = ncol(X)
  N = nrow(X)

  alpha = (y-1/2)*n

  Z = colSums(X*alpha) + P0 %*% m0;
  ## PsiToBeta = solve(t(X) %*% X) %*% t(X);

  w = rep(0,N)
  ## w = w.known;
  beta = rep(0.0, p)

  output <- list(w = matrix(nrow=samp, ncol=N),
                 beta = matrix(nrow=samp, ncol=p)
                 )

  ## c_k = (1:200-1/2)^2 * pi^2 * 4;

  ## Timing
  start.time = proc.time()
  
  ## Sample
  for ( j in 1:(samp+burn) )
  {
    if (j==burn+1) start.ess = proc.time();
    
    ## draw w
    psi = drop(X%*%beta)
    ## Sum of gamma: poor approximation when psi is large!  Causes crash.
    ## w = rpg.gamma(N, n, psi)
    ## Devroye is faster anyway.
    w = rpg.devroye(N, n, psi);

    ## draw beta - Joint Sample.
    PP = t(X) %*% (X * w) + P0;
    ## U = chol(PP);
    ## m = backsolve(U, Z, transpose=TRUE);
    ## m = backsolve(U, m);
    ## beta = m + backsolve(U, rnorm(p))
    S = chol2inv(chol(PP));
    m = S %*% as.vector(Z);
    beta = m + t(chol(S)) %*% rnorm(p);

    # Record if we are past burn-in.
    if (j>burn) {
        output$w[j-burn,] <- w
        output$beta[j-burn,] <- beta
    }

    if (j %% verbose == 0) { print(paste("LogitPG: Iteration", j)); }
  }

  end.time = proc.time()
  output$total.time = end.time - start.time
  output$ess.time   = end.time - start.ess

  ## Add new data to output.
  output$"y" = y;
  output$"X" = X;
  output$"n" = n;

  output
} ## logit.gibbs.R