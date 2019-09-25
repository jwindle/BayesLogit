## Copyright 2013 Nick Polson, James Scott, and Jesse Windle.

## This file is part of BayesLogit, distributed under the GNU General Public
## License version 3 or later and without ANY warranty, implied or otherwise.


## h14 = scan("~/Projects/RPackage/BayesLogit/Code/R/h1to4.txt")
## t14 = scan("~/Projects/RPackage/BayesLogit/Code/R/t1to4.txt")
## d14 = scan("~/Projects/RPackage/BayesLogit/Code/R/d1to4.txt")
## e14 = 1/d14

t14 = c(0.64, 0.68, 0.72, 0.75, 0.78, 0.8, 0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95, 0.96, 0.98, 1, 1.01, 1.03, 1.04, 1.06, 1.07, 1.09, 1.1, 1.12, 1.13, 1.15, 1.16, 1.17, 1.19, 1.2, 1.21, 1.23, 1.24, 1.25, 1.26, 1.28, 1.29, 1.3, 1.32, 1.33, 1.34, 1.35, 1.36, 1.38, 1.39, 1.4, 1.41, 1.42, 1.44, 1.45, 1.46, 1.47, 1.48, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.65, 1.66, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.9, 1.91, 1.92, 1.93, 1.95, 1.96, 1.97, 1.98, 1.99, 2, 2.01, 2.02, 2.03, 2.04, 2.05, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.51, 2.52, 2.53, 2.54, 2.55, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.76, 2.77, 2.78, 2.79, 2.8, 2.81, 2.82, 2.83, 2.84, 2.85, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3, 3.01, 3.02, 3.03, 3.04, 3.06, 3.07, 3.08, 3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2, 3.21, 3.22, 3.23, 3.24, 3.25, 3.27, 3.28, 3.29, 3.3, 3.31, 3.32, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.49, 3.5, 3.51, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.58, 3.59, 3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.8, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.95, 3.96, 3.97, 3.98, 3.99, 4, 4.01, 4.02, 4.03, 4.04, 4.05, 4.06, 4.07, 4.08, 4.09, 4.1, 4.11, 4.12, 4.13)


################################################################################
                                ##    J^*    ##
################################################################################

p.ch.left <- function(t, h)
{
  2^h * pgamma(1/t, shape=0.5, rate=0.5*h^2, lower.tail=FALSE)
}

p.ch.right <- function(t, h)
{
  (4/pi)^h * pgamma(t, shape=h, rate=pi^2 * 0.125, lower.tail=FALSE)
}

a.coef.alt <- function(n, x, h)
{
  ## a_n(x,h).
  ## You could precompute lgamma(h), log(2).
  d.n = (2 * n + h)
  l.out = h * log(2) - lgamma(h) + lgamma(n+h) - lgamma(n+1) + log(d.n) -
    0.5 * log(2 * pi * x^3) - 0.5 * d.n^2 / x;
  out = exp(l.out)
  out
}

c.2.coef <- function(n, x)
{
  c.n = (n+1/2) * pi
  out = (1 - 1/(x*c.n^2)) * c.n^2 * x * exp(-c.n^2 * x / 2)
  if (n > 0) out = 2 * out
  out
}

pg.a.coef.alt <- function(n, x, h, z=0)
{
  cosh(z/2)^h * exp(-0.5 * z^2 * x) * 4 * a.coef.alt(n, 4 * x, h)
}

jj.m1 <- function(b,z)
{
  if (z > 1e-12)
    b * tanh(z) / z
  else
    b * (1 - (1/3) * z^2 + (2/15) * z^4 - (17/315) * z^6)
}

jj.m2 <- function(b, z)
{
 
  if (z > 1e-12)
    (b+1) * b * (tanh(z)/z)^2 + b * ((tanh(z)-z)/z^3)
  else
    (b+1) * b * (1 - (1/3) * z^3 + (2/15) * z^4 - (17/315) * z^6)^2 +
      b * ((-1/3) + (2/15) * z - (17/315) * z^3);
}

pg.m1 <- function(b,z)
{
  jj.m1(b,z/2) / 4
}

pg.m2 <- function(b,z)
{
  jj.m2(b,z/2) / 16
}

##------------------------------------------------------------------------------
                          ## SAMPLE TRUNCATED GAMMA ##
##------------------------------------------------------------------------------

## A RELATIVELY EASY METHOD TO IMPLEMENT
rltgamma.dagpunar.1 <- function(shape=1, rate=1, trnc=1)
{
  ## y ~ Ga(shape, rate, trnc)
  ## x = y/t
  ## x ~ Ga(shape, rate t, 1)
  a = shape
  b = rate * trnc

  if (shape <  1) return(NA)
  if (shape == 1) return(rexp(1) / rate + trnc);
  
  d1 = b-a
  d3 = a-1
  c0 = 0.5 * (d1 + sqrt(d1^2 + 4 * b)) / b

  x = 0
  accept = FALSE
  
  while (!accept) {
    x  = b + rexp(1) / c0
    u  = runif(1)

    l.rho = d3 * log(x) - x * (1-c0);
    l.M   = d3 * log(d3 / (1-c0)) - d3

    accept = log(u) <= (l.rho - l.M)
  }

  x = x / b
  y = trnc * x
  y
}

rltgamma.dagpunar <- function(num=1, shape=1, rate=1, trnc=1)
{
  shape = array(shape, dim=num)
  rate  = array(rate , dim=num)
  trnc = array(trnc, dim=num)

  y = rep(0, num)
  for (i in 1:num) y[i] = rltgamma.dagpunar.1(shape[i], rate[i], trnc[i]);

  y
}

rrtinvch2.ch.1 <- function(h, trnc)
{
  R = 1 / (trnc * h^2)
  E = rexp(2)
  while ( (E[1]^2) > (2 * E[2] / R)) {
    ## cat("E", E[1], E[2], E[1]^2, 2*E[2] / R, "\n")
    E = rexp(2)
  }
  ## cat("E", E[1], E[2], "\n")
  ## W^2 = (1 + R*E[1])^2 / R is left truncated chi^2(1) I_{(R,\infty)}.
  ## So X is right truncated inverse chi^2(1).
  X = R / (1 + R*E[1])^2
  X = h^2 * X
  X
}

rrtinvch2.ch <- function(num, h, trnc)
{
  out = rep(0, num)
  for (i in 1:num)
    out[i] = rrtinvch2.ch.1(h, trnc)
  out
}

##------------------------------------------------------------------------------

a.coef.alt.1 <- function(n,x, trnc)
{
  if ( x>trnc )
    pi * (n+0.5) * exp( -(n+0.5)^2*pi^2*x/2 )
  else
    (2/pi/x)^1.5 * pi * (n+0.5) * exp( -2*(n+0.5)^2/x )
}

g.tilde <- function(x, h, trnc)
{
  if (x > trnc)
    (0.5 * pi)^h * x^(h-1) / gamma(h) * exp(-pi^2 / 8 * x)
  else
    2^h * h * (2 * pi * x^3)^(-0.5) * exp(-0.5 * h^2 / x)
}

rch.1 <- function(h)
{
  if (h < 1) return(NA);
  
  rate  = pi^2 / 8
  idx   = floor((h-1) * 100) + 1;
  trnc  = t14[idx];
  
  num.trials = 0;
  total.iter = 0;

  p.l = p.ch.left (trnc, h);
  p.r = p.ch.right(trnc, h);
  p   = p.r / (p.l + p.r);

  max.inner = 100
      
  while (TRUE)
    {
      num.trials = num.trials + 1;

      if ( runif(1) < p ) {
        ## Left truncated gamma
        X = rltgamma.dagpunar.1(shape=h, rate=rate, trnc=trnc)
      }
      else {
        ## Right truncated inverse Chi^2
        X = rrtinvch2.ch.1(h, trnc)
      }

      ## C = cosh(Z) * exp( -0.5 * Z^2 * X )
      ## Don't need to multiply everything by C, since it cancels in inequality.

      S = a.coef.alt(0,X,h)
      ## B = a.coef.alt.1(0, X, trnc)
      D = g.tilde(X, h, trnc)
      Y = runif(1) * D
      n = 0

      ## cat("B,C,left", B, C, X<trnc, "\n")

      a.n = S
      decreasing = FALSE

      while (n < max.inner)
        {
          n = n + 1
          total.iter = total.iter + 1;
          
          a.prev = a.n
          a.n = a.coef.alt(n, X, h)
          ## b.n = a.coef.alt.1(n, X, trnc)
          decreasing = a.n < a.prev
          ## if (!decreasing) cat("n:", n, "; ");
          
          if ( n %% 2 == 1 )
            {
              S = S - a.n
              ## B = B - b.n
              if ( Y<=S && decreasing) break
            }
          else
            {
              S = S + a.n
              ## B = B + b.n
              if ( Y>S && decreasing) break
            }
        }

      ## cat("S,B =", S, B, "\n")

      if ( Y<=S ) break     
    }

  out = list("x"=X, "n"=num.trials, "total.iter"=total.iter)
  out
}

rch <- function(num, h)
{
  h = array(h, num)
  x = rep(0, num)
  for (i in 1:num) {
    out = rch.1(h[i])
    x[i] = out$x
  }
  x
}


################################################################################
                            ## PLOTTING DENSITIES ##
################################################################################

if (FALSE)
{

  ## a_n is absolutely summable.  When h = 2 sum |a_n| = 0.5 as x -> infty.
  
  ## source("Ch.R")
  dx    = 0.1
  xgrid = seq(dx, 10, dx)
  y1    = xgrid
  y2    = xgrid
  n     = length(xgrid)

  for (i in 1:n) {
    y1[i] = 0
    y2[i] = 0
    for (j in 0:200) {
      y1[i] = y1[i] + (-1)^j * a.coef.alt(j,xgrid[i],2)
      y2[i] = y2[i] + c.2.coef(j,xgrid[i])
    }
    
  }

  plot(xgrid, y1)
  points(xgrid, y2, col=2)
  
}

##------------------------------------------------------------------------------
                           ## PLOTTING PG DENSITY ##
##------------------------------------------------------------------------------

if (FALSE)
{
  ## source("~/Projects/RPackage/BayesLogit/Code/R/Ch.R")
  dx    = 0.01
  xgrid = seq(dx, 3, dx)
  y1    = xgrid
  y2    = xgrid
  y3    = xgrid
  y4    = xgrid
  y5    = xgrid
  n     = length(xgrid)

  for (i in 1:n) {
    y1[i] = 0
    y2[i] = 0
    y3[i] = 0
    y4[i] = 0
    y5[i] = 0
    for (j in 0:200) {
      y1[i] = y1[i] + (-1)^j * pg.a.coef.alt(j,xgrid[i],1)
      y2[i] = y2[i] + (-1)^j * pg.a.coef.alt(j,xgrid[i],2)
      y3[i] = y3[i] + (-1)^j * pg.a.coef.alt(j,xgrid[i],3)
      y4[i] = y4[i] + (-1)^j * pg.a.coef.alt(j,xgrid[i],1, 2.0)
      y5[i] = y5[i] + (-1)^j * pg.a.coef.alt(j,xgrid[i],1, 4.0)
    }
    
  }

  ## png(filename="pg-dens.png", width=800, height=400)
  
  par(mfrow=c(1,2))
  ymax = max(c(y1,y2,y3))
  plot(xgrid, y1, type="l", ylim=c(0,ymax), main="Density of PG(b,0)", xlab="x", ylab="f(x|b,0)")
  lines(xgrid, y2, type="l", col=2, lty=2)
  lines(xgrid, y3, type="l", col=4, lty=4)

  legend("topright", legend=c("b=1", "b=2", "b=3"), col=c(1,2,4), lty=c(1,2,4))

  ## hist(rpg.devroye(10000, 1, 0), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(0, 0, 0, 16, maxColorValue=255))
  ## hist(rpg.devroye(10000, 2, 0), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(255, 0, 0, 16, maxColorValue=255))
  ## hist(rpg.devroye(10000, 3, 0), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(0, 0, 255, 16, maxColorValue=255))

  ymax = max(c(y1,y4,y5))
  plot(xgrid, y1, type="l", ylim=c(0,ymax), xlim=c(0,1),
       main=expression(paste("Density of PG(1,", psi, ")", sep="")),
       xlab="x",
       ylab=expression(paste("f(x|1,", psi, ")", sep="")))
  lines(xgrid, y4, type="l", col=2, lty=2)
  lines(xgrid, y5, type="l", col=4, lty=4)

  legend("topright",
         legend=c(expression(paste(psi, "=0", sep="")),
                  expression(paste(psi, "=2", sep="")),
                  expression(paste(psi, "=4", sep=""))),
         col=c(1,2,4), lty=c(1,2,4))

  ## hist(rpg.devroye(10000, 1, 0), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(0, 0, 0, 16, maxColorValue=255))
  ## hist(rpg.devroye(10000, 1, 2), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(255, 0, 0, 16, maxColorValue=255))
  ## hist(rpg.devroye(10000, 1, 4), add=TRUE, prob=TRUE, breaks=100,
  ##      col=rgb(0, 0, 255, 16, maxColorValue=255))

  dev.off()
  
}

################################################################################
                                ## TILTED J^* ##
################################################################################

## pigauss - cumulative distribution function for Inv-Gauss(mu, lambda).
##------------------------------------------------------------------------------
pigauss <- function(x, Z=1, lambda=1)
{
  ## I believe this works when Z = 0
  ## Z = 1/mu
  b = sqrt(lambda / x) * (x * Z - 1);
  a = sqrt(lambda / x) * (x * Z + 1) * -1.0;
  y = exp(pnorm(b, log.p=TRUE)) + exp(2 * lambda * Z + pnorm(a, log.p=TRUE));
                                        # y2 = 2 * pnorm(-1.0 / sqrt(x));
  y
}

p.tilt.left <- function(trnc, h, z)
{
  out = 0
  if (z == 0) {
    out = p.ch.left(trnc, h)
  }
  else {
    out = (2^h * exp(-z*h)) * pigauss(trnc, Z=z/h, h^2)
  }
  out
}

p.tilt.right <- function(trcn, h, z)
{
  ## Note: this works when z=0
  lambda.z = pi^2/8 + z^2/2
  (pi/2/lambda.z)^h * pgamma(trcn, shape=h, rate=lambda.z, lower.tail=FALSE)
}

## rigauss - sample from Inv-Gauss(mu, lambda).
##------------------------------------------------------------------------------
rigauss <- function(mu, lambda)
{
  nu = rnorm(1);
  y  = nu^2;
  x  = mu + 0.5 * mu^2 * y / lambda -
    0.5 * mu / lambda * sqrt(4 * mu * lambda * y + (mu*y)^2);
  if (runif(1) > mu / (mu + x)) {
    x = mu^2 / x;
  }
  x
}

rrtigauss.ch.1 <- function(h, z, trnc=1)
{
  ## trnc is truncation point
  z = abs(z);
  mu = h/z;
  X = trnc + 1;
  if (mu > trnc) {
    alpha = 0.0;
    while (runif(1) > alpha) {
      X = rrtinvch2.ch.1(h, trnc)
      alpha = exp(-0.5 * z^2 * X);
    }
    ## cat("rtigauss.ch, part i:", X, "\n");
  }
  else {
    while (X > trnc) {
      lambda = h^2;
      Y = rnorm(1)^2;
      X = mu + 0.5 * mu^2 / lambda * Y -
        0.5 * mu / lambda * sqrt(4 * mu * lambda * Y + (mu * Y)^2);
      if ( runif(1) > mu / (mu + X) ) {
        X = mu^2 / X;
      }
    }
    ## cat("rtiguass, part ii:", X, "\n");
  }
  X;
}

rrtigauss.ch <- function(num, h, z, trnc=1)
{
  x = rep(0, num)
  for (i in 1:num)
    x[i] = rrtigauss.ch.1(h, z, trnc)
  x
}

rpg.alt.1 <- function(h, z)
{
  z = z/2
  if (h < 1 || h > 4) return(NA);
  
  rate  = pi^2 / 8 + z^2 / 2
  idx   = floor((h-1) * 100) + 1;
  trnc  = t14[idx];
  
  p.l = p.tilt.left (trnc, h, z);
  p.r = p.tilt.right(trnc, h, z);
  p   = p.r / (p.l + p.r);

  ## cat("prob.right:", p, "\n")

  num.trials = 0;
  total.iter = 0;
  
  max.outer = 1000
  max.inner = 1000
      
  while (num.trials < max.outer)
    {
      num.trials = num.trials + 1;

      uu = runif(1)
      if ( uu < p ) {
        ## Left truncated gamma
        X = rltgamma.dagpunar.1(shape=h, rate=rate, trnc=trnc)
      }
      else {
        ## Right truncated inverse Chi^2
        ## Note: this sampler works when z=0.
        X = rrtigauss.ch.1(h, z, trnc)
      }

      ## C = cosh(Z) * exp( -0.5 * Z^2 * X )
      ## Don't need to multiply everything by C, since it cancels in inequality.

      S = a.coef.alt(0,X,h)
      ## S = a.coef.alt.1(0, X, trnc)
      D = g.tilde(X, h, trnc)
      Y = runif(1) * D
      n = 0

      ## cat("B,C,left", B, C, X<trnc, "\n")
      ## cat("test gt:", g.tilde(trnc * 0.1, h, trnc), "\n");
      ## cat("X, Y, S, gt:", X, Y, S, D,"\n");
      
      a.n = S
      decreasing = FALSE
      go = TRUE

      while (go && n < max.inner)
        {
          n = n + 1
          total.iter = total.iter + 1;
          
          a.prev = a.n
          a.n = a.coef.alt(n, X, h)
          ## a.n = a.coef.alt.1(n, X, trnc)
          decreasing = a.n <= a.prev
          ## if (!decreasing) cat("n:", n, "; ");
          
          if ( n %% 2 == 1 )
            {
              S = S - a.n
              ## B = B - b.n
              if ( Y<=S && decreasing) break
              ## if ( Y<=S && decreasing) return(0.25 * X)
            }
          else
            {
              S = S + a.n
              ## B = B + b.n
              if ( Y>S && decreasing) go = break
            }
        }

      ## cat("S,B =", S, B, "\n")

      ## Need to check max.outer
      if ( Y<=S ) break     
    }

  X = 0.25 * X
  out = list("x"=X, "num.trials"=num.trials, "total.iter"=total.iter)

  out
}

## rpg.alt.4 <- function(num, h, z)
## {
##   if (h < 1) return(NA)

##   z = array(z, num)
##   h = array(h, num)

##   out = list(
##     draw = rep(0, num),
##     num.trials = rep(0, num),
##     total.iter = rep(0, num)
##     )
  
##   for (i in 1:num) {
##     draw = rpg.alt.1(h[i], z[i])
##     out$draw[i] = draw$x
##     out$num.trials[i] = draw$num.trials
##     out$total.iter[i] = draw$total.iter
##   }

##   out
## }


rpg.alt.R <- function(num, h, z)
{
    if (any(h < 1)) { return(NA) }
    
    h = array(h, num)
    z = array(z, num)

    n = floor( (h-1.) / 4. )
    remain = h - 4. * n
    
    for (i in 1:num) {
        x = 0.0
        for (j in 1:n[i]) {
            draw = rpg.alt.1(4.0, z[i])
            x = x + draw$x
        }
        if (remain[i] > 4.0) {
            draw = rpg.alt.1(0.5 * remain[i], z[i]) + rpg.alt.1(0.5 * remain[i], z[i]) 
            x = x + draw$x
        } else {
            draw = rpg.alt.1(remain[i], z[i])
            x = x + draw$x
        }
        out$draw[i] = x
        out$num.trials[i] = draw$num.trials
        out$total.iter[i] = draw$total.iter
    }

    out$draw
}


################################################################################
                                ## CHECK rpg ##
################################################################################

if (FALSE)
{

  ## source("Ch.R")
  h = 2.3
  z = 1.1

  num = 20000

  samp.1 = rpg.4(num, h, z)
  samp.2 = rpg.gamma(num, h, z)

  mean(samp.1$draw)
  mean(samp.2)

  mean(samp.1$total.iter)

  hist(samp.1$draw, prob=TRUE, breaks=100)
  hist(samp.2, prob=TRUE, add=TRUE, col="#22000022", breaks=100)
  
}

if (FALSE) {

  ## source("Ch.R")
  source("ManualLoad.R")

  nsamp = 10000
  n = 1
  z = 0

  seed = sample.int(10000, 1)
  ## seed = 8922

  set.seed(seed)
  samp.a = rpg.alt(nsamp, n, z)
  ## samp.d = rpg.devroye(nsamp, n, z)
  set.seed(seed)
  samp.4 = rpg.4(nsamp, n, z)
  
  mean(samp.a)
  ## mean(samp.d)
  mean(samp.4$draw)

  hist(samp.a, prob=TRUE, breaks=100)
  hist(samp.4$draw, prob=TRUE, breaks=100, col="#99000022", add=TRUE)
  
  h = 1.5
  z = 0

  set.seed(seed)
  samp.a = rpg.alt(nsamp, h, z)
  ## samp.g = rpg.gamma(nsamp, h, z)
  set.seed(seed)
  samp.4 = rpg.4(nsamp, h, z)

  mean(samp.a)
  ## mean(samp.g)
  mean(samp.4$draw)

  hist(samp.a, prob=TRUE, breaks=100)
  hist(samp.4$draw, prob=TRUE, breaks=100, col="#99000022", add=TRUE)
  
}

if (FALSE)
{

  ## source("Ch.R")
  source("ManualLoad.R")

  reps   = 10
  nsamp = 10000
  n = 4
  z = 0

  seed = sample.int(10000, 1)
  ## seed = 8922

  set.seed(seed)
  start.a = proc.time()
  for (i in 1:reps) {
    samp.a = rpg(nsamp, n, z)
  }
  end.a = proc.time();
  diff.a = end.a - start.a

  set.seed(seed)
  start.d = proc.time()
  for (i in 1:reps) {
    samp.d = rpg.devroye(nsamp, n, z)
  }
  end.d = proc.time()
  diff.d = end.d - start.d

  mean(samp.a)
  mean(samp.d)

  diff.a
  diff.d

  ## hist(samp.a, prob=TRUE, breaks=100)
  ## hist(samp.4$draw, prob=TRUE, breaks=100, col="#99000022", add=TRUE)
  
  h = 3.5
  z = 0

 set.seed(seed)
  start.a = proc.time()
  for (i in 1:reps) {
    samp.a = rpg(nsamp, h, z)
  }
  end.a = proc.time();
  diff.a = end.a - start.a

  set.seed(seed)
  start.g = proc.time()
  for (i in 1:reps) {
    samp.g = rpg.gamma(nsamp, h, z)
  }
  end.g = proc.time()
  diff.g = end.g - start.g

  mean(samp.a)
  mean(samp.g)

  diff.a
  diff.g

  ## hist(samp.a, prob=TRUE, breaks=100)
  ## hist(samp.4$draw, prob=TRUE, breaks=100, col="#99000022", add=TRUE)
  
}

################################################################################

################################################################################

if (FALSE)
{
  ## Preliminary: approximating using normal.
  
  ## source("Ch.R")
  dx    = 0.1
  xgrid = seq(dx, 30, dx)
  y1    = xgrid
  y2    = xgrid
  n     = length(xgrid)
  h     = 20
  
  for (i in 1:n) {
    y1[i] = 0
    y2[i] = 0
    for (j in 0:200) {
      y1[i] = y1[i] + (-1)^j * a.coef.alt(j,xgrid[i],h)
    }
  }

  m1 = jj.m1(h, 0)
  m2 = jj.m2(h, 0)
  V  = m2 - m1^2;
  y2 = dnorm(xgrid, m1, sqrt(V));
  
  ## png(filename="pg-dens.png", width=800, height=400)
  par(mfrow=c(1,2))
  ymax = max(c(y1,y2))
  plot(xgrid, y1, type="l", ylim=c(0,ymax), main="Density of PG(b,0)", xlab="x", ylab="f(x|b,0)")
  lines(xgrid, y2, type="l", col=2, lty=2)

  
}
