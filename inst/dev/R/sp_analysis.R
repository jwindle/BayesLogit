## Analysis and plotting code for the saddlepoint approximation to PG(h, z).
## Not part of the package — source() this file directly for exploration.
##
## Requires the package to be loaded: library(BayesLogit)
## Several functions below also need SaddlePointApprox.R internals accessible,
## so either load the package or source SaddlePointApprox.R first.


################################################################################
## Density series coefficient for PG(h, z)
################################################################################

a.coef.sp <- function(n, x, h, z)
{
  d.n   = (2 * n + h)
  l.out = h*log(2) - lgamma(h) + lgamma(n+h) - lgamma(n+1) + log(d.n) -
          0.5*log(2*pi*x^3) - 0.5*d.n^2/x - 0.5*z^2*x
  cosh(z)^h * exp(l.out)
}


################################################################################
## Additional density approximations in u- and x-space
################################################################################

log_sin_rt <- function(u)
{
  r   = sqrt(abs(u))
  out = log(sinh(r))
  out[u > 0] = log(sin(r[u > 0]))
  out
}

u.dens <- function(u, n=1, z=0)
{
  x = BayesLogit:::utox.mgf(u)
  0.5 * cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x)/u)^0.5 *
    exp(-n * BayesLogit:::log_cos_rt(u) - n*0.5*(u+z^2)*x)
}

u.dens.2 <- function(u, n=1, z=0)
{
  x = BayesLogit:::utox.mgf(u)
  0.5 * cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x)/u)^0.5 *
    exp(-n*(log_sin_rt(u) - log(sqrt(abs(u))) - log(x)) - n*0.5*(u+z^2)*x)
}

u.a1 <- function(u, n=1, z=0)
{
  x = 1 + 1/3*u + 3/15*u^2
  0.5 * cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x)/u)^0.5 *
    exp(-n * BayesLogit:::log_cos_rt(u) - n*0.5*(u+z^2)*x)
}

u.a2 <- function(u, n=1, z=0)
{
  x = 1 + 1/3*u + 3/15*u^2
  0.5 * cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x)/u)^0.5 *
    exp(-n*(log(pi/2) + log(x)) - n*0.5*(u+z^2)*x)
}

x.dens.1 <- function(u, n=1, z=0)
{
  x = BayesLogit:::utox.mgf(u)
  cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x)/u)^(-0.5) *
    exp(-n * BayesLogit:::log_cos_rt(u) - n*0.5*(u+z^2)*x)
}

x.dens.2 <- function(u, n=1, z=0)
{
  x = BayesLogit:::utox.mgf(u)
  cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x)/u)^(-0.5) *
    exp(-n*(log_sin_rt(u) - log(sqrt(abs(u))) - log(x)) - n*0.5*(u+z^2)*x)
}

x.a1 <- function(u, n=1, z=0)
{
  x = BayesLogit:::utox.mgf(u)
  cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x)/u)^(-0.5) *
    exp(-n*(log(1-u/6) - log(x)) - n*0.5*(u+z^2)*x)
}

approx1 <- function(x, h, z)
{
  (h/2/pi)^0.5 * x^(pi/2-1)
}


################################################################################
## Grid-based inversion of utox.mgf (u -> x map)
################################################################################

invert.x.1 <- function(x, xgrid, ugrid)
{
  out = list("u"=NA)
  if (x > max(xgrid)) return(out)
  dxgrid = x - xgrid
  idx    = which.min(abs(dxgrid))
  inc    = 0
  if (dxgrid[idx] < 0) inc =  1
  if (dxgrid[idx] > 0) inc = -1
  dx = xgrid[idx+inc] - xgrid[idx]
  du = ugrid[idx+inc] - ugrid[idx]
  ap = ugrid[idx]
  if (inc != 0) ap = du/dx * (x - xgrid[idx]) + ap
  out$u   = ap
  out$idx = idx
  out
}

invert.x <- function(x, xgrid, ugrid)
{
  ap = rep(0, length(x))
  for (i in seq_along(x)) ap[i] = invert.x.1(x[i], xgrid, ugrid)$u
  ap
}


################################################################################
## Plot 1: SP approximation vs. exact PG density
################################################################################

if (FALSE) {
  ## source("inst/dev/R/sp_analysis.R")
  z = 0.0
  n = 10.0

  umin  = -100
  umax  = (pi/2)^2 - 0.5
  du    = 0.1
  ugrid = seq(umin, umax, du)

  tgrid  = 0.5 * (ugrid + z^2)
  xgrid  = BayesLogit:::utox.mgf(ugrid)
  kgrid  = BayesLogit:::k.mgf(tgrid, z)
  term1  = (0.5*n/pi)^0.5 * (xgrid^2 + (1-xgrid)/ugrid)^(-0.5)
  term2  = exp(n * (kgrid - tgrid*xgrid))
  sp.approx = term1 * term2

  y1 = xgrid
  N  = length(xgrid)
  for (i in 1:N) {
    y1[i] = 0
    for (j in 0:200)
      y1[i] = y1[i] + (-1)^j * a.coef.sp(j, xgrid[i]*n, n, z) * n
  }

  m1 = jj.m1(n, z) / n
  m2 = jj.m2(n, z) / n^2
  y2 = dnorm(xgrid, m1, sqrt(m2 - m1^2))

  ymax = max(sp.approx, na.rm=TRUE)
  par(mfrow=c(1,2))
  plot(xgrid, y1, type="l", ylim=c(0,ymax),
       main="Density of JJ(b,z)", xlab="x", ylab="f(x|b,0)")
  lines(xgrid, sp.approx, col=4, lty=4)
  legend("topright", legend=c("J*","S.P."), col=c(1,4), lty=c(1,4))

  plot(xgrid/4*n, 4*y1/n, type="l", ylim=c(0,4*ymax/n),
       main=paste("Density of PG(", n, ")", sep=""), xlab="x", ylab="f")
  lines(xgrid/4*n, 4*sp.approx/n, col=4, lty=4)
  legend("topright", legend=c("PG","S.P."), col=c(1,4), lty=c(1,4))
}


################################################################################
## Plot 2: eta envelope construction
################################################################################

if (FALSE) {
  z = 0.0

  umin  = -10
  umax  = (pi/2)^2 - 0.5
  du    = 0.01
  ugrid = seq(umin, umax, du)
  xgrid = BayesLogit:::utox.mgf(ugrid)
  tgrid = 0.5 * (ugrid + z^2)
  kgrid = BayesLogit:::k.mgf(tgrid, z)

  deltaxgrid = sapply(xgrid, function(x) BayesLogit:::delta.val(x)$val)

  x.l = 0.75
  x.r = 4/3
  l.u = invert.x(x.l, xgrid, ugrid)
  u.r = invert.x(x.r, xgrid, ugrid)
  t.l = 0.5 * l.u
  t.r = 0.5 * u.r

  left.slope  = -t.l - 0.5/x.l^2
  right.slope = -t.r - 1/x.r
  l.int = BayesLogit:::k.mgf(t.l, z) - t.l*x.l - 0.5*(1-1/x.l)
  r.int = BayesLogit:::k.mgf(t.r, z) - t.r*x.r - log(x.r)

  left.line  = (xgrid - x.l) * left.slope  + l.int
  right.line = (xgrid - x.r) * right.slope + r.int

  pw.line = -deltaxgrid
  lc = which.min(abs(xgrid-1))
  pw.line[1:lc]        = left.line[1:lc]
  pw.line[lc:length(pw.line)] = right.line[lc:length(pw.line)]

  equigrid  = seq(0, max(xgrid), 0.1)
  zerogrid  = rep(0, length(equigrid))

  par(mfrow=c(1,2))
  plot(xgrid, (kgrid - tgrid*xgrid) - deltaxgrid, col=1, type="l",
       xlab="x", ylab="", main="eta envelope")
  lines(xgrid, -deltaxgrid, lty=2)
  lines(xgrid, kgrid - tgrid*xgrid, col=2, lty=3)
  lines(equigrid, zerogrid, col=2, lty=4)
  lines(xgrid, pw.line, col=3)
  legend("bottom",
         legend=c("eta-phi(1)","phi-phi(1)","-delta","envelope"),
         col=c(1,2,1,3), lty=c(1,3,2,1))

  plot(xgrid, kgrid - tgrid*xgrid, col=2, type="l", lty=1,
       xlab="x", ylab="", main="phi envelope")
  lines(xgrid, pw.line + deltaxgrid, col=3)
  lines(equigrid, zerogrid, col=2, lty=4)
  legend("bottom", legend=c("phi-phi(m)","phi envelope"), col=c(2,3), lty=c(3,1))
}


################################################################################
## Plot 3: u-space density approximation comparisons
################################################################################

if (FALSE) {
  z = 0.0; n = 10.0
  ugrid = seq(-100, (pi/2)^2 - 0.5, 0.1)
  xgrid = BayesLogit:::utox.mgf(ugrid)
  tgrid = 0.5 * (ugrid + z^2)
  kgrid = BayesLogit:::k.mgf(tgrid, z)

  udens1 = 0.5*cosh(z)^n * (0.5*n/pi)^0.5 * (xgrid^2 + (1-xgrid)/ugrid)^0.5 *
             exp(n*(kgrid - tgrid*xgrid))
  udens2 = u.dens(ugrid, n, z)
  udens3 = u.a1(ugrid, n, z)
  udens4 = u.a2(ugrid, n, z)

  plot(ugrid, udens1, type="l", main="u-space density", xlab="u", ylab="")
  lines(ugrid, udens2, col=2)
  lines(ugrid, udens3, col=3)
  lines(ugrid, udens4, col=4)
  legend("topright", legend=c("exact","u.dens","u.a1","u.a2"),
         col=1:4, lty=1)
}
