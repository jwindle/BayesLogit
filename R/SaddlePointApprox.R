# (C) Nicholas Polson, James Scott, Jesse Windle, 2012-2019

# This file is part of BayesLogit.

# BayesLogit is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# BayesLogit.  If not, see <https://www.gnu.org/licenses/>.



## HERE I AM USING log laplace TRANSFORM.

## source("Ch.R")

a.coef.sp <- function(n, x, h, z)
{
  ## a_n(x,h).
  ## You could precompute lgamma(h), log(2).
  d.n = (2 * n + h)
  l.out = h * log(2) - lgamma(h) + lgamma(n+h) - lgamma(n+1) + log(d.n) -
    0.5 * log(2 * pi * x^3) - 0.5 * d.n^2 / x - 0.5 * z^2  * x;
  out = cosh(z)^h * exp(l.out)
  out
}

utox.laplace <- function(u)
{
  ## tanh(sqrt(u))/sqrt(u))
  out = u
  gt.idc = u >  1e-6
  lt.idc = u < -1e-6
  md.idc = u <= 1e-6 & u >= -1e-6

  r = sqrt(abs(u))
  gt.val = r[gt.idc]
  out[gt.idc] = tanh(gt.val)/gt.val
  lt.val = r[lt.idc]
  out[lt.idc] = tan(lt.val)/ lt.val
  md.val = r[md.idc]
  out[md.idc] = 1 - (1/3) * md.val^2 + (2/15) * md.val^4 - (17/315) * md.val^6
  out
}

utox.mgf <- function(u)
{
  utox.laplace(-1*u)
}

k.laplace <- function(t,z=0)
{
  s = 2*t + z^2;
  u = sqrt(abs(s));
  out = log(cosh(u))
  out[s<0] = log(cos(u[s<0]));
  ## out = ifelse(s >= 0, log(cosh(u)), log(cos(u)))
  out = log(cosh(z)) - out
  out
}

k1.laplace <- function(t,z=0)
{
  u = 2*t+z^2
  -1 * utox.laplace(u)^2
}

k2.laplace <- function(t,z=0)
{
  u = 2*t+z^2
  a1 = k2.laplace(t,z);
  a2 = (1-utox.laplace(u)) / u
  a1 - a2
}

k.mgf <- function(t,z=0)
{
  k.laplace(-1*t,z)
}

k1.mgf <- function(t, z=0)
{
  -1* k1.laplace(-1*t,z)
}

k2.mgf <- function(t,z=0)
{
  k2.laplace(-1*t,z)
}

################################################################################

log.cos.rt <- function(u)
{
  r   = sqrt(abs(u))
  out = log(cosh(r))
  out[u>0] = log(cos(r[u>0]));
  out
}

log.sin.rt <- function(u)
{
  r   = sqrt(abs(u))
  out = log(sinh(r))
  out[u>0] = log(sin(r[u>0]));
  out
}

u.dens <- function(u, n=1, z=0)
{
  x = utox.mgf(u)
  out = 0.5 * cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x) / u)^0.5 *
    exp(-n * log.cos.rt(u) - n * 0.5 * (u+z^2) * x);
  out
}

u.dens.2 <- function(u, n=1, z=0)
{
  x = utox.mgf(u)
  out = 0.5 * cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x) / u)^0.5 *
    exp(-n * (log.sin.rt(u) - log(sqrt(abs(u))) - log(x)) - n * 0.5 * (u+z^2) * x);
  out
}

u.a1 <- function(u, n=1, z=0)
{
  x = 1 + 1/3 * u + 3/15 * u^2
  out = 0.5 * cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x) / u)^0.5 *
    exp(-n * log.cos.rt(u) - n * 0.5 * (u+z^2) * x);
  out
}

u.a2 <- function(u, n=1, z=0)
{
  x = 1 + 1/3 * u + 3/15 * u^2
  out = 0.5 * cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x) / u)^0.5 *
    exp(-n * (log(pi/2) + log(x)) - n * 0.5 * (u+z^2) * x);
  out
}

x.dens.1 <- function(u, n=1, z=0)
{
  x = utox.mgf(u)
  out = cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x) / u)^(-0.5) *
    exp(-n * log.cos.rt(u) - n * 0.5 * (u+z^2) * x);
  out
}

x.dens.2 <- function(u, n=1, z=0)
{
  x = utox.mgf(u)
  out = cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x) / u)^(-0.5) *
    exp(-n * (log.sin.rt(u) - log(sqrt(abs(u))) - log(x))
        - n * 0.5 * (u+z^2) * x);
  out
}

x.a1 <- function(u, n=1, z=0)
{
  x = utox.mgf(u)
  out = cosh(z)^n * (0.5*n/pi)^0.5 * (x^2 + (1-x) / u)^(-0.5) *
    exp(-n * (log(1-u/6) - log(x))
        - n * 0.5 * (u+z^2) * x);
  out
}

delta.val <- function(x)
{
  ifelse(x < 1, 0.5 * (1-1/x), log(x))
}


invert.x.1 <- function(x, xgrid, ugrid)
{
  out = list("u"=NA)
  if (x > max(xgrid)) return(out)
  dxgrid = x - xgrid
  idx = which.min(abs(dxgrid))
  inc = 0
  if (dxgrid[idx] < 0) inc = 1
  if (dxgrid[idx] > 0) inc = -1
  dx = xgrid[idx+inc] - xgrid[idx]
  du = ugrid[idx+inc] - ugrid[idx]
  ap = ugrid[idx]
  if (inc!=0) ap = du / dx * (x - xgrid[idx]) + ap
  out$u   = ap
  out$idx = idx
  out
}

invert.x <- function(x, xgrid, ugrid)
{
  N = length(x)
  ap = rep(0, N)
  for (i in 1:N) {
    out = invert.x.1(x[i], xgrid, ugrid)
    ap[i] = out$u
  }
  ap
}

################################################################################

approx1 <- function(x, h, z)
{
  ## cosh(z)^h * (h/2/pi)^0.5 * (x^2 * (2 - x))^(-0.5) * cosh(1/x)^(-h) * exp(-0.5 * h / x) * exp(-0.5 * z^2 * h * x)
  (h/2/pi)^0.5 * x^(pi/2-1) ## * exp(-h*x*arccos(2/pi/x)^2) # *  cosh(z)^h * exp(-0.5 * z^2 * h * x)
}

################################################################################

if (FALSE)
{
  par(mfrow=c(1,3))
  
  umin = -10
  umax = (pi/2)^2 - 0.5
  du = 0.01
  ugrid = seq(umin, umax, du)
  xgrid = utox.mgf(ugrid)
  plot(xgrid, ugrid)

  tgrid = seq(-10,10, 0.1)
  kgrid = k(tgrid, 0.0)
  ## plot(tgrid, kgrid)

  z = 0.0
  n = 100.0
  ugrid = seq(umin, umax, du)
  ## ugrid = 2*tgrid - z^2
  tgrid = 0.5 * (ugrid + z^2)
  xgrid = utox.mgf(ugrid);
  kgrid = k.mgf(tgrid, z)
  term1 = (0.5 * n / pi)^(0.5) * (xgrid^2 + (1-xgrid) / ugrid)^(-0.5)
  term2 = exp(n * (kgrid - tgrid * xgrid));
  ygrid = term1 * term2

  plot(xgrid, ygrid)
  plot(xgrid, log(ygrid))

  udens =  (0.5 * n / pi)^(0.5) * (xgrid^2 + (1-xgrid) / ugrid)^(0.5) *
    exp(n * (kgrid - tgrid * xgrid));
  plot(ugrid, udens)
  plot(ugrid, log(udens))

}

if (FALSE)
{
  ## Preliminary: approximating using normal.

  ## source("SaddlePointApprox.R")
  ## png("pg-sp-dens.png", width=600, height=300); par(mfrow=c(1,2))
  
  z = 0.0
  n = 10.0

  umin = -100
  umax = (pi/2)^2 - 0.5
  du = 0.1
  ugrid = seq(umin, umax, du)
  
  tgrid = 0.5 * (ugrid + z^2)
  xgrid = utox.mgf(ugrid);
  kgrid = k.mgf(tgrid, z)
  term1 = (0.5 * n / pi)^(0.5) * (xgrid^2 + (1-xgrid) / ugrid)^(-0.5)
  term2 = exp(n * (kgrid - tgrid * xgrid));
  sp.approx = term1 * term2
  ## plot(xgrid, sp.approx)
  
  y1    = xgrid
  y2    = xgrid
  N     = length(xgrid)

  for (i in 1:N) {
    y1[i] = 0
    y2[i] = 0
    for (j in 0:200) {
      y1[i] = y1[i] + (-1)^j * a.coef.sp(j,xgrid[i] * n ,n, z) * n
    }
  }

  m1 = jj.m1(n, z) / n
  m2 = jj.m2(n, z) / n^2
  V  = m2 - m1^2;
  sv = sqrt(V)
  y2 = dnorm(xgrid, m1, sv);
  y3 =  dt((xgrid - m1) / sv, 6) / sv

  ## y4 = approx1(xgrid, n, z)
  
  ## png(filename="pg-dens.png", width=800, height=400)
  ## par(mfrow=c(1,2))
  ymax = max(c(sp.approx), na.rm=TRUE)
  plot(xgrid, y1, type="l", ylim=c(0,ymax), main="Density of JJ(b,z)", xlab="x", ylab="f(x|b,0)")
  ## lines(xgrid, y2, col=2, lty=2)
  ## lines(xgrid, y3, col=3, lty=3)
  ## lines(xgrid, y3, type="l", col=2, lty=2)
  lines(xgrid, sp.approx, col=4, lty=4)
  legend("topright", legend=c("J*", "S.P."), col=c(1,4), lty=c(1,4))

  plot(xgrid, log(y1), ylim=c(log(min(sp.approx, na.rm=TRUE)),log(ymax)), type="l")  
  lines(xgrid, log(sp.approx), col=4, lty=4)

  ## a0 = x.dens.2(ugrid, n, z)
  ## lines(xgrid, log(a0), col=2, lty=2)

  a1 = x.a1(ugrid, n, z)
  lines(xgrid, log(a1), col=2, lty=2)

  ## a2 = x.a2(ugrid, n, z)
  ## lines(xgrid, log(a2), col=3, lty=3)

  ## ----------
  
  ymax = max(c(sp.approx), na.rm=TRUE)
  plot(xgrid/4*n, 4*y1/n, type="l", ylim=c(0,4*ymax/n),
       main=paste("Density of PG(", n, ")", sep=""), xlab="x", ylab="f")
  ## lines(xgrid, y2, col=2, lty=2)
  ## lines(xgrid, y3, col=3, lty=3)
  ## lines(xgrid, y3, type="l", col=2, lty=2)
  lines(xgrid/4*n, 4*sp.approx/n, col=4, lty=4)
  legend("topright", legend=c("PG", "S.P."), col=c(1,4), lty=c(1,4))

  ##----------------------------------------------------------------------------

  ## png("eta-phi-envelope.png", width=800, height=400)
  par(mfrow=c(1,2))
  
  equigrid = seq(0, max(xgrid), 0.1)
  zerogrid = rep(0, length(equigrid))
  deltaxgrid = delta.val(xgrid)
  
  plot(xgrid, (kgrid - tgrid * xgrid) - deltaxgrid, col=1, type="l",
       xlab="x", ylab="phi(x)-phi(1) scale", main="eta envelope")
  lines(xgrid, -1*deltaxgrid, col=1, lty=2)

  lines(xgrid, (kgrid - tgrid * xgrid), col=2, lty=3)
  lines(equigrid, zerogrid, col=2, lty=4)

  x.l = .75
  x.r = 4/3
  l.u = invert.x(x.l, xgrid, ugrid)
  u.r = invert.x(x.r, xgrid, ugrid)
  t.l = 0.5 * l.u
  t.r = 0.5 * u.r

  left.slope = -t.l - 0.5 / x.l^2
  right.slope = -t.r - 1 / x.r

  l.int = k.mgf(t.l, z) - t.l * x.l - 0.5 * (1-1/x.l)
  r.int = k.mgf(t.r, z) - t.r * x.r - log(x.r)
  
  left.line = (xgrid - x.l) * left.slope + l.int
  right.line = (xgrid - x.r) * right.slope + r.int

  pw.line = -deltaxgrid
  ## left.cross = which.min(abs(left.line+deltaxgrid))
  left.cross = which.min(abs(xgrid-1))
  left.idx = 1:left.cross
  pw.line[left.idx] = left.line[1:left.cross]
  ## right.cross = which.min(abs(right.line+deltaxgrid))
  right.cross = which.min(abs(xgrid-1))
  right.idx = right.cross:length(pw.line)
  pw.line[right.idx] = right.line[right.idx]
 
  
  ## lines(xgrid, left.line, col=3)
  ## lines(xgrid, right.line, col=3)
  lines(xgrid, pw.line, col=3)

  legend("bottom", legend=c("eta(x)-phi(1)", "phi(x)-phi(1)", "-delta(x)", "eta envelope"),
         col=c(1,2,1,3), lty=c(1,3,2,1))
  
  plot(xgrid, (kgrid - tgrid * xgrid), col=2, type="l", lty=1,
       xlab="x", ylab="phi(x)-phi(1) scale", main="phi envelope")
  lines(xgrid, pw.line + deltaxgrid, col=3)
  lines(equigrid, zerogrid, col=2, lty=4)

  legend("bottom", legend=c("phi(x)-phi(m)","phi envelope"), col=c(2,3), lty=c(3,1))
  
  ##----------------------------------------------------------------------------

  ## Key here for approximation.  Use second order approx below.  I think second
  ## is better because we then have x^n for large x.
  plot(xgrid, log.cos.rt(ugrid) + 0.5)
  plot(xgrid, log.sin.rt(ugrid) - 0.5 * log(abs(ugrid)))

  b1 = log.cos.rt(ugrid) + 0.5 * ugrid * xgrid
  b2 = log.sin.rt(ugrid) -  0.5 * log(abs(ugrid)) - log(xgrid) + 0.5 * ugrid * xgrid

  plot(xgrid, b1)
  plot(xgrid, b2)
  lines(xgrid, -log(xgrid), col=2)

  ## z = 0.0
  b3 = log.sin.rt(ugrid) - 0.5 * log(abs(ugrid)) + 0.5 * ugrid * xgrid + 0.5 * z^2 * xgrid
  b4 = (xgrid^2 + (1-xgrid) / ugrid)
  
  plot(xgrid, b3)
  plot(xgrid, b4^-0.5)
  lines(xgrid, 1/xgrid, col=2)
  
  plot(xgrid, exp(b3))
  plot(xgrid, b3)
  
  a = 1.09 + z^2 / 2
  b = 0.175
  iggrid = (a * xgrid + b / xgrid) - 1.3
  lines(xgrid, iggrid, col=2)

  b5 = log.sin.rt(ugrid) - 0.5 * log(abs(ugrid))
  plot(xgrid, b5)
  plot(xgrid, ugrid * xgrid)

  plot(xgrid, b4 / xgrid^2)
  plot(xgrid, b4 / xgrid^(3))

  plot(xgrid, b4 / xgrid^(3))
  plot(xgrid, b4)
  lines(xgrid, xgrid^3, col=2)

  plot(xgrid, -(1-xgrid)/ugrid)
  plot(xgrid, xgrid^3 - b4)

  ## Plotting in u.
  plot(ugrid, b4 / xgrid^2)
  plot(ugrid, b4 / xgrid^3)
  
  ## k3
  b6 = 2 * xgrid * b4 - 0.5 * (1-xgrid) / tgrid^2 - 0.5 * b4 / tgrid
  b7 = b6 * xgrid

  plot(xgrid, b6)
  plot(xgrid, b7)

  b8 = -2 * b4 / xgrid^3 + b6 / xgrid^2 / b4
  b9 = -3 * b4 / xgrid^4 + b6 / xgrid^3 / b4

  plot(xgrid, b8)
  plot(xgrid, b9)

  b10 = -2 * b4^2 / xgrid^3 + b6 / xgrid^2 
  b11 = -3 * b4^2 / xgrid^4 + b6 / xgrid^3 

  plot(tgrid, b10)
  plot(tgrid, b11)

  b12 = cot(sqrt(2*tgrid))
  
  ##----------------------------------------------------------------------------
  
  ## plot(xgrid, sp.approx, col=4, type="l")
  udens1 =  0.5 * cosh(z)^n * (0.5 * n / pi)^(0.5) * (xgrid^2 + (1-xgrid) / ugrid)^(0.5) *
    exp(n * (kgrid - tgrid * xgrid));
  
  udens2 = u.dens(ugrid, n, z)
  udens3 = u.a1(ugrid, n, z)
  udens4 = u.a2(ugrid, n, z)
  udens5 = u.dens.2(ugrid, n, z)
  plot(ugrid, udens1, type="l")
  lines(ugrid, udens2, col=2)
  lines(ugrid, udens3, col=3)
  lines(ugrid, udens4, col=4)
  lines(ugrid, udens5, col=5)

  blah1 = - log.cos.rt(ugrid)
  blah2 = log(pi/2) + log(xgrid)

  ##----------------------------------------------------------------------------

  N = 10000
  c.n = ((1:N)-1/2)^2 * pi^2 / 2
  
  lw.fact <- function(u, n=1)
    {
      a = log(1 - outer(u, c.n, "/"));
      n * apply(a, 1, sum)
    }

  lwf.1 = lw.fact(ugrid, n)
  lwf.2 = n * log.cos.rt(ugrid)

  plot(ugrid, lwf.1)
  plot(ugrid, lwf.2)

  calc.psi <- function(c.n)
    {
      N = length(c.n)
      psi.0 = psi.1 = psi.2 = psi.3 = rep(0, N)
      psi.alt = psi.0
      for (i in 1:N) {
        a.i = 1-c.n[i]/c.n[-i]
        b.i = 1 / (c.n[i] * a.i)
        psi.0[i] = log(c.n[i]) - sum(log(abs(a.i)))
        psi.1[i] = sum(b.i)
        psi.2[i] = sum(b.i^2)
        psi.3[i] = sum(b.i^3)
        psi.alt[i] = c.n[i] * prod(1/a.i)
      }

      phi.0 = exp(psi.0)
      phi.1 = phi.0 * psi.1
      phi.2 = phi.0 * psi.1^2 + phi.0 * psi.2
      phi.3 = phi.0 * psi.1^3 + 3 * phi.0 * psi.1 * psi.2 + phi.0 * psi.3

      out <- list(psi.0=psi.0, psi.1=psi.1, psi.2=psi.2, psi.3=psi.3, psi.alt=psi.alt,
                  phi.0=phi.0, phi.1=phi.1, phi.2=phi.2, phi.3=phi.3)
      out
    }

  the.psi = calc.psi(c.n)
  plot(c.n, exp(the.psi$psi.0))
  
}
