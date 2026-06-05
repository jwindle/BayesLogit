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


################################################################################
## CGF utilities
##
## Convention: k.laplace(t, z) = log E[e^{-tX}] for X ~ PG(1, z).
##             k.mgf   (t, z) = log E[e^{ tX}] for X ~ PG(1, z).
## utox.laplace(u): solves the saddlepoint equation x = tanh(sqrt(u))/sqrt(u)
##                 (Laplace side: u > 0 => x in (0,1)).
################################################################################

utox.laplace <- function(u)
{
  out = u
  r   = sqrt(abs(u))
  pos = u >  1e-6
  neg = u < -1e-6
  mid = !pos & !neg
  out[pos] = tanh(r[pos]) / r[pos]
  out[neg] = tan(r[neg])  / r[neg]
  out[mid] = 1 - (1/3)*r[mid]^2 + (2/15)*r[mid]^4 - (17/315)*r[mid]^6
  out
}

utox.mgf <- function(u) utox.laplace(-u)

## log(cosh(sqrt(|u|))) for u <= 0; log(cos(sqrt(u))) for u > 0.
log_cos_rt <- function(u)
{
  r   = sqrt(abs(u))
  out = log(cosh(r))
  out[u > 0] = log(cos(r[u > 0]))
  out
}

k.laplace <- function(t, z=0)
{
  s   = 2*t + z^2
  u   = sqrt(abs(s))
  out = log(cosh(u))
  out[s < 0] = log(cos(u[s < 0]))
  log(cosh(z)) - out
}

## K'(t) = -utox.laplace(2t + z^2)
k1.laplace <- function(t, z=0)
{
  -utox.laplace(2*t + z^2)
}

## K''(t) = x^2 - (1-x)/s  where s = 2t+z^2, x = utox.laplace(s)
k2.laplace <- function(t, z=0)
{
  s = 2*t + z^2
  x = utox.laplace(s)
  ifelse(abs(s) >= 1e-6,
         x^2 - (1-x)/s,
         x^2 - 1/3 + (2/15)*s)
}

k.mgf  <- function(t, z=0)  k.laplace(-t, z)
k1.mgf <- function(t, z=0) -k1.laplace(-t, z)
k2.mgf <- function(t, z=0)  k2.laplace(-t, z)


################################################################################
## Saddlepoint sampler helpers
################################################################################

## Scalar inverse of utox.laplace: find s s.t. utox.laplace(s) = x.
##   x in (0,1]: s > 0 (tanh branch).  x > 1: s < 0 (tan branch).
## Mirrors C++ v_eval / InvertY, where the sign convention is v_C++ = -s_R.
v.eval <- function(x)
{
  if (abs(x - 1) < 1e-10) return(0.0)
  if (x < 1)
    uniroot(function(s) utox.laplace(s) - x, lower=1e-10,          upper=1000)$root
  else
    uniroot(function(s) utox.laplace(s) - x, lower=-(pi/2)^2+1e-6, upper=-1e-10)$root
}

## Piecewise delta function with midpoint mid (generalises the mid=1 case).
## Returns list(val, der) matching C++ delta_func.
delta.val <- function(x, mid=1)
{
  if (x >= mid)
    list(val = log(x) - log(mid), der = 1/x)
  else
    list(val = 0.5*(1 - 1/x) - 0.5*(1 - 1/mid), der = 0.5/x^2)
}

## Tangent line to eta = phi - delta at x.
## z is the internally halved value (z_internal = abs(z_original)/2).
## Returns list(slope, icept).
tangent.to.eta <- function(x, z, mid)
{
  s       = v.eval(x)
  t       = 0.5 * (z^2 - s)
  phi.val = log(cosh(z)) - log_cos_rt(-s) - t * x
  phi.der = -t

  d       = delta.val(x, mid)
  eta.der = phi.der - d$der
  eta.val = phi.val - d$val

  list(slope = eta.der,
       icept  = eta.val - eta.der * x)
}

## Saddlepoint density approximation. z is the internally halved value.
sp.approx.R <- function(x, n, z)
{
  s   = v.eval(x)
  t   = 0.5 * (z^2 - s)
  phi = log(cosh(z)) - log_cos_rt(-s) - t * x
  K2  = if (abs(s) >= 1e-6) x^2 - (1-x)/s else x^2 - 1/3 + (2/15)*s
  exp(0.5*log(0.5*n/pi) - 0.5*log(K2) + n*phi)
}

## Left-truncated Gamma(shape, rate) on (trunc, Inf) via inverse-CDF.
ltgamma.sp <- function(shape, rate, trunc)
{
  p = pgamma(trunc, shape=shape, rate=rate)
  qgamma(p + runif(1) * (1 - p), shape=shape, rate=rate)
}

## Right-truncated inverse-chi^2(lambda) on (0, trunc):
##   X = 1/Y,  Y ~ chi^2(lambda) on (1/trunc, Inf).
rtinvchi2.sp <- function(lambda, trunc)
{
  1 / ltgamma.sp(lambda/2, rate=0.5, trunc=1/trunc)
}

## Truncated IG(mu, lambda) on (0, trunc). Mirrors C++ PolyaGammaApproxSP::rtigauss.
rtigauss.sp <- function(mu, lambda, trunc)
{
  X = trunc + 1
  if (mu > trunc) {
    alpha = 0.0
    while (runif(1) > alpha) {
      X     = rtinvchi2.sp(lambda, trunc)
      alpha = exp(-0.5 * lambda / mu^2 * X)
    }
  } else {
    while (X > trunc) X = rigauss(mu, lambda)
  }
  X
}


################################################################################
## R saddlepoint sampler for PG(h, z) 
################################################################################

rpg.sp.1 <- function(h, z, max.iter=200)
{
  z = 0.5 * abs(z)

  xl  = utox.laplace(z^2)   # mode of phi
  md  = xl * 1.1             # midpoint
  xr  = xl * 1.2             # right point

  # K2 at midpoint
  s.md = v.eval(md)
  K2md = if (abs(s.md) >= 1e-6) md^2 - (1-md)/s.md else md^2 - 1/3 + (2/15)*s.md

  m2 = md^2
  al = m2 * md / K2md
  ar = m2      / K2md

  # Tangent lines
  ll = tangent.to.eta(xl, z, md)
  lr = tangent.to.eta(xr, z, md)
  rl = -ll$slope;  il = ll$icept
  rr = -lr$slope;  ir = lr$icept

  lcn   = 0.5 * log(0.5 * h / pi)
  rt2rl = sqrt(2 * rl)

  # Proposal weights
  wl = exp(0.5*log(al) - h*rt2rl + h*il + 0.5*h/md) *
       pigauss(md, 1/rt2rl, h)
  wr = exp(0.5*log(ar) + lcn +
           (-h*log(h*rr) + h*ir - h*log(md) + lgamma(h))) *
       pgamma(md, shape=h, rate=h*rr, lower.tail=FALSE)

  pl = wl / (wl + wr)

  go   = TRUE
  iter = 0

  while (go && iter < max.iter) {
    iter = iter + 1

    if (runif(1) < pl) {
      X      = rtigauss.sp(1/rt2rl, h, md)
      phi.ev = h*(il - rl*X) + 0.5*h*((1 - 1/X) - (1 - 1/md))
      F      = exp(0.5*log(al) + lcn - 1.5*log(X) + phi.ev)
    } else {
      X      = ltgamma.sp(h, h*rr, md)
      phi.ev = h*(ir - rr*X) + h*(log(X) - log(md))
      F      = exp(0.5*log(ar) + lcn + phi.ev) / X
    }

    if (F * runif(1) < sp.approx.R(X, h, z)) go = FALSE
  }

  h * 0.25 * X
}

rpg.sp.R <- function(num=1, h=1, z=0.0)
{
  if (any(h < 1)) stop("rpg.sp.R: h must be >= 1.")

  if (length(h) != num) h = array(h, num)
  if (length(z) != num) z = array(z, num)

  x = rep(0.0, num)
  for (i in 1:num) x[i] = rpg.sp.1(h[i], z[i])
  x
}
