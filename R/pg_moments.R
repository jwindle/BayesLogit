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
## Moments of the J*(b, z) and PG(b, z) distributions
##
## PG(b, z) = (1/4) * J*(b, z/2), so pg.m1/m2 scale jj.m1/m2 accordingly.
## Taylor series branches handle the z -> 0 limit stably.
################################################################################

jj.m1 <- function(b, z)
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
      b * ((-1/3) + (2/15) * z - (17/315) * z^3)
}

## First moment: E[X] for X ~ PG(b, z)
pg.m1 <- function(b, z)
{
  jj.m1(b, z/2) / 4
}

## Second moment: E[X^2] for X ~ PG(b, z)
pg.m2 <- function(b, z)
{
  jj.m2(b, z/2) / 16
}

## Variance: Var[X] for X ~ PG(b, z)
pg.var <- function(b, z)
{
  pg.m2(b, z) - pg.m1(b, z)^2
}
