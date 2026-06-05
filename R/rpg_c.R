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
                               ## POLYAGAMMA ##
################################################################################

## Draw PG(h, z)
##------------------------------------------------------------------------------
rpg.gamma <- function(num=1, h=1, z=0.0, trunc=200)
{
    ## Check Parameters.
    if (sum(h<0)!=0) {
        stop("rpg.gamma: h must be greater than zero.");
    }
    if (trunc < 1) {
        stop("rpg.gamma: xtrunc must be > 0.");
    }

    x = rep(0, num);

    if (length(h) != num) { h = array(h, num); }
    if (length(z) != num) { z = array(z, num); }

    OUT = .C("rpg_gamma", x, h, z, as.integer(num), as.integer(trunc), PACKAGE="BayesLogit");

    OUT[[1]]
}

## Draw PG(n, z) where n is a natural number.
##------------------------------------------------------------------------------
rpg.devroye <- function(num=1, h=1, z=0.0)
{
    n = h
    
    ## Check Parameters.
    if (any(n<0)) {
      stop("rpg.devroye: h must be greater than zero.");
    }

    x = rep(0, num);

    if (length(n) != num) { n = array(n, num); }
    if (length(z) != num) { z = array(z, num); }

    OUT = .C("rpg_devroye", x, as.integer(n), z, as.integer(num), PACKAGE="BayesLogit");

    OUT[[1]]
}

## Draw PG(h, z) where h is \geq 1.
##------------------------------------------------------------------------------
rpg.alt <- function(num=1, h=1, z=0.0)
{
    ## Check Parameters.
    if (any(h<1)) {
      stop("rpg.alt: h must be >= 1.");
    }

    x = rep(0, num);

    if (length(h) != num) { h = array(h, num); }
    if (length(z) != num) { z = array(z, num); }

    OUT = .C("rpg_alt", x, h, z, as.integer(num), PACKAGE="BayesLogit");

    OUT[[1]]
}


## Draw PG(h, z) using SP approx where h is \geq 1.
##------------------------------------------------------------------------------
rpg.sp <- function(num=1, h=1, z=0.0)
{
    ## Check Parameters.
    if (any(h<1)) {
      stop("rpg.sp: h must be >= 1.");
    }

    x = rep(0, num);
    iter = rep(0, num);

    if (length(h) != num) { h = array(h, num); }
    if (length(z) != num) { z = array(z, num); }

    ## Faster if we do not track iter.
    OUT = .C("rpg_sp", x, h, z, as.integer(num), as.integer(iter), PACKAGE="BayesLogit");

    ## ## Tracking total iterations
    ## out = list() 
    ## if (!track.iter)
    ##   out = OUT[[1]]
    ## else
    ##   out = list(samp=OUT[[1]], iter=OUT[[5]])
    ## out

    OUT[[1]]
}

## Draw PG(n, z)
##------------------------------------------------------------------------------
rpg <- function(num=1, h=1, z=0.0)
{
    ## Check Parameters.
    if (any(h<=0)) {
      stop("rpg: h must be > 0.");
    }

    x = rep(0, num);

    if (length(h) != num) { h = array(h, num); }
    if (length(z) != num) { z = array(z, num); }

    ## Faster if we do not track iter.
    OUT = .C("rpg_hybrid", x, h, z, as.integer(num), PACKAGE="BayesLogit");

    OUT[[1]]
}

