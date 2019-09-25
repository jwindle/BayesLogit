

# "Shape" parameter and method used:
# Devroye: 1 or 2
# Alt: 1-13, but not 1 or 2
# SP: 13-170
# Normal: 170 and above

library("BayesLogit")

rpg.sp.R(1, 100, 12)

rpg.sp(1, 100, 12)


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


nsamp = 2000

out = list()

bseq = c(1, 3, 12, 36, 100)
zseq = c(0., 1., 2., 4., 12.)

lenb = length(bseq)
lenz = length(zseq)

# zseq = 0.0

count = 0
for (j in 1:lenz) {
  for (i in 1:lenb) {
    count = count + 1
    
    b = bseq[i]
    z = zseq[j]

    cat("Working on b =", b, ", z =", z, "\n")

    samp = list()

    samp[["gamma.R"]]   = rpg.gamma.R(nsamp, b, z, trunc=1000)

    samp[["devroye.R"]] = rpg.devroye.R(nsamp, b, z)

    ## samp[["alt.R"]]     = rpg.alt.R(nsamp, b, z)

    samp[["sp.R"]]      = rpg.sp.R(nsamp, b, z)

    samp[["gamma.C"]]   = rpg.gamma(nsamp, b, z, trunc=1000)
    
    samp[["devroye.C"]] = rpg.devroye(nsamp, b, z)
    
    ## samp[["alt.C"]]     = rpg.alt(nsamp, b, z)
    
    samp[["sp.C"]] = rpg.sp(nsamp, b, z)

    samp[["rpg.C"]] = rpg(nsamp, b, z)

    m1 = pg.m1(b,z)
    m2 = pg.m2(b,z)

    m1.samp = sapply(samp, mean)
    m2.samp = sapply(samp, function(x){mean(x^2)})

    out[[count]] = list()
    out[[count]]$param = c(b=b, z=z)
    out[[count]]$moments = c(m1=m1, m2=m2)
    out[[count]]$moments.samp = data.frame(m1=m1.samp, m2=m2.samp)
    out[[count]]$samp  = samp
    out[[count]]$error = data.frame(m1=m1-m1.samp, m2=m2-m2.samp)
    out[[count]]$rel.error = data.frame(m1=(m1-m1.samp)/m1, m2=(m2-m2.samp)/m2)
    
  }
}


# Compare 
for (i in 1:(lenz*lenb)) {
    print(out[[i]]$param)
    print(out[[i]]$moments)
    ## print(out[[i]]$error)
    print(out[[i]]$rel.error)
}


# Plot qq plots
par(mfrow=c(5,6))
par(mar=c(2, 2, 1, 0))

j = 1
for (i in 1:lenb) {
    k = (j-1)*lenb + i
    with(out[[k]]$samp, qqplot(gamma.R, devroye.R))
    with(out[[k]], title(main=paste("devroye.R vs. gamma.R", "b:", param[1], "z:", param[2]), cex.main=0.5))
    with(out[[k]]$samp, qqplot(gamma.R, sp.R))
    with(out[[k]], title(main=paste("sp.R vs. gamma.R", "b:", param[1], "z:", param[2]), cex.main=0.5))
    with(out[[k]]$samp, qqplot(gamma.R, gamma.C))
    with(out[[k]], title(main=paste("gamma.C vs. gamma.R", "b:", param[1], "z:", param[2]), cex.main=0.5))
    with(out[[k]]$samp, qqplot(gamma.R, devroye.C))
    with(out[[k]], title(main=paste("devroye.C vs. gamma.R", "b:", param[1], "z:", param[2]), cex.main=0.5))
    with(out[[k]]$samp, qqplot(gamma.R, sp.C))
    with(out[[k]], title(main=paste("sp.C vs. gamma.R", "b:", param[1], "z:", param[2]), cex.main=0.5))
    with(out[[k]]$samp, qqplot(gamma.R, rpg.C))
    with(out[[k]], title(main=paste("rpg.C vs. gamma.R", "b:", param[1], "z:", param[2]), cex.main=0.5))
}

