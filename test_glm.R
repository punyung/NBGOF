# try to test gene-wise estimation accuracy
rm=list(ls())
source("R/nb.gof.m.R")
library(doMC)

seed = 539
sim = 999
conf.env = 0.95

m = 1000 # 1000 genes
n = 6    # 6 samples
s = 1e6  # lib.sizes
offset = log(s)  # make sure offset for NBP functions should take the log

## simulate coefficients:
beta = matrix(0, m, 1)  # beta0 only
set.seed(seed)
beta[,1] = rnorm(m, 5, 2) - offset   # beta0 ~ normal (mean and a based on real RNA-Seq data)

## design matrix (single group only):
x = matrix(rep(1,n))

## specify mean levels:
mu = round(t(s * exp(x %*% t(beta))))
pi = mu/s


alpha1 = -0.2  # NB1.8
phi0 = 0.02
alpha0 = log(phi0)
phi.nbp = phi0 * pi^alpha1
range(phi.nbp)
cbind(mu[,1], phi.nbp[,1])

## add noise:
a = 0.1
set.seed(seed)
phi.noi.vec = phi.nbp[ ,1] * exp(matrix(runif(m, -a, a), nr=m, nc=1))
phi.noi = matrix(phi.noi.vec, nr=m, nc=n)
range(phi.noi)
cbind(mu[,1], phi.noi[,1])  # make sure phi's are in reasonable range


## generate NBP response with added noise:
set.seed(seed)
y = rnbinom(m * n, mu=mu, size=1/phi.noi)
dim(y) = dim(mu)
rownames(y) = paste("g", seq(1,m), sep="")
colnames(y) = paste("s", seq(1,n), sep="")
plot(mu, phi.noi, log="xy")

# evaluate dispersion model
ftgc.noip = nb.gof.m(counts=y, x=x, sim=sim, model="Tagwise-Common", method="CoxReid", seed=1, ncores=NULL)
