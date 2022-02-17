library(TMB)


# I used the guide here as a starting point:
# https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecm.1470&file=ecm1470-sup-0001-AppendixS1.pdf

nT = 200

x = matrix(0, nrow = 2, ncol = nT + 1)

# Standard deviation of the process variation
sdp <- 0.1
# Set the seed, so we can reproduce the results
set.seed(553)
# For-loop that simulates the state through time,
for(t in 1:nT){
  # This is the process equation
  x[, t+1] <- x[, t] + rnorm(nrow(x), 0, sdp)
  # Note that this index is shifted compared to equation in text,
  # because we assume the first value to be at time 0
}

plot(0:nT, x[1,],
     pch = 19, cex = 0.7, col="red", ty = "o",
     xlab = "t", ylab = expression(x[t]), las = 1)
points(0:nT, x[2,],
     pch = 19, cex = 0.7, col="green", ty = "o",
     xlab = "t", ylab = expression(x[t]), las = 1)

# Create a vector that will keep track of the observations
y <- matrix(ncol = nT, nrow = 5)
# Standard deviation of the observation error
sdo <- 0.1

y[1,] <- x[1,2:(nT+1)] + rnorm(nT, 0, sdo)
y[2,] <- x[1,2:(nT+1)] + rnorm(nT, 0, sdo)
y[3,] <- .5*x[1,2:(nT+1)] + .5*x[2,2:(nT+1)] + rnorm(nT, 0, sdo)
y[4,] <- x[2,2:(nT+1)] + rnorm(nT, 0, sdo)
y[5,] <- x[2,2:(nT+1)] + rnorm(nT, 0, sdo)


compile("simpleDFA.cpp")
dyn.load(dynlib("simpleDFA"))

dataTmb <- list(y =y)

nfac = 2 # Number of factors
ny = nrow(y) # Number of time series

Zinit = matrix(rnorm(ny * nfac), ncol = nfac)
# Set constrained elements to zero
constrInd = rep(1:nfac, each = ny) > rep(1:ny,  nfac)
Zinit[constrInd] = 0
Zinit

tmbPar =  list(logSdO = 0, Z = Zinit,
                 x=matrix(c(rep(0, nfac), rnorm(nfac * nT)), ncol = nT+1, nrow = nfac))

# Set up parameter constraints. Elements set to NA will be fixed and not estimated.
Zmap = matrix(ncol = nfac, nrow = ny)
Zmap[constrInd] = NA
Zmap[!constrInd] = 1:sum(!constrInd)
xmap = matrix(ncol = nT + 1, nrow = nfac)
xmap[,1] = NA
xmap[(nfac + 1) : length(tmbPar$x)] = 1:(length(tmbPar$x) - nfac)
tmbMap = list(Z = as.factor(Zmap), x = as.factor(xmap))


tmbObj = MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= "x", DLL= "simpleDFA")

names(tmbObj)

tmbOpt = nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))

tmbOpt$message

sdRep <- summary(sdreport(tmbObj))

sdRep[grepl('Z|sdo', rownames(sdRep)),]


# Get point estimates
#x_hat <- matrix(sdRep[rownames(sdRep)=="x",1], nrow = nfac)

x_hat = (tmbObj$env$parList)()$x

Z_hat = (tmbObj$env$parList)(par=tmbOpt$par)$Z

matplot(t(y), pch =20)

matpoints(t(Z_hat %*% x_hat), type = 'l', lwd = 3)

matplot(t(x_hat), type = 'l')

# Compute AIC

# Only works if obj has been already optimized
AIC.tmb = function(obj, tol = 0.01) {
  # Simple convergence check
  stopifnot(max(abs(obj$gr(obj$env$last.par.best[obj$env$lfixed()]))) < tol)
  as.numeric(2 * obj$env$value.best + 2*sum(obj$env$lfixed()))
}

AIC.tmb(tmbObj)

#aic = 2 * as.numeric(tmbObj$fn(tmbOpt$par)) + 2 * length(tmbOpt$par)
#aic


# Plot rotated trends, see https://atsa-es.github.io/atsa-labs/sec-dfa-rotating-loadings.html
matplot(t(solve(varimax(Z_hat)$rotmat) %*% x_hat), type = 'l')

Z_hat %*% varimax(Z_hat)$rotmat
