require(R2jags)
require(mvtnorm)

rm(list=ls())

debug <- FALSE
plot <- TRUE
print <- TRUE
report.DIC <- TRUE
full.data <- FALSE

gamma.rate.resid <- 1.0
gamma.shape.resid <- 1.0
gamma.rate.species <- 1.0
gamma.shape.species <-  1.0
beta.par <- 6
max.r <- 0.4

model.file="traits-environment.txt"

if (debug) {
  n.chains <- 1
  n.burnin <- 500
  n.iter <- 1000
  n.thin <- 1
  ## to allow replication of results across runs in JAGS
  ##
  set.seed(1)
} else {
  n.chains <- 5
  n.burnin <- 5000
  n.iter <- 25000
  n.thin <- 25
}

standardize <- function(x) {
  if (is.numeric(x)) {
    y <- (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
  } else {
    y <- x
  }
  y
}

get.mean.vector <- function(x) {
  x.mean <- apply(x, 1, mean, na.rm=TRUE)
}

drop.levels <- function(dat) {
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}

likelihood <- function(y, mu, Sigma) {
  llike <- 0.0
  for (i in 1:nrow(y)) {
    llike <- llike + dmvnorm(y[i,], mu[i,], Sigma[,], log=TRUE)
  }
  llike
}

Dbar <- function(y, mu, Sigma) {
  dbar <- numeric(0)
  n.rep <- dim(mu)[1]
  for (i in 1:n.rep) {
    if ((i %% 10) == 0) {
      cat(".", sep="")
      flush.console()
    }
    if ((i %% 500) == 0) {
      cat(i, "\n")
      flush.console()
    }
    dbar[i] <- -2.0*likelihood(y, mu[i,,], Sigma[i,,])
  }
  mean(dbar)
}

Dhat <- function(y, mu.mean, Sigma.mean) {
  -2.0*likelihood(y, mu.mean, Sigma.mean)
}

## read in appropriate data set
##
if (full.data) {
  combined <- read.csv("traits-environment.csv", na.strings=".", header=TRUE)
} else {
  combined <- read.csv("combined.csv", na.strings=".", header=TRUE)
}

## prepare the data for JAGS
##
if (full.data) {
  species <- as.numeric(combined$source_pop)
  sla <- standardize(combined$SLA)
  area <- standardize(combined$leaf_area)
  sd <- standardize(combined$stomatal_density)
  lwr <- standardize(combined$LWR)
  spi <- standardize(combined$SPI)
  map <- standardize(combined$MAP)
  mat <- standardize(combined$MAT)
  ratio <- standardize(combined$rain_DecJanFeb)
  ## set up temporary data frame for cleaning and manipulation
  ##
  tmp <- data.frame(species=species,
                    sla=sla,
                    area=area,
                    sd=sd,
                    lwr=lwr,
                    spi=spi,
                    map=map,
                    mat=mat,
                    ratio=ratio)
  ## remove lines for which all response variables are missing
  ## if sla is missing, all response variables are
  ##
  tmp <- subset(tmp, !is.na(sla), drop=TRUE)
  ## pull out lines with all response variables present
  ## if sd is there, all response variables are
  ##
  complete <- subset(tmp, !is.na(sd), drop=TRUE)
  ## pull out lines with sd and spi missing
  ## if sd is missing, so is spi
  ##
  incomplete <- subset(tmp, is.na(sd), drop=FALSE)
} else {
  species <- as.numeric(combined$species)
  sla <- standardize(combined$SLA)
  area <- standardize(combined$Leaf_area)
  sd <- standardize(combined$SD)
  lwr <- standardize(combined$Lwratio)
  spi <- standardize(combined$SPI)
  map <- standardize(combined$MAP)
  mat <- standardize(combined$MAT)
  ratio <- standardize(combined$ratio)
  n.samp <- nrow(combined)
}

n.species <- max(species, na.rm=TRUE)
n.dim <- 5
n.species.dim <- n.species*n.dim

## structured random effect matrix to reflect co-ancestry of species
##
Ginv <- diag(x=1.0, nrow=n.species, ncol=n.species)

## construct response matrix
##
y <- as.matrix(data.frame(sla,
                          area,
                          sd,
                          lwr,
                          spi))

## parameters for Wishart prior
##
## nrow = ncol = # of parameters, i.e., ncol(y) == 5
## nu <- nrow + 2 makes it as vague as possible
## Note: nu > nrow + 1 required for distribution to be
##       non-degenerate
##
Omega <- diag(x=1.0, nrow=ncol(y), ncol=ncol(y))
nu <- nrow(Omega) + 2

## prior precision on regression coefficients
##
tau <- 0.1


jags.data <- c("species",
               "y",
               "map",
               "mat",
               "ratio",
               "n.samp",
               "n.species",
               "n.dim",
               "n.species.dim",
               "Ginv",
               "tau",
               "gamma.rate.resid",
               "gamma.shape.resid",
               "gamma.rate.species",
               "gamma.shape.species",
               "beta.par",
               "max.r")
jags.par <- c("beta.map",
               "beta.mat",
               "beta.ratio",
               "mu",
               "rho.resid",
               "Sigma.resid",
               "rho.species",
               "Sigma.species")

fit <- jags(data=jags.data,
            inits=NULL,
            parameters=jags.par,
            model.file=model.file,
            n.chains=n.chains,
            n.burnin=n.burnin,
            n.iter=n.iter,
            n.thin=n.thin,
            DIC=TRUE,
            working.directory=".")

mu.mean <- fit$BUGSoutput$mean$mu
Sigma.mean <- fit$BUGSoutput$mean$Sigma.resid
mu <- fit$BUGSoutput$sims.list$mu
Sigma <- fit$BUGSoutput$sims.list$Sigma.resid
summary<-fit$BUGSoutput$summary

if (print) {
  opt.old <- options(width=120)
  filename <- paste("results-",
                    gsub(":", "-",
                         gsub(" ", "-", Sys.time())),
                    ".txt",
                    sep="")
  if (!debug) {
    sink(filename, split=TRUE)
  }

  cat("beta.par:            ", beta.par, "\n")
  cat("max.r:               ", max.r, "\n\n")
  print(fit, digits.summary=3)
  if (!debug) {
    sink()
  }
  if (report.DIC) {
    dbar <- Dbar(y, mu, Sigma)
    dhat <- Dhat(y, mu.mean, Sigma.mean)
    pD <- dbar - dhat
    DIC <- dbar + pD
    if (!debug) {
      sink(filename, append=TRUE, split=TRUE)
    }
    cat("\n",
        "Dbar: ", dbar, "\n",
        "Dhat: ", dhat, "\n",
        "pD:   ", pD, "\n",
        "DIC:  ", DIC, "\n")
    if (!debug) {
      sink()
    }
  }
  options(opt.old)
}


if (plot) {
  dev.new()
  old.par <- par(mfrow=c(2,3))
  for (i in 1:5) {
    x.plot <- y[,i]
    y.plot <- mu.mean[,i]
    plot(x.plot, y.plot, xlab="Observed", ylab="Predicted",
         pch=16, cex=0.75,
         main=c("LMA", "Area", "LWR", "SPI", "SD")[i])
  }
  par(old.par)
  dev.new()
  old.par <- par(mfrow=c(2,3))
  for (i in 1:5) {
    x.plot <- mu.mean[,i]
    y.plot <- y[,i]-mu.mean[,i]
    plot(x.plot, y.plot, xlab="Predicted", ylab="Residual",
         pch=16, cex=0.75,
         main=c("LMA", "Area", "LWR", "SPI", "SD")[i])
  }
  par(old.par)
}

filename <- paste("results-",
                  gsub(":", "-",
                       gsub(" ", "-", Sys.time())),
                  ".Rsave",
                  sep="")
save(fit, file=filename)
