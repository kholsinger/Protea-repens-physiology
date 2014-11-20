require(R2jags)
require(mvtnorm)

rm(list=ls())

use.PCA <- TRUE

debug <- FALSE
## prior parameters for covariance matrix
##
gamma.rate.resid <- 1.0
gamma.shape.resid <- 1.0
gamma.rate.species <- 1.0
gamma.shape.species <-  1.0
beta.par <- 10
max.r <- 0.6
## prior precision on regression coefficients and year effect
##
tau <- 0.1

if (use.PCA) {
  model.file <- "traits-environment-pca.txt"
} else {
  model.file <- "traits-environment.txt"
}

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

if (use.PCA) {
  combined <- read.csv("traits-environment-pca.csv", na.strings=".", header=TRUE)
} else {
  combined <- read.csv("traits-environment.csv", na.strings=".", header=TRUE)
}

## prepare the data for JAGS
##
species <- as.numeric(combined$source_pop)
sla <- standardize(combined$SLA)
area <- standardize(combined$leaf_area)
sd <- standardize(combined$stomatal_density)
lwr <- standardize(combined$LWR)
spi <- standardize(combined$SPI)
if (use.PCA) {
  pca1 <- standardize(combined$Prin1_temp)
  pca2 <- standardize(combined$Prin2_dry)
  pca3 <- standardize(combined$Prin3_map)
} else {
  map <- standardize(combined$MAP)
  mat <- standardize(combined$MAT)
  ratio <- standardize(combined$rain_DecJanFeb)
}
## year read in as integer: convert to factor and back to numeric
## to make it (1,2) instead of (2013,2014)
##
year <- as.numeric(as.factor(combined$year))

## set up temporary data frame for cleaning and manipulation
##
if (use.PCA) {
  tmp <- data.frame(species=species,
                    sla=sla,
                    area=area,
                    lwr=lwr,
                    sd=sd,
                    spi=spi,
                    pca1=pca1,
                    pca2=pca2,
                    pca3=pca3,
                    year=year)
} else {
  tmp <- data.frame(species=species,
                    sla=sla,
                    area=area,
                    lwr=lwr,
                    sd=sd,
                    spi=spi,
                    map=map,
                    mat=mat,
                    ratio=ratio,
                    year=year)
}
## remove lines for which all response variables are missing
## if sla is missing, all response variables are
##
tmp <- subset(tmp, !is.na(sla), drop=TRUE)

## pull out remaining lines with all response variables present
## if sd is there, all response variables are
##
complete <- subset(tmp, !is.na(sd), drop=TRUE)

## pull out remaining lines with sd and spi missing
## if sd is missing, so is spi
##
incomplete <- subset(tmp, is.na(sd), drop=FALSE)

## re-extract vectors for JAGS
##
## complete data
##
species <- complete$species
sla <- complete$sla
area <- complete$area
sd <- complete$sd
lwr <- complete$lwr
spi <- complete$spi
if (use.PCA) {
  pca1 <- complete$pca1
  pca2 <- complete$pca2
  pca3 <- complete$pca3
} else {
  map <- complete$map
  mat <- complete$mat
  ratio <- complete$ratio
}
year <- complete$year
n.samp <- nrow(complete)
## missing sd and spi
##
species.inc <- incomplete$species
sla.inc <- incomplete$sla
area.inc <- incomplete$area
lwr.inc <- incomplete$lwr
if (use.PCA) {
  pca1.inc <- incomplete$pca1
  pca2.inc <- incomplete$pca2
  pca3.inc <- incomplete$pca3
} else {
  map.inc <- incomplete$map
  mat.inc <- incomplete$mat
  ratio.inc <- incomplete$ratio
}
year.inc <- incomplete$year
n.samp.inc <- nrow(incomplete)

n.species <- max(species, na.rm=TRUE)
n.dim <- 5
n.dim.inc <- 3
n.species.dim <- n.species*n.dim

## structured random effect matrix to reflect co-ancestry of species
##
Ginv <- diag(x=1.0, nrow=n.species, ncol=n.species)

## construct response matrix
##
## IMPORTANT NOTE: As currently written this code assumes that the columns with
## missing response data come immediately after those without missing data.
## So SLA, AREA, and LWR must precede SD and SPI
##
y <- as.matrix(data.frame(sla,
                          area,
                          lwr,
                          sd,
                          spi))
z <- as.matrix(data.frame(sla.inc,
                          area.inc,
                          lwr.inc))
if (use.PCA) {
  jags.data <- c("species",
                 "species.inc",
                 "y",
                 "z",
                 "pca1",
                 "pca2",
                 "pca3",
                 "year",
                 "pca1.inc",
                 "pca2.inc",
                 "pca3.inc",
                 "year.inc",
                 "n.samp",
                 "n.samp.inc",
                 "n.species",
                 "n.dim",
                 "n.dim.inc",
                 "n.species.dim",
                 "Ginv",
                 "tau",
                 "gamma.rate.resid",
                 "gamma.shape.resid",
                 "gamma.rate.species",
                 "gamma.shape.species",
                 "beta.par",
                 "max.r")
  jags.par <- c("beta.pca1",
                "beta.pca2",
                "beta.pca3",
                "rho.resid",
                "Sigma.resid",
                "rho.species",
                "Sigma.species",
                "yr")

} else {
  jags.data <- c("species",
                 "species.inc",
                 "y",
                 "z",
                 "map",
                 "mat",
                 "ratio",
                 "year",
                 "map.inc",
                 "mat.inc",
                 "ratio.inc",
                 "year.inc",
                 "n.samp",
                 "n.samp.inc",
                 "n.species",
                 "n.dim",
                 "n.dim.inc",
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
                "rho.resid",
                "Sigma.resid",
                "rho.species",
                "Sigma.species",
                "yr")
}

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

## print the results
##
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
options(opt.old)

if (!debug) {
  filename <- paste("results-",
                    gsub(":", "-",
                         gsub(" ", "-", Sys.time())),
                    ".Rsave",
                    sep="")
  save(fit, file=filename)
}
