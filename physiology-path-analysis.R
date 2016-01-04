require(R2jags)

rm(list=ls())

debug <- FALSE

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

fit <- function(dat, label, model) {
  model.file <- paste("physiology-path-analysis-model-", model, ".txt", sep="")

  ## set up multi-response vectors
  ##
  leaf <- as.matrix(data.frame(dat$sla,
                               dat$lwr,
                               dat$lfarea,
                               dat$spi,
                               dat$sdi))
  if (model == 1) {
    phys <- as.matrix(data.frame(dat$lftemp,
                                 dat$fluor,
                                 dat$trans,
                                 dat$cond))
  } else if (model == 2) {
    phys <- as.matrix(data.frame(dat$lftemp,
                                 dat$fluor,
                                 dat$trans))
  } else if (model == 3) {
    phys <- as.matrix(data.frame(dat$lftemp,
                                 dat$fluor,
                                 dat$trans))
  } else if (model == 4) {
    phys <- as.matrix(data.frame(dat$lftemp,
                                 dat$fluor,
                                 dat$cond,
                                 dat$trans))
  } else if (model == 5) {
    phys <- as.matrix(data.frame(dat$lftemp,
                                 dat$fluor,
                                 dat$cond,
                                 dat$trans))
  }
  ## set up individual covariates
  ##
  photo <- dat$photo
  wue <- dat$wue
  cond <- dat$cond
  humid <- dat$humid
  temp <- dat$temp
  plant <- dat$plant
  n.plant <- max(plant)
  n.samp <- nrow(leaf)
  ## parameters for Wishart prior
  ##
  ## nrow = ncol = # of parameters, i.e., ncol(y) == 4
  ## nu <- nrow + 2 makes it as vague as possible
  ## Note: nu > nrow + 1 required for distribution to be
  ##       non-degenerate
  ##
  Omega.lf <- diag(x=1.0, nrow=ncol(leaf), ncol=ncol(leaf))
  nu.lf <- nrow(Omega.lf) + 2
  Omega.ph <- diag(x=1.0, nrow=ncol(phys), ncol=ncol(phys))
  nu.ph <- nrow(Omega.ph) + 2

  ## prior precision on path coefficients (and plant effect)
  ##
  tau <- 0.1
  ## parameter for gamma prior on residual variance for photo
  ##
  nu <- 1

  if (model == 1) {
    jags.data <- c("leaf",
                   "phys",
                   "photo",
                   "humid",
                   "temp",
                   "plant",
                   "Omega.lf",
                   "Omega.ph",
                   "nu.lf",
                   "nu.ph",
                   "nu",
                   "tau",
                   "n.plant",
                   "n.samp")
  } else if (model == 2) {
    jags.data <- c("leaf",
                   "phys",
                   "photo",
                   "cond",
                   "humid",
                   "temp",
                   "plant",
                   "Omega.lf",
                   "Omega.ph",
                   "nu.lf",
                   "nu.ph",
                   "nu",
                   "tau",
                   "n.plant",
                   "n.samp")
  } else if (model == 3) {
    jags.data <- c("leaf",
                   "phys",
                   "photo",
                   "cond",
                   "humid",
                   "temp",
                   "plant",
                   "Omega.lf",
                   "Omega.ph",
                   "nu.lf",
                   "nu.ph",
                   "nu",
                   "tau",
                   "n.plant",
                   "n.samp")
  } else if (model == 4) {
    jags.data <- c("leaf",
                   "phys",
                   "photo",
                   "humid",
                   "temp",
                   "plant",
                   "Omega.lf",
                   "Omega.ph",
                   "nu.lf",
                   "nu.ph",
                   "nu",
                   "tau",
                   "n.plant",
                   "n.samp")
  } else if (model == 5) {
    jags.data <- c("leaf",
                   "phys",
                   "wue",
                   "humid",
                   "temp",
                   "plant",
                   "Omega.lf",
                   "Omega.ph",
                   "nu.lf",
                   "nu.ph",
                   "nu",
                   "tau",
                   "n.plant",
                   "n.samp")
  }
  if (model == 1) {
    jags.par <- c("phx.lf",
                  "phx.ph",
                  "rho.lf",
                  "rho.ph",
                  "sigmasq.photo",
                  "sigmasq.indiv.photo",
                  "sigmasq.indiv.lftemp",
                  "sigmasq.indiv.fluor",
                  "sigmasq.indiv.trans",
                  "sigmasq.indiv.cond",
                  "beta.lftemp.sla",
                  "beta.lftemp.lwr",
                  "beta.lftemp.lfarea",
                  "beta.lftemp.spi",
                  "beta.lftemp.sdi",
                  "beta.lftemp.humid",
                  "beta.lftemp.temp",
                  "beta.fluor.sla",
                  "beta.fluor.lwr",
                  "beta.fluor.lfarea",
                  "beta.fluor.spi",
                  "beta.fluor.sdi",
                  "beta.fluor.humid",
                  "beta.fluor.temp",
                  "beta.trans.sla",
                  "beta.trans.lwr",
                  "beta.trans.lfarea",
                  "beta.trans.spi",
                  "beta.trans.sdi",
                  "beta.trans.humid",
                  "beta.trans.temp",
                  "beta.cond.sla",
                  "beta.cond.lwr",
                  "beta.cond.lfarea",
                  "beta.cond.spi",
                  "beta.cond.sdi",
                  "beta.cond.humid",
                  "beta.cond.temp",
                  "beta.photo.lftemp",
                  "beta.photo.fluor",
                  "beta.photo.trans",
                  "beta.photo.cond",
                  "beta.photo.humid",
                  "mu.lf",
                  "mu.ph",
                  "mu.photo")
  } else if (model == 2) {
    jags.par <- c("phx.lf",
                  "phx.ph",
                  "rho.lf",
                  "rho.ph",
                  "sigmasq.photo",
                  "sigmasq.cond",
                  "sigmasq.indiv.lftemp",
                  "sigmasq.indiv.fluor",
                  "sigmasq.indiv.trans",
                  "sigmasq.indiv.cond",
                  "beta.lftemp.sla",
                  "beta.lftemp.lwr",
                  "beta.lftemp.lfarea",
                  "beta.lftemp.spi",
                  "beta.lftemp.sdi",
                  "beta.lftemp.humid",
                  "beta.lftemp.temp",
                  "beta.fluor.sla",
                  "beta.fluor.lwr",
                  "beta.fluor.lfarea",
                  "beta.fluor.spi",
                  "beta.fluor.sdi",
                  "beta.fluor.humid",
                  "beta.fluor.temp",
                  "beta.trans.sla",
                  "beta.trans.lwr",
                  "beta.trans.lfarea",
                  "beta.trans.spi",
                  "beta.trans.sdi",
                  "beta.trans.humid",
                  "beta.trans.temp",
                  "beta.cond.sla",
                  "beta.cond.lwr",
                  "beta.cond.lfarea",
                  "beta.cond.spi",
                  "beta.cond.sdi",
                  "beta.cond.trans",
                  "beta.cond.lftemp",
                  "beta.photo.lftemp",
                  "beta.photo.fluor",
                  "beta.photo.cond",
                  "beta.photo.humid",
                  "mu.lf",
                  "mu.ph",
                  "mu.photo")
  } else if (model == 3) {
    jags.par <- c("phx.lf",
                  "phx.ph",
                  "rho.lf",
                  "rho.ph",
                  "sigmasq.photo",
                  "sigmasq.cond",
                  "sigmasq.indiv.lftemp",
                  "sigmasq.indiv.fluor",
                  "sigmasq.indiv.trans",
                  "sigmasq.indiv.cond",
                  "beta.lftemp.sla",
                  "beta.lftemp.lwr",
                  "beta.lftemp.lfarea",
                  "beta.lftemp.spi",
                  "beta.lftemp.sdi",
                  "beta.lftemp.humid",
                  "beta.lftemp.temp",
                  "beta.fluor.sla",
                  "beta.fluor.lwr",
                  "beta.fluor.lfarea",
                  "beta.fluor.spi",
                  "beta.fluor.sdi",
                  "beta.fluor.humid",
                  "beta.fluor.temp",
                  "beta.trans.sla",
                  "beta.trans.lwr",
                  "beta.trans.lfarea",
                  "beta.trans.spi",
                  "beta.trans.sdi",
                  "beta.trans.humid",
                  "beta.trans.temp",
                  "beta.cond.trans",
                  "beta.cond.lftemp",
                  "beta.photo.lftemp",
                  "beta.photo.fluor",
                  "beta.photo.cond",
                  "beta.photo.humid",
                  "beta.photo.temp",
                  "mu.lf",
                  "mu.ph",
                  "mu.photo")
  } else if (model == 4) {
    jags.par <- c("phx.lf",
                  "phx.ph",
                  "rho.lf",
                  "rho.ph",
                  "sigmasq.photo",
                  "sigmasq.indiv.lftemp",
                  "sigmasq.indiv.fluor",
                  "sigmasq.indiv.cond",
                  "beta.lftemp.sla",
                  "beta.lftemp.lwr",
                  "beta.lftemp.lfarea",
                  "beta.lftemp.spi",
                  "beta.lftemp.sdi",
                  "beta.lftemp.humid",
                  "beta.lftemp.temp",
                  "beta.fluor.sla",
                  "beta.fluor.lwr",
                  "beta.fluor.lfarea",
                  "beta.fluor.spi",
                  "beta.fluor.sdi",
                  "beta.fluor.humid",
                  "beta.fluor.temp",
                  "beta.cond.sla",
                  "beta.cond.lwr",
                  "beta.cond.lfarea",
                  "beta.cond.spi",
                  "beta.cond.sdi",
                  "beta.cond.humid",
                  "beta.cond.temp",
                  "beta.trans.sla",
                  "beta.trans.lwr",
                  "beta.trans.lfarea",
                  "beta.trans.spi",
                  "beta.trans.sdi",
                  "beta.trans.humid",
                  "beta.trans.temp",
                  "beta.photo.lftemp",
                  "beta.photo.fluor",
                  "beta.photo.cond",
                  "mu.lf",
                  "mu.ph",
                  "mu.photo")
  } else if (model == 5) {
    jags.par <- c("phx.lf",
                  "phx.ph",
                  "rho.lf",
                  "rho.ph",
                  "sigmasq.wue",
                  "sigmasq.indiv.lftemp",
                  "sigmasq.indiv.fluor",
                  "sigmasq.indiv.cond",
                  "beta.lftemp.sla",
                  "beta.lftemp.lwr",
                  "beta.lftemp.lfarea",
                  "beta.lftemp.spi",
                  "beta.lftemp.sdi",
                  "beta.lftemp.humid",
                  "beta.lftemp.temp",
                  "beta.fluor.sla",
                  "beta.fluor.lwr",
                  "beta.fluor.lfarea",
                  "beta.fluor.spi",
                  "beta.fluor.sdi",
                  "beta.fluor.humid",
                  "beta.fluor.temp",
                  "beta.cond.sla",
                  "beta.cond.lwr",
                  "beta.cond.lfarea",
                  "beta.cond.spi",
                  "beta.cond.sdi",
                  "beta.cond.humid",
                  "beta.cond.temp",
                  "beta.trans.sla",
                  "beta.trans.lwr",
                  "beta.trans.lfarea",
                  "beta.trans.spi",
                  "beta.trans.sdi",
                  "beta.trans.humid",
                  "beta.trans.temp",
                  "beta.wue.lftemp",
                  "beta.wue.fluor",
                  "beta.wue.cond",
                  "mu.lf",
                  "mu.ph",
                  "mu.wue")
  }

  dat.fit <- jags(data=jags.data,
                  inits=NULL,
                  parameters=jags.par,
                  model.file=model.file,
                  n.chains=n.chains,
                  n.burnin=n.burnin,
                  n.iter=n.iter,
                  n.thin=n.thin,
                  DIC=TRUE,
                  working.directory=".")
  opt.old <- options(width=120)
  filename <- paste("results-path-model-",
                    model,
                    "-",
                    gsub(":", "-",
                         gsub(" ", "-", Sys.time())),
                    ".txt",
                    sep="")
  if (!debug) {
    sink(filename, split=TRUE)
  }
  cat(label, "\n", model.file, "\n\n")
  print(dat.fit, digits.summary=3)
  if (!debug) {
    sink()
  }
  dat.fit
}

analysis <- function(dehoop, kleinm, model) {
  dehoop.fit <- fit(dehoop, "Path analysis results at De Hoop", model)
  kleinm.fit <- fit(kleinm, "Path analysis results at Kleinmond", model)

  filename <- paste("results-path-",
                    gsub(":", "-",
                         gsub(" ", "-", Sys.time())),
                    ".Rsave",
                    sep="")
  save(dehoop,
       dehoop.fit,
       kleinm,
       kleinm.fit,
       file=filename)
}

raw <- read.csv("physiology-path-analysis.csv",
                header=TRUE,
                na.strings=".")
clean <- data.frame(sla=standardize(raw$SLA),
                    lwr=standardize(raw$LWR),
                    lfarea=standardize(raw$area_cm2),
                    spi=standardize(raw$SPI),
                    sdi=standardize(raw$stomata_per_mm2),
                    lftemp=standardize(raw$leaf_surface_temp),
                    fluor=standardize(raw$Y_avg),
                    cond=standardize(raw$Cond),
                    trans=standardize(raw$transp),
                    photo=standardize(raw$photosynthesis),
                    wue=standardize(raw$photosynthesis/raw$Cond),
                    humid=standardize(raw$relative_humidity),
                    temp=standardize(raw$air_temp),
                    site=raw$site,
                    plant=raw$plant)

## extract data from each site for separate analyses
##
dehoop <- droplevels(subset(clean, site=="De_Hoop"))
kleinm <- droplevels(subset(clean, site=="Kleinmond"))
## make plant a factor in each site (do this separately so that
## indices are consecutive in each site)
##
dehoop$plant <- as.numeric(as.factor(dehoop$plant))
kleinm$plant <- as.numeric(as.factor(kleinm$plant))

