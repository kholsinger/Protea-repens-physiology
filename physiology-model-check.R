require(R2jags)
require(ggplot2)

rm(list=ls())

load(file="results-path-2014-08-23-19-50-40.Rsave")

check.fit <- function(dat, fit, label, suffix) {
  plot.data <- data.frame(var="LFTEMP",
                          obs=dat$lftemp,
                          pred=fit$BUGSoutput$mean$mu.ph[,1])
  plot.data <- rbind(plot.data,
                     data.frame(var="FLUOR",
                                obs=dat$fluor,
                                pred=fit$BUGSoutput$mean$mu.ph[,2]))
  plot.data <- rbind(plot.data,
                     data.frame(var="TRANS",
                                obs=dat$trans,
                                pred=fit$BUGSoutput$mean$mu.ph[,3]))
  plot.data <- rbind(plot.data,
                     data.frame(var="PHOTO",
                                obs=dat$photo,
                                pred=fit$BUGSoutput$mean$mu.photo))
  plot.data$resid <- plot.data$obs - plot.data$pred
  dev.new()
  p <- ggplot(plot.data, aes(x=obs, y=pred)) +
       geom_point() +
       xlab("Observed") +
       ylab("Predicted") +
       facet_wrap(~ var) +
       ggtitle(paste("Observed vs. Predicted for", label))
  print(p)
  filename <- paste("obs-vs-pred-", suffix, ".pdf", sep="")
  ggsave(filename)
  dev.new()
  p <- ggplot(plot.data, aes(x=pred, y=resid)) +
       geom_point() +
       xlab("Predicted") +
       ylab("Residual") +
       facet_wrap(~ var) +
       ggtitle(paste("Predicted vs. Residual for", label))
  print(p)
  filename <- paste("pred-vs-resid-", suffix, ".pdf", sep="")
  ggsave(filename)
}

check.fit(dehoop, dehoop.fit, "De Hoop", "dehoop")
check.fit(kleinm, kleinm.fit, "Kleinmond", "kleinm")
