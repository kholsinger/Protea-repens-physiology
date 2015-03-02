require(plyr)

combined <- read.csv("traits-environment.csv", na.strings=".", header=TRUE)

summary <- ddply(combined,
                 c("source_population"),
                 summarise,
                 sla = mean(SLA, na.rm=TRUE),
                 leaf.area = mean(leaf_area, na.rm=TRUE),
                 stomatal.density = mean(stomatal_density, na.rm=TRUE),
                 lwr = mean(LWR, na.rm=TRUE),
                 stomatal.pore.index=mean(SPI, na.rm=TRUE))

write.csv(summary,
          file="population-means.csv",
          row.names=FALSE)
