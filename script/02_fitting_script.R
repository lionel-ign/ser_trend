######################################
##
## Companion Rscript to the manuscript:
## Turning point in the productivity of western European forests
## associated with a climate change footprint
## Author: Lionel Hertzog
## contact: lionel.hertzog@ign.fr
## Date: 12/04/2024
## Aim: this scripts fits the models used in the manuscript
##
######################################

# set working directory where the repository has been cloned
setwd("~/LIF/Croissance/ser_trend/")

set.seed(20241124)
# load libraries
library(rstan)
library(brms)
options(mc.cores = parallel::detectCores())
library(tidyverse)

# load data (check the 01_create_data.R script)
d_all <- read.csv("data/dat_dendro.csv")
clim5 <- read.csv("data/dat_clim5.csv")
d2 <- inner_join(d_all, clim5, by = c("ser", "year" = "win"))


# some fine tuning of the data to ease model fit
d_all$mqd_std <- scale(d_all$mqd)
d_all$v_std <- scale(d_all$V)
d_all$year_std <- scale(d_all$year)
d_all$year2.1 <- poly(d_all$year, 2)[,1]
d_all$year2.2 <- poly(d_all$year, 2)[,2]
d_all$year3.1 <- poly(d_all$year, 3)[,1]
d_all$year3.2 <- poly(d_all$year, 3)[,2]
d_all$year3.3 <- poly(d_all$year, 3)[,3]
d_all$grecorand <- factor(substr(d_all$ser, 1, 1))


# a linear model of year effect
mm <- brm(pvv ~ mqd_std + year_std + (1 | greco / ser) +
            (0 + year_std | greco / ser),
          data = d_all)

# save this
saveRDS(mm, "model/model_mm.rds")

# quadratic effect
mm2 <- brm(pvv ~ mqd_std + year2.1 + year2.2 + (1 | grecorand / ser) +
            (0 + year2.1 | grecorand / ser) + (0 + year2.2 | grecorand / ser),
          data = d_all)

# save this
saveRDS(mm2, "model/model_mm2.rds")

# cubic effects
mm3 <- brm(pvv ~ mqd_std + year3.1 + year3.2 + year3.3 + (1 | grecorand / ser) +
             (0 + year3.1 | grecorand / ser) + (0 + year3.2 | grecorand / ser) +
             (0 + year3.3 | grecorand / ser),
           data = d_all)

# save this
saveRDS(mm3, "model/model_mm3.rds")

# compare between models
c_loo <- loo(mm, mm2, mm3)


## Climatic model with tmoy,  tmoy² and deth
pp <- poly(d2$tmoy, 2)
d2$tmoy.1 <- pp[,1] / sd(pp[,1])
d2$tmoy.2 <- pp[,2] / sd(pp[,2])
d2$deths <- scale(d2$deth)[,1]
d2$mqd_std <- scale(d2$mqd)[,1]
d2$grecorand <- factor(substr(d2$ser, 1, 1))

m_c <- brm(pvv ~ mqd_std + tmoy.1 + tmoy.2 + 
             deths + 
             (1 | grecorand / ser) +
             (0 + tmoy.1 | grecorand / ser) +
             (0 + tmoy.2 | grecorand /ser) +
             (0 + deths | grecorand / ser),
           data = d2)

saveRDS(m_c, "model/m_climate.rds")

