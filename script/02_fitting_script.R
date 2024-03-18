## model fitting script for the analysis of the temporal trend of prodctivity
## at the ser level

# load libraries
library(rstan)
library(brms)
options(mc.cores = parallel::detectCores())
library(tidyverse)
#library(pls)

# load data (created in the 01_create_data.R script)
d_all <- readRDS("LIF/Croissance/ser_trend/data/d_allobj.rds")

# some fine tuning of the data to ease model fit
d_all$pvv <- with(d_all, (PV / V) * 100)
d_all$mqd_std <- scale(d_all$mqd)
d_all$year_std <- scale(d_all$year)
d_all$year2.1 <- poly(d_all$year, 2)[,1]
d_all$year2.2 <- poly(d_all$year, 2)[,2]
d_all$year3.1 <- poly(d_all$year, 3)[,1]
d_all$year3.2 <- poly(d_all$year, 3)[,2]
d_all$year3.3 <- poly(d_all$year, 3)[,3]
d_all$grecorand <- factor(d_all$greco,
                          labels = LETTERS[1:11])


# a linear model of year effect
mm <- brm(pvv ~ mqd_std + year_std + (1 | greco / ser) +
            (0 + year_std | greco / ser),
          data = d_all)

# save this
saveRDS(mm, "LIF/Croissance/ser_trend/model/model_mm.rds")

# quadratic effect
mm2 <- brm(pvv ~ mqd_std + year2.1 + year2.2 + (1 | grecorand / ser) +
            (0 + year2.1 | grecorand / ser) + (0 + year2.2 | grecorand / ser),
          data = d_all)

# save this
saveRDS(mm2, "LIF/Croissance/ser_trend/model/model_mm2.rds")

# cubic effects
mm3 <- brm(pvv ~ mqd_std + year3.1 + year3.2 + year3.3 + (1 | grecorand / ser) +
             (0 + year3.1 | grecorand / ser) + (0 + year3.2 | grecorand / ser) +
             (0 + year3.3 | grecorand / ser),
           data = d_all)

# save this
saveRDS(mm3, "LIF/Croissance/ser_trend/model/model_mm3.rds")

# compare between models
c_loo <- loo(mm, mm2, mm3)


# climatic model with first PC axis on difference to 30-year reference
meteo_ser <- read.csv("LIF/Croissance/ser_trend/data/meteo_ser_ref_pied.csv")
pcc <- prcomp(meteo_ser[,3:14], scale. = TRUE)

meteo_ser$pc1 <- pcc$x[,1]

meteo_ser %>%
  rename(year = win) %>%
  select(ser, year, pc1, dethspring, dethsummer) %>%
  inner_join(d_all[,c("year", "ser", "grecorand", "V",
                      "PV", "pvv", "mqd_std")],
             by = c("ser", "year")) %>%
  mutate(dethsus = scale(dethsummer),
         dethsps = scale(dethspring)) -> d2

d2$pc1.1 <- poly(d2$pc1, 2)[,1]
d2$pc1.2 <- poly(d2$pc1, 2)[,2]

write.csv(d2, "LIF/Croissance/ser_trend/data/d2.csv", row.names = FALSE)

m_d <- brm(pvv ~ mqd_std + pc1.1 + pc1.2 + dethsus + 
             dethsps + pc1.1:dethsus + pc1.1:dethsps +
              (1 | grecorand / ser) +
             (0 + pc1.1 | grecorand / ser) + (0 + pc1.2 | grecorand / ser) +
             (0 + dethsus | grecorand / ser) + 
             (0 + dethsps | grecorand / ser) +
             (0 + pc1.1:dethsus | grecorand / ser) +
             (0 + pc1.1:dethsps | grecorand / ser),
           data = d2)

# save this
saveRDS(m_d, "LIF/Croissance/ser_trend/model/model_md.rds")

# fit the reduced models
m_1 <- brm(pvv ~ mqd_std + pc1.1 + 
             (1 | grecorand / ser) +
             (0 + pc1.1 | grecorand / ser),
           data = d2)

m_2 <- brm(pvv ~ mqd_std + dethsps +
             dethsus +
             (1 | grecorand / ser) +
             (0 + dethsps | grecorand / ser) + 
             (0 + dethsus | grecorand / ser),
           data = d2)

m_3 <- brm(pvv ~ mqd_std + pc1.1 + pc1.2 +
             (1 | grecorand / ser) +
             (0 + pc1.1 | grecorand / ser) +
             (0 + pc1.2 | grecorand / ser),
           data = d2)

m_4 <- brm(pvv ~ mqd_std + dethsps +
             dethsus + pc1.1 +
             (1 | grecorand / ser) +
             (0 + pc1.1 | grecorand / ser) +
             (0 + dethsps | grecorand / ser) + 
             (0 + dethsus | grecorand / ser),
           data = d2)



############# Code for SI analysis ##########


