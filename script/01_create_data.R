######################################
##
## Companion Rscript to the manuscript:
## Turning point in forest productivity revealed from 40 years
## of national forest inventory data
## Author: Lionel Hertzog
## contact: lionel.hertzog@ign.fr
## Date: 12/04/2024
## Aim: this script present the steps to create the datasets used in the analysis
## Note that the raw forest inventory data are available for internal use only
## the lines 1-80 are only for documentation purposes.
## The climatic data are not given in the repository but can be asked from:
## christian.piedallu@agroparistech.fr
##
############################

# load libraries
library(data.table)
library(inventR)
library(tidyverse)
library(terra)

## Derive inventory estimations

# create data with temporal window
yy <- data.table(year_st = seq(2005, 2020, by = 5),
                 year_fi = seq(2009, 2024, by = 5))
yy[4,"year_fi"] <- 2022
# from 2005 onwards
vpv_nm <- yy[, ocre_nm(ans = year_st:year_fi, vent_p = c("ser_86"), vars = c("V", "PV", "GTOT", "NT"), moy_ha = TRUE),
             by = .(year_st, year_fi)]
# from 1978 to 2004
vpv_am <- dep_cyc_anref[base == "dendro",
                        ocre_am(dep = dep, cyc = cyc, vent_p = "ser_86", 
                                vars = c("V", "PV", "GTOT", "NT"),
                                moy_ha = TRUE),
                        by = .(dep, cyc)]


# format results
vpv_am %>%
  inner_join(dep_cyc_anref, by = c("dep", "cyc")) %>% # add year info
  filter(anref < 2005, as.numeric(nbpoints) > 99) %>% # remove lines from after 2005 and with less than 100 pts
  mutate(year = case_when(anref < 1985 ~ 1982,
                              anref < 1990 ~ 1987,
                              anref < 1995 ~ 1992,
                              anref < 2000 ~ 1997,
                              anref < 2005 ~ 2002)) %>% # create year category per 5-years windows
  group_by(ser_86.ser_86, vars, year) %>%
  summarise(moy_ha = sum(as.numeric(nbpoints) * moy_ha, na.rm = TRUE) /
              sum(as.numeric(nbpoints), na.rm = TRUE),
            nbpoint = sum(as.numeric(nbpoints), na.rm=TRUE)) %>% # compute weighted means
  pivot_wider(names_from = vars, values_from = moy_ha) %>% # put in tab format
  rename(ser = ser_86.ser_86) %>%
  mutate(greco = substr(ser, 1, 1)) -> d_am


vpv_nm %>%
  mutate(nbpoints = as.numeric(nbpoints)) %>%
  filter(visite == 1,
         nbpoints > 29) %>% # only keep lines with more than 30 pts
  mutate(year = year_st + 2,
         ser = ser_86.ser_86,
         greco = substr(ser_86.ser_86, 1, 1)) %>%
  pivot_wider(id_cols = c(greco, ser, year),
              names_from = vars,
              values_from = moy_ha) -> d_nm

# put am and nm together
d_all <- rbind(d_am[,c("greco", "ser", "year", "V", "PV", "GTOT", "NT")],
               d_nm[,c("greco", "ser", "year", "V", "PV", "GTOT", "NT")])
# compute mqd from basal area and number of trees
d_all$mqd <- with(d_all, sqrt((GTOT*4) / (NT * pi)))

# change greco label
d_all$greco <- factor(d_all$greco, labels = c("Grand ouest", "Centre nord",
                                              "Grand est", "Vosges", "Jura",
                                              "Sud ouest", "Massif central",
                                              "Alpes", "Pyrénnés", "Provence",
                                              "Corse"))
# rate of growth
d_all$pvv <- with(d_all, PV / V)

# save this
write.csv(d_all, "LIF/Croissance/ser_trend/data/d_all.csv", row.names = FALSE)

## compute climatic conditions in the ser
# set some vars
out_5year <- list() # an object to put the 5-year ser data
out_year <- NULL # an object to put the yearly ser data
ser <- vect("LIF/IFN_stuff/data/geodata/ser_27572.gpkg") # shapefile of SER
ser_2154 <- project(ser, "epsg:2154")


## loop through the variables
for(var in c("tmin", "tmax", "tmoy", "deth")){
  ## read the data
  ll <- list.files(paste0("LIF/Meteo/DIGICLIM/", var), pattern = "au|hi|_et|pr")
  ll_yr <- plyr::laply(str_split(ll, "_"), "[", 2)
  ras <- rast(lapply(ll[which(ll_yr %in% as.character(1970:2020))],
                     function(x) rast(paste0("LIF/Meteo/DIGICLIM/", var, "/", x))))
  ## rename the layers
  names(ras) <- paste(rep(1970:2020, each = 4),
                      rep(c("autumn", "summer", "winter", "spring"), times = 51),
                      sep = "_")
  
  ## compute 30-year reference value (1980-2009)
  ras_ref <- tapp(ras[[41:160]], rep(1:4, times = 30), mean)
  ## save the mean reference value per ser
  ser_ref <-  terra::extract(ras_ref, ser_2154, mean, na.rm=TRUE)
  serr <- data.frame(type = var,
                     ser = ser_2154$codeser,
                     value = apply(ser_ref[,-1], 1, mean))
  write.csv(serr, paste0("LIF/Croissance/ser_trend/data/", var, "_serref.csv"),
            row.names = FALSE)
  
  ## compute yearly difference to reference value
  ras_diff <- (ras - ras_ref)

  ## average per ser
  ser_diff <- terra::extract(ras_diff, ser_2154, median, na.rm=TRUE)
  out_year <- rbind(out_year, ser_diff)

  ## average per 5-year window
  ser_diff %>%
    pivot_longer(-1) %>%
    separate(name, c("year", "season")) %>%
    filter(year > 1979) %>%
    mutate(win = case_when(year < 1985 ~ 1982,
                           year < 1990 ~ 1987,
                           year < 1995 ~ 1992,
                           year < 2000 ~ 1997,
                           year < 2005 ~ 2002,
                           year < 2010 ~ 2007,
                           year < 2015 ~ 2012,
                           year < 2020 ~ 2017,
                           year < 2023 ~ 2022)) %>%
    group_by(ID, win, season) %>%
    summarise(dev = mean(value, na.rm = TRUE)) %>%
    pivot_wider(names_from = season, values_from = dev,
                names_prefix = var) %>%
    mutate(ser = ser$codeser[ID]) -> ress

  ## put in one object with the other variables
  out_5year[[length(out_5year) + 1]] <- ress
}

# recuperate yearly data
out_year %>%
  mutate(var = rep(c("tmin", "tmax", "tmoy", "deth"), each = 86)) %>%
  pivot_longer('1970_autumn':'2020_spring') %>%
  separate(name, c("year", "season")) %>%
  filter(year > 1979) %>%
  mutate(ser = ser$codeser[ID]) %>%
  pivot_wider(names_from = c(var, season), values_from = value,
              names_sep = "") -> meteo_year

# save this
write.csv(meteo_year[,c(3, 2, 4:19)],
          "LIF/Croissance/ser_trend/data/meteo_ser_ref_pied_year.csv",
          row.names = FALSE)

# gather 5-year data
outt <- out[[1]][,-1]
outt %>%
  inner_join(out[[2]][,-1], by = c("ser", "win")) %>%
  inner_join(out[[3]][,-1], by = c("ser", "win")) %>%
  inner_join(out[[4]][,-1], by = c("ser", "win")) -> outt

# save
write.csv(outt[,c(6, 1, 2:5, 7:18)], "LIF/Croissance/ser_trend/data/meteo_ser_ref_pied.csv", 
          row.names = FALSE)

## create the climatic data
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

write.csv(d2, "LIF/Croissance/ser_trend/data/d_climate.csv", row.names = FALSE)



## get per ser the tree species making up 80% of the volume
spe <- ocre_nm(ans = 2005:2022, vent_p = "ser_86", vent_a = "ess",
               vars = "V", moy_ha = TRUE)
# species code
ess <- read.csv2("LIF/IFN_stuff/data/pommier/essence_code_clean.csv")
ess$ess <- sprintf("%02d", ess$ess)

spe %>%
  filter(visite == 1) %>%
  rename(ser = ser_86.ser_86,
         ess = ess.ess) %>%
  select(ser, ess, tot) %>%
  group_by(ser) %>%
  mutate(prop = tot / sum(tot)) %>%
  arrange(ser, desc(prop)) %>%
  mutate(cump = cumsum(prop)) %>%
  filter(cump <= 0.8 | prop >= 0.8) %>%
  inner_join(ess, by = "ess") -> dd

write.csv(dd,
          "LIF/Croissance/ser_trend/output/ser_species.csv",
          row.names = FALSE)