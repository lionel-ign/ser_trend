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
yy[4,"year_fi"] <- 2023
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
# compute quadratic mean diameter  from basal area and number of trees
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
write.csv(d_all, "LIF/Croissance/ser_trend/data/dat_dendro.csv", row.names = FALSE)

## compute climatic conditions in the ser
# shapefile of SER
ser <- vect("LIF/IFN_stuff/data/geodata/ser_27572.gpkg") 
ser_2154 <- project(ser, "epsg:2154")

## anomalies in mean temperature and water deficit of growing season
ltemp <- list.files("LIF/Meteo/DIGICLIM/tmoy", pattern = "_13", full.names = TRUE)[1:60]
ras <- rast(lapply(ltemp, function(x) rast(x)))
names(ras) <- 1961:2020
# compute 30-years reference (1973-2002)
ras_ref <- app(ras[[13:42]], mean)
## the difference
ras_diff <- ras - ras_ref
## average per ser
ser_diff <- terra::extract(ras_diff, ser_2154, median, na.rm=TRUE)
ser_diff$ser <- ser_2154$codeser
ser_diff %>%
  pivot_longer(as.character(1961:2020)) %>%
  mutate(year = as.numeric(name)) %>%
  rename(tmoy = value) %>%
  select(ser, year, tmoy) -> tmoy
# average per 5-year window
tmoy %>%
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
  group_by(ser, win) %>%
  summarise(tmoy = mean(tmoy)) -> tmoy5

## now for water deficit of growing season (march-august)
lwat <- list.files("LIF/Meteo/DIGICLIM/deth", pattern = "_et|pr", full.names = TRUE)[1:120]
ras <- rast(lapply(lwat, function(x) rast(x)))
## sum deficit over spring and summer
ras <- tapp(ras, rep(1:60, each = 2), sum)
names(ras) <- 1961:2020
# compute 30-years reference (1973-2002)
ras_ref <- app(ras[[13:42]], mean)
## the difference
ras_diff <- ras - ras_ref
## average per ser
ser_diff <- terra::extract(ras_diff, ser_2154, median, na.rm=TRUE)
ser_diff$ser <- ser_2154$codeser
ser_diff %>%
  pivot_longer(as.character(1961:2020)) %>%
  mutate(year = as.numeric(name)) %>%
  rename(deth = value) %>%
  select(ser, year, deth) -> deth
# per 5-year window
deth %>%
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
  group_by(ser, win) %>%
  summarise(deth = mean(deth)) -> deth5

# combine the metrics
clim <- inner_join(tmoy, deth, by = c("ser", "year"))
clim5 <- inner_join(tmoy5, deth5, by = c("ser", "win"))

write_csv(clim, "~/LIF/Croissance/ser_trend/data/dat_clim.csv")
write_csv(clim5, "~/LIF/Croissance/ser_trend/data/dat_clim5.csv")


## get per ser the tree species making up 80% of the volume
spe <- ocre_nm(ans = 2005:2023, vent_p = "ser_86", vent_a = "ess",
               vars = "V", moy_ha = TRUE)
# species code
ess <- inventR::ListDataMod("ess")

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
  inner_join(ess[,c("mode", "libelle")], by = c("ess" = "mode")) -> dd

write.csv(dd,
          "LIF/Croissance/ser_trend/output/ser_species.csv",
          row.names = FALSE)

## compute volume and productivity of dead and harvested trees

# compute pv and v for dead and harvested trees
conn <- inventR::connect_db()
g3a <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "g3arbre"))
e2p <- tbl(conn, dbplyr::in_schema("inv_exp_nm", "e2point"))

# first grab values for living trees
g3a %>%
  select(npp, w, pv, v) %>%
  inner_join(e2p, by = "npp") %>%
  group_by(incref, greco, npp, poids) %>%
  summarise(v = sum(w * v, na.rm = TRUE),
            pv = sum(w * pv, na.rm = TRUE)) %>%
  collect() %>%
  mutate(vinit = v - pv) %>%
  group_by(greco, incref) %>%
  summarise(v.alive = sum(v * poids) / sum(poids),
            pv.alive = sum(pv * poids) / sum(poids)) -> vpv_alive

# then grab v of dead or harvested trees
g3a %>%
  select(npp, w, v, pv, veget5) %>%
  filter(veget5 %in% c("5", "6", "M", "7", "2", "T", "C")) %>%
  inner_join(e2p, by = "npp") %>%
  group_by(incref, greco, npp, poids) %>%
  summarise(vd = sum(w * v, na.rm = TRUE),
            pvd = sum(w * pv, na.rm = TRUE)) %>%
  collect() %>%
  group_by(incref, greco) %>%
  summarise(v.deathp = sum(vd * poids) / sum(poids),
            pv.deathp = sum(pvd * poids) / sum(poids)) %>%
  mutate(incref = incref + 5) -> vpv_dis

# then grab pv and v of dead trees
g3a %>%
  select(npp, w, v, pv, veget5) %>%
  filter(veget5 %in% c("5", "M", "2")) %>%
  inner_join(e2p, by = "npp") %>%
  group_by(incref, greco, npp, poids) %>%
  summarise(vd = sum(w * v, na.rm = TRUE),
            pvd = sum(w * pv, na.rm = TRUE)) %>%
  collect() %>%
  group_by(incref, greco) %>%
  summarise(v.death = sum(vd * poids) / sum(poids),
            pv.death = sum(pvd * poids) / sum(poids)) %>%
  mutate(incref = incref + 5) -> vpv_dis2

# merge the two
vpv_alive %>%
  full_join(vpv_dis, by = c("greco", "incref")) %>%
  full_join(vpv_dis2, by = c("greco", "incref")) %>%
  mutate(alive = pv.alive / v.alive,
         alive_dead = ifelse(!is.na(v.death), (pv.death + pv.alive) /
                               (v.death + v.alive), alive),
         alive_dead_harvested = ifelse(!is.na(v.deathp), (pv.deathp + pv.alive) /
                                         (v.deathp + v.alive), alive)) %>%
  pivot_longer(alive:alive_dead_harvested) %>%
  filter(!is.na(greco),
         incref > 4) -> dd

dd$name <- factor(dd$name, 
                  labels = c("alive", "alive + dead", "alive + dead +\nharvested"))

dd %>%
  mutate(year = incref + 2005) %>%
  select(greco, year, name, value) %>%
  filter(year > 2014) -> ddd

ddd$greco <- factor(ddd$greco, labels = c("Grand ouest", "Centre nord",
                                          "Grand est", "Vosges", "Jura",
                                          "Sud ouest", "Massif central",
                                          "Alpes", "Pyrénées", "Méditerranée",
                                          "Corse"))

write.csv(ddd, "data/dat_deadharvested.csv", row.names = FALSE)
