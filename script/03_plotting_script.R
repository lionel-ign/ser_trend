######################################
##
## Companion Rscript to the manuscript:
## Turning point in forest productivity revealed from 40 years
## of national forest inventory data
## Author: Lionel Hertzog
## contact: lionel.hertzog@ign.fr
## Date: 12/04/2024
## Aim: this scripts creates the plots shown in the manuscript
##
######################################

## set working directory where the repository has been cloned
setwd("~/LIF/Croissance/ser_trend/")

## load libraries
library(brms)
library(tidyverse)
library(sf)
library(patchwork)

## load data
d_all <- read.csv("data/dat_dendro.csv")
clim5 <- read.csv("data/dat_clim5.csv")
d_climate <- inner_join(d_all, clim5, by = c("ser", "year" = "win"))
ddd <- read.csv("data/dat_deadharvested.csv")

## load model (output from script 02_fitting_script.R)
mm2 <- readRDS("model/model_mm2.rds")
m_c <- readRDS("model/m_climate.rds")


## load ser shape
ser_sf <- st_read("~/LIF/IFN_stuff/data/geodata/ser_27572.gpkg")

##### Extract fitted model coefficient per SER ######

## extract model coefficient from model mm2
n_ser <- length(unique(d_all$ser))
n_grecoser <- d_all %>%
  group_by(greco) %>%
  summarise(N = length(unique(ser)))
## get linear part of the equation
dd_year1 <- as_draws_df(mm2, variable = "year2.1", regex = TRUE)
## compute per ser the coefficient
pp_lin <- dd_year1[,rep(1, n_ser)] +
  as.matrix(dd_year1[,4:14][, rep(1:11, times = n_grecoser$N)]) +
  dd_year1[,15:99]
pp_lin_df <- as.data.frame(posterior_summary(pp_lin, probs = c(0.1, 0.9)))
# add ser name  
ll <- names(dd_year1)[15:99]
pp_lin_df$id <- substr(ll, 17, nchar(ll))
pp_lin_df <- separate(pp_lin_df, id, into = c("greco", "ser", "drop"), sep = "_|,")
pp_lin_df$greco <- factor(pp_lin_df$greco,
                          labels = c("Grand ouest", "Centre nord",
                                     "Grand est", "Vosges", "Jura",
                                     "Sud ouest", "Massif central",
                                     "Alpes", "Pyrénnés", "Méditerranée",
                                     "Corse"))
## posterior probability of negative coefficient
pp_lin_df$prob <- apply(pp_lin, 2, function(x) sum(x<0)) / 4000


### get the quadratic term of the model
dd_year2 <- as_draws_df(mm2, variable = "year2.2", regex = TRUE)
# posterior summary per ser
pp_quad <- dd_year2[,rep(1, n_ser)] +
  as.matrix(dd_year2[,4:14][, rep(1:11, times = n_grecoser$N)]) +
  dd_year2[,15:99]
pp_quad_df <- as.data.frame(posterior_summary(pp_quad, probs = c(0.1, 0.9)))
# add ser name  
ll <- names(dd_year2)[15:99]
pp_quad_df$id <- substr(ll, 17, nchar(ll))
pp_quad_df <- separate(pp_quad_df, id, into = c("greco", "ser", "drop"), sep = "_|,")
pp_quad_df$greco <- factor(pp_quad_df$greco,
                           labels = c("Grand ouest", "Centre nord",
                                      "Grand est", "Vosges", "Jura",
                                      "Sud ouest", "Massif central",
                                      "Alpes", "Pyrénnés", "Méditerranée",
                                      "Corse"))
pp_quad_df$prob <- apply(pp_quad, 2, function(x) sum(x<0)) / 4000



##### Figures in main document #####

### Figure 1
## trend per ser with greco/ser map

## part 1: derive predictions plus
## some data wraggling to format the data
y_seq <- min(d_all$year):max(d_all$year)
ppo <- poly(d_all$year, 2)
ppo_pred <- predict(ppo, y_seq)

yy <- data.frame(year = y_seq,
                 year2.1 = ppo_pred[,1],
                 year2.2 = ppo_pred[,2],
                 mqd_std = 0)

#yy$mqd_std <- 0
yy_ser <- yy[rep(1:nrow(yy), each = length(unique(d_all$ser))),]
yy_ser$ser <- rep(unique(d_all$ser), times = nrow(yy))
yy_ser$grecorand <- substr(yy_ser$ser, 1, 1)
yy_ser$greco <- factor(yy_ser$grecorand,
                       labels = c("Grand ouest", "Centre nord",
                                  "Grand est", "Vosges", "Jura",
                                  "Sud ouest", "Massif central",
                                  "Alpes", "Pyrénées", "Méditerranée",
                                  "Corse"))

## the fitted temporal trend
pp_ser <- fitted(mm2, newdata = yy_ser, probs = c(0.1, 0.9))
pp_ser <- cbind(pp_ser, yy_ser)
pp_ser$coll <- ifelse(pp_quad_df$prob >= 0.9, "quad. trend",
                     ifelse(pp_lin_df$prob >= 0.9, "lin. decline", "no trend"))

## part 2: derive the maps
# a map for the greco/ser
ser_sf %>%
  mutate(greco = substr(codeser, 1, 1)) %>%
  group_by(greco) %>%
  summarise(geom = st_union(geom)) -> greco

greco$greco <- factor(greco$greco, labels = c("Grand ouest", "Centre nord",
                                              "Grand est", "Vosges", "Jura",
                                              "Sud ouest", "Massif central",
                                              "Alpes", "Pyrénées", "Méditerranée",
                                              "Corse"))
ser_sf$clim <- rep(c("Oceanic", "Oceanic", "Semi-continental",
                     "Mountain", "Mountain", "Oceanic",
                     "Mountain", "Mountain", "Mountain",
                     "Mediterranean", "Mediterranean"),
                   times = c(6, 21, 8, 2, 2, 12, 14, 6, 5, 7, 3))

gmap <- ggplot() +
  geom_sf(data=ser_sf, color="white", aes(fill=clim)) +
  geom_sf(data=greco, color="black", fill=NA, size = 2) +
  geom_sf_label(data=greco, aes(label = greco)) +
  scale_fill_discrete(name = "Climatic zone:") +
  guides(fill = guide_legend(direction = "horizontal",
                             position = "bottom", nrow = 2)) +
  theme_void()

# grab the legen
gleg <- ggpubr::get_legend(gmap)

# remove the legend
gmap2 <- gmap + guides(fill="none")

# define color scheme
cols <- c("lin. decline" = "#F8766D",
          "no trend" = "#00BA38",
          "quad. trend" = "#619CFF")
# make the plots per SER
gpred <- plyr::dlply(pp_ser, plyr::as.quoted("greco"),
                     function(g){ggplot(g, aes(x=year,
                                               y=Estimate,
                                               ymin=Q10, ymax=Q90,
                                               group = ser)) +
                         geom_ribbon(alpha = 0.1, aes(fill = coll)) +
                         geom_line(aes(color = coll)) +
                         theme_classic() +
                         scale_color_manual(values = cols) +
                         scale_fill_manual(values = cols) +
                         scale_x_continuous(breaks = c(1980, 2000, 2020)) +
                         labs(x = "",
                              y = "",
                              title = unique(g$greco)) +
                         guides(fill="none", color="none")})

# grab the trend legend
gfoo <- ggplot(pp_ser, aes(x=year, y=Estimate, ymin=Q10, ymax=Q90)) +
  geom_line(aes(color=coll)) +
  geom_ribbon(aes(fill=coll), alpha=0.1) +
  scale_color_discrete(name="Trend shape") +
  scale_fill_discrete(name="Trend shape")
gpleg <- ggpubr::get_legend(gfoo)

# define grid
grr <- "
#ABC#
DEEEF
GEEEH
IEEEJ
#KLM#"

# add the x and y-axis label to specific plots
gmc <- gpred[[7]] + labs(y = "Predicted PR (%)")
gme <- gpred[[10]] + labs(x = "Year")

# and the plot
gf <- ggpubr::as_ggplot(gpleg) + gpred[[2]] + gpred[[3]] + gpred[[1]] +
  gmap + gpred[[4]] +
  gmc + gpred[[5]] +  gpred[[6]] +  gpred[[8]] +
  gpred[[9]] + gme + gpred[[11]] + plot_layout(design = grr)

ggsave("figures/fig1.png", gf, width=8, height=8)

### Figure 2

# merge pre-2005 data per biogeographical region
d_all %>%
  group_by(year, grecorand) %>%
  summarise(v = mean(V),
            pv = mean(PV),
            pvv = pv/v) %>%
  mutate(name = "alive only") %>%
  rename(value = pvv) %>%
  select(grecorand, year, name, value) -> df
df$greco <- factor(df$grecorand, labels = c("Grand ouest", "Centre nord",
                                            "Grand est", "Vosges", "Jura",
                                            "Sud ouest", "Massif central",
                                            "Alpes", "Pyrénées", "Méditerranée",
                                            "Corse"))


ddf <- rbind(ddd, df[,names(ddd)])
ddf$name <- factor(ddf$name, levels = c("alive only", "alive + dead",
                                        "alive + dead +\nharvested",
                                        "alive"))


gt <- ggplot(subset(ddf, name != "alive"), aes(x=year, y=value*100, color=name)) +
  geom_point(alpha=0.5) +
  facet_wrap(vars(greco), scales="free") +
  stat_smooth(se = FALSE, formula = y ~ poly(x, 2), method = "lm") +
  labs(x = "Inventory year",
       y = "Rate of productivity (%)") +
  scale_color_discrete(name = "Trees considered:") +
  scale_x_continuous(breaks = c(1980, 2000, 2020))

# plot trend in stocks of living trees per ser
d_all %>%
  select(ser, year, V) %>%
  nest(data = -ser) %>%
  mutate(mm = map(data, ~lm(V~year, data = .x)),
         tm = map(mm, broom::tidy)) %>%
  unnest(tm) %>%
  filter(term == "year") ->tt

gl <- ggplot(tt, aes(x=fct_reorder(ser, estimate, min),
                     y=estimate,
                     ymin=estimate-1.96*std.error,
                     ymax = estimate+1.96*std.error)) +
  geom_hline(yintercept = 0, linetype="dashed", color="red", linewidth=1.25) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  labs(x = "Forest region",
       y = "Trend in stocks of living trees") +
  theme(axis.text.y = element_text(size = 5))

ggsave("figures/vtrend.png", gt, width = 6, height = 9)

ga <- gl / gt
ggsave("figures/fig2.png", ga, width = 12, height = 12)


### Figure 3 cor from temp and clim model
n1 <- arrange(d_all, ser, year)
n1$mqd_std <- 0
tt <- fitted(mm2, newdata = n1, summary=FALSE)

nn <- arrange(d2, ser, year)
nn$mqd_std <- 0
nn$grecorand <- substr(nn$ser, 1, 1)

cc <- fitted(m_c, newdata = nn, summary=FALSE)

cctt <- cbind(pivot_longer(as.data.frame(tt), 1:637),
              pivot_longer(as.data.frame(cc), 1:637))
cctt$ser <- n1$ser[as.numeric(substr(cctt$name, 2, nchar(cctt$name)))]
cctt$iter <- rep(1:4000, each = 637)
names(cctt)[1:4] <- c("n1", "temp", "n2", "clim")
cctt %>%
  group_by(iter, ser) %>%
  summarise(r = cor(temp, clim)) -> co

co %>%
  group_by(ser) %>%
  summarise(cm = mean(r),
            cl = quantile(r, probs = 0.1),
            ch = quantile(r, probs = 0.9)) -> dc

dc %>%
  inner_join(ser_sf, by = c("ser" = "codeser")) -> dcf

g3 <- ggplot() +
  geom_sf(data=dcf, aes(fill = cm,
                        geometry = geom),
          color = "grey50") +
  scale_fill_continuous(type="viridis", direction = -1,
                        name = "Correlation\ncoefficient") +
  theme_void()

gd <- ggplot(dc, aes(x=cm)) +
  geom_density() +
  theme_classic() +
  labs(x = "Correlation coefficient",
       y = "Density") +
  theme(text = element_text(size = 10))

ga <- g3 +inset_element(gd, right = 1, top = 1, left = 0.75,
                        bottom = 0.75, align_to = "full")

ggsave("figures/fig3.png", ga, width=8, height=8)


### Figure 4
## model estimates from the climate model
bb <- coef(m_c, summary = FALSE)
bbd <- plyr::adply(bb$grecorand, c(2, 3), function(x) data.frame(m = mean(x),
                                                                 lci = quantile(x, probs = 0.05),
                                                                 uci = quantile(x, probs = 0.95)))
bbd <- subset(bbd, X2 %in% c("tmoy.1", "tmoy.2","deths"))
bbd$X1 <- factor(bbd$X1, labels = c("Grand ouest", "Centre nord",
                                    "Grand est", "Vosges", "Jura",
                                    "Sud ouest", "Massif central",
                                    "Alpes", "Pyrénées", "Méditerranée",
                                    "Corse"))
bbd$X2 <- factor(bbd$X2, labels = c("Temperature anomalies\nlinear",
                                    "Temperature anomalies\nquadratic",
                                    "Water deficit\nanomalies (growing season)"))
ff <- as.data.frame(fixef(m_c, probs = c(0.05, 0.95)))
ff$X2 <- rownames(ff)
names(ff)[1] <- "m"
ffs <- subset(ff, X2 %in% c("tmoy.1", "tmoy.2","deths"))
ffs$X2 <- factor(ffs$X2, labels = c("Water deficit\nanomalies (growing season)",
                                    "Temperature anomalies\nlinear",
                                    "Temperature anomalies\nquadratic"))

ge <- ggplot(bbd, aes(x = X1, y = m, ymin=lci, ymax=uci)) +
  geom_linerange() +
  geom_point() +
  geom_hline(data=ffs, aes(yintercept = m), linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(vars(X2), scales = "free", ncol = 2, nrow = 3) +
  coord_flip() +
  labs(x = "",
       y = "Estimate (90% CrI)")

ggsave("figures/fig4.png", ge)


##### Figures in SI #####

### Figure S1-3 trend in stand growing stock, productivity and quadratic mean diameter

## compute overall mean
d_all %>%
  group_by(year) %>%
  summarise(V = mean(V),
            PV = mean(PV),
            mqd = mean(mqd)) -> da

d_all$greco <- factor(d_all$grecorand, labels = c("Grand ouest", "Centre nord",
                                              "Grand est", "Vosges", "Jura",
                                              "Sud ouest", "Massif central",
                                              "Alpes", "Pyrénées", "Méditerranée",
                                              "Corse"))


g_v <- ggplot(d_all, aes(x=year, y=V)) +
  geom_jitter() +
  stat_smooth() +
  geom_line(data=da, color="red", linetype = "dashed", linewidth = 1.25) +
  facet_wrap(~greco,scales = "free") +
  labs(x = "Year",
       y = "Volume (m^3/ha)")

ggsave("LIF/Croissance/figures/figS1.png", g_v)


g_pv <- ggplot(d_all, aes(x=year, y=PV)) +
  geom_jitter() +
  geom_smooth() +
  geom_line(data=da, color="red", linetype = "dashed", size = 1.25)+
  facet_wrap(~greco,scales = "free") +
  labs(x = "Year",
       y = "Volume productivity (m^3 / ha / an)")

ggsave("LIF/Croissance/figures/figS2.png", g_pv)

g_mqd <- ggplot(d_all, aes(x=year, y=mqd)) +
  geom_jitter() +
  stat_smooth() +
  geom_line(data=da, color="red", linetype = "dashed", size = 1.25)+
  facet_wrap(~greco,scales = "free") +
  labs(x = "Year",
       y = "Mean quadratic diameter (m / ha)")

ggsave("LIF/Croissance/figures/figS3.png", g_mqd)

### Figure S4 relation productivity volume
gp <- ggplot(d_all, aes(x=V, y=PV)) +
  geom_point() +
  stat_smooth() +
  labs(x = "Stand growing stocks (m^3/ha)",
       y = "Productivity (m^3/ha/year)")

ggsave("figures/figS4.png", gp)

### Figure S5 trend in volume of living, dead and harvested trees
## cannot be run without inventR package
dt <- data.table(year = 2005:2018)
d1 <- dt[1:5, ocre_nm(ans=year, vars="V", vent_a = "veget5",
                      vent_p = "greco", moy_ha = TRUE), by = year]
d2 <- dt[6:7, ocre_nm(ans=year, vars="V", vent_a = "veget5",
                      vent_p = "greco", moy_ha = TRUE), by = year]
d3 <- dt[8:9, ocre_nm(ans=year, vars="V", vent_a = "veget5",
                      vent_p = "greco", moy_ha = TRUE), by = year]
d4 <- dt[10:14, ocre_nm(ans=year, vars="V", vent_a = "veget5",
                        vent_p = "greco", moy_ha = TRUE), by = year]
da <- rbind(d1, d2, d3, d4, use.names = FALSE)
da %>%
  filter(visite == "1") %>%
  mutate(status = case_when(veget5.veget5 %in% c("0", "(null)", "A", "N", "1") ~ "alive",
                            veget5.veget5 %in% c("2", "T", "M", "5", "C") ~ "dead",
                            veget5.veget5 %in% c("6", "7") ~ "harvested")) %>%
  group_by(year, greco.greco, status) %>%
  summarise(v = sum(moy_ha)) -> mm

mm$greco.greco <- factor(mm$greco.greco, labels = c("Grand ouest", "Centre nord",
                                                    "Grand est", "Vosges", "Jura",
                                                    "Sud ouest", "Massif central",
                                                    "Alpes", "Pyrénnés", "Méditerranée",
                                                    "Corse"))

g1 <- ggplot(subset(mm, year > 2009), aes(x=year, y=v, color=status)) +
  geom_point() +
  stat_smooth(method="lm") +
  facet_wrap(vars(greco.greco), scales="free_y") +
  labs(y = "Stem volume (m^3/ha)")

ggsave("~/LIF/Croissance/ser_trend/figures/figS5.png", g1)

## Figure S6
### get residuals
res_mm2 <- residuals(mm2)
### get predicted values
pred_mm2 <- posterior_summary(posterior_predict(mm2))
### residuals plot 1, res ~ fitted
res_pred <- as.data.frame(cbind(res_mm2, pred_mm2))
names(res_pred) <- paste(rep(c("res", "pred"), each = 4),
                         rep(c("ess", "err", "lci", "uci"), times = 2),
                         sep = "_")
res_pred$greco <- d_all$greco

res1 <- ggplot(res_pred, aes(x=pred_ess, xmin=pred_lci, xmax=pred_uci,
                             y=res_ess, ymin=res_lci,ymax=res_uci)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(alpha = 0.3) +
  geom_errorbarh(alpha = 0.3) +
  geom_point(alpha = 0.5) +
  facet_wrap(vars(greco), scales = "free") +
  labs(x = "Fitted values (95% CrI)",
       y = "Reiduals (95% CrI)")

ggsave("figures/figS6.png", res1)

### Figure S7
### residual plot 2, res ~ x
res_pred$year <- d_all$year
res_pred$mqd <- d_all$mqd
res_pred$v <- d_all$V

res2.1 <- ggplot(res_pred, aes(x=year, y=res_ess, group = year)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(x = "Years",
       y = "Residuals")

res2.2 <- ggplot(res_pred, aes(x=mqd, y=res_ess)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(x = "Mean quadratic diameter (m / ha)",
       y = "Residuals")

res2.3 <- ggplot(res_pred, aes(x=v, y=res_ess)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(x = "Volume (m^3 / ha)",
       y = "Residuals")

res2.4 <- pp_check(mm2, ndraws = 100) +
  labs(x = "Productivity (m^3/ha/an)",
       y = "Density") +
  scale_color_discrete(labels = c("Observation", "Prediction"),
                       name = "") +
  theme_gray()

resa <- ggpubr::ggarrange(res2.1, res2.2, res2.3, res2.4, nrow = 2, ncol = 2,
                          labels = "AUTO")
ggsave("figures/figS7.png", resa)

## figure S8 relation between productivity and quadratic mean diameter
# plot of mqd vs pvv
# gather fitted effect
cc <- conditional_effects(mm2, "mqd_std")$mqd_std
cc$mqd <- cc$mqd_std * sd(d_all$mqd) + mean(d_all$mqd)
cc$pvv <- cc$estimate__ * 100
cc$low <- cc$lower__ * 100
cc$high <- cc$upper__ * 100
d_all$pvv <- with(d_all, (PV/V) * 100)

gmqd <- ggplot(cc, aes(x=mqd, y=pvv)) +
  geom_ribbon(aes(ymin=low, ymax=high), alpha = 0.25) +
  geom_line() +
  geom_point(data=d_all) +
  labs(x = "Quadratic mean diameter (m/ha)",
       y = "Productivity rate rate (%/year)")

ggsave("LIF/Croissance/ser_trend/figures/figS8.png", gmqd)


### Figure S9 map of the quadratic trend coefficient
## join the ser shape to the coefs
pp_ser %>%
  filter(coll == "quad. trend") %>%
  group_by(ser) %>%
  slice_max(Estimate) -> dd

## add geoms info
ser_sf %>%
  left_join(dd, by = c("codeser" = "ser")) -> ssf

gg_osf <- ggplot() +
  geom_sf(data=ssf, aes(fill=year)) +
  viridis::scale_fill_viridis(direction=-1)

ggsave("figures/figS9.png", gg_osf)


### Figure S10 compute where the ser ended up at the end of the period
newd <- data.frame(ser = rep(unique(d_all$ser), each = 2),
                   year = rep(c(1982, 2022), times = 85),
                   year2.1 = rep(c(min(d_all$year2.1), max(d_all$year2.1)), times = 85),
                   year2.2 = rep(c(0.06030563, 0.05852646), times = 85),
                   mqd_std = 0)
newd$grecorand <- substr(newd$ser, 1, 1)
tp <- fitted(mm2, newdata = newd, summary = FALSE)
tpp <- as.data.frame(posterior_summary(tp[,seq(2, 170, 2)] - tp[,seq(1, 170, 2)],
                                       probs = c(0.1, 0.9)))
tpp$ser <- unique(d_all$ser)
tpp$quad <- pp_quad_df$prob

## indicate sig level
tpp$sig <- case_when(tpp$Q10 < 0 & tpp$Q90 < 0 ~ "sig",
                     tpp$Q10 > 0 & tpp$Q90 > 0 ~ "sig",
                     TRUE ~ "not sig.")

ser_sf %>%
  inner_join(tpp, by = c("codeser" = "ser")) %>%
  mutate(a = ifelse(sig == "sig", Estimate, NA)) -> ss

gs <- ggplot(ss) +
  geom_sf(aes(fill=a)) +
  scale_fill_gradient2(name = "Changes in\nproductivity\nrate")

ggsave("figures/figS10.png", gs, width = 8, height = 8)


### Figure S11
## average normal temperature per SER against predicted trend form
## load the climate normal data (all parameters)
ppp <- do.call(rbind,
               lapply(list.files(path = "data/",
                                 pattern = "serref", full.names = TRUE),
                      read.csv))
ppp <- subset(ppp, ser %in% pp_quad_df$ser)
ppp$greco <- substr(ppp$ser, 1, 1)
# create groups of ser
ppp %>%
  inner_join(pp_quad_df[,c("ser", "prob")], by = "ser") %>%
  inner_join(pp_lin_df[,c("ser", "prob")], by = "ser") %>%
  mutate(coll = case_when(prob.x >= 0.9 ~ "quad.trend",
                          prob.y >= 0.9 ~ "lin. decline",
                          .default = "no trend")) %>%
  pivot_wider(names_from = type,
              values_from = value) -> ppw
## separate lowland from mountain
ppw$moutain <- ifelse(ppw$greco%in%c("D", "E", "G", "H", "I", "K"), "mountain", "lowland")


# first plot with tmoy
hl <- data.frame(moutain = "lowland",
                 tmoy = 11.85)
g1 <- ggplot(ppw, aes(x=coll, y=tmoy/10)) + 
  geom_jitter(width = 0.1)+
  stat_summary(fun.data = "mean_cl_boot", color = "red") +
  geom_hline(data=hl, aes(yintercept=tmoy)) +
  facet_wrap(vars(moutain), scales="free") +
  labs(x = "", y = "Mean temperature (°C)")

# then with deth
hl <- data.frame(moutain = "mountain",
                 deth = 25)
g2 <- ggplot(ppw, aes(x=coll, y=deth/10)) + 
  geom_jitter(width = 0.1)+
  stat_summary(fun.data = "mean_cl_boot", color = "red") +
  geom_hline(data=hl, aes(yintercept=deth)) +
  facet_wrap(vars(moutain), scales="free") +
  labs(x = "Trend shape", y = "Mean water deficit (mm)")


ga <- g1 / g2
ggsave("figures/figS11.png", ga, width=8, height=8)









