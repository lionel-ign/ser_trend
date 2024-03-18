## plotting script for the analysis of the temporal trends in productivity
## at the ser level

library(brms)
library(tidyverse)
library(sf)
library(patchwork)

## load data
d_all <- readRDS("LIF/Croissance/ser_trend/data/d_allobj.rds")
d2 <- read.csv("LIF/Croissance/ser_trend/data/d2.csv")

## load model (output from script 02_fitting_script.R)
mm2 <- readRDS("LIF/Croissance/ser_trend/model/model_mm2.rds")
m_d <- readRDS("LIF/Croissance/ser_trend/model/model_md.rds")
m_clim <- readRDS("LIF/Croissance/ser_trend/model/clim_models.rds")

## load ser shape
ser_sf <- st_read("LIF/IFN_stuff/data/geodata/ser_27572.gpkg")

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
                                  "Alpes", "Pyrénnés", "Méditerranée",
                                  "Corse"))

## the fitted temporal trend
pp_ser <- fitted(mm2, newdata = yy_ser, probs = c(0.1, 0.9))
pp_ser <- cbind(pp_ser, yy_ser)
pp_ser <- inner_join(pp_ser, ppw[,c("ser", "coll")], by = "ser")

## part 2: derive the maps
# a map for the greco/ser
ser_sf %>%
  mutate(greco = substr(codeser, 1, 1)) %>%
  group_by(greco) %>%
  summarise(geom = st_union(geom)) -> greco

greco$greco <- factor(greco$greco, labels = c("Grand ouest", "Centre nord",
                                              "Grand est", "Vosges", "Jura",
                                              "Sud ouest", "Massif central",
                                              "Alpes", "Pyrénnés", "Méditerranée",
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
  scale_fill_discrete(name = "Climatic zone") +
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
ABCD#
EFFFG
HFFFI
JFFFK
#LMN#"

# add the x and y-axis label to specific plots
gmc <- gpred[[7]] + labs(y = "Predicted AGR (%)")
gme <- gpred[[10]] + labs(x = "Year")

# and the plot
gf <- ggpubr::as_ggplot(gpleg) + gpred[[2]] + gpred[[3]] + gpred[[4]] + gpred[[1]] +
  gmap2 + ggpubr::as_ggplot(gleg) +
  gmc + gpred[[5]] +  gpred[[6]] +  gpred[[8]] +
  gpred[[9]] + gme + gpred[[11]] + plot_layout(design = grr)

ggsave("LIF/Croissance/ser_trend/figures/fig1.png", gf)

### Figure 2
## map of the quadratic trend coefficient
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

ggsave("LIF/Croissance/ser_trend/figures/fig2.png", gg_osf)


### Figure 3
## check correlation between trend shape and temperature anomalies
# an helper function to go through the models
get_corrs <- function(models, mm2, newdata, yy_ser, group="greco"){
  # potentially restrict the temporal model to the available climatic data
  yy_ser %>%
    arrange(ser, year) %>%
    inner_join(newdata[,c("year", "ser")], by = c("ser", "year")) -> yy_ser2
  
  # derive temporal model prediction
  pp_ser <- posterior_epred(mm2, newdata = yy_ser2)
  
  # derive model predictions
  mpred <- plyr::llply(models, function(model) posterior_epred(model,
                                                               newdata=newdata))
  
  # go through the posterior iterations and grab the correlation coeff
  mmat <- plyr::llply(mpred, function(mat) cbind(mat, pp_ser))
  out1 <- NULL
  out2 <- NULL
  if(group == "greco"){
    for(g in unique(newdata$grecorand)){
      id <- newdata$grecorand == g
      nb <- sum(id)
      corr <- plyr::laply(mmat, function(mat) apply(mat[,id], 1,
                                                    function(x) cor(x[1:nb], x[(nb + 1): (nb * 2)],
                                                                    method = "spearman")))
      
      # summarise
      corr_df <- plyr::adply(corr, 1, function(cor) data.frame(cor_m = mean(cor),
                                                               cor_l = quantile(cor, probs = 0.1),
                                                               cor_h = quantile(cor, probs = 0.9),
                                                               greco = g))
      corr_df$label <- c("T", "CWD", "T +\nT²", "T +\nCWD", "T +\nT² +\nCWD")
      out1 <- rbind(out1, corr_df)
      
      # get probs
      out2 <- rbind(out2, data.frame(greco = rep(g, 3),
                                     test = c("T+CWD > CWD", "T+CWD > T+T²",
                                              "T+CWD > T+T²+CWD"),
                                     probs = c(sum(corr[2,] < corr[4,]) / 4000,
                                               sum(corr[3,] < corr[4,]) / 4000,
                                               sum(corr[5,] < corr[4,]) / 4000)))
    }
  }
  if(group == "ser"){
    for(s in unique(newdata$ser)){
      id <- newdata$ser == s
      nb <- sum(id)
      corr <- plyr::laply(mmat, function(mat) apply(mat[,id], 1,
                                                    function(x) cor(x[1:nb], x[(nb + 1): (nb * 2)],
                                                                    method = "spearman")))
      
      # summarise
      corr_df <- plyr::adply(corr, 1, function(cor) data.frame(cor_m = mean(cor),
                                                               cor_l = quantile(cor, probs = 0.1),
                                                               cor_h = quantile(cor, probs = 0.9),
                                                               ser = s))
      corr_df$label <- c("T", "CWD", "T +\nT²", "T +\nCWD", "T +\nT² +\nCWD")
      out1 <- rbind(out1, corr_df)
      
      # get probs
      out2 <- rbind(out2, data.frame(ser = rep(s, 4),
                                     test = c("T > CWD","T+CWD > CWD",
                                              "T+CWD > T+T²",
                                              "T+CWD > T+T²+CWD"),
                                     probs = c(sum(corr[1,] > corr[2,]) / 4000,
                                               sum(corr[2,] < corr[4,]) / 4000,
                                               sum(corr[3,] < corr[4,]) / 4000,
                                               sum(corr[5,] < corr[4,]) / 4000)))
    }
  }
  
  return(list(out1, out2))
}

newdata <- d2
newdata$mqd_std <- 0
m_clim[[5]] <- m_d

cc_ser <- get_corrs(m_clim, mm2, newdata, yy_ser, group="ser")

# select ser based on optimum or not
cc2 <- inner_join(cc_ser[[1]], pp_quad_df[,c("ser", "prob")], by = "ser")
cc2$type <- ifelse(cc2$prob >= 0.9, "quadratic trend", "linear trend")
cc2$grec <- substr(cc2$ser, 1, 1)

g_c2 <- ggplot(cc2, aes(x=label, y=cor_m,  color=type)) +
  geom_point(position = position_dodge2(width=0.5), alpha=0.25) +
  stat_summary(fun.data = "mean_cl_boot", 
               position = position_dodge2(width=0.5),
               size=1.5, shape=18) +
  facet_wrap(vars(type)) +
  labs(x = "Climatic model",
       y = "Correlation coefficient between trend and climatic model") +
  guides(color="none")

ggsave("LIF/Croissance/ser_trend/figures/fig3.png", g_c2)


### Figure 4
## model estimates from the climate model
bb <- coef(m_d, summary = FALSE)
bbd <- plyr::adply(bb$grecorand, c(2, 3), function(x) data.frame(m = mean(x),
                                                                 lci = quantile(x, probs = 0.05),
                                                                 uci = quantile(x, probs = 0.95)))
bbd <- subset(bbd, X2 %in% c("pc1.1", "pc1.2",
                             "dethsus","dethsps",
                             "pc1.1:dethsus",
                             "pc1.1:dethsps"))
bbd$X1 <- factor(bbd$X1, labels = c("Grand ouest", "Centre nord",
                                    "Grand est", "Vosges", "Jura",
                                    "Sud ouest", "Massif central",
                                    "Alpes", "Pyrénées", "Méditerranée",
                                    "Corse"))
bbd$X2 <- factor(bbd$X2, labels = c("Temperature anomalies\nlinear",
                                    "Temperature anomalies\nquadratic",
                                    "Summer water\ndeficit anomalies",
                                    "Spring water\ndeficit anomalies"))
ff <- as.data.frame(fixef(m_d, probs = c(0.05, 0.95)))
ff$X2 <- rownames(ff)
names(ff)[1] <- "m"
ffs <- subset(ff, X2 %in% c("pc1.1", "pc1.2",
                            "dethsps", "dethsus"))
ffs$X2 <- factor(ffs$X2, labels = c("Spring water\ndeficit anomalies",
                                    "Summer water\ndeficit anomalies",
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

ggsave("LIF/Croissance/ser_trend/figures/fig4.png", ge)


##### Figures in SI #####

## Figure S5
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

ggsave("LIF/Croissance/ser_trend/figures/res1_mm2.png", res1)

### Figure S6
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
ggsave("LIF/Croissance/ser_trend/figures/res2_mm2.png", resa)

### Figure S7
## compute where the ser ended up at the end of the period
newd <- data.frame(ser = rep(unique(d_all$ser), each = 2),
                   year = rep(c(1982, 2021), times = 85),
                   year2.1 = rep(c(-0.06421230,  0.05267874), times = 85),
                   year2.2 = rep(c(0.06025386, 0.04951907), times = 85),
                   mqd_std = 0)
newd$grecorand <- substr(newd$ser, 1, 1)
tp <- fitted(mm2, newdata = newd, summary = FALSE)
tpp <- as.data.frame(posterior_summary(tp[,seq(2, 170, 2)] - tp[,seq(1, 170, 2)],
                                       probs = c(0.1, 0.9)))
tpp$ser <- unique(d_all$ser)
tpp$quad <- pp_quad_df$prob

gg_t <- ggplot(tpp, aes(x=forcats::fct_reorder(ser, Estimate, min),
                        y = Estimate, ymin=Q10, ymax=Q90,
                        color=quad<0.9)) +
  geom_linerange() +
  geom_point() +
  coord_flip() +
  geom_hline(yintercept = 0, linetype="dashed") +
  labs(x = "",
       y = "Predicted difference between prediction at the end and at the beginning")

ggsave("LIF/Croissance/ser_trend/figures/diff_end_beg.png", gg_t)

# new version as a map
## indicate sig level
tpp$sig <- case_when(tpp$Q10 < 0 & tpp$Q90 < 0 ~ "sig",
                     tpp$Q10 > 0 & tpp$Q90 > 0 ~ "sig",
                     TRUE ~ "not sig.")

ser_sf %>%
  inner_join(tpp, by = c("codeser" = "ser")) %>%
  mutate(a = ifelse(sig == "sig", Estimate, NA)) -> ss

gs <- ggplot(ss) +
  geom_sf(aes(fill=a)) +
  scale_fill_gradient2(name = "Changes in\nabsolute growth\nrate (%)")

ggsave("LIF/Croissance/ser_trend/figures/ser_changes.png", gs, width = 8, height = 8)


### Figure S8
## average normal temperature per SER against predicted trend form
## load the climate normal data (all parameters)
ppp <- do.call(rbind,
               lapply(list.files(path = "LIF/Croissance/ser_trend/data/",
                                 pattern = "serref", full.names = TRUE),
                      read.csv))
ppp <- subset(ppp, ser %in% pp_quad_df$ser)
ppp$greco <- substr(ppp$ser, 1, 1)
# create groups of ser
ppp$coll <- rep(ifelse(pp_quad_df$prob >= 0.9, "quad. trend",
                       ifelse(pp_lin_df$prob >= 0.9, "lin. decline", "no trend")),
                4)
## widden it
ppp %>%
  pivot_wider(names_from = type,
              values_from = value) -> ppw

gt <- ggplot(ppw, aes(x=coll, y=tmoy/10)) +
  geom_jitter(width=0.1) +
  geom_hline(yintercept = c(7.7, 12)) +
  labs(x="Trend shape",
       y="Average temperature (°C)")
ggsave("LIF/Croissance/ser_trend/figures/trend_temp_ser.png", gt)


### Figure S9
## correlation between temporal trend and climate implied temporal trend
# plot the temporal changes implied by the model
nn <- d2
nn$mqd_std <- 0
nn <- cbind(nn, fitted(m_d, newdata=nn, probs=c(0.1, 0.9)))

nn$greco <- factor(nn$grecorand, labels = c("Grand ouest", "Centre nord",
                                            "Grand est", "Vosges", "Jura",
                                            "Sud ouest", "Massif central",
                                            "Alpes", "Pyrénnés", "Méditerranée",
                                            "Corse"))


# plot the correlation between fitted trend and fitted climate
nn %>%
  select(greco, ser, year, Estimate, Q10, Q90) %>%
  rename(cm = Estimate, cl = Q10, ch=Q90) %>%
  inner_join(pp_ser[,c("ser", "year", "Estimate", "Q10", "Q90")],
             by = c("ser", "year")) -> dd

gg_tc <- ggplot(dd, aes(x=cm, xmin=cl,xmax=ch,
                        y=Estimate, ymin=Q10, ymax=Q90)) +
  geom_linerange(alpha=0.25) +
  geom_errorbarh(alpha=0.25) +
  geom_point(alpha = 0.25) +
  geom_abline(slope=1, intercept = 0) +
  facet_wrap(vars(greco), scales = "free") +
  labs(x = "Prediction climate model",
       y = "Prediction trend model")

ggsave("LIF/Croissance/ser_trend/figures/climate_trend_ref_pied.png", gg_tc)