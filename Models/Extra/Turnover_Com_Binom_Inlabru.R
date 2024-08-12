library(tidyverse)
library(INLA)
library(fmesher)
library(inlabru)
library(sf)

#### Data Setup ####
vegGeo <- read_rds("./Vegetation/vegGeo.rds")

vegINLA <- vegGeo %>% mutate(fBlock = as.factor(Block),
                             Plot_code = paste0(Plot, Block),
                             whiteStick = as.logical(whiteStick),
                             dungD = as.logical(dungD),
                             treat = factor(Plot,
                                            levels = c("III", "I", "II", "IV"))) %>% 
  filter(Year > 2002)#%>% filter(!is.na(nvcM1)) 

# anything nonzero has at least one NVC turn over 
vegINLA <- vegINLA %>% mutate(nvcTurnBin = ifelse(nvcTurn != 1, 0, 1))

# first, the mesh
mesh_TurnComBin <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_TurnComBin) + theme_bw()
#mesh_TurnComBin$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_TurnComBin <- inla.spde2.matern(mesh = mesh_TurnComBin)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_TurnComBin <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)

#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_TurnComBin <- ~ #Intercept(1) + 
  # interaction terms don't change sheepOff but kill the signal in cowOff due to construction
  sheepOff + cowOff + 
  #treat +
  # humidity, precipitation, and pressure are too correlated to include in a single model
  # simply by effect size and significance, pressure is best; harder to interpret though
  # frost days have no effect
  windSpeed + meanT + solarRad + precip + 
  # dung, bareGround, and wood have almost no observations
  # water, rock, and mud also fairly low
  # flowers and tussock take away from the NVC community variables (less positive)
  litter + rock + #graze +
  # flower and tussock are mainly M and U
  elev +
  isW + isH + isU + isM + 
  dungD + whiteStick

inlaTempC_TurnComBin <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  dungD + whiteStick + 
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

inlaSpatTempAR1C_TurnComBin <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  dungD + whiteStick + 
  field(geometry, 
        model = aMatSPDE_TurnComBin, 
        group = Year, 
        group_mapper = yearMap_TurnComBin,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

# inlaSpatTempAR2C_TurnComBin <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_TurnComBin, 
#         group = Year, 
#         group_mapper = yearMap_TurnComBin,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 2))
# 
# inlaSpatTempAR3C_TurnComBin <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_TurnComBin, 
#         group = Year, 
#         group_mapper = yearMap_TurnComBin,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_TurnComBin <- like(
  formula = nvcTurnBin ~ .,
  # full distribution list: inla.list.models()
  family = "binomial",
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_TurnComBin)
)

likeTemp_TurnComBin <- like(
  formula = nvcTurnBin ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "binomial", 
  #beta.censor.value = 0.05,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_TurnComBin)
)

likeSpatTempAR_TurnComBin <- like(
  formula = nvcTurnBin ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "binomial", 
  #beta.truncation = 0.05,
  data=vegINLA,
  domain = list(geometry = mesh_TurnComBin)
)

# run the actual fit
# null model
fitNull_TurnComBin <- bru(components = inlaNullC_TurnComBin,
                           likeNull_TurnComBin,
                           options = list(
                             # add linear combinations
                             lincomb = c(lc1, lc2),
                             #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                             control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                             control.predictor = list(compute=TRUE),
                             #control.inla = list(int.strategy = "eb"),
                             verbose = FALSE))

# temporal model
fitTemp_TurnComBin <- bru(components = inlaTempC_TurnComBin,
                           likeTemp_TurnComBin,
                           options = list(
                             # add linear combinations
                             lincomb = c(lc1, lc2),
                             #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                             control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                             control.predictor = list(compute=TRUE),
                             #control.inla = list(int.strategy = "eb"),
                             verbose = FALSE))

# spatiotemporal model
# fairly memory-intensive
fitSpatTempAR1_TurnComBin <- bru(components = inlaSpatTempAR1C_TurnComBin,
                                  likeSpatTempAR_TurnComBin,
                                  options = list(
                                    # add linear combinations
                                    lincomb = c(lc1, lc2),
                                    #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                    control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                    control.predictor = list(compute=TRUE),
                                    #control.inla = list(int.strategy = "eb"),
                                    verbose = FALSE))

# AR2/3 untested; likely do not converge


#### Model Output ####
# WAIC, DIC, CPO
# CPO needs no NAs to work perfectly; na.rm is good otherwise
cbind(fitNull_TurnComBin$waic$waic[1], fitNull_TurnComBin$dic$dic[1], -sum(log(fitNull_TurnComBin$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_TurnComBin$waic$waic[1], fitTemp_TurnComBin$dic$dic[1], -sum(log(fitTemp_TurnComBin$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_TurnComBin$waic$waic[1], fitSpatTempAR1_TurnComBin$dic$dic[1], -sum(log(fitSpatTempAR1_TurnComBin$cpo$cpo), na.rm = TRUE))
#[1,] 35045.41 35044.69 17522.7
#[1,] 34330.06 34331.62 17165.06
#[1,] 30417.98 30836.95 15240.68

# dung and white stick
#[1,] 323375.6 323306.6 161729.8
#[1,] 323365.6 323300.6 161725.3

# no benefits to using AR2 for height
# cbind(fitSpatTempAR2_TurnComBin$waic$waic[1], fitSpatTempAR2_TurnComBin$dic$dic[1], -sum(log(fitSpatTempAR2_TurnComBin$cpo$cpo), na.rm = TRUE))
# cbind(fitSpatTempAR3_TurnComBin$waic$waic[1], fitSpatTempAR3_TurnComBin$dic$dic[1], -sum(log(fitSpatTempAR3_TurnComBin$cpo$cpo), na.rm = TRUE))


# summary
summary(fitSpatTempAR1_TurnComBin)

#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_TurnComBin)
brinla::bri.fixed.plot(fitTemp_TurnComBin)
brinla::bri.fixed.plot(fitSpatTempAR1_TurnComBin)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTempAR1_TurnComBin$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTempAR1_TurnComBin$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
inla.pmarginal(q = 0,
                   marginal= fitSpatTempAR1_TurnComBin$marginals.fixed$sheepOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_TurnComBin$marginals.fixed$cowOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_TurnComBin$marginals.lincomb.derived$offtake)


#### Marginal distribution ####

plots_TurnComBin <- list()

# Null model
inlaNullPvalT_TurnComBin<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_TurnComBin[i]<-inla.pmarginal(q=vegINLA$nvcTurnBin[i],
                                               marginal=fitNull_TurnComBin$marginals.fitted.values[[i]])
}

postProb1_TurnComBin <- data.frame(predicted = fitNull_TurnComBin$summary.fitted.values$mean[1:length(vegINLA$nvcTurnBin)],
                                    observed = vegINLA$nvcTurnBin,
                                    lower = fitNull_TurnComBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnBin)],
                                    upper = fitNull_TurnComBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnBin)],
                                    p.value = inlaNullPvalT_TurnComBin) 

min1_TurnComBin <- min(postProb1_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)
max1_TurnComBin <- max(postProb1_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnComBin[[1]] <- ggplot(postProb1_TurnComBin, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnComBin[[2]] <- ggplot(postProb1_TurnComBin, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_TurnComBin, max1_TurnComBin), y = c(min1_TurnComBin, max1_TurnComBin)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_TurnComBin<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_TurnComBin[i]<-inla.pmarginal(q=vegINLA$nvcTurnBin[i],
                                               # lots of infinity values, so this will not work
                                               marginal=fitTemp_TurnComBin$marginals.fitted.values[[i]])
}
postProb2_TurnComBin <- data.frame(predicted = fitTemp_TurnComBin$summary.fitted.values$mean[1:length(vegINLA$nvcTurnBin)],
                                    observed = vegINLA$nvcTurnBin,
                                    lower = fitTemp_TurnComBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnBin)],
                                    upper = fitTemp_TurnComBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnBin)],
                                    p.value = inlaTempPvalT_TurnComBin) 

min2_TurnComBin <- min(postProb2_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)
max2_TurnComBin <- max(postProb2_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnComBin[[3]] <- ggplot(postProb2_TurnComBin, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnComBin[[4]] <- ggplot(postProb2_TurnComBin, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_TurnComBin, max2_TurnComBin), y = c(min2_TurnComBin, max2_TurnComBin)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_TurnComBin<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_TurnComBin[i]<-inla.pmarginal(q=vegINLA$nvcTurnBin[i],
                                                      marginal=fitSpatTempAR1_TurnComBin$marginals.fitted.values[[i]])
}

postProb3_TurnComBin <- data.frame(predicted = fitSpatTempAR1_TurnComBin$summary.fitted.values$mean[1:length(vegINLA$nvcTurnBin)],
                                    observed = vegINLA$nvcTurnBin,
                                    lower = fitSpatTempAR1_TurnComBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnBin)],
                                    upper = fitSpatTempAR1_TurnComBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnBin)],
                                    p.value = inlaSpatTempAR1PvalT_TurnComBin) 

min3_TurnComBin <- min(postProb3_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)
max3_TurnComBin <- max(postProb3_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnComBin[[5]] <- ggplot(postProb3_TurnComBin, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnComBin[[6]] <- ggplot(postProb3_TurnComBin, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_TurnComBin, max3_TurnComBin), y = c(min3_TurnComBin, max3_TurnComBin)) + 
  theme_bw()
 
cowplot::plot_grid(plotlist = plots_TurnComBin, nrow = 3)

library(pROC)
rocobj <- roc(postProb3_TurnComBin$observed, postProb3_TurnComBin$predicted)
auc <- round(auc(postProb3_TurnComBin$observed, postProb3_TurnComBin$predicted),4)
ggroc(rocobj, colour = 'steelblue', size = 2) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()


# inlaSpatTempAR2PvalT_TurnComBin<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR2PvalT_TurnComBin[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR2_TurnComBin$marginals.fitted.values[[i]])
# }
# 
# postProb4_TurnComBin <- data.frame(predicted = fitSpatTempAR2_TurnComBin$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR2_TurnComBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR2_TurnComBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR2PvalT_TurnComBin) 
# 
# min4_TurnComBin <- min(postProb4_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)
# max4_TurnComBin <- max(postProb4_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_TurnComBin[[7]] <- ggplot(postProb4_TurnComBin, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_TurnComBin[[8]] <- ggplot(postProb4_TurnComBin, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min4_TurnComBin, max4_TurnComBin), y = c(min4_TurnComBin, max4_TurnComBin)) + 
#   theme_bw()
# 
# inlaSpatTempAR3PvalT_TurnComBin<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR3PvalT_TurnComBin[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR3_TurnComBin$marginals.fitted.values[[i]])
# }
# 
# postProb5_TurnComBin <- data.frame(predicted = fitSpatTempAR3_TurnComBin$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR3_TurnComBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR3_TurnComBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR3PvalT_TurnComBin) 
# 
# min5_TurnComBin <- min(postProb5_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)
# max5_TurnComBin <- max(postProb5_TurnComBin[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_TurnComBin[[9]] <- ggplot(postProb5_TurnComBin, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_TurnComBin[[10]] <- ggplot(postProb5_TurnComBin, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min5_TurnComBin, max5_TurnComBin), y = c(min5_TurnComBin, max5_TurnComBin)) + 
#   theme_bw()
# 
# cowplot::plot_grid(plotlist = plots_TurnComBin, nrow = 6)



#### Plotting ####
# currently non-functional

#### Exports ####
write_rds(fitNull_TurnComBin, "./Models/Objects/fitNull_TurnComBin.rds")
write_rds(fitTemp_TurnComBin, "./Models/Objects/fitTemp_TurnComBin.rds")
write_rds(fitSpatTempAR1_TurnComBin, "./Models/Objects/fitSpatTempAR1_TurnComBin.rds")
#write_rds(fitSpatTempAR2_TurnComBin, "./Models/Objects/fitSpatTempAR2_TurnComBin.rds")
#write_rds(fitSpatTempAR3_TurnComBin, "./Models/Objects/fitSpatTempAR3_TurnComBin.rds")
