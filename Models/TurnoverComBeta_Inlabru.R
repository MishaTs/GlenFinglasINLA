library(tidyverse)
library(INLA)
library(fmesher)
library(inlabru)
library(sf)

#### Data Setup ####
vegGeo <- read_rds("./Vegetation/vegGeo.rds")

vegINLA <- vegGeo %>% mutate(fBlock = as.factor(Block),
                             Plot_code = paste0(Plot, Block),
                             # as.numeric converts 0, 1 to 1, 2 for some reason
                             whiteStick = as.numeric(whiteStick) - 1,
                             dungD = as.numeric(dungD),
                             treat = factor(Plot,
                                            levels = c("III", "I", "II", "IV")),
                             #nvcEffort = nvcVars - nvcSumComs,
                             sameObserver = as.numeric(Recorder == lag(Recorder))) %>% 
  filter(Year > 2002)#%>% filter(!is.na(nvcM1)) 

# generate a "squeeze" term to eliminate discrete 0s and 1s from the data
# we want this to be a consistent distance from both 0 and 1 despite much more granularity at upper intervals rather than lower ones
squeeze_TurnComBeta <- (1 - sort(unique(vegINLA$nvcTurn))[n_distinct(vegINLA$nvcTurn) - 2])/10
vegINLA <- vegINLA %>% mutate(nvcTurnSq = ifelse(nvcTurn == 1, 1 - squeeze_TurnComBeta,
                                                 ifelse(nvcTurn == 0, squeeze_TurnComBeta,
                                                        nvcTurn)))

# t(apply(table(vegINLA$isW, vegINLA$flower),1, function(x) x/sum(x)))

# first, the mesh
mesh_TurnComBeta <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_TurnComBeta) + theme_bw()
#mesh_TurnComBeta$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_TurnComBeta <- inla.spde2.matern(mesh = mesh_TurnComBeta)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_TurnComBeta <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)

#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_TurnComBeta <- ~ #Intercept(1) + 
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
  litter + #rock + 
  #graze +
  # flower and tussock are mainly M and U
  elev +
  isW + isH + isU + isM +
  #dungD + 
  whiteStick +
  #nvcEffort + 
  sameObserver

inlaTempC_TurnComBeta <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + #rock + 
  #graze +
  elev +
  isW + isH + isU + isM +
  #dungD + 
  whiteStick +
  #nvcEffort +
  sameObserver + 
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

inlaSpatTempAR1C_TurnComBeta <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + #rock + 
  #graze +
  elev +
  isW + isH + isU + isM +
  #dungD + 
  whiteStick + 
  #nvcEffort +
  sameObserver +
  field(geometry, 
        model = aMatSPDE_TurnBeta, 
        group = Year, 
        group_mapper = yearMap_TurnBeta,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

# hidden alternative AR term formulae since models did not converge or were not used
# inlaSpatTempAR2C_TurnComBeta <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_TurnComBeta, 
#         group = Year, 
#         group_mapper = yearMap_TurnComBeta,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 2))
# 
# inlaSpatTempAR3C_TurnComBeta <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_TurnComBeta, 
#         group = Year, 
#         group_mapper = yearMap_TurnComBeta,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_TurnComBeta <- like(
  formula = nvcTurnSq ~ .,
  # full distribution list: inla.list.models()
  family = "beta",
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_TurnComBeta)
)

likeTemp_TurnComBeta <- like(
  formula = nvcTurnSq ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "beta", 
  #beta.censor.value = 0.05,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_TurnComBeta)
)

likeSpatTempAR_TurnComBeta <- like(
  formula = nvcTurnSq ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "beta", 
  #beta.truncation = 0.05,
  data=vegINLA,
  domain = list(geometry = mesh_TurnComBeta)
)

# run the actual fit
# null model
# does not accept factor variables
# Error in ibm_jacobian.bru_mapper_linear(mapper[["mappers"]][[x]], input = input[[x]],  : 
#   The input to a bru_mapper_linear evaluation must be numeric or logical.
fitNull_TurnComBeta <- bru(components = inlaNullC_TurnComBeta,
                        likeNull_TurnComBeta,
                        options = list(
                          # add linear combinations
                          lincomb = c(lc1, lc2),
                          #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                          control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                          control.predictor = list(compute=TRUE),
                          #control.inla = list(int.strategy = "eb"),
                          verbose = FALSE))

# temporal model
fitTemp_TurnComBeta <- bru(components = inlaTempC_TurnComBeta,
                        likeTemp_TurnComBeta,
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
fitSpatTempAR1_TurnComBeta <- bru(components = inlaSpatTempAR1C_TurnComBeta,
                               likeSpatTempAR_TurnComBeta,
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
cbind(fitNull_TurnComBeta$waic$waic[1], fitNull_TurnComBeta$dic$dic[1], -sum(log(fitNull_TurnComBeta$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_TurnComBeta$waic$waic[1], fitTemp_TurnComBeta$dic$dic[1], -sum(log(fitTemp_TurnComBeta$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_TurnComBeta$waic$waic[1], fitSpatTempAR1_TurnComBeta$dic$dic[1], -sum(log(fitSpatTempAR1_TurnComBeta$cpo$cpo), na.rm = TRUE))
#[1,] -80359.34 -80352.17 -40179.67
#[1,] -80555.73 -80497.38 -40277.87
#[1,] -81639.26 -81419.14 -40819.61

# no nvcEffort
#[1,] -80360.39 -80354.02 -40180.19
#[1,] -80550.86 -80497.04 -40275.43
#[1,] -81643.02 -81434.54 -40821.49

# no benefits to using AR2 
# cbind(fitSpatTempAR2_TurnComBeta$waic$waic[1], fitSpatTempAR2_TurnComBeta$dic$dic[1], -sum(log(fitSpatTempAR2_TurnComBeta$cpo$cpo), na.rm = TRUE))
# cbind(fitSpatTempAR3_TurnComBeta$waic$waic[1], fitSpatTempAR3_TurnComBeta$dic$dic[1], -sum(log(fitSpatTempAR3_TurnComBeta$cpo$cpo), na.rm = TRUE))


# summary
summary(fitSpatTempAR1_TurnComBeta)

#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_TurnComBeta)
brinla::bri.fixed.plot(fitTemp_TurnComBeta)
brinla::bri.fixed.plot(fitSpatTempAR1_TurnComBeta)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTempAR1_TurnComBeta$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTempAR1_TurnComBeta$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
1 - inla.pmarginal(q = 0,
                   marginal= fitSpatTempAR1_TurnComBeta$marginals.fixed$sheepOff)
inla.pmarginal(q = 0,
                   marginal= fitSpatTempAR1_TurnComBeta$marginals.fixed$cowOff)
inla.pmarginal(q = 0,
                   marginal= fitSpatTempAR1_TurnComBeta$marginals.lincomb.derived$offtake)


#### Marginal distribution ####
plots_TurnComBeta <- list()

# Null model
inlaNullPvalT_TurnComBeta<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_TurnComBeta[i]<-inla.pmarginal(q=vegINLA$nvcTurnSq[i],
                                            marginal=fitNull_TurnComBeta$marginals.fitted.values[[i]])
}

postProb1_TurnComBeta <- data.frame(predicted = fitNull_TurnComBeta$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSq)],
                                 observed = vegINLA$nvcTurnSq,
                                 lower = fitNull_TurnComBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSq)],
                                 upper = fitNull_TurnComBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSq)],
                                 p.value = inlaNullPvalT_TurnComBeta) 

min1_TurnComBeta <- min(postProb1_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)
max1_TurnComBeta <- max(postProb1_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnComBeta[[1]] <- ggplot(postProb1_TurnComBeta, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnComBeta[[2]] <- ggplot(postProb1_TurnComBeta, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_TurnComBeta, max1_TurnComBeta), y = c(min1_TurnComBeta, max1_TurnComBeta)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_TurnComBeta<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_TurnComBeta[i]<-inla.pmarginal(q=vegINLA$nvcTurnSq[i],
                                            marginal=fitTemp_TurnComBeta$marginals.fitted.values[[i]])
}
postProb2_TurnComBeta <- data.frame(predicted = fitTemp_TurnComBeta$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSq)],
                                 observed = vegINLA$nvcTurnSq,
                                 lower = fitTemp_TurnComBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSq)],
                                 upper = fitTemp_TurnComBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSq)],
                                 p.value = inlaTempPvalT_TurnComBeta) 

min2_TurnComBeta <- min(postProb2_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)
max2_TurnComBeta <- max(postProb2_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnComBeta[[3]] <- ggplot(postProb2_TurnComBeta, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnComBeta[[4]] <- ggplot(postProb2_TurnComBeta, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_TurnComBeta, max2_TurnComBeta), y = c(min2_TurnComBeta, max2_TurnComBeta)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_TurnComBeta<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_TurnComBeta[i]<-inla.pmarginal(q=vegINLA$nvcTurnSq[i],
                                                   marginal=fitSpatTempAR1_TurnComBeta$marginals.fitted.values[[i]])
}

postProb3_TurnComBeta <- data.frame(predicted = fitSpatTempAR1_TurnComBeta$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSq)],
                                 observed = vegINLA$nvcTurnSq,
                                 lower = fitSpatTempAR1_TurnComBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSq)],
                                 upper = fitSpatTempAR1_TurnComBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSq)],
                                 p.value = inlaSpatTempAR1PvalT_TurnComBeta) 

min3_TurnComBeta <- min(postProb3_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)
max3_TurnComBeta <- max(postProb3_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnComBeta[[5]] <- ggplot(postProb3_TurnComBeta, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnComBeta[[6]] <- ggplot(postProb3_TurnComBeta, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_TurnComBeta, max3_TurnComBeta), y = c(min3_TurnComBeta, max3_TurnComBeta)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_TurnComBeta, nrow = 3)

# hidden alternative AR term checks since models did not converge or were not used
# inlaSpatTempAR2PvalT_TurnComBeta<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR2PvalT_TurnComBeta[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR2_TurnComBeta$marginals.fitted.values[[i]])
# }
# 
# postProb4_TurnComBeta <- data.frame(predicted = fitSpatTempAR2_TurnComBeta$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR2_TurnComBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR2_TurnComBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR2PvalT_TurnComBeta) 
# 
# min4_TurnComBeta <- min(postProb4_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)
# max4_TurnComBeta <- max(postProb4_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_TurnComBeta[[7]] <- ggplot(postProb4_TurnComBeta, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_TurnComBeta[[8]] <- ggplot(postProb4_TurnComBeta, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min4_TurnComBeta, max4_TurnComBeta), y = c(min4_TurnComBeta, max4_TurnComBeta)) + 
#   theme_bw()
# 
# inlaSpatTempAR3PvalT_TurnComBeta<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR3PvalT_TurnComBeta[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR3_TurnComBeta$marginals.fitted.values[[i]])
# }
# 
# postProb5_TurnComBeta <- data.frame(predicted = fitSpatTempAR3_TurnComBeta$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR3_TurnComBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR3_TurnComBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR3PvalT_TurnComBeta) 
# 
# min5_TurnComBeta <- min(postProb5_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)
# max5_TurnComBeta <- max(postProb5_TurnComBeta[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_TurnComBeta[[9]] <- ggplot(postProb5_TurnComBeta, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_TurnComBeta[[10]] <- ggplot(postProb5_TurnComBeta, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min5_TurnComBeta, max5_TurnComBeta), y = c(min5_TurnComBeta, max5_TurnComBeta)) + 
#   theme_bw()
# 
# cowplot::plot_grid(plotlist = plots_TurnComBeta, nrow = 6)

#### Exports ####
write_rds(fitNull_TurnComBeta, "./Models/Objects/fitNull_TurnComBeta.rds")
write_rds(fitTemp_TurnComBeta, "./Models/Objects/fitTemp_TurnComBeta.rds")
write_rds(fitSpatTempAR1_TurnComBeta, "./Models/Objects/fitSpatTempAR1_TurnComBeta.rds")
#write_rds(fitSpatTempAR2_TurnComBeta, "./Models/Objects/fitSpatTempAR2_TurnComBeta.rds")
#write_rds(fitSpatTempAR3_TurnComBeta, "./Models/Objects/fitSpatTempAR3_TurnComBeta.rds")
