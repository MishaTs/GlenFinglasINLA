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
squeeze_TurnBeta <- (1 - sort(unique(vegINLA$nvcTurnSub))[n_distinct(vegINLA$nvcTurnSub) - 2])/10
vegINLA <- vegINLA %>% mutate(nvcTurnSubSq = ifelse(nvcTurnSub == 1, 1 - squeeze_TurnBeta, 
                                                     ifelse(nvcTurnSub == 0, squeeze_TurnBeta,
                                                            nvcTurnSub)))

# t(apply(table(vegINLA$isW, vegINLA$flower),1, function(x) x/sum(x)))

# first, the mesh
mesh_TurnBeta <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_TurnBeta) + theme_bw()
#mesh_TurnBeta$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_TurnBeta <- inla.spde2.matern(mesh = mesh_TurnBeta)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_TurnBeta <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)

#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_TurnBeta <- ~ #Intercept(1) + 
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

inlaTempC_TurnBeta <- ~ #Intercept(1) + 
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

inlaSpatTempAR1C_TurnBeta <- ~ #Intercept(1) + 
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

# inlaSpatTempAR2C_TurnBeta <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_TurnBeta, 
#         group = Year, 
#         group_mapper = yearMap_TurnBeta,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 2))
# 
# inlaSpatTempAR3C_TurnBeta <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_TurnBeta, 
#         group = Year, 
#         group_mapper = yearMap_TurnBeta,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_TurnBeta <- like(
  formula = nvcTurnSubSq ~ .,
  # full distribution list: inla.list.models()
  family = "beta",
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_TurnBeta)
)

likeTemp_TurnBeta <- like(
  formula = nvcTurnSubSq ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "beta", 
  #beta.censor.value = 0.05,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_TurnBeta)
)

likeSpatTempAR_TurnBeta <- like(
  formula = nvcTurnSubSq ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "beta", 
  #beta.truncation = 0.05,
  data=vegINLA,
  domain = list(geometry = mesh_TurnBeta)
)

# run the actual fit
# null model
# does not accept factor variables
fitNull_TurnBeta <- bru(components = inlaNullC_TurnBeta,
                       likeNull_TurnBeta,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# temporal model
fitTemp_TurnBeta <- bru(components = inlaTempC_TurnBeta,
                       likeTemp_TurnBeta,
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
fitSpatTempAR1_TurnBeta <- bru(components = inlaSpatTempAR1C_TurnBeta,
                              likeSpatTempAR_TurnBeta,
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
cbind(fitNull_TurnBeta$waic$waic[1], fitNull_TurnBeta$dic$dic[1], -sum(log(fitNull_TurnBeta$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_TurnBeta$waic$waic[1], fitTemp_TurnBeta$dic$dic[1], -sum(log(fitTemp_TurnBeta$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_TurnBeta$waic$waic[1], fitSpatTempAR1_TurnBeta$dic$dic[1], -sum(log(fitSpatTempAR1_TurnBeta$cpo$cpo), na.rm = TRUE))
# with everything: dungD, whiteStick, nvcEffort
#[1,] -240638.3 -240630.9 -120319.2
#[1,] -241972.3 -241816.1 -120986.1
#[1,] -243478.9 -243241.2 -121739.4

# no dungD
#[1,] -240639 -240631.8 -120319.5
#[1,] -241971.9 -241818.2 -120985.9
#[1,] -243482.1 -243240.7 -121741

# no dungD, no nvcEffort
#[1,] -240634.7 -240628 -120317.3
#[1,] -241976.1 -241822.1 -120988
#[1,] -243473.1 -243240.6 -121736.5

# no dungD, no litter
# no improvement in precip coefficients either
#[1,] -240630.1 -240623.4 -120315.1
#[1,] -241969.7 -241808.7 -120984.8
#[1,] -243457.5 -243232.8 -121728.7

# no dungD, yes graze
#[1,] -240639.2 -240631.7 -120319.6
#[1,] -241976.1 -241820.2 -120988.1
#[1,] -243454.7 -243230.3 -121727.3

# no dungD, yes sameObserver
#[1,] -240912.9 -240905.4 -120456.5
#[1,] -241988.2 -241846.6 -120994.1
#[1,] -243472.8 -243248.3 -121736.3

# no nvcEffort, dungD, yes sameObserver
#[1,] -240909.7 -240902.6 -120454.9
#[1,] -241983 -241845.6 -120991.5
#[1,] -243480.4 -243251.6 -121740.1


# summary
summary(fitSpatTempAR1_TurnBeta)

#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, the intercept is meaningless and they should be evaluated for "relative effect"
brinla::bri.fixed.plot(fitNull_TurnBeta)
brinla::bri.fixed.plot(fitTemp_TurnBeta)
brinla::bri.fixed.plot(fitSpatTempAR1_TurnBeta)
# get custom quantiles
# this tests the 0 level for a variable
1 - inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_TurnBeta$marginals.fixed$sheepOff)
1 - inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_TurnBeta$marginals.fixed$cowOff)
1 - inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_TurnBeta$marginals.lincomb.derived$offtake)


#### Marginal distribution ####

plots_TurnBeta <- list()

# Null model
inlaNullPvalT_TurnBeta<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_TurnBeta[i]<-inla.pmarginal(q=vegINLA$nvcTurnSubSq[i],
                                           marginal=fitNull_TurnBeta$marginals.fitted.values[[i]])
}

postProb1_TurnBeta <- data.frame(predicted = fitNull_TurnBeta$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSubSq)],
                                observed = vegINLA$nvcTurnSubSq,
                                lower = fitNull_TurnBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSubSq)],
                                upper = fitNull_TurnBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSubSq)],
                                p.value = inlaNullPvalT_TurnBeta) 

min1_TurnBeta <- min(postProb1_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)
max1_TurnBeta <- max(postProb1_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnBeta[[1]] <- ggplot(postProb1_TurnBeta, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnBeta[[2]] <- ggplot(postProb1_TurnBeta, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_TurnBeta, max1_TurnBeta), y = c(min1_TurnBeta, max1_TurnBeta)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_TurnBeta<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_TurnBeta[i]<-inla.pmarginal(q=vegINLA$nvcTurnSubSq[i],
                                           # lots of infinity values, so this will not work
                                           marginal=fitTemp_TurnBeta$marginals.fitted.values[[i]])
}
postProb2_TurnBeta <- data.frame(predicted = fitTemp_TurnBeta$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSubSq)],
                                observed = vegINLA$nvcTurnSubSq,
                                lower = fitTemp_TurnBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSubSq)],
                                upper = fitTemp_TurnBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSubSq)],
                                p.value = inlaTempPvalT_TurnBeta) 

min2_TurnBeta <- min(postProb2_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)
max2_TurnBeta <- max(postProb2_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnBeta[[3]] <- ggplot(postProb2_TurnBeta, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnBeta[[4]] <- ggplot(postProb2_TurnBeta, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_TurnBeta, max2_TurnBeta), y = c(min2_TurnBeta, max2_TurnBeta)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_TurnBeta<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_TurnBeta[i]<-inla.pmarginal(q=vegINLA$nvcTurnSubSq[i],
                                                  marginal=fitSpatTempAR1_TurnBeta$marginals.fitted.values[[i]])
}

postProb3_TurnBeta <- data.frame(predicted = fitSpatTempAR1_TurnBeta$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSubSq)],
                                observed = vegINLA$nvcTurnSubSq,
                                lower = fitSpatTempAR1_TurnBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSubSq)],
                                upper = fitSpatTempAR1_TurnBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSubSq)],
                                p.value = inlaSpatTempAR1PvalT_TurnBeta) 

min3_TurnBeta <- min(postProb3_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)
max3_TurnBeta <- max(postProb3_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnBeta[[5]] <- ggplot(postProb3_TurnBeta, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnBeta[[6]] <- ggplot(postProb3_TurnBeta, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_TurnBeta, max3_TurnBeta), y = c(min3_TurnBeta, max3_TurnBeta)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_TurnBeta, nrow = 3)

# inlaSpatTempAR2PvalT_TurnBeta<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR2PvalT_TurnBeta[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR2_TurnBeta$marginals.fitted.values[[i]])
# }
# 
# postProb4_TurnBeta <- data.frame(predicted = fitSpatTempAR2_TurnBeta$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR2_TurnBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR2_TurnBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR2PvalT_TurnBeta) 
# 
# min4_TurnBeta <- min(postProb4_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)
# max4_TurnBeta <- max(postProb4_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_TurnBeta[[7]] <- ggplot(postProb4_TurnBeta, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_TurnBeta[[8]] <- ggplot(postProb4_TurnBeta, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min4_TurnBeta, max4_TurnBeta), y = c(min4_TurnBeta, max4_TurnBeta)) + 
#   theme_bw()
# 
# inlaSpatTempAR3PvalT_TurnBeta<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR3PvalT_TurnBeta[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR3_TurnBeta$marginals.fitted.values[[i]])
# }
# 
# postProb5_TurnBeta <- data.frame(predicted = fitSpatTempAR3_TurnBeta$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR3_TurnBeta$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR3_TurnBeta$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR3PvalT_TurnBeta) 
# 
# min5_TurnBeta <- min(postProb5_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)
# max5_TurnBeta <- max(postProb5_TurnBeta[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_TurnBeta[[9]] <- ggplot(postProb5_TurnBeta, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_TurnBeta[[10]] <- ggplot(postProb5_TurnBeta, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min5_TurnBeta, max5_TurnBeta), y = c(min5_TurnBeta, max5_TurnBeta)) + 
#   theme_bw()
# 
# cowplot::plot_grid(plotlist = plots_TurnBeta, nrow = 6)



#### Plotting ####
# currently non-functional

#### Exports ####
write_rds(fitNull_TurnBeta, "./Models/Objects/fitNull_TurnBeta.rds")
write_rds(fitTemp_TurnBeta, "./Models/Objects/fitTemp_TurnBeta.rds")
write_rds(fitSpatTempAR1_TurnBeta, "./Models/Objects/fitSpatTempAR1_TurnBeta.rds")
#write_rds(fitSpatTempAR2_TurnBeta, "./Models/Objects/fitSpatTempAR2_TurnBeta.rds")
#write_rds(fitSpatTempAR3_TurnBeta, "./Models/Objects/fitSpatTempAR3_TurnBeta.rds")
