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
vegINLA <- vegINLA %>% mutate(nvcTurnSubBin = ifelse(nvcTurnSub != 1, 0, 1))

# t(apply(table(vegINLA$isW, vegINLA$flower),1, function(x) x/sum(x)))

# first, the mesh
mesh_TurnBin <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_TurnBin) + theme_bw()
#mesh_TurnBin$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_TurnBin <- inla.spde2.matern(mesh = mesh_TurnBin)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_TurnBin <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)

#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_TurnBin <- ~ #Intercept(1) + 
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
  litter + #rock + #graze +
  # flower and tussock are mainly M and U
  elev +
  isW + isH + isU + isM# +
#dungD + whiteStick

inlaTempC_TurnBin <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + #rock + #graze +
  elev +
  isW + isH + isU + isM +
  #dungD + whiteStick + 
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

inlaSpatTempAR1C_TurnBin <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + #rock + #graze +
  elev +
  isW + isH + isU + isM +
  #dungD + whiteStick + 
  field(geometry, 
        model = aMatSPDE_TurnBin, 
        group = Year, 
        group_mapper = yearMap_TurnBin,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

# inlaSpatTempAR2C_TurnBin <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_TurnBin, 
#         group = Year, 
#         group_mapper = yearMap_TurnBin,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 2))
# 
# inlaSpatTempAR3C_TurnBin <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_TurnBin, 
#         group = Year, 
#         group_mapper = yearMap_TurnBin,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_TurnBin <- like(
  formula = nvcTurnSubBin ~ .,
  # full distribution list: inla.list.models()
  family = "binomial",
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_TurnBin)
)

likeTemp_TurnBin <- like(
  formula = nvcTurnSubBin ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "binomial", 
  #beta.censor.value = 0.05,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_TurnBin)
)

likeSpatTempAR_TurnBin <- like(
  formula = nvcTurnSubBin ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "binomial", 
  #beta.truncation = 0.05,
  data=vegINLA,
  domain = list(geometry = mesh_TurnBin)
)

# run the actual fit
# null model
# does not accept factor variables
# Error in ibm_jacobian.bru_mapper_linear(mapper[["mappers"]][[x]], input = input[[x]],  : 
#   The input to a bru_mapper_linear evaluation must be numeric or logical.
fitNull_TurnBin <- bru(components = inlaNullC_TurnBin,
                    likeNull_TurnBin,
                    options = list(
                      # add linear combinations
                      lincomb = c(lc1, lc2),
                      #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                      control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                      control.predictor = list(compute=TRUE),
                      #control.inla = list(int.strategy = "eb"),
                      verbose = FALSE))

# temporal model
# does not converge
fitTemp_TurnBin <- bru(components = inlaTempC_TurnBin,
                    likeTemp_TurnBin,
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
# does not converge
fitSpatTempAR1_TurnBin <- bru(components = inlaSpatTempAR1C_TurnBin,
                           likeSpatTempAR_TurnBin,
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
cbind(fitNull_TurnBin$waic$waic[1], fitNull_TurnBin$dic$dic[1], -sum(log(fitNull_TurnBin$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_TurnBin$waic$waic[1], fitTemp_TurnBin$dic$dic[1], -sum(log(fitTemp_TurnBin$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_TurnBin$waic$waic[1], fitSpatTempAR1_TurnBin$dic$dic[1], -sum(log(fitSpatTempAR1_TurnBin$cpo$cpo), na.rm = TRUE))
#[1,] 46455.65 46455.67 23227.83
#[1,] 43624.23 43630.55 21812.34
#[1,] 40575.79 40647.1 20295.05

# dung and white stick
#[1,] 323375.6 323306.6 161729.8
#[1,] 323365.6 323300.6 161725.3

# no benefits to using AR2 for height
# cbind(fitSpatTempAR2_TurnBin$waic$waic[1], fitSpatTempAR2_TurnBin$dic$dic[1], -sum(log(fitSpatTempAR2_TurnBin$cpo$cpo), na.rm = TRUE))
# cbind(fitSpatTempAR3_TurnBin$waic$waic[1], fitSpatTempAR3_TurnBin$dic$dic[1], -sum(log(fitSpatTempAR3_TurnBin$cpo$cpo), na.rm = TRUE))


# summary
summary(fitSpatTempAR1_TurnBin)

#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_TurnBin)
brinla::bri.fixed.plot(fitTemp_TurnBin)
brinla::bri.fixed.plot(fitSpatTempAR1_TurnBin)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTempAR1_TurnBin$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTempAR1_TurnBin$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
1 - inla.pmarginal(q = 0,
                   marginal= fitSpatTempAR1_TurnBin$marginals.fixed$sheepOff)
1 - inla.pmarginal(q = 0,
                   marginal= fitSpatTempAR1_TurnBin$marginals.fixed$cowOff)
1 - inla.pmarginal(q = 0,
                   marginal= fitSpatTempAR1_TurnBin$marginals.lincomb.derived$offtake)


#### Marginal distribution ####

plots_TurnBin <- list()

# Null model
inlaNullPvalT_TurnBin<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_TurnBin[i]<-inla.pmarginal(q=vegINLA$nvcTurnSubBin[i],
                                        marginal=fitNull_TurnBin$marginals.fitted.values[[i]])
}

postProb1_TurnBin <- data.frame(predicted = fitNull_TurnBin$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSubBin)],
                             observed = vegINLA$nvcTurnSubBin,
                             lower = fitNull_TurnBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSubBin)],
                             upper = fitNull_TurnBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSubBin)],
                             p.value = inlaNullPvalT_TurnBin) 

min1_TurnBin <- min(postProb1_TurnBin[, c('upper', 'observed')], na.rm = TRUE)
max1_TurnBin <- max(postProb1_TurnBin[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnBin[[1]] <- ggplot(postProb1_TurnBin, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnBin[[2]] <- ggplot(postProb1_TurnBin, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_TurnBin, max1_TurnBin), y = c(min1_TurnBin, max1_TurnBin)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_TurnBin<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_TurnBin[i]<-inla.pmarginal(q=vegINLA$nvcTurnSubBin[i],
                                        # lots of infinity values, so this will not work
                                        marginal=fitTemp_TurnBin$marginals.fitted.values[[i]])
}
postProb2_TurnBin <- data.frame(predicted = fitTemp_TurnBin$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSubBin)],
                             observed = vegINLA$nvcTurnSubBin,
                             lower = fitTemp_TurnBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSubBin)],
                             upper = fitTemp_TurnBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSubBin)],
                             p.value = inlaTempPvalT_TurnBin) 

min2_TurnBin <- min(postProb2_TurnBin[, c('upper', 'observed')], na.rm = TRUE)
max2_TurnBin <- max(postProb2_TurnBin[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnBin[[3]] <- ggplot(postProb2_TurnBin, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnBin[[4]] <- ggplot(postProb2_TurnBin, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_TurnBin, max2_TurnBin), y = c(min2_TurnBin, max2_TurnBin)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_TurnBin<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_TurnBin[i]<-inla.pmarginal(q=vegINLA$nvcTurnSubBin[i],
                                               marginal=fitSpatTempAR1_TurnBin$marginals.fitted.values[[i]])
}

postProb3_TurnBin <- data.frame(predicted = fitSpatTempAR1_TurnBin$summary.fitted.values$mean[1:length(vegINLA$nvcTurnSubBin)],
                             observed = vegINLA$nvcTurnSubBin,
                             lower = fitSpatTempAR1_TurnBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$nvcTurnSubBin)],
                             upper = fitSpatTempAR1_TurnBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$nvcTurnSubBin)],
                             p.value = inlaSpatTempAR1PvalT_TurnBin) 

min3_TurnBin <- min(postProb3_TurnBin[, c('upper', 'observed')], na.rm = TRUE)
max3_TurnBin <- max(postProb3_TurnBin[, c('upper', 'observed')], na.rm = TRUE)

plots_TurnBin[[5]] <- ggplot(postProb3_TurnBin, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_TurnBin[[6]] <- ggplot(postProb3_TurnBin, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_TurnBin, max3_TurnBin), y = c(min3_TurnBin, max3_TurnBin)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_TurnBin, nrow = 3)

# tests for alternative AR structures hidden
# inlaSpatTempAR2PvalT_TurnBin<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR2PvalT_TurnBin[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR2_TurnBin$marginals.fitted.values[[i]])
# }
# 
# postProb4_TurnBin <- data.frame(predicted = fitSpatTempAR2_TurnBin$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR2_TurnBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR2_TurnBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR2PvalT_TurnBin) 
# 
# min4_TurnBin <- min(postProb4_TurnBin[, c('upper', 'observed')], na.rm = TRUE)
# max4_TurnBin <- max(postProb4_TurnBin[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_TurnBin[[7]] <- ggplot(postProb4_TurnBin, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_TurnBin[[8]] <- ggplot(postProb4_TurnBin, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min4_TurnBin, max4_TurnBin), y = c(min4_TurnBin, max4_TurnBin)) + 
#   theme_bw()
# 
# inlaSpatTempAR3PvalT_TurnBin<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR3PvalT_TurnBin[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR3_TurnBin$marginals.fitted.values[[i]])
# }
# 
# postProb5_TurnBin <- data.frame(predicted = fitSpatTempAR3_TurnBin$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR3_TurnBin$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR3_TurnBin$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR3PvalT_TurnBin) 
# 
# min5_TurnBin <- min(postProb5_TurnBin[, c('upper', 'observed')], na.rm = TRUE)
# max5_TurnBin <- max(postProb5_TurnBin[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_TurnBin[[9]] <- ggplot(postProb5_TurnBin, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_TurnBin[[10]] <- ggplot(postProb5_TurnBin, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min5_TurnBin, max5_TurnBin), y = c(min5_TurnBin, max5_TurnBin)) + 
#   theme_bw()
# 
# cowplot::plot_grid(plotlist = plots_TurnBin, nrow = 6)

#### Exports ####
write_rds(fitNull_TurnBin, "./Models/Objects/fitNull_TurnBin.rds")
write_rds(fitTemp_TurnBin, "./Models/Objects/fitTemp_TurnBin.rds")
write_rds(fitSpatTempAR1_TurnBin, "./Models/Objects/fitSpatTempAR1_TurnBin.rds")
#write_rds(fitSpatTempAR2_TurnBin, "./Models/Objects/fitSpatTempAR2_TurnBin.rds")
#write_rds(fitSpatTempAR3_TurnBin, "./Models/Objects/fitSpatTempAR3_TurnBin.rds")