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
                             sameObserver = as.numeric(Recorder == lag(Recorder))) #%>% filter(!is.na(nvcM1)) 

# t(apply(table(vegINLA$isW, vegINLA$flower),1, function(x) x/sum(x)))

# first, the mesh
mesh_Height <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_Height) + theme_bw()
#mesh_Height$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_Height <- inla.spde2.matern(mesh = mesh_Height)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_Height <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)

#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_Height <- ~ #Intercept(1) + 
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
  dungD + whiteStick +
  sameObserver

inlaTempC_Height <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  dungD + whiteStick + 
  sameObserver +
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

inlaSpatTempAR1C_Height <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  dungD + whiteStick + 
  sameObserver +
  field(geometry, 
        model = aMatSPDE_Height, 
        group = Year, 
        group_mapper = yearMap_Height,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

# inlaSpatTempAR2C_Height <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_Height, 
#         group = Year, 
#         group_mapper = yearMap_Height,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 2))
# 
# inlaSpatTempAR3C_Height <- ~ #Intercept(1) + 
#   sheepOff + cowOff + 
#   #treat +
#   windSpeed + meanT + solarRad + precip + 
#   litter + rock + #graze +
#   elev +
#   isW + isH + isU + isM +
#   dungD + whiteStick + 
#   field(geometry, 
#         model = aMatSPDE_Height, 
#         group = Year, 
#         group_mapper = yearMap_Height,
#         control.group = list(model = "ar",
#                              # order <= 10
#                              order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_Height <- like(
  formula = htAvg ~ .,
  # full distribution list: inla.list.models()
  # possible distributions: tweedie, lognormal, normal
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_Height)
)

likeTemp_Height <- like(
  formula = htAvg ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_Height)
)

likeSpatTempAR_Height <- like(
  formula = htAvg ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  domain = list(geometry = mesh_Height)
)

# run the actual fit
# null model
# does not accept factor variables
# Error in ibm_jacobian.bru_mapper_linear(mapper[["mappers"]][[x]], input = input[[x]],  : 
#   The input to a bru_mapper_linear evaluation must be numeric or logical.
fitNull_Height <- bru(components = inlaNullC_Height,
                       likeNull_Height,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# temporal model
# for gaussian family:
#Assertion failed: (arg->Q->a[0] >= 0.0), function GMRFLib_tabulate_Qfunc_core, file tabulate-Qfunc.c, line 147.
#sh: line 1: 15505 Segmentation fault: 11  '/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/INLA/bin/mac.arm64/inla.run' -s -v -t16:1 -P compact '/private/var/folders/m4/xyhttgk14tl61_q2f5y9yj7r0000gn/T/RtmpDkSlB6/file161666a5d6bdf/Model.ini' > '/private/var/folders/m4/xyhttgk14tl61_q2f5y9yj7r0000gn/T/RtmpDkSlB6/file161666a5d6bdf/Logfile.txt'
# *** inla.core.safe:  The inla program failed, but will rerun in case better initial values may help. try=1/1 
# *** inla.core.safe:  rerun with improved initial values 
fitTemp_Height <- bru(components = inlaTempC_Height,
                       likeTemp_Height,
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
fitSpatTempAR1_Height <- bru(components = inlaSpatTempAR1C_Height,
                              likeSpatTempAR_Height,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

# AR2 does not converge
# fitSpatTempAR2_Height <- bru(components = inlaSpatTempAR2C_Height,
#                              likeSpatTempAR_Height,
#                              options = list(
#                                # add linear combinations
#                                lincomb = c(lc1, lc2),
#                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
#                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
#                                control.predictor = list(compute=TRUE),
#                                #control.inla = list(int.strategy = "eb"),
#                                verbose = FALSE))
# 
# # AR3 not tested
# fitSpatTempAR3_Height <- bru(components = inlaSpatTempAR3C_Height,
#                              likeSpatTempAR_Height,
#                              options = list(
#                                # add linear combinations
#                                lincomb = c(lc1, lc2),
#                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
#                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
#                                control.predictor = list(compute=TRUE),
#                                #control.inla = list(int.strategy = "eb"),
#                                verbose = FALSE))


#### Model Output ####
# WAIC, DIC, CPO
# CPO needs no NAs to work perfectly; na.rm is good otherwise
cbind(fitNull_Height$waic$waic[1], fitNull_Height$dic$dic[1], -sum(log(fitNull_Height$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_Height$waic$waic[1], fitTemp_Height$dic$dic[1], -sum(log(fitTemp_Height$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_Height$waic$waic[1], fitSpatTempAR1_Height$dic$dic[1], -sum(log(fitSpatTempAR1_Height$cpo$cpo), na.rm = TRUE))
#[1,] 344111.1 344108.4 172055.5
#[1,] 338394.4 338397.8 169198.1
#[1,] 323419.9 323347.4 161751.6

# how about sameObserver?
#[1,] 344020.4 344018.2 172010.2
#[1,] 338091.9 338094.6 169046.8
#[1,] 323382.8 323318 161731.4

# no benefits to using AR2 for height
# cbind(fitSpatTempAR2_Height$waic$waic[1], fitSpatTempAR2_Height$dic$dic[1], -sum(log(fitSpatTempAR2_Height$cpo$cpo), na.rm = TRUE))
# cbind(fitSpatTempAR3_Height$waic$waic[1], fitSpatTempAR3_Height$dic$dic[1], -sum(log(fitSpatTempAR3_Height$cpo$cpo), na.rm = TRUE))


# summary
summary(fitSpatTempAR1_Height)

#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_Height)
brinla::bri.fixed.plot(fitTemp_Height)
brinla::bri.fixed.plot(fitSpatTempAR1_Height)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTempAR1_Height$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTempAR1_Height$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_Height$marginals.fixed$sheepOff)
1 - inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_Height$marginals.fixed$cowOff)
1 - inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_Height$marginals.lincomb.derived$offtake)


#### Marginal distribution ####

plots_Height <- list()

# Null model
inlaNullPvalT_Height<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_Height[i]<-inla.pmarginal(q=vegINLA$htAvg[i],
                                           marginal=fitNull_Height$marginals.fitted.values[[i]])
}

postProb1_Height <- data.frame(predicted = fitNull_Height$summary.fitted.values$mean[1:length(vegINLA$htAvg)],
                                observed = vegINLA$htAvg,
                                lower = fitNull_Height$summary.fitted.values$`0.025quant`[1:length(vegINLA$htAvg)],
                                upper = fitNull_Height$summary.fitted.values$`0.975quant`[1:length(vegINLA$htAvg)],
                                p.value = inlaNullPvalT_Height) 

min1_Height <- min(postProb1_Height[, c('upper', 'observed')], na.rm = TRUE)
max1_Height <- max(postProb1_Height[, c('upper', 'observed')], na.rm = TRUE)

plots_Height[[1]] <- ggplot(postProb1_Height, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_Height[[2]] <- ggplot(postProb1_Height, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_Height, max1_Height), y = c(min1_Height, max1_Height)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_Height<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_Height[i]<-inla.pmarginal(q=vegINLA$htAvg[i],
                                           # lots of infinity values, so this will not work
                                           marginal=fitTemp_Height$marginals.fitted.values[[i]])
}
postProb2_Height <- data.frame(predicted = fitTemp_Height$summary.fitted.values$mean[1:length(vegINLA$htAvg)],
                                observed = vegINLA$htAvg,
                                lower = fitTemp_Height$summary.fitted.values$`0.025quant`[1:length(vegINLA$htAvg)],
                                upper = fitTemp_Height$summary.fitted.values$`0.975quant`[1:length(vegINLA$htAvg)],
                                p.value = inlaTempPvalT_Height) 

min2_Height <- min(postProb2_Height[, c('upper', 'observed')], na.rm = TRUE)
max2_Height <- max(postProb2_Height[, c('upper', 'observed')], na.rm = TRUE)

plots_Height[[3]] <- ggplot(postProb2_Height, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_Height[[4]] <- ggplot(postProb2_Height, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_Height, max2_Height), y = c(min2_Height, max2_Height)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_Height<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_Height[i]<-inla.pmarginal(q=vegINLA$htAvg[i],
                                                  marginal=fitSpatTempAR1_Height$marginals.fitted.values[[i]])
}

postProb3_Height <- data.frame(predicted = fitSpatTempAR1_Height$summary.fitted.values$mean[1:length(vegINLA$htAvg)],
                                observed = vegINLA$htAvg,
                                lower = fitSpatTempAR1_Height$summary.fitted.values$`0.025quant`[1:length(vegINLA$htAvg)],
                                upper = fitSpatTempAR1_Height$summary.fitted.values$`0.975quant`[1:length(vegINLA$htAvg)],
                                p.value = inlaSpatTempAR1PvalT_Height) 

min3_Height <- min(postProb3_Height[, c('upper', 'observed')], na.rm = TRUE)
max3_Height <- max(postProb3_Height[, c('upper', 'observed')], na.rm = TRUE)

plots_Height[[5]] <- ggplot(postProb3_Height, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_Height[[6]] <- ggplot(postProb3_Height, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_Height, max3_Height), y = c(min3_Height, max3_Height)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_Height, nrow = 3)

# inlaSpatTempAR2PvalT_Height<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR2PvalT_Height[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR2_Height$marginals.fitted.values[[i]])
# }
# 
# postProb4_Height <- data.frame(predicted = fitSpatTempAR2_Height$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR2_Height$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR2_Height$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR2PvalT_Height) 
# 
# min4_Height <- min(postProb4_Height[, c('upper', 'observed')], na.rm = TRUE)
# max4_Height <- max(postProb4_Height[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_Height[[7]] <- ggplot(postProb4_Height, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_Height[[8]] <- ggplot(postProb4_Height, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min4_Height, max4_Height), y = c(min4_Height, max4_Height)) + 
#   theme_bw()
# 
# inlaSpatTempAR3PvalT_Height<-rep(NA, nrow=(vegINLA))
# for(i in 1:nrow(vegINLA)){
#   inlaSpatTempAR3PvalT_Height[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
#                                                   marginal=fitSpatTempAR3_Height$marginals.fitted.values[[i]])
# }
# 
# postProb5_Height <- data.frame(predicted = fitSpatTempAR3_Height$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
#                                 observed = vegINLA$hillBryo1,
#                                 lower = fitSpatTempAR3_Height$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
#                                 upper = fitSpatTempAR3_Height$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
#                                 p.value = inlaSpatTempAR3PvalT_Height) 
# 
# min5_Height <- min(postProb5_Height[, c('upper', 'observed')], na.rm = TRUE)
# max5_Height <- max(postProb5_Height[, c('upper', 'observed')], na.rm = TRUE)
# 
# plots_Height[[9]] <- ggplot(postProb5_Height, aes(x = p.value)) + 
#   geom_histogram() +
#   labs(y = "Count", x = "Posterior probability") +
#   theme_bw()
# 
# plots_Height[[10]] <- ggplot(postProb5_Height, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(y = "Observed", x = "Fitted") +
#   lims(x = c(min5_Height, max5_Height), y = c(min5_Height, max5_Height)) + 
#   theme_bw()
# 
# cowplot::plot_grid(plotlist = plots_Height, nrow = 6)



#### Plotting ####
# currently non-functional
pix_Height <- fm_pixels(mesh_Height)
#Warning in handle_problems(e_input) :
#The input evaluation 'Year' for 'field.group' failed. Perhaps the data object doesn't contain the needed variables? Falling back to '1'.
pred_Height <- predict(
  fitSpatTempAR1_Height, pix_Height,
  # it looks like we need to provide actual values for covariates (and year) for predictions to work
  ~ field #+ Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + rock + elev + isW + isH + isU + isM 
)
samp <- generate(fitSpatTempAR1_Height, pix_Height,
                 ~ field + Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + rock + elev + isW + isH + isU + isM,
                 n.samples = 1
)
pred_Height$sample <- samp[, 1]

#### Exports ####
write_rds(fitNull_Height, "./Models/Objects/fitNull_Height.rds")
write_rds(fitTemp_Height, "./Models/Objects/fitTemp_Height.rds")
write_rds(fitSpatTempAR1_Height, "./Models/Objects/fitSpatTempAR1_Height.rds")
write_rds(fitSpatTempAR2_Height, "./Models/Objects/fitSpatTempAR2_Height.rds")
write_rds(fitSpatTempAR3_Height, "./Models/Objects/fitSpatTempAR3_Height.rds")
