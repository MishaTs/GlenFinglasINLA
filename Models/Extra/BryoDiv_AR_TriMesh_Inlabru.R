library(tidyverse)
library(INLA)
library(fmesher)
library(inlabru)

#### Data Setup ####
vegGeo <- read_rds("./Vegetation/vegGeo.rds")

vegINLA <- vegGeo %>% mutate(fBlock = as.factor(Block),
                             Plot_code = paste0(Plot, Block),
                             treat = factor(Plot,
                                            levels = c("III", "I", "II", "IV"))) #%>% filter(!is.na(nvcM1))

# t(apply(table(vegINLA$isW, vegINLA$flower),1, function(x) x/sum(x)))

# first, the mesh
mesh_BryoDiv <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_BryoDiv) + theme_bw()
#mesh_BryoDiv$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_BryoDiv <- inla.spde2.matern(mesh = mesh_BryoDiv)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_BryoDiv <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)
yearMapTri_BryoDiv <- bru_mapper(fm_mesh_1d(sort(as.vector(unique(st_drop_geometry(vegINLA) %>% 
                                                                    filter(!is.na(hillVasc0)) %>% 
                                                                    select(Year)))$Year)), indexed = TRUE)



#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_BryoDiv <- ~ #Intercept(1) + 
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
  isW + isH + isU + isM

inlaTempC_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

inlaSpatTempAR1C_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMap_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

inlaSpatTempAR2C_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMap_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 2))

inlaSpatTempAR3C_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMap_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 3))

inlaSpatTempAR1TriC_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMapTri_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

inlaSpatTempAR2TriC_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMapTri_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 2))

inlaSpatTempAR3TriC_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoDiv, 
        group = Year, 
        group_mapper = yearMapTri_BryoDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_BryoDiv <- like(
  formula = hillBryo1 ~ .,
  # possible distributions: gamma and tweedie
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  # dtweedie: Reached upper limit. Increase TWEEDIE_MAX_IDX or scale your data!
  scale = 0.5,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_BryoDiv)
)

likeTemp_BryoDiv <- like(
  formula = hillBryo1 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_BryoDiv)
)

likeSpatTempAR_BryoDiv <- like(
  formula = hillBryo1 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  domain = list(geometry = mesh_BryoDiv)
)

# run the actual fit
# null model
fitNull_BryoDiv <- bru(components = inlaNullC_BryoDiv,
                       likeNull_BryoDiv,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# temporal model
fitTemp_BryoDiv <- bru(components = inlaTempC_BryoDiv,
                       likeTemp_BryoDiv,
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
fitSpatTempAR1_BryoDiv <- bru(components = inlaSpatTempAR1C_BryoDiv,
                           likeSpatTempAR_BryoDiv,
                           options = list(
                             # add linear combinations
                             lincomb = c(lc1, lc2),
                             #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                             control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                             control.predictor = list(compute=TRUE),
                             #control.inla = list(int.strategy = "eb"),
                             verbose = FALSE))

fitSpatTempAR2_BryoDiv <- bru(components = inlaSpatTempAR2C_BryoDiv,
                              likeSpatTempAR_BryoDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR3_BryoDiv <- bru(components = inlaSpatTempAR3C_BryoDiv,
                              likeSpatTempAR_BryoDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR1Tri_BryoDiv <- bru(components = inlaSpatTempAR1TriC_BryoDiv,
                              likeSpatTempAR_BryoDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR2Tri_BryoDiv <- bru(components = inlaSpatTempAR2TriC_BryoDiv,
                              likeSpatTempAR_BryoDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR3Tri_BryoDiv <- bru(components = inlaSpatTempAR3TriC_BryoDiv,
                              likeSpatTempAR_BryoDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

#### Model Output ####
# WAIC, DIC, CPO
# CPO needs no NAs to work perfectly; na.rm is good otherwise
cbind(fitNull_BryoDiv$waic$waic[1], fitNull_BryoDiv$dic$dic[1], -sum(log(fitNull_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_BryoDiv$waic$waic[1], fitTemp_BryoDiv$dic$dic[1], -sum(log(fitTemp_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_BryoDiv$waic$waic[1], fitSpatTempAR1_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR1_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR2_BryoDiv$waic$waic[1], fitSpatTempAR2_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR2_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR3_BryoDiv$waic$waic[1], fitSpatTempAR3_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR3_BryoDiv$cpo$cpo), na.rm = TRUE))
# ar1
#[1,] 12025.41 11994.13 6041.836
#[1,] 12024.58 11993.68 6041.598
# ar2
#[1,] 12003.54 11987.77 6038.287
#[1,] 12007.56 11974.84 6039.013
# ar3
#[1,] 11993.83 11949.38 6030.014
#[1,] 11989.41 11953.99 6027.315
# ar4 -- worse DIC and CPO and only slightly better WAIC (maybe)
#[1,] 11992.9 11967.33 6032.676

# summary
summary(fitSpatTempAR1_BryoDiv)

# info on priors
#fitSpatTempAR1_BryoDiv$bru_info$model$effects$field
#inla.set.control.fixed.default()
#names(inla.models()$latent$ar)
#names(fitSpatTempAR1_BryoDiv$all.hyper$predictor$hyper$theta)


#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_BryoDiv)
brinla::bri.fixed.plot(fitTemp_BryoDiv)
brinla::bri.fixed.plot(fitSpatTempAR1_BryoDiv)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTempAR1_BryoDiv$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTempAR1_BryoDiv$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoDiv$marginals.fixed$sheepOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoDiv$marginals.fixed$cowOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoDiv$marginals.lincomb.derived$offtake)


#### Marginal distribution ####

plots_BryoDiv <- list()

# Null model
inlaNullPvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                           marginal=fitNull_BryoDiv$marginals.fitted.values[[i]])
}

postProb1_BryoDiv <- data.frame(predicted = fitNull_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitNull_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitNull_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaNullPvalT_BryoDiv) 

min1_BryoDiv <- min(postProb1_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max1_BryoDiv <- max(postProb1_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[1]] <- ggplot(postProb1_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[2]] <- ggplot(postProb1_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_BryoDiv, max1_BryoDiv), y = c(min1_BryoDiv, max1_BryoDiv)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                           # lots of infinity values, so this will not work
                                           marginal=fitTemp_BryoDiv$marginals.fitted.values[[i]])
}
postProb2_BryoDiv <- data.frame(predicted = fitTemp_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitTemp_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitTemp_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaTempPvalT_BryoDiv) 

min2_BryoDiv <- min(postProb2_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max2_BryoDiv <- max(postProb2_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[3]] <- ggplot(postProb2_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[4]] <- ggplot(postProb2_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_BryoDiv, max2_BryoDiv), y = c(min2_BryoDiv, max2_BryoDiv)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                               marginal=fitSpatTempAR1_BryoDiv$marginals.fitted.values[[i]])
}

postProb3_BryoDiv <- data.frame(predicted = fitSpatTempAR1_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR1_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR1_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR1PvalT_BryoDiv) 

min3_BryoDiv <- min(postProb3_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max3_BryoDiv <- max(postProb3_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[5]] <- ggplot(postProb3_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[6]] <- ggplot(postProb3_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_BryoDiv, max3_BryoDiv), y = c(min3_BryoDiv, max3_BryoDiv)) + 
  theme_bw()

inlaSpatTempAR2PvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR2PvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                  marginal=fitSpatTempAR2_BryoDiv$marginals.fitted.values[[i]])
}

postProb4_BryoDiv <- data.frame(predicted = fitSpatTempAR2_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR2_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR2_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR2PvalT_BryoDiv) 

min4_BryoDiv <- min(postProb4_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max4_BryoDiv <- max(postProb4_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[7]] <- ggplot(postProb4_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[8]] <- ggplot(postProb4_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min4_BryoDiv, max4_BryoDiv), y = c(min4_BryoDiv, max4_BryoDiv)) + 
  theme_bw()

inlaSpatTempAR3PvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR3PvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                  marginal=fitSpatTempAR3_BryoDiv$marginals.fitted.values[[i]])
}

postProb5_BryoDiv <- data.frame(predicted = fitSpatTempAR3_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR3_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR3_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR3PvalT_BryoDiv) 

min5_BryoDiv <- min(postProb5_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max5_BryoDiv <- max(postProb5_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv[[9]] <- ggplot(postProb5_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv[[10]] <- ggplot(postProb5_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min5_BryoDiv, max5_BryoDiv), y = c(min5_BryoDiv, max5_BryoDiv)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_BryoDiv, nrow = 6)

cbind(fitSpatTempAR1Tri_BryoDiv$waic$waic[1], fitSpatTempAR1Tri_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR1Tri_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR2Tri_BryoDiv$waic$waic[1], fitSpatTempAR2Tri_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR2Tri_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR3Tri_BryoDiv$waic$waic[1], fitSpatTempAR3Tri_BryoDiv$dic$dic[1], -sum(log(fitSpatTempAR3Tri_BryoDiv$cpo$cpo), na.rm = TRUE))

plots_BryoDivAlt <- plots_BryoDivAlt

# Spatial temporal model
inlaSpatTempAR1TriPvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1TriPvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                     marginal=fitSpatTempAR1Tri_BryoDiv$marginals.fitted.values[[i]])
}

postProb6_BryoDiv <- data.frame(predicted = fitSpatTempAR1Tri_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR1Tri_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR1Tri_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR1TriPvalT_BryoDiv) 

min6_BryoDiv <- min(postProb6_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max6_BryoDiv <- max(postProb6_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDivAlt[[11]] <- ggplot(postProb6_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDivAlt[[12]] <- ggplot(postProb6_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min6_BryoDiv, max6_BryoDiv), y = c(min6_BryoDiv, max6_BryoDiv)) + 
  theme_bw()

inlaSpatTempAR2TriPvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR2TriPvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                     marginal=fitSpatTempAR2Tri_BryoDiv$marginals.fitted.values[[i]])
}

postProb7_BryoDiv <- data.frame(predicted = fitSpatTempAR2Tri_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR2Tri_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR2Tri_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR2TriPvalT_BryoDiv) 

min7_BryoDiv <- min(postProb7_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max7_BryoDiv <- max(postProb7_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDivAlt[[13]] <- ggplot(postProb7_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDivAlt[[14]] <- ggplot(postProb7_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min7_BryoDiv, max7_BryoDiv), y = c(min7_BryoDiv, max7_BryoDiv)) + 
  theme_bw()

inlaSpatTempAR3TriPvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR3TriPvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                     marginal=fitSpatTempAR3Tri_BryoDiv$marginals.fitted.values[[i]])
}

postProb8_BryoDiv <- data.frame(predicted = fitSpatTempAR3Tri_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempAR3Tri_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempAR3Tri_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempAR3TriPvalT_BryoDiv) 

min8_BryoDiv <- min(postProb8_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max8_BryoDiv <- max(postProb8_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDivAlt[[15]] <- ggplot(postProb8_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDivAlt[[16]] <- ggplot(postProb8_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min8_BryoDiv, max8_BryoDiv), y = c(min8_BryoDiv, max8_BryoDiv)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_BryoDivAlt, nrow = 8)


#### Plotting ####
pix_BryoDiv <- fm_pixels(mesh_BryoDiv)
pred_BryoDiv <- predict(
  fitSpatTempAR1_BryoDiv, pix_BryoDiv,
  # it looks like we need to provide actual values for covariates (and year) for predictions to work
  ~ field + Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + rock + elev + isW + isH + isU + isM 
)
samp <- generate(fitSpatTempAR1_BryoDiv, pix_BryoDiv,
                 ~ field + Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + rock + elev + isW + isH + isU + isM,
                 n.samples = 1
)
pred_BryoDiv$sample <- samp[, 1]


