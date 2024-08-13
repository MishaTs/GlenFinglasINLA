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
mesh_BryoRich <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_BryoRich) + theme_bw()
#mesh_BryoRich$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_BryoRich <- inla.spde2.matern(mesh = mesh_BryoRich)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_BryoRich <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)

#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_BryoRich <- ~ #Intercept(1) + 
  # interaction terms don't change sheepOff but kill the signal in cowOff due to construction
  sheepOff + cowOff + 
  #treat +
  # humidity, precipitation, and pressure are too correlated to include in a single model
  # simply by effect size and significance, pressure is best; harder to interpret though
  # frost days have no effect
  windSpeed + #meanT + 
  solarRad + precip + 
  # specific inlabru interaction term syntax
  mix(~ -1 + meanT:elev, model = "fixed") +
  # dung, bareGround, and wood have almost no observations
  # water, rock, and mud also fairly low
  # flowers and tussock take away from the NVC community variables (less positive)
  litter + rock + #graze +
  # flower and tussock are mainly M and U
  elev +
  isW + isH + isU + isM

inlaTempC_BryoRich <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + #meanT + 
  solarRad + precip + 
  # specific inlabru interaction term syntax
  mix(~ -1 + meanT:elev, model = "fixed") +
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

inlaSpatTempAR1C_BryoRich <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + #meanT + 
  solarRad + precip + 
  # specific inlabru interaction term syntax
  mix(~ -1 + meanT:elev, model = "fixed") +
  litter + rock + #graze +
  #elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoRich, 
        group = Year, 
        group_mapper = yearMap_BryoRich,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

inlaSpatTempAR2C_BryoRich <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoRich, 
        group = Year, 
        group_mapper = yearMap_BryoRich,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 2))

inlaSpatTempAR3C_BryoRich <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_BryoRich, 
        group = Year, 
        group_mapper = yearMap_BryoRich,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_BryoRich <- like(
  formula = hillBryo0 ~ .,
  # full distribution list: inla.list.models()
  # possible distributions: nbionimal and poisson
  family = "nbinomial", 
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_BryoRich)
)

likeTemp_BryoRich <- like(
  formula = hillBryo0 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "nbinomial", 
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_BryoRich)
)

likeSpatTempAR_BryoRich <- like(
  formula = hillBryo0 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "nbinomial", 
  data=vegINLA,
  domain = list(geometry = mesh_BryoRich)
)

# run the actual fit
# null model
fitNull_BryoRich <- bru(components = inlaNullC_BryoRich,
                       likeNull_BryoRich,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         # might want to consider to setting ``control.predictor=list(link=...)''
                         # otherwise fitted.values in the result-object computed w/ identity for NAs
                         # https://rdrr.io/github/inbo/INLA/man/control.predictor.html
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# temporal model
fitTemp_BryoRich <- bru(components = inlaTempC_BryoRich,
                       likeTemp_BryoRich,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         # might want to consider to setting ``control.predictor=list(link=...)''
                         # otherwise fitted.values in the result-object computed w/ identity for NAs
                         # https://rdrr.io/github/inbo/INLA/man/control.predictor.html
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# spatiotemporal model
# fairly memory-intensive
fitSpatTempAR1_BryoRich <- bru(components = inlaSpatTempAR1C_BryoRich,
                              likeSpatTempAR_BryoRich,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                # might want to consider to setting ``control.predictor=list(link=...)''
                                # otherwise fitted.values in the result-object computed w/ identity for NAs
                                # https://rdrr.io/github/inbo/INLA/man/control.predictor.html
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR2_BryoRich <- bru(components = inlaSpatTempAR2C_BryoRich,
                              likeSpatTempAR_BryoRich,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR3_BryoRich <- bru(components = inlaSpatTempAR3C_BryoRich,
                              likeSpatTempAR_BryoRich,
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
cbind(fitNull_BryoRich$waic$waic[1], fitNull_BryoRich$dic$dic[1], -sum(log(fitNull_BryoRich$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_BryoRich$waic$waic[1], fitTemp_BryoRich$dic$dic[1], -sum(log(fitTemp_BryoRich$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_BryoRich$waic$waic[1], fitSpatTempAR1_BryoRich$dic$dic[1], -sum(log(fitSpatTempAR1_BryoRich$cpo$cpo), na.rm = TRUE))
#[1,] 15679.96 15685.27 7839.985
#[1,] 15545.16 15553.15 7772.586
#[1,] 15374.49 15393.96 7687.258

# no major improvement via AR2-AR3
cbind(fitSpatTempAR2_BryoRich$waic$waic[1], fitSpatTempAR2_BryoRich$dic$dic[1], -sum(log(fitSpatTempAR2_BryoRich$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR3_BryoRich$waic$waic[1], fitSpatTempAR3_BryoRich$dic$dic[1], -sum(log(fitSpatTempAR3_BryoRich$cpo$cpo), na.rm = TRUE))



# summary
summary(fitSpatTempAR1_BryoRich)

#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_BryoRich)
brinla::bri.fixed.plot(fitTemp_BryoRich)
brinla::bri.fixed.plot(fitSpatTempAR1_BryoRich)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTempAR1_BryoRich$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTempAR1_BryoRich$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoRich$marginals.fixed$sheepOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoRich$marginals.fixed$cowOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTempAR1_BryoRich$marginals.lincomb.derived$offtake)


#### Marginal distribution ####

plots_BryoRich <- list()

# Null model
inlaNullPvalT_BryoRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_BryoRich[i]<-inla.pmarginal(q=vegINLA$hillBryo0[i],
                                           marginal=fitNull_BryoRich$marginals.fitted.values[[i]])
}

postProb1_BryoRich <- data.frame(predicted = fitNull_BryoRich$summary.fitted.values$mean[1:length(vegINLA$hillBryo0)],
                                observed = vegINLA$hillBryo0,
                                lower = fitNull_BryoRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo0)],
                                upper = fitNull_BryoRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo0)],
                                p.value = inlaNullPvalT_BryoRich) 

min1_BryoRich <- min(postProb1_BryoRich[, c('upper', 'observed')], na.rm = TRUE)
max1_BryoRich <- max(postProb1_BryoRich[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoRich[[1]] <- ggplot(postProb1_BryoRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoRich[[2]] <- ggplot(postProb1_BryoRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_BryoRich, max1_BryoRich), y = c(min1_BryoRich, max1_BryoRich)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_BryoRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_BryoRich[i]<-inla.pmarginal(q=vegINLA$hillBryo0[i],
                                           # lots of infinity values, so this will not work
                                           marginal=fitTemp_BryoRich$marginals.fitted.values[[i]])
}
postProb2_BryoRich <- data.frame(predicted = fitTemp_BryoRich$summary.fitted.values$mean[1:length(vegINLA$hillBryo0)],
                                observed = vegINLA$hillBryo0,
                                lower = fitTemp_BryoRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo0)],
                                upper = fitTemp_BryoRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo0)],
                                p.value = inlaTempPvalT_BryoRich) 

min2_BryoRich <- min(postProb2_BryoRich[, c('upper', 'observed')], na.rm = TRUE)
max2_BryoRich <- max(postProb2_BryoRich[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoRich[[3]] <- ggplot(postProb2_BryoRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoRich[[4]] <- ggplot(postProb2_BryoRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_BryoRich, max2_BryoRich), y = c(min2_BryoRich, max2_BryoRich)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_BryoRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_BryoRich[i]<-inla.pmarginal(q=vegINLA$hillBryo0[i],
                                                  marginal=fitSpatTempAR1_BryoRich$marginals.fitted.values[[i]])
}

postProb3_BryoRich <- data.frame(predicted = fitSpatTempAR1_BryoRich$summary.fitted.values$mean[1:length(vegINLA$hillBryo0)],
                                observed = vegINLA$hillBryo0,
                                lower = fitSpatTempAR1_BryoRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo0)],
                                upper = fitSpatTempAR1_BryoRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo0)],
                                p.value = inlaSpatTempAR1PvalT_BryoRich) 

min3_BryoRich <- min(postProb3_BryoRich[, c('upper', 'observed')], na.rm = TRUE)
max3_BryoRich <- max(postProb3_BryoRich[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoRich[[5]] <- ggplot(postProb3_BryoRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoRich[[6]] <- ggplot(postProb3_BryoRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_BryoRich, max3_BryoRich), y = c(min3_BryoRich, max3_BryoRich)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_BryoRich, nrow = 3)

inlaSpatTempAR2PvalT_BryoRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR2PvalT_BryoRich[i]<-inla.pmarginal(q=vegINLA$hillBryo0[i],
                                                  marginal=fitSpatTempAR2_BryoRich$marginals.fitted.values[[i]])
}

postProb4_BryoRich <- data.frame(predicted = fitSpatTempAR2_BryoRich$summary.fitted.values$mean[1:length(vegINLA$hillBryo0)],
                                observed = vegINLA$hillBryo0,
                                lower = fitSpatTempAR2_BryoRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo0)],
                                upper = fitSpatTempAR2_BryoRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo0)],
                                p.value = inlaSpatTempAR2PvalT_BryoRich) 

min4_BryoRich <- min(postProb4_BryoRich[, c('upper', 'observed')], na.rm = TRUE)
max4_BryoRich <- max(postProb4_BryoRich[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoRich[[7]] <- ggplot(postProb4_BryoRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoRich[[8]] <- ggplot(postProb4_BryoRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min4_BryoRich, max4_BryoRich), y = c(min4_BryoRich, max4_BryoRich)) + 
  theme_bw()

inlaSpatTempAR3PvalT_BryoRich<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR3PvalT_BryoRich[i]<-inla.pmarginal(q=vegINLA$hillBryo0[i],
                                                  marginal=fitSpatTempAR3_BryoRich$marginals.fitted.values[[i]])
}

postProb5_BryoRich <- data.frame(predicted = fitSpatTempAR3_BryoRich$summary.fitted.values$mean[1:length(vegINLA$hillBryo0)],
                                observed = vegINLA$hillBryo0,
                                lower = fitSpatTempAR3_BryoRich$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo0)],
                                upper = fitSpatTempAR3_BryoRich$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo0)],
                                p.value = inlaSpatTempAR3PvalT_BryoRich) 

min5_BryoRich <- min(postProb5_BryoRich[, c('upper', 'observed')], na.rm = TRUE)
max5_BryoRich <- max(postProb5_BryoRich[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoRich[[9]] <- ggplot(postProb5_BryoRich, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoRich[[10]] <- ggplot(postProb5_BryoRich, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min5_BryoRich, max5_BryoRich), y = c(min5_BryoRich, max5_BryoRich)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_BryoRich, nrow = 5)

#### Exports ####
write_rds(fitNull_BryoRich, "./Models/Objects/fitNull_BryoRich.rds")
write_rds(fitTemp_BryoRich, "./Models/Objects/fitTemp_BryoRich.rds")
write_rds(fitSpatTempAR1_BryoRich, "./Models/Objects/fitSpatTempAR1_BryoRich.rds")
write_rds(fitSpatTempAR2_BryoRich, "./Models/Objects/fitSpatTempAR2_BryoRich.rds")
write_rds(fitSpatTempAR3_BryoRich, "./Models/Objects/fitSpatTempAR3_BryoRich.rds")