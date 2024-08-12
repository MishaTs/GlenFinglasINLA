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
mesh_VascDiv <- fm_mesh_2d_inla(
  loc = vegINLA$geometry, 
  # doesn't look like any changes where max.edge > 0.05 have much impact
  # much finer mesh possible otherwise, but left as-is for computing benefits for now
  # outer edge set to 1 but doesn't matter
  max.edge = c(0.05,1), # map units inside and outside
  # cutoff is min edge
  #cutoff = 50, offset = c(100, 300),
  crs = fm_crs(vegINLA)
)
ggplot() + gg(mesh_VascDiv) + theme_bw()
#mesh_VascDiv$n # 2044 points

# Making SPDE
# Quantify distance between points and define Matern correlation on the mesh
#inla.spde2.pcmatern allows for penalised complexity priors
aMatSPDE_VascDiv <- inla.spde2.matern(mesh = mesh_VascDiv)#, 
#prior.range = c(10, 0.5), 
#prior.sigma = c(.5, .5))#,
#alpha = 2) 
# code from slide 10 of Finn's ppt
# https://www.maths.ed.ac.uk/~flindgre/talks/Lindgren_RSS2023.pdf

# map across years for the temporal component
yearMap_VascDiv <- bru_mapper(fm_mesh_1d(sort(unique(vegINLA$Year))), indexed = TRUE)


#### Modelling ####
# first define the model components
# high-leverage variables end up being quite impactful (e.g., mostly 0 with a few rows of 1)
inlaNullC_VascDiv <- ~ #Intercept(1) + 
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

inlaTempC_VascDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field1(Year, model='ar1', replicate = Point) +
  field2(fBlock, model = "iid") + 
  field3(Plot_code, model = "iid")

inlaSpatTempAR1C_VascDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_VascDiv, 
        group = Year, 
        group_mapper = yearMap_VascDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 1))

inlaSpatTempAR2C_VascDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_VascDiv, 
        group = Year, 
        group_mapper = yearMap_VascDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 2))

inlaSpatTempAR3C_VascDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  field(geometry, 
        model = aMatSPDE_VascDiv, 
        group = Year, 
        group_mapper = yearMap_VascDiv,
        control.group = list(model = "ar",
                             # order <= 10
                             order = 3))

# define linear combinations
lc1 <- inla.make.lincomb(sheepOff = 1, cowOff = 1)
names(lc1) <- "offtake"
lc2 <- inla.make.lincomb(isW = 1, isH = 1, isU = 1, isM = 1)
names(lc2) <- "nvc"

# build the likelihoods
likeNull_VascDiv <- like(
  formula = hillVasc1 ~ .,
  # possible distributions: gamma and tweedie
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  # dtweedie: Reached upper limit. Increase TWEEDIE_MAX_IDX or scale your data!
  scale = 0.5,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_VascDiv),
)

likeTemp_VascDiv <- like(
  formula = hillVasc1 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  # domain optional for non-spatial model
  domain = list(geometry = mesh_VascDiv),
)

likeSpatTempAR_VascDiv <- like(
  formula = hillVasc1 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  domain = list(geometry = mesh_VascDiv),
)

# run the actual fit
# null model
fitNull_VascDiv <- bru(components = inlaNullC_VascDiv,
                       likeNull_VascDiv,
                       options = list(
                         # add linear combinations
                         lincomb = c(lc1, lc2),
                         #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                         control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                         control.predictor = list(compute=TRUE),
                         #control.inla = list(int.strategy = "eb"),
                         verbose = FALSE))

# temporal model
fitTemp_VascDiv <- bru(components = inlaTempC_VascDiv,
                       likeTemp_VascDiv,
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
fitSpatTempAR1_VascDiv <- bru(components = inlaSpatTempAR1C_VascDiv,
                              likeSpatTempAR_VascDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))


fitSpatTempAR2_VascDiv <- bru(components = inlaSpatTempAR2C_VascDiv,
                              likeSpatTempAR_VascDiv,
                              options = list(
                                # add linear combinations
                                lincomb = c(lc1, lc2),
                                #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                control.predictor = list(compute=TRUE),
                                #control.inla = list(int.strategy = "eb"),
                                verbose = FALSE))

fitSpatTempAR3_VascDiv <- bru(components = inlaSpatTempAR3C_VascDiv,
                              likeSpatTempAR_VascDiv,
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
cbind(fitNull_VascDiv$waic$waic[1], fitNull_VascDiv$dic$dic[1], -sum(log(fitNull_VascDiv$cpo$cpo), na.rm = TRUE))
cbind(fitTemp_VascDiv$waic$waic[1], fitTemp_VascDiv$dic$dic[1], -sum(log(fitTemp_VascDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR1_VascDiv$waic$waic[1], fitSpatTempAR1_VascDiv$dic$dic[1], -sum(log(fitSpatTempAR1_VascDiv$cpo$cpo), na.rm = TRUE))
# [1,] 18355.02 18354.58 9177.527
# [1,] 17886.87 17888.42 8943.589
# [1,] 16999.89 17001.69 8541.799
cbind(fitSpatTempAR2_VascDiv$waic$waic[1], fitSpatTempAR2_VascDiv$dic$dic[1], -sum(log(fitSpatTempAR2_VascDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempAR3_VascDiv$waic$waic[1], fitSpatTempAR3_VascDiv$dic$dic[1], -sum(log(fitSpatTempAR3_VascDiv$cpo$cpo), na.rm = TRUE))
# [1,] 16998.35 16997.57 8541.369
# [1,] 16996.2 16995.93 8540.013

# summary
summary(fitSpatTempAR3_VascDiv)




#### Fixed effects ####
# note that there are no observations without at least one of H, M, U, and W community classifications
# since these are non-exclusive, all are slightly positively linked to diversity
brinla::bri.fixed.plot(fitNull_VascDiv)
brinla::bri.fixed.plot(fitTemp_VascDiv)
brinla::bri.fixed.plot(fitSpatTempAR1_VascDiv)
brinla::bri.fixed.plot(fitSpatTempAR2_VascDiv)
brinla::bri.fixed.plot(fitSpatTempAR3_VascDiv)
# get custom quantiles
# doesn't work due to Inf and -Inf values
#inla.pmarginal(c(0.05, 0.95), fitSpatTemp_VascDiv$marginals.fitted.values)
#inla.qmarginal(p = c(0.12, 0.88),
#               marginal= fitSpatTemp_VascDiv$marginals.fixed$sheepOff)
# this tests the 0 level for a variable
inla.pmarginal(q = 0,
               marginal= fitSpatTemp_VascDiv$marginals.fixed$sheepOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTemp_VascDiv$marginals.fixed$cowOff)
inla.pmarginal(q = 0,
               marginal= fitSpatTemp_VascDiv$marginals.lincomb.derived$offtake)


#### Marginal distribution ####
plots_VascDiv <- list()

# Null model
inlaNullPvalT_VascDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaNullPvalT_VascDiv[i]<-inla.pmarginal(q=vegINLA$hillVasc1[i],
                                           marginal=fitNull_VascDiv$marginals.fitted.values[[i]])
}

postProb1_VascDiv <- data.frame(predicted = fitNull_VascDiv$summary.fitted.values$mean[1:length(vegINLA$hillVasc1)],
                                observed = vegINLA$hillVasc1,
                                lower = fitNull_VascDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc1)],
                                upper = fitNull_VascDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc1)],
                                p.value = inlaNullPvalT_VascDiv) 

min1_VascDiv <- min(postProb1_VascDiv[, c('upper', 'observed')], na.rm = TRUE)
max1_VascDiv <- max(postProb1_VascDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_VascDiv[[1]] <- ggplot(postProb1_VascDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascDiv[[2]] <- ggplot(postProb1_VascDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min1_VascDiv, max1_VascDiv), y = c(min1_VascDiv, max1_VascDiv)) + 
  theme_bw()

# Temporal model
inlaTempPvalT_VascDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaTempPvalT_VascDiv[i]<-inla.pmarginal(q=vegINLA$hillVasc1[i],
                                           # lots of infinity values, so this will not work
                                           marginal=fitTemp_VascDiv$marginals.fitted.values[[i]])
}
postProb2_VascDiv <- data.frame(predicted = fitTemp_VascDiv$summary.fitted.values$mean[1:length(vegINLA$hillVasc1)],
                                observed = vegINLA$hillVasc1,
                                lower = fitTemp_VascDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc1)],
                                upper = fitTemp_VascDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc1)],
                                p.value = inlaTempPvalT_VascDiv) 

min2_VascDiv <- min(postProb2_VascDiv[, c('upper', 'observed')], na.rm = TRUE)
max2_VascDiv <- max(postProb2_VascDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_VascDiv[[3]] <- ggplot(postProb2_VascDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascDiv[[4]] <- ggplot(postProb2_VascDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min2_VascDiv, max2_VascDiv), y = c(min2_VascDiv, max2_VascDiv)) + 
  theme_bw()

# Spatial temporal model
inlaSpatTempAR1PvalT_VascDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR1PvalT_VascDiv[i]<-inla.pmarginal(q=vegINLA$hillVasc1[i],
                                               marginal=fitSpatTempAR1_VascDiv$marginals.fitted.values[[i]])
}

postProb3_VascDiv <- data.frame(predicted = fitSpatTempAR1_VascDiv$summary.fitted.values$mean[1:length(vegINLA$hillVasc1)],
                                observed = vegINLA$hillVasc1,
                                lower = fitSpatTempAR1_VascDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc1)],
                                upper = fitSpatTempAR1_VascDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc1)],
                                p.value = inlaSpatTempAR1PvalT_VascDiv) 

min3_VascDiv <- min(postProb3_VascDiv[, c('upper', 'observed')], na.rm = TRUE)
max3_VascDiv <- max(postProb3_VascDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_VascDiv[[5]] <- ggplot(postProb3_VascDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascDiv[[6]] <- ggplot(postProb3_VascDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min3_VascDiv, max3_VascDiv), y = c(min3_VascDiv, max3_VascDiv)) + 
  theme_bw()

#cowplot::plot_grid(plotlist = plots_VascDiv, nrow = 3)

inlaSpatTempAR2PvalT_VascDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR2PvalT_VascDiv[i]<-inla.pmarginal(q=vegINLA$hillVasc1[i],
                                                  marginal=fitSpatTempAR2_VascDiv$marginals.fitted.values[[i]])
}

postProb4_VascDiv <- data.frame(predicted = fitSpatTempAR2_VascDiv$summary.fitted.values$mean[1:length(vegINLA$hillVasc1)],
                                observed = vegINLA$hillVasc1,
                                lower = fitSpatTempAR2_VascDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc1)],
                                upper = fitSpatTempAR2_VascDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc1)],
                                p.value = inlaSpatTempAR2PvalT_VascDiv) 

min4_VascDiv <- min(postProb4_VascDiv[, c('upper', 'observed')], na.rm = TRUE)
max4_VascDiv <- max(postProb4_VascDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_VascDiv[[7]] <- ggplot(postProb4_VascDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascDiv[[8]] <- ggplot(postProb4_VascDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min4_VascDiv, max4_VascDiv), y = c(min4_VascDiv, max4_VascDiv)) + 
  theme_bw()

inlaSpatTempAR3PvalT_VascDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempAR3PvalT_VascDiv[i]<-inla.pmarginal(q=vegINLA$hillVasc1[i],
                                                  marginal=fitSpatTempAR3_VascDiv$marginals.fitted.values[[i]])
}

postProb5_VascDiv <- data.frame(predicted = fitSpatTempAR3_VascDiv$summary.fitted.values$mean[1:length(vegINLA$hillVasc1)],
                                observed = vegINLA$hillVasc1,
                                lower = fitSpatTempAR3_VascDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillVasc1)],
                                upper = fitSpatTempAR3_VascDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillVasc1)],
                                p.value = inlaSpatTempAR3PvalT_VascDiv) 

min5_VascDiv <- min(postProb5_VascDiv[, c('upper', 'observed')], na.rm = TRUE)
max5_VascDiv <- max(postProb5_VascDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_VascDiv[[9]] <- ggplot(postProb5_VascDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_VascDiv[[10]] <- ggplot(postProb5_VascDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min5_VascDiv, max5_VascDiv), y = c(min5_VascDiv, max5_VascDiv)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_VascDiv, nrow = 5)

#### Plotting ####
# hyperparameter plots
hyper_VascDiv <- data.frame(p = fitSpatTempAR3_VascDiv$marginals.hyperpar$`p parameter for Tweedie`,
                            dispersion = fitSpatTempAR3_VascDiv$marginals.hyperpar$`Dispersion parameter for Tweedie`)
# initial conditions don't seem to affect things too much; lots of information about the parametrisation
# mean and CI for hyperparameters

# non-negative continuous if 1>p>2 w/ zero inflation
# p = 1; quasipoisson; p = 2 is a gamma
ggplot(hyper_VascDiv, aes(x =p.x, y = p.y)) + geom_line() + theme_bw()
ggplot(hyper_VascDiv, aes(x =dispersion.x, y = dispersion.y)) + geom_line() + theme_bw()

# spatial prediction plots

# create data for 2023 predictions and block group labels for plotting later
vegPred <- st_drop_geometry(vegINLA) %>% ungroup() %>% filter(Year == 2023) %>% 
  mutate(blockGroup = factor(Block,
                             levels = LETTERS[1:6],
                             labels = c("AB", "AB", "CD", "CD", "EF", "EF")))

# create mesh only within the study area and not intermediate mesh
ppxl_VascDiv <- fm_pixels(mesh_VascDiv, mask = plotBounds)
# add the raw data for predictions
ppxl_all_VascDiv <- fm_cprod(ppxl, vegPred %>% select(c(Year, sheepOff, cowOff, windSpeed, meanT, solarRad, precip, litter, rock,
                                                        elev, isW, isH, isU, isM, hillBryo1)))

# run the built-in inlabru function to get predictions
spatPred_VascDiv <- predict(
  fitSpatTempAR3_VascDiv, ppxl_all,
  # use the full formula with all covariates
  ~ field + Intercept + sheepOff + cowOff + windSpeed + meanT + solarRad + precip + litter + rock + elev + isW + isH + isU + isM
)

# modify data to get it to the response scale for mean and lower/upper CIs  
spatPred_VascDiv <- spatPred_VascDiv %>% mutate(linkMean = exp(mean),
                                                linkLower = exp(q0.025),
                                                linkUpper = exp(q0.975)) %>% 
  # roundabout way using spatial checks to categorise each observation to the block group
  mutate(ab = ifelse(st_intersects(geometry, plotBounds %>% filter(blockGroup == "A&B")),"A&B",""),
         cd = ifelse(st_intersects(geometry, plotBounds %>% filter(blockGroup == "C&D")),"C&D",""),
         ef = ifelse(st_intersects(geometry, plotBounds %>% filter(blockGroup == "E&F")),"E&F",""),
         blockGroup = paste0(ab, cd, ef))

# plot each block group separately bc geom_sf and gg don't support facet_wrapping
# specific plot elements are "sausage making" and not really relevant beyond visibility
spatPlots <- list()

spatPlots[[1]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NAC&DNA") %>% relocate(linkLower), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(breaks = seq(-4.453, -4.447, by = 0.002),
                     minor_breaks = seq(-4.453, -4.446, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.275, 56.281, by = 0.002),
                     minor_breaks = seq(56.274, 56.281, by = 0.001)) 

spatPlots[[2]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NANAE&F") %>% relocate(linkLower), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(minor_breaks = seq(-4.408, -4.396, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.295, 56.301, by = 0.002),
                     minor_breaks = seq(56.295, 56.301, by = 0.001)) 

spatPlots[[3]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "A&BNANA") %>% relocate(linkLower), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        panel.grid.minor = element_line(colour = "#CCCCCC",
                                        linetype = 2)) +
  labs(fill='Pred Bryo\nDiversity') +
  scale_x_continuous(breaks = seq(-4.380, -4.374, by = 0.002),
                     minor_breaks = seq(-4.381, -4.374, by = 0.001)) +
  scale_y_continuous(minor_breaks = seq(56.262, 56.270, by = 0.001)) 

spatPlots[[4]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NAC&DNA") %>% relocate(linkMean), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(breaks = seq(-4.453, -4.447, by = 0.002),
                     minor_breaks = seq(-4.453, -4.446, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.275, 56.281, by = 0.002),
                     minor_breaks = seq(56.274, 56.281, by = 0.001)) 

spatPlots[[5]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NANAE&F") %>% relocate(linkMean), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(minor_breaks = seq(-4.408, -4.396, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.295, 56.301, by = 0.002),
                     minor_breaks = seq(56.295, 56.301, by = 0.001)) 

spatPlots[[6]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "A&BNANA") %>% relocate(linkMean), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        panel.grid.minor = element_line(colour = "#CCCCCC",
                                        linetype = 2)) +
  labs(fill='Pred Bryo\nDiversity') +
  scale_x_continuous(breaks = seq(-4.380, -4.374, by = 0.002),
                     minor_breaks = seq(-4.381, -4.374, by = 0.001)) +
  scale_y_continuous(minor_breaks = seq(56.262, 56.270, by = 0.001)) 

spatPlots[[7]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NAC&DNA") %>% relocate(linkUpper), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(breaks = seq(-4.453, -4.447, by = 0.002),
                     minor_breaks = seq(-4.453, -4.446, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.275, 56.281, by = 0.002),
                     minor_breaks = seq(56.274, 56.281, by = 0.001)) 

spatPlots[[8]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "NANAE&F") %>% relocate(linkUpper), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        legend.position = 'none') +
  scale_x_continuous(minor_breaks = seq(-4.408, -4.396, by = 0.001)) +
  scale_y_continuous(breaks = seq(56.295, 56.301, by = 0.002),
                     minor_breaks = seq(56.295, 56.301, by = 0.001)) 

spatPlots[[9]] <- ggplot() +
  gg(lambda1 %>% filter(blockGroup == "A&BNANA") %>% relocate(linkUpper), geom = "tile") +
  scale_fill_viridis_c(begin = 0, end = 1, limits=c(min(lambda1$linkLower), max(lambda1$linkUpper))) + 
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_line(colour = "#ABB0B8",
                                  linetype = 2),
        panel.grid.minor = element_line(colour = "#CCCCCC",
                                        linetype = 2)) +
  labs(fill='Pred Bryo\nDiversity') +
  scale_x_continuous(breaks = seq(-4.380, -4.374, by = 0.002),
                     minor_breaks = seq(-4.381, -4.374, by = 0.001)) +
  scale_y_continuous(minor_breaks = seq(56.262, 56.270, by = 0.001)) 

# plot everything together using cowplot() to get it in a fancy grid
# 3 rows for comparing all sub-categories
cowplot::plot_grid(plotlist = spatPlots, nrow = 3,
                   # add labels for figure caption too
                   labels = c('A','','','B','','','C','',''), label_size = 10,
                   align = "h", axis = "b")
# can align widths to make things nicer, but not needed for now
#, rel_widths = c(0.9495,1.5,1.0105))

# test the fit with krieged diversity -- quite bad, but not sure which is "correct
# ggplot(lambda1, aes(x = linkMean, y = hillBryo1)) + geom_point() +
#        geom_abline(slope = 1, intercept = 0) +
#        labs(y = "Observed", x = "Fitted") +
#        lims(x = c(0, 7), y = c(0, 7)) + 
#        theme_bw()

#### Exports ####
write_rds(fitNull_VascDiv, "./Models/Objects/fitNull_VascDiv.rds")
write_rds(fitTemp_VascDiv, "./Models/Objects/fitTemp_VascDiv.rds")
write_rds(fitSpatTempAR1_VascDiv, "./Models/Objects/fitSpatTempAR1_VascDiv.rds")
write_rds(fitSpatTempAR2_VascDiv, "./Models/Objects/fitSpatTempAR2_VascDiv.rds")
write_rds(fitSpatTempAR3_VascDiv, "./Models/Objects/fitSpatTempAR3_VascDiv.rds")

write_rds(spatPred_VascDiv, "./Models/Objects/vascSpatPreds.rds")