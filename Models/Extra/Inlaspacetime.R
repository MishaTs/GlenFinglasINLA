# this is a rudimentary tester script meant to be used in tandem with BryoDivInlabru.R
# models WILL NOT run if run without BryoDivInlabru.R lines 1-33 loaded into the namespace
# posterior checks (lines 60+ here) also compare against the separable spacetime model, so it's good to have those loaded too
library(INLAspacetime)

# extend mesh outside of data range to avoid boundary effects
years_BryoDiv <- seq(min(vegINLA$Year) - 1, max(vegINLA$Year) + 1, by = 1)
meshYr_BryoDiv <- fm_mesh_1d(years_BryoDiv, boundary = "free")
# visualise mesh
ggplot() + gg(meshYr_BryoDiv)

# define model
spdeST_BryoDiv102 <- stModel.define(smesh = mesh_BryoDiv,
                                    tmesh = meshYr_BryoDiv,
                                    # 102, 121, 202 and 220
                                    # with the priors below, separable AR > 202 > 102
                                    # 121 doesn't converge with non-zero a for prt
                                    model = "121",
                                    # priors are required for these models
                                    # prior.rs, prior.rt and prior.sigma as vectors with length two (U, a)
                                    # define the corresponding PC-prior such that
                                    # P(r_s<U)=a, P(r_t<U)=a or P(sigma>U)=a
                                    # If a=0 then U is taken to be the fixed value of the parameter
                                    control.priors = list(
                                      prs = c(0.01, 0.05), ## P(spatial range < 0.001) = 0.1
                                      prt = c(0.05, 0.05), ## P(temporal range < 1) = 0.1
                                      psigma = c(5, 0.05) ## P(sigma > 5) = 0.1
                                      ))

# build likelihoods
inlaSpatTempCNonSep_BryoDiv <- ~ #Intercept(1) + 
  sheepOff + cowOff + 
  #treat +
  windSpeed + meanT + solarRad + precip + 
  litter + rock + #graze +
  elev +
  isW + isH + isU + isM +
  
  field(list(space = geometry, 
             time = Year),
        model = spdeST_BryoDiv102)

likeSpatTempNonSep_BryoDiv <- like(
  formula = hillBryo1 ~ .,
  # tweedie is a superset of gaussian, gamma, poisson, and inverse gaussians
  family = "tweedie", 
  scale = 0.5,
  data=vegINLA,
  #domain = list(geometry = mesh_BryoDiv),
)

# fit the model
fitSpatTempNonSep_BryoDiv <- bru(components = inlaSpatTempCNonSep_BryoDiv,
                                 likeSpatTempNonSep_BryoDiv,
                                 options = list(
                                   # add linear combinations
                                   lincomb = c(lc1, lc2),
                                   #quantiles = c(0.025, 0.05, 0.5, 0.95, 0.975),
                                   control.compute = list(waic=TRUE, cpo=TRUE, dic = TRUE, return.marginals.predictor=TRUE),
                                   control.predictor = list(compute=TRUE),
                                   #control.inla = list(int.strategy = "eb"),
                                   verbose = FALSE))

# test fit metrics
cbind(fitSpatTemp_BryoDiv$waic$waic[1], fitSpatTemp_BryoDiv$dic$dic[1], -sum(log(fitSpatTemp_BryoDiv$cpo$cpo), na.rm = TRUE))
cbind(fitSpatTempNonSep_BryoDiv$waic$waic[1], fitSpatTempNonSep_BryoDiv$dic$dic[1], -sum(log(fitSpatTempNonSep_BryoDiv$cpo$cpo), na.rm = TRUE))
# 202
#[1,] 12449.78 12446.08 6225.153
# 102
#[1,] 12451.33 12447.77 6225.877
# 121, 220
#no convergence

# plot building off existing old plots from original Bryophyte diversity
plots_BryoDiv2 <- list()
plots_BryoDiv2[[1]] <- plots_BryoDiv[[1]]
plots_BryoDiv2[[2]] <- plots_BryoDiv[[2]]
plots_BryoDiv2[[3]] <- plots_BryoDiv[[3]]
plots_BryoDiv2[[4]] <- plots_BryoDiv[[4]]
plots_BryoDiv2[[5]] <- plots_BryoDiv[[5]]
plots_BryoDiv2[[6]] <- plots_BryoDiv[[6]]


# Spatial temporal model
# 202 GOF is way worse
inlaSpatTempNonSepPvalT_BryoDiv<-rep(NA, nrow=(vegINLA))
for(i in 1:nrow(vegINLA)){
  inlaSpatTempNonSepPvalT_BryoDiv[i]<-inla.pmarginal(q=vegINLA$hillBryo1[i],
                                                     marginal=fitSpatTempNonSep_BryoDiv$marginals.fitted.values[[i]])
}

postProb4_BryoDiv <- data.frame(predicted = fitSpatTempNonSep_BryoDiv$summary.fitted.values$mean[1:length(vegINLA$hillBryo1)],
                                observed = vegINLA$hillBryo1,
                                lower = fitSpatTempNonSep_BryoDiv$summary.fitted.values$`0.025quant`[1:length(vegINLA$hillBryo1)],
                                upper = fitSpatTempNonSep_BryoDiv$summary.fitted.values$`0.975quant`[1:length(vegINLA$hillBryo1)],
                                p.value = inlaSpatTempNonSepPvalT_BryoDiv) 

min4_BryoDiv <- min(postProb4_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)
max4_BryoDiv <- max(postProb4_BryoDiv[, c('upper', 'observed')], na.rm = TRUE)

plots_BryoDiv2[[7]] <- ggplot(postProb4_BryoDiv, aes(x = p.value)) + 
  geom_histogram() +
  labs(y = "Count", x = "Posterior probability") +
  theme_bw()

plots_BryoDiv2[[8]] <- ggplot(postProb4_BryoDiv, aes(x = predicted, y = observed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(y = "Observed", x = "Fitted") +
  lims(x = c(min4_BryoDiv, max4_BryoDiv), y = c(min4_BryoDiv, max4_BryoDiv)) + 
  theme_bw()

cowplot::plot_grid(plotlist = plots_BryoDiv2, nrow = 4)

