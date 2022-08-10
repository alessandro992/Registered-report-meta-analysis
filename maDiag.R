# Outlier diagnostics -----------------------------------------------------

#'# For "being in nature" studies

#+eval = FALSE
# Initial outlier diagnostics
# Univariate MA
maUniNature <- dataNature %$% rma(yi = yi, vi = vi, method = "REML", slab = result)

#+eval = FALSE
# MA diagnostics
baujat(maUniNature, symbol = "slab")

#+eval = FALSE
#fit FE model to all possible subsets
goshPlotNature <- gosh(maUniNature, progbar = TRUE, subsets = 1000, parallel = "multicore")
plot(goshPlotNature, out = 45, breaks = 50) # Testing the influence of single outliers

#+eval = FALSE
# Influence diagnostics
infNature <- influence(maUniNature, progbar = T)

#+eval = FALSE
### Plot the influence diagnostics
plot(infNature)

#'# For "emotional social support" studies

#+eval = FALSE
# Initial outlier diagnostics
# Univariate MA
maUniSocial <- dataSocial %$% rma(yi = yi, vi = vi, method = "REML", slab = result)

#+eval = FALSE
# MA diagnostics
baujat(maUniSocial, symbol = "slab")

#+eval = FALSE
#fit FE model to all possible subsets
goshPlotSocial <- gosh(maUniSocial, progbar = TRUE, subsets = 1000, parallel = "multicore")
plot(goshPlotSocial, out = 12, breaks = 50) # Testing the influence of single outliers

#+eval = FALSE
# Influence diagnostics
infSocial <- influence(maUniSocial, progbar = T)

#+eval = FALSE
### Plot the influence diagnostics
plot(infSocial)

#+eval = TRUE
# Outlier removal in case of a need
# Excluding improbably big effect sizes or ES with improbably small SE, i.e. excerting a big influence on the MA model due to combination of huge ES and small variance.
# Sensitivity analysis with the outlying ESs included will be reported as well.
# dat[c(),] <- NA

# Missing data ------------------------------------------------------------

# For being in nature studies
dat %>% filter(strategy == 1) %$% table(is.na(.$yi))
dat %>% filter(strategy == 1) %>% missmap(rank.order = TRUE, margins = c(10, 0), legend = F)

# For social support studies
dat %>% filter(strategy == 2) %$% table(is.na(.$yi))
dat %>% filter(strategy == 2) %>% missmap(rank.order = TRUE, margins = c(10, 0), legend = F)

#'### Percentage of missing data overall
# dat %>% filter(strategy == 2) %$% paste(round(sum(is.na(.))/prod(dim(.))*100, 3), "%", sep = "") # insert column numbers