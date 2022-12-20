## ----include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(knitr)
library(kableExtra)
opts_chunk$set(concordance=TRUE)
opts_chunk$set(include=FALSE)

## ----preliminaries, echo=FALSE, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------
options(
  prompt = "R> ",
	continue = "+  ",
	width = 70,
	digits = 4,
	useFancyQuotes = FALSE
)
opts.default<- options()

set.seed(648634)
library("starve")
library("raster")
library("rnaturalearth")
library("tmap")
library("stars")
library("DHARMa")
library("parallel")


## ----bird_fit, include=TRUE, message=FALSE, warning=FALSE, results="hide", cache=TRUE-------------------------------------------------------------------------------------------------------------
set.seed(30795)
bird_fit<- strv_prepare(
  cnt ~ time(year),
  bird_survey,
  distribution = "poisson",
  fit = TRUE
)


## ----include=TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bird_fit


## ----include=TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
space_parameters(bird_fit)


## ----pois_simulation, include=TRUE, message=FALSE, warning=FALSE, cache=TRUE----------------------------------------------------------------------------------------------------------------------
pois_sim<- strv_simulate(
  bird_fit,
  conditional = TRUE
)
pois_sim

## ----negbin_simulation, include=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------
negbin_sim<- bird_fit
response_distribution(negbin_sim)<- "negative binomial"
response_parameters(negbin_sim)$cnt["overdispersion", "par"]<- 4
negbin_sim<-strv_simulate(
  negbin_sim,
  conditional = TRUE
)

## ----compois_simulation, include=FALSE, cache=TRUE------------------------------------------------------------------------------------------------------------------------------------------------
compois_sim<- bird_fit
response_distribution(compois_sim)<- "compois"
response_parameters(compois_sim)$cnt["dispersion", "par"]<- 0.2
compois_sim<- strv_simulate(
  compois_sim,
  conditional = TRUE
)


## ----pois_CDF, include=TRUE, warning=F, error=F, message=F, cache=TRUE----------------------------------------------------------------------------------------------------------------------------
# Simulate 100 new datasets and combine the new simulated observations
# into a matrix where each column is a different simulation.
#
# These will be used to create the bootstrap CDF
pois_CDF<- do.call(
  cbind,
  mclapply(
    seq(100),
    function(i) {
      sim<- strv_simulate(
        bird_fit,
        conditional = TRUE
      )
      return(dat(sim)$cnt)
    },
    mc.cores = 8
  )
)


## ----bird_dharma, include=TRUE, warning=FALSE, error=FALSE, message=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------
# Compute the quantile residuals for the wren
# dataset relative to the fitted Poisson model
library(DHARMa)
bird_dharma<- createDHARMa(
  simulated = pois_CDF,
  observed = dat(bird_fit)$cnt,
  integer = TRUE
)

## ----pois_dharma, include=FALSE, cache=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------
pois_dharma<- createDHARMa(
  simulated = pois_CDF,
  observed = dat(pois_sim)$cnt,
  integer = TRUE
)

## ----negbin_dharma, include=FALSE, cache=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
negbin_dharma<- createDHARMa(
  simulated = pois_CDF,
  observed = dat(negbin_sim)$cnt,
  integer = TRUE
)

## ----compois_dharma, include=FALSE, cache=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------
compois_dharma<- createDHARMa(
  simulated = pois_CDF,
  observed = dat(compois_sim)$cnt,
  integer = TRUE
)


## ----qqplot, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bird_resid<- sort(bird_dharma$scaledResiduals)
underdispersed_resid<- sort(compois_dharma$scaledResiduals)
true_resid<- sort(pois_dharma$scaledResiduals)
overdispersed_resid<- sort(negbin_dharma$scaledResiduals)
qqunif<- seq(0, 1, length = length(bird_resid))

png("qqplot.png", width = 600, height = 600, pointsize = 20)
plot(x = c(0,1),
     y = c(0,1),
     main = "Q-Q Plot",
     ylab = "Observed",
     xlab = "Expected (Relative to Poisson)",
     type = "n")
lines(x = qqunif, y = bird_resid, lty = 2)
lines(x = qqunif, y = underdispersed_resid, col = "blue")
lines(x = qqunif, y = true_resid)
lines(x = qqunif, y = overdispersed_resid, col = "red")
legend(
  "topleft",
  legend = c("Bird Data", "Under-dispersed", "Poisson", "Over-dispersed"),
  col = c("black", "blue", "black", "red"),
  lty = c(2, 1, 1, 1)
)
dev.off()


## ----compois, include=TRUE, message=FALSE, warning=FALSE, cache=TRUE------------------------------------------------------------------------------------------------------------------------------
# Copy the existing model
bird_compois<- bird_fit
# Set the new response distribution
response_distribution(bird_compois)<- "compois"
# Fit the model. `silent = TRUE' suppresses the optimizer tracing
bird_compois<- strv_fit(
  bird_compois,
  silent = TRUE
)

## ----include=TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Print the convergence message
convergence(bird_compois)
space_parameters(bird_compois)
# Print the estimated parameters for the
# Conway-Maxwell-Poisson distribution
response_parameters(bird_compois)

## ----include=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
p<- response_parameters(bird_compois)$cnt["dispersion", "par"]
CI<- with(
  response_parameters(bird_compois)$cnt,
  par[[1]] + c(-1, +1) * 1.96 * se[[1]]
)


## ----bird_raster, include=TRUE, message=FALSE, warning=FALSE, results="hide", cache=TRUE----------------------------------------------------------------------------------------------------------
missouri<- subset(
  ne_states(iso_a2 = "US", returnclass = "sf"),
	name == "Missouri",
	select = "name"
)
st_crs(missouri)<- 4326
raster_to_pred<- rasterize(
  missouri,
  raster(missouri, nrow = 20, ncol = 20),
  getCover = TRUE
)
raster_to_pred[raster_to_pred == 0]<- NA

bird_stars<- strv_predict(
  bird_fit,
  raster_to_pred,
  time = 1994:2018
)

## ----include=TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bird_stars


## ----include=TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tmap)
missouri_tm<- tm_shape(missouri, is.master = TRUE) + tm_borders()
count_tm<- tm_shape(bird_survey) +
  tm_dots(col = "cnt", title = "Count", size = 0.7) +
  tm_facets(by = "year", free.coords = FALSE)
intensity_tm<- tm_shape(bird_stars["response"]) + tm_raster(title = "Intensity")
stderr_tm<- tm_shape(bird_stars["response_se"]) + tm_raster(title = "Std. Error")

legend_tm<- tm_layout(
  legend.outside.size = 0.15,
  panel.label.size = 1.4
)

count_map<- count_tm + missouri_tm + legend_tm
intensity_map<- intensity_tm + missouri_tm + legend_tm
stderr_map<- stderr_tm + missouri_tm + legend_tm


## ----include=TRUE, echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
  bird_pars<- rbind(
    time_parameters(bird_fit)$cnt[, c("par", "se")],
    space_parameters(bird_fit)$cnt[c("sd", "range"), c("par", "se")]
  )
  rownames(bird_pars)<- c(
    "(log) global mean $\\mu$",
    "AR(1) correlation $\\phi$",
    "temporal std. dev. $\\sigma$",
    "spatial std. dev. $\\tau$",
    "spatial range $\\rho$"
  )
  colnames(bird_pars)<- c("Estimate", "Standard Error")
  knitr::kable(bird_pars, "latex", escape = FALSE)


## ----echo=FALSE, include=TRUE, fig.height=5.5, fig.width=7, res=144-------------------------------------------------------------------------------------------------------------------------------
    count_map


## ----echo=FALSE, include=TRUE, fig.height=5, fig.width=7, res=144---------------------------------------------------------------------------------------------------------------------------------
    intensity_map


## ----echo=FALSE, include=TRUE, fig.height=5, fig.width=7, res=144---------------------------------------------------------------------------------------------------------------------------------
    stderr_map


## ----include=TRUE, eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------
## stack<- inla.stack(
##   data = list(cnt = bird_survey$cnt),
##   A = list(1, A_matrix),
##   effects = list(
##     Intercept = rep(1, nrow(bird_survey)),
##     index = w
##   )
## )
## fit<- inla(
##   cnt ~ -1 + Intercept +
##     f(w, model = spde, group = w.group, control.group = list(model = "ar1", hyper = ar1.spec)),
##   family = "poisson",
##   data = inla.stack.data(stack),
##   control.predictor = list(A = inla.stack.A(stack)),
##   control.compute = list(config = TRUE)
## )

