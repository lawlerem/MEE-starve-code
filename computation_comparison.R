#
# starve code to fit the bird_survey data with varying n_neighbours
#

library(starve)
library(tmap)
library(rnaturalearth)
library(rnaturalearthhires)
library(raster)

missouri<- subset(ne_states(iso_a2 = "US",returnclass="sf"),
	name == "Missouri",
	select = "name")
st_crs(missouri)<- 4326
pred_locations<- raster(missouri,nrow=20,ncol=20)



# n_neighbours = 3
# start mem = 823112K; max mem = 1061008; mem usage = 237896 K = 0.238 GB

# 3.798 seconds
fit_timing<- system.time({
  bird_fit<- strv_prepare(
    cnt ~ time(year),
    bird_survey,
    distribution = "poisson",
    fit = TRUE,
    n_neighbours = 3
  )
})

# 1.691 seconds
pred_timing<- system.time({
  bird_pred3<- strv_predict(
    bird_fit,
    pred_locations,
    time = 2013:2014
  )
})


# n_neighbours = 10
# start mem = 1118604; max mem = 1667188 ; mem usage = 548584 MiB = 0.55 GB
#
# 4.528 seconds
fit_timing<- system.time({
  bird_fit<- strv_prepare(
    cnt ~ time(year),
    bird_survey,
    distribution = "poisson",
    fit = TRUE
  )
})

# 3.294 seconds
pred_timing<- system.time({
  bird_pred10<- strv_predict(
    bird_fit,
    pred_locations,
    time = 2013:2014
  )
})




# n_neighbours = 30
# start mem = 1339968; max mem = 4032484; mem usage = 2692516 MiB = 2.7 GB
#
# 11.237 seconds
fit_timing<- system.time({
  bird_fit<- strv_prepare(
    cnt ~ time(year),
    bird_survey,
    distribution = "poisson",
    fit = TRUE,
    n_neighbours = 30
  )
})

# 17.238 seconds
pred_timing<- system.time({
  bird_pred30<- strv_predict(
    bird_fit,
    pred_locations,
    time = 2013:2014
  )
})


legend_tm<- tm_layout(
  legend.outside.size = 0.15,
  panel.label.size = 0.8
)

breaks<- c(0,5,10,15,20,25,30,35)
pred3map<- tm_shape(bird_pred3["response"]) + tm_raster(title = "3 neighbours", breaks = breaks) + legend_tm
pred10map<- tm_shape(bird_pred30["response"]) + tm_raster(title = "10 neighbours", breaks = breaks) + legend_tm
pred30map<- tm_shape(bird_pred30["response"]) + tm_raster(title = "30 neighbours", breaks = breaks) + legend_tm

png("starve_preds.png", height = 1000, width = 800, res = 144)
tmap_arrange(pred3map, pred10map, pred30map, nrow = 3)
dev.off()

breaks<- c(0, 2, 4, 6, 8, 10 , 12, 14, 16)
pred3semap<- tm_shape(bird_pred3["response_se"]) + tm_raster(title = "3 neighbours", breaks = breaks) + legend_tm
pred10semap<- tm_shape(bird_pred30["response_se"]) + tm_raster(title = "10 neighbours", breaks = breaks) + legend_tm
pred30semap<- tm_shape(bird_pred30["response_se"]) + tm_raster(title = "30 neighbours", breaks = breaks) + legend_tm

png("starve_ses.png", height = 1000, width = 800, res = 144)
tmap_arrange(pred3semap, pred10semap, pred30semap, nrow = 3)
dev.off()








#
# INLA code to fit bird_survey data
#

library(INLA)
library(sf)
data("bird_survey",package="starve")

parallel<- FALSE

locations<- st_coordinates(bird_survey)
year<- bird_survey$year - min(bird_survey$year)+1
nyear<- length(unique(year))

if( parallel ) {inla.setOption(num.threads = 8)} else {inla.setOption(num.threads=1)}
# w/ parallel
# start mem = 221072 K; max mem = 875012 K; mem usage = 653940 K = 0.66 GB
# 19.710 seconds]
#
# no parallel
# start mem = 713876 K; max mem = 893524 K; mem usage = 64228 K = 0.18 GB
# 26.106 seconds
#
fit_timing<- system.time({
  mesh<- inla.mesh.2d(locations, max.n = 5*nrow(bird_survey))
  A_matrix<- inla.spde.make.A(
    mesh,
    loc=locations,
    group=year,
    n.group=nyear
  )
  spde<- inla.spde2.pcmatern(mesh = mesh, prior.range=c(0.1,0.9), prior.sigma = c(0.5,0.5))
  ar1.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))
  w<- inla.spde.make.index(name="w", n.spde=spde$n.spde, n.group=nyear)
  stack<- inla.stack(
    data = list(cnt = bird_survey$cnt),
    A = list(1,A_matrix),
    effects = list(
      Intercept = rep(1,nrow(bird_survey)),
      index = w
    )
  )
  fit<- inla(
    cnt~ -1 + Intercept + f(w, model=spde, group=w.group, control.group=list(model="ar1",hyper=ar1.spec)),
    family = "poisson",
    data = inla.stack.data(stack),
    control.predictor = list(A = inla.stack.A(stack)),
    control.compute = list(config=TRUE)
  )
})



#
# If seed=1 then simulations (inla.posterior.sample) aren't in parallel
# If seed=0 then simulations are in parallel
#
# w/ parallel : 2.992 seconds
# no parallel : 6.975 seconds
#
pred_timing<- system.time({
  nsample<- 100
  s <- inla.posterior.sample(n = nsample, fit, intern = TRUE, seed=ifelse(parallel,0,1), add.names = FALSE)
  sw<- do.call(cbind,lapply(s,`[[`,"latent"))
  sw<- sw[rownames(sw) %in% paste0("w:",seq(length(w$w))),]

  nxy<- c(20,20)
  projgrid <- inla.mesh.projector(
    mesh,
    xlim = range(locations[, 1]),
    ylim = range(locations[, 2]),
    dims = nxy
  )

  xmean<- list()
  i<- 0
  for(j in (2013:2014 - min(bird_survey$year)+1) ) {
    i<- i+1
    xmean[[i]] <- inla.mesh.project(
      projgrid,
      rowMeans(sw[w$w.group==j,])
    )
  }

  xsd<- list()
  i<- 0
  for(j in (2013:2014 - min(bird_survey$year)+1) ) {
    i<- i+1
    xmean[[i]] <- inla.mesh.project(
      projgrid,
      apply(sw[w$w.group==j,],1,sd,na.rm=TRUE)
    )
  }
})














#
# Code used for novice TMB code
#

library(TMB)
library(starve)
data(bird_survey)
library(units)
library(raster)
library(rnaturalearth)
library(rnaturalearthhires)

compile("full_separable/full_separable.cpp")
dyn.load(dynlib("full_separable/full_separable"))

#
# 131.101 seconds
#
fit_timing<- system.time({
data<- list(
  y = bird_survey$cnt,
  ds = drop_units(set_units(st_distance(bird_survey),"km")),
  dt = as.matrix(dist(bird_survey$year)),
  pred_ds = matrix(0,nrow=0,ncol=nrow(bird_survey)),
  pred_dt = matrix(0,nrow=0,ncol=nrow(bird_survey))
)
para<- list(
  mu = 0,
  working_st_sd = 0,
  working_s_range = 3,
  working_t_ar1 = 0,
  re = numeric(nrow(bird_survey)),
  pred_re = numeric(0)
)
rand<- c("re","pred_re")

obj<- MakeADFun(
  data = data,
  para = para,
  random = rand,
  DLL = "full_separable"
)
opt<- nlminb(
  obj$par,
  obj$fn,
  obj$gr
)
sdr<- sdreport(obj)
})


#
# Predictions max out memory at 36 locations x 1 year, with no std. errors
# (16 GB of RAM)
#
# 37.984 seconds
#
missouri<- subset(ne_states(iso_a2 = "US",returnclass="sf"),
	name == "Missouri",
	select = "name")
st_crs(missouri)<- 4326
pred_locations<- st_as_sf(rasterToPoints(raster(missouri,nrow=6,ncol=6),spatial=TRUE))
pred_locations$year<- 2014

pred_timing<- system.time({
  data$pred_ds = drop_units(set_units(st_distance(pred_locations,bird_survey),"km"))
  data$pred_dt = as.matrix(dist(c(pred_locations$year,bird_survey$year)))[1:nrow(pred_locations),1:nrow(bird_survey)+nrow(pred_locations)]
  para$pred_re = numeric(nrow(pred_locations))
  pred_obj<- MakeADFun(
    data = data,
    para = para,
    random = rand,
    DLL = "full_separable"
  )
  obj$fn(opt$par)
  # pred_sdr<- sdreport(pred_obj)
})
