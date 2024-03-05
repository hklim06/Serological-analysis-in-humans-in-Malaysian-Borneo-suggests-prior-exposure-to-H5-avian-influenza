###################################################################################################################
############################## Example models of species distributions and serology results ################################
###################################################################################################################
rm(list=ls())
library(INLA)
library(fields)
library(brinla)
library(ggplot2)
library(sp)
library(tidyverse)
library(bestglm)
library(MASS)
library(wesanderson)
library(sjPlot)
library(ggstatsplot)
library(ape)
library(lme4)
library(sp)
library(rgdal)
library(plyr)
library(spdep)
library(sf)

## Import data
# R:  ---> shorebird contact
# Z: --> poultry ownership
# Y:  --> H5N1 ELISA responses

## Import your input file
df <- read.csv("/users/R/input_file.csv")

## Load prediction data
pr <- read.csv("/users/R/prediction_layer.csv")
## Prepare covariates
pr$shorebirds <- NA
pr$shorebirds_n <- 1
pr$hr_per_house <- NA
pr$sampled_per_house <- 1
pr$poultry_n <- 1 
pr$poultry_nn <- NA 
pr$link <- 1
colnames(pr)

## Variable selection - stepAIC
mf <- glm(cbind(shorebirds, shorebirds_n)~ ., family="binomial", data=df)
summary(mf)
ms <- mf %>% stepAIC(direction = "both")
logLik(ms)
summary(ms)

## Assess spatial correlation with Moran's I
morans_shorebird_data <- df[c(24, 25, 1:23)]
morans_shorebird_data <- morans_shorebird_data %>% drop_na(shorebirds)
morans_shorebird_data <- morans_shorebird_data %>% drop_na(shorebirds_n)
inc.dists <- as.matrix(dist(cbind(morans_shorebird_data$GPS_XSS_X, morans_shorebird_data$GPS_XSS_Y)))
inc.dists.inv <- 1/inc.dists
diag(inc.dists.inv) <- 0
inc.dists.inv[1:5, 1:5]
inc.dists.inv[is.infinite(inc.dists.inv)] <- 0
morans_shorebird_data$output <- (morans_shorebird_data$shorebirds)/(morans_shorebird_data$shorebirds_n)
Moran.I(morans_shorebird_data$output, inc.dists.inv, scaled = FALSE, na.rm = TRUE)

############################## Model preparation ############################## 
z.n <- df$poultry_n 
z <- df$poultry_nn 
r <- df$shorebirds
r.n <- df$shorebirds_n
y <- df$hr_per_house
y.n <- df$sampled_per_house

## Set intercepts
z.b0 <- rep(1, length(z)) 
r.b0 <- rep(1, length(r))
y.b0 <- rep(1, length(y))
            
## Mean-centre and scale covariates
en_df <- data.frame(scale(df[3:23]), center=TRUE, scale=TRUE) 
en_df$center <- NULL; en_df$scale <- NULL
en <- en_df

## Extract covariates
z.ndvi <- en$ndvi
r.ndvi <- en$ndvi
y.ndvi <- en$ndvi

z.roads <- en$roads
r.roads <- en$roads
y.roads <- en$roads

z.pop <- en$pop
r.pop <- en$pop
y.pop <- en$pop

z.slp <- en$slp
r.slp <- en$slp
y.slp <- en$slp

z.asp <- en$asp
r.asp <- en$asp
y.asp <- en$asp

z.sea <- en$sea
r.sea <- en$sea
y.sea <- en$sea

z.avgt <- en$avgt
r.avgt <- en$avgt
y.avgt <- en$avgt

z.mdr <- en$mdr
r.mdr <- en$mdr
y.mdr <- en$mdr

z.maxt <- en$maxt
r.maxt <- en$maxt
y.maxt <- en$maxt

z.mint <- en$mint
r.mint <- en$mint
y.mint <- en$mint

z.bf <- en$bf
r.bf <- en$bf
y.bf <- en$bf

z.mg <- en$mg
r.mg <- en$mg
y.mg <- en$mg

z.rb <- en$rb
r.rb <- en$rb
y.rb <- en$rb

z.ag <- en$ag
r.ag <- en$ag
y.ag <- en$ag

z.op <- en$op
r.op <- en$op
y.op <- en$op

z.ir <- en$ir
r.ir <- en$ir
y.ir <- en$ir

z.of <- en$of
r.of <- en$of
y.of <- en$of

z.dem <- en$dem
r.dem <- en$dem
y.dem <- en$dem

z.prec <- en$prec
r.prec <- en$prec
y.prec <- en$prec

z.seas <- en$seas
r.seas <- en$seas
y.seas <- en$seas

## Scale coordinates
xmin0 <- min(df$GPS_XSS_X)
ymin0 <- min(df$GPS_XSS_Y)
df$x <- df$GPS_XSS_X - xmin0
df$y <- df$GPS_XSS_Y - ymin0
df$x <- df$x/1000
df$y <- df$y/1000

## Set domain
xmin <- min(df$x)
ymin <- min(df$y)
xmax <- max(df$x)
ymax <- max(df$y)
long <- c(xmin, xmin, xmax, xmax)
lat <- c(ymin, ymax, ymin, ymax)
boundary <- cbind(long,lat)

## Create mesh
pts <- cbind(df$x, df$y)
mesh1 <- inla.mesh.create.helper(points=pts, max.edge=c(1,2), cut=0.5)
plot(mesh1)
points(pts[,1], pts[,2], pch=19, cex=.5, col="red")

## Create observation matrix
A <- inla.spde.make.A(mesh1, loc=pts)
ind <- inla.spde.make.index('s', mesh1$n)

# Priors: precision = 1/ sigma ^2 - weakly informative priors
fixed.priors <- list(mean.intercept = 0, prec.intercept=1/100, mean = list(ndvi=0, roads=0, pop=0, dem=0, slp=0, asp=0, sea=0, lc=0, avgt=0, mdr=0, maxt=0, prec=0, 
                                                                           seas=0, bf=0, mg=0, rb=0, ag=0, op=0, ir=0, of=0), 
                     prec=list(ndvi=1/10, roads=1/10, pop=1/10, dem=1/10, slp=1/10, asp=1/10, sea=1/10, lc=1/10, avgt=1/10, mdr=1/10, maxt=1/10, prec=1/10, 
                               seas=1/10, bf=1/10, mg=1/10, rb=1/10, ag=1/10, op=1/10, ir=1/10, of=1/10))

spde <- inla.spde2.matern(mesh1, alpha = 2)

############################## Model shorebird contact (r) ###############################
# Data
dat.r <- list(r.b0=r.b0, r.ndvi=r.ndvi, r.roads=r.roads, r.pop=r.pop, r.dem=r.dem, r.sea=r.sea, y.lc=y.lc, y.avgt=y.avgt,
              r.maxt=r.maxt, r.prec=r.prec, r.seas=r.seas, r.bf=r.bf, r.op=r.op, r.ir=r.ir, r.of=r.of)
# Priors
fixed.priors.r <- list(mean.intercept = 0, prec.intercept=1/100,
                       mean = list(r.ndvi=0, r.roads=0, r.pop=0, r.dem=0, r.sea=0, r.lc=0, r.avgt=0,
                                   r.maxt=0, r.prec=0, r.seas=0, r.bf=0, r.op=0, r.ir=0, r.of=0),
                       prec=list(r.b0 = 1/10, r.ndvi=1/10, r.roads=1/10, r.pop=1/10, r.dem=1/10, r.sea=1/10, r.lc=1/10, r.avgt=1/10,
                                 r.maxt=1/10, r.prec=1/10, r.seas=1/10, r.bf=1/10, r.op=1/10, r.ir=1/10, r.of=1/10))
## Create data stack
stk.r <- inla.stack(data=list(r=r, y=cbind(NA, r)), A=list(A,1), tag="est.r", 
                    effects=list(i.r=1:spde$n.spde, list(data.frame(dat.r))))
## Model 0: no spatial effect
form0.r <- r ~ 0 + r.b0 + r.roads + r.dem + r.sea + r.maxt + r.prec + r.bf + r.of
res0.r <- inla(form0.r, family = "binomial", data=inla.stack.data(stk.r), Ntrials=r.n,
                      control.predictor = list(A=inla.stack.A(stk.r), compute=TRUE),
                      control.fixed=fixed.priors.r, control.family=list(link="logit"),
                      control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
                      control.inla = list(strategy = "gaussian", int.strategy = "eb"))
summary(res0.r)
bri.fixed.plot(res0.r)

# Model fit
str(sdat <- inla.stack.index(stk.r, 'est.r')$data)
fitted.values0 <- res0.r$summary.fitted.values$mean[sdat]
observed.values <- r/r.n
plot(fitted.values0, observed.values)
qplot(df$GPS_XSS_X, df$GPS_XSS_Y, colour=fitted.values0)+scale_color_gradient(low="blue", high="red")

## AUC
ROC_auc <- performance(prediction(fitted.values0,observed.values),"auc")
ROC_auc@y.values[[1]] # AUC


## Model 1: spatial effect
form1.r <- r ~ 0 + r.b0 + r.roads + r.dem + r.sea + r.maxt + r.prec + r.bf + r.of + f(i.r, model=spde)
res1.r <- inla(form1.r, family = "binomial", data=inla.stack.data(stk.r), Ntrials=r.n,
               control.predictor = list(A=inla.stack.A(stk.r), compute=TRUE),
               control.fixed=fixed.priors.r, control.family=list(link="logit"),
               control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
               control.inla = list(strategy = "gaussian", int.strategy = "eb"))
summary(res1.r)
bri.fixed.plot(res1.r)

# Model fit
str(sdat <- inla.stack.index(stk.r, 'est.r')$data)
fitted.values1 <- res1.r$summary.fitted.values$mean[sdat]
observed.values <- r/r.n

## AUC
ROC_auc <- performance(prediction(fitted.values1, observed.values),"auc")
ROC_auc@y.values[[1]] # AUC

# Calculate spatial range
spde.est <- inla.spde2.result(inla = res1.r, name = "i.r", spde = spde, do.transf = TRUE)
inla.zmarginal(spde.est$marginals.range.nominal[[1]])

############################## Model H5 contact (y) ###############################
# Data
dat.y <- list(y.b0=y.b0, y.ndvi=y.ndvi, y.dem=y.dem, y.slp=y.slp, y.asp=y.asp, y.sea=y.sea, y.avgt=y.avgt,
              y.mdr=y.mdr, y.maxt=y.maxt, y.mint=y.mint, y.seas=y.seas, y.bf=y.bf, 
              y.mg=y.mg, y.ag=y.ag, y.op=y.op,
              y.ir=y.ir, y.of=y.of, y.prec=y.prec)
# Priors
fixed.priors.y <- list(mean.intercept = 0, prec.intercept=1/100,
                       mean = list(y.ndvi=0, y.dem=0, y.slp=0, y.asp=0, y.sea=0, y.avgt=0,
                                   y.mdr=0, y.maxt=0, y.mint=0, y.seas=0, y.bf=0,
                                   y.mg=0, y.ag=0, y.op=0,
                                   y.ir=0, y.of=0, y.prec=0),
                       prec=list(y.b0 = 1/10, y.ndvi=1/10, y.dem=1/10, y.slp=1/10, y.asp=1/10,
                                 y.sea=1/10, y.avgt=1/10,
                                 y.mdr=1/10, y.maxt=1/10, y.mint=1/10, y.seas=1/10, y.bf=1/10,
                                 y.mg=1/10, y.ag=1/10, y.op=1/10,
                                 y.ir=1/10, y.of=1/10, y.prec=1/10))
## Create data stack
stk.y <- inla.stack(data=list(y=y, s=cbind(NA, y)), A=list(A,1), tag="est.y", 
                    effects=list(i.s=1:spde$n.spde, list(data.frame(dat.y))))
## Model 0: no spatial effect
form0.y <- y ~ 0 + y.b0 + y.ndvi + y.dem + y.sea + y.mdr + y.mint + y.prec + y.seas + y.ir
res0.y <- inla(form0.y, family = "binomial", data=inla.stack.data(stk.y), Ntrials=y.n,
                     control.predictor = list(A=inla.stack.A(stk.y), compute=TRUE),
                     control.fixed=fixed.priors.y, control.family=list(link="logit"),
                     control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
                     control.inla = list(strategy = "gaussian", int.strategy = "eb"))
summary(res0.y)
bri.fixed.plot(res0.y)

# Model fit
str(sdat <- inla.stack.index(stk.y, 'est.y')$data)
fitted.values0 <- res0.y$summary.fitted.values$mean[sdat]
observed.values <- y/y.n

## AUC
ROC_auc <- performance(prediction(fitted.values0,observed.values),"auc")
ROC_auc@y.values[[1]] # AUC

## Model 1: spatial effect
form.y <- y ~ 0 + y.b0 + y.ndvi + y.sea + y.seas + y.ir + y.mdr + y.mint + y.dem + f(i.s, model=spde)
res.y <- inla(form.y, family = "binomial", data=inla.stack.data(stk.y), Ntrials=y.n,
              control.predictor = list(A=inla.stack.A(stk.y), compute=TRUE),
              control.fixed=fixed.priors.y, control.family=list(link="logit"),
              control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
              control.inla = list(strategy = "gaussian", int.strategy = "eb"))
summary(res.y)
bri.fixed.plot(res.y)

# Model fit
str(sdat <- inla.stack.index(stk.y, 'est.y')$data)
fitted.values1 <- res.y$summary.fitted.values$mean[sdat]
observed.values <- y/y.n

## AUC
ROC_auc <- performance(prediction(fitted.values1, observed.values),"auc")
ROC_auc@y.values[[1]] # AUC

# Calculate spatial range
spde.est <- inla.spde2.result(inla = res.y, name = "i.s", spde = spde, do.transf = TRUE)
inla.zmarginal(spde.est$marginals.range.nominal[[1]])

############################## Joint models (y+r) ############################## 
## Joint data stack
stk.ry <- inla.stack(stk.y, stk.r)

## Joint priors
fixed.priors <- list(mean.intercept = 0, prec.intercept=1/100,
                     mean = list(y.ndvi=0, y.dem=0, y.sea=0, y.mdr=0, y.mint=0, y.prec=0, y.seas=0, y.ir=0,
                                 r.roads=0, r.dem=0, r.sea=0, r.maxt=0,
                                 r.prec=0, r.bf=0, r.of=0),
                     prec=list(y.b0 = 1/10, y.ndvi=1/10, y.dem=1/10, y.sea=1/10, y.mdr=1/10, 
                               y.mint=1/10, y.prec=1/10, y.seas=1/10, y.ir=1/10, r.b0=1/10,
                               r.roads=1/10, r.dem=1/10, r.sea=1/10, r.maxt=1/10,
                               r.prec=1/10, r.bf=1/10, r.of=1/10))
## Joint non-spatial model, s and y
form.ry0 <- y ~ 0 + y.b0 + y.ndvi + y.sea + y.seas + y.ir + y.mdr + y.mint + y.dem + 
  r.b0 + r.roads + r.dem + r.sea + r.maxt + r.prec + r.bf + r.of
res.ry0 <- inla(form.ry0, family=c('binomial', 'binomial'), data=inla.stack.data(stk.ry), 
                Ntrials = c(y.n, r.n),
                control.fixed = fixed.priors, 
                control.predictor=list(A=inla.stack.A(stk.ry),compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
                control.inla = list(strategy = "gaussian"))
summary(res.ry0) 
bri.fixed.plot(res.ry0)

## Joint model
form.ry <- y ~ 0 + y.b0 + y.ndvi + y.sea + y.seas + y.ir + y.mdr + y.mint + y.dem + 
  r.b0 + r.roads + r.dem + r.sea + r.maxt + r.prec + r.bf + r.of +
  f(i.r, model=spde) + f(i.s, model=spde)
res.ry <- inla(form.ry, family=c('binomial', 'binomial'), data=inla.stack.data(stk.ry), 
               Ntrials = c(y.n, r.n),
               control.fixed = fixed.priors, 
               control.predictor=list(A=inla.stack.A(stk.ry),compute=TRUE),
               control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
               control.inla = list(strategy = "gaussian", int.strategy = "eb"))
summary(res.ry)
spde.est.ry <- inla.spde2.result(inla = res.ry, name = "i.s", spde = spde, do.transf = TRUE)
inla.zmarginal(spde.est.ry$marginals.range.nominal[[1]])

## Joint model with shared spatial component
stk.s <- inla.stack(data=list(y=y, s=cbind(y, NA)), A=list(A,1), tag="est.s",
                    effects=list(list(i.s=1:spde$n.spde, i.sc=1:spde$n.spde),
                                 list(data.frame(dat.y))))
## Joint data stack
stk.sr <- inla.stack(stk.s, stk.r)
form.ryc <-  y ~ 0 + y.b0 + y.ndvi + y.sea + y.seas + y.ir + y.mdr + y.mint + y.dem + 
  r.b0 + r.roads + r.dem + r.sea + r.maxt + r.prec + r.bf + r.of +
  f(i.r, model=spde) + f(i.s, model=spde) + f(i.sc, copy="i.s", fixed=FALSE)
res.ryc <- inla(form.ryc, family=c('binomial', 'binomial'), data=inla.stack.data(stk.sr),
                Ntrials = c(y.n, r.n),
                #list(E=NULL, E=df$pop), #what does this mean
                control.fixed = fixed.priors, 
                control.predictor=list(A=inla.stack.A(stk.sr),compute=TRUE),
                control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE),
                control.inla = list(strategy = "gaussian", int.strategy = "eb"))
summary(res.ryc) 
spde.est.ryc <- inla.spde2.result(inla = res.ryc, name = "i.s", spde = spde, do.transf = TRUE)
inla.zmarginal(spde.est.ryc$marginals.range.nominal[[1]])



