rm(list=ls())

library(tidyverse)
library(rstan)
library(RColorBrewer)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd('~/dropbox/working/taylor_host_parasite/github/')
#################################################
## DATA #########################################
#################################################
d  <- read_csv('../data/data.csv')
dd <- read.csv('../data/replicates.csv') %>% group_by(date) %>% summarise(sdI=sd(I)/sqrt(2),sdP=sd(P)/sqrt(2))
d  <- merge(d,dd,by='date')

d$sdI[is.na(d$sdI)] <- median(d$sdI,na.rm=TRUE)
d$sdP[is.na(d$sdP)] <- median(d$sdP,na.rm=TRUE)
d$sdP[d$sdP==0]     <- min(d$sdP[d$sdP!=0])

d <- d[order(d$days),]

X1 <- d %>% filter(days < 120)
X2 <- d %>% filter(days > 120 & days < 170)
X3 <- d %>% filter(days > 170)

XX <- list()
XX[[1]] <- X1
XX[[2]] <- X2
XX[[3]] <- X3

#####################################################
## STAN #############################################
#####################################################
mod_events <- stan_model(file='src/events_h2.stan')

dt       <- 0.1
ndays    <- 25
sig_perc <- 0.1
rs       <- c(0.21,0.24,0.35)
as       <- c(4.40E-7, 5.58E-8, 1.34E-8)
hs       <- c(2.96,2.79,2.46)
es       <- c(600,550,40)
ms       <- c(0.48,1.02,0.26)

initf  <- function(){list(
  H0 = runif(1,0,2.0*12080.0),
  I0 = fit_dat$I[1], 
  P0 = fit_dat$P[1],
  r = mean(rs),
  h1 = mean(hs),
  h2 = 0.5,
  a = mean(as),
  m = mean(ms),
  e = mean(es))}

tmod <- seq(0,ndays,length.out=ndays/dt)

POST_events <- list()

for(i in 1:3){
  fit_dat  <- XX[[i]]
  fit_time <- fit_dat$days - fit_dat$days[1]

  sigma_I <- rep(mean(fit_dat$sdI),nrow(fit_dat))
  sigma_P <- rep(mean(fit_dat$sdP),nrow(fit_dat))
  
  dat <- list(N=nrow(fit_dat),
            I=fit_dat$I,
            P=fit_dat$P,
            dt=dt,
            time_points=fit_time/dt + 1,
            sigma_I=sigma_I,
            sigma_P=sigma_P,
            maxt = ndays/dt,
            K = 1E7)

  mcmc <- sampling(mod_events,data=dat,init=initf)
  
  post                <- extract(mcmc)
  POST_events[[i]]    <- post
}

post <- POST_events[[1]]
flux_df_event1 <- data.frame(time_days=seq(0,ndays,length.out=dat$maxt),
                             growH=colMeans(post$growH*dt),growH_sd=apply(post$growH*dt,2,sd),
                             lossH=colMeans(post$lossH*dt),lossH_sd=apply(post$lossH*dt,2,sd),
                             hand =colMeans(post$hand*dt), hand_sd=apply(post$hand*dt,  2,sd),
                             growP=colMeans(post$growP*dt),growP_sd=apply(post$growP*dt,2,sd),
                             lossP=colMeans(post$lossP*dt),lossP_sd=apply(post$lossP*dt,2,sd))
post <- POST_events[[2]]
flux_df_event2 <- data.frame(time_days=seq(0,ndays,length.out=dat$maxt),
                             growH=colMeans(post$growH*dt),growH_sd=apply(post$growH*dt,2,sd),
                             lossH=colMeans(post$lossH*dt),lossH_sd=apply(post$lossH*dt,2,sd),
                             hand =colMeans(post$hand*dt), hand_sd=apply(post$hand*dt,  2,sd),
                             growP=colMeans(post$growP*dt),growP_sd=apply(post$growP*dt,2,sd),
                             lossP=colMeans(post$lossP*dt),lossP_sd=apply(post$lossP*dt,2,sd))
post <- POST_events[[3]]
flux_df_event3 <- data.frame(time_days=seq(0,ndays,length.out=dat$maxt),
                             growH=colMeans(post$growH*dt),growH_sd=apply(post$growH*dt,2,sd),
                             lossH=colMeans(post$lossH*dt),lossH_sd=apply(post$lossH*dt,2,sd),
                             hand =colMeans(post$hand*dt), hand_sd=apply(post$hand*dt,  2,sd),
                             growP=colMeans(post$growP*dt),growP_sd=apply(post$growP*dt,2,sd),
                             lossP=colMeans(post$lossP*dt),lossP_sd=apply(post$lossP*dt,2,sd))

write.csv(file='~/dropbox/working/taylor_host_parasite/results/flux_df_event1_h2.csv',flux_df_event1,row.names=FALSE)
write.csv(file='~/dropbox/working/taylor_host_parasite/results/flux_df_event2_h2.csv',flux_df_event2,row.names=FALSE)
write.csv(file='~/dropbox/working/taylor_host_parasite/results/flux_df_event3_h2.csv',flux_df_event3,row.names=FALSE)

##--SAVE RESULTS--#######################################
save(file='../results/POST_events.rdata',POST_events)



