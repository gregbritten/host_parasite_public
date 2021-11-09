rm(list=ls())

library(tidyverse)
library(rstan)
library(RColorBrewer)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#setwd('') #set your working directory
#################################################
## DATA #########################################
#################################################
dh <- read.csv('../data/lab_heterocapsa.csv')
ds <- read.csv('../data/lab_scrippsiella.csv')
ddh <- dh %>% group_by(days) %>% summarise(Hmean=mean(H), Imean=mean(I), Pmean=mean(P), sdH=sd(H)/sqrt(3), sdI=sd(I)/sqrt(3),sdP=sd(P)/sqrt(3))
dds <- ds %>% group_by(days) %>% summarise(Hmean=mean(H), Imean=mean(I), Pmean=mean(P), sdH=sd(H)/sqrt(3), sdI=sd(I)/sqrt(3),sdP=sd(P)/sqrt(3))

DAT      <- list()
DAT[[1]] <- ddh
DAT[[2]] <- dds
#####################################################
## STAN #############################################
#####################################################
mod    <- stan_model(file='src/lab_h2.stan')

dt    <- 0.01   #divide by 100 to match time variable resolution in lab data
ndays <- 8

rs <- c(0.21,0.24,0.35)
as <- c(4.40E-7, 5.58E-8, 1.34E-8)
hs <- c(2.96,2.79,2.46)
es <- c(600,550,40)
ms <- c(0.48,1.02,0.26)

initf  <- function(){list(
  H0 = fit_dat$Hmean[1],
  I0 = fit_dat$Imean[1], 
  P0 = fit_dat$Pmean[1],
  r = mean(rs),
  h1 = mean(hs),
  h2 = 0.5,
  a = mean(as),
  m = mean(ms),
  e = mean(es))}

tmod <- seq(0,ndays,length.out=ndays/dt)

##--SET DATASET--###########################
POST_lab <- list()
for(i in 1:2){
  fit_dat  <- DAT[[i]]   #units of 1/100 day to match lab time resolution
  fit_time <- (fit_dat$days - fit_dat$days[1])*100
  
  sigma_H <- rep(mean(fit_dat$sdH),nrow(fit_dat))
  sigma_I <- rep(mean(fit_dat$sdI),nrow(fit_dat))
  sigma_P <- rep(mean(fit_dat$sdP),nrow(fit_dat))

  dat <- list(N=nrow(fit_dat),
              H=fit_dat$Hmean,
              I=fit_dat$Imean,
              P=fit_dat$Pmean,
              dt=dt,
              time_points=as.integer(fit_time + 1),
              sigma_H=sigma_H,
              sigma_I=sigma_I,
              sigma_P=sigma_P,
              maxt = ndays/dt,
              K = 1E7)

  mcmc    <- sampling(mod,data=dat,init=initf)
  post    <- extract(mcmc)

  POST_lab[[i]]    <- post
}

save(file='../results/POST_lab.rdata',POST_lab)


##--TABLES--#############################################
df1 <- data.frame(mean = c(mean(POST_lab[[1]]$r), mean(POST_lab[[1]]$a/(1+POST_lab[[1]]$h1*POST_lab[[1]]$a)), mean(POST_lab[[1]]$h2), mean(POST_lab[[1]]$m), mean(POST_lab[[1]]$e)),
                  sd   = c(sd(POST_lab[[1]]$r), sd(POST_lab[[1]]$a/(1+POST_lab[[1]]$h1*POST_lab[[1]]$a)), sd(POST_lab[[1]]$h2), sd(POST_lab[[1]]$m), sd(POST_lab[[1]]$e)))
df1$cv <- df1$sd/df1$mean
df1

df2 <- data.frame(mean = c(mean(POST_lab[[2]]$r), mean(POST_lab[[2]]$a/(1+POST_lab[[2]]$h1*POST_lab[[2]]$a)), mean(POST_lab[[2]]$h2), mean(POST_lab[[2]]$m), mean(POST_lab[[2]]$e)),
                  sd   = c(sd(POST_lab[[2]]$r), sd(POST_lab[[2]]$a/(1+POST_lab[[2]]$h1*POST_lab[[2]]$a)), sd(POST_lab[[2]]$h2), sd(POST_lab[[2]]$m), sd(POST_lab[[2]]$e)))
df2$cv <- df2$sd/df2$mean
df2

##--FLUXES--#############################################
#post <- POST[[1]]
post <- POST_lab[[1]]
flux_df_heterocapsa <- data.frame(time_days=seq(0,ndays,length.out=dat$maxt),
                             growH=colMeans(post$growH*dt),growH_sd=apply(post$growH*dt,2,sd),
                             lossH=colMeans(post$lossH*dt),lossH_sd=apply(post$lossH*dt,2,sd),
                             hand =colMeans(post$hand*dt), hand_sd=apply(post$hand*dt,2,sd),
                             growP=colMeans(post$growP*dt),growP_sd=apply(post$growP*dt,2,sd),
                             lossP=colMeans(post$lossP*dt),lossP_sd=apply(post$lossP*dt,2,sd))

#post <- POST[[2]]
post <- POST_lab[[2]]
flux_df_scrippsiella <- data.frame(time_days=seq(0,ndays,length.out=dat$maxt),
                                  growH=colMeans(post$growH*dt),growH_sd=apply(post$growH*dt,2,sd),
                                  lossH=colMeans(post$lossH*dt),lossH_sd=apply(post$lossH*dt,2,sd),
                                  hand =colMeans(post$hand*dt), hand_sd=apply(post$hand*dt,2,sd),
                                  growP=colMeans(post$growP*dt),growP_sd=apply(post$growP*dt,2,sd),
                                  lossP=colMeans(post$lossP*dt),lossP_sd=apply(post$lossP*dt,2,sd))

write.csv(file='~/dropbox/working/taylor_host_parasite/results/flux_df_heterocapsa_h2.csv',flux_df_heterocapsa,row.names=FALSE)
write.csv(file='~/dropbox/working/taylor_host_parasite/results/flux_df_scrippsiella_h2.csv',flux_df_scrippsiella,row.names=FALSE)

