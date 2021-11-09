rm(list=ls())
library(tidyverse)
library(RColorBrewer)

#setwd('~/dropbox/working/taylor_host_parasite/github/') #set your working directory

###############################################################################
##--FIELD EVENTS--#############################################################
###############################################################################
##--LOAD DATA--##############################
d  <- read_csv('../data/data.csv')
dd <- read.csv('../data/replicates.csv') %>% group_by(date) %>% summarise(sdI=sd(I)/sqrt(2),sdP=sd(P)/sqrt(2))
d  <- merge(d,dd,by='date')

##--ESTIMATE OF OBSERVATION ERROR--####################
d$sdI[is.na(d$sdI)] <- median(d$sdI,na.rm=TRUE)
d$sdP[is.na(d$sdP)] <- median(d$sdP,na.rm=TRUE)
d$sdP[d$sdP==0] <- min(d$sdP[d$sdP!=0])

d <- d[order(d$days),]

##--SPLIT UP INTO EVENTS--################
XX <- list()
XX[[1]] <- d %>% filter(days < 120)
XX[[2]] <- d %>% filter(days > 120 & days < 170)
XX[[3]] <- d %>% filter(days > 170)

##--TIME-STEPPING PARAMETERS FOR MODEL--########################
dt    <- 0.1
ndays <- 25
tmod  <- seq(0,ndays,length.out=ndays/dt)

##--LOAD EVENTS FITTING RESULTS FOR PLOTTING--###################
load('../results/POST_events.rdata')

###############################################################################
##--LAB EVENTS--###############################################################
###############################################################################
##--LOAD LAB DATA--########
dh <- read.csv('../data/lab_heterocapsa.csv')
ds <- read.csv('../data/lab_scrippsiella.csv')

DAT      <- list()
DAT[[1]] <- dh %>% group_by(days) %>% summarise(Hmean=mean(H), Imean=mean(I), Pmean=mean(P), sdH=sd(H)/sqrt(3), sdI=sd(I)/sqrt(3),sdP=sd(P)/sqrt(3))
DAT[[2]] <- ds %>% group_by(days) %>% summarise(Hmean=mean(H), Imean=mean(I), Pmean=mean(P), sdH=sd(H)/sqrt(3), sdI=sd(I)/sqrt(3),sdP=sd(P)/sqrt(3))

dt       <- 0.01   #divide by 100 to match time variable resolution in lab data
ndays    <- 8
tmod_lab <- seq(0,ndays,length.out=ndays/dt)

##--LOAD FITTING RESULTS FOR LAB EVENTS--##########################
load('../results/POST_lab.rdata')

#################################################################################
##--PLOTS--######################################################################
#################################################################################

##--FIELD EVENT TIME SERIES--########################
pdf('../plots/time_series_fits_events.pdf',height=6,width=8.5)
par(mfcol=c(3,3),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:3){
  post <- POST_events[[i]]
  X    <- XX[[i]]
  
  Hmod   <- colMeans(post$x[,,1])
  Hmodsd <- apply(post$x[,,1],2,sd) 
  Imod   <- colMeans(post$x[,,2])
  Imodsd <- apply(post$x[,,2],2,sd)
  Pmod   <- colMeans(post$x[,,3])
  Pmodsd <- apply(post$x[,,3],2,sd)
  
  plot(tmod,Hmod,type='l',lty=1,xlab='',ylab='',ylim=c(0,max(Hmod + 2*Hmodsd)),lwd=2)
    lines(tmod,Hmod - 2*Hmodsd,type='l',lwd=1,lty=2)
    lines(tmod,Hmod + 2*Hmodsd,type='l',lwd=1,lty=2)
    mtext(paste('Event #',i,sep=''),adj=0)
    if(i==1) mtext(side=2,'Host Abundance',line=2.5)
  plot(tmod,Imod,type='l',lty=1,xlab='',ylab='',ylim=c(0,max(Imod + 2*Imodsd, X$I)),lwd=2)
    lines(tmod,Imod - 2*Imodsd,type='l',lwd=1,lty=2)
    lines(tmod,Imod + 2*Imodsd,type='l',lwd=1,lty=2)
    points(X$days-X$days[1], X$I,pch=8,col='red')
    if(i==1) mtext(side=2,'Infected Host Abundance',line=2.5)
  plot(tmod,Pmod,type='l',lty=1,xlab='',ylab='',ylim=c(0,max(Pmod + 2*Pmodsd, X$P)),lwd=2)
    lines(tmod,Pmod - 2*Pmodsd,type='l',lwd=1,lty=2)
    lines(tmod,Pmod + 2*Pmodsd,type='l',lwd=1,lty=2)
    points(X$days-X$days[1], X$P,pch=8,col='red')
    if(i==1) mtext(side=2,'Spore Abundance',line=2.5)
}  
mtext(side=1,outer=TRUE,'Days',line=0.5)
dev.off()

##--LAB EVENT TIME SERIES--############################
pdf('../plots/time_series_fits_lab.pdf',height=6,width=8.5*(2/3))
par(mfcol=c(3,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:2){
  post <- POST_lab[[i]]
  dat  <- DAT[[i]]

  Hmod   <- colMeans(post$x[,,1])
  Hmodsd <- apply(post$x[,,1],2,sd) 
  Imod   <- colMeans(post$x[,,2])
  Imodsd <- apply(post$x[,,2],2,sd)
  Pmod   <- colMeans(post$x[,,3])
  Pmodsd <- apply(post$x[,,3],2,sd)

  plot(tmod_lab,Hmod,type='l',lty=1,xlab='',ylab='',ylim=c(0,max(Hmod + 2*Hmodsd, dat$Hmean)),lwd=2)
  lines(tmod_lab,Hmod - 2*Hmodsd,type='l',lwd=1,lty=2)
  lines(tmod_lab,Hmod + 2*Hmodsd,type='l',lwd=1,lty=2)
  if(i==1) mtext(side=2,'Host Abundance',line=2.5,cex=0.8)
  points(dat$days-dat$days[1], dat$Hmean,pch=8,col='red')
  if(i==1) mtext(expression(italic('heterocapsa')),adj=0)
  if(i==2) mtext(expression(italic('scrippsiella')),adj=0)

  plot(tmod_lab,Imod,type='l',lty=1,xlab='',ylab='',ylim=c(0,max(Imod + 2*Imodsd, dat$Imean)),lwd=2)
  lines(tmod_lab,Imod - 2*Imodsd,type='l',lwd=1,lty=2)
  lines(tmod_lab,Imod + 2*Imodsd,type='l',lwd=1,lty=2)
  points(dat$days-dat$days[1], dat$Imean,pch=8,col='red')
  if(i==1) mtext(side=2,'Infected Host Abundance',line=2.5,cex=0.8)

  plot(tmod_lab,Pmod,type='l',lty=1,xlab='',ylab='',ylim=c(0,max(Pmod + 2*Pmodsd, dat$Pmean)),lwd=2)
  lines(tmod_lab,Pmod - 2*Pmodsd,type='l',lwd=1,lty=2)
  lines(tmod_lab,Pmod + 2*Pmodsd,type='l',lwd=1,lty=2)
  points(dat$days-dat$days[1], dat$Pmean,pch=8,col='red')
  if(i==1) mtext(side=2,'Spore Abundance',line=2.5,cex=0.8)
}
mtext(side=1,outer=TRUE,'Days',line=0.5)
dev.off()


###########################################################################################
## POSTERIOR PARAMETER DISTRIBUTIONS ACROSS EVENTS ########################################
###########################################################################################

post1_events <- POST_events[[1]]
post2_events <- POST_events[[2]]
post3_events <- POST_events[[3]]
post_scrip   <- POST_lab[[2]]

cols <- brewer.pal(n = 8, name = "Dark2")
nbrks <- 60

hist2 <- function(x,xlim,col,breaks,ylim,add){
  hist(x,xlim=xlim,col=col,breaks=breaks,main='',ylim=ylim,add=add,border=TRUE)
  abline(v=quantile(x,probs=c(0.025,0.5,0.975)),lty=c(2,1,2),col=col,lwd=0.8)
}

##--PARAMETER DISTRIBUTIONS--###############################
pdf('../plots/parameter_distributions_11_05_2021.pdf',height=5,width=8)
par(mfrow=c(2,3),mar=c(2,2,2,2),oma=c(2,2,2,2))
uplim <- max(c(post1_events$r,post3_events$r,post_scrip$r))

hist2(x=post1_events$r, xlim=c(0,uplim), col=cols[1], breaks=seq(-0.1,uplim,length.out=nbrks),ylim=c(0,1500),add=FALSE)
hist2(x=post3_events$r, xlim=c(0,1), col=cols[3], breaks=seq(-0.1,uplim,length.out=nbrks),ylim=c(0,1500),add=TRUE)
  mtext(expression(italic('r')),side=1,line=2.5)
hist2(x=post_scrip$r, xlim=c(0,uplim), col=cols[2], breaks=seq(-0.1,uplim,length.out=nbrks),ylim=c(0,1500),add=TRUE)
  legend('topright',col=cols[c(1,3,2)],legend=c('Event #1','Event #3',expression(italic('Scrippsiella'))),bty='n',pch=15)

uplim <- log10(max(c(post1_events$a/(1+post1_events$a*post1_events$h1),post3_events$a/(1+post3_events$a*post3_events$h1),post_scrip$a/(1+post_scrip$a*post_scrip$h1))))
hist2(x=log10(post1_events$a/(1+post1_events$a*post1_events$h1)), xlim=c(-8.3,uplim),col=cols[1], breaks=seq(-8.3,uplim,length.out=nbrks),ylim=c(0,1500),add=FALSE)
hist2(x=log10(post3_events$a/(1+post3_events$a*post3_events$h1)), xlim=c(-8.3,uplim),col=cols[3], breaks=seq(-8.3,uplim,length.out=nbrks),ylim=c(0,1500),add=TRUE)
hist2(x=log10(post_scrip$a/(1+post_scrip$a*post_scrip$h1)), xlim=c(-8.3,uplim),col=cols[2], breaks=seq(-8.3,uplim,length.out=nbrks),ylim=c(0,1500),add=TRUE)
mtext(expression('log10{'~italic('a/(1+ah'[1]*')')~'}'),side=1,line=2.5)
  
uplim <- max(c(post1_events$h2,post3_events$h2,post_scrip$h2))
hist2(x=post1_events$h2, xlim=c(0,uplim), col=cols[1], breaks=seq(-0.1,uplim,length.out=nbrks), ylim=c(0,1500), add=FALSE)
hist2(x=post3_events$h2, xlim=c(0,uplim), col=cols[3], breaks=seq(-0.1,uplim,length.out=nbrks), ylim=c(0,1500), add=TRUE)
hist2(x=post_scrip$h2, xlim=c(0,uplim), col=cols[2], breaks=seq(-0.1,uplim,length.out=nbrks), ylim=c(0,1500), add=TRUE)
  mtext(expression(italic('h'[2])),side=1,line=2.5)

uplim <- max(c(post1_events$m,post3_events$m,post_scrip$m))
hist2(x=post1_events$m, xlim=c(0,uplim),col=cols[1],   breaks=seq(0,uplim,length.out=nbrks),ylim=c(0,1500),add=FALSE)
hist2(x=post3_events$m, xlim=c(0,uplim),col=cols[3],   breaks=seq(0,uplim,length.out=nbrks),ylim=c(0,1500),add=TRUE)
hist2(x=post_scrip$m, xlim=c(0,uplim),col=cols[2],   breaks=seq(0,uplim,length.out=nbrks),ylim=c(0,1500),add=TRUE)
  mtext(expression(italic('m')),side=1,line=2.5)

uplim <- max(c(post1_events$e,post3_events$e,post_scrip$e))
hist2(x=post1_events$e,xlim=c(0,uplim),col=cols[1],   breaks=seq(0,uplim,length.out=nbrks*2),ylim=c(0,1500),add=FALSE)
hist2(x=post3_events$e,xlim=c(0,uplim),col=cols[3],   breaks=seq(0,uplim,length.out=nbrks*2),ylim=c(0,1500),add=TRUE)
hist2(x=post_scrip$e,xlim=c(0,uplim),col=cols[2],   breaks=seq(0,uplim,length.out=nbrks*2),ylim=c(0,1500),add=TRUE)
  mtext(expression(italic('e')),side=1,line=2.5)
dev.off()


##--PARAMETER DISTRIBUTIONS FOR H1 and H2--##########################################
pdf('../plots/parameter_distributions.pdf',height=4,width=8)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(2,2,2,2),cex.axis=0.8)
  uplim <- max(c(post1_events$h1,post3_events$h1,post_scrip$h1))
  hist2(x=post1_events$h1, xlim=c(0,uplim),col=cols[1], breaks=seq(0,uplim,length.out=nbrks),ylim=c(0,600),add=FALSE)
  hist2(x=post3_events$h1, xlim=c(0,uplim),col=cols[3], breaks=seq(0,uplim,length.out=nbrks),ylim=c(0,600),add=TRUE)
  hist2(x=post_scrip$h1, xlim=c(0,uplim),col=cols[2], breaks=seq(0,uplim,length.out=nbrks),ylim=c(0,600),add=TRUE)
    mtext(expression(italic('h'[1])),side=1,line=2.5)
    legend('topleft',col=cols[c(1,3,2)],legend=c('Event #1','Event #3',expression(italic('Scrippsiella'))),bty='n',pch=15)
  
  uplim <- max(c(post1_events$h2,post3_events$h2,post_scrip$h2))
  hist2(x=post1_events$h2, xlim=c(0,uplim), col=cols[1], breaks=seq(-0.1,uplim,length.out=nbrks), ylim=c(0,600), add=FALSE)
  hist2(x=post3_events$h2, xlim=c(0,uplim), col=cols[3], breaks=seq(-0.1,uplim,length.out=nbrks), ylim=c(0,600), add=TRUE)
  hist2(x=post_scrip$h2, xlim=c(0,uplim), col=cols[2], breaks=seq(-0.1,uplim,length.out=nbrks), ylim=c(0,600), add=TRUE)
    mtext(expression(italic('h'[2])),side=1,line=2.5)
dev.off()

#############################################################################
##--TABLE OF RATES--#########################################################
#############################################################################S
df1 <- data.frame(mean = c(mean(post1_events$r), mean(post1_events$a/(1+post1_events$h1*post1_events$a)), mean(post1_events$h2), mean(post1_events$m), mean(post1_events$e)),
           sd   = c(sd(post1_events$r), sd(post1_events$a/(1+post1_events$h1*post1_events$a)), sd(post1_events$h2), sd(post1_events$m), sd(post1_events$e)))
df1$cv <- df1$sd/df1$mean
df1

df2 <- data.frame(mean = c(mean(post2_events$r), mean(post2_events$a/(1+post2_events$h1*post2_events$a)), mean(post2_events$h2), mean(post2_events$m), mean(post2_events$e)),
                  sd   = c(sd(post2_events$r), sd(post2_events$a/(1+post2_events$h1*post2_events$a)), sd(post2_events$h2), sd(post2_events$m), sd(post2_events$e)))
df2$cv <- df2$sd/df2$mean
df2

df3 <- data.frame(mean = c(mean(post3_events$r), mean(post3_events$a/(1+post3_events$h1*post3_events$a)), mean(post3_events$h2), mean(post3_events$m), mean(post3_events$e)),
                  sd   = c(sd(post3_events$r), sd(post3_events$a/(1+post3_events$h1*post3_events$a)), sd(post3_events$h2), sd(post3_events$m), sd(post3_events$e)))
df3$cv <- df3$sd/df3$mean
df3
