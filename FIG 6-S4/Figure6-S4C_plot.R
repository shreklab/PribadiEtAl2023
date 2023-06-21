# _                                
# platform       x86_64-w64-mingw32               
# arch           x86_64                           
# os             mingw32                          
# crt            ucrt                             
# system         x86_64, mingw32                  
# status                                          
# major          4                                
# minor          2.1                              
# year           2022                             
# month          06                               
# day            23                               
# svn rev        82513                            
# language       R                                
# version.string R version 4.2.1 (2022-06-23 ucrt)
# nickname       Funny-Looking Kid

# Mike Rieger, Update 04/06/2023

# Figure6-S4C_plot: Generate plot for Fig 6-S4C  WormWatcher data P/C delta with bootstrap intervals & sig asterisks.
rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 6-S4/"
prefix = "Figure6-S4A-C"
genotype=c("wt","cat2")
frameRate=15 # 15 frames/hr
lastFrame=301
firstFrame=1
alpha=0.05
loProb=alpha/2
hiProb=1-alpha/2
limitsOverage=1.5

######### Load Data #########
setwd(wd)
delta=vector("list",length=length(genotype))
names(delta)=genotype
for(i in genotype){
  load(paste0(prefix,"_",i,".Rdata"))
  delta[[i]] = DistPCbootstrap$sav2 - DistPCbootstrap$sav1
}

source("../alwaysLoad.R")

# Reset Prefix
prefix = "Figure6-S4C"
######## Generate Plot #######

xlim=c(0,ceiling(lastFrame*(1/frameRate)))
xtick=seq(from=xlim[1],to=xlim[2],by=4)
ylim=bestupperlim(limitsOverage*range(as.numeric(delta[[1]]),as.numeric(delta[[2]])),mod=10,digits=1)
ytick=seq(from=ylim[1],to=ylim[2],by=0.5)
t=((firstFrame:lastFrame)-1)*1/frameRate


pdf(paste0(prefix,"_deltas.pdf"))
plot.new()
plot.window(xlim=xlim,ylim=ylim)

lines(xlim,rep(0,2))

#### Determine Conf Intervals for Plotting for each condition
lo1=apply(delta[[1]],1,quantile,probs=alpha/2)
hi1=apply(delta[[1]],1,quantile,probs=(1-alpha/2))
polygon(c(t,rev(t)),c(lo1,rev(hi1)),col = adjustcolor("#433d6a",alpha.f = 0.5),border=NA)
lo2=apply(delta[[2]],1,quantile,probs=alpha/2)
hi2=apply(delta[[2]],1,quantile,probs=(1-alpha/2))
polygon(c(t,rev(t)),c(lo2,rev(hi2)),col = adjustcolor("#cf9b50",alpha.f = 0.5),border=NA)



#### Plot Bootstrap Mean Estimate lines ##########
mn1=rowMeans(delta[[1]])
lines(t,mn1,col="#433d6a",lwd=2)
mn2=rowMeans(delta[[2]])
lines(t,mn2,col="#cf9b50",lwd=2)

###### Determine Confidence Interval Overlap == 0 #######
confOverlap=rep(NA,length(t))
for(i in 1:length(confOverlap)){
  x2=c(lo2[i],hi2[i])
  x1=c(lo1[i],hi1[i])
  confOverlap[i]=computeConfIntOverlap(x1,x2)
}
z=which(confOverlap==0)
tz=t[which(confOverlap==0)]
###### Plot asterisks ###############################
points(tz,rep(-0.5,length(tz)),pch="*",lwd=2) # 0.5 is a low value, can move these around in Illustrator.


axis(1,at=xtick)
axis(2,at=ytick,las=2)
title(main=paste0(" | Nsim = \n",DistPCbootstrap.report$Nsim,
                  " | %0 overlap = ",round(mean(confOverlap==0)*100,digits=2)),ylab = "Delta",xlab="hrs")
dev.off()
