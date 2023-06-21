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

# Figure6-S4B_plot: Generate plot for Fig 6-S4B  WormWatcher data with bootstrap intervals & sig asterisks, cat2.
rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 6-S4/"
prefix = "Figure6-S4A-C"
genotype="cat2"
frameRate=15 # 15 frames/hr
lastFrame=301
firstFrame=1
alpha=0.05
loProb=alpha/2
hiProb=1-alpha/2
limitsOverage=1.4

######### Load Data #########
setwd(wd)
load(paste0(prefix,"_",genotype,".Rdata"))
source("../alwaysLoad.R")

# reset prefix
prefix="Figure6-S4B"

######## Generate Plot #######

lawnRad=mean(d$Lawn_Radius)*d$PxScale[1]
xlim=c(0,ceiling(lastFrame*(1/frameRate)))
xtick=seq(from=xlim[1],to=xlim[2],by=4)
#ylim=bestupperlim(c(1/limitsOverage,limitsOverage)*range(d$DistCenterMm),mod=10,digits=2)
ylim=c(0,3.5)
ytick=seq(from=ylim[1],to=ylim[2],by=0.5)
t=((firstFrame:lastFrame)-1)*1/frameRate

d=d[d$Genotype == genotype,]
pdf(paste0(prefix,"_",genotype,".pdf"))
plot.new()
plot.window(xlim=xlim,ylim=ylim)
polygon(c(xlim,rev(xlim)),c(0,0,rep(lawnRad,2)),col = "#EFD99E",border = NA)
lines(xlim,rep(lawnRad,2))
#### Plot individual plate lines, adjust the alpha manually in Illustrator
for(plate in unique(d$PlateID)){
  y=d[d$PlateID==plate,]
  y=y[order(y$Frame,decreasing = FALSE),]
  if(y$Condition[1]=="predator"){
    lines(t,y$DistCenterMm,col="#FFAE97")
  }else if(y$Condition[1]=="control"){
    lines(t,y$DistCenterMm,col="#B3EDFF")
  }
}
#### Determine Conf Intervals for Plotting for each condition
lo1=apply(DistPCbootstrap$sav1,1,quantile,probs=0.025)
hi1=apply(DistPCbootstrap$sav1,1,quantile,probs=0.975)
polygon(c(t,rev(t)),c(lo1,rev(hi1)),col = adjustcolor("#5CB3E5",alpha.f = 0.5),border=NA)
lo2=apply(DistPCbootstrap$sav2,1,quantile,probs=0.025)
hi2=apply(DistPCbootstrap$sav2,1,quantile,probs=0.975)
polygon(c(t,rev(t)),c(lo2,rev(hi2)),col = adjustcolor("#D36228",alpha.f = 0.5),border=NA)

#### Plot Bootstrap Mean Estimate lines ##########
mn1=rowMeans(DistPCbootstrap$sav1)
lines(t,mn1,col="#5CB3E5",lwd=2)
mn2=rowMeans(DistPCbootstrap$sav2)
lines(t,mn2,col="#D36228",lwd=2)

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
points(tz,rep(0.5,length(tz)),pch="*",lwd=2) # 0.5 is a low value, can move these around in Illustrator.


axis(1,at=xtick)
axis(2,at=ytick,las=2)
title(main=paste0("pKS(overlap)~ \n",prettyPValue(DistPCbootstrap$poverlap),
                  " | Nsim = \n",DistPCbootstrap.report$Nsim,
                  " | %0 overlap = ",round(mean(confOverlap==0)*100,digits=2),
                  " | t first zero =", round(tz[1],digits=1), "hours"),ylab = "dist from center(mm)",xlab="hrs")
dev.off()
