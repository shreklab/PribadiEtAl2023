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

# Figure 9D_plot: Generate plot for Fig 9D Poff in control vs. predator in different dop receptor triple mutants.
rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 9/"
prefix = "Figure9D-F"
mdl = "pOff"
plotjitter=0.08
limitsOverage=1.5


setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports

load(paste0(prefix,".RData"))

# reset prefix for output
prefix = "Figure9D"
# atab: anova tab
# d: raw data
# dd: processed data <- data points in scatter.
# exemplarData: <- means & confidence estimates
# ph: posthoc comparisons

# All we really need is: dd for data points, and exemplardata

cls=c("#5CB3E5","#D36228")
names(cls)=levels(dd$Condition)
dd$cols = cls[as.character(dd$Condition)]

exemplarData = exemplarData[[mdl]] # overwrite with just the table we need
exemplarData$cols = cls[as.character(exemplarData$Condition)]



################# Make Plots of Condition X Time #########################
### 4D displays the STREAKY condition "data.predator" vs. the NORMAL condition in control (mock) arenas.
data.predator = dd[dd$Condition == "predator",] # raw data for the predator.
data.control = dd[dd$Condition == "control",] # raw data for the control.

exemplar.predator = exemplarData[exemplarData$Condition == "predator",]
exemplar.control = exemplarData[exemplarData$Condition == "control",]

data.predator$x = as.numeric(data.predator$Genotype) + 2*plotjitter + runif(nrow(data.predator),-plotjitter,plotjitter)
exemplar.predator$x = as.numeric(exemplar.predator$Genotype)

data.control$x = as.numeric(data.control$Genotype) - 2*plotjitter + runif(nrow(data.control),-plotjitter,plotjitter)
exemplar.control$x = as.numeric(exemplar.control$Genotype)

### Plot limits ##
xlims=c(-1,1)+range(as.numeric(dd$Genotype))
xticks=sort(unique(as.numeric(dd$Genotype)))
xlabs=levels(dd$Genotype)
ylims=c(0,100)
yticks=seq(ylims[1],ylims[2],by=20)

#Print as PDF
pdf(paste0(prefix,".pdf"))
plot.new()
plot.window(xlim=xlims,ylim=ylims)

# predator data points
points(data.predator$x,
       data.predator$pct.off,col=data.predator$cols,pch=19)

# control data points
points(data.control$x,
       data.control$pct.off,col=data.control$cols,pch=19)

axis(1,at=xticks,labels=xlabs,cex.axis=0.6)
axis(2,at=yticks,las=2,cex.axis=0.6)
title(ylab=mdl,main=prefix)

for(i in 1:nrow(exemplar.control)){
  # Estimate point
  points(exemplar.control[i,"x"]-plotjitter,
         exemplar.control[i,"Prob"]*100,pch=19)
  # Up and Down CI Line
  lines(rep(exemplar.control[i,"x"],2)-plotjitter,
        c(exemplar.control[i,"Prob.lwr"],exemplar.control[i,"Prob.upr"])*100,lwd=2)
  lines(rep(exemplar.control[i,"x"],2)-plotjitter+c(-1,1)*plotjitter,
        rep(exemplar.control[i,"Prob.lwr"],2)*100,lwd=2)
  lines(rep(exemplar.control[i,"x"],2)-plotjitter+c(-1,1)*plotjitter,
        rep(exemplar.control[i,"Prob.upr"],2)*100,lwd=2)
}


for(i in 1:nrow(exemplar.predator)){
  # Estimate point
  points(exemplar.predator[i,"x"]+4*plotjitter,
         exemplar.predator[i,"Prob"]*100,pch=19)
  # Up and Down CI Line
  lines(rep(exemplar.predator[i,"x"],2)+4*plotjitter,
        c(exemplar.predator[i,"Prob.lwr"],exemplar.predator[i,"Prob.upr"])*100,lwd=2)
  lines(rep(exemplar.predator[i,"x"],2)+4*plotjitter+c(-1,1)*plotjitter,
        rep(exemplar.predator[i,"Prob.lwr"],2)*100,lwd=2)
  lines(rep(exemplar.predator[i,"x"],2)+4*plotjitter+c(-1,1)*plotjitter,
        rep(exemplar.predator[i,"Prob.upr"],2)*100,lwd=2)
}

dev.off()