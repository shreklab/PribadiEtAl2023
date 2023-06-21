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

# Figure 4F_plot: Generate plot for Fig 4F Poff in predator/normal vs. streaky/control
rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 4/"
prefix = "Figure4D-G"
mdl = "pOff"
plotjitter=0.08
limitsOverage=1.5


setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports


load(paste0(prefix,".RData"))
# reset prefix for plots:

prefix = "Figure4F"
# atab: anova tab
# d: raw data
# dd: processed data <- data points in scatter.
# exemplarData: <- means & confidence estimates
# ph: posthoc comparisons

# All we really need is: dd for data points, and exemplardata

cls=c("#5CB3E5","#728baa","#D36228","#8a3537")
names(cls)=levels(interaction(dd$Lawn,dd$Condition,sep="_"))
dd$cols = cls[as.character(interaction(dd$Lawn,dd$Condition,sep="_"))]

exemplarData = exemplarData[[mdl]] # overwrite with just the table we need
exemplarData$cols = cls[as.character(interaction(exemplarData$Lawn,exemplarData$Condition,sep="_"))]



################# Make Plots of Condition X Time #########################
### 4D displays the STREAKY condition "data.predator" vs. the NORMAL condition in control (mock) arenas.
data.predator = dd[dd$Lawn == "norm" & dd$Condition == "predator",] # raw data for the predator.
data.control = dd[dd$Lawn == "streaky" & dd$Condition == "control",] # raw data for the control.

exemplar.predator = exemplarData[exemplarData$Lawn == "norm" & exemplarData$Condition == "predator",]
exemplar.control = exemplarData[exemplarData$Lawn == "streaky" & exemplarData$Condition == "control",]

#data.predator$x = data.predator$t + 2*plotjitter + runif(nrow(data.predator),-plotjitter,plotjitter)
data.predator$x = data.predator$t + 2*plotjitter + data.predator$plotjitter
exemplar.predator$x = exemplar.predator$t

#data.control$x = data.control$t - 2*plotjitter + runif(nrow(data.control),-plotjitter,plotjitter)
data.control$x = data.control$t - 2*plotjitter + data.control$plotjitter

exemplar.control$x = exemplar.control$t

### Plot limits ##
xlims=c(-1,1)+range(d$t)
xticks=sort(unique(d$t))
ylims=c(0,100)
yticks=seq(ylims[1],ylims[2],by=20)

#Print as PDF
pdf(paste0(prefix,".pdf"))
plot.new()
plot.window(xlim=xlims,ylim=ylims)

# predator fit shading
polygon(x=c(exemplar.predator$x,rev(exemplar.predator$x)),
        y=c(exemplar.predator$Prob.lwr*100,rev(exemplar.predator$Prob.upr*100)),
        col = adjustcolor(rep(exemplar.predator$cols,2),alpha.f = 0.5),border = NA)
# control fit shading
polygon(x=c(exemplar.control$x,rev(exemplar.control$x)),
        y=c(exemplar.control$Prob.lwr*100,rev(exemplar.control$Prob.upr*100)),
        col = adjustcolor(rep(exemplar.control$cols,2),alpha.f = 0.5),border = NA)


# predator fit line
lines(exemplar.predator$x, exemplar.predator$Prob*100,col=exemplar.predator$cols,lwd=2)
# control fit line
lines(exemplar.control$x, exemplar.control$Prob*100,col=exemplar.control$cols,lwd=2)

# predator data points
points(data.predator$x,
       data.predator$pct.off,col=data.predator$cols,pch=19)

# control data points
points(data.control$x,
       data.control$pct.off,col=data.control$cols,pch=19)

axis(1,at=xticks,cex.axis=0.6)
axis(2,at=yticks,las=2,cex.axis=0.6)
title(ylab=mdl,main=paste0(prefix))

dev.off()