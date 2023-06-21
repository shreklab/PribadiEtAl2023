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

# Figure 4B_plot: Generate plot for Fig 4B (CV (Dispersion) from Center in filed arena)
rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 4/"
prefix = "Figure4B-C"
mdl = "CVdist"
plotjitter=0.08
limitsOverage=1.5


setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports


load(paste0(prefix,".RData"))
# reset prefix for plots:

prefix = "Figure4C"
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
data.predator = dd[dd$Condition == "predator",] # raw data for the predator.
data.control = dd[dd$Condition == "control",] # raw data for the control.

exemplar.predator = exemplarData[exemplarData$Condition == "predator",]
exemplar.control = exemplarData[exemplarData$Condition == "control",]

data.predator$x = data.predator$t + 2*plotjitter + runif(nrow(data.predator),-plotjitter,plotjitter)
exemplar.predator$x = exemplar.predator$t

data.control$x = data.control$t - 2*plotjitter + runif(nrow(data.control),-plotjitter,plotjitter)
exemplar.control$x = exemplar.control$t

### Plot limits ##
xlims=c(-1,1)+range(d$t)
xticks=sort(unique(d$t))
ylims=bestupperlim(c(1/limitsOverage,limitsOverage)*range(dd$CVdist),mod=10,digits=2)
yticks=seq(ylims[1],ylims[2],by=0.1)

#Print as PDF
pdf(paste0(prefix,".pdf"))
plot.new()
plot.window(xlim=xlims,ylim=ylims)

# predator fit shading
polygon(x=c(exemplar.predator$x,rev(exemplar.predator$x)),
        y=c(exemplar.predator$lwr,rev(exemplar.predator$upr)),
        col = adjustcolor(rep(exemplar.predator$cols,2),alpha.f = 0.5),border = NA)
# control fit shading
polygon(x=c(exemplar.control$x,rev(exemplar.control$x)),
        y=c(exemplar.control$lwr,rev(exemplar.control$upr)),
        col = adjustcolor(rep(exemplar.control$cols,2),alpha.f = 0.5),border = NA)


# predator fit line
lines(exemplar.predator$x, exemplar.predator$estimate,col=exemplar.predator$cols,lwd=2)
# control fit line
lines(exemplar.control$x, exemplar.control$estimate,col=exemplar.control$cols,lwd=2)

# predator data points
points(data.predator$x,
       data.predator$CVdist,col=data.predator$cols,pch=19)

# control data points
points(data.control$x,
       data.control$CVdist,col=data.control$cols,pch=19)

axis(1,at=xticks,cex.axis=0.6)
axis(2,at=yticks,las=2,cex.axis=0.6)
title(ylab=mdl,main=paste0(prefix))

dev.off()