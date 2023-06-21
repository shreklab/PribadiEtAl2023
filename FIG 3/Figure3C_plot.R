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

# Figure 3C_plot: Generate plot for Fig 3C (POff in smeared lawns or predator)
rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 3/"
prefix = "Figure3B-C"
mdl = "pOff"
plotjitter=0.08


setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports


load(paste0(prefix,".RData"))
# reset prefix for plots:

prefix = "Figure3C"
# atab: anova tab
# d: raw data
# dd: processed data <- data points in scatter.
# exemplarData: <- means & confidence estimates
# ph: posthoc comparisons

# All we really need is: dd for data points, and exemplardata

cls=c("#5CB3E5","#818181","#D36228")
names(cls)=levels(dd$Condition)
# This column is called Interaction for robustness with other Figures, but 
# this particular dataset is a single factor term.
dd$x = as.numeric(interaction(dd$Condition)) # just for robustness across code, other figures have interaction term labels
dd$cols=cls[dd$x]
exemplarData = exemplarData[[mdl]] # overwrite with just the table we need
exemplarData$x = 1:length(levels(interaction(dd$Condition))) #"interaction" will already be the rownames


xlims=c(-1,1)+range(dd$x)
xticks=exemplarData$x
xlabs=levels(interaction(dd$Condition))
ylims=c(0,100)
yticks=seq(0,100,by=20)

#Print as EPS & as PDF, one is always a little better than the other
setEPS()
postscript(paste0(prefix,".eps"))
plot.new()
plot.window(xlim=xlims,ylim=ylims)
points(dd$x+runif(nrow(dd),min=-plotjitter,max=plotjitter),
       dd$pct.off,col=dd$cols,pch=19)
axis(1,at=xticks,labels=xlabs,cex.axis=0.6)
axis(2,at=yticks,las=2,cex.axis=0.6)
title(ylab=mdl,main=prefix)
for(i in 1:nrow(exemplarData)){
  # Estimate point
  points(exemplarData[i,"x"]+2*plotjitter,
         exemplarData[i,"Prob"]*100,pch=19)
  # Up and Down CI Line
  lines(rep(exemplarData[i,"x"],2)+2*plotjitter,
        c(exemplarData[i,"Prob.lwr"],exemplarData[i,"Prob.upr"])*100,lwd=2)
  lines(rep(exemplarData[i,"x"],2)+2*plotjitter+c(-1,1)*plotjitter,
        rep(exemplarData[i,"Prob.lwr"],2)*100,lwd=2)
  lines(rep(exemplarData[i,"x"],2)+2*plotjitter+c(-1,1)*plotjitter,
        rep(exemplarData[i,"Prob.upr"],2)*100,lwd=2)
}
dev.off()

pdf(paste0(prefix,".pdf"))
plot.new()
plot.window(xlim=xlims,ylim=ylims)
points(dd$x+runif(nrow(dd),min=-plotjitter,max=plotjitter),
       dd$pct.off,col=dd$cols,pch=19)
axis(1,at=xticks,labels=xlabs,cex.axis=0.6)
axis(2,at=yticks,las=2,cex.axis=0.6)
title(ylab=mdl,main=prefix)
for(i in 1:nrow(exemplarData)){
  # Estimate point
  points(exemplarData[i,"x"]+2*plotjitter,
         exemplarData[i,"Prob"]*100,pch=19)
  # Up and Down CI Line
  lines(rep(exemplarData[i,"x"],2)+2*plotjitter,
        c(exemplarData[i,"Prob.lwr"],exemplarData[i,"Prob.upr"])*100,lwd=2)
  lines(rep(exemplarData[i,"x"],2)+2*plotjitter+c(-1,1)*plotjitter,
        rep(exemplarData[i,"Prob.lwr"],2)*100,lwd=2)
  lines(rep(exemplarData[i,"x"],2)+2*plotjitter+c(-1,1)*plotjitter,
        rep(exemplarData[i,"Prob.upr"],2)*100,lwd=2)
}
dev.off()