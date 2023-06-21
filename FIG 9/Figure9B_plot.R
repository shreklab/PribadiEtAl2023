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

# Figure 9B_plot: Generate effects contrast plot for changes to the baseline off-lawn laying probability.

rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 9/"
prefix = "Figure9A-C"
mdl = "pOff"
plotscale = "log2"
plotjitter=0.08
limitsOverage=1.25


setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports

load(paste0(prefix,".RData"))

# reset prefix for output
prefix = "Figure9B"
# atab: anova tab
# d: raw data
# dd: processed data <- data points in scatter.
# exemplarData: <- means & confidence estimates
# ph: posthoc comparisons

# All we use are the post hocs 
# In this case, we care about posthocs (fold changes) in linear scale between predator & mock at 1HR, 2HRS, and 24HRS respectively.

ph = ph[[mdl]]

phPC = ph[c("predator_wt - control_wt",
          "predator_d12 - control_d12",
          "predator_d13 - control_d13",
          "predator_d14 - control_d14",
          "predator_d23 - control_d23",
          "predator_d24 - control_d24",
          "predator_d34 - control_d34"),]

# For presentation, sort based on effects similar to WT
k = order(phPC$estimate[-1],decreasing = TRUE)

ph = ph[c("control_d12 - control_wt",
          "control_d13 - control_wt",
          "control_d14 - control_wt",
          "control_d23 - control_wt",
          "control_d24 - control_wt",
          "control_d34 - control_wt"
),]
ph = ph[k,]


if(plotscale == "linear"){
  deltas = ph[,c("linear.FCOR","linear.lwr","linear.upr")]
}else if(plotscale == "log"){
  deltas =  ph[,c("estimate","lwr","upr")]
}else if(plotscale == "log2"){
  # change base from e to 2:
  deltas =  ph[,c("estimate","lwr","upr")]
  deltas = deltas*log2(exp(1))
}


### Plot limits ##
if(plotscale == "linear"){
  xlims=c(-2,bestupperlim(limitsOverage*max(deltas),mod=5,digits = 1)) # -2 leaves space for labels.
  xticks=seq(0,xlims[2],by=2) # 1 is the fold change line.
}else if(plotscale == "log" | plotscale == "log2"){
  xlims=bestupperlim(c(-2,1)*limitsOverage+range(deltas),mod=5,digits=1)
  xticks=seq(xlims[1],xlims[2],by=0.5) # leaves space for labels, should include 0. 0 is the no change predator/control line.
}
### Y-axis is the label axis ###

ylims=c(-1,1)+range(as.numeric(dd$Genotype))
yticks=sort(unique(as.numeric(dd$Genotype)))[-1]-1 # compared to wt
#ylabs=levels(dd$Genotype)[-1] # compared to wt
ylabs=levels(dd$Genotype)[-1][k]


#Print as PDF
pdf(paste0(prefix,".pdf"))
plot.new()
plot.window(xlim=xlims,ylim=ylims)
# Confidence bounds
for(i in 1:nrow(deltas)){
  points(deltas[i,1],i,pch=19)
  lines(deltas[i,c(2,3)],rep(yticks[i],2),lwd=2)
  lines(rep(deltas[i,2],2),i+c(-plotjitter,plotjitter),lwd=2)
  lines(rep(deltas[i,3],2),i+c(-plotjitter,plotjitter),lwd=2)
}
# add the no change line
axis(1,at=xticks,cex.axis=0.6)
if(plotscale == "linear"){
  nochange=1
}else{
  nochange=0
}
lines(rep(nochange,2),ylims)
text(rep(xlims[1],length(ylabs)),yticks,labels = ylabs,adj=0)
axis(2,at=yticks,labels = rep("",length(yticks)),las=2,cex.axis=0.6)
if(plotscale=="linear"){
  xlab="Fold change in the Odds Ratio (P(off lawn)/P(on lawn)) in control conditions vs. WT"
}else if(plotscale == "log"){
  xlab="Log fold change in the Odds Ratio (P(off lawn)/P(on lawn)) in control conditions vs. WT"
}else if(plotscale == "log2"){
  xlab="Log2 fold change in the Odds Ratio (P(off lawn)/P(on lawn)) in control conditions vs. WT"
}

title(xlab=xlab,main=prefix)
box(which = "plot",lty="solid")

dev.off()