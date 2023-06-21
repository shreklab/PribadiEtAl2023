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

# Figure 1-S1D_plot: Generate plot for Figure 1-S1D nlp-29 fluorescence over time.
rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 1-S1/"
prefix = "Figure1-S1D"
mdl = "norm_log2ratio"
plotjitter=0.08
limitsOverage=1.5


setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports


load(paste0(prefix,".RData"))
# atab: anova tab
# d: raw data
# exemplarData: <- means & confidence estimates
# ph: posthoc comparisons

# All we really need is: dd for data points, and exemplardata

cls=c("#5CB3E5","black","#E6A126","#059E74","#0774B4","#D36228","#CC79A7")
names(cls)=levels(d$Predator)
d$cols=cls[as.character(d$Predator)]
#### For Figure 1-S1D, make a series of plots that plots each curve of nlp-29 fluorescence
# w.r.t time in a predator vs. control. Each time point should be off-set by 2*jitter for the points


# Grab Exemplar Data (fits)
exemplarData = exemplarData$norm_log2ratio
exemplarData$cols = cls[as.character(exemplarData$Predator)]


for(i in levels(d$Predator)[-1]){
  data.predator = d[d$Predator == i,] # raw data for the predator.
  data.control = d[d$Predator == "Control",] # raw data for the control.
  
  exemplar.predator = exemplarData[exemplarData$Predator == i,]
  exemplar.control = exemplarData[exemplarData$Predator == "Control",]
  
  
  # Set plotting 
  # t in 1:6
  # set control to be t - 2*plotjitter and predator to be t + 2*plotjitter for the scatter plot.
  
  
  data.predator$x = data.predator$t + 2*plotjitter + runif(nrow(data.predator),-plotjitter,plotjitter)
  exemplar.predator$x = exemplar.predator$t
  
  data.control$x = data.control$t - 2*plotjitter + runif(nrow(data.control),-plotjitter,plotjitter)
  exemplar.control$x = exemplar.control$t
  
  xlims=c(-1,1)+range(d$t)
  xticks=sort(unique(d$t))
  ylims=bestupperlim(limitsOverage*range(d$norm_log2ratio),mod=10,digits=2)
  yticks=round(seq(ylims[1],ylims[2],by=0.5),digits=1)
  
  pdf(paste0(prefix,"_PredatorVsControl_",i,".pdf"))
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
         data.predator$norm_log2ratio,col=data.predator$cols,pch=19)
  
  # control data points
  points(data.control$x,
         data.control$norm_log2ratio,col=data.control$cols,pch=19)
  
  axis(1,at=xticks,cex.axis=0.6)
  axis(2,at=yticks,las=2,cex.axis=0.6)
  title(ylab=mdl,main=paste0(prefix,"|",i,"| Control"))
  
  dev.off()
}