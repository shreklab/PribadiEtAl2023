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

# Figure 1-S3: Analysis of egg distributions in presence of either control or P. uniformis conditioned media

rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 1-S3/"
dataSheet = "Figure1-S3"
rawdata = "../../SourceData/RawData.xlsx"
#plannedposthoc = "../../SourceData/PostHocComparisons.xlsx"
# this dataset only has two groups, no complex post hoc comparisons are needed.
#posthocSheet = "Figure1-S2A-B"
alpha=0.05

library(readxl) # reads excel files
library(car)    # tests for LMs/GLMs
library(multcomp) # multiple adjustment for post hoc, for LMs & GLMs
library(stringr)
## Multcomp provides simple framework using the multivariate t- or Z-test adjustment of LMs or GLMs/LMMs

setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports


######### READ DATA ######### 
d = as.data.frame(read_excel(rawdata,sheet=dataSheet))

######### SPECIFY FACTOR STRUCTURE #########
d$Condition = factor(
  d$Condition, 
  levels=c("control","predator")
)

######### DATA PROCESSING ######### 
## Binomial Aggregate Data: # eggs off, on, and total by arena
dd = aggregate(off ~ Arena*Condition,d,sum) # compute the number of eggs off lawn by arena
dd$on = aggregate(1-off ~ Arena*Condition,d,sum)$`1 - off` # compute the number of eggs on lawn by arena
dd$total = dd$off + dd$on
dd$pct.off = (dd$off/dd$total)*100

## Average Eggs per C. elegans in arena
dd$nCel = 6
dd$eggsPerWorm = dd$total/dd$nCel

## Distributional Parameters of Egg Distance from center of arena
## Mean Distance From Edge
dd$MeanDistance = aggregate(DistCenterMm ~ Arena*Condition,d,mean)$DistCenterMm # Mean distance from center of the arena
## Lower 25% Quartile
dd$LowerDistance = aggregate(DistCenterMm ~ Arena*Condition,d,quantile,probs=0.25)$DistCenterMm
## Upper 75% Quartile
dd$UpperDistance = aggregate(DistCenterMm ~ Arena*Condition,d,quantile,probs=0.75)$DistCenterMm 
## Coefficient of Variation, Distance From Edge
dd$SDdist = aggregate(DistCenterMm ~ Arena*Condition,d,sd)$DistCenterMm # Standard Deviation
dd$CVdist = dd$SDdist/dd$MeanDistance # Coefficient of Variation

#Write Processed Data to File#
filename = paste0(dataSheet,"_processedData.csv")
write.csv(dd,file=filename,row.names = FALSE)

######### MODELS #########
mdls = list(
  pOff = glm(cbind(off,on) ~ Condition,dd,family=binomial),
  eggsPerWorm = lm(eggsPerWorm ~ Condition,dd),
  MeanDistance = lm(MeanDistance ~ Condition,dd),
  LowerDistance = lm(LowerDistance ~ Condition,dd),
  UpperDistance = lm(UpperDistance ~ Condition,dd),
  CVdist = lm(CVdist ~ Condition,dd)
)

######### Omnibus Tests ######### 
atab = vector("list",length(mdls))
names(atab) = names(mdls)
for(i in names(atab)){
  atab[[i]] = Anova(mdls[[i]]) ## car Anova performs Type II ANOVA F tests for outputs of lm() and Type II Likelihood Ratio for outputs of glm()
}
# where names include "-" symbols, a "printHypothesis()" error occurs but does not affect output.

# Generate a report from the omnibus tests and write to disk
filename=paste0(dataSheet,"_omnibusReports.txt")
for(i in 1:length(atab)){
  if(i == 1){
    capture.output(paste0("######",names(atab)[i],"######"),file = filename)
    capture.output(atab[[i]],file = filename,append = TRUE)
  }else{
    capture.output(paste0("######",names(atab)[i],"######"),file = filename, append=TRUE)
    capture.output(atab[[i]],file = filename,append = TRUE)
  }
}


######### PostHoc Comparisons Tests #########
# Data frame for descriptive stats and model matrix
exemplarData = expand.grid(Condition = levels(dd$Condition))

# Model Matrix
X=model.matrix(~Condition,exemplarData)
rownames(X) = levels(dd$Condition)
# Planned Contrast Design
# In this dataset, there are only two groups. P value should be derived from omnibus tests. However,
# pull the estimate for the predator - control with associate confidence interval:
cX = t(X["predator",]-X["control",])
rownames(cX)="predator - control"
# Compute posthocs where omnibus returns significant effects
ph=vector("list",length=length(atab))
names(ph)=names(atab)
for(i in names(ph)){
  a.tmp = atab[[i]]
  p.tmp = a.tmp[,grepl("Pr",colnames(a.tmp))]
  if(sum(p.tmp<alpha,na.rm = TRUE)>0){
    ht = glht(mdls[[i]],linfct = cX)
    sht = summary(ht)
    ciht = confint(ht)
    if(i=="pOff"){ # compute the inverse logOR for reporting
    ph[[i]]=data.frame(
      estimate = sht$test$coefficients,
      se = sht$test$sigma,
      lwr = ciht$confint[,"lwr"],
      upr = ciht$confint[,"upr"],
      z = sht$test$tstat,
      p.raw = 2*pnorm(abs(sht$test$tstat),mean=0,sd=1,lower.tail=FALSE), # two-tailed unadjusted for comparison before & after correction.
      p.adj = sht$test$pvalues,
      linear.FCOR = exp(sht$test$coefficients),
      linear.lwr = exp(ciht$confint[,"lwr"]),
      linear.upr = exp(ciht$confint[,"upr"]),
      p.adj.tidy = prettyPValue(sht$test$pvalues),
      p.sig = sigmarker(sht$test$pvalues)
    )
    }else {# report LM tests based on multivariate t-test corrections
      ph[[i]]=data.frame(
      estimate = sht$test$coefficients,
      se = sht$test$sigma,
      lwr = ciht$confint[,"lwr"],
      upr = ciht$confint[,"upr"],
      t = sht$test$tstat,
      p.raw = 2*pt(abs(sht$test$tstat),df = mdls[[i]]$df.residual,lower.tail = FALSE), # two-tailed unadjusted for comparison before & after correction.
      p.adj = sht$test$pvalues,
      p.adj.tidy = prettyPValue(sht$test$pvalues),
      p.sig = sigmarker(sht$test$pvalues)
      )
    }
  }else{ # no sig effects, do not run posthocs.
    ph[[i]] = "No significant effects detected at the omnibus level."
  }
}

### Write PostHoc Output To Report ###
filename=paste0(dataSheet,"_postHocReports.txt")
for(i in 1:length(ph)){
  if(i == 1){
    capture.output(paste0("######",names(ph)[i],"######"),file = filename)
    capture.output(ph[[i]],file = filename,append = TRUE)
  }else{
    capture.output(paste0("######",names(ph)[i],"######"),file = filename, append=TRUE)
    capture.output(ph[[i]],file = filename,append = TRUE)
  }
}

######### Group Descriptive Stats Report #########

# Overwrite exemplarData as a list, like the others:
exemplarData = vector("list",length=length(mdls))
names(exemplarData) = names(mdls)
# Use multcomp functions as a shortcut to generate the estimates, ses, and CIs quickly, 
# but without "p adjustment"
for(i in names(exemplarData)){
  ht = glht(mdls[[i]],linfct = X)
  sht = summary(ht,test=adjusted(type="none"))
  ciht = confint(ht,calpha=univariate_calpha())
  if(i=="pOff"){ # also compute the linear scale inverse logit
    # For GLM, [lwr - upr] based on z-intervals
    exemplarData[[i]] = data.frame(
      estimate = sht$test$coefficients,
      se = sht$test$sigma,
      lwr = ciht$confint[,"lwr"],
      upr = ciht$confint[,"upr"],
      Prob = invlogor(sht$test$coefficients),
      Prob.lwr = invlogor(ciht$confint[,"lwr"]),
      Prob.upr = invlogor(ciht$confint[,"upr"]),
      Narenas = as.numeric(table(dd$Condition))
    )
  }else{ # omit invlogor step
    # For LM, [lwr - upr] based on t-intervals
    exemplarData[[i]] = data.frame(
      estimate = sht$test$coefficients,
      se = sht$test$sigma,
      lwr = ciht$confint[,"lwr"],
      upr = ciht$confint[,"upr"]
    ) 
  }
}

### Write Descriptive Stats Output To Report ###
filename=paste0(dataSheet,"_descriptiveStatsReports.txt")
for(i in 1:length(exemplarData)){
  if(i == 1){
    capture.output(paste0("######",names(exemplarData)[i],"######"),file = filename)
    capture.output(exemplarData[[i]],file = filename,append = TRUE)
  }else{
    capture.output(paste0("######",names(exemplarData)[i],"######"),file = filename, append=TRUE)
    capture.output(exemplarData[[i]],file = filename,append = TRUE)
  }
}

### Write Output R Structure ###
filename=paste0(dataSheet,".RData")
save(d,dd,atab,exemplarData,mdls,atab,ph,file = filename)