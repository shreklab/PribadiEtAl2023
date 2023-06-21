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

# Figure 1-S1Danalysis: Analysis of GFP/dsRed ratio (nlp-29 GFP experiments) over 2-6 hours

rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 1-S1/"
dataSheet = "Figure1-S1D"
rawdata = "../../SourceData/RawData.xlsx"
plannedposthoc = "../../SourceData/PostHocComparisons.xlsx"
posthocSheet = "Figure1-S1D"
alpha=0.05

library(readxl) # reads excel files
library(car)    # tests for LMs/GLMs
library(multcomp) # multiple adjustment for post hoc, for LMs & GLMs
## Multcomp provides simple framework using the multivariate t- or Z-test adjustment of LMs or GLMs/LMMs

setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports


######### READ DATA ######### 
d = as.data.frame(read_excel(rawdata,sheet=dataSheet))

######### SPECIFY FACTOR STRUCTURE #########
d$Predator = factor(
  d$Predator, 
  levels=c("Control","stab","eud-1","PS312","RS5194","JU1051M","JU1051F"
  )
)

# mean center the time variable for fitting.
d$t.fit = d$t - mean(range(d$t))

######### MODELS #########
mdls = list(
  norm_log2ratio = lm(norm_log2ratio ~ Predator*poly(t.fit,degree=2,raw=TRUE),d)
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
exemplarData = expand.grid(Predator = levels(d$Predator),t=sort(unique(d$t)))
exemplarData$t.fit = exemplarData$t - mean(range(exemplarData$t))

# Model Matrix
X=model.matrix(~Predator*poly(t.fit,degree = 2,raw=TRUE),exemplarData)
rownames(X) = levels(interaction(exemplarData$Predator,exemplarData$t,sep = "_"))
# These are the expressions for the points. For Estimations of the slopes & the quad. slopes, we need additional lines to X:
X2=model.matrix(~Predator*poly(t.fit,degree = 2,raw=TRUE),exemplarData[exemplarData$t==2,]) # t doesn't matter, just need an extracted set of predators.
X2[,!grepl("t\\.fit",colnames(X2))]=0 # 0 out intercept columns
X2.1 = X2
X2.2 = X2
# linear slope vectors
tval=2-mean(range(exemplarData$t))
tval2 = tval^2
X2.1[X2.1==tval]=1
X2.1[X2.1==tval2]=0
rownames(X2.1)=paste0(levels(exemplarData$Predator),"_linSlope")
# X2.1 is the linear slope term
X2.2[X2.2==tval]=0
X2.2[X2.2==tval2]=1
rownames(X2.2)=paste0(levels(exemplarData$Predator),"_quadSlope")
# X2.2 is the quadratic slope term
X=rbind(X,X2.1,X2.2)
# X is now combined: Pairwise comparisons at each time point
# The overal changes to the linear and quadratic slope terms.

# Planned Contrast Design
cXDesign = as.data.frame(read_excel(plannedposthoc,sheet=posthocSheet))
# Planned contrasts are pulled from the posthocSheet
# These can include individual group comparisons, as well as more complex comparisons
# such as the comparisons between 2 deltas.

cXlevels=sort(unique(cXDesign$level))
# Level 1: contrasts between rows of X
# Level 2: contrasts between rows of cX (contrasts of contrasts)
# Level 3: contrasts between rows of contrasts of contrasts
# There are no 4 way interactions in these data.
cX=vector("list",length=length(cXlevels))

for(lev in cXlevels){
  cX.tmp = cXDesign[cXDesign$level == lev,]
  if(lev==1){ #Use X
    cX[[lev]] = X[cX.tmp$row1,] - X[cX.tmp$row2,]
    if(nrow(cX.tmp)==1){
      cX[[lev]]=t(as.matrix(cX[[lev]]))
    }
    rownames(cX[[lev]])=cX.tmp$label
  }else if(lev>1){# Use previous contrast matrix
    cX[[lev]] = cX[[lev-1]][cX.tmp$row1,] - cX[[lev-1]][cX.tmp$row2,]
    if(nrow(cX.tmp)==1){
      cX[[lev]]=t(as.matrix(cX[[lev]]))
    }
    rownames(cX[[lev]])=cX.tmp$label
  }
}
cX=do.call(rbind,cX)

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
  }else{ # no sig effects, do not run posthocs.
    ph[[i]] = "No significant effects detected at the omnibus level."
  }
}

### Write PostHoc Output To Report ###
options(max.print=999999) #temp reset max.print to fit everything on report.
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
options(max.print=1000)
######### Group Descriptive Stats Report #########

# Overwrite exemplarData as a list, like the others:
# Reconfigure X for prediction:
X=model.matrix(~Predator*poly(t.fit,degree = 2,raw=TRUE),exemplarData)
rownames(X) = levels(interaction(exemplarData$Predator,exemplarData$t,sep = "_"))
XSlopes = rbind(X2.1,X2.2)
exemplarData.old=exemplarData
slopeData.old = expand.grid(Predator=levels(d$Predator),slopeParameter=c("linear","quadratic"))

exemplarData = vector("list",length=length(mdls))
names(exemplarData) = names(mdls)

# Generate slopeData for summary stats on the curvature terms.
slopeData = vector("list",length=length(mdls))
names(slopeData) = names(mdls)
# Use multcomp functions as a shortcut to generate the estimates, ses, and CIs quickly, 
# but without "p adjustment"
for(i in names(exemplarData)){
  ht = glht(mdls[[i]],linfct = X)
  sht = summary(ht,test=adjusted(type="none"))
  ciht = confint(ht,calpha=univariate_calpha())
  # For LM, [lwr - upr] based on t-intervals
  exemplarData[[i]] = data.frame(
    Predator = exemplarData.old$Predator,
    t = exemplarData.old$t,
    t.fit = exemplarData.old$t.fit,
    estimate = sht$test$coefficients,
    se = sht$test$sigma,
    lwr = ciht$confint[,"lwr"],
    upr = ciht$confint[,"upr"],
    N = as.numeric(table(interaction(d$Predator,d$t)))
    )
}


for(i in names(slopeData)){
  ht = glht(mdls[[i]],linfct = XSlopes)
  sht = summary(ht,test=adjusted(type="none"))
  ciht = confint(ht,calpha=univariate_calpha())
  # For LM, [lwr - upr] based on t-intervals
  slopeData[[i]] = data.frame(
    Predator = slopeData.old$Predator,
    slopeParameter = slopeData.old$slopeParameter,
    estimate = sht$test$coefficients,
    se = sht$test$sigma,
    lwr = ciht$confint[,"lwr"],
    upr = ciht$confint[,"upr"],
    p = sht$test$pvalues,
    p.tidy = prettyPValue(sht$test$pvalues),
    p.sig = sigmarker(sht$test$pvalues)
  )
}


### Write Descriptive Stats Output To Report ###
filename=paste0(dataSheet,"_descriptiveStatsReports.txt")
for(i in 1:length(exemplarData)){
  if(i == 1){
    capture.output(paste0("######",names(exemplarData)[i],"######"),file = filename)
    capture.output(exemplarData[[i]],file = filename,append = TRUE)
    capture.output(slopeData[[i]],file = filename, append = TRUE)
  }else{
    capture.output(paste0("######",names(exemplarData)[i],"######"),file = filename, append=TRUE)
    capture.output(exemplarData[[i]],file = filename,append = TRUE)
    capture.output(slopeData[[i]],file = filename, append = TRUE)
  }
}
### Write Output R Structure ###
filename=paste0(dataSheet,".RData")
save(d,atab,exemplarData,mdls,atab,ph,file = filename)