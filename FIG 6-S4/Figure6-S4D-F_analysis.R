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

# Figure 6-S4D-F: Bootstrap Worm Watcher Data, dat1 mutants.

rm(list=ls())  # for console work
graphics.off() # for console work

wd="~/../Dropbox/chalasanilabsync/mrieger/Manuscripts/PribadiEtAl-2022/FINALv2/Figures/FIG 6-S4/"
dataSheet = "Figure6-S4D-F"
rawdata = "../../SourceData/RawData.xlsx"
alpha=0.05
setwd(wd)
source("../alwaysLoad.R") # base set of custom functions for reports

########## Parameters for Bootstrap ###############
seedval=2023
Nsim=10^5
frameRate=15 # 15 frames/hr
loProb=alpha/2
hiProb=1-alpha/2
lastFrame=301
firstFrame=1

########### Packages #############################
library(readxl) # reads excel files
library(stringr)
#source("computeKSBootstrap.R") # Function for bootstrap analysis
# Compute KS Bootstrap Comparing two samples
computeKSBootstrap = function(formula,data,contrast,contrastvar,samplevar,alpha=0.05,Nsim=1000,seed=1000){
  vars=all.vars(formula)
  response=vars[1]
  x=vars[2]
  #(x must be uniform)
  sample1=data[data[,contrastvar]==contrast[1],]
  sample2=data[data[,contrastvar]==contrast[2],]
  
  set.seed(seed)
  #initialize output
  output=list(
    poverlap=NA,
    n1=length(unique(sample1[,samplevar])),
    n2=length(unique(sample2[,samplevar])),
    sav1=matrix(NA,nrow=length(unique(data[,x])),ncol=Nsim),
    sav2=matrix(NA,nrow=length(unique(data[,x])),ncol=Nsim),
    dalt=rep(NA,length=Nsim),
    dnull=rep(NA,length=Nsim)
  )
  # Simulation loop
  print("Running Sims")
  for(i in 1:Nsim){
    # Null samples
    s1names=sample(unique(data[,samplevar]),size=output$n1,replace=TRUE)
    s1=matrix(NA,nrow=length(unique(data[,x])),ncol=length(s1names))
    for(s in 1:length(s1names)){
      tmp=data[data[,samplevar]==s1names[s],]
      tmp=tmp[order(tmp[,x],decreasing = FALSE),]
      s1[,s]=tmp[,response]
    }
    
    s2names=sample(unique(data[,samplevar]),size=output$n2,replace=TRUE)
    s2=matrix(NA,nrow=length(unique(data[,x])),ncol=length(s2names))
    for(s in 1:length(s2names)){
      tmp=data[data[,samplevar]==s2names[s],]
      tmp=tmp[order(tmp[,x],decreasing = FALSE),]
      s2[,s]=tmp[,response]
    }
    
    s1=rowMeans(s1)
    s2=rowMeans(s2)
    output$dnull[i] = max(abs(s2-s1))
    
    #Alt Samples
    s1names=sample(unique(sample1[,samplevar]),size=output$n1,replace=TRUE)
    s1=matrix(NA,nrow=length(unique(sample1[,x])),ncol=length(s1names))
    for(s in 1:length(s1names)){
      tmp=data[data[,samplevar]==s1names[s],]
      tmp=tmp[order(tmp[,x],decreasing = FALSE),]
      s1[,s]=tmp[,response]
    }
    
    s2names=sample(unique(sample2[,samplevar]),size=output$n2,replace=TRUE)
    s2=matrix(NA,nrow=length(unique(sample2[,x])),ncol=length(s2names))
    for(s in 1:length(s2names)){
      tmp=data[data[,samplevar]==s2names[s],]
      tmp=tmp[order(tmp[,x],decreasing = FALSE),]
      s2[,s]=tmp[,response]
    }
    s1=rowMeans(s1)
    s2=rowMeans(s2)
    output$sav1[,i]=s1
    output$sav2[,i]=s2
    
    output$dalt[i] = max(abs(s2-s1))
  }
  print("Finished Running Sims")
  ## Compute pOverlap: #dnull samples that exist within the 1-alpha confidence bounds of dalt
  lwr=quantile(output$dalt,probs=alpha/2)
  upr=quantile(output$dalt,probs=(1-alpha/2))
  overlap = sum(output$dnull >= lwr & output$dnull <= upr)
  output$poverlap = overlap/Nsim
  ## Compute two-tailed p-peak: number of samples that are as extreme or greater than the average of dalt.
  pk=mean(output$dalt)
  output$ppeak=((sum(output$dnull>=abs(pk))+sum(output$dnull<=-abs(pk)))/Nsim)
  return(output)
}


######### Load Data ################################
d = as.data.frame(read_excel(rawdata,sheet=dataSheet))
d=d[d$Frame <=301,] # There are sometimes spurious additional +1 frame at the end of some recordings, trim to exactly 301 frames

d$Genotype = factor(d$Genotype,levels=c("wt","dat1"))
for( i in levels(d$Genotype)){
######### Run Bootstrap ###########################
DistPCbootstrap = computeKSBootstrap(DistCenterMm ~ Frame, data=d[d$Genotype==i,],contrast = c("control","predator"),
                                     contrastvar = "Condition",samplevar = "PlateID",Nsim = Nsim,seed = seedval, alpha=alpha)

###### Report #####################################
DistPCbootstrap.report = c(list(Nsim=Nsim),DistPCbootstrap[c("n1","n2","dnull","dalt","poverlap","ppeak")])
DistPCbootstrap.report$dalt = quantile(DistPCbootstrap.report$dalt,probs = c(loProb,hiProb))
DistPCbootstrap.report$dnull = quantile(DistPCbootstrap.report$dnull,probs = c(loProb,hiProb))

###### Write Report to Text File ##############################
capture.output(DistPCbootstrap.report,file=paste0(dataSheet,"_",i,"_BootstrapReport.txt"))

### Write Output R Structure ###
filename=paste0(dataSheet,"_",i,".RData")
save(d,DistPCbootstrap,DistPCbootstrap.report,file = filename)
}