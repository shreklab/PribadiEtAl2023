source("~/../Dropbox/chalasanilabsync/mrieger/Code/salkcolorpalette.R")
pdf("salkcolorpalette_reference.pdf")
par(mfrow=c(4,4))
for(i in names(salkcol)){
  barplot(1,col=salkcol[i],yaxt="n",xlab=i)
}
par(mfrow=c(1,1))
dev.off()