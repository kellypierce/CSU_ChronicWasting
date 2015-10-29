# ------------------------------------------------------------------
# -- CIs for Simulated Data
# ------------------------------------------------------------------

#rm(list=ls())

library(gplots)
library(ggplot2)
library(RColorBrewer)
#library(rPython)

# -- KL CI -- better with extreme proportions
KL.CI<-function(proportion, num.pos, sample.size, data){
  z.sq<-1.96**2
  delta<-((z.sq/3)+(1/6))*((1-(2*proportion))/sample.size)
  nu.proportion<-(proportion*(1-proportion))/(sample.size-1)
  delta.sq<-delta**2
  ci.u<-proportion+delta+sqrt((z.sq*nu.proportion)+delta.sq)
  ci.l<-proportion+delta-sqrt((z.sq*nu.proportion)+delta.sq)
  updated.data<-cbind(data,ci.u,ci.l)
  return(updated.data)}

# -- plot the CIs

plot.prev.CI<-function(data.set, xlab, ylab, name, main){
  
  ## -- read in data, add prevalence, and calculate CI
  data.set[4]<-data.set[2]/data.set[3]
  data.CI<-KL.CI(data.set[4], data.set[2], data.set[3], data.set)
  names(data.CI)<-c('animal', 'k', 'n', 'p','ci.u', 'ci.l')
  
  ## -- get colors; add extra because the lightest colors are hard to see; reverse the series so the darkest are used first
  p<-rev(brewer.pal((length(unique(data.set[,1]))+2), 'Blues'))
  
  ## -- build and save the pdf
  pdf(file=name, width=15, height=10)
  plotCI(x=c(1:length(data.CI[,1])), y=data.CI[,4], 
         sfrac=0.005, asp=FALSE, lwd=2, pch=16, cex.main=2, 
         cex.lab=2, cex.axis=2, ui=data.CI$ci.u, li=data.CI$ci.l, 
         col=p[data.set[,1]], barcol=p[data.set[,1]], 
         gap=0, xlab='', ylab='', xaxt='n', main=paste("Confidence Intervals", main))
  legend(x="topleft", legend=unique(data.set[,1]), col=p, cex=1.3, pch=16, bty='n')
  dev.off()
}

# -- notes on color indexing:
# how does this color business work? well, data.set[,1] contains factors (animal name), 
# and when used as indexes evaluate to integers... each unique animal name 
# evaluates to a unique number
# so... we can get the anmial names to actually index a vector 
# of colors without explicitly converting them to integers

# -- read in data
#data<-read.csv('/Users/kellypierce/Dropbox/ModelFitting/FutureProof/SimulationScripts/OnePref_TwoTrans_Simulated_lite/OnePref_TwoTrans_SimDataLITE_forMCMC.csv', header=FALSE)

# -- make plot

#plot.prev.CI(data, 'Simulated Dataset', 'Prevalence', 'Sim_7_1_5.pdf', 'phi=0.7, rhoTD=0.1, rhoDT=0.5')


