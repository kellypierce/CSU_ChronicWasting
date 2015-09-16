# DNA-B sample quantification summary statistics

dna.quant<-read.csv('~/Dropbox/ColoState/Projects/Outputs//DNA_Quantification.csv')
head(dna.quant)

par(mfrow=c(1,1))
conc<-hist(dna.quant$Conc...ng.ul., breaks=seq(0.0, 10, 0.5), main="DNA-B Concentration", xlab="ng/ul", labels=TRUE, ylim=c(0,50))
abline(v=(2+(2/3)), col='red', lwd=5, lty=3)
quant<-hist(dna.quant$Mass, main="DNA-B Estimated Mass in 375ul", breaks=seq(0,3500,250), xlab="ng", labels=TRUE, ylim=c(0,25))
abline(v=1000, col='red', lwd=5, lty=3)

png(file='~/Dropbox/ColoState/Projects/Outputs/DNA_for_Digest.png', height=10, width=20, unit='cm', res=300)
hist(dna.quant$Approx..ng.DNA.avail..For.final.digest, breaks=seq(0, 3000, 100), 
     xlab="Approx. ng DNA avail for digest", main="", labels=TRUE, ylim=c(0,35))
dev.off()

min(dna.quant[complete.cases(dna.quant),]$Approx..ng.DNA.avail..For.final.digest)

exhaust.sample<-dna.quant[which(dna.quant$Approx..ng.DNA.avail..For.final.digest<500),]
head(exhaust.sample)
exhaust.DNA.A<-exhaust.sample[which(exhaust.sample$DNA.A.or.DNA.B.== 'A'),]
exhaust.DNA.A

# Buffy coat extraction summary statistics

# lazy entry of recovery stats (no excel sheet yet)
twenty<-c(850, 90.5, 1279.2, 1606.8, 1021.8, 35.1, 579.15, 573.5, 873.6, 657.15, 430.95)
forty<-c(1131, 795.6, 1154.4, 1053, 2277.6, 943.8, 1567.8, 626.3, 842.4, 343.2, 1583.4)

png(file='~/Dropbox/ColoState/Projects/Outputs/BuffyCoatExt.png', height=10, width=12, units='cm', res=300)
par(mfrow=c(1,1), mar=c(4,5,2,2))
boxplot(recovery, las=1, ylab='', main='Buffy Coat Extractions', 
        names=c(expression(paste('20', mu, 'l')), expression(paste('40', mu, 'l'))))
mtext(side=2, text='Mass DNA (ng)', line=4)
dev.off()
