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
