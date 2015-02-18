#Number of reads required for different estimates of phi

library(magicaxis)

#######################################################
# Bolnick formula for calculating required read depth #
#######################################################

d=30 #coverage depth
fl=550 #average fragment length
x=20000 #desired number of SNPs
n=270 #number of samples

phi=seq(0,0.01,by=0.001)

nr=n*d*(x/phi)/fl
nl=nr/200000000

png(file='~/Dropbox/ColoState//Projects/Outputs/LanesNeeded.png', height=10, width=10,
    res=300, unit='cm')
par(mar=c(5,6,6,2))
plot(phi, nl, xlab="Per Site SNP Probability",  
     ylab="Lanes of Illumina HiSeq\n(est 200 million reads/lane)",
     main="20,000 SNPs\n30X Coverage Depth\n250 Samples", las=1)
dev.off()

#######################################################################
# How many SNPs do we expect to find for a given fragment size range? #
#######################################################################

gs=3000000000 #genome size
dr=0.01 #percent in size range, expressed as a decimal
it=200 #illumina type; 2x100=200
br=gs*dr/fl*it #no. bases to read per sample
ps=br*seq(0, 0.01, 0.001) #possible number of SNPs
png(file='~/Dropbox/ColoState/Projects/Outputs/Expected_SNPs.png', height=10, width=10, res=300, unit='cm')
plot(seq(0, 0.01, 0.001), ps, axes=FALSE, ylab="", xlab="Per Base SNP Probability", main="SNPs from 1% of Deer Genome")
magaxis(side=1:2, las=1, ylab="Expected Number of SNPs", mtline=3.5)
dev.off()
