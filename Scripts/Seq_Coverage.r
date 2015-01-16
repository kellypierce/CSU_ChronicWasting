#Number of reads required for different estimates of phi

d=30 #coverage depth
fl=300 #average fragment length
x=20000 #desired number of SNPs
n=250 #number of samples

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
