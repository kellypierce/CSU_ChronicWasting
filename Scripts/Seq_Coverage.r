#Number of reads required for different estimates of phi

d=30 #coverage depth
fl=300 #average fragment length
x=20000 #desired number of SNPs
n=250 #number of samples

phi=seq(0,0.01,by=0.001)

nr=n*d*(x/phi)/fl
nl=nr/200000000
par(mar=c(5,5,5,2))
png(file='~/Dropbox/ColoState//Projects/Outputs/LanesNeeded.png', height=10, width=10,
    res=300, unit='cm')
plot(phi, nl, xlab="Per Site SNP Probability",  
     ylab="Lanes of Illumina HiSeq\n(est 200 million reads/lane)",
     main="30X Coverage Depth/20,000 SNPs/250 Samples", las=1)
dev.off()
