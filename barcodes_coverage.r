bc<-read.csv('~/Dropbox/ColoState/Projects/Outputs/BarcodesCoverage.csv')

num.bc<-as.numeric(as.character(bc$Sequence.Count))
cc.bc<-num.bc[complete.cases(num.bc)]
nz.bc<-cc.bc[cc.bc>0]
png('~/Dropbox/ColoState/Projects/Outputs/round1_lib_read_counts.png', height=10, width=10, units='cm', res=300)
boxplot(nz.bc, axes=F, ylab='Number Reads per Individual (Millions)')
axis(side=2, las=1, labels=seq(0, 300, 100), at=seq(0, 3000000, 1000000))
dev.off()

length(nz.bc)
