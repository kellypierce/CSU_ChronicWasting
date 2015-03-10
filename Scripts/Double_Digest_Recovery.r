# DNA recovery

post.dig<-read.csv('~/Dropbox/ColoState/Projects/Outputs/PilotLibraryConstruction.csv')

hist(post.dig$X..Recovery)
boxplot(post.dig$X..Recovery)

mass.lost<-post.dig$Mass.Digested-post.dig$Mass.Recovered..ng..assuming.35ul.elution.post.quant.
boxplot(mass.lost)

high.yield<-post.dig[which(post.dig$Mass.Recovered..ng..assuming.35ul.elution.post.quant.>350),]
png(file='~/Dropbox/ColoState/Projects/Outputs/Post_Digest_DNA_Recovery.png', height=10, width=15, res=300, units='cm')
par(mar=c(5,5,5,2))
plot(post.dig$Mass.Digested, post.dig$Mass.Recovered..ng..assuming.35ul.elution.post.quant., pch=1, las=1, xlab='Mass Digested (ng)',
     ylab='Mass Recovered (ng)', main="Pre- and Post-Digest DNA Mass", xlim=c(997,1005), ylim=c(200,700), cex.lab=1.5, cex.main=1.5, axes=FALSE)
points(high.yield$Mass.Digested, high.yield$Mass.Recovered..ng..assuming.35ul.elution.post.quant., pch=16, col='red')
axis(side=1, at=seq(997,1005,1))
axis(side=2, at=seq(200,700,100), las=1)
dev.off()