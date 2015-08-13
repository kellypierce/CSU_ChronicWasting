## White-tailed deer WGS assembly info

# tab-delim. file with info on contig length
wtd.wgs<-read.table('~/Downloads/AEGZ01.txt')
names(wtd.wgs)<-c('#', 'Accession', 'Name', 'Has annotation',
                  'Length', '# proteins')

# how big are the contigs?
max(wtd.wgs$Length)
min(wtd.wgs$Length)
hist(wtd.wgs$Length) # too wide a range/too highly skewed to be informative
hist(wtd.wgs$Length, breaks=100, axes=FALSE)
axis(side=1, at=seq(0,20000,1000))

big.contig<-subset(wtd.wgs, Length>1000)
small.contig<-subset(wtd.wgs, Length<=1000)

# 
hist(big.contig$Length, breaks=40, xlim=c(500, 20000))
hist(small.contig$Length)
