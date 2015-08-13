# ddRADseq library QC for first batch of samples processed

data<-read.csv('~/Dropbox/ColoState/Projects/Outputs/DNA_Prep_QC.csv')
#why are there so many extra columns??
data<-data[,1:16] #that's better

#duplicated samples each get their own row... remove them
unique.samples<-unique(data$Sample.ID)
u.data<-c()
for(i in 1:length(unique.samples)){
  #print(data[data$Sample.ID==unique.samples[i],][1,])
  u.data<-rbind(u.data, data[data$Sample.ID==unique.samples[i],][1,])
}

# plot successes and failures
par(mfrow=c(1,1), mar=c(3,3,2,2))
barplot(table(u.data$Pass.Fail.Final.QC), ylim=c(0,200))

# count successes and failures
failed<-subset(u.data, Pass.Fail.Final.QC=='fail')
passed<-subset(u.data, Pass.Fail.Final.QC=='pass')

# count and plot reasons for failures
outcomes<-table(u.data$Fail.reason)
png('~/Dropbox/ColoState/Projects/Outputs/first_pass_library_success.png', height=10, width=15, 
    units='cm', res=300)
bar<-barplot(outcomes, ylim=c(0,225))
text(x=bar, y=outcomes+20, label=outcomes)
dev.off()

# some failed samples are in libraries with successful samples
# visualize the composition of libraries in terms of failed and successful samples
out.lib<-table(u.data$Fail.reason, u.data$Library.Number)
split.lib<-out.lib[,4:6]
split.lib<-rbind(split.lib[1,], split.lib[4,])
png('~/Dropbox/ColoState/Projects/Outputs/first_pass_library_success_T4_vs_low_conc.png', height=12, width=18, 
    units='cm', res=300)
barplot(split.lib, beside=T, ylim=c(0,10), legend.text=c('low', 'T4'))
dev.off()
