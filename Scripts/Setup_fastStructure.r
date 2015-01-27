# extract sample info from structure output (via stacks)
# make a .fam file for input into fastStructure
# .fam file should have sample name as found in structure file,
# family info (dummy info since we don't know that), and sex info

# This script find which samples were included in the final library and
# pair them with their respective sex and dummy-family info.

library(stringr)

# what the hell does the VCF file look like?
vcf.data<-read.table("~/Dropbox/ddRADseq/pseudoref_raccoon_counts_maf0.1_minmeanDP20_minGQ25.recode.vcf")
dim(vcf.data)
vcf.data[1,]

# read in the structure data to get the sample IDs
stru.data<-read.table("~/Dropbox//ddRADseq/pseudoref_raccoon_counts_maf0.1_minmeanDP20_minGQ25_structure.stru", skip=1)
in.library<-stru.data[,1]

# read in a separate file that contains sample IDs and sex info for the entire dataset
sample.info<-read.csv("~/Dropbox/ddRADseq/CLNWR_tickID_sex_data.csv")

# separate barcodes and sample IDs
barcodes<-str_match(in.library, pattern='_.*$')
IDs<-str_match(in.library, pattern='[O|R][0123456789]+T[0123456789]+')
bc.dict<-data.frame(cbind(IDs, barcodes))
names(bc.dict)<-c('tick.id', 'barcode')
head(bc.dict)

# merge sample info and barcode data
merged.IDs<-merge(bc.dict, sample.info)

# properly code the sex of the tick (1=M, 2=F)
sex.code<-sub('M', 1, merged.IDs$sex, ignore.case=TRUE)
sex.code<-sub('F', 2, sex.code, ignore.case=TRUE)
merged.IDs<-cbind(merged.IDs, sex.code)
merged.unique<-unique(merged.IDs)

# make full fam file with numeric IDs for all samples
# first 4 columns are numeric ID
# fifth column is sex code
# sixth column is -9 (missing phenotype data)

fam.file<-cbind(seq(0,length(merged.unique$sex.code)-1),
                seq(0,length(merged.unique$sex.code)-1),
                seq(0,length(merged.unique$sex.code)-1),
                seq(0,length(merged.unique$sex.code)-1),
                merged.unique$sex.code,
                rep(-9, times=length(merged.unique$sex.code)))

#### CLUMSY WAY TO MERGE ####

# find the IDs that preceed underscores. IDs are four or five characters long
# after finding the ID sequences, remove the underscores and barcode info coded in the structure data
four<-str_match(in.library, pattern='^...._')
four<-str_match(four, pattern='^....')
five<-str_match(in.library, pattern='^....._')
five<-str_match(five, pattern='^.....')
all.IDs<-c(four, five)
all.IDs<-data.frame(all.IDs[complete.cases(all.IDs)])
names(all.IDs)<-'tick.id'

# merge the two tables and remove duplicates
in.library.sex<-merge(all.IDs, sample.info)
head(in.library.sex)
no.duplicates<-unique(in.library.sex)