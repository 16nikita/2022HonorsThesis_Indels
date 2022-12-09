{
#READ DSBs & PREP FOR CROSSMAP:
DSB_v3<-read.delim("3126_peaks_information_2k_midpoint_v3_extend_2K.bed",sep='\t', header = FALSE)
colnames(DSB_v3)<-c('chrom','start','end')
DSB_v3$midpoint<-(DSB_v3$start+DSB_v3$end)/2
DSB_v3$length<-DSB_v3$end-DSB_v3$midpoint
DSB_v3$start<-DSB_v3$midpoint
DSB_v3$end<-DSB_v3$midpoint +1
DSB_v3$start<-as.integer(round(DSB_v3$start))
DSB_v3$end<-as.integer(round(DSB_v3$end))
write.table(DSB_v3,"DSB_v3.bed", row.names = F,quote=FALSE, sep = '\t')

#READ CROSSMAP OUTPUT:
DSB<-read.delim("DSB_v4.bed",sep='\t', header = FALSE)
colnames(DSB)<-c('chrom','start','end','midpoint','length')
DSB$end<-DSB$start+DSB$length
DSB$start<-DSB$start-DSB$length
DSB$start<-as.integer(round(DSB$start))
DSB$end<-as.integer(round(DSB$end))
DSB$length<-DSB$end-DSB$start
DSB<-DSB[!grepl("B73", DSB$chrom),]
write.table(DSB,"DSB_V4.bed", row.names = F,quote=FALSE, sep = '\t')
}
{
##READ INDELS (up to 50kb)
large_indels<-read.delim("edit_hufford_indels_NAM.txt",sep='\t', header = FALSE)
colnames(large_indels)<-c("type","chrom","start","end","length")
large_indels$chrom<-gsub("^chr","",large_indels$chrom)
indelsLarge_all<-large_indels[which(large_indels$length<=50000),]
indelsLarge_all<-indelsLarge_all[which(indelsLarge_all$type=="del" | indelsLarge_all$type=="ins"),]
indelsLarge_all<-indelsLarge_all[-c(1)]
indelsLarge_all$midpoint<-(indelsLarge_all$start+indelsLarge_all$end)/2
indelsLarge_all$length<-indelsLarge_all$end-indelsLarge_all$midpoint
indelsLarge_all$start<-indelsLarge_all$midpoint
indelsLarge_all$end<-indelsLarge_all$midpoint +1
indelsLarge_all$start<-as.integer(round(indelsLarge_all$start))
indelsLarge_all$end<-as.integer(round(indelsLarge_all$end))
write.table(indelsLarge_all,"indelsLarge_all_v5.bed", row.names = F,quote=FALSE, sep = '\t')

#AFTER CROSSMAP
indelsLarge_all<-read.delim("indelsLarge_all_v4.bed",sep='\t', header = FALSE)
colnames(indelsLarge_all)<-c("chrom","start","end","length","midpoint")
indelsLarge_all$start<-as.integer(round(indelsLarge_all$start))
indelsLarge_all$end<-as.integer(round(indelsLarge_all$end))
indelsLarge_all$end<-indelsLarge_all$start+indelsLarge_all$length
indelsLarge_all$start<-indelsLarge_all$start-indelsLarge_all$length
indelsLarge_all$length<-indelsLarge_all$end-indelsLarge_all$start
indelsLarge_all<-indelsLarge_all[!grepl("B73", indelsLarge_all$chrom),]
indelsLarge_all$start<-as.integer(round(indelsLarge_all$start))
indelsLarge_all$end<-as.integer(round(indelsLarge_all$end))

indelsLarge_10kb<-indelsLarge_all[!(indelsLarge_all$length>10000),]
indelsLarge_20kb<-subset(indelsLarge_all, length <= 20000 & length > 10000)
indelsLarge_30kb<-subset(indelsLarge_all, length <= 30000 & length > 20000)
indelsLarge_40kb<-subset(indelsLarge_all, length <= 40000 & length > 30000)
indelsLarge_50kb<-subset(indelsLarge_all, length <= 50000 & length > 40000)
write.table(indelsLarge_all[-c(4:5)],"indelsLarge.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indelsLarge_10kb[-c(4:5)],"indelsLarge_10kb.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indelsLarge_20kb[-c(4:5)],"indelsLarge_20kb.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indelsLarge_30kb[-c(4:5)],"indelsLarge_30kb.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indelsLarge_40kb[-c(4:5)],"indelsLarge_40kb.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indelsLarge_50kb[-c(4:5)],"indelsLarge_50kb.bed", row.names = F,quote=FALSE, sep = '\t')
}

#READ INDELS 100BP-500BP
{
  indelsLarge_all<-read.delim("indelsLarge_all_v4.bed",sep='\t', header = FALSE)
  colnames(indelsLarge_all)<-c("chrom","start","end","length","midpoint")
  indelsLarge_all$start<-as.integer(round(indelsLarge_all$start))
  indelsLarge_all$end<-as.integer(round(indelsLarge_all$end))
  indelsLarge_all$end<-indelsLarge_all$start+indelsLarge_all$length
  indelsLarge_all$start<-indelsLarge_all$start-indelsLarge_all$length
  indelsLarge_all$length<-indelsLarge_all$end-indelsLarge_all$start
  indelsLarge_all<-indelsLarge_all[!grepl("B73", indelsLarge_all$chrom),]
  indelsLarge_all$start<-as.integer(round(indelsLarge_all$start))
  indelsLarge_all$end<-as.integer(round(indelsLarge_all$end))
  
  indelsLarge_all<-indelsLarge_all[!(indelsLarge_all$length>500),]
  indelsLarge_all<-indelsLarge_all %>% distinct()
  indelsLarge_100bp<-indelsLarge_all[!(indelsLarge_all$length>150),]
  indelsLarge_200bp<-subset(indelsLarge_all, length <= 200 & length > 150)
  indelsLarge_300bp<-subset(indelsLarge_all, length <= 300 & length > 200)
  indelsLarge_400bp<-subset(indelsLarge_all, length <= 400 & length > 300)
  indelsLarge_500bp<-subset(indelsLarge_all, length <= 500 & length > 400)
  write.table(indelsLarge_all[-c(4:5)],"indelsLarge_to500bp.bed", row.names = F,quote=FALSE, sep = '\t')
  write.table(indelsLarge_100bp[-c(4:5)],"indelsLarge_100bp.bed", row.names = F,quote=FALSE, sep = '\t')
  write.table(indelsLarge_200bp[-c(4:5)],"indelsLarge_200bp.bed", row.names = F,quote=FALSE, sep = '\t')
  write.table(indelsLarge_300bp[-c(4:5)],"indelsLarge_300bp.bed", row.names = F,quote=FALSE, sep = '\t')
  write.table(indelsLarge_400bp[-c(4:5)],"indelsLarge_400bp.bed", row.names = F,quote=FALSE, sep = '\t')
  write.table(indelsLarge_500bp[-c(4:5)],"indelsLarge_500bp.bed", row.names = F,quote=FALSE, sep = '\t')

}
#ASSIGN INSERTION LENGTHS:
library(stringr)
insertions_60bp<-read.delim("insertions_NAM_60bp_2.bed",sep='\t', header = FALSE)
insertions_60bp<-insertions_60bp[c(1:8)]
insertions_60bp<-insertions_60bp[which(insertions_60bp$V8=="PASS"),]
insertions_60bp$length<-str_length(insertions_60bp$V7)
insertions_60bp<-insertions_60bp[-c(4:8)]
colnames(insertions_60bp)<-c('chrom','start','end','length')
insertions_all<-insertions_60bp[!(insertions_60bp$length>50),]

#PREPPING INSERTION DATASET FOR CROSSMAP:
insertions_all$midpoint<-(insertions_all$start+insertions_all$end)/2
insertions_all$start<-insertions_all$midpoint
insertions_all$end<-insertions_all$midpoint +1
insertions_all$start<-as.integer(round(insertions_all$start))
insertions_all$end<-as.integer(round(insertions_all$end))
write.table(insertions_all,"insertions_50bp_v5.bed", row.names = F,quote=FALSE, sep = '\t')

#READING CROSSMAP OUTPUT (INSERTIONS):
insertions_all<-read.delim("insertions_50bp_v4.bed",sep='\t', header = FALSE)
colnames(insertions_all)<-c("chrom","start","end","length","midpoint")
insertions_all$start<-as.integer(round(insertions_all$start))
insertions_all$end<-insertions_all$start+1
insertions_all$chrom<-gsub("^chr","",insertions_all$chrom)
insertions_10bp<-insertions_all[!(insertions_all$length>10),]
insertions_20bp<-subset(insertions_all, length <= 20 & length > 10)
insertions_30bp<-subset(insertions_all, length <= 30 & length > 20)
insertions_40bp<-subset(insertions_all, length <= 40 & length > 30)
insertions_50bp<-subset(insertions_all, length <= 50 & length > 40)

write.table(insertions_all[-c(4:5)],"insertions_v4.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(insertions_10bp[-c(4:5)],"insertions_10bp.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(insertions_20bp[-c(4:5)],"insertions_20bp.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(insertions_30bp[-c(4:5)],"insertions_30bp.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(insertions_40bp[-c(4:5)],"insertions_40bp.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(insertions_50bp[-c(4:5)],"insertions_50bp.bed", row.names = F,quote=FALSE, sep = '\t')


#MAKE DELETION INDEL FILE:
library(stringr)
indels_NAM<-read.delim("deletions_NAM_60bp.bed",sep='\t', header = FALSE)
indels_NAM<-indels_NAM[c(1:8)]
indels_NAM<-indels_NAM[which(indels_NAM$V8=="PASS"),]
indels_NAM<-indels_NAM[-c(4:8)]
colnames(indels_NAM)<-c('chrom','start','end')
indels_NAM$length<-indels_NAM$end-indels_NAM$start
indels_NAM<-indels_NAM[!(indels_NAM$length>50),]

#PREPPING INDEL DATASET FOR CROSSMAP:
indels_NAM$midpoint<-(indels_NAM$end+indels_NAM$start)/2
indels_NAM$length<-indels_NAM$end-indels_NAM$midpoint
indels_NAM$start<-indels_NAM$midpoint
indels_NAM$end<-indels_NAM$midpoint +1
indels_NAM$start<-as.integer(round(indels_NAM$start))
indels_NAM$end<-as.integer(round(indels_NAM$end))
indels_NAM$chrom<-gsub("^chr","",indels_NAM$chrom)
write.table(indels_NAM,"deletions_50bp_v5.bed", row.names = F,quote=FALSE, sep = '\t')

#READING CROSSMAP OUTPUT:
indels_NAMv4<-read.delim("deletions_50bp_v4.bed",sep='\t', header = FALSE)
colnames(indels_NAMv4)<-c("chrom","start","end","length","midpoint")
indels_NAMv4<-indels_NAMv4[!grepl("B73", indels_NAMv4$chrom),]
indels_NAMv4$start<-as.integer(round(indels_NAMv4$start))
indels_NAMv4$end<-as.integer(round(indels_NAMv4$end))
indels_NAMv4$end<-indels_NAMv4$start+indels_NAMv4$length
indels_NAMv4$start<-indels_NAMv4$start-indels_NAMv4$length
indels_NAMv4$start<-as.integer(round(indels_NAMv4$start))
indels_NAMv4$end<-as.integer(round(indels_NAMv4$end))
indels_NAMv4$length<-indels_NAMv4$end-indels_NAMv4$start
indels_10bp<-indels_NAMv4[!(indels_NAMv4$length>10),]
indels_20bp<-subset(indels_NAMv4, length <= 20 & length > 10)
indels_30bp<-subset(indels_NAMv4, length <= 30 & length > 20)
indels_40bp<-subset(indels_NAMv4, length <= 40 & length > 30)
indels_50bp<-subset(indels_NAMv4, length <= 50 & length > 40)
write.table(indels_NAMv4,"indels_NAMv4.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indels_10bp,"indels_10bp.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indels_20bp,"indels_20bp.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indels_30bp,"indels_30bp.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indels_40bp,"indels_40bp.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indels_50bp,"indels_50bp.bed", row.names = F,quote=FALSE, sep = '\t')

# insertions_NAM_50bp<-indels_NAM[(indels_NAM$type=="insertions"),]
# insertions_NAM_50bp<-insertions_NAM_50bp[-c(4:6)]
# deletions_NAM_50bp<-indels_NAM[(indels_NAM$type=="deletions"),]
# deletions_NAM_50bp<-deletions_NAM_50bp[-c(4:6)]
# write.table(insertions_NAM_50bp,"insertions_NAM_50bp.bed", row.names = F,quote=FALSE, sep = '\t')
# write.table(deletions_NAM_50bp,"deletions_NAM_50bp.bed", row.names = F,quote=FALSE, sep = '\t')

#DISTRIBUTION/DENSITY PLOTS OF EACH FILE

indels_NAM<-read.delim("indelsLarge_to500bp.bed",sep='\t', header = TRUE)[c(1:3)]
indels_NAM<-indels_10bp
indels_NAM<-indels_20bp
indels_NAM<-indels_30bp
indels_NAM<-indels_40bp
indels_NAM<-indels_50bp

indels_NAM<-read.delim("indels_all.bed",sep='\t', header = TRUE)[c(1:3)]
indels_NAM<-read.delim("indelsLarge_to500bp.bed",sep='\t', header = TRUE)[c(1:3)]
indels_NAM<-read.delim("indelsLarge.bed",sep='\t', header = TRUE)[c(1:3)]


colnames(indels_NAM)<-c('chrom','start','end')
indels_NAM$midpoint<-(indels_NAM$start+indels_NAM$end)/2
ref_indels1 <- indels_NAM[which(indels_NAM$chrom == 1),]
ref_indels2 <- indels_NAM[which(indels_NAM$chrom == 2),]
ref_indels3 <- indels_NAM[which(indels_NAM$chrom == 3),]
ref_indels4 <- indels_NAM[which(indels_NAM$chrom == 4),]
ref_indels5 <- indels_NAM[which(indels_NAM$chrom == 5),]
ref_indels6 <- indels_NAM[which(indels_NAM$chrom == 6),]
ref_indels7 <- indels_NAM[which(indels_NAM$chrom == 7),]
ref_indels8 <- indels_NAM[which(indels_NAM$chrom == 8),]
ref_indels9 <- indels_NAM[which(indels_NAM$chrom == 9),]
ref_indels10 <- indels_NAM[which(indels_NAM$chrom == 10),]

library(Rcpp)
library(ggplot2)
library(dplyr) 
library(dlookr)
library(tidyverse)
library(OneR)


#Bin indels on each chromosome into equal (~1 Mb) intervals, create dataframe with bin info (position, freq, rate, gene density)  
#Variables foo.X1 and foo.X2 represent interval start and end positions, respectively
fix_table<-function(bin){
  bin <- within(bin, foo<-data.frame(do.call('rbind', strsplit(as.character(bin$levels), ',', fixed=TRUE))))
  bin <- do.call(data.frame, bin)
  bin <- bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
  bin <- bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
  bin[1,4]<-substring(bin$levels[1],unlist(gregexpr("\\[", bin$levels[1]))[1]+1,unlist(gregexpr(',', bin$levels[1]))[1]-1)
  bin$foo.X1 <- as.numeric(bin$foo.X1)
  bin$length <- as.numeric(bin$foo.X2)-as.numeric(bin$foo.X1)
  bin$density <- (bin$freq/(bin$length/1000000))
  return(as.data.frame(bin))
}
remove_outliers<-function(bin){
  quartiles <- quantile(bin$density, probs=c(.25, .75), na.rm = FALSE)
  bin_no_outlier <- subset(bin,bin$density < quartiles[2] + 1.5*IQR(bin$density))
  return(as.data.frame(bin_no_outlier))
}
{
indels_bin1 <- as.data.frame(summary(binning(ref_indels1$midpoint, nbins = round(max(ref_indels1$start)/1000000), type = "kmeans")))
indels_bin1 <- fix_table(indels_bin1)
indels_bin1 <- remove_outliers(indels_bin1)
indel_1_plot<- ggplot(indels_bin1, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 1", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin1, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)
# #indel_1_plot<- ggplot(indels_bin1, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
#   theme_minimal() + labs(title="Indel Density Along Chromosome 1", x="Position (Mb)", y = "Indel Density (Indels/Mb)")

indel_1_plot

indels_bin2 <- as.data.frame(summary(binning(ref_indels2$midpoint, nbins = round(max(ref_indels2$start)/1000000), type = "kmeans")))
indels_bin2 <- fix_table(indels_bin2)
indels_bin2 <- remove_outliers(indels_bin2)
indel_2_plot<- ggplot(indels_bin2, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 2", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin2, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)

indels_bin3 <- as.data.frame(summary(binning(ref_indels3$midpoint, nbins = round(max(ref_indels3$start)/1000000), type = "kmeans")))
indels_bin3 <- fix_table(indels_bin3)
indels_bin3 <- remove_outliers(indels_bin3)
indel_3_plot<-ggplot(indels_bin3, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 3", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin3, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)

indels_bin4 <- as.data.frame(summary(binning(ref_indels4$midpoint, nbins = round(max(ref_indels4$start)/1000000), type = "kmeans")))
indels_bin4 <- fix_table(indels_bin4)
indels_bin4 <- remove_outliers(indels_bin4)
indel_4_plot<-ggplot(indels_bin4, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 4", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin4, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)

indels_bin5 <- as.data.frame(summary(binning(ref_indels5$midpoint, nbins = round(max(ref_indels5$start)/1000000), type = "kmeans")))
indels_bin5 <- fix_table(indels_bin5)
indels_bin5 <- remove_outliers(indels_bin5)
indel_5_plot<-ggplot(indels_bin5, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 5", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin5, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)

indels_bin6 <- as.data.frame(summary(binning(ref_indels6$midpoint, nbins = round(max(ref_indels6$start)/1000000), type = "kmeans")))
indels_bin6 <- fix_table(indels_bin6)
indels_bin6 <- remove_outliers(indels_bin6)
indel_6_plot<-ggplot(indels_bin6, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 6", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin6, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)

indels_bin7 <- as.data.frame(summary(binning(ref_indels7$midpoint, nbins = round(max(ref_indels7$start)/1000000), type = "kmeans")))
indels_bin7 <- fix_table(indels_bin7)
indels_bin7 <- remove_outliers(indels_bin7)
indel_7_plot<-ggplot(indels_bin7, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 7", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin7, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)

indels_bin8 <- as.data.frame(summary(binning(ref_indels8$midpoint, nbins = round(max(ref_indels8$start)/1000000), type = "kmeans")))
indels_bin8 <- fix_table(indels_bin8)
indels_bin8 <- remove_outliers(indels_bin8)
indel_8_plot<-ggplot(indels_bin8, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 8", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin8, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)

indels_bin9 <- as.data.frame(summary(binning(ref_indels9$midpoint, nbins = round(max(ref_indels9$start)/1000000), type = "kmeans")))
indels_bin9 <- fix_table(indels_bin9)
indels_bin9 <- remove_outliers(indels_bin9)
indel_9_plot<-ggplot(indels_bin9, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 9", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin9, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)

indels_bin10 <- as.data.frame(summary(binning(ref_indels10$midpoint, nbins = round(max(ref_indels10$start)/1000000), type = "kmeans")))
indels_bin10 <- fix_table(indels_bin10)
indels_bin10 <- remove_outliers(indels_bin10)
indel_10_plot<-ggplot(indels_bin10, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="Indel Density Along Chromosome 10", 
                         x="Position (Mb)", y = "Indel Density (Indels/Mb)")+ geom_smooth(data = indels_bin10, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)
indels_bin1<-cbind(chrom=1,indels_bin1)
indels_bin2<-cbind(chrom=2,indels_bin2)
indels_bin3<-cbind(chrom=3,indels_bin3)
indels_bin4<-cbind(chrom=4,indels_bin4)
indels_bin5<-cbind(chrom=5,indels_bin5)
indels_bin6<-cbind(chrom=6,indels_bin6)
indels_bin7<-cbind(chrom=7,indels_bin7)
indels_bin8<-cbind(chrom=8,indels_bin8)
indels_bin9<-cbind(chrom=9,indels_bin9)
indels_bin10<-cbind(chrom=10,indels_bin10)

annotation_file = rbind(indels_bin1,indels_bin2,indels_bin3,indels_bin4,indels_bin5,indels_bin6,indels_bin7,indels_bin8,indels_bin9,indels_bin10)[-c(2,3,4,7)]
colnames(annotation_file)=c('chrom','start','end','indel_per_Mb')
rownames(annotation_file) <- NULL
annotation_file=cbind(ID=rownames(annotation_file),annotation_file)

# maize_genome<-read.delim("maize.genome",header=FALSE)
# colnames(chrom_file)<-c('chrom','start','end','centromere')
# chrom_file$chrom<-c(1:10)
# cent<-c(121+133,92.2+101.6,87.4+89.4,71.8+93,94.8+118.4,31.8+32.7,44.9+54.1,55.8+80.7,34.1+34.2,37.1+46)
# cent<-cent/2
# chrom_file$centromere<-cent*1000000
# chrom_file$end<-maize_genome$V3
#save(chrom_file, file = "chrom_file.RData")

#install.packages("chromoMap")
library(chromoMap)
library('RColorBrewer')
display.brewer.all()
chromoMap(ch.files = list(chrom_file), data.files = list(annotation_file), title = "Large Indel Density Across Each Chromosome",data_based_color_map = T,
          data_type = "numeric", legend = T,lg_y= 500,n_win.factor = 1, export.options = T, plots = "bar", plot_color="purple",data_colors = list(brewer.pal(n = 3, name = "PuRd")))
#without chromosomes
# chromoMap(ch.files = list(chrom_file), data.files = list(annotation_file), title = "Large Indel Density Across Each Chromosome",data_based_color_map = T,
#           data_type = "numeric", legend = T,display.chr = c(F), lg_y= 500,n_win.factor = 1, export.options = T, plots = "bar", plot_color="purple",data_colors = list(brewer.pal(n = 3, name = "PuRd")))


# saveRDS(indel_1_plot,file = "50bpdeletion1.rds")
# saveRDS(indel_2_plot,file = "50bpdeletion2.rds")
# saveRDS(indel_3_plot,file = "50bpdeletion3.rds")
# saveRDS(indel_4_plot,file = "50bpdeletion4.rds")
# saveRDS(indel_5_plot,file = "50bpdeletion5.rds")
# saveRDS(indel_6_plot,file = "50bpdeletion6.rds")
# saveRDS(indel_7_plot,file = "50bpdeletion7.rds")
# saveRDS(indel_8_plot,file = "50bpdeletion8.rds")
# saveRDS(indel_9_plot,file = "50bpdeletion9.rds")
# saveRDS(indel_10_plot,file = "50bpdeletion10.rds")

}
indel_distribution <- ggarrange(indel_1_plot,indel_2_plot,indel_3_plot,indel_4_plot,indel_5_plot, indel_6_plot,indel_7_plot,indel_8_plot,indel_9_plot,indel_10_plot,
                                common.legend=TRUE,
                                labels = c("A", "B","C", "D","E","F","G","H","I","J"),
                                ncol = 5, nrow = 2)
annotate_figure(indel_distribution, top = text_grob("Indels (1-250 bp)", 
                                      color = "black", face = "bold", size = 14))
indel_distribution
#Read in NAM CO data
NAM <- read.table("NAM_US_COs_v4.txt", header = TRUE)
colnames(NAM) <- c("Chr", "CO Start", "CO End")
NAM<-NAM[!grepl("chrB73", NAM$Chr),]
NAM <- NAM[order(NAM$Chr,NAM$`CO Start`),]
NAM$Chr<-gsub("^chr","",NAM$Chr)

chr1_CO <- NAM[ which(NAM$Chr == 1),]
chr1_CO$midpoint <- (chr1_CO$`CO Start`+ chr1_CO$`CO End`)/2

chr2_CO <- NAM[ which(NAM$Chr == 2),]
chr2_CO$midpoint <- (chr2_CO$`CO Start`+ chr2_CO$`CO End`)/2

chr3_CO <- NAM[ which(NAM$Chr == 3),]
chr3_CO$midpoint <- (chr3_CO$`CO Start`+ chr3_CO$`CO End`)/2

chr4_CO <- NAM[ which(NAM$Chr == 4),]
chr4_CO$midpoint <- (chr4_CO$`CO Start`+ chr4_CO$`CO End`)/2

chr5_CO <- NAM[ which(NAM$Chr == 5),]
chr5_CO$midpoint <- (chr5_CO$`CO Start`+ chr5_CO$`CO End`)/2

chr6_CO <- NAM[ which(NAM$Chr == 6),]
chr6_CO$midpoint <- (chr6_CO$`CO Start`+ chr6_CO$`CO End`)/2

chr7_CO <- NAM[ which(NAM$Chr == 7),]
chr7_CO$midpoint <- (chr7_CO$`CO Start`+ chr7_CO$`CO End`)/2

chr8_CO <- NAM[ which(NAM$Chr == 8),]
chr8_CO$midpoint <- (chr8_CO$`CO Start`+ chr8_CO$`CO End`)/2

chr9_CO <- NAM[ which(NAM$Chr == 9),]
chr9_CO$midpoint <- (chr9_CO$`CO Start`+ chr9_CO$`CO End`)/2

chr10_CO <- NAM[ which(NAM$Chr == 10),]
chr10_CO$midpoint <- (chr10_CO$`CO Start`+ chr10_CO$`CO End`)/2

#bin crossovers into 1 Mb bins
fix_table<-function(bin){
  bin <- within(bin, foo<-data.frame(do.call('rbind', strsplit(as.character(bin$levels), ',', fixed=TRUE))))
  bin <- do.call(data.frame, bin)
  bin <- bin %>% dplyr::mutate(foo.X1 = as.numeric(gsub("\\(", "", foo.X1)))
  bin <- bin %>% dplyr::mutate(foo.X2 = as.numeric(gsub("]", "", foo.X2)))
  bin[1,4]<-substring(bin$levels[1],unlist(gregexpr("\\[", bin$levels[1]))[1]+1,unlist(gregexpr(',', bin$levels[1]))[1]-1)
  bin$foo.X1 <- as.numeric(bin$foo.X1)
  bin$length <- as.numeric(bin$foo.X2)-as.numeric(bin$foo.X1)
  bin$density <- (bin$freq/(bin$length/1000000))
  return(bin)
}

chr1_CO <- chr1_CO[order(chr1_CO$`CO Start`),]
chr1_bin <- as.data.frame(summary(binning(chr1_CO$midpoint, nbins = max(chr1_CO$`CO Start`)/200000, type = "kmeans")))
#transforming data; making bin interval into 2 columns
chr1_bin<-fix_table(chr1_bin)
chr1_CO_plot<- ggplot(chr1_no_outlier, aes(x=foo.X1/1000000, y=density)) + geom_bar(stat="identity", width=0.7, fill="steelblue")+
  theme_minimal() + labs(title="CO Density Along Chromosome 1", 
                         x="Position (Mb)", y = "CO Density (Indels/Mb)")+ geom_smooth(data = chr1_no_outlier, aes(x=foo.X1/1000000, y=density), colour="black",size=0.5,se = FALSE)


chr2_CO <- chr2_CO[order(chr2_CO$`CO Start`),]
chr2_bin <- as.data.frame(summary(binning(chr2_CO$midpoint, nbins = max(chr2_CO$`CO Start`)/200000, type = "kmeans")))
chr2_bin<-fix_table(chr2_bin)

chr3_CO <- chr3_CO[order(chr3_CO$`CO Start`),]
chr3_bin <- as.data.frame(summary(binning(chr3_CO$midpoint, nbins = max(chr3_CO$`CO Start`)/200000, type = "kmeans")))
chr3_bin<-fix_table(chr3_bin)

chr4_CO <- chr4_CO[order(chr4_CO$`CO Start`),]
chr4_bin <- as.data.frame(summary(binning(chr4_CO$midpoint, nbins = max(chr4_CO$`CO Start`)/200000, type = "kmeans")))
chr4_bin<-fix_table(chr4_bin)

chr5_CO <- chr5_CO[order(chr5_CO$`CO Start`),]
chr5_bin <- as.data.frame(summary(binning(chr5_CO$midpoint, nbins = max(chr5_CO$`CO Start`)/200000, type = "kmeans")))
chr5_bin<-fix_table(chr5_bin)

chr6_CO <- chr6_CO[order(chr6_CO$`CO Start`),]
chr6_bin <- as.data.frame(summary(binning(chr6_CO$midpoint, nbins = max(chr6_CO$`CO Start`)/200000, type = "kmeans")))
chr6_bin<-fix_table(chr6_bin)

chr7_CO <- chr7_CO[order(chr7_CO$`CO Start`),]
chr7_bin <- as.data.frame(summary(binning(chr7_CO$midpoint, nbins = max(chr7_CO$`CO Start`)/200000, type = "kmeans")))
chr7_bin<-fix_table(chr7_bin)

chr8_CO <- chr8_CO[order(chr8_CO$`CO Start`),]
chr8_bin <- as.data.frame(summary(binning(chr8_CO$midpoint, nbins = max(chr8_CO$`CO Start`)/200000, type = "kmeans")))
chr8_bin<-fix_table(chr8_bin)

chr9_CO <- chr9_CO[order(chr9_CO$`CO Start`),]
chr9_bin <- as.data.frame(summary(binning(chr9_CO$midpoint, nbins = max(chr9_CO$`CO Start`)/200000, type = "kmeans")))
chr9_bin<-fix_table(chr9_bin)

chr10_CO <- chr10_CO[order(chr10_CO$`CO Start`),]
chr10_bin <- as.data.frame(summary(binning(chr10_CO$midpoint, nbins = max(chr10_CO$`CO Start`)/200000, type = "kmeans")))
chr10_bin<-fix_table(chr10_bin)


#Spearmen correlation test between indel density and CO density
correlate_density<-function(indel_bin,chr_bin){
  indel_bin$'CO density' <-NA
  for(i in 1:nrow(indel_bin)){
    for(k in 1:nrow(chr_bin)){
      if(isTRUE((chr_bin$foo.X1[k]>=indel_bin$foo.X1[i] ) && (chr_bin$foo.X2[k] <= indel_bin$foo.X2[i]))){
        indel_bin$'CO density'[i] <- chr_bin$density[k]
      }
    }
  }
  return(indel_bin)
}
indels_bin1_test<-correlate_density(indels_bin1,chr1_bin)
indels_bin1_test<-na.omit(indels_bin1_test)
cor.test(indels_bin1_test$density, indels_bin1_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin2_test<-correlate_density(indels_bin2[-c(8:9)],chr2_bin)
indels_bin2_test<-na.omit(indels_bin2_test)
cor.test(indels_bin2_test$density, indels_bin2_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin3_test<-correlate_density(indels_bin3,chr3_bin)
indels_bin3_test<-na.omit(indels_bin3_test)
cor.test(indels_bin3_test$density, indels_bin3_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin4_test<-correlate_density(indels_bin4,chr4_bin)
indels_bin4_test<-na.omit(indels_bin4_test)
cor.test(indels_bin4_test$density, indels_bin4_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin5_test<-correlate_density(indels_bin5,chr5_bin)
indels_bin5_test<-na.omit(indels_bin5_test)
cor.test(indels_bin5_test$density, indels_bin5_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin6_test<-correlate_density(indels_bin6,chr6_bin)
indels_bin6_test<-na.omit(indels_bin6_test)
cor.test(indels_bin6_test$density, indels_bin6_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin7_test<-correlate_density(indels_bin7,chr7_bin)
indels_bin7_test<-na.omit(indels_bin7_test)
cor.test(indels_bin7_test$density, indels_bin7_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin8_test<-correlate_density(indels_bin8,chr8_bin)
indels_bin8_test<-na.omit(indels_bin8_test)
cor.test(indels_bin8_test$density, indels_bin8_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin9_test<-correlate_density(indels_bin9,chr9_bin)
indels_bin9_test<-na.omit(indels_bin9_test)
cor.test(indels_bin9_test$density, indels_bin9_test$'CO density',  method = "spearman", alternative = "greater")

indels_bin10_test<-correlate_density(indels_bin10,chr10_bin)
indels_bin10_test<-na.omit(indels_bin10_test)
cor.test(indels_bin10_test$density, indels_bin10_test$'CO density',  method = "spearman", alternative = "greater")

total_cor<-rbind(indels_bin1_test,indels_bin2_test,indels_bin3_test,indels_bin4_test,indels_bin5_test,indels_bin6_test,indels_bin7_test,indels_bin8_test,indels_bin9_test,indels_bin10_test)
cor.test(total_cor$density, total_cor$'CO density',  method = "spearman", alternative = "greater")
#   rho  0.4944223 

correlations <- as.data.frame(matrix(ncol=2, nrow=11))
colnames(correlations)<-c("chromosome","rho")
correlations$chromosome<-c("1","2","3","4","5","6","7","8","9","10","total")
correlations$rho<-c(0.4568343,0.6801134,0.3164979,0.4581461,0.559593,0.4500214,0.3526726,0.5666199,0.275576,0.6389917,0.4944223)
