library(dplyr)
#QUICK RELOAD:
DSB<-read.delim("DSB_V4.bed")
DSB_2kb_up<-read.delim("DSB_2kb_up.bed")
DSB_2kb_down<-read.delim("DSB_2kb_down.bed")
DSB_coldspots<-read.delim("DSB_coldspotsv4.bed")
DSB_cold_2kb_up<-read.delim("DSB_cold_2kb_up.bed")
DSB_cold_2kb_down<-read.delim("DSB_cold_2kb_down.bed")

##READ IN INTERSECT DATA:
# INTERSECT WITH PREDICTED CO REGIONS
{
  intersect_pred<-read.delim("intersect_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_pred$type<-"all deletions"
  intersect_pred$intersect<-"intersect"
  intersect_pred<-intersect_pred[-c(7:8)]
  colnames(intersect_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_pred<-read.delim("intersect_10bp_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_10bp_pred$type<-"1-10bp"
  intersect_10bp_pred$intersect<-"intersect"
  intersect_10bp_pred<-intersect_10bp_pred[-c(7:8)]
  colnames(intersect_10bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_pred<-read.delim("intersect_20bp_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_20bp_pred$type<-"11-20bp"
  intersect_20bp_pred$intersect<-"intersect"
  intersect_20bp_pred<-intersect_20bp_pred[-c(7:8)]
  colnames(intersect_20bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  #intersect_30bp_pred<-as.data.frame(matrix(data=0,ncol=7,nrow=1))
  intersect_30bp_pred<-read.delim("intersect_30bp_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_30bp_pred$type<-"21-30bp"
  intersect_30bp_pred$intersect<-"intersect"
  intersect_30bp_pred<-intersect_30bp_pred[-c(7:8)]
  colnames(intersect_30bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_pred<-read.delim("intersect_40bp_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_40bp_pred$type<-"31-40bp"
  intersect_40bp_pred$intersect<-"intersect"
  intersect_40bp_pred<-intersect_40bp_pred[-c(7:8)]
  colnames(intersect_40bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_pred<-read.delim("intersect_50bp_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_50bp_pred$type<-"41-50bp"
  intersect_50bp_pred$intersect<-"intersect"
  intersect_50bp_pred<-intersect_50bp_pred[-c(7:8)]
  colnames(intersect_50bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_pred<-read.delim("intersect_2kbup_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_2kbup_pred$type<-"all deletions"
  intersect_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_pred<-read.delim("intersect_10bp_2kbup_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_pred$type<-"1-10bp"
  intersect_10bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_pred<-read.delim("intersect_20bp_2kbup_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_pred$type<-"11-20bp"
  intersect_20bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_pred<-read.delim("intersect_30bp_2kbup_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_pred$type<-"21-30bp"
  intersect_30bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_pred<-read.delim("intersect_40bp_2kbup_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_pred$type<-"31-40bp"
  intersect_40bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_pred<-read.delim("intersect_50bp_2kbup_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_pred$type<-"41-50bp"
  intersect_50bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_pred<-read.delim("intersect_2kbdown_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_2kbdown_pred$type<-"all deletions"
  intersect_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_pred<-read.delim("intersect_10bp_2kbdown_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_pred<-read.delim("intersect_20bp_2kbdown_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_pred<-read.delim("intersect_30bp_2kbdown_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_pred<-read.delim("intersect_40bp_2kbdown_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_pred<-read.delim("intersect_50bp_2kbdown_aL_DSB.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}
# INTERSECT WITH RANDOM CO DESERT REGIONS
{
  intersect_cold_pred<-read.delim("intersect_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_cold_pred$type<-"all deletions"
  intersect_cold_pred$intersect<-"intersect"
  colnames(intersect_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_cold_pred<-read.delim("intersect_10bp_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_10bp_cold_pred$type<-"1-10bp"
  intersect_10bp_cold_pred$intersect<-"intersect"
  colnames(intersect_10bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_cold_pred<-read.delim("intersect_20bp_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_20bp_cold_pred$type<-"11-20bp"
  intersect_20bp_cold_pred$intersect<-"intersect"
  colnames(intersect_20bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_cold_pred<-read.delim("intersect_30bp_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_30bp_cold_pred$type<-"21-30bp"
  intersect_30bp_cold_pred$intersect<-"intersect"
  colnames(intersect_30bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  #no data?
  #intersect_40bp_cold_pred<-read.delim("intersect_40bp_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_40bp_cold_pred<-as.data.frame(matrix(data=0,ncol=7,nrow=1))
  intersect_40bp_cold_pred$type<-"31-40bp"
  intersect_40bp_cold_pred$intersect<-"intersect"
  colnames(intersect_40bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_cold_pred<-read.delim("intersect_50bp_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_50bp_cold_pred$type<-"41-50bp"
  intersect_50bp_cold_pred$intersect<-"intersect"
  colnames(intersect_50bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_cold_pred<-read.delim("intersect_2kbup_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_2kbup_cold_pred$type<-"all deletions"
  intersect_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_cold_pred<-read.delim("intersect_10bp_2kbup_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_cold_pred<-read.delim("intersect_20bp_2kbup_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_cold_pred<-read.delim("intersect_30bp_2kbup_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  #intersect_40bp_2kbup_cold_pred<-read.delim("intersect_40bp_2kbup_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_cold_pred<-as.data.frame(matrix(data=0,ncol=7,nrow=1))
  intersect_40bp_2kbup_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_cold_pred<-read.delim("intersect_50bp_2kbup_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_cold_pred<-read.delim("intersect_2kbdown_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_2kbdown_cold_pred$type<-"all deletions"
  intersect_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_cold_pred<-read.delim("intersect_10bp_2kbdown_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_cold_pred<-read.delim("intersect_20bp_2kbdown_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_cold_pred<-read.delim("intersect_30bp_2kbdown_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_cold_pred<-read.delim("intersect_40bp_2kbdown_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_cold_pred<-read.delim("intersect_50bp_2kbdown_aL_randDSB.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}

##STATISTIC 1: INDEL OVERLAP
#Calculate percentage of Indels that overlap with DSBs
indels_NAMv4<-read.delim(file="indelsLarge.bed")
indels_10bp<-read.delim(file="indelsLarge_10kb.bed")
indels_20bp<-read.delim(file="indelsLarge_20kb.bed")
indels_30bp<-read.delim(file="indelsLarge_30kb.bed")
indels_40bp<-read.delim(file="indelsLarge_40kb.bed")
indels_50bp<-read.delim(file="indelsLarge_50kb.bed")
indels_NAMv4<-indels_NAMv4[which(indels_NAMv4$end-indels_NAMv4$start>500),]
indels_10bp<-indels_NAMv4[!(indels_NAMv4$end-indels_NAMv4$start>10000),]
write.table(indels_NAMv4,"indelsLarge.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(indels_10bp,"indelsLarge_10kb.bed", row.names = F,quote=FALSE, sep = '\t')
#saved old files as indelsLarge2 and indelsLarge2_10kb

percentage_overlap<-function(intersect,indels,type){
  percentage_per_chrom<-c()
  intersect<-intersect %>% distinct()
  intersect<-intersect %>% distinct(A_start, .keep_all = TRUE)
  if(type=="insertions" || type=="deletions"){
    for(i in 1:10){
      percentage_per_chrom[i]<- (nrow(intersect[which(intersect$`A_chrom` == i),])/nrow(indels[which(indels$chrom == i)&&which(indels$type == type),]))*100
    }
  }else{
    for(i in 1:10){
      percentage_per_chrom[i]<- (nrow(intersect[which(intersect$`A_chrom` == i),])/nrow(indels[which(indels$`chrom` == i),]))*100
    }
  }
  print(percentage_per_chrom)
}
{
intersect_overlap<-percentage_overlap(intersect_pred,indels_NAMv4,"else")
intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred,indels_10bp,"else")
intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred,indels_20bp,"else")
intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred,indels_30bp,"else")
intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred,indels_40bp,"else")
intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred,indels_50bp,"else")

intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred,indels_NAMv4,"else")
intersect_10bp_2kbup_overlap<-percentage_overlap(intersect_10bp_2kbup_pred,indels_10bp,"else")
intersect_20bp_2kbup_overlap<-percentage_overlap(intersect_20bp_2kbup_pred,indels_20bp,"else")
intersect_30bp_2kbup_overlap<-percentage_overlap(intersect_30bp_2kbup_pred,indels_30bp,"else")
intersect_40bp_2kbup_overlap<-percentage_overlap(intersect_40bp_2kbup_pred,indels_40bp,"else")
intersect_50bp_2kbup_overlap<-percentage_overlap(intersect_50bp_2kbup_pred,indels_50bp,"else")

intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred,indels_NAMv4,"else")
intersect_10bp_2kbdown_overlap<-percentage_overlap(intersect_10bp_2kbdown_pred,indels_10bp,"else")
intersect_20bp_2kbdown_overlap<-percentage_overlap(intersect_20bp_2kbdown_pred,indels_20bp,"else")
intersect_30bp_2kbdown_overlap<-percentage_overlap(intersect_30bp_2kbdown_pred,indels_30bp,"else")
intersect_40bp_2kbdown_overlap<-percentage_overlap(intersect_40bp_2kbdown_pred,indels_40bp,"else")
intersect_50bp_2kbdown_overlap<-percentage_overlap(intersect_50bp_2kbdown_pred,indels_50bp,"else")

intersect_cold_overlap<-percentage_overlap(intersect_cold_pred,indels_NAMv4,"else")
intersect_10bp_cold_overlap<-percentage_overlap(intersect_10bp_cold_pred,indels_10bp,"else")
intersect_20bp_cold_overlap<-percentage_overlap(intersect_20bp_cold_pred,indels_20bp,"else")
intersect_30bp_cold_overlap<-percentage_overlap(intersect_30bp_cold_pred,indels_30bp,"else")
intersect_40bp_cold_overlap<-percentage_overlap(intersect_40bp_cold_pred,indels_40bp,"else")
intersect_50bp_cold_overlap<-percentage_overlap(intersect_50bp_cold_pred,indels_50bp,"else")

intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred,indels_NAMv4,"else")
intersect_10bp_2kbup_cold_overlap<-percentage_overlap(intersect_10bp_2kbup_cold_pred,indels_10bp,"else")
intersect_20bp_2kbup_cold_overlap<-percentage_overlap(intersect_20bp_2kbup_cold_pred,indels_20bp,"else")
intersect_30bp_2kbup_cold_overlap<-percentage_overlap(intersect_30bp_2kbup_cold_pred,indels_30bp,"else")
intersect_40bp_2kbup_cold_overlap<-percentage_overlap(intersect_40bp_2kbup_cold_pred,indels_40bp,"else")
intersect_50bp_2kbup_cold_overlap<-percentage_overlap(intersect_50bp_2kbup_cold_pred,indels_50bp,"else")

intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred,indels_NAMv4,"else")
intersect_10bp_2kbdown_cold_overlap<-percentage_overlap(intersect_10bp_2kbdown_cold_pred,indels_10bp,"else")
intersect_20bp_2kbdown_cold_overlap<-percentage_overlap(intersect_20bp_2kbdown_cold_pred,indels_20bp,"else")
intersect_30bp_2kbdown_cold_overlap<-percentage_overlap(intersect_30bp_2kbdown_cold_pred,indels_30bp,"else")
intersect_40bp_2kbdown_cold_overlap<-percentage_overlap(intersect_40bp_2kbdown_cold_pred,indels_40bp,"else")
intersect_50bp_2kbdown_cold_overlap<-percentage_overlap(intersect_50bp_2kbdown_cold_pred,indels_50bp,"else")
}

#T-TEST
{
library(matrixStats)
# stat1<-as.data.frame(cbind( predicted=c(intersect_overlap,intersect_2kbup_overlap,intersect_2kbdown_overlap), cold=c(intersect_cold_overlap,intersect_2kbup_cold_overlap,intersect_2kbdown_cold_overlap) ))
# t.test(stat1$predicted , stat1$cold, paired=TRUE, var.equal=TRUE)
# model<-lm(predicted~cold,data=stat1)
# plot(model)

stat1<-as.data.frame(cbind(intersect_overlap,intersect_cold_overlap,intersect_2kbup_overlap,intersect_2kbup_cold_overlap,intersect_2kbdown_overlap,intersect_2kbdown_cold_overlap))
stat1[nrow(stat1) + 1,] <- colMeans(as.matrix(stat1))
stat1[nrow(stat1) + 1,] <- colVars(as.matrix(stat1))
stat1[nrow(stat1) + 1,] <-NA
rownames(stat1)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","mean","variance","p-value")
options(scipen=999)
for(i in seq(from=1, to=ncol(stat1)-1, by=2)){
  ttest<-t.test(stat1[1:10,i] , stat1[1:10,i+1], paired=TRUE, var.equal=TRUE)
  stat1["p-value",i]<-ttest$p.value
}
write.table(stat1,file="stat1_1_indel_DSB.txt",row.names = T,quote=FALSE, sep = '\t')

stat1_2<-as.data.frame(cbind( intersect_10bp_overlap,intersect_10bp_cold_overlap,intersect_20bp_overlap,intersect_20bp_cold_overlap,intersect_30bp_overlap,intersect_30bp_cold_overlap,intersect_40bp_overlap, intersect_40bp_cold_overlap,intersect_50bp_overlap,intersect_50bp_cold_overlap))
stat1_2[nrow(stat1_2) + 1,] <- colMeans(as.matrix(stat1_2))
stat1_2[nrow(stat1_2) + 1,] <- colVars(as.matrix(stat1_2))
stat1_2[nrow(stat1_2) + 1,] <-NA
rownames(stat1_2)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","mean","variance","p-value")
options(scipen=999)
for(i in seq(from=1, to=ncol(stat1_2)-1, by=2)){
  ttest<-t.test(stat1_2[1:10,i] , stat1_2[1:10,i+1], paired=TRUE, var.equal=TRUE)
  stat1_2["p-value",i]<-ttest$p.value
}
write.table(stat1_2,file="stat1_2_indel_DSB.txt",row.names = T,quote=FALSE, sep = '\t')

stat1_3<-as.data.frame(cbind( intersect_10bp_2kbup_overlap,intersect_10bp_2kbup_cold_overlap,intersect_20bp_2kbup_overlap,intersect_20bp_2kbup_cold_overlap,intersect_30bp_2kbup_overlap,intersect_30bp_2kbup_cold_overlap,intersect_40bp_2kbup_overlap, intersect_40bp_2kbup_cold_overlap,intersect_50bp_2kbup_overlap,intersect_50bp_2kbup_cold_overlap))
stat1_3[nrow(stat1_3) + 1,] <- colMeans(as.matrix(stat1_3))
stat1_3[nrow(stat1_3) + 1,] <- colVars(as.matrix(stat1_3))
stat1_3[nrow(stat1_3) + 1,] <-NA
rownames(stat1_3)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","mean","variance","p-value")
options(scipen=999)
for(i in seq(from=1, to=ncol(stat1_3)-1, by=2)){
  ttest<-t.test(stat1_3[1:10,i] , stat1_3[1:10,i+1], paired=TRUE, var.equal=TRUE)
  stat1_3["p-value",i]<-ttest$p.value
}
write.table(stat1_3,file="stat1_3_indel_DSB.txt",row.names = T,quote=FALSE, sep = '\t')

stat1_4<-as.data.frame(cbind( intersect_10bp_2kbdown_overlap,intersect_10bp_2kbdown_cold_overlap,intersect_20bp_2kbdown_overlap,intersect_20bp_2kbdown_cold_overlap,intersect_30bp_2kbdown_overlap,intersect_30bp_2kbdown_cold_overlap,intersect_40bp_2kbdown_overlap, intersect_40bp_2kbdown_cold_overlap,intersect_50bp_2kbdown_overlap,intersect_50bp_2kbdown_cold_overlap))
stat1_4[nrow(stat1_4) + 1,] <- colMeans(as.matrix(stat1_4))
stat1_4[nrow(stat1_4) + 1,] <- colVars(as.matrix(stat1_4))
stat1_4[nrow(stat1_4) + 1,] <-NA
rownames(stat1_4)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","mean","variance","p-value")
options(scipen=999)
for(i in seq(from=1, to=ncol(stat1_4)-1, by=2)){
  ttest<-t.test(stat1_4[1:10,i] , stat1_4[1:10,i+1], paired=TRUE, var.equal=TRUE)
  stat1_4["p-value",i]<-ttest$p.value
}
write.table(stat1_4,file="stat1_4_indel_DSB.txt",row.names = T,quote=FALSE, sep = '\t')
}

##STATISTIC 2: CO OVERLAP
#Calculate percentage of COs that overlap with indels
percentage_overlap<-function(intersect,CO){
  percentage_per_chrom<-c()
  intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
  for(i in 1:10){
    percentage_per_chrom[i]<- (nrow(intersect[which(intersect$`A_chrom` == i),])/nrow(CO[which(CO$`chrom` == i),]))*100 
  }
  print(percentage_per_chrom)
}
{
intersect_overlap<-percentage_overlap(intersect_pred,DSB)
intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred,DSB)
intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred,DSB)
intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred,DSB)
intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred,DSB)
intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred,DSB)

intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred,DSB)
intersect_10bp_2kbup_overlap<-percentage_overlap(intersect_10bp_2kbup_pred,DSB)
intersect_20bp_2kbup_overlap<-percentage_overlap(intersect_20bp_2kbup_pred,DSB)
intersect_30bp_2kbup_overlap<-percentage_overlap(intersect_30bp_2kbup_pred,DSB)
intersect_40bp_2kbup_overlap<-percentage_overlap(intersect_40bp_2kbup_pred,DSB)
intersect_50bp_2kbup_overlap<-percentage_overlap(intersect_50bp_2kbup_pred,DSB)

intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred,DSB)
intersect_10bp_2kbdown_overlap<-percentage_overlap(intersect_10bp_2kbdown_pred,DSB)
intersect_20bp_2kbdown_overlap<-percentage_overlap(intersect_20bp_2kbdown_pred,DSB)
intersect_30bp_2kbdown_overlap<-percentage_overlap(intersect_30bp_2kbdown_pred,DSB)
intersect_40bp_2kbdown_overlap<-percentage_overlap(intersect_40bp_2kbdown_pred,DSB)
intersect_50bp_2kbdown_overlap<-percentage_overlap(intersect_50bp_2kbdown_pred,DSB)

intersect_cold_overlap<-percentage_overlap(intersect_cold_pred,DSB_coldspots)
intersect_10bp_cold_overlap<-percentage_overlap(intersect_10bp_cold_pred,DSB_coldspots)
intersect_20bp_cold_overlap<-percentage_overlap(intersect_20bp_cold_pred,DSB_coldspots)
intersect_30bp_cold_overlap<-percentage_overlap(intersect_30bp_cold_pred,DSB_coldspots)
intersect_40bp_cold_overlap<-percentage_overlap(intersect_40bp_cold_pred,DSB_coldspots)
intersect_50bp_cold_overlap<-percentage_overlap(intersect_50bp_cold_pred,DSB_coldspots)

intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred,DSB_cold_2kb_up)
intersect_10bp_2kbup_cold_overlap<-percentage_overlap(intersect_10bp_2kbup_cold_pred,DSB_cold_2kb_up)
intersect_20bp_2kbup_cold_overlap<-percentage_overlap(intersect_20bp_2kbup_cold_pred,DSB_cold_2kb_up)
intersect_30bp_2kbup_cold_overlap<-percentage_overlap(intersect_30bp_2kbup_cold_pred,DSB_cold_2kb_up)
intersect_40bp_2kbup_cold_overlap<-percentage_overlap(intersect_40bp_2kbup_cold_pred,DSB_cold_2kb_up)
intersect_50bp_2kbup_cold_overlap<-percentage_overlap(intersect_50bp_2kbup_cold_pred,DSB_cold_2kb_up)

intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred,DSB_cold_2kb_down)
intersect_10bp_2kbdown_cold_overlap<-percentage_overlap(intersect_10bp_2kbdown_cold_pred,DSB_cold_2kb_down)
intersect_20bp_2kbdown_cold_overlap<-percentage_overlap(intersect_20bp_2kbdown_cold_pred,DSB_cold_2kb_down)
intersect_30bp_2kbdown_cold_overlap<-percentage_overlap(intersect_30bp_2kbdown_cold_pred,DSB_cold_2kb_down)
intersect_40bp_2kbdown_cold_overlap<-percentage_overlap(intersect_40bp_2kbdown_cold_pred,DSB_cold_2kb_down)
intersect_50bp_2kbdown_cold_overlap<-percentage_overlap(intersect_50bp_2kbdown_cold_pred,DSB_cold_2kb_down)
}
#overlap % stats:
{
  {
    hot<-cbind(intersect_overlap, c(rep("hot",10)))
    cold<-cbind(intersect_cold_overlap, c(rep("cold",10)))
    stat1<-as.data.frame(rbind(hot,cold))
    colnames(stat1)<-c("y","type")
    data1 <- data.frame(y=stat1$y,group=factor(stat1$type))
    fit = lm(y ~ group, data1)
    one.way <- aov(fit)
    summary(one.way)
    
  }
  #final frame:
  {
    c1<-mean(as.numeric(stat1[which(stat1$type=="hot"),]$y))
    c2<-mean(as.numeric(stat1[which(stat1$type=="cold"),]$y))
    c1
    c2
  }
}
##STATISTIC 3: indel density
library(dplyr)
#mean indel density on each CO
percentage_overlap<-function(intersect){
  colnames(intersect) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect<-intersect %>% distinct(A_start, .keep_all = TRUE)
  duplicated<-as.data.frame(table(intersect$'B_start'))
  intersect$length<-intersect$'B_end'-intersect$'B_start'
  intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
  intersect<-intersect[order(intersect$B_start),]
  intersect$density<-duplicated$Freq
  intersect$percent<- (intersect$'density'/intersect$length)
  percentage_per_chrom<-c()
  for(i in 1:10){
    intersect_temp<-intersect[which(intersect$`A_chrom` == i),]
    if(is.na(mean(intersect_temp$percent))){
      percentage_per_chrom[i]<-0
    }else{
      percentage_per_chrom[i]<- mean(intersect_temp$percent)
    }
    
  }
  print(percentage_per_chrom)
}
{
intersect_overlap<-percentage_overlap(intersect_pred)
intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred)
intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred)
intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred)
intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred)
intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred)

intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred)
intersect_10bp_2kbup_overlap<-percentage_overlap(intersect_10bp_2kbup_pred)
intersect_20bp_2kbup_overlap<-percentage_overlap(intersect_20bp_2kbup_pred)
intersect_30bp_2kbup_overlap<-percentage_overlap(intersect_30bp_2kbup_pred)
intersect_40bp_2kbup_overlap<-percentage_overlap(intersect_40bp_2kbup_pred)
intersect_50bp_2kbup_overlap<-percentage_overlap(intersect_50bp_2kbup_pred)

intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred)
intersect_10bp_2kbdown_overlap<-percentage_overlap(intersect_10bp_2kbdown_pred)
intersect_20bp_2kbdown_overlap<-percentage_overlap(intersect_20bp_2kbdown_pred)
intersect_30bp_2kbdown_overlap<-percentage_overlap(intersect_30bp_2kbdown_pred)
intersect_40bp_2kbdown_overlap<-percentage_overlap(intersect_40bp_2kbdown_pred)
intersect_50bp_2kbdown_overlap<-percentage_overlap(intersect_50bp_2kbdown_pred)

intersect_cold_overlap<-percentage_overlap(intersect_cold_pred)
intersect_10bp_cold_overlap<-percentage_overlap(intersect_10bp_cold_pred)
intersect_20bp_cold_overlap<-percentage_overlap(intersect_20bp_cold_pred)
intersect_30bp_cold_overlap<-percentage_overlap(intersect_30bp_cold_pred)
intersect_40bp_cold_overlap<-percentage_overlap(intersect_40bp_cold_pred)
intersect_50bp_cold_overlap<-percentage_overlap(intersect_50bp_cold_pred)

intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred)
intersect_10bp_2kbup_cold_overlap<-percentage_overlap(intersect_10bp_2kbup_cold_pred)
intersect_20bp_2kbup_cold_overlap<-percentage_overlap(intersect_20bp_2kbup_cold_pred)
intersect_30bp_2kbup_cold_overlap<-percentage_overlap(intersect_30bp_2kbup_cold_pred)
intersect_40bp_2kbup_cold_overlap<-percentage_overlap(intersect_40bp_2kbup_cold_pred)
intersect_50bp_2kbup_cold_overlap<-percentage_overlap(intersect_50bp_2kbup_cold_pred)

intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred)
intersect_10bp_2kbdown_cold_overlap<-percentage_overlap(intersect_10bp_2kbdown_cold_pred)
intersect_20bp_2kbdown_cold_overlap<-percentage_overlap(intersect_20bp_2kbdown_cold_pred)
intersect_30bp_2kbdown_cold_overlap<-percentage_overlap(intersect_30bp_2kbdown_cold_pred)
intersect_40bp_2kbdown_cold_overlap<-percentage_overlap(intersect_40bp_2kbdown_cold_pred)
intersect_50bp_2kbdown_cold_overlap<-percentage_overlap(intersect_50bp_2kbdown_cold_pred)
}

#Stats: ANOVA, tukey/bonferonni correction
#ANOVA test for significance
{
  percentage_overlap<-function(intersect,dataset){
    intersect<-intersect %>% distinct(A_start, .keep_all = TRUE)
    duplicated<-as.data.frame(table(intersect$'B_start'))
    intersect$length<-intersect$'B_end'-intersect$'B_start'
    intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
    intersect<-intersect[order(intersect$B_start),]
    intersect$count<-duplicated$Freq
    intersect$percent<- (intersect$count/intersect$length)
    intersect$dataset<-dataset
    return(intersect)
  }
  {
    intersect_overlap<-percentage_overlap(intersect_pred,"hot")
    intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred,"hot")
    intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred,"hot")
    intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred,"hot")
    intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred,"hot")
    intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred,"hot")
    
    intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred,"hot")
    intersect_10bp_2kbup_overlap<-percentage_overlap(intersect_10bp_2kbup_pred,"hot")
    intersect_20bp_2kbup_overlap<-percentage_overlap(intersect_20bp_2kbup_pred,"hot")
    intersect_30bp_2kbup_overlap<-percentage_overlap(intersect_30bp_2kbup_pred,"hot")
    intersect_40bp_2kbup_overlap<-percentage_overlap(intersect_40bp_2kbup_pred,"hot")
    intersect_50bp_2kbup_overlap<-percentage_overlap(intersect_50bp_2kbup_pred,"hot")
    
    intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred,"hot")
    intersect_10bp_2kbdown_overlap<-percentage_overlap(intersect_10bp_2kbdown_pred,"hot")
    intersect_20bp_2kbdown_overlap<-percentage_overlap(intersect_20bp_2kbdown_pred,"hot")
    intersect_30bp_2kbdown_overlap<-percentage_overlap(intersect_30bp_2kbdown_pred,"hot")
    intersect_40bp_2kbdown_overlap<-percentage_overlap(intersect_40bp_2kbdown_pred,"hot")
    intersect_50bp_2kbdown_overlap<-percentage_overlap(intersect_50bp_2kbdown_pred,"hot")
    
    intersect_cold_overlap<-percentage_overlap(intersect_cold_pred,"cold")
    intersect_10bp_cold_overlap<-percentage_overlap(intersect_10bp_cold_pred,"cold")
    intersect_20bp_cold_overlap<-percentage_overlap(intersect_20bp_cold_pred,"cold")
    intersect_30bp_cold_overlap<-percentage_overlap(intersect_30bp_cold_pred,"cold")
    intersect_40bp_cold_overlap<-percentage_overlap(intersect_40bp_cold_pred,"cold")
    intersect_50bp_cold_overlap<-percentage_overlap(intersect_50bp_cold_pred,"cold")
    
    intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred,"cold")
    intersect_10bp_2kbup_cold_overlap<-percentage_overlap(intersect_10bp_2kbup_cold_pred,"cold")
    intersect_20bp_2kbup_cold_overlap<-percentage_overlap(intersect_20bp_2kbup_cold_pred,"cold")
    intersect_30bp_2kbup_cold_overlap<-percentage_overlap(intersect_30bp_2kbup_cold_pred,"cold")
    intersect_40bp_2kbup_cold_overlap<-percentage_overlap(intersect_40bp_2kbup_cold_pred,"cold")
    intersect_50bp_2kbup_cold_overlap<-percentage_overlap(intersect_50bp_2kbup_cold_pred,"cold")
    
    intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred,"cold")
    intersect_10bp_2kbdown_cold_overlap<-percentage_overlap(intersect_10bp_2kbdown_cold_pred,"cold")
    intersect_20bp_2kbdown_cold_overlap<-percentage_overlap(intersect_20bp_2kbdown_cold_pred,"cold")
    intersect_30bp_2kbdown_cold_overlap<-percentage_overlap(intersect_30bp_2kbdown_cold_pred,"cold")
    intersect_40bp_2kbdown_cold_overlap<-percentage_overlap(intersect_40bp_2kbdown_cold_pred,"cold")
    intersect_50bp_2kbdown_cold_overlap<-percentage_overlap(intersect_50bp_2kbdown_cold_pred,"cold")
  }
  
  #STATISTIC 1
  stat2<-as.data.frame(rbind(intersect_overlap,intersect_10bp_overlap,intersect_20bp_overlap,intersect_30bp_overlap,intersect_40bp_overlap,intersect_50bp_overlap,intersect_cold_overlap,intersect_10bp_cold_overlap,intersect_20bp_cold_overlap,intersect_30bp_cold_overlap,intersect_40bp_cold_overlap,intersect_50bp_cold_overlap))
  data1 <- data.frame(y=(stat2$count),group=factor(stat2$dataset),group2=factor(stat2$type))
  data1<- subset(data1, group2!="all deletions")
  fit= glm(y~group * group2, data1, family = poisson(link = "log"))
  #* interaction effect- impact of hot/cold depends on indel size
  two.way <- aov(fit)
  summary(two.way)
  
  #bonferroni for non normal distribution
  data1$group3<-paste(data1$group,data1$group2)
  pairwise.t.test(data1$y,data1$group, p.adjust.method="bonferroni")
  pairwise.t.test(data1$y,data1$group2, p.adjust.method="bonferroni")
  
  ggplot(data1, aes(x=factor(group2), y=y, fill=group)) +
    geom_boxplot()+theme_classic2()
  ggboxplot(data1, x = "group2", y = "y",
            color = "group", palette = "jco")+
    stat_compare_means(method = "anova")+theme_cleveland()
  
  
  tukey.test <- TukeyHSD(two.way)
  tukey.test
  plot(tukey.test)
  #ANOVA test indicated significant diff (F=0.15, p=0.861)
  
  #STATISTIC 2
  stat2<-as.data.frame(rbind(intersect_overlap,intersect_2kbup_overlap,intersect_2kbdown_overlap,intersect_cold_overlap,intersect_2kbup_cold_overlap,intersect_2kbdown_cold_overlap))
  data1 <- data.frame(y=(stat2$count),group=factor(stat2$dataset),group2=factor(stat2$intersection))
  fit= glm(y~group * group2, data1, family = poisson(link = "log"))
  #* interaction effect- impact of hot/cold depends on indel size
  two.way <- aov(fit)
  summary(two.way)
  
  #bonferroni for non normal distribution
  data1$group3<-paste(data1$group,data1$group2)
  pairwise.t.test(data1$y,data1$group, p.adjust.method="bonferroni")
  pairwise.t.test(data1$y,data1$group2, p.adjust.method="bonferroni")
  pairwise.t.test(data1$y,data1$group3, p.adjust.method="bonferroni")
  
  tukey.test <- TukeyHSD(two.way)
  tukey.test
  plot(tukey.test)
  
  library(multcompView)
  cld <- multcompLetters4(two.way, tukey.test)
  
  library(car)
  leveneTest(fit)
  #Levene's test indicated equal variances (F = 2.5291, p = 0.03962)
  
  #anderson-darling normality test
  library(nortest)
  ad.test(data1$y) #departure from normality
  
}
#BARGRAPHS - for each: change TITLE, INDEL SIZES, LEGEND, AXIS TITLE POSITION
{
  #create dataframes:
  {
    #DATAFRAME OF CO PREDICTED
    type <- c(rep(c("all indels") ,3),rep(c("1-10bp") ,3),rep(c("11-20bp") ,3),rep(c("21-3bp") ,3),rep(c("31-40bp") ,3),rep(c("41-50bp") ,3))
    intersect<-c(c("Interval"),c("2kb upstream") ,c("2kb downstream"))
    #intersect<-c(c("Interval"),c("2kb upstream") ,c("2kb downstream"),c("Interval"),c("2kb upstream") ,c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10))
    value <- c(mean(intersect_overlap),mean(intersect_2kbup_overlap),mean(intersect_2kbdown_overlap),mean(intersect_10bp_overlap),mean(intersect_10bp_2kbup_overlap),mean(intersect_10bp_2kbdown_overlap),mean(intersect_20bp_overlap),mean(intersect_20bp_2kbup_overlap),mean(intersect_20bp_2kbdown_overlap),mean(intersect_30bp_overlap),mean(intersect_30bp_2kbup_overlap),mean(intersect_30bp_2kbdown_overlap),mean(intersect_40bp_overlap),mean(intersect_40bp_2kbup_overlap),mean(intersect_40bp_2kbdown_overlap),mean(intersect_50bp_overlap),mean(intersect_50bp_2kbup_overlap),mean(intersect_50bp_2kbdown_overlap))
    sd <- c(sd(intersect_overlap)/sqrt(length((intersect_overlap))),sd(intersect_2kbup_overlap)/sqrt(length((intersect_2kbup_overlap))),sd(intersect_2kbdown_overlap)/sqrt(length((intersect_2kbdown_overlap))),sd(intersect_10bp_overlap)/sqrt(length((intersect_10bp_overlap))),sd(intersect_10bp_2kbup_overlap)/sqrt(length((intersect_10bp_2kbup_overlap))),sd(intersect_10bp_2kbdown_overlap)/sqrt(length((intersect_10bp_2kbdown_overlap))),sd(intersect_20bp_overlap)/sqrt(length((intersect_20bp_overlap))),sd(intersect_20bp_2kbup_overlap)/sqrt(length((intersect_20bp_2kbup_overlap))),sd(intersect_20bp_2kbdown_overlap)/sqrt(length((intersect_20bp_2kbdown_overlap))),sd(intersect_30bp_overlap)/sqrt(length((intersect_30bp_overlap))),sd(intersect_30bp_2kbup_overlap)/sqrt(length((intersect_30bp_2kbup_overlap))),sd(intersect_30bp_2kbdown_overlap)/sqrt(length((intersect_30bp_2kbdown_overlap))),sd(intersect_40bp_overlap)/sqrt(length((intersect_40bp_overlap))),sd(intersect_40bp_2kbup_overlap)/sqrt(length((intersect_40bp_2kbup_overlap))),sd(intersect_40bp_2kbdown_overlap)/sqrt(length((intersect_40bp_2kbdown_overlap))),sd(intersect_50bp_overlap)/sqrt(length((intersect_50bp_overlap))),sd(intersect_50bp_2kbup_overlap)/sqrt(length((intersect_50bp_2kbup_overlap))),sd(intersect_50bp_2kbdown_overlap)/sqrt(length((intersect_50bp_2kbdown_overlap))))
    data <- data.frame(type,value,intersect,sd)
    data$dataset<-"hot"
    
    #DATAFRAME OF CO DESERT
    value <- c(mean(intersect_cold_overlap),mean(intersect_2kbup_cold_overlap),mean(intersect_2kbdown_cold_overlap),mean(intersect_10bp_cold_overlap),mean(intersect_10bp_2kbup_cold_overlap),mean(intersect_10bp_2kbdown_cold_overlap),mean(intersect_20bp_cold_overlap),mean(intersect_20bp_2kbup_cold_overlap),mean(intersect_20bp_2kbdown_cold_overlap),mean(intersect_30bp_cold_overlap),mean(intersect_30bp_2kbup_cold_overlap),mean(intersect_30bp_2kbdown_cold_overlap),mean(intersect_40bp_cold_overlap),mean(intersect_40bp_2kbup_cold_overlap),mean(intersect_40bp_2kbdown_cold_overlap),mean(intersect_50bp_cold_overlap),mean(intersect_50bp_2kbup_cold_overlap),mean(intersect_50bp_2kbdown_cold_overlap))
    sd <- c(sd(intersect_cold_overlap)/sqrt(length((intersect_cold_overlap))),sd(intersect_2kbup_cold_overlap)/sqrt(length((intersect_2kbup_cold_overlap))),sd(intersect_2kbdown_cold_overlap)/sqrt(length((intersect_2kbdown_cold_overlap))),sd(intersect_10bp_cold_overlap)/sqrt(length((intersect_10bp_cold_overlap))),sd(intersect_10bp_2kbup_cold_overlap)/sqrt(length((intersect_10bp_2kbup_cold_overlap))),sd(intersect_10bp_2kbdown_cold_overlap)/sqrt(length((intersect_10bp_2kbdown_cold_overlap))),sd(intersect_20bp_cold_overlap)/sqrt(length((intersect_20bp_cold_overlap))),sd(intersect_20bp_2kbup_cold_overlap)/sqrt(length((intersect_20bp_2kbup_cold_overlap))),sd(intersect_20bp_2kbdown_cold_overlap)/sqrt(length((intersect_20bp_2kbdown_cold_overlap))),sd(intersect_30bp_cold_overlap)/sqrt(length((intersect_30bp_cold_overlap))),sd(intersect_30bp_2kbup_cold_overlap)/sqrt(length((intersect_30bp_2kbup_cold_overlap))),sd(intersect_30bp_2kbdown_cold_overlap)/sqrt(length((intersect_30bp_2kbdown_cold_overlap))),sd(intersect_40bp_cold_overlap)/sqrt(length((intersect_40bp_cold_overlap))),sd(intersect_40bp_2kbup_cold_overlap)/sqrt(length((intersect_40bp_2kbup_cold_overlap))),sd(intersect_40bp_2kbdown_cold_overlap)/sqrt(length((intersect_40bp_2kbdown_cold_overlap))),sd(intersect_50bp_cold_overlap)/sqrt(length((intersect_50bp_cold_overlap))),sd(intersect_50bp_2kbup_cold_overlap)/sqrt(length((intersect_50bp_2kbup_cold_overlap))),sd(intersect_50bp_2kbdown_cold_overlap)/sqrt(length((intersect_50bp_2kbdown_cold_overlap))))
    data_cold <- data.frame(type,value,intersect,sd)
    data_cold[is.na(data_cold)]<-0
    data_cold$dataset<-"cold"
    full_data<-rbind(data,data_cold)
    full_data<-full_data[order(full_data$type,full_data$dataset),]
    full_data$levels<-interaction(full_data$type,full_data$dataset)
    full_data$levels <- factor(full_data$levels , levels = unique(full_data$levels ))
    full_data$color<-paste(full_data$intersect,full_data$dataset)
  }
  library(ggplot2)
  library(viridis)
  library(hrbrthemes)
  dev.off()
  full_data1<-full_data[which(full_data$intersect=="Interval"),]
  full_data1<- subset(full_data1, type!="all indels")
  full_data1$posthoc<-c("","","","","","","","","","")
  #plot 1: indel density 
  perc<-ggplot(full_data1, aes(fill=color,x=levels,y=value)) +
    geom_bar(stat="identity",position="dodge") +
    scale_fill_manual(labels = c("DSB",
                                 "Negative Control Region"),
                      values = c("Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("Large Indel Density on DSB") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(),width=.2)+
    geom_text(aes(label=posthoc, y=value+sd), position=position_dodge(.5), vjust=-.5) +
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(0,max(full_data1$value + full_data1$sd))) +
    annotate("text", c(1.5,3.5,5.5,7.5,9.5), y = -.00005, label = c("500bp-10 kb","10-20 kb","20-30 kb","30-40 kb","40-50 kb")) +
    annotate("text", x = 5.5, y = -.0001,
             label = c("Dataset and Indel Size"),fontface =2) +
    ylab("Indel Density (indels/bp)")+
    guides(fill=guide_legend("Dataset"))+
    theme_classic() +
    theme(plot.margin = unit(c(3, 1, 5, 1), "lines"),
          plot.title = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face="bold"),
          legend.position='top',legend.spacing.x= unit(.5, 'cm'))
  #plot margins unit(c(top, right, bottom, left),
  perc
  g2 <- ggplot_gtable(ggplot_build(perc))
  g2$layout$clip[g2$layout$name == "panel"] <- "off"
  grid.draw(g2)
  
  full_data2<-full_data[which(full_data$type=="all indels"),]
  full_data2$posthoc<-c("","","","","","")
  perc2<-ggplot(full_data2, aes(fill=color,x=levels,y=value)) +
    geom_bar(stat="identity",position="dodge", width = .5) +
    scale_fill_manual(labels = c("2kb downstream", "2kb downstream",
                                 "2kb upstream","2kb upstream",
                                 "DSB",
                                 "Negative Control"),
                      values = c("2kb downstream hot"= "darkgoldenrod1",
                                 "2kb downstream cold" = "deepskyblue",
                                 "2kb upstream hot" = "darkorange1",
                                 "2kb upstream cold" = "dodgerblue2",
                                 "Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("Large Indel Density on DSB") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(.5),width=.2)+
    geom_text(aes(label=posthoc, y=value+sd), position=position_dodge(.5), vjust=-.5) +
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(0,max(full_data2$value + full_data2$sd))) +
    annotate("text", x = 1:2, y = -.00005,
             label = rep(c("Negative Control Region", "DSB Region"), 1)) +
    annotate("text", x = 1.5, y = -.0001,
             label = c("Dataset"),fontface =2) +
    ylab("Indel Density (indels/bp)")+
    guides(fill=guide_legend("Dataset"))+
    theme_classic() +
    theme(plot.margin = unit(c(3, 1, 5, 1), "lines"),
          plot.title = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_text(face="bold"),
          legend.position='top',legend.spacing.x= unit(.5, 'cm'))
  #plot margins unit(c(top, right, bottom, left),
  perc2
  g2 <- ggplot_gtable(ggplot_build(perc2))
  g2$layout$clip[g2$layout$name == "panel"] <- "off"
  grid.draw(g2)
}


#SUMMARY TABLE OF ALL RESULTS
test <- matrix(nrow=6)
rownames(test)<-c("All indels","1-10bp","11-20bp","21-30bp","31-40bp","41-50bp")
test<-cbind(test,"Mean Indel / DSB Overlap"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Indel / Random Site Overlap"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test<-cbind(test,"Mean DSB Overlap"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Random Site Overlap"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test<-cbind(test,"Mean Indel Density on DSB (indels/bp)"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Indel Density on Random Site (indels/bp)"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test_2<-as.data.frame(test)[-c(1)]
write.table(test_2,file="DSBs_aL_output.txt",row.names = T,quote=FALSE, sep = '\t')


#bargraphs of data
library(ggplot2)
library(viridis)
library(hrbrthemes)

perc<-ggplot(data, aes(fill=intersect,x=reorder(type, value),y=value)) +
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual(values = c("Interval" = "tan1",
                               "2kb upstream" = "sienna1",
                               "2kb downstream" = "brown1"))+
  facet_wrap(~chrom,ncol=5) +
  ggtitle("Indel Density \nOn DSBs Per Chromosome",subtitle = "Average Density: 0.26") +
  theme(plot.title = element_text(face="bold", size=20, hjust=0, color="#555555")) +
  geom_text(aes(label = round(value,digits=4)),position = position_stack(vjust = 0.5),size=2.5)+
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.position='top',legend.spacing.x= unit(1.0, 'cm')) +
  ylab("Percentage Overlap")+
  xlab("Type")+guides(fill=guide_legend("Type of intersection"))
perc

perc_cold<-ggplot(data_cold, aes(fill=intersect,x=reorder(type, value),y=value)) +
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual(values = c("Interval" = "turquoise2",
                               "2kb upstream" = "deepskyblue",
                               "2kb downstream" = "blue"))+
  facet_wrap(~chrom,ncol=5) +
  ggtitle("Indel Density \nOn Random Sites Per Chromosome",subtitle = "Average Density: 0.002") +
  theme(plot.title = element_text(face="bold", size=20, hjust=0, color="#555555")) +
  geom_text(aes(label = round(value,digits=4)),position = position_stack(vjust = 0.5),size=2.5)+
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.position='top',legend.spacing.x= unit(1.0, 'cm')) +
  ylab("Percentage Overlap")+
  xlab("Type")+guides(fill=guide_legend("Type of intersection"))
perc_cold

##STATISTIC 4: indel density at specific sites
library(dplyr)
#mean indel density on each CO
percentage_overlap<-function(intersect){
  intersect<-intersect %>% distinct()
  duplicated<-as.data.frame(table(intersect$'B_start'))
  intersect$length<-intersect$'B_end'-intersect$'B_start'
  intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
  intersect$density<-duplicated$Freq
  intersect$percent<- (intersect$'density'/intersect$length)
  return(intersect)
}
{
  intersect_overlap<-percentage_overlap(intersect_pred)
  intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred)
  intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred)
  all<-rbind(intersect_overlap,intersect_2kbup_overlap,intersect_2kbdown_overlap)
  all<-all[order(all$B_chrom,all$B_start),]
  all<-all[c(1:6,9,12)]
  insert_DSB<-function(all){
    for(i in seq(from=1, to=nrow(all)-1, by=2)){
      if(all$intersection=="2kb up"){
        
      }
    }
  }
  all_1<-all[ which(all$A_chrom==1),]
  all_2<-all[ which(all$A_chrom==2),]  
  all_2$A_start<-min(all_2$A_start)
  all_2$A_end<-max(all_2$A_end)
  M=matrix((all_2$B_start+all_2$B_end)/2- all_2$A_start,ncol=1)
  
  all_2=cbind(all_2,M)
  
  #all_2=all_2[all_2$M<10000,]
  
  tapply(all_2$percent, cut(all_2$M, seq(0, max(all_2$B_end), by=100)), mean)->mean_all_2
  #mean_all_2<-na.omit(mean_all_2)
  plot (smooth.spline(seq(-5000,4999,100),mean_all_2,spar=0.3),col="red",lwd=3,type="l",xlab ="Distance from midpoint of DSBs", ylab="Indel Density (indels/bp)")

  # plot (smooth.spline(seq(-5000,4999,100),mean_all_2,spar=0.35),col="red",lwd=3,type="l",ylim=c(0,120),xlab ="Distance from midpoint of COs", ylab="DNA methylation (%)")
  # 
  intersect_cold_overlap<-percentage_overlap(intersect_cold_pred)
  intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred)
  intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred)
 
}
