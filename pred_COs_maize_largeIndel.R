
#QUICK RELOAD
CO_predicted<-read.delim("CO_predicted.bed")
CO_pred_cold_2kb_up<-read.delim("CO_pred_cold_2kb_up.bed")
CO_pred_cold_2kb_down<-read.delim("CO_pred_cold_2kb_down.bed")
CO_pred_coldspots<-read.delim("CO_pred_coldspots.bed")
CO_pred_2kb_up<-read.delim("CO_pred_2kb_up.bed")
CO_pred_2kb_down<-read.delim("CO_pred_2kb_down.bed")

#clustered datasets
CO_predicted<-read.delim("CO_predicted_3.bed",header=FALSE)
colnames(CO_predicted)<-c('chrom','start','end')
CO_pred_cold_2kb_up<-read.delim("CO_pred_cold_2kb_up_2.bed")
CO_pred_cold_2kb_down<-read.delim("CO_pred_cold_2kb_down_2.bed")
CO_pred_coldspots<-read.delim("CO_pred_coldspots_2.bed")
CO_pred_2kb_up<-read.delim("CO_pred_2kb_up_2.bed")
CO_pred_2kb_down<-read.delim("CO_pred_2kb_down_2.bed")

#Read in overlap data, calculate percentage of hotspots that overlap
##READ IN INTERSECT DATA:
# INTERSECT WITH PREDICTED CO REGIONS
{
intersect_pred<-read.delim("intersect_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_pred$type<-"all deletions"
intersect_pred$intersect<-"intersect"
colnames(intersect_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_10bp_pred<-read.delim("intersect_10bp_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_10bp_pred$type<-"1-10bp"
intersect_10bp_pred$intersect<-"intersect"
colnames(intersect_10bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_20bp_pred<-read.delim("intersect_20bp_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_20bp_pred$type<-"11-20bp"
intersect_20bp_pred$intersect<-"intersect"
colnames(intersect_20bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_30bp_pred<-read.delim("intersect_30bp_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_30bp_pred$type<-"21-30bp"
intersect_30bp_pred$intersect<-"intersect"
colnames(intersect_30bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_40bp_pred<-read.delim("intersect_40bp_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_40bp_pred$type<-"31-40bp"
intersect_40bp_pred$intersect<-"intersect"
colnames(intersect_40bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_50bp_pred<-read.delim("intersect_50bp_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_50bp_pred$type<-"41-50bp"
intersect_50bp_pred$intersect<-"intersect"
colnames(intersect_50bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')

intersect_2kbup_pred<-read.delim("intersect_2kbup_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_2kbup_pred$type<-"all deletions"
intersect_2kbup_pred$intersect<-"2kb up"
colnames(intersect_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_10bp_2kbup_pred<-read.delim("intersect_10bp_2kbup_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_10bp_2kbup_pred$type<-"1-10bp"
intersect_10bp_2kbup_pred$intersect<-"2kb up"
colnames(intersect_10bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_20bp_2kbup_pred<-read.delim("intersect_20bp_2kbup_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_20bp_2kbup_pred$type<-"11-20bp"
intersect_20bp_2kbup_pred$intersect<-"2kb up"
colnames(intersect_20bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_30bp_2kbup_pred<-read.delim("intersect_30bp_2kbup_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_30bp_2kbup_pred$type<-"21-30bp"
intersect_30bp_2kbup_pred$intersect<-"2kb up"
colnames(intersect_30bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_40bp_2kbup_pred<-read.delim("intersect_40bp_2kbup_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_40bp_2kbup_pred$type<-"31-40bp"
intersect_40bp_2kbup_pred$intersect<-"2kb up"
colnames(intersect_40bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_50bp_2kbup_pred<-read.delim("intersect_50bp_2kbup_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_50bp_2kbup_pred$type<-"41-50bp"
intersect_50bp_2kbup_pred$intersect<-"2kb up"
colnames(intersect_50bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')

intersect_2kbdown_pred<-read.delim("intersect_2kbdown_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_2kbdown_pred$type<-"all deletions"
intersect_2kbdown_pred$intersect<-"2kb down"
colnames(intersect_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_10bp_2kbdown_pred<-read.delim("intersect_10bp_2kbdown_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_10bp_2kbdown_pred$type<-"1-10bp"
intersect_10bp_2kbdown_pred$intersect<-"2kb down"
colnames(intersect_10bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_20bp_2kbdown_pred<-read.delim("intersect_20bp_2kbdown_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_20bp_2kbdown_pred$type<-"11-20bp"
intersect_20bp_2kbdown_pred$intersect<-"2kb down"
colnames(intersect_20bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_30bp_2kbdown_pred<-read.delim("intersect_30bp_2kbdown_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_30bp_2kbdown_pred$type<-"21-30bp"
intersect_30bp_2kbdown_pred$intersect<-"2kb down"
colnames(intersect_30bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_40bp_2kbdown_pred<-read.delim("intersect_40bp_2kbdown_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_40bp_2kbdown_pred$type<-"31-40bp"
intersect_40bp_2kbdown_pred$intersect<-"2kb down"
colnames(intersect_40bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_50bp_2kbdown_pred<-read.delim("intersect_50bp_2kbdown_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_50bp_2kbdown_pred$type<-"41-50bp"
intersect_50bp_2kbdown_pred$intersect<-"2kb down"
colnames(intersect_50bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}
# INTERSECT WITH RANDOM CO DESERT REGIONS
{
intersect_cold_pred<-read.delim("intersect_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_cold_pred$type<-"all deletions"
intersect_cold_pred$intersect<-"intersect"
intersect_cold_pred<-intersect_cold_pred[-c(7:8)]
colnames(intersect_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_10bp_cold_pred<-read.delim("intersect_10bp_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_10bp_cold_pred$type<-"1-10bp"
intersect_10bp_cold_pred$intersect<-"intersect"
intersect_10bp_cold_pred<-intersect_10bp_cold_pred[-c(7:8)]
colnames(intersect_10bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_20bp_cold_pred<-read.delim("intersect_20bp_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_20bp_cold_pred$type<-"11-20bp"
intersect_20bp_cold_pred$intersect<-"intersect"
intersect_20bp_cold_pred<-intersect_20bp_cold_pred[-c(7:8)]
colnames(intersect_20bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_30bp_cold_pred<-read.delim("intersect_30bp_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_30bp_cold_pred$type<-"21-30bp"
intersect_30bp_cold_pred$intersect<-"intersect"
intersect_30bp_cold_pred<-intersect_30bp_cold_pred[-c(7:8)]
colnames(intersect_30bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_40bp_cold_pred<-read.delim("intersect_40bp_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_40bp_cold_pred$type<-"31-40bp"
intersect_40bp_cold_pred$intersect<-"intersect"
intersect_40bp_cold_pred<-intersect_40bp_cold_pred[-c(7:8)]
colnames(intersect_40bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_50bp_cold_pred<-read.delim("intersect_50bp_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_50bp_cold_pred$type<-"41-50bp"
intersect_50bp_cold_pred$intersect<-"intersect"
intersect_50bp_cold_pred<-intersect_50bp_cold_pred[-c(7:8)]
colnames(intersect_50bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')

intersect_2kbup_cold_pred<-read.delim("intersect_2kbup_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_2kbup_cold_pred$type<-"all deletions"
intersect_2kbup_cold_pred$intersect<-"2kb up"
colnames(intersect_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_10bp_2kbup_cold_pred<-read.delim("intersect_10bp_2kbup_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_10bp_2kbup_cold_pred$type<-"1-10bp"
intersect_10bp_2kbup_cold_pred$intersect<-"2kb up"
colnames(intersect_10bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_20bp_2kbup_cold_pred<-read.delim("intersect_20bp_2kbup_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_20bp_2kbup_cold_pred$type<-"11-20bp"
intersect_20bp_2kbup_cold_pred$intersect<-"2kb up"
colnames(intersect_20bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_30bp_2kbup_cold_pred<-read.delim("intersect_30bp_2kbup_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_30bp_2kbup_cold_pred$type<-"21-30bp"
intersect_30bp_2kbup_cold_pred$intersect<-"2kb up"
colnames(intersect_30bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_40bp_2kbup_cold_pred<-read.delim("intersect_40bp_2kbup_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_40bp_2kbup_cold_pred$type<-"31-40bp"
intersect_40bp_2kbup_cold_pred$intersect<-"2kb up"
colnames(intersect_40bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_50bp_2kbup_cold_pred<-read.delim("intersect_50bp_2kbup_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_50bp_2kbup_cold_pred$type<-"41-50bp"
intersect_50bp_2kbup_cold_pred$intersect<-"2kb up"
colnames(intersect_50bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')

intersect_2kbdown_cold_pred<-read.delim("intersect_2kbdown_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_2kbdown_cold_pred$type<-"all deletions"
intersect_2kbdown_cold_pred$intersect<-"2kb down"
colnames(intersect_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_10bp_2kbdown_cold_pred<-read.delim("intersect_10bp_2kbdown_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_10bp_2kbdown_cold_pred$type<-"1-10bp"
intersect_10bp_2kbdown_cold_pred$intersect<-"2kb down"
colnames(intersect_10bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_20bp_2kbdown_cold_pred<-read.delim("intersect_20bp_2kbdown_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_20bp_2kbdown_cold_pred$type<-"11-20bp"
intersect_20bp_2kbdown_cold_pred$intersect<-"2kb down"
colnames(intersect_20bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_30bp_2kbdown_cold_pred<-read.delim("intersect_30bp_2kbdown_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_30bp_2kbdown_cold_pred$type<-"21-30bp"
intersect_30bp_2kbdown_cold_pred$intersect<-"2kb down"
colnames(intersect_30bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_40bp_2kbdown_cold_pred<-read.delim("intersect_40bp_2kbdown_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_40bp_2kbdown_cold_pred$type<-"31-40bp"
intersect_40bp_2kbdown_cold_pred$intersect<-"2kb down"
colnames(intersect_40bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
intersect_50bp_2kbdown_cold_pred<-read.delim("intersect_50bp_2kbdown_cold_clustLarge_pred.txt",sep='\t', header = FALSE)
intersect_50bp_2kbdown_cold_pred$type<-"41-50bp"
intersect_50bp_2kbdown_cold_pred$intersect<-"2kb down"
colnames(intersect_50bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}
##STATISTIC 1: INDEL OVERLAP
#Calculate percentage of Indels that overlap with COs
indels_NAMv4<-read.delim(file="indelsLarge.bed")
indels_10bp<-read.delim(file="indelsLarge_10kb.bed")
indels_20bp<-read.delim(file="indelsLarge_20kb.bed")
indels_30bp<-read.delim(file="indelsLarge_30kb.bed")
indels_40bp<-read.delim(file="indelsLarge_40kb.bed")
indels_50bp<-read.delim(file="indelsLarge_50kb.bed")
library(dplyr)
percentage_overlap<-function(intersect,indels,type){
  intersect<-intersect %>% distinct()
  intersect<-intersect %>% distinct(A_start, .keep_all = TRUE)
  percentage_per_chrom<-c()
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
library(matrixStats)
# stat1<-as.data.frame(cbind( predicted=c(intersect_overlap,intersect_2kbup_overlap,intersect_2kbdown_overlap), cold=c(intersect_cold_overlap,intersect_2kbup_cold_overlap,intersect_2kbdown_cold_overlap) ))
# t.test(stat1$predicted , stat1$cold, paired=TRUE, var.equal=TRUE)
# model<-lm(predicted~cold,data=stat1)
# plot(model)

#T-TESTS
{
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
write.table(stat1,file="stat1_1_indel_predCO.txt",row.names = T,quote=FALSE, sep = '\t')

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
write.table(stat1_2,file="stat1_2_indel_predCO.txt",row.names = T,quote=FALSE, sep = '\t')

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
write.table(stat1_3,file="stat1_3_indel_predCO.txt",row.names = T,quote=FALSE, sep = '\t')

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
write.table(stat1_4,file="stat1_4_indel_predCO.txt",row.names = T,quote=FALSE, sep = '\t')
}

# stat1_2<-as.data.frame(cbind( predicted=c(intersect_10bp_overlap,intersect_20bp_overlap,intersect_30bp_overlap,intersect_40bp_overlap,intersect_50bp_overlap), cold =c(intersect_10bp_cold_overlap,intersect_20bp_cold_overlap,intersect_30bp_cold_overlap,intersect_40bp_cold_overlap,intersect_50bp_cold_overlap)))
# t.test(stat1_2$predicted , stat1_2$cold, paired=TRUE, var.equal=TRUE)
# model<-lm(predicted~cold,data=stat1_2)
# plot(model)

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
CO_pred_2kb_up<-read.delim(file="CO_pred_2kb_up.bed")
CO_pred_2kb_down<-read.delim(file="CO_pred_2kb_down.bed")
CO_pred_cold_2kb_up<-read.delim(file="CO_pred_cold_2kb_up.bed")
CO_pred_cold_2kb_down<-read.delim(file="CO_pred_cold_2kb_down.bed")

intersect_overlap<-percentage_overlap(intersect_pred,CO_predicted)
intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred,CO_predicted)
intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred,CO_predicted)
intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred,CO_predicted)
intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred,CO_predicted)
intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred,CO_predicted)

intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred,CO_pred_2kb_up)
intersect_10bp_2kbup_overlap<-percentage_overlap(intersect_10bp_2kbup_pred,CO_pred_2kb_up)
intersect_20bp_2kbup_overlap<-percentage_overlap(intersect_20bp_2kbup_pred,CO_pred_2kb_up)
intersect_30bp_2kbup_overlap<-percentage_overlap(intersect_30bp_2kbup_pred,CO_pred_2kb_up)
intersect_40bp_2kbup_overlap<-percentage_overlap(intersect_40bp_2kbup_pred,CO_pred_2kb_up)
intersect_50bp_2kbup_overlap<-percentage_overlap(intersect_50bp_2kbup_pred,CO_pred_2kb_up)

intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred,CO_pred_2kb_down)
intersect_10bp_2kbdown_overlap<-percentage_overlap(intersect_10bp_2kbdown_pred,CO_pred_2kb_down)
intersect_20bp_2kbdown_overlap<-percentage_overlap(intersect_20bp_2kbdown_pred,CO_pred_2kb_down)
intersect_30bp_2kbdown_overlap<-percentage_overlap(intersect_30bp_2kbdown_pred,CO_pred_2kb_down)
intersect_40bp_2kbdown_overlap<-percentage_overlap(intersect_40bp_2kbdown_pred,CO_pred_2kb_down)
intersect_50bp_2kbdown_overlap<-percentage_overlap(intersect_50bp_2kbdown_pred,CO_pred_2kb_down)

intersect_cold_overlap<-percentage_overlap(intersect_cold_pred,CO_pred_coldspots)
intersect_10bp_cold_overlap<-percentage_overlap(intersect_10bp_cold_pred,CO_pred_coldspots)
intersect_20bp_cold_overlap<-percentage_overlap(intersect_20bp_cold_pred,CO_pred_coldspots)
intersect_30bp_cold_overlap<-percentage_overlap(intersect_30bp_cold_pred,CO_pred_coldspots)
intersect_40bp_cold_overlap<-percentage_overlap(intersect_40bp_cold_pred,CO_pred_coldspots)
intersect_50bp_cold_overlap<-percentage_overlap(intersect_50bp_cold_pred,CO_pred_coldspots)

intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred,CO_pred_cold_2kb_up)
intersect_10bp_2kbup_cold_overlap<-percentage_overlap(intersect_10bp_2kbup_cold_pred,CO_pred_cold_2kb_up)
intersect_20bp_2kbup_cold_overlap<-percentage_overlap(intersect_20bp_2kbup_cold_pred,CO_pred_cold_2kb_up)
intersect_30bp_2kbup_cold_overlap<-percentage_overlap(intersect_30bp_2kbup_cold_pred,CO_pred_cold_2kb_up)
intersect_40bp_2kbup_cold_overlap<-percentage_overlap(intersect_40bp_2kbup_cold_pred,CO_pred_cold_2kb_up)
intersect_50bp_2kbup_cold_overlap<-percentage_overlap(intersect_50bp_2kbup_cold_pred,CO_pred_cold_2kb_up)

intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred,CO_pred_cold_2kb_down)
intersect_10bp_2kbdown_cold_overlap<-percentage_overlap(intersect_10bp_2kbdown_cold_pred,CO_pred_cold_2kb_down)
intersect_20bp_2kbdown_cold_overlap<-percentage_overlap(intersect_20bp_2kbdown_cold_pred,CO_pred_cold_2kb_down)
intersect_30bp_2kbdown_cold_overlap<-percentage_overlap(intersect_30bp_2kbdown_cold_pred,CO_pred_cold_2kb_down)
intersect_40bp_2kbdown_cold_overlap<-percentage_overlap(intersect_40bp_2kbdown_cold_pred,CO_pred_cold_2kb_down)
intersect_50bp_2kbdown_cold_overlap<-percentage_overlap(intersect_50bp_2kbdown_cold_pred,CO_pred_cold_2kb_down)
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

#T-TEST:
{
  library(matrixStats)
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
  write.table(stat1,file="stat3_1_indel_predCO.txt",row.names = T,quote=FALSE, sep = '\t')
  
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
  write.table(stat1_2,file="stat3_2_indel_predCO.txt",row.names = T,quote=FALSE, sep = '\t')
  
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
  write.table(stat1_3,file="stat3_3_indel_predCO.txt",row.names = T,quote=FALSE, sep = '\t')
  
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
  write.table(stat1_4,file="stat3_4_indel_predCO.txt",row.names = T,quote=FALSE, sep = '\t')
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
  full_data1$posthoc<-c("Ab","ab","AB","aB","AB","aB","AB","aB","AB","aB")
  #plot 1: indel density 
  perc<-ggplot(full_data1, aes(fill=color,x=levels,y=value)) +
    geom_bar(stat="identity",position="dodge") +
    scale_fill_manual(labels = c("Predicted CO",
                                 "Negative Control Region"),
                      values = c("Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("Large Indel Density on Predicted COs") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(),width=.2)+
    geom_text(aes(label=posthoc, y=value+sd), position=position_dodge(.5), vjust=-.5) +
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(0,max(full_data1$value + full_data1$sd))) +
    annotate("text", c(1.5,3.5,5.5,7.5,9.5), y = -.00004, label = c("500bp-10 kb","10-20 kb","20-30 kb","30-40 kb","40-50 kb")) +
    annotate("text", x = 5.5, y = -.00008,
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
  full_data2$posthoc<-c("AB","AB","AB","aB","aB","ab")
  perc2<-ggplot(full_data2, aes(fill=color,x=levels,y=value)) +
    geom_bar(stat="identity",position="dodge", width = .5) +
    scale_fill_manual(labels = c("2kb downstream", "2kb downstream",
                                 "2kb upstream","2kb upstream",
                                 "Predicted CO",
                                 "Negative Control"),
                      values = c("2kb downstream hot"= "darkgoldenrod1",
                                 "2kb downstream cold" = "deepskyblue",
                                 "2kb upstream hot" = "darkorange1",
                                 "2kb upstream cold" = "dodgerblue2",
                                 "Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("Large Indel Density on Predicted COs") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(.5),width=.2)+
    geom_text(aes(label=posthoc, y=value+sd), position=position_dodge(.5), vjust=-.5) +
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(0,max(full_data2$value + full_data2$sd))) +
    annotate("text", x = 1:2, y = -.0001,
             label = rep(c("Negative Control Region", "Predicted CO Region"), 1)) +
    annotate("text", x = 1.5, y = -.0002,
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
test<-cbind(test,"Mean Indel / Pred CO Overlap"=c(mean(data[which(data$type == "all insertions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Indel / Random Site Overlap"=c(mean(data_cold[which(data_cold$type == "all insertions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test<-cbind(test,"Mean Predicted CO Overlap"=c(mean(data[which(data$type == "all insertions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Predicted Random Site Overlap"=c(mean(data_cold[which(data_cold$type == "all insertions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test<-cbind(test,"Mean Indel Density on Predicted CO (indels/bp)"=c(mean(data[which(data$type == "all insertions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Indel Density on Predicted Random Site (indels/bp)"=c(mean(data_cold[which(data_cold$type == "all insertions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test_2<-as.data.frame(test)[-c(1)]
test_2 <- test_2[, c(5,6,1,2,3,4)]
write.table(test_2,file="pred_COs_aL_output.txt",row.names = T,quote=FALSE, sep = '\t')


