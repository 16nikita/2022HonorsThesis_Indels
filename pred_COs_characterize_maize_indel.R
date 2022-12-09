library(dplyr)
#QUICK RELOAD
CO_predicted<-read.delim("CO_predicted_3.bed",header=FALSE)
colnames(CO_predicted)<-c('chrom','start','end')
CO_pred_cold_2kb_up<-read.delim("CO_pred_cold_2kb_up_2.bed")
CO_pred_cold_2kb_down<-read.delim("CO_pred_cold_2kb_down_2.bed")
CO_pred_coldspots<-read.delim("CO_pred_coldspots_2.bed")
CO_pred_2kb_up<-read.delim("CO_pred_2kb_up_2.bed")
CO_pred_2kb_down<-read.delim("CO_pred_2kb_down_2.bed")
chrom_feat<-read.delim("chrom_feat.bed",sep='\t', header = TRUE)

#Where are hotspots:
CO_pred_loc<-read.delim("where_are_CO_pred.txt",sep='\t', header = FALSE)
colnames(CO_pred_loc) <- c('chrom', 'start', 'end','B_chrom', 'B_start', 'B_end','type','overlap in bp')
CO_pred_cold_loc<-read.delim("where_are_CO_pred_cold.txt",sep='\t', header = FALSE)
colnames(CO_pred_cold_loc) <- c('chrom', 'start', 'end','A_midpoint','A_length','B_chrom', 'B_start', 'B_end','type','overlap in bp')
gene_CO_pred<-CO_pred_loc[which(CO_pred_loc$type=='gene'),]
TSS_CO_pred<-CO_pred_loc[which(CO_pred_loc$type=='TSS'),]
TTS_CO_pred<-CO_pred_loc[which(CO_pred_loc$type=='TTS'),]
gene_CO_pred_cold<-CO_pred_cold_loc[which(CO_pred_cold_loc$type=='gene'),]
TSS_CO_pred_cold<-CO_pred_cold_loc[which(CO_pred_cold_loc$type=='TSS'),]
TTS_CO_pred_cold<-CO_pred_cold_loc[which(CO_pred_cold_loc$type=='TTS'),]
write.table(gene_CO_pred[c(1:3)],file="gene_CO_pred.bed",row.names = F,quote=FALSE, sep = '\t')
write.table(TSS_CO_pred[c(1:3)],file="TSS_CO_pred.bed",row.names = F,quote=FALSE, sep = '\t')
write.table(TTS_CO_pred[c(1:3)],file="TTS_CO_pred.bed",row.names = F,quote=FALSE, sep = '\t')
write.table(gene_CO_pred_cold[c(1:3)],file="gene_CO_pred_cold.bed",row.names = F,quote=FALSE, sep = '\t')
write.table(TSS_CO_pred_cold[c(1:3)],file="TSS_CO_pred_cold.bed",row.names = F,quote=FALSE, sep = '\t')
write.table(TTS_CO_pred_cold[c(1:3)],file="TTS_CO_pred_cold.bed",row.names = F,quote=FALSE, sep = '\t')

##READ IN INTERSECT DATA:
# INTERSECT WITH PREDICTED CO REGIONS
{
  intersect_pred<-read.delim("intersect_gene_CO_pred.txt",sep='\t', header = FALSE)
  intersect_pred$type<-"all deletions"
  intersect_pred$intersect<-"gene"
  colnames(intersect_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_pred<-read.delim("intersect_10bp_gene_CO_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_pred$type<-"1-10bp"
  intersect_10bp_pred$intersect<-"gene"
  colnames(intersect_10bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_pred<-read.delim("intersect_20bp_gene_CO_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_pred$type<-"11-20bp"
  intersect_20bp_pred$intersect<-"gene"
  colnames(intersect_20bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_pred<-read.delim("intersect_30bp_gene_CO_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_pred$type<-"21-30bp"
  intersect_30bp_pred$intersect<-"gene"
  colnames(intersect_30bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_pred<-read.delim("intersect_40bp_gene_CO_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_pred$type<-"31-40bp"
  intersect_40bp_pred$intersect<-"gene"
  colnames(intersect_40bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_pred<-read.delim("intersect_50bp_gene_CO_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_pred$type<-"41-50bp"
  intersect_50bp_pred$intersect<-"gene"
  colnames(intersect_50bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_pred<-read.delim("intersect_TSS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_2kbup_pred$type<-"all deletions"
  intersect_2kbup_pred$intersect<-"TSS"
  colnames(intersect_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_pred<-read.delim("intersect_10bp_TSS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_pred$type<-"1-10bp"
  intersect_10bp_2kbup_pred$intersect<-"TSS"
  colnames(intersect_10bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_pred<-read.delim("intersect_20bp_TSS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_pred$type<-"11-20bp"
  intersect_20bp_2kbup_pred$intersect<-"TSS"
  colnames(intersect_20bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_pred<-read.delim("intersect_30bp_TSS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_pred$type<-"21-30bp"
  intersect_30bp_2kbup_pred$intersect<-"TSS"
  colnames(intersect_30bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_pred<-read.delim("intersect_40bp_TSS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_pred$type<-"31-40bp"
  intersect_40bp_2kbup_pred$intersect<-"TSS"
  colnames(intersect_40bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_pred<-read.delim("intersect_50bp_TSS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_pred$type<-"41-50bp"
  intersect_50bp_2kbup_pred$intersect<-"TSS"
  colnames(intersect_50bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_pred<-read.delim("intersect_TTS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_2kbdown_pred$type<-"all deletions"
  intersect_2kbdown_pred$intersect<-"TTS"
  colnames(intersect_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_pred<-read.delim("intersect_10bp_TTS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_pred$intersect<-"TTS"
  colnames(intersect_10bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_pred<-read.delim("intersect_20bp_TTS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_pred$intersect<-"TTS"
  colnames(intersect_20bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_pred<-read.delim("intersect_30bp_TTS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_pred$intersect<-"TTS"
  colnames(intersect_30bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_pred<-read.delim("intersect_40bp_TTS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_pred$intersect<-"TTS"
  colnames(intersect_40bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_pred<-read.delim("intersect_40bp_TTS_CO_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_pred$intersect<-"TTS"
  colnames(intersect_50bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}
# INTERSECT WITH RANDOM REGIONS
{
  intersect_cold_pred<-read.delim("intersect_gene_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_cold_pred$type<-"all deletions"
  intersect_cold_pred$intersect<-"gene"
  colnames(intersect_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_cold_pred<-read.delim("intersect_10bp_gene_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_10bp_cold_pred$type<-"1-10bp"
  intersect_10bp_cold_pred$intersect<-"gene"
  colnames(intersect_10bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_cold_pred<-read.delim("intersect_20bp_gene_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_20bp_cold_pred$type<-"11-20bp"
  intersect_20bp_cold_pred$intersect<-"gene"
  colnames(intersect_20bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_cold_pred<-read.delim("intersect_30bp_gene_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_30bp_cold_pred$type<-"21-30bp"
  intersect_30bp_cold_pred$intersect<-"gene"
  colnames(intersect_30bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_cold_pred<-read.delim("intersect_40bp_gene_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_40bp_cold_pred$type<-"31-40bp"
  intersect_40bp_cold_pred$intersect<-"gene"
  colnames(intersect_40bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_cold_pred<-read.delim("intersect_50bp_gene_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_50bp_cold_pred$type<-"41-50bp"
  intersect_50bp_cold_pred$intersect<-"gene"
  colnames(intersect_50bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_cold_pred<-read.delim("intersect_TSS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_2kbup_cold_pred$type<-"all deletions"
  intersect_2kbup_cold_pred$intersect<-"TSS"
  colnames(intersect_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_cold_pred<-read.delim("intersect_10bp_TSS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbup_cold_pred$intersect<-"TSS"
  colnames(intersect_10bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_cold_pred<-read.delim("intersect_20bp_TSS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbup_cold_pred$intersect<-"TSS"
  colnames(intersect_20bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_cold_pred<-read.delim("intersect_30bp_TSS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbup_cold_pred$intersect<-"TSS"
  colnames(intersect_30bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_cold_pred<-read.delim("intersect_40bp_TSS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbup_cold_pred$intersect<-"TSS"
  colnames(intersect_40bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_cold_pred<-read.delim("intersect_50bp_TSS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbup_cold_pred$intersect<-"TSS"
  colnames(intersect_50bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_cold_pred<-read.delim("intersect_TTS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_2kbdown_cold_pred$type<-"all deletions"
  intersect_2kbdown_cold_pred$intersect<-"TTS"
  colnames(intersect_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_cold_pred<-read.delim("intersect_10bp_TTS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_cold_pred$intersect<-"TTS"
  colnames(intersect_10bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_cold_pred<-read.delim("intersect_20bp_TTS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_cold_pred$intersect<-"TTS"
  colnames(intersect_20bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_cold_pred<-read.delim("intersect_30bp_TTS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_cold_pred$intersect<-"TTS"
  colnames(intersect_30bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_cold_pred<-read.delim("intersect_40bp_TTS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_cold_pred$intersect<-"TTS"
  colnames(intersect_40bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_cold_pred<-read.delim("intersect_40bp_TTS_CO_pred_cold.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_cold_pred$intersect<-"TTS"
  colnames(intersect_50bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}

#calculate distance of closest gene from predicted COs:
genes<-chrom_feat[which(chrom_feat$type=='gene'),]
find_closest<-function(genes,CO){
  CO$closest<-NA
  for(i in 1:nrow(CO)){
    closest<-genes[which.min(abs(CO$start[i]-genes$start)),]
    CO$closest[i]<-abs(CO$start[i]-closest$start)
  }
  return(CO)
}
CO_predicted_2<-find_closest(genes,CO_predicted)
write.table(CO_predicted_2,file="Pred_COs_distance.bed",row.names = F,quote=FALSE, sep = '\t')

#ASSOCIATION BETWEEN DISTANCE TO GENE & INDEL DENSITY?:
DISTANCE_pred<-read.delim("distance_CO_pred.txt",sep='\t', header = FALSE)
colnames(DISTANCE_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','distance to gene','overlap in bp')
percentage_overlap<-function(intersect){
  intersect<-intersect %>% distinct()
  duplicated<-as.data.frame(table(intersect$'B_start'))
  intersect$length<-intersect$'B_end'-intersect$'B_start'
  intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
  intersect$density<-duplicated$Freq
  intersect$percent<- (intersect$'density'/intersect$length)
  return(intersect)
}
distance_overlap<-percentage_overlap(DISTANCE_pred)
spearmen<-cor.test(distance_overlap$percent, distance_overlap$'distance to gene',  method = "spearman", alternative = "greater")

##STATISTIC 1: INDEL OVERLAP
#Calculate percentage of Indels that overlap with COs
indels_NAMv4<-read.delim(file="indels_all.bed")
indels_10bp<-read.delim(file="indels_10bp_all.bed")
indels_20bp<-read.delim(file="indels_20bp_all.bed")
indels_30bp<-read.delim(file="indels_30bp_all.bed")
indels_40bp<-read.delim(file="indels_40bp_all.bed")
indels_50bp<-read.delim(file="indels_50bp_all.bed")
percentage_overlap<-function(intersect,indels,type){
  intersect<-intersect %>% distinct()
  intersect<-intersect %>% distinct(A_start, .keep_all = TRUE)
  indels<-indels %>% distinct()
  
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
write.table(stat1,file="stat1_1_indel_hotspot.txt",row.names = T,quote=FALSE, sep = '\t')

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
write.table(stat1_2,file="stat1_2_indel_hotspot.txt",row.names = T,quote=FALSE, sep = '\t')

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
write.table(stat1_3,file="stat1_3_indel_hotspot.txt",row.names = T,quote=FALSE, sep = '\t')

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
write.table(stat1_4,file="stat1_4_indel_hotspot.txt",row.names = T,quote=FALSE, sep = '\t')
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
intersect_overlap<-percentage_overlap(intersect_pred,hotspots_2)
intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred,hotspots_2)
intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred,hotspots_2)
intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred,hotspots_2)
intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred,hotspots_2)
intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred,hotspots_2)

intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred,hotspots_2_2kbup)
intersect_10bp_2kbup_overlap<-percentage_overlap(intersect_10bp_2kbup_pred,hotspots_2_2kbup)
intersect_20bp_2kbup_overlap<-percentage_overlap(intersect_20bp_2kbup_pred,hotspots_2_2kbup)
intersect_30bp_2kbup_overlap<-percentage_overlap(intersect_30bp_2kbup_pred,hotspots_2_2kbup)
intersect_40bp_2kbup_overlap<-percentage_overlap(intersect_40bp_2kbup_pred,hotspots_2_2kbup)
intersect_50bp_2kbup_overlap<-percentage_overlap(intersect_50bp_2kbup_pred,hotspots_2_2kbup)

intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred,hotspots_2_2kbdown)
intersect_10bp_2kbdown_overlap<-percentage_overlap(intersect_10bp_2kbdown_pred,hotspots_2_2kbdown)
intersect_20bp_2kbdown_overlap<-percentage_overlap(intersect_20bp_2kbdown_pred,hotspots_2_2kbdown)
intersect_30bp_2kbdown_overlap<-percentage_overlap(intersect_30bp_2kbdown_pred,hotspots_2_2kbdown)
intersect_40bp_2kbdown_overlap<-percentage_overlap(intersect_40bp_2kbdown_pred,hotspots_2_2kbdown)
intersect_50bp_2kbdown_overlap<-percentage_overlap(intersect_50bp_2kbdown_pred,hotspots_2_2kbdown)

intersect_cold_overlap<-percentage_overlap(intersect_cold_pred,hotspot_random)
intersect_10bp_cold_overlap<-percentage_overlap(intersect_10bp_cold_pred,hotspot_random)
intersect_20bp_cold_overlap<-percentage_overlap(intersect_20bp_cold_pred,hotspot_random)
intersect_30bp_cold_overlap<-percentage_overlap(intersect_30bp_cold_pred,hotspot_random)
intersect_40bp_cold_overlap<-percentage_overlap(intersect_40bp_cold_pred,hotspot_random)
intersect_50bp_cold_overlap<-percentage_overlap(intersect_50bp_cold_pred,hotspot_random)

intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred,hotspot_random_2kb_up)
intersect_10bp_2kbup_cold_overlap<-percentage_overlap(intersect_10bp_2kbup_cold_pred,hotspot_random_2kb_up)
intersect_20bp_2kbup_cold_overlap<-percentage_overlap(intersect_20bp_2kbup_cold_pred,hotspot_random_2kb_up)
intersect_30bp_2kbup_cold_overlap<-percentage_overlap(intersect_30bp_2kbup_cold_pred,hotspot_random_2kb_up)
intersect_40bp_2kbup_cold_overlap<-percentage_overlap(intersect_40bp_2kbup_cold_pred,hotspot_random_2kb_up)
intersect_50bp_2kbup_cold_overlap<-percentage_overlap(intersect_50bp_2kbup_cold_pred,hotspot_random_2kb_up)

intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred,hotspot_random_2kb_down)
intersect_10bp_2kbdown_cold_overlap<-percentage_overlap(intersect_10bp_2kbdown_cold_pred,hotspot_random_2kb_down)
intersect_20bp_2kbdown_cold_overlap<-percentage_overlap(intersect_20bp_2kbdown_cold_pred,hotspot_random_2kb_down)
intersect_30bp_2kbdown_cold_overlap<-percentage_overlap(intersect_30bp_2kbdown_cold_pred,hotspot_random_2kb_down)
intersect_40bp_2kbdown_cold_overlap<-percentage_overlap(intersect_40bp_2kbdown_cold_pred,hotspot_random_2kb_down)
intersect_50bp_2kbdown_cold_overlap<-percentage_overlap(intersect_50bp_2kbdown_cold_pred,hotspot_random_2kb_down)
}

##STATISTIC 3: mean indel density
library(dplyr)
#mean indel density on each CO
percentage_overlap<-function(intersect){
  intersect<-intersect %>% distinct(A_start, .keep_all = TRUE)
  duplicated<-as.data.frame(table(intersect$'B_start'))
  intersect$length<-intersect$'B_end'-intersect$'B_start'
  intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
  intersect$density<-duplicated$Freq
  intersect$percent<- (intersect$'density'/intersect$length)
  percentage_per_chrom<-c()
  for(i in 1:10){
    intersect_temp<-intersect[which(intersect$`A_chrom` == i),]
    percentage_per_chrom[i]<- mean(intersect_temp$percent)
    if(is.na(percentage_per_chrom[i])){
      percentage_per_chrom[i]<-0
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

}

#BARGRAPHS- for characterizing each dataset intersect w gene/tss/tts
{
  #ANOVA test for significance
  {
    intersect<-c(rep(c("Gene"),10),rep(c("TSS"),10),rep(c("TTS"),10))
    value <- c((intersect_overlap),(intersect_2kbup_overlap),(intersect_2kbdown_overlap))
    data1 <- data.frame(y=value,group=factor(intersect))
    fit = lm(y ~ group, data1)
    #p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
    shapiro.test(value)
    
    #kruskal-wallis- no normal distribution, equal variance
    
    library(devtools)
    library(qqplotr)
    gg <- ggplot(data = data1, mapping = aes(sample = value)) +
      geom_qq_band(bandType = "ks", mapping = aes(fill = "KS"), alpha = 0.5) +
      geom_qq_band(bandType = "ts", mapping = aes(fill = "TS"), alpha = 0.5) +
      geom_qq_band(bandType = "pointwise", mapping = aes(fill = "Normal"), alpha = 0.5) +
      geom_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5) +
      stat_qq_line() +
      stat_qq_point() +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      scale_fill_discrete("Bandtype")
    gg
    
    library(car)
    leveneTest(y ~ group, data = data1)
    #Levene's test indicated equal variances (F = .0094, p = .9907)
    one.way <- aov(fit)
    summary(one.way)
    #ANOVA test indicated no significant diff (F=0.505, p=0.609)
  }
  
  #create dataframes:
  {
    #DATAFRAME OF CO PREDICTED
    type <- c(rep(c("all indels") ,3))
    intersect<-c(c("Gene"),c("TSS") ,c("TTS"))
    value <- c(mean(intersect_overlap),mean(intersect_2kbup_overlap),mean(intersect_2kbdown_overlap))
    se <- c(sd(intersect_overlap)/sqrt(length((intersect_overlap))),sd(intersect_2kbup_overlap)/sqrt(length((intersect_2kbup_overlap))),sd(intersect_2kbdown_overlap)/sqrt(length((intersect_2kbdown_overlap))))
    data <- data.frame(type,value,intersect,sd)
    #data$dataset<-"hot"
    
  }
  library(ggplot2)
  library(viridis)
  library(hrbrthemes)
  dev.off()
  full_data<-data
  #plot 1: indel density 
  perc<-ggplot(full_data, aes(fill=intersect,x=intersect,y=value)) +
    geom_bar(stat="identity",position="dodge",width=.5) +
    scale_fill_manual(values = c("Gene" = "sandybrown",
                                 "TSS" = "deepskyblue",
                                 "TTS" = 'palegreen'))+
    ggtitle("Small Indel Density on Predicted Crossovers in Gene/TSS/TTS") +
    geom_errorbar(aes(ymax = value + se, ymin = value - se), position = position_dodge(),width=.2)+
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    xlab("Predicted CO Feature Location")+
    ylab("Indel Density (indels/bp)")+
    #annotate("label", x = 4,y =max(full_data$value), label = "F=0.15, P=0.861")+
    #guides(fill=guide_legend("Type of intersection"))+
    coord_cartesian(ylim = c(0, max(full_data$value+se)),
                    clip = 'off')+
    theme_classic() +
    theme(plot.margin = unit(c(3, 3, 5, 1), "lines"),
          plot.title = element_text(face="bold"),
          axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          legend.position="none")
  #plot margins unit(c(top, right, bottom, left),
  perc
  
  
}

##STATISTIC 4: Looking closer at indel density
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
}

{
intersect_overlap<-intersect_overlap[order(intersect_overlap$`B_chrom`,intersect_overlap$`B_start`),]
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==1),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 1 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==2),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 2 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==3),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 3 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==4),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 4 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==5),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 5 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==6),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 6 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==7),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 7 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==8),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 8 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==9),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 9 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")
intersect_overlap_plot<-intersect_overlap[which(intersect_overlap$B_chrom==10),]
plot(intersect_overlap_plot$B_start/1000000,intersect_overlap_plot$percent,type="p",main="Chromosome 10 Hotspots & Indel Densities",xlab="Hotspot Position (Mb)",ylab="Density (indels/bp)")

#find outliers
out <- boxplot.stats(intersect_overlap$percent)$out
out_ind <- which(intersect_overlap$percent %in% c(out))
stat_outliers<- intersect_overlap[out_ind, ]

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==1),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==1),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=3, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 1 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+  expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_1.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==2),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==2),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=3, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 2 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+  expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_2.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==3),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$'B_chrom'==3),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=2, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 3 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+ expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_3.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==4),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==4),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=3, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 4 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+ expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_4.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==5),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==5),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=3, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 5 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+ expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_5.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==6),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==6),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=3, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 6 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+ expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_6.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==7),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==7),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=3, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 7 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+ expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_7.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==8),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==8),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=3, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 8 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+ expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_8.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==9),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==9),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=3, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 9 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+ expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_9.rds")

max_temp<-max(indels_NAMv4[which(indels_NAMv4$chrom==10),]$start)/1000000
intersect_overlap_plot<-stat_outliers[which(stat_outliers$B_chrom==10),]
intersect_overlap_plot<-ggplot(intersect_overlap_plot, aes(x=B_start/1000000, y=percent)) + geom_bar(stat="identity", width=2, fill="steelblue")+
  theme_minimal() + labs(title="Chromosome 10 Indel Densities at Outlier Hotspots", 
                         x="Hotspot Position (Mb)", y = "Deletion Density (Indels/bp)")+ expand_limits(x = c(0,max_temp))
saveRDS(intersect_overlap_plot,file = "outlier_hotspot_10.rds")

}

##STATISTIC 5: metaplots of indel density
library(dplyr)
#mean indel density on each CO
# percentage_overlap<-function(intersect){
#   intersect<-intersect %>% distinct()
#   duplicated<-as.data.frame(table(intersect$'B_start'))
#   intersect$length<-intersect$'B_end'-intersect$'B_start'
#   intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
#   intersect$density<-duplicated$Freq
#   intersect$percent<- (intersect$'density'/intersect$length)
#   return(intersect)
# }
percentage_overlap<-function(intersect){
  intersect<-intersect %>% distinct()
  #intersect<-intersect %>% distinct(A_start,keep.all=TRUE)
  intersect<-intersect[order(intersect$B_chrom,intersect$B_start),]
  duplicated<-as.data.frame(table(intersect$'B_start'))
  colnames(duplicated)<-c("B_start","Freq")
  merged <- merge(intersect, duplicated, by = "B_start", sort = F, all.x = T)
  merged$length<-merged$'B_end'-merged$'B_start'
  merged$percent<-merged$Freq/merged$length
  intersect<-merged
  return(intersect)
}
{
  metaplot<-function(all, all_cold){
    intersect<-all[ which(all$intersection=="intersect"),]  
    intersect_cold<-all_cold[ which(all_cold$intersection=="intersect"),]  
    
    M=matrix(((intersect$B_start+intersect$B_end)/2)-intersect$A_start,ncol=1)
    intersect=cbind(intersect,M)
    intersect=intersect[intersect$M<5000,]
    tapply(intersect$percent, cut(intersect$M, seq(-5000,5000, by=10)), mean)->mean_intersect
    test<-as.data.frame(mean_intersect)
    test$mean_intersect[is.na(test$mean_intersect)]<-0
    test<-unlist(test)
    max<-max(test)
    
    M=matrix(((intersect_cold$B_start+intersect_cold$B_end)/2)- intersect_cold$A_start,ncol=1)
    intersect_cold=cbind(intersect_cold,M)
    intersect_cold=intersect_cold[intersect_cold$M<10000,]
    tapply(intersect_cold$percent, cut(intersect_cold$M, seq(-5000,5000, by=10)), mean)->mean_intersect_cold
    test2<-as.data.frame(mean_intersect_cold)
    test2$mean_intersect_cold[is.na(test2$mean_intersect_cold)]<-0
    test2<-unlist(test2)
    plot(smooth.spline(seq(-5000,4999,10),test,spar=0.5),col="red",lwd=1,type="l",ylim=c(0,max),main="Indel Densities at Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
    lines(smooth.spline(seq(-5000,4999,10),test2,spar=0.5),col="blue",lwd=1,type="l",main="Indel Densities at Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
    # saveRDS(intersect_plot, file="intersect_meta.rds")
    
  }
  metaplot_2<-function(all, all_cold){
    #Get max value for plot
    intersect<-all[ which(all$intersection=="intersect"),]  
    intersect_cold<-all_cold[ which(all_cold$intersection=="intersect"),]  
    
    M=matrix(((intersect$B_start+intersect$B_end)/2)-intersect$A_start,ncol=1)
    intersect=cbind(intersect,M)
    intersect=intersect[intersect$M<5000,]
    tapply(intersect$percent, cut(intersect$M, seq(-5000,5000, by=10)), mean)->mean_intersect
    test<-as.data.frame(mean_intersect)
    test$mean_intersect[is.na(test$mean_intersect)]<-0
    test<-unlist(test)
    max<-max(test)
    
    down2kb<-all[ which(all$intersection=="2kb down"),]  
    up2kb<-all[ which(all$intersection=="2kb up"),]
    down2kb_cold<-all_cold[ which(all_cold$intersection=="2kb down"),]  
    up2kb_cold<-all_cold[ which(all_cold$intersection=="2kb up"),]
    
    M=matrix((down2kb$B_start+down2kb$B_end)/2- down2kb$A_start,ncol=1)
    down2kb=cbind(down2kb,M)
    down2kb=down2kb[down2kb$M<20000,]
    tapply(down2kb$percent, cut(down2kb$M, seq(-10000, 10000, by=100)), mean)->mean_down2kb
    test<-as.data.frame(mean_down2kb)
    test$mean_down2kb[is.na(test$mean_down2kb)]<-0
    test<-unlist(test)
    
    M=matrix((down2kb_cold$B_start+down2kb_cold$B_end)/2- down2kb_cold$A_start,ncol=1)
    down2kb_cold=cbind(down2kb_cold,M)
    down2kb_cold=down2kb_cold[down2kb_cold$M<20000,]
    tapply(down2kb_cold$percent, cut(down2kb_cold$M, seq(-10000, 10000, by=100)), mean)->mean_down2kb_cold
    test2<-as.data.frame(mean_down2kb_cold)
    test2$mean_down2kb_cold[is.na(test2$mean_down2kb_cold)]<-0
    test2<-unlist(test2)
    plot(smooth.spline(seq(-10000,9999,100),test,spar=0.5),col="red",lwd=1,type="l",ylim=c(0,max),main="Indel Densities 20kb down from Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
    lines(smooth.spline(seq(-10000,9999,100),test2,spar=0.5),col="blue",lwd=1,type="l",main="Indel Densities 20kb down from Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
    
    M=matrix((up2kb$B_start+up2kb$B_end)/2- up2kb$A_start,ncol=1)
    up2kb=cbind(up2kb,M)
    up2kb=up2kb[up2kb$M<20000,]
    tapply(up2kb$percent, cut(up2kb$M, seq(-10000, 10000, by=100)), mean)->mean_up2kb
    test<-as.data.frame(mean_up2kb)
    test$mean_up2kb[is.na(test$mean_up2kb)]<-0
    test<-unlist(test)
    
    M=matrix((up2kb_cold$B_start+up2kb_cold$B_end)/2- up2kb_cold$A_start,ncol=1)
    up2kb_cold=cbind(up2kb_cold,M)
    up2kb_cold=up2kb_cold[up2kb_cold$M<20000,]
    tapply(up2kb_cold$percent, cut(up2kb_cold$M, seq(-10000, 10000, by=100)), mean)->mean_up2kb_cold
    test2<-as.data.frame(mean_up2kb_cold)
    test2$mean_up2kb_cold[is.na(test2$mean_up2kb_cold)]<-0
    test2<-unlist(test2)
    plot(smooth.spline(seq(-10000,9999,100),test,spar=0.5),col="red",lwd=1,type="l",ylim=c(0,max),main="Indel Densities 20kb up from Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
    lines(smooth.spline(seq(-10000,9999,100),test2,spar=0.5),col="blue",lwd=1,type="l",main="Indel Densities 20kb up from Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
  
  }
  intersect_overlap<-percentage_overlap(intersect_pred)
  intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred)
  intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred)
  all<-rbind(intersect_overlap,intersect_2kbup_overlap,intersect_2kbdown_overlap)
  all<-all[order(all$B_chrom,all$B_start),]
  all<-all[c(1:6,9,12)]
  intersect_cold_overlap<-percentage_overlap(intersect_cold_pred)
  intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred)
  intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred)
  all_cold<-rbind(intersect_cold_overlap,intersect_2kbup_cold_overlap,intersect_2kbdown_cold_overlap)
  all_cold<-all_cold[order(all_cold$B_chrom,all_cold$B_start),]
  all_cold<-all_cold[c(1:6,9,12)]
  metaplot(all, all_cold)
  metaplot_2(all, all_cold)
  
}