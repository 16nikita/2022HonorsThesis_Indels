library(dplyr)
#QUICK RELOAD
genes<-read.delim("genes_maize.bed",sep='\t', header = TRUE)
genes_2kb_up<-read.delim("gene_2kb_up.bed",sep='\t', header = TRUE)
genes_2kb_down<-read.delim("gene_2kb_down.bed",sep='\t', header = TRUE)
genes_cold<-read.delim("gene_coldspots.bed",sep='\t', header = TRUE)
genes_cold_2kb_up<-read.delim("gene_cold_2kb_up.bed",sep='\t', header = TRUE)
genes_cold_2kb_down<-read.delim("gene_cold_2kb_down.bed",sep='\t', header = TRUE)

#Where are genes:
genes<-read.delim("chrom_feat.bed",sep='\t', header = TRUE)
genes<-genes[which(genes$type=="gene"),]
DSB<-genes
#write.table(genes,"genes_maize.bed", row.names = F,quote=FALSE, sep = '\t')
genes_2kb_up <- CO_up(data.frame(matrix(ncol = 3,nrow = nrow(genes))),genes)
genes_2kb_down <-CO_down(data.frame(matrix(ncol = 3,nrow = nrow(genes))),genes)
#write.table(genes_2kb_up,"gene_2kb_up.bed", row.names = F,quote=FALSE, sep = '\t')
#write.table(genes_2kb_down,"gene_2kb_down.bed", row.names = F,quote=FALSE, sep = '\t')



#CREATE PREDICTED GENES COLD DATASET
{
  distance<-function(DSB){
    DSB$distance<-NA
    for(i in 1:nrow(DSB)-1){
      DSB$distance[i]<-DSB$start[i+1]-DSB$end[i]
    }
    DSB <- DSB[order(DSB$distance, decreasing = TRUE), ]
    DSB <- data.table(DSB)
    DSB <- DSB[ , head(.SD, 10)]
    DSB <- DSB[order(DSB$chrom, DSB$start), ]
    return(DSB)
  }
  library(data.table)
  DSB<-DSB[order(DSB$chrom,DSB$start),]
  DSB_1<-DSB[which(DSB$chrom==1),]
  DSB_2<-DSB[which(DSB$chrom==2),]
  DSB_3<-DSB[which(DSB$chrom==3),]
  DSB_4<-DSB[which(DSB$chrom==4),]
  DSB_5<-DSB[which(DSB$chrom==5),]
  DSB_6<-DSB[which(DSB$chrom==6),]
  DSB_7<-DSB[which(DSB$chrom==7),]
  DSB_8<-DSB[which(DSB$chrom==8),]
  DSB_9<-DSB[which(DSB$chrom==9),]
  DSB_10<-DSB[which(DSB$chrom==10),]
  DSB_cold_1<-distance(DSB_1)
  DSB_cold_2<-distance(DSB_2)
  DSB_cold_3<-distance(DSB_3)
  DSB_cold_4<-distance(DSB_4)
  DSB_cold_5<-distance(DSB_5)
  DSB_cold_6<-distance(DSB_6)
  DSB_cold_7<-distance(DSB_7)
  DSB_cold_8<-distance(DSB_8)
  DSB_cold_9<-distance(DSB_9)
  DSB_cold_10<-distance(DSB_10)
  DSB_cold_regions<-rbind(DSB_cold_1,DSB_cold_2,DSB_cold_3,DSB_cold_4,DSB_cold_5,DSB_cold_6,DSB_cold_7,DSB_cold_8,DSB_cold_9,DSB_cold_10)
  write.table(DSB_cold_regions,"gene_cold_regions.bed", row.names = F,quote=FALSE, sep = '\t')
  
  CO_up<-function(up,hotspots){
    colnames(up)<- c('chrom','start','end')
    up$chrom<-hotspots$chrom
    up$start<-hotspots$start-2000
    up$end<-hotspots$start
    up<-up[which(up$start>0),]
    return(up)
  }
  CO_down<-function(down,hotspots){
    colnames(down)<- c('chrom','start','end')
    down$chrom<-hotspots$chrom
    down$start<-hotspots$end
    down$end<-hotspots$end +2000
    return(down)
  }
  CO_coldspots<-read.delim("gene_cold_regions.bed",sep='\t', header = TRUE)
  library(dplyr)
  library(reshape)
  create_coldspot<-function(rows,sample,chrom,size){
    set.seed(200)
    c<-data.frame(matrix(nrow=rows,ncol=nrow(sample)))
    for(j in 1:ncol(c)){
      for(k in 1:nrow(c)){
        c[k,j]<-sample(sample$start[j]:sample$end[j],size=1)
      }
    }
    coldspots<-data.frame(matrix(nrow=(nrow(c)*ncol(c)),ncol=3))
    colnames(coldspots)<-c('chrom','start','end')
    coldspots$start<- unlist(c)
    if(nrow(coldspots)>rows){
      copy<-coldspots$start
      coldspots <- slice(coldspots, c(1:rows))
      coldspots$start<-sample(copy,size=rows)
    }
    coldspots<-coldspots %>% distinct(start, .keep_all = TRUE)
    coldspots<-coldspots[order(coldspots$start),]
    for(i in 1:(nrow(coldspots)-1)){
      if(coldspots$start[i]+size>=coldspots$start[i+1]){
        diff<-(coldspots$start[i]+size)-(coldspots$start[i+1])
        coldspots$start[i+1]<-coldspots$start[i+1]+diff
      }
    }
    coldspots$chrom<-chrom
    coldspots$end<-coldspots$start+size
    options(scipen=999999)
    return(coldspots)
  }
  #store row counts
  chrom1<- nrow(DSB[which(DSB$chrom == 1),])
  chrom2<- nrow(DSB[which(DSB$chrom == 2),])
  chrom3<- nrow(DSB[which(DSB$chrom == 3),])
  chrom4<- nrow(DSB[which(DSB$chrom == 4),])
  chrom5<- nrow(DSB[which(DSB$chrom == 5),])
  chrom6<- nrow(DSB[which(DSB$chrom == 6),])
  chrom7<- nrow(DSB[which(DSB$chrom == 7),])
  chrom8<- nrow(DSB[which(DSB$chrom == 8),])
  chrom9<- nrow(DSB[which(DSB$chrom == 9),])
  chrom10<- nrow(DSB[which(DSB$chrom == 10),])
  
  install.packages("valr")
  library(valr)
  
  CO_coldspots_1<-create_coldspot(chrom1,CO_coldspots[which(CO_coldspots$chrom == 1),],1,3677)
  CO_coldspots_2<-create_coldspot(chrom2,CO_coldspots[which(CO_coldspots$chrom == 2),],2,3677)
  CO_coldspots_3<-create_coldspot(chrom3,CO_coldspots[which(CO_coldspots$chrom == 3),],3,3677)
  CO_coldspots_4<-create_coldspot(chrom4,CO_coldspots[which(CO_coldspots$chrom == 4),],4,3677)
  CO_coldspots_5<-create_coldspot(chrom5,CO_coldspots[which(CO_coldspots$chrom == 5),],5,3677)
  CO_coldspots_6<-create_coldspot(chrom6,CO_coldspots[which(CO_coldspots$chrom == 6),],6,3677)
  CO_coldspots_7<-create_coldspot(chrom7,CO_coldspots[which(CO_coldspots$chrom == 7),],7,3677)
  CO_coldspots_8<-create_coldspot(chrom8,CO_coldspots[which(CO_coldspots$chrom == 8),],8,3677)
  CO_coldspots_9<-create_coldspot(chrom9,CO_coldspots[which(CO_coldspots$chrom == 9),],9,3677)
  CO_coldspots_10<-create_coldspot(chrom10,CO_coldspots[which(CO_coldspots$chrom == 10),],10,3677)
  
  #Bind to create dataframe with all random regions
  CO_random<-rbind(CO_coldspots_1,CO_coldspots_2,CO_coldspots_3,CO_coldspots_4,CO_coldspots_5,CO_coldspots_6,CO_coldspots_7,CO_coldspots_8,CO_coldspots_9,CO_coldspots_10)
  write.table(CO_random,"gene_coldspots.bed", row.names = F,quote=FALSE, sep = '\t')
  
  #Creating datasets of 2kb upstream and 2kb downstream from CO "coldspot" intervals
  CO_random_2kb_up <- CO_up(data.frame(matrix(ncol = 3,nrow = nrow(CO_random))),CO_random)
  CO_random_2kb_down <-CO_down(data.frame(matrix(ncol = 3,nrow = nrow(CO_random))),CO_random)
  write.table(CO_random_2kb_up,"gene_cold_2kb_up.bed", row.names = F,quote=FALSE, sep = '\t')
  write.table(CO_random_2kb_down,"gene_cold_2kb_down.bed", row.names = F,quote=FALSE, sep = '\t')
  
}


##READ IN INTERSECT DATA:
# INTERSECT WITH PREDICTED CO REGIONS
{
  intersect_pred<-read.delim("intersect_a_pred.txt",sep='\t', header = FALSE)
  intersect_pred$type<-"all deletions"
  intersect_pred$intersect<-"intersect"
  colnames(intersect_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_pred<-read.delim("intersect_10bp_a_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_pred$type<-"1-10bp"
  intersect_10bp_pred$intersect<-"intersect"
  colnames(intersect_10bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_pred<-read.delim("intersect_20bp_a_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_pred$type<-"11-20bp"
  intersect_20bp_pred$intersect<-"intersect"
  colnames(intersect_20bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_pred<-read.delim("intersect_30bp_a_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_pred$type<-"21-30bp"
  intersect_30bp_pred$intersect<-"intersect"
  colnames(intersect_30bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_pred<-read.delim("intersect_40bp_a_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_pred$type<-"31-40bp"
  intersect_40bp_pred$intersect<-"intersect"
  colnames(intersect_40bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_pred<-read.delim("intersect_50bp_a_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_pred$type<-"41-50bp"
  intersect_50bp_pred$intersect<-"intersect"
  colnames(intersect_50bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_pred<-read.delim("intersect_2kbup_a_pred.txt",sep='\t', header = FALSE)
  intersect_2kbup_pred$type<-"all deletions"
  intersect_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_pred<-read.delim("intersect_10bp_2kbup_a_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_pred$type<-"1-10bp"
  intersect_10bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_pred<-read.delim("intersect_20bp_2kbup_a_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_pred$type<-"11-20bp"
  intersect_20bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_pred<-read.delim("intersect_30bp_2kbup_a_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_pred$type<-"21-30bp"
  intersect_30bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_pred<-read.delim("intersect_40bp_2kbup_a_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_pred$type<-"31-40bp"
  intersect_40bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_pred<-read.delim("intersect_50bp_2kbup_a_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_pred$type<-"41-50bp"
  intersect_50bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_pred<-read.delim("intersect_2kbdown_a_pred.txt",sep='\t', header = FALSE)
  intersect_2kbdown_pred$type<-"all deletions"
  intersect_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_pred<-read.delim("intersect_10bp_2kbdown_a_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_pred<-read.delim("intersect_20bp_2kbdown_a_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_pred<-read.delim("intersect_30bp_2kbdown_a_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_pred<-read.delim("intersect_40bp_2kbdown_a_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_pred<-read.delim("intersect_50bp_2kbdown_a_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}
# INTERSECT WITH RANDOM CO DESERT REGIONS
{
  intersect_cold_pred<-read.delim("intersect_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_cold_pred$type<-"all deletions"
  intersect_cold_pred$intersect<-"intersect"
  intersect_cold_pred<-intersect_cold_pred[-c(7:8)]
  colnames(intersect_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_cold_pred<-read.delim("intersect_10bp_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_cold_pred$type<-"1-10bp"
  intersect_10bp_cold_pred$intersect<-"intersect"
  intersect_10bp_cold_pred<-intersect_10bp_cold_pred[-c(7:8)]
  colnames(intersect_10bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_cold_pred<-read.delim("intersect_20bp_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_cold_pred$type<-"11-20bp"
  intersect_20bp_cold_pred$intersect<-"intersect"
  intersect_20bp_cold_pred<-intersect_20bp_cold_pred[-c(7:8)]
  colnames(intersect_20bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_cold_pred<-read.delim("intersect_30bp_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_cold_pred$type<-"21-30bp"
  intersect_30bp_cold_pred$intersect<-"intersect"
  intersect_30bp_cold_pred<-intersect_30bp_cold_pred[-c(7:8)]
  colnames(intersect_30bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_cold_pred<-read.delim("intersect_40bp_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_cold_pred$type<-"31-40bp"
  intersect_40bp_cold_pred$intersect<-"intersect"
  intersect_40bp_cold_pred<-intersect_40bp_cold_pred[-c(7:8)]
  colnames(intersect_40bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_cold_pred<-read.delim("intersect_50bp_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_cold_pred$type<-"41-50bp"
  intersect_50bp_cold_pred$intersect<-"intersect"
  intersect_50bp_cold_pred<-intersect_50bp_cold_pred[-c(7:8)]
  colnames(intersect_50bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_cold_pred<-read.delim("intersect_2kbup_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_2kbup_cold_pred$type<-"all deletions"
  intersect_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_cold_pred<-read.delim("intersect_10bp_2kbup_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_cold_pred<-read.delim("intersect_20bp_2kbup_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_cold_pred<-read.delim("intersect_30bp_2kbup_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_cold_pred<-read.delim("intersect_40bp_2kbup_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_cold_pred<-read.delim("intersect_50bp_2kbup_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_cold_pred<-read.delim("intersect_2kbdown_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_2kbdown_cold_pred$type<-"all deletions"
  intersect_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_cold_pred<-read.delim("intersect_10bp_2kbdown_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_cold_pred<-read.delim("intersect_20bp_2kbdown_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_cold_pred<-read.delim("intersect_30bp_2kbdown_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_cold_pred<-read.delim("intersect_40bp_2kbdown_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_cold_pred<-read.delim("intersect_50bp_2kbdown_cold_a_pred.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}


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
intersect_overlap<-percentage_overlap(intersect_pred,genes)
intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred,genes)
intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred,genes)
intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred,genes)
intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred,genes)
intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred,genes)

intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred,genes_2kb_up)
intersect_10bp_2kbup_overlap<-percentage_overlap(intersect_10bp_2kbup_pred,genes_2kb_up)
intersect_20bp_2kbup_overlap<-percentage_overlap(intersect_20bp_2kbup_pred,genes_2kb_up)
intersect_30bp_2kbup_overlap<-percentage_overlap(intersect_30bp_2kbup_pred,genes_2kb_up)
intersect_40bp_2kbup_overlap<-percentage_overlap(intersect_40bp_2kbup_pred,genes_2kb_up)
intersect_50bp_2kbup_overlap<-percentage_overlap(intersect_50bp_2kbup_pred,genes_2kb_up)

intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred,genes_2kb_down)
intersect_10bp_2kbdown_overlap<-percentage_overlap(intersect_10bp_2kbdown_pred,genes_2kb_down)
intersect_20bp_2kbdown_overlap<-percentage_overlap(intersect_20bp_2kbdown_pred,genes_2kb_down)
intersect_30bp_2kbdown_overlap<-percentage_overlap(intersect_30bp_2kbdown_pred,genes_2kb_down)
intersect_40bp_2kbdown_overlap<-percentage_overlap(intersect_40bp_2kbdown_pred,genes_2kb_down)
intersect_50bp_2kbdown_overlap<-percentage_overlap(intersect_50bp_2kbdown_pred,genes_2kb_down)

intersect_cold_overlap<-percentage_overlap(intersect_cold_pred,genes_cold)
intersect_10bp_cold_overlap<-percentage_overlap(intersect_10bp_cold_pred,genes_cold)
intersect_20bp_cold_overlap<-percentage_overlap(intersect_20bp_cold_pred,genes_cold)
intersect_30bp_cold_overlap<-percentage_overlap(intersect_30bp_cold_pred,genes_cold)
intersect_40bp_cold_overlap<-percentage_overlap(intersect_40bp_cold_pred,genes_cold)
intersect_50bp_cold_overlap<-percentage_overlap(intersect_50bp_cold_pred,genes_cold)

intersect_2kbup_cold_overlap<-percentage_overlap(intersect_2kbup_cold_pred,genes_cold_2kb_up)
intersect_10bp_2kbup_cold_overlap<-percentage_overlap(intersect_10bp_2kbup_cold_pred,genes_cold_2kb_up)
intersect_20bp_2kbup_cold_overlap<-percentage_overlap(intersect_20bp_2kbup_cold_pred,genes_cold_2kb_up)
intersect_30bp_2kbup_cold_overlap<-percentage_overlap(intersect_30bp_2kbup_cold_pred,genes_cold_2kb_up)
intersect_40bp_2kbup_cold_overlap<-percentage_overlap(intersect_40bp_2kbup_cold_pred,genes_cold_2kb_up)
intersect_50bp_2kbup_cold_overlap<-percentage_overlap(intersect_50bp_2kbup_cold_pred,genes_cold_2kb_up)

intersect_2kbdown_cold_overlap<-percentage_overlap(intersect_2kbdown_cold_pred,genes_cold_2kb_down)
intersect_10bp_2kbdown_cold_overlap<-percentage_overlap(intersect_10bp_2kbdown_cold_pred,genes_cold_2kb_down)
intersect_20bp_2kbdown_cold_overlap<-percentage_overlap(intersect_20bp_2kbdown_cold_pred,genes_cold_2kb_down)
intersect_30bp_2kbdown_cold_overlap<-percentage_overlap(intersect_30bp_2kbdown_cold_pred,genes_cold_2kb_down)
intersect_40bp_2kbdown_cold_overlap<-percentage_overlap(intersect_40bp_2kbdown_cold_pred,genes_cold_2kb_down)
intersect_50bp_2kbdown_cold_overlap<-percentage_overlap(intersect_50bp_2kbdown_cold_pred,genes_cold_2kb_down)
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
  write.table(stat1,file="stat3_1_indel_hotspots.txt",row.names = T,quote=FALSE, sep = '\t')
  
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
  write.table(stat1_2,file="stat3_2_indel_hotspots.txt",row.names = T,quote=FALSE, sep = '\t')
  
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
  write.table(stat1_3,file="stat3_3_indel_hotspots.txt",row.names = T,quote=FALSE, sep = '\t')
  
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
  write.table(stat1_4,file="stat3_4_indel_hotspots.txt",row.names = T,quote=FALSE, sep = '\t')
}

{
#DATAFRAME OF CO PREDICTED
chrom <- c(rep(1:10,18) )
chrom<-paste("Chromosome",chrom,sep=" ")
type <- c(rep(c("all indels") ,30),rep(c("1-10bp") ,30),rep(c("11-20bp") ,30),rep(c("21-30bp") ,30),rep(c("31-40bp") ,30),rep(c("41-50bp") ,30))
intersect<-c(rep(c("Gene"),10),rep(c("TSS") ,10),rep(c("TTS"),10),rep(c("Gene"),10),rep(c("TSS") ,10),rep(c("TTS"),10),rep(c("Gene"),10),rep(c("TSS") ,10),rep(c("TTS"),10),rep(c("Gene"),10),rep(c("TSS") ,10),rep(c("TTS"),10),rep(c("Gene"),10),rep(c("TSS") ,10),rep(c("TTS"),10),rep(c("Gene"),10),rep(c("TSS") ,10),rep(c("TTS"),10))
value <- c(intersect_overlap,intersect_2kbup_overlap,intersect_2kbdown_overlap,intersect_10bp_overlap,intersect_10bp_2kbup_overlap,intersect_10bp_2kbdown_overlap,intersect_20bp_overlap,intersect_20bp_2kbup_overlap,intersect_20bp_2kbdown_overlap,intersect_30bp_overlap,intersect_30bp_2kbup_overlap,intersect_30bp_2kbdown_overlap,intersect_40bp_overlap,intersect_40bp_2kbup_overlap,intersect_40bp_2kbdown_overlap,intersect_50bp_overlap,intersect_50bp_2kbup_overlap,intersect_50bp_2kbdown_overlap)
data <- data.frame(chrom,type,value,intersect)
store<-mean(data$value)
data[data==0] <- NA
data<-na.omit(data)

#DATAFRAME OF CO DESERT
value <- c(intersect_cold_overlap,intersect_2kbup_cold_overlap,intersect_2kbdown_cold_overlap,intersect_10bp_cold_overlap,intersect_10bp_2kbup_cold_overlap,intersect_10bp_2kbdown_cold_overlap,intersect_20bp_cold_overlap,intersect_20bp_2kbup_cold_overlap,intersect_20bp_2kbdown_cold_overlap,intersect_30bp_cold_overlap,intersect_30bp_2kbup_cold_overlap,intersect_30bp_2kbdown_cold_overlap,intersect_40bp_cold_overlap,intersect_40bp_2kbup_cold_overlap,intersect_40bp_2kbdown_cold_overlap,intersect_50bp_cold_overlap,intersect_50bp_2kbup_cold_overlap,intersect_50bp_2kbdown_cold_overlap)
data_cold <- data.frame(chrom,type,value,intersect)
options(scipen = 999)
data_cold[is.na(data_cold)]<-0
store2<-mean(data_cold$value)
data_cold[data_cold==0] <- NA
data_cold<-na.omit(data_cold)
}

#SUMMARY TABLE OF ALL RESULTS
test <- matrix(nrow=6)
rownames(test)<-c("All indels","1-10bp","11-20bp","21-30bp","31-40bp","41-50bp")
test<-cbind(test,"Mean Indel / Gene Overlap"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Indel / Random Site Overlap"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test<-cbind(test,"Mean Gene Overlap"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Random Site Overlap"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test<-cbind(test,"Mean Indel Density on Genes (indels/bp)"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Indel Density on Random Site (indels/bp)"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test_2<-as.data.frame(test)[-c(1)]
write.table(test_2,file="hotspot_a_output.txt",row.names = T,quote=FALSE, sep = '\t')

#bargraphs of data
library(ggplot2)
library(viridis)
library(hrbrthemes)

perc<-ggplot(data, aes(fill=intersect,x=reorder(type, value),y=value)) +
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual(values = c("Gene" = "tan1",
                               "TSS" = "sienna1",
                               "TTS" = "brown1"))+
  facet_wrap(~chrom,ncol=5) +
  ggtitle("Indel Density (indels/bp) \nOn Genes Per Chromosome") +
  theme(plot.title = element_text(face="bold", size=20, hjust=0, color="#555555")) +
  geom_text(aes(label = round(value,digits=4)),position = position_stack(vjust = 0.5),size=2.5)+
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.position='top',legend.spacing.x= unit(1.0, 'cm')) +
  ylab("Percentage Overlap")+
  xlab("Type")+guides(fill=guide_legend("Type of intersection"))
perc
max<-max(data[which(data$intersect=="Gene"),]$value)+max(data[which(data$intersect=="TSS"),]$value)+max(data[which(data$intersect=="TTS"),]$value)

ggsave("indel_density_hotspot_feat.pdf", width = 20, height = 10)
perc_cold<-ggplot(data_cold, aes(fill=intersect,x=reorder(type, value),y=value)) +
  geom_bar(stat="identity",position="stack") +
  scale_fill_manual(values = c("Gene" = "turquoise2",
                               "TSS" = "deepskyblue",
                               "TTS" = "blue"))+
  facet_wrap(~chrom,ncol=5) +
  ggtitle("Indel Density (indels/bp) \nOn Random Sites Per Chromosome",subtitle = "Average Density: .00076") +
  theme(plot.title = element_text(face="bold", size=20, hjust=0, color="#555555")) +
  geom_text(aes(label = round(value,digits=4)),position = position_stack(vjust = 0.5),size=2.5)+
  theme(axis.text.x = element_text(angle=90)) +
  theme(legend.position='top',legend.spacing.x= unit(1.0, 'cm')) +
  ylab("Percentage Overlap")+
  xlab("Type")+guides(fill=guide_legend("Type of intersection"))+
  ylim(0,max)
perc_cold
ggsave("indel_density_coldspot_feat.pdf", width = 20, height = 10)


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