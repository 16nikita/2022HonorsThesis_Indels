library(dplyr)
#QUICK RELOAD
hotspots<-read.delim("hotspot.bed")
hotspot_2kb_up<-read.delim("hotspot_2kb_up.bed")
hotspot_2kb_down<-read.delim("hotspot_2kb_down.bed")
hotspot_random<-read.delim("hotspot_random.bed")
hotspot_random_2kb_up<-read.delim("hotspotRand_2kb_up.bed")
hotspot_random_2kb_down<-read.delim("hotspotRand_2kb_down.bed")

hotspots<-read.delim("hotspot_2.bed")
hotspot_2kb_up<-read.delim("hotspots_2_2kbup.bed")
hotspot_2kb_down<-read.delim("hotspots_2_2kbdown.bed")
hotspot_random<-read.delim("hotspot_2_random.bed")
hotspot_random_2kb_up<-read.delim("hotspotRand_2_2kb_up.bed")
hotspot_random_2kb_down<-read.delim("hotspotRand_2_2kb_down.bed")

hotspot_2kb_lowMeth<-read.delim("hotspots_2kb_lowMethylation.bed",header=FALSE)
colnames(hotspot_2kb_lowMeth)<-c('chrom','start','end','n_CO','CpG','CHG','Mnase')
hotspot_2kblowMeth_2kbUP<-CO_up(data.frame(matrix(ncol = 3,nrow = nrow(hotspot_2kb_lowMeth))),hotspot_2kb_lowMeth)
hotspot_2kblowMeth_2kbDOWN<-CO_down(data.frame(matrix(ncol = 3,nrow = nrow(hotspot_2kb_lowMeth))),hotspot_2kb_lowMeth)
write.table(hotspot_2kb_lowMeth,"hotspot_2kb_lowMeth.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(hotspot_2kblowMeth_2kbUP,"hotspot_2kblowMeth_UP.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(hotspot_2kblowMeth_2kbDOWN,"hotspot_2kblowMeth_DOWN.bed", row.names = F,quote=FALSE, sep = '\t')

#FOR METAPLOTS: extending all 10kb
{
hotspots_2<-hotspots[which((hotspots$end-hotspots$start)<=10000),]
hotspots_2$midpoint<-(hotspots_2$end+hotspots_2$start)/2
# for(i in 1:nrow(hotspots_2)){
#   if((hotspots_2$end[i]-hotspots_2$start[i])<10000){
#     hotspots_2$start[i]<-hotspots_2$midpoint[i]-5000
#     hotspots_2$end[i]<-hotspots_2$midpoint[i]+5000
#   }
# }
for(i in 1:nrow(hotspots_2)){
    hotspots_2$start[i]<-hotspots_2$midpoint[i]-5000
    hotspots_2$end[i]<-hotspots_2$midpoint[i]+5000
}
CO_up<-function(up,hotspots){
  colnames(up)<- c('chrom','start','end')
  up$chrom<-hotspots$chrom
  up$start<-hotspots$start-20000
  up$end<-hotspots$start
  up<-up[which(up$start>0),]
  return(up)
}
CO_down<-function(down,hotspots){
  colnames(down)<- c('chrom','start','end')
  down$chrom<-hotspots$chrom
  down$start<-hotspots$end
  down$end<-hotspots$end +20000
  return(down)
}
hotspots_2_2kbup<-CO_up(data.frame(matrix(ncol = 3,nrow = nrow(hotspots_2))),hotspots_2)
hotspots_2_2kbdown<-CO_down(data.frame(matrix(ncol = 3,nrow = nrow(hotspots_2))),hotspots_2)
# write.table(hotspots_2[-c(4)],"hotspot_2.bed", row.names = F,quote=FALSE, sep = '\t')
# write.table(hotspots_2_2kbup,"hotspots_2_2kbup.bed", row.names = F,quote=FALSE, sep = '\t')
# write.table(hotspots_2_2kbdown,"hotspots_2_2kbdown.bed", row.names = F,quote=FALSE, sep = '\t')
}
#Random data for Predicted COs- negative control, unique regions that do not overlap with hotspots
set.seed(200)
#Create dataframe of "coldspots" that do not overlap with one another, excluding hotspots, sample from list of CO desert regions
{
#READ IN DESERT REGIONS:
CO_coldspots<-read.delim("genome_length_without_centramere_sort_without_CO_distance_coldspot.bed",sep='\t', header = FALSE)
colnames(CO_coldspots)<-c('chrom','start','end','length')
library(dplyr)
library(reshape)
create_coldspot<-function(rows,sample,chrom,size){
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
    coldspots$start<-sample(coldspots$start,size=rows)
  }
  coldspots<-coldspots %>% distinct(start, .keep_all = TRUE)
  coldspots<-coldspots[order(coldspots$start),]
  for(i in 1:(nrow(coldspots)-1)){
    if(coldspots$start[i]+size>=coldspots$start[i+1]){
      diff<-size-abs(coldspots$start[i+1]-coldspots$start[i])
      coldspots$start[i+1]<-coldspots$start[i+1]+diff
    }
  }
  coldspots<-coldspots[order(coldspots$start),]
  coldspots$chrom<-chrom
  coldspots$end<-coldspots$start+size
  return(coldspots)
}
#store row counts
chrom1<- nrow(hotspots_2[which(hotspots_2$chrom == 1),])
chrom2<- nrow(hotspots_2[which(hotspots_2$chrom == 2),])
chrom3<- nrow(hotspots_2[which(hotspots_2$chrom == 3),])
chrom4<- nrow(hotspots_2[which(hotspots_2$chrom == 4),])
chrom5<- nrow(hotspots_2[which(hotspots_2$chrom == 5),])
chrom6<- nrow(hotspots_2[which(hotspots_2$chrom == 6),])
chrom7<- nrow(hotspots_2[which(hotspots_2$chrom == 7),])
chrom8<- nrow(hotspots_2[which(hotspots_2$chrom == 8),])
chrom9<- nrow(hotspots_2[which(hotspots_2$chrom == 9),])
chrom10<- nrow(hotspots_2[which(hotspots_2$chrom == 10),])

CO_coldspots_1<-create_coldspot(chrom1,CO_coldspots[which(CO_coldspots$chrom == 1),],1,10000)
CO_coldspots_2<-create_coldspot(chrom2,CO_coldspots[which(CO_coldspots$chrom == 2),],2,10000)
CO_coldspots_3<-create_coldspot(chrom3,CO_coldspots[which(CO_coldspots$chrom == 3),],3,10000)
CO_coldspots_4<-create_coldspot(chrom4,CO_coldspots[which(CO_coldspots$chrom == 4),],4,10000)
CO_coldspots_5<-create_coldspot(chrom5,CO_coldspots[which(CO_coldspots$chrom == 5),],5,10000)
CO_coldspots_6<-create_coldspot(chrom6,CO_coldspots[which(CO_coldspots$chrom == 6),],6,10000)
CO_coldspots_7<-create_coldspot(chrom7,CO_coldspots[which(CO_coldspots$chrom == 7),],7,10000)
CO_coldspots_8<-create_coldspot(chrom8,CO_coldspots[which(CO_coldspots$chrom == 8),],8,10000)
CO_coldspots_9<-create_coldspot(chrom9,CO_coldspots[which(CO_coldspots$chrom == 9),],9,10000)
CO_coldspots_10<-create_coldspot(chrom10,CO_coldspots[which(CO_coldspots$chrom == 10),],10,10000)

#Bind to create dataframe with all random regions
hotspot_random<-rbind(CO_coldspots_1,CO_coldspots_2,CO_coldspots_3,CO_coldspots_4,CO_coldspots_5,CO_coldspots_6,CO_coldspots_7,CO_coldspots_8,CO_coldspots_9,CO_coldspots_10)

#PREP COLDSPOT DATAFRAME FOR CROSSMAP: (convert v2 to v4)
hotspot_random$midpoint<-as.integer((hotspot_random$start+hotspot_random$end)/2)
hotspot_random$length<-hotspot_random$end-hotspot_random$midpoint
hotspot_random$start<-hotspot_random$midpoint
hotspot_random$end<-hotspot_random$midpoint+1
write.table(hotspot_random,"hotspot_2_random_v2.bed", row.names = F,quote=FALSE, sep = '\t')

#READING CO COlDSPOTS CROSSMAP OUTPUT:
hotspot_random<-read.delim("hotspot_2_random_v4.bed",sep='\t', header = FALSE)
colnames(hotspot_random)<-c("chrom","start","end","midpoint","length")
hotspot_random$end<-hotspot_random$start+hotspot_random$length
hotspot_random$start<-hotspot_random$start-hotspot_random$length
hotspot_random$length<-hotspot_random$end-hotspot_random$start
hotspot_random<-hotspot_random[!grepl("B73", hotspot_random$chrom),]
hotspot_random<-hotspot_random[which(hotspot_random$start>0),]
write.table(hotspot_random,"hotspot_2_random.bed", row.names = F,quote=FALSE, sep = '\t')

#Creating datasets of 2kb upstream and 2kb downstream from CO "coldspot" intervals
hotspot_random_2kb_up <- CO_up(data.frame(matrix(ncol = 3,nrow = nrow(hotspot_random))),hotspot_random)
hotspot_random_2kb_down <-CO_down(data.frame(matrix(ncol = 3,nrow = nrow(hotspot_random))),hotspot_random)
write.table(hotspot_random_2kb_up,"hotspotRand_2_2kb_up.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(hotspot_random_2kb_down,"hotspotRand_2_2kb_down.bed", row.names = F,quote=FALSE, sep = '\t')
}

##READ IN INTERSECT DATA:
# INTERSECT WITH PREDICTED CO REGIONS
{
  intersect_pred<-read.delim("intersect_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_pred$type<-"all indels"
  intersect_pred$intersect<-"intersect"
  colnames(intersect_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_pred<-read.delim("intersect_10bp_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_pred$type<-"1-10bp"
  intersect_10bp_pred$intersect<-"intersect"
  colnames(intersect_10bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_pred<-read.delim("intersect_20bp_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_pred$type<-"11-20bp"
  intersect_20bp_pred$intersect<-"intersect"
  colnames(intersect_20bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_pred<-read.delim("intersect_30bp_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_pred$type<-"21-30bp"
  intersect_30bp_pred$intersect<-"intersect"
  colnames(intersect_30bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_pred<-read.delim("intersect_40bp_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_pred$type<-"31-40bp"
  intersect_40bp_pred$intersect<-"intersect"
  colnames(intersect_40bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_pred<-read.delim("intersect_50bp_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_pred$type<-"41-50bp"
  intersect_50bp_pred$intersect<-"intersect"
  colnames(intersect_50bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_pred<-read.delim("intersect_2kbup_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_2kbup_pred$type<-"all deletions"
  intersect_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_pred<-read.delim("intersect_10bp_2kbup_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_pred$type<-"1-10bp"
  intersect_10bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_pred<-read.delim("intersect_20bp_2kbup_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_pred$type<-"11-20bp"
  intersect_20bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_pred<-read.delim("intersect_30bp_2kbup_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_pred$type<-"21-30bp"
  intersect_30bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_pred<-read.delim("intersect_40bp_2kbup_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_pred$type<-"31-40bp"
  intersect_40bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_pred<-read.delim("intersect_50bp_2kbup_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_pred$type<-"41-50bp"
  intersect_50bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_pred<-read.delim("intersect_2kbdown_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_2kbdown_pred$type<-"all deletions"
  intersect_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_pred<-read.delim("intersect_10bp_2kbdown_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_pred<-read.delim("intersect_20bp_2kbdown_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_pred<-read.delim("intersect_30bp_2kbdown_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_pred<-read.delim("intersect_40bp_2kbdown_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_pred<-read.delim("intersect_50bp_2kbdown_a_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}
# INTERSECT WITH RANDOM CO DESERT REGIONS
{
  intersect_cold_pred<-read.delim("intersect_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_cold_pred$type<-"all deletions"
  intersect_cold_pred$intersect<-"intersect"
  intersect_cold_pred<-intersect_cold_pred[-c(7:8)]
  colnames(intersect_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_cold_pred<-read.delim("intersect_10bp_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_10bp_cold_pred$type<-"1-10bp"
  intersect_10bp_cold_pred$intersect<-"intersect"
  intersect_10bp_cold_pred<-intersect_10bp_cold_pred[-c(7:8)]
  colnames(intersect_10bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_cold_pred<-read.delim("intersect_20bp_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_20bp_cold_pred$type<-"11-20bp"
  intersect_20bp_cold_pred$intersect<-"intersect"
  intersect_20bp_cold_pred<-intersect_20bp_cold_pred[-c(7:8)]
  colnames(intersect_20bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_cold_pred<-read.delim("intersect_30bp_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_30bp_cold_pred$type<-"21-30bp"
  intersect_30bp_cold_pred$intersect<-"intersect"
  intersect_30bp_cold_pred<-intersect_30bp_cold_pred[-c(7:8)]
  colnames(intersect_30bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_cold_pred<-read.delim("intersect_40bp_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_40bp_cold_pred$type<-"31-40bp"
  intersect_40bp_cold_pred$intersect<-"intersect"
  intersect_40bp_cold_pred<-intersect_40bp_cold_pred[-c(7:8)]
  colnames(intersect_40bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_cold_pred<-read.delim("intersect_50bp_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_50bp_cold_pred$type<-"41-50bp"
  intersect_50bp_cold_pred$intersect<-"intersect"
  intersect_50bp_cold_pred<-intersect_50bp_cold_pred[-c(7:8)]
  colnames(intersect_50bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_cold_pred<-read.delim("intersect_2kbup_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_2kbup_cold_pred$type<-"all deletions"
  intersect_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_cold_pred<-read.delim("intersect_10bp_2kbup_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_cold_pred<-read.delim("intersect_20bp_2kbup_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_cold_pred<-read.delim("intersect_30bp_2kbup_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_cold_pred<-read.delim("intersect_40bp_2kbup_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_cold_pred<-read.delim("intersect_50bp_2kbup_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_cold_pred<-read.delim("intersect_2kbdown_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_2kbdown_cold_pred$type<-"all deletions"
  intersect_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_cold_pred<-read.delim("intersect_10bp_2kbdown_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_cold_pred<-read.delim("intersect_20bp_2kbdown_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_cold_pred<-read.delim("intersect_30bp_2kbdown_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_cold_pred<-read.delim("intersect_40bp_2kbdown_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_cold_pred<-read.delim("intersect_50bp_2kbdown_a_coldspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}

##READ IN INTERSECT DATA FOR SUBSET OF HOTSPOTS:
# INTERSECT WITH PREDICTED CO REGIONS
{
  intersect_pred<-read.delim("intersect_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_pred$type<-"all deletions"
  intersect_pred$intersect<-"intersect"
  colnames(intersect_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_pred<-read.delim("intersect_10bp_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_pred$type<-"1-10bp"
  intersect_10bp_pred$intersect<-"intersect"
  colnames(intersect_10bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_pred<-read.delim("intersect_20bp_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_pred$type<-"11-20bp"
  intersect_20bp_pred$intersect<-"intersect"
  colnames(intersect_20bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_pred<-read.delim("intersect_30bp_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_pred$type<-"21-30bp"
  intersect_30bp_pred$intersect<-"intersect"
  colnames(intersect_30bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_pred<-read.delim("intersect_40bp_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_pred$type<-"31-40bp"
  intersect_40bp_pred$intersect<-"intersect"
  colnames(intersect_40bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_pred<-read.delim("intersect_50bp_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_pred$type<-"41-50bp"
  intersect_50bp_pred$intersect<-"intersect"
  colnames(intersect_50bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_pred<-read.delim("intersect_2kbup_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_2kbup_pred$type<-"all deletions"
  intersect_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_pred<-read.delim("intersect_10bp_2kbup_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_pred$type<-"1-10bp"
  intersect_10bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_pred<-read.delim("intersect_20bp_2kbup_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_pred$type<-"11-20bp"
  intersect_20bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_pred<-read.delim("intersect_30bp_2kbup_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_pred$type<-"21-30bp"
  intersect_30bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_pred<-read.delim("intersect_40bp_2kbup_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_pred$type<-"31-40bp"
  intersect_40bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_pred<-read.delim("intersect_50bp_2kbup_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_pred$type<-"41-50bp"
  intersect_50bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_pred<-read.delim("intersect_2kbdown_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_2kbdown_pred$type<-"all deletions"
  intersect_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_pred<-read.delim("intersect_10bp_2kbdown_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_pred<-read.delim("intersect_20bp_2kbdown_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_pred<-read.delim("intersect_30bp_2kbdown_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_pred<-read.delim("intersect_40bp_2kbdown_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_pred<-read.delim("intersect_50bp_2kbdown_sub_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}
# INTERSECT WITH RANDOM CO DESERT REGIONS
{
  intersect_cold_pred<-read.delim("intersect_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_cold_pred$type<-"all deletions"
  intersect_cold_pred$intersect<-"intersect"
  intersect_cold_pred<-intersect_cold_pred[-c(7:8)]
  colnames(intersect_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_cold_pred<-read.delim("intersect_10bp_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_10bp_cold_pred$type<-"1-10bp"
  intersect_10bp_cold_pred$intersect<-"intersect"
  intersect_10bp_cold_pred<-intersect_10bp_cold_pred[-c(7:8)]
  colnames(intersect_10bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_cold_pred<-read.delim("intersect_20bp_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_20bp_cold_pred$type<-"11-20bp"
  intersect_20bp_cold_pred$intersect<-"intersect"
  intersect_20bp_cold_pred<-intersect_20bp_cold_pred[-c(7:8)]
  colnames(intersect_20bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_cold_pred<-read.delim("intersect_30bp_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_30bp_cold_pred$type<-"21-30bp"
  intersect_30bp_cold_pred$intersect<-"intersect"
  intersect_30bp_cold_pred<-intersect_30bp_cold_pred[-c(7:8)]
  colnames(intersect_30bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_cold_pred<-read.delim("intersect_40bp_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_40bp_cold_pred$type<-"31-40bp"
  intersect_40bp_cold_pred$intersect<-"intersect"
  intersect_40bp_cold_pred<-intersect_40bp_cold_pred[-c(7:8)]
  colnames(intersect_40bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_cold_pred<-read.delim("intersect_50bp_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_50bp_cold_pred$type<-"41-50bp"
  intersect_50bp_cold_pred$intersect<-"intersect"
  intersect_50bp_cold_pred<-intersect_50bp_cold_pred[-c(7:8)]
  colnames(intersect_50bp_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbup_cold_pred<-read.delim("intersect_2kbup_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_2kbup_cold_pred$type<-"all deletions"
  intersect_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_cold_pred<-read.delim("intersect_10bp_2kbup_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_cold_pred<-read.delim("intersect_20bp_2kbup_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_cold_pred<-read.delim("intersect_30bp_2kbup_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_cold_pred<-read.delim("intersect_40bp_2kbup_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_cold_pred<-read.delim("intersect_50bp_2kbup_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbup_cold_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_cold_pred<-read.delim("intersect_2kbdown_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_2kbdown_cold_pred$type<-"all deletions"
  intersect_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_cold_pred<-read.delim("intersect_10bp_2kbdown_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_cold_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_cold_pred<-read.delim("intersect_20bp_2kbdown_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_cold_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_cold_pred<-read.delim("intersect_30bp_2kbdown_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_cold_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_cold_pred<-read.delim("intersect_40bp_2kbdown_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_cold_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_cold_pred<-read.delim("intersect_50bp_2kbdown_sub_coldspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_cold_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_cold_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_cold_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}

##READ IN INTERSECT DATA FOR LOWEST METHYLATED 2KB HOTSPOTS:
# INTERSECT WITH PREDICTED CO REGIONS
{

  intersect_pred<-read.delim("intersect_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_pred$type<-"all deletions"
  intersect_pred$intersect<-"intersect"
  colnames(intersect_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','n_CO','CpG','ChG','Mnase','overlap in bp','type','intersection')
  intersect_10bp_pred<-read.delim("intersect_10bp_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_pred$type<-"1-10bp"
  intersect_10bp_pred$intersect<-"intersect"
  colnames(intersect_10bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','n_CO','CpG','ChG','Mnase','overlap in bp','type','intersection')
  intersect_20bp_pred<-read.delim("intersect_20bp_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_pred$type<-"11-20bp"
  intersect_20bp_pred$intersect<-"intersect"
  colnames(intersect_20bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','n_CO','CpG','ChG','Mnase','overlap in bp','type','intersection')
  intersect_30bp_pred<-read.delim("intersect_30bp_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_pred$type<-"21-30bp"
  intersect_30bp_pred$intersect<-"intersect"
  colnames(intersect_30bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','n_CO','CpG','ChG','Mnase','overlap in bp','type','intersection')
  intersect_40bp_pred<-read.delim("intersect_40bp_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_pred$type<-"31-40bp"
  intersect_40bp_pred$intersect<-"intersect"
  colnames(intersect_40bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','n_CO','CpG','ChG','Mnase','overlap in bp','type','intersection')
  intersect_50bp_pred<-read.delim("intersect_50bp_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_pred$type<-"41-50bp"
  intersect_50bp_pred$intersect<-"intersect"
  colnames(intersect_50bp_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','n_CO','CpG','ChG','Mnase','overlap in bp','type','intersection')
  
  intersect_2kbup_pred<-read.delim("intersect_2kbup_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_2kbup_pred$type<-"all deletions"
  intersect_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbup_pred<-read.delim("intersect_10bp_2kbup_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbup_pred$type<-"1-10bp"
  intersect_10bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_10bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbup_pred<-read.delim("intersect_20bp_2kbup_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbup_pred$type<-"11-20bp"
  intersect_20bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_20bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbup_pred<-read.delim("intersect_30bp_2kbup_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbup_pred$type<-"21-30bp"
  intersect_30bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_30bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbup_pred<-read.delim("intersect_40bp_2kbup_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbup_pred$type<-"31-40bp"
  intersect_40bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_40bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbup_pred<-read.delim("intersect_50bp_2kbup_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbup_pred$type<-"41-50bp"
  intersect_50bp_2kbup_pred$intersect<-"2kb up"
  colnames(intersect_50bp_2kbup_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  
  intersect_2kbdown_pred<-read.delim("intersect_2kbdown_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_2kbdown_pred$type<-"all deletions"
  intersect_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_10bp_2kbdown_pred<-read.delim("intersect_10bp_2kbdown_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_10bp_2kbdown_pred$type<-"1-10bp"
  intersect_10bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_10bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_20bp_2kbdown_pred<-read.delim("intersect_20bp_2kbdown_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_20bp_2kbdown_pred$type<-"11-20bp"
  intersect_20bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_20bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_30bp_2kbdown_pred<-read.delim("intersect_30bp_2kbdown_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_30bp_2kbdown_pred$type<-"21-30bp"
  intersect_30bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_30bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_40bp_2kbdown_pred<-read.delim("intersect_40bp_2kbdown_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_40bp_2kbdown_pred$type<-"31-40bp"
  intersect_40bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_40bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
  intersect_50bp_2kbdown_pred<-read.delim("intersect_50bp_2kbdown_lowMETH_hotspot.txt",sep='\t', header = FALSE)
  intersect_50bp_2kbdown_pred$type<-"41-50bp"
  intersect_50bp_2kbdown_pred$intersect<-"2kb down"
  colnames(intersect_50bp_2kbdown_pred) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type','intersection')
}

indels_hotspots<-intersect_pred[c(1:3)]
colnames(indels_hotspots)<-c("chrom","start","end")
indels_coldhotspots<-intersect_cold_pred[c(1:3)]
colnames(indels_coldhotspots)<-c("chrom","start","end")
write.table(indels_hotspots,file="indels_hotspots.bed",row.names = F,quote=FALSE, sep = '\t')
write.table(indels_coldhotspots,file="indels_coldhotspots.bed",row.names = F,quote=FALSE, sep = '\t')


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
intersect_overlap<-percentage_overlap(intersect_pred,hotspots)
intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred,hotspots)
intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred,hotspots)
intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred,hotspots)
intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred,hotspots)
intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred,hotspots)

intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred,hotspots_2kbup)
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
##STATISTIC 3: mean indel density
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
# lin_regression <- lm(as.numeric(n_CO) ~ as.numeric(Mnase) + as.numeric(percent), data = intersect_overlap)
# summary(lin_regression)
#write.table(intersect_overlap,file="hotspot_2kb_lowMeth_indel_density.bed",row.names = T,quote=FALSE, sep = '\t')
intersect_10bp_overlap<-percentage_overlap(intersect_10bp_pred)
intersect_20bp_overlap<-percentage_overlap(intersect_20bp_pred)
intersect_30bp_overlap<-percentage_overlap(intersect_30bp_pred)
intersect_40bp_overlap<-percentage_overlap(intersect_40bp_pred)
intersect_50bp_overlap<-percentage_overlap(intersect_50bp_pred)

intersect_2kbup_overlap<-percentage_overlap(intersect_2kbup_pred)
#write.table(intersect_2kbup_overlap,file="hotspot_2kb_lowMeth_2kbUP_indel_density.bed",row.names = T,quote=FALSE, sep = '\t')
intersect_10bp_2kbup_overlap<-percentage_overlap(intersect_10bp_2kbup_pred)
intersect_20bp_2kbup_overlap<-percentage_overlap(intersect_20bp_2kbup_pred)
intersect_30bp_2kbup_overlap<-percentage_overlap(intersect_30bp_2kbup_pred)
intersect_40bp_2kbup_overlap<-percentage_overlap(intersect_40bp_2kbup_pred)
intersect_50bp_2kbup_overlap<-percentage_overlap(intersect_50bp_2kbup_pred)

intersect_2kbdown_overlap<-percentage_overlap(intersect_2kbdown_pred)
#write.table(intersect_2kbdown_overlap,file="hotspot_2kb_lowMeth_2kbDOWN_indel_density.bed",row.names = T,quote=FALSE, sep = '\t')
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
type <- c(rep(c("all deletions") ,30),rep(c("1-10bp") ,30),rep(c("11-20bp") ,30),rep(c("21-30bp") ,30),rep(c("31-40bp") ,30),rep(c("41-50bp") ,30))
intersect<-c(rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10))
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
test<-cbind(test,"Mean Indel / Hotspot Overlap"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Indel / Random Site Overlap"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test<-cbind(test,"Mean Hotspot Overlap"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Random Site Overlap"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test<-cbind(test,"Mean Indel Density on Hotspots (indels/bp)"=c(mean(data[which(data$type == "all deletions"),]$value),mean(data[which(data$type == "1-10bp"),]$value),mean(data[which(data$type == "11-20bp"),]$value),mean(data[which(data$type == "21-30bp"),]$value),mean(data[which(data$type == "31-40bp"),]$value),mean(data[which(data$type == "41-50bp"),]$value)))
test<-cbind(test,"Mean Indel Density on Random Site (indels/bp)"=c(mean(data_cold[which(data_cold$type == "all deletions"),]$value),mean(data_cold[which(data_cold$type == "1-10bp"),]$value),mean(data_cold[which(data_cold$type == "11-20bp"),]$value),mean(data_cold[which(data_cold$type == "21-30bp"),]$value),mean(data_cold[which(data_cold$type == "31-40bp"),]$value),mean(data_cold[which(data_cold$type == "41-50bp"),]$value)))

test_2<-as.data.frame(test)[-c(1)]
write.table(test_2,file="hotspot_a_output.txt",row.names = T,quote=FALSE, sep = '\t')

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
  full_data1$posthoc<-c("AB","aB","Ab","ab","Ab","ab","Ab","ab","Ab","ab")
  #plot 1: indel density 
  perc<-ggplot(full_data1, aes(fill=color,x=levels,y=value)) +
    geom_bar(stat="identity",position="dodge") +
    scale_fill_manual(labels = c("CO Hotspot",
                                 "Negative Control Region"),
                      values = c("Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("Small Indel Density on CO Hotspots") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(),width=.2)+
    geom_text(aes(label=posthoc, y=value+sd), position=position_dodge(.5), vjust=-.5) +
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(min(full_data1$value - full_data1$sd),max(full_data1$value + full_data1$sd))) +
    annotate("text", c(1.5,3.5,5.5,7.5,9.5), y = -.002, label = c("1-10 bp","11-20 bp","21-30 bp","31-40 bp","41-50 bp")) +
    annotate("text", x = 5.5, y = -.004,
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
  full_data2$posthoc<-c("AB","Ab","Ab","aC","ac","ac")
  perc2<-ggplot(full_data2, aes(fill=color,x=levels,y=value)) +
    geom_bar(stat="identity",position="dodge", width = .5) +
    scale_fill_manual(labels = c("2kb downstream", "2kb downstream",
                                 "2kb upstream","2kb upstream",
                                 "CO Hotspot",
                                 "Negative Control"),
                      values = c("2kb downstream hot"= "darkgoldenrod1",
                                 "2kb downstream cold" = "deepskyblue",
                                 "2kb upstream hot" = "darkorange1",
                                 "2kb upstream cold" = "dodgerblue2",
                                 "Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("Small Indel Density on CO Hotspots") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(.5),width=.2)+
    geom_text(aes(label=posthoc, y=value+sd), position=position_dodge(.5), vjust=-.5) +
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(0,max(full_data2$value + full_data2$sd))) +
    annotate("text", x = 1:2, y = -.002,
             label = rep(c("Negative Control Region", "CO Hotspot Region"), 1)) +
    annotate("text", x = 1.5, y = -.003,
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
  write.table(intersect_overlap,file="hotspot_indel_density.bed",row.names = T,quote=FALSE, sep = '\t')
  
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
    interval<-plot(smooth.spline(seq(-5000,4999,10),test,spar=0.5),col="brown1",lwd=2,type="l",ylim=c(0,max),main="Indel Densities at Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")+
      lines(smooth.spline(seq(-5000,4999,10),test2,spar=0.5),col="royalblue3",lwd=2,type="l",main="Indel Densities at Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
    par(mar=c(5, 4, 4, 8), xpd=TRUE)
    legend("topright", inset=c(-.35, 0), c("Interval hot", "Interval cold"), title="Dataset",
           col=c("brown1", "royalblue3"), lty=1,lwd=2,cex = .8)

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
    twokbdown<-plot(smooth.spline(seq(-10000,9999,100),test,spar=0.5),col="brown1",lwd=2,type="l",ylim=c(0,max),main="Indel Densities 20kb down from Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")+
      lines(smooth.spline(seq(-10000,9999,100),test2,spar=0.5),col="royalblue3",lwd=2,type="l",main="Indel Densities 20kb down from Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
    par(mar=c(5, 4, 4, 8), xpd=TRUE)
    legend("topright", inset=c(-.35, 0), c("Interval hot", "Interval cold"), title="Dataset",
           col=c("brown1", "royalblue3"), lty=1,lwd=2,cex = .8)
    
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
    twokbup<-plot(smooth.spline(seq(-10000,9999,100),test,spar=0.5),col="brown1",lwd=2,type="l",ylim=c(0,max),main="Indel Densities 20kb up from Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")+
      lines(smooth.spline(seq(-10000,9999,100),test2,spar=0.5),col="royalblue3",lwd=2,type="l",main="Indel Densities 20kb up from Eli's Hotspots",xlab ="Distance from midpoint", ylab="Indel Density (indels/bp)")
    par(mar=c(5, 4, 4, 8), xpd=TRUE)
    legend("topright", inset=c(-.35, 0), c("Interval hot", "Interval cold"), title="Dataset",
           col=c("brown1", "royalblue3"), lty=1,lwd=2,cex = .8)
    #return(twokbdown)
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