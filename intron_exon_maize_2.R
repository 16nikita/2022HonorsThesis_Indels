library(dplyr)
#QUICK RELOAD
introns<-read.delim("intron_v4.bed",sep='\t', header = TRUE)
exons<-read.delim("exon_v4.bed",sep='\t', header = TRUE)


##READ IN INTERSECT DATA: small indels
{
all_chrom_intron_feat<-read.delim("all_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
all_chrom_intron_feat$type<-"intron"
chrom_intron_feat_10bp<-read.delim("10bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
chrom_intron_feat_10bp$type<-"intron"
chrom_intron_feat_20bp<-read.delim("20bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
chrom_intron_feat_20bp$type<-"intron"
chrom_intron_feat_30bp<-read.delim("30bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
chrom_intron_feat_30bp$type<-"intron"
chrom_intron_feat_40bp<-read.delim("40bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
chrom_intron_feat_40bp$type<-"intron"
chrom_intron_feat_50bp<-read.delim("50bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
chrom_intron_feat_50bp$type<-"intron"

all_chrom_exon_feat<-read.delim("all_chrom_exon_feat.txt",sep='\t', header = FALSE)
all_chrom_exon_feat$type<-"exon"
chrom_exon_feat_10bp<-read.delim("10bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
chrom_exon_feat_10bp$type<-"exon"
chrom_exon_feat_20bp<-read.delim("20bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
chrom_exon_feat_20bp$type<-"exon"
chrom_exon_feat_30bp<-read.delim("30bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
chrom_exon_feat_30bp$type<-"exon"
chrom_exon_feat_40bp<-read.delim("40bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
chrom_exon_feat_40bp$type<-"exon"
chrom_exon_feat_50bp<-read.delim("50bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
chrom_exon_feat_50bp$type<-"exon"
}
##READ IN INTERSECT DATA: medium indels
{
  all_M_chrom_intron_feat<-read.delim("all_M_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  all_M_chrom_intron_feat$type<-"intron"
  chrom_intron_feat_100bp<-read.delim("100bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_100bp$type<-"intron"
  chrom_intron_feat_200bp<-read.delim("200bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_200bp$type<-"intron"
  chrom_intron_feat_300bp<-read.delim("300bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_300bp$type<-"intron"
  chrom_intron_feat_400bp<-read.delim("400bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_400bp$type<-"intron"
  chrom_intron_feat_500bp<-read.delim("500bp_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_500bp$type<-"intron"
  
  all_M_chrom_exon_feat<-read.delim("all_M_chrom_exon_feat.txt",sep='\t', header = FALSE)
  all_M_chrom_exon_feat$type<-"exon"
  chrom_exon_feat_100bp<-read.delim("100bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_100bp$type<-"exon"
  chrom_exon_feat_200bp<-read.delim("200bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_200bp$type<-"exon"
  chrom_exon_feat_300bp<-read.delim("300bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_300bp$type<-"exon"
  chrom_exon_feat_400bp<-read.delim("400bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_400bp$type<-"exon"
  chrom_exon_feat_500bp<-read.delim("500bp_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_500bp$type<-"exon"
}
##READ IN INTERSECT DATA: large indels
{
  all_L_chrom_intron_feat<-read.delim("all_L_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  all_L_chrom_intron_feat$type<-"intron"
  chrom_intron_feat_10kb<-read.delim("10kb_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_10kb$type<-"intron"
  chrom_intron_feat_20kb<-read.delim("20kb_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_20kb$type<-"intron"
  chrom_intron_feat_30kb<-read.delim("30kb_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_30kb$type<-"intron"
  chrom_intron_feat_40kb<-read.delim("40kb_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_40kb$type<-"intron"
  chrom_intron_feat_50kb<-read.delim("50kb_chrom_intron_feat.txt",sep='\t', header = FALSE)[-c(7:8)]
  chrom_intron_feat_50kb$type<-"intron"
  
  all_L_chrom_exon_feat<-read.delim("all_L_chrom_exon_feat.txt",sep='\t', header = FALSE)
  all_L_chrom_exon_feat$type<-"exon"
  chrom_exon_feat_10kb<-read.delim("10kb_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_10kb$type<-"exon"
  chrom_exon_feat_20kb<-read.delim("20kb_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_20kb$type<-"exon"
  chrom_exon_feat_30kb<-read.delim("30kb_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_30kb$type<-"exon"
  chrom_exon_feat_40kb<-read.delim("40kb_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_40kb$type<-"exon"
  chrom_exon_feat_50kb<-read.delim("50kb_chrom_exon_feat.txt",sep='\t', header = FALSE)
  chrom_exon_feat_50kb$type<-"exon"
}

##STATISTIC 1: INDEL OVERLAP
#Calculate percentage of Indels that overlap with COs
indels_NAMv4<-read.delim(file="indels_all.bed")
indels_10bp<-read.delim(file="indels_10bp_all.bed")
indels_20bp<-read.delim(file="indels_20bp_all.bed")
indels_30bp<-read.delim(file="indels_30bp_all.bed")
indels_40bp<-read.delim(file="indels_40bp_all.bed")
indels_50bp<-read.delim(file="indels_50bp_all.bed")

indels_NAMv4<-read.delim(file="indelsLarge_to500bp.bed")
indels_10bp<-read.delim(file="indelsLarge_100bp.bed")
indels_20bp<-read.delim(file="indelsLarge_200bp.bed")
indels_30bp<-read.delim(file="indelsLarge_300bp.bed")
indels_40bp<-read.delim(file="indelsLarge_400bp.bed")
indels_50bp<-read.delim(file="indelsLarge_500bp.bed")

indels_NAMv4<-read.delim(file="indelsLarge.bed")
indels_10bp<-read.delim(file="indelsLarge_10kb.bed")
indels_20bp<-read.delim(file="indelsLarge_20kb.bed")
indels_30bp<-read.delim(file="indelsLarge_30kb.bed")
indels_40bp<-read.delim(file="indelsLarge_40kb.bed")
indels_50bp<-read.delim(file="indelsLarge_50kb.bed")

percentage_overlap<-function(intersect,indels,type){
  colnames(intersect) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type')
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
intersect_exon_overlap<-percentage_overlap(all_chrom_exon_feat,indels_NAMv4,"else")
intersect_10bp_exon_overlap<-percentage_overlap(chrom_exon_feat_10bp,indels_10bp,"else")
intersect_20bp_exon_overlap<-percentage_overlap(chrom_exon_feat_20bp,indels_20bp,"else")
intersect_30bp_exon_overlap<-percentage_overlap(chrom_exon_feat_30bp,indels_30bp,"else")
intersect_40bp_exon_overlap<-percentage_overlap(chrom_exon_feat_40bp,indels_40bp,"else")
intersect_50bp_exon_overlap<-percentage_overlap(chrom_exon_feat_50bp,indels_50bp,"else")

intersect_intron_overlap<-percentage_overlap(all_chrom_intron_feat,indels_NAMv4,"else")
intersect_10bp_intron_overlap<-percentage_overlap(chrom_intron_feat_10bp,indels_10bp,"else")
intersect_20bp_intron_overlap<-percentage_overlap(chrom_intron_feat_20bp,indels_20bp,"else")
intersect_30bp_intron_overlap<-percentage_overlap(chrom_intron_feat_30bp,indels_30bp,"else")
intersect_40bp_intron_overlap<-percentage_overlap(chrom_intron_feat_40bp,indels_40bp,"else")
intersect_50bp_intron_overlap<-percentage_overlap(chrom_intron_feat_50bp,indels_50bp,"else")

}
library(matrixStats)
{
stat1<-as.data.frame(cbind(intersect_exon_overlap,intersect_intron_overlap))
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

stat1_2<-as.data.frame(cbind( intersect_10bp_exon_overlap,intersect_10bp_intron_overlap,intersect_20bp_exon_overlap,intersect_20bp_intron_overlap,intersect_30bp_exon_overlap,intersect_30bp_intron_overlap,intersect_40bp_exon_overlap, intersect_40bp_intron_overlap,intersect_50bp_exon_overlap,intersect_50bp_intron_overlap))
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
}

##STATISTIC 2: mean indel density
library(dplyr)
#mean indel density on each feature
percentage_overlap<-function(intersect){
  colnames(intersect) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type')
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
    percentage_per_chrom[i]<- mean(intersect_temp$percent)
  }
  print(percentage_per_chrom)
}

{
  intersect_exon_overlap<-percentage_overlap(all_chrom_exon_feat)
  intersect_10bp_exon_overlap<-percentage_overlap(chrom_exon_feat_10bp)
  intersect_20bp_exon_overlap<-percentage_overlap(chrom_exon_feat_20bp)
  intersect_30bp_exon_overlap<-percentage_overlap(chrom_exon_feat_30bp)
  intersect_40bp_exon_overlap<-percentage_overlap(chrom_exon_feat_40bp)
  intersect_50bp_exon_overlap<-percentage_overlap(chrom_exon_feat_50bp)
  
  intersect_intron_overlap<-percentage_overlap(all_chrom_intron_feat)
  intersect_10bp_intron_overlap<-percentage_overlap(chrom_intron_feat_10bp)
  intersect_20bp_intron_overlap<-percentage_overlap(chrom_intron_feat_20bp)
  intersect_30bp_intron_overlap<-percentage_overlap(chrom_intron_feat_30bp)
  intersect_40bp_intron_overlap<-percentage_overlap(chrom_intron_feat_40bp)
  intersect_50bp_intron_overlap<-percentage_overlap(chrom_intron_feat_50bp)
  
  intersect_M_exon_overlap<-percentage_overlap(all_M_chrom_exon_feat)
  intersect_100bp_exon_overlap<-percentage_overlap(chrom_exon_feat_100bp)
  intersect_200bp_exon_overlap<-percentage_overlap(chrom_exon_feat_200bp)
  intersect_300bp_exon_overlap<-percentage_overlap(chrom_exon_feat_300bp)
  intersect_400bp_exon_overlap<-percentage_overlap(chrom_exon_feat_400bp)
  intersect_500bp_exon_overlap<-percentage_overlap(chrom_exon_feat_500bp)
  
  intersect_M_intron_overlap<-percentage_overlap(all_M_chrom_intron_feat)
  intersect_100bp_intron_overlap<-percentage_overlap(chrom_intron_feat_100bp)
  intersect_200bp_intron_overlap<-percentage_overlap(chrom_intron_feat_200bp)
  intersect_300bp_intron_overlap<-percentage_overlap(chrom_intron_feat_300bp)
  intersect_400bp_intron_overlap<-percentage_overlap(chrom_intron_feat_400bp)
  intersect_500bp_intron_overlap<-percentage_overlap(chrom_intron_feat_500bp)
  
  intersect_L_exon_overlap<-percentage_overlap(all_L_chrom_exon_feat)
  intersect_10kb_exon_overlap<-percentage_overlap(chrom_exon_feat_10kb)
  intersect_20kb_exon_overlap<-percentage_overlap(chrom_exon_feat_20kb)
  intersect_30kb_exon_overlap<-percentage_overlap(chrom_exon_feat_30kb)
  intersect_40kb_exon_overlap<-percentage_overlap(chrom_exon_feat_40kb)
  intersect_50kb_exon_overlap<-percentage_overlap(chrom_exon_feat_50kb)
  
  intersect_L_intron_overlap<-percentage_overlap(all_L_chrom_intron_feat)
  intersect_10kb_intron_overlap<-percentage_overlap(chrom_intron_feat_10kb)
  intersect_20kb_intron_overlap<-percentage_overlap(chrom_intron_feat_20kb)
  intersect_30kb_intron_overlap<-percentage_overlap(chrom_intron_feat_30kb)
  intersect_40kb_intron_overlap<-percentage_overlap(chrom_intron_feat_40kb)
  intersect_50kb_intron_overlap<-percentage_overlap(chrom_intron_feat_50kb)
  
}
#T-TEST:
{
  library(matrixStats)
  stat1<-as.data.frame(cbind(intersect_exon_overlap,intersect_intron_overlap,intersect_M_exon_overlap,intersect_M_intron_overlap,intersect_L_exon_overlap,intersect_L_intron_overlap))
  stat1[nrow(stat1) + 1,] <- colMeans(as.matrix(stat1))
  stat1[nrow(stat1) + 1,] <- colVars(as.matrix(stat1))
  stat1[nrow(stat1) + 1,] <- colSds(as.matrix(stat1))/sqrt(10)
  stat1[nrow(stat1) + 1,] <-NA
  rownames(stat1)<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","mean","variance","se","p-value")
  options(scipen=999)
  for(i in seq(from=1, to=ncol(stat1)-1, by=2)){
    ttest<-t.test(stat1[1:10,i] , stat1[1:10,i+1], paired=TRUE, var.equal=TRUE)
    stat1["p-value",i]<-ttest$p.value
  }
  
  type<-rep(c("Exon","Intron"),3)
  size<-c(rep(c("Small Indels"),2),rep(c("Medium Indels"),2),rep(c("Large Indels"),2))
  value<-unname(unlist(stat1[11,]))
  se<-unname(unlist(stat1[12,]))
  data <- data.frame(type,size,value,se)
  data$posthoc<-c("a","a","A","A","A","A")
  
  exon_analysis<- ggplot(data, aes(fill=size,x=type,y=value)) +
    geom_bar(stat="identity",position="dodge", width = .5)+
    geom_text(aes(label=posthoc, y=value+se), position=position_dodge(.5), vjust=-.5) +
    #geom_errorbar(aes(ymax = value + se, ymin = value - se), position = position_dodge(),width=.2)+
    scale_fill_manual(values = c("Large Indels"= "olivedrab3",
                                 "Medium Indels" = "steelblue1",
                                 "Small Indels" = "sandybrown"))+
    ggtitle("Comparison of Indel Density in Intronic and Exonic Regions")+
    ylab("Indel Density (indels/bp)")+theme_classic()+xlab("Region")
  exon_analysis
  
  #ANOVA test for significance
  {
    percentage_overlap<-function(intersect){
      colnames(intersect) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type')
      intersect<-intersect %>% distinct(A_start, .keep_all = TRUE)
      duplicated<-as.data.frame(table(intersect$'B_start'))
      intersect$length<-intersect$'B_end'-intersect$'B_start'
      intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
      intersect<-intersect[order(intersect$B_start),]
      intersect$count<-duplicated$Freq
      intersect$percent<- (intersect$count/intersect$length)
      return(intersect)
      # percentage_per_chrom<-c()
      # for(i in 1:10){
      #   intersect_temp<-intersect[which(intersect$`A_chrom` == i),]
      #   percentage_per_chrom[i]<- mean(intersect_temp$percent)
      # }
      # print(percentage_per_chrom)
    }
    {
      intersect_exon_overlap<-percentage_overlap(all_chrom_exon_feat)
      intersect_10bp_exon_overlap<-percentage_overlap(chrom_exon_feat_10bp)
      intersect_20bp_exon_overlap<-percentage_overlap(chrom_exon_feat_20bp)
      intersect_30bp_exon_overlap<-percentage_overlap(chrom_exon_feat_30bp)
      intersect_40bp_exon_overlap<-percentage_overlap(chrom_exon_feat_40bp)
      intersect_50bp_exon_overlap<-percentage_overlap(chrom_exon_feat_50bp)
      intersect_exon_overlap$size<-"s"
      intersect_10bp_exon_overlap$size<-"s"
      intersect_20bp_exon_overlap$size<-"s"
      intersect_30bp_exon_overlap$size<-"s"
      intersect_40bp_exon_overlap$size<-"s"
      intersect_50bp_exon_overlap$size<-"s"
      
      intersect_intron_overlap<-percentage_overlap(all_chrom_intron_feat)
      intersect_10bp_intron_overlap<-percentage_overlap(chrom_intron_feat_10bp)
      intersect_20bp_intron_overlap<-percentage_overlap(chrom_intron_feat_20bp)
      intersect_30bp_intron_overlap<-percentage_overlap(chrom_intron_feat_30bp)
      intersect_40bp_intron_overlap<-percentage_overlap(chrom_intron_feat_40bp)
      intersect_50bp_intron_overlap<-percentage_overlap(chrom_intron_feat_50bp)
      intersect_intron_overlap$size<-"s"
      intersect_10bp_intron_overlap$size<-"s"
      intersect_20bp_intron_overlap$size<-"s"
      intersect_30bp_intron_overlap$size<-"s"
      intersect_40bp_intron_overlap$size<-"s"
      intersect_50bp_intron_overlap$size<-"s"
      
      intersect_M_exon_overlap<-percentage_overlap(all_M_chrom_exon_feat)
      intersect_100bp_exon_overlap<-percentage_overlap(chrom_exon_feat_100bp)
      intersect_200bp_exon_overlap<-percentage_overlap(chrom_exon_feat_200bp)
      intersect_300bp_exon_overlap<-percentage_overlap(chrom_exon_feat_300bp)
      intersect_400bp_exon_overlap<-percentage_overlap(chrom_exon_feat_400bp)
      intersect_500bp_exon_overlap<-percentage_overlap(chrom_exon_feat_500bp)
      intersect_M_exon_overlap$size<-"m"
      intersect_100bp_exon_overlap$size<-"m"
      intersect_200bp_exon_overlap$size<-"m"
      intersect_300bp_exon_overlap$size<-"m"
      intersect_400bp_exon_overlap$size<-"m"
      intersect_500bp_exon_overlap$size<-"m"
      
      intersect_M_intron_overlap<-percentage_overlap(all_M_chrom_intron_feat)
      intersect_100bp_intron_overlap<-percentage_overlap(chrom_intron_feat_100bp)
      intersect_200bp_intron_overlap<-percentage_overlap(chrom_intron_feat_200bp)
      intersect_300bp_intron_overlap<-percentage_overlap(chrom_intron_feat_300bp)
      intersect_400bp_intron_overlap<-percentage_overlap(chrom_intron_feat_400bp)
      intersect_500bp_intron_overlap<-percentage_overlap(chrom_intron_feat_500bp)
      intersect_M_intron_overlap$size<-"m"
      intersect_100bp_intron_overlap$size<-"m"
      intersect_200bp_intron_overlap$size<-"m"
      intersect_300bp_intron_overlap$size<-"m"
      intersect_400bp_intron_overlap$size<-"m"
      intersect_500bp_intron_overlap$size<-"m"
      
      
      intersect_L_exon_overlap<-percentage_overlap(all_L_chrom_exon_feat)
      intersect_10kb_exon_overlap<-percentage_overlap(chrom_exon_feat_10kb)
      intersect_20kb_exon_overlap<-percentage_overlap(chrom_exon_feat_20kb)
      intersect_30kb_exon_overlap<-percentage_overlap(chrom_exon_feat_30kb)
      intersect_40kb_exon_overlap<-percentage_overlap(chrom_exon_feat_40kb)
      intersect_50kb_exon_overlap<-percentage_overlap(chrom_exon_feat_50kb)
      intersect_L_exon_overlap$size<-"L"
      intersect_10kb_exon_overlap$size<-"L"
      intersect_20kb_exon_overlap$size<-"L"
      intersect_30kb_exon_overlap$size<-"L"
      intersect_40kb_exon_overlap$size<-"L"
      intersect_50kb_exon_overlap$size<-"L"
      
      intersect_L_intron_overlap<-percentage_overlap(all_L_chrom_intron_feat)
      intersect_10kb_intron_overlap<-percentage_overlap(chrom_intron_feat_10kb)
      intersect_20kb_intron_overlap<-percentage_overlap(chrom_intron_feat_20kb)
      intersect_30kb_intron_overlap<-percentage_overlap(chrom_intron_feat_30kb)
      intersect_40kb_intron_overlap<-percentage_overlap(chrom_intron_feat_40kb)
      intersect_50kb_intron_overlap<-percentage_overlap(chrom_intron_feat_50kb)
      intersect_L_intron_overlap$size<-"L"
      intersect_10kb_intron_overlap$size<-"L"
      intersect_20kb_intron_overlap$size<-"L"
      intersect_30kb_intron_overlap$size<-"L"
      intersect_40kb_intron_overlap$size<-"L"
      intersect_50kb_intron_overlap$size<-"L"
      
    }
    stat2<-as.data.frame(rbind(intersect_exon_overlap,intersect_intron_overlap,intersect_M_exon_overlap,intersect_M_intron_overlap,intersect_L_exon_overlap,intersect_L_intron_overlap))
    data1 <- data.frame(y=(stat2$count),group=factor(stat2$size),group2=factor(stat2$type),offset=(stat2$length))
    #fit = lm(y ~ group*group2, data1)
    fit= glm(y~group * group2 ,data1,  family = poisson(link = "log"),offset(log(offset)))
    #* interaction effect- impact of hot/cold depends on indel size
    
    library(devtools)
    library(qqplotr)
    gg <- ggplot(data = data1, mapping = aes(sample = y)) +
      geom_qq_band(bandType = "ks", mapping = aes(fill = "KS"), alpha = 0.5) +
      geom_qq_band(bandType = "ts", mapping = aes(fill = "TS"), alpha = 0.5) +
      geom_qq_band(bandType = "pointwise", mapping = aes(fill = "Normal"), alpha = 0.5) +
      geom_qq_band(bandType = "boot", mapping = aes(fill = "Bootstrap"), alpha = 0.5) +
      stat_qq_line() +
      stat_qq_point() +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      scale_fill_discrete("Bandtype")
    gg
    
    two.way <- aov(fit)
    summary(two.way)
    
    #bonferroni for non normal distribution
    data1$group3<-paste(data1$group,data1$group2)
    pairwise.t.test(data1$y,data1$group3, p.adjust.method="bonferroni")

    tukey.test <- TukeyHSD(two.way)
    tukey.test
    plot(tukey.test)
    #ANOVA test indicated significant diff (F=0.15, p=0.861)
    
    library(car)
    leveneTest(fit)
    #Levene's test indicated equal variances (F = 2.5291, p = 0.03962)
    
    #anderson-darling normality test
    library(nortest)
    ad.test(data1$y) #departure from normality
    
    #Shapiro Wilk Test: p-value > 0.05 = normal 
    aov_residuals <- residuals(object = two.way)
    shapiro.test(data1$y)
    plot(two.way, 1) #check residual v fitted for outliers
    plot(two.way, 2) #check if residuals are along straight line
    #normality condition violated!! p<0.05
    
  }
  
}



TE_randPred<-read.delim("intersect_TE_randPred.txt",sep='\t', header = FALSE)[-c(4:5,9:10)]
TE_randPred$type<-"pred co"
TE_randPCO<-read.delim("intersect_TE_randPCO.txt",sep='\t', header = FALSE)[-c(4:5,9:10)]
TE_randPCO$type<-"PCO"
TE_coldspot<-read.delim("intersect_TE_coldspot.txt",sep='\t', header = FALSE)[-c(4:5,9:10)]
TE_coldspot$type<-"eli"
TE_randDSB<-read.delim("intersect_TE_randDSB.txt",sep='\t', header = FALSE)[-c(4:5)]
TE_randDSB$type<-"emp dsb"

randPredDSB<-read.delim("pred_DSB_coldspotsv4.bed",sep='\t', header = TRUE)
randRCO<-read.delim("RCO_random_2.bed",sep='\t', header = TRUE)
randPred<-read.delim("CO_pred_coldspots_2.bed",sep='\t', header = TRUE)
randNAM<-read.delim("CO_random.bed",sep='\t', header = TRUE)
coldspot<-read.delim("hotspot_random.bed",sep='\t', header = TRUE)
randDSB<-read.delim("DSB_coldspotsv4.bed",sep='\t', header = TRUE)


TE_Pred<-read.delim("intersect_TE_Pred.txt",sep='\t', header = FALSE)[-c(4:5)]
TE_Pred$type<-"pred co"
TE_PCO<-read.delim("intersect_TE_PCO.txt",sep='\t', header = FALSE)[-c(4:5,9:10)]
TE_PCO$type<-"PCO"
TE_hotspot<-read.delim("intersect_TE_hotspot.txt",sep='\t', header = FALSE)[-c(4:5)]
TE_hotspot$type<-"eli"
TE_DSB<-read.delim("intersect_TE_DSB.txt",sep='\t', header = FALSE)[-c(4:5,9:10)]
TE_DSB$type<-"emp dsb"

PredDSB<-read.delim("pred_DSB_v4_3.bed",sep='\t', header = FALSE)
RCO<-read.delim("Ruth_COs_2.bed",sep='\t', header = FALSE)
Pred<-read.delim("CO_predicted_3.bed",sep='\t', header = FALSE)
NAM<-read.delim("NAM_COs.bed",sep='\t', header = TRUE)
hotspot<-read.delim("hotspot_2.bed",sep='\t', header = TRUE)
DSB<-read.delim("DSB_V4.bed",sep='\t', header = TRUE)[c(1:3)]


##STATISTIC 3: % of cold regions that overlap with TEs
library(dplyr)
percentage_overlap<-function(intersect,CO){
  colnames(intersect) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','overlap in bp','type')
  colnames(CO)<-c('chrom','start','end')
  intersect<-intersect %>% distinct()
  intersect<-intersect %>% distinct(B_start, .keep_all = TRUE)
  percentage_per_chrom<-c()
  for(i in 1:10){
    percentage_per_chrom[i]<- (nrow(intersect[which(intersect$`A_chrom` == i),])/nrow(CO[which(CO$`chrom` == i),]))*100 
  }
  print(percentage_per_chrom)
}
#coldspots overlap
perc_randPredDSB<-percentage_overlap(TE_randPredDSB,randPredDSB)
perc_randRCO<-percentage_overlap(TE_randRCO,randRCO)
perc_randPred<-percentage_overlap(TE_randPred,randPred)
perc_randNAM<-percentage_overlap(TE_randNAM,randNAM)
perc_coldspot<-percentage_overlap(TE_coldspot,coldspot)
perc_randDSB<-percentage_overlap(TE_randDSB,randDSB)
perc_coldwTE<-rbind(perc_randPredDSB,perc_randRCO,perc_randPred,perc_randNAM,perc_coldspot,perc_randDSB)
colnames(perc_coldwTE)<-c(1,2,3,4,5,6,7,8,9,10)
perc_coldwTE<-cbind(perc_coldwTE,mean=rowMeans(perc_coldwTE),sd=rowSds(perc_coldwTE), median=rowMedians(perc_coldwTE))

#TE density
perc_randPred<-percentage_overlap(TE_randPred)
perc_randPCO<-percentage_overlap(TE_randPCO)
perc_coldspot<-percentage_overlap(TE_coldspot)
perc_randDSB<-percentage_overlap(TE_randDSB)
perc_TE_cold_density<-rbind(perc_randPred,perc_randPCO,perc_coldspot,perc_randDSB)
colnames(perc_TE_cold_density)<-c(1,2,3,4,5,6,7,8,9,10)
perc_TE_cold_density<-as.data.frame(cbind(perc_TE_cold_density,mean=rowMeans(perc_TE_cold_density),sd=rowSds(perc_TE_cold_density), median=rowMedians(perc_TE_cold_density)))[-c(1:10)]
perc_Pred<-percentage_overlap(TE_Pred)
perc_PCO<-percentage_overlap(TE_PCO)
perc_hotspot<-percentage_overlap(TE_hotspot)
perc_DSB<-percentage_overlap(TE_DSB)
perc_hot_TEdensity<-rbind(perc_Pred,perc_PCO,perc_hotspot,perc_DSB)
colnames(perc_hot_TEdensity)<-c(1,2,3,4,5,6,7,8,9,10)
perc_hot_TEdensity<-as.data.frame(cbind(perc_hot_TEdensity,mean=rowMeans(perc_hot_TEdensity),sd=rowSds(perc_hot_TEdensity), median=rowMedians(perc_hot_TEdensity)))[-c(1:10)]

perc_TEdensity<-as.data.frame(cbind(perc_TE_cold_density,perc_hot_TEdensity))
row.names(perc_TEdensity)<-c('Predicted Crossovers',"Penny's Crossovers","Eli's Hotspots","DSBs" )
colnames(perc_TEdensity)<-c("Mean TE Density (TEs/bp) on Hot Region","SD of TE Density (TEs/bp) on Hot Region","Median TE Density (TEs/bp) on Hot Region","Mean TE Density (TEs/bp) on Cold Region","SD TE Density (TEs/bp) on Cold Region","Median TE Density (TEs/bp) on Cold Region")
write.csv(perc_TEdensity,file="perc_TEdensity.csv",row.names = T,quote=FALSE)
#ANOVA test for significance
{
  intersect<-c(rep(c("Cold"),10),rep(c("Hot"),10))
  value <- c((perc_randPCO),(perc_PCO))
  data1 <- data.frame(y=value,group=factor(intersect))
  fit = lm(y ~ group, data1)
  #p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.
  shapiro.test(value)
  
  #kruskal-wallis- no normal distribution, equal variance
  kruskal.test(y ~ group, data = data1)
  
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
  #Levene's test indicated unequal variances (F = 0.0116, p = .9885)
  one.way <- aov(fit)
  summary(one.way)
  #ANOVA test indicated no significant diff (F=0.15, p=0.861)
}


#hotspots
perc_PredDSB<-percentage_overlap(TE_PredDSB,PredDSB)
perc_RCO<-percentage_overlap(TE_RCO,RCO)
perc_Pred<-percentage_overlap(TE_Pred,Pred)
perc_NAM<-percentage_overlap(TE_NAM,NAM)
perc_coldspot<-percentage_overlap(TE_hotspot,hotspot)
perc_DSB<-percentage_overlap(TE_DSB,DSB)
perc_hotwTE<-rbind(perc_PredDSB,perc_RCO,perc_Pred,perc_NAM,perc_coldspot,perc_DSB)
colnames(perc_hotwTE)<-c(1,2,3,4,5,6,7,8,9,10)
perc_hotwTE<-cbind(perc_hotwTE,mean=rowMeans(perc_hotwTE),sd=rowSds(perc_hotwTE), median=rowMedians(perc_hotwTE))

#TE density

write.csv(perc_hotwTE,file="perc_hotwTE.csv",row.names = T,quote=FALSE)
write.csv(perc_TEdensity,file="perc_TEdensity.csv",row.names = T,quote=FALSE)

#check indel density on open chromatin
intersect_open<-read.delim("intersect_open_pred.txt",sep='\t', header = FALSE)
intersect_10bp_open<-read.delim("intersect_10bp_open_pred.txt",sep='\t', header = FALSE)
intersect_20bp_open<-read.delim("intersect_20bp_open_pred.txt",sep='\t', header = FALSE)
intersect_30bp_open<-read.delim("intersect_30bp_open_pred.txt",sep='\t', header = FALSE)
intersect_40bp_open<-read.delim("intersect_40bp_open_pred.txt",sep='\t', header = FALSE)
intersect_50bp_open<-read.delim("intersect_50bp_open_pred.txt",sep='\t', header = FALSE)

indels_open<-percentage_overlap(intersect_open)
#cor.test(indels_open$percent, indels_open$occupancy,  method = "spearman", alternative = "greater")
lin_regression <- lm(log(percent) ~ occupancy, data = indels_open)
summary(lin_regression)

indels_10bp_open<-percentage_overlap(intersect_10bp_open)
lin_regression2 <- lm(percent ~ occupancy, data = indels_10bp_open)
summary(lin_regression2)

indels_20bp_open<-percentage_overlap(intersect_20bp_open)
indels_30bp_open<-percentage_overlap(intersect_30bp_open)
indels_40bp_open<-percentage_overlap(intersect_40bp_open)
indels_50bp_open<-percentage_overlap(intersect_50bp_open)
open_chrom_density<-rbind(indels_open,indels_10bp_open,indels_20bp_open,indels_30bp_open,indels_40bp_open,indels_50bp_open)
colnames(open_chrom_density)<-c(1,2,3,4,5,6,7,8,9,10)
open_chrom_density<-cbind(open_chrom_density,mean=rowMeans(open_chrom_density),sd=rowSds(open_chrom_density), median=rowMedians(open_chrom_density))
##STATISTIC 4:  mean indel density on open chromatin
library(dplyr)
percentage_overlap<-function(intersect){
  colnames(intersect) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','peak #','occupancy','overlap in bp')
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
  return(intersect)
}

