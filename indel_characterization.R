#indel characterization
TTS<-read.delim("TTS_v4_region.bed",header=FALSE)[-c(4:5)]
TSS<-read.delim("TSS_v4_region.bed",header=FALSE)[-c(4:5)]
genes<-read.delim("gene_v4_region.bed",header=FALSE)[-c(4:5)]
TTS$type <-"TTS"
TSS$type <-"TSS"
genes$type <-"gene"

chrom_feat<-rbind(TTS,TSS,genes)
colnames(chrom_feat)<-c("chrom","start","end","type")
chrom_feat<-chrom_feat[order(chrom_feat$chrom,chrom_feat$start),]
chrom_feat<-chrom_feat[!grepl("B73", chrom_feat$chrom),]
chrom_feat<-chrom_feat[!grepl("Pt", chrom_feat$chrom),]
chrom_feat<-chrom_feat[!grepl("Mt", chrom_feat$chrom),]
write.table(chrom_feat,file="chrom_feat.bed",row.names = F,quote=FALSE, sep = '\t')
chrom_feat$length<-chrom_feat$end-chrom_feat$start

all_chrom_feat<-read.delim("all_chrom_feat.txt",sep='\t', header = FALSE)
chrom_feat_10bp<-read.delim("10bp_chrom_feat.txt",sep='\t', header = FALSE)
chrom_feat_20bp<-read.delim("20bp_chrom_feat.txt",sep='\t', header = FALSE)
chrom_feat_30bp<-read.delim("30bp_chrom_feat.txt",sep='\t', header = FALSE)
chrom_feat_40bp<-read.delim("40bp_chrom_feat.txt",sep='\t', header = FALSE)
chrom_feat_50bp<-read.delim("50bp_chrom_feat.txt",sep='\t', header = FALSE)

all_chrom_feat_medium<-read.delim("all_chrom_medium_feat.txt",sep='\t', header = FALSE)
chrom_feat_10bp_medium<-read.delim("10bp_chrom_medium_feat.txt",sep='\t', header = FALSE)
chrom_feat_20bp_medium<-read.delim("20bp_chrom_medium_feat.txt",sep='\t', header = FALSE)
chrom_feat_30bp_medium<-read.delim("30bp_chrom_medium_feat.txt",sep='\t', header = FALSE)
chrom_feat_40bp_medium<-read.delim("40bp_chrom_medium_feat.txt",sep='\t', header = FALSE)
chrom_feat_50bp_medium<-read.delim("50bp_chrom_medium_feat.txt",sep='\t', header = FALSE)

all_chrom_feat_large<-read.delim("all_chrom_large_large_feat.txt",sep='\t', header = FALSE)
chrom_feat_10kb<-read.delim("10bp_chrom_large_large_feat.txt",sep='\t', header = FALSE)
chrom_feat_20kb<-read.delim("20bp_chrom_large_large_feat.txt",sep='\t', header = FALSE)
chrom_feat_30kb<-read.delim("30bp_chrom_large_large_feat.txt",sep='\t', header = FALSE)
chrom_feat_40kb<-read.delim("40bp_chrom_large_large_feat.txt",sep='\t', header = FALSE)
chrom_feat_50kb<-read.delim("50bp_chrom_large_large_feat.txt",sep='\t', header = FALSE)

indels_NAMv4_large<-read.delim(file="indelsLarge.bed")
indels_10bp_large<-read.delim(file="indelsLarge_10kb.bed")
indels_20bp_large<-read.delim(file="indelsLarge_20kb.bed")
indels_30bp_large<-read.delim(file="indelsLarge_30kb.bed")
indels_40bp_large<-read.delim(file="indelsLarge_40kb.bed")
indels_50bp_large<-read.delim(file="indelsLarge_50kb.bed")

indels_NAMv4<-read.delim(file="indels_all.bed")
indels_10bp<-read.delim(file="indels_10bp_all.bed")
indels_20bp<-read.delim(file="indels_20bp_all.bed")
indels_30bp<-read.delim(file="indels_30bp_all.bed")
indels_40bp<-read.delim(file="indels_40bp_all.bed")
indels_50bp<-read.delim(file="indels_50bp_all.bed")

indels_NAMv4_medium<-read.delim(file="indelsLarge_to500bp.bed")
indels_10bp_medium<-read.delim(file="indelsLarge_100bp.bed")
indels_20bp_medium<-read.delim(file="indelsLarge_200bp.bed")
indels_30bp_medium<-read.delim(file="indelsLarge_300bp.bed")
indels_40bp_medium<-read.delim(file="indelsLarge_400bp.bed")
indels_50bp_medium<-read.delim(file="indelsLarge_500bp.bed")

chrom_feat<-read.delim(file="chrom_feat.bed")
TSS_count<-nrow(chrom_feat[which(chrom_feat$type=="TSS"),])
TTS_count<-nrow(chrom_feat[which(chrom_feat$type=="TTS"),])
gene_count<-nrow(chrom_feat[which(chrom_feat$type=="gene"),])

where<-function(chrom,indels){
  colnames(chrom) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','type','overlap in bp')
  table<-as.data.frame(matrix(data=NA,nrow=1,ncol=4))
  colnames(table)<-c("% indels TSS","% indels TTS","% indels Genic","% indels Non-Genic")
  chrom<-chrom %>% distinct()
  chrom_2<-chrom %>% distinct(A_start, .keep_all = TRUE)
  TSS_count_1<-nrow(chrom_2[which(chrom_2$type=="TSS"),])
  TTS_count_1<-nrow(chrom_2[which(chrom_2$type=="TTS"),])
  gene_count_1<-nrow(chrom_2[which(chrom_2$type=="gene"),])
  
  indels<-indels %>% distinct()
  indel_count<-nrow(indels)
  
  table$`% indels TSS`<-(TSS_count_1/indel_count)*100
  table$`% indels TTS`<-(TTS_count_1/indel_count)*100
  table$`% indels Genic`<-(gene_count_1/indel_count)*100
  table$`% indels Non-Genic`<-100-(table$`% indels TSS`+table$`% indels TTS`+table$`% indels Genic`)
    return(table)
}


where_all_indels<-where(all_chrom_feat,indels_NAMv4)
where_10bp_indels<-where(chrom_feat_10bp,indels_10bp)
where_20bp_indels<-where(chrom_feat_20bp,indels_20bp)
where_30bp_indels<-where(chrom_feat_30bp,indels_30bp)
where_40bp_indels<-where(chrom_feat_40bp,indels_40bp)
where_50bp_indels<-where(chrom_feat_50bp,indels_50bp)

where_are_they<-rbind(where_all_indels,where_10bp_indels,where_20bp_indels,where_30bp_indels,where_40bp_indels,where_50bp_indels)
rownames(where_are_they)<-c("all indels","10bp indels","20bp indels","30bp indels","40bp indels","50bp indels")

where_all_indels_medium<-where(all_chrom_feat_medium,indels_NAMv4_medium)
where_10bp_indels_medium<-where(chrom_feat_10bp_medium,indels_10bp_medium)
where_20bp_indels_medium<-where(chrom_feat_20bp_medium,indels_20bp_medium)
where_30bp_indels_medium<-where(chrom_feat_30bp_medium,indels_30bp_medium)
where_40bp_indels_medium<-where(chrom_feat_40bp_medium,indels_40bp_medium)
where_50bp_indels_medium<-where(chrom_feat_50bp_medium,indels_50bp_medium)

where_are_they_medium<-rbind(where_all_indels_medium,where_10bp_indels_medium,where_20bp_indels_medium,where_30bp_indels_medium,where_40bp_indels_medium,where_50bp_indels_medium)
rownames(where_are_they_medium)<-c("all indels","100-150bp indels","151-200bp indels","201-300bp indels","301-400bp indels","401-500bp indels")

where_all_Large_indels<-where(all_chrom_feat_large,indels_NAMv4_large)
where_10kb_indels<-where(chrom_feat_10kb,indels_10bp_large)
where_20kb_indels<-where(chrom_feat_20kb,indels_20bp_large)
where_30kb_indels<-where(chrom_feat_30kb,indels_30bp_large)
where_40kb_indels<-where(chrom_feat_40kb,indels_40bp_large)
where_50kb_indels<-where(chrom_feat_50kb,indels_50bp_large)

where_are_they_large<-rbind(where_all_Large_indels,where_10kb_indels,where_20kb_indels,where_30kb_indels,where_40kb_indels,where_50kb_indels)
rownames(where_are_they_large)<-c("all indels","500bp-10kb indels","20kb indels","30kb indels","40kb indels","50kb indels")

write.csv(where_are_they,file="where_are_they.csv",row.names = T,quote=FALSE)
write.csv(where_are_they_large,file="where_are_they_large.csv",row.names = T,quote=FALSE)
write.csv(where_are_they_medium,file="where_are_they_medium.csv",row.names = T,quote=FALSE)
