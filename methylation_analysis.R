
#QUICK RELOAD
{
  #'cold' datasets
  DSB_coldspotsv4<-read.delim("DSB_coldspotsv4.bed")
  CO_random_2<-read.delim("CO_random_2.bed")
  CO_pred_coldspots<-read.delim("CO_pred_coldspots.bed")
  hotspot_2_random<-read.delim("hotspot_2_random.bed")
  RCO_random<-read.delim("RCO_random.bed")
  cold_pred_DSB_v4<-read.delim("pred_DSB_coldspotsv4.bed")
  
  #'warm' datasets
  DSB_V4<-read.delim("DSB_V4.bed")
  NAM_COs_2<-read.delim("NAM_COs_2.bed")
  CO_predicted<-read.delim("CO_predicted.bed")
  hotspot_2<-read.delim("hotspot_2.bed")
  Ruth_COs<-read.delim("Ruth_COs.bed")
  pred_DSB_v4<-read.delim("pred_DSB_v4.bed")
}
#EXTEND DATASETS FOR METAPLOT:
{
extend_10kb<-function(data){
  data$midpoint<-(data$start+data$end)/2
  for(i in 1:nrow(data)){
    data$start[i]<-data$midpoint[i]-5000
    data$end[i]<-data$midpoint[i]+5000
  }
  data<-data[which(data$start>0),]
  data$length<-(data$end-data$start)
  return(data)
}
extend_20kb<-function(data){
  data$midpoint<-(data$start+data$end)/2
  for(i in 1:nrow(data)){
    data$start[i]<-data$midpoint[i]-10000
    data$end[i]<-data$midpoint[i]+10000
  }
  data<-data[which(data$start>0),]
  return(data)
}

Predicted_DSB<-extend_10kb(pred_DSB_v4)
Predicted_DSB_Cold<-extend_10kb(cold_pred_DSB_v4)

DSB<-extend_10kb(DSB_V4)
DSB_cold<-extend_10kb(DSB_coldspotsv4)

NAM_CO<-extend_20kb(pred_DSB_v4)
NAM_CO_cold<-extend_20kb(pred_DSB_v4)

Eli_Hotspot<-extend_20kb(hotspot_2)
Eli_Hotspot_cold<-extend_20kb(hotspot_2_random)

Ruth_Hotspot<-extend_10kb(Ruth_COs)
Ruth_Hotspot_cold<-extend_10kb(RCO_random)

Pred_CO<-extend_10kb(CO_predicted)
Pred_CO_cold<-extend_10kb(CO_pred_coldspots)

write.table(Predicted_DSB,"Predicted_DSB_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(Predicted_DSB_Cold,"Predicted_DSB_Cold_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(DSB,"DSB_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(DSB_cold,"DSB_cold_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(NAM_CO,"NAM_CO_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(NAM_CO_cold,"NAM_CO_cold_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(Eli_Hotspot,"Eli_hotspot_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(Eli_Hotspot_cold,"Eli_hotspot_cold_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(Ruth_Hotspot,"RCO_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(Ruth_Hotspot_cold,"RCO_cold_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(Pred_CO,"Pred_CO_extend.bed", row.names = F,quote=FALSE, sep = '\t')
write.table(Pred_CO_cold,"Pred_CO_cold_extend.bed", row.names = F,quote=FALSE, sep = '\t')
}

CHG<-read.delim("CHG_2.bed")
CHH<-read.delim("CHH_2.bed")
CpG<-read.delim("CpG_2.bed")



no_ind_dsb<-read.delim("no_indel_DSB.bed",sep='\t', header = FALSE)
no_ind_hotspot<-read.delim("no_indel_hotspot.bed",sep='\t', header = FALSE)
no_ind_predCO<-read.delim("no_indel_predCO.bed",sep='\t', header = FALSE)
no_ind_PCO<-read.delim("no_indel_empCO.bed",sep='\t', header = FALSE)

##READ IN INTERSECT DATA:
{  
  CHG_methylation_DSB_hot<-read.delim("CHGmethylation_DSB_hot.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_DSB_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHGmethylation_PredCO_hot<-read.delim("CHGmethylation_PredCO_hot.txt",sep='\t', header = FALSE)
  colnames(CHGmethylation_PredCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHG_methylation_hotspot_hot<-read.delim("CHGmethylation_hotspot_hot.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_hotspot_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHG_methylation_PCO_hot<-read.delim("CHGmethylation_PCO_hot.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_PCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','B_midpoint','B_length','overlap_in_bp')
  
  CHH_methylation_DSB_hot<-read.delim("CHHmethylation_DSB_hot.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_DSB_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHHmethylation_PredCO_hot<-read.delim("CHHmethylation_PredCO_hot.txt",sep='\t', header = FALSE)
  colnames(CHHmethylation_PredCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHH_methylation_hotspot_hot<-read.delim("CHHmethylation_hotspot_hot.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_hotspot_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHH_methylation_PCO_hot<-read.delim("CHHmethylation_PCO_hot.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_PCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','B_midpoint','B_length','overlap_in_bp')
  
  CpG_methylation_DSB_hot<-read.delim("CpGmethylation_DSB_hot.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_DSB_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpGmethylation_PredCO_hot<-read.delim("CpGmethylation_PredCO_hot.txt",sep='\t', header = FALSE)
  colnames(CpGmethylation_PredCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpG_methylation_hotspot_hot<-read.delim("CpGmethylation_hotspot_hot.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_hotspot_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpG_methylation_PCO_hot<-read.delim("CpGmethylation_PCO_hot.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_PCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','B_midpoint','B_length','overlap_in_bp')
  
  CpG_methylation_DSB<-read.delim("CpGmethylation_DSB.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_DSB) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','B_midpoint','B_length','overlap_in_bp')
  CpGmethylation_PredCO<-read.delim("CpGmethylation_PredCO.txt",sep='\t', header = FALSE)
  colnames(CpGmethylation_PredCO) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpG_methylation_hotspot<-read.delim("CpGmethylation_hotspot.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_hotspot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpG_methylation_PCO<-read.delim("CpGmethylation_PCO.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_PCO) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  
  CHG_methylation_DSB<-read.delim("CHGmethylation_DSB.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_DSB) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','B_midpoint','B_length','overlap_in_bp')
  CHGmethylation_PredCO<-read.delim("CHGmethylation_PredCO.txt",sep='\t', header = FALSE)
  colnames(CHGmethylation_PredCO) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHG_methylation_hotspot<-read.delim("CHGmethylation_hotspot.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_hotspot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHG_methylation_PCO<-read.delim("CHGmethylation_PCO.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_PCO) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  
  CHH_methylation_DSB<-read.delim("CHHmethylation_DSB.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_DSB) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','B_midpoint','B_length','overlap_in_bp')
  CHHmethylation_PredCO<-read.delim("CHHmethylation_PredCO.txt",sep='\t', header = FALSE)
  colnames(CHHmethylation_PredCO) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHH_methylation_hotspot<-read.delim("CHHmethylation_hotspot.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_hotspot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHH_methylation_PCO<-read.delim("CHHmethylation_PCO.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_PCO) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  
}

##READ IN INTERSECT DATA for methylation of indels:
{  
  CHG_methylation_DSB_hot<-read.delim("CHGmethylation_DSB_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_DSB_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHGmethylation_NAM_hot<-read.delim("CHGmethylation_NAM_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHGmethylation_NAM_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHGmethylation_PredCO_hot<-read.delim("CHGmethylation_PredCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHGmethylation_PredCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHG_methylation_hotspot_hot<-read.delim("CHGmethylation_hotspot_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_hotspot_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHG_methylation_RCO_hot<-read.delim("CHGmethylation_RCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_RCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHG_methylation_PCO_hot<-read.delim("CHGmethylation_PCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHG_methylation_PCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  
   CHH_methylation_DSB_hot<-read.delim("CHHmethylation_DSB_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_DSB_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHHmethylation_NAM_hot<-read.delim("CHHmethylation_NAM_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHHmethylation_NAM_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHHmethylation_PredCO_hot<-read.delim("CHHmethylation_PredCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHHmethylation_PredCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHH_methylation_hotspot_hot<-read.delim("CHHmethylation_hotspot_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_hotspot_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHH_methylation_RCO_hot<-read.delim("CHHmethylation_RCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_RCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CHH_methylation_PCO_hot<-read.delim("CHHmethylation_PCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CHH_methylation_PCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  
  CpG_methylation_DSB_hot<-read.delim("CpGmethylation_DSB_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_DSB_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpGmethylation_NAM_hot<-read.delim("CpGmethylation_NAM_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CpGmethylation_NAM_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpGmethylation_PredCO_hot<-read.delim("CpGmethylation_PredCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CpGmethylation_PredCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpG_methylation_hotspot_hot<-read.delim("CpGmethylation_hotspot_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_hotspot_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpG_methylation_RCO_hot<-read.delim("CpGmethylation_RCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_RCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  CpG_methylation_PCO_hot<-read.delim("CpGmethylation_PCO_hot_indels.txt",sep='\t', header = FALSE)
  colnames(CpG_methylation_PCO_hot) <- c('A_chrom', 'A_start', 'A_end','% methylated','B_chrom', 'B_start', 'B_end','overlap_in_bp')
  
}

##STATISTIC 5: metaplots of indel density
library(dplyr)
#mean indel density on each CO
{
  
  metaplot_10kb<-function(all, all_cold){
    all<-all[order(all$B_chrom,all$B_start),]
    all_cold<-all_cold[order(all_cold$B_chrom,all_cold$B_start),]
    
    M=matrix(((all$B_start+all$B_end)/2)-all$A_start,ncol=1)
    all=cbind(all,M)
    all=all[all$M<5000,]
    tapply(all$'% methylated', cut(all$M, seq(-5000,5000, by=100)), mean)->mean_all
    test<-as.data.frame(mean_all)
    test$mean_all[is.na(test$mean_all)]<-0
    test<-unlist(test)
    min<-min(test)
    
    M=matrix(((all_cold$B_start+all_cold$B_end)/2)- all_cold$A_start,ncol=1)
    all_cold=cbind(all_cold,M)
    all_cold=all_cold[all_cold$M<5000,]
    tapply(all_cold$'% methylated', cut(all_cold$M, seq(-5000,5000, by=100)), mean)->mean_all_cold
    test2<-as.data.frame(mean_all_cold)
    test2$mean_all_cold[is.na(test2$mean_all_cold)]<-0
    test2<-unlist(test2)
    max<-max(test2)
    plot(smooth.spline(seq(-5000,4999,100),test,spar = 0.5),col="red",lwd=1,type="l",ylim=c(min,max),main="% CHG Methylation at Predicted COs",xlab ="Distance from midpoint", ylab="Methylation Level (%)")
    lines(smooth.spline(seq(-5000,4999,100),test2,spar = 0.5),col="blue",lwd=1,type="l",main="% CHG Methylation at Predicted COs",xlab ="Distance from midpoint", ylab="Methylation Level (%)")

  }
  metaplot_20kb<-function(all, all_cold){
    all<-all[order(all$B_chrom,all$B_start),]
    all_cold<-all_cold[order(all_cold$B_chrom,all_cold$B_start),]
    
    M=matrix(((all$B_start+all$B_end)/2)-all$A_start,ncol=1)
    all=cbind(all,M)
    all=all[all$M<10000,]
    tapply(all$'% methylated', cut(all$M, seq(-10000,10000, by=100)), mean)->mean_all
    test<-as.data.frame(mean_all)
    test$mean_all[is.na(test$mean_all)]<-0
    test<-unlist(test)
    min<-min(test)
    
    M=matrix(((all_cold$B_start+all_cold$B_end)/2)- all_cold$A_start,ncol=1)
    all_cold=cbind(all_cold,M)
    all_cold=all_cold[all_cold$M<10000,]
    tapply(all_cold$'% methylated', cut(all_cold$M, seq(-10000,10000, by=100)), mean)->mean_all_cold
    test2<-as.data.frame(mean_all_cold)
    test2$mean_all_cold[is.na(test2$mean_all_cold)]<-0
    test2<-unlist(test2)
    max<-max(test2)
    plot(smooth.spline(seq(-10000,9999,100),test,spar = 0.5),col="red",lwd=1,type="l",ylim=c(min,max),main="% CHG Methylation at Eli's Hotspots",xlab ="Distance from midpoint", ylab="Methylation Level (%)")
    lines(smooth.spline(seq(-10000,9999,100),test2,spar = 0.5),col="blue",lwd=1,type="l",main="% CHG Methylation at Eli's Hotspots",xlab ="Distance from midpoint", ylab="Methylation Level (%)")
    
  }
  metaplot_2kb<-function(all, all_cold){
    all<-all[order(all$B_chrom,all$A_start),]
    all_cold<-all_cold[order(all_cold$B_chrom,all_cold$A_start),]
    min<-min(all_cold$'% methylated')
    max<-max(all$'% methylated')
    #,ylim=c(min,max)
    plot(smooth.spline(all$A_start/1000000,all$'% methylated',spar = 1),col="red",lwd=2,type="l",ylim=c(min,max),main="% CHG Methylation at Indels Intersecting with Predicted Crossovers",xlab ="Distance from chromosome start", ylab="Methylation Level (%)")
    lines(smooth.spline(all_cold$A_start/1000000,all_cold$'% methylated',spar = 1),col="blue",lwd=2,type="l",main="% CHG Methylation at Indels Intersecting with Predicted Crossovers",xlab ="Distance from chromosome start", ylab="Methylation Level (%)")
    par(mar=c(5, 4, 4, 10), xpd=TRUE)
    legend("topright", inset=c(-.35, 0), c("'Hot' Indels", "'Cold' Indels"), title="Dataset",
           col=c("brown1", "royalblue3"), lty=1,lwd=2,cex = .8)
    
  }

  metaplot_2kb(CHG_methylation_DSB_hot, CHG_methylation_DSB)
  metaplot_2kb(CHG_methylation_hotspot_hot, CHG_methylation_hotspot)
  metaplot_2kb(CHG_methylation_PCO_hot, CHG_methylation_PCO)
  metaplot_2kb(CHGmethylation_PredCO_hot, CHGmethylation_PredCO)
  
  metaplot_10kb(CHGmethylation_PredCO_hot, CHGmethylation_PredCO)
  metaplot_10kb(CHG_methylation_RCO_hot, CHG_methylation_RCO)
  metaplot_10kb(CHG_methylation_PredDSB_hot, CHG_methylation_PredDSB)
  metaplot_20kb(CHGmethylation_NAM_hot, CHGmethylation_NAM)
  metaplot_20kb(CHG_methylation_hotspot_hot, CHG_methylation_hotspot)
  
  #Summary Table of % Methylation of Indels 
  test <- as.data.frame(matrix(nrow=8,ncol=3))
  colnames(test)<-c("Mean % CHG Methylation","Mean % CHH Methylation","Mean % CpG Methylation")
  rownames(test)<-c("DSBs","Random (DSB)","Eli's Hotspots","Random (Eli)","Penny's Crossovers","Random (PCO)","Predicted Crossovers","Random (PredCO)")
  
  test$`Mean % CHG Methylation`<-c(mean(CHG_methylation_DSB_hot$'% methylated'),mean(CHG_methylation_DSB$'% methylated'),mean(CHG_methylation_hotspot_hot$'% methylated'),mean(CHG_methylation_hotspot$'% methylated'),mean(CHG_methylation_PCO_hot$'% methylated'),mean(CHG_methylation_PCO$'% methylated'),mean(CHGmethylation_PredCO_hot$'% methylated'),mean(CHGmethylation_PredCO$'% methylated'))
  
  test$`Mean % CHH Methylation`<-c(mean( CHH_methylation_DSB_hot$'% methylated'),mean( CHH_methylation_DSB$'% methylated'),mean( CHH_methylation_hotspot_hot$'% methylated'),mean( CHH_methylation_hotspot$'% methylated'),mean( CHH_methylation_PCO_hot$'% methylated'),mean( CHH_methylation_PCO$'% methylated'),mean( CHHmethylation_PredCO_hot$'% methylated'),mean( CHHmethylation_PredCO$'% methylated'))
  
  test$`Mean % CpG Methylation`<-c(mean(CpG_methylation_DSB_hot$'% methylated'),mean(CpG_methylation_DSB$'% methylated'),mean(CpG_methylation_hotspot_hot$'% methylated'),mean(CpG_methylation_hotspot$'% methylated'),mean(CpG_methylation_PCO_hot$'% methylated'),mean(CpG_methylation_PCO$'% methylated'),mean(CpGmethylation_PredCO_hot$'% methylated'),mean(CpGmethylation_PredCO$'% methylated'))
  
  #modifying tables for analysis
  {
  CHG_methylation_DSB_hot$type<-"hot"
  CHG_methylation_DSB_hot$dataset<-"DSB"
  CHG_methylation_DSB$type<-"cold"
  CHG_methylation_DSB$dataset<-"DSB"
  
  CHG_methylation_hotspot_hot$type<-"hot"
  CHG_methylation_hotspot_hot$dataset<-"hotspot"
  CHG_methylation_hotspot$type<-"cold"
  CHG_methylation_hotspot$dataset<-"hotspot"
  
  CHG_methylation_PCO_hot$type<-"hot"
  CHG_methylation_PCO_hot$dataset<-"empCO"
  CHG_methylation_PCO$type<-"cold"
  CHG_methylation_PCO$dataset<-"empCO"

  CHGmethylation_PredCO_hot$type<-"hot"
  CHGmethylation_PredCO_hot$dataset<-"predCO"
  CHGmethylation_PredCO$type<-"cold"
  CHGmethylation_PredCO$dataset<-"predCO"
  
  CHH_methylation_DSB_hot$type<-"hot"
  CHH_methylation_DSB_hot$dataset<-"DSB"
  CHH_methylation_DSB$type<-"cold"
  CHH_methylation_DSB$dataset<-"DSB"
  
  CHH_methylation_hotspot_hot$type<-"hot"
  CHH_methylation_hotspot_hot$dataset<-"hotspot"
  CHH_methylation_hotspot$type<-"cold"
  CHH_methylation_hotspot$dataset<-"hotspot"
  
  CHH_methylation_PCO_hot$type<-"hot"
  CHH_methylation_PCO_hot$dataset<-"empCO"
  CHH_methylation_PCO$type<-"cold"
  CHH_methylation_PCO$dataset<-"empCO"
  
  CHHmethylation_PredCO_hot$type<-"hot"
  CHHmethylation_PredCO_hot$dataset<-"predCO"
  CHHmethylation_PredCO$type<-"cold"
  CHHmethylation_PredCO$dataset<-"predCO"
  
  CpG_methylation_DSB_hot$type<-"hot"
  CpG_methylation_DSB_hot$dataset<-"DSB"
  CpG_methylation_DSB$type<-"cold"
  CpG_methylation_DSB$dataset<-"DSB"
  
  CpG_methylation_hotspot_hot$type<-"hot"
  CpG_methylation_hotspot_hot$dataset<-"hotspot"
  CpG_methylation_hotspot$type<-"cold"
  CpG_methylation_hotspot$dataset<-"hotspot"
  
  CpG_methylation_PCO_hot$type<-"hot"
  CpG_methylation_PCO_hot$dataset<-"empCO"
  CpG_methylation_PCO$type<-"cold"
  CpG_methylation_PCO$dataset<-"empCO"
  
  CpGmethylation_PredCO_hot$type<-"hot"
  CpGmethylation_PredCO_hot$dataset<-"predCO"
  CpGmethylation_PredCO$type<-"cold"
  CpGmethylation_PredCO$dataset<-"predCO"
}
  
  CHG_methylation<-rbind(CHG_methylation_DSB_hot,CHG_methylation_DSB,CHG_methylation_hotspot_hot,CHG_methylation_hotspot,CHG_methylation_PCO_hot,CHG_methylation_PCO,CHGmethylation_PredCO_hot,CHGmethylation_PredCO)
  CpG_methylation<-rbind(CpG_methylation_DSB_hot,CpG_methylation_DSB,CpG_methylation_hotspot_hot,CpG_methylation_hotspot,CpG_methylation_PCO_hot,CpG_methylation_PCO,CpGmethylation_PredCO_hot,CpGmethylation_PredCO)
  CHH_methylation<-rbind(CHH_methylation_DSB_hot,CHH_methylation_DSB,CHH_methylation_hotspot_hot,CHH_methylation_hotspot,CHH_methylation_PCO_hot,CHH_methylation_PCO,CHHmethylation_PredCO_hot,CHHmethylation_PredCO)
  
  data1 <- data.frame(y=CpG_methylation$'% methylated',group=factor(CpG_methylation$type),group2=factor(CpG_methylation$dataset),group3=factor(paste(CpG_methylation$type,CpG_methylation$dataset)))
  ad.test(data1$y) #departure from normality
  fit= lm(y~group2/ factor(group) ,data1)
  two.way <- aov(fit)
  summary(two.way)
  ggplot(data1, aes(x=factor(group2), y=y, fill=group)) +
    geom_boxplot()
  ggboxplot(data1, x = "group2", y = "y",
            color = "group", palette = "jco")+
    stat_compare_means(method = "anova")
  data1$group3<-paste(data1$group,data1$group2)
  pairwise.t.test(data1$y,data1$group3, p.adjust.method="bonferroni")
  
  tukey.test <- TukeyHSD(two.way)
  tukey.test
  plot(tukey.test)
  #ANOVA test indicated significant diff (F=0.15, p=0.861)
  
  library(car)
  leveneTest(fit)
  #Levene's test indicated equal variances (F = 2.5291, p = 0.03962)
  
 
  

  
  test_2 <- as.data.frame(matrix(nrow=12,ncol=3))
  colnames(test_2)<-c("% CHG Methylation","% CHH Methylation","% CpG Methylation")
  rownames(test_2)<-c("NAM CO","Random (NAM)","DSBs","Random (DSB)","Eli's Hotspots","Random (Eli)","Predicted DSB","Random (PDSB)","Ruth's Hotspots","Random (RCO)","Predicted CO","Random (PCO)")
  
  test_2$`% CHG Methylation`<-c(median(CHGmethylation_NAM_hot$'% methylated'),median(CHGmethylation_NAM$'% methylated'),median(CHG_methylation_DSB_hot$'% methylated'),median(CHG_methylation_DSB$'% methylated'),median(CHG_methylation_hotspot_hot$'% methylated'),median(CHG_methylation_hotspot$'% methylated'),median(CHG_methylation_PredDSB_hot$'% methylated'),median(CHG_methylation_PredDSB$'% methylated'),median(CHG_methylation_RCO_hot$'% methylated'),median(CHG_methylation_RCO$'% methylated'),median(CHGmethylation_PredCO_hot$'% methylated'),median(CHGmethylation_PredCO$'% methylated'))
  
  test_2$`% CHH Methylation`<-c(median(CHHmethylation_NAM_hot$'% methylated'),median(CHHmethylation_NAM$'% methylated'),median(CHH_methylation_DSB_hot$'% methylated'),median(CHH_methylation_DSB$'% methylated'),median(CHH_methylation_hotspot_hot$'% methylated'),median(CHH_methylation_hotspot$'% methylated'),median(CHH_methylation_PredDSB_hot$'% methylated'),median(CHH_methylation_PredDSB$'% methylated'),median(CHH_methylation_RCO_hot$'% methylated'),median(CHH_methylation_RCO$'% methylated'),median(CHHmethylation_PredCO_hot$'% methylated'),median(CHHmethylation_PredCO$'% methylated'))
  
  test_2$`% CpG Methylation`<-c(median(CpGmethylation_NAM_hot$'% methylated'),median(CpGmethylation_NAM$'% methylated'),median(CpG_methylation_DSB_hot$'% methylated'),median(CpG_methylation_DSB$'% methylated'),median(CpG_methylation_hotspot_hot$'% methylated'),median(CpG_methylation_hotspot$'% methylated'),median(CpG_methylation_PredDSB_hot$'% methylated'),median(CpG_methylation_PredDSB$'% methylated'),median(CpG_methylation_RCO_hot$'% methylated'),median(CpG_methylation_RCO$'% methylated'),median(CpGmethylation_PredCO_hot$'% methylated'),median(CpGmethylation_PredCO$'% methylated'))
  
  write.table(test,"methylation_indels.txt", row.names = T,quote=FALSE, sep = '\t')
  write.table(test_2,"methylation_median_indels.txt", row.names = T,quote=FALSE, sep = '\t')
  
}

##stat 6: cold data
{
  #modifying tables for analysis
  {
    CHG_methylation_DSB_hot$type<-"hot"
    CHG_methylation_DSB_hot$dataset<-"DSB"
    CHG_methylation_DSB$type<-"cold"
    CHG_methylation_DSB$dataset<-"DSB"
    
    CHG_methylation_hotspot_hot$type<-"hot"
    CHG_methylation_hotspot_hot$dataset<-"hotspot"
    CHG_methylation_hotspot$type<-"cold"
    CHG_methylation_hotspot$dataset<-"hotspot"
    
    CHG_methylation_PCO_hot$type<-"hot"
    CHG_methylation_PCO_hot$dataset<-"empCO"
    CHG_methylation_PCO$type<-"cold"
    CHG_methylation_PCO$dataset<-"empCO"
    
    CHGmethylation_PredCO_hot$type<-"hot"
    CHGmethylation_PredCO_hot$dataset<-"predCO"
    CHGmethylation_PredCO$type<-"cold"
    CHGmethylation_PredCO$dataset<-"predCO"
    
    CHH_methylation_DSB_hot$type<-"hot"
    CHH_methylation_DSB_hot$dataset<-"DSB"
    CHH_methylation_DSB$type<-"cold"
    CHH_methylation_DSB$dataset<-"DSB"
    
    CHH_methylation_hotspot_hot$type<-"hot"
    CHH_methylation_hotspot_hot$dataset<-"hotspot"
    CHH_methylation_hotspot$type<-"cold"
    CHH_methylation_hotspot$dataset<-"hotspot"
    
    CHH_methylation_PCO_hot$type<-"hot"
    CHH_methylation_PCO_hot$dataset<-"empCO"
    CHH_methylation_PCO$type<-"cold"
    CHH_methylation_PCO$dataset<-"empCO"
    
    CHHmethylation_PredCO_hot$type<-"hot"
    CHHmethylation_PredCO_hot$dataset<-"predCO"
    CHHmethylation_PredCO$type<-"cold"
    CHHmethylation_PredCO$dataset<-"predCO"
    
    CpG_methylation_DSB_hot$type<-"hot"
    CpG_methylation_DSB_hot$dataset<-"DSB"
    CpG_methylation_DSB$type<-"cold"
    CpG_methylation_DSB$dataset<-"DSB"
    
    CpG_methylation_hotspot_hot$type<-"hot"
    CpG_methylation_hotspot_hot$dataset<-"hotspot"
    CpG_methylation_hotspot$type<-"cold"
    CpG_methylation_hotspot$dataset<-"hotspot"
    
    CpG_methylation_PCO_hot$type<-"hot"
    CpG_methylation_PCO_hot$dataset<-"empCO"
    CpG_methylation_PCO$type<-"cold"
    CpG_methylation_PCO$dataset<-"empCO"
    
    CpGmethylation_PredCO_hot$type<-"hot"
    CpGmethylation_PredCO_hot$dataset<-"predCO"
    CpGmethylation_PredCO$type<-"cold"
    CpGmethylation_PredCO$dataset<-"predCO"
  }
  
  CHG_methylation<-rbind(CHG_methylation_DSB[-c(8:9)],CHG_methylation_hotspot,CHG_methylation_PCO[-c(9:10)],CHGmethylation_PredCO,CHG_methylation_DSB_hot,CHG_methylation_hotspot_hot,CHG_methylation_PCO_hot,CHGmethylation_PredCO_hot)
  CpG_methylation<-rbind(CpG_methylation_DSB[-c(8:9)],CpG_methylation_hotspot,CpG_methylation_PCO[-c(9:10)],CpGmethylation_PredCO,CpG_methylation_DSB_hot,CpG_methylation_hotspot_hot,CpG_methylation_PCO_hot,CpGmethylation_PredCO_hot)
  CHH_methylation<-rbind(CHH_methylation_DSB[-c(8:9)],CHH_methylation_hotspot,CHH_methylation_PCO[-c(9:10)],CHHmethylation_PredCO,CHH_methylation_DSB_hot,CHH_methylation_hotspot_hot,CHH_methylation_PCO_hot,CHHmethylation_PredCO_hot)
  
  data1 <- data.frame(y=CHH_methylation$'% methylated',group=factor(CHH_methylation$type),group2=factor(CHH_methylation$dataset),group3=factor(paste(CHH_methylation$type,CHH_methylation$dataset)))
  ad.test(data1$y) #departure from normality
  fit= lm(y~group2/ factor(group) ,data1)
  two.way <- aov(fit)
  summary(two.way)
  ggplot(data1, aes(x=factor(group2), y=y, fill=group)) +
    geom_boxplot()
  ggboxplot(data1, x = "group2", y = "y",
            color = "group", palette = "jco")+
    stat_compare_means(method = "anova")
  data1$group3<-paste(data1$group,data1$group2)
  pairwise.t.test(data1$y,data1$group3, p.adjust.method="bonferroni")
  
  tukey.test <- TukeyHSD(two.way)
  tukey.test
  plot(tukey.test)
  #ANOVA test indicated significant diff (F=0.15, p=0.861)
  
  library(car)
  leveneTest(fit)
  #Levene's test indicated equal variances (F = 2.5291, p = 0.03962)
  
}