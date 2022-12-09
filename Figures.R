
#chromoMap
#indel datasets- small, medium, large
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


#mnase and indel density:
{
small_mnase<-read.delim("intersect_open_pred.txt",sep='\t', header = FALSE)[-c(7)]
med_mnase<-read.delim("intersect_med_open_pred.txt",sep='\t', header = FALSE)[-c(7)]
large_mnase<-read.delim("intersect_large_open_pred.txt",sep='\t', header = FALSE)[-c(7)]

percentage_overlap<-function(intersect){
  colnames(intersect) <- c('A_chrom', 'A_start', 'A_end','B_chrom', 'B_start', 'B_end','occupancy','overlap in bp')
  intersect<-intersect %>% distinct(A_start, .keep_all = TRUE)
  duplicated<-as.data.frame(table(intersect$'B_start'))
  intersect$length<-intersect$'B_end'-intersect$'B_start'
  intersect <- intersect %>% distinct(B_start, .keep_all = TRUE)
  intersect<-intersect[order(intersect$B_start),]
  intersect$density<-duplicated$Freq
  intersect$percent<- (intersect$'density'/intersect$length)
  return(intersect)
}

chr1<-percentage_overlap(small_mnase)
chr1<-percentage_overlap(med_mnase)
chr1<-percentage_overlap(large_mnase)
cor.test(log(chr1$percent), log(chr1$occupancy),  method = "spearman", alternative = "greater")
chr1<-chr1[order(chr1$percent),]
values <- loess(log(chr1$occupancy)~ log(chr1$percent))
ggscatter(chr1, x = "percent", y = "occupancy", 
          add = "reg.line",add.params = list(color = "red"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Indel Density (indels/bp)", ylab = "Nucleosome occupancy")
plot(log(chr1$percent), log(chr1$occupancy))
lines(log(chr1$percent),predict(values),
      col = "blue",
      lwd = 2)
}
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
  #annotation_file$indel_per_Mb<-annotation_file$indel_per_Mb
  
  maize_genome<-read.delim("maize.genome",header=FALSE)
  colnames(chrom_file)<-c('chrom','start','end','centromere')
  chrom_file$chrom<-c(1:10)
  cent<-c(121+133,92.2+101.6,87.4+89.4,71.8+93,94.8+118.4,31.8+32.7,44.9+54.1,55.8+80.7,34.1+34.2,37.1+46)
  cent<-cent/2
  chrom_file$centromere<-cent*1000000
  chrom_file$end<-maize_genome$V3
  #save(chrom_file, file = "chrom_file.RData")
  load("chrom_file.RData")
  #install.packages("chromoMap")
  library(chromoMap)
  library('RColorBrewer')
  display.brewer.all()
  chromoMap(ch.files = list(chrom_file), data.files = list(annotation_file), title = "Medium Indel Density Across Each Chromosome",data_based_color_map = T,
            data_type = "numeric",  plot_y_domain = list(c(0,.03)),legend = T,lg_x=100,lg_y= 500,n_win.factor = 1, export.options = T, plots = "bar", plot_color="purple",data_colors = list(brewer.pal(n = 3, name = "PuRd")))
  #without chromosomes
  # chromoMap(ch.files = list(chrom_file), data.files = list(annotation_file), title = "Large Indel Density Across Each Chromosome",data_based_color_map = T,
  #           data_type = "numeric", legend = T,display.chr = c(F), lg_y= 500,n_win.factor = 1, export.options = T, plots = "bar", plot_color="purple",data_colors = list(brewer.pal(n = 3, name = "PuRd")))
  
}



#no data
intersect_40bp_cold_pred<-as.data.frame(matrix(data=0,ncol=7,nrow=1))

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


#calculating INDEL DENSITY:
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

#calculating other feature density:
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
    scale_fill_manual(labels = c("_____",
                                 "Negative Control Region"),
                      values = c("Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("____ Indel Density on _____") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(),width=.2)+
    geom_text(aes(label=posthoc, y=value+sd), position=position_dodge(.5), vjust=-.5) +
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(min(full_data1$value - full_data1$sd),max(full_data1$value + full_data1$sd))) +
    annotate("text", c(1.5,3.5,5.5,7.5,9.5), y = -.0005, label = c("1-10 bp","11-20 bp","21-30 bp","31-40 bp","41-50 bp")) +
    #annotate("text", c(1.5,3.5,5.5,7.5,9.5), y = -.000015, label = c("100-150 bp","151-200 bp","201-300 bp","301-400 bp","401-500 bp")) +
    #annotate("text", c(1.5,3.5,5.5,7.5,9.5), y = -.00003, label = c("500bp-10 kb","10-20 kb","20-30 kb","30-40 kb","40-50 kb")) +
    annotate("text", x = 5.5, y = -.0008,
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
                                 "____",
                                 "Negative Control"),
                      values = c("2kb downstream hot"= "darkgoldenrod1",
                                 "2kb downstream cold" = "deepskyblue",
                                 "2kb upstream hot" = "darkorange1",
                                 "2kb upstream cold" = "dodgerblue2",
                                 "Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("____ Indel Density on ____") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(.5),width=.2)+
    geom_text(aes(label=posthoc, y=value+sd), position=position_dodge(.5), vjust=-.5) +
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(0,max(full_data2$value + full_data2$sd))) +
    annotate("text", x = 1:2, y = -.0005,
             label = rep(c("Negative Control Region", "____ Region"), 1)) +
    annotate("text", x = 1.5, y = -.001,
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

#BARGRAPHS - for each: change TITLE, INDEL SIZES, LEGEND, AXIS TITLE POSITION
{
  #create dataframes:
  {
    #DATAFRAME OF CO PREDICTED
    type <- c(rep(c("all indels") ,3),rep(c("1-10bp") ,3),rep(c("11-20bp") ,3),rep(c("21-3bp") ,3),rep(c("31-40bp") ,3),rep(c("41-50bp") ,3))
    intersect<-c(c("Interval"),c("2kb upstream") ,c("2kb downstream"))
    #intersect<-c(c("Interval"),c("2kb upstream") ,c("2kb downstream"),c("Interval"),c("2kb upstream") ,c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10),rep(c("Interval"),10),rep(c("2kb upstream") ,10),rep(c("2kb downstream"),10))
    value <- c(mean(intersect_overlap),mean(intersect_2kbup_overlap),mean(intersect_2kbdown_overlap),mean(intersect_10bp_overlap),mean(intersect_10bp_2kbup_overlap),mean(intersect_10bp_2kbdown_overlap),mean(intersect_20bp_overlap),mean(intersect_20bp_2kbup_overlap),mean(intersect_20bp_2kbdown_overlap),mean(intersect_30bp_overlap),mean(intersect_30bp_2kbup_overlap),mean(intersect_30bp_2kbdown_overlap),mean(intersect_40bp_overlap),mean(intersect_40bp_2kbup_overlap),mean(intersect_40bp_2kbdown_overlap),mean(intersect_50bp_overlap),mean(intersect_50bp_2kbup_overlap),mean(intersect_50bp_2kbdown_overlap))
    sd <- c(sd(intersect_overlap),sd(intersect_2kbup_overlap),sd(intersect_2kbdown_overlap),sd(intersect_10bp_overlap),sd(intersect_10bp_2kbup_overlap),sd(intersect_10bp_2kbdown_overlap),sd(intersect_20bp_overlap),sd(intersect_20bp_2kbup_overlap),sd(intersect_20bp_2kbdown_overlap),sd(intersect_30bp_overlap),sd(intersect_30bp_2kbup_overlap),sd(intersect_30bp_2kbdown_overlap),sd(intersect_40bp_overlap),sd(intersect_40bp_2kbup_overlap),sd(intersect_40bp_2kbdown_overlap),sd(intersect_50bp_overlap),sd(intersect_50bp_2kbup_overlap),sd(intersect_50bp_2kbdown_overlap))
    data <- data.frame(type,value,intersect,sd)
    data$dataset<-"hot"
    
    #DATAFRAME OF CO DESERT
    value <- c(mean(intersect_cold_overlap),mean(intersect_2kbup_cold_overlap),mean(intersect_2kbdown_cold_overlap),mean(intersect_10bp_cold_overlap),mean(intersect_10bp_2kbup_cold_overlap),mean(intersect_10bp_2kbdown_cold_overlap),mean(intersect_20bp_cold_overlap),mean(intersect_20bp_2kbup_cold_overlap),mean(intersect_20bp_2kbdown_cold_overlap),mean(intersect_30bp_cold_overlap),mean(intersect_30bp_2kbup_cold_overlap),mean(intersect_30bp_2kbdown_cold_overlap),mean(intersect_40bp_cold_overlap),mean(intersect_40bp_2kbup_cold_overlap),mean(intersect_40bp_2kbdown_cold_overlap),mean(intersect_50bp_cold_overlap),mean(intersect_50bp_2kbup_cold_overlap),mean(intersect_50bp_2kbdown_cold_overlap))
    sd <- c(sd(intersect_cold_overlap),sd(intersect_2kbup_cold_overlap),sd(intersect_2kbdown_cold_overlap),sd(intersect_10bp_cold_overlap),sd(intersect_10bp_2kbup_cold_overlap),sd(intersect_10bp_2kbdown_cold_overlap),sd(intersect_20bp_cold_overlap),sd(intersect_20bp_2kbup_cold_overlap),sd(intersect_20bp_2kbdown_cold_overlap),sd(intersect_30bp_cold_overlap),sd(intersect_30bp_2kbup_cold_overlap),sd(intersect_30bp_2kbdown_cold_overlap),sd(intersect_40bp_cold_overlap),sd(intersect_40bp_2kbup_cold_overlap),sd(intersect_40bp_2kbdown_cold_overlap),sd(intersect_50bp_cold_overlap),sd(intersect_50bp_2kbup_cold_overlap),sd(intersect_50bp_2kbdown_cold_overlap))
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
  #plot 1: indel density 
  perc<-ggplot(full_data1, aes(fill=color,x=levels,y=value)) +
    geom_bar(stat="identity",position="dodge") +
    scale_fill_manual(labels = c("____",
                                 "Coldspot"),
                      values = c("Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("Small Indel Density on ____") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(width=1))+
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(min(full_data1$value - full_data1$sd),max(full_data1$value + full_data1$sd))) +
    annotate("text", x = 1:12, y = -.0002,
             label = rep(c("Cold", "Hot"), 6)) +
    annotate("text", c(1.5,3.5,5.5,7.5,9.5,11.5), y = -.0008, label = c("1-10 bp","11-20 bp","21-30 bp","31-40 bp","41-50 bp","All Indels")) +
    #annotate("text", c(1.5,3.5,5.5,7.5,9.5), y = -.000015, label = c("100-150 bp","151-200 bp","201-300 bp","301-400 bp","401-500 bp")) +
    #annotate("text", c(1.5,3.5,5.5,7.5,9.5), y = -.00003, label = c("500bp-10 kb","10-20 kb","20-30 kb","30-40 kb","40-50 kb")) +
    annotate("text", x = 6.5, y = -.0003,
             label = c("Dataset (cold/hot) and Indel Size"),fontface =2) +
    ylab("Indel Density (indels/bp)")+
    guides(fill=guide_legend("Type of intersection"))+
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
  perc2<-ggplot(full_data2, aes(fill=color,x=levels,y=value)) +
    geom_bar(stat="identity",position="dodge", width = 1) +
    scale_fill_manual(labels = c("2kb downstream", "2kb downstream",
                                 "2kb upstream","2kb upstream",
                                 "_____",
                                 "Coldspot"),
                      values = c("2kb downstream hot"= "darkgoldenrod1",
                                 "2kb downstream cold" = "deepskyblue",
                                 "2kb upstream hot" = "darkorange1",
                                 "2kb upstream cold" = "dodgerblue2",
                                 "Interval hot" = "brown1",
                                 "Interval cold" = "royalblue3"))+
    ggtitle("Small Indel Density on _____") +
    geom_errorbar(aes(ymax = value + sd, ymin = value - sd), position = position_dodge(width=1))+
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    coord_cartesian(ylim = c(0,max(full_data2$value + full_data2$sd))) +
    annotate("text", x = 1:2, y = -.00008,
             label = rep(c("Cold", "Hot"), 1)) +
    #annotate("text", c(1.5,3.5,5.5,7.5,9.5,11.5), y = -.00013, label = c("1-10 bp","11-20 bp","21-30 bp","31-40 bp","41-50 bp","All Indels")) +
    annotate("text", x = 1.5, y = -.00013,
             label = c("Dataset"),fontface =2) +
    ylab("Indel Density (indels/bp)")+
    guides(fill=guide_legend("Type of intersection"))+
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
    #Levene's test indicated equal variances (F = 0.0116, p = .9885)
    one.way <- aov(fit)
    summary(one.way)
    #ANOVA test indicated no significant diff (F=0.15, p=0.861)
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
    ggtitle("Small Indel Density on ____ in Gene/TSS/TTS") +
    geom_errorbar(aes(ymax = value + se, ymin = value - se), position = position_dodge(),width=.2)+
    theme(plot.title = element_text(face="bold", size=15, hjust=.5)) +
    xlab("_____ Feature Location")+
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

#TESTING: exon/intron
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
  
  exon_analysis<- ggplot(data, aes(fill=size,x=type,y=value)) +
    geom_bar(stat="identity",position="dodge", width = .5)+
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
    
    ggplot(data1, aes(x=factor(group2), y=y, fill=group)) +
      geom_boxplot()
    ggboxplot(data1, x = "group2", y = "y",
              color = "group", palette = "jco")+
      stat_compare_means(method = "anova")
    
    
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