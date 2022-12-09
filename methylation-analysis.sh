#!/bin/sh
echo running bedtools intersect
bedtools intersect -a CHG_2.bed -b DSB_V4.bed -wo > CHGmethylation_DSB_hot.txt
bedtools intersect -a CHG_2.bed -b CO_predicted_3.bed -wo > CHGmethylation_PredCO_hot.txt
bedtools intersect -a CHG_2.bed -b hotspot_2.bed -wo > CHGmethylation_hotspot_hot.txt
bedtools intersect -a CHG_2.bed -b penny_CO_v4.bed -wo > CHGmethylation_PCO_hot.txt

bedtools intersect -a CHH_2.bed -b DSB_V4.bed -wo > CHHmethylation_DSB_hot.txt
bedtools intersect -a CHH_2.bed -b CO_predicted_3.bed -wo > CHHmethylation_PredCO_hot.txt
bedtools intersect -a CHH_2.bed -b hotspot_2.bed -wo > CHHmethylation_hotspot_hot.txt
bedtools intersect -a CHH_2.bed -b penny_CO_v4.bed -wo > CHHmethylation_PCO_hot.txt

bedtools intersect -a CpG_2.bed -b DSB_V4.bed -wo > CpGmethylation_DSB_hot.txt
bedtools intersect -a CpG_2.bed -b CO_predicted_3.bed -wo > CpGmethylation_PredCO_hot.txt
bedtools intersect -a CpG_2.bed -b hotspot_2.bed -wo > CpGmethylation_hotspot_hot.txt
bedtools intersect -a CpG_2.bed -b penny_CO_v4.bed -wo > CpGmethylation_PCO_hot.txt

bedtools intersect -a CHH_2.bed -b no_indel_DSB.bed -wo > CHHmethylation_DSB.txt
bedtools intersect -a CHH_2.bed -b no_indel_predCO.bed -wo > CHHmethylation_PredCO.txt
bedtools intersect -a CHH_2.bed -b no_indel_hotspot.bed -wo > CHHmethylation_hotspot.txt
bedtools intersect -a CHH_2.bed -b no_indel_empCO.bed -wo > CHHmethylation_PCO.txt

bedtools intersect -a CHG_2.bed -b no_indel_DSB.bed -wo > CHGmethylation_DSB.txt
bedtools intersect -a CHG_2.bed -b no_indel_predCO.bed -wo > CHGmethylation_PredCO.txt
bedtools intersect -a CHG_2.bed -b no_indel_hotspot.bed -wo > CHGmethylation_hotspot.txt
bedtools intersect -a CHG_2.bed -b no_indel_empCO.bed -wo > CHGmethylation_PCO.txt

bedtools intersect -a CpG_2.bed -b no_indel_DSB.bed -wo > CpGmethylation_DSB.txt
bedtools intersect -a CpG_2.bed -b no_indel_predCO.bed -wo > CpGmethylation_PredCO.txt
bedtools intersect -a CpG_2.bed -b no_indel_hotspot.bed -wo > CpGmethylation_hotspot.txt
bedtools intersect -a CpG_2.bed -b no_indel_empCO.bed -wo > CpGmethylation_PCO.txt
