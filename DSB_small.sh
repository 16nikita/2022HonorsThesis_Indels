#!/bin/sh
echo running bedtools intersect

bedtools intersect -a indels_all.bed -b DSB_V4.bed -wo > intersect_a_DSB.txt
bedtools intersect -a indels_10bp_all.bed -b DSB_V4.bed -wo > intersect_10bp_a_DSB.txt
bedtools intersect -a indels_20bp_all.bed -b DSB_V4.bed -wo > intersect_20bp_a_DSB.txt
bedtools intersect -a indels_30bp_all.bed -b DSB_V4.bed -wo > intersect_30bp_a_DSB.txt
bedtools intersect -a indels_40bp_all.bed -b DSB_V4.bed -wo > intersect_40bp_a_DSB.txt
bedtools intersect -a indels_50bp_all.bed -b DSB_V4.bed -wo > intersect_50bp_a_DSB.txt

bedtools intersect -a indels_all.bed -b DSB_2kb_up.bed -wo > intersect_2kbup_a_DSB.txt
bedtools intersect -a indels_10bp_all.bed -b DSB_2kb_up.bed -wo > intersect_10bp_2kbup_a_DSB.txt
bedtools intersect -a indels_20bp_all.bed -b DSB_2kb_up.bed -wo > intersect_20bp_2kbup_a_DSB.txt
bedtools intersect -a indels_30bp_all.bed -b DSB_2kb_up.bed -wo > intersect_30bp_2kbup_a_DSB.txt
bedtools intersect -a indels_40bp_all.bed -b DSB_2kb_up.bed -wo > intersect_40bp_2kbup_a_DSB.txt
bedtools intersect -a indels_50bp_all.bed -b DSB_2kb_up.bed -wo > intersect_50bp_2kbup_a_DSB.txt

bedtools intersect -a indels_all.bed -b DSB_2kb_down.bed -wo > intersect_2kbdown_a_DSB.txt
bedtools intersect -a indels_10bp_all.bed -b DSB_2kb_down.bed -wo > intersect_10bp_2kbdown_a_DSB.txt
bedtools intersect -a indels_20bp_all.bed -b DSB_2kb_down.bed -wo > intersect_20bp_2kbdown_a_DSB.txt
bedtools intersect -a indels_30bp_all.bed -b DSB_2kb_down.bed -wo > intersect_30bp_2kbdown_a_DSB.txt
bedtools intersect -a indels_40bp_all.bed -b DSB_2kb_down.bed -wo > intersect_40bp_2kbdown_a_DSB.txt
bedtools intersect -a indels_50bp_all.bed -b DSB_2kb_down.bed -wo > intersect_50bp_2kbdown_a_DSB.txt

bedtools intersect -a indels_all.bed -b DSB_coldspotsv4.bed -wo > intersect_a_randDSB.txt
bedtools intersect -a indels_10bp_all.bed -b DSB_coldspotsv4.bed -wo > intersect_10bp_a_randDSB.txt
bedtools intersect -a indels_20bp_all.bed -b DSB_coldspotsv4.bed -wo > intersect_20bp_a_randDSB.txt
bedtools intersect -a indels_30bp_all.bed -b DSB_coldspotsv4.bed -wo > intersect_30bp_a_randDSB.txt
bedtools intersect -a indels_40bp_all.bed -b DSB_coldspotsv4.bed -wo > intersect_40bp_a_randDSB.txt
bedtools intersect -a indels_50bp_all.bed -b DSB_coldspotsv4.bed -wo > intersect_50bp_a_randDSB.txt

bedtools intersect -a indels_all.bed -b DSB_cold_2kb_down.bed -wo > intersect_2kbdown_a_randDSB.txt
bedtools intersect -a indels_10bp_all.bed -b DSB_cold_2kb_down.bed -wo > intersect_10bp_2kbdown_a_randDSB.txt
bedtools intersect -a indels_20bp_all.bed -b DSB_cold_2kb_down.bed -wo > intersect_20bp_2kbdown_a_randDSB.txt
bedtools intersect -a indels_30bp_all.bed -b DSB_cold_2kb_down.bed -wo > intersect_30bp_2kbdown_a_randDSB.txt
bedtools intersect -a indels_40bp_all.bed -b DSB_cold_2kb_down.bed -wo > intersect_40bp_2kbdown_a_randDSB.txt
bedtools intersect -a indels_50bp_all.bed -b DSB_cold_2kb_down.bed -wo > intersect_50bp_2kbdown_a_randDSB.txt

bedtools intersect -a indels_all.bed -b DSB_cold_2kb_up.bed -wo > intersect_2kbup_a_randDSB.txt
bedtools intersect -a indels_10bp_all.bed -b DSB_cold_2kb_up.bed -wo > intersect_10bp_2kbup_a_randDSB.txt
bedtools intersect -a indels_20bp_all.bed -b DSB_cold_2kb_up.bed -wo > intersect_20bp_2kbup_a_randDSB.txt
bedtools intersect -a indels_30bp_all.bed -b DSB_cold_2kb_up.bed -wo > intersect_30bp_2kbup_a_randDSB.txt
bedtools intersect -a indels_40bp_all.bed -b DSB_cold_2kb_up.bed -wo > intersect_40bp_2kbup_a_randDSB.txt
bedtools intersect -a indels_50bp_all.bed -b DSB_cold_2kb_up.bed -wo > intersect_50bp_2kbup_a_randDSB.txt
