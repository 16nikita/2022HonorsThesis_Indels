#!/bin/sh
echo running bedtools intersect

bedtools intersect -a indels_all.bed -b penny_CO_v4.bed -wo > intersect_a_PCO.txt
bedtools intersect -a indels_10bp_all.bed -b penny_CO_v4.bed -wo > intersect_10bp_a_PCO.txt
bedtools intersect -a indels_20bp_all.bed -b penny_CO_v4.bed -wo > intersect_20bp_a_PCO.txt
bedtools intersect -a indels_30bp_all.bed -b penny_CO_v4.bed -wo > intersect_30bp_a_PCO.txt
bedtools intersect -a indels_40bp_all.bed -b penny_CO_v4.bed -wo > intersect_40bp_a_PCO.txt
bedtools intersect -a indels_50bp_all.bed -b penny_CO_v4.bed -wo > intersect_50bp_a_PCO.txt

bedtools intersect -a indels_all.bed -b PCO_2kb_up.bed -wo > intersect_2kbup_a_PCO.txt
bedtools intersect -a indels_10bp_all.bed -b PCO_2kb_up.bed -wo > intersect_10bp_2kbup_a_PCO.txt
bedtools intersect -a indels_20bp_all.bed -b PCO_2kb_up.bed -wo > intersect_20bp_2kbup_a_PCO.txt
bedtools intersect -a indels_30bp_all.bed -b PCO_2kb_up.bed -wo > intersect_30bp_2kbup_a_PCO.txt
bedtools intersect -a indels_40bp_all.bed -b PCO_2kb_up.bed -wo > intersect_40bp_2kbup_a_PCO.txt
bedtools intersect -a indels_50bp_all.bed -b PCO_2kb_up.bed -wo > intersect_50bp_2kbup_a_PCO.txt

bedtools intersect -a indels_all.bed -b PCO_2kb_down.bed -wo > intersect_2kbdown_a_PCO.txt
bedtools intersect -a indels_10bp_all.bed -b PCO_2kb_down.bed -wo > intersect_10bp_2kbdown_a_PCO.txt
bedtools intersect -a indels_20bp_all.bed -b PCO_2kb_down.bed -wo > intersect_20bp_2kbdown_a_PCO.txt
bedtools intersect -a indels_30bp_all.bed -b PCO_2kb_down.bed -wo > intersect_30bp_2kbdown_a_PCO.txt
bedtools intersect -a indels_40bp_all.bed -b PCO_2kb_down.bed -wo > intersect_40bp_2kbdown_a_PCO.txt
bedtools intersect -a indels_50bp_all.bed -b PCO_2kb_down.bed -wo > intersect_50bp_2kbdown_a_PCO.txt

bedtools intersect -a indels_all.bed -b PCO_random.bed -wo > intersect_a_randPCO.txt
bedtools intersect -a indels_10bp_all.bed -b PCO_random.bed -wo > intersect_10bp_a_randPCO.txt
bedtools intersect -a indels_20bp_all.bed -b PCO_random.bed -wo > intersect_20bp_a_randPCO.txt
bedtools intersect -a indels_30bp_all.bed -b PCO_random.bed -wo > intersect_30bp_a_randPCO.txt
bedtools intersect -a indels_40bp_all.bed -b PCO_random.bed -wo > intersect_40bp_a_randPCO.txt
bedtools intersect -a indels_50bp_all.bed -b PCO_random.bed -wo > intersect_50bp_a_randPCO.txt

bedtools intersect -a indels_all.bed -b PCO_cold_2kb_down.bed -wo > intersect_2kbdown_a_randPCO.txt
bedtools intersect -a indels_10bp_all.bed -b PCO_cold_2kb_down.bed -wo > intersect_10bp_2kbdown_a_randPCO.txt
bedtools intersect -a indels_20bp_all.bed -b PCO_cold_2kb_down.bed -wo > intersect_20bp_2kbdown_a_randPCO.txt
bedtools intersect -a indels_30bp_all.bed -b PCO_cold_2kb_down.bed -wo > intersect_30bp_2kbdown_a_randPCO.txt
bedtools intersect -a indels_40bp_all.bed -b PCO_cold_2kb_down.bed -wo > intersect_40bp_2kbdown_a_randPCO.txt
bedtools intersect -a indels_50bp_all.bed -b PCO_cold_2kb_down.bed -wo > intersect_50bp_2kbdown_a_randPCO.txt

bedtools intersect -a indels_all.bed -b PCO_cold_2kb_up.bed -wo > intersect_2kbup_a_randPCO.txt
bedtools intersect -a indels_10bp_all.bed -b PCO_cold_2kb_up.bed -wo > intersect_10bp_2kbup_a_randPCO.txt
bedtools intersect -a indels_20bp_all.bed -b PCO_cold_2kb_up.bed -wo > intersect_20bp_2kbup_a_randPCO.txt
bedtools intersect -a indels_30bp_all.bed -b PCO_cold_2kb_up.bed -wo > intersect_30bp_2kbup_a_randPCO.txt
bedtools intersect -a indels_40bp_all.bed -b PCO_cold_2kb_up.bed -wo > intersect_40bp_2kbup_a_randPCO.txt
bedtools intersect -a indels_50bp_all.bed -b PCO_cold_2kb_up.bed -wo > intersect_50bp_2kbup_a_randPCO.txt
