#!/bin/sh
echo running bedtools intersect
#INTERSECTION WITH PREDICTED CROSSOVERS after clustering (small indels, gradient)
bedtools intersect -a indels_all.bed -b CO_predicted_3.bed -wo > intersect_a_pred.txt
bedtools intersect -a indels_10bp_all.bed -b CO_predicted_3.bed -wo > intersect_10bp_a_pred.txt
bedtools intersect -a indels_20bp_all.bed -b CO_predicted_3.bed -wo > intersect_20bp_a_pred.txt
bedtools intersect -a indels_30bp_all.bed -b CO_predicted_3.bed -wo > intersect_30bp_a_pred.txt
bedtools intersect -a indels_40bp_all.bed -b CO_predicted_3.bed -wo > intersect_40bp_a_pred.txt
bedtools intersect -a indels_50bp_all.bed -b CO_predicted_3.bed -wo > intersect_50bp_a_pred.txt

bedtools intersect -a indels_all.bed -b CO_pred_2kb_up_2.bed -wo > intersect_2kbup_a_pred.txt
bedtools intersect -a indels_10bp_all.bed -b CO_pred_2kb_up_2.bed -wo > intersect_10bp_2kbup_a_pred.txt
bedtools intersect -a indels_20bp_all.bed -b CO_pred_2kb_up_2.bed -wo > intersect_20bp_2kbup_a_pred.txt
bedtools intersect -a indels_30bp_all.bed -b CO_pred_2kb_up_2.bed -wo > intersect_30bp_2kbup_a_pred.txt
bedtools intersect -a indels_40bp_all.bed -b CO_pred_2kb_up_2.bed -wo > intersect_40bp_2kbup_a_pred.txt
bedtools intersect -a indels_50bp_all.bed -b CO_pred_2kb_up_2.bed -wo > intersect_50bp_2kbup_a_pred.txt

bedtools intersect -a indels_all.bed -b CO_pred_2kb_down_2.bed -wo > intersect_2kbdown_a_pred.txt
bedtools intersect -a indels_10bp_all.bed -b CO_pred_2kb_down_2.bed -wo > intersect_10bp_2kbdown_a_pred.txt
bedtools intersect -a indels_20bp_all.bed -b CO_pred_2kb_down_2.bed -wo > intersect_20bp_2kbdown_a_pred.txt
bedtools intersect -a indels_30bp_all.bed -b CO_pred_2kb_down_2.bed -wo > intersect_30bp_2kbdown_a_pred.txt
bedtools intersect -a indels_40bp_all.bed -b CO_pred_2kb_down_2.bed -wo > intersect_40bp_2kbdown_a_pred.txt
bedtools intersect -a indels_50bp_all.bed -b CO_pred_2kb_down_2.bed -wo > intersect_50bp_2kbdown_a_pred.txt

bedtools intersect -a indels_all.bed -b CO_pred_coldspots_2.bed -wo > intersect_cold_a_pred.txt
bedtools intersect -a indels_10bp_all.bed -b CO_pred_coldspots_2.bed -wo > intersect_10bp_cold_a_pred.txt
bedtools intersect -a indels_20bp_all.bed -b CO_pred_coldspots_2.bed -wo > intersect_20bp_cold_a_pred.txt
bedtools intersect -a indels_30bp_all.bed -b CO_pred_coldspots_2.bed -wo > intersect_30bp_cold_a_pred.txt
bedtools intersect -a indels_40bp_all.bed -b CO_pred_coldspots_2.bed -wo > intersect_40bp_cold_a_pred.txt
bedtools intersect -a indels_50bp_all.bed -b CO_pred_coldspots_2.bed -wo > intersect_50bp_cold_a_pred.txt

bedtools intersect -a indels_all.bed -b CO_pred_cold_2kb_down_2.bed -wo > intersect_2kbdown_cold_a_pred.txt
bedtools intersect -a indels_10bp_all.bed -b CO_pred_cold_2kb_down_2.bed -wo > intersect_10bp_2kbdown_cold_a_pred.txt
bedtools intersect -a indels_20bp_all.bed -b CO_pred_cold_2kb_down_2.bed -wo > intersect_20bp_2kbdown_cold_a_pred.txt
bedtools intersect -a indels_30bp_all.bed -b CO_pred_cold_2kb_down_2.bed -wo > intersect_30bp_2kbdown_cold_a_pred.txt
bedtools intersect -a indels_40bp_all.bed -b CO_pred_cold_2kb_down_2.bed -wo > intersect_40bp_2kbdown_cold_a_pred.txt
bedtools intersect -a indels_50bp_all.bed -b CO_pred_cold_2kb_down_2.bed -wo > intersect_50bp_2kbdown_cold_a_pred.txt

bedtools intersect -a indels_all.bed -b CO_pred_cold_2kb_up_2.bed -wo > intersect_2kbup_cold_a_pred.txt
bedtools intersect -a indels_10bp_all.bed -b CO_pred_cold_2kb_up_2.bed -wo > intersect_10bp_2kbup_cold_a_pred.txt
bedtools intersect -a indels_20bp_all.bed -b CO_pred_cold_2kb_up_2.bed -wo > intersect_20bp_2kbup_cold_a_pred.txt
bedtools intersect -a indels_30bp_all.bed -b CO_pred_cold_2kb_up_2.bed -wo > intersect_30bp_2kbup_cold_a_pred.txt
bedtools intersect -a indels_40bp_all.bed -b CO_pred_cold_2kb_up_2.bed -wo > intersect_40bp_2kbup_cold_a_pred.txt
bedtools intersect -a indels_50bp_all.bed -b CO_pred_cold_2kb_up_2.bed -wo > intersect_50bp_2kbup_cold_a_pred.txt
