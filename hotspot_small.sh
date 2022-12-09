#!/bin/sh
echo running bedtools intersect
bedtools intersect -a indels_all.bed -b hotspot_2.bed -wo > intersect_sub_hotspot.txt
bedtools intersect -a indels_10bp_all.bed -b hotspot_2.bed -wo > intersect_10bp_sub_hotspot.txt
bedtools intersect -a indels_20bp_all.bed -b hotspot_2.bed -wo > intersect_20bp_sub_hotspot.txt
bedtools intersect -a indels_30bp_all.bed -b hotspot_2.bed -wo > intersect_30bp_sub_hotspot.txt
bedtools intersect -a indels_40bp_all.bed -b hotspot_2.bed -wo > intersect_40bp_sub_hotspot.txt
bedtools intersect -a indels_50bp_all.bed -b hotspot_2.bed -wo > intersect_50bp_sub_hotspot.txt

bedtools intersect -a indels_all.bed -b hotspots_2_2kbup.bed -wo > intersect_2kbup_sub_hotspot.txt
bedtools intersect -a indels_10bp_all.bed -b hotspots_2_2kbup.bed -wo > intersect_10bp_2kbup_sub_hotspot.txt
bedtools intersect -a indels_20bp_all.bed -b hotspots_2_2kbup.bed -wo > intersect_20bp_2kbup_sub_hotspot.txt
bedtools intersect -a indels_30bp_all.bed -b hotspots_2_2kbup.bed -wo > intersect_30bp_2kbup_sub_hotspot.txt
bedtools intersect -a indels_40bp_all.bed -b hotspots_2_2kbup.bed -wo > intersect_40bp_2kbup_sub_hotspot.txt
bedtools intersect -a indels_50bp_all.bed -b hotspots_2_2kbup.bed -wo > intersect_50bp_2kbup_sub_hotspot.txt

bedtools intersect -a indels_all.bed -b hotspots_2_2kbdown.bed -wo > intersect_2kbdown_sub_hotspot.txt
bedtools intersect -a indels_10bp_all.bed -b hotspots_2_2kbdown.bed -wo > intersect_10bp_2kbdown_sub_hotspot.txt
bedtools intersect -a indels_20bp_all.bed -b hotspots_2_2kbdown.bed -wo > intersect_20bp_2kbdown_sub_hotspot.txt
bedtools intersect -a indels_30bp_all.bed -b hotspots_2_2kbdown.bed -wo > intersect_30bp_2kbdown_sub_hotspot.txt
bedtools intersect -a indels_40bp_all.bed -b hotspots_2_2kbdown.bed -wo > intersect_40bp_2kbdown_sub_hotspot.txt
bedtools intersect -a indels_50bp_all.bed -b hotspots_2_2kbdown.bed -wo > intersect_50bp_2kbdown_sub_hotspot.txt

bedtools intersect -a indels_all.bed -b hotspot_2_random.bed -wo > intersect_sub_coldspot.txt
bedtools intersect -a indels_10bp_all.bed -b hotspot_2_random.bed -wo > intersect_10bp_sub_coldspot.txt
bedtools intersect -a indels_20bp_all.bed -b hotspot_2_random.bed -wo > intersect_20bp_sub_coldspot.txt
bedtools intersect -a indels_30bp_all.bed -b hotspot_2_random.bed -wo > intersect_30bp_sub_coldspot.txt
bedtools intersect -a indels_40bp_all.bed -b hotspot_2_random.bed -wo > intersect_40bp_sub_coldspot.txt
bedtools intersect -a indels_50bp_all.bed -b hotspot_2_random.bed -wo > intersect_50bp_sub_coldspot.txt

bedtools intersect -a indels_all.bed -b hotspotRand_2_2kb_down.bed -wo > intersect_2kbdown_sub_coldspot.txt
bedtools intersect -a indels_10bp_all.bed -b hotspotRand_2_2kb_down.bed -wo > intersect_10bp_2kbdown_sub_coldspot.txt
bedtools intersect -a indels_20bp_all.bed -b hotspotRand_2_2kb_down.bed -wo > intersect_20bp_2kbdown_sub_coldspot.txt
bedtools intersect -a indels_30bp_all.bed -b hotspotRand_2_2kb_down.bed -wo > intersect_30bp_2kbdown_sub_coldspot.txt
bedtools intersect -a indels_40bp_all.bed -b hotspotRand_2_2kb_down.bed -wo > intersect_40bp_2kbdown_sub_coldspot.txt
bedtools intersect -a indels_50bp_all.bed -b hotspotRand_2_2kb_down.bed -wo > intersect_50bp_2kbdown_sub_coldspot.txt

bedtools intersect -a indels_all.bed -b hotspotRand_2_2kb_up.bed -wo > intersect_2kbup_sub_coldspot.txt
bedtools intersect -a indels_10bp_all.bed -b hotspotRand_2_2kb_up.bed -wo > intersect_10bp_2kbup_sub_coldspot.txt
bedtools intersect -a indels_20bp_all.bed -b hotspotRand_2_2kb_up.bed -wo > intersect_20bp_2kbup_sub_coldspot.txt
bedtools intersect -a indels_30bp_all.bed -b hotspotRand_2_2kb_up.bed -wo > intersect_30bp_2kbup_sub_coldspot.txt
bedtools intersect -a indels_40bp_all.bed -b hotspotRand_2_2kb_up.bed -wo > intersect_40bp_2kbup_sub_coldspot.txt
bedtools intersect -a indels_50bp_all.bed -b hotspotRand_2_2kb_up.bed -wo > intersect_50bp_2kbup_sub_coldspot.txt
