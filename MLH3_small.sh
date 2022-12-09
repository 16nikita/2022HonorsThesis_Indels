#!/bin/sh
echo running bedtools intersect
bedtools intersect -a indels_all.bed -b mlh3_peaks.bed -wo > intersect_a_mlh3.txt
bedtools intersect -a indels_10bp_all.bed -b mlh3_peaks.bed -wo > intersect_10bp_a_mlh3.txt
bedtools intersect -a indels_20bp_all.bed -b mlh3_peaks.bed -wo > intersect_20bp_a_mlh3.txt
bedtools intersect -a indels_30bp_all.bed -b mlh3_peaks.bed -wo > intersect_30bp_a_mlh3.txt
bedtools intersect -a indels_40bp_all.bed -b mlh3_peaks.bed -wo > intersect_40bp_a_mlh3.txt
bedtools intersect -a indels_50bp_all.bed -b mlh3_peaks.bed -wo > intersect_50bp_a_mlh3.txt
