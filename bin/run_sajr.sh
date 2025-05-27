#!/bin/bash -e

path=`dirname $0`

echo "All parameters: $@"

bam=$1
out=$2
ann=$3
strand=$4
genome_pos=$5
max_reads=$6
use_bam_dedupl=$7

cmd="java -Xmx20G -jar ${path}/sajr.ss.jar \
    count_reads ${path}/sajr.config \
        -batch_in=${bam} \
        -batch_out=${out} \
        -ann_in=${ann} \
        -stranded=${strand} \
        -genome_pos=${genome_pos} \
        -count_only_border_reads=true \
        -max_reads=${max_reads} \
        -use_bam_dedupl=${use_bam_dedupl}"

echo $cmd
$cmd