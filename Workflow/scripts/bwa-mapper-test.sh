#!/bin/bash

header=$(zcat $1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
echo "Read Group @RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"
path_bwaindex_genome=$3
threads=$4
output=$5
bwa mem \
-M \
-t $threads \
-v 3 \
"$path_bwaindex_genome" \
$1 $2 | samblaster -M | samtools fixmate - - | samtools sort -@ $threads -O bam -o $output
