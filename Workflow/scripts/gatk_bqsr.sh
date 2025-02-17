#!/bin/bash

# Script to perform BQSR using GATK
# Usage: bash run_bqsr.sh <input_bam> <reference_fasta> <known_sites_vcf> <output_bam>

# Check input arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: bash $0 <input_bam> <reference_fasta> <known_sites_vcf> <output_bam>"
    exit 1
fi

# Assign input arguments to variables
INPUT_BAM=$1
REFERENCE_FASTA=$2
KNOWN_SITES_VCF=$3
OUTPUT_BAM=$4


# Path to GATK (modify this path if GATK is not in your PATH)
GATK_CMD="workflow/scripts/gatk-4.2.3.0/gatk --java-options -Xmx200g"
echo GATK_CMD
# Step 1: Generate recalibration table
echo "Step 1: Generating recalibration table..."
$GATK_CMD BaseRecalibrator \
    -I $INPUT_BAM \
    -R $REFERENCE_FASTA \
    --known-sites $KNOWN_SITES_VCF \
    -O recal_data.table

if [ $? -ne 0 ]; then
    echo "Error in BaseRecalibrator step. Exiting..."
    exit 1
fi

# Step 2: Apply recalibration
echo "Step 2: Applying recalibration..."
$GATK_CMD ApplyBQSR \
    -I $INPUT_BAM \
    -R $REFERENCE_FASTA \
    --bqsr-recal-file recal_data.table \
    -O $OUTPUT_BAM

if [ $? -ne 0 ]; then
    echo "Error in ApplyBQSR step. Exiting..."
    exit 1
fi

echo "BQSR completed successfully. Recalibrated BAM file: $OUTPUT_BAM"

