## Integrative Genomic Analysis of High-Grade Serous Ovarian Cancer Using Liquid-Based Cytology and Tumor Samples – modelling genomic instability

This project focuses on high-grade serous ovarian cancer (HGSC), an aggressive and lethal gynecologic malignancy driven by genomic instability. Previous studies from our group have assessed genomic aberrations, particularly copy number alterations (CNAs), to classify tumors and liquid-based cervical cytology samples from women with and without HGSC. However, existing CNA-based approaches are inherently subjective and limited in capturing focal genomic events with molecular specificity, providing only a broad genomic overview. Consequently, these methods constrain our deeper understanding of the mechanisms driving tumorigenesis and metastasis progression.

To overcome these limitations, we are developing an advanced framework that integrates shallow whole-genome sequencing (sWGS) data with SRIQ, a machine-learning method for identifying significant and shared focal and long CNA regions across diverse sample types, including cervical cytology samples, blood, plasma, and tumor tissues. To enhance specificity, CNA regions overlapping with normal control samples are excluded. The identified CNAs will be validated using an external cohort dataset with similar sample types. Additionally, given the genomic similarities between HGSC and triple-negative breast cancer (TNBC), particularly basal-like subtypes, we will further validate these CNA regions using TNBC data from the SCAN-B. Finally, these CNAs will be mapped to genes and regulatory regions to assess their functional relevance. The role/s of the involved genes will further be evaluated in relation to tumor and patient-related variables, including degree of differentiation, metastatic spread and immune cell markers (as examples) using in-house tissue microarrays (n=140 HGSC, n=231 TNBC, n=yyy SEC). 

Taken together by integrating advanced genomic analysis, machine learning, and bioinformatics, this study aims to uncover novel genomic aberrations and their molecular functions across diverse tumor and sample types, especially liquid-based. These insights will enhance our understanding of genomic instability in tumorigenesis and metastasis, and key oncogenic drivers, ultimately contributing to the development of more effective targeted therapies.

#### COGIP (Characterization of Oncogenic Genomic Instability analysis Pipeline):
The version of tools and packages to be used will be specified in each step (see Chapter 3). The scripts within the pipeline are based on Snakemake (6.15.1), Python (v3.11.6), R (v4.3.2), and Java (JDK 11).

(1) **Preprocessing**. This step includes quality assessment and quality trimming on the raw reads. (Fastp will be used for QC and trimming, together with fastqc and multiQC to generate the QC reports.)

(2) **Alignment**. The human reference genome will be indexed. The reads will be mapped to the reference genome. (BWA will be used for both indexing and alignment.)

(3) **Clean-up**. After alignment, the SAM files will be sorted and the PCR duplicates will be marked and removed. Also, the .sorted.deduplicated.sam will be converted to BAM files. The BAM files will be indexed for later analysis. (Picard will be used for sorting SAM, marking duplicates, removing duplicates and converting SAM to BAM. samtools will be used for generating the clean_up stats and for indexing the BAM files.)

(4) **Relative copy number profile**. The BAM files will be analyzed through fixed-size binning, filtering, correction, and normalization to generate the read counts per bin. This data will then be used for the segmentation of bins and for generating the relative copy number profile. (QDNAseq will be used for this step.)

(5) **Ploidy and cellularity solutions**. The output file from QDNAseq contains a relative copy number, and we need to estimate ploidy and cellularity in our samples to generate our final absolute copy number profile for comparison. (Rascal will be used for this step to find the solutions that best fit our study samples.)

(6) **Absolute copy number profile**. We will further use other information (such as TP53 allele frequency) to infer the tumour fraction to select the best ploidy and cellularity solution. We apply this best solution to our relative copy number profile and generate the final absolute copy number profile for each sample. (Rascal will be used for this step.)

#### (1) Preprocessing, (2) Alignment, (3) Clean-up, (4) Relative copy number profile, (5) Ploidy and cellularity solutions, (6) Absolute copy number profile,
```bash
snakemake --use-conda  --configfile config/config.yaml --cores 30  --snakefile workflow/Snakemake_HCsig_Pipeline.smk
snakemake --use-conda  --configfile config/config.yaml --cores 30  --snakefile workflow/Snakemake_QDNASeq_RASCAL_CN_matrix.smk
```
