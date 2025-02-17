__authors__ = "Srinivas Veerla, Guyuan TANG"
__copyright__ = "Copyright 2024, Srinivas Veerla"
__email__ = "srinivas.veerla@med.lu.se"
__license__ = "GPL-3"

import pandas as pd

#Configuration file
configfile: "config/config.yaml"

#Read sample information file
sample_info = (pd.read_table(config['samples'],     
    dtype={'sample_name':str, 'patient':str, 'type':str,'fastq_1':str, 'fastq_2':str})
    .set_index('sample_name', drop=False))

#Working directory location
work_dir = config['workdir']

#Result data files storage location
results = config['outputdir']

# specify the bin size used for annotation
binsize = str(config['QDNAseq']['binsize'])
print(binsize)

rule all:	
    #input:  expand(results + '{sample}/02_alignment/{sample}.sorted.bam', sample=sample_info.sample_name)	
    #input: expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.seg.tsv' , sample=sample_info.sample_name)
    #input: expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_df.sample_name)
    input:  expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_info.sample_name)
    

########### 1 sWGS raw data preprocessing ##########
rule fastp:
    input:
        R1 = lambda wildcards: sample_info.loc[wildcards.sample, 'fastq_1'],
        R2 = lambda wildcards: sample_info.loc[wildcards.sample, 'fastq_2']
    output:
        R1 = results + '{sample}/01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        html = results + '{sample}/01_preprocess/html/{sample}_fastp.html',
        R2 = results + '{sample}/01_preprocess/reads/{sample}_R2_preprocess.fastq.gz' 
    log: 'log/fastp/{sample}_fastp.log'
    threads: 20
    params: json = results + '{sample}/01_preprocess/html/{sample}_fastp.json'
    #conda: "envs/preprocess_env.yaml"
    shell: """
        fastp --detect_adapter_for_pe \
        --correction --cut_right --thread {threads} \
        --html {output.html} --json {params.json} \
        --in1 {input.R1} --in2 {input.R2} \
        --out1 {output.R1} --out2 {output.R2} \
        2>{log}
    
    rm {params.json}
    """
 # quality assessment of preprocessed reads with fastqc
rule fastqc:
    ### we will use fastqc to generate the quality control stats from the outputs of fastp
    input:
        R1_seq = results + '{sample}/01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        R2_seq = results + '{sample}/01_preprocess/reads/{sample}_R2_preprocess.fastq.gz'
    output:
        R1_html = results + '{sample}/01_preprocess/html/{sample}_R1_preprocess_fastqc.html',
        R1_qc = results + '{sample}/01_preprocess/reports/{sample}_R1_preprocess_fastqc.zip',
        R2_html = results + '{sample}/01_preprocess/html/{sample}_R2_preprocess_fastqc.html',
        R2_qc = results + '{sample}/01_preprocess/reports/{sample}_R2_preprocess_fastqc.zip'
    log: 'log/fastqc/{sample}.fastqc.log'
    params: 
        outdir = results + '{sample}/01_preprocess/reports/',
        out_html = results + '{sample}/01_preprocess/html/'
    threads: 20
    #conda: 'envs/preprocess_env.yaml'
    shell: """
    fastqc -o {params.outdir} {input.R1_seq} {input.R2_seq} 2>{log}
    mv {params.outdir}*_fastqc.html {params.out_html}
    """

########## 2 Alignment ####################

# 2 mapping the reads to the indexed reference genome
rule map_reads:
    ### use bwa again for alignment
    input: 
        #idx = rules.bwa_index.output,
        link_up = rules.fastqc.input,
        R1 = results + '{sample}/01_preprocess/reads/{sample}_R1_preprocess.fastq.gz',
        R2 = results + '{sample}/01_preprocess/reads/{sample}_R2_preprocess.fastq.gz'
    output:
        results + '{sample}/02_alignment/{sample}.sorted.bam'
    log: 'log/bwa_mapping/{sample}/{sample}.log'
    params:
        index_ref = 'resources/genome/hg38.fa'
    #conda: 'envs/alignment.yaml'
    threads: config['bwa_mapping']['threads']
    shell: '''
    workflow/scripts/bwa-mapper.sh {input.R1} {input.R2} {params.index_ref} {threads} {output} 2>{log}
    samtools index -@ {threads} {output}
    '''

# 3 MarkDuplicates and de-duplication
rule de_duplicate:
    ### using Picard to remove PCR duplicates
    input: 
        link_up = rules.map_reads.output,
        bamFile= results + '{sample}/02_alignment/{sample}.sorted.bam'
    output:
        sorted_bam = results + '{sample}/03_clean_up/{sample}.sorted.bam',
        tmp_dir = directory(results + '{sample}/tmp/'),
        dedup_bam = results + '{sample}/03_clean_up/{sample}.sorted_dedup.bam'        
    log: 'log/de_duplicate/{sample}.log'
    threads:config['bwa_mapping']['threads']
    params:
        metrix_file = results + '{sample}/03_clean_up/{sample}.metrics.txt'
    #conda: 'envs/clean_up.yaml'
    shell: """
       mkdir -p {output.tmp_dir}

    samtools index -@ {threads} {input.bamFile}

    java -Xmx200G -Djava.io.tmpdir={output.tmp_dir} -jar workflow/scripts/picard.jar \
        SortSam \
        I={input.bamFile} \
        O={output.sorted_bam} \
        SORT_ORDER=coordinate \
        TMP_DIR={output.tmp_dir} \
        2>{log}


    #rm -rf {output.tmp_dir}

   java -Xmx200G -Djava.io.tmpdir={output.tmp_dir} -jar workflow/scripts/picard.jar MarkDuplicates \
        INPUT={output.sorted_bam} \
        OUTPUT={output.dedup_bam} \
        METRICS_FILE={params.metrix_file} \
        REMOVE_DUPLICATES=true \
        ASSUME_SORT_ORDER=coordinate \
        CLEAR_DT=false \
        2>{log}
    
    #samtools view -bo {output.dedup_bam} {input}
    #qualimap bamqc -bam {output.dedup_bam} --java-mem-size=100G   
    #rm -f {output.sorted_bam}
    """

########## 3 GATK BaseRecalibration ##########

rule GATK_bqsr:
    input:
        link_up = rules.de_duplicate.output,
        sorted_bam = results + '{sample}/03_clean_up/{sample}.sorted.bam',
    	bamfile =  results + '{sample}/03_clean_up/{sample}.sorted_dedup.bam'  
    output:
         results + '{sample}/03_clean_up/{sample}.sorted_dedup_GATK_bqsr.bam'
    log: 'log/GATK_bqsr_mapping/{sample}/{sample}.log'
    params:
        index_ref = 'resources/genome/hg38.fa',
        known_sites = 'resources/hg38/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz'
    #conda: 'envs/alignment.yaml'
    threads: config['bwa_mapping']['threads']
    shell: '''
     workflow/scripts/gatk_bqsr.sh {input.bamfile} {params.index_ref} {params.known_sites} {output}  2>{log}
     samtools index -@ {threads} {output}
     rm -f {input.sorted_bam}
     rm -f {input.bamfile}
    '''

########## 4 Generating relative CN profile ###########
rule QDNAseq:
    ### QDNAseq will be applied to generate relative copy number profile stored in RDS and tsv files for later analyses    
    input: 
    	link_up = rules.GATK_bqsr.output,
    	bamfile = results + '{sample}/03_clean_up/{sample}.sorted_dedup_GATK_bqsr.bam'
    output:
        rds = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds',
        igv = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.igv',
        seg_tsv = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.seg.tsv'    
    params:
        sample = '{sample}',
        binsize = config['QDNAseq']['binsize'],
        outdir = results + '{sample}/04_relative_CN/' + binsize + 'kb/',
        maxSize = config['QDNAseq']['maxSize']
    threads: config['QDNAseq']['threads']   
    script: 'scripts/runQDNAseq_custom.R' 

########## 5 Ploidy and cellularity solution ################

rule CN_solution: # the rule is the same with rule all at this moment
    input:
        expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_info.sample_name)


########## 6 Calculate the optimal solutions (ploidy and cellularity) of the samples ##########
rule rascal_solution:
    ### Rascal will be applied to calculate the optimal solutions for further deriving the absolute copy numbers
    input:
        link_up = rules.QDNAseq.output,
        rds = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds'
    output:
        solution = results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv'
    params:
        output_prefix = results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb',
        min_cellularity = config['Rascal']['min_cellularity'],
        script = 'workflow/scripts/fit_CN_solution.R'
    shell: '''    
    echo "Its here"
    Rscript {params.script} -i {input.rds} -o {params.output_prefix} --min-cellularity {params.min_cellularity} 
    '''


