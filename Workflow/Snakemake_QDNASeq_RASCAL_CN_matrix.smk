__authors__ = "Srinivas Veerla, Guyuan TANG"
__copyright__ = "Copyright 2024, Srinivas Veerla"
__email__ = "srinivas.veerla@med.lu.se"
__license__ = "GPL-3"

import pandas as pd

# specify the configuration file
configfile: "config/config.yaml"

# specify the samples and their groups
sample_df = (pd.read_csv(config['samples'],     
    dtype={'sample_name':str, 'path':str})
    .set_index('sample_name', drop=False))

# specify working directory
wdir = config['workdir']

# specify the results location (output directory)
results = config['outputdir']

# specify the bin size used for annotation
binsize = str(config['QDNAseq']['binsize'])

rule all:
	input:  expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_df.sample_name)	

rule bamFiles:
    input: expand(wdir + '{sample}.bam', sample=sample_df.sample_name)
    
rule relative_CN:
    input:
        rds = expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds', sample=sample_df.sample_name),
        tsv = expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.seg.tsv', sample=sample_df.sample_name)

rule abs_seg_CN:
    input:
        expand(results + '{sample}/05_absolute_CN/{sample}_' + binsize + 'kb_seg.tsv', sample=sample_df.sample_name)

rule CN_solution: # the rule is the same with rule all at this moment
    input:
        expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_df.sample_name)

        
########## Generating relative CN profile ###########
rule QDNAseq:
    ### QDNAseq will be applied to generate relative copy number profile stored in RDS and tsv files for later analyses    
    input: 
    	link_up = rules.bamFiles.input,
    	bamfile = wdir + '{sample}.bam'
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
    script: 'scripts/runQDNAseq.R'
    
########## Ploidy and cellularity solution ################
"""
The final output for this step would be the solutions of ploidy and tumour purity for each sample. It is also the last step in the first snakemake file (Snakefile_solution.smk)
"""
rule rascal_env:
    output:
        "log/rascal_settle_info.txt"
    script: 'scripts/rascal_env.R'

########## Calculate the optimal solutions (ploidy and cellularity) of the samples
rule rascal_solution:
    ### Rascal will be applied to calculate the optimal solutions for further deriving the absolute copy numbers
    input:
        link_up = rules.relative_CN.input,
        env_set = 'log/rascal_settle_info.txt',
        rds = results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds'
    output:
        solution = results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv'
    params:
        output_prefix = results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb',
        min_cellularity = config['Rascal']['min_cellularity'],
        script = 'scripts/fit_CN_solution.R'
    threads: 5
    shell: '''
    Rscript {params.script} -i {input.rds} -o {params.output_prefix} --min-cellularity {params.min_cellularity}
    '''

########## CN Signatures and sample-by-component matrix ####################
"""
The final output for this step would be the matrix files (including a matrix txt file, a matrix object RDS, a simple heatmap for both sample-by-component matrix and sample-by-signature matrix) for each group. The matrices containing sample-by-signature information.
"""
# Check all the files are well prepared
rule cn_sig_check:
    input:
        link_up = rules.abs_seg_CN.input,
        scripts = rules.cn_sig_git.output
    output:
        'log/cn_sig_settle.txt'
    shell: """
    echo 'Finished preparation for signature validation.' > {output}
    """

# Validate the CN signatures in our samples
rule CN_signature:
    input:
        # check whether the essential git repository have been downloaded
        main = 'scripts/cnsignatures/main_functions.R',
        helper = 'scripts/cnsignatures/helper_functions.R',
        # check all the absolute copy number profiles have successfully been generated
        check_data = rules.cn_sig_check.output,
        # the input segment files
        seg_files = rules.abs_seg_CN.input
    output:
        CN_sig = results + 'signatures/CN_sig/CN_sig.SSmatrix.rds'
    params:
        sample_info = config['samples'],
        def_SC = 'scripts/cnsignatures/data/feat_sig_mat.rds',
        indir = results,
        binsize = binsize,
        outdir = results + 'signatures/CN_sig/'
    threads: 20
    script: 'scripts/CN_sig.R'

