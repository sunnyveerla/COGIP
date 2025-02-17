__authors__ = "Srinivas Veerla, Guyuan TANG"
__copyright__ = "Copyright 2024, Srinivas Veerla"
__email__ = "srinivas.veerla@med.lu.se"
__license__ = "GPL-3"

import pandas as pd

# specify the configuration file
configfile: "config/config.yaml"

# specify the samples and their groups
sample_info = (pd.read_table(config['samples'],     
    dtype={'sample_name':str, 'path':str})
    .set_index('sample_name', drop=False))

# specify working directory
wdir = config['workdir']

# specify the results location (output directory)
results = config['outputdir']

# specify the bin size used for annotation
binsize = str(config['QDNAseq']['binsize'])

rule all:
	input:  expand(results + '{sample}/05_absolute_CN/solutions/{sample}_' + binsize + 'kb.solution.csv', sample=sample_info.sample_name)	

rule relative_CN:
    input:
        rds = expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.rds', sample=sample_info.sample_name),
        tsv = expand(results + '{sample}/04_relative_CN/' + binsize + 'kb/{sample}_' + binsize + 'kb.seg.tsv', sample=sample_info.sample_name)

rule bamFiles:
    input: expand(wdir + '{sample}.bam', sample=sample_info.sample_name)
        
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



