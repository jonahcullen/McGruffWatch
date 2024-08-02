import re
import pandas as pd
# needed to list the lists we have from the list directory in a system-agnostic way
import os

# snakemake documentation implies that setting this casuses it to be the default config file, but it will be overriden if the user uses --configfile "" as an arg
configfile:"./af_config.yaml"

singularity: config["sif_path"]

# MAKE THIS A PART OF THE AF CONTAINER - NO NEED TO REGEN EVERYTIME
# read in the dog.regions.txt file as a table using pandas, and a regions variable that's a list of those chrom values
chrmsDF = pd.read_table("dog.regions.txt")
regions = list(chrmsDF['chrom'])

# get the lists that exist in the provided list directory and check that they have a .list extension, make a list of their paths
lists = os.listdir(config["lists_directory"])
list_paths = list()

# generate a header based on the group names we have already
new_header = "chrm,pos,ref,alt,major,minor,maf.all,"
for filename in lists:
	if filename[-5:] == ".list":
		group_name = str(filename.split('.')[0])
		new_header = new_header + f"{group_name}.count.homref,{group_name}.count.het,{group_name}.count.homvar,{group_name}.count.nocall,{group_name}.maf,"
		list_paths.append(config["lists_directory"] + "/" + filename)
# finish off the new header
new_header = new_header + "all.count.homref,all.count.het,all.count.homvar,all.count.nocall,Allele,Consequence,IMPACT,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,DISTANCE,STRAND,VARIANT_CLASS,ENSP,SOURCE,SIFT_pred,SIFT_score,PhyloP_score"


# for some reason, this needs to not be in the input section of the process_variants thing
af_args = " ".join(list_paths)

rule all:
    input:
        "af_analysis.no_intergenic.csv"  

# ADD MORE STRUCTURE TO DIRS LIKE split/{region}/{region}.vcf.gz
# takes in a vcf file and outputs many processed vcf files using the options provided and bcftools
rule select_variants_chrom:
    input:
        vcf = config["vcf_path"]
    output:
        vcf = "results/{region}/{region}.split.vcf.gz"
    params:
        reg = lambda wildcards, output: re.sub(r'\.(?=[^.]*$)', ':', wildcards.region)
    threads: 12
    resources:
        time = 600,
        mem_mb = 40000
    shell:
        '''
            bcftools view \
                -o {output.vcf} \
                -r {params.reg} \
                --threads {threads} \
                -Oz \
                {input.vcf}
        '''

rule process_variants:
    input:
        vcf = "results/{region}/{region}.split.vcf.gz"
    output:
        table = "results/{region}/{region}.csv"
    threads: 1
    resources:
        time = 1440,
        mem_mb = 20000
    shell:
        '''
            ./af_analysis.mv.py \
                {input.vcf} \
                {af_args} > {output.table}
        '''

rule concat_files:
    input:
        expand("results/{region}/{region}.csv", region=regions)
    output:
        "af_analysis.csv" 
    params:
        header = str(new_header)
    threads: 1
    resources:
        time = 600,
        mem_mb = 20000
    shell:
        '''
            echo {params.header} > {output}
            cat {input} >> {output}
        '''

# takes the generated above af_analysis.csv and removes specific areas from it (intergenic ones), and generates a new output without intergenic regions
# this localrules means this isn't submitted as its own job, and is just run on the host node since it's so easy/short (I THINK???)
localrules: remove_intergenic
rule remove_intergenic:
    input:
        "af_analysis.csv" 
    output:
        "af_analysis.no_intergenic.csv" 
    shell:
        '''
            awk 'NR == 1 || !/intergenic/' {input} > {output}
        '''
