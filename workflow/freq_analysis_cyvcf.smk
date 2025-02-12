
import pandas as pd

regions = list(pd.read_table("/share/stern/mwvandew/joanna_cats/ref/cat/felcat126/GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.chunks.txt", header = None)[0])

rule all:
    input:
        "output/final/af_analysis.no_intergenic.csv"
 
rule process_variants:
    input:
        list_dir = "lists",
	vcf = "../colony_cats.vep.filter.noJim.vcf.gz" 
    output:
        table = "output/split/{region}.csv"
    params:
        region = "{region}"
    threads: 1
    conda: "env/hts.yaml"
    resources:
        time = 1440,
        mem_mb = 20000
    shell:
        '''
            ./calc_af_cyvcf.py \
                -v {input.vcf} \
                -d {input.list_dir} \
                -r {params.region} \
		> {output.table}
        '''

rule concat_files:
    input:
        expand("output/split/{region}.csv", region=regions)
    output:
        "output/final/af_analysis.csv"  
    threads: 1
    resources:
        time = 60,
        mem_mb = 20000
    shell:
        '''
            cat {input} |awk 'NR==1 || !/chrm/' > {output} 
        '''

localrules: remove_intergenic
rule remove_intergenic:
    input:
        "output/final/af_analysis.csv" 
    output:
        "output/final/af_analysis.no_intergenic.csv" 
    shell:
        '''
            awk 'NR == 1 || !/intergenic/' {input} > {output}
        '''

