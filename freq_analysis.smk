import re
import pandas as pd

singularity: '/panfs/jay/groups/0/fried255/shared/gatk4_workflow/rescuer/AFREQ/CalcFix/afreq.sif'

# MAKE THIS A PART OF THE AF CONTAINER - NO NEED TO REGEN EVERYTIME
chrmsDF = pd.read_table("dog.regions.txt")
regions = list(chrmsDF['chrom'])

rule all:
    input:
       #expand(
       #   #"results/{region}/{region}.split.vcf.gz",
       #    "results/{region}/{region}.csv",
       #    region=regions
       #   #region='chr27.31108325-46662488'
       #),
        "af_analysis.no_intergenic.csv"  

# ADD MORE STRUCTURE TO DIRS LIKE split/{region}/{region}.vcf.gz
rule select_variants_chrom:
    input:
        vcf = "/panfs/jay/groups/0/fried255/fried255/working/pipeline/UU_Cfam_GSD_1.0_ROSY/20230809/joint_call.UU_Cfam_GSD_1.0_ROSY.20230809.vep.vcf.gz" 
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
        "Toy/Steve/Lists/stpd.list",
        "Toy/Steve/Lists/pwdg.list",
        "Toy/Steve/Lists/oodle.list",
        "Toy/Steve/Lists/ctrlbrd.list",
        vcf = "results/{region}/{region}.split.vcf.gz"
       #vcf = "chr27_snp.vcf"
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
                {input[0]} \
                {input[1]} \
                {input[2]} \
                {input[3]} > {output.table}
        '''

rule concat_files:
    input:
        expand("results/{region}/{region}.csv", region=regions)
    output:
        "af_analysis.csv" 
    params:
       #header = "chrm,pos,ref,alt,major,minor,maf.all,aff.count.homref,aff.count.het,aff.count.homvar,aff.count.nocall,aff.maf,risk.count.homref,risk.count.het,risk.count.homvar,risk.count.nocall,risk.maf,ctl.count.homref,ctl.count.het,ctl.count.homvar,ctl.count.nocall,ctl.maf,all.count.homref,all.count.het,all.count.homvar,all.count.nocall,allele,consequence,impact,gene,transcript,biotype,exon,hgvsc,hgvsp,cdna_position,cds_position,protein_position,amino_acid,codons,variant_class,protein"
        header = (
            "chrm,pos,ref,alt,major,minor,maf.all,"
            "stpd.count.homref,stpd.count.het,stpd.count.homvar,stpd.count.nocall,stpd.maf,"
            "pwdg.count.homref,pwdg.count.het,pwdg.count.homvar,pwdg.count.nocall,pwdg.maf,"
            "oodles.count.homref,oodles.count.het,oodles.count.homvar,oodles.count.nocall,oodles.maf,"
            "ctrl.count.homref,ctrl.count.het,ctrl.count.homvar,ctrl.count.nocall,ctrl.maf,"
            "all.count.homref,all.count.het,all.count.homvar,all.count.nocall,Allele,"
            "Consequence,IMPACT,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,cDNA_position,"
            "CDS_position,Protein_position,Amino_acids,Codons,DISTANCE,STRAND,VARIANT_CLASS,"
            "ENSP,SOURCE,SIFT_pred,SIFT_score,PhyloP_score"
        )
    threads: 1
    resources:
        time = 600,
        mem_mb = 20000
    shell:
        '''
            echo {params.header} > {output}
            cat {input} >> {output}
        '''

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

