import os
import sys
from collections import Counter
import vcf
from pathlib import Path

def reduceCSQ(csqs):
    mainCols = [0, 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 18, 19, 21, 24, 25, 27, 28, 30]
    updated = []
    for record in csqs:
        record = record.split('|')
        if '&' in record[-1]:
            record[-1] = "&".join(sorted(set(record[-1].split('&')), key=record[-1].split('&').index))
        result = ','.join([record[i] for i in mainCols])
        updated.append(result)

    return updated

def calcAF(groupList, minorAlt=True):
    #calculate MAF
    alleles = [allele for genotype in groupList for allele in genotype.split('/') if allele != '.']
    minor_count = len(list(filter(lambda x: x != '0', alleles))) if minorAlt else len(list(filter(lambda x: x == '0', alleles)))
    try:
        maf = round(minor_count/len(alleles), 6)
    except ZeroDivisionError:
        maf = 'NA'

    genotypes = Counter(groupList)
    homRef = genotypes.get('0/0', 0)
    miss = genotypes.get('./.', genotypes.get('.',0))
    het = sum(count for genos, count in genotypes.items() 
              if genos not in ['0/0', './.', '.'] and genos.split('/')[0] != genos.split('/')[1])
    homVar = sum(count for genos, count in genotypes.items() 
                 if genos not in ['0/0', './.', '.'] and genos.split('/')[0] == genos.split('/')[1])
    
    return [homRef, het, homVar, miss, maf]

# Parse dog lists dynamically
dogs = {}
group_order = []
for list_file in sys.argv[2:]:
    group = Path(list_file).stem
    group_order.append(group)
    with open(list_file) as f:
        for sample in f:
            dogs[sample.strip()] = group

vcfR = vcf.Reader(filename=sys.argv[1])

# Generate header based on groups
base_cols = ['chrm', 'pos', 'ref', 'alt', 'major', 'minor', 'maf.all']
group_cols = []
for group in group_order:
    group_cols.extend([f'{group}.count.homref', f'{group}.count.het', 
                      f'{group}.count.homvar', f'{group}.count.nocall', f'{group}.maf'])
all_cols = ','.join([*base_cols, *group_cols, 'all.count.homref', 'all.count.het', 
                     'all.count.homvar', 'all.count.nocall'])

print(f"{all_cols},Allele,Consequence,IMPACT,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON," \
      f"cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,DISTANCE,STRAND," \
      f"VARIANT_CLASS,ENSP,SOURCE,SIFT_pred,SIFT_score,PhyloP_score")

for line in vcfR:
    chrm, pos, ref, alt = line.CHROM, line.POS, line.REF, '/'.join(map(str, line.ALT))
    AF = sum(line.INFO.get('AF', [0]))
    minor_alt = AF < 0.5
    major, minor, maf = (ref, alt, AF) if minor_alt else (alt, ref, 1-AF)
    first = ','.join(map(str, [chrm, pos, ref, alt, major, minor, round(maf, 6)]))
    
    # process groups dynamically
    main = []
    group_data = {group: [] for group in group_order}
    
    for sample in line.samples:
        gt = sample['GT'].replace('|', '/')
        main.append(gt)
        group = dogs.get(sample.sample)
        if group:
            group_data[group].append(gt)
    
    # calculate stats per group
    group_stats = []
    for group in group_order:
        stats = calcAF(group_data[group], minor_alt)
        group_stats.append(','.join(map(str, stats)))
    
    # calculate all stats
    all_stats = ','.join(map(str, calcAF(main, minor_alt)[:4]))
    
    # process consequences
    csq = reduceCSQ(line.INFO.get('CSQ', []))
    for consequence in csq:
        print(f"{first},{','.join(group_stats)},{all_stats},{consequence}")

