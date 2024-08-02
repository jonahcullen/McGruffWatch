#!/usr/bin/env python

import os
import sys
from collections import Counter

import vcf

# ask if this mainCols var needs to be generalized at some point? - H
def reduceCSQ(csqs):
    mainCols = [0, 1, 2, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 18, 19, 21, 24, 25, 27, 28, 30]
    updated = []
    for record in csqs:
        record = record.split('|')
        if '&' in record[-1]:
            record[-1] = "&".join(sorted(set(record[-1].split('&')), key=record[-1].split('&').index))
        result = ','.join([record[i] for i in mainCols])
        updated.append(result)
    
    return(updated)

def majority_element(nums):
    return '1' if sum(int(num) for num in nums) > len(nums) / 2 else '0'

def calcAF(groupList, minorAlt=True):
    #calculate MAF
    alleles = [allele for genotype in groupList for allele in genotype.split('/') if allele != '.']
    # where allele freq < 0.5
    minor_count = len(list(filter(lambda x: x != '0', alleles)))
    if not minorAlt:
        # where allele freq > 0.5
        minor_count = len(list(filter(lambda x: x == '0', alleles)))
    try:
        maf = round(minor_count/len(alleles), 6)
       #maf = round(alleles.count(alt_count)/len(alleles), 6)
    except ZeroDivisionError:
        maf = 'NA'

    #Count Genotypes
    genotypes = Counter(groupList)
    homRef = genotypes.get('0/0', 0)
    miss = genotypes.get('./.', genotypes.get('.',0))
    het = 0
    homVar = 0
    for genos, count in genotypes.items():
        if genos != '0/0' and genos != './.' and genos != '.':
            geno = genos.split('/')
            if geno[0] == geno[1]:
                homVar += count
            else:
                het += count
    
    return([homRef, het, homVar, miss, maf]) 

#parse dog list  
dogs = {}
# make a seperate list to hold group names (I think this is inefficient but I want to KISS here)
group_names = list()
for i in sys.argv[2:]:
    # this line gets the name of the file (or actually, a string of anything before the first . in the final name of the path passed)
    group = os.path.basename(i).split('.')[0]
    group_names.append(group)
    with open(i) as f:
        for sample in f:
            sample = sample.strip()
            dogs[sample] = group

#read vcf
vcfR = vcf.Reader(filename = sys.argv[1])

for line in vcfR:
    chrm, pos, ref, alt = line.CHROM, line.POS, line.REF, '/'.join(map(str, line.ALT))
    AF = line.INFO.get('AF')
    AF = sum(AF)
    minor_alt = True
    if AF < 0.5:
        major, minor, maf = ref, alt, AF
    else:
        major, minor, maf = alt, ref, 1-AF
        minor_alt = False
    maf = round(maf, 6)
    first = ','.join(map(str, [chrm, pos, ref, alt, major, minor, maf]))
    # make a dictionary to store the groups with the name of the group as the key value
    groups_dict = {}
    for group_name in group_names:
        groups_dict[group_name] = []
    # still need main list!
    main = []
    for i in line.samples:
        sample, gt = i.sample, i['GT'].replace('|', '/')
        main.append(gt)
        group = os.path.basename(dogs.get(sample, "")) if dogs.get(sample) is not None else None
        groups_dict[group].append(gt)

    counts = list()
    # for each item in the dictionary of groups, run the usual calcAF and build a new list of them
    for key, value in groups_dict.items():
        counts.append(','.join(map(str, calcAF(groups_dict[key], minor_alt))))
    countsAll = ','.join(map(str,calcAF(main, minor_alt)[:4]))
    csq = reduceCSQ(line.INFO.get('CSQ', []))
    for consequence in csq:
        print(f"{first}," + str(",".join(counts)) + f",{countsAll},{consequence}")
