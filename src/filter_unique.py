#!/usr/bin/env python3
import sys
import re

infile = sys.argv[1]
gff3 = sys.argv[2]
high_coverage_gene = sys.argv[3]

coverage = 0.90

# find start codon
def start_codon_index(seq1,gene):
    if gene == 'rps19':
        start_codon = 'V'
    elif gene == 'rpl2':
        start_codon = 'T'
    else:
        start_codon = 'M'
    if seq1[0] == start_codon:
        return 0
    else:
        return -1

high_cov_gene = {}
with open(high_coverage_gene) as f:
    for line in f:
        high_cov_gene[line.strip()] = ''

# input pep sequence
dc = {}
with open(infile) as fa:
    for line in fa:
        if line.startswith('>'):
            ID = line.strip().split(':')[1].split(' ')[0]
            dc[ID] = ''
        elif not line.startswith('#'):
            dc[ID] += line.strip('\n')

# filter no translation
trans_geneid = {}
for i in dc:
    if 'transcript:'+i in high_cov_gene:
        seq = dc[i]
        geneid = i.split('_')[0]
        if '.' not in seq:
            start_index = start_codon_index(seq,geneid)
            if start_index == 0:
                if geneid not in trans_geneid:
                    trans_geneid[geneid] = [i,len(seq)]
                else:
                    if len(seq) > trans_geneid[geneid][1]:
                        trans_geneid[geneid] = [i,len(seq)]

#
with open(gff3) as f:
    for i in f:
        i = i.strip()
        try:
            geneid = i.split('\t')[8].split(';')[0]
            geneid_info = re.split('[=:]',geneid)
            geneid = geneid_info[len(geneid_info)-1].split('.')[0]
            ref_geneid = geneid.split('_')[0]
            if geneid == trans_geneid[ref_geneid][0]:
                i = i.replace(geneid,ref_geneid)
                print (i)
        except KeyError:
            continue
