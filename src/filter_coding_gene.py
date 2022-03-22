#!/usr/bin/env python3
import sys
import re

infile = sys.argv[1]
gff3 = sys.argv[2]
high_coverage_gene = sys.argv[3]
pt_fai = sys.argv[4]

coverage = 0.90

#
high_cov_gene = {}
with open(high_coverage_gene) as f:
    for line in f:
        high_cov_gene[line.strip()] = ''

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
# find end codon
def end_codon_index(seq1):
    lt = seq1.split('.',1)
    if len(lt) == 1:
        return len(seq1)-1
    else:
        return len(lt[0])-1

# gene length
pt_dc = {}
with open(pt_fai) as fa:
    for i in fa:
        i = i.strip().split('\t')
        pt_dc[i[0]] = int(i[1])

# input pep sequence
dc = {}
with open(infile) as fa:
    for line in fa:
        if line.startswith('>'):
            ID = line.strip().split('>')[1].split(' ')[0]
            dc[ID] = ''
        elif not line.startswith('#'):
            dc[ID] += line.strip('\n')

# filter no translation
trans_geneid = {}
for i in dc:
    if i in high_cov_gene:
        seq = dc[i]
        geneid = i.split('_')[2]
        start_index = start_codon_index(seq,geneid)
        if start_index == 0:
            end_index = end_codon_index(seq)
            genelegth = pt_dc[geneid]
            new_length = end_index-start_index+1
            gene_coverage = float(new_length)/genelegth
            if gene_coverage >= coverage:
                gene_id = i.split(':')[1]
                trans_geneid[gene_id] = ''

#
#
with open(gff3) as f:
    for i in f:
        i = i.strip().split('\t')
        if i[2] == 'gene' or i[2] == 'mRNA':
            geneid = i[8].split(';')[0]
            geneid_info = re.split('[=:]',geneid)
            geneid = geneid_info[len(geneid_info)-1]
            info = '\t'.join(i)
            if geneid in trans_geneid:
                print (info)
            else:
                info = info.replace(geneid,geneid+'_pseudogene')
                print(info)
        else:
            geneid = i[8].split(';')[0]
            geneid_info = re.split('[=:.]',geneid)
            geneid = geneid_info[len(geneid_info)-2]
            info = '\t'.join(i)
            if geneid in trans_geneid:
                print ('\t'.join(i))
            else:
                info = info.replace(geneid,geneid+'_pseudogene')
                print(info)
