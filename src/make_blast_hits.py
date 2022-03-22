#!/usr/bin/env python3

import sys
import math

inputfile1 = sys.argv[1]
name1 = inputfile1.split('/')
name1 = name1[len(name1)-1].split('.')[0]
dc1 = {}
with open(inputfile1) as f:
	for i in f:
		info = i.strip().split('\t')
		gene = info[3].split('_')[2]
		dc1.setdefault(info[4],[]).append(gene)
#
inputfile2 = sys.argv[2]
name2 = inputfile2.split('/')
name2 = name2[len(name2)-1].split('.')[0]
dc2 = {}
with open(inputfile2) as f:
	for i in f:
		info = i.strip().split('\t')
		gene = info[3].split('_')[2]
		dc2.setdefault(info[4],[]).append(gene)


for key1 in dc1.keys():
	info1 =  dc1[key1]
	for key2 in dc2.keys():
		info2 = dc2[key2]
		overlap = list(set(info1).intersection(set(info2)))
		length1 = len(info1)
		length2 = len(info2)
		overlap_length = float(len(overlap))
		overlap_ratio1 = overlap_length/length1
		overlap_ratio2 = overlap_length/length2
		
		if len(overlap) > 0:
			print('\t'.join([name1+'.'+key1,name2+'.'+key2,'100.000','100','0','0','1','100','1','100','1.01e-10','100']))

