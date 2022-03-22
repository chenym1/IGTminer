#!/usr/bin/env python3

import sys

total_list = []

with open(sys.argv[1]) as f:

	tmp = f.readline()
	while tmp:
		if tmp.startswith('# seqname source feature start end score strand frame attributes'):
			tmp = f.readline()
			tmp = f.readline()
			group = []
			while not tmp.startswith('#'):
				group.append(tmp.strip())
				tmp = f.readline()
			total_list.append(group)
		else:
			tmp = f.readline()

#
name = sys.argv[2]

gene_num = 1
for i in range(len(total_list)):
	info = total_list[i]
	geneid = ''
	cds_num = 1
	exon_num = 1
	for j in range(len(info)):
		tmp = info[j].split('\t')
		label = 'HGTminer'
		if tmp[2] == 'gene':
			geneid = name+'_'+str(gene_num)+'_'+tmp[8].split(' ')[4]
			geneinfo = [tmp[0],label,tmp[2],tmp[3],tmp[4],'.',tmp[6],tmp[7],'ID='+geneid]
			mRNAinfo = [tmp[0],label,'mRNA',tmp[3],tmp[4],'.',tmp[6],tmp[7],'ID=transcript:'+geneid+';Parent='+geneid+';']
			print ('\t'.join(geneinfo))
			print ('\t'.join(mRNAinfo))
		elif tmp[2] == 'cds':
			cdsinfo = [tmp[0],label,'CDS',tmp[3],tmp[4],'.',tmp[6],tmp[7],'ID=CDS:'+geneid+'.'+str(cds_num)+';Parent=transcript:'+geneid+';']
			cds_num = cds_num +1
			print ('\t'.join(cdsinfo))
		elif tmp[2] == 'exon':
			exoninfo = [tmp[0],label,'exon',tmp[3],tmp[4],'.',tmp[6],tmp[7],'ID=exon:'+geneid+'.'+str(exon_num)+';Parent=transcript:'+geneid+';']
			exon_num = exon_num+1
			print ('\t'.join(exoninfo))
	gene_num = gene_num +1
