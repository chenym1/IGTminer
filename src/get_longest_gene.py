#!/usr/bin/env python3

import re
def longestfasta(infile, gff, name):
	
	dc = {}
	with open(infile) as fa:
		for line in fa:
			if line.startswith('>'):
				ID = line.strip().split('>')[1].split(' ')[0]
				dc[ID] = ''
			else:
				dc[ID] += line
	
	# transcript ID to gene ID
	dc_t2g = {}
	with open(gff) as FILE:
		for j in FILE:
			j = j.strip().split('\t')
			try:
				if j[2] == 'mRNA' or j[2] == 'transcript':
					info = j[8].split(';')
					for info_num in range(len(info)):
						info_tmp = info[info_num]
						if info_tmp.startswith('ID'):
							tid = info_tmp.split('=')[1]
						elif info_tmp.startswith('Parent') or info_tmp.startswith('geneID'):
							gid = info_tmp.split('=')[1]
					dc_t2g[tid] = gid
			except IndexError:
				continue
	# longest seq
	longest_dc = {}
	for ID in dc_t2g.keys():
		try:
			seq = dc[ID]
			gene_ID = dc_t2g[ID]
			if gene_ID not in longest_dc:
				longest_dc[gene_ID] = [ID,seq]
			else:
				length = len(seq)
				if length > len(longest_dc[gene_ID][1]):
					longest_dc[gene_ID] = [ID,seq]
		except KeyError:
			continue

	longest_t2g = {}
	out_seq = open(name,'w')
	for i in longest_dc.keys():
		longest_t2g[longest_dc[i][0]] = i
		seqname = i.replace('gene:','').replace('transcript:','')
		out_seq.write('>'+seqname+'\n')
		out_seq.write(longest_dc[i][1])
	
#
from optparse import OptionParser
# ===========================================
def main():
    usage = "Usage: %prog -i pep.fa -f gff\n" \
            "Description: find the longest protein sequence in genome"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file", metavar="FILE")
    parser.add_option("-f", dest="gff",
		  help="gff", metavar="FILE")
    parser.add_option("-n", dest="name",
                  help="prefix name", metavar="FILE")
    (options, args) = parser.parse_args()
    longestfasta (options.infile,options.gff,options.name)


# ===========================================
if __name__ == "__main__":
    main()
