#!/usr/bin/env python3
import re
#
def Filter (blastfile,beda,bedb,matchpair):
	# match pair
	matchpair = matchpair.split(',')
	matchA = matchpair[0]
	matchB = matchpair[1]
	## A chr pos
	with open(beda) as ainfo:
		adc = {}
		for i in ainfo:
			i = i.strip().split('\t')
			gene = i[3]
			CHR = i[0]
			if gene not in adc:
				adc[gene] = CHR
	### B chr pos
	with open(bedb) as binfo:
		bdc = {}
		for i in binfo:
			i = i.strip().split('\t')
			gene = i[3]
			CHR = i[0]
			if gene not in bdc:
				bdc[gene] = CHR
	###
	with open(blastfile) as blast:
		for i in blast:
			try:
				i = i.strip().split('\t')
				chra1 = adc[i[0]]
				chrb1 = bdc[i[1]]
				chra = re.sub('[0-9][0-9]|[0-9]','N',chra1)
				chrb = re.sub('[0-9][0-9]|[0-9]','N',chrb1)
				if re.search(matchA,chra) and re.search(matchB,chrb):
					print ('\t'.join(i))
			except KeyError:
				continue
###
from optparse import OptionParser
def main():
        usage = "Usage: %prog -i blast_A_to_B -l A.bed -f B.bed -m chrNA,TuN\n" \
                "Description: Match hit by same chromosome\n"
        parser = OptionParser(usage)
        parser.add_option("-i", dest="blastfile",
                  help="blast_A_to_B", metavar="FILE")
        parser.add_option("-l", dest="beda",
                  help="A.bed", metavar="FILE")
        parser.add_option("-f", dest="bedb",
                  help="B.bed", metavar="FILE")
        parser.add_option("-m", dest="matchpair",
                  help="chrNA,TuN", metavar="FILE")
        (options, args) = parser.parse_args()
        Filter(options.blastfile,options.beda,options.bedb,options.matchpair)
###
if __name__ == "__main__":
        main()
