#!/usr/bin/env python3
import re
#
def Filter (bed,matchchr):
	###
	with open(bed) as FILE:
		for i in FILE:
			i = i.strip().split('\t')
			chr1 = i[0]
			chra = re.sub('[0-9][0-9]|[0-9]','N',chr1)
			if re.search(matchchr,chra):
				print ('\t'.join(i))
###
from optparse import OptionParser
def main():
        usage = "Usage: %prog -i blast_A_to_B -l A.bed -m chrNA\n" \
                "Description: Match hit by same chromosome\n"
        parser = OptionParser(usage)
        parser.add_option("-l", dest="bed",
                  help="A.bed", metavar="FILE")
        parser.add_option("-m", dest="matchchr",
                  help="chrNA", metavar="FILE")
        (options, args) = parser.parse_args()
        Filter(options.bed,options.matchchr)
###
if __name__ == "__main__":
        main()
