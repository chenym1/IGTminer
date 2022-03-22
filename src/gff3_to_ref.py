#!/usr/bin/env python3
import sys
with open(sys.argv[1]) as f:
	for i in f:
		i = i.strip().split('\t')
		info = i[0].split('__')
		CHR = info[0]
		start_ref = int(info[1])
		start = int(i[3])
		end = int(i[4])
		i[0] = CHR
		i[3] = str(start_ref+start-1)
		i[4] = str(start_ref+end-1)
		print('\t'.join(i))
