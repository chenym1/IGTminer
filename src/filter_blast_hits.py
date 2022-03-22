#!/usr/bin/env python3

import sys

chrsize = sys.argv[1]

if sys.argv[2]:
    IN = open(sys.argv[2])
else:
    IN = sys.stdin

chrdc = {}
with open(chrsize) as f:
	for i in f:
		i = i.strip().split('\t')
		chrdc[i[0]] = float(i[1])

f = IN.readlines()
for i in f:
	try:
		i = i.strip().split('\t')
		pos1 = float(i[8])
		pos2 = float(i[9])
		if float(i[2]) >= 75:
			if pos1 > pos2:
				start = pos2
				end = pos1
			else:
				start = pos1
				end = pos2
			#
			if start - 40001 < 0:
					start = 0
			else:
				start = start - 40001
			if end + 40000 > chrdc[i[1]]:
				end = chrdc[i[1]]
			else:
				end = end + 40000
			#
			start = int(start)
			end = int(end)
			print ('\t'.join([i[1],str(start),str(end)]))
	except KeyError:
		continue
if sys.argv[2]:
	IN.close()
