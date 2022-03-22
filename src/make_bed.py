#!/usr/bin/env python3
import sys
import numpy as np

dc = {}
name = sys.argv[1].split('/')
name = name[len(name)-1].split('.')[0]
with open(sys.argv[1]) as f:
	for i in f:
		i = i.strip().split('\t')
		dc.setdefault(i[4],[]).append([i[0],i[1],i[2]])
#
for i in dc:
	info = dc[i]
	CHR = ''
	start = []
	end = []
	for j in range(len(info)):
		CHR = info[j][0]
		start.append(int(info[j][1]))
		end.append(int(info[j][2]))
	left = str(int(np.median(start)))
	right = str(int(np.median(end)))
	print('\t'.join([CHR,left,right,name+'.'+i,'.','+']))
