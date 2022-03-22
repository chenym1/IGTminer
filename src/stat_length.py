#!/usr/bin/env python3
import sys
dc = {}
with open(sys.argv[1]) as fa:
    for line in fa:
        if line.startswith('>'):
            ID = line.strip().split('>')[1].split(' ')[0]
            dc[ID] = ''
        elif not line.startswith('#'):
            dc[ID] += line.strip('\n')

for i in dc:
	ID = i
	length = len(dc[i])
	print('\t'.join([ID,str(length)]))
