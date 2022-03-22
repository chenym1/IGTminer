#!/usr/bin/env python3

import sys

spec_list = []
with open(sys.argv[1]) as f:
	for i in f:
		i = i.strip()
		spec_list.append(i)

with open(sys.argv[2]) as f:
	for i in f:
		i = i.strip().split('\t')
		out = []
		for j in range(len(spec_list)):
			spec = spec_list[j]
			matchs = [tmp.split('.')[1] for tmp in i if tmp.startswith(spec)]
			if len(matchs)>0:
				out.append(','.join(matchs))
			else:
				out.append('.')
		print ('\t'.join(out))
