#!/usr/bin/env python

import sys

dc = {}
with open(sys.argv[1]) as f:
    for i in f:
        dc[i.strip()] = ''

with open(sys.argv[2]) as f:
    for i in f:
        i = i.strip().split('\t')
        if i[0] not in dc and i[1] not in dc:
            print ('\t'.join(i))