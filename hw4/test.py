#!/usr/bin/python

import sys

result = True
old = -float('inf')
now = 0
for i in xrange(int(sys.argv[1])):
	with open('output%02d.txt'%i,'r') as f:
		for l in f:
			now = int(l.strip())
			if now<old:
				result = False
			old =  now
if result:
	print 'TRUE'
else:
	print 'FALSE'
