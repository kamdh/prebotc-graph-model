#!/usr/bin/env python
from subprocess import call
import genER

n = 300
k = 3
numgraphs = 10


for i in range(numgraphs):
    fn = "testER_%d.gml" % i
    #args = [n
    output = genER.main(['', str(n), str(k), fn, '-d'])
