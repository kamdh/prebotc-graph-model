#!/usr/bin/env python
import numpy as np
import genER
import os

n = 300
ks = np.arange(1, 4.5, 0.5)
pIs = np.arange(0, 0.31, 0.1)
numReplicas = range(8)
baseDir = os.path.abspath("sweep_random")
print "base directory: " + baseDir
paramFn = os.path.abspath("../model/param_files/model2.json")
print "param filename: " + paramFn
oputDir = os.path.abspath( os.path.join( os.environ['HOME'],
                                         'prebotc', 'data',
                                         'sweep_random' ) )
print "output directory: " + oputDir
cmdFn = os.path.join(baseDir, "sweepfile")
print "sweep filename: " + cmdFn
tf = 180000

try:
    os.makedirs(baseDir)
except OSError:
    pass
try:
    os.remove(cmdFn)
except OSError:
    pass
cmdFile = open(cmdFn, 'w')

## loop through parameters
for k in ks:
    for pI in pIs:
        for rep in numReplicas:
            baseName = "er_n%d_k%0.1f_deg_pI%0.2f_rep%d" % (n, k, pI, rep)
            graphFn = os.path.join( baseDir, baseName + ".gml" )
            print graphFn
            #output = genER.main(['', str(n), str(k), graphFn, '--deg', '-pI', str(pI)])
            simOutFn = os.path.join( oputDir, baseName + ".mat" )
            print simOutFn
            cmd = "./runmodel.py " + paramFn + " " + graphFn + \
                " " + simOutFn + " -tf " + str(tf) + " -q\n"
            cmdFile.write(cmd)
cmdFile.close()
