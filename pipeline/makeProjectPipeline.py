#!/usr/bin/env python
import numpy as np
import os
from itertools import product

#### modify these:
n = 300
ks = np.arange(1.0, 6.5, 0.5)
pIs = np.arange(0.0, 1.05, 0.05)
gEs = np.arange(2.0, 6.0, 1.0)
gIs = np.arange(2.0, 6.0, 1.0)
reps = range(8)
tf = 200000
projName = "random_fine_g_sweep"
####

print "setting up project '" + projName + "'"
cwd = os.getcwd()
dataDir = os.path.join(os.environ['HOME'], 'work', 'prebotc', 'data', projName)
graphDir = os.path.join(dataDir, 'graphs')
oputDir = os.path.join(dataDir, 'output')
postDir = os.path.join(dataDir, 'post')
plotDir = os.path.join(dataDir, 'plots')
srcDir = os.path.join(os.environ['HOME'], 'work', 'prebotc', 'src')
paramFn = os.path.join(srcDir, 'model', 'param_files', 'BPR_syn.json')
cmdFn0 = os.path.join(srcDir, 'pipeline', projName + "_graphs")
cmdFn1 = os.path.join(srcDir, 'pipeline', projName + "_sweepfile")
cmdFn2 = os.path.join(srcDir, 'pipeline', projName + "_post")
cmdFn3 = os.path.join(srcDir, 'pipeline', projName + "_collect_control")
graphOpts = "--deg"
modelOpts = "-q --save_spikes"
postOpts = "-f 60 --bin 40"

try:
    os.makedirs(graphDir)
    os.makedirs(oputDir)
    os.makedirs(postDir)
    os.makedirs(plotDir)
except OSError:
    pass
cmd0 = open(cmdFn0, 'w')
cmd1 = open(cmdFn1, 'w')
cmd2 = open(cmdFn2, 'w')
cmd3 = open(cmdFn3, 'w')

## loop through parameters
for (rep, pI, k, gE, gI) in product(reps, pIs, ks, gEs, gIs):
    baseName = "er_n%d_k%0.1f_deg_pI%0.2f_rep%d_gE%0.1f_gI%0.1f" % (n, k, pI, 
                                                                    rep, gE, gI)
    graphFn = os.path.join(graphDir, baseName + ".gml")
    simOutFn = os.path.join(oputDir, baseName + ".mat")
    postOutFn = os.path.join(postDir, baseName + "_post.mat")
    cmd = ' '.join(
        [os.path.join(srcDir, 'graphs', 'genER.py'),
         str(n), str(k), graphFn, '-pI', str(pI), '-gE', str(gE),
         '-gI', str(gI), graphOpts]) + '\n'
    cmd0.write(cmd)
    cmd = ' '.join(
        [os.path.join(srcDir, 'model', "runmodel.py"),
         paramFn, graphFn, simOutFn, '-tf', str(tf), modelOpts]) + "\n"
    cmd1.write(cmd)
    cmd = ' '.join(
        [os.path.join(srcDir, 'postprocessing', "doPost.py"),
         simOutFn, postOutFn, postOpts]) + "\n"
    cmd2.write(cmd)
    iput = " ".join([graphFn, simOutFn, postOutFn, str(k), str(pI), str(rep),
                     str(gE), str(gI)]) + "\n"
    cmd3.write(iput)
cmd0.close()
cmd1.close()
cmd2.close()
cmd3.close()


