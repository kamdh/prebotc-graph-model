#!/usr/bin/env python
import numpy as np
import os
import re

projName = "g_sweep_fix_all"
#projName = "random_fine_2"
srcDir = os.path.join(os.environ['HOME'], 'work', 'prebotc', 'src')

control_fn = projName + "_collect_control"
cmd_fn = projName + "_modify_gs"
f = open(control_fn,"r")
fout = open(cmd_fn,"w")
for line in f:
    graphFn, simOutFn, postOutFn, k, pI, rep, gE, gI = line.split(" ")
    k = float(k)
    pI = float(pI)
    rep = int(rep)
    gE = float(gE)
    gI = float(gI)
    tmp = re.sub(r'gE\d+\.\d+','gE0.0',graphFn)
    baseGraphFn = re.sub(r'gI\d+\.\d+','gI0.0',tmp)
    if (gE != 0.0) or (gI != 0.0):
        newline =  ' '.join(
            [os.path.join(srcDir, 'graphs', 'modify_g.py'),
             baseGraphFn, str(gE), str(gI), graphFn]) + '\n'
        fout.write(newline)
fout.close()
f.close()
