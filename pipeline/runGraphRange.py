#!/usr/bin/env python
import numpy as np
import subprocess

n = 300
k = 2
#pIs = np.arange(0, 0.3, 0.01)
pIs = [0.0, 0.14, 0.29]

for pI in pIs:
    graphfn = "graphs/er_300_2_deg_%0.2f.gml" % pI
    outputfn = "model/output/er_300_2_deg_%0.2f_b" % pI
    print graphfn
    subprocess.call(["model/runmodel.py", "model/param_files/model2_b.json",
                     graphfn, outputfn])
    
