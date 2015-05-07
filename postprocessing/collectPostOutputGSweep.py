#!/usr/bin/env python
"""
collectPostOutputGSweep.py

collect output from different runs & reps into single data structures
to be used with the newer sweeps across gE, gI
"""

import sys
import os
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import progressbar

#### config here
projName = "random_extend"
outFn = os.path.join(os.environ['HOME'], 'work', 'prebotc', 
                      'data', projName, 'post/collected.mat')

## setup variables
srcDir = os.path.join(os.environ['HOME'], 'work', 'prebotc', 'src')
iputFn = os.path.join(srcDir, 'pipeline', projName + "_collect_control")
postFn = os.path.join(srcDir, 'pipeline', projName + "_post")
errFn = os.path.join(srcDir, 'pipeline', projName + "_post_err")
## load
f = open(iputFn, 'r')
lines = f.readlines()
splitLines = np.array([ line.split() for line in lines ])
fp = open(postFn, 'r')
postLines = fp.readlines()
## casting as numpy array allows slicing
f.close()
fErr = open(errFn, 'w')
## get the parameters for which we will collect
ks = np.array(np.unique(splitLines[:,3]), dtype=np.float)
ks = np.sort(ks)
print "ks: " + str(ks)
pIs = np.array(np.unique(splitLines[:,4]), dtype=np.float)
pIs = np.sort(pIs)
print "pIs: " + str(pIs)
gEs = np.array(np.unique(splitLines[:,6]), dtype=np.float)
gEs = np.sort(gEs)
print "gEs: " + str(gEs)
gIs = np.array(np.unique(splitLines[:,7]), dtype=np.float)
gIs = np.sort(gIs)
print "gIs: " + str(gIs)
reps = np.array(np.unique(splitLines[:,5]), dtype=np.float)
reps = np.sort(reps)
numk = len(ks)
numpI = len(pIs)
numgE = len(gEs)
numgI = len(gIs)
numRep = len(reps)
print "num k: " + str(numk)
print "num pI: " + str(numpI)
print "num gE: " + str(numgE)
print "num gI: " + str(numgI)
print "num rep: " + str(numRep)

## setup the collected arrays
chiArray = np.zeros((numk, numpI, numgE, numgI, numRep), dtype=np.float)
fMax = np.zeros((numk, numpI, numgE, numgI, numRep), dtype=np.float)
lag = np.zeros((numk, numpI, numgE, numgI, numRep), dtype=np.float)
op_angle_mean = np.zeros((numk, numpI, numgE, numgI, numRep), dtype=np.float)
op_angle_std = np.zeros((numk, numpI, numgE, numgI, numRep), dtype=np.float)
num_expir = np.zeros((numk, numpI, numgE, numgI, numRep), dtype=np.float)
avg_firing_rate = np.zeros((numk, numpI, numgE, numgI, numRep), dtype=np.float)
amplitude_irregularity = np.zeros((numk, numpI, numgE, numgI, numRep),
                                  dtype=np.float)
ibi_irregularity = np.zeros((numk, numpI, numgE, numgI, numRep),
                            dtype=np.float)

print("Looping over all postprocessing output....")
bar_updates = 100
widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()]
bar = progressbar.ProgressBar(maxval=bar_updates, widgets=widgets)
bar.start()
bar_i = 0
nstep=len(splitLines)
for i in range(len(splitLines)):
    run = splitLines[i,:]
    postFile = run[2]
    k = float(run[3])
    pI = float(run[4])
    rep = float(run[5])
    gE = float(run[6])
    gI = float(run[7])
    idx = (np.where(ks == k), np.where(pIs == pI), np.where(gEs == gE),
           np.where(gIs == gI), np.where(rep == reps))
    try:
        M = scipy.io.loadmat(postFile)
        chiArray[idx] = float(M['chi'])
        fMax[idx] = float(M['peak_freq'])
        lag[idx] = float(M['peak_lag'])
        op_angle_mean[idx] = float(M['op_angle_mean'])
        op_angle_std[idx] = float(M['op_angle_std'])
        num_expir[idx] = float(M['num_expir'])
        avg_firing_rate[idx] = float(M['avg_firing_rate'])
        amplitude_irregularity[idx] = float(M['amplitude_irregularity'])
        ibi_irregularity[idx] = float(M['ibi_irregularity'])
    except (IOError, KeyError):
        # simOutFn = run[1]
        # cmd = "./doPost.py " + simOutFn + " " + postFile + "\n"
        print postFile + " is missing"
        print "Writing command:"
        cmd = postLines[i]
        print cmd
        fErr.write(cmd)
    except:
        print "Unexpected error in " + postFile
        print "Writing command:"
        cmd = postLines[i]
        print cmd
        fErr.write(cmd)
    if ( i % np.floor(nstep/bar_updates) ) == 0:
        bar.update(bar_i)
        bar_i+=1
bar.finish()

# means over reps
chiArray_mean=np.mean(chiArray,axis=4)
fMax_mean=np.mean(fMax, axis=4)
lag_mean=np.mean(lag, axis=4)
op_angle_mean_mean=np.mean(op_angle_mean, axis=4)
op_angle_std_mean=np.mean(op_angle_std, axis=4)
num_expir_mean=np.mean(num_expir,axis=4)
avg_firing_rate_mean=np.mean(avg_firing_rate,axis=4)
amplitude_irregularity_mean=np.mean(amplitude_irregularity,axis=4)
ibi_irregularity_mean=np.mean(ibi_irregularity,axis=4)
# standard deviations over reps
chiArray_std=np.std(chiArray,axis=4)
fMax_std=np.std(fMax, axis=4)
lag_std=np.std(lag, axis=4)
op_angle_mean_std=np.std(op_angle_mean, axis=4)
op_angle_std_std=np.std(op_angle_std, axis=4)
num_expir_std=np.std(num_expir,axis=4)
avg_firing_rate_std=np.std(avg_firing_rate,axis=4)
amplitude_irregularity_std=np.std(amplitude_irregularity,axis=4)
ibi_irregularity_std=np.std(ibi_irregularity,axis=4)

X = np.transpose(np.tile(ks,(numpI,1)))
Y = np.tile(pIs,(numk,1))

Xg= np.transpose(np.tile(gEs,(numgI,1)))
Yg= np.tile(gIs,(numgE,1))

fErr.close()

scipy.io.savemat(outFn,
                 mdict={'X':X,
                        'Y':Y,
                        'Xg': Xg,
                        'Yg': Yg,
                        'chiArray':chiArray_mean,
                        'fMax': fMax_mean,
                        'lag': lag_mean,
                        'op_angle_mean': op_angle_mean_mean,
                        'op_angle_std': op_angle_std_mean,
                        'num_expir': num_expir_mean,
                        'avg_firing_rate': avg_firing_rate_mean,
                        'amplitude_irregularity':amplitude_irregularity_mean,
                        'ibi_irregularity':ibi_irregularity_mean,
                        'chiArray_std':chiArray_std,
                        'fMax_std': fMax_std,
                        'lag_std': lag_std,
                        'op_angle_mean_std': op_angle_mean_std,
                        'op_angle_std_std': op_angle_std_std,
                        'num_expir_std': num_expir_std,
                        'avg_firing_rate_std': avg_firing_rate_std,
                        'amplitude_irregularity_std':amplitude_irregularity_std,
                        'ibi_irregularity_std':ibi_irregularity_std,
                        'ks': ks,
                        'pIs': pIs,
                        'gEs': gEs,
                        'gIs': gIs,
                        'reps': reps})
