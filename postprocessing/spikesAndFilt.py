#!/usr/bin/env python
#
# spikesAndFilt.py
#
# python version of postprocessing script for identifying spikes and
# filtering the output.
#


import os
import sys
import numpy as np
import scipy
import scipy.signal
import scipy.io
import argparse



def parse_args(argv):
    # defaults
    transient = 20000 # ms
    spikeThresh = -20 # mV
    filtWin = 20 # ms
    # parsing
    parser = argparse.ArgumentParser(prog="spikesAndFilt",
                                     description='filter model output')
    parser.add_argument('sim', help='model output (.mat) file')
    parser.add_argument('output', help='output (.mat) filename')
    parser.add_argument('--transient', '-t', help='transient time, ms (default: %(default))', 
                        type=float, default=transient)
    parser.add_argument('--sec', '-s', action='store_true',
                    help='time units are in seconds (default: ms)')
    parser.add_argument('--thresh', '-th', help='spike threshold, mV (default: %(default))',
                        type=float, default=spikeThresh)
    parser.add_argument('--fwin', '-f', help='filter window, ms (default: %(default))',
                        type=float, default=filtWin)
    args = parser.parse_args(argv[1:])
    return args.sim, args.output, args.transient, args.sec, args.thresh, args.fwin

def chop_transient(Y, transient, dt, scalet):
    firstIdx = int( np.ceil( transient / (dt*scalet) ) )
    return Y[:,firstIdx:]

def find_spikes(Y, threshold):
    indices = scipy.signal.argrelmax(Y, axis=1) # list w/ 1st and 2nd coords of maxima
    mask = np.where( Y[indices] > threshold )
    newIndices = ( indices[0][mask],
                   indices[1][mask] )
    return newIndices

def spikes_of_neuron(spikes, neuron):
    return spikes[1][ np.where( spikes[0] == neuron ) ]

def filt_window(std = 20, width = None, normalize = 1):
    if width is None:
        width = std*4+1
    scipy.signal.gaussian(width, std)
    w = scipy.signal.gaussian(width, std)
    if not normalize == 0:
        w = normalize * w / sum(w)
    return w

def main(argv = None):
    if argv is None:
        argv = sys.argv
    simFn, outFn, trans, secFlag, spikeThresh, fWidth = parse_args(argv)
    if secFlag:
        scalet = 1e-3
    else:
        scalet = 1
    # this will load our simulation output
    # X contains Y, vTypes, dt, t0, tf, paramFn, graphFn, absErr, relErr
    # assumes no --save_full
    X = scipy.io.loadmat(simFn)
    V = chop_transient( X['Y'], trans, X['dt'], scalet )
    numNeur = np.shape(V)[0]
    tMax = np.shape(V)[1]
    spikes = find_spikes( V, spikeThresh ) # sparse format (coords of spikes)
    # filter spikes with gaussian window
    spikeMat = np.zeros( np.shape( V ) ) # dense format
    spikeMat[spikes] = 1
    w = filt_window()
    spikeFil = scipy.signal.fftconvolve( spikeMat, w[ np.newaxis, : ]  ) # fft is fast
    
    return spikeFil, spikeMat

# run the main stuff
if __name__ == '__main__':
    status = main()
    sys.exit(status)
