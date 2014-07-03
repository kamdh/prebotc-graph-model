#!/usr/bin/env python
#
# spikesAndFilt.py
#
# python version of postprocessing script for identifying spikes and
# filtering the output.
#

import sys
import numpy as np
import scipy.signal
import scipy.io
import argparse
#from matplotlib.pyplot import acorr, psd
from scikits.talkbox.tools.correlations import acorr

def parse_args(argv):
    # defaults
    transient = 40000 # ms
    spikeThresh = -20 # mV
    fSigma = 20 # ms
    binWidth = 40 # ms
    cutoff = 0.6
    # parsing
    parser = argparse.ArgumentParser(prog="doPost",
                                     description='do postprocessing of model output')
    parser.add_argument('sim', help='model output (.mat) file')
    parser.add_argument('output', help='output (.mat) filename')
    parser.add_argument('--transient', '-t', help='transient time, ms (default: %(default))', 
                        type=float, default=transient)
    parser.add_argument('--sec', '-s', action='store_true',
                        help='time units are in seconds (default: ms)')
    parser.add_argument('--thresh', '-th', 
                        help='spike threshold, mV (default: %(default))',
                        type=float, default=spikeThresh)
    parser.add_argument('--fsig', '-f', 
                        help='filter standard deviation, ms (default: %(default))',
                        type=float, default=fSigma)
    parser.add_argument('--bin', '-b', 
                        help='bin width, ms (default: %(default))',
                        type=float, default=binWidth)
    parser.add_argument('--cut', '-c', 
                        help='burst cutoff parameter (default: %(default))',
                        type=float, default=cutoff)
    args = parser.parse_args(argv[1:])
    return args.sim, args.output, args.transient, args.sec, args.thresh, \
        args.fsig, args.bin, args.cut

def chop_transient(Y, transient, dt):
    firstIdx = int( np.ceil( transient / dt ) - 1 )
    return Y[:,firstIdx:]

def find_spikes(Y, threshold):
    indices = scipy.signal.argrelmax(Y, axis=1) # list w/ 1st and 2nd coords of maxima
    mask = np.where( Y[indices] > threshold )
    newIndices = ( indices[0][mask],
                   indices[1][mask] )
    spikeMat = np.zeros( np.shape( Y ), dtype=np.int ) # dense format
    spikeMat[newIndices] = 1
    return newIndices, spikeMat

def spikes_of_neuron(spikes, neuron):
    return spikes[1][ np.where( spikes[0] == neuron ) ]

def spikes_filt(spikeMat, sampFreq, fSigma):
    """
    Filter the spike timeseries. Returns both neuron-by-neuron timeseries
    filtered with a gaussian kernel and the population data filtered
    with a butterworth filter.

    Parameters
    ----------
    spikeMat: the numneuron x time matrix of spikes
    sampFreq: sample frequency
    fSigma:   variance of gaussian

    Returns
    -------
    F: gaussian filtered matrix same shape as spikeMat
    B: butterworth filtered population timeseries
    """
    def filt_window_gauss(sampFreq, std = 20, width = None, normalize = 1):
        if width is None:
            width = std*4+1
        width /= sampFreq
        std /= sampFreq
        scipy.signal.gaussian(width, std)
        w = scipy.signal.gaussian(width, std)
        if not normalize == 0:
            w = normalize * w / sum(w)
        return w
    def filt_gauss(spikeMat, sampFreq, fSigma=20):
        w = filt_window_gauss(sampFreq, std=fSigma, normalize=0)
        spikeFil = scipy.signal.fftconvolve( spikeMat, w[ np.newaxis, : ], mode='same' )
        return spikeFil
    def filt_butter(Y, sampFreq):
        sampFreq *= 1e-3 # ms -> s, filter defn in terms of Hz
        order = 2
        Ny = 0.5 / sampFreq # Nyquist frequency
        coF = 4
        coF1 = coF / Ny
        b, a = scipy.signal.butter(order, coF1, btype='low')
        filtNs = scipy.signal.filtfilt(b, a, Y)
        return filtNs
    spikeFil = filt_gauss(spikeMat, sampFreq, fSigma=fSigma)
    integratedSignal = filt_butter( np.sum(spikeMat, axis=0), sampFreq )
    return spikeFil, integratedSignal

def bin_spikes(spikeMat, binWidth, dt):
    numNeurons= np.shape(spikeMat)[0]
    dataPts = np.shape(spikeMat)[1]
    stride = int( np.ceil( binWidth / dt ) )
    bins = range(0, dataPts, stride)
    whichBins = np.digitize( range(0, dataPts), bins )
    numBins = len( bins )
    binnedSpikes = np.zeros( (numNeurons, numBins), dtype=np.int )
    for i in range(numBins):
        binMask = np.where( whichBins == i )[0] # mask of data belonging to bin i, tuple
        binData = spikeMat[:, binMask]
        binnedSpikes[:,i] = np.sum(binData, axis=1).flatten()
    return bins, binnedSpikes

def synchrony_stats(V, dt, maxlags=3000):
    """
    Synchrony measures

    Parameters
    ----------
        V: numneuron x time
        maxlags: maximal lag for autocorrelation, default = 3000 ms
    
    Returns
    -------
        chi: synchrony measure
        autocorr: autocorrelation of population avg \bar{V}(t)
    """
    Vpop = np.mean( V, axis=0 ) # pop avg
    sigmaV = np.mean( np.square(Vpop) ) - np.square( np.mean(Vpop) )
    sigmaVi = np.mean( np.square(V), axis=1 ) - np.square( np.mean(V, axis=1) )
    sigmaViMean = np.mean( sigmaVi )
    chisq = sigmaV / sigmaViMean
    chi = np.sqrt( chisq )
    autocorr = acorr( Vpop - np.mean(Vpop), onesided=True, scale='coeff' )
    return chi, autocorr

def peak_freq_welch(B, dt):
    """
    Compute the Welch periodogram (psd) and return the peak frequency
    """
    f, Pxx = scipy.signal.welch( B - np.mean(B), fs = 1000/dt,
                                 return_onesided = True )
    idx = np.argmax(Pxx)
    fMax = f[idx]
    lag = 1/fMax
    return lag, fMax, f, Pxx

def ISI(raster):
    whenSpiking = np.nonzero( raster )[0]
    ISIs = np.diff( whenSpiking )
    return ISIs

def burst_lens(raster):
    ## apply not to raster since we count 0s
    newraster = np.array( np.logical_not(raster), dtype=np.int )
    w = np.hstack( (1, newraster, 1) ) # pad with 1s
    runs_zeros = np.nonzero( np.diff(w) == 1 )[0] - \
        np.nonzero( np.diff(w) == -1 )[0]
    return runs_zeros

def burst_stats(B, cutoff, dt):
    """
    Estimate when the population is bursting by comparing filtered
    activity B with a threshold = cutoff*(max(B) - min(B)) + min(B).

    Parameters
    ----------
        B: butterworth filtered signal
        cutoff: fraction of variation to define bursting
    
    Returns
    -------
        dutyCycle
        muIBI: mean of IBI distribution
        cvIBI: cv of IBI distribution
        muB: mean burst duration
        cvB: cv of burst durations
    """
    if cutoff <= 0 or cutoff > 1:
        raise Exception("cutoff out of range")
    minB = np.min(B[10:-10])
    maxB = np.max(B[10:-10])
    thresh = cutoff * (maxB - minB) + minB
    bursting = B > thresh
    dutyCycle = np.float( np.sum( bursting ) ) / bursting.shape[0]
    IBIs = ISI(bursting) * dt
    burstLengths = burst_lens(bursting) * dt
    muIBI = np.mean( IBIs )
    cvIBI = np.std( IBIs ) / muIBI
    muB = np.mean( burstLengths )
    cvB = np.std( burstLengths ) / muB
    return dutyCycle, muIBI, cvIBI, muB, cvB, IBIs, burstLengths

def main(argv = None):
    if argv is None:
        argv = sys.argv
    simFn, outFn, trans, secFlag, spikeThresh, fSigma, binWidth, cutoff\
        = parse_args(argv)
    if secFlag:
        scalet = 1e3
    else:
        scalet = 1
    ## load simulation output
    ## X contains Y, vTypes, dt, t0, tf, paramFn, graphFn, absErr, relErr
    ## assumes no --save_full
    X = scipy.io.loadmat(simFn)
    ## begin postprocessing
    dt = float( X['dt'] ) * scalet
    V = chop_transient( X['Y'], trans, dt )
    numNeur = np.shape(V)[0]
    tMax = np.shape(V)[1]
    ## identifying spikes
    spikes, spikeMat = find_spikes( V, spikeThresh )
    bins, spikeMatBin = bin_spikes( spikeMat, binWidth, dt )
    ## filtering
    spikeFil, butterInt = spikes_filt( spikeMat, dt, fSigma )
    spikeFilBin, butterIntBin = spikes_filt( spikeMatBin, dt*binWidth, fSigma )
    psthBin = np.sum( spikeMatBin, axis=0 )
    ## synchrony measures
    chi, autocorr = synchrony_stats(V, dt)
    lag, fMax, f, Pxx = peak_freq_welch( butterIntBin, dt*binWidth )
    ## burst measures
    dutyCycle, muIBI, cvIBI, muB, cvB, IBIs, burstLengths = \
        burst_stats( butterIntBin, cutoff, dt*binWidth )
    ## save output
    scipy.io.savemat( outFn,
                      mdict = {'V': V,
                               'spikes': spikes,
                               #'spikeMat': spikeMat,
                               'bins': bins,
                               'spikeMatBin': spikeMatBin,
                               #'spikeFil': spikeFil,
                               #'butterInt': butterInt,
                               'spikeFilBin': spikeFilBin,
                               'butterIntBin': butterIntBin,
                               'psthBin': psthBin,
                               'chi': chi,
                               'autocorr': autocorr,
                               'lag': lag,
                               'fMax': fMax,
                               'dutyCycle': dutyCycle,
                               'muIBI': muIBI,
                               'cvIBI': cvIBI,
                               'muB': muB,
                               'cvB': cvB,
                               'IBIs': IBIs,
                               'burstLengths': burstLengths
                               } 
                      )

    

# run the main stuff
if __name__ == '__main__':
    status = main()
    sys.exit(status)
