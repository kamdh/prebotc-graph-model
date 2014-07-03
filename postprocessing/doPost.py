#!/usr/bin/env python
'''
doPost.py

python version of postprocessing script for identifying spikes and filtering the
output.
'''
import sys
import numpy as np
import scipy.signal
import scipy.io
import argparse
import networkx as nx
#from matplotlib.pyplot import acorr, psd
from scikits.talkbox.tools.correlations import acorr
from IPython import embed
from sklearn.decomposition import NMF

maxorder=20
eta_norm_pts = 10

def parse_args(argv):
    # defaults
    transient = 40000 # ms
    spike_thresh = -20 # mV
    f_sigma = 20 # ms
    butter_high = 10 # Hz
    butter_low = 1 # Hz
    bin_width = 20 # ms
    cutoff = 0.5
    # parsing
    parser = argparse.ArgumentParser(prog="doPost",
                                     description=('Postprocessing of' 
                                                  ' model output'))
    parser.add_argument('sim', help='model output (.mat) file')
    parser.add_argument('output', help='output (.mat) filename')
    parser.add_argument('--transient', '-t', 
                        help='transient time, ms (default: %(default)s)', 
                        type=float, default=transient)
    parser.add_argument('--sec', '-s', action='store_true',
                        help='time units are in seconds (default: ms)')
    parser.add_argument('--volt', '-V', action='store_true',
                        help=('file contains voltage traces '
                              '(default: sparse spike trains)'))
    parser.add_argument('--thresh', 
                        help='spike threshold, mV (default: %(default)s)',
                        type=float, default=spike_thresh) 
    parser.add_argument('--fsig', '-f', 
                        help=('filter standard deviation, ms '
                        '(default: %(default)s)'),
                        type=float, default=f_sigma)
    parser.add_argument('--butter_high', 
                        help=('Butterworth filter upper cutoff frequency, Hz '
                              '(default: %(default)s)'),
                        type=float, default=butter_high)
    parser.add_argument('--butter_low', 
                        help=('Butterworth filter lower cutoff frequency, Hz '
                              '(default: %(default)s)'),
                        type=float, default=butter_low)
    parser.add_argument('--bin_width', '-b', 
                        help='bin width, ms (default: %(default)s)',
                        type=float, default=bin_width)
    parser.add_argument('--cut', '-c', 
                        help='burst cutoff parameter (default: %(default)s)',
                        type=float, default=cutoff)
    args = parser.parse_args(argv[1:])
    return args.sim, args.output, args.transient, args.sec, args.thresh, \
        args.fsig, args.butter_low, args.butter_high, args.bin_width,\
        args.cut, args.volt

def chop_transient(data, transient, dt):
    '''
    Remove a transient from the data
    '''
    firstIdx = int(np.ceil(transient / dt) - 1)
    return data[:,firstIdx:]

def find_spikes(data, threshold):
    '''
    Find spikes in voltage data by taking relative maxima
    '''
    indices = scipy.signal.argrelmax(data, axis=1) # 1st and 2nd coords of maxima
    mask = np.where(data[indices] > threshold)
    new_indices = (indices[0][mask],
                   indices[1][mask])
    spike_mat = np.zeros(np.shape(data), dtype=np.int) # dense format
    spike_mat[new_indices] = 1
    return new_indices, spike_mat

def spikes_of_neuron(spikes, neuron):
    '''
    Return time indices of spiking of a given neuron
    '''
    return spikes[1][np.where(spikes[0] == neuron)]

def spikes_filt(spike_mat, samp_freq, f_sigma, butter_freq):
    '''
    Filter the spike timeseries. Returns both neuron-by-neuron timeseries
    filtered with a gaussian kernel and the population data filtered
    with a butterworth filter.

    Parameters
    ==========
    spike_mat: the numneuron x time matrix of spikes
    samp_freq: sample frequency
    f_sigma:   variance of gaussian
    butter_freq: butterworth filter cutoff frequency(s)

    Returns
    =======
    spike_fil: gaussian filtered matrix, same shape as spike_mat
    int_signal: butterworth filtered population timeseries
    spike_fil_butter: butterworth filtered matrix, same shape as spike_mat
    '''
    def filt_window_gauss(samp_freq, std = 20, width = None, normalize = 1):
        if width is None:
            width = std*4+1
        width /= samp_freq
        std /= samp_freq
        w = scipy.signal.gaussian(width, std)
        if not normalize == 0:
            w = normalize * w / sum(w)
        return w
    def filt_gauss(spike_mat, samp_freq, f_sigma=20):
        w = filt_window_gauss(samp_freq, std=f_sigma, normalize=1)
        spike_fil = scipy.signal.fftconvolve(spike_mat, w[ np.newaxis, : ], 
                                             mode='same')
        #spike_fil = scipy.signal.convolve(spike_mat, w[ np.newaxis, : ], 
        #                                  mode='same')
        return spike_fil
    def filt_butter(data, samp_freq, butter_freq, axis=-1):
        samp_freq *= 1e-3 # ms -> s, filter defn in terms of Hz
        order = 2
        ny = 0.5 / samp_freq # Nyquist frequency
        cof = butter_freq # cutoff freq, in Hz
        cof1 = cof / ny
        b, a = scipy.signal.butter(order, cof1, btype='low')
        filtNs = scipy.signal.filtfilt(b, a, data, axis=axis)
        filtNs /= samp_freq
        return filtNs
    spike_fil = filt_gauss(spike_mat, samp_freq, f_sigma=f_sigma) 
    int_signal = filt_butter(np.mean(spike_mat, axis=0), 
                             samp_freq, butter_freq)
    spike_fil_butter = filt_butter(spike_fil, samp_freq, butter_freq, axis=1)
    return spike_fil, int_signal, spike_fil_butter

def bin_spikes(spike_mat, bin_width, dt):
    '''
    Bin spikes

    Parameters
    ==========
      spike_mat: matrix of spikes, (num_neuron x num_time)
      bin_width: bin width in time units
      dt: sampling frequency in spike mat

    Returns
    =======
      bins: an array of the bin locations in time units
      binned_spikes: a new matrix (num_neuron x num_bins)
    '''
    num_neurons= np.shape(spike_mat)[0]
    num_times = np.shape(spike_mat)[1]
    stride = int(np.ceil(bin_width / dt))
    bins = np.arange(0, num_times, stride, dtype=np.float)
    which_bins = np.digitize(range(0, num_times), bins)
    num_bins = len(bins)
    binned_spikes = np.zeros((num_neurons, num_bins), dtype=np.int)
    for i in range(num_bins):
        bin_mask = np.where(which_bins == i)[0] # mask data in bin i, tuple
        bin_data = spike_mat[:,bin_mask]
        binned_spikes[:,i] = np.sum(bin_data, axis=1).flatten()
    return bins, binned_spikes

def synchrony_stats(data, dt, maxlags=3000):
    '''
    Synchrony measures

    Parameters
    ==========
      data: numneuron x time
      maxlags: maximal lag for autocorrelation, default = 3000 ms
    
    Returns
    =======
      chi: synchrony measure
      autocorr: autocorrelation of population avg \bar{data}(t)
    '''
    data_pop = np.mean(data, axis=0) # pop avg
    sigma_pop = np.mean(np.square(data_pop)) - np.square(np.mean(data_pop))
    sigma = np.mean(np.square(data), axis=1) - np.square(np.mean(data, axis=1))
    sigma_mean = np.mean(sigma)
    chisq = sigma_pop / sigma_mean
    chi = np.sqrt(chisq)
    autocorr = acorr(data_pop - np.mean(data_pop), onesided=True, scale='coeff')
    return chi, autocorr

def peak_freq_welch(data, dt):
    '''
    Compute the Welch periodogram (psd) and return the peak frequency
    '''
    # f, Pxx = scipy.signal.welch(data, fs = 1000/dt, detrend='constant',
    #                              return_onesided = True, nperseg=2**14)
    freq, power = scipy.signal.periodogram(data, fs = 1000/dt,
                                       return_onesided=True, 
                                       detrend='constant')
    idx = np.argmax(power)
    peak_freq = freq[idx]
    peak_lag = 1/peak_freq
    return peak_lag, peak_freq, freq, power

def isi(raster):
    '''
    Finds the inter-event (spike)-interval of a raster (0-1) array,
    where 1 indicates an event.

    Parameters
    ==========
      raster: 0-1 array

    Returns
    =======
      isi_vec: vector of inter-event-intervals
    '''
    whenSpiking = np.nonzero(raster)[0]
    isi_vec = np.diff(whenSpiking)
    isi_vec = isi_vec[isi_vec != 1]
    return isi_vec

def burst_lens(raster):
    '''
    Burst lengths in a raster; 0 indicates not bursting, 1 indicated bursting.
    Counts length of consecutive strings of 1s in the raster.

    Parameters
    ==========
      raster: 0-1 array

    Returns
    =======
      runs_zeros: length of consecutive strings of 1s 
    '''
    ## apply not to raster since we count 0s
    new_raster = np.array(np.logical_not(raster), dtype=np.int)
    w = np.hstack((1, new_raster, 1)) # pad with 1s
    runs_zeros = np.nonzero(np.diff(w) == 1)[0] - np.nonzero(np.diff(w) == -1)[0]
    return runs_zeros

def burst_starts(raster):
    '''
    Find the starting points of each burst, where raster entries jump
    from 0 to 1.
    '''
    return np.where(np.diff(raster) == 1)[0] + 1

def burst_stats(data, cutoff, dt):
    '''
    Estimate when the population is bursting by comparing filtered
    activity data with a threshold = cutoff*(max(data) - min(data)) + min(data).

    Parameters
    ==========
      data: butterworth filtered signal
      cutoff: fraction of variation to define bursting
      dt: sampling period of butterworth
    
    Returns
    =======
      dutyCycle
      muIBI: mean of IBI distribution
      cvIBI: cv of IBI distribution
      muB: mean burst duration
      cvB: cv of burst durations
    '''
    if cutoff <= 0: #or cutoff > 1:
        raise Exception("cutoff out of range")
    # mindata = np.min(data[10:-10])
    # maxdata = np.max(data[10:-10])
    # thresh = cutoff * (maxdata - mindata) + mindata
    mean_data = np.mean(data[20:-20])
    std_data = np.std(data[20:-20])
    thresh = mean_data + std_data*cutoff
    bursting = np.array(data > thresh, dtype=np.float)
    duty_cycle = np.sum(bursting) / bursting.shape[0]
    ibi_vec = isi(bursting) * dt
    burst_start_locs = burst_starts(bursting)
    burst_lengths = burst_lens(bursting)
    burst_peak_locs = np.zeros(burst_start_locs.shape)
    burst_peaks = np.zeros(burst_start_locs.shape)
    bad_bursts=0
    for i in range(len(burst_start_locs)):
        # find peak of burst i
        burst_index = burst_start_locs[i] + range(burst_lengths[i])
        tmp = data[ burst_index ]
        # peakInd = np.argmax(tmp)
        peak_index = scipy.signal.argrelmax(tmp)[0]
        if len(peak_index) > 1 or not peak_index:
            ## more than one peak found or empty list
            #print "more than one local max found in burst " + str(i)
            bad_bursts += 1
            peak_index = np.argmax(tmp)
        # add 1 for Matlab 1-based indexing
        burst_peak_locs[i] = burst_index[ np.int(peak_index) ] + 1
        burst_peaks[i] = data[ burst_index[ np.int(peak_index) ] ]
    ibi_mean = np.mean(ibi_vec)
    ibi_cv = np.std(ibi_vec) / ibi_mean
    burst_length_mean = np.mean(burst_lengths * dt)
    burst_length_cv = np.std(burst_lengths * dt) / burst_length_mean
    burst_start_locs += 1 # for Matlab
    return (duty_cycle, ibi_mean, ibi_cv, burst_length_mean, burst_length_cv, 
            ibi_vec, burst_lengths, burst_start_locs, burst_peak_locs, 
            burst_peaks, bursting, bad_bursts)

def graph_attributes(graph_fn):
    '''
    Load vertex attributes so we can access them from MATLAB.
    Note that if the attribute does not exist, this will be an empty array.
    '''
    g = nx.read_gml(graph_fn)
    vertex_types = np.array(nx.get_node_attributes(g, 'type').values(),
                           dtype=np.int)
    vertex_inh = np.array(nx.get_node_attributes(g, 'inh').values(), 
                         dtype=np.int)
    vertex_respir_area = np.array(
        nx.get_node_attributes(g,'respir_area').values(),
        dtype=np.int)
    graph_adj = nx.adjacency_matrix(g, weight='gsyn')
    return vertex_types, vertex_inh, vertex_respir_area, graph_adj

def event_trig_avg(events, data, normalize=False, pts=10):
    '''
    Compute an event-triggered average.

    Parameters
    ==========
    events, ndarray
      Array of event indices.
    data, ndarray, ndim=2
      Array to be averaged along dim 1 relative to the events.
    normalize, bool, optional
      Whether to normalize to [-.5, .5]
    '''
    breakpts = np.array(
        np.hstack((0, (events[0:-1] + events[1:]) / 2., data.shape[1]-1)),
        dtype=np.int)
    if normalize:
        from scipy.interpolate import griddata
        max_interval = 2*pts
        xgrid1 = np.linspace(-.5, 0, num=pts, endpoint=False)
        xgrid2 = np.linspace(0, .5, num=pts)
    else:
        max_interval = 2*np.max(np.hstack((events-breakpts[0:-1],
                                          breakpts[1:]-events)))
    midpt = np.floor(max_interval / 2)
    numevents = events.shape[0]
    eta = np.zeros((data.shape[0], max_interval))
    for i in range(numevents):
        timeidx = range(int(breakpts[i]), int(breakpts[i+1]))
        thisevent = events[i] 
        center = np.where(timeidx==thisevent)[0].astype(int)
        if normalize:
            xs1 = (np.array(timeidx[:center], dtype=np.float) 
                   - timeidx[center])
            xs1 /= xs1[0]*(-2)
            xs2 = (np.array(timeidx[center:], dtype=np.float) 
                   - timeidx[center])
            xs2 /= xs2[-1]*2
            fun1 = lambda x: scipy.interpolate.griddata(xs1, x, xgrid1)
            fun2 = lambda x: scipy.interpolate.griddata(xs2, x, xgrid2)
            toadd =np.hstack((
                np.apply_along_axis(fun1,
                                    1, data[:,timeidx[:center]]),
                np.apply_along_axis(fun2,
                                    1, data[:,timeidx[center:]])))
            eta += toadd
        else:
            lpad = midpt - center
            rpad = max_interval - (len(timeidx)+lpad)
            eta += np.lib.pad(data[:, timeidx], ((0,0), (lpad,rpad)), 
                              'constant', constant_values=(0,0))
    eta /= float(numevents)
    return eta

def nmf_error(eta):
    '''
    Decompose the ETAs using nonnegative matrix factorization with 1 component.
    The reconstruction error is an estimate of the non-inspiratory activity.
    '''
    nmf = NMF(n_components=1)
    nmf.fit(eta)
    reconstruction_err = nmf.reconstruction_err_
    return float(reconstruction_err)

def main(argv = None):
    if argv is None:
        argv = sys.argv
    (simFn, outFn, trans, sec_flag, spike_thresh, f_sigma, butter_low, 
     butter_high, bin_width, cutoff, are_volts) = parse_args(argv)
    butter_freq = np.array([butter_low, butter_high])
    if sec_flag:
        scalet = 1e3
    else:
        scalet = 1
    ## load simulation output
    ## assumes no --save_full
    sim_output = scipy.io.loadmat(simFn)
    ## begin postprocessing
    graph_fn = str(sim_output['graphFn'][0])
    dt = float(sim_output['dt']) * scalet
    data = chop_transient(sim_output['Y'], trans, dt)
    num_neurons = np.shape(data)[0]
    tMax = np.shape(data)[1]
    ## spike trains
    if are_volts:
        spikes, spike_mat = find_spikes(data, spike_thresh)
    else:
        spike_mat = data.todense()
        spikes = data.nonzero()
    bins, spike_mat_bin = bin_spikes(spike_mat, bin_width, dt)
    ## filtering
    # spike_fil, butter_int = spikes_filt(spike_mat, dt, f_sigma, butter_freq)
    spike_fil_bin, butter_int_bin, spike_fil_butter = spikes_filt(spike_mat_bin, 
                                                                  dt*bin_width, 
                                                                  f_sigma, 
                                                                  butter_freq)
    psth_bin = np.sum(spike_mat_bin, axis=0)
    ## synchrony measures
    if are_volts:
        chi, autocorr = synchrony_stats(data, dt)
    else:
        chi, autocorr = synchrony_stats(spike_fil_bin, dt*bin_width)
    peak_lag, peak_freq, freq, power = peak_freq_welch(psth_bin, dt*bin_width)
    ## burst measures
    (duty_cycle, ibi_mean, ibi_cv, burst_length_mean, burst_length_cv, ibi_vec,
     burst_lengths, burst_start_locs, burst_peak_locs, burst_peaks, bursting, 
     bad_bursts) = burst_stats(butter_int_bin, cutoff, dt*bin_width)
    pop_burst_peak = scipy.signal.argrelmax(butter_int_bin, order=maxorder)[0]
    # plt.plot(bins, butter_int_bin)
    # plt.plot(bins[which_bins_max], butter_int_bin[which_bins_max], 'r*')
    # eta = event_trig_avg(pop_burst_peak, spike_mat_bin)
    eta = event_trig_avg(pop_burst_peak, spike_fil_bin)
    eta_t = (np.arange(eta.shape[1])-eta.shape[1]/2)*bin_width
    eta_norm = event_trig_avg(pop_burst_peak, spike_fil_bin, normalize=True,
                              pts=eta_norm_pts)
    eta_t_norm = np.linspace(-0.5, 0.5, 2*eta_norm_pts)
    eta[eta < 0] = 0
    eta_norm[eta_norm < 0] = 0
    eta_nmf_err = nmf_error(eta)
    eta_norm_nmf_err = nmf_error(eta_norm)
    ## graph data
    graph_fn='../../data/random_fine_2/graphs/er_n300_k2.0_deg_pI0.20_rep0.gml'
    (vertex_types, vertex_inh, 
     vertex_respir_area, graph_adj) = graph_attributes(graph_fn)
    embed()
    ## save output
    scipy.io.savemat(outFn,
                      mdict = {'bins': bins,
                               'spike_mat_bin': spike_mat_bin,
                               # 'thetas': thetas,
                               # 'spike_fil': spike_fil,
                               # 'butter_int': butter_int,
                               'spike_fil_bin': spike_fil_bin,
                               'butter_int_bin': butter_int_bin,
                               'psth_bin': psth_bin,
                               'chi': chi,
                               'autocorr': autocorr,
                               'peak_lag': peak_lag,
                               'peak_freq': peak_freq,
                               'duty_cycle': duty_cycle,
                               'ibi_mean': ibi_mean,
                               'ibi_cv': ibi_cv,
                               'burst_length_mean': burst_length_mean,
                               'burst_length_cv': burst_length_cv,
                               'ibi_vec': ibi_vec,
                               'burst_lengths': burst_lengths,
                               'burst_start_locs': burst_start_locs,
                               'burst_peak_locs': burst_peak_locs,
                               'burst_peaks': burst_peaks,
                               'bursting': bursting,
                               'bad_bursts': bad_bursts,
                               'eta': eta,
                               'eta_t': eta_t,
                               'eta_norm': eta_norm,
                               'eta_t_norm': eta_t_norm,
                               'do_post_argv': " ".join(argv),
                               'vertex_types': vertex_types,
                               'vertex_inh': vertex_inh,
                               'vertex_respir_area': vertex_respir_area,
                               'graph_adj': graph_adj
                               },
                      oned_as='column'
                     )

# run the main stuff
if __name__ == '__main__':
    main()

