#!/usr/bin/env python
'''
This file does all the post processing for a given mat file at once. This includes: 

1) Deleting Transient
2)Binning the spikes 
3)Filter spikes using gaussian distribution 
4)Using butterworth filter to remove high frequency signals to smooth 
5)Finds Phase Lag and Population Correlation
'''
import sys
import numpy as np
import scipy.signal
import scipy.io
import argparse
import networkx as nx
import matplotlib.pyplot as plt

maxorder=20
eta_norm_pts = 10


def parse_args(argv):
    # defaults
    transient = 10000 # ms
    spike_thresh = -20 # mV
    f_sigma = 20 # ms
    butter_high = 4 # Hz
    butter_low = -np.inf # Hz
    bin_width = 20 # ms
    cutoff = 0.5
    # parsing
    parser = argparse.ArgumentParser(prog="doPost",
                                     description=('Postprocessing of' 
                                                  ' model output'))
    parser.add_argument('sim', help='model output (.mat) file')
    parser.add_argument('output', help='output (.jpg) filename')
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


'''
    This method chops of the transient stage of the data for better processing

    parameters: data-Data being passed in to chop
                transient- time and which you want to chop till
                dt-the change in time of the model

     return: The modified data excluding transient stage 
'''
def chop_transient(data, transient, dt):
    
    firstIdx = int(np.ceil(transient / dt) - 1)
    return data[:,firstIdx:]

'''
    Find spikes in voltage data by taking relative maxima

    parameters: data- self-explanatory
                threshhold- The voltage at which you start to count data as a spike
    return: new_indices-location of the maxima
            spike_mat-dense matrix containing 1 or 0 based on if a spike is present
'''
def find_spikes(data, threshold):
    
    indices = scipy.signal.argrelmax(data, axis=1) # 1st and 2nd coords of maxima
    mask = np.where(data[indices] > threshold)
    new_indices = (indices[0][mask],
                   indices[1][mask])
    spike_mat = np.zeros(np.shape(data), dtype=np.int) # dense format
    spike_mat[new_indices] = 1
    return new_indices, spike_mat

'''
    Return time indices of spiking of a given neuron
'''
def spikes_of_neuron(spikes, neuron):
    
    return spikes[1][np.where(spikes[0] == neuron)]

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

def spikes_filt(spike_mat, samp_freq, f_sigma, butter_freq):
    '''
    Filter the spike timeseries. Returns both neuron-by-neuron timeseries
    filtered with a gaussian kernel and the population data filtered
    with a butterworth filter.

    Parameters
    ==========
    spike_mat: the numneuron x time matrix of spikes
    samp_freq: period (in ms) between measurements in spike_mat
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
        '''
        Filter data with a 2nd order butterworth filter.
        
        Parameters
        ==========
          data: ndarray
          samp_freq: sampling period (s)
          butter_freq: [cutoff_low, cutoff_high] (Hz), can be infinite
          axis (optional): axis along which to filter, default = -1
        Returns
        =======
          filtNs: filtered version of data
        '''
        order = 2
        ny = 0.5 / samp_freq # Nyquist frequency
        cof = butter_freq / ny # normalized cutoff freq
        if np.isneginf(cof[0]) and np.isfinite(cof[1]):
            # lowpass
            cof1 = cof[1]
            b, a = scipy.signal.butter(order, cof1, btype='low')
            filtNs = scipy.signal.filtfilt(b, a, data, axis=axis)
        elif np.isfinite(cof[0]) and np.isinf(cof[1]):
            # highpass
            cof1 = cof[0]
            b, a = scipy.signal.butter(order, cof1, btype='high')
            filtNs = scipy.signal.filtfilt(b, a, data, axis=axis)
        elif np.isfinite(cof[0]) and np.isfinite(cof[1]):
            # bandpass
            b, a = scipy.signal.butter(order, cof, btype='band')
            filtNs = scipy.signal.filtfilt(b, a, data, axis=axis)
        else:
            raise Exception('filt_butter called with bad cutoff frequency')
        filtNs /= samp_freq # normalize to rate
        return filtNs
    spike_fil = filt_gauss(spike_mat, samp_freq, f_sigma=f_sigma) 
    int_signal = filt_butter(np.mean(spike_mat, axis=0), 
                             samp_freq*1e-3, butter_freq)
    spike_fil_butter = filt_butter(spike_fil, samp_freq*1e-3, 
                                   butter_freq, axis=1)
    return spike_fil, int_signal, spike_fil_butter

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
def bin_spikes(spike_mat, bin_width, dt):
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


##Implementation of the mean to make sure all code is running correctly
def find_mean(data):
    sum = 0;
    for x in data:
        sum = sum + x
    return sum/len(data)


'''def xcorr(signal_1,signal_2,taolen):
    xcorrl = np.empty(len(signal_1)/taolen)
    index = 0
    
    signal_1 = np.asarray(signal_1)
    signal_2 = np.asarray(signal_2)	
    ## m1,m2 are the means of signals 1 and 2, respectivly
    m1 = np.mean(signal_1)
    m2 = np.mean(signal_2)
    
    #loops through all values of tao 
    for tao in range(0,(len(signal_1)/taolen)):
        sum = 0
        if tao == 0:
            average = np.mean((signal_1 - m1)*(signal_2-m2))
        else:
            average = np.mean((signal_1[0:-tao] - m1)*(signal_2[tao:] - m2))
        xcorrl[index] = average
        index = index + 1
    return xcorrl
'''


'''
    This is computes the cross correlation for two signals
    
    paramters:
    signal_1: first signal you want to use
    signal_2: second signal you want to use
    taolen: this number determines how much of the tao to use
    returns:
    values of the cross correlation
'''
def xcorr(signal_1,signal_2):
    
    signal_1 = np.asarray(signal_1)
    signal_2 = np.asarray(signal_2)

    #Centering the data, giving it a zero mean to reduce variance and not worry about peak differences
    m1 = np.mean(signal_1)
    m2 = np.mean(signal_2)
    
    signal_1_centered = (signal_1 - m1) / (np.std(signal_1) * len(signal_1))
    signal_2_centered = (signal_2 - m2) / np.std(signal_2)

    xcorr = scipy.signal.correlate(signal_1_centered,signal_2_centered)
    return xcorr

'''
    Gets info from the graph to be used for plotting
'''
def get_graphinfo(graph_fn):
    graph = nx.read_gml(graph_fn)

    cells_inhib = np.array(nx.get_node_attributes(graph, 'inh').values(), 
                         dtype=np.int)
    graph_edges = nx.edges(graph)
    number_of_nodes = nx.number_of_nodes(graph)
    degree_histogram = nx.degree_histogram(graph)
    return cells_inhib, graph_edges,number_of_nodes,degree_histogram


'''
    This method gets the time at which the peak occurs for a signal
    goes through the given peak_times and finds at which point the signal is the
    strongest

    peak_times: the times at which a peak in the signal occurs
    signal: The signal that you want to find the max of
'''
def find_max_time(peak_times,signal):
    max_time = np.nan
    for t in peak_times:
        if np.isnan(max_time):
            max_time = t
        elif signal[t] > signal[max_time]:
            max_time = t
    return max_time

'''
    This method finds the phase lag and population correlation for the data given. Use the max peak from the autocorrelation and then the cross correlation peak in the middle of those two 
    
    input:
        xcorr-The cross correlation signal 
        autocorr-An autocorrelations signal to be ran against
        
    output:
        phase-The phase lag/difference between two populations 
        pop_corr-The correlation between both signals
    
'''
def find_metrics(xcorr,autocorr):
    max_time_cross = np.nan;
    peak_auto = scipy.signal.argrelmax(autocorr)[0].tolist()
    peak_cross = scipy.signal.argrelmax(xcorr)[0].tolist()
    max_time = find_max_time(peak_auto,autocorr)
    for i in range(peak_auto.index(max_time)+1,len(peak_auto)):
        if autocorr[peak_auto[i]] > 0:
            max_time_next = peak_auto[i]
            break
    for x in peak_cross:
        if x >= max_time and x <= max_time_next:
            max_time_cross = x
            break
    auto_period = max_time_next - max_time
    auto_cross_perioid = max_time_cross - max_time
    phase = float(auto_cross_perioid)/float(auto_period)
    return phase, xcorr[max_time_cross]


'''
    This method is adapted from the old main methods of the code, this method will do all the post processing 
    and allow for it to be ran indpendently of main to allow for passing in of dictionaries withouth saving and loading them to the hard disc to avoid excess memory usage
    
    Output: 
        mdict - The dictionary of final variables and results. Can either be saved or used as is.
    
'''
def run(sim_output,trans,sec_flag,spike_thresh,f_sigma,butter_low,butter_high,bin_width,cutoff,are_volts):

    butter_freq = np.array([butter_low,butter_high])
    if sec_flag:
        scalet=1e3
    else:
        scalet = 1

    graph_fn = ''
    if isinstance(sim_output['graphFn'],np.ndarray):
        graph_fn = str(sim_output['graphFn'][0])
    else:
        graph_fn = sim_output['graphFn']

    #Retrieve parameters from dictionary and data
    dt = float(sim_output['dt'])*scalet
    data = chop_transient(sim_output['Y'],trans,dt)
    num_neurons = np.shape(data)[0]
    tmax = np.shape(data)[1]

    #Generate spike trains from the data and bin the spikes
    if are_volts:
        spikes, spike_mat = find_spikes(data, spike_thresh)
    else:
        data = scipy.sparse.csc.csc_matrix(data)
        spike_mat= data.todense()
        spikes = data.nonzero()
    bins, spike_mat_bin = bin_spikes(spike_mat, bin_width, dt)

    #Get the different versions of the filtered data
    spike_fil_bin, butter_int_bin, spike_fil_butter = spikes_filt(spike_mat_bin[:num_neurons/2],
                                                              dt*bin_width,
                                                              f_sigma,
                                                              butter_freq)
    spike_fil_bin2, butter_int_bin2, spike_fil_butter2 = spikes_filt(spike_mat_bin[num_neurons/2:],
                                                                     dt*bin_width,
                                                                     f_sigma,
                                                                     butter_freq)

    cross_correlation = xcorr(butter_int_bin2,butter_int_bin)
    auto_cross_correlation1 = xcorr(butter_int_bin,butter_int_bin)
    auto_cross_correlation2 = xcorr(butter_int_bin2,butter_int_bin2)

    phase_lag,pop_corr = find_metrics(cross_correlation,auto_cross_correlation1)

    ##graph attributes
    cells_inhib,graph_edges,number_of_nodes,degree_histogram = get_graphinfo(graph_fn)
    
    mdict = {'bins':bins,
             'spike_mat':spike_mat,
             'spike_mat_bin':spike_mat_bin,
             'spike_fil_bin':spike_fil_bin,
             'spike_fil_bin':spike_fil_bin2,
             'butter_int_bin': butter_int_bin,
             'butter_int_bin2': butter_int_bin2,
             'cross_correlation': cross_correlation,
             'auto_cross_correlation1':auto_cross_correlation1,
             'auto_cross_correlation2':auto_cross_correlation2,
             'cells_inhib': cells_inhib,
             'graph_edges':graph_edges,
             'number_of_nodes':number_of_nodes,
             'degree_histogram':degree_histogram,
             'phase_lag': phase_lag,
             'pop_correlation': pop_corr,
             'time': sim_output['tf'],
             'bin_width': bin_width
            }

    return mdict




def main(argv=None):
    
    should_save = True
    if argv is None:
        argv = sys.argv
    else:
        should_save = False

    (simFn, outFn, trans, sec_flag, spike_thresh, f_sigma, butter_low,
     butter_high, bin_width, cutoff, are_volts) = parse_args(argv)

    sim_output = scipy.io.loadmat(simFn)
    
    post_dict = run(sim_output,trans,sec_flag,spike_thresh,f_sigma,butter_low,butter_high,bin_width,cutoff,are_volts)

    if should_save:
        scipy.io.savemat(outFn,post_dict,oned_as ='column')
    else:
        return post_dict

if __name__ == '__main__':
    status = main()
    sys.exit(status) 
