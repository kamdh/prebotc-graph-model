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
import cmath
import math

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
    peak_order = 30
    peak_percentile = 75
    eta_norm_pts=8
    op_abs_thresh=0.2
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
        args.cut, args.volt,peak_order,peak_percentile,eta_norm_pts,op_abs_thresh,


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
        if x > max_time and x < max_time_next and xcorr[x] > 0:
            max_time_cross = x
            break
    auto_period = max_time_next - max_time
    auto_cross_perioid = max_time_cross - max_time
    phase = float(auto_cross_perioid)/float(auto_period)
    return phase, xcorr[max_time_cross]

'''
    This method finds the population burst peaks for a given signal, uses a percentile filter to elimnate finding noisy peaks 
    
    input:
        signal-This is the signal you want to find the peaks of
        peak_order-The number of points of comparison for each peak on each side of the current value
        peak_percentile-The percentage threshold the peak must meet 
        dt-The time step
    
    output:
        pop_burst_peak-Peaks of the signal that pass the given criteria for a peak
'''
def burst_stats(signal,peak_order,peak_percentile,dt):
    pop_burst_peak=scipy.signal.argrelmax(signal, order=peak_order)[0]
    pop_burst_peak=pop_burst_peak[signal[pop_burst_peak] >
                                    np.percentile(signal,peak_percentile)]
    return pop_burst_peak

'''
    This method is used to get the phi (phase differences) between signals, 
    here we use a moving window to find the refrence perioid to calculate phi
    
    input:
        pop_burst_peak1(2)-This is the time for the peaks from signal 1 or signal 2 
        bins-This is the bin info for the signals after some post processing 
        
    output: 
        phis-List of phis for the signals
'''
def get_phis(pop_burst_peak1,pop_burst_peak2,bins):
    phis = []
    windowStartIndex = 0
    windowEndIndex = 1
    while windowEndIndex < len(pop_burst_peak1):
        windowStart = pop_burst_peak1[windowStartIndex]
        windowEnd = pop_burst_peak1[windowEndIndex]
        peaksInWindow = [i for i in pop_burst_peak2 if i >= windowStart and i <= windowEnd]
        for peak in peaksInWindow:
            phi = (bins[peak] - bins[windowStart]) / (bins[windowEnd] - bins[windowStart])
            phis.append(phi)
        windowStartIndex = windowEndIndex
        windowEndIndex = windowEndIndex + 1
    return phis

'''
    Map phi values to a circle to accuratley take mean and std of the values 
    
    input:
        phis- Phi values that are in [0,1]
    output:
        phis- Phi values that are now mapped to [0,2pi] represents radians
'''
def map_phi_to_complex(phis):
    complex = []
    for i in range(len(phis)):
        radians = 2*np.pi*phis[i]
        complex.append(cmath.rect(1,radians))
    return complex

'''
    This will get the mean phi and variance using circular statistics 
    
    input:
        complex_values- This is a list of complex values that are gotten from the phi values 
        
    output:
        mean_angle- This is the mean angle of the phi values, represents what the average phase is (can be converted back)
        variance_circular- This is the variance of the angles, 0 represents all phi values are the same.
'''
def get_circular_statistics(complex_values):
    mean_resultant = np.mean(complex_values)
    mean_angle = cmath.phase(mean_resultant)
    variance_circular = abs(mean_resultant)
    return mean_angle,variance_circular

'''
    This converts the mean angle back to the standard phi values which lies in [0,1]
    
    input:
        mean_angle- This is the mean angle that was calculated from the list of phis
    
    output: 
        This is the converted average phi values that now consisted with other metrics
'''
def get_normalized_phi(mean_angle):
    if mean_angle < 0:
        return (2*math.pi + mean_angle) / (2*math.pi)
    else:
        return mean_angle / (2*math.pi)

def synchrony_stats(data, dt, maxlags=3000):
    '''
        Synchrony measures
        
        Parameters
        ==========
        data: numneuron x time
        dt: time spacing
        maxlags: maximal lag for autocorrelation, default=3000 ms
        
        Returns
        =======
        chi: synchrony measure
        autocorr: autocorrelation of population avg \bar{data}(t)
        '''
    data_pop=np.mean(data, axis=0) # pop avg
    sigma_pop=np.mean(np.square(data_pop)) - np.square(np.mean(data_pop))
    sigma=np.mean(np.square(data), axis=1) - np.square(np.mean(data, axis=1))
    sigma_mean=np.mean(sigma)
    chisq=sigma_pop / sigma_mean
    chi=np.sqrt(chisq)
    mean_subtract=data_pop - np.mean(data_pop)
    autocorr=scipy.signal.correlate(mean_subtract, mean_subtract,
                                    mode='valid')
    return chi, autocorr

def order_param(eta_norm, eta_t_norm, op_abs_thresh):
    '''
        Compute the order parameter for the normalized (phase) ETAs.
        
        Parameters
        ==========
        eta_norm: normalized ETA array
        eta_t_norm: [-.5, .5] phases corresponding to second axis of array
        op_abs_thresh: float
        
        Returns
        =======
        ops: array of complex valued order parameters, np.nan if undefined
        op_abs: magnitudes
        op_angle: angles
        op_mask: mask of ops with magnitude above threshold
        op_angle_mean: mean angle of significant ops
        op_angle_std: standard deviation of significant ops
        '''
    assert op_abs_thresh < 0.5 and op_abs_thresh >= 0.0,\
        'op_abs_thresh out of range'
    num_neurons=eta_norm.shape[0]
    num_bins=eta_norm.shape[1]
    dtheta=np.min(np.diff(eta_t_norm))
    # below will generate NaNs if the normalization is 0
    density_eta=eta_norm/np.tile(np.sum(eta_norm, axis=1),(num_bins,1)).T
    ops=np.sum(density_eta*
               np.exp(1.0j*
                      np.tile(eta_t_norm,(num_neurons,1))*
                      (2*np.pi)),
               axis=1)
    op_angle=np.angle(ops)/(2*np.pi)
    op_abs=np.abs(ops)
    op_mask=op_abs > op_abs_thresh
    op_angle_mean=np.nanmean(op_angle[op_mask])
    op_angle_std=np.nanstd(op_angle[op_mask])
    return (ops,op_abs,op_angle,op_mask,op_angle_mean,op_angle_std)

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
        Whether to normalize to phase variable
        '''
    breakpts=np.array(
        np.hstack((0, (events[0:-1] + events[1:]) / 2., data.shape[1]-1)),
        dtype=np.int)
    if normalize:
        from scipy.interpolate import griddata
        max_interval=2*pts
        fullrange=np.linspace(-.5, .5, num=max_interval)
        xgrid1=fullrange[0:pts]
        xgrid2=fullrange[pts:]
    else:
        max_interval=2*np.max(np.hstack((events-breakpts[0:-1],
                                        breakpts[1:]-events)))
    midpt=int(np.floor(max_interval / 2))
    numevents=events.shape[0]-2 # don't use 1st and last due to boundary
    eta=np.zeros((data.shape[0], max_interval))
    for j in range(numevents):
        i=j+1
        timeidx=np.arange(int(breakpts[i]), int(breakpts[i+1]), dtype=np.int)
        thisevent=events[i]
        center=int(np.where(timeidx==thisevent)[0].astype(int))
        if normalize:
            xs1=np.array(timeidx[:center] - timeidx[center], dtype=np.float)
            xs1 /= xs1[0]*(-2.0)
            xs2=np.array(timeidx[center+1:] - timeidx[center], dtype=np.float)
            xs2 /= xs2[-1]*2.0
            xs=np.hstack((xs1, xs2))
            toadd=np.apply_along_axis(lambda x:
                                      scipy.interpolate.griddata(
                                          xs, x, fullrange),
                                      1, data[:,timeidx])
            eta += toadd
        else:
            lpad=midpt - center
            rpad=max_interval - (len(timeidx)+lpad)
            eta += np.pad(data[:, timeidx], ((0,0), (lpad,rpad)), 
                          'constant', constant_values=(0,0))
    eta /= float(numevents)
    eta[eta < 0] = 0
    return eta




'''
    This method is adapted from the old main methods of the code, this method will do all the post processing 
    and allow for it to be ran indpendently of main to allow for passing in of dictionaries withouth saving and loading them to the hard disc to avoid excess memory usage
    
    Output: 
        mdict - The dictionary of final variables and results. Can either be saved or used as is.
    
'''
def run(sim_output,trans,sec_flag,spike_thresh,f_sigma,butter_low,butter_high,bin_width,cutoff,are_volts,peak_order,peak_percentile,eta_norm_pts,op_abs_thresh):

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

    #Calculate Correlation Values
    cross_correlation = xcorr(butter_int_bin2,butter_int_bin)
    auto_cross_correlation1 = xcorr(butter_int_bin,butter_int_bin)
    auto_cross_correlation2 = xcorr(butter_int_bin2,butter_int_bin2)

    #phase_lag,pop_corr = find_metrics(cross_correlation,auto_cross_correlation1)

    #graph attributes
    cells_inhib,graph_edges,number_of_nodes,degree_histogram = get_graphinfo(graph_fn)

    #Calculating Values for Circle Map
    pop_burst_peak1 = burst_stats(butter_int_bin,peak_order,peak_percentile,dt*bin_width/1000.)
    pop_burst_peak2 = burst_stats(butter_int_bin2,peak_order,peak_percentile,dt*bin_width/1000.)
    phis = get_phis(pop_burst_peak1,pop_burst_peak2,bins)
    complex_phis = map_phi_to_complex(phis)
    mean_angle,variance_angle = get_circular_statistics(complex_phis)
    mean_phi = get_normalized_phi(mean_angle)
    #std_phi = np.std(phis)

    #Get Synchrony Values for each signal
    chi1,chi1_auto = synchrony_stats(spike_fil_bin,dt*bin_width/1000.)
    chi2,chi2_auto = synchrony_stats(spike_fil_bin2,dt*bin_width/1000.)

    '''##Compute event triggered averages and get individual cell statistics
    ##Population 1
    ##Normalize time to phase variable [-.5,.5]
    eta1_norm = event_trig_avg(pop_burst_peak1,spike_fil_bin,normalize=True,pts=eta_norm_pts)
    eta1_t_norm = np.linspace(-0.5, 0.5, 2*eta_norm_pts)
    ##Order Parameters
    (ops1,op_abs1,op_angle1,op_mask1,
     op_angle_mean1,op_angle_std1)=order_param(eta1_norm,eta1_t_norm,op_abs_thresh)
    ##Population 2
    ##Normalize time to phase variable [-.5,.5]
    eta2_norm = event_trig_avg(pop_burst_peak2,spike_fil_bin2,normalize=True,pts=eta_norm_pts)
    eta2_t_norm = np.linspace(-0.5, 0.5, 2*eta_norm_pts)
    ##Order Parameters
    (ops2,op_abs2,op_angle2,op_mask2,
     op_angle_mean2,op_angle_std2)=order_param(eta2_norm,eta2_t_norm,op_abs_thresh)'''
    


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
             #'phase_lag': phase_lag,
             #'pop_correlation': pop_corr,
             'time': sim_output['tf'],
             'bin_width': bin_width,
             'phis' : phis,
             'mean_phi': mean_phi,
             'variance_angle' : variance_angle,
             'chi1' : chi1,
             'chi2' : chi2,
	     'pop_burst_peak1': pop_burst_peak1,
	     'pop_burst_peak2': pop_burst_peak2
             #'op_abs1' : op_abs1,
             #'op_angle1' : op_angle1,
             #'op_angle_mean1' : op_angle_mean1,
             #'op_angle_std1' : op_angle_std1,
             #'op_abs2' : op_abs2,
             #'op_angle2' : op_angle2,
             #'op_angle_mean2' : op_angle_mean2,
             #'op_angle_std2' : op_angle_std2
            }

    return mdict




def main(argv=None):
    
    should_save = True
    if argv is None:
        argv = sys.argv
    else:
        
        should_save = False

    (simFn, outFn, trans, sec_flag, spike_thresh, f_sigma, butter_low,
     butter_high, bin_width, cutoff, are_volts,peak_order,peak_percentile,eta_norm_pts,op_abs_thresh) = parse_args(argv)

    sim_output = scipy.io.loadmat(simFn)
    
    post_dict = run(sim_output,trans,sec_flag,spike_thresh,f_sigma,butter_low,butter_high,bin_width,cutoff,are_volts,peak_order,peak_percentile,eta_norm_pts,op_abs_thresh)

    if should_save:
        scipy.io.savemat(outFn,post_dict,oned_as ='column')
    else:
        return post_dict

if __name__ == '__main__':
    status = main()
    sys.exit(status) 
