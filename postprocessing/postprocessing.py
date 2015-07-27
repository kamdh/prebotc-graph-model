'''
postprocessing.py

A library of functions used for postprocessing the preBotC model output.
'''

import numpy as np
import scipy.signal
from IPython import embed

def chop_transient(data, transient, dt):
    '''
    Remove a transient from the data
    '''
    if transient > 0:
        firstIdx=int(np.ceil(transient / dt) - 1)
        return data[:,firstIdx:]
    else:
        return data

def find_spikes(data, threshold):
    '''
    Find spikes in voltage data by taking relative maxima
    '''
    indices=scipy.signal.argrelmax(data, axis=1) # coords 1-2 of maxima
    mask=np.where(data[indices] > threshold)
    new_indices=(indices[0][mask],
                 indices[1][mask])
    spike_mat=np.zeros(np.shape(data), dtype=np.int) # dense format
    spike_mat[new_indices]=1
    return new_indices, spike_mat

def spikes_of_neuron(spikes, neuron):
    '''
    Return time indices of spiking of a given neuron
    '''
    return spikes[1][np.where(spikes[0] == neuron)]

def filter_spikes(spike_mat, samp_dt, f_sigma, butter_freq):
    '''
    Filter the spike timeseries. Returns both neuron-by-neuron timeseries
    filtered with a gaussian kernel and the population data filtered
    with a butterworth filter.

    Parameters
    ==========
    spike_mat: the numneuron x time matrix of spikes
    samp_dt: period (in s) between measurements in spike_mat
    f_sigma: variance of gaussian (ms)
    butter_freq: butterworth filter cutoff frequency(s)

    Returns
    =======
    spike_fil: gaussian filtered matrix, same shape as spike_mat (spikes/s)
    int_signal: butterworth filtered population timeseries (spikes/neuron/s)
    '''
    def filt_window_gauss(samp_dt, std=20, width=None, normalize=1):
        if width is None:
            width=std*4+1
        width /= (1000.0*samp_dt)
        w=scipy.signal.gaussian(width, std)
        if not normalize == 0:
            w=normalize * w / sum(w)
        return w
    
    def filt_gauss(spike_mat, samp_dt, f_sigma=20):
        w=filt_window_gauss(samp_dt, std=f_sigma, normalize=1)
        spike_fil=scipy.signal.fftconvolve(spike_mat, w[ np.newaxis, : ], 
                                           mode='same')
        #spike_fil=scipy.signal.convolve(spike_mat, w[ np.newaxis, : ], 
                    #                                  mode='same')
        spike_fil/=samp_dt
        spike_fil[np.abs(spike_fil) < 1e-10]=0.0
        return spike_fil
    
    def filt_butter(data, samp_dt, butter_freq, axis=-1):
        '''
        Filter data with a 2nd order butterworth filter.
        
        Parameters
        ==========
          data: ndarray
          samp_dt: sampling period (s)
          butter_freq: [cutoff_low, cutoff_high] (Hz), can be infinite
          axis (optional): axis along which to filter, default=-1
        Returns
        =======
          data_butter: filtered version of data
        '''
        order=2
        ny=0.5 / samp_dt # Nyquist frequency
        cof=butter_freq / ny # normalized cutoff freq
        if np.isneginf(cof[0]) and np.isfinite(cof[1]):
            # lowpass
            cof1=cof[1]
            b, a=scipy.signal.butter(order, cof1, btype='low')
            data_butter=scipy.signal.filtfilt(b, a, data, axis=axis)
        elif np.isfinite(cof[0]) and np.isinf(cof[1]):
            # highpass
            cof1=cof[0]
            b, a=scipy.signal.butter(order, cof1, btype='high')
            data_butter=scipy.signal.filtfilt(b, a, data, axis=axis)
        elif np.isfinite(cof[0]) and np.isfinite(cof[1]):
            # bandpass
            b, a=scipy.signal.butter(order, cof, btype='band')
            data_butter=scipy.signal.filtfilt(b, a, data, axis=axis)
        else:
            raise Exception('filt_butter called with bad cutoff frequency')
        data_butter /= samp_dt # normalize to rate
        return data_butter
    
    spike_fil=filt_gauss(spike_mat, samp_dt, f_sigma=f_sigma) 
    int_signal=filt_butter(np.mean(spike_mat, axis=0), 
                           samp_dt, butter_freq)
    ## removed below because it is a large matrix for high samp_dt
    # spike_fil_butter=filt_butter(spike_fil, samp_dt,
    #                                butter_freq, axis=1)
    int_signal[int_signal < 0]=0.0
    return spike_fil, int_signal

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
    num_times=np.shape(spike_mat)[1]
    stride=int(np.ceil(bin_width / dt))
    bins=np.arange(0, num_times, stride, dtype=np.float)
    which_bins=np.digitize(range(0, num_times), bins)
    num_bins=len(bins)
    binned_spikes=np.zeros((num_neurons, num_bins), dtype=np.float)
    for i in range(num_bins):
        bin_mask=np.where(which_bins == i)[0] # mask data in bin i, tuple
        bin_data=spike_mat[:,bin_mask]
        binned_spikes[:,i]=np.sum(bin_data, axis=1).flatten()
    return bins, binned_spikes

def bin_subsamp(data, bins):
    '''
    Subsample timeseries into the bins. Accepts the same paramaters as
    bin_spikes but does not return the bins themselves.

    Parameters
    ==========
      data: matrix of timeseries, (num_neuron x num_time)
      bins: indices of the bin endpoints

    Returns
    =======
      binned_data: a new matrix (num_neuron x num_bins)
    '''
    binned_data=data[:,np.array(bins, dtype=np.int)]
    # num_neurons= np.shape(data)[0]
    # num_times=np.shape(data)[1]
    # stride=int(np.ceil(bin_width / dt))
    # #bins=np.arange(0, num_times, stride, dtype=np.float)
    # bins=range(0, num_times, stride)
    # binned_data=data[:,bins] # subsample points
    # # which_bins=np.digitize(range(0, num_times), bins)
    # # num_bins=len(bins)
    # # binned_data=np.zeros((num_neurons, num_bins), dtype=np.int)
    # # # costly but better than subsampling points
    # # for i in range(num_bins):
    # #     bin_mask=np.where(which_bins == i)[0] # mask data in bin i, tuple
    # #     bin_data=data[:,bin_mask]
    # #     binned_data[:,i]=np.mean(bin_data, axis=1).flatten()
    return binned_data

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
    # autocorr=acorr(data_pop - np.mean(data_pop), onesided=True, 
    #                  scale='coeff')
    mean_subtract=data_pop - np.mean(data_pop)
    autocorr=scipy.signal.correlate(mean_subtract, mean_subtract, 
                                    mode='valid')
    return chi, autocorr

def peak_freq_welch(data, dt):
    '''
    Compute the Welch periodogram (psd) and return the peak frequency

    Parameters
    ==========
      data: array of 
      dt: time interval (s) between data points

    Returns
    =======
      peak_lag: dominant period in psd
      peak_freq: dominant frequency in psd
      freq: frequencies in spectrogram
      power: power in spectrogram
    '''
    # f, Pxx=scipy.signal.welch(data, fs=1000/dt, detrend='constant',
    #                              return_onesided=True, nperseg=2**14)
    data=np.array(data, dtype=np.float)
    freq, power=scipy.signal.periodogram(data, fs=1.0/dt,
                                         return_onesided=True, 
                                         detrend='constant')
    idx=np.argmax(power)
    peak_freq=freq[idx]
    peak_lag=1/peak_freq
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
    whenSpiking=np.nonzero(raster)[0]
    isi_vec=np.diff(whenSpiking)
    isi_vec=isi_vec[isi_vec != 1]
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
    new_raster=np.array(np.logical_not(raster), dtype=np.int)
    w=np.hstack((1, new_raster, 1)) # pad with 1s
    runs_zeros=np.nonzero(np.diff(w) == 1)[0] - np.nonzero(np.diff(w) == -1)[0]
    return runs_zeros

def burst_starts(raster):
    '''
    Find the starting points of each burst, where raster entries jump
    from 0 to 1.
    '''
    return np.where(np.diff(raster) == 1)[0] + 1

def burst_stats_old(data, cutoff, dt):
    '''
    Estimate when the population is bursting by comparing filtered
    activity data with a threshold=cutoff*(max(data) - min(data)) + min(data).

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
    mean_data=np.mean(data[20:-20])
    std_data=np.std(data[20:-20])
    thresh=mean_data + std_data*cutoff
    bursting=np.array(data > thresh, dtype=np.float)
    duty_cycle=np.sum(bursting) / bursting.shape[0]
    ibi_vec=isi(bursting) * dt
    burst_start_locs=burst_starts(bursting)
    burst_lengths=burst_lens(bursting)
    burst_peak_locs=np.zeros(burst_start_locs.shape)
    burst_peaks=np.zeros(burst_start_locs.shape)
    bad_bursts=0
    for i in range(len(burst_start_locs)):
        # find peak of burst i
        burst_index=burst_start_locs[i] + range(burst_lengths[i])
        tmp=data[ burst_index ]
        # peakInd=np.argmax(tmp)
        peak_index=scipy.signal.argrelmax(tmp)[0]
        if len(peak_index) > 1 or not peak_index:
            ## more than one peak found or empty list
            #print "more than one local max found in burst " + str(i)
            bad_bursts += 1
            peak_index=np.argmax(tmp)
        # add 1 for Matlab 1-based indexing
        burst_peak_locs[i]=burst_index[ np.int(peak_index) ] + 1
        burst_peaks[i]=data[ burst_index[ np.int(peak_index) ] ]
    ibi_mean=np.mean(ibi_vec)
    ibi_cv=np.std(ibi_vec) / ibi_mean
    burst_length_mean=np.mean(burst_lengths * dt)
    burst_length_cv=np.std(burst_lengths * dt) / burst_length_mean
    burst_start_locs += 1 # for Matlab
    return (duty_cycle, ibi_mean, ibi_cv, burst_length_mean, burst_length_cv, 
            ibi_vec, burst_lengths, burst_start_locs, burst_peak_locs, 
            burst_peaks, bursting, bad_bursts)

def graph_attributes(graph_fn):
    '''
    Load vertex attributes so we can access them from MATLAB.
    Note that if the attribute does not exist, this will be an empty array.
    '''

    import networkx as nx
    g=nx.read_gml(graph_fn)
    vertex_types=np.array(nx.get_node_attributes(g, 'type').values(),
                          dtype=np.int)
    vertex_inh=np.array(nx.get_node_attributes(g, 'inh').values(), 
                        dtype=np.int)
    vertex_respir_area=np.array(
        nx.get_node_attributes(g,'respir_area').values(),
        dtype=np.int)
    graph_adj=nx.adjacency_matrix(g, weight='gsyn')
    bin_adj=np.matrix(graph_adj > 0, dtype=np.float)
    adj_inh=np.multiply(np.tile(vertex_inh,(g.number_of_nodes(),1)).T,
                        bin_adj)
    adj_exc=np.multiply(np.tile(1-vertex_inh,(g.number_of_nodes(),1)).T,
                        bin_adj)
    return (vertex_types, vertex_inh, vertex_respir_area, graph_adj,
            bin_adj, adj_exc, adj_inh)

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

def eta_vertex_inputs(eta,bin_adj,vertex_inh):
    '''
    Compute the event triggered average summed over the neighbors
    of each cell. This characterizes their input.

    Parameters
    ==========
    eta, event triggered average
    bin_adj, adjacency matrix
    vertex_inh, logical array indicating inhibitory type

    Returns
    =======
    exc_input, excitatory input to each cell
    inh_input, inhibitory input to each cell    
    '''
    m=eta.shape[0]
    n=eta.shape[1]
    inh_input=np.zeros((m,n))
    exc_input=np.zeros((m,n))
    for v in range(m):
        neighbors=np.array(np.where(bin_adj[:,v])[0]).flatten()
        for u in neighbors:
            if vertex_inh[u]:
                inh_input[v,:]+=eta[u,:]
            else:
                exc_input[v,:]+=eta[u,:]
    return exc_input, inh_input

def nmf_error(eta):
    '''
    Decompose the ETAs using nonnegative matrix factorization with 1 component.
    The reconstruction error is an estimate of the non-inspiratory activity.
    '''
    from sklearn.decomposition import NMF
    nmf=NMF(n_components=1)
    nmf.fit(eta)
    reconstruction_err=nmf.reconstruction_err_
    return float(reconstruction_err)

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

def predict_ops_olm(ops, predictors):
    import statsmodels.api as sm
    model=sm.OLS(ops, predictors)
    results=model.fit()
    return results

def irregularity_score(ts):
    return np.mean(np.abs(np.diff(ts))/ts[0:-1])
    
def burst_stats(signal,peak_order,peak_percentile,dt):
    pop_burst_peak=scipy.signal.argrelmax(signal, order=peak_order)[0]
    pop_burst_peak=pop_burst_peak[signal[pop_burst_peak] >
                                  np.percentile(signal,peak_percentile)]
    pop_burst_trough=scipy.signal.argrelmin(signal, order=peak_order)[0]
    pop_burst_trough=pop_burst_trough[signal[pop_burst_trough] <
                                      np.percentile(signal,100.-peak_percentile)]
    ibi_vec=np.diff(pop_burst_peak)*dt/1000.0
    ibi_mean=np.mean(ibi_vec)
    ibi_cv=np.std(ibi_vec)/ibi_mean
    ibi_irregularity=irregularity_score(ibi_vec)
    amplitude_irregularity=irregularity_score(signal[pop_burst_peak])
    amplitude_cv=np.std(signal[pop_burst_peak])/np.mean(signal[pop_burst_peak])
    peak_to_trough=(signal[pop_burst_peak].mean() -
                    signal[pop_burst_trough].mean())/signal.mean()
    return (pop_burst_peak,pop_burst_trough,ibi_vec,ibi_mean,ibi_cv,
        ibi_irregularity, amplitude_irregularity, amplitude_cv, peak_to_trough)
