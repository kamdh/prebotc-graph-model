#!/usr/bin/env python
'''
doPost.py

python version of postprocessing script for identifying spikes and filtering the
output.
'''

import warnings
import sys
import numpy as np
import scipy.signal
import scipy.io
import argparse
#from IPython import embed
from postprocessing import *
from classify_phase import fit_MRF_pseudolikelihood

def parse_args(argv):
    # defaults
    transient=20000 # ms
    spike_thresh=-20.0 # mV
    f_sigma=60.0 # ms
    butter_high=4.0 # Hz
    butter_low=-np.inf # Hz
    bin_width=50 # ms
    cutoff=0.5
    peak_spacing=1.0 # s
    eta_norm_pts=8
    op_abs_thresh=0.2
    silent_firing_rate=0.1 # Hz
    peak_percentile=75
    # parsing
    parser=argparse.ArgumentParser(prog="doPost",
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
                        help='file contains voltage traces ' + \
                        '(default: sparse spike trains)')
    parser.add_argument('--thresh', 
                        help='spike threshold, mV (default: %(default)s)',
                        type=float, default=spike_thresh) 
    parser.add_argument('--fsig', '-f', 
                        help='filter standard deviation, ms ' + \
                        '(default: %(default)s)',
                        type=float, default=f_sigma)
    parser.add_argument('--butter_high', 
                        help='Butterworth filter upper cutoff frequency, Hz ' +\
                        '(default: %(default)s)', 
                        type=float, default=butter_high)
    parser.add_argument('--butter_low', 
                        help='Butterworth filter lower cutoff frequency, Hz '+\
                        '(default: %(default)s)',
                        type=float, default=butter_low)
    parser.add_argument('--bin_width', '-b', 
                        help='bin width, ms (default: %(default)s)',
                        type=float, default=bin_width)
    parser.add_argument('--cut', '-c', 
                        help='burst cutoff parameter (default: %(default)s)',
                        type=float, default=cutoff)
    parser.add_argument('--peak_spacing', 
                        help='maximum order for defining peak number of bins'+\
                        ', seconds (default: %(default)s)',
                        type=float, default=peak_spacing)
    parser.add_argument('--eta_norm_pts',
                        help='half the number of points in [-.5, .5] onto '+\
                        'which to interpolate (default: %(default)s)',
                        type=int, default=eta_norm_pts)
    parser.add_argument('--op_abs_thresh', 
                        help='threshold order parameter magnitude above '+\
                        'which to calculate statistics ' +\
                        '(default: %(default)s)',
                        type=float, default=op_abs_thresh)
    parser.add_argument('--silent_firing_rate',
                        help='threshold firing rate below which to classify '+\
                        'neuron as silent, Hz (default: %(default)s)',
                        type=float, default=silent_firing_rate)
    parser.add_argument('--peak_percentile',
                        help='percentile threshold for peak detection ' +\
                        '(default: %(default)s)',
                        type=float, default=peak_percentile)
    args=parser.parse_args(argv[1:])
    return (args.sim, args.output, args.transient, args.sec, args.thresh,
            args.fsig, args.butter_low, args.butter_high, args.bin_width,
            args.cut, args.volt, args.peak_spacing, args.peak_percentile,
            args.eta_norm_pts,
            args.op_abs_thresh, args.silent_firing_rate)
    
def main(argv=None):
    '''
    main function
    '''
    ## Setup parameters for postprocessing
    if argv is None:
        argv=sys.argv
        (simFn, outFn, trans, sec_flag,
         spike_thresh, f_sigma, butter_low,
         butter_high, bin_width, cutoff, are_volts, 
         peak_spacing, peak_percentile, eta_norm_pts, op_abs_thresh,
         silent_firing_rate)=parse_args(argv)

    butter_freq=np.array([butter_low, butter_high])
    peak_order = int(np.round(peak_spacing / bin_width * 1e3))

    if sec_flag:
        scalet=1e3
    else:
        scalet=1
        
    ## Load simulation output
    sim_output=scipy.io.loadmat(simFn)

    ## Setup a few more variables, assumes no --save_full
    if sim_output['saveStr'][0] == 'full':
        raise Exception('output should not be from --save_full option')

    graph_fn=str(sim_output['graphFn'][0])
    dt=float(sim_output['dt']) * scalet
    data=chop_transient(sim_output['Y'], trans, dt)
    num_neurons=np.shape(data)[0]

    ## Load in graph data
    (vertex_types, vertex_inh, vertex_respir_area, graph_adj, bin_adj,
     adj_exc, adj_inh)=graph_attributes(graph_fn)
    
    ## Begin postprocessing
    ## spike trains
    if are_volts:
        spikes, spike_mat=find_spikes(data,spike_thresh)
    else:
        spike_mat=data.todense()
        spikes=data.nonzero()

    ## Bin spikes
    bins, spike_mat_bin=bin_spikes(spike_mat, bin_width, dt)
    
    ## Filter spike raster for integrated activity, filtered spike trains
    spike_fil, butter_int=filter_spikes(spike_mat, dt/1000.0, 
                                        f_sigma, butter_freq)
    elec_neur=np.bitwise_or(vertex_types==2, # TS
                            vertex_types==3) # Q
    spike_fil_elec, butter_int_elec= \
      filter_spikes(spike_mat[elec_neur,], dt/1000.0, f_sigma, butter_freq)
    spike_fil_bin=bin_subsamp(spike_fil, bins)
    spike_fil_bin_elec=bin_subsamp(spike_fil_elec, bins)
    butter_int_bin=bin_subsamp(butter_int, bins).flatten()
    butter_int_bin_elec=bin_subsamp(butter_int_elec, bins).flatten()
    # spike_fil_bin,butter_int_bin=filter_spikes(spike_mat_bin, 
    #                                            dt*bin_width/1000.0, 
    #                                            f_sigma, 
    #                                            butter_freq)

    ## Peri-neuron time histogram
    psth_bin=np.sum(spike_mat_bin,axis=0)
    ## Synchrony measures from autocorrelation
    if are_volts:
        chi, autocorr=synchrony_stats(data, dt)
    else:
        chi, autocorr=synchrony_stats(spike_fil_bin, dt*bin_width/1000.)
    chi_elec, _=synchrony_stats(spike_fil_bin_elec, dt*bin_width/1000.)
    
    ## peak freqs
    peak_lag, peak_freq, freq, power=peak_freq_welch(psth_bin, 
                                                     dt*bin_width/1000.0)

    ## Old burst stats, remove??
    # (duty_cycle, ibi_mean, ibi_cv, burst_length_mean, burst_length_cv, ibi_vec,
    #  burst_lengths, burst_start_locs, burst_peak_locs, burst_peaks, bursting, 
    #  bad_bursts)=burst_stats(butter_int_bin, cutoff, dt*bin_width)
    ## New burst stats
    firing_rates=np.sum(spike_mat_bin,axis=1)/((bins[-1]-bins[0])/1000.0)
    avg_firing_rate=np.mean(firing_rates)

    ## Compute the population burst peaks
    pop_burst_peak,pop_burst_trough,ibi_vec,ibi_mean,ibi_cv, \
      ibi_irregularity, amplitude_mean, \
      amplitude_irregularity, amplitude_cv, peak_to_trough = \
      burst_stats(butter_int_bin,peak_order,peak_percentile,dt*bin_width/1000.)
    print "ptt: " + str(peak_to_trough)
    pop_burst_peak_elec,pop_burst_trough_elec,ibi_vec_elec,ibi_mean_elec,\
      ibi_cv_elec, ibi_irregularity_elec, \
      amplitude_mean_elec, amplitude_irregularity_elec, \
      amplitude_cv_elec, peak_to_trough_elec \
       = burst_stats(butter_int_bin_elec,peak_order,peak_percentile,
                     dt*bin_width/1000.)
    print "ptt elec: " + str(peak_to_trough_elec)
        
    ## Compute event triggered averages
    ## First, using absolute time
    eta=event_trig_avg(pop_burst_peak, spike_fil_bin)
    eta_t=(np.arange(eta.shape[1])-eta.shape[1]/2)*bin_width
    ## Second, normalize time to be a phase variable [-.5, .5]
    eta_norm=event_trig_avg(pop_burst_peak, spike_fil_bin, normalize=True,
                            pts=eta_norm_pts)
    eta_t_norm=np.linspace(-0.5, 0.5, 2*eta_norm_pts)

    ## Order parameters   
    (ops,op_abs,op_angle,op_mask,
     op_angle_mean,op_angle_std)=order_param(eta_norm,eta_t_norm,op_abs_thresh)

    ## Classify cells as inspiratory/expiratory/tonic/silent
    mask_expir=np.bitwise_and(op_mask,
        np.bitwise_or(op_angle<-0.25,op_angle>0.25))
    mask_inspir=np.bitwise_and(op_mask,
        np.bitwise_and(op_angle<=0.25,op_angle>=-0.25))
    mask_silent=firing_rates<silent_firing_rate
    mask_tonic=np.bitwise_not(mask_expir | mask_inspir | mask_silent)
    num_silent=np.sum(mask_silent)
    num_expir=np.sum(mask_expir)
    num_inspir=np.sum(mask_inspir)
    num_tonic=np.sum(mask_tonic)
    (exc_input,inh_input)=eta_vertex_inputs(eta_norm,bin_adj,vertex_inh)

    ## Average input those cells receive
    avg_inspir=np.mean(np.hstack((inh_input[mask_inspir,:],
                                  exc_input[mask_inspir,:])),axis=0)
    avg_expir=np.mean(np.hstack((inh_input[mask_expir,:],
                                 exc_input[mask_expir,:])),axis=0)
    avg_tonic=np.mean(np.hstack((inh_input[mask_tonic,:],
                                 exc_input[mask_tonic,:])),axis=0)
    avg_silent=np.mean(np.hstack((inh_input[mask_silent,:],
                                  exc_input[mask_silent,:])),axis=0)
    mrf_coeff=fit_MRF_pseudolikelihood(adj_exc[op_mask][:,op_mask],
                                       adj_inh[op_mask][:,op_mask],
                                       np.array(mask_expir[op_mask],
                                                dtype='float'))
    
    ## Save output
    scipy.io.savemat(outFn,
                     mdict={'bins': bins,
                            'spike_mat_bin': spike_mat_bin,
                            'spike_fil_bin': spike_fil_bin,
                            'butter_int_bin': butter_int_bin,
                            'butter_int_bin_elec': butter_int_bin_elec,
                            'psth_bin': psth_bin,
                            'chi': chi,
                            'chi_elec': chi_elec,
                            'autocorr': autocorr,
                            'peak_lag': peak_lag,
                            'peak_freq': peak_freq,
                            'pop_burst_peak': pop_burst_peak,
                            'pop_burst_trough': pop_burst_trough,
                            'ibi_vec': ibi_vec,
                            'ibi_mean': ibi_mean,
                            'ibi_cv': ibi_cv,
                            'ibi_irregularity' : ibi_irregularity,
                            'amplitude_mean' : amplitude_mean,
                            'amplitude_irregularity' : amplitude_irregularity,
                            'amplitude_cv' : amplitude_cv,
                            'peak_to_trough': peak_to_trough,
                            'pop_burst_peak_elec': pop_burst_peak_elec,
                            'pop_burst_trough_elec': pop_burst_trough_elec,
                            'ibi_vec_elec': ibi_vec_elec,
                            'ibi_mean_elec': ibi_mean_elec,
                            'ibi_cv_elec': ibi_cv_elec,
                            'ibi_irregularity_elec' : ibi_irregularity_elec,
                            'amplitude_mean_elec' : amplitude_mean_elec,
                            'amplitude_irregularity_elec' : \
                            amplitude_irregularity_elec,
                            'amplitude_cv_elec' : amplitude_cv_elec,
                            'peak_to_trough_elec':peak_to_trough_elec,
                            'eta': eta,
                            'eta_t': eta_t,
                            'eta_norm': eta_norm,
                            'eta_t_norm': eta_t_norm,
                            'do_post_argv': " ".join(argv),
                            'vertex_types': vertex_types,
                            'vertex_inh': vertex_inh,
                            'vertex_respir_area': vertex_respir_area,
                            'graph_adj': graph_adj,
                            'ops': ops,
                            'op_angle_mean': op_angle_mean,
                            'op_angle_std': op_angle_std,
                            'avg_firing_rate': avg_firing_rate,
                            'num_silent': num_silent,
                            'num_expir': num_expir,
                            'num_inspir': num_inspir,
                            'num_tonic': num_tonic,
                            'eta_avg_inspir': avg_inspir,
                            'eta_avg_expir': avg_expir,
                            'eta_avg_tonic': avg_tonic,
                            'eta_avg_silent': avg_silent,
                            'mask_inspir': np.array(mask_inspir,dtype=np.int),
                            'mask_expir': np.array(mask_expir,dtype=np.int),
                            'mask_silent': np.array(mask_silent,dtype=np.int),
                            'mask_tonic': np.array(mask_tonic,dtype=np.int),
                            'num_silent': num_silent,
                            'num_expir': num_expir,
                            'num_inspir': num_inspir,
                            'num_tonic': num_tonic,
                            'adj_exc': adj_exc,
                            'adj_inh': adj_inh,
                            'bin_adj': bin_adj,
                            'mrf_coeff': mrf_coeff
                        },
                     oned_as='column', do_compression=True)

# run the main stuff
if __name__ == '__main__':
    main()
