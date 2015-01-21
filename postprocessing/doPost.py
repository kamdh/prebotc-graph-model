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
from IPython import embed
from postprocessing import *

def parse_args(argv):
    # defaults
    transient=20000 # ms
    spike_thresh=-20.0 # mV
    f_sigma=20.0 # ms
    butter_high=4.0 # Hz
    butter_low=-np.inf # Hz
    bin_width=20 # ms
    cutoff=0.5
    peak_order=20
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
    parser.add_argument('--peak_order', 
                        help='maximum order for defining peak number of bins'+\
                        ', bin count (default: %(default)s)',
                        type=int, default=peak_order)
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
            args.cut, args.volt, args.peak_order, args.peak_percentile,
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
         peak_order, peak_percentile, eta_norm_pts, op_abs_thresh,
         silent_firing_rate)=parse_args(argv)

    butter_freq=np.array([butter_low, butter_high])

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
    spike_fil_bin=bin_subsamp(spike_fil, bins)
    butter_int_bin=bin_subsamp(butter_int, bins).flatten()
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
        chi, autocorr=synchrony_stats(spike_fil_bin, dt*bin_width)
    
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
    pop_burst_peak=scipy.signal.argrelmax(butter_int_bin, order=peak_order)[0]
    pop_burst_peak=pop_burst_peak[butter_int_bin[pop_burst_peak] >
                                  np.percentile(butter_int_bin,peak_percentile)]
    ibi_vec=np.diff(pop_burst_peak)*bin_width*dt/1000.0
    ibi_mean=np.mean(ibi_vec)
    ibi_cv=np.std(ibi_vec)/ibi_mean

    ## Compute event triggered averages
    ## First, using absolute time
    eta=event_trig_avg(pop_burst_peak, spike_fil_bin)
    eta_t=(np.arange(eta.shape[1])-eta.shape[1]/2)*bin_width
    ## Second, normalize time to be a phase variable [-.5, .5]
    eta_norm=event_trig_avg(pop_burst_peak, spike_fil_bin, normalize=True,
                            pts=eta_norm_pts)
    eta_t_norm=np.linspace(-0.5, 0.5, 2*eta_norm_pts)

    ## Nonnegative matrix factorizations as expiratory measure (old)
    # eta_nmf_err=nmf_error(eta) 
    # eta_norm_nmf_err=nmf_error(eta_norm)

    ## Load in graph data for matlab to use later, so far unused
    (vertex_types, vertex_inh, vertex_respir_area, graph_adj, bin_adj,
     adj_exc, adj_inh)=graph_attributes(graph_fn)
    ## Order parameters   
    (ops,op_abs,op_angle,op_mask,
     op_angle_mean,op_angle_std)=order_param(eta_norm,eta_t_norm,op_abs_thresh)

    ## Classify cells as inspiratory/expiratory/tonic/silent
    inh_in_deg=np.sum(np.multiply(bin_adj,
                                  np.tile(vertex_inh,#*firing_rates,
                                          (300,1)).T), 0)
    exc_in_deg=np.sum(np.multiply(bin_adj,
                                  np.tile((1-vertex_inh),#*firing_rates,
                                          (300,1)).T), 0)
    inh_in_deg=np.array(inh_in_deg).flatten()
    exc_in_deg=np.array(exc_in_deg).flatten()
    expir_mask=np.bitwise_and(op_mask,
        np.bitwise_or(op_angle<-0.25,op_angle>0.25))
    inspir_mask=np.bitwise_and(op_mask,
        np.bitwise_and(op_angle<=0.25,op_angle>=-0.25))
    silent_mask=firing_rates<silent_firing_rate
    tonic_mask=np.bitwise_not(expir_mask | inspir_mask | silent_mask)
    num_silent=np.sum(silent_mask)
    num_expir=np.sum(expir_mask)
    num_inspir=np.sum(inspir_mask)
    num_tonic=np.sum(tonic_mask)
    (exc_input,inh_input)=eta_vertex_inputs(eta_norm,bin_adj,vertex_inh)

    ## Average input those cells receive
    avg_inspir=np.mean(np.hstack((inh_input[inspir_mask,:],
                                  exc_input[inspir_mask,:])),axis=0)
    avg_expir=np.mean(np.hstack((inh_input[expir_mask,:],
                                 exc_input[expir_mask,:])),axis=0)
    avg_tonic=np.mean(np.hstack((inh_input[tonic_mask,:],
                                 exc_input[tonic_mask,:])),axis=0)
    avg_silent=np.mean(np.hstack((inh_input[silent_mask,:],
                                  exc_input[silent_mask,:])),axis=0)
    inh_inspir_in_deg=np.array(np.sum(np.multiply(bin_adj,
        np.tile((vertex_inh & inspir_mask), #*firing_rates,
                (300,1)).T),
        0)).flatten()
    inh_expir_in_deg=np.array(np.sum(np.multiply(bin_adj,
        np.tile((vertex_inh & expir_mask), #*firing_rates,
                (300,1)).T),
        0)).flatten()
    exc_inspir_in_deg=np.array(np.sum(np.multiply(bin_adj,
        np.tile(((1-vertex_inh)&inspir_mask), #*firing_rates,
                (300,1)).T),
        0)).flatten()
    exc_expir_in_deg=np.array(np.sum(np.multiply(bin_adj,
        np.tile(((1-vertex_inh)&expir_mask), #*firing_rates,
                (300,1)).T),
        0)).flatten()
    # fit_ols=predict_ops_olm(np.abs(op_angle[op_mask]),
    #                         np.column_stack((vertex_types[op_mask],
    #                                          inh_inspir_in_deg[op_mask],
    #                                          inh_expir_in_deg[op_mask],
    #                                          exc_inspir_in_deg[op_mask],
    #                                          exc_expir_in_deg[op_mask])))
    Ys=np.array(expir_mask[op_mask],dtype='float')
    Xs=np.column_stack((#vertex_types[op_mask],
                       # np.ones((op_mask.sum(),)), # only for statsmodels
                        inh_inspir_in_deg[op_mask],
                        inh_expir_in_deg[op_mask],
                        exc_inspir_in_deg[op_mask],
                        exc_expir_in_deg[op_mask]))
    from sklearn.linear_model import LogisticRegression
    from sklearn import cross_validation
    logmodel=LogisticRegression(penalty='l1')
    scores=cross_validation.cross_val_score(logmodel,Xs,Ys,
                                            cv=10,scoring='f1')
    # print 'F1 scores for logit'
    # print scores
    print("F1 score summary: %0.2f (+/- %0.2f)" % \
          (scores.mean(), scores.std() * 2))
    logmodel.fit(Xs,Ys)
    print("Coefficients for full fit: %s") % str(logmodel.raw_coef_)
    print("Sign of coefficients: %s") % str(np.sign(logmodel.coef_))
    #fit_log=predict_ops_logit(Ys,Xs)
    # print fit_log.summary()
    # print 'params'
    # print fit_log.params
    # print 'odds ratios'
    # print np.exp(fit_log.params)    


    # ## some plotting
    # import matplotlib.pyplot as plt
    # plt.ion()
    # plt.plot(bins*dt/1000.0,butter_int_bin)
    # plt.hold(True)
    # plt.plot(bins[pop_burst_peak]*dt/1000.0,
    #          butter_int_bin[pop_burst_peak],'ko')
    
    # plt.figure()
    # plt.hold(True)
    # plt.plot(avg_inspir)
    # plt.plot(avg_expir)
    # plt.plot(avg_tonic)
    # plt.legend(('inspir','expir','tonic'),'upper right')

    # plt.figure()
    # plt.imshow(spike_mat_bin,cmap='gray')


    # plt.figure()
    # xs=np.linspace(butter_int_bin.min(), butter_int_bin.max(),500)
    # import scipy.stats as ss
    # fit_a,fit_loc,fit_b=ss.gamma.fit(butter_int_bin)
    # plt.hist(butter_int_bin, normed=True,bins=100)
    # plt.ion()
    # plt.plot(xs, ss.gamma.pdf(xs,fit_a,fit_loc,fit_b))

    mrf_coeff=fit_MRF_pseudolikelihood(adj_exc[op_mask][:,op_mask],
                                       adj_inh[op_mask][:,op_mask],
                                       Ys)
    
    embed()

    # import pymc
    
    # def pymc_simple(data):
    #     '''
    #     Single gamma model
    #     '''
    #     alpha=pymc.Gamma("alpha",alpha=1,beta=1)
    #     beta=pymc.Gamma("beta",alpha=1,beta=1)
    #     y=pymc.Gamma("y",alpha=alpha, beta=beta, value=data, observed=True)
    #     return locals()
    
    # def pymc_model(data):
    #     '''
    #     Mixed Gamma model
    #     '''
    #     p=pymc.Uniform("p",0,1)
    #     ber=pymc.Bernoulli("ber",p=p,size=len(data))
    #     alpha1=pymc.Gamma("alpha1",alpha=1,beta=1)
    #     alpha2=pymc.Gamma("alpha2",alpha=1,beta=1)
    #     beta1=pymc.Gamma("beta1",alpha=1,beta=1)
    #     beta2=pymc.Gamma("beta2",alpha=1,beta=1)
    #     @pymc.deterministic
    #     def par_alpha(ber=ber,alpha1=alpha1,alpha2=alpha2):
    #         return ber*alpha1 + (1-ber)*alpha2

    #     @pymc.deterministic
    #     def par_beta(ber=ber,beta1=beta1,beta2=beta2):
    #         return ber*beta1 + (1-ber)*beta2

    #     y=pymc.Gamma("y",alpha=par_alpha, beta=par_beta,
    #                  value=data, observed=True)
    #     return locals()

    # simple_fit=pymc.MCMC(pymc_simple(butter_int_bin))
    # sipmle_fit.
    # alpha_simp=simple_fit.stats()['alpha']['mean']
    # model_fit=pymc.MCMC(pymc_model(butter_int_bin))

    ## Save output
    scipy.io.savemat(outFn,
                     mdict={'bins': bins,
                            'spike_mat_bin': spike_mat_bin,
                            # 'spike_fil': spike_fil,
                            # 'butter_int': butter_int,
                            'spike_fil_bin': spike_fil_bin,
                            'butter_int_bin': butter_int_bin,
                            'psth_bin': psth_bin,
                            'chi': chi,
                            'autocorr': autocorr,
                            'peak_lag': peak_lag,
                            'peak_freq': peak_freq,
                            # 'duty_cycle': duty_cycle,
                            'ibi_mean': ibi_mean,
                            'ibi_cv': ibi_cv,
                            # 'burst_length_mean': burst_length_mean,
                            # 'burst_length_cv': burst_length_cv,
                            'ibi_vec': ibi_vec,
                            # 'burst_lengths': burst_lengths,
                            # 'burst_start_locs': burst_start_locs,
                            # 'burst_peak_locs': burst_peak_locs,
                            # 'burst_peaks': burst_peaks,
                            # 'bursting': bursting,
                            # 'bad_bursts': bad_bursts,
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
                            'pop_burst_peak': pop_burst_peak,
                            'avg_firing_rate': avg_firing_rate,
                            'num_silent': num_silent,
                            'num_expir': num_expir,
                            'num_inspir': num_inspir,
                            'num_tonic': num_tonic,
                            'eta_avg_inspir': avg_inspir,
                            'eta_avg_expir': avg_expir,
                            'eta_avg_tonic': avg_tonic,
                            'eta_avg_silent': avg_silent,
                            'mask_inspir': inspir_mask,
                            'mask_expir': expir_mask,
                            'mask_silent': silent_mask,
                            'mask_tonic': tonic_mask
                        },

                     oned_as='column')
    

# run the main stuff
if __name__ == '__main__':
    main()
