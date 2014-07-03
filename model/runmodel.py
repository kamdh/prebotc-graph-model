#!/usr/bin/env python
import os
import sys
import prebotc_BPR as prebotc
#import prebotc_weave as prebotc
import numpy as np
import scipy.sparse
import scipy.integrate
import time
import argparse
import scipy.io
import progressbar

def parse_args(argv):
    # # defaults for original model with seconds as time units
    # dt = 1e-3
    # tf = 30
    # defaults for BPR model
    dt = 1.0 # ms
    tf = 30.0e3 # ms
    # shared defaults
    t0 = 0.0 # ms
    abs_error = 1e-6
    rel_error = 1e-4
    spike_thresh = -15.0 # mV
    refractory = 6.0 # ms
    parser = argparse.ArgumentParser(prog="runmodel",
                                     description='run the preBotC model')
    parser.add_argument('-t0', type=float, default=t0,
                        help='initial time (default: %(default)s ms)')
    parser.add_argument('-tf', type=float, default=tf,
                        help='final time (default: %(default)s ms)')
    parser.add_argument('-dt', type=float, default=dt,
                        help='time step (default: %(default)s ms)')
    parser.add_argument('param', help='parameter pkl file')
    parser.add_argument('graph', help='graph gml file')
    parser.add_argument('output', help='output (.mat) filename')
    parser.add_argument('--abs_err', '-a', type=float, 
                        help='absolute error (default: %(default)s)',
                        default=abs_error)
    parser.add_argument('--rel_err', '-r', type=float, 
                        help='relative error (default: %(default)s)',
                        default=rel_error)
    parser.add_argument('--save_full', '-F', action='store_true',
                        help=('save all state variables '
                              '(default: store membrane potentials)'))
    parser.add_argument('--save_spikes', '-S', action='store_true',
                        help=('save just spike times '
                              '(default: store membrane potentials); '
                              'you should use dt < 10 ms'''))
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='silence output (for running in batch mode)')
    parser.add_argument('--spike_thresh', type=float,
                        default=spike_thresh,
                        help='spike threshold (default:%(default)s mV)')
    parser.add_argument('--refractory', type=float, default=refractory,
                        help='refractory period (default: %(default)s ms)')
    args = parser.parse_args(argv[1:])
    assert not ( args.save_spikes and args.save_full ), \
        "only one of --save_spikes and --save_full can be set"
    if args.save_spikes and args.dt > 1:
        print 'Warning, dt is possibly too large to resolve spikes'
    return args.t0, args.tf, args.dt, args.param, args.graph, \
        args.output, args.abs_err, args.rel_err, args.save_full, \
        args.save_spikes, args.quiet, args.spike_thresh, args.refractory

def main(argv=None):
    if argv is None:
        argv = sys.argv

    (t0, tf, dt, paramFn, graphFn, outFn, abs_error, rel_error, save_full, 
     save_spikes, quiet, spike_thresh, refractory) = parse_args(argv)
    # compute the number of steps required
    Nstep = np.ceil(tf/dt)
    if not quiet:
        print "Loading parameters, graph, and setting up IC's"
    my_params = prebotc.params(paramFn)
    num_vertices, num_edges, graph_params = prebotc.graph(graphFn)
    y, N = prebotc.ics(num_vertices, num_edges)
    #y, N = prebotc.ics(num_vertices, num_edges, random=False)
    # rhs of ODE with parameters evaluated
    # f is the rhs with parameters evaluated
    f = lambda t, y: prebotc.rhs(t, y, 
                                 graph_params,
                                 my_params)
    # data structure to output, timeseries or sparse raster
    if save_full:
        # all state variables
        save_state = np.zeros( (N, Nstep+1) )
    elif save_spikes:
        # just spike times, as sparse matrix
        # not boolean because type conversion bug in matlab output
        save_state = scipy.sparse.dok_matrix( (num_vertices, Nstep+1) )
        last_spike = np.ones( num_vertices ) * (-np.inf)
    else:
        # timeseries of spikes
        save_state = np.zeros( (num_vertices, Nstep+1) ) 
    r = scipy.integrate.ode(f)
    r.set_initial_value(y, t0)
    # # other integration methods
    # # method 1: BDF
    # r.set_integrator(
    #     'vode', 
    #     method='bdf', 
    #     with_jacobian = False,
    #     order=3,
    #     rtol= rel_error,
    #     atol= abs_error
    #     )
    # # method 2: Dormand-Price
    # r.set_integrator(
    #     'dopri5', 
    #     rtol = rel_error,
    #     atol = abs_error
    #     )
    # # method 3: VODE
    # r.set_integrator('vode',
    #                  rtol = rel_error,
    #                  atol = abs_error)
    if not quiet:
        print "Running integration loop...."
        t = time.time()
        bar_updates = 100
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()]
        bar = progressbar.ProgressBar(maxval=bar_updates, widgets=widgets)
        bar.start()
        j = 0
    i = 0
    while r.successful() and r.t < tf:
        r.integrate(r.t + dt)
        y = r.y.copy()
        if save_full:
            save_state[:, i] = y
        elif save_spikes:
            spikers = prebotc.spiking(y, num_vertices, spike_thresh)
            for neur in spikers:
                # only count if the new trigger occurs after reasonable delay
                if dt*( float(i) - last_spike[neur] ) >  refractory:
                    save_state[neur, i] = 1
                    last_spike[neur] = i
        else:
            save_state[:, i] = prebotc.voltages(y, num_vertices)
        i += 1
        if not quiet:
            if ( i % np.floor(Nstep/bar_updates) ) == 0:
                bar.update(j)
                j += 1
    if not save_spikes:
        save_state = save_state[:, 0:(i-1)]
    else:
        save_state.resize( (num_vertices, i-1) )
    if not quiet:
        bar.finish()
        elapsed = time.time() - t
        print "Done!\nElapsed: %1.2fs" % elapsed
        # Time saving
        t = time.time()
        print "Saving output...."
    if save_full:
        save_str = 'full'
    elif save_spikes:
        save_str = 'spikes'
    else:
        save_str = 'V'
    # save output
    scipy.io.savemat(outFn, 
                     mdict={'Y': save_state,
                            'dt': dt,
                            't0': t0,
                            'tf': tf,
                            'paramFn': os.path.abspath(paramFn),
                            'graphFn': os.path.abspath(graphFn),
                            'absErr': abs_error,
                            'relErr': rel_error,
                            'saveStr': save_str,
                            'finalState': y
                            },
                     oned_as = 'column')
    if not quiet:
        elapsed = time.time() - t
        print "Done!\nSave time: %1.2fs" % elapsed

# run the main stuff
if __name__ == '__main__':
    status = main()
    sys.exit(status)
