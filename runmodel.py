#!/usr/bin/env python
import sys
import prebotc_weave as prebotc
import numpy as np
import graph_tool as gt
import scipy.io
import scipy.integrate
import time
import argparse
import pickle

## constants
toolbar_width = 40
num_eqns_per_vertex = 7 #V, Na m, Na h, K n, hp Nap, Ca Can, Na pump
num_eqns_per_edge = 1

def main(argv=None):
    if argv is None:
        argv = sys.argv
    # input argument defaults
    dt = 1e-4
    t0 = 0
    tf = 30
    abs_error = 1e-10
    rel_error = 1e-9
    # parse arguments
    parser = argparse.ArgumentParser(prog="runmodel",
                                     description='run the preBotC model')
    parser.add_argument('-t0', type=float, default=t0,
                        help='initial time (default: %(default)s)')
    parser.add_argument('-tf', type=float, default=tf,
                        help='final time (default: %(default)s)')
    parser.add_argument('-dt', type=float, default=dt,
                        help='time step (default: %(default)s)')
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
                        help='save all state variables (default: just membrane potentials)')
    args = parser.parse_args(argv[1:])
    # store in usual variables
    t0 = args.t0
    tf = args.tf
    dt = args.dt
    paramFn = args.param
    graphFn = args.graph
    outFn = args.output
    abs_error = args.abs_err
    rel_error = args.rel_err
    save_full = args.save_full
    # compute the number of steps required
    Nstep = np.ceil(tf/dt)
    # load model parameters
    f = open(paramFn, 'r')
    my_params = pickle.load(f)
    f.close()
    # load graph topology
    print "Loading and processing graph"
    g = gt.load_graph(graphFn)
    g.reindex_edges()
    num_vertices = g.num_vertices()
    num_edges = g.num_edges()
    # store vertex types
    vertex_types = np.array( g.vertex_properties["type"].get_array(), 
                             dtype=np.int )
    # construct an edge list
    edge_list = np.zeros( (num_edges, 3) )
    # also a lookup table for in-edges
    # this requires a degree list
    in_degrees = np.array( g.degree_property_map("in").get_array(),
                           dtype=np.int )
    max_degree = np.max( in_degrees )
    # "ragged" array of in-edges
    if num_edges > 0:
        in_edges = np.zeros( (num_vertices, max_degree), dtype=np.int )
        gsyn_props = g.edge_properties["gsyn"]
    else:
        in_edges = np.zeros( (num_vertices, max_degree), dtype=np.int )
        gsyn_props = []
    # for looping
    in_edge_ct = np.zeros( (num_vertices,), dtype=np.int )
    i = 0
    for e in g.edges():
        source_index = int( e.source() )
        target_index = int( e.target() )
        edge_list[i,...] = [source_index, 
                            target_index,
                            gsyn_props[e]]
        in_edges[ target_index, in_edge_ct[target_index] ] = i
        # increment indices
        in_edge_ct[ target_index ] += 1
        i += 1
    ## setup initial conditions
    # state will contain vertex variables & edge
    # variables in a 1d array
    N = num_vertices*num_eqns_per_vertex +\
        num_edges*num_eqns_per_edge
    # state vector y encodes vertex and edge data
    y = np.zeros(N)
    for i in range( num_vertices ):
        # vertex data in 0:num_eqns_per_vertex*num_vertices-1
        j = range(i*num_eqns_per_vertex, (i+1)*num_eqns_per_vertex)
        y[j] = [
            -0.026185387764343,
             0.318012107836673,
             0.760361103277830,
             0.681987892188221,
             0.025686471226045,
             0.050058183820371,
             4.998888741335261
             ]
    offset = num_vertices*num_eqns_per_vertex
    for i in range( num_edges ):
        j = range(offset + i*num_eqns_per_edge,
                  offset + (i+1)*num_eqns_per_edge)
        y[j] = 0.000001090946631
    
    # f is the rhs with parameters evaluated
    def f(t, y):
        dydt = prebotc.rhs(t, y, 
                           vertex_types = vertex_types,
                           edge_list = edge_list, 
                           in_degrees = in_edge_ct,
                           in_edges = in_edges,
                           **my_params)
        return dydt

  
    #### integration
    ## scipy.odeint
    # tRange = np.arange(t0, tf+dt, dt)
    # save_state = scipy.integrate.odeint(f, y, tRange,
    #                                     rtol = rel_error,
    #                                     atol = abs_error)
    ## scipy.integrate.ode
    # output vector of states
    if save_full:
        save_state = np.zeros( (N, Nstep+1) )
    else:
        save_state = np.zeros( (num_vertices, Nstep+1) ) 
    r = scipy.integrate.ode(f)
    # method 1: BDF
    # r.set_integrator(
    #     'vode', 
    #     method='bdf', 
    #     with_jacobian = False,
    #     order=3,
    #     rtol= rel_error,
    #     atol= abs_error
    #     )
    # method 2: Dormand-Price
    # r.set_integrator(
    #     'dop853', 
    #     rtol = rel_error,
    #     atol = abs_error
    #     )
    # method 3: LSODE
    # r.set_integrator('lsode',
    #                  rtol = rel_error,
    #                  atol = abs_error)

    ## integrate the ODE
    print "Running integration loop...."
    # start timing integration
    t = time.time()
    # setup simple toolbar
    # from http://stackoverflow.com/questions/3160699/python-progress-bar
    global toolbar_width
    sys.stdout.write("[%s]" % (" " * toolbar_width ))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to 
    r.set_initial_value(y, t0)
    i = 0
    while r.successful() and r.t < tf:
        r.integrate(r.t + dt)
        y = r.y.copy()
        if save_full:
            save_state[:, i] = y
        else:
            save_state[:, i] = y[ 0:(num_vertices*num_eqns_per_vertex):num_eqns_per_vertex ]
        i += 1
        if ( i % (Nstep/toolbar_width) )==0 :
            sys.stdout.write('#')
            sys.stdout.flush()
        # if ( (i)%report_every ) == 0:
        #     print "%1.2f" % r.t
    save_state = save_state[:, 0:(i-1)]
    elapsed = time.time() - t
    print "\nDone!\nElapsed: %1.2fs" % elapsed
    # Time saving
    t = time.time()
    print "Saving output...."
    scipy.io.savemat(outFn, 
                     mdict={'Y': save_state,
                            'vTypes': vertex_types},
                     oned_as = 'column')
    elapsed = time.time() - t
    print "Done!\nSave time: %1.2fs" % elapsed

# run the main stuff
if __name__ == '__main__':
    status = main()
    sys.exit(status)
