#!/usr/bin/env python

import sys
#import prebotc_pure as prebotc
#import prebotc_cython as prebotc
import prebotc_weave as prebotc
import numpy as np
import graph_tool as gt
import scipy.io
import scipy.integrate
import pickle

paramFn = 'param_files/test.pkl'
outFn = 'output/test.mat'
graphFn = '../graphs/test.gml'
dt = 1e-4
t0 = 0.0
tf = 5
Nstep = int(round(tf/dt))
report_every = 1000
num_eqns_per_vertex = 7 #V, Na m, Na h, K n, hp Nap, Ca Can, Na pump
num_eqns_per_edge = 1
abs_error = 1e-9
rel_error = 1e-8

def main(argv=None):
    # parse arguments (not used yet)
    if argv is None:
        argv = sys.argv
    # load parameters
    f = open(paramFn, 'r')
    my_params = pickle.load(f)
    f.close()
    # load graph topology
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
    if num_edges > 0:
        # "ragged" array of in-edges
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
        #print(j)
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
        #print(j)
        y[j] = 0.000001090946631
    #print(N)
    print y
    
    # f is the rhs with parameters evaluated
    def f(t, y):
        dydt = prebotc.rhs(t, y, 
                           vertex_types,
                           edge_list, 
                           in_edge_ct,
                           in_edges,
                           my_params)
        return dydt
    
    # output vector of states
    save_state = np.zeros( (N, Nstep) ) 

    ## hard-coded Euler method
    t = t0;
    for i in range(Nstep):
        dydt = f(t, y)
        y = y + dydt * dt # fwd Euler
        #save_state[:, i] = y[ 0:(num_vertices*num_eqns_per_vertex):num_eqns_per_vertex ] # just voltages
        save_state[:, i] = y; # all vars
        t = t + dt;
        if ( (i+1)%report_every ) == 0:
            print t
            
    scipy.io.savemat(outFn, mdict={'Y': save_state},
                     oned_as = 'col')
    
# run the main stuff
if __name__ == '__main__':
    status = main()
    sys.exit(status)
