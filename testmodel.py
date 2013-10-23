#!/usr/bin/env python

import sys
import prebotc_pure as prebotc
#import prebotc_cython as prebotc
import numpy as np
import graph_tool as gt
import scipy.io
import scipy.integrate

outfn = 'test.mat'
dt = 5e-5
t0 = 0.0
tf = 1.0
Nstep = int(round(tf/dt))
report_every = 100
num_eqns_per_vertex = 7 #V, Na m, Na h, K n, hp Nap, Ca Can, Na pump
num_eqns_per_edge = 1
abs_error = 1e-6
rel_error = 1e-4
my_params = dict(
    #gP = 4.0,
    EL = -0.0625,
    gCaNS = 10.0,
    gPS = 0.5,
    ELS = -0.06,
    gCaNI = 0.0,
    gPI = 4.0,
    ELI = -0.0625,
    gCaNTS = 0.0,
    gPTS = 5.,
    ELTS = -0.062,
    gCaNSil = 0.0,
    gPSil = 0.6,
    ELSil = -0.0605,
    alpha = 6.6e-2,
    Cab = 0.05,
    Cm = 0.045,
    EK = -0.075,
    eCa = 0.0007,
    ehp = 0.001,
    ECaN = 0.0,
    Esyn = 0.0,
    gK = 30.0,
    gL = 3.0,
    gNa = 160.0,
    Iapp = 0.0,
    kIP3 = 1.2e+6,
    ks = 1.0,
    kNa = 10.0,
    kCa = 22.5e+3,
    kCAN = 0.9,
    Nab = 5.0,
    rp = 0.2,
    siCAN = -0.05,
    sih = 0.005,
    sihp = 0.006,
    sim = -0.0085,
    simp = -0.006,
    siN = -0.005,
    sis = -0.003,
    tauh = 0.015,
    tauhp = 0.001,
    taum = 0.001,
    taun = 0.030,
    taus = 0.015,
    Qh = -0.030,
    Qhp = -0.048,
    Qm = -0.036,
    Qmp = -0.040,
    Qn = -0.030,
    Qs = 0.015,
    ENa = 0.045,
    gsyn = 2.5
    )

def main(argv=None):
    # parse arguments (not used yet)
    if argv is None:
        argv = sys.argv

    # load T Dashevskiy's graph topology
    g = gt.load_graph("Dashevskiy.gml")
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
    in_edges = np.zeros( (num_vertices, max_degree), dtype=np.int )
    gsyn_props = g.edge_properties["gsyn"]
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
    
    # f is the rhs with parameters evaluated
    def f(t, y):
        dydt = prebotc.rhs(t, y, 
                           vertex_types = vertex_types,
                           edge_list = edge_list, 
                           in_degrees = in_edge_ct,
                           in_edges = in_edges,
                           **my_params)
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
            
    scipy.io.savemat(outfn, mdict={'Y': save_state},
                     oned_as = 'col')
    
# run the main stuff
if __name__ == '__main__':
    status = main()
    sys.exit(status)
