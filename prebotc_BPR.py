#!/usr/bin/env python
# Kameron Decker Harris
#
# pre-BotC model ODEs from Park & Rubin, 2013
# so-called Butera, Park, Rubin model

import numpy as np
import networkx as nx
import json
import os
from scipy import weave
from scipy.weave import converters
import math  # allow functions in math module

# constants
num_eqns_per_vertex = 5
num_eqns_per_edge = 1

def params(paramFn):
    def parse_I_aps(var, paramFn):
        if isinstance(var, (int, float, long)):
            # then it is constant
            fun = lambda t: var
        elif isinstance(var, basestring):
            # then it should be python code for a function
            code = "lambda t: " + var
            # below is not safe, but we trust the param file
            bytecode = compile( code, paramFn, 'eval' )
            fun = eval( bytecode )
        else:
            raise Exception("wrong type for I_aps in " + paramFn)
        return fun
    with open(paramFn) as f:
        my_params = json.load(f)
    my_params["I_apsE"] = parse_I_aps( my_params["I_apsE"], paramFn )
    my_params["I_apsI"] = parse_I_aps( my_params["I_apsI"], paramFn )
    return my_params

def graph(graphFn):
    g = nx.read_gml(graphFn)
    #g.reindex_edges()
    num_vertices = g.number_of_nodes()
    num_edges = g.number_of_edges()
    # store vertex types
    vertex_types = np.array(nx.get_node_attributes(g, 'type').values(),
                            dtype=np.int)
    vertex_inh = np.array(nx.get_node_attributes(g, 'inh').values(),
                          dtype=np.int )
    vertex_respir_area = np.array(
        nx.get_node_attributes(g, 'respir_area').values(), dtype=np.int)
    # construct an edge list
    edge_list = np.zeros( (num_edges, 3) )
    # also a lookup table for in-edges
    # this requires a degree list
    in_degrees = np.array( g.in_degree().values(), dtype=np.int )
    max_degree = np.max( in_degrees )
    # "ragged" array of in-edges
    in_edges = np.zeros( (num_vertices, max_degree), dtype=np.int )
    if num_edges > 0:
        gsyn_props = nx.get_edge_attributes(g, 'gsyn')
    else:
        gsyn_props = []
    # for looping
    in_edge_ct = np.zeros( (num_vertices,), dtype=np.int )
    i = 0
    for e in g.edges():
        source_index = int( e[0] )
        target_index = int( e[1] )
        is_inh = vertex_inh[ source_index ]
        if is_inh == 1:
            inh_mult = -1
        else:
            inh_mult = 1
        edge_list[i,...] = [ source_index, 
                             target_index,
                             inh_mult * abs(gsyn_props[e]) ]
        in_edges[ target_index, in_edge_ct[target_index] ] = i
        # increment indices
        in_edge_ct[ target_index ] += 1
        i += 1
    graph_params = (vertex_types, vertex_inh, vertex_respir_area,
                    edge_list, in_edge_ct, in_edges)
    return num_vertices, num_edges, graph_params

def ics(num_vertices, num_edges, random=True):
    # state will contain vertex variables & edge
    # variables in a 1d array
    N = num_vertices*num_eqns_per_vertex + num_edges*num_eqns_per_edge
    # state vector y encodes vertex and edge data
    y = np.zeros(N)
    if random:
        for i in range( num_vertices ):
            # vertex data in 0:num_eqns_per_vertex*num_vertices-1
            j = range(i*num_eqns_per_vertex, (i+1)*num_eqns_per_vertex)
            # ranges below come from range of variables in sample output run
            y[j] = [
                (7.06 + 53.6) * np.random.random_sample() - 53.6,
                (0.95 - 0.00) * np.random.random_sample() + 0.00,
                (0.64 - 0.21) * np.random.random_sample() + 0.21,
                (0.95 - 0.02) * np.random.random_sample() + 0.02,
                (0.93 - 0.46) * np.random.random_sample() + 0.46
            ]
        offset = num_vertices*num_eqns_per_vertex
        for i in range( num_edges ):
            j = range(offset + i*num_eqns_per_edge,
                      offset + (i+1)*num_eqns_per_edge)
            y[j] = (0.0025 - 0) * np.random.random_sample() + 0
    else:
        for i in range( num_vertices ):
            # vertex data in 0:num_eqns_per_vertex*num_vertices-1
            j = range(i*num_eqns_per_vertex, (i+1)*num_eqns_per_vertex)
            y[j] = [
                -50,
                 0.004,
                 0.33,
                 0.03,
                 0.93
                 ]
        offset = num_vertices*num_eqns_per_vertex
        for i in range( num_edges ):
            j = range(offset + i*num_eqns_per_edge,
                      offset + (i+1)*num_eqns_per_edge)
            y[j] = 0.0000011
    return y, N

def voltages(y, num_vertices):
    V = y[ 0:(num_vertices*num_eqns_per_vertex):num_eqns_per_vertex ]
    return V

def spiking(y, num_vertices, thresh):
    V = voltages(y, num_vertices)
    spiking_neurons = np.where( V > thresh )[0]
    return spiking_neurons

def rhs( t, y,
         graph_params,
         params
         ):
    # unpack all of the parameters
    (vertex_types, vertex_inh, vertex_respir_area, edge_list,
     in_degrees, in_edges) = graph_params
    Cms = params['Cms']
    I_apsE = params['I_apsE']
    I_apsI = params['I_apsI']
    vna = params['vna']
    vk = params['vk']
    vleaks = params['vleaks']
    vsyn = params['vsyn']
    vsynE = params['vsynE']
    vsynI = params['vsynI']
    vm = params['vm']
    vn = params['vn']
    vmp = params['vmp']
    vh = params['vh']
    sm = params['sm']
    sn = params['sn']
    smp = params['smp']
    sh = params['sh']
    ssyn = params['ssyn']
    taunb = params['taunb']
    tauhb = params['tauhb']
    tausyn = params['tausyn']
    gk = params['gk']
    gna = params['gna']
    gnap = params['gnap']
    gl_CS = params['gl_CS']
    gl_CI = params['gl_CI']
    gl_TS = params['gl_TS']
    gl_Q = params['gl_Q']
    Kcan = params['Kcan']
    ncan = params['ncan']
    gcan_CS = params['gcan_CS']
    gcan_CI = params['gcan_CI']
    gcan_TS = params['gcan_TS']
    gcan_Q = params['gcan_Q']
    I = params['I']
    Ct = params['Ct']
    fi = params['fi']
    LL = params['LL']
    P = params['P']
    Ki = params['Ki']
    Ka = params['Ka']
    Ve = params['Ve']
    Ke = params['Ke']
    A = params['A']
    Kd = params['Kd']
    sigma = params['sigma']
    ksyn = params['ksyn']
    # initialize vector field
    dydt = np.zeros(y.shape[0])
    # evaluate time-varying input currents
    I_apsE_t = I_apsE(t)
    I_apsI_t = I_apsI(t)
    code = """
int num_vertices = Nvertex_types[0];
int num_edges = Nedge_list[0];
int offset = num_vertices*num_eqns_per_vertex;
int i,j,k;
double gcan_array[] = {gcan_CS, gcan_CI, gcan_TS, gcan_Q};
double gl_array[]   = {gl_CS, gl_CI, gl_TS, gl_Q};
//// vertex variables
for (i=0; i<num_vertices; i++) {
  double minf, ninf, minfp, hinf, taun, tauh, I_na, I_k, I_nap, I_l,
    caninf, I_can, J_ER_in, J_ER_out, Ce, I_syn, gcan, gl, I_aps;
  int type_idx;
  int is_inh;
  j = i*num_eqns_per_vertex; // index of first (V) variable
  type_idx = (int) vertex_types(i); // neuron type
  is_inh = (int) vertex_inh(i); // whether or not it's inhibitory
  //// type-dependent values
  gcan = gcan_array[type_idx];
  gl = gl_array[type_idx];
  if (is_inh == 0) {
    I_aps = (double) I_apsE_t;
  } 
  else if (is_inh == 1) {
    I_aps = (double) I_apsI_t;
  } 
  else {
    throw \"is_inh returned nonsensical result\";
  }
  //// activation/inactivation variables
  minf  = 1/(1+exp(( y(j) - vm )/sm));
  ninf  = 1/(1+exp(( y(j) - vn )/sn));
  minfp = 1/(1+exp(( y(j) - vmp )/smp));
  hinf  = 1/(1+exp(( y(j) - vh )/sh));
  //// time constants
  taun = taunb/cosh(( y(j) - vn )/(2*sn));
  tauh = tauhb/cosh(( y(j) - vh )/(2*sh));
  //// calculate currents
  I_na  = gna * pow(minf, 3) * (1 - y(j+1)) * (y(j) - vna);
  I_k   = gk * pow(y(j+1), 4) * (y(j) - vk);
  I_nap = gnap * minfp * y(j+2) * (y(j) - vna);
  I_l   = gl * (y(j) - vleaks);
  //// CaN current
  caninf = 1/(1+pow(Kcan/y(j+3), ncan));
  I_can = gcan * caninf * (y(j) - vna);
  //// Fluxes in and out of ER
  //// l is fraction of open IP3 channels
  Ce = (Ct - y(j+3))/sigma;
  J_ER_in=(LL + P * pow( (I*y(j+3)*y(j+4))/( (I+Ki)*(y(j+3)+Ka) ), 3) )
             *(Ce - y(j+3));
  J_ER_out=Ve * pow(y(j+3), 2) / ( pow(Ke, 2) + pow(y(j+3), 2) );
  //// calculate synaptic current I_syn
  I_syn = 0.0;
  for (k=0; k < (int)in_degrees(i); k++) {
    int edgeid = (int) in_edges(i,k);
    double gSyn = edge_list( edgeid, 2 );
    double synVar = y( offset + edgeid );
    if ( gSyn < 0.0 ) {
      // inhibitory
      I_syn += abs(gSyn) * synVar * ( y(j) - vsynI );
    } 
    else {
      // excitatory
      I_syn += abs(gSyn) * synVar * ( y(j) - vsynE );
    }
  }
  if ( (int)in_degrees(i) == 0 ) {
    I_syn = 0.0;
  }
  //// set the derivatives
  // v'
  dydt(j) = -( I_k + I_na + I_nap + I_l - I_aps + I_can + I_syn ) / Cms;
  // n'
  dydt(j+1) = (ninf - y(j+1))/taun;
  // h'
  dydt(j+2) = (hinf - y(j+2))/tauh;
  // C'
  dydt(j+3) = fi * (J_ER_in - J_ER_out);
  // l'
  dydt(j+4) = A * (Kd - (y(j+3) + Kd) * y(j+4));
 }
// edge variables
for (i=0; i<num_edges; i++) {
  double v_source, v_target, msyninf;
  int sourceid = (int) edge_list(i, 0);
  double gsyn = edge_list(i, 2);
  v_source = y( sourceid * num_eqns_per_vertex );
  //v_target = edge_list(i, 1);
  j = offset + i;
  msyninf = 1/( 1 + exp( ( v_source - vsyn )/ssyn ) );
  dydt(j) = ( (1 - y(j)) * msyninf - ksyn * y(j) ) / tausyn;
}
"""
    weave.inline(code, 
                 ['t', 'y', 'vertex_types', 'vertex_inh', 'edge_list', 
                  'in_degrees', 'in_edges', 'Cms', 'I_apsE_t', 'I_apsI_t', 
                  'vna', 'vk', 'vleaks', 'vsynE', 'vsynI', 'vsyn', 'vm', 'vn',
                  'vmp', 'vh', 'sm', 'sn', 'smp', 'sh', 'ssyn', 'taunb', 
                  'tauhb', 'tausyn', 'gk', 'gna', 'gnap', 
                  'gl_CS', 'gl_CI', 'gl_TS', 'gl_Q', 'Kcan', 'ncan', 
                  'gcan_CS', 'gcan_CI', 'gcan_TS', 'gcan_Q', 'I', 'Ct', 'fi',
                  'LL', 'P', 'Ki', 'Ka', 'Ve', 'Ke', 
                  'A', 'Kd', 'sigma', 'ksyn', 'num_eqns_per_vertex', 'dydt'],
                 verbose = 1,
                 type_converters = converters.blitz, 
                 compiler='gcc',
                 headers=['<math.h>']
                 )
    return dydt

# extra business
def _test():
    return 0

if __name__ == "__main__":
    _test()

