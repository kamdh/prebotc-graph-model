#!/usr/bin/env python
# Kameron Decker Harris
# pre-BotC model ODEs

import numpy as np
import networkx as nx
import json
import os
from scipy import weave
from scipy.weave import converters

# constants
num_eqns_per_vertex = 7 # V, Na m, Na h, K n, hp Nap, Ca Can, Na pump
num_eqns_per_edge = 1   # not necessary at the moment

def params(paramFn):
    with open(paramFn) as f:
        my_params = json.load(f)
    return my_params

def graph(graphFn):
    g = nx.read_gml(graphFn)
    #g.reindex_edges()
    num_vertices = g.number_of_nodes()
    num_edges = g.number_of_edges()
    # store vertex types
    vertex_types = np.array( nx.get_node_attributes(g, 'type').values(),
                             dtype=np.int )
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
        edge_list[i,...] = [source_index, 
                            target_index,
                            gsyn_props[e]]
        in_edges[ target_index, in_edge_ct[target_index] ] = i
        # increment indices
        in_edge_ct[ target_index ] += 1
        i += 1
    
    return num_vertices, num_edges, vertex_types, edge_list, \
        in_edge_ct, in_edges

def ics(num_vertices, num_edges):
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
    return y, N

def voltages(y, num_vertices):
    V = y[ 0:(num_vertices*num_eqns_per_vertex):num_eqns_per_vertex ]
    return V

def rhs(
    t, # time
    y, # state variables
    vertex_types, 
    edge_list, 
    in_degrees,
    in_edges,
    params
    ):

    EL = params['EL']
    gCaNS = params['gCaNS']
    gPS = params['gPS']
    ELS = params['ELS']
    gCaNI = params['gCaNI']
    gPI = params['gPI']
    ELI = params['ELI']
    gCaNTS = params['gCaNTS']
    gPTS = params['gPTS']
    ELTS = params['ELTS']
    gCaNSil = params['gCaNSil']
    gPSil = params['gPSil']
    ELSil = params['ELSil']
    alpha = params['alpha']
    Cab = params['Cab']
    Cm = params['Cm']
    EK = params['EK']
    eCa = params['eCa']
    ehp = params['ehp']
    ECaN = params['ECaN']
    Esyn = params['Esyn']
    gK = params['gK']
    gL = params['gL']
    gNa = params['gNa']
    Iapp = params['Iapp']
    kIP3 = params['kIP3']
    ks = params['ks']
    kNa = params['kNa']
    kCa = params['kCa']
    kCAN = params['kCAN']
    Nab = params['Nab']
    siCAN = params['siCAN']
    sih = params['sih']
    sihp = params['sihp']
    sim = params['sim']
    simp = params['simp']
    siN = params['siN']
    sis = params['sis']
    tauh = params['tauh']
    tauhp = params['tauhp']
    taum = params['taum']
    taun = params['taun']
    taus = params['taus']
    Qh = params['Qh']
    Qhp = params['Qhp']
    Qm = params['Qm']
    Qmp = params['Qmp']
    Qn = params['Qn']
    Qs = params['Qs']
    ENa = params['ENa']
    gsyn = params['gsyn']

    dydt = np.zeros(y.shape[0]) # initialize vector field
    
    support_code = '''
extern "C"{
double infHN( double A, double B, double V );
double tau( double A, double B, double C, double V );
double phi( double x, double kNa );
}

double infHN( double A, double B, double V ) {
  return 1.0/(1.0 + exp( (V - A)/B ));
}

double tau( double A, double B, double C, double V ) {
  return A / cosh( 0.5 * (V - B)/C );
}
    
double phi( double x, double kNa ) {
  return pow(x,3) / (pow(x, 3) + pow(kNa, 3));
}
'''
    
    code = """
//// edit

int num_vertices = Nvertex_types[0];
int num_edges = Nedge_list[0];
int offset = num_vertices*num_eqns_per_vertex;
int i,j,k;
double gCaN_array[] = {gCaNS, gCaNI, gCaNTS, gCaNSil};
double gP_array[]   = {gPS, gPI, gPTS, gPSil};
double EL_array[]   = {ELS, ELI, ELTS, ELSil};
//// vertex variables
for (i=0; i<num_vertices; i++) {
  double EL_cur, gP_cur, gCaN_cur, I_L, I_Na, I_K,
    I_NaP, I_CaN, I_pump, I_syn;
  int type_idx;
  j = i*num_eqns_per_vertex; // index of first (V) variable
  type_idx = (int) vertex_types(i); // neuron type
  //// type-dependent values
  EL_cur = EL_array[type_idx]; // leak potential
  gP_cur = gP_array[type_idx]; // Na p conductance
  gCaN_cur = gCaN_array[type_idx]; // Ca n conductance
  //// calculate currents
  // I_L(V)
  I_L = gL * ( y(j) - EL_cur );
  // I_Na(V, h, m)
  I_Na = gNa * pow(y(j+1),3) * y(j+2) * (y(j) - ENa);
  // I_K(V, n)
  I_K = gK * pow(y(j+3),4) *(y(j) - EK);
  // I_NaP(V, hp)
  I_NaP = gP_cur * infHN( Qmp, simp, y(j) ) * y(j+4) * (y(j) - ENa);
  // I_CaN(V, Ca)
  I_CaN = gCaN_cur * (y(j) - ECaN) * infHN( kCAN, siCAN, y(j+5) );
  // I_pump(Na)
  I_pump = 0.2 * ( phi( y(j+6), kNa ) - phi( Nab, kNa ) );
  //// calculate synaptic current I_syn
  double gSynTot = 0.0;
  double synVarMean = 0.0;
  for (k=0; k < (int)in_degrees(i); k++) {
    int edgeid = (int) in_edges(i,k);
    gSynTot += y(offset + edgeid) *
      edge_list( edgeid, 2 );
    synVarMean += y(offset + edgeid);
  }
  if ( in_degrees(i) > 0 ) {
    // Ohm's law
    I_syn = gSynTot * (y(j) - Esyn);
    synVarMean = synVarMean / in_degrees(i);
  } else {
    I_syn = 0.0;
    synVarMean = 0.0;
  }
  //// set the derivatives
  // voltage V
  dydt(j) = -( Iapp + 
               I_L +
               I_Na + 
               I_K + 
               I_NaP + 
               I_CaN + 
               I_pump +
               I_syn ) / Cm;
  // Na m
  dydt(j+1) = ( infHN(Qm,sim,y(j)) - y(j+1) ) / tau(taum, Qm, sim, y(j));
  // Na h
  dydt(j+2) = ( infHN(Qh,sih,y(j)) - y(j+2) ) / tau(tauh, Qh, sih,y(j));
  // K n
  dydt(j+3) = ( infHN(Qn,siN,y(j)) - y(j+3) ) / tau(taun, Qn, siN,y(j));
  // hp Nap
  dydt(j+4) = ehp * ( infHN( Qhp, sihp, y(j) ) - y(j+4) ) /
    tau( tauhp, Qhp, sihp, y(j) );
  // Ca Can
  dydt(j+5) = eCa * ( kIP3 * synVarMean -
                      kCa * ( y(j+5) - Cab ) );
  // Na pump
  dydt(j+6) = alpha * ( - I_CaN - I_pump );
 }
// edge variables
for (i=0; i<num_edges; i++) {
  double v_source, v_target;
  int sourceid = (int) edge_list(i,0);
  v_source = y( sourceid * num_eqns_per_vertex );
  //v_target = edge_list(i, 1);
  j = offset + i;
  dydt(j) = ( (1 - y(j)) *
              infHN( Qs, sis, v_source ) -
              ks * y(j) ) / taus;
}

"""


    weave.inline(code, ['t', 'y', 'vertex_types','edge_list',
                        'in_degrees', 'in_edges', 'EL', 'gCaNS', 
                        'gPS', 'ELS', 'gCaNI', 'gPI', 'ELI', 
                        'gCaNTS', 'gPTS', 'ELTS', 'gCaNSil', 
                        'gPSil', 'ELSil', 'alpha', 'Cab', 'Cm', 
                        'EK', 'eCa', 'ehp', 'ECaN', 'Esyn', 
                        'gK', 'gL', 'gNa', 'Iapp', 'kIP3', 'ks',
                        'kNa', 'kCa', 'kCAN', 'Nab', 'siCAN', 
                        'sih', 'sihp', 'sim', 'simp', 'siN', 'sis', 
                        'tauh', 'tauhp', 'taum', 'taun', 'taus', 
                        'Qh', 'Qhp', 'Qm', 'Qmp', 'Qn', 'Qs', 
                        'ENa', 'gsyn', 'num_eqns_per_vertex',
                        'dydt'],
                 support_code = support_code,
                 verbose = 1,
                 type_converters = converters.blitz, 
                 compiler='gcc',
                 # library_dirs = [os.getcwd()],
                 # libraries=['helpers'],
                 # runtime_library_dirs = [os.getcwd()],
                 headers=['<math.h>']
                 )

    
    return dydt

# extra business
def _test():
    return 0

if __name__ == "__main__":
    _test()
