#!/usr/bin/env python
# Kameron Decker Harris
# pre-BotC model ODEs

import numpy as np

# constants
num_eqns_per_vertex = 7 #V, Na m, Na h, K n, hp Nap, Ca Can, Na pump
num_eqns_per_edge = 1

def rhs(
    t, # time
    y, # state variables
    vertex_types, 
    edge_list, 
    in_degrees,
    in_edges,
    EL, gCaNS, 
    gPS, ELS, gCaNI, gPI,
    ELI, gCaNTS, gPTS,
    ELTS, gCaNSil, gPSil,
    ELSil, alpha, Cab, 
    Cm, EK, eCa, ehp, ECaN,
    Esyn, gK, gL, gNa, Iapp, 
    kIP3, ks, kNa, kCa, 
    kCAN, Nab,
    siCAN, sih, sihp, sim, 
    simp, siN, sis, tauh, 
    tauhp, taum, taun,
    taus, Qh, Qhp, Qm,
    Qmp, Qn, Qs, ENa, gsyn
    ):

    ## setup arrays
    # Ca n conductances
    gCaN_array = (gCaNS, gCaNI, gCaNTS, gCaNSil);
    # Na p conductances
    gP_array   = (gPS, gPI, gPTS, gPSil);
    # leak potentials
    EL_array   = (ELS, ELI, ELTS, ELSil);

    ## helper functions
    def infHN(A, B, V):
        return 1.0/(1.0 + np.exp( (V - A)/B ))

    def tau( A, B, C, V ):
        return A / np.cosh( 0.5 * (V-B)/C )
    
    def phi( x ):
        return x**3 / (x**3 + kNa**3)

    ## main routine
    dydt = np.zeros(y.shape[0]) # initialize vector field
    num_vertices = vertex_types.shape[0]
    num_edges = edge_list.shape[0]
    offset = num_vertices*num_eqns_per_vertex # for indexing edges
    # vertex variables
    for i in range( num_vertices ):
        j = i*num_eqns_per_vertex # index of first (V) variable
        type_idx = int(vertex_types[i]) # neuron type
        ## type-dependent values
        EL_cur = EL_array[type_idx] # leak potential
        gP_cur = gP_array[type_idx] # Na p conductance
        gCaN_cur = gCaN_array[type_idx] # Ca n conductance
        ## calculate currents
        # I_L(V)
        I_L = gL * ( y[j] - EL_cur )
        # I_Na(V, h, m)
        I_Na = gNa * y[j+1]**3 * y[j+2] * (y[j] - ENa)
        # I_K(V, n)
        I_K = gK * y[j+3]**4 *(y[j] - EK)
        # I_NaP(V, hp)
        I_NaP = gP_cur * infHN( Qmp, simp, y[j] ) * y[j+4] * (y[j] - ENa)
        # I_CaN(V, Ca)
        I_CaN = gCaN_cur * (y[j] - ECaN) * infHN( kCAN, siCAN, y[j+5] )
        # I_pump(Na)
        I_pump = 0.2 * ( phi( y[j+6] ) - phi( Nab ) )
        ## calculate synaptic current I_syn
        if in_degrees[i] > 0:
            # first get in-edges
            these_in_edges = in_edges[i, 0:(in_degrees[i]-1)]
            # find the synaptic variables for those edges
            syn_variables = y[offset + these_in_edges]
            # and the conductances
            syn_conductances = edge_list[ these_in_edges, 2 ]
            # Ohm's law
            I_syn = np.dot(syn_variables, syn_conductances) * (y[j] - Esyn)
        else:
            I_syn = 0
            syn_variables = (0)
        ## set the derivatives
        # voltage V
        dydt[j] = -( Iapp + 
                     I_L +
                     I_Na + 
                     I_K + 
                     I_NaP + 
                     I_CaN + 
                     I_pump +
                     I_syn ) / Cm
        # Na m
        dydt[j+1] = ( infHN(Qm,sim,y[j]) - y[j+1] ) / tau(taum, Qm, sim, y[j])
        # Na h
        dydt[j+2] = ( infHN(Qh,sih,y[j]) - y[j+2] ) / tau(tauh, Qh, sih,y[j])
        # K n
        dydt[j+3] = ( infHN(Qn,siN,y[j]) - y[j+3] ) / tau(taun, Qn, siN,y[j])
        # hp Nap
        dydt[j+4] = ehp * ( infHN(Qhp, sihp, y[j]) - y[j+4] ) /\
            tau(tauhp, Qhp, sihp, y[j])
        # Ca Can
        # dydt[j+5] = eCa * ( kIP3 * syn_variables[0] - \
        #                         kCa * ( y[j+5] - Cab ) )
        dydt[j+5] = eCa * ( kIP3 * np.mean(syn_variables) - \
                                kCa * ( y[j+5] - Cab ) )
        # Na pump
        dydt[j+6] = alpha * ( - I_CaN - I_pump )
    # edge variables
    for i in range( num_edges ):
        this_edge = edge_list[i,:]
        v_source = this_edge[0]
        v_target = this_edge[1]
        j = range(offset + i*num_eqns_per_edge,
                  offset + (i+1)*num_eqns_per_edge)
        dydt[j] = ( (1 - y[j]) * \
                        infHN( Qs, sis, y[ v_source*num_eqns_per_vertex ] ) - \
                        ks * y[j] ) / taus
    return dydt

# extra business
def _test():
    return 0

if __name__ == "__main__":
    _test()
