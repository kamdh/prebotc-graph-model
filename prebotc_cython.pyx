#!/usr/bin/env python
# Kameron Decker Harris
# pre-BotC model ODEs

import numpy as np
cimport numpy as np
cimport cython


# cython datatypes
DTYPE = np.float
ctypedef np.float_t DTYPE_t
INT = np.int
ctypedef np.int_t INT_t

@cython.boundscheck(False) # turn of bounds-checking for entire function
def rhs(
    np.float t, # time
    np.ndarray[DTYPE_t, ndim=1] y, # state variables
    np.ndarray[INT_t, ndim=1] vertex_types, 
    np.ndarray[DTYPE_t, ndim=2] edge_list, 
    np.ndarray[INT_t, ndim=1] in_degrees,
    np.ndarray[INT_t, ndim=2] in_edges, # graph
    np.float EL, np.float gCaNS, 
    np.float gPS, np.float ELS, np.float gCaNI, np.float gPI,
    np.float ELI, np.float gCaNTS, np.float gPTS,
    np.float ELTS, np.float gCaNSil, np.float gPSil,
    np.float ELSil, np.float alpha, np.float Cab, 
    np.float Cm, np.float EK, np.float eCa, np.float ehp, np.float ECaN,
    np.float Esyn, np.float gK, np.float gL, np.float gNa, np.float Iapp, 
    np.float kIP3, np.float ks, np.float kNa, np.float kCa, 
    np.float kCAN, np.float Nab, 
    np.float siCAN, np.float sih, np.float sihp, np.float sim, 
    np.float simp, np.float siN, np.float sis, np.float tauh, 
    np.float tauhp, np.float taum, np.float taun,
    np.float taus, np.float Qh, np.float Qhp, np.float Qm,
    np.float Qmp, np.float Qn, np.float Qs, np.float ENa, np.float gsyn
    ):
    ## check types
    # assert y.dtype == np.float
    # assert vertex_types.dtype == long
    # assert edge_list.dtype == np.float
    # assert in_degrees.dtype == long
    # assert in_edges.dtype == long

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
    cdef np.ndarray dydt = \
        np.zeros(y.shape[0], dtype=DTYPE) # initialize vector field
    # constants
    cdef unsigned long num_eqns_per_vertex, num_vertices, num_edges
    num_vertices = vertex_types.shape[0]
    num_edges = edge_list.shape[0]
    num_eqns_per_vertex = 7 #V, Na m, Na h, K n, hp Nap, Ca Can, Na pump
    cdef unsigned long offset
    offset = num_vertices*num_eqns_per_vertex # for indexing edges
    # vertex variables
    cdef unsigned long i, j, k, type_idx
    for i in range( num_vertices ):
        j = i*num_eqns_per_vertex # index of V(i)
        type_idx = <unsigned int>(vertex_types[i]) # neuron type
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
        # first get in-edges
        these_in_edges = in_edges[i, 0:(in_degrees[i]-1)]
        # find the synaptic variables for those edges
        syn_variables = y[offset + these_in_edges]
        # and the conductances
        syn_conductances = edge_list[ these_in_edges, 2 ]
        # Ohm's law
        I_syn = np.dot(syn_variables, syn_conductances) * (y[j] - Esyn)
        ## set the derivatives
        # voltage V
        dydt[j] = - ( Iapp + 
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
        dydt[j+2] = ( infHN(Qh,sih,y[j]) - y[j+2] ) / tau(tauh, Qh, sih, y[j])
        # K n
        dydt[j+3] = ( infHN(Qn,siN,y[j]) - y[j+3] ) / tau(taun, Qn, siN, y[j])
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
        v_source = y[ this_edge[0] * num_eqns_per_vertex ]
        # v_target = y[ this_edge[1] * num_eqns_per_vertex ]
        j = offset + i
        dydt[j] = \
            ( (1 - y[j]) * \
                  infHN( Qs, sis, v_source ) - \
                  ks * y[j] ) / taus
    return dydt

# extra business
def _test():
    return 0

if __name__ == "__main__":
    _test()
