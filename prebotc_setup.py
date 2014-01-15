import pickle
import numpy as np
import graph_tool as gt

def params(paramFn):
    f = open(paramFn, 'r')
    my_params = pickle.load(f)
    f.close()
    return my_params

def graph(graphFn):
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
    
    return num_vertices, num_edges, vertex_types, edge_list, \
        in_edge_ct, in_edges

def ics(num_vertices, num_edges, num_eqns_per_vertex, num_eqns_per_edge):
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

