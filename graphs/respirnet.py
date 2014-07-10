import networkx as nx
import numpy as np
import random
import itertools

def assign_type(pTypes):
    numTypes = len(pTypes)
    # node is type 0, 1, 2, 3 wp given in pTypes
    return np.random.choice(range(numTypes), p=pTypes)

def assign_inh(pI):
    # node is inhibitory wp pI
    if np.random.rand() < pI:
        inh = 1
    else:
        inh = 0
    return inh
        
def assign_gsyn(graph, e, gE, gI):
    # assign excitatory or inhibitory synapses based on whether projecting neuron
    # is inhibitory or excitatory
    if nx.get_node_attributes(graph, 'inh')[ e[0] ] == 1:
        # inhibitory
        gsyn = gI
    else:
        # excitatory
        gsyn = gE
    return gsyn

def er_prebot(n, p, pTypes, pI, gE, gI):
    assert isinstance(n, int), 'n should be integer'
    assert 0 <= p and p <= 1, 'p out of range'
    assert 0 <= pI and pI <= 1, 'pI out of range'
    assert sum(pTypes) == 1, 'pTypes must sum to 1'
    # generate graph
    graph = nx.fast_gnp_random_graph(n, p, directed=True)
    # assign node types (uniformly at random)
    nx.set_node_attributes(graph, 'type', 
                           {x: assign_type(pTypes) for x in graph.nodes()})
    nx.set_node_attributes(graph, 'inh',
                           {x: assign_inh(pI) for x in graph.nodes()})
    nx.set_edge_attributes(graph, 'gsyn', 
                           {x: assign_gsyn(graph, x, gE, gI)
                            for x in graph.edges()})
    return graph

def er_prebot_bot(n, p_mat_I, p_mat_E, pTypes, pI, gE, gI):
    '''
    Generate a model for preBot and Bot coupled respiratory groups.
    This is a signed (E-I) version of a 2-group stochastic block model.

    Parameters
    ==========
    n: int
        Number of nodes in each complex, preBot and Bot
    p_mat_I: 2x2 ndarray
        Probability of inhibitory projections
    p_mat_E: 2x2 ndarray
        Probability of excitatory projections
    pTypes: float array
        Probability of each node type 0, 1, ..., numtypes
    pI: float
        Probability a node is inhibitory
    gE: float
        Max conductance of excitatory connections
    gI: float
        Max conductance of inhibitory connections
    '''
    assert isinstance(n, int), 'n should be integer'
    assert p_mat_I.ndim == 2 and \
        np.all( np.array(p_mat_I.shape) == 2 ), 'p_mat_I wrong shape'
    assert np.all( p_mat_I <= 1 ) and np.all( p_mat_I >= 0 ), \
        'p_mat_I out of range'
    assert p_mat_E.ndim == 2 and \
        np.all( np.array(p_mat_E.shape) == 2 ), 'p_mat_E wrong shape'
    assert np.all( p_mat_E <= 1 ) and np.all( p_mat_E >= 0 ), \
        'p_mat_E out of range'
    assert 0 <= pI and pI <= 1, 'pI_inter out of range'
    assert sum(pTypes) == 1, 'pTypes must sum to 1'
    assert isinstance(gE, float), 'gE should be float'
    assert isinstance(gI, float), 'gI should be float'

    # setup nodes
    graph = nx.empty_graph(2*n)
    graph = nx.DiGraph(graph)
    nx.set_node_attributes(graph, 'respir_area', {x: 0 if x < n else 1
                                                  for x in graph.nodes()})
    nx.set_node_attributes(graph, 'type', 
                           {x: assign_type(pTypes) for x in graph.nodes()})
    nx.set_node_attributes(graph, 'inh',
                           {x: assign_inh(pI) for x in graph.nodes()})
    # now add edges
    edges = itertools.permutations(range(2*n),2)
    for e in edges:
        # draw edges with differing probabilities depending on types of endpoints
        u = graph.node[ e[0] ]
        v = graph.node[ e[1] ]
        area1 = u['respir_area']
        area2 = v['respir_area']
        if u['inh'] == 1:
            # inhibitory projection
            p = p_mat_I[ area1, area2 ]
        elif u['inh'] == 0:
            # excitatory
            p = p_mat_E[ area1, area2 ]
        if random.random() < p:
            graph.add_edge(*e)
    # update conductances... node types should be inherited from before
    nx.set_edge_attributes(graph, 'gsyn', 
                           {x: assign_gsyn(graph, x, gE, gI)
                            for x in graph.edges()})
    graph.name = "er_prebot_bot"
    return graph
