#!/usr/bin/env python
import sys
import numpy as np
import networkx as nx


G = nx.DiGraph()
G.add_node(0, { 'type':1, 'inh':0 })
G.add_node(1, { 'type':3, 'inh':0 })
G.add_node(2, { 'type':1, 'inh':1 })
G.add_node(3, { 'type':3, 'inh':0 })
G.add_edge(0,1, {'gsyn':2.5})
G.add_edge(2,3, {'gsyn':-2.5})
nx.write_gml(G, 'test_syn.gml')
