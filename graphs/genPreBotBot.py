#!/usr/bin/env python

import respirnet
import numpy as np
import networkx as nx


n = 300
pMatE = np.array([ (2.95/(n-1), 0.05/(n-1)), 
                   (0.05/(n-1), 2.95/(n-1)) ])
pMatI = np.array([ (1.0/(n-1), 2.0/(n-1)), 
                   (2.0/(n-1), 1.0/(n-1)) ])
pI = 0.2
gE = 2.5
gI = 2.5
pTypes = [0, 0.25, 0.45, 0.3]

g = respirnet.er_prebot_bot(n, pMatI, pMatE, pTypes, pI, gE, gI)
nx.write_gml(g, 'test5.gml')
