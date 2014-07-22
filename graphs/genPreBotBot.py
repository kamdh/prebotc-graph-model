#!/usr/bin/env python

import respirnet
import numpy as np
import networkx as nx


n = 200
pI = 0.2
gE = 2.5
gI = 2.5
pTypes = [0, 0.25, 0.45, 0.3]
a = 2.0
b = 1.0
c = 2.0

pMatE = np.array([(a, b),
                  (b, a)]) / (n-1)
pMatI = np.array([(a, c),
                  (c, a)]) / (n-1)
g = respirnet.er_prebot_bot(n, pMatI, pMatE, pTypes, pI, gE, gI)
nx.write_gml(g, 'test_prebot_bot.gml')
