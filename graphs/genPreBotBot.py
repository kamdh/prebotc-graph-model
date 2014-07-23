#!/usr/bin/env python

import respirnet
import numpy as np
import networkx as nx


n0 = 200
n1 = 100
pI = 0.2
gE = 2.5
gI = 2.5
pTypes = [0, 0.25, 0.45, 0.3]
a = 3.0
b = 0.5
c = 1.0
d = 3.0

pMatE = np.array([(a/(n0-1), b/(n1-1)),
                  (b/(n0-1), a/(n1-1))])
pMatI = np.array([(c/(n0-1), d/(n1-1)),
                  (d/(n0-1), c/(n1-1))])
g = respirnet.er_prebot_bot(n0, n1, pMatI, pMatE, pTypes, pI, gE, gI)
nx.write_gml(g, 'test_prebot_bot.gml')
