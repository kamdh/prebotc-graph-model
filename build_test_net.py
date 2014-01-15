""" Generate disconnected testrespiratory network """

import graph_tool.all as gt
import pickle

f = open('param_files/param_dashevskiy.pkl')
p = pickle.load(f)
f.close()

# setup "enum" types
CS = 0
CI = 1
TS = 2
Sil = 3

gCaN_array = (p['gCaNS'], p['gCaNI'], p['gCaNTS'], p['gCaNSil']);
gP_array   = (p['gPS'], p['gPI'], p['gPTS'], p['gPSil']);
EL_array   = (p['ELS'], p['ELI'], p['ELTS'], p['ELSil']);

g = gt.Graph()
g.add_vertex(4)
edge_gsyn = g.new_edge_property("double")
vertex_type = g.new_vertex_property("int")
# vertex_gCaN = g.new_vertex_property("double")
# vertex_gP = g.new_vertex_property("double")
# vertex_EL = g.new_vertex_property("double")

# CS neuron
v = g.vertex(0)
vertex_type[v] = CS
# vertex_gCaN[v] = gCaN_array[CS]
# vertex_gP[v] = gP_array[CS]
# vertex_EL[v] = EL_array[CS]

# CI neuron
v = g.vertex(1)
vertex_type[v] = CI

# TS neuron
v = g.vertex(2)
vertex_type[v] = TS

# SIL neuron
v = g.vertex(3)
vertex_type[v] = Sil

g.edge_properties["gsyn"] = edge_gsyn
g.vertex_properties["type"] = vertex_type
g.save("../graphs/test.gml", fmt="gml")
