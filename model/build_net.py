""" Generate Tatiana Dashevskiy's respiratory network """

import graph_tool.all as gt

# load parameters
gP = 4.0
EL = -0.0625
gCaNS = 10.0
gPS = 0.5
ELS = -0.06
gCaNI = 0.0
gPI = 4.0
ELI = -0.0625
gCaNTS = 0.0
gPTS = 5.
ELTS = -0.062
gCaNSil = 0.0
gPSil = 0.6
ELSil = -0.0605
alfa = 6.6e-2
Cab = 0.05
Cm = 0.045
EK = -0.075
eCa = 0.0007
ehp = 0.001
ECaN = 0.0
Esyn = 0.0
gK = 30.0
gL = 3.0
gNa = 160.0
Iapp = 0.0
kIP3 = 1.2e+6
ks = 1.0
kNa = 10.0
kCa = 22.5e+3
kCAN = 0.9
Nab = 5.0
rp = 0.2
siCAN = -0.05
sih = 0.005
sihp = 0.006
sim = -0.0085
simp = -0.006
siN = -0.005
sis = -0.003
tauh = 0.015
tauhp = 0.001
taum = 0.001
taun = 0.030
taus = 0.015
Qh = -0.030
Qhp = -0.048
Qm = -0.036
Qmp = -0.040
Qn = -0.030
Qs = 0.015
ENa = 0.045
gsyn = 2.5

# setup "enum" types
CS = 0
CI = 1
TS = 2
Sil = 3

gCaN_array = (gCaNS, gCaNI, gCaNTS, gCaNSil);
gP_array   = (gPS, gPI, gPTS, gPSil);
EL_array   = (ELS, ELI, ELTS, ELSil);


neurTypes = (CS, CI, TS, TS, CI, Sil, Sil,
             CS, Sil, TS, CS, CI, TS, TS, 
             CI, Sil, Sil, CS, Sil, TS, CS,
             CI, TS, TS, CI, Sil, Sil, CS, Sil, TS)

gsyn1 = (gsyn, gsyn, -15.0, -15.0, gsyn, gsyn, 
         gsyn, gsyn, 0.0, gsyn, gsyn, gsyn, -15.0, 
         -15.0, gsyn, gsyn, gsyn, gsyn, 0.0, gsyn, 
         gsyn, gsyn, -15.0, -15.0, 0.0, gsyn,
         gsyn, gsyn, 0.0, gsyn)

gsyn2 = (0.0, gsyn, 0.0, 0.0, 0.0, gsyn,
         gsyn, 0.0, gsyn, gsyn, 0.0, 0.0,
         0.0, 0.0, 0.0, gsyn, gsyn, 0.0, gsyn, 
         0.0, 0.0, gsyn, 0.0, 0.0, 0.0, 
         gsyn, gsyn, 0.0, gsyn, gsyn)

syn1 = ( 7, 0, 0, 2, 25, 4, 4, 0, 8, 7, 17, 10, 10, 12, 14, 14, 14, 
         10, 18, 17, 27, 20, 20, 22, 24, 24, 24, 20, 28, 27 )
syn2 = ( 0, 10, 2, 3, 4, 3, 7, 7, 9, 16, 10, 11, 12, 13, 5, 13, 17, 
         17, 19, 19, 20, 0, 22, 23, 24, 23, 27, 27, 29, 6 )

g = gt.Graph()
g.add_vertex(30)
edge_gsyn = g.new_edge_property("double")
vertex_type = g.new_vertex_property("int")
vertex_gCaN = g.new_vertex_property("double")
vertex_gP = g.new_vertex_property("double")
vertex_EL = g.new_vertex_property("double")

for vidx in range(30):
    v = g.vertex(vidx)
    vertex_type[v] = neurTypes[vidx]
    type_idx = neurTypes[vidx]
    vertex_gCaN[v] = gCaN_array[type_idx]
    vertex_gP[v] = gP_array[type_idx]
    vertex_EL[v] = EL_array[type_idx]
    s = g.vertex( syn1[vidx] )
    e = g.add_edge( s, v )
    edge_gsyn[e] = gsyn1[vidx]
    s = g.vertex( syn2[vidx] )
    e = g.add_edge( s, v )
    edge_gsyn[e] = gsyn2[vidx]

g.edge_properties["gsyn"] = edge_gsyn
g.vertex_properties["type"] = vertex_type
g.save("Dashevskiy.gml", fmt="gml")
