#!/usr/bin/env python
import pickle

f = open('param_files/param_dashevskiy.pkl', 'w')
my_params = dict(
    #gP = 4.0,
    EL = -0.0625,
    # CS params
    gCaNS = 10.0,
    gPS = 0.5,
    ELS = -0.06,
    # CI params
    gCaNI = 0.0,
    gPI = 4.0,
    ELI = -0.0625,
    # TS params
    gCaNTS = 0.0,
    gPTS = 5.,
    ELTS = -0.062,
    # Sil / quiescent params
    gCaNSil = 0.0,
    gPSil = 0.6,
    ELSil = -0.0605,
    # shared params
    alpha = 6.6e-2,
    Cab = 0.05,
    Cm = 0.045, # nF
    EK = -0.075,
    eCa = 0.0007, # .0007, epsilon in Rubin
    ehp = 0.001,
    ECaN = 0.0,
    Esyn = 0.0,
    gK = 30.0,
    gL = 3.0,  # 2-4, units: nS
    gNa = 160.0,
    Iapp = 0.0,
    kIP3 = 1.2e+6, # microM / s
    ks = 1.0,
    kNa = 10.0,
    kCa = 22.5e+3,
    kCAN = 0.9,
    Nab = 5.0,
    siCAN = -0.05,
    sih = 0.005,
    sihp = 0.006,
    sim = -0.0085,
    simp = -0.006,
    siN = -0.005,
    sis = -0.003,
    tauh = 0.015,
    tauhp = 0.001,
    taum = 0.001,
    taun = 0.030,
    taus = 0.015,
    Qh = -0.030,
    Qhp = -0.048,
    Qm = -0.036,
    Qmp = -0.040,
    Qn = -0.030,
    Qs = 0.015,
    ENa = 0.045,
    gsyn = 2.5
    )
pickle.dump(my_params, f)
f.close()
