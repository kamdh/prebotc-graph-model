#!/usr/bin/env python
import sys
import argparse
import numpy as np
import networkx as nx
from respirnet import er_prebot

## constants
ntypes = 4

def main(argv=None):
    if argv is None:
        argv = sys.argv
    # parse arguments
    parser = argparse.ArgumentParser(prog="genER",
        description='generate random graphs from G(n,p)'
        )
    parser.add_argument('n', type=int, help='number of nodes')
    parser.add_argument('p', type=float, help='probability of an edge')
    parser.add_argument('output', help='output filename')
    parser.add_argument('--deg', action='store_true',
                 help='interpret p as expected degree, i.e. p <- p/(n-1)')
    parser.add_argument('-pI', type=float, 
                        help='probability of an inhibitory neuron',
                        default=0.0)
    parser.add_argument('-gE', type=float,
                        help='conductance of excitatory (E) synapses, nS',
                        default = 2.5)
    parser.add_argument('-gI', type=float,
                        help='conductance of inhibitory (I) synapses, nS',
                        default = 2.5)
    parser.add_argument('--balanced', action='store_true',
                        help='set gI for a balanced network, given gE and pI')
    parser.add_argument('-pCS', type=float,
                        help='probability of CS neuron (default: %(default)s)',
                        default = 0.00)
    parser.add_argument('-pCI', type=float,
                        help='probability of CI neuron (default: %(default)s)',
                        default = 0.25)
    parser.add_argument('-pTS', type=float,
                        help='probability of TS neuron (default: %(default)s)',
                        default = 0.45)
    args = parser.parse_args(argv[1:])
    n = args.n
    if args.deg:
        p = 1.0 * args.p / (2*(n-1)) # compute p for kavg = args.p
    else:
        p = args.p
    pI = args.pI
    assert 0 <= pI and pI <= 1, 'pI out of range'
    gE = args.gE
    gI = args.gI
    if args.balanced and pI != 0:
        gI = (pI-1)*gE/gI
    assert args.pCS >= 0 and args.pCS <= 1, 'pCS out of range'
    assert args.pCI >= 0 and args.pCI <= 1, 'pCI out of range'
    assert args.pTS >= 0 and args.pTS <= 1, 'pTS out of range'
    sumPTypes = args.pCS + args.pCI + args.pTS
    assert sumPTypes <= 1, \
        'probabilities do not sum to <= 1'
    pTypes = [args.pCS, args.pCI, args.pTS, 1-sumPTypes]
    graph = er_prebot(n, p, pTypes, pI, gE, gI)
    nx.write_gml(graph, args.output)

# run the main stuff
if __name__ == '__main__':
    status = main()
    sys.exit(status)
