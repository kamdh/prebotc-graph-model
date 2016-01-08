#!/usr/bin/env python
import sys
import argparse
import networkx as nx
from respirnet import assign_gsyn

def main(argv=None):
    if argv is None:
        argv = sys.argv
    # parse arguments
    parser = argparse.ArgumentParser(prog="modify_g",
        description='modify gE and gI'
        )
    parser.add_argument('input', help='input filename')
    parser.add_argument('gE', type=float, help='excitatory conductance, nS')
    parser.add_argument('gI', type=float, help='inhibitory conductance, nS')
    parser.add_argument('output', help='output filename')
    args = parser.parse_args(argv[1:])
    gE = args.gE
    gI = args.gI
    graph=nx.read_gml(args.input)
    nx.set_edge_attributes(graph, 'gsyn',
                           {x: assign_gsyn(graph, x, gE, gI)
                            for x in graph.edges()})
    nx.write_gml(graph, args.output)

# run the main stuff
if __name__ == '__main__':
    status = main()
    sys.exit(status)
