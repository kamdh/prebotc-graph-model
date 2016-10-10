#!/usr/bin/env python

import respirnet
import numpy as np
import networkx as nx
import sys
import argparse

'''
    This class outputs a file (.gml) that contains a network for the pre-bot and botzinger complex based on the desired network parameters. This class allows you to vary the architecture based on the amount of intra population inhibition. Recommended range is gamma between 0 and 1. Zero will translate to a model similar to the half center oscillator and one will be unifrom inhibition across the network
    
    INPUTS:
    n0 - The number of nodes that are in population 1 (pre-bot)
    n1 - The number of nodes that are in populatuon 2 (bot)
    intra - The amount of intra population inhibition
    inter - The amound of inter population inhibition
    #####OPTIONAL Inputs####
    
    pI - The probability of an inhibitory neuron being created
    gE - The conductance for excitatory snyapses in nS
    gI - The conductance for inhibitory synapses in nS
    
    Return:
    output - file name for the .gml file the graph will be saved to
    '''

def main(argv = None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(prog="genPreBotBot_gamma",
                                     description = 'Generates Graph based on Block Model with varying amounts of intra inhibition')
    parser.add_argument('n0', type = int, help='number of nodes in pop1')
    parser.add_argument('n1', type = int, help='number of nodes in pop2')
    parser.add_argument('intraDegree', type = float, help='the average degree for a population to itself')
    parser.add_argument('interDegree',type = float, help='the average degree for a population to the other')
    parser.add_argument('output', help='output filename')
    parser.add_argument('-pI', type = float,
                        help = 'probability of inhibitory neuron',
                        default = .2)
    parser.add_argument('-gE', type=float,
                        help='conductance of excitatory (E) synapses, nS',
                        default = 2.5)
    parser.add_argument('-gI', type=float,
                        help='conductance of inhibitory (I) synapses, nS',
                        default = 2.5)
        
    args = parser.parse_args(argv[1:])
    n0 = args.n0
    n1 = args.n1
    intraDegree = args.intraDegree
    interDegree = args.interDegree
    output = args.output
    gE = args.gE
    gI = args.gI
    pI = args.pI
                        
    #Each column of a matrix is responsible for the probability that index 1 connects to index 2
                        
    #Excitatory connection probability matrix
    pMatE = np.array([ (3.0/(n0-1), 0.00/n1),
                     (0.00/n0, 3.0/(n1-1)) ])
                        
                        
    #Inhibitory connection probaility matrix
    pMatI = np.array([ (intraDegree/(n0-1), interDegree/n1),
                        (interDegree/n0, intraDegree/(n1-1)) ])
                        
    #Probability distribution for neuron types (?,Bursting,Tonic,Quiescent)
    pTypes = [0, 0.25, 0.45, 0.3]
        
    g = respirnet.er_prebot_bot(n0,n1,pMatI, pMatE, pTypes, pI, gE, gI)
    nx.write_gml(g, output)

if __name__ == '__main__':
    status = main()
    sys.exit(status)

