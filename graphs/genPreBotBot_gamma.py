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
        gamma - The amount of intra population inhibition
        
        #####OPTIONAL Inputs####
        
        pI - The probability of an inhibitory neuron being created 
        gE - The conductance for excitatory snyapses in nS
        gI - The conductance for inhibitory synapses in nS
        
    Return:
        output - file name for the .gml file the graph will be saved to
        
    TO-DO:
        VERIFY INPUTS MATCH DESIRED (CHECK IF type= x takes care of it)
'''

def main(argv = None):
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(prog="genPreBotBot_gamma",
            description = 'Generates Graph based on Block Model with varying amounts of intra inhibition')
    parser.add_argument('n0', type = int, help='number of nodes in pop1')
    parser.add_argument('n1', type = int, help='number of nodes in pop2')
    parser.add_argument('gamma', type = float, help='the percentage of unifrom inhbiiton')
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
    #Multiply by a factor of 3/2 to maintain an average degree of 6 for each neuron
    gamma = 3*float(args.gamma) / 2
    output = args.output
    gE = args.gE
    gI = args.gI
    pI = args.pI
    
    #Each column of a matrix is responsible for the probability that index 1 connects to index 2

    #Excitatory connection probability matrix
    pMatE = np.array([ (3.0/(n0-1), 0.05/n1),
                       (0.05/n0, 3.0/(n1-1)) ])


    #Inhibitory connection probaility matrix
    pMatI = np.array([ (gamma/(n0-1), (3.0-gamma)/n1),
                       ((3.0-gamma)/n0, gamma/(n1-1)) ])

    #Probability distribution for neuron types (?,Bursting,Tonic,Quiescent)
    pTypes = [0, 0.25, 0.45, 0.3]
    
    print pMatE
    print pMatI
		
    g = respirnet.er_prebot_bot(n0,n1,pMatI, pMatE, pTypes, pI, gE, gI)
    nx.write_gml(g, output + '.gml')

if __name__ == '__main__':
    status = main()
    sys.exit(status)
