import unittest
import networkx as nx
import os
import numpy as np
import sys
import math

'''
    This class allows for testing of graph generation using unit tests and saves them to a log file so it can be ran from the command line 
    efficiently
'''
    
class TestGraph(unittest.TestCase):
    
    '''
       Constructor that takes the desired properties that are going to be used for testing
       
       Parameters:
        testname-This is the name of the testmethod you are calling, this allows for efficient testing and adding using cmd line arguments
        graph-This is the graph name that you are looking at
        gamma-This is the amount of intra inhibitory
       
       NOTE: Must use the full path name for the grap or it will not be able to find it
    '''
    def __init__(self,testname,graph_file,gamma):
        super(TestGraph,self).__init__(testname)
        self.graph_file = graph_file
        self.gamma = 3 * float(gamma) / 2

    def setUp(self):
        pass

    '''
        This method tests the connection distribution of the graph generation methods. Since it is impossible to eliminate all the
        variance in the distribution due this test is used for graphs with large amounts of neurons in each population >= 10000
        this allows for the probability distribution to be much more accurate and allow for accurate comparisons
    '''
    def test_connectionDistribution(self):
        graph = nx.read_gml(self.graph_file)
        
        #Gets the probability distribution for connections
        n0 = len([d for n,d in graph.nodes_iter(data=True) if d['respir_area'] == 0])
        n1 = len([d for n,d in graph.nodes_iter(data=True) if d['respir_area'] == 1])
        
        area_size = [n0,n1]

        pMatE = np.array([ (3.0/(n0-1), .05/n1),
                          (.05/n0, 3.0/(n1-1)) ])
        pMatI = np.array([ (self.gamma/(n0-1), (3.0-self.gamma)/n1),
                          ((3.0-self.gamma)/n0, self.gamma/(n1-1)) ])
                          
        #First level list represents area number
        #Second level list represents type of neuron index 0 is exict and index 1 is inh

        excit_inh_node_count = [[0,0],[0,0]]

        #Gets totals for the types of neuron in each population
        for i in range(n0 + n1):
            if i < n0:
                if graph.node[i]['inh'] == 0:
                    excit_inh_node_count[0][0] += 1
                else:
                    excit_inh_node_count[0][1] += 1
            else:
                if graph.node[i]['inh'] == 0:
                    excit_inh_node_count[1][0] += 1
                else:
                    excit_inh_node_count[1][1] += 1
    
        excit_sampled = [[0,0],[0,0]]
        inh_sampled = [[0,0],[0,0]]

        #Calculates actual observed connectivity
        
        edges = graph.edges()

        for e in edges:
            source = e[0]
            target = e[1]
            source_index = 0
            target_index = 0
            
            if source < n0:
                source_index = 0
            else:
                source_index = 1
            if target < n0:
                target_index = 0
            else:
                target_index = 1
            
            if graph.node[source]['inh'] == 0:
                excit_sampled[source_index][target_index] += 1
            else:
                inh_sampled[source_index][target_index] += 1

        excit_prob_sampled = [[-1,-1],[-1,-1]]
        inh_prob_sampled = [[-1,-1],[-1,-1]]

        #Calculates sample probability for the given connections
        for i in range(2):
            for j in range(2):
                nodes_in_j = area_size[j]
                if i == j:
                    nodes_in_j -= 1
                if excit_inh_node_count[i][0] != 0:
                    excit_prob_sampled[i][j] = float(excit_sampled[i][j]) / (excit_inh_node_count[i][0] * nodes_in_j)
                if excit_inh_node_count[i][1] != 0:
                    inh_prob_sampled[i][j] = float(inh_sampled[i][j]) / (excit_inh_node_count[i][1] * nodes_in_j)


        #Calculates the percentage difference between two populations
        percentage_diff_excit = (abs(pMatE-excit_prob_sampled) / pMatE) * 100
        percentage_diff_inh = (abs(pMatI-inh_prob_sampled) / pMatI) * 100

        self.assertTrue(percentage_diff_excit[0][0] <= 5 and percentage_diff_excit[0][1] <= 5 and percentage_diff_excit[1][0] <= 5 and percentage_diff_excit[1][1] <= 5)
        self.assertTrue(percentage_diff_inh[0][0] <= 5 and percentage_diff_inh[0][1] <= 5 and percentage_diff_inh[1][0] <= 5 and percentage_diff_inh[1][1] <= 5)


if __name__ == '__main__':
    graph = sys.argv[1]
    gamma = sys.argv[2]
    output = sys.argv[3]
    log_file = output
    f = open(log_file,"w")
    suite = unittest.TestSuite()
    suite.addTest(TestGraph("test_connectionDistribution",graph,gamma))
    unittest.TextTestRunner(f,verbosity=2).run(suite)
    f.close()











