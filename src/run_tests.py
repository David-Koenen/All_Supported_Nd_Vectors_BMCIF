#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 13:13:41 2023

@author: david
"""
from find_all_supported import find_all_supported
from graphgenerator import bicrit_test, best_test_n, worst_test_n
from graphgenerator import netgen_random_mcfp, netgen_random_mcfp_from_instance






def run_instance_from_file(nodes,arcs,number,m1,m2,m3,m4):
        """
        Runs an instance of a bi-objective minimum cost flow problem from a file given in the folder instances.
        Results will be saved 
        
        Note: Instnace must be in the folder. If you want to run instance_new_50_100_3 (50 nodes, 100 arcs, instance number 3)
        call run_instance_from_file(50,100,3,....)
        
        The other parameters are set, note that these can be easily modified. 
        
       Parameter:
          nodes (int): Number of nodes in the graph.
          arcs (int): Number of arcs in the graph.
          number (int): Instance number (used in file naming).
          m1 (bool): Flag to activate the first method (e.g., All-Optimal-Method).
          m2 (bool): Flag to activate the second method (e.g., All-Distinct-Method).
          m3 (bool): Flag to activate the third method (e.g., Epsilon-Method).
          m4 (bool): Flag to activate the fourth method (e.g., New-Epsilon-Method).
       """
        G = None
        name = 'instance_new'+ str(nodes) + '_'+str(arcs)+'_'+str(number)
        G, costs, secondary_costs = netgen_random_mcfp_from_instance(number, name, nodes, arcs, 1, 10) 
        if G != None:
            find_all_supported(name, G, costs, secondary_costs,  m1, m2, m3, m4, m4, '../tests/Tests_Results.csv')


def create_and_run_instance(nodes,arcs,seed,m1,m2,m3,m4):
    
        """
        Creates and runs an instance of a bi-objective minimum cost flow problem with the given nodes and arcs.
        The other parameter of the nwtwork generator as mincost, .... are fixed, but can be easily modified.  

      Parameter:
      nodes (int): Number of nodes in the graph.
      arcs (int): Number of arcs in the graph.
      seed (int): Random seed for instance generation.
      m1 (bool): Flag to activate the first method.
      m2 (bool): Flag to activate the second method.
      m3 (bool): Flag to activate the third method.
      m4 (bool): Flag to activate the fourth method.
      """    
    
        G = None
        name = 'instance_'+ str(nodes) + '_'+str(arcs)+'_'+str(seed)
        G, costs, secondary_costs = netgen_random_mcfp(nodes, seed, 2, 2, arcs, 1, 10, 50, 0, 50, name, 1, 10)
        if G != None:
            find_all_supported(name, G, costs, secondary_costs,  m1, m2, m3, m4, True, 'tests/Tests_Results.csv')
    



"""
Create best and worst case examples from the paper.

"""


def run_best_tests():
    """
   Runs the "best case" test instances, creating graphs with specific properties.
   """
    for i in range(5,10):
        G,costs,secondary_costs = best_test_n(i,10,5)
        name = 'best_' + str(i) + '_M=' + str(10) + '_l=' + str(5) 
        find_all_supported(name, G, costs, secondary_costs,  True, True, True, True, True, 'tests/test_best.csv')
    
    G, costs, secondary_costs =  best_test_n(20,1,3)
    name = 'best_' + str(20) + '_M=' + str(1) + '_l=' + str(3)
    find_all_supported(name, G, costs, secondary_costs,  True, True, True, True, True, 'tests/test_best.csv')
    #786432 




def run_worst_tests():
    """
    Runs the "worst case" test instances, creating graphs with specific properties.
    """
    for i in range(5,11):
        G,costs,secondary_costs = worst_test_n(i,5)
        name = 'worst_' + str(i) + '_l=' + str(5)
        find_all_supported(name, G, costs, secondary_costs,  True, True, True, True, True, 'tests/test_worst.csv')

    
    G, costs, secondary_costs = worst_test_n(20,5)
    name =  'worst_' + str(20) + '_l=' + str(3)
    find_all_supported(name, G, costs, secondary_costs,  True, True, True, True, True, 'tests/test_worst.csv')
    #33649



       
                
