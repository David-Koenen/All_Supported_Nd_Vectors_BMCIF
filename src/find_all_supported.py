#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:40:04 2023

@author: David Könen

This module implements functions to find all supported nondominated vectors or efficient solutions
in a bi-objective minimum cost flow problem using  different methods. It includes functionality
for dichotomic search, all optimal flows, distinct values, epsilon, and the new epsilon method.

"""


import numpy as np
import pandas as pd
from all_best import all_optimal_flows
from all_distinct_values import get_all_distinct_per_face
from dichotomic_search import dichotomic_search


from epsilon import eps
from new_epsilon import new_eps
from timeout_decorator import timeout,TimeoutError
import networkx as nx 
from my_network import net_simplex



@timeout(500)
def allSupportedSolutions(HG):
        """
        Determine all Supported Efficient Solutions. Goes through one nondominated face. 
        """
        try:
            t, n = all_optimal_flows(HG)
            return t,n
        except:
            print('Fehler all Supported')
            return None,None
    
    
    
@timeout(500)
def allDistinct(HG2,optimal_flow,mst):
       """
       Determine all Supported NonDominated Vectors with All Distinct Method. Goes through on nondominated face. 
       """
       t, n = get_all_distinct_per_face(HG2,optimal_flow,mst)
       return t,n
    
@timeout(500)
def Epsilon(HG3,flow_cost):
        """
        Determine all Supported NonDominated Vectors with Epsilon Method on one nondominated face. 
        """
        t,n = eps(HG3,flow_cost)
        return t, n
    
    
@timeout(500)
def NewEpsilon(HG4,tree_arcs,optimal_flow,mst):
        """
        Determine all Supported NonDominated Vectors with New-Epsilon Method on one nondominated face. 
        """

        t,n = new_eps(HG4,tree_arcs,optimal_flow,mst)
        return t,n


@timeout(100)
def updatePiandRC(G,w, mst):
    """
    Update the Duals and reduced_costs as in Ahuja Chapter 11.4 (p.412)
    Returns: duals as pi, and reduced_cost as rc both as dictinoaries. 
    
    
    Parameters:
        G (nx.DiGraph): The directed graph.
        w (str): The weight attribute to use ('weight', 'new_cost', or 'new_cost2').
        mst (nx.Graph): Minimum spanning tree of the graph.

    Returns:
        tuple: Dictionaries of dual variables (pi) and reduced costs (rc).
    
    """
    
    
    ###########################################################################
    #  Determine potentials 
    ###########################################################################
    
    pred = nx.dfs_predecessors(mst,source=1)
    pre = list(nx.dfs_preorder_nodes(mst,source=1))
 
    pi = {n: {} for n in G.nodes}
    
    pi[1] = 0
    l = 1 
    j = pre[l]
    while l != len(G.nodes)-1:
        i = pred[j]
        if G.has_edge(i,j): 
            pi[j] = pi[i] - G[i][j][w]
        else: 
            pi[j] = pi[i] + G[j][i][w]
        l = l+1
        j = pre[l]
        
    i = pred[j]
    if G.has_edge(i,j): 
            pi[j] = pi[i] - G[i][j][w]
    else: 
            pi[j] = pi[i] + G[j][i][w]            

    
     ###########################################################################
     # Determine reduced costs. 
     ###########################################################################    


    rc = {}
    for u,v in G.edges():
        rc[(u,v)]= G[u][v][w] -pi[u] + pi[v]
    
    return pi, rc



def find_all_supported(name,G,costs,secondary_costs,m1,m2,m3,m4,save,save_name):
    r"""Find all supported nondominated vectors or efficient solutions in  
    a bi-objective minimum cost flow problem using different methods. 

    G is a digraph with edge costs and capacities and in which nodes
    have demand, i.e., they want to send or receive some amount of
    flow. A negative demand means that the node wants to send flow, a
    positive demand means that the node want to receive flow. A flow on
    the digraph G satisfies all demand if the net flow into each node
    is equal to the demand of that node.
    
    The Function prints all Numbers, Execution Time for each activated method 
    and save all results to a excel file if the variable save equals true. 
    

    Parameters
    ----------
    name: String
        Instance Name of the given bi-objective minimum cost-flow problem. 
    
    G : NetworkX graph
        DiGraph on which a minimum cost flow satisfying all demands is
        to be found.

    costs : Integer Vector for each arc.
        First objective of G.
    
    new_costs : Integer Vector for each arc.
        Second objective of G. 
        
    m1-m4 : Boolean
        Activate method if m_x equals True.
    
    save: Boolean
        Activate the Save-Function if save equals True.
        
    save-name: String 
        Name of the excel-file, where the results should be loaded to. 
    """
    
    if G != None:
        
     
        instance = name
        extreme_vectors, execution_time, lambdas  = dichotomic_search(G)
        
        
        print('Start:', instance)
        print('Start determining all supported vectors:')
        
        print('Step 1: Deterining all extrem supported vectors with dichotomic search:')
        
        print('Number of nondominated extrem points = ', len(extreme_vectors))
        print(f"Complete Execution time: {execution_time} seconds")
        
        print('')
        
        
        
        costs = np.array(costs)
        secondary_costs = np.array(secondary_costs)

         
        print('Step2: Determining all supported vectors using the four different Methods')
        print('')
        
        
        if m1 is True:
            execution_time_1 =  0
            number_1 = 0 
        else:
            execution_time_1 =  None
            number_1 = None
       
        if m2 is True:
            execution_time_2 =  0
            number_2 = 0 
        else:
            execution_time_2 =  None
            number_2 = None
        
        if m3 is True:
            execution_time_3 =  0
            number_3 = 0 
        else:
            execution_time_3 =  None
            number_3 = None
        
        if m4 is True:
            execution_time_4 =  0
            number_4 = 0 
        else:
            execution_time_4 =  None
            number_4 = None
        
        
        
        
        ###########################################################################
        # For each phase determine all supported nondominated vectors 
        ###########################################################################
        
        
                
        for i in range(0,len(lambdas)):
            new_costs =  np.array([])
            new_costs = lambdas[i][0]*costs + lambdas[i][1]*secondary_costs
                    #print(new_costs)
            for i,e in enumerate(G.edges()):
                G.edges[e]['new_cost'] = new_costs[i]
                G.edges[e]['new_cost2'] = 0.99*(new_costs[i])+0.01*(G.edges[e]['weight'])
                
                    
             ###########################################################################
             # Solves optimal solution with new_cost. Optimal solution for the given face
             ###########################################################################  
            try:
                optimal_flow, flowCost, tree_arcs, upper_cap, lower_caps, reduced_costs = net_simplex(G, 30, 'new_cost2')  
                if optimal_flow == None:
                    return 
            
                mst = nx.Graph()
                for e in tree_arcs:
                    mst.add_edge(*e)
                
                pi,reduced_costs = updatePiandRC(G, 'new_cost' ,  mst)
    
                # Determine reduced costs
                for  (u,v) in G.edges:
                    flow = optimal_flow[(u,v)]
                    cost = G.edges[(u,v)]['new_cost']
                    weight = G.edges[(u,v)]['weight']
                    lower_cap = 0
                    capacity = G.edges[(u,v)]['capacity']
                    
                    reduced_cost = reduced_costs[(u,v)]
                    G.add_edge(u, v, new_cost = cost, capacity = capacity, lower_cap = lower_cap, flow = flow, reduced_cost = reduced_cost,weight=weight)
                        
                        
                
                    
                ###########################################################################
                # Creates reduced graph, delete all arcs with reduced_cost unequal to zero 
                # Creates the solution space where each solution maps on a vector one the 
                # given face. 
                ###########################################################################  
                
                
                HG= G.copy()
                
                for u,v in G.edges():
                    if  reduced_costs[(u,v)] != 0 and (u,v) not in tree_arcs:
                        HG.nodes[u]['demand'] =  HG.nodes[u]['demand'] + HG.edges[(u,v)]['flow']
                        HG.nodes[v]['demand'] =   HG.nodes[v]['demand'] - HG.edges[(u,v)]['flow']
                        HG.remove_edge(u,v)
                        
                
                HG1 = HG.copy()
                HG2 = HG.copy()
                HG3 = HG.copy()
                HG4 = HG.copy()
               
               
            
                ###########################################################################
                # Method1: Determine All Supoorted Efficient Solutions
                ###########################################################################
                
                if m1 is True:
                    
                    
                    try:
                        t,n = allSupportedSolutions(HG1)
                        execution_time_1 = execution_time_1 + t
                        number_1 = number_1 + n 
                    except TimeoutError: 
                        print(name, "Function timed out: Method1") 
                        execution_time_1 = None 
                        number_1= None
                        m1 = False
                    except: 
                        print(name,'Fehler in Instance oder Code')
                        execution_time_1 = None 
                        number_1= None
                        m1 = False 
                        
                        
                        
                ###########################################################################
                # Method2: Determine All Supoorted NonDominated Vectors
                # All Distinct Method 
                ###########################################################################
                
                
                
                
                if m2 is True:
                    
                    
                    try:
                        t,n = allDistinct(HG2,optimal_flow, mst)
                        execution_time_2 = execution_time_2 + t
                        number_2 = number_2 + n 
                    except TimeoutError: 
                        print(name, "Function timed out: Method2") 
                        execution_time_2 = None 
                        number_2= None
                        m2 = False
                    except: 
                        print(name,'Fehler in Instance oder Code')
                        execution_time_2 = None 
                        number_2 = None
                        m2 = False 
                        
                
                
              
                
                ###########################################################################
                # Method3: Determine All Supoorted NonDominated Vectors
                # Epsilon Method
                ###########################################################################
                        
                if m3 is True:
                    
                    
                    try:
                        t,n = Epsilon(HG3,optimal_flow)
                        execution_time_3 = execution_time_3 + t
                        number_3 = number_3 + n 
                    except TimeoutError: 
                        print(name, "Function timed out: Method3") 
                        execution_time_3 = None 
                        number_3= None
                        m3 = False
                    except: 
                        print(name,'Fehler in Instance oder Code')
                        execution_time_3 = None 
                        number_3= None
                        m3 = False 
           
                # ###########################################################################
                # # Method4: Determine All Supoorted NonDominated Vectors
                # # New-Epsilon Method
                # ###########################################################################
        
                if m4 is True:
                    
                    
                    try:
                        t,n = NewEpsilon(HG3,tree_arcs,optimal_flow,mst)
                        execution_time_4 = execution_time_4 + t
                        number_4 = number_4 + n 
                    except TimeoutError: 
                        print(name, "Function timed out: Method4") 
                        execution_time_4 = None 
                        number_4 = None
                        m4 = False
                    except: 
                        print(name,'Fehler in Instance oder Code')
                        execution_time_4 = None 
                        number_4= None
                        m4 = False 
           
          
           
            except TimeoutError: 
                print(name, "Function timed out because of cycling")
          
            
            ###########################################################################
            # Print Number of Solutions and save them in the test file.
            ###########################################################################
            
            
        if number_1 is not None:
                    print("Method 1: All-Optimal-Method:") 
                    print('Number of  all supported solutions = ',number_1-len(lambdas)+1,f"Complete Execution time: {execution_time_1} seconds")
                    number1 =  number_1-len(lambdas)+1
        else: 
                number1 = None
            
            
            
        if number_2 is not None:
                    print("Method 2: Get distinct values:")  
                    print('Number of branches in get:all_distinct = ',number_2,f"Complete Execution time: {execution_time_2} seconds")
                    number2 =  number_2
        else: 
                    number2 = None
                
                
        if number_3 is not None:
                    
                    print("Method 3: Epsilon-Method:")  
                    print('Number of  supported vectors= ',number_3-len(lambdas)+1,f"Complete Execution time: {execution_time_3} seconds")
                    number3 = number_3-len(lambdas)+1
        else: 
                    number3 = None
                    
          
                
        if number_4 is not None:
                    print("Method 4: New-Epsilon-Method:")  
                    print('Number of  supported vectors= ',number_4-len(lambdas)+1,f"Complete Execution time: {execution_time_4} seconds")
                    number4 =  number_4-len(lambdas)+1
        else: 
                      number4 = None
        
        if save is True:
                    df = pd.read_csv(save_name)
                    new_row_series = pd.Series([instance, G.number_of_nodes(), G.number_of_edges(), len(extreme_vectors), execution_time, number3, number1, execution_time_1, number2, execution_time_2, execution_time_3, execution_time_4, number3, number4], index=['InstanceName', 'NumberNodes', 'NumberArcs', 'ExtremePoints', 'TimeExtremeVectors', 'SupportedPoints', 'NumberOptimalSolutions', 'TimeAllOptimal','NumberBranches2','TimeAllDistinct','TimeEpsilon','TimeNewEpsilon','SupEpsilon','SupNEpsilon'])
                                   
                    # Convert the pd.Series to a DataFrame
                    new_row_df = new_row_series.to_frame().T
            
                    # Append the new row to the original DataFrame
                    df = pd.concat([df, new_row_df], ignore_index=True)
                    file_path = save_name
                    df.to_csv(file_path, index=False)  # Set index=False to exclude the DataFrame index in the CSV
            
                    print(f'DataFrame saved to {file_path}')
                    print('')
                    print('')
