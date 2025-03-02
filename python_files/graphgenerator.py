#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 09:14:55 2023

@author: David Könen
"""
import pynetgen as py 
import networkx as nx
import time
import random



def reformulate_graph(G):
    """
    Reformulates the graph to remove anti-parallel arcs. If an edge (u, v) exists
    along with (v, u), it introduces an intermediate node to resolve the conflict.
    Additionally, removes edges with zero capacity and ensures the graph remains connected.
    
    Parameters:
    G : networkx.DiGraph
        Input directed graph.
    
    Returns:
    networkx.DiGraph
        Reformulated graph.
    """
    HG = G.copy()
    k = 0 
    for u,v in HG.edges():
        if G.has_edge(u, v):
            if G.has_edge(v,u):
                k += 1
                number = G.number_of_nodes()+1
                G.add_node(number,demand = 0)
                G.add_edge(v,number, weight = 0, capacity = G[v][u]['capacity'], flow = 0, reduced_cost = 0, secondary_cost = 0 )
                G.add_edge(number,u, weight = G[v][u]['weight'], capacity = G[v][u]['capacity'], flow = 0, reduced_cost = 0, secondary_cost = G[v][u]['secondary_cost'])
                G.remove_edge(v,u)
    print(f"Reformulation: Added {k} nodes")
    
    ###########################################################################
    # If upper capacity equals lower capacity delete arc 
    ###########################################################################    
    
    
    l = 0 
    
    HG2 = G.copy()
    for u,v in HG2.edges():
        if G[u][v]['capacity'] == 0 :
            G.remove_edge(u, v)
            l = 1
    
    ###########################################################################
    # Tests if Graph is still connected
    ###########################################################################
    
    if l == 1:
        G2 = nx.Graph()
        for n in G.nodes():
            G2.add_node(n)
        for u,v in G.edges():
            G2.add_edge(u,v)
       
    
       
        
        if(len(list(nx.connected_components(G2))) == 1):
               print('Reformulation: G is still connected')
       
        else: 
              print(' After Reformulation: G is not connected anymore')
              raise ValueError("G is not connected")
            
    return G


def generate_random_numbers(n, min_val, max_val,seed):
    """
    Generate a list of random numbers within a given range.
    
    Parameters:
    n : int
        Number of random numbers to generate.
    min_val : int
        Minimum possible value.
    max_val : int
        Maximum possible value.
    seed : int
        Seed for reproducibility.
    
    Returns:
    list
        List of generated random numbers.
    """
    random.seed(seed)
    return [random.randint(min_val, max_val) for _ in range(n)]

def bicrit_test():
    """
    Constructs a test graph based on an example from the paper.
    Reads edge and node data from the file `instances/bicrit`.
    
    Returns:
    tuple
        A directed graph, list of edge costs, and list of secondary costs.
    """
    number_of_nodes= 5
    #py.netgen_generate( seed = 4, nodes = number_of_nodes, fname= 'instance_1')

    G = nx.DiGraph()
    for  i in range(1,number_of_nodes+1):
         G.add_node(i,demand= 0)


    demands = {node: 0 for node in G.nodes()}

    datei = open('instances/bicrit','r')
    #print(datei.read())

    costs = []

    for z in datei:
        if z[0] == 'n':
            nodes = z.split()
            G.nodes[int(nodes[1])]['demand'] = -int(nodes[2])
            demands[int(nodes[1])] = int(nodes[2])
        if z[0] == 'a':
            arcs = z.split()
            costs.append(int(arcs[5]))
            G.add_edge(int(arcs[1]), int(arcs[2]), weight = int(arcs[5]), capacity = int(arcs[4]), flow = 0, reduced_cost = 0, secondary_cost = 0 )
    #print(demands) 
    datei.close()

    secondary_costs = [5,1,5,9,7,2,4]
    for i,e in enumerate(G.edges()):
        G.edges[e]['secondary_cost'] = secondary_costs[i]
    return G,costs,secondary_costs
    
    
    
    
def netgen_random_mcfp(number_of_nodes, seed, sources, sinks, density, mincost, maxcost, supply, mincap, maxcap, fname,minsecondcost,maxsecondcost):
    
    """
    Creates an netgen Random Minimum Cost Flow Problem with a feasible solution using all
    the above parameters. 
    
    Returns:
    tuple
         DiGraph H, list of edge costs, and list of secondary costs.
    
    
    
    """
    
    try:
        py.netgen_generate( seed = seed, nodes = number_of_nodes, fname= "instances/"+fname, sources = sources, sinks = sinks, density = density, mincost= mincost, maxcost = maxcost, supply= supply, mincap=mincap, maxcap=maxcap)
    
        G = nx.DiGraph()
        for  i in range(1,number_of_nodes+1):
             G.add_node(i,demand= 0)
    
    
        demands = {node: 0 for node in G.nodes()}
    
        datei = open( "instances/"+fname,'r')
      
    
        costs = []
    
        for z in datei:
            if z[0] == 'n':
                nodes = z.split()
                G.nodes[int(nodes[1])]['demand'] = -int(nodes[2])
                demands[int(nodes[1])] = int(nodes[2])
            if z[0] == 'a':
                arcs = z.split()
                costs.append(int(arcs[5]))
                G.add_edge(int(arcs[1]), int(arcs[2]), weight = int(arcs[5]), capacity = int(arcs[4]), flow = 0, reduced_cost = 0, secondary_cost = 0 )
        #print(demands) 
        datei.close()
        
    
        
        secondary_costs = generate_random_numbers(len(G.edges), minsecondcost, maxsecondcost,seed)
        #print(secondary_costs)
        for i,e in enumerate(G.edges()):
            G.edges[e]['secondary_cost'] = secondary_costs[i]
            #print( G.edges[e]['secondary_cost'])
        H = reformulate_graph(G)
        c = nx.get_edge_attributes(G, "weight")
        s = nx.get_edge_attributes(G, "secondary_cost")
        secondary_costs = []
        costs = []
        for (u,v) in c:
            costs.append( c[(u,v)] )
            secondary_costs.append(s[(u,v)])
        return H,costs,secondary_costs
    except:  
        print('Instance Error during pygen Network Creation. Try another seed or other values') 
        return None,None,None
            
    
    
def netgen_random_mcfp_from_instance(seed,fname,number_of_nodes,number_of_arcs,minsecondcost,maxsecondcost):
        """
        Creates an DiGraph with all parameters and second cost function of the instance 
        from the file `tests/fname`. 
        Parameter seed gives the Seed used for the creation of the second cost function.
    
        Returns:
        tuple
            DiGraph H, list of edge costs, and list of secondary costs.
        """
        
        try:
            #py.netgen_generate( seed = seed, nodes = number_of_nodes, fname= "instances/"+fname, sources = sources, sinks = sinks, density = density, mincost= mincost, maxcost = maxcost, supply= supply, mincap=mincap, maxcap=maxcap)
        
            G = nx.DiGraph()
            for  i in range(1,number_of_nodes+1):
                 G.add_node(i,demand= 0)
        
        
            demands = {node: 0 for node in G.nodes()}
        
            datei = open( "instances/"+ number_of_nodes + "_"+ number_of_arcs + "/" +fname,'r')
          
        
            costs = []
        
            for z in datei:
                if z[0] == 'n':
                    nodes = z.split()
                    G.nodes[int(nodes[1])]['demand'] = -int(nodes[2])
                    demands[int(nodes[1])] = int(nodes[2])
                if z[0] == 'a':
                    arcs = z.split()
                    costs.append(int(arcs[5]))
                    G.add_edge(int(arcs[1]), int(arcs[2]), weight = int(arcs[5]), capacity = int(arcs[4]), flow = 0, reduced_cost = 0, secondary_cost = 0 )
            #print(demands) 
            datei.close()
            
        
            
            secondary_costs = generate_random_numbers(len(G.edges), minsecondcost, maxsecondcost,seed)
            #print(secondary_costs)
            for i,e in enumerate(G.edges()):
                G.edges[e]['secondary_cost'] = secondary_costs[i]
                #print( G.edges[e]['secondary_cost'])
            H = reformulate_graph(G)
            c = nx.get_edge_attributes(G, "weight")
            s = nx.get_edge_attributes(G, "secondary_cost")
            secondary_costs = []
            costs = []
            for (u,v) in c:
                costs.append( c[(u,v)] )
                secondary_costs.append(s[(u,v)])
            return H,costs,secondary_costs
        except:  
            print('Instance Error during pygen Network Creation. Try another seed or other values') 
            return None,None,None



def best_test_n(n,M,l):
    
    """
    Creates the best-case example from the paper. 
    
    Parameters
    ----------
    
    n = Number of Nodes
    M, l capacity of the given arcs. 
    
   Returns:
   tuple
       A directed graph, list of edge costs, and list of secondary costs.
    
    """
    
    G = nx.DiGraph()
    G.add_node(1,demand = -2)
    for i in range(2,n):
            G.add_node(i,demand = 0)
    G.add_node(n, demand = 2 )
    
    for i in range(1,n):
        G.add_edge(i, i+1, weight = 0, capacity = (n-3)*M+(l+2) , flow = 0, reduced_cost = 0, secondary_cost = 0 ) 
    
    for i in range(1,n-2):
        G.add_edge(n-i,n-i-2, weight = 0 , capacity = M, flow = 0 , reduced_cost = 0, secondary_cost = 0)
    
    G.add_edge(n,n-2, weight = 1, capacity = l, flow= 0, reduced_cost = 0, secondary_cost = -1 )
    
    
    secondary_costs = []
    costs =  []
    for i,e in enumerate(G.edges()):
        secondary_costs.append(G.edges[e]['secondary_cost'])
        costs.append(G.edges[e]['weight'])
    return G,costs,secondary_costs
    


def worst_test_n(n,l):
    
    """
    Creates the worst-case example from the paper. 
    
    Parameters
    ----------
    
    n = Number of Nodes
    l = capacity of the given arc. 
    
    Returns:
    tuple
        A directed graph, list of edge costs, and list of secondary costs.
    """
    
    
    G = nx.DiGraph()

    G.add_node(1,demand = -2)
    for i in range(2,n):
            G.add_node(i,demand = 0)
    G.add_node(n, demand = 2 )
    
    
    for i in range(1,n):
        G.add_edge(i, i+1, weight = 0, capacity = 2+l , flow = 0, reduced_cost = 0, secondary_cost = 0 ) 
    
    for i in range(2,n):
        G.add_edge(n,n-i, weight = 1 , capacity = l, flow = 0 , reduced_cost = 0, secondary_cost = -1)
    

    
    
    secondary_costs = []
    costs =  []
    for i,e in enumerate(G.edges()):
        secondary_costs.append(G.edges[e]['secondary_cost'])
        costs.append(G.edges[e]['weight'])
    return G,costs,secondary_costs




