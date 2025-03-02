#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:58:58 2023

@author: David Könen
"""


"""
This module implements the network simplex algorithm for solving minimum cost flow problems on directed connected graphs.
"""
import networkx as nx

class network_simplex:
    """
    Class implementing the network simplex algorithm.
    """
    
    def __init__(self,G):
        """
        Initializes the network simplex solver.

        Parameter:
            G (nx.DiGraph): The directed graph.
        """
        self.tree_arc = []  #list of tree_arcs
        self.lower_cap = [] # list of arcs on lower capacity
        self.upper_cap = [] # list of arcs in upper capacity
        self.edge_flow = {}  #dict of edge flow for each edge
        self.pi = {}  # dict of node potentials for each node
        self.rc = {} # dict of reduced_cost for each edge
        self.pre = [] # list ot preorder for current tree
        self.pred = {} # dictonary of predeccesor for each node of current tree
        self.T = nx.Graph() # Tree - Graph 
        self.entry_arc = None
        self.leaving_arc = None
        self.direct = None
        self.cycle = []
        self.lamda = None
     
        


    def initial_tree(self,G,HG,w):
        """
        Initializes the tree structure and flow based on the given graph and weight attribute.

        Parameter:
            G (nx.DiGraph): The original directed graph.
            HG (nx.DiGraph): A copy of the graph with potentially modified demands.
            w (str): The weight attribute to use ('weight', 'new_cost', or 'new_cost2').
        """
        if w == 'weight':
            for u in G.nodes:
                if u != 1:
                    b = HG.nodes[u]['demand'] 
                    if b >= 0:
                        if HG.has_edge(u,1) and HG.edges[(u,1)]['capacity'] > b:
                            self.tree_arc.append((u,1))
                            self.edge_flow[(u,1)] = b
                        elif HG.has_edge(u,1) or HG.has_edge(1,u):
                            number = HG.number_of_nodes() +1
                            HG.add_node(number,demand = 0)
                            HG.add_edge(u,number, weight = 0, capacity = 1000)
                            HG.add_edge(number,1, weight = 1000, capacity = 1000)
                            self.edge_flow[(u,number)] = b 
                            self.edge_flow[(number,1)] = b 
                            self.tree_arc.append((u,number))
                            self.tree_arc.append((number,1))
    
                        else:
                            HG.add_edge(u,1, weight = 1000, capacity = 1000)
                            self.edge_flow[(u,1)] = b 
                            self.tree_arc.append((u,1))
                    if b < 0:
                       if HG.has_edge(1,u) and HG.edges[(1,u)]['capacity'] > -b:
                           self.tree_arc.append((1,u))
                           self.edge_flow[(1,u)] = -b
                       elif HG.has_edge(1,u) or HG.has_edge(u,1):
                           number = HG.number_of_nodes() +1
                           HG.add_node(number,demand = 0)
                           HG.add_edge(1,number, weight = 0, capacity = 1000)
                           HG.add_edge(number,u, weight = 1000, capacity = 1000)
                           self.edge_flow[(1,number)] = -b 
                           self.edge_flow[(number,u)] = -b 
                           self.tree_arc.append((1,number))
                           self.tree_arc.append((number,u))    
                       else:
                           HG.add_edge(1,u, weight = 1000, capacity = 1000)
                           self.edge_flow[(1,u)] = -b 
                           self.tree_arc.append((1,u))
        
        elif w == 'new_cost':
                for u in G.nodes:
                    if u != 1:
                        b = HG.nodes[u]['demand'] 
                        if b >= 0:
                            if HG.has_edge(u,1) and HG.edges[(u,1)]['capacity'] > b:
                                self.tree_arc.append((u,1))
                                self.edge_flow[(u,1)] = b
                            elif HG.has_edge(u,1) or HG.has_edge(1,u):
                                number = HG.number_of_nodes() +1
                                HG.add_node(number,demand = 0)
                                HG.add_edge(u,number, new_cost = 0, capacity = 1000)
                                HG.add_edge(number,1, new_cost = 10000, capacity = 1000)
                                self.edge_flow[(u,number)] = b 
                                self.edge_flow[(number,1)] = b 
                                self.tree_arc.append((u,number))
                                self.tree_arc.append((number,1))
        
                            else:
                                HG.add_edge(u,1, new_cost = 10000, capacity = 1000)
                                self.edge_flow[(u,1)] = b 
                                self.tree_arc.append((u,1))
                        if b < 0:
                           if HG.has_edge(1,u) and HG.edges[(1,u)]['capacity'] > -b:
                               self.tree_arc.append((1,u))
                               self.edge_flow[(1,u)] = -b
                           elif HG.has_edge(1,u) or HG.has_edge(u,1):
                               number = HG.number_of_nodes() +1
                               HG.add_node(number,demand = 0)
                               HG.add_edge(1,number, new_cost = 0, capacity = 10000)
                               HG.add_edge(number,u, new_cost = 10000, capacity = 10000)
                               self.edge_flow[(1,number)] = -b 
                               self.edge_flow[(number,u)] = -b 
                               self.tree_arc.append((1,number))
                               self.tree_arc.append((number,u))    
                           else:
                               HG.add_edge(1,u, new_cost = 10000, capacity = 10000)
                               self.edge_flow[(1,u)] = -b 
                               self.tree_arc.append((1,u))
        
        elif w== 'new_cost2':
            for u in G.nodes:
                if u != 1:
                    b = HG.nodes[u]['demand'] 
                    if b >= 0:
                        if HG.has_edge(u,1) and HG.edges[(u,1)]['capacity'] > b:
                            self.tree_arc.append((u,1))
                            self.edge_flow[(u,1)] = b
                        elif HG.has_edge(u,1) or HG.has_edge(1,u):
                            number = HG.number_of_nodes() +1
                            HG.add_node(number,demand = 0)
                            HG.add_edge(u,number, new_cost2 = 0, capacity = 100000)
                            HG.add_edge(number,1, new_cost2 = 1000000, capacity = 100000)
                            self.edge_flow[(u,number)] = b 
                            self.edge_flow[(number,1)] = b 
                            self.tree_arc.append((u,number))
                            self.tree_arc.append((number,1))
    
                        else:
                            HG.add_edge(u,1, new_cost2 = 1000000, capacity = 100000)
                            self.edge_flow[(u,1)] = b 
                            self.tree_arc.append((u,1))
                    if b < 0:
                       if HG.has_edge(1,u) and HG.edges[(1,u)]['capacity'] > -b:
                           self.tree_arc.append((1,u))
                           self.edge_flow[(1,u)] = -b
                       elif HG.has_edge(1,u) or HG.has_edge(u,1):
                           number = HG.number_of_nodes() +1
                           HG.add_node(number,demand = 0)
                           HG.add_edge(1,number, new_cost2 = 0, capacity = 100000)
                           HG.add_edge(number,u, new_cost2 = 1000000, capacity = 100000)
                           self.edge_flow[(1,number)] = -b 
                           self.edge_flow[(number,u)] = -b 
                           self.tree_arc.append((1,number))
                           self.tree_arc.append((number,u))    
                       else:
                           HG.add_edge(1,u, new_cost2 = 1000000, capacity = 100000)
                           self.edge_flow[(1,u)] = -b 
                           self.tree_arc.append((1,u))
                        
        
        for e in HG.edges:
            if e not in self.tree_arc:
                self.lower_cap.append(e)
                self.edge_flow[e] = 0
        
        for n in HG.nodes():
            self.T.add_node(n)
        for e in self.tree_arc:
            self.T.add_edge(*e)
    
        self.pred = nx.dfs_predecessors(self.T,source=1)
        self.pre = list(nx.dfs_preorder_nodes(self.T,source=1))
     
        #print('pred:', self.pred)
        #print('pre', self.pre)
        self.updatePiandRC(HG,w)
        
        
    def updatePiandRC(self,HG,w):
        """
        Update the Duals and reduced_costs as in Ahuja Chapter 11.4 (p.412)
        Returns: duals as pi, and reduced_cost as rc both as dictinoaries. 
        """
        
        
        ###########################################################################
        #  Determine potentials 
        ###########################################################################
    
        
        
        self.pi[1] = 0
        l = 1 
        j = self.pre[l]
        while l != len(HG.nodes)-1:
            i = self.pred[j]
            if HG.has_edge(i,j): 
                self.pi[j] = self.pi[i] - HG[i][j][w]
            else: 
                self.pi[j] = self.pi[i] + HG[j][i][w]
            l = l+1
            j = self.pre[l]
            
        i = self.pred[j]
        if HG.has_edge(i,j): 
                self.pi[j] = self.pi[i] - HG[i][j][w]
        else: 
                self.pi[j] = self.pi[i] + HG[j][i][w]            

        
         ###########################################################################
         # Determine reduced costs. 
         ###########################################################################    



        for u,v in HG.edges():
            self.rc[(u,v)]= (HG[u][v][w] - self.pi[u] + self.pi[v])
        
        
    def identify_arc(self,HG,c):
        i = 0
       
        self.entry_arc = None
        min_cost = 0
        self.direct = 0
        
        if i <= c: 
            for e in self.lower_cap:
                if self.rc[e] < min_cost:
                    self.entry_arc = e
                    min_cost = self.rc[e]
                    self.direct = 1
                    i += 1
        
            for e in self.upper_cap:
                if -self.rc[e] < min_cost:
                    self.entry_arc = e
                    min_cost = -self.rc[e]
                    self.direct = -1 
                    i += 1
                    
                    
    def identify_cycle(self,HG):
     
        self.cycle = {}
        self.lamda = None
           
           ###########################################################################
           #  Search for leaving arc and change the arc 
           ###########################################################################
                
        if self.entry_arc != None:
            
            self.cycle[self.entry_arc] = 1
            self.leaving_arc = None
            
            ###########################################################################
            #  Creates undirectd tree to afterwards determine the lowest common ancestor.
            ###########################################################################
            
            
            T2 = nx.DiGraph()
            for n in self.T.nodes:
                T2.add_node(n)
            for e in self.T.edges:
                if self.pred[e[1]] == e[0]:
                    T2.add_edge(*e)
                else:
                    T2.add_edge(e[1],e[0])
           
            
            #Forward arc 
            if self.direct == 1:

                u = self.entry_arc[0] 
                
                v = self.entry_arc[1]
                
                min_capacity = HG[u][v]['capacity']
                self.leaving_arc = (u,v)
                a = nx.lowest_common_ancestor(T2, u, v)
                
                i = u     
                while i != a :
                    j = self.pred[i]
                    if HG.has_edge(j,i):
                        self.cycle[(j,i)] = 1 
                        #print('cycle1',self.cycle)
                        capacity = HG[j][i]['capacity'] -self.edge_flow[(j,i)]
                        if capacity < min_capacity:
                            min_capacity = capacity
                            self.leaving_arc = (j,i)
                    else:
                        self.cycle[(i,j)] = -1
                        #print('cycle2',self.cycle)
                        capacity = self.edge_flow[(i,j)]
                        if capacity < min_capacity:
                            min_capacity = capacity
                            self.leaving_arc = (i,j)
                    i = self.pred[i]
            
                
            
                i = v
                while i != a:
    
                    j = self.pred[i]
    
                    if HG.has_edge(i,j):
                          self.cycle[(i,j)] = 1
                          #print('cycle3',self.cycle)
                          capacity = HG[i][j]['capacity'] - self.edge_flow[(i,j)]
                          if capacity <= min_capacity:
                              min_capacity = capacity
                              self.leaving_arc = (i,j)
                    else:
                            self.cycle[(j,i)] = -1
                            #print('cycle4',self.cycle)
                            capacity = self.edge_flow[(j,i)]
                            if capacity <= min_capacity:
                                min_capacity = capacity
                                self.leaving_arc = (j,i)
                    i = self.pred[i]
            
            
            # Backward cycle         
            elif self.direct == -1:
                u = self.entry_arc[0] 
                v = self.entry_arc[1] 
                
                min_capacity = self.edge_flow[(u,v)]
                self.leaving_arc = (u,v)
                self.cycle[(u,v)] = - 1 
            
                a = nx.lowest_common_ancestor(T2, u, v)

                i = v 
                while i != a :
                       
                        j = self.pred[i]
                        if HG.has_edge(i,j):
                            self.cycle[(i,j)] = - 1 
                            #print('cycle5',self.cycle)
                            capacity = self.edge_flow[(i,j)]
                            if capacity < min_capacity:
                                min_capacity = capacity
                                self.leaving_arc = (i,j)
                        else:
                            self.cycle[(j,i)] = 1
                            #print('cycle6',self.cycle)
                            capacity = HG[j][i]['capacity'] - self.edge_flow[(j,i)]
                            if capacity < min_capacity:
                                min_capacity = capacity
                                self.leaving_arc = (j,i)
                        i =self.pred[i]
             
             
                i = u
                while i != a:
    
                    j = self.pred[i]

                    if HG.has_edge(i,j):
                        self.cycle[(i,j)] = 1
                        #print('cycle7',self.cycle)
                        capacity = HG[i][j]['capacity'] - self.edge_flow[(i,j)]
                        if capacity <= min_capacity:
                            min_capacity = capacity
                            self.leaving_arc = (i,j)
                          
                    else:
                      self.cycle[(j,i)] = -1 
                      #print('cycle8',self.cycle)
                      capacity = self.edge_flow[(j,i)]
                      if capacity <= min_capacity:
                          min_capacity = capacity
                          self.leaving_arc = (j,i)
                    i = self.pred[i]  
            
            
            self.lamda = min_capacity
    

    def augment_flow(self,HG):
        for (i,j) in self.cycle:
            self.edge_flow[(i,j)] += self.cycle[(i,j)]*self.lamda         
            
    def update(self,HG,w):
        self.tree_arc.append(self.entry_arc)
        self.tree_arc.remove(self.leaving_arc)
        self.T.add_edge(*self.entry_arc)
        self.T.remove_edge(*self.leaving_arc)
        
        
        if self.cycle[self.entry_arc] == 1:
            self.lower_cap.remove(self.entry_arc)
        else:
            self.upper_cap.remove(self.entry_arc)
        if self.edge_flow[self.leaving_arc] == 0:
            self.lower_cap.append(self.leaving_arc)
        else:
            self.upper_cap.append(self.leaving_arc)
        
        
        self.pred = nx.dfs_predecessors(self.T,source=1)
        self.pre = list(nx.dfs_preorder_nodes(self.T,source=1))
        
        self.updatePiandRC(HG,w)
        
        
    
    def network_simplex_run(self,G,c,w):
        HG = G.copy()
        for u in HG.nodes():
            HG.nodes[u]['demand'] = -   HG.nodes[u]['demand']
        self.initial_tree(G,HG,w)
        
        #print('Tree:', self.tree_arc, 'Lower', self.lower_cap, 'Upper', self.upper_cap, 'Flow', self.edge_flow)                
        #print('rc', self.rc)
        
        finish = 0
        while finish != 1:
           
            self.identify_arc(HG, c)
            if self.entry_arc != None:
                #print(self.entry_arc)
                self.identify_cycle(HG)
                #print('cycle',self.cycle)
                #print('lambda', self.lamda)        
                self.augment_flow(HG)
                #print('Flow',self.edge_flow)
                self.update(HG,w)
                #print(self.tree_arc)
            else: finish = 1
        
        
        optimal_flow = {}
        optimal_value = 0
        flow_dict = None
        rc2 = {}
        tree_arc2 = []
        lower_cap2 = []
        upper_cap2 = []
        
      
        
        for e in self.tree_arc:
            if e in G.edges:
                tree_arc2.append(e)
        
        for e in self.lower_cap:
            if e in G.edges:
                lower_cap2.append(e)
                
        for e in self.upper_cap:
            if e in G.edges:
                upper_cap2.append(e)
                
        
        for e in G.edges():
            optimal_flow[e] = self.edge_flow[e]
            optimal_value += self.edge_flow[e]*G.edges[e][w]
            rc2[e] = self.rc[e]          
        
        return optimal_flow, optimal_value, tree_arc2, upper_cap2, lower_cap2, rc2 
        
        
        
def net_simplex(G,c,w):
    """
   Solves the minimum cost flow problem using the network simplex algorithm.

   This function creates an instance of the `network_simplex` class,
   initializes it with the given graph `G`, and then runs the
   network simplex algorithm with the specified parameters `c` and `w`.

   Parameter:
       G (nx.DiGraph): The directed graph representing the network.
       c (int): A parameter controlling the number of arcs considered for entry into the tree.
                It influences the selection of the entering arc.
       w (str): The weight attribute to use for cost calculations.
                It should be one of 'weight', 'new_cost', or 'new_cost2'.

   Returns:
       tuple: A tuple containing the optimal flow, optimal value, tree arcs,
              upper capacity arcs, lower capacity arcs, and reduced costs.
              Specifically:
              - optimal_flow (dict): A dictionary representing the optimal flow for each edge.
              - optimal_value (float): The total cost of the optimal flow.
              - tree_arc2 (list): A list of arcs that are part of the optimal tree.
              - upper_cap2 (list): A list of arcs at their upper capacity in the optimal solution.
              - lower_cap2 (list): A list of arcs at their lower capacity in the optimal solution.
              - rc2 (dict): A dictionary of reduced costs for each edge in the graph.
   """
    solver = network_simplex(G)
    return solver.network_simplex_run(G,c,w)