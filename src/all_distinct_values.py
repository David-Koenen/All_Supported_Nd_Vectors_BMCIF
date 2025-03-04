#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 14:48:08 2023

@author: david
"""

import networkx as nx
import time
from min_cost_solver import min_cost_solve




def updatePiandRC(G,w, mst):
   """
   Update the Duals and reduced_costs as in Ahuja Chapter 11.4 (p.412)
   Returns: duals as pi, and reduced_cost as rc both as dictinoaries. 
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


class all_distint_values:
    ''' Determine all distinct cost flows of a graph and 
    determines the number of branches needed. It retuns the number of branches and the execution time'''
    
    def __init__(self,rG):
        # Graph as input
        self.number =  0 
        self.min_edge = {}
        
    def solve_G(self,G):
        
        ###########################################################################
        # Solves first minimum cost flow problem 
        ###########################################################################

        
        return min_cost_solve(G, 1) 



    def create_residual_graph(self,graph):
        
        ###########################################################################
        # Create residual graph with reduced costs 
        ###########################################################################

        
        residual_graph = nx.DiGraph()

        for u, v, attr in graph.edges(data=True):
            capacity = attr['capacity']
            lower_cap = attr['lower_cap']
            flow = attr.get('flow', 0)
            reduced_cost = attr['reduced_cost']
            n_cost = attr['cost']
            residual_capacity = capacity - flow
            if residual_capacity > 0:
                residual_graph.add_edge(u, v, capacity=residual_capacity, cost= reduced_cost, ncost= n_cost)

            if flow > lower_cap:
                residual_graph.add_edge(v, u, capacity=flow-lower_cap, cost = -reduced_cost, ncost = -n_cost)

        return residual_graph
    
    
    
    def determine_minimal_cycle(self,rG,shortest_path,predecessors):
        
        
        ###########################################################################
        # Determine minimal cycle!
        # Warning min_value = 1000 could be too small if costs are too large! 
        ###########################################################################

        
        cycle = []
        self.min_edge = {}
        min_value = 1000
        for u,v in rG.edges():
            if rG.edges[(u,v)]['cost'] != 0:
                value = shortest_path[v][u] + rG.edges[(u,v)]['cost']
                if min_value > value:
                    min_value = value
                    self.min_edge = (u,v)
        if self.min_edge:
            path = nx.reconstruct_path(self.min_edge[1], self.min_edge[0], predecessors)
        
            length = len(path)
            for i in range (0,length-1):
                edge = (path[i], path[i+1])
                cycle.append(edge)
            edge = (self.min_edge[0], self.min_edge[1])
            cycle.append(edge)       
        return self.min_edge, cycle
    
    def update_node_potential(self,rG):
        
        ###########################################################################
        # Update Knotenpotentiale. 
        ###########################################################################
        
        edges = []
        for node1 in rG.nodes:
            for node2 in rG.nodes():
                if not rG.has_edge(node1, node2):
                    e = (node1,node2)
                    edges.append(e)
                    rG.add_edge(node1, node2, capacity = 1, cost= 1000, ncost = 1000)
        lenghts = nx.single_source_bellman_ford_path_length(rG,1,weight='ncost')
        # for e in edges:
        #     rG.remove_edge(*e)
        
        
        return lenghts
            
    
    def determine_second_best(self,G):
        
        
        ###########################################################################
        # Determine second best flow
        ###########################################################################
        
        rG = self.create_residual_graph(G)
        l = self.update_node_potential(rG)

        
        
        for u,v in G.edges:
            reduced_cost = G.edges[(u,v)]['cost']+l[u]-l[v]
            G.edges[(u,v)]['reduced_cost'] = reduced_cost
       
        rG  = self.create_residual_graph(G)
        
       
        predecessors, shortest_path = nx.floyd_warshall_predecessor_and_distance(rG, weight = 'cost')
        
        diffedge, cycle = self.determine_minimal_cycle(rG, shortest_path, predecessors)
        if cycle: 
            

            if G.has_edge(*diffedge):
                diffedge = diffedge
            else:
                diffedge = (diffedge[1],diffedge[0])
            
            
            #flow  =  nx.get_edge_attributes(HG, "flow")
            #new_flow = nx.get_edge_attributes(HG, "flow")
            
            newG = G.copy()
            
            
            for u,v in cycle:
                if G.has_edge(u,v):
                    newG.edges[(u,v)]['flow'] = newG.edges[(u,v)]['flow'] +1 
                    
                else:
                    newG.edges[(v,u)]['flow'] = newG.edges[(v,u)]['flow'] - 1 
                 
            #print(nx.get_edge_attributes(G, "flow"))
            cost = 0 
            for u,v in newG.edges():
                cost += newG.edges[(u,v)]['cost']*newG.edges[(u,v)]['flow']
            #print(nx.get_edge_attributes(newG, "flow"),'cost: ', cost) 
            self.number += 1
                 
            if newG.edges[diffedge]['flow'] > G.edges[diffedge]['flow']:
                newG.edges[diffedge]['lower_cap'] = G.edges[diffedge]['flow'] + 1 
                self.determine_second_best(newG)
                G.edges[diffedge]['capacity'] = G.edges[diffedge]['flow'] 
                self.determine_second_best(G)
            else:
                newG.edges[diffedge]['capacity'] = G.edges[diffedge]['flow'] -1
                self.determine_second_best(newG)
                G.edges[diffedge]['lower_cap'] = G.edges[diffedge]['flow'] 
                self.determine_second_best(G)
            
            
        
        
    def determine_all_distinct(self,G,optimal_flow,mst):
        
        start_time = time.time()
        self.number += 1
        # optimal_flow, flow_Dict, flowCost, tree_arcs, upper_cap, lower_caps, reduced_costs, duals = self.solve_G(G)

    
        
        pi,reduced_costs = updatePiandRC(G, 'weight' ,  mst)
        #optimal_flow, flowDict, flowCost, tree_arcs, upper_cap, lower_caps, reduced_costs, duals = min_cost_solve(HG, 1)
        
        # Determine reduced costs
        for  (u,v) in G.edges:
            flow = optimal_flow[(u,v)]
            cost = G.edges[(u,v)]['weight']
            weight = G.edges[(u,v)]['weight']
            lower_cap = 0
            capacity = G.edges[(u,v)]['capacity']
            
            reduced_cost = reduced_costs[(u,v)]
            G.add_edge(u, v, cost = cost, capacity = capacity, lower_cap = lower_cap, flow = flow, reduced_cost = reduced_cost,weight=weight)
      
                
        rG = self.create_residual_graph(G)
        predecessors, shortest_path = nx.floyd_warshall_predecessor_and_distance(rG, weight = 'cost')
        
        diffedge, cycle = self.determine_minimal_cycle(rG, shortest_path, predecessors)
        if cycle: 
            

            if G.has_edge(*diffedge):
                diffedge = diffedge
            else:
                diffedge = (diffedge[1],diffedge[0])
            
            
            newG = G.copy()
            
            
            for u,v in cycle:
                if G.has_edge(u,v):
                    newG.edges[(u,v)]['flow'] = newG.edges[(u,v)]['flow'] +1 
                    
                else:
                    newG.edges[(v,u)]['flow'] = newG.edges[(v,u)]['flow'] - 1 
                 
            
            cost = 0 
            for u,v in G.edges():
                cost += G.edges[(u,v)]['cost']*G.edges[(u,v)]['flow']
            #print(nx.get_edge_attributes(G, "flow"), 'cost: ', cost)
            cost = 0 
            for u,v in newG.edges():
                cost += newG.edges[(u,v)]['cost']*newG.edges[(u,v)]['flow']
            #print(nx.get_edge_attributes(newG, "flow"),'cost: ',    cost) 
            
            self.number += 1
                 
            if newG.edges[diffedge]['flow'] > G.edges[diffedge]['flow']:
                newG.edges[diffedge]['lower_cap'] = G.edges[diffedge]['flow'] + 1 
                self.determine_second_best(newG)
                G.edges[diffedge]['capacity'] = G.edges[diffedge]['flow'] 
                self.determine_second_best(G)
            else:
                newG.edges[diffedge]['capacity'] = G.edges[diffedge]['flow'] -1
                self.determine_second_best(newG)
                G.edges[diffedge]['lower_cap'] = G.edges[diffedge]['flow'] 
                self.determine_second_best(G)
        end_time = time.time()
        execution_time = end_time - start_time
        return execution_time, self.number
        

    
    
    
def get_all_distinct_per_face(HG,optimal_flow,mst):
    
    ###########################################################################
    # Determine supported vectors on a face. 
    ###########################################################################
    
    dist = all_distint_values(HG)
    execution_time, number = dist.determine_all_distinct(HG,optimal_flow,mst)  
    return execution_time,number

