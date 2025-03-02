#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 13:29:35 2023

@author: david
"""

import networkx as nx
import time
from min_cost_solver import min_cost_solve






class Directed_DFS:
    ''' Determine the dfs-tree and returns a zero-cycle if one exist and 
         an arc whehre the new flow differs from the old one'''
    
    def __init__(self,rG):
        # Graph as input
        # DFS intern variables
        self.time = 0
        self.visited = {node: False for node in rG.nodes()}
        self.start_time = {node: 0 for node in rG.nodes()}
        self.end_time = {node: 0 for node in rG.nodes()}
        self.parent =  {node: node for node in rG.nodes()}
        self.backward = []
        self.cross = []
        self.forward = []
        self.tree = []
        self.SBAlow = {node: 0 for node in rG.nodes()}
        self.cycle = []
        self.diffedge = {}

        
    def dfs(self, rG, node):
        self.visited[node] = True
        self.start_time[node] = self.time
        
        edge = (node, self.parent[node])
        if rG.has_edge(*edge):
            self.SBAlow[node] = self.SBAlow[self.parent[node]]
        else:
            self.SBAlow[node] = self.start_time[node]
        
        self.time += 1

                
        for neighbor in rG.neighbors(node):
            if not self.visited[neighbor]:
               self.parent[neighbor] = node
               edge = (node,neighbor)
               self.tree.append(edge)
               self.dfs(rG, neighbor)
            else:
               # when the parent node is traversed after the neighbour node
               if self.start_time[node] > self.start_time[neighbor] and self.end_time[neighbor] == 0:
                   edge = (node,neighbor)
                   self.backward.append(edge)
                   #print('Back edge')
               # when the neighbour node is a descendant but not a part of tree
               elif self.start_time[node] < self.end_time[neighbor]:
                   edge = (node,neighbor)
                   self.forward.append(edge)
                   #print('forward edge')
               elif self.start_time[node] > self.end_time[neighbor]:
                   edge = (node,neighbor)
                   self.cross.append(edge)
                   #print('cross edge')
        self.end_time[node] = self.time
        self.time += 1
        
        

     
    def directedDFS(self,rG):
        """ Perform depth-first search on the given graph """
        for node in rG.nodes():
            if not self.visited[node]:
                self.dfs( rG, node)

    def print_all_SBAlow_values(self):
        print(self.SBAlow)
        print(self.start_time)
        print(self.backward)
        print(self.tree)
        print(self.forward)
        print(self.parent)
        print(self.cross)
        
    def get_tree(self,rG):
        
        
        ###########################################################################
        # Determine tree
        ###########################################################################
        
        T = nx.DiGraph()
        for node in rG.nodes():
            T.add_node(node)
        for e in self.tree:
            T.add_edge(*e)
        return T
    
    def find_lca(self,T,node1,node2):
        
        
        ###########################################################################
        # Find Lowest Common Ancestor
        ###########################################################################
        
        return nx.lowest_common_ancestor(T, node1, node2)
        
    
    
    def Find_0_cycle(self,rG):
        
        ###########################################################################
        # Determine zero-cylce. 
        ###########################################################################
        self.directedDFS(rG)
        for u,v in self.backward:
            edge = (u,v)
            self.diffedge = edge 
            if v != self.parent[u]:
                edge = (u,v)
                self.cycle.append(edge)
                while v  != u:
                    edge = (self.parent[u],u)
                    self.cycle.append(edge)
                    u = self.parent[u]
                return self.diffedge, self.cycle
        for u,v in self.forward:
            edge = (u,v)
            self.diffedge = edge 
            if self.SBAlow[v] <= self.start_time[u]:
                edge = (u,v)
                self.cycle.append(edge)
                while v != u :
                    edge = (v,self.parent[v])
                    self.cycle.append(edge)
                    v = self.parent[v]
                return self.diffedge, self.cycle
        T = self.get_tree(rG)
        for u,v in self.cross:
            edge = (u,v)
            self.diffedge = edge 
            a = self.find_lca(T, u, v)
            if a:
                if self.SBAlow[v] <=  self.start_time[a]:
                    edge = (u,v)
                    self.cycle.append(edge)
                    while v != a:
                        edge=(v,self.parent[v])
                        self.cycle.append(edge)
                        v = self.parent[v]
                    while a != u :
                        edge = (self.parent[u],u)
                        self.cycle.append(edge)
                        u = self.parent[u]
                    return self.diffedge, self.cycle
        return self.diffedge,self.cycle    
    
    
    
    
class All_Flows:
    ''' Returns the number of all optimal flows and 
    the execution time for a given graph (MCFP)'''
    def __init__(self,G):
        self.number = 0
    
    def Get_another_optimal_flow(self,HG):
        
        
        rG = self.create_residual_graph(HG)
    
    
        t = Directed_DFS(rG)
        diffedge, cycle = t.Find_0_cycle(rG)
        
        if cycle: 
    
            if HG.has_edge(*diffedge):
                diffedge = diffedge
            else:
                diffedge = (diffedge[1],diffedge[0])
        
        
            newHG = HG.copy()
        
            for u,v in cycle:
                if HG.has_edge(u,v):
                    newHG.edges[(u,v)]['flow'] = newHG.edges[(u,v)]['flow'] +1 
                    
                else:
                    newHG.edges[(v,u)]['flow'] = newHG.edges[(v,u)]['flow'] - 1 
                 
            #print(nx.get_edge_attributes(newHG, "    G['costs'] = new_costsflow"))
            self.number += 1
                 
            if newHG.edges[diffedge]['flow'] > HG.edges[diffedge]['flow']:
                newHG.edges[diffedge]['lower_cap'] = HG.edges[diffedge]['flow'] + 1 
                self.Get_another_optimal_flow(newHG)
                HG.edges[diffedge]['capacity'] = HG.edges[diffedge]['flow'] 
                self.Get_another_optimal_flow(HG)
            else:
                newHG.edges[diffedge]['capacity'] = HG.edges[diffedge]['flow'] -1
                self.Get_another_optimal_flow(newHG)
                HG.edges[diffedge]['lower_cap'] = HG.edges[diffedge]['flow'] 
                self.Get_another_optimal_flow(HG)
                
                
                
                
    def solve_G(self,G):
        return min_cost_solve(G, 1) 


    def create_residual_graph(self,graph):
        
        
        ###########################################################################
        # Create residual graph 
        ###########################################################################
        
        residual_graph = nx.DiGraph()

        for u, v, attr in graph.edges(data=True):
            capacity = attr['capacity']
            lower_cap = attr['lower_cap']
            flow = attr.get('flow', 0)
            residual_capacity = capacity - flow
            if residual_capacity > 0:
                residual_graph.add_edge(u, v, capacity=residual_capacity, flow=0)

            if flow > lower_cap:
                residual_graph.add_edge(v, u, capacity=flow-lower_cap, flow=0)

        return residual_graph
    
    
    def Get_all_optimal_flows(self,HG):
        
            start_time = time.time()
        
            self.number += 1
           
                    
            rG = self.create_residual_graph(HG)

              ###########################################################################
              # Determine tree and find cycle
              ###########################################################################

            t = Directed_DFS(rG)
            diffedge, cycle = t.Find_0_cycle(rG)

            if cycle: 
                

                if HG.has_edge(*diffedge):
                    diffedge = diffedge
                else:
                    diffedge = (diffedge[1],diffedge[0])
                
                
                flow  =  nx.get_edge_attributes(HG, "flow")
                new_flow = nx.get_edge_attributes(HG, "flow")
                
                newHG = HG.copy()
                
                
                for u,v in cycle:
                    if HG.has_edge(u,v):
                        newHG.edges[(u,v)]['flow'] = newHG.edges[(u,v)]['flow'] +1 
                        
                    else:
                        newHG.edges[(v,u)]['flow'] = newHG.edges[(v,u)]['flow'] - 1 
                     
                #print(nx.get_edge_attributes(HG, "flow"))
                #print(nx.get_edge_attributes(newHG, "flow")) 
                self.number += 1
                
                ###########################################################################
                # Branch
                ###########################################################################
                     
                if newHG.edges[diffedge]['flow'] > HG.edges[diffedge]['flow']:
                    newHG.edges[diffedge]['lower_cap'] = HG.edges[diffedge]['flow'] + 1 
                    self.Get_another_optimal_flow(newHG)
                    HG.edges[diffedge]['capacity'] = HG.edges[diffedge]['flow'] 
                    self.Get_another_optimal_flow(HG)
                else:
                    newHG.edges[diffedge]['capacity'] = HG.edges[diffedge]['flow'] -1
                    self.Get_another_optimal_flow(newHG)
                    HG.edges[diffedge]['lower_cap'] = HG.edges[diffedge]['flow'] 
                    self.Get_another_optimal_flow(HG)
            end_time = time.time()
            execution_time = end_time - start_time
            #print('Number of optimal Flows = ', self.number,f"Execution time: {execution_time} seconds")
            
            return execution_time, self.number
   
def all_optimal_flows(G):
    all = All_Flows(G)
    execution_time, number = all.Get_all_optimal_flows(G)
    return execution_time, number








