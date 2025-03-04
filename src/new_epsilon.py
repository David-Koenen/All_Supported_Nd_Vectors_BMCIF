#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 14:36:48 2023

@author: David Könen
"""
import networkx as nx
import time

from min_cost_solver import min_cost_solve 
import gurobipy as gp


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



class new_epsilon:
    def __init__(self,G):
        # Graph as input
        self.number =  0 
        self.arc_dict =   {}
        self.arc_direction={}
        
    def cycle(self, G, a ,T):
        
        ###########################################################################
        # Determines the cylce for the given nontree arc 
        ###########################################################################  
        
        if  G.edges[a]['flow']== 0:
            self.arc_direction[a] = 1
            
            path = nx.shortest_path(T, a[1], a[0])
     
            for i,n in enumerate(path[:-1]):
                n1 = path[i]
                n2 = path[i+1]
                if G.has_edge(n1,n2):
                    self.arc_dict[(n1,n2)].update({a: 1})
             
                else:
        
                    self.arc_dict[(n2,n1)].update({a: -1})
            #return cycle
        else: 
            self.arc_direction[a] = -1
        

            path = nx.shortest_path(T, a[1], a[0])
 
            for i,n in enumerate(path[:-1]):
                n1 = path[i]
                n2 = path[i+1]
                if G.has_edge(n1,n2):
                    self.arc_dict[(n1,n2)].update({a: -1})
                else:
                    self.arc_dict[(n2,n1)].update({a: 1})
                    
  
    
        
    def epsilion_scalarization(self,HG,tree_arcs,optimal_flow,mst):
    #try:   
        
        
        pi,reduced_costs = updatePiandRC(HG, 'weight' ,  mst)
        #optimal_flow, flowDict, flowCost, tree_arcs, upper_cap, lower_caps, reduced_costs, duals = min_cost_solve(HG, 1)
        
        # Determine reduced costs
        for  (u,v) in HG.edges:
            flow = optimal_flow[(u,v)]
            cost = HG.edges[(u,v)]['new_cost']
            weight = HG.edges[(u,v)]['weight']
            lower_cap = 0
            capacity = HG.edges[(u,v)]['capacity']
            
            reduced_cost = reduced_costs[(u,v)]
            HG.add_edge(u, v, new_cost = cost, capacity = capacity, lower_cap = lower_cap, flow = flow, reduced_cost = reduced_cost,weight=weight)
      
        
        ###########################################################################
        # Creates the cycles for each nontree arc
        ###########################################################################  

                
                
        non_tree_edges = []
        for e in HG.edges:
            if e not in tree_arcs:
                non_tree_edges.append(e)
        
        self.arc_dict = {a: {} for a in tree_arcs}
        self.arc_direction = {a: {} for a in non_tree_edges} 
         
        
        T = nx.Graph()
        for n in HG.nodes():
            T.add_node(n)
        for e in tree_arcs:
            T.add_edge(*e)

       
        for e in non_tree_edges:
          self.cycle(HG,e,T)
          
          
          
        ###########################################################################
        # Creates the new-epsilon model that should be solved in gurobi
        ###########################################################################  
  
         
         
        # Create a new Gurobi model
        model = gp.Model("Minimum_Cost_Flow")
        
        # Decision variables - flow on each edge (non-negative integers).
        lambdas = {}
        for (u, v) in non_tree_edges:
            lambdas[(u, v)] = model.addVar(vtype=gp.GRB.INTEGER, name=f"Lambda_{u}_{v}", lb=0, ub = HG[u][v]['capacity'])
        
        # Objective function - minimize the total cost
        obj1 = gp.quicksum(HG[u][v]['reduced_cost'] * lambdas[(u, v)] * self.arc_direction[(u, v)] for (u, v) in non_tree_edges)
        model.setObjective(obj1, gp.GRB.MINIMIZE)
        
        
        for (u,v) in tree_arcs:
            expr = HG[u][v]['flow'] + gp.quicksum(lambdas[(a, b)] * self.arc_dict[(u, v)][(a, b)] for a, b in self.arc_dict[(u, v)])
            model.addConstr(expr <= HG[u][v]['capacity'])
            model.addConstr(expr >= 0)
        
        
        # Solve the ILP.
        self.number += 1
    
    
        
        opt_val1 = 0
        
        
        # Adding the constraint to the Gurobi model
        constraint_target = model.addConstr(gp.quicksum(HG[u][v]['reduced_cost'] * lambdas[(u, v)] * self.arc_direction[(u, v)] for (u, v) in non_tree_edges) >= opt_val1 + 1, "Constraint_target")

        
       
        # Solve the MIP (Here as LP because of the intger property)
        model.setParam(gp.GRB.Param.Threads, 1)
        model.Params.OutputFlag = 0  # Suppress Gurobi output
        model.optimize()
        
        # Check if the optimization was successful
        if model.status == gp.GRB.OPTIMAL:
            model.remove(constraint_target)
            self.number += 1
            opt_val1 = model.ObjVal
    
        
        ###########################################################################
        # add the target constraint and solves the model until no new solution 
        # can be found. 
        ###########################################################################  
  
    
        while model.status == gp.GRB.OPTIMAL:
            
            constraint_target = model.addConstr(gp.quicksum(HG[u][v]['reduced_cost'] * lambdas[(u, v)] * self.arc_direction[(u, v)] for (u,v) in non_tree_edges) >= opt_val1 +1, "Constraint_target")
            model.setParam(gp.GRB.Param.Threads, 1)
            model.Params.OutputFlag = 0  # Suppress Gurobi output
            model.optimize()  
            if model.status == gp.GRB.OPTIMAL:
                model.remove(constraint_target)
                self.number += 1
                opt_val1 = model.ObjVal 
                
        model.dispose()
        return self.number
     #except: print('Fehler New-Epsilon') 
     #return 0
        
        
def new_eps(G,tree_arcs,optimal_flow,mst):
    
    """
    Uses the new-epsilon method to get all supported nondominated vectors 
    on one face. 
    
    Parameters
    ----------
    
    G : NetworkX graph
        DiGraph on which all supported nondominated vectors are to found 
        
    Returns
    -------
    number: the number of supported nondominated vectos found on the given face
    execution_time: time needed to determine these vectors. 
    """
    
    start_time = time.time() 
    epsilon2 = new_epsilon(G)
    number = epsilon2.epsilion_scalarization(G,tree_arcs,optimal_flow,mst)
    end_time= time.time()
    execution_time = end_time - start_time
    return execution_time,number
