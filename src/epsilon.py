#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 13:00:59 2023

@author: David Könen
"""

import networkx as nx
import time

from min_cost_solver import min_cost_solve 
import gurobipy as gp



class epsilon:
    def __init__(self,G):
        # Graph as input
        self.number =  0 
        
        
        
    def epsilion_scalarization(self,HG,optimal_flow):
    
        
        
        
       # Create a new Gurobi model
       model = gp.Model("Minimum_Cost_Flow")
        
    
            
        
       # Decision variables - flow on each edge (non-negative integers).
       flow = {}
       for (u, v) in HG.edges():
            flow[(u, v)] = model.addVar(vtype=gp.GRB.INTEGER, name=f"Flow_{u}_{v}", lb=0, ub = HG[u][v]['capacity'])
        
       # Objective function - minimize the total cost.
       obj1 = gp.LinExpr(
        [HG[u][v]['weight'] for (u, v) in HG.edges()],
        [flow[(u, v)] for (u, v) in HG.edges()]
        )
        
        
       model.setObjective(obj1, gp.GRB.MINIMIZE)
        
       # Constraints - flow conservation constraint, capacity constraints already in upper bound. 
       for u in HG.nodes():
            list_of_edges = []
            for e in HG.edges():
                         if e[1] == u :
                             list_of_edges.append(e)
                    
            flow_out = gp.quicksum(flow[(u,v)] for v in HG[u])
            flow_in = gp.quicksum(flow[e] for e in list_of_edges)
            demand = HG.nodes[u]['demand']
            model.addConstr( flow_in - flow_out == demand, name=f"FlowConservation_{u}")

        
       # Solve the MIP (Here as LP because of the intger property)
       model.setParam(gp.GRB.Param.Threads, 1)
       model.Params.OutputFlag = 0  # Suppress Gurobi output
       #model.optimize()
       
       # Check if the optimization was successful
       #if model.status == gp.GRB.OPTIMAL:
       #   self.number += 1
       #   opt_val1 = model.ObjVal
       #   print('opt',opt_val1)
       #self.number += 1
       opt_val1 = 0
       self.number+= 1
       for e in HG.edges:
           opt_val1 += optimal_flow[e] * HG.edges[e]['weight']
  
       
       constraint_target = model.addConstr(gp.quicksum(HG[u][v]['weight'] * flow[(u, v)] for (u,v) in HG.edges()) >= opt_val1 + 1, "Constraint_target")
       model.setParam(gp.GRB.Param.Threads, 1)
       model.Params.OutputFlag = 0  # Suppress Gurobi output
       model.optimize()  
       if model.status == gp.GRB.OPTIMAL:
          model.remove(constraint_target)
          self.number += 1
          opt_val1 = model.ObjVal 
         ###########################################################################
         # add the target constraint and solves the model until no new solution 
         # can be found. 
         ###########################################################################  
        
       # Constraints - flow conservation and capacity constraints.
       while model.status == gp.GRB.OPTIMAL:
            constraint_target = model.addConstr(gp.quicksum(HG[u][v]['weight'] * flow[(u, v)] for (u,v) in HG.edges()) >= opt_val1 + 1, "Constraint_target")
            model.setParam(gp.GRB.Param.Threads, 1)
            model.Params.OutputFlag = 0  # Suppress Gurobi output
            model.optimize()  
            if model.status == gp.GRB.OPTIMAL:
                model.remove(constraint_target)
                self.number += 1
                opt_val1 = model.ObjVal 
    
       model.dispose() 
       return self.number
        
def eps(G,optimal_flow):
    
    """
    Uses the Epsilon method to get all supported nondominated vectors 
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
    epsilon2 = epsilon(G)
    number = epsilon2.epsilion_scalarization(G,optimal_flow)
    end_time= time.time()
    execution_time = end_time - start_time
    return execution_time,number
    