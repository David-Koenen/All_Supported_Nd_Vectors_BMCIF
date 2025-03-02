#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 10:07:26 2023

@author: David Könen

This module implements functions to solve minimum cost flow problems using Gurobi.
It includes functionality to find optimal tree structures and update dual variables.
"""

import gurobipy as gp
import networkx as nx
from timeout_decorator import timeout,TimeoutError





    

def change_tree(G, tree_arcs, upper_cap, lower_cap, rc, mst,optimal_flow):
    """
    Due to Gurobipy issues, we might not get an optimal tree structure (T,L,U)
    Change_tree will obtain such an optimal tree structure!
    !!WARNING!! Cycle?
    
    Finds min_arc <-- (entering arc) and the change_arc <-- (leaving_arc)
    
    Parameters:
      G (nx.DiGraph): The original directed graph.
      tree_arcs (list): List of arcs in the current tree structure.
      upper_cap (list): List of arcs at their upper capacity.
      lower_cap (list): List of arcs at their lower capacity.
      rc (dict): Dictionary of reduced costs for each arc.
      mst (nx.Graph): Minimum spanning tree of the graph.
      optimal_flow (dict): Dictionary of optimal flow values for each arc.

  Returns:
      tuple: Updated tree_arcs, upper_cap, lower_cap, mst, and a flag indicating if a change occurred.

    
    """
    
    min_arc = None
    min_cost = 0
    direct = 0
    set_change = 0 
    
    pred = nx.dfs_predecessors(mst,source=1)
    
    ###########################################################################
    #  Creates undirectd tree to afterwards determine the lowest common ancestor.
    ###########################################################################
    
    
    T = nx.DiGraph()
    for n in mst.nodes:
        T.add_node(n)
    for e in mst.edges:
        if pred[e[1]] == e[0]:
            T.add_edge(*e)
        else:
            T.add_edge(e[1],e[0])

    
    ###########################################################################
    #  Determines entering arc and its direction. 
    ###########################################################################

    for e in lower_cap:
        if rc[e] < min_cost:
            min_arc = e
            min_cost = rc[e]
            direct = 1
    
    for e in upper_cap:
        if -rc[e] < min_cost:
            min_arc = e
            min_cost = -rc[e]
            direct = -1 
            
            
       ###########################################################################
       #  Search for leaving arc and change the arc 
       ###########################################################################
            
    if min_cost != 0 :
        set_change = 1
        change_arc = None
        pred = nx.dfs_predecessors(mst,source=1)
        #Forward arc 
        if direct == 1:

            u = min_arc[0] 

            v = min_arc[1]

            a = nx.lowest_common_ancestor(T, u, v)

            i = v
            while i != a:

                j = pred[i]

                if G.has_edge(i,j):
                        if optimal_flow[(i,j)] == G[i][j]['capacity']:
                            change_arc = (i,j)
                            upper_cap.append(change_arc)
                            break
                else:
                        if optimal_flow[(j,i)] == 0:
                            change_arc = (j,i)
                            lower_cap.append(change_arc)
                            break 
                i = pred[i]
            i = u     
            while i != a :
                j = pred[i]
                if G.has_edge(j,i):
                            if optimal_flow[(j,i)] == G[j][i]['capacity']:
                                change_arc = (j,i)
                                upper_cap.append(change_arc)
                                break
                else:
                            if optimal_flow[(i,j)] == 0:
                                change_arc = (i,j)
                                lower_cap.append(change_arc)
                                break 
                i = pred[i]
            
        # Backward cycle         
        elif direct == -1:
            u = min_arc[0] 
            v = min_arc[1] 
        
            a = nx.lowest_common_ancestor(T, u, v)
    
            i = u
            while i != a:
                

                j = pred[i]

                if G.has_edge(i,j):
                       
                        if optimal_flow[(i,j)] == G[i][j]['capacity']:
                            change_arc = (i,j)
                            upper_cap.append(change_arc)
                            break
                else:
                      
                        if optimal_flow[(j,i)] == 0:
                            change_arc = (j,i)
                            lower_cap.append(change_arc)
                            break 
                i = pred[i]  
            i = v 
            while i != a :
                   
                    j = pred[i]
                    if G.has_edge(j,i):
                               
                                if optimal_flow[(j,i)] == G[j][i]['capacity']:
                                    change_arc = (j,i)
                                    upper_cap.append(change_arc)
                                    break
                    else:
                                
                                if optimal_flow[(i,j)] == 0:
                                    change_arc = (i,j)
                                    lower_cap.append(change_arc)
                                    break 
                    i = pred[i]
        
        
        ###########################################################################
        #  Change the old structure (T,L,U) to the new one (T*,L*,U*)
        ###########################################################################

        if change_arc != None:
            
            mst.remove_edge(*change_arc)
            mst.add_edge(*min_arc)
            tree_arcs.remove(change_arc)
            tree_arcs.append(min_arc)
            if direct == 1:
                lower_cap.remove(min_arc)
            else: 
                upper_cap.remove(min_arc)
            
    
        
    return  tree_arcs, upper_cap, lower_cap, mst, set_change 




def updatePiandRC(G,w, mst):
    """
    Updates dual variables (pi) and reduced costs (rc) as in Ahuja Chapter 11.4 (p.412).

    Parameter:
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
            


@timeout(4000)
def min_cost_solve(G, weight):
    
    
    """
    Uses a gurobi-python version to get an optimal tree soluiton
    It also determines an optimal TreeStructure (T,L,U)
    
    Parameters
    ----------
    
    G : NetworkX graph
        DiGraph on which a minimum cost flow satisfying all demands is
        to be found.
        
    weight: If weight = 1, G[u][v]['weight'] is used if weight = 0 G[u][v]['newcost']
            is used. 
        
    Returns
    -------
    Tree-Struture: tree_arcs, upper_cap, lower_cap 
    Reduced_costs for each arc and duals variables for each node 
        
    
    """
    
    
    ###########################################################################
    # Create complete Gurobi Model. 
    ###########################################################################    
    
    
    # Create a new Gurobi model
    model = gp.Model("Minimum_Cost_Flow")
    
    
    # Checking for which objective should be searched. 
    if weight == 1:
        w = 'weight'
    elif weight == 0:
        w = 'new_cost'
    else:
        w = 'new_cost2'
        
    
    # Decision variables - flow on each edge (non-negative integers).
    flow = {}
    for (u, v) in G.edges():
        flow[(u, v)] = model.addVar(name=f"Flow_{u}_{v}", lb=0, ub = G[u][v]['capacity'])
    
    # Objective function - minimize the total cost.
    obj1 = gp.LinExpr(
    [G[u][v][w] for (u, v) in G.edges()],
    [flow[(u, v)] for (u, v) in G.edges()]
    )
    
    
    model.setObjective(obj1, gp.GRB.MINIMIZE)
    
    # Constraints - flow conservation constraint, capacity constraints already in upper bound. 
    for u in G.nodes():
        list_of_edges = []
        for e in G.edges():
                     if e[1] == u :
                         list_of_edges.append(e)
                
        flow_out = gp.quicksum(flow[(u,v)] for v in G[u])
        flow_in = gp.quicksum(flow[e] for e in list_of_edges)
        demand = G.nodes[u]['demand']
        model.addConstr( flow_in - flow_out == demand, name=f"FlowConservation_{u}")
    
    
    # Solve the MIP (Here as LP because of the intger property)
    model.Params.OutputFlag = 0  # Suppress Gurobi output
    model.optimize()
    
    
    ###########################################################################
    # Get optimal_flow, flow_dict, duals, and reduced costs
    ###########################################################################    
    
    # Get the optimal flow values
    optimal_flow = {(u, v): flow[(u, v)].x for (u, v) in G.edges()}
    
    # Get the optimal flow value
    optimal_value = model.objVal
    
    #Create the flow_dict with f[u][v] you get the complete flow. 
    flow_dict = {}
    for n in G.nodes():
        flow_dict[n]= {v : flow[u,v].x for (u,v) in G.edges()}
    
    
    # print(optimal_flow)
        
    # Extract dual variables (Pi) for flow conservation constraints
    duals = {i+1: constr.Pi for i,constr in enumerate(model.getConstrs()) if "FlowConservation" in constr.ConstrName}
    
     # Extract reduced costs (RC) for decision variables
    reduced_costs = {(u,v): flow[(u,v)].RC for (u,v) in G.edges()}

    
 
    
    # Extract vbasis and cbasis
    vbasis = model.getAttr('VBasis', model.getVars())
    #cbasis = model.getAttr('CBasis', model.getConstrs())
    
    
    variables = model.getVars()
    

    basic_variables = []
    
    list_of_edges = list(G.edges)
    for i, var in enumerate(variables):
         if vbasis[i] == gp.GRB.BASIC:
             basic_variables.append(list_of_edges[i])
    
        
    ###########################################################################
    # If we do not have an optimal tree structure due to 
    # many slack variables we use function change_tree to get such an optimal
    # tree stucture. 
    ###################################################################min_cost_solve########   
    try:
        if len(basic_variables) != len(G.nodes())-1: 
            T = nx.Graph()
            for node in G.nodes():
                T.add_node(n)
            for arc in G.edges():
                if arc in basic_variables:
                    T.add_edge(*arc, weight = 0 )
                elif reduced_costs[arc] == 0: 
                    T.add_edge(*arc, weight = 5 )
                else:
                    T.add_edge(*arc, weight = 20 )
            
            mst = nx.minimum_spanning_tree(T)
        
        
            tree_arcs = []
            non_tree_arcs = []
            lower_cap = []
            upper_cap = []
            
            
            for e in mst.edges():
                if G.has_edge(*e):
                    tree_arcs.append(e)
                else:
                    tree_arcs.append((e[1],e[0]))
            
            for i,e in enumerate(G.edges()):
                if e not in tree_arcs:
                    non_tree_arcs.append(e)
                    if vbasis[i] == -1:
                        lower_cap.append(e)
                    else:
                        upper_cap.append(e)
        
                
          
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
        
            
            
         
            rc = {}
            for u,v in G.edges():
                rc[(u,v)]= G[u][v][w] -pi[u] + pi[v]
            
        
            set_change = 1
            
        
            while set_change != 0:
                tree_arcs, upper_cap, lower_cap, mst, set_change = change_tree(G, tree_arcs, upper_cap, lower_cap, rc, mst, optimal_flow) 
                pi,rc = updatePiandRC(G,w, mst)
                
                
            model.dispose()
        
            return optimal_flow, flow_dict,  optimal_value , tree_arcs, upper_cap, lower_cap, rc, pi 

    ###########################################################################
    # If we get already an optimal tree structure we use the given reduced_costs 
    # and duals. 
    ####### WARNING !!!! ######
    # Are the reduced_costs and duals, the ones we want?! Tests this. 
    # If problems determine node potentials and reduced_costs and return them. 
    ###########################################################################   
    
        else:
            tree_arcs = basic_variables
            upper_cap = []
            lower_cap = []
            for i,e in enumerate(G.edges()):
                if e not in basic_variables:
                    if vbasis[i] == -1:
                        lower_cap.append(e)
                    else:
                        upper_cap.append(e)
            model.dispose()
    
            return optimal_flow, flow_dict,  optimal_value , tree_arcs, upper_cap, lower_cap , reduced_costs, duals 
  
    except TimeoutError: 
        print("Function timed out, Warning Cycling may occured here!")
        return None,None, None, None, None, None, None, None 


def min_cost_solve_without_T(G, weight):
    
    
    """
    Uses a gurobi-python version to get an optimal soluiton

    
    Parameters
    ----------
    
    G : NetworkX graph
        DiGraph on which a minimum cost flow satisfying all demands is
        to be found.
        
    weight: If weight = 1, G[u][v]['weight'] is used if weight = 0 G[u][v]['newcost']
            is used. 
        
    Returns
    -------
    optimal_flow
        
    
    """
    
    
    ###########################################################################
    # Create complete Gurobi Model. 
    ###########################################################################    
    
    
    # Create a new Gurobi model
    model = gp.Model("Minimum_Cost_Flow")
    
    
    # Checking for which objective should be searched. 
    if weight == 1:
        w = 'weight'
    else:
        w = 'new_cost'
        
    
    # Decision variables - flow on each edge (non-negative integers).
    flow = {}
    for (u, v) in G.edges():
        flow[(u, v)] = model.addVar(name=f"Flow_{u}_{v}", lb=0, ub = G[u][v]['capacity'])
    
    # Objective function - minimize the total cost.
    obj1 = gp.LinExpr(
    [G[u][v][w] for (u, v) in G.edges()],
    [flow[(u, v)] for (u, v) in G.edges()]
    )
    
    
    model.setObjective(obj1, gp.GRB.MINIMIZE)
    
    # Constraints - flow conservation constraint, capacity constraints already in upper bound. 
    for u in G.nodes():
        list_of_edges = []
        for e in G.edges():
                     if e[1] == u :
                         list_of_edges.append(e)
                
        flow_out = gp.quicksum(flow[(u,v)] for v in G[u])
        flow_in = gp.quicksum(flow[e] for e in list_of_edges)
        demand = G.nodes[u]['demand']
        model.addConstr( flow_in - flow_out == demand, name=f"FlowConservation_{u}")
    
    
    # Solve the MIP (Here as LP because of the intger property)
    model.Params.OutputFlag = 0  # Suppress Gurobi output
    model.optimize()
    
    
    ###########################################################################
    # Get optimal_flow, flow_dict, duals, and reduced costs
    ###########################################################################    
    
    # Get the optimal flow values
    optimal_flow = {(u, v): flow[(u, v)].x for (u, v) in G.edges()}
    
    # Get the optimal flow value
    optimal_value = model.objVal
    
    #Create the flow_dict with f[u][v] you get the complete flow. 
    flow_dict = {}
    for n in G.nodes():
        flow_dict[n]= {v : flow[u,v].x for (u,v) in G.edges()}
    
    
    return optimal_flow, flow_dict,  optimal_value




