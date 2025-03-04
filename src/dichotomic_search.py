#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 09:09:52 2023

@author: david

This script implements a dichotomic search algorithm to find extreme nondominated vectors
in a multi-objective optimization problem using network flow methods.
"""



import time

from min_cost_solver import min_cost_solve_without_T

class dichotom:
    
    def __init__(self,G):
        """
        Initializes the Dichotom class.
       
        Parameters:
        G : networkx.Graph
           The graph on which the dichotomic search is performed.
       """
        self.extrem_vec = []
        
        
        
    def dichotomic(self,G):
        
        """
        Perform a dichotomic search to find extreme nondominated vectors.
        
        Parameters:
        G : networkx.Graph
            The graph with weighted edges and multiple objectives.
        
        Returns:
        tuple: (List of extreme nondominated vectors, execution time, lambda weights per face)
        """
        
        ###########################################################################
        # Starts dichotomic search with first finding the two  lexicographic best 
        # vectors corresponding to first and then second objective. 
        ###########################################################################    
        
        start_time = time.time()
        
        for e in G.edges():
            G.edges[e]['new_cost'] = (0.99 *  G.edges[e]['weight'] + 0.01 * G.edges[e]['secondary_cost'])


        #flowCost, flowDict, edge_flow, node_potentials = network_simplex(G,demand='demand',capacity='capacity',weight= 'new_cost')  
        optimal_flow, flow_Dict, flowCost = min_cost_solve_without_T(G, 0)       

            
        y1_1 = 0
        y1_2 = 0
   
        for i,e in enumerate(G.edges()):
            y1_1 += G.edges[e]['weight']*optimal_flow[(e[0],e[1])]
            y1_2 += G.edges[e]['secondary_cost']*optimal_flow[(e[0],e[1])]

        for e in G.edges():
            G.edges[e]['new_cost'] = (0.01 *  G.edges[e]['weight'] + 0.99 * G.edges[e]['secondary_cost'])


        optimal_flow, flow_Dict, flowCost = min_cost_solve_without_T(G, 0) 
        y2_1 = 0
        y2_2 = 0

        for i,e in enumerate(G.edges()):
            y2_1 += G.edges[e]['weight']*optimal_flow[(e[0],e[1])]
            y2_2 += G.edges[e]['secondary_cost']*optimal_flow[(e[0],e[1])]

        y1 = [y1_1, y1_2]
        y2 = [y2_1,y2_2]
        
        
        
        if y1 == y2:
            self.extrem_vec.append(y1)
            return self.extrem_vec, 1 , []

        #print(y1,y2)
        self.extrem_vec.append(y1)
        self.extrem_vec.append(y2)
        self.dichtomic_search(G, y1, y2)
        
        self.extrem_vec.sort()
        end_time = time.time()
        execution_time = end_time- start_time
        lambdas = self.determine_obj_per_face(self.extrem_vec)
        #print('Extrem nondominated vectors:', self.extrem_vec, 'Weights:', lambdas)
        return self.extrem_vec, execution_time, lambdas
    
    def dichtomic_search(self,G,ys,yr):
        
        """
        Recursively perform a dichotomic search to find intermediate extreme vectors.
        
        Parameters:
        G : networkx.Graph
            The graph with weighted edges and multiple objectives.
        ys : list
            One extreme vector.
        yr : list
            Another extreme vector.
        """
        
        ###########################################################################
        # For two extreme vectors, determine another extreme point inbetween 
        # if one exists and recusively does this for newly founded extreme points. 
        ###########################################################################    
        
        for e in G.edges():
            G.edges[e]['new_cost'] = ((ys[1]-yr[1]) *  G.edges[e]['weight'] + (yr[0]-ys[0]) * G.edges[e]['secondary_cost'])
        optimal_flow, flow_Dict, flowCost = min_cost_solve_without_T(G, 0)
      
        
        ynew_1 = 0
        ynew_2 = 0
        
        for i,e in enumerate(G.edges()):
            ynew_1 += G.edges[e]['weight']*optimal_flow[(e[0],e[1])]
            ynew_2 += G.edges[e]['secondary_cost']*optimal_flow[(e[0],e[1])]
        ynew= [ynew_1,ynew_2]
        #print((ys[1]-yr[1])*ys[0] +  (yr[0]-ys[0])*ys[1])
        #print(((ys[1]-yr[1])*ynew_1 +  (yr[0]-ys[0])*ynew_2))
        if  ((ys[1]-yr[1])*ys[0] +  (yr[0]-ys[0])*ys[1]) != ((ys[1]-yr[1])*ynew_1 +  (yr[0]-ys[0])*ynew_2):
            #print('hey')
            self.extrem_vec.append(ynew)
            self.dichtomic_search(G, ys, ynew)
            self.dichtomic_search(G, ynew, yr)
        
            
    def determine_obj_per_face(self,extreme_vectors):
        
        """
        Determine the weights for each consecutive extreme vector pair (supported nondominated faces).
        
        Parameters:
        extreme_vectors : list
            List of extreme nondominated vectors.
        
        Returns:
        list: Weights (lambda values) for each face.
        """
        
        ###########################################################################
        # Determine the weights for each consecutive extreme vectors, i.e, each 
        # supported nondominated face
        ###########################################################################  
        
        lambdas = []
        for i in range(0,len(extreme_vectors)-1):
            lambda_1 =  extreme_vectors[i][1]- extreme_vectors[i+1][1]
            lambda_2 = extreme_vectors[i+1][0] - extreme_vectors[i][0]
            lambdas.append([lambda_1,lambda_2])
        return lambdas  
    
    
    
def dichotomic_search(G):
   """
   Wrapper function to perform dichotomic search on a given graph.
   
   Parameters:
   G : networkx.Graph
       The graph with weighted edges and multiple objectives.
   
   Returns:
   tuple: (List of extreme nondominated vectors, execution time, lambda weights per face)
   """
   dich = dichotom(G)
   extreme_vectors, execution_time, lambdas  = dich.dichotomic(G)
   return extreme_vectors, execution_time, lambdas