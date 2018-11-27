#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math

from ortools.constraint_solver import pywrapcp
from ortools.constraint_solver import routing_enums_pb2
from timeit import default_timer as timer
import pandas as pd

def euclid_distance(x1, y1, x2, y2):
    # Euclidean distance between points.
    dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    return dist

def create_distance_matrix(locations):
# Create the distance matrix.
    size = len(locations)
    dist_matrix = {}
    for from_node in range(size):
        dist_matrix[from_node] = {}
    
    for from_node in range(size):
        for to_node in range(from_node,size):
            x1 = locations[from_node][0]
            y1 = locations[from_node][1]
            x2 = locations[to_node][0]
            y2 = locations[to_node][1]
            d = euclid_distance(x1, y1, x2, y2)
            if d < 400:
                dist_matrix[from_node][to_node] = d
    print("Created distance matrix")
    return dist_matrix

def create_distance_callback(dist_matrix):
  # Create the distance callback.

    def distance_callback(from_node, to_node):
        if to_node < from_node:
            from_node, to_node = to_node, from_node
            
        if to_node in dist_matrix[from_node]:
            return dist_matrix[to_node][from_node]
        else:
            return 10000

    return distance_callback

def main():
    # Create the data.
    locations = create_data_array()
    dist_matrix = create_distance_matrix(locations)
    dist_callback = create_distance_callback(dist_matrix)
    tsp_size = len(locations)
    num_routes = 1
    depot = 0

    # Create routing model.
    if tsp_size > 0:
        routing = pywrapcp.RoutingModel(tsp_size, num_routes, depot)
        search_parameters = pywrapcp.RoutingModel.DefaultSearchParameters()
        search_parameters.local_search_metaheuristic = (
            routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH)
        search_parameters.time_limit_ms = 600*1000
        search_parameters.first_solution_strategy = (
                                    routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)
        routing.SetArcCostEvaluatorOfAllVehicles(dist_callback)
        # Solve the problem.
        assignment = routing.SolveWithParameters(search_parameters)
        if assignment:

            # Solution cost.
            print("Total distance: " + str(assignment.ObjectiveValue()) + "\n")

            # Inspect solution.
            # Only one route here; otherwise iterate from 0 to routing.vehicles() - 1.
            route_number = 0
            node = routing.Start(route_number)
            start_node = node
            route = ''

            while not routing.IsEnd(node):
                route += str(node) + ' -> '
                node = assignment.Value(routing.NextVar(node))
            route += '0'
#             print("Route:\n\n" + route)
        else:
            print('No solution found.')
    else:
        print('Specify an instance greater than 0.')
def create_data_array():
    cities = pd.read_csv('cities.csv')
    locations = []
    for i in range(len(cities)):
        locations.append([cities.X[i],cities.Y[i]])
    
    print("Length of locations: ", len(locations))
    return locations

if __name__ == '__main__':
    start = timer()
    main()
    end = timer()
    print("Time: ", end-start)

