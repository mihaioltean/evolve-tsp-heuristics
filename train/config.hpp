//
//  config.hpp
//  mep_evolve
//
//  Created by Virginia Niculescu on 23/02/16.
//  Copyright © 2016 Virginia Niculescu. All rights reserved.
//

#ifndef config_hpp
#define config_hpp

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <thread>
#include <mutex>
#include <float.h>

#define num_operators 6

// +   -1
// -   -2
// *   -3
// /   -4

extern char operators_string[num_operators][10];

// variables
// some indecses in the variables array

#define min_distance_in_matrix 0     // 0 - min distance in the distance matrix
#define max_distance_in_matrix 1     // 1 - max distance in the distance matrix
#define average_distance_in_matrix 2 // 2 - average distance in the matrix
#define length_so_far 3              // 3 - length so far
#define distance_to_next_node 4      // 4 - distance to the next node
#define distance_to_nearest_node 5   // 5 - distance to the nearest node
#define distance_to_fartest_node 6        // 6 - distance to the fartest node
#define average_distance_to_unvisited 7
#define num_visited_nodes 8          // number of visisted nodes so far
#define num_total_nodes 9            // total number of nodes of that graph

//---------------------------------------------------------------------------
struct t_graph{
    int num_nodes;			// number of nodes
    double **distance;      // distance matrix
    double min_distance, max_distance, average_distance;
    double optimal_length;  // known optimal solution
};



//---------------------------------------------------------------------------
struct t_code3{
    int op;				// either a variable, operator or constant;
    // variables are indexed from 0: 0,1,2,...;
    // constants are indexed from num_variables
    // operators are -1, -2, -3...
    int adr1, adr2;    // pointers to arguments
};
//---------------------------------------------------------------------------
struct t_chromosome{
    t_code3 *prg;        // the program - a string of genes
    double *constants; // an array of constants
    
    double fitness;        // the fitness (or the error)
};
//---------------------------------------------------------------------------
struct t_parameters{
    int code_length;             // number of instructions in a t_chromosome
    int num_generations;
    int num_sub_populations;       // number of subpopulations
    int sub_population_size;                // subpopulation size
    double mutation_probability, crossover_probability;
    int num_constants;
    double constants_min, constants_max;   // the array for constants
    double variables_probability, operators_probability, constants_probability;
    
   
};



//---------------------------------------------------------------------------
void allocate_chromosome(t_chromosome &c, t_parameters &params);
//---------------------------------------------------------------------------
void delete_chromosome(t_chromosome &c);
//---------------------------------------------------------------------------
void allocate_training_data(double **&data, double *&target, int num_training_data, int num_variables);
//---------------------------------------------------------------------------
void allocate_partial_expression_values(double ***&expression_value, int num_training_data, int code_length, int num_threads);
//---------------------------------------------------------------------------
void delete_partial_expression_values(double ***&expression_value, int code_length, int num_threads);
//---------------------------------------------------------------------------


#endif /* config_hpp */
