//
//  genetic_mep.hpp
//  mep_evolve
//
//  Created by Virginia Niculescu on 23/02/16.
//  Copyright © 2016 Virginia Niculescu. All rights reserved.
//

#ifndef genetic_mep_hpp
#define genetic_mep_hpp


#include "config.hpp"



//---------------------------------------------------------------------------
void copy_individual(t_chromosome& dest, const t_chromosome& source, t_parameters &params);
//---------------------------------------------------------------------------
void generate_random_chromosome(t_chromosome &a, t_parameters &params, int num_variables); // randomly initializes the individuals

//---------------------------------------------------------------------------
void mutation(t_chromosome &a_chromosome, t_parameters params, int num_variables); // mutate the individual

//---------------------------------------------------------------------------
void one_cut_point_crossover(const t_chromosome &parent1, const t_chromosome &parent2, t_parameters &params, t_chromosome &offspring1, t_chromosome &offspring2);
//---------------------------------------------------------------------------
void uniform_crossover(const t_chromosome &parent1, const t_chromosome &parent2, t_parameters &params, t_chromosome &offspring1, t_chromosome &offspring2);
//---------------------------------------------------------------------------
int sort_function(const void *a, const void *b);
//comparator for quick sort
//---------------------------------------------------------------------------
int tournament_selection(t_chromosome *a_sub_pop, int sub_pop_size, int tournament_size)     // Size is the size of the tournament
;
//---------------------------------------------------------------------------
double evaluate(t_chromosome &a_t_chromosome, int code_length, int num_variables, double *vars_values, double *partial_values_array);
//---------------------------------------------------------------------------
void compute_local_variables(t_graph &graph, int num_visited, int current_node, int *node_visited, double *vars_values);
//--------------------------------------------------------------------
void fitness(t_chromosome &individual, int code_length, t_graph *training_graphs, int num_training_graphs, int num_variables, double * vars_values, double *partial_values_array);
//-----------------------------------------------------------------

void evolve_one_subpopulation(int *current_subpop_index, std::mutex* mutex, t_chromosome ** sub_populations, int generation_index, t_parameters *params, t_graph *training_graphs, int num_training_graphs, int num_variables, double* vars_values);
//--------------------------------------------------------------------
void compute_global_variables(t_graph *training_graphs, int num_training_graphs);
//--------------------------------------------------------------------


#endif /* genetic_mep_hpp */
