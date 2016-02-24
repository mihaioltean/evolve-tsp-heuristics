//
//  genetic_mep.cpp
//  mep_evolve
//
//  Created by Virginia Niculescu on 23/02/16.
//  Copyright © 2016 Virginia Niculescu. All rights reserved.
//

#include "genetic_mep.hpp"
#include "config.hpp"
#include <math.h>


//---------------------------------------------------------------------------
void copy_individual(t_chromosome& dest, const t_chromosome& source, t_parameters &params)
{
    for (int i = 0; i < params.code_length; i++)
        dest.prg[i] = source.prg[i];
    for (int i = 0; i < params.num_constants; i++)
        dest.constants[i] = source.constants[i];
    dest.fitness = source.fitness;
}
//---------------------------------------------------------------------------
void generate_random_chromosome(t_chromosome &a, t_parameters &params, int num_variables) // randomly initializes the individuals
{
    // generate constants first
    for (int c = 0; c < params.num_constants; c++)
        a.constants[c] = rand() / double(RAND_MAX) * (params.constants_max - params.constants_min) + params.constants_min;
    
    // on the first position we can have only a variable or a constant
    double sum = params.variables_probability + params.constants_probability;
    double p = rand() / (double)RAND_MAX * sum;
    
    if (p <= params.variables_probability)
        a.prg[0].op = rand() % num_variables;
    else
        a.prg[0].op = num_variables + rand() % params.num_constants;
    
    // for all other genes we put either an operator, variable or constant
    for (int i = 1; i < params.code_length; i++) {
        double p = rand() / (double)RAND_MAX;
        
        if (p <= params.operators_probability)
            a.prg[i].op = -rand() % num_operators - 1;        // an operator
        else
            if (p <= params.operators_probability + params.variables_probability)
                a.prg[i].op = rand() % num_variables;     // a variable
            else
                a.prg[i].op = num_variables + rand() % params.num_constants; // index of a constant
        
        a.prg[i].adr1 = rand() % i;
        a.prg[i].adr2 = rand() % i;
    }
}
//---------------------------------------------------------------------------
void mutation(t_chromosome &a_chromosome, t_parameters params, int num_variables) // mutate the individual
{
    // mutate each symbol with the given probability
    // first gene must be a variable or constant
    double p = rand() / (double)RAND_MAX;
    if (p < params.mutation_probability) {
        double sum = params.variables_probability + params.constants_probability;
        double p = rand() / (double)RAND_MAX * sum;
        
        if (p <= params.variables_probability)
            a_chromosome.prg[0].op = rand() % num_variables;
        else
            a_chromosome.prg[0].op = num_variables + rand() % params.num_constants;
    }
    // other genes
    for (int i = 1; i < params.code_length; i++) {
        p = rand() / (double)RAND_MAX;      // mutate the operator
        if (p < params.mutation_probability) {
            // we mutate it, but we have to decide what we put here
            p = rand() / (double)RAND_MAX;
            
            if (p <= params.operators_probability)
                a_chromosome.prg[i].op = -rand() % num_operators - 1;
            else
                if (p <= params.operators_probability + params.variables_probability)
                    a_chromosome.prg[i].op = rand() % num_variables;
                else
                    a_chromosome.prg[i].op = num_variables + rand() % params.num_constants; // index of a constant
        }
        
        p = rand() / (double)RAND_MAX;      // mutate the first address  (adr1)
        if (p < params.mutation_probability)
            a_chromosome.prg[i].adr1 = rand() % i;
        
        p = rand() / (double)RAND_MAX;      // mutate the second address   (adr2)
        if (p < params.mutation_probability)
            a_chromosome.prg[i].adr2 = rand() % i;
    }
    // mutate the constants
    for (int c = 0; c < params.num_constants; c++) {
        p = rand() / (double)RAND_MAX;
        if (p < params.mutation_probability)
            a_chromosome.constants[c] = rand() / double(RAND_MAX) * (params.constants_max - params.constants_min) + params.constants_min;
    }
    
}
//---------------------------------------------------------------------------
void one_cut_point_crossover(const t_chromosome &parent1, const t_chromosome &parent2, t_parameters &params, t_chromosome &offspring1, t_chromosome &offspring2)
{
    int cutting_pct = rand() % params.code_length;
    for (int i = 0; i < cutting_pct; i++) {
        offspring1.prg[i] = parent1.prg[i];
        offspring2.prg[i] = parent2.prg[i];
    }
    for (int i = cutting_pct; i < params.code_length; i++) {
        offspring1.prg[i] = parent2.prg[i];
        offspring2.prg[i] = parent1.prg[i];
    }
    // now the constants
    if (params.num_constants) {
        cutting_pct = rand() % params.num_constants;
        for (int i = 0; i < cutting_pct; i++) {
            offspring1.constants[i] = parent1.constants[i];
            offspring2.constants[i] = parent2.constants[i];
        }
        for (int i = cutting_pct; i < params.num_constants; i++) {
            offspring1.constants[i] = parent2.constants[i];
            offspring2.constants[i] = parent1.constants[i];
        }
    }
}
//---------------------------------------------------------------------------
void uniform_crossover(const t_chromosome &parent1, const t_chromosome &parent2, t_parameters &params, t_chromosome &offspring1, t_chromosome &offspring2)
{
    for (int i = 0; i < params.code_length; i++)
        if (rand() % 2) {
            offspring1.prg[i] = parent1.prg[i];
            offspring2.prg[i] = parent2.prg[i];
        }
        else {
            offspring1.prg[i] = parent2.prg[i];
            offspring2.prg[i] = parent1.prg[i];
        }
    
    // constants
    for (int i = 0; i < params.num_constants; i++)
        if (rand() % 2) {
            offspring1.constants[i] = parent1.constants[i];
            offspring2.constants[i] = parent2.constants[i];
        }
        else {
            offspring1.constants[i] = parent2.constants[i];
            offspring2.constants[i] = parent1.constants[i];
        }
}
//---------------------------------------------------------------------------
int sort_function(const void *a, const void *b)
{// comparator for quick sort
    if (((t_chromosome *)a)->fitness > ((t_chromosome *)b)->fitness)
        return 1;
    else
        if (((t_chromosome *)a)->fitness < ((t_chromosome *)b)->fitness)
            return -1;
        else
            return 0;
}
//---------------------------------------------------------------------------
int tournament_selection(t_chromosome *a_sub_pop, int sub_pop_size, int tournament_size)     // Size is the size of the tournament
{
    int r, p;
    p = rand() % sub_pop_size;
    for (int i = 1; i < tournament_size; i++) {
        r = rand() % sub_pop_size;
        p = a_sub_pop[r].fitness < a_sub_pop[p].fitness ? r : p;
    }
    return p;
}
//---------------------------------------------------------------------------
double evaluate(t_chromosome &a_t_chromosome, int code_length, int num_variables, double *vars_values, double *partial_values_array)
{
    for (int i = 0; i < code_length; i++)
        switch (a_t_chromosome.prg[i].op) {
            case -1:// +
                partial_values_array[i] = partial_values_array[a_t_chromosome.prg[i].adr1] + partial_values_array[a_t_chromosome.prg[i].adr2];
                break;
            case -2:// -
                partial_values_array[i] = partial_values_array[a_t_chromosome.prg[i].adr1] - partial_values_array[a_t_chromosome.prg[i].adr2];
                break;
            case -3:// *
                partial_values_array[i] = partial_values_array[a_t_chromosome.prg[i].adr1] * partial_values_array[a_t_chromosome.prg[i].adr2];
                break;
            case -4:// max
                partial_values_array[i] = partial_values_array[a_t_chromosome.prg[i].adr1] > partial_values_array[a_t_chromosome.prg[i].adr2] ? partial_values_array[a_t_chromosome.prg[i].adr1] : partial_values_array[a_t_chromosome.prg[i].adr2];
                break;
            case -5:// min
                partial_values_array[i] = partial_values_array[a_t_chromosome.prg[i].adr1] < partial_values_array[a_t_chromosome.prg[i].adr2] ? partial_values_array[a_t_chromosome.prg[i].adr1] : partial_values_array[a_t_chromosome.prg[i].adr2];
                break;
            default:
                if (a_t_chromosome.prg[i].op < num_variables)
                    partial_values_array[i] = vars_values[a_t_chromosome.prg[i].op];
                else
                    partial_values_array[i] = a_t_chromosome.constants[a_t_chromosome.prg[i].op - num_variables];
        }
    
    return partial_values_array[code_length - 1]; // last gene is the one providing the output
}

//---------------------------------------------------------------------------
void compute_local_variables(t_graph &graph, int num_visited, int current_node, int *node_visited, double *vars_values)
{
    vars_values[distance_to_nearest_node] = 1E307; // min distance to unvisited
    vars_values[distance_to_fartest_node] = 0;     // max_distance to unvisited
    vars_values[average_distance_to_unvisited] = 0;  // average distance from current to all unvisited
    
    for (int i = 0; i < graph.num_nodes; i++)
        if (!node_visited[i]) {
            if (vars_values[distance_to_nearest_node] > graph.distance[current_node][i])
                vars_values[distance_to_nearest_node] = graph.distance[current_node][i];
            if (vars_values[distance_to_fartest_node] < graph.distance[current_node][i])
                vars_values[distance_to_fartest_node] = graph.distance[current_node][i];
            
            vars_values[average_distance_to_unvisited] += graph.distance[current_node][i];
        }
    
    vars_values[average_distance_to_unvisited] /= (double)(graph.num_nodes - num_visited);
}
//--------------------------------------------------------------------
void fitness(t_chromosome &individual, int code_length, t_graph *training_graphs, int num_training_graphs, int num_variables, double * vars_values, double *partial_values_array)
{
    // fitness is the sum of errors over all training graphs.
    // error is the distance from the optimum
    
    // we cannot send the entire graph as a parameter to the function
    // so we extract only some summarized information from it
    // we fill the "vars_values" array with this summarized information
    // a function having vars_values as a parameter will be evolved
    
    individual.fitness = 0;
    
    for (int k = 0; k < num_training_graphs; k++) {
        
        // set the value of variables
        vars_values[min_distance_in_matrix] = training_graphs[k].min_distance;
        vars_values[max_distance_in_matrix] = training_graphs[k].max_distance;
        vars_values[average_distance_in_matrix] = training_graphs[k].average_distance;
        vars_values[num_total_nodes] = training_graphs[k].num_nodes;
        
        // initialize memory for the path
        int *tsp_path = new int[training_graphs[k].num_nodes];
        int *node_visited = new int[training_graphs[k].num_nodes];
        for (int i = 0; i < training_graphs[k].num_nodes; i++)
            node_visited[i] = 0;
        
        // init path
        int count_nodes = 0;
        tsp_path[0] = 0; // first node in the path is 0;
        node_visited[0] = 1;
        count_nodes++;
        double path_length = 0;
        
        while (count_nodes < training_graphs[k].num_nodes) {
            // fill the path node by node
            double min_eval = DBL_MAX;
            int best_node = -1;
            // compute other statistics
            vars_values[length_so_far] = path_length;
            vars_values[num_visited_nodes] = count_nodes;
            compute_local_variables(training_graphs[k], count_nodes, tsp_path[count_nodes - 1], node_visited, vars_values);
            
            for (int nod = 0; nod < training_graphs[k].num_nodes; nod++)
                // consider each unvisited node
                if (!node_visited[nod]) {// not visited yet
                    vars_values[distance_to_next_node] = training_graphs[k].distance[tsp_path[count_nodes - 1]][nod];
                    double eval = evaluate(individual, code_length, num_variables, vars_values, partial_values_array);
                    if (eval < min_eval) {
                        best_node = nod; // keep the one with minimal evaluation
                        min_eval = eval;
                    }
                }
            // add the best next node to the path
            path_length += training_graphs[k].distance[tsp_path[count_nodes - 1]][best_node];
            node_visited[best_node] = 1;
            
            tsp_path[count_nodes] = best_node;
            count_nodes++;
        }
        // connect the last with the first in the path
        path_length += training_graphs[k].distance[tsp_path[count_nodes - 1]][tsp_path[0]];
        individual.fitness += (path_length - training_graphs[k].optimal_length) / training_graphs[k].optimal_length * 100;
        // keep it in percent from the optimal solution, otherwise we have scalling problems
        
        delete[] tsp_path;
        delete[] node_visited;
    }
    individual.fitness /= (double)num_training_graphs; // average over the number of training graphs
}

//-----------------------------------------------------------------
int nextPopulation(int *current_subpop_index, std::mutex* mutex){
    std::lock_guard<std::mutex> lock(*mutex);
    int pop_index = *current_subpop_index;
    (*current_subpop_index)++;
    return pop_index;
}

//-----------------------------------------------------------------

void evolve_one_subpopulation(int *current_subpop_index, std::mutex* mutex, t_chromosome ** sub_populations, int generation_index, t_parameters *params, t_graph *training_graphs, int num_training_graphs, int num_variables, double* vars_values)


{
    
    int pop_index = 0;
    while (*current_subpop_index < params->num_sub_populations) {// still more subpopulations to evolve?
        
        /*
        {
        //
        while (!mutex->try_lock()) {}// create a lock so that multiple threads will not evolve the same sub population
        
        pop_index = *current_subpop_index;
        (*current_subpop_index)++;
       
        mutex->unlock();
        }*/
        
        pop_index = nextPopulation(current_subpop_index, mutex);
        
        // pop_index is the index of the subpopulation evolved by the current thread
        
        t_chromosome *a_sub_population = sub_populations[pop_index];
        
        t_chromosome offspring1, offspring2;
        allocate_chromosome(offspring1, *params);
        allocate_chromosome(offspring2, *params);
        
        double *partial_values_array = new double[params->code_length];
        
        if (generation_index == 0) {
            for (int i = 0; i < params->sub_population_size; i++) {
                generate_random_chromosome(a_sub_population[i], *params, num_variables);
                
                fitness(a_sub_population[i], params->code_length, training_graphs, num_training_graphs, num_variables, vars_values, partial_values_array);
                
            }
            // sort ascendingly by fitness inside this population
            qsort((void *)a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]), sort_function);
        }
        else // next generations
            for (int k = 0; k < params->sub_population_size; k += 2) {
                // we increase by 2 because at each step we create 2 offspring
                
                // choose the parents using binary tournament
                int r1 = tournament_selection(a_sub_population, params->sub_population_size, 2);
                int r2 = tournament_selection(a_sub_population, params->sub_population_size, 2);
                // crossover
                double p_0_1 = rand() / double(RAND_MAX); // a random number between 0 and 1
                if (p_0_1 < params->crossover_probability)
                    one_cut_point_crossover(a_sub_population[r1], a_sub_population[r2], *params, offspring1, offspring2);
                else {// no crossover so the offspring are a copy of the parents
                    copy_individual(offspring1, a_sub_population[r1], *params);
                    copy_individual(offspring2, a_sub_population[r2], *params);
                }
                // mutate the result and compute fitness
                mutation(offspring1, *params, num_variables);
                fitness(offspring1, params->code_length, training_graphs, num_training_graphs, num_variables, vars_values, partial_values_array);
                
                // mutate the other offspring too
                mutation(offspring2, *params, num_variables);
                fitness(offspring2, params->code_length, training_graphs, num_training_graphs, num_variables, vars_values, partial_values_array);
                
                // replace the worst in the population
                if (offspring1.fitness < a_sub_population[params->sub_population_size - 1].fitness) {
                    copy_individual(a_sub_population[params->sub_population_size - 1], offspring1, *params);
                    qsort((void *)a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]), sort_function);
                }
                if (offspring2.fitness < a_sub_population[params->sub_population_size - 1].fitness) {
                    copy_individual(a_sub_population[params->sub_population_size - 1], offspring2, *params);
                    qsort((void *)a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]), sort_function);
                }
            }
        
        delete_chromosome(offspring1);
        delete_chromosome(offspring2);
    }
}
//--------------------------------------------------------------------
void compute_global_variables(t_graph *training_graphs, int num_training_graphs)
{
    // here we compute the min, max and average distance in a given graph
    for (int k = 0; k < num_training_graphs; k++) {
        training_graphs[k].max_distance = training_graphs[k].distance[0][1];
        training_graphs[k].min_distance = training_graphs[k].distance[0][1];
        training_graphs[k].average_distance = 0;
        for (int i = 0; i < training_graphs[k].num_nodes; i++)
            for (int j = 0; j < training_graphs[k].num_nodes; j++)
                if (i != j) {
                    training_graphs[k].average_distance += training_graphs[k].distance[i][j];
                    if (training_graphs[k].max_distance < training_graphs[k].distance[i][j])
                        training_graphs[k].max_distance = training_graphs[k].distance[i][j];
                    if (training_graphs[k].min_distance > training_graphs[k].distance[i][j])
                        training_graphs[k].min_distance = training_graphs[k].distance[i][j];
                }
        training_graphs[k].average_distance /= (training_graphs[k].num_nodes * training_graphs[k].num_nodes) - training_graphs[k].num_nodes;
    }
}
//--------------------------------------------------------------------

