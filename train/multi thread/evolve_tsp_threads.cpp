//---------------------------------------------------------------------------
//   Multi Expression Programming Software - with multiple subpopulations and threads
//   Copyright Mihai Oltean  (mihai.oltean@gmail.com)
//   Version 2022.01.27.0

//   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//---------------------------------------------------------------------------
//   Each subpopulation is evolved in a separate thread
//   The subpopulations have a circular structure
//   At the end of each generation we move, from each subpopulation, some individuals in the next one

//   I recommend to check the basic variant first (without subpopulations and threads)

//   Compiled with Microsoft Visual C++ 2019
//   Requires C++11 or newer (for thread support)

//   how to use it: 
//   just create a console project and copy-paste the content this file in the main file of the project

//   More info at:  https://www.mepx.org

//   Please reports any sugestions and/or bugs to:     mihai.oltean@gmail.com

//   Training data file must have the following format (see building1.txt and cancer1.txt):
//   building1 and cancer1 data are taken from PROBEN1

//   x11 x12 ... x1n f1
//   x21 x22 ....x2n f2
//   .............
//   xm1 xm2 ... xmn fm

//   where m is the number of training data
//   and n is the number of variables.

//--------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <thread>
#include <mutex>

#include "graphs.h"
#include "mep_rands.h"

#define NUM_OPERATORS 6

#define OP_ADD -1
#define OP_SUB -2
#define OP_MUL -3
#define OP_MIN -4
#define OP_MAX -5
#define OP_IFABCD -6

// +   -1
// -   -2
// *   -3
// /   -4

const char operators_string[NUM_OPERATORS][10] = { "+", "-", "*", "min", "max", "ifabcd"};

// variables
// some indexes in the variables array

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

#define num_variables 10

//---------------------------------------------------------------------------
struct t_code3{
	int op;				// either a variable, operator or constant; 
	// variables are indexed from 0: 0,1,2,...; 
	// constants are indexed from num_variables
	// operators are -1, -2, -3...
	int addr1, addr2, addr3, addr4;    // pointers to arguments
};
//---------------------------------------------------------------------------
struct t_mep_chromosome{
	t_code3 *prg;        // the program - a string of genes
	double *constants; // an array of constants

	double fitness;        // the fitness (or the error)
	int best_index;
};
//---------------------------------------------------------------------------
struct t_mep_parameters{
	unsigned int code_length;             // number of instructions in a t_mep_chromosome
	unsigned int num_generations;
	unsigned int num_sub_populations;       // number of subpopulations
	unsigned int sub_population_size;                // subpopulation size
	double mutation_probability, crossover_probability;
	unsigned int num_constants;
	double constants_min, constants_max;   // the array for constants
	double variables_probability, operators_probability, constants_probability;

	unsigned int num_threads; // num threads. 
	//for best performances the number of subpopulations should be multiple of num_threads.
	// num_thread should no exceed the number of processor cores.
};
//---------------------------------------------------------------------------
void allocate_chromosome(t_mep_chromosome &c, const t_mep_parameters &params)
{
	c.prg = new t_code3[params.code_length];
	if (params.num_constants)
		c.constants = new double[params.num_constants];
	else
		c.constants = NULL;
}
//---------------------------------------------------------------------------
void delete_chromosome(t_mep_chromosome &c)
{
	if (c.prg) {
		delete[] c.prg;
		c.prg = NULL;
	}
	if (c.constants) {
		delete[] c.constants;
		c.constants = NULL;
	}
}
//---------------------------------------------------------------------------
void copy_individual(t_mep_chromosome& dest, const t_mep_chromosome& source, const t_mep_parameters &params)
{
	for (unsigned int i = 0; i < params.code_length; i++)
		dest.prg[i] = source.prg[i];
	for (unsigned int i = 0; i < params.num_constants; i++)
		dest.constants[i] = source.constants[i];
	dest.fitness = source.fitness;
	dest.best_index = source.best_index;
}
//---------------------------------------------------------------------------
void generate_random_chromosome(t_mep_chromosome &a, const t_mep_parameters &params, t_seed& seed) // randomly initializes the individuals
{
	// generate constants first
	for (unsigned int c = 0; c < params.num_constants; c++)
		a.constants[c] = mep_real_rand(seed, params.constants_min, params.constants_max);

	// on the first position we can have only a variable or a constant
	double sum = params.variables_probability + params.constants_probability;
	double p = mep_real_rand(seed, 0, sum);

	if (p <= params.variables_probability)
		a.prg[0].op = mep_unsigned_int_rand(seed, 0, num_variables - 1);
	else
		a.prg[0].op = num_variables + mep_unsigned_int_rand(seed, 0, params.num_constants - 1);

	// for all other genes we put either an operator, variable or constant
	for (unsigned int i = 1; i < params.code_length; i++) {
		double p = mep_real_rand(seed, 0, 1);

		if (p <= params.operators_probability)
			a.prg[i].op = -1 -mep_unsigned_int_rand(seed, 0, NUM_OPERATORS - 1);       // an operator
		else {
			if (p <= params.operators_probability + params.variables_probability)
				a.prg[i].op = mep_unsigned_int_rand(seed, 0, num_variables - 1);     // a variable
			else
				a.prg[i].op = num_variables + mep_unsigned_int_rand(seed, 0, params.num_constants - 1); // index of a constant
		}
		a.prg[i].addr1 = mep_unsigned_int_rand(seed, 0, i - 1);
		a.prg[i].addr2 = mep_unsigned_int_rand(seed, 0, i - 1);
		a.prg[i].addr3 = mep_unsigned_int_rand(seed, 0, i - 1);
		a.prg[i].addr4 = mep_unsigned_int_rand(seed, 0, i - 1);
	}
}
//---------------------------------------------------------------------------
void mutation(t_mep_chromosome &a_chromosome, const t_mep_parameters &params, t_seed& seed) // mutate the individual
{
	// mutate each symbol with the given probability
	// first gene must be a variable or constant
	double p = mep_real_rand(seed, 0, 1);
	if (p < params.mutation_probability) {
		double sum = params.variables_probability + params.constants_probability;
		double p = mep_real_rand(seed, 0, sum);

		if (p <= params.variables_probability)
			a_chromosome.prg[0].op = mep_unsigned_int_rand(seed, 0, num_variables - 1);
		else
			a_chromosome.prg[0].op = num_variables + mep_unsigned_int_rand(seed, 0, params.num_constants - 1);
	}
	// other genes
	for (unsigned int i = 1; i < params.code_length; i++) {
		p = mep_real_rand(seed, 0, 1);      // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = mep_real_rand(seed, 0, 1);

			if (p <= params.operators_probability)
				a_chromosome.prg[i].op = -1 - mep_unsigned_int_rand(seed, 0, NUM_OPERATORS - 1);
			else
				if (p <= params.operators_probability + params.variables_probability)
					a_chromosome.prg[i].op = mep_unsigned_int_rand(seed, 0, num_variables - 1);
				else
					a_chromosome.prg[i].op = num_variables + mep_unsigned_int_rand(seed, 0, params.num_constants - 1); // index of a constant
		}

		p = mep_real_rand(seed, 0, 1);      // mutate the first address  (addr1)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr1 = mep_unsigned_int_rand(seed, 0, i - 1);

		p = mep_real_rand(seed, 0, 1);      // mutate the second address   (addr2)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr2 = mep_unsigned_int_rand(seed, 0, i - 1);
		p = mep_real_rand(seed, 0, 1);      // mutate the second address   (addr3)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr3 = mep_unsigned_int_rand(seed, 0, i - 1);
		p = mep_real_rand(seed, 0, 1);      // mutate the second address   (addr4)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr4 = mep_unsigned_int_rand(seed, 0, i - 1);
	}
	// mutate the constants
	for (unsigned int c = 0; c < params.num_constants; c++) {
		p = mep_real_rand(seed, 0, 1);
		if (p < params.mutation_probability)
			a_chromosome.constants[c] = mep_real_rand(seed, params.constants_min, params.constants_max);
	}
}
//---------------------------------------------------------------------------
void one_cut_point_crossover(const t_mep_chromosome &parent1, const t_mep_chromosome &parent2, 
	const t_mep_parameters &params, 
	t_mep_chromosome &offspring1, t_mep_chromosome &offspring2, t_seed& seed)
{
	unsigned int cutting_pct = 1 + mep_unsigned_int_rand(seed, 0, params.code_length - 3);
	for (unsigned int i = 0; i < cutting_pct; i++) {
		offspring1.prg[i] = parent1.prg[i];
		offspring2.prg[i] = parent2.prg[i];
	}
	for (unsigned int i = cutting_pct; i < params.code_length; i++) {
		offspring1.prg[i] = parent2.prg[i];
		offspring2.prg[i] = parent1.prg[i];
	}
	// now the constants
	if (params.num_constants) {
		cutting_pct = 1 + mep_unsigned_int_rand(seed, 0, params.num_constants - 3);
		for (unsigned int i = 0; i < cutting_pct; i++) {
			offspring1.constants[i] = parent1.constants[i];
			offspring2.constants[i] = parent2.constants[i];
		}
		for (unsigned int i = cutting_pct; i < params.num_constants; i++) {
			offspring1.constants[i] = parent2.constants[i];
			offspring2.constants[i] = parent1.constants[i];
		}
	}
}
//---------------------------------------------------------------------------
void uniform_crossover(const t_mep_chromosome &parent1, const t_mep_chromosome &parent2, 
	const t_mep_parameters &params, 
	t_mep_chromosome &offspring1, t_mep_chromosome &offspring2, t_seed& seed)
{
	for (unsigned int i = 0; i < params.code_length; i++)
		if (rand_int_01(seed)) {
			offspring1.prg[i] = parent1.prg[i];
			offspring2.prg[i] = parent2.prg[i];
		}
		else {
			offspring1.prg[i] = parent2.prg[i];
			offspring2.prg[i] = parent1.prg[i];
		}

	// constants
	for (unsigned int i = 0; i < params.num_constants; i++)
		if (rand_int_01(seed)) {
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
	if (((t_mep_chromosome *)a)->fitness > ((t_mep_chromosome *)b)->fitness)
		return 1;
	else
		if (((t_mep_chromosome *)a)->fitness < ((t_mep_chromosome *)b)->fitness)
			return -1;
		else
			return 0;
}
//---------------------------------------------------------------------------
unsigned int tournament_selection(const t_mep_chromosome *a_sub_pop, int sub_pop_size, unsigned int tournament_size, t_seed& seed)     // Size is the size of the tournament
{
	unsigned int p = mep_unsigned_int_rand(seed, 0, sub_pop_size - 1);
	for (unsigned int i = 1; i < tournament_size; i++) {
		int r = mep_unsigned_int_rand(seed, 0, sub_pop_size - 1);
		p = a_sub_pop[r].fitness < a_sub_pop[p].fitness ? r : p;
	}
	return p;
}
//---------------------------------------------------------------------------
double evaluate(const t_mep_chromosome &a_t_mep_chromosome, int code_length, double *vars_values, double *partial_values_array, int output_gene_index)
{
	for (int i = 0; i <= output_gene_index; i++)
		switch (a_t_mep_chromosome.prg[i].op) {
		case OP_ADD:// +
			partial_values_array[i] = partial_values_array[a_t_mep_chromosome.prg[i].addr1] + partial_values_array[a_t_mep_chromosome.prg[i].addr2];
			break;
		case OP_SUB:// -
			partial_values_array[i] = partial_values_array[a_t_mep_chromosome.prg[i].addr1] - partial_values_array[a_t_mep_chromosome.prg[i].addr2];
			break;
		case OP_MUL:// *
			partial_values_array[i] = partial_values_array[a_t_mep_chromosome.prg[i].addr1] * partial_values_array[a_t_mep_chromosome.prg[i].addr2];
			break;
		case OP_MAX:// max
			partial_values_array[i] = partial_values_array[a_t_mep_chromosome.prg[i].addr1] > partial_values_array[a_t_mep_chromosome.prg[i].addr2] ? partial_values_array[a_t_mep_chromosome.prg[i].addr1] : partial_values_array[a_t_mep_chromosome.prg[i].addr2];
			break;
		case OP_MIN:// min
			partial_values_array[i] = partial_values_array[a_t_mep_chromosome.prg[i].addr1] < partial_values_array[a_t_mep_chromosome.prg[i].addr2] ? partial_values_array[a_t_mep_chromosome.prg[i].addr1] : partial_values_array[a_t_mep_chromosome.prg[i].addr2];
			break;
		case OP_IFABCD:// a<b?c:d
			partial_values_array[i] = partial_values_array[a_t_mep_chromosome.prg[i].addr1] < 
										partial_values_array[a_t_mep_chromosome.prg[i].addr2] ? 
										partial_values_array[a_t_mep_chromosome.prg[i].addr3] : 
										partial_values_array[a_t_mep_chromosome.prg[i].addr4];
			break;
		default:
			if (a_t_mep_chromosome.prg[i].op < num_variables)
				partial_values_array[i] = vars_values[a_t_mep_chromosome.prg[i].op];
			else
				partial_values_array[i] = a_t_mep_chromosome.constants[a_t_mep_chromosome.prg[i].op - num_variables];
	}

	return partial_values_array[output_gene_index]; 
}
//---------------------------------------------------------------------------
void print_mep_chromosome(const t_mep_chromosome& a, const t_mep_parameters &params)
{
	printf("The t_mep_chromosome is:\n");

	for (unsigned int i = 0; i < params.num_constants; i++)
		printf("constants[%d] = %lf\n", i, a.constants[i]);

	for (unsigned int i = 0; i < params.code_length; i++)
		if (a.prg[i].op < 0) {
			if (a.prg[i].op == OP_IFABCD)
				printf("%d: %s %d %d %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].addr1, a.prg[i].addr2, a.prg[i].addr3, a.prg[i].addr4);
			else
				printf("%d: %s %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].addr1, a.prg[i].addr2);
		}
		else
			if (a.prg[i].op < num_variables)
				printf("%d: inputs[%d]\n", i, a.prg[i].op);
			else
				printf("%d: constants[%d]\n", i, a.prg[i].op - num_variables);

	printf("Fitness = %lf\n", a.fitness);
	printf("Best gene = %d\n", a.best_index);
}
//---------------------------------------------------------------------------
void print_chromosome_to_file(FILE* f, const t_mep_chromosome& a, const t_mep_parameters& params)
{
	for (unsigned int i = 0; i < params.code_length; i++)
		if (a.prg[i].op < 0) {
			if (a.prg[i].op == OP_IFABCD)
				fprintf(f, "%d: %d %d %d %d %d\n", i, a.prg[i].op, a.prg[i].addr1, a.prg[i].addr2, a.prg[i].addr3, a.prg[i].addr4);
			else
				fprintf(f, "%d: %d %d %d\n", i, a.prg[i].op, a.prg[i].addr1, a.prg[i].addr2);
		}
		else
			if (a.prg[i].op < num_variables)
				fprintf(f, "%d: %d\n", i, a.prg[i].op);
			else
				fprintf(f, "%d: %d %lf\n", i, a.prg[i].op, a.constants[a.prg[i].op - num_variables]);

	fprintf(f, "%d\n", a.best_index);
}
//---------------------------------------------------------------------------
void compute_local_variables(const t_graph &graph, int num_visited, int current_node, const int *node_visited, double *vars_values)
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
void fitness(t_mep_chromosome &individual, int code_length, 
		const t_graph *training_graphs, int num_training_graphs, 
		double * vars_values, double *partial_values_array)
{
	// fitness is the sum of errors over all training graphs.
	// error is the distance from the optimum

	// we cannot send the entire graph as a parameter to the function
	// so we extract only some summarized information from it
	// we fill the "vars_values" array with this summarized information
	// a function having vars_values as a parameter will be evolved

	individual.fitness = DBL_MAX;
	for (int g = 0; g < code_length; g++) {
		double error = 0;
		for (int k = 0; k < num_training_graphs; k++) {

			// set the value of variables
			vars_values[min_distance_in_matrix] = training_graphs[k].min_distance;
			vars_values[max_distance_in_matrix] = training_graphs[k].max_distance;
			vars_values[average_distance_in_matrix] = training_graphs[k].average_distance;
			vars_values[num_total_nodes] = training_graphs[k].num_nodes;

			// initialize memory for the path
			int* tsp_path = new int[training_graphs[k].num_nodes];
			int* node_visited = new int[training_graphs[k].num_nodes];
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
						double eval = evaluate(individual, code_length, vars_values, partial_values_array, g);
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
			error += (path_length - training_graphs[k].optimal_length) / training_graphs[k].optimal_length * 100;
			// keep it in percent from the optimal solution, otherwise we have scalling problems

			delete[] tsp_path;
			delete[] node_visited;
		}
		error /= (double)num_training_graphs; // average over the number of training graphs
		if (individual.fitness > error) {
			individual.fitness = error;
			individual.best_index = g;
		}
	}
}
//-----------------------------------------------------------------

void evolve_one_subpopulation(unsigned int *current_subpop_index, std::mutex* mutex,
		t_mep_chromosome ** sub_populations, unsigned int generation_index, 
	const t_mep_parameters &params, const t_graph *training_graphs, int num_training_graphs, 
	double* vars_values, t_seed* seed)
{
	unsigned int pop_index = 0;
	while (*current_subpop_index < params.num_sub_populations) {// still more subpopulations to evolve?

		while (!mutex->try_lock()) {}// create a lock so that multiple threads will not evolve the same sub population
		pop_index = *current_subpop_index;
		(*current_subpop_index)++;
		mutex->unlock();

		// pop_index is the index of the subpopulation evolved by the current thread

		t_mep_chromosome *a_sub_population = sub_populations[pop_index];

		t_mep_chromosome offspring1, offspring2;
		allocate_chromosome(offspring1, params);
		allocate_chromosome(offspring2, params);

		double *partial_values_array = new double[params.code_length];

		if (generation_index == 0) {
			for (unsigned int i = 0; i < params.sub_population_size; i++) {
				generate_random_chromosome(a_sub_population[i], params, seed[pop_index]);
				
				fitness(a_sub_population[i], params.code_length, training_graphs, num_training_graphs, vars_values, partial_values_array);

			}
			// sort ascendingly by fitness inside this population
			qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
		}
		else // next generations
			for (unsigned int k = 0; k < params.sub_population_size; k += 2) {
				// we increase by 2 because at each step we create 2 offspring

				// choose the parents using binary tournament
				unsigned int r1 = tournament_selection(a_sub_population, params.sub_population_size, 2, seed[pop_index]);
				unsigned int r2 = tournament_selection(a_sub_population, params.sub_population_size, 2, seed[pop_index]);
				// crossover
				double p_0_1 = mep_real_rand(seed[pop_index], 0, 1); // a random number between 0 and 1
				if (p_0_1 < params.crossover_probability)
					one_cut_point_crossover(a_sub_population[r1], a_sub_population[r2], params, offspring1, offspring2, seed[pop_index]);
				else {// no crossover so the offspring are a copy of the parents
					copy_individual(offspring1, a_sub_population[r1], params);
					copy_individual(offspring2, a_sub_population[r2], params);
				}
				// mutate the result and compute fitness
				mutation(offspring1, params, seed[pop_index]);
				fitness(offspring1, params.code_length, training_graphs, num_training_graphs, vars_values, partial_values_array);

				// mutate the other offspring too
				mutation(offspring2, params, seed[pop_index]);
				fitness(offspring2, params.code_length, training_graphs, num_training_graphs, vars_values, partial_values_array);

				// replace the worst in the population
				if (offspring1.fitness < a_sub_population[params.sub_population_size - 1].fitness) {
					copy_individual(a_sub_population[params.sub_population_size - 1], offspring1, params);
					qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
				}
				if (offspring2.fitness < a_sub_population[params.sub_population_size - 1].fitness) {
					copy_individual(a_sub_population[params.sub_population_size - 1], offspring2, params);
					qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
				}
			}

		delete_chromosome(offspring1);
		delete_chromosome(offspring2);
	}
}
//---------------------------------------------------------------------------
void start_steady_state(const t_mep_parameters &params, 
	const t_graph *training_graphs, int num_training_graphs,
	FILE *f_out)
{
	// a steady state model - 
	// Newly created inviduals replace the worst ones (if the offspring are better) in the same (sub) population.

	// allocate memory for all sub populations
	t_mep_chromosome **sub_populations; // an array of sub populations
	sub_populations = new t_mep_chromosome*[params.num_sub_populations];
	for (unsigned int p = 0; p < params.num_sub_populations; p++) {
		sub_populations[p] = new t_mep_chromosome[params.sub_population_size];
		for (unsigned int i = 0; i < params.sub_population_size; i++)
			allocate_chromosome(sub_populations[p][i], params); // allocate each individual in the subpopulation 
	}

	// allocate memory
	double** vars_values = new double*[params.num_threads];
	for (unsigned int t = 0; t < params.num_threads; t++)
		vars_values[t] = new double[num_variables];

	t_seed* seeds = new t_seed[params.num_sub_populations];
	for (unsigned int p = 0; p < params.num_sub_populations; p++)
		seeds[p].init(p);

	// an array of threads. Each sub population is evolved by a thread
	std::thread **mep_threads = new std::thread*[params.num_threads];
	// we create a fixed number of threads and each thread will take and evolve one subpopulation, then it will take another one
	std::mutex mutex;
	// we need a mutex to make sure that the same subpopulation will not be evolved twice by different threads

	unsigned int best_individual_subpop_index;
	// evolve for a fixed number of generations
	for (unsigned int generation = 0; generation < params.num_generations; generation++) { // for each generation

		//current_subpop_index = 0;
		unsigned int current_subpop_index = 0;
		for (unsigned int t = 0; t < params.num_threads; t++)
			mep_threads[t] = new std::thread(evolve_one_subpopulation, &current_subpop_index, &mutex, sub_populations, generation, params, training_graphs, num_training_graphs, vars_values[t], seeds);

		for (unsigned int t = 0; t < params.num_threads; t++) {
			mep_threads[t]->join();
			delete mep_threads[t];
		}

		// find the best individual
		best_individual_subpop_index = 0; // the index of the subpopulation containing the best invidual
		for (unsigned int p = 1; p < params.num_sub_populations; p++)
			if (sub_populations[p][0].fitness < sub_populations[best_individual_subpop_index][0].fitness)
				best_individual_subpop_index = p;
		// compute average
		double average_fitness = 0;
		for (unsigned int p = 0; p < params.num_sub_populations; p++)
			for (unsigned int i = 0; i < params.sub_population_size; i++)
				average_fitness += sub_populations[p][i].fitness;

		average_fitness /= (params.num_sub_populations * params.sub_population_size);

		printf("generation %d, best fitness = %lf; average fitness = %lf\n", generation, sub_populations[best_individual_subpop_index][0].fitness, average_fitness);
		fprintf(f_out, "%d %lf %lf\n", generation, sub_populations[best_individual_subpop_index][0].fitness, average_fitness);
		print_chromosome_to_file(f_out, sub_populations[best_individual_subpop_index][0], params);

		// now copy one individual from one population to the next one.
		// the copied invidual will replace the worst in the next one (if is better)

		for (unsigned int p = 0; p < params.num_sub_populations; p++) {
			unsigned int  k = mep_unsigned_int_rand(seeds[0], 0, params.sub_population_size - 1);// the individual to be copied
			// replace the worst in the next population (p + 1) - only if is better
			unsigned int index_next_pop = (p + 1) % params.num_sub_populations; // index of the next subpopulation (taken in circular order)
			if (sub_populations[p][k].fitness < sub_populations[index_next_pop][params.sub_population_size - 1].fitness) {
				copy_individual(sub_populations[index_next_pop][params.sub_population_size - 1], sub_populations[p][k], params);
				qsort((void *)sub_populations[index_next_pop], params.sub_population_size, sizeof(sub_populations[0][0]), sort_function);
			}
		}
	}
	delete[] seeds;
	delete[] mep_threads;

	// print best t_mep_chromosome

	print_mep_chromosome(sub_populations[best_individual_subpop_index][0], params);
	// free memory

	for (unsigned int p = 0; p < params.num_sub_populations; p++) {
		for (unsigned int i = 0; i < params.sub_population_size; i++)
			delete_chromosome(sub_populations[p][i]);
		delete[] sub_populations[p];
	}
	delete[] sub_populations;

	for (unsigned int t = 0; t < params.num_threads; t++)
		delete[] vars_values[t];
	delete[] vars_values;
}
//--------------------------------------------------------------------

int main(void)
{
	t_mep_parameters params;
	params.num_sub_populations = 20;
	params.sub_population_size = 10;			// the number of individuals in population  (must be an even number!)
	params.code_length = 50;
	params.num_generations = 30;				// the number of generations
	params.mutation_probability = 0.02;          // mutation probability
	params.crossover_probability = 0.9;         // crossover probability

	params.variables_probability = 0.4;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 10; // use 5 constants from -1 ... +1 interval
	params.constants_min = -1;
	params.constants_max = 1;

	params.num_threads = 20;

	t_graph *training_graphs = NULL;
	int num_training_graphs = 0;

	const char* path_str = "c:/Mihai/work/evolve-tsp-heuristics/data/";

	if (!read_training_data(path_str, training_graphs, num_training_graphs)) {
		printf("Cannot find input file(s) in path %s! Please specify the correct path!", path_str);
		getchar();
		return 1;
	}

	compute_global_variables(training_graphs, num_training_graphs);

	srand(0); 

	printf("evolving...\n");
	// store start time
	std::chrono::time_point<std::chrono::system_clock> start;
	start = std::chrono::system_clock::now();

	FILE* f_out = fopen("c:/temp/tsp_output.txt", "w");

	if (!f_out) {
		printf("Cannot create file for output! Please specify a valid path!\n");
		getchar();
		return 1;
	}

	start_steady_state(params, training_graphs, num_training_graphs, f_out);

	fclose(f_out);

	std::chrono::time_point<std::chrono::system_clock> end;
	end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	printf("run time = %lf\n", elapsed_seconds.count());

	delete_training_graphs(training_graphs, num_training_graphs);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------