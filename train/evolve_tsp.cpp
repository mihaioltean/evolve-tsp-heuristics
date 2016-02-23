//---------------------------------------------------------------------------
//   Designing Travelling Salesman Heuristics with Evolutionary Algorithms
//   Copyright (C) 2002-2015, Mihai Oltean  (mihai.oltean@gmail.com)
//   Version 2015.08.23

//   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//---------------------------------------------------------------------------

//   Compiled with Microsoft Visual C++ 2013
//   Also compiled with XCode 5.

//   New versions of this program will be available at  www.cs.ubbcluj.ro/~moltean

//   Please reports any sugestions and/or bugs to       mihai.oltean@gmail.com

//   The purpose of this program is to evolve a heuristic better than Nearest Neighbour
//   The evolved heuristic construct the path node by node, like NN 
//   but chooses the next node based on an evolved function

//   paper to read:
//   Oltean Mihai, Dumitrescu, D., Evolving TSP Heuristics using Multi Expression Programming, 
//   International Conference on Computational Sciences, ICCS'04, 6-9 June, Poland, 
//   Edited by M. Bubak, G.van Albada, P. Sloot, and J. Dongarra, Vol II, pp. 670-673, 
//   Springer-Verlag, Berlin, 2004

//   this source code is a little bit different compared to the source code used in the published paper
//   here we trained the evolved function against some instances in TSPLIB
//   originally we have trained the heuristic against some random TSP instances

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define num_operators 5

// +   -1
// -   -2
// *   -3
// max   -4
// min   -5

char operators_string[num_operators][10] = { "+", "-", "*", "min", "max" };

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
struct s_graph{
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
	int pop_size;
double mutation_probability, crossover_probability;
	int num_constants;
	double constants_min, constants_max;   // the array for constants
	double variables_probability, operators_probability, constants_probability;
};
//---------------------------------------------------------------------------
void allocate_chromosome(t_chromosome &c, t_parameters &params)
{
	c.prg = new t_code3[params.code_length];
	if (params.num_constants)
		c.constants = (double*)malloc(params.num_constants * sizeof(double));
	else
		c.constants = NULL;
}
//---------------------------------------------------------------------------
void delete_chromosome(t_chromosome &c)
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
void allocate_partial_expression_values(double **&expression_value, int num_training_data, int code_length)
{
	expression_value = new double*[code_length];
	for (int i = 0; i < code_length; i++)
		expression_value[i] = new double[num_training_data];
}
//---------------------------------------------------------------------------
void delete_partial_expression_values(double **&expression_value, int code_length)
{
	if (expression_value) {
		for (int i = 0; i < code_length; i++)
			delete[] expression_value[i];
		delete[] expression_value;
	}
}
//---------------------------------------------------------------------------
bool read_training_data(s_graph *&training_graphs, int &num_training_graphs)
{
	// training is done on 4 graphs
	num_training_graphs = 4;
	training_graphs = new s_graph[num_training_graphs];

	int k = 0; // count the graphs
	{
		FILE* f = fopen("data\\bayg29.tsp", "r");
		if (!f)
			return false;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double*[training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];
		// now read the data
		for (int i = 0; i < training_graphs[k].num_nodes - 1; i++)
			for (int j = i; j < training_graphs[k].num_nodes; j++)
				if (i != j) {
					fscanf(f, "%lf", &training_graphs[k].distance[i][j]);
					training_graphs[k].distance[j][i] = training_graphs[k].distance[i][j];
				}
				else
					training_graphs[k].distance[i][i] = 0;

		// now read the length of the shortest path
		fscanf(f, "%lf", &training_graphs[k].optimal_length);

		fclose(f);
	}
	{
		k++;
		FILE* f = fopen("data\\a280.tsp", "r");
		if (!f)
			return false;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double*[training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double *x = new double[training_graphs[k].num_nodes];
		double *y = new double[training_graphs[k].num_nodes];
		int index;
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			fscanf(f, "%d%lf%lf", &index, &x[i], &y[i]);
		// now read the length of the shortest path
		fscanf(f, "%lf", &training_graphs[k].optimal_length);
		fclose(f);

		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			for (int j = 0; j < training_graphs[k].num_nodes; j++)
				training_graphs[k].distance[i][j] = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));

		delete[] x;
		delete[] y;
	}
	{
		k++;
		FILE* f = fopen("data\\berlin52.tsp", "r");
		if (!f)
			return false;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double*[training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double *x = new double[training_graphs[k].num_nodes];
		double *y = new double[training_graphs[k].num_nodes];
		int index;
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			fscanf(f, "%d%lf%lf", &index, &x[i], &y[i]);
		// now read the length of the shortest path
		fscanf(f, "%lf", &training_graphs[k].optimal_length);
		fclose(f);

		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			for (int j = 0; j < training_graphs[k].num_nodes; j++)
				training_graphs[k].distance[i][j] = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));

		delete[] x;
		delete[] y;
	}
	{
		k++;
		FILE* f = fopen("data\\bier127.tsp", "r");
		if (!f)
			return false;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double*[training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double *x = new double[training_graphs[k].num_nodes];
		double *y = new double[training_graphs[k].num_nodes];
		int index;
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			fscanf(f, "%d%lf%lf", &index, &x[i], &y[i]);
		// now read the length of the shortest path
		fscanf(f, "%lf", &training_graphs[k].optimal_length);
		fclose(f);

		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			for (int j = 0; j < training_graphs[k].num_nodes; j++)
				training_graphs[k].distance[i][j] = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));

		delete[] x;
		delete[] y;
	}	return true;
}
//---------------------------------------------------------------------------
void delete_training_graphs(s_graph *&training_graphs, int num_training_graphs)
{
	if (training_graphs)
		for (int i = 0; i < num_training_graphs; i++) {
			for (int j = 0; j < training_graphs[i].num_nodes; j++)
				delete[] training_graphs[i].distance[j];
			delete[] training_graphs[i].distance;
		}
	delete[] training_graphs;
}
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
// evaluate Individual
void mutation(t_chromosome &a, t_parameters params, int num_variables) // mutate the individual
{
	// mutate each symbol with the given probability
	double p = rand() / (double)RAND_MAX;
	if (p < params.mutation_probability) {
		double sum = params.variables_probability + params.constants_probability;
		double p = rand() / (double)RAND_MAX * sum;

		if (p <= params.variables_probability)
			a.prg[0].op = rand() % num_variables;
		else
			a.prg[0].op = num_variables + rand() % params.num_constants;
	}

	for (int i = 1; i < params.code_length; i++) {
		p = rand() / (double)RAND_MAX;      // mutate the operator
		if (p < params.mutation_probability) {
			p = rand() / (double)RAND_MAX;

			if (p <= params.operators_probability)
				a.prg[i].op = -rand() % num_operators - 1;
			else
				if (p <= params.operators_probability + params.variables_probability)
					a.prg[i].op = rand() % num_variables;
				else
					a.prg[i].op = num_variables + rand() % params.num_constants; // index of a constant
		}

		p = rand() / (double)RAND_MAX;      // mutate the first address  (adr1)
		if (p < params.mutation_probability)
			a.prg[i].adr1 = rand() % i;

		p = rand() / (double)RAND_MAX;      // mutate the second address   (adr2)
		if (p < params.mutation_probability)
			a.prg[i].adr2 = rand() % i;
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
void print_chromosome(t_chromosome& a, t_parameters &params, int num_variables)
{
	printf("The t_chromosome is:\n");

	for (int i = 0; i < params.num_constants; i++)
		printf("constants[%d] = %lf\n", i, a.constants[i]);

	for (int i = 0; i < params.code_length; i++)
		if (a.prg[i].op < 0)
			printf("%d: %s %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].adr1, a.prg[i].adr2);
		else
			if (a.prg[i].op < num_variables)
				printf("%d: inputs[%d]\n", i, a.prg[i].op);
			else
				printf("%d: constants[%d]\n", i, a.prg[i].op - num_variables);

	printf("Fitness = %lf\n", a.fitness);
}
//---------------------------------------------------------------------------
int tournament_selection(t_chromosome *pop, int pop_size, int tournament_size)     // Size is the size of the tournament
{
	// tournament selection, 
	// returns the best out of tournament_size 

	int r, p;
	p = rand() % pop_size;
	for (int i = 1; i < tournament_size; i++) {
		r = rand() % pop_size;
		p = pop[r].fitness < pop[p].fitness ? r : p;
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
void compute_local_variables(s_graph &graph, int num_visited, int current_node, int *node_visited, double *vars_values)
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
void fitness(t_chromosome &individual, int code_length, s_graph *training_graphs, int num_training_graphs, int num_variables, double * vars_values, double *partial_values_array)
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
void start_steady_state_mep(t_parameters &params, s_graph *training_graphs, int &num_training_graphs, int num_variables, double* vars_values)       
// Steady-State evolutionary algorithm
// we work with one population and the newly created individuals will replace the existing worst individuals
{

	// allocate memory

	t_chromosome *population;
	population = new t_chromosome[params.pop_size];
	for (int i = 0; i < params.pop_size; i++)
		allocate_chromosome(population[i], params);

	t_chromosome offspring1, offspring2;
	allocate_chromosome(offspring1, params);
	allocate_chromosome(offspring2, params);

	double *partial_values_array = new double[params.code_length];

	// initialize
	for (int i = 0; i < params.pop_size; i++) {
		generate_random_chromosome(population[i], params, num_variables);
		fitness(population[i], params.code_length, training_graphs, num_training_graphs, num_variables, vars_values, partial_values_array);
	}
	// sort ascendingly by fitness
	qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);
	// fitness is computed as percent from the optimum
	printf("generation %d, best fitness = %lf percent from optimum\n", 0, population[0].fitness);

	for (int g = 1; g < params.num_generations; g++) {
		for (int k = 0; k < params.pop_size; k += 2) {
			// choose the parents using binary tournament
			int r1 = tournament_selection(population, params.pop_size, 2);
			int r2 = tournament_selection(population, params.pop_size, 2);
			// crossover
			double p = rand() / double(RAND_MAX);
			if (p < params.crossover_probability)
				one_cut_point_crossover(population[r1], population[r2], params, offspring1, offspring2);
			else {// no crossover so the offspring are a copy of the parents
				copy_individual(offspring1, population[r1], params);
				copy_individual(offspring2, population[r2], params);
			}
			// mutate the result
			mutation(offspring1, params, num_variables);
			fitness(offspring1, params.code_length, training_graphs, num_training_graphs, num_variables, vars_values, partial_values_array);

			mutation(offspring2, params, num_variables);
			fitness(offspring2, params.code_length, training_graphs, num_training_graphs, num_variables, vars_values, partial_values_array);

			// replace the worst in the population
			if (offspring1.fitness < population[params.pop_size - 1].fitness) {
				copy_individual(population[params.pop_size - 1], offspring1, params);
				qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);
			}
			if (offspring2.fitness < population[params.pop_size - 1].fitness) {
				copy_individual(population[params.pop_size - 1], offspring2, params);
				qsort((void *)population, params.pop_size, sizeof(population[0]), sort_function);
			}
		}
		printf("generation %d, best fitness = %lf percent from optimum\n", g, population[0].fitness);
	}
	// print best t_chromosome
	print_chromosome(population[0], params, num_variables);

	// free memory
	delete_chromosome(offspring1);
	delete_chromosome(offspring2);

	for (int i = 0; i < params.pop_size; i++)
		delete_chromosome(population[i]);
	delete[] population;

	delete_training_graphs(training_graphs, num_training_graphs);
	delete[] partial_values_array;
}
//--------------------------------------------------------------------
void compute_global_variables(s_graph *training_graphs, int num_training_graphs)
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
int main(void)
{
	t_parameters params;

	params.pop_size = 50;						    // the number of individuals in population  (must be an even number!)
	params.code_length = 50;
	params.num_generations = 100;					// the number of generations
	params.mutation_probability = 0.01;              // mutation probability
	params.crossover_probability = 0.9;             // crossover probability

	params.variables_probability = 0.4;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 3; // use 3 constants from -1 ... +1 interval
	params.constants_min = -1;
	params.constants_max = 1;

	s_graph *training_graphs = NULL;
	int num_training_graphs = 0;

	if (!read_training_data(training_graphs, num_training_graphs)) {
		printf("Cannot find input file(s)! Please specify the full path!");
		getchar();
		return 1;
	}

	compute_global_variables(training_graphs, num_training_graphs);

	int num_variables = 10;
	double* vars_values = new double[num_variables];

	srand(1);

	printf("evolving...\n");
	start_steady_state_mep(params, training_graphs, num_training_graphs, num_variables, vars_values);

	delete[] vars_values;

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------