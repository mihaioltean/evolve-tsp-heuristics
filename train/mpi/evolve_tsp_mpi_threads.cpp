//---------------------------------------------------------------------------
//   Multi Expression Programming Software - with multiple subpopulations and threads
//   Copyright Mihai Oltean  (mihai.oltean@gmail.com)
//   Version 2016.04.02.2 // year.month.day.build#

//   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//---------------------------------------------------------------------------
//   Each subpopulation is evolved in a separate thread
//   The subpopulations have a circular structure
//   At the end of each generation we move, from each subpopulation, some individuals in the next one

//   I recommend to check the basic variant first (without subpopulations and threads)

//   Compiled with Microsoft Visual C++ 2013
//   Also compiled with XCode 5.
//   Requires C++11 or newer (for thread support)

//   how to use it: 
//   just create a console project and copy-paste the content this file in the main file of the project

//   More info at:  http://www.mep.cs.ubbcluj.ro

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

//#define USE_THREADS
#define USE_MPI

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef USE_THREADS
#include <thread>
#include <mutex>
#endif

#include <float.h>

//#include <vld.h> // for detecting memory leaks in VC++



#ifdef USE_MPI

#include "mpi.h"

#endif


#define num_operators 7

// +   -1
// -   -2
// *   -3
// min   -4
// max   -5

char operators_string[num_operators][10] = { "+", "-", "*", "min", "max", "sin", "cos"};

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

	//------------------------------------------------------------------------------
	void to_string(char * s_dest, int code_length, int num_constants)
	{
		char tmp_s[100];
		s_dest[0] = 0;
		for (int i = 0; i < code_length; i++) {
			sprintf(tmp_s, "%d ", prg[i].op);
			strcat(s_dest, tmp_s);
			sprintf(tmp_s, "%d ", prg[i].adr1);
			strcat(s_dest, tmp_s);
			sprintf(tmp_s, "%d ", prg[i].adr2);
			strcat(s_dest, tmp_s);
		}
		for (int i = 0; i < num_constants; i++) {
			sprintf(tmp_s, "%lg ", constants[i]);
			strcat(s_dest, tmp_s);
		}
		sprintf(tmp_s, "%lg ", fitness);
		strcat(s_dest, tmp_s);
	}
	//------------------------------------------------------------------------------
	void from_string(char* s_source, int code_length, int num_constants)
	{
		int num_consumed = 0;
		for (int i = 0; i < code_length; i++) {
			sscanf(s_source, "%d%d%d%n", &prg[i].op, &prg[i].adr1, &prg[i].adr2, &num_consumed);
			s_source += num_consumed;
		}
		for (int i = 0; i < num_constants; i++) {
			sscanf(s_source, "%lf%n", &constants[i], &num_consumed);
			s_source += num_consumed;
		}
		sscanf(s_source, "%lf", &fitness);
	}
	//------------------------------------------------------------------------------

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
#ifdef USE_THREADS
	int num_threads; // num threads.
#endif
	//for best performances the number of subpopulations should be multiple of num_threads.
	// num_thread should no exceed the number of processor cores.
};
//---------------------------------------------------------------------------
void allocate_chromosome(t_chromosome &c, t_parameters &params)
{
	c.prg = new t_code3[params.code_length];
	if (params.num_constants)
		c.constants = new double[params.num_constants];
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
void allocate_training_data(double **&data, double *&target, int num_training_data, int num_variables)
{
	target = new double[num_training_data];
	data = new double*[num_training_data];
	for (int i = 0; i < num_training_data; i++)
		data[i] = new double[num_variables];
}
//---------------------------------------------------------------------------
void allocate_partial_expression_values(double ***&expression_value, int num_training_data, int code_length, int num_threads)
{
	// partial values are stored in a matrix of size: code_length x num_training_data
	// for each thread we have a separate matrix
	expression_value = new double**[num_threads];
	for (int t = 0; t < num_threads; t++) {
		expression_value[t] = new double*[code_length];
		for (int i = 0; i < code_length; i++)
			expression_value[t][i] = new double[num_training_data];
	}
}
//---------------------------------------------------------------------------
void delete_partial_expression_values(double ***&expression_value, int code_length, int num_threads)
{
	if (expression_value) {
		for (int t = 0; t < num_threads; t++) {
			for (int i = 0; i < code_length; i++)
				delete[] expression_value[t][i];
			delete[] expression_value[t];
		}
		delete[] expression_value;
	}
}
//---------------------------------------------------------------------------
bool read_training_data(t_graph *&training_graphs, int &num_training_graphs)
{
	// training is done on 4 graphs
	num_training_graphs = 4;
	training_graphs = new t_graph[num_training_graphs];

	int k = 0; // count the graphs
	{
		FILE* f = fopen("data//bayg29.tsp", "r");
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
		FILE* f = fopen("data//a280.tsp", "r");
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
		FILE* f = fopen("data//berlin52.tsp", "r");
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
		FILE* f = fopen("data//bier127.tsp", "r");
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
void delete_training_graphs(t_graph *&training_graphs, int num_training_graphs)
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
		case -6:// sin
			partial_values_array[i] = sin(partial_values_array[a_t_chromosome.prg[i].adr1]);
			break;
		case -7:// cos
			partial_values_array[i] = cos(partial_values_array[a_t_chromosome.prg[i].adr1]);
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

			for (int node = 0; node < training_graphs[k].num_nodes; node++)
				// consider each unvisited node
				if (!node_visited[node]) {// not visited yet
					vars_values[distance_to_next_node] = training_graphs[k].distance[tsp_path[count_nodes - 1]][node];
					double eval = evaluate(individual, code_length, num_variables, vars_values, partial_values_array);
					if (eval < min_eval) {
						best_node = node; // keep the one with minimal evaluation
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
#ifdef USE_THREADS
void evolve_one_subpopulation(int *current_subpop_index, std::mutex* mutex, t_chromosome ** sub_populations, int generation_index, t_parameters *params, t_graph *training_graphs, int num_training_graphs, int num_variables, double* vars_values)
#else
void evolve_one_subpopulation(int *current_subpop_index, t_chromosome ** sub_populations, int generation_index, t_parameters *params, t_graph *training_graphs, int num_training_graphs, int num_variables, double* vars_values)
#endif
{
	int pop_index = 0;
	while (*current_subpop_index < params->num_sub_populations) {// still more subpopulations to evolve?
#ifdef USE_THREADS
		while (!mutex->try_lock()) {}// create a lock so that multiple threads will not evolve the same sub population
		pop_index = *current_subpop_index;
		(*current_subpop_index)++;
		mutex->unlock();
#else
        pop_index = *current_subpop_index;
        (*current_subpop_index)++;
#endif
        
		// pop_index is the index of the subpopulation evolved by the current thread
		if (pop_index < params->num_sub_populations) {
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
}
//---------------------------------------------------------------------------
void start_steady_state(t_parameters &params, t_graph *training_graphs, int num_training_graphs, int num_variables, int num_procs, int current_proc_id)
{

	int size_to_send = params.code_length * sizeof(t_code3) + params.num_constants * 20 + 20;
	char *s_source = new char[size_to_send];
	char *s_dest = new char[size_to_send];

	t_chromosome receive_chromosome;
	allocate_chromosome(receive_chromosome, params);


	// a steady state model - 
	// Newly created inviduals replace the worst ones (if the offspring are better) in the same (sub) population.

	// allocate memory for all sub populations
	t_chromosome **sub_populations; // an array of sub populations
	sub_populations = new t_chromosome*[params.num_sub_populations];
	for (int p = 0; p < params.num_sub_populations; p++) {
		sub_populations[p] = new t_chromosome[params.sub_population_size];
		for (int i = 0; i < params.sub_population_size; i++)
			allocate_chromosome(sub_populations[p][i], params); // allocate each individual in the subpopulation 
	}
#ifdef USE_THREADS
	// allocate memory
	double** vars_values = new double*[params.num_threads];
	for (int t = 0; t < params.num_threads; t++)
		vars_values[t] = new double[num_variables];

	// an array of threads. Each sub population is evolved by a thread
	std::thread **mep_threads = new std::thread*[params.num_threads];
	// we create a fixed number of threads and each thread will take and evolve one subpopulation, then it will take another one
	std::mutex mutex;
	// we need a mutex to make sure that the same subpopulation will not be evolved twice by different threads
#else
    double* vars_values = new double[num_variables];
#endif
	
	// evolve for a fixed number of generations
	for (int generation = 0; generation < params.num_generations; generation++) { // for each generation

#ifdef USE_THREADS
		int current_subpop_index = 0;
		for (int t = 0; t < params.num_threads; t++)
			mep_threads[t] = new std::thread(evolve_one_subpopulation, &current_subpop_index, &mutex, sub_populations, generation, &params, training_graphs, num_training_graphs, num_variables, vars_values[t]);

		for (int t = 0; t < params.num_threads; t++) {
			mep_threads[t]->join();
			delete mep_threads[t];
		}
#else
       //for (int p = 0; p < params.num_sub_populations; p++)
        int p = 0;
               evolve_one_subpopulation(&p, sub_populations, generation, &params, training_graphs, num_training_graphs, num_variables, vars_values);
           
#endif
		// find the best individual
		int best_individual_subpop_index = 0; // the index of the subpopulation containing the best invidual
		
		for (int p = 1; p < params.num_sub_populations; p++) 
			if (sub_populations[p][0].fitness < sub_populations[best_individual_subpop_index][0].fitness)
				best_individual_subpop_index = p;

		printf("proc_id=%d, generation=%d, best=%lf\n", current_proc_id, generation, sub_populations[best_individual_subpop_index][0].fitness);

		if (current_proc_id == 0) {
			FILE* f = fopen("tst_log.txt", "a");
			fprintf(f, "proc_id=%d, generation=%d, best=%lf\n", current_proc_id, generation, sub_populations[best_individual_subpop_index][0].fitness);
			fclose(f);
		}

		// now copy one individual from one population to the next one.
		// the copied invidual will replace the worst in the next one (if is better)

		for (int p = 0; p < params.num_sub_populations; p++) {
			int  k = rand() % params.sub_population_size;// the individual to be copied
			// replace the worst in the next population (p + 1) - only if is better
			int index_next_pop = (p + 1) % params.num_sub_populations; // index of the next subpopulation (taken in circular order)
			if (sub_populations[p][k].fitness < sub_populations[index_next_pop][params.sub_population_size - 1].fitness) {
				copy_individual(sub_populations[index_next_pop][params.sub_population_size - 1], sub_populations[p][k], params);
				qsort((void *)sub_populations[index_next_pop], params.sub_population_size, sizeof(sub_populations[0][0]), sort_function);
			}
		}
        
#ifdef USE_MPI
		// here I have to copy few individuals from one process to another process
		for (int i = 0; i < 1; i++) {
			int source_sub_population_index = rand() % params.num_sub_populations;
			int chromosome_index = rand() % params.sub_population_size;
			int tag = 0;

			sub_populations[source_sub_population_index][chromosome_index].to_string(s_dest, params.code_length, params.num_constants);
			MPI_Send((void*)s_dest, size_to_send, MPI_CHAR, (current_proc_id + 1) % num_procs, tag, MPI_COMM_WORLD);

			MPI_Status status;
			int flag;
			MPI_Iprobe(!current_proc_id ? num_procs - 1 : current_proc_id - 1, tag, MPI_COMM_WORLD, &flag, &status);
			if (flag) {
				MPI_Recv(s_source, size_to_send, MPI_CHAR, !current_proc_id ? num_procs - 1 : current_proc_id - 1, tag, MPI_COMM_WORLD, &status);
				receive_chromosome.from_string(s_source, params.code_length, params.num_constants);

				int dest_sub_population_index = rand() % params.num_sub_populations;
				if (receive_chromosome.fitness < sub_populations[dest_sub_population_index][params.sub_population_size - 1].fitness) {
					copy_individual(sub_populations[dest_sub_population_index][params.sub_population_size - 1], receive_chromosome, params);
					qsort((void *)sub_populations[dest_sub_population_index], params.sub_population_size, sizeof(sub_populations[0][0]), sort_function);
				}
			}
		}
	

#endif
    }
    
    #ifdef USE_THREADS
	delete[] mep_threads;
#endif
		// print best t_chromosome
	// any of them can be printed because if we allow enough generations, the populations will become identical
	print_chromosome(sub_populations[0][0], params, num_variables);
	// free memory

	for (int p = 0; p < params.num_sub_populations; p++) {
		for (int i = 0; i < params.sub_population_size; i++)
			delete_chromosome(sub_populations[p][i]);
		delete[] sub_populations[p];
	}
	delete[] sub_populations;

    #ifdef USE_THREADS
	for (int t = 0; t < params.num_threads; t++)
		delete[] vars_values[t];
	delete[] vars_values;
#endif
    
	delete[] s_dest;
	delete[] s_source;

	delete_chromosome(receive_chromosome);
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
int main(int argc, char* argv[])
{
	t_parameters params;
	params.num_sub_populations = 1;
	params.sub_population_size = 50;						    // the number of individuals in population  (must be an even number!)
	params.code_length = 50;
	params.num_generations = 10000;					// the number of generations
	params.mutation_probability = 0.01;              // mutation probability
	params.crossover_probability = 0.9;             // crossover probability

	params.variables_probability = 0.4;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 10; // use 3 constants from -1 ... +1 interval
	params.constants_min = -1;
	params.constants_max = 1;

    #ifdef USE_THREADS
	params.num_threads = 4;
#endif
    
	t_graph *training_graphs = NULL;
	int num_training_graphs = 0;

	if (!read_training_data(training_graphs, num_training_graphs)) {
		printf("Cannot find input file(s)! Please specify the full path!\n");
        printf("Press Enter ...");
		getchar();
		return 1;
	}

	compute_global_variables(training_graphs, num_training_graphs);

	int num_variables = 10;

	int current_proc_id = 0;

	int num_procs = 0;
#ifdef USE_MPI

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &current_proc_id);
#endif

	srand(current_proc_id); // we run each process with a different seed


	printf("Evolving... proc_id = %d\n", current_proc_id);
	start_steady_state(params, training_graphs, num_training_graphs, num_variables, num_procs, current_proc_id);

	delete_training_graphs(training_graphs, num_training_graphs);

#ifdef USE_MPI
	if (current_proc_id == 0) {


		printf("Press Enter ...");
		getchar();
	}

	MPI_Finalize();
#endif

	return 0;
}
//--------------------------------------------------------------------