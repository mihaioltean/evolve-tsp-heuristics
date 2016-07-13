//---------------------------------------------------------------------------
//   Multi Expression Programming Software - with multiple subpopulations and threads and MPI
//   Copyright Mihai Oltean  (mihai.oltean@gmail.com), Virginia Niculescu (vniculescu@cs.ubbcluj.ro)
//   Version 2016.06.26.0 // year.month.day.build#

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
//   just create an empty console project and add this file to the project

//   More info at:  http://www.mep.cs.ubbcluj.ro

//   Please reports any sugestions and/or bugs to:     mihai.oltean@gmail.com or to vniculescu@cs.ubbcluj.ro


//--------------------------------------------------------------------

//#define USE_THREADS
#define USE_MPI

#define SIMPLIFY

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <time.h>

#ifdef USE_THREADS
#include <thread>
#endif

#include <float.h>

//#include <vld.h> // for detecting memory leaks in VC++

#ifdef USE_MPI
#include "mpi.h"
#endif

#define max_file_name 50

#define num_operators 7

// +   -1
// -   -2
// *   -3
// min   -4
// max   -5

#define O_ADDITION -1
#define O_SUBTRACTION -2
#define O_MULTIPLICATION -3
#define O_MIN -4
#define O_MAX -5
#define O_SIN -6
#define O_COS -7

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
// problem parameters
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
	int num_migrations ;
	int num_training_graphs;
};
//---------------------------------------------------------------------------
//execution parameters
//---------------------------------------------------------------------------
struct run_parameters{
	char run_config_file[max_file_name];
	char config_file[max_file_name];
	char log_file[max_file_name];
	char result_file[max_file_name];
	int num_threads; // num threads.
	int num_procs; // num processes
	int current_id; // rank
	int recv_no; // no of receivings  for one process
};
//---------------------------------------------------------------------------
// related to chromosomes
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
	t_code3 *simplified_prg;
	double *constants; // an array of constants

	double fitness;        // the fitness (or the error)
	int num_utilized_instructions; // num_utilized_instructions by the simplified program

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
	void to_string_simplified(char * s_dest, int num_constants)
	{
		char tmp_s[100];
		s_dest[0] = 0;
		sprintf(tmp_s, "%d ", num_utilized_instructions);
		strcat(s_dest, tmp_s);
		for (int i = 0; i < num_utilized_instructions; i++) {
			sprintf(tmp_s, "%d ", simplified_prg[i].op);
			strcat(s_dest, tmp_s);
			sprintf(tmp_s, "%d ", simplified_prg[i].adr1);
			strcat(s_dest, tmp_s);
			sprintf(tmp_s, "%d ", simplified_prg[i].adr2);
			strcat(s_dest, tmp_s);
		}
		sprintf(tmp_s, "%d ", num_constants);
		strcat(s_dest, tmp_s);
		for (int i = 0; i < num_constants; i++) {
			sprintf(tmp_s, "%lg ", constants[i]);
			strcat(s_dest, tmp_s);
		}
		sprintf(tmp_s, "%lg ", fitness);
		strcat(s_dest, tmp_s);
	}
	//------------------------------------------------------------------------------
	bool to_file_simplified(char* file_name, int code_length, int num_constants)
	{
		FILE *f = fopen(file_name, "w");
		if (!f)
			return false;

		simplify(code_length);

		char s[1000];

		to_string_simplified(s, num_constants);
		fputs(s, f);

		fclose(f);
		return true;
	}
	//------------------------------------------------------------------------------
	void mark(int k, bool* marked)
	// mark all utilized instructions
	{
		if ((prg[k].op < 0) && !marked[k]) {
			mark(prg[k].adr1, marked);

			switch (prg[k].op) {
				case O_ADDITION:
					mark(prg[k].adr2, marked);
					break;
				case O_SUBTRACTION:
					mark(prg[k].adr2, marked);
					break;
				case O_MULTIPLICATION:
					mark(prg[k].adr2, marked);
					break;
				case O_MIN:
					mark(prg[k].adr2, marked);
					break;
				case O_MAX:
					mark(prg[k].adr2, marked);
					break;
					// sin, cos are already marked because they have only 1 parameter
			}
		}
		marked[k] = true;
	}
	//---------------------------------------------------------------------------
	void simplify(int code_length)
	{
#ifdef SIMPLIFY
		bool *marked = new bool[code_length];
		for (int i = 0; i < code_length; marked[i++] = false);
		mark(code_length - 1, marked);

		// how many are skipped until a given instruction
		int *skipped = new int[code_length];
		if (!marked[0])
			skipped[0] = 1;
		else
			skipped[0] = 0;
		for (int i = 1; i < code_length; i++)
			if (!marked[i])
				skipped[i] = skipped[i - 1] + 1;
			else
				skipped[i] = skipped[i - 1];

		if (simplified_prg)
			delete[] simplified_prg;
		simplified_prg = new t_code3[code_length];

		num_utilized_instructions = 0;
		for (int i = 0; i < code_length; i++)
			if (marked[i]) {
				simplified_prg[num_utilized_instructions] = prg[i];
				if (prg[i].op < 0) {// operator
					simplified_prg[num_utilized_instructions].adr1 -= skipped[prg[i].adr1];
					simplified_prg[num_utilized_instructions].adr2 -= skipped[prg[i].adr2];
				}
				num_utilized_instructions++;
			}

		delete[] skipped;
		delete[] marked;
#else
		if (simplified_prg)
			delete[] simplified_prg;
		simplified_prg = new t_code3[code_length];
		num_utilized_instructions = code_length;
		for (int i = 0; i < code_length; i++)
			simplified_prg[i] = prg[i];
#endif
	}
	//---------------------------------------------------------------------------

};
//---------------------------------------------------------------------------
void allocate_chromosome(t_chromosome &c, t_parameters &params)
{
	c.prg = new t_code3[params.code_length];
	c.simplified_prg = NULL;
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
	if (c.simplified_prg) {
		delete[] c.simplified_prg;
		c.simplified_prg = NULL;
	}

	if (c.constants) {
		delete[] c.constants;
		c.constants = NULL;
	}
}
//---------------------------------------------------------------------------
void copy_individual(t_chromosome& dest, const t_chromosome& source, t_parameters &params)
{
	for (int i = 0; i < params.code_length; i++)
		dest.prg[i] = source.prg[i];
	for (int i = 0; i < params.num_constants; i++)
		dest.constants[i] = source.constants[i];
	dest.fitness = source.fitness;
	dest.num_utilized_instructions = source.num_utilized_instructions;
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

	a.simplify(params.code_length);
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
	a_chromosome.simplify(params.code_length);
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
	offspring1.simplify(params.code_length);
	offspring2.simplify(params.code_length);
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
	offspring1.simplify(params.code_length);
	offspring2.simplify(params.code_length);

}
//---------------------------------------------------------------------------
int comp_function(const void *a, const void *b)
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
double evaluate(t_code3 *prg, int head_index, int num_variables, double *vars_values, double *partial_values_array, double *constants)
{
	for (int i = 0; i <= head_index; i++)
		switch (prg[i].op) {
			case O_ADDITION:// +
				partial_values_array[i] = partial_values_array[prg[i].adr1] + partial_values_array[prg[i].adr2];
				break;
			case O_SUBTRACTION:// -
				partial_values_array[i] = partial_values_array[prg[i].adr1] - partial_values_array[prg[i].adr2];
				break;
			case O_MULTIPLICATION:// *
				partial_values_array[i] = partial_values_array[prg[i].adr1] * partial_values_array[prg[i].adr2];
				break;
			case O_MIN:// min
				partial_values_array[i] = partial_values_array[prg[i].adr1] < partial_values_array[prg[i].adr2] ? partial_values_array[prg[i].adr1] : partial_values_array[prg[i].adr2];
				break;
			case O_MAX:// max
				partial_values_array[i] = partial_values_array[prg[i].adr1] > partial_values_array[prg[i].adr2] ? partial_values_array[prg[i].adr1] : partial_values_array[prg[i].adr2];
				break;
			case O_SIN:// sin
				partial_values_array[i] = sin(partial_values_array[prg[i].adr1]);
				break;
			case O_COS:// cos
				partial_values_array[i] = cos(partial_values_array[prg[i].adr1]);
				break;
			default:
				if (prg[i].op < num_variables)
					partial_values_array[i] = vars_values[prg[i].op];
				else
					partial_values_array[i] = constants[prg[i].op - num_variables];
		}

	return partial_values_array[head_index]; // last gene is the one providing the output
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
bool file_print_chromosome(t_chromosome& a, t_parameters &params, int num_variables, char* filename)
{

	FILE* f = fopen(filename, "a");
	if (!f)
		return false;

	fprintf(f, "The t_chromosome is:\n");

	for (int i = 0; i < params.num_constants; i++)
		fprintf(f,"constants[%d] = %lf\n", i, a.constants[i]);

	for (int i = 0; i < params.code_length; i++)
		if (a.prg[i].op < 0)
			fprintf(f,"%d: %s %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].adr1, a.prg[i].adr2);
		else
		if (a.prg[i].op < num_variables)
			fprintf(f,"%d: inputs[%d]\n", i, a.prg[i].op);
		else
			fprintf(f,"%d: constants[%d]\n", i, a.prg[i].op - num_variables);

	fprintf(f,"Fitness = %lf\n", a.fitness);

	fprintf(f, "simplified = ");
	a.simplify(params.code_length);

	char s[1000];

	a.to_string_simplified(s, params.num_constants);
	fputs(s, f);

	fclose(f);
	return true;
}
//--------------------------------------------------------------------
//related to training graphs
//--------------------------------------------------------------------

struct t_graph{
	int num_nodes;			// number of nodes
	double **distance;      // distance matrix
	double min_distance, max_distance, average_distance;
	double optimal_length;  // known optimal solution
};
//---------------------------------------------------------------------------
bool allocate_training_graphs(t_graph *&training_graphs, int num_training_graphs){
	training_graphs = new t_graph[num_training_graphs];

	return true;
}
//---------------------------------------------------------------------------

int read_training_data(t_graph *training_graphs)
{
	// training is done on 4 graphs
	//4 graphs are read
	int k = 0; // count the graphs
	{
		FILE* f = fopen("data//bayg29.tsp", "r");
		if (!f)
			return 0;

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
			return 0;

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
			return 0;

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
			return 0;

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
	}	return 1;
}
//---------------------------------------------------------------------------
void delete_training_graphs(t_graph *&training_graphs, int num_training_graphs)
{
	if (training_graphs) {
		for (int i = 0; i < num_training_graphs; i++) {
			for (int j = 0; j < training_graphs[i].num_nodes; j++)
				delete[] training_graphs[i].distance[j];
			delete[] training_graphs[i].distance;
		}
		delete[] training_graphs;
	}
}

//--------------------------------------------------------------------
// here we compute the min, max and average distance in a given graph
//--------------------------------------------------------------------
void compute_global_variables(t_graph *training_graphs, int num_training_graphs) {

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
		training_graphs[k].average_distance /=
				(training_graphs[k].num_nodes * training_graphs[k].num_nodes) - training_graphs[k].num_nodes;
	}
}

//---------------------------------------------------------------------------
void calibrate_vars_values(double * vars_values, t_graph *training_graphs, int k){
	vars_values[min_distance_in_matrix]     /=  training_graphs[k].max_distance; //in order to reduce to [0,1] range
	vars_values[max_distance_in_matrix]     /=  training_graphs[k].max_distance; //in order to reduce to [0,1] range
	vars_values[average_distance_in_matrix] /=  training_graphs[k].max_distance;//in order to reduce to [0,1] range

	vars_values[num_total_nodes]            /=  training_graphs[k].num_nodes ; //in order to reduce to [0,1] range

	vars_values[length_so_far]              /=  training_graphs[k].max_distance;
	vars_values[num_visited_nodes]          /=  training_graphs[k].num_nodes;
	vars_values[distance_to_next_node]      /=  training_graphs[k].max_distance;

}
//---------------------------------------------------------------------------

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
//---------------------------------------------------------------------------

void fitness(t_chromosome &individual, t_graph *training_graphs, int num_training_graphs,
			 int num_variables, double * vars_values, double *partial_values_array)

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

					calibrate_vars_values(vars_values,training_graphs, k);

					double eval = evaluate(individual.simplified_prg, individual.num_utilized_instructions - 1, num_variables, vars_values, partial_values_array, individual.constants);
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
// this function evolves in one execution thread a number o sub_population
// inside the range defined by [ start, stop)
//-----------------------------------------------------------------
void evolve_subpopulation_range( int generation_index, t_chromosome ** sub_populations,
								 t_parameters *params, t_graph *training_graphs,
								 int num_variables, double* vars_values, int start, int stop){

	for(int pop_index=start; pop_index< stop;pop_index++){
		t_chromosome *a_sub_population = sub_populations[pop_index];
		//printf("inside evolve range index= %d\n", pop_index);
		t_chromosome offspring1, offspring2;

		allocate_chromosome(offspring1, *params);
		allocate_chromosome(offspring2, *params);
		double *partial_values_array = new double[params->code_length];

		if (generation_index == 0) {
			for (int i = 0; i < params->sub_population_size; i++) {
				generate_random_chromosome(a_sub_population[i], *params, num_variables);
				fitness(a_sub_population[i], training_graphs, params->num_training_graphs, num_variables,
						vars_values, partial_values_array);
			}
			// sort ascendingly by fitness inside this population
			qsort((void *) a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]), comp_function);
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
				fitness(offspring1, training_graphs, params->num_training_graphs, num_variables,
						vars_values, partial_values_array);

				// mutate the other offspring too
				mutation(offspring2, *params, num_variables);
				fitness(offspring2, training_graphs, params->num_training_graphs, num_variables,
						vars_values, partial_values_array);

				// replace the worst in the population
				if (offspring1.fitness < a_sub_population[params->sub_population_size - 1].fitness) {
					copy_individual(a_sub_population[params->sub_population_size - 1], offspring1, *params);
					qsort((void *) a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]),
						  comp_function);
				}
				if (offspring2.fitness < a_sub_population[params->sub_population_size - 1].fitness) {
					copy_individual(a_sub_population[params->sub_population_size - 1], offspring2, *params);
					qsort((void *) a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]),
						  comp_function);
				}
			}
		delete_chromosome(offspring1);
		delete_chromosome(offspring2);

	}
}
//-----------------------------------------------------------------
// all the sub_populations are evolved
// this function that works either in case of using threads or in case when threads are not used
//-----------------------------------------------------------------
void evolve_subpopulations( run_parameters & r_params, t_chromosome ** sub_populations, int generation_index,
							t_parameters *params, t_graph *training_graphs,
							int num_variables, double* vars_values) {

#ifdef USE_THREADS
	//printf("debug proc %d  before starting generation %d after threads array alloc\n", r_params.current_id, generation_index);
	std::thread **t = new std::thread*[r_params.num_threads];

	int start = 0, stop, step, rest;


	step = params->num_sub_populations / r_params.num_threads;
	rest = params->num_sub_populations % r_params.num_threads;
	stop = (rest > 0) ? start + step + 1 : start + step;
	//printf("debug proc %d  before starting generation %d\n", r_params.current_id, generation_index);

	for (int i = 0; i < r_params.num_threads; i++) {
		t[i] = new std::thread(
				evolve_subpopulation_range,
				generation_index, sub_populations, params, training_graphs,
				num_variables, vars_values, start, stop
		);

		rest--;
		start = stop;
		stop = (rest > 0) ? start + step + 1 : start + step;
		if (stop>params->num_sub_populations) stop = params->num_sub_populations;
	}
	for (int i = 0; i < r_params.num_threads; i++) {
		t[i]->join();
	}

//	printf("debug proc %d t generation %d final\n", r_params.current_id, generation_index);
	for (int i = 0; i < r_params.num_threads; i++) delete t[i];
	delete[] t;
//	printf("debug proc %d t generation %d final after deleting threads\n", r_params.current_id, generation_index);
#else
	evolve_subpopulation_range( generation_index, sub_populations, params, training_graphs,
								num_variables, vars_values, 0, params->num_sub_populations);
#endif
}

//--------------------------------------------------------------------
// MPI interchange
// this defines the interchange between subpopulations on different MPI processes
//--------------------------------------------------------------------
void interchange(run_parameters& r_params, t_parameters&  params,t_chromosome **sub_populations, t_chromosome & receive_chromosome){

#ifdef USE_MPI
	MPI_Request send_request= MPI_REQUEST_NULL;
	MPI_Request recv_request= MPI_REQUEST_NULL;
	MPI_Status status;

	int size_to_send = params.code_length * sizeof(t_code3) + params.num_constants * 20 + 20;
	char *s_source = new char[size_to_send];
	char *s_dest = new char[size_to_send];
	for (int i = 0; i < params.num_migrations; i++) {
		int source_sub_population_index = rand() % params.num_sub_populations;
		int chromosome_index = rand() % params.sub_population_size;
		int tag = 0;

		/*if(generation>0) {
            if (send_request != MPI_REQUEST_NULL) MPI_Cancel(&send_request);//MPI_Request_free(&send_request);//
            //MPI_Wait(&send_request, &status);
        }*/
		sub_populations[source_sub_population_index][chromosome_index].to_string(s_source, params.code_length, params.num_constants);

		//MPI_Isend(s_source, size_to_send, MPI_CHAR, (r_params.current_id + 1) % r_params.num_procs, tag, MPI_COMM_WORLD,&send_request);
		MPI_Send((void*)s_source, size_to_send, MPI_CHAR, (r_params.current_id + 1) % r_params.num_procs, tag, MPI_COMM_WORLD);


		//int flag;
		//MPI_Iprobe(!r_params.current_id ? r_params.num_procs - 1 : r_params.current_id - 1, tag, MPI_COMM_WORLD, &flag, &status);
		//if (flag)
		{
			MPI_Recv(s_dest, size_to_send, MPI_CHAR, !r_params.current_id? r_params.num_procs - 1 : r_params.current_id - 1, tag, MPI_COMM_WORLD, &status);

			//if (recv_request != MPI_REQUEST_NULL) MPI_Request_free(&recv_request);
			//MPI_Irecv(s_dest, size_to_send, MPI_CHAR,  !current_proc_id ? num_procs - 1 :  current_proc_id - 1,
			//		  tag, MPI_COMM_WORLD,  &recv_request);

			//MPI_Test(&recv_request, &flag, &status);
			//if (flag)
			{
				r_params.recv_no++;
				receive_chromosome.from_string(s_dest, params.code_length, params.num_constants);
				int dest_sub_population_index = rand() % params.num_sub_populations;
				if (receive_chromosome.fitness < sub_populations[dest_sub_population_index][params.sub_population_size - 1].fitness) {
					copy_individual(sub_populations[dest_sub_population_index][params.sub_population_size - 1],
									receive_chromosome, params);
					qsort((void *) sub_populations[dest_sub_population_index], params.sub_population_size,
						  sizeof(sub_populations[0][0]), comp_function);

				}
			}
		}
	}
#endif
}
//---------------------------------------------------------------------------
//print to log file if USE_MPI
//---------------------------------------------------------------------------
int print_to_log_file(t_parameters &params, run_parameters& r_params, int generation, t_chromosome **sub_populations, int best_individual_subpop_index) {
/*
 * the function return an error code
 * ->0 means no errors
 */

    #ifdef USE_MPI
	char message[1040];
	char expression[1000];

	sub_populations[best_individual_subpop_index][0].simplify(params.code_length);

	sub_populations[best_individual_subpop_index][0].to_string_simplified(expression, params.num_constants);

	sprintf(message, "pid= %d gen %d best_fitness=%f simplified = %s\n", r_params.current_id, generation,
			sub_populations[best_individual_subpop_index][0].fitness, expression );
	MPI_File file;
	MPI_Status status;

	int err;
	//while (err)
	{
		err = MPI_File_open(MPI_COMM_SELF,
							r_params.log_file,
							MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &file);


	}
	//int err=MPI_File_seek(*file, MPI_SEEK_END,0);

	if (err) {
		printf("try to open the file %s in proc %d with error %d\n", r_params.log_file, r_params.current_id, err);
		MPI_Abort(MPI_COMM_WORLD, err);
	}
	else {
		err = MPI_File_write(file, message, (int) strlen(message), MPI_CHAR, &status);

		if (err)
			printf("write on file %s with error= %d in proc %d message %s of length %d\n", r_params.log_file, err,
				   r_params.current_id, message, (int) strlen(message));
	}
	MPI_File_close(&file);

	return  err;
#endif
}
//---------------------------------------------------------------------------
// the main evolution function
//	-  a steady state model -
// Newly created inviduals replace the worst ones (if the offspring are better) in the same (sub) population.
//---------------------------------------------------------------------------

void start_steady_state(t_parameters &params, t_graph *training_graphs,	int num_variables, run_parameters &r_params)
{
	t_chromosome receive_chromosome;
	allocate_chromosome(receive_chromosome, params);

	//populations creation

	// allocate memory for all sub populations
	t_chromosome **sub_populations; // an array of sub populations

	sub_populations = new t_chromosome*[params.num_sub_populations];
	for (int p = 0; p < params.num_sub_populations; p++) {
		sub_populations[p] = new t_chromosome[params.sub_population_size];
		for (int i = 0; i < params.sub_population_size; i++)
			allocate_chromosome(sub_populations[p][i], params); // allocate each individual in the subpopulation
	}

	// allocate memory for variables
	double* vars_values = new double[num_variables];

	// evolve for a fixed number of generations
	for (int generation = 0; generation < params.num_generations; generation++) { // for each generation
		//int p = 0;
		//for(int p = 0;p<params.sub_population_size; p++)
		evolve_subpopulations( r_params, sub_populations, generation, &params, training_graphs,
							   num_variables, vars_values);

		// find the best individual
		int best_individual_subpop_index = 0; // the index of the subpopulation containing the best invidual
		for (int p = 1; p < params.num_sub_populations; p++)
			if (sub_populations[p][0].fitness < sub_populations[best_individual_subpop_index][0].fitness)
				best_individual_subpop_index = p;

//print to log file the intermediary results
#ifndef USE_MPI
		FILE* f = fopen("tst_log.txt", "a");
			fprintf(f, "generation=%d, best=%lf\n", generation, sub_populations[best_individual_subpop_index][0].fitness);
			fclose(f);
#else
		if (generation%10==0) {
			int err = print_to_log_file(params, r_params, generation, sub_populations, best_individual_subpop_index);
			if (err) printf("error while writing to log file inside process %d \n", r_params.current_id);
		}
#endif

// now copy an individual from one population to the next one.
// the copied invidual will replace the worst in the next one (if is better)
		for (int p = 0; p < params.num_sub_populations; p++) {
			int  k = rand() % params.sub_population_size;// the individual to be copied
			// replace the worst in the next population (p + 1) - only if is better
			int index_next_pop = (p + 1) % params.num_sub_populations; // index of the next subpopulation (taken in circular order)
			if (sub_populations[p][k].fitness < sub_populations[index_next_pop][params.sub_population_size - 1].fitness) {

				copy_individual(sub_populations[index_next_pop][params.sub_population_size - 1], sub_populations[p][k], params);

				qsort((void *) sub_populations[index_next_pop], params.sub_population_size,
					  sizeof(sub_populations[0][0]), comp_function);
			}
		}

#ifdef USE_MPI
		// here we have to copy few individuals from one process to the next process
///////////////////////////////////////////////////////////////////////////////////////
		interchange(r_params, params,sub_populations, receive_chromosome );
		///////////////////////////////////////////////////////////////////////////////////////
#endif
	}
	// print best t_chromosome
	// any of them can be printed because if we allow enough generations, the populations will become identical
	//print_chromosome(sub_populations[0][0], params, num_variables);
#ifdef USE_MPI
	if (r_params.current_id==0)
#endif
		file_print_chromosome(sub_populations[0][0], params, num_variables,r_params.result_file);

	// free memory - for allocated populations
	for (int p = 0; p < params.num_sub_populations; p++) {
		for (int i = 0; i < params.sub_population_size; i++)
			delete_chromosome(sub_populations[p][i]);
		delete[] sub_populations[p];
	}
	delete[] sub_populations;
	delete_chromosome(receive_chromosome);
}


//---------------------------------------------------------------------------
//print to log file if USE_MPI the no of elements thathas been received from previous process
//---------------------------------------------------------------------------
void print_num_migrations(run_parameters &r_params){
#ifdef USE_MPI
	char * message;
	message = new char[40];

	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(message,"in proc id = %d recv_no=%d\n", r_params.current_id, r_params.recv_no);

	//sprintf(message, "pid= %d gen %d best_fitness=%f\n",  r_params.current_id,  generation,
	//		sub_populations[best_individual_subpop_index][0].fitness);
	MPI_File file;
	MPI_Status status;
	int err = MPI_File_open(MPI_COMM_SELF,
							r_params.log_file,
							MPI_MODE_CREATE| MPI_MODE_WRONLY|MPI_MODE_APPEND,	MPI_INFO_NULL, &file);


	if (err) printf("try to open file %s in proc %d error %d\n",r_params.log_file,r_params.current_id ,err);
	else {
		err =   MPI_File_write(file, message, (int) strlen(message), MPI_CHAR, &status);
		if(err)
		printf(" tried in proc %d print to file %s the message %s with status(no of char written) = %lu\n", r_params.current_id, r_params.log_file, message, status._ucount);
	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_File_sync(file);

	//printf("write on file in proc %d message %s of length %d error %d\n",r_params.current_id ,message, (int)strlen(message),err);

	MPI_File_close(&file);
	//MPI_Barrier(MPI_COMM_WORLD);
	//MPI_File_sync(file);
#endif
}

//--------------------------------------------------------------------
void init_run_params(run_parameters& r_params) {
	strcpy(r_params.log_file, "log");
	strcpy(r_params.result_file, "result");

	r_params.current_id = 0;
	r_params.num_procs = 1;
	r_params.num_threads = 1;// it is recommended to be a number that divide params.num_sub_populations
	r_params.recv_no = 0;
}
//--------------------------------------------------------------------
void init_params(t_parameters& params) {
	params.num_sub_populations = 8;
	params.sub_population_size = 30;                            // the number of individuals in population  (must be an even number!)
	params.code_length = 50;
	params.num_generations = 5;                    // the number of generations
	params.mutation_probability = 0.1;              // mutation probability
	params.crossover_probability = 0.9;             // crossover probability

	params.variables_probability = 0.4;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability;
	// / sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 3; // use 3 constants from -1 ... +1 interval
	params.constants_min = -1;
	params.constants_max = 1;
	params.num_migrations =2;
	params.num_training_graphs = 4;
}


//--------------------------------------------------------------------
// set the names of the output files depending on the params
void set_name_files(t_parameters& t_params, run_parameters &r_params){

	char np[20];

#ifndef USE_MPI

	#ifdef USE_THREADS
	sprintf(np, "_g%d_t%d.txt", t_params.num_generations, r_params.num_threads);
#else
	sprintf(np,"_g%d.txt",t_params.num_generations);
#endif

#else

#ifdef USE_THREADS
	sprintf(np, "_g%d_t%d_p%d.txt", t_params.num_generations, r_params.num_threads,r_params.num_procs);
#else
	sprintf(np,"_g%d_p%d.txt",t_params.num_generations,r_params.num_procs);
#endif
#endif
	strcat(r_params.result_file, np);
	strcat(r_params.log_file, np);
}
//--------------------------------------------------------------------
//init the evolve parameters from a config file (the file name is in r_params.config_file
//--------------------------------------------------------------------
int init_params_config_file(run_parameters & r_params, t_parameters& params) {

	FILE* f = fopen(r_params.config_file,"r");
	if (!f){
		return 1; //error - config file does not exist r couldn't be opened
	}

	fscanf(f, "%d",&params.num_sub_populations );
	if (params.num_sub_populations<1) params.num_sub_populations = 10;

	fscanf(f, "%d",&params.sub_population_size );     // the number of individuals in population  (must be an even number!)
	if (params.sub_population_size<1) params.sub_population_size = 10;

	fscanf(f, "%d",&params.code_length );
	if (params.code_length<1) params.code_length = 10;

	fscanf(f, "%d",&params.num_generations);// the number of generations
	if(params.num_generations<1) params.num_generations=10;

	fscanf(f, "%lf",&params.mutation_probability );   // mutation probability
	if (params.mutation_probability<0 || params.mutation_probability>1) params.mutation_probability = 0.1;

	fscanf(f, "%lf",&params.crossover_probability );   // crossover probability
	if (params.crossover_probability<0 || params.crossover_probability>1) params.crossover_probability = 0.9;

	fscanf(f, "%lf",&params.variables_probability  );   // crossover probability
	if (params.variables_probability<0 || params.variables_probability >1) params.variables_probability  = 0.4;

	fscanf(f, "%lf",&params.operators_probability );   // crossover probability
	if (params.operators_probability<0 || params.operators_probability>1) params.operators_probability = 0.5;

	// / sum of variables_prob + operators_prob + constants_prob MUST BE 1 !
	params.constants_probability = 1 - params.variables_probability - params.operators_probability;


	fscanf(f, "%d",&params.num_constants );   // use a no. of constants from (constants_min...constants_max )interval
	if( params.num_constants<1)
		params.num_constants = 3;
	fscanf(f, "%lf",&params.constants_min );
	fscanf(f, "%lf",&params.constants_max);
	if( params.constants_min >=	params.constants_max){
		params.constants_min = -1;
		params.constants_max = 1;
	}

	fscanf(f, "%d",&params.num_migrations );   // use a no. of individuals that are moved between super-populations
	if( params.num_migrations<1) params.num_migrations = 1;


	fclose(f);
	params.num_training_graphs = 4;

	return 0;
}
//--------------------------------------------------------------------
//init the run parameters from a run_config file
//--------------------------------------------------------------------
int init_run_params_config_file(run_parameters& r_params) {

	r_params.recv_no = 0;

	FILE* f = fopen(r_params.run_config_file,"r");
	if (!f){
		return 1; //error - config file does not exist or couldn't be opened
	}
	fscanf(f, "%s",r_params.log_file );         // the name of the log file
	fscanf(f, "%s",r_params.result_file );      // the name of the result file

#ifdef USE_THREADS
	fscanf(f, "%d",&r_params.num_threads);
	if (r_params.num_threads == 0) r_params.num_threads = 1;
#endif

	fclose(f);

	return 0;
}
//--------------------------------------------------------------------
//init the names of the configuration files
//--------------------------------------------------------------------
int init_run_params_command_line(run_parameters& r_params, int argc, char* argv[]) {

	r_params.recv_no = 0;

//	take the names of the config file
	if (argc < 2)
		strcpy(r_params.run_config_file, "run_config.txt");
	else
	{
		strcpy(r_params.run_config_file, argv[1]);

	}
	if (argc < 3)
		strcpy(r_params.config_file,"config.txt");
	else{
		strcpy(r_params.config_file, argv[2]);

	}
	return 0;
}

//--------------------------------------------------------------------
// write the parameters on the result file and on the log file
//--------------------------------------------------------------------
int init_files(t_parameters & t_params, run_parameters &r_params){
	if (r_params.current_id==0)
	{
		//init result file
		FILE* f = fopen(r_params.result_file,"a");

		if (!f) return 1;

		fprintf(f,"_____________________________________________________________\n");
		fprintf(f,"_____________________________________________________________\n");
		fprintf(f,"num_sub_populations %d \n", t_params.num_sub_populations );
		fprintf(f,"sub_population_size %d \n", t_params.sub_population_size );     // the number of individuals in population  (must be an even number!)
		fprintf(f,"code_length %d \n", t_params.code_length );
		fprintf(f,"num_generations %d\n", t_params.num_generations );                    // the number of generations
		fprintf(f,"mutation_probability %f \n",t_params.mutation_probability );              // mutation probability
		fprintf(f,"crossover_probability %f \n",t_params.crossover_probability );             // crossover probability

		fprintf(f,"variables_probability %f\n", t_params.variables_probability );
		fprintf(f,"operators_probability %f\n", t_params.operators_probability );
		fprintf(f,"constants_probability %f \n", t_params.constants_probability );
		// / sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

		fprintf(f,"num_constants %d \n", t_params.num_constants); // use 3 constants from -1 ... +1 interval
		fprintf(f,"constants_min %f\n",t_params.constants_min );
		fprintf(f,"constants_max %f \n", t_params.constants_max );

		fprintf(f,"num_training_graphs = %d\n", t_params.num_training_graphs);
		fprintf(f,"_______________________________________________________________\n");

		fclose(f);

//init log file
		f = fopen(r_params.log_file,"a");

		if (!f) return 2;

		fprintf(f,"_____________________________________________________________\n");
		fprintf(f,"_____________________________________________________________\n");
		fprintf(f,"num_sub_populations %d \n", t_params.num_sub_populations );
		fprintf(f,"sub_population_size %d \n", t_params.sub_population_size );     // the number of individuals in population  (must be an even number!)
		fprintf(f,"code_length %d \n", t_params.code_length );
		fprintf(f,"num_generations %d\n", t_params.num_generations );                    // the number of generations
		fprintf(f,"mutation_probability %f \n",t_params.mutation_probability );              // mutation probability
		fprintf(f,"crossover_probability %f \n",t_params.crossover_probability );             // crossover probability

		fprintf(f,"variables_probability %f\n", t_params.variables_probability );
		fprintf(f,"operators_probability %f\n", t_params.operators_probability );
		fprintf(f,"constants_probability %f \n", t_params.constants_probability );
		// / sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

		fprintf(f,"num_constants %d \n", t_params.num_constants); // use 3 constants from -1 ... +1 interval
		fprintf(f,"constants_min %f\n",t_params.constants_min );
		fprintf(f,"constants_max %f \n", t_params.constants_max );

		fprintf(f,"num_training_graphs = %d\n", t_params.num_training_graphs);
		fprintf(f,"___________________________________________________\n");
		fclose(f);



	}
	return 0;
}


//--------------------------------------------------------------------
// main function for non MPI case
//--------------------------------------------------------------------
int  main_no_mpi(t_parameters& t_params,run_parameters& r_params,t_graph *training_graphs,int argc, char* argv[]){

	init_run_params_command_line(r_params,  argc, argv);

	if  (!init_run_params_config_file(r_params) )
		init_run_params(r_params);

	if (!init_params_config_file( r_params, t_params))
		init_params(t_params);

	set_name_files(t_params, r_params);

	if (allocate_training_graphs(training_graphs, t_params.num_training_graphs)) {
		int read_sum = read_training_data(training_graphs);
		if (read_sum < 1) {
			printf("Cannot find input file(s)! Please specify the full path!\n");
			printf("Press Enter ...");
			getchar();
			return 1;
		}

		compute_global_variables(training_graphs, t_params.num_training_graphs);
		int num_variables = 10;

		//srand(current_proc_id); // we run each process with a different seed

		printf("Evolving. ..\n");

		start_steady_state(t_params, training_graphs, num_variables, r_params);
		delete_training_graphs(training_graphs, t_params.num_training_graphs);

		printf("Press enter ...");
		getchar();

		return 0;
	}
	else
		return 1;

}
//--------------------------------------------------------------------
//--------------------------------------------------------------------
int main(int argc, char* argv[])
{

	t_parameters t_params;
	run_parameters r_params;

	t_graph *training_graphs = NULL;

#ifndef USE_MPI
	return main_no_mpi( t_params, r_params, training_graphs,  argc, argv);
#else

	double starttime, endtime, computetime, compute_sum=0, compute_max=0;


#ifndef  USE_THREADS

	MPI_Init(&argc, &argv);
#else


	int provided, flag,claimed, errs;
		MPI_Init_thread( &argc, &argv,MPI_THREAD_FUNNELED, &provided );

		MPI_Is_thread_main( &flag );
		if (!flag) {
			errs++;
			printf( "This thread called init_thread but Is_thread_main gave false\n" );fflush(stdout);
		}
		else  printf( "This thread called init_thread with thread level %d\n", provided );fflush(stdout);
		MPI_Query_thread( &claimed );
		if (claimed != provided) {
			errs++;
			printf( "Query thread gave thread level %d but Init_thread gave %d\n", claimed, provided );fflush(stdout);
		}
		else ;//printf( "Query thread gave thread level %d but Init_thread gave %d\n", claimed, provided );fflush(stdout);
#endif


	MPI_Comm_size(MPI_COMM_WORLD, &r_params.num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &r_params.current_id);


	init_run_params_command_line(r_params,   argc, argv);

	if  (init_run_params_config_file(r_params) )
		init_run_params(r_params);

	if (init_params_config_file(r_params, t_params))
		init_params(t_params);

//adjust the output filenames with params
	set_name_files(t_params, r_params);


// write the parameters values on  the result file
	if (init_files(t_params, r_params)) {
		printf("result and log files could not be opend");
		return 1;
	};


//all processes will read the graphs
	bool alloc_signal = allocate_training_graphs(training_graphs, t_params.num_training_graphs);
	if (alloc_signal) {

		//reading graphs
		int read_sum = 0, verif_read_sum = 0;
		read_sum = read_training_data(training_graphs);
		MPI_Allreduce(&read_sum, &verif_read_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if (verif_read_sum == r_params.num_procs) {

			//an attempt to improve writing into the log file - no to be open and closed each time
			//	MPI_File file;
			//	int err = MPI_File_open(MPI_COMM_WORLD,
			//							"logfile.txt",
			//						MPI_MODE_CREATE| MPI_MODE_WRONLY|MPI_MODE_SEQUENTIAL,	MPI_INFO_NULL, &file);

			//	printf(" open file ierr=%d ", err);

			starttime = MPI_Wtime();

			compute_global_variables(training_graphs, t_params.num_training_graphs);
			int num_variables = 10;

			srand(r_params.current_id+r_params.num_procs); // we run each process with a different seed

			printf("Evolving. proc ID=%d ..\n", r_params.current_id);
			//start the evolution
			start_steady_state(t_params, training_graphs, num_variables,  r_params);

			//	MPI_File_close(&file);

			endtime = MPI_Wtime();
			computetime = endtime - starttime;



// write the execution times  on  the result file
			MPI_Reduce(&computetime, &compute_sum, 1, MPI_DOUBLE, MPI_SUM, 0,  MPI_COMM_WORLD);
			MPI_Reduce(&computetime, &compute_max, 1, MPI_DOUBLE, MPI_MAX, 0,  MPI_COMM_WORLD);
			if (r_params.current_id==0)
			{
				FILE* f = fopen(r_params.result_file,"a");
				fprintf(f,"___________________________________________________________________\n");
				fprintf(f, "number of processes = %d \n", r_params.num_procs);
				fprintf(f, "average computation time per process is %f seconds\n", compute_sum/r_params.num_procs);
				fprintf(f, "maximum  computation time of processes is %f seconds\n", compute_max);
				fprintf(f,"___________________________________________________________________\n");
				fprintf(f,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
				fclose(f);
			}

			print_num_migrations(r_params);

			MPI_Barrier(MPI_COMM_WORLD);
			if (r_params.current_id==0) {
				FILE *f = fopen(r_params.log_file, "a");
				fprintf(f, "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
				fclose(f);
			}

		}
		delete_training_graphs(training_graphs, t_params.num_training_graphs);
	}
	MPI_Finalize();

#endif

	return 0;
}
//--------------------------------------------------------------------
//--------------------------------------------------------------------